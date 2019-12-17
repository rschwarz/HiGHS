/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HDualMulti.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "lp_data/HConst.h"
#include "simplex/HDual.h"
#include "simplex/HPrimal.h"
#include "simplex/SimplexTimer.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <set>
#include <stdexcept>

using std::cout;
using std::endl;

void HDual::iterateMulti() {
  slice_PRICE = 1;

  // Report candidate
  majorChooseRow();
  minorChooseRow();
  if (rowOut == -1) {
    invertHint = INVERT_HINT_POSSIBLY_OPTIMAL;
    return;
  }

  // Assign the slice_row_ep, skip if possible
  if (1.0 * multi_finish[multi_nFinish].row_ep->count / solver_num_row < 0.01)
    slice_PRICE = 0;

  if (slice_PRICE) {
#pragma omp parallel
#pragma omp single
    chooseColumnSlice(multi_finish[multi_nFinish].row_ep);
  } else {
    chooseColumn(multi_finish[multi_nFinish].row_ep);
  }
  // If we failed.
  if (invertHint) {
    majorUpdate();
    return;
  }

  minorUpdate();
  majorUpdate();
}

void HDual::majorChooseRow() {
  /**
   * 0. Initial check to see if we need to do it again
   */
  if (workHMO.simplex_info_.update_count == 0) multi_chooseAgain = 1;
  if (!multi_chooseAgain) return;
  multi_chooseAgain = 0;
  multi_iteration++;

  /**
   * Major loop:
   *     repeat 1-5, until we found a good sets of choices
   */
  int* choiceIndex = new int[multi_num];
  for (;;) {
    // 1. Multiple CHUZR
    int initialCount = 0;
    
    dualRHS.choose_multi_HGauto(&choiceIndex[0], &initialCount, multi_num);
    //        dualRHS.choose_multi_global(&choiceIndex[0], &initialCount,
    //        multi_num);
    if (initialCount == 0 && dualRHS.workCutoff == 0) {
      // OPTIMAL
      return;
    }

    // 2. Shrink the size by cutoff
    int choiceCount = 0;
    for (int i = 0; i < initialCount; i++) {
      int iRow = choiceIndex[i];
      if (dualRHS.work_infeasibility[iRow] / dualRHS.workEdWt[iRow] >=
          dualRHS.workCutoff) {
        choiceIndex[choiceCount++] = iRow;
      }
    }

    if (initialCount == 0 || choiceCount <= initialCount / 3) {
      // Need to do the list again
      dualRHS.create_infeasList(columnDensity);
      continue;
    }

    // 3. Store the choiceIndex to buffer
    for (int ich = 0; ich < multi_num; ich++) multi_choice[ich].rowOut = -1;
    for (int ich = 0; ich < choiceCount; ich++)
      multi_choice[ich].rowOut = choiceIndex[ich];

    // 4. Parallel BTRAN and compute weight
    majorChooseRowBtran();

    // 5. Update row densities
    for (int ich = 0; ich < multi_num; ich++) {
      if (multi_choice[ich].rowOut >= 0) {
        row_epDensity *= 0.95;
        row_epDensity += 0.05 * multi_choice[ich].row_ep.count / solver_num_row;
      }
    }

    if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
      // 6. Check updated and computed weight - just for dual steepest edge
      int countWrongEdWt = 0;
      for (int i = 0; i < multi_num; i++) {
	const int iRow = multi_choice[i].rowOut;
	if (iRow < 0) continue;
	double updated_edge_weight = dualRHS.workEdWt[iRow];
	computed_edge_weight = dualRHS.workEdWt[iRow] = multi_choice[i].infeasEdWt;
	//      if (updated_edge_weight < 0.25 * computed_edge_weight) {
	if (!acceptDualSteepestEdgeWeight(updated_edge_weight)) {
	  multi_choice[i].rowOut = -1;
	  countWrongEdWt++;
	}
      }
      if (countWrongEdWt <= choiceCount / 3) break;
    } else {
      // No checking required if not using dual steepest edge so break
      break;
    }
  }
  delete[] choiceIndex;

  // 6. Take other info associated with choices
  double pami_cutoff = 0.95;
  for (int i = 0; i < multi_num; i++) {
    const int iRow = multi_choice[i].rowOut;
    if (iRow < 0) continue;
    // Other info
    multi_choice[i].baseValue = baseValue[iRow];
    multi_choice[i].baseLower = baseLower[iRow];
    multi_choice[i].baseUpper = baseUpper[iRow];
    multi_choice[i].infeasValue = dualRHS.work_infeasibility[iRow];
    multi_choice[i].infeasEdWt = dualRHS.workEdWt[iRow];
    multi_choice[i].infeasLimit =
        dualRHS.work_infeasibility[iRow] / dualRHS.workEdWt[iRow];
    multi_choice[i].infeasLimit *= pami_cutoff;
  }

  // 6. Finish count
  multi_nFinish = 0;
}

void HDual::majorChooseRowBtran() {
  HighsTimer& timer = workHMO.timer_;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  timer.start(simplex_info.clock_[BtranClock]);

  // 4.1. Prepare BTRAN buffer
  int multi_ntasks = 0;
  int multi_iRow[HIGHS_THREAD_LIMIT];
  int multi_iwhich[HIGHS_THREAD_LIMIT];
  double multi_EdWt[HIGHS_THREAD_LIMIT];
  HVector_ptr multi_vector[HIGHS_THREAD_LIMIT];
  for (int ich = 0; ich < multi_num; ich++) {
    if (multi_choice[ich].rowOut >= 0) {
      multi_iRow[multi_ntasks] = multi_choice[ich].rowOut;
      multi_vector[multi_ntasks] = &multi_choice[ich].row_ep;
      multi_iwhich[multi_ntasks] = ich;
      multi_ntasks++;
    }
  }

  // 4.2 Perform BTRAN
#pragma omp parallel for schedule(static, 1)
  for (int i = 0; i < multi_ntasks; i++) {
    const int iRow = multi_iRow[i];
    HVector_ptr work_ep = multi_vector[i];
    work_ep->clear();
    work_ep->count = 1;
    work_ep->index[0] = iRow;
    work_ep->array[iRow] = 1;
    work_ep->packFlag = true;
    factor->btran(*work_ep, row_epDensity);
    if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
      // For Dual steepest edge we know the exact weight as the 2-norm of work_ep
      multi_EdWt[i] = work_ep->norm2();
    } else {
      // For Devex (and Dantzig) we take the updated edge weight
      multi_EdWt[i] = dualRHS.workEdWt[iRow];
    }
  }

  // 4.3 Put back edge weights: the edge weights for the chosen rows
  // are stored in multi_choice[*].infeasEdWt
  for (int i = 0; i < multi_ntasks; i++)
    multi_choice[multi_iwhich[i]].infeasEdWt = multi_EdWt[i];

  timer.stop(simplex_info.clock_[BtranClock]);
}

void HDual::minorChooseRow() {
  /**
   * 1. Find which to go out
   *        Because we had other checking code
   *        We know current best is OK to be used
   */
  multi_iChoice = -1;
  double bestMerit = 0;
  for (int ich = 0; ich < multi_num; ich++) {
    const int iRow = multi_choice[ich].rowOut;
    if (iRow < 0) continue;
    double infeasValue = multi_choice[ich].infeasValue;
    double infeasEdWt = multi_choice[ich].infeasEdWt;
    double infeasMerit = infeasValue / infeasEdWt;
    if (bestMerit < infeasMerit) {
      bestMerit = infeasMerit;
      multi_iChoice = ich;
    }
  }

  /**
   * 2. Obtain other info for current sub-optimization choice
   */
  rowOut = -1;
  if (multi_iChoice != -1) {
    MChoice* workChoice = &multi_choice[multi_iChoice];

    // Assign useful variables
    rowOut = workChoice->rowOut;
    columnOut = workHMO.simplex_basis_.basicIndex_[rowOut];
    double valueOut = workChoice->baseValue;
    double lowerOut = workChoice->baseLower;
    double upperOut = workChoice->baseUpper;
    deltaPrimal = valueOut - (valueOut < lowerOut ? lowerOut : upperOut);
    sourceOut = deltaPrimal < 0 ? -1 : 1;

    // Assign buffers
    MFinish* finish = &multi_finish[multi_nFinish];
    finish->rowOut = rowOut;
    finish->columnOut = columnOut;
    finish->row_ep = &workChoice->row_ep;
    finish->column = &workChoice->column;
    finish->columnBFRT = &workChoice->columnBFRT;
    // Save the edge weight - over-written later when using Devex
    finish->EdWt = workChoice->infeasEdWt;

    // Disable current row
    workChoice->rowOut = -1;
  }
}

void HDual::minorUpdate() {

  if (dual_edge_weight_mode == DualEdgeWeightMode::DEVEX && !new_devex_framework && rowOut>=0) {
    const double updated_edge_weight = dualRHS.workEdWt[rowOut];
    new_devex_framework = newDevexFramework(updated_edge_weight);
    minor_new_devex_framework = new_devex_framework;
  }

  // Minor update - store roll back data
  MFinish* finish = &multi_finish[multi_nFinish];
  finish->moveIn = workHMO.simplex_basis_.nonbasicMove_[columnIn];
  finish->shiftOut = workHMO.simplex_info_.workShift_[columnOut];
  finish->flipList.clear();
  for (int i = 0; i < dualRow.workCount; i++)
    finish->flipList.push_back(dualRow.workData[i].first);

  // Minor update - key parts
  minorUpdateDual();
  minorUpdatePrimal();
  minorUpdatePivots();
  minorUpdateRows();
  if (minor_new_devex_framework) {
    printf("Iter %7d (Major %7d): Minor new Devex framework\n",
	   workHMO.scaled_solution_params_.simplex_iteration_count,
	   multi_iteration);
    minorInitialiseDevexFramework();
  }
  multi_nFinish++;

  // Minor update - check for the next iteration
  int countRemain = 0;
  for (int i = 0; i < multi_num; i++) {
    int iRow = multi_choice[i].rowOut;
    if (iRow < 0) continue;
    double myInfeas = multi_choice[i].infeasValue;
    double myWeight = multi_choice[i].infeasEdWt;
    countRemain += (myInfeas / myWeight > multi_choice[i].infeasLimit);
  }
  if (countRemain == 0) multi_chooseAgain = 1;
  //    if (multi_nFinish + 1 == multi_num)
  //        multi_chooseAgain = 1;

}

void HDual::minorUpdateDual() {
  /**
   * 1. Update the dual solution
   *    XXX Data parallel (depends on the ap partition before)
   */
  if (thetaDual == 0) {
    shift_cost(workHMO, columnIn, -workDual[columnIn]);  // model->shiftCost(columnIn,
                                                         // -workDual[columnIn]);
  } else {
    dualRow.update_dual(thetaDual);//, columnOut);
    if (slice_PRICE) {
      for (int i = 0; i < slice_num; i++)
        slice_dualRow[i].update_dual(thetaDual);//, columnOut);
    }
  }
  workDual[columnIn] = 0;
  workDual[columnOut] = -thetaDual;
  shift_back(workHMO, columnOut);  // model->shiftBack(columnOut);

  /**
   * 2. Apply global bound flip
   */
  dualRow.update_flip(multi_finish[multi_nFinish].columnBFRT);

  /**
   * 3. Apply local bound flips
   */
  for (int ich = 0; ich < multi_num; ich++) {
    if (ich == multi_iChoice || multi_choice[ich].rowOut >= 0) {
      HVector* this_ep = &multi_choice[ich].row_ep;
      for (int i = 0; i < dualRow.workCount; i++) {
        double dot = matrix->compute_dot(*this_ep, dualRow.workData[i].first);
        multi_choice[ich].baseValue -= dualRow.workData[i].second * dot;
      }
    }
  }
}

void HDual::minorUpdatePrimal() {
  MChoice* choice = &multi_choice[multi_iChoice];
  MFinish* finish = &multi_finish[multi_nFinish];
  double valueOut = choice->baseValue;
  double lowerOut = choice->baseLower;
  double upperOut = choice->baseUpper;
  if (deltaPrimal < 0) {
    thetaPrimal = (valueOut - lowerOut) / alphaRow;
    finish->basicBound = lowerOut;
  }
  if (deltaPrimal > 0) {
    thetaPrimal = (valueOut - upperOut) / alphaRow;
    finish->basicBound = upperOut;
  }
  finish->thetaPrimal = thetaPrimal;

  // With Devex, the edge weight for the pivotal row is computed in
  // chooseColumn (if not computed in chooseColumnSlice), so can't be stored in the PAMI data structure.
  //
  // Restore the weight to multi_finish so that it's picked up in
  // minorUpdatePrimal();
  if (dual_edge_weight_mode == DualEdgeWeightMode::DEVEX) finish->EdWt = computed_edge_weight;



  /**
   * 5. Update the other primal value
   *    By the pivot (thetaPrimal)
   */
  if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
    const double dl_EdWt = fabs(choice->infeasEdWt - finish->EdWt);
    if (dl_EdWt > 0) {
      printf("minorUpdatePrimal: Delta Edge Weight = |%11.4g-%11.4g| = %11.4g\n",
	     choice->infeasEdWt, finish->EdWt, dl_EdWt);
    }
  }

  if (dual_edge_weight_mode == DualEdgeWeightMode::DEVEX) {
    // Transform the edge weight of the pivotal row according to the
    // simplex update
    double updated_pivotal_edge_weight = finish->EdWt / (alphaRow * alphaRow);
    updated_pivotal_edge_weight = max(1.0, updated_pivotal_edge_weight);
    // Store the Devex weight of the leaving row now - OK since it's
    // stored in finish->EdWt and the updated weights are stored in
    // multi_choice[*].infeasEdWt
    finish->EdWt = updated_pivotal_edge_weight;
  }
  for (int ich = 0; ich < multi_num; ich++) {
    if (multi_choice[ich].rowOut >= 0) {
      HVector* this_ep = &multi_choice[ich].row_ep;
      double dot = matrix->compute_dot(*this_ep, columnIn);
      multi_choice[ich].baseValue -= thetaPrimal * dot;
      double value = multi_choice[ich].baseValue;
      double lower = multi_choice[ich].baseLower;
      double upper = multi_choice[ich].baseUpper;
      double infeas = 0;
      if (value < lower - Tp) infeas = value - lower;
      if (value > upper + Tp) infeas = value - upper;
      infeas *= infeas;
      multi_choice[ich].infeasValue = infeas;
      if (dual_edge_weight_mode == DualEdgeWeightMode::DEVEX) {
	const double updated_pivotal_edge_weight = finish->EdWt;
	double aa_iRow = dot;
	multi_choice[ich].infeasEdWt = max(multi_choice[ich].infeasEdWt, updated_pivotal_edge_weight * aa_iRow * aa_iRow);
      }
    }
  }
}
void HDual::minorUpdatePivots() {
  MFinish* finish = &multi_finish[multi_nFinish];
  update_pivots(workHMO, columnIn, rowOut, sourceOut);
  if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
    // Transform the edge weight of the pivotal row according to the
    // simplex update
    finish->EdWt /= (alphaRow * alphaRow);
  }
  finish->basicValue = workHMO.simplex_info_.workValue_[columnIn] + thetaPrimal;
  update_matrix(workHMO, columnIn,
                columnOut);  // model->updateMatrix(columnIn, columnOut);
  finish->columnIn = columnIn;
  finish->alphaRow = alphaRow;
  // Move this to Simplex class once it's created
  // simplex_method.record_pivots(columnIn, columnOut, alphaRow);
  workHMO.scaled_solution_params_.simplex_iteration_count++;
}

void HDual::minorUpdateRows() {
  HighsTimer& timer = workHMO.timer_;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  timer.start(simplex_info.clock_[UpdateRowClock]);
  const HVector* Row = multi_finish[multi_nFinish].row_ep;
  int updateRows_inDense =
      (Row->count < 0) || (Row->count > 0.1 * solver_num_row);
  if (updateRows_inDense) {
    int multi_nTasks = 0;
    int multi_iwhich[HIGHS_THREAD_LIMIT];
    double multi_xpivot[HIGHS_THREAD_LIMIT];
    HVector_ptr multi_vector[HIGHS_THREAD_LIMIT];

    /*
     * Dense mode
     *  1. Find which ones to do and the pivotX
     *  2. Do all of them in task parallel
     */

    // Collect tasks
    for (int ich = 0; ich < multi_num; ich++) {
      if (multi_choice[ich].rowOut >= 0) {
        HVector* next_ep = &multi_choice[ich].row_ep;
        double pivotX = matrix->compute_dot(*next_ep, columnIn);
        if (fabs(pivotX) < HIGHS_CONST_TINY) continue;
        multi_vector[multi_nTasks] = next_ep;
        multi_xpivot[multi_nTasks] = -pivotX / alphaRow;
        multi_iwhich[multi_nTasks] = ich;
        multi_nTasks++;
      }
    }

    // Perform tasks
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < multi_nTasks; i++) {
      HVector_ptr nextEp = multi_vector[i];
      const double xpivot = multi_xpivot[i];
      nextEp->saxpy(xpivot, Row);
      nextEp->tight();
      if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
	multi_xpivot[i] = nextEp->norm2();
      }
    }

    // Put weight back
    if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
      for (int i = 0; i < multi_nTasks; i++)
	multi_choice[multi_iwhich[i]].infeasEdWt = multi_xpivot[i];
    }
  } else {
    // Sparse mode: just do it sequentially
    for (int ich = 0; ich < multi_num; ich++) {
      if (multi_choice[ich].rowOut >= 0) {
        HVector* next_ep = &multi_choice[ich].row_ep;
        double pivotX = matrix->compute_dot(*next_ep, columnIn);
        if (fabs(pivotX) < HIGHS_CONST_TINY) continue;
        next_ep->saxpy(-pivotX / alphaRow, Row);
        next_ep->tight();
	if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
	  multi_choice[ich].infeasEdWt = next_ep->norm2();
	}
      }
    }
  }
  timer.stop(simplex_info.clock_[UpdateRowClock]);
}

void HDual::minorInitialiseDevexFramework() {
  // Set the local Devex weights to 1
  for (int i = 0; i < multi_num; i++) {
    multi_choice[i].infeasEdWt = 1.0;
  }
  minor_new_devex_framework = false;
}

void HDual::majorUpdate() {
  /**
   * 0. See if it's ready to perform a major update
   */
  if (invertHint) multi_chooseAgain = 1;
  if (!multi_chooseAgain) return;

  // Major update - FTRANs
  majorUpdateFtranPrepare();
  majorUpdateFtranParallel();
  majorUpdateFtranFinal();

  // Major update - check for roll back
  for (int iFn = 0; iFn < multi_nFinish; iFn++) {
    MFinish* iFinish = &multi_finish[iFn];
    HVector* iColumn = iFinish->column;
    int iRowOut = iFinish->rowOut;
    double alphaC = fabs(iColumn->array[iRowOut]);
    double alphaR = fabs(iFinish->alphaRow);
    double compare = min(alphaC, alphaR);
    double alphaDiff = fabs(alphaC - alphaR);
    // int startUpdate = workHMO.simplex_info_.update_count - multi_nFinish;
    if (alphaDiff / compare > 1e-8 && workHMO.simplex_info_.update_count > 0) {
      cout << "REPORT " << workHMO.simplex_lp_.model_name_
           << " NEED-ROLL-BACK   ";
      cout << workHMO.scaled_solution_params_.simplex_iteration_count << " alpha = " << alphaC
           << " alphaR = " << alphaR << " diff = " << alphaDiff / compare
           << "  multi_nFinish = " << multi_nFinish << endl;
      invertHint = INVERT_HINT_POSSIBLY_SINGULAR_BASIS;
      // if (startUpdate > 0) {
      majorRollback();
      return;
      // }
    }
  }

  // Major update - primal and factor
  majorUpdatePrimal();
  majorUpdateFactor();
  if (new_devex_framework) {
    //    printf("Iter %7d: New Devex framework\n", workHMO.scaled_solution_params_.simplex_iteration_count);
    const bool parallel = true;
    initialiseDevexFramework(parallel);
  }
}

void HDual::majorUpdateFtranPrepare() {
  // Prepare FTRAN BFRT buffer
  columnBFRT.clear();
  for (int iFn = 0; iFn < multi_nFinish; iFn++) {
    MFinish* finish = &multi_finish[iFn];
    HVector* Vec = finish->columnBFRT;
    matrix->collect_aj(*Vec, finish->columnIn, finish->thetaPrimal);

    // Update this buffer by previous Row_ep
    for (int jFn = iFn - 1; jFn >= 0; jFn--) {
      MFinish* jFinish = &multi_finish[jFn];
      double* jRow_epArray = &jFinish->row_ep->array[0];
      double pivotX = 0;
      for (int k = 0; k < Vec->count; k++) {
        int iRow = Vec->index[k];
        pivotX += Vec->array[iRow] * jRow_epArray[iRow];
      }
      if (fabs(pivotX) > HIGHS_CONST_TINY) {
        pivotX /= jFinish->alphaRow;
        matrix->collect_aj(*Vec, jFinish->columnIn, -pivotX);
        matrix->collect_aj(*Vec, jFinish->columnOut, pivotX);
      }
    }
    columnBFRT.saxpy(1, Vec);
  }

  // Prepare regular FTRAN buffer
  for (int iFn = 0; iFn < multi_nFinish; iFn++) {
    MFinish* iFinish = &multi_finish[iFn];
    HVector* iColumn = iFinish->column;
    iColumn->clear();
    iColumn->packFlag = true;
    matrix->collect_aj(*iColumn, iFinish->columnIn, 1);
  }
}

void HDual::majorUpdateFtranParallel() {
  HighsTimer& timer = workHMO.timer_;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  timer.start(simplex_info.clock_[FtranMixParClock]);

  // Prepare buffers
  int multi_ntasks = 0;
  double multi_density[HIGHS_THREAD_LIMIT * 2 + 1];
  HVector_ptr multi_vector[HIGHS_THREAD_LIMIT * 2 + 1];
  // BFRT first
  multi_density[multi_ntasks] = columnDensity;
  multi_vector[multi_ntasks] = &columnBFRT;
  multi_ntasks++;
  if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
    // Then DSE
    for (int iFn = 0; iFn < multi_nFinish; iFn++) {
      multi_density[multi_ntasks] = rowdseDensity;
      multi_vector[multi_ntasks] = multi_finish[iFn].row_ep;
      multi_ntasks++;
    }
  }
  // Then Column
  for (int iFn = 0; iFn < multi_nFinish; iFn++) {
    multi_density[multi_ntasks] = columnDensity;
    multi_vector[multi_ntasks] = multi_finish[iFn].column;
    multi_ntasks++;
  }

  // Perform FTRAN
#pragma omp parallel for schedule(dynamic, 1)
  for (int i = 0; i < multi_ntasks; i++) {
    HVector_ptr rhs = multi_vector[i];
    double density = multi_density[i];
    factor->ftran(*rhs, density);
  }

  // Update ticks
  for (int iFn = 0; iFn < multi_nFinish; iFn++) {
    MFinish* finish = &multi_finish[iFn];
    HVector* Col = finish->column;
    HVector* Row = finish->row_ep;
    total_FT_inc_TICK += Col->syntheticTick;  // Was .pseudoTick
    total_FT_inc_TICK += Row->syntheticTick;  // Was .pseudoTick
  }

  // Update rates
  for (int iFn = 0; iFn < multi_nFinish; iFn++) {
    MFinish* finish = &multi_finish[iFn];
    HVector* Col = finish->column;
    HVector* Row = finish->row_ep;
    columnDensity = 0.95 * columnDensity + 0.05 * Col->count / solver_num_row;
    rowdseDensity = 0.95 * rowdseDensity + 0.05 * Row->count / solver_num_row;
  }
  timer.stop(simplex_info.clock_[FtranMixParClock]);
}

void HDual::majorUpdateFtranFinal() {
  HighsTimer& timer = workHMO.timer_;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  timer.start(simplex_info.clock_[FtranMixFinalClock]);
  int updateFTRAN_inDense = dualRHS.workCount < 0;
  if (updateFTRAN_inDense) {
    for (int iFn = 0; iFn < multi_nFinish; iFn++) {
      multi_finish[iFn].column->count = -1;
      multi_finish[iFn].row_ep->count = -1;
      double* myCol = &multi_finish[iFn].column->array[0];
      double* myRow = &multi_finish[iFn].row_ep->array[0];
      for (int jFn = 0; jFn < iFn; jFn++) {
        int pivotRow = multi_finish[jFn].rowOut;
        const double pivotAlpha = multi_finish[jFn].alphaRow;
        const double* pivotArray = &multi_finish[jFn].column->array[0];
        double pivotX1 = myCol[pivotRow];
        double pivotX2 = myRow[pivotRow];

        // The FTRAN regular buffer
        if (fabs(pivotX1) > HIGHS_CONST_TINY) {
          const double pivot = pivotX1 / pivotAlpha;
#pragma omp parallel for
          for (int i = 0; i < solver_num_row; i++)
            myCol[i] -= pivot * pivotArray[i];
          myCol[pivotRow] = pivot;
        }
        // The FTRAN-DSE buffer
        if (fabs(pivotX2) > HIGHS_CONST_TINY) {
          const double pivot = pivotX2 / pivotAlpha;
#pragma omp parallel for
          for (int i = 0; i < solver_num_row; i++)
            myRow[i] -= pivot * pivotArray[i];
          myRow[pivotRow] = pivot;
        }
      }
    }
  } else {
    for (int iFn = 0; iFn < multi_nFinish; iFn++) {
      MFinish* finish = &multi_finish[iFn];
      HVector* Col = finish->column;
      HVector* Row = finish->row_ep;
      for (int jFn = 0; jFn < iFn; jFn++) {
        MFinish* jFinish = &multi_finish[jFn];
        int pivotRow = jFinish->rowOut;
        double pivotX1 = Col->array[pivotRow];
        // The FTRAN regular buffer
        if (fabs(pivotX1) > HIGHS_CONST_TINY) {
          pivotX1 /= jFinish->alphaRow;
          Col->saxpy(-pivotX1, jFinish->column);
          Col->array[pivotRow] = pivotX1;
        }
        // The FTRAN-DSE buffer
        double pivotX2 = Row->array[pivotRow];
        if (fabs(pivotX2) > HIGHS_CONST_TINY) {
          pivotX2 /= jFinish->alphaRow;
          Row->saxpy(-pivotX2, jFinish->column);
          Row->array[pivotRow] = pivotX2;
        }
      }
    }
  }
  timer.stop(simplex_info.clock_[FtranMixFinalClock]);
}

void HDual::majorUpdatePrimal() {
  const bool updatePrimal_inDense = dualRHS.workCount < 0;
  if (updatePrimal_inDense) {
    // Update the RHS in dense
    const double* mixArray = &columnBFRT.array[0];
    double* local_work_infeasibility = &dualRHS.work_infeasibility[0];
#pragma omp parallel for schedule(static)
    for (int iRow = 0; iRow < solver_num_row; iRow++) {
      baseValue[iRow] -= mixArray[iRow];
      const double value = baseValue[iRow];
      const double less = baseLower[iRow] - value;
      const double more = value - baseUpper[iRow];
      double infeas = less > Tp ? less : (more > Tp ? more : 0);
      //      local_work_infeasibility[iRow] = infeas * infeas;
      if (workHMO.simplex_info_.store_squared_primal_infeasibility) 
	local_work_infeasibility[iRow] = infeas * infeas;
      else
	local_work_infeasibility[iRow] = fabs(infeas);
    }

    if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE ||
	(dual_edge_weight_mode == DualEdgeWeightMode::DEVEX && !new_devex_framework)) {
      // Update any edge weights (except weights for pivotal rows) in dense
      for (int iFn = 0; iFn < multi_nFinish; iFn++) {
	// multi_finish[iFn].EdWt has already been transformed to correspond to the new basis
	const double updated_pivotal_edge_weight = multi_finish[iFn].EdWt;
	const double* colArray = &multi_finish[iFn].column->array[0];
	double* EdWt = &dualRHS.workEdWt[0];
	if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
	  // Update steepest edge weights 
	  const double* dseArray = &multi_finish[iFn].row_ep->array[0];
	  const double Kai = -2 / multi_finish[iFn].alphaRow;
#pragma omp parallel for schedule(static)
	  for (int iRow = 0; iRow < solver_num_row; iRow++) {
	    const double aa_iRow = colArray[iRow];
	    EdWt[iRow] += aa_iRow * (updated_pivotal_edge_weight * aa_iRow + Kai * dseArray[iRow]);
	    if (EdWt[iRow] < 1e-4) EdWt[iRow] = 1e-4;
	  }
	} else {
	  for (int iRow = 0; iRow < solver_num_row; iRow++) {
	    const double aa_iRow = colArray[iRow];
	    EdWt[iRow] = max(EdWt[iRow], updated_pivotal_edge_weight * aa_iRow * aa_iRow);
	  }
	  num_devex_iterations++;
	}
      }
    }
  } else {
    // Update primal and pivots
    dualRHS.update_primal(&columnBFRT, 1);
    dualRHS.update_infeasList(&columnBFRT);

    // Update any edge weights (except weights for pivotal rows) and infeasList
    for (int iFn = 0; iFn < multi_nFinish; iFn++) {
      MFinish* finish = &multi_finish[iFn];
      HVector* Col = finish->column;
      const double updated_pivotal_edge_weight = finish->EdWt;
      if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
	// Update steepest edge weights 
	HVector* Row = finish->row_ep;
	double Kai = -2 / finish->alphaRow;
	dualRHS.updateWeightDualSteepestEdge(Col, updated_pivotal_edge_weight, Kai, &Row->array[0]);
      } else if (dual_edge_weight_mode == DualEdgeWeightMode::DEVEX && !new_devex_framework) {
	// Update rest of weights
	dualRHS.updateWeightDevex(Col, updated_pivotal_edge_weight);
	num_devex_iterations++;
      }
      dualRHS.update_infeasList(Col);
    }
  }

  // Update primal value for the pivots
  for (int iFn = 0; iFn < multi_nFinish; iFn++) {
    MFinish* finish = &multi_finish[iFn];
    int iRow = finish->rowOut;
    double value = baseValue[iRow] - finish->basicBound + finish->basicValue;
    dualRHS.update_pivots(iRow, value);
  }

  if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
    // Update weight value for the pivots
    for (int iFn = 0; iFn < multi_nFinish; iFn++) {
      const int iRow = multi_finish[iFn].rowOut;
      const double updated_pivotal_edge_weight = multi_finish[iFn].EdWt;
      const double* colArray = &multi_finish[iFn].column->array[0];
      const double* dseArray = &multi_finish[iFn].row_ep->array[0];
      double Kai = -2 / multi_finish[iFn].alphaRow;
      for (int jFn = 0; jFn < iFn; jFn++) {
	int jRow = multi_finish[jFn].rowOut;
	double value = colArray[jRow];
	double EdWt = dualRHS.workEdWt[jRow];
	EdWt += value * (updated_pivotal_edge_weight * value + Kai * dseArray[jRow]);
	if (EdWt < min_dual_steepest_edge_weight) EdWt = min_dual_steepest_edge_weight;
	dualRHS.workEdWt[jRow] = EdWt;
      }
      dualRHS.workEdWt[iRow] = updated_pivotal_edge_weight;
    }
  } else if (dual_edge_weight_mode == DualEdgeWeightMode::DEVEX && !new_devex_framework) {
    for (int iFn = 0; iFn < multi_nFinish; iFn++) {
      const int iRow = multi_finish[iFn].rowOut;
      const double updated_pivotal_edge_weight = multi_finish[iFn].EdWt;
      const double* colArray = &multi_finish[iFn].column->array[0];
      // Pivotal row is for the current basis: weights are required for
      // the next basis so have to divide the current (exact) weight by
      // the pivotal value
      for (int jFn = 0; jFn < iFn; jFn++) {
	int jRow = multi_finish[jFn].rowOut;
	const double aa_iRow = colArray[iRow];
	double EdWt = dualRHS.workEdWt[jRow];
	EdWt = max(EdWt, updated_pivotal_edge_weight * aa_iRow * aa_iRow);
	//	  dualRHS.workEdWt[jRow] = EdWt;
      }
      dualRHS.workEdWt[iRow] = updated_pivotal_edge_weight;
      num_devex_iterations++;
    }
  }
  if (dual_edge_weight_mode == DualEdgeWeightMode::DEVEX) new_devex_framework = true;
  checkNonUnitWeightError("999");
}

void HDual::majorUpdateFactor() {
  /**
   * 9. Update the factor by CFT
   */
  int* iRows = new int[multi_nFinish];
  for (int iCh = 0; iCh < multi_nFinish - 1; iCh++) {
    multi_finish[iCh].row_ep->next = multi_finish[iCh + 1].row_ep;
    multi_finish[iCh].column->next = multi_finish[iCh + 1].column;
    iRows[iCh] = multi_finish[iCh].rowOut;
  }
  iRows[multi_nFinish - 1] = multi_finish[multi_nFinish - 1].rowOut;
  if (multi_nFinish > 0)
    update_factor(workHMO, multi_finish[0].column, multi_finish[0].row_ep, iRows, &invertHint);
  if (total_FT_inc_TICK > total_INVERT_TICK * 1.5 &&
      workHMO.simplex_info_.update_count > 200)
    invertHint = INVERT_HINT_SYNTHETIC_CLOCK_SAYS_INVERT;
  delete[] iRows;
}

void HDual::majorRollback() {
  for (int iFn = multi_nFinish - 1; iFn >= 0; iFn--) {
    MFinish* finish = &multi_finish[iFn];

    // 1. Roll back pivot
    workHMO.simplex_basis_.nonbasicMove_[finish->columnIn] = finish->moveIn;
    workHMO.simplex_basis_.nonbasicFlag_[finish->columnIn] = 1;
    workHMO.simplex_basis_.nonbasicMove_[finish->columnOut] = 0;
    workHMO.simplex_basis_.nonbasicFlag_[finish->columnOut] = 0;
    workHMO.simplex_basis_.basicIndex_[finish->rowOut] = finish->columnOut;

    // 2. Roll back matrix
    update_matrix(workHMO, finish->columnOut, finish->columnIn);

    // 3. Roll back flips
    for (unsigned i = 0; i < finish->flipList.size(); i++) {
      flip_bound(workHMO, finish->flipList[i]);
    }

    // 4. Roll back cost
    workHMO.simplex_info_.workShift_[finish->columnIn] = 0;
    workHMO.simplex_info_.workShift_[finish->columnOut] = finish->shiftOut;

    // 5. The iteration count
    workHMO.scaled_solution_params_.simplex_iteration_count--;
  }
}

bool HDual::checkNonUnitWeightError(std::string message) {
  bool error_found = false;
  if (dual_edge_weight_mode == DualEdgeWeightMode::DANTZIG) {
    double unit_wt_error = 0;
    for (int iRow = 0; iRow < solver_num_row; iRow++) {
      unit_wt_error += fabs(dualRHS.workEdWt[iRow]-1.0);
    }
    error_found = unit_wt_error>1e-4;
    if (error_found) printf("Non-unit Edge weight error of %g: %s\n", unit_wt_error, message.c_str());
  }
  return error_found;
}
