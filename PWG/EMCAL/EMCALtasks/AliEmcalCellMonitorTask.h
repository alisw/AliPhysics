#ifndef ALIEMCALCELLMONITOR_H_
#define ALIEMCALCELLMONITOR_H_
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"

class THistManager;
class AliEMCALGeometry;

/**
 * @class AliEmcalCellMonitorTask
 * @brief Simple monitoring task for cell related quantities
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @since May 19, 2016
 * @ingroup EMCALFWTASKS
 *
 * This task monitors simple cell quantities like
 * - Amplitude distribution
 * - Time distribution
 * - Count rate in col-row space
 * - Integrated amplitude in col-row space
 *
 * It can be added to the train using the add macro
 * ~~~
 * $ALICE_PHYSICS/PWG/EMCAL/AddEmcalCellMonitorTask.C
 * ~~~
 */
class AliEmcalCellMonitorTask : public AliAnalysisTaskSE {
public:

  /**
   * Dummy constructor, for ROOT I/O only
   */
  AliEmcalCellMonitorTask();

  /**
   * Default constructor, setting also the name and defining the output type
   * @param name Name of the task
   */
  AliEmcalCellMonitorTask(const char *name);

  /**
   * Destructor
   */
  virtual ~AliEmcalCellMonitorTask();

protected:

  /**
   * Create the output histograms
   *
   * For all supermodules the following histograms will be created:
   * - cellAmplitude with the cell amplitude distribution
   * - cellTime with the cell time distribution
   *
   * For each supermodule the followign histograms will be created:
   * - cellAmpSM with the integrated amplitude for cells within a supermodule in col and row
   * - cellCountSM with the count rate for cells within a supermodule in col and row
   */
  virtual void UserCreateOutputObjects();

  /**
   * Event loop.
   *
   * Running over all cells and filling the cell-related
   * quantities specified in the class documentation.
   *
   * @param[in]
   */
  virtual void UserExec(Option_t *);

private:
  THistManager                        *fHistManager;        //!<! Histogram handler
  AliEMCALGeometry                    *fGeometry;           //!<! EMCAL geometry

  AliEmcalCellMonitorTask(const AliEmcalCellMonitorTask &ref);
  AliEmcalCellMonitorTask &operator=(const AliEmcalCellMonitorTask &ref);

  /// \cond CLASSIMP
  ClassDef(AliEmcalCellMonitorTask, 1);
  /// \endcond
};

#endif /* ALIEMCALCELLMONITORTASK_H_ */
