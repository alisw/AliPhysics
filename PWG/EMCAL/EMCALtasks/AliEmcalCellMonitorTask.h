#ifndef ALIEMCALCELLMONITOR_H_
#define ALIEMCALCELLMONITOR_H_
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"
#include <TCustomBinning.h>
#include <TString.h>
#include <vector>

class TArrayD;
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

  /**
   * Set the minimum amplitude requested for a cell in order to use the
   * cell time in the cell time plot
   * @param minAmp Min. cell amplitude
   */
  void SetMinCellAmplitudeForCellTime(Double_t minAmp) { fMinCellAmplitude = minAmp; }

  /**
   * Define trigger selection. It can consist of trigger bits and a trigger string
   * @param triggerbits Trigger bit selection
   * @param triggerstring Trigger string (i.e. EG1, EG2, DG1, ...)
   */
  void SetRequestTrigger(ULong_t triggerbits, TString triggerstring = "") {
    fRequestTrigger = triggerbits;
    fTriggerString = triggerstring;
  }

  /**
   * Set number of Cells in order to match run2 geometry
   */
  void SetRun2() { fNumberOfCells = 17664; }

  /**
   * Mark cell as bad (i.e) exclude from running
   * @param cellId Cell to be masked
   */
  void SetBadCell(Int_t cellId);

protected:

  /**
   * @class AliEmcalCellMonitorAmplitudeBinning
   * @brief Defining binning in amplitude direction
   *
   * Create binning for the amplitude:
   * 0 - 1 GeV: 100 MeV bins
   * 1 - 10 GeV: 500 MeV bins
   * 10 - 20 GeV: 1 GeV bins
   * 20 - 50 GeV: 2 GeV bins
   * 50 - 100 GeV: 5 GeV bins
   * 100 - 200 GeV: 10 GeV bins
   */
  class AliEmcalCellMonitorAmplitudeBinning : public TCustomBinning {
  public:

    /**
     * Constructor, defining minimum and ranges
     */
    AliEmcalCellMonitorAmplitudeBinning();

    /**
     * Destructor
     */
    virtual ~AliEmcalCellMonitorAmplitudeBinning() {}
  };

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

  /**
   * Check whether cell with a given ID is masked by the user
   * @param cellId Cell ID to be checked
   * @return True if the cell is masked, false otherwise
   */
  bool IsCellMasked(Int_t cellId) const;

private:
  THistManager                        *fHistManager;        //!<! Histogram handler
  AliEMCALGeometry                    *fGeometry;           //!<! EMCAL geometry

  Double_t                            fMinCellAmplitude;    ///< Min. cell amplitude requested for cell time and frequency
  ULong_t                             fRequestTrigger;      ///< Trigger selection
  TString                             fTriggerString;       ///< Trigger string in addition to trigger selection
  Int_t                               fNumberOfCells;       ///< Number of cells

  std::vector<Int_t>                  fMaskedCells;         ///< Vector of masked cells

  AliEmcalCellMonitorTask(const AliEmcalCellMonitorTask &ref);
  AliEmcalCellMonitorTask &operator=(const AliEmcalCellMonitorTask &ref);

  /// \cond CLASSIMP
  ClassDef(AliEmcalCellMonitorTask, 1);
  /// \endcond
};

#endif /* ALIEMCALCELLMONITORTASK_H_ */
