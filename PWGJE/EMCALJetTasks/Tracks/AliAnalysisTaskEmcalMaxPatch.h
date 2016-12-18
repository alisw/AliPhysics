#ifndef ALIANALYSISTASKEMCALMAXPATCH_H
#define ALIANALYSISTASKEMCALMAXPATCH_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcal.h"
#include <TString.h>

class THistManager;

namespace EMCalTriggerPtAnalysis {

class AliEMCalTriggerWeightHandler;

/**
 * @class AliAnalysisTaskEmcalMaxPatch
 * @brief Simple task monitoring the energy spectrum of the maximum patch in the event.
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @since Aug 30th, 2016
 *
 * With the help of the maximum trigger patch energy spectrum trigger rejection factors
 * can be calculated as function of the trigger threshold.
 */
class AliAnalysisTaskEmcalMaxPatch: public AliAnalysisTaskEmcal {
public:

  /**
   * Dummy (I/O) constructor
   */
  AliAnalysisTaskEmcalMaxPatch();

  /**
   * Main constructor, intializing the task
   * @param[in] name Name of the task
   */
  AliAnalysisTaskEmcalMaxPatch(const char *name);

  /**
   * Destructor, cleaning up
   */
  virtual ~AliAnalysisTaskEmcalMaxPatch();

  /**
   * Specify trigger to be selected. By default INT7 is used
   * as trigger.
   * @param[in] trigger Trigger (after physics selection) to be selected
   */
  void SetTriggerSelection(ULong_t trigger) { fSelectTrigger = trigger; };

  /**
   * Apply additional event selection using the trigger string
   * (for EMCAL level1 triggers)
   * @param[in] pattern Trigger string pattern to be requested
   */
  void SetTriggerPattern(const TString &pattern) { fTriggerPattern = pattern; }

protected:

  /**
   * Creating output histograms. Also initializing the
   * analysis utils in case they are not provided from
   * outside.
   */
  virtual void               UserCreateOutputObjects();

  /**
   * Perform event selection. In this task the standard
   * ALICE event selection for pPb events is used
   * @return True if the event is selected, false otherwise
   */
  virtual Bool_t             IsEventSelected();

  /**
   * Select the highest EGA and EJE patch per event and
   * fill the histogram for the corresponding spectrum with
   * the energy measurement.
   * @return Always true
   */
  virtual Bool_t             Run();

private:
  const AliEMCalTriggerWeightHandler    *fWeightHandler;    ///< Weight handler (optional)
  THistManager                          *fHistos;           //!<! Local Histogram handler. Histograms are stored in the AliEmcalList later
  ULong_t                               fSelectTrigger;     ///< Trigger bit selection
  TString                               fTriggerPattern;    ///< Trigger pattern string (i.e. EG1)

  AliAnalysisTaskEmcalMaxPatch(const AliAnalysisTaskEmcalMaxPatch &);
  AliAnalysisTaskEmcalMaxPatch &operator=(const AliAnalysisTaskEmcalMaxPatch &);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalMaxPatch, 1);
  /// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIANALYSISTASKEMCALMAXPATCH_H */
