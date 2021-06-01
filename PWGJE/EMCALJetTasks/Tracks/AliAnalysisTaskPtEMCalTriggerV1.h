/**
 * \file AliAnalysisTaskPtEMCalTriggerV1.h
 * \brief Declaration of the re-structured analysis task of high-\f$ p_{t} \f$ tracks
 * in triggered events. Task only behaves as steering task for analysis components.
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Dec 12, 2014
 */
#ifndef ALIANALYSISTASKPTEMCALTRIGGERV1_H
#define ALIANALYSISTASKPTEMCALTRIGGERV1_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcalJet.h"
#include "AliEMCalTriggerAnaClassManager.h"
#include <TObjArray.h>
#include <TString.h>

class TArrayD;

namespace PWGJE {

namespace EMCALJetTasks {

class AliEMCalTriggerAnaTriggerClass;
class AliEMCalTriggerAnaTriggerDecisionConfig;
class AliEMCalTriggerBinningComponent;
class AliEMCalTriggerEventData;
class AliEMCalTriggerTaskGroup;

/**
 * \class AliAnalysisTaskPtEMCalTriggerV1
 * \brief Re-structured analysis task of high-\f$ p_{t} \f$ tracks in triggered events.
 *
 * Re-structured analysis task of the pt analysis on EMCal-triggered events:
 * Analysis steps are moved to analysis components, which are grouped by a common
 * event selection. The analysis task steers the event builder, runs each group,
 * and collects the output of all groups.
 */
class AliAnalysisTaskPtEMCalTriggerV1: public AliAnalysisTaskEmcalJet {
public:
  AliAnalysisTaskPtEMCalTriggerV1();
  AliAnalysisTaskPtEMCalTriggerV1(const char *name);
  virtual ~AliAnalysisTaskPtEMCalTriggerV1();

  void AddAnalysisGroup(AliEMCalTriggerTaskGroup *taskGroup);
  virtual void UserCreateOutputObjects();
  virtual Bool_t Run();

  void SetBinning(const char *dimname, int nbins, double *binning);
  void SetBinning(const char *dimname, const TArrayD &binning);

  /**
   * Add trigger class to the list of trigger classes
   * \param triggerclass Trigger class to be added
   */
  void AddTriggerClass(AliEMCalTriggerAnaTriggerClass * triggerclass) { fTriggerClassManager->AddTriggerClass(triggerclass); }

  /**
   * Set the name of the jet container for generator level jets
   * \param name
   */
  void SetMCJetContainerName(const char *name)          { fMCJetContainer = name; }

  /**
   * Set the name of the jet container for jets in data
   * \param name
   */
  void SetDataJetContainerName(const char *name)        { fDataJetContainer = name; }

  /**
   * Set trigger selection into debug mode
   * \param doDebug If true we run the trigger selection in debug mode.
   */
  void SetTriggerDebug(Bool_t doDebug = kTRUE) { fDoTriggerDebug = doDebug; }

  /**
   * Set configuration for the trigger decision generation
   * \param config configuration of the trigger selection
   */
  void SetTriggerDecisionConfig(AliEMCalTriggerAnaTriggerDecisionConfig *config) { fTriggerDecisionConfig = config; }

protected:
  AliEMCalTriggerEventData *BuildEvent();
  void FixTrackInputEvent(AliVTrack *trk);

  TObjArray                           *fTaskGroups;                       ///< grouped analysis components
  AliEMCalTriggerBinningComponent     *fBinning;                          ///< Global binning component
  AliEMCalTriggerAnaTriggerDecisionConfig *fTriggerDecisionConfig;		    ///< Configuration for the trigger decision handling
  AliEMCalTriggerAnaClassManager      *fTriggerClassManager;              ///< Manager for trigger classes
  TString                              fMCJetContainer;                   ///< Name of the Monte-Carlo jet container
  TString                              fDataJetContainer;                 ///< Data jet container name
  Bool_t                               fSwapTriggerThresholds;            ///< Swap thresholds of the low and high threshold trigger
  Bool_t                               fDoTriggerDebug;                   ///< Debug trigger decision creator

private:
  AliAnalysisTaskPtEMCalTriggerV1(const AliAnalysisTaskPtEMCalTriggerV1 &);
  AliAnalysisTaskPtEMCalTriggerV1 &operator=(const AliAnalysisTaskPtEMCalTriggerV1 &);

  ClassDef(AliAnalysisTaskPtEMCalTriggerV1, 1);
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIANALYSISTASKPTEMCALTRIGGERV1_H */
