#ifndef ALIANALYSISTASKEMCALTRIGGERMULTIPLICITY_H
#define ALIANALYSISTASKEMCALTRIGGERMULTIPLICITY_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcalTriggerBase.h"

class AliEmcalTrackSelection;

namespace PWGJE {
  
namespace EMCALJetTasks {

/**
 * @class AliAnalysisTaskEmcalTriggerMultiplicity
 * @brief Study of multiplicity distrubtions in EMCAL triggered events
 * @ingroup PWGJETASKS
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Nov 10, 2016
 *
 * Simple study of the multiplicity distributions of various probes (VZERO amplitude,
 * SPD tracklets, global tracks, EMCAL clusters) in min. bias and EMCAL-triggered events.
 */
class AliAnalysisTaskEmcalTriggerMultiplicity : public AliAnalysisTaskEmcalTriggerBase {
public:

  /***
   * Dummy constructor
   */
  AliAnalysisTaskEmcalTriggerMultiplicity();

  /**
   * Named constructor: Default constructor for users
   * @param[in] name Name of the task.
   */
  AliAnalysisTaskEmcalTriggerMultiplicity(const char *name);

  /**
   * Destructor
   */
  virtual ~AliAnalysisTaskEmcalTriggerMultiplicity();

  /**
   * Enable Sumw2 when creating the histograms. Attention: Enabling Sumw2
   * will increase memory consumption significantly. Option should only be
   * used in case histograms are filled with a weight.
   * @param[in] doEnable If true Sumw2 is enabled for all histograms
   */
  void EnableSumw2(Bool_t doEnable) { fEnableSumw2 = doEnable; }

  /**
   * Configure internal virtual track selection used for the
   * determination of the track multiplicity according to
   * predefined cut configurations.
   * @param[in] cutname Name of the cut configuration
   * @param[in] isAOD Switch between ESD and AOD mode
   */
  void InitializeTrackCuts(const TString &cutname, bool isAOD);

  /**
   * Set the virtual track selection, used to determine the
   * track multiplicity
   * @param[in] sel Track selection object
   */
  void SetEmcalTrackSelection(AliEmcalTrackSelection *sel) { fTrackSel = sel; }


  /**
   * Create and configure trigger mutiplicity task.
   * @param[in] nclusters Name of the cluster container (default: "usedefault" - auto-configures the task)
   * @param[in] suffix Container name suffix
   * @return Fully configured task
   */
  static AliAnalysisTaskEmcalTriggerMultiplicity *AddTaskEmcalTriggerMultiplicity(const TString &nclusters = "usedefault", const TString &suffix = "");


protected:

  /**
   * Creating histograms for the multiplicity distributions
   */
  virtual void CreateUserHistos();

  /**
   * No runtime user objects needed in this task
   */
  virtual void CreateUserObjects();

  /**
   * Filling event coutner histograms for the different trigger classes
   */
  virtual void UserFillHistosAfterEventSelection();

  /**
   * Determined multiplicities (V0A/C, tracklet, track, EMCAL clusters)
   * and fill histograms
   * @return always true
   */
  virtual bool Run();

private:
  AliEmcalTrackSelection     *fTrackSel;         ///< EMCAL virtual track selection
  Bool_t                                  fEnableSumw2;       ///< Setter for enabling Sumw2

  AliAnalysisTaskEmcalTriggerMultiplicity(const AliAnalysisTaskEmcalTriggerMultiplicity &);
  AliAnalysisTaskEmcalTriggerMultiplicity &operator=(const AliAnalysisTaskEmcalTriggerMultiplicity &);

  ClassDef(AliAnalysisTaskEmcalTriggerMultiplicity, 1);
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIANALYSISTASKEMCALTRIGGERMULTIPLICITY_H */
