#ifndef ALIANALYSISTASKEMCALTRIGGERBASE_H
#define ALIANALYSISTASKEMCALTRIGGERBASE_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TString.h>
#include <vector>

#include "AliAnalysisTaskEmcal.h"
#include "AliCutValueRange.h"
#include "AliEmcalTriggerOfflineSelection.h"

class TClonesArray;
class THistManager;
class AliOADBContainer;

namespace EMCalTriggerPtAnalysis {

/**
 * @class AliAnalysisTaskEmcalTriggerBase
 * @brief Base class for analyses using EMCAL triggers
 * @ingroup PWGJETASKS
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @date Oct.  6, 2016
 *
 * Base class for analyses using EMCAL-triggered events. This class provides several
 * extra funcitonalities shared among different tasks
 * - Handling of downscale weights for downscaled triggers
 * - Handling of online trigger cleanup (requirement of recalc patch above threshold to select event)
 * - Trigger classification, also for exclusive classes (classes not containing lower classes)
 * - Histogram management via THistManager
 *
 * This class is abstract. Users have to implement at least the following two functions
 * ~~~{.cxx}
 * void CreateUserHistos();     // Instanciation of histograms within the histogram manager
 * void CreateUserObjects();    // Handling of other objects the user might need in the analysis
 * ~~~
 *
 * The task runs a common (default) event selection. Users can implement additional event selection
 * ~~~{.cxx}
 * Bool_t IsUserEventSelected();
 * ~~~
 * This function is running only in case the common event selection is passed.
 *
 * Two functions implement monitoring
 * ~~~{.cxx}
 * void UserFillHistosBeforeEventSelection();   // Filling distribuions before event selection
 * void UserFillHistosAfterEventSelection();    // Filling distribuions after event selection
 * ~~~
 *
 * Otherwise the class uses functionality of the AliAnalysisTaskEmcal. This includes the main
 * function
 * ~~~{.cxx}
 * Bool_t Run();
 * ~~~
 * which contains the main event loop.
 */
class AliAnalysisTaskEmcalTriggerBase : public AliAnalysisTaskEmcal{
public:

  /**
   * Dummy I/O constructor
   */
  AliAnalysisTaskEmcalTriggerBase();

  /**
   * Main Constructor: Initializes the task and sets the name
   * @param[in] name Name of the task
   */
  AliAnalysisTaskEmcalTriggerBase(const char *name);

  /**
   * Destructor
   */
  virtual ~AliAnalysisTaskEmcalTriggerBase();

  /**
   * Set the name of the OADB container with the downscale factors.
   * Once it is available, downscale weights can be obtained via
   * ~~~{.cxx}
   * // Get the weight of the EG2 trigger for the given run
   * this->GetTriggerWeight("EG2");
   * ~~~
   * @param[in] oadbname Name of the OADB container with the trigger downscale factors
   */
  void SetDownscaleOADB(TString oadbname) { fNameDownscaleOADB = oadbname; }

  /**
   * If true then noise events (events without recalc trigger patch above threshold) are
   * excluded from the analysis.
   * @param[in] doExclude If true then noise events are excluded from the analysis
   */
  void SetExcludeNoiseEvents(Bool_t doExclude = true) { fRejectNoiseEvents = doExclude; }

  /**
   * If true then noise events (events without recalc trigger patch above threshold) are
   * explicitly selected for the analysis.
   * @param[in] doSelect If true only noise events are used in the analysis
   */
  void SetSelectNoiseEvents(Bool_t doSelect = true) { fSelectNoiseEvents = doSelect; }

  /**
   * Add absolute ID of a FastOR to be masked (excluded from trigger patches)
   * @param[in] fastorID Absolute ID of a fastor to be masked
   */
  void AddMaskedFastor(int fastorID) { fMaskedFastors.push_back(fastorID); }

  /**
   * Set the name of the file with the OADB container containing the masked FastORs
   * @param[in] oadbname Name of the OADB container file
   */
  void SetMaskedFastorOADB(TString oadbname) { fNameMaskedFastorOADB = oadbname; }

  /**
   * Set an offline trigger selection (selection of trigger according to the presence
   * of an trigger patch from cells above energy)
   * @param[in] sel Trigger offline selection
   */
  void SetOfflineTriggerSelection(AliEmcalTriggerOfflineSelection *sel) { fTriggerSelection = sel; }

  /**
   * Set z-range of the primary vertex which is selected
   * @param[in] zmin Min. allowed z-value
   * @param[in] zmax Max. allowed z-value
   */
  void SetVertexCut(double zmin, double zmax) { fVertexCut.SetLimits(zmin, zmax); }

  /**
   * Specify whether the trigger decision should be done from trigger patches
   * @param doUse If true the trigger string is rebuilt from recalc patches
   */
  void UseTriggerPatches(Bool_t doUse) { fTriggerStringFromPatches = doUse; }

  /**
   * Setting trigger threshold for online trigger selection
   * @param[in] triggerclass Name of the trigger class
   * @param[in] threshold Online trigger threshold
   */
  void SetOnlineTriggerThreshold(const TString &triggerclass, Int_t threshold);

  /**
   * Set the name of teh OADB container with trigger acceptance maps. Trigger
   * acceptance maps will be used in the trigger offline selection to mimic
   * the acceptance observed in data. Only useful in simulations and the mimicing
   * of the trigger in min. bias data.
   * @param[in] nameAcceptanceOADB Location of the OADB container with the acceptance maps
   */
  void SetTiggerAcceptanceOADB(const TString &nameAcceptanceOADB) { fNameAcceptanceOADB = nameAcceptanceOADB; }

protected:
  /**
   * Steering of object creation. Handles general objects (histogram manager
   * and analysis utils) and steers user objects
   * - CreateUserHistos
   * - CreateUserObjects
   */
  virtual void UserCreateOutputObjects();

  /**
   * Run default event selection
   * - Vertex-z cut
   * - Pileup cut
   * Default event selection contains also trigger selection. In
   * case a user event selection is provided on top, it will be
   * executed after the default event selection.
   * @return True if the event is selected
   */
  virtual bool IsEventSelected();

  /**
   * New framework function: can be used by the user to implement
   * an event selection which extends the common event selection.
   * @return True if the event is selected or then function is not implemented
   */
  virtual bool IsUserEventSelected() { return true; }

  /**
   * New framework function: Implemented by users to create
   * histograms within the common histogram handler. Called
   * in UserCreateOutputObjects
   */
  virtual void CreateUserHistos() = 0;

  /**
   * New framework function: Create or initialize objects needed
   * by the user which are not handled by the EMCAL framework.
   */
  virtual void CreateUserObjects() = 0;

  /**
   * Perform gloabl initializations
   * Used for the moment for
   * - Initialize OADB container with downscaling factors
   * - Defining STU online trigger thresholds
   */
  virtual void ExecOnce();

  /**
   * Run change method. Called when the run number of the new event
   * is different compared to the run number of the previous event.
   * Used for loading of the downscale factor for a given
   * run from the downscale OADB.
   * @param[in] runnumber Number of the new run.
   */
  virtual void RunChanged(Int_t runnuber);

  /**
   * New framework function: Can be used by the user to fill histograms
   * before event selection
   */
  virtual void UserFillHistosBeforeEventSelection() { }

  /**
   * New framework function: Can be used by the user to full histograms
   * after event selection.
   */
  virtual void UserFillHistosAfterEventSelection() { }

  /**
   * Creates a list of trigger classes supported by this framework.
   * It can be used by the users when creating or filling histograms
   * according to trigger classes.
   * @return List of supported trigger classes
   */
  std::vector<TString> GetSupportedTriggers();

  /**
   * Get a trigger class dependent event weight. The weight
   * is defined as 1/downscalefactor. The downscale factor
   * is taken from the OADB. For triggers which are not downscaled
   * the weight is always 1.
   * @param[in] triggerclass Class for which to obtain the trigger.
   * @return Downscale facror for the trigger class (1 if trigger is not downscaled or no OADB container is available)
   */
  Double_t GetTriggerWeight(const TString &triggerclass) const;

  /**
   * Apply trigger selection using offline patches and trigger thresholds based on offline ADC Amplitude
   * @param[in] triggerpatches Trigger patches found by the trigger maker
   * @return String with EMCAL trigger decision
   */
  TString GetFiredTriggerClassesFromPatches(const TClonesArray* triggerpatches) const;

  /**
   * Steering of the trigger selection:
   * Combines the selection of triggers from event trigger string,
   * offline energy selection, and online noise rejection / selection.
   * Also handles exclusive trigger classes (classes which do not contain
   * triggers from lower classes).
   */
  void TriggerSelection();

  /**
   * Second approach: We assume masked fastors are already handled in the
   * trigger maker. In this case we treat masked fastors similarly to online
   * masked fastors and apply online cuts on the recalc ADC value. Implemented
   * for the moment only for L1 triggers.
   * @return True if the event has at least 1 recalc patch above threshold
   */
  bool SelectOnlineTrigger(AliEmcalTriggerOfflineSelection::EmcalTriggerClass trigger) const;

  /**
   * Checks whether online trigger thresholds are initialized. All trigger
   * classes are required to be set for this.
   * @return True if thresholds are initialized
   */
  bool OnlineThresholdsInitialized() const;

  /**
   * Get STU online trigger threshold by the index in AliEmcalTriggerOfflineSelection
   * @param[in] trg Index of the trigger class
   * @return Online trigger threshold
   */
  Int_t GetOnlineTriggerThresholdByIndex(AliEmcalTriggerOfflineSelection::EmcalTriggerClass trg) const;

  /**
   * Get STU online trigger threshold by the name of the online trigger class
   * @param[in] name Name of the trigger class
   * @return Online trigger threshold
   */
  Int_t GetOnlineTriggerThresholdByName(const TString &name) const;

  /**
   * Select trigger patches firing the trigger for patches above threshold for
   * a given trigger class.
   * @param[in] triggerclass Name of the trigger class from which to apply the threshold
   * @param[in] adc ADC value of the trigger patch at Level1
   * @return True if the patch is selected
   */
  Bool_t SelectFiredPatch(const TString &triggerclass, Int_t adc) const;


  THistManager                    *fHistos;                   ///< Task Histogram container

  AliEmcalTriggerOfflineSelection *fTriggerSelection;         ///< Offline trigger selection
  Bool_t                          fTriggerStringFromPatches;  ///< Do rebuild the trigger string from trigger patches
  std::vector<TString>            fSelectedTriggers;          //!<! Triggers selected for given event

  AliCutValueRange<double>        fVertexCut;                 ///< Cut on the z-position of the primary vertex

  TString                         fNameDownscaleOADB;         ///< Name of the downscale OADB container
  AliOADBContainer                *fDownscaleOADB;            //!<! Container with downscale factors for different triggers
  TObjArray                       *fDownscaleFactors;         //!<! Downscalfactors for given run
  TString                         fNameMaskedFastorOADB;      ///< Name of the masked fastor OADB container
  AliOADBContainer                *fMaskedFastorOADB;         //!<! Container with masked fastors
  std::vector<int>                fMaskedFastors;             ///< List of masked fastors
  TObjArray                       fOnlineTriggerThresholds;   ///< Trigger thresholds applied at online level
  TString                         fNameAcceptanceOADB;        ///< Name of the OADB container with the trigger acceptance

  Bool_t                          fSelectNoiseEvents;         ///< Explicitly select events triggered only by noisy fastors
  Bool_t                          fRejectNoiseEvents;         ///< Reject events triggered by noisy fastors

private:
  AliAnalysisTaskEmcalTriggerBase(const AliAnalysisTaskEmcalTriggerBase &);
  AliAnalysisTaskEmcalTriggerBase &operator=(const AliAnalysisTaskEmcalTriggerBase &);
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIANALYSISTASKEMCALTRIGGERBASE_H */
