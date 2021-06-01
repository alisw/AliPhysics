#ifndef ALIANALYSISTASKEMCALTRIGGERBASE_H
#define ALIANALYSISTASKEMCALTRIGGERBASE_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TString.h>
#include <string>
#include <vector>

#include "AliAnalysisTaskEmcal.h"
#include "AliCutValueRange.h"
#include "AliEmcalTriggerOfflineSelection.h"

class TClonesArray;
class THistManager;
class AliOADBContainer;

namespace PWG { namespace EMCAL { class AliEmcalTriggerDecisionContainer; } }

namespace PWGJE {

namespace EMCALJetTasks {

/**
 * @class AliAnalysisTaskEmcalTriggerBase
 * @brief Base class for analyses using EMCAL triggers
 * @ingroup PWGJETASKS
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @date Oct.  6, 2016
 *
 * Base class for analyses using EMCAL-triggered events. This class provides several
 * extra functionalities shared among different tasks
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
 * void UserFillHistosBeforeEventSelection();   // Filling distributions before event selection
 * void UserFillHistosAfterEventSelection();    // Filling distributions after event selection
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
   * @brief Dummy I/O constructor
   */
  AliAnalysisTaskEmcalTriggerBase();

  /**
   * @brief Main Constructor
   *
   * Initializes the task and sets the name
   *
   * @param[in] name Name of the task
   */
  AliAnalysisTaskEmcalTriggerBase(const char *name);

  /**
   * @brief Destructor
   */
  virtual ~AliAnalysisTaskEmcalTriggerBase();

  /**
   * @brief Enable / Disable histograms for the DCAL triggers.
   *
   * If true, DCAL trigger classes are added to the list of
   * supported triggers and the online/offline trigger decision
   * for DCAL triggers is evaluated.
   *
   * @param[in] doEnable If true DCAL triggers are enabled
   */
  void EnableDCALTriggers(Bool_t doEnable) { fEnableDCALTriggers = doEnable; }

  /**
   * @brief Enable EMCAL/DCAL combined triggers (OR of EMCAL and DCAL triggers at same threshold)
   * @param doEnable If true EMCAL and DCAL combined triggers are enabled
   */
  void EnableEDCombinedTriggers(Bool_t doEnable) { fEnableEDCombinedTriggers = doEnable; }

  /**
   * @brief Enable T0-based (INT8, EMC8, DMC8) trigger suite (Default: Off)
   * @param doEnable If true T0-based triggers are enabled
   */
  void EnableT0Triggers(Bool_t doEnable) { fEnableT0Triggers = doEnable; }

  /**
   * @brief Enable VZERO-based (INT7, EMC7, DMC7) trigger suite (Default: On)
   * @param doEnable If true VZERO-based triggers are enabled
   */
  void EnableVZEROTriggers(Bool_t doEnable) { fEnableV0Triggers = doEnable; }

  /**
   * @brief Enable EMCAL triggers without coincidence with INT triggers
   * 
   * Exotic case, only of relevance when EMCAL is in PHYSICS_2 or in 
   * quiet beam runs where VZERO is not in readout
   * 
   * @param doEnable If true also EMCAL triggers without INT triggers are enabled
   */
  void EnableNoINTTriggers(Bool_t doEnable) { fEnableNoINTTriggers = doEnable; }

  /**
   * @brief Enable centrality (CENT/SEMICENT) triggers (only relevant for Pb-Pb)
   * @param doEnable If true centrality triggers are enabled
   */
  void EnableCentralityTriggers(Bool_t doEnable) { fEnableCentralityTriggers = doEnable; }

  /**
   * @brief Switch on selection of centrality triggers for PbPb 2018
   * 
   * Selection done based on trigger string using trigger classes for centrality triggers in 2018.
   * Triggers are not yet supported by the physics selection.
   * 
   * @param doSelect If true trigger selection is enabled
   */
  void SetSelectCentralityTriggers2018(Bool_t doSelect) { fSelectCentralityTriggers2018 = doSelect; }

  /**
   * @brief Set the name of the OADB container with the downscale factors.
   *
   * Once it is available, downscale weights can be obtained via
   * ~~~{.cxx}
   * // Get the weight of the EG2 trigger for the given run
   * this->GetTriggerWeight("EG2");
   * ~~~
   *
   * @param[in] oadbname Name of the OADB container with the trigger downscale factors
   */
  void SetDownscaleOADB(EMCAL_STRINGVIEW oadbname) { fNameDownscaleOADB = oadbname.data(); }

  /**
   * @brief Load the downscale factors run-by-run from the OCDB
   * @param[in] doApply If true downscale factors are loaded from the OCDB
   */
  void SetApplyDownscaleCorrectionFromOCDB(Bool_t doApply) { fUseDownscaleCorrectionFormOCDB = doApply; }

  /**
   * @brief Set an offline trigger selection.
   *
   * The offline trigger selection will select events
   * based on the presence of a trigger patch from
   * cells above energy.
   *
   * @param[in] sel Trigger offline selection
   */
  void SetOfflineTriggerSelection(AliEmcalTriggerOfflineSelection *sel) { fTriggerSelection = sel; }

  /**
   * @brief Providing access to the offline trigger selection.
   *
   * Note that the trigger offline selection need to be
   * defined from outside the task.
   *
   * @return Trigger offline selection
   */
  AliEmcalTriggerOfflineSelection *GetOfflineTriggerSelection() const { return fTriggerSelection; }

  /**
   * @brief Switch on/off vertex cuts
   * @param[in] doUse If true vertex cuts will be applied
   */
  void SetApplyVertexCuts(bool doUse) { fApplyVertexCuts = doUse; }

  /**
   * @brief Set z-range of the primary vertex which is selected
   * @param[in] zmin Min. allowed z-value
   * @param[in] zmax Max. allowed z-value
   */
  void SetVertexCut(double zmin, double zmax) { fVertexCut.SetLimits(zmin, zmax); }

  /**
   * @brief Use vertex from SPD
   * 
   * For productions without TPC (i.e. CALOFAST cluster, muon_calo_pas1)
   * @param[in] doUse If true the SPD vertex will be used instread of the vertex from tracks
   */
  void SetUseSPDVertex(bool doUse) { fUseSPDVertex = doUse; }

  /**
   * @brief Defining whether to require trigger bits.
   *
   * Attention: Relies on the presence of the physics selection - to be
   * switched off in case of new data. By default trigger bits are required.
   *
   * @param[in] doUse If true trigger bits are required, if false not
   */
  void SetUseTriggerBits(Bool_t doUse) { fUseTriggerBits = doUse; }

  /**
   * @brief Defining whether to require bunch crossing events.
   *
   * To be switch off in case of noise studies (i.e. beam-empty
   * or empty-empty events). By default bunch crossing is required.
   *
   * @param[in] doRequire If true bunch-bunch crossing is required.
   */
  void SetRequireBunchCrossing(Bool_t doRequire) { fRequireBunchCrossing = doRequire; }

  /**
   * @brief Define whether cuts in AliAnalysisUtils are used in the event selection.
   *
   * Not to be applied on new data or data from muon_calo_pass1
   * @param[in] doRequire If true AliAnalysisUtils are used in addition for the event selection
   */
  void SetRequireAnalysisUtils(Bool_t doRequire) { fRequireAnalysisUtils = doRequire; }

  /**
   * @brief Run event loop only on min. bias events.
   *
   * In this case EMCAL triggers are ignored, and the trigger selection code is not run.
   * Also the calo trigger patch object container is not linked to the task / not required.
   * @param[in] exclusivemb If true only min. bias events are analyzed
   */
  void SetExclusiveMinBias(Bool_t exclusivemb) { fExclusiveMinBias = exclusivemb; SetCaloTriggerPatchInfoName(""); }

  /**
   * @brief Use trigger selection container in addition to trigger string
   * 
   * @param[in] doUse If true results from the trigger decision container are used in addition 
   */
  void SetUseTriggerSelectionContainer(Bool_t doUse) { fUseTriggerSelectionContainer = doUse; }

  /**
   * @brief Set the name of the trigger decision container
   * 
   * @param[in] nameCont Name of the trigger decision container 
   */
  void SetNameTriggerSelectionContainer(EMCAL_STRINGVIEW nameCont) { fNameTriggerSelectionContainer = nameCont.data();}

protected:

  /**
   * @brief Steering of object creation.
   *
   * Handles general objects (histogram manager
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
   * @param useExclusiveTriggers If true also exclusive triggers are added
   * @return List of supported trigger classes
   */
  std::vector<TString> GetSupportedTriggers(Bool_t useExclusiveTriggers = true) const;

  /**
   * Get a trigger class dependent event weight. The weight
   * is defined as 1/downscalefactor. The downscale factor
   * is taken from the OADB. For triggers which are not downscaled
   * the weight is always 1.
   * @param[in] triggerclass Class for which to obtain the trigger.
   * @return Downscale facror for the trigger class (1 if trigger is not downscaled or no OADB container is available)
   */
  Double_t GetTriggerWeight(EMCAL_STRINGVIEW triggerclass) const;

  /**
   * Steering of the trigger selection:
   * Combines the selection of triggers from event trigger string,
   * offline energy selection, and online noise rejection / selection.
   * Also handles exclusive trigger classes (classes which do not contain
   * triggers from lower classes).
   */
  void TriggerSelection();

  /**
   * @brief Match trigger pattern
   * 
   * Trigger pattern can
   * - start with = - full class name must be found
   * - contain strings separated with one or more | - any of the classes to be found 
   * 
   * @param pattern Pattern used for checking
   * @param triggerclass Trigger class to be checked
   * @return true Pattern found
   * @return false Pattern not found
   */
  bool MatchTriggerFromPattern(EMCAL_STRINGVIEW pattern, EMCAL_STRINGVIEW triggerclass) const;

  /**
   * @brief Matching triggers in pattern with entry in trigger decision container
   * 
   * Trigger will ignore =. | will be used to split several triggers to be checked
   * 
   * @param pattern Pattern to be checked
   * @param trgcont Trigger container with selected classes from patches
   * @return true Pattern found in the trigger decision container 
   * @return false Pattern not found in the trigger decision container
   */
  bool MatchTriggerFromContainer(EMCAL_STRINGVIEW pattern, const PWG::EMCAL::AliEmcalTriggerDecisionContainer *trgcont) const;

  /**
   * Define name of the cluster container used to read EMCAL cluster information
   * from
   * @param[in] clustercontname Name of the cluster container
   */
  void SetClusterContainer(EMCAL_STRINGVIEW clustercontname) { fNameClusterContainer = clustercontname.data(); }

  /**
   * @brief Read the downscale factors from the OCDB
   */
  void PrepareDownscaleFactorsFormOCDB();


  THistManager                    *fHistos;                   ///< Task Histogram container

  Bool_t                          fUseTriggerBits;            ///< Switch whether using trigger bits (relies on physics selection)
  Bool_t                          fRequireBunchCrossing;      ///< Require bunch-bunch events (tag -B- in trigger string)
  Bool_t                          fUseDownscaleCorrectionFormOCDB; ///< Use downscale factors from OCDB
  AliEmcalTriggerOfflineSelection *fTriggerSelection;         ///< Offline trigger selection
  std::vector<TString>            fSelectedTriggers;          //!<! Triggers selected for given event
  TString                         fNameClusterContainer;      ///< Name of the cluster container in the event

  Bool_t                          fRequireAnalysisUtils;      ///< Switch whether to require event selection in AliAnalysisUtils
  Bool_t                          fUseSPDVertex;              ///< Use SPD vertex (for productions without TPC)
  Bool_t                          fApplyVertexCuts;           ///< Apply vertex cuts (default: True)
  AliCutValueRange<double>        fVertexCut;                 ///< Cut on the z-position of the primary vertex

  TString                         fNameDownscaleOADB;         ///< Name of the downscale OADB container
  AliOADBContainer                *fDownscaleOADB;            //!<! Container with downscale factors for different triggers
  TObjArray                       *fDownscaleFactors;         //!<! Downscalfactors for given run
  TString                         fNameTriggerSelectionContainer; ///< Name of the trigger selection container

  Bool_t                          fEnableDCALTriggers;        ///< Enable / Disable event selection for DCAL trigger classes
  Bool_t                          fEnableEDCombinedTriggers;  ///< Enable OR combination of EMCAL and DCAL combined triggers
  Bool_t                          fEnableV0Triggers;          ///< Enable VZERO-based triggers (default)
  Bool_t                          fEnableT0Triggers;          ///< Enable triggers depending on T0 (INT8, EMC8, EMC8EGA, EMC8EJE) - default off
  Bool_t                          fEnableNoINTTriggers;       ///< Process EMCAL triggers without coincidence with INT triggers - exotic case - default off
  Bool_t                          fEnableCentralityTriggers;  ///< Enable central / semi-central trigger
  Bool_t                          fExclusiveMinBias;          ///< Only look at Min. Bias trigger
  Bool_t                          fUseTriggerSelectionContainer;    ///< Use trigger decision in trigger selection container
  Bool_t                          fSelectCentralityTriggers2018;    ///< Select centrality triggers 2018 based on trigger string (missing support by physics selection yet)

private:
  AliAnalysisTaskEmcalTriggerBase(const AliAnalysisTaskEmcalTriggerBase &);
  AliAnalysisTaskEmcalTriggerBase &operator=(const AliAnalysisTaskEmcalTriggerBase &);
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIANALYSISTASKEMCALTRIGGERBASE_H */
