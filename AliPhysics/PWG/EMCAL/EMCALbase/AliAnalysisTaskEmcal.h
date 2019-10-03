/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#ifndef ALIANALYSISTASKEMCAL_H
#define ALIANALYSISTASKEMCAL_H

class TClonesArray;
class TString;
class AliMCParticle;
class AliVCluster;
class AliVTrack;
class AliVParticle;
class AliVCaloCells;
class TH1;
class TProfile;
class AliEMCALGeometry;
class AliGenPythiaEventHeader;
class AliGenHerwigEventHeader;
class AliVCaloTrigger;
class AliAnalysisUtils;
class AliEMCALTriggerPatchInfo;
class AliAODTrack;
class AliEmcalPythiaInfo;
class AliAODInputHandler;
class AliESDInputHandler;

#include "Rtypes.h"
#include "TArrayI.h"

#include "AliParticleContainer.h"
#include "AliMCParticleContainer.h"
#include "AliTrackContainer.h"
#include "AliClusterContainer.h"
#include "AliEmcalList.h"
#include "AliEventCuts.h"
#include "AliEmcalStringView.h"


#include "AliAnalysisTaskSE.h"
/**
 * @class AliAnalysisTaskEmcal
 * @brief Base task in the EMCAL framework
 * @ingroup EMCALCOREFW
 * @author Marta Verweij
 * @author Salvatore Aiola
 *
 * # Extra features of the AliAnalysisTaskEmcal
 * 
 * The AliAnalysisTaskEmcal adds common steps needed in the work flow 
 * of EMCAL-related or jet analyses. It inherits itself from AliAnalysisTaskSE,
 * and user tasks can inherit from it in order to profit from the extra 
 * features - the task is fully compatible with the ALICE Analysis Framework.
 * 
 * ## Virtual functions to be implemented in user tasks
 * 
 * In contrast to the normal AliAnalysisTaskSE, the main event
 * loop function to be implemented by the user is called Run.
 * This function is only called in case the event was selected
 * previously.
 * 
 * ~~~{.cxx}
 * bool MyUserTask::Run() {
 *   // User code here ...
 * 
 *   return kTRUE;
 * }
 * ~~~
 * 
 * The function is called after the event selection for every event passing the 
 * event selection. AliAnalysisTaskEmcal offers a second function called for each
 * event called FillHistograms. The functions focus is on filling histograms and
 * is only called in case the Run returns true
 * 
 * ~~~{.cxx}
 * bool MyUserTask::FillHistograms() {
 *   // Fill user histograms
 *   myEventCounter->Fill(1.);
 * }
 * ~~~
 * 
 * Implementation of FillHistograms is not mandatory, user histograms can also be
 * filled in Run. 
 * 
 * User output objects are created in UserCreateOutputObjects. As the AliAnalysisTaskEmcal
 * performs initalizations in this function itsel, the user tasks must call the parent
 * UserCreateOutputObjects. AliAnalysisTaskEmcal provides a combined output container fOutput
 * where users can add their histograms to. The main constructor of a user task should construct
 * the AliAnalysisTaskEmcal using its named constructor with the second argument true.
 * 
 * ~~~{.cxx}
 * void MyUserTask::UserCreateOutputObjects() {
 *   AliAnalysisTaskEmcal::UserCreateOutputObjects();
 * 
 *   // User initialization here
 *   fEventCounter = new TH1F("fEventCounter", "Event counter", 1, 0.5, 1.5);
 * 
 *   // Add histograms to the combined output container
 *   fOutput->Add(fEventCounter);
 *   PostData(1, fOutput);
 * }
 * ~~~
 * 
 * AliAnalysisTaskEmcal provides several optional interface functions which the user
 * can implement. 
 * 
 * | Function name     | Default implementation | Purpose                                            | 
 * |-------------------| -----------------------|----------------------------------------------------|
 * | IsEventSelected   | Yes                    | Perform event selection (default: AliEventCuts)    |
 * | IsTriggerSelected | Yes                    | Perform trigger selection (default: INT7)          |
 * | RunChanged        | No                     | Initializations when a new run starts (i.e. OADB)  |
 * | UserExecOnce      | No                     | Initializations which need knowledge about dataset |
 * | UserFileChanged   | No                     | Tasks to be performed when an input file changes   |
 * 
 * ## Connecting EMCAL container
 *
 * A key feature is the handling of EMCAL containers (cluster/
 * particle/track). Users can create containers and attach it
 * to their analysis task.
 *
 * ~~~{.cxx}
 * AliParticleContainer *cont = task->AddParticleContainer("MCParticlesSelected");
 * cont->SetEtaLimits(-0.7, 0.7);
 * ~~~
 *
 * Note that the name in the container must match the name of the
 * TClonesArray in the input event connected to the container.
 * Containers can be accessed inside the task either by their name
 * or by their index as they were added to the task. The indices
 * of cluster-, particle- or track-containers are not mixed. Containers
 * provide an easy access to content created by other tasks and
 * attached to the event as a TClonesArray.
 * 
 * ## Handling of jet-jet productions in pt-hard bins
 * 
 * Productions in pt-hard bins require special treatent. User histogrmas 
 * need to be weighted by the cross section / number of trials  in order 
 * to be physically meaningfull. In addition productions can contain outlier
 * events, where generated tracks/clusters/jets contain a pt much larger
 * than the hard process of the event. Those events need to be removed,
 * otherwise they distort the distributions. The AliAnalysisTaskEmcal
 * provides helper to handle these special cases.
 * 
 * ### Automatic weighting of user histograms 
 * 
 * User results are stored in the list fOutput (of type AliEmcalList).
 * When merging the AliEmcalList users can switch on scaled merging.
 * In this case the AliEmcalList determines automatically when different
 * pt-hard bins are added and scales the user histograms before. In order
 * to use the feature the dataset must be declared as a pt-hard production.
 * In case the pt-hard binning is different from the default pt-hard binning
 * (10 pt-hard bins) it need to be specified. The following example setup the
 * AliAnalysisTaskEmcal for a 20-bin sample
 * 
 * ~~~{.cxx}
 * // In dataset configuration
 * const int knPthardBins=21;
 * TArrayD kPtHardBinning;
 * kPtHardBinning.Set(22);
 * Int_t binning[]={0,5, 7, 9, 12, 16, 21, 28, 36, 45, 57, 70, 85, 99, 115, 132, 150, 169, 190, 212, 235,1000};
 * for(Int_t bin=0;bin<22;bin++){ kPtHardBinning[bin]=binning[bin]; }
 * 
 * // In wagon configuration
 * __R_ADDTASK__->SetIsPythia(kTRUE);              // Mark production as pt-hard production
 * __R_ADDTASK__->SetMakeGeneralHistograms(kTRUE); // Create weighting histograms
 * __R_ADDTASK__->SetUsePtHardBinScaling(kTRUE); // Switch on pt-hard bin scaling of the output histograms
 * // Set the pt-hard binning of the production
 * __R_ADDTASK__->SetNumberOfPtHardBins(knPthardBins);
 * if(knPthardBins!=11) __R_ADDTASK__-> SetUserPtHardBinning(kPtHardBinning);
 * ~~~
 * 
 * In case the production is a production which does not yet support the pt-hard
 * structure (i.e. LHC11a1) but the pt-hard bins are produced as separate samples, 
 * the production must be marked as old production. In this case the pt-hard bin
 * cannot be determined automatically, so the cross sections and weights are put
 * into the first bin for all pt-hard bins. The following example marks the production
 * as old pt-hard production:
 * 
 * ~~~{.cxx}
 * __R_ADDTASK__->SetIsPythia(kTRUE);
 * __R_ADDTASK__->SetGetPtHardBinFromPath(kFALSE);
 * ~~~
 * 
 * Note that for these samples the automatic weighting by the EMCAL list does not
 * work and must be switched off. Users must weight the output manually.
 * 
 * ### Rejection of outlierts
 * 
 * Outlier events can distort distributions after weighting. The reason is that single
 * entries with a pt much larger than the pt-hard get a weight which might be an order 
 * of magnitude or more larger than the weight of entries in the proper pt-hard bin. 
 * Outlier rejection in the AliAnalysisTaskEmcal is performed rejecting events which have
 * at least one jet / track / cluster with a pt > event pt-hard * scale factor. In order
 * to run the outlier rejection it has to be switched on and the scaling factor needs to
 * be set. The following example enables the outlier rejection with a scaling factor of 3
 * for jets.
 * 
 * ~~~{.cxx}
 * __R_ADDTASK__->SetMCFilter();
 * __R_ADDTASK__->SetJetPtFactor(3.);
 * ~~~
 * 
 * This means that the event is rejected if it has at least a jet with a pt > 3 * event pt-hard.
 *
 * For more information refer to \subpage EMCALAnalysisTask
 */
class AliAnalysisTaskEmcal : public AliAnalysisTaskSE {
 public:

  /**
   * @enum EDataType_t
   * @brief Switch for the data type
   */
  enum EDataType_t {
    kUnknownDataType,
    kESD,
    kAOD
  };

  /**
   * @enum BeamType
   * @brief Switch for the beam type
   */
  enum BeamType {
    kNA       = -1,//!< Undefined
    kpp       = 0, //!< Proton-Proton
    kAA       = 1, //!< Nucleus-Nucleus
    kpA       = 2  //!< Proton-Nucleus
  };

  /**
   * @enum TriggerType
   * @brief Switch for EMCAL trigger types
   */
  enum TriggerType {
    kND       = -1,//!< Undefined
    kJ1       = 0, //!< EMCAL Level1 jet trigger, high threshold
    kJ2       = 1, //!< EMCAL Level1 jet trigger, low threshold
    kG1	      = 2, //!< EMCAL Level1 gamma trigger, high threshold
    kG2       = 3, //!< EMCAL Level1 gamma trigger, low threshold
    kL0       = 4  //!< EMCAL Level0 trigger
  };

  /**
   * @enum TriggerCategory
   * @brief Online trigger categories
   */
  enum TriggerCategory {
    kTriggerLevel0      = 0,   //!< Level0 trigger patch
    kTriggerLevel1Jet   = 1,   //!< Level1 jet trigger patch
    kTriggerLevel1Gamma = 2,   //!< Level1 gamma trigger patch
    kTriggerRecalcJet   = 3,   //!< Recalculated jet trigger patch; does not need to be above trigger threshold
    kTriggerRecalcGamma = 4    //!< kRecalculated gamma trigger patch; does not need to be above trigger threshold
  };

  /**
   * @enum EMCalTriggerMode_t
   * @brief Handling of the EMCAL trigger thresholds
   */
  enum EMCalTriggerMode_t {
    kNoSpecialTreatment,       //!< No special treatment for EMCal triggers
    kOverlapWithLowThreshold   //!< The overlap between low and high threshold trigger is assigned to the lower threshold only
  };

  /**
   * @brief Default constructor.
   */
  AliAnalysisTaskEmcal();

  /**
   * @brief Standard constructor. Should be used by the user.
   *
   * Note: This constructor also handles the general histograms. In
   * case the second parameter is true, then general histograms (see
   * UserCreateOutputObjects and FillHistograms) are created and filled
   * by the task, and a container is provided handling the user histograms.
   * @param[in] name Name of the task
   * @param[in] histo If true then general histograms are filled by the task
   */
  AliAnalysisTaskEmcal(const char *name, Bool_t histo=kFALSE);

  /**
   * @brief Destructor
   */
  virtual ~AliAnalysisTaskEmcal();

  /**
   * @brief Create new particle container and attach it to the task.
   *
   * The name provided to this function must match the name of the array attached
   * to the new container inside the input event.
   * @param[in] n Name of the container and the array the container points to
   * @return Pointer to the new particle container
   */
  AliParticleContainer       *AddParticleContainer(const char *n);

  /**
   * @brief Create new track container and attach it to the task.
   *
   * The name provided to this function must match the name of the array
   * attached to the new container inside the input event.
   * @param[in] n Name of the container and the array the container points to
   * @return Pointer to the new track container
   */
  AliTrackContainer          *AddTrackContainer(const char *n);

  /**
   * @brief Create new container for MC particles and attach it to the task.
   *
   * The name provided to this function must match the name of the array attached
   * to the new container inside the input event.
   * @param[in] n Name of the container and the array the container points to
   * @return Pointer to the new container for MC particles
   */
  AliMCParticleContainer     *AddMCParticleContainer(const char *n);

  /**
   * @brief Create new cluster container and attach it to the task.
   *
   * The name provided to this function must match the name of the array
   * attached to the new container inside the input event.
   * @param[in] n Name of the container and the array the container points to
   * @return Pointer to the new cluster container
   */
  AliClusterContainer        *AddClusterContainer(const char *n);

  void                        AdoptParticleContainer(AliParticleContainer* cont)    { fParticleCollArray.Add(cont)                        ; }
  void                        AdoptTrackContainer(AliTrackContainer* cont)          { AdoptParticleContainer(cont)                        ; }
  void                        AdoptMCParticleContainer(AliMCParticleContainer* cont){ AdoptParticleContainer(cont)                        ; }
  void                        AdoptClusterContainer(AliClusterContainer* cont)      { fClusterCollArray.Add(cont)                         ; }

  /**
   * @brief Get \f$ i^{th} \f$ particle container attached to this task
   * @param[in] i Index of the particle container
   * @return Particle container found for the given index (NULL if no particle container exists for that index)
   */
  AliParticleContainer       *GetParticleContainer(Int_t i=0)         const;

  /**
   * @brief Find particle container attached to this task according to its name
   * @param[in] name Name of the particle container
   * @return Particle container found under the given name
   */
  AliParticleContainer       *GetParticleContainer(const char* name)  const;

  /**
   * @brief Get \f$ i^{th} \f$ cluster container attached to this task
   * @param[in] i Index of the cluster container
   * @return Cluster container found for the given index (NULL if no cluster container exists for that index)
   */
  AliClusterContainer        *GetClusterContainer(Int_t i=0)          const;

  /**
   * @brief Find cluster container attached to this task according to its name
   * @param[in] name Name of the cluster container
   * @return Cluster container found under the given name
   */
  AliClusterContainer        *GetClusterContainer(const char* name)   const;
  AliMCParticleContainer     *GetMCParticleContainer(Int_t i=0)               const { return dynamic_cast<AliMCParticleContainer*>(GetParticleContainer(i))   ; }
  AliMCParticleContainer     *GetMCParticleContainer(const char* name)        const { return dynamic_cast<AliMCParticleContainer*>(GetParticleContainer(name)); }
  AliTrackContainer          *GetTrackContainer(Int_t i=0)                    const { return dynamic_cast<AliTrackContainer*>(GetParticleContainer(i))        ; }
  AliTrackContainer          *GetTrackContainer(const char* name)             const { return dynamic_cast<AliTrackContainer*>(GetParticleContainer(name))     ; }
  AliEventCuts               &GetEventCuts()                                        { return fAliEventCuts; }
  void                        RemoveParticleContainer(Int_t i=0)                    { fParticleCollArray.RemoveAt(i)                      ; } 
  void                        RemoveClusterContainer(Int_t i=0)                     { fClusterCollArray.RemoveAt(i)                       ; } 
  void                        SetCaloCellsName(const char *n)                       { fCaloCellsName     = n                              ; }
  void                        SetCaloTriggerPatchInfoName(const char *n)            { fCaloTriggerPatchInfoName = n                       ; }
  void                        SetCaloTriggersName(const char *n)                    { fCaloTriggersName  = n                              ; }
  void                        SetCentRange(Double_t min, Double_t max)              { fMinCent           = min  ; fMaxCent = max          ; }
  void                        SetCentralityEstimator(const char *c)                 { fCentEst           = c                              ; }

  /**
   * @brief Apply cut on \f$ p_{t} \f$ for all clusters in container with
   * index c
   * @param[in] cut \f$ p_{t} \f$-cut to be applied
   * @param[in] c Index of the cluster container affected by the cut
   */
  void                        SetClusPtCut(Double_t cut, Int_t c=0);

  /**
   * @brief Apply cut on cluster time for clusters in container with
   * index c
   * @param[in] min Min. cluster time
   * @param[in] max Max. cluster time
   * @param[in] c Index of the cluster container affected by the cut
   */
  void                        SetClusTimeCut(Double_t min, Double_t max, Int_t c=0);

  void                        SetEventPlaneVsEmcal(Double_t ep)                     { fEventPlaneVsEmcal = ep                             ; }
  void                        SetForceBeamType(BeamType f)                          { fForceBeamType     = f                              ; }
  void                        SetHistoBins(Int_t nbins, Double_t min, Double_t max) { fNbins = nbins; fMinBinPt = min; fMaxBinPt = max    ; }
  /**
   * @brief Enables internal event selection in embedding by recycling unused events.
   * @param[in] b Enables recycled unused embedded events.
   * @deprecated: This is now handled automatically by the embedding helper, so the option is redundant.
   */
  void                        SetRecycleUnusedEmbeddedEventsMode(Bool_t b)          { AliWarning("Enabling recycling of unused embedded events is deprecated. It is now automatically controlled in the embedding helper. You can remove this call."); }
  void                        SetIsEmbedded(Bool_t i)                               { fIsEmbedded        = i                              ; }

  /**
   * @brief Define production as pythia pt-hard production
   * 
   * In this case the scaling histograms (cross section and number of trials) 
   * are created in case the general histograms are enabled (AliAnalysisTaskEmcal::MakeGeneralHistograms(true)).
   * The cross section and number of trials are read from the associated cross section file.
   * 
   * @param i If true the production is handled as a pt-hard production
   */
  void                        SetIsPythia(Bool_t i)                                 { fIsPythia          = i                              ; }

  /**
   * @brief Define production as herwig (6) pt-hard production
   * 
   * In this case the scaling histograms (cross section and number of trials) 
   * are created in case the general histograms are enabled (AliAnalysisTaskEmcal::MakeGeneralHistograms(true)).
   * The cross section and number of trials are read from the associated cross section file.
   * 
   * @param i If true the production is handled as a pt-hard production
   */
  void                        SetIsHerwig(Bool_t i)                                 { fIsHerwig          = i                              ; }

  /**
   * @brief Enable general histograms
   * 
   * Among general histograms are the QA histograms (vertex distribution, rejection reason), normalization
   * histograms and weighting histograms (in case the production is marked as pt-hard production).
   * 
   * @param g If true general histograms of the AliAnalysisTaskEmcal are created and filled.
   */
  void                        SetMakeGeneralHistograms(Bool_t g)                    { fGeneralHistograms = g                              ; }

  /**
   * @brief Switch on/off getting \f$ p_{t,hard}\f$ bin from the file path.
   *
   * By default the \f$ p_{t,hard}\f$ bin is obtained from the file path in case the
   * task is defined to run on a \f$ p_{t,hard}\f$ production. New samples decode
   * the \f$ p_{t,hard}\f$ bin number in the file path. However this is not supported
   * for
   * @param[in] docheck
   */
  void                        SetGetPtHardBinFromPath(Bool_t docheck)               { fGetPtHardBinFromName = docheck; }

  /**
   * @brief Set the number of \f$ p_{t}\f$-hard bins
   * @param[in] nbins Number of \f$ p_{t}\f$-hard bins
   */
  void                        SetNumberOfPtHardBins(Int_t nbins)                    { fNPtHardBins = nbins; }

  /**
   * @brief Set a non-standard \f$ p_{t}\f$-hard binning
   *
   * The array reflects the bin limits, therefore the size must
   * match the number of \f$ p_{t}\f$-hard bins + 1, otherwise
   * the binning is not used in order to create the bin labels
   *
   * @param[in] binning Non-standard binning to be applied
   */
  void                        SetUserPtHardBinning(const TArrayI &binning)          { fPtHardBinning = binning; }

  void                        SetMCLabelShift(Int_t s)                              { fMCLabelShift      = s                              ; }
  void                        SetMinMCLabel(Int_t s)                                { fMinMCLabel        = s                              ; }
  void                        SetMinNTrack(Int_t min)                               { fMinNTrack         = min                            ; }
  void                        SetZvertexDiffValue(Double_t cut)                     { fZvertexDiff  = cut                            ; }
  void                        SetMinPtTrackInEmcal(Double_t min)                    { fMinPtTrackInEmcal = min                            ; }
  virtual void                SetNCentBins(Int_t n)                                 { fNcentBins         = n                              ; }
  void                        SetNeedEmcalGeom(Bool_t n)                            { fNeedEmcalGeom     = n                              ; }
  void                        SetCountDownscaleCorrectedEvents(Bool_t d)            { fCountDownscaleCorrectedEvents =  d                 ; }
  void                        SetOffTrigger(UInt_t t)                               { fOffTrigger        = t                              ; }

  /**
   * @brief Apply cut on the pseudorapidity \f$ \eta \f$ of the all tracks in the
   * track container with index c
   * @param[in] min Minimum of the allowed track \f$ \eta \f$ range
   * @param[in] max Maximum of the allowed track \f$ \eta \f$ range
   * @param[in] c Index of the particle container affected by the cut
   */
  void                        SetTrackEtaLimits(Double_t min, Double_t max, Int_t c=0);

  /**
   * @brief Apply cut on azimuthal angle \f$ \phi \f$ of the all tracks in the
   * track container with index c
   * @param[in] min Minimum of the allowed track \f$ \phi \f$ range
   * @param[in] max Maximum of the allowed track \f$ \phi \f$ range
   * @param[in] c Index of the particle container affected by the cut
   */
  void                        SetTrackPhiLimits(Double_t min, Double_t max, Int_t c=0);

  /**
   * @brief Apply cut on the transverse momentum \f$ p_{t} \f$ of all tracks in
   * the track container with index c.
   * @param[in] cut \f$ p_{t} \f$-cut to be applied
   * @param[in] c Index of the particle container affected by the cut
   */
  void                        SetTrackPtCut(Double_t cut, Int_t c=0);

  void                        SetTrigClass(const char *n)                           { fTrigClass         = n                              ; }
  void                        SetMinBiasTriggerClassName(const char *n)             { fMinBiasRefTrigger = n                              ; }
  void                        SetTriggerTypeSel(TriggerType t)                      { fTriggerTypeSel    = t                              ; } 
  void                        SetUseAliAnaUtils(Bool_t b, Bool_t bRejPilup = kTRUE) { fUseAliAnaUtils    = b ; fRejectPileup = bRejPilup  ; }

  /**
   * @brief Use internal (old) event selection
   * @param[in] doUse It true use the old internal event selection instead of AliEventCuts
   * @deprecated: Will be removed soon due to a naming clash with AliAnalysisTaskEmbeddingHelper, 
   * please use SetUseBuiltinEventSelection instead
   */
  void                        SetUseInternalEventSelection(Bool_t doUse)            { fUseBuiltinEventSelection = doUse                  ; }

  /**
   * @brief Use internal (old) event selection
   * @param[in] doUse It true use the old internal event selection instead of AliEventCuts
   */
  void                        SetUseBuiltinEventSelection(Bool_t doUse)            { fUseBuiltinEventSelection = doUse                  ; }
  
  /**
   * @brief Set pre-configured event cut object
   * 
   * In order to allow custon configuration set a pre-configured AliEventCuts object. 
   * Attention: The object must be constructed with argument false in order to have the 
   * histograms stored correctly. This sets the usage of the interal event selection to
   * false automatically.
   * 
   * @param[in] cuts Cut object to be used
   */
  void                        SetVzRange(Double_t min, Double_t max)                { fMinVz             = min  ; fMaxVz   = max          ; }
  void                        SetMinVertexContrib(Int_t min)                        { fMinVertexContrib = min                             ; }
  void                        SetUseSPDTrackletVsClusterBG(Bool_t b)                { fTklVsClusSPDCut   = b                              ; }
  void                        SetEMCalTriggerMode(EMCalTriggerMode_t m)             { fEMCalTriggerMode  = m                              ; }
  void                        SetUseNewCentralityEstimation(Bool_t b)               { fUseNewCentralityEstimation = b                     ; }
  void                        SetGeneratePythiaInfoObject(Bool_t b)                 { fGeneratePythiaInfoObject = b                       ; }
  void                        SetPythiaInfoName(const char *n)                      { fPythiaInfoName    = n                              ; }
  const TString&              GetPythiaInfoName()                             const { return fPythiaInfoName                              ; }
  const AliEmcalPythiaInfo   *GetPythiaInfo()                                 const { return fPythiaInfo                                  ; }

  /**
   * @brief Switch on pt-hard bin scaling
   * 
   * If enabled the output histograms are automatically scaled in the 
   * merging prodcess. This requires a two step merging:
   * - 1. All runs merged together (scaling not applied)
   * - 2. All pt-hard bins are merged (scaling applied)
   * 
   * Needs the production to be defined as pt-hard production 
   * (AliAnalysisTaskEmcal::SetIsPythia(kTRUE)).
   * 
   * @param b If true pt-hard bin scaling is enabled.
   */
  void                        SetUsePtHardBinScaling(Bool_t b)                      { fUsePtHardBinScaling = b                            ; }

  /**
   * @brief Switch on MC outlier rejection
   * 
   * Outlier rejection should be only used for productions in pt-hard
   * bins. Events are rejected as outliers in case at least one jet / cluster / track
   * in the event has a pt larger outlier fraction * event pt-hard.
   * 
   * The outlier rejection is only performed if at least either of the jet / cluster / track
   * pt factor is set.
   */
  void                        SetMCFilter()                                         { fMCRejectFilter = kTRUE                             ; }

  /**
   * @brief Switch off MC outlier rejection
   */
  void                        ResetMCFilter()                                       { fMCRejectFilter = kFALSE                            ; }

  /**
   * @brief Set the jet pt factor for the outlier rejection
   * 
   * Events are rejected in case they contain at least one jet with pt > event pt-hard * jet pt factor
   * @param f Jet pt factor
   */
  void                        SetJetPtFactor(Float_t f)                             { fPtHardAndJetPtFactor = f                           ; }

  /**
   * @brief Get the jet pt factor for the outlier rejection
   * @return Jet pt factor
   */
  Float_t                     JetPtFactor()                                         { return fPtHardAndJetPtFactor                        ; }

  /**
   * @brief Set the cluster pt factor for the outlier rejection
   * 
   * Events are rejected in case they contain at least one EMCAL cluster with pt > event pt-hard * cluster pt factor
   * @param f Cluster pt factor
   */
  void                        SetClusterPtFactor(Float_t f)                         { fPtHardAndClusterPtFactor = f                       ; }

  /**
   * @brief Get the cluster pt factor for the outlier rejection
   * @return Cluster pt factor
   */
  Float_t                     ClusterPtFactor()                                     { return fPtHardAndClusterPtFactor                    ; }

  /**
   * @brief Set the track pt factor for the outlier rejection
   * 
   * Events are rejected in case they contain at least one track with pt > event pt-hard * track pt factor
   * @param f Track pt factor
   */
  void                        SetTrackPtFactor(Float_t f)                           { fPtHardAndTrackPtFactor = f                         ; }

  /**
   * @brief Get the track pt factor for the outlier rejection
   * @return Track pt factor
   */
  Float_t                     TrackPtFactor()                                       { return fPtHardAndTrackPtFactor                      ; }

  // Static Utilities
  /**
   * @brief Add an AOD handler to the analysis manager
   * @return pointer to the new AOD handler
   */
  static AliAODInputHandler*  AddAODHandler();

  /**
   * @brief Add a ESD handler to the analysis manager
   * @return pointer to the new ESD handler
   */
  static AliESDInputHandler*  AddESDHandler();

 protected:
  /**
   * @brief Load parton info
   * @param event
   */
  void                        LoadPythiaInfo(AliVEvent *event);

  void                        SetRejectionReasonLabels(TAxis* axis);

  /**
   * @brief Cluster selection
   *
   * Check whether cluster is accepted under the cluster
   * selection specified in the cluster container stored
   * under index c
   * @param[in] clus EMCAL cluster to be accepted
   * @param[in] c Index of the cluster container
   * @return true if cluster is accepted
   * @deprecated Use GetClusterContainer(c)->AcceptCluster(clust) instead
   */
  Bool_t                      AcceptCluster(AliVCluster *clus, Int_t c = 0)      const;

  /**
   * Check whether track is accepted under the track
   * selection specified in the track container stored
   * under index c
   * @param[in] track Track to be checked
   * @param[in] c Index of the track container
   * @return true if track is accepted
   * @deprecated Use GetParticleContainer(c)->AcceptParticle(track) instead
   */
  Bool_t                      AcceptTrack(AliVParticle *track, Int_t c = 0)      const;

  /**
   * @brief Add object to event
   * @param[in] obj Object to be added
   * @param[in] attempt If true don't handle error
   */
  void                        AddObjectToEvent(TObject *obj, Bool_t attempt = kFALSE);

  /**
   * @brief Get particle p if accepted from  container with index c
   * If particle not accepted return 0
   * @param[in] p Index of the particle inside the Particle container
   * @param[in] c Index of the particle container
   * @return Accepted particle (NULL if no accepted particle with the index is found inside the container)
   */
  AliVParticle               *GetAcceptParticleFromArray(Int_t p, Int_t c=0)     const;

  /**
   * @brief Get cluster cl if accepted from  container c
   * If particle not accepted return 0
   * @param[in] cl Index of the cluster inside the container
   * @param[in] c Indexd of the cluster container
   * @return Accepted cluster (NULL if no accepted cluster with index is found)
   */
  AliVCluster                *GetAcceptClusterFromArray(Int_t cl, Int_t c=0)     const;


  /**
   * @brief Read a TClonesArray from event.
   *
   * Attention: Both the name of the array and the name of the object
   * stored inside must match.
   * @param[in] name Name of the array to be read in
   * @param[in] clname Name of the type of the objects stored in the array
   * @return Pointer to the TClonesArray (NULL if not found)
   */
  TClonesArray               *GetArrayFromEvent(const char *name, const char *clname=0);

  /**
   * @brief Get beam type
   *
   * Possible beam types defined in @ref BeamType
   * ESDs have it directly, AODs get it from hardcoded run number ranges
   * @return Beam type of the run.
   */
  BeamType                    GetBeamType()                                      const;

  /**
   * @brief Get \f$ i^{th} \f$ TClonesArray with particles
   * @param[in] i Index of the particle container inside the task
   * @return TClonesArray connected to the container (NULL if not found)
   */
  TClonesArray               *GetParticleArray(Int_t i=0)                        const;

  /**
   * @brief Get \f$ i^{th} \f$ TClonesArray with EMCAL clusters
   * @param[in] i Index of the cluster container inside the task
   * @return TClonesArray connected to the container (NULL if not found)
   */
  TClonesArray               *GetClusterArray(Int_t i=0)                         const;

  /**
   * @brief Get number of particles in container attached to this task with index i
   * @param[in] i Index of then particle container
   * @return Number of particles in container
   */
  Int_t                       GetNParticles(Int_t i=0)                           const;

  /**
   * @brief Get number of clusters in the cluster container attached to this task with
   * index i
   * @param[in] i Index of the container
   * @return Number of clusters inside the container
   */
  Int_t                       GetNClusters(Int_t i=0)                            const;

  /**
   * @brief Get main trigger match
   *
   * For the selection of the main trigger patch, high and low threshold triggers of a given category are grouped
   * If there are more than 1 main patch of a given trigger category (i.e. different high and low threshold patches),
   * the highest one according to the ADC value is taken. In case doSimpleOffline is true, then only the patches from
   * the simple offline trigger are used.
   * @param[in] triggersel Type of the EMCAL Level0/Level1 trigger
   * @param[in] doSimpleOffline If true simple offline patches are used
   * @return Main trigger patch for the trigger category (defined as highest energy patch)
   */
  AliEMCALTriggerPatchInfo   *GetMainTriggerPatch(TriggerCategory triggersel = kTriggerLevel1Jet, Bool_t doSimpleOffline = kFALSE);

  /**
   * @brief Check if event has a given trigger type
   *
   * Trigger types defined in @ref TriggerType.
   * @param triggersel Trigger type to check
   * @return True fo the trigger type is found in the event,
   * false otherwise.
   */
  Bool_t		                  HasTriggerType(TriggerType triggersel);

  /**
   * @brief Get list of selected triggers of the given event.
   *
   * The trigger selection is performed using EMCAL
   * patches found for the given event
   * @return A bitmap encoding the trigger found for the event
   */
  ULong_t 		                GetTriggerList();

  /**
   * @brief Loading PYTHIA information from external cross section file into the task
   *
   * Get the cross section and the trails either from pyxsec.root or from pysec_hists.root
   * Get the pt hard bin from the file path
   * This is to called in Notify and should provide the path to the AOD/ESD file
   * (Partially copied from AliAnalysisHelperJetTasks)
   * @param[in] currFile Name of the current ESD/AOD file
   * @param[out] fXsec Cross section calculated by PYTHIA
   * @param[out] fTrials Number of trials needed by PYTHIA
   * @param[out] pthard \f$ p_{t} \f$-hard bin, extracted from path name
   * @return True if parameters were obtained successfully, false otherwise
   */
  Bool_t                      PythiaInfoFromFile(const char* currFile, Float_t &fXsec, Float_t &fTrials, Int_t &pthard);
  /**
   * @brief Determines if a track is inside the EMCal acceptance.
   *
   * Selection is done using \f$\eta\f$/\f$\phi\f$ at the vertex (no propagation).
   * Includes +/- edges. Useful to determine whether track propagation should be attempted.
   * @param[in] part Particle to check
   * @param[in] edges Size of the edges in \f$\phi\f$ excluded from the EMCAL acceptance
   * @return True if a particle is inside the EMCAL acceptance, false otherwise
   */
  Bool_t                      IsTrackInEmcalAcceptance(AliVParticle* part, Double_t edges=0.9) const;

  /**
   * @brief Copy some information about the Pythia event in a PythaInfo object
   * @param[in] mcEvent Monte Carlo event from which the information is obtained
   */
  void                        GeneratePythiaInfoObject(AliMCEvent* mcEvent);

  /**
   * @brief Filter the mc tails in pt-hard distributions
   *
   * See https://twiki.cern.ch/twiki/bin/view/ALICE/JetMCProductionsCrossSections#How_to_reject_tails_in_the_pT_ha
   * @return kTRUE if it is not a MC outlier
   */
  virtual Bool_t              CheckMCOutliers();

  /**
   * @brief Main initialization function on the worker
   *
   * Performing run-independent initialization. This consists of
   * - Determining data type (ESD/AOD)
   * - Creating general QA histograms
   *
   * Attention: Histograms are only created in case the task is
   * configured for this (second argument in the named constructor).
   * In this case the container fOuput is created which can be used
   * by the users to handle and store their histograms. In this case
   * the users must overwrite this function in their tasks and call
   * this function right at the beginning of their function.
   *
   * The general QA histograms monitor event related observables like
   * the z-position of the primary vertex before and after event selection,
   * the trigger classes selecting the event and the event rejection
   * reason, but also Monte-Carlo related observables like the cross
   * section, the number of trials and the \f$ p_{t} \f$-hard bin in
   * case of a corresponding production.
   *
   * Important for user tasks: User tasks need to implement the function
   * UserCreateOutputObjects themselves. In case the output objects uses
   * the common output container fOutput, the function UserCreateOutputOjects
   * from AliAnalysisTaskEmcal needs to be called at the beginning. A
   * typical function in the user task could look like this:
   *
   * ~~~{.cxx}
   * void AliAnalysisTaskInheritingFromEmcal::UserCreateOutputObjects() {
   *   // Call UserCreateOutputObjects from AliAnalysisTaskEmcal in
   *   // order to initialize the common output container and the
   *   // general histograms
   *   AliAnalysisTaskEmcal::UserCreateOutputObjects();
   *
   *   fTestHisto = new TH1F("fTestHisto", "Test Histogram", 100, 0., 100);
   *   fOutput->Add(fTestHisto);
   * }
   * ~~~
   */
  void                        UserCreateOutputObjects();

  /**
   * @brief Event loop, called for each event.
   *
   * The function consists of three
   * steps:
   * -# Event selection
   * -# Running the user code
   * -# Filling general (QA) histograms
   * The event selection steps are documented in the function IsEventSelected.
   *
   * Users must not overwrite this function. Instead the virtual function Run
   * should be user and implemented by the user. The return value of the Run
   * function decides on whether general histograms are filled.
   *
   * In case the task is not yet initialized, which is the case for the first
   * event, the UserExec performs several basic initialization steps, documented
   * in the functions ExecOnce. Note that this is only done for the first event
   * and only for properties which need the presence of an input event.
   *
   * Function also steers events triggered by signals
   * - if a run changed
   * - if a file changed
   * Users can implement own actions when overloading the functions UserRunChanged and
   * UserFileChanged.
   *
   * @param[in] option Not used
   */
  void                        UserExec(Option_t *option);

  /**
   * @brief Notifying the user that the input data file has
   * changed and performing steps needed to be done.
   *
   * Only set the signal that the file has changed. The steps
   * for file changed are handled by the function FileChanged,
   * called by UserExec once it sees the signal fFileChanged.
   *
   * @return Always true
   */
  Bool_t                      UserNotify();

  /**
   * @brief  Steps to be executed when a few file is loaded into the
   * input handler
   *
   * The method is called in UserExec once it sees the FileChanged signal,
   * which is set in UserNotify() (as this handled asynchronously). As it
   * is called within UserExec one can expect that the event is fully initialized.
   *
   * This function is of relevance for analysis of  Monte-Carlo productions in
   * \f$ p_{t} \f$-hard bins as it reads the pythia cross section and the number
   * of  trials from the file pyxsec.root and fills the relevant distributions with
   * the values obtained.
   *
   * Function UserFileChanged is provided in order to allow for user code in
   * derived classes to be executed when the file changed.
   *
   * @return False if the data tree or the data file
   * doesn't exist, true otherwise
   */
  Bool_t					            FileChanged();

  /**
   * @brief Perform steps needed to initialize the analysis.
   *
   * This function relies on the presence of an input
   * event (ESD or AOD event). Consequently it is called
   * internally by UserExec for the first event.
   *
   * This function connects all containers attached to
   * this task to the corresponding arrays in the
   * input event. Furthermore it initializes the geometry.
   */
  virtual void                ExecOnce();

  /**
   * @brief Filling general histograms.
   *
   * Among the general histograms that are filled only in case of
   * running over MC productions are
   * - \f$ p_{t} \f$-hard bin
   * - Cross section after event selection
   * - Number of trials after event selection
   * - Number of events after event selection
   * In any case the vertex distribution is filled as general
   * histograms. For heavy ion collisions also the centrality
   * distribution and the event plane distribution are filled.
   * @return Always true
   */
  virtual Bool_t              FillGeneralHistograms();

  /**
   * @brief Performing event selection.
   *
   * By default the event selection is delegated to the 
   * AliEventCuts. Users can specify to use a builtin 
   * event select (old method) by calling SetUseInternalEventSelection.
   *
   * The Run function is only called in case the event
   * is selected by this function, providing a default
   * event selection. In case the user wants to apply
   * a different event selection, the function should
   * be implemented by the user.
   *
   * @return True if the event is selected.
   */
  virtual Bool_t              IsEventSelected();

  /**
   * @brief Selection of a hardware trigger
   * 
   * The function is used in the default implementation of 
   * IsEventSelected in order to perform the trigger selection.
   * Users can reimplement the function in order to perform the
   * trigger selection of their choise. The default implementation
   * checks for the trigger bits specified in the task configuration.
   * 
   * @return True if the event is selected as triggered event, false otherwise 
   */
  virtual Bool_t              IsTriggerSelected();

  /**
   * @brief Retrieve common objects from event.
   *
   * Several object used for the analysis are handled by the
   * AliAnalysisTaskEmcal, and provided as member variables to
   * the task itself and to derived task. Among them are
   * event global objects
   * - beam type
   * - EMCAL cells
   * - centrality percentile
   * - event plane angle
   * - triggers
   * - vertex position (both from global tracks and from SPD tracklets)
   *
   * as well as the arrays for
   * - EMCAL clusters
   * - particles
   *
   * In case of MC simulation the MC headers are loaded as well and in
   * case the task is flagged to run over a pythia or herwig production
   *
   * @return Always true
   */
  virtual Bool_t              RetrieveEventObjects();

  /**
   * @brief Process tasks relevant when a file with a different run number is processed
   *
   * Method exclusively called when the run is changed (new run number differing
   * from old run number). Can be used for run-dependent initializations (i.e.
   * setting parameters from the OADB)
   */
  virtual void                RunChanged(Int_t /*newrun*/)          {}

  /**
   * @brief Task initializations handled in user tasks
   *
   * Interface for user code executed when the first event is called.
   * At this step we know run number and data type and can therefore
   * do proper initializations.
   */
  virtual void                UserExecOnce()                    {}

  /**
   * @brief Virtual method for user code to be executed when a file changed
   *
   * The event is expected to be fully configured at that stage
   */
  virtual void				        UserFileChanged()					{}

  /**
   * @brief Function filling histograms
   *
   * This function optionally fills histograms created by the users. Can
   * access data previously handled by the user Run function.
   * @return
   */
  virtual Bool_t              FillHistograms()                  { return kTRUE                 ; }

  /**
   * @brief Run function. This is the core function of the analysis and
   * contains the user code. Therefore users have to implement this
   * function.
   *
   * @return True if event is selected, false otherwise
   */
  virtual Bool_t              Run()                             { return kTRUE                 ; }
    
  // Static utilities
  /**
   * @brief Calculate \f$\phi\f$ and \f$\eta\f$ difference between a track (t) and a cluster (c).
   *
   * The position of the track is obtained on the EMCAL surface
   * @param[in] t Track to check
   * @param[in] v Cluster to check
   * @param[out] phidiff Distance in \f$\phi\f$ between cluster and track
   * @param[out] etadiff Distance in \f$\eta\f$ between cluster and track
   */
  static void                 GetEtaPhiDiff(const AliVTrack *t, const AliVCluster *v, Double_t &phidiff, Double_t &etadiff);

  /**
   * @brief Get track type encoded from bits 20 and 21.
   * @param[in] t Track to check
   * @return
   */
  static Byte_t               GetTrackType(const AliVTrack *t);

  /**
   * @brief Decode track type
   *
   * Return track type: 0 = filterBit1, 1 = filterBit2 && ITS, 2 = filterBit2 && !ITS.
   * Returns 3 if filterBit1 and filterBit2 do not test.
   * WARNING: only works with AOD tracks and AOD filter bits must be provided. Otherwise will always return 0.
   * @param[in] aodTrack Track to check
   * @param filterBit1 First filter bit
   * @param filterBit2 Second filter bit
   * @return See explanation
   */
  static Byte_t               GetTrackType(const AliAODTrack *aodTrack, UInt_t filterBit1, UInt_t filterBit2);

  /**
   * Calculate Delta Phi.
   * @param[in] phia \f$ \phi \f$ of the first particle
   * @param[in] phib \f$ \phi \f$ of the second particle
   * @param[in] rangeMin Minimum \f$ \phi \f$ range
   * @param[in] rangeMax Maximum \f$ \phi \f$ range
   * @return Difference in \f$ \phi \f$
   */
  static Double_t             DeltaPhi(Double_t phia, Double_t phib, Double_t rMin = -TMath::Pi()/2, Double_t rMax = 3*TMath::Pi()/2);

  /**
   * Generate array with fixed binning within min and max with n bins. The array containing the bin
   * edges set will be created by this function. Attention, this function does not take care about
   * memory it allocates - the array needs to be deleted outside of this function
   * @param[in] n Number of bins
   * @param[in] min Minimum value for the binning
   * @param[in] max Maximum value for the binning
   * @return Array containing the bin edges created bu this function
   */
  static Double_t*            GenerateFixedBinArray(Int_t n, Double_t min, Double_t max);

  /**
   * Generate array with fixed binning within min and max with n bins. The parameter array
   * will contain the bin edges set by this function. Attention, the array needs to be
   * provided from outside with a size of n+1
   * @param[in] n Number of bins
   * @param[in] min Minimum value for the binning
   * @param[in] max Maximum value for the binning
   * @param[out] array Array containing the bin edges
   */
  static void                 GenerateFixedBinArray(Int_t n, Double_t min, Double_t max, Double_t* array);

  /**
   * @brief Perform event selection (old method)
   * 
   * This contains
   * - Selection of the trigger class
   * - Selection according to the centrality class
   * - Selection of event with good vertex quality
   * - Selection of the event plane orientation
   * - Selection of the multiplicity (including
   *   above minimum \f$ p_{t} \f$ and tracks in the
   *   EMCAL acceptance
   *
   * Note that for the vertex selection both the usage
   * of the analysis util and the range of the z-position
   * of the primary vertex need to be specified.
   *
   * In case the event is rejected, a histogram
   * monitoring the rejection reason is filled with
   * the bin corresponding to the source of the rejection
   * of the current event.
   * 
   * @return True if the event is selected, false if it is rejected
   */
  Bool_t                      IsEventSelectedInternal();

  /**
   * @brief Calculates the fraction of momentum z of part 1 w.r.t. part 2 in the direction of part 2.
   * @param[in] part1 Momentum vector for which the relative fraction is calculated
   * @param[in] part2 Reference momentum vector for the calculation
   * @return Relative fraction of momentum of particle 1 with respect to particle 2
   */
  static Double_t             GetParallelFraction(AliVParticle* part1, AliVParticle* part2);

  /**
   * @brief Calculates the fraction of momentum z of vect 1 w.r.t. part 2 in the direction of part 2.
   * @param[in] vect1 Momentum vector for which the relative fraction is calculated
   * @param[in] part2 Reference momentum vector for the calculation
   * @return Relative fraction of momentum of particle 1 with respect to particle 2
   */
  static Double_t             GetParallelFraction(const TVector3& vect1, AliVParticle* part2);

  static Double_t             fgkEMCalDCalPhiDivide;       ///<  phi value used to distinguish between DCal and EMCal

  // Task configuration
  TString                     fPythiaInfoName;             ///< name of pythia info object
  BeamType                    fForceBeamType;              ///< forced beam type
  Bool_t                      fGeneralHistograms;          ///< whether or not it should fill some general histograms
  Bool_t                      fLocalInitialized;           ///< whether or not the task has been already initialized
  Bool_t					            fFileChanged;				         //!<! Signal triggered when the file has changed
  Bool_t                      fCreateHisto;                ///< whether or not create histograms
  TString                     fCaloCellsName;              ///< name of calo cell collection
  TString                     fCaloTriggersName;           ///< name of calo triggers collection
  TString                     fCaloTriggerPatchInfoName;   ///< trigger patch info array name
  Double_t                    fMinCent;                    ///< min centrality for event selection
  Double_t                    fMaxCent;                    ///< max centrality for event selection
  Double_t                    fMinVz;                      ///< min vertex for event selection
  Double_t                    fMaxVz;                      ///< max vertex for event selection
  Int_t                       fMinVertexContrib;           ///< Min. number of vertex contributors
  Double_t                    fTrackPtCut;                 ///< cut on track pt in event selection
  Int_t                       fMinNTrack;                  ///< minimum nr of tracks in event with pT>fTrackPtCut
  Double_t                    fZvertexDiff;                ///< upper limit for distance between primary and SPD vertex
  Bool_t                      fUseAliAnaUtils;             ///< used for LHC13* data: z-vtx, Ncontributors, z-vtx resolution cuts
  Bool_t                      fRejectPileup;               ///< Reject pilup using function AliAnalysisUtils::IsPileUpEvent()
  Bool_t                      fTklVsClusSPDCut;            ///< Apply tracklet-vs-cluster SPD cut to reject background events in pp
  UInt_t                      fOffTrigger;                 ///< offline trigger for event selection
  TString                     fTrigClass;                  ///< trigger class name for event selection
  TString                     fMinBiasRefTrigger;          ///< Name of the minmum bias reference trigger, used in the calculation of downscale-corrected event numbers
  TriggerType                 fTriggerTypeSel;             ///< trigger type to select based on trigger patches
  Int_t                       fNbins;                      ///< no. of pt bins
  Double_t                    fMinBinPt;                   ///< min pt in histograms
  Double_t                    fMaxBinPt;                   ///< max pt in histograms
  Double_t                    fMinPtTrackInEmcal;          ///< min pt track in emcal
  Double_t                    fEventPlaneVsEmcal;          ///< select events which have a certain event plane wrt the emcal
  Double_t                    fMinEventPlane;              ///< minimum event plane value
  Double_t                    fMaxEventPlane;              ///< maximum event plane value
  TString                     fCentEst;                    ///< name of V0 centrality estimator
  Bool_t                      fIsEmbedded;                 ///< trigger, embedded signal
  Bool_t                      fIsPythia;                   ///< trigger, if it is a PYTHIA production
  Bool_t                      fIsHerwig;                   ///< trigger, if it is a HERWIG production
  Bool_t                      fGetPtHardBinFromName;       ///< Obtain pt-hard bin from file path
  Int_t                       fSelectPtHardBin;            ///< select one pt hard bin for analysis
  Int_t                       fMinMCLabel;                 ///< minimum MC label value for the tracks/clusters being considered MC particles
  Int_t                       fMCLabelShift;               ///< if MC label > fMCLabelShift, MC label -= fMCLabelShift
  Int_t                       fNcentBins;                  ///< how many centrality bins
  Bool_t                      fNeedEmcalGeom;              ///< whether or not the task needs the emcal geometry
  TObjArray                   fParticleCollArray;          ///< particle/track collection array
  TObjArray                   fClusterCollArray;           ///< cluster collection array
  ULong_t                     fTriggers;                   ///< list of fired triggers
  EMCalTriggerMode_t          fEMCalTriggerMode;           ///< EMCal trigger selection mode
  Bool_t                      fUseNewCentralityEstimation; ///< Use new centrality estimation (for 2015 data)
  Bool_t                      fGeneratePythiaInfoObject;   ///< Generate Pythia info object
  Bool_t                      fUsePtHardBinScaling;        ///< Use \f$ p_{t}\f$-hard bin scaling in merging
  Bool_t                      fUseXsecFromHeader;          //!<! Use cross section from header instead of pyxsec.root (purely transient)
  Bool_t                      fMCRejectFilter;             ///< enable the filtering of events by tail rejection
  Bool_t                      fCountDownscaleCorrectedEvents; ///< Count event number corrected for downscaling
  Bool_t                      fUseBuiltinEventSelection;   ///< Use builtin event selection of the AliAnalysisTaskEmcal instead of AliEventCuts
  Float_t                     fPtHardAndJetPtFactor;       ///< Factor between ptHard and jet pT to reject/accept event.
  Float_t                     fPtHardAndClusterPtFactor;   ///< Factor between ptHard and cluster pT to reject/accept event.
  Float_t                     fPtHardAndTrackPtFactor;     ///< Factor between ptHard and track pT to reject/accept event.

  // Service fields
  Int_t                       fRunNumber;                  //!<!run number (triggering RunChanged())
  AliEventCuts                fAliEventCuts;               ///< Event cuts (run2 defaults)
  AliAnalysisUtils           *fAliAnalysisUtils;           //!<!vertex selection (optional)
  Bool_t                      fIsEsd;                      //!<!whether it's an ESD analysis
  AliEMCALGeometry           *fGeom;                       //!<!emcal geometry
  TClonesArray               *fTracks;                     //!<!tracks
  TClonesArray               *fCaloClusters;               //!<!clusters
  AliVCaloCells              *fCaloCells;                  //!<!cells
  AliVCaloTrigger            *fCaloTriggers;               //!<!calo triggers
  TClonesArray               *fTriggerPatchInfo;           //!<!trigger patch info array
  Double_t                    fCent;                       //!<!event centrality
  Int_t                       fCentBin;                    //!<!event centrality bin
  Double_t                    fEPV0;                       //!<!event plane V0
  Double_t                    fEPV0A;                      //!<!event plane V0A
  Double_t                    fEPV0C;                      //!<!event plane V0C
  Double_t                    fVertex[3];                  //!<!event vertex
  Double_t                    fVertexSPD[3];               //!<!event Svertex
  Int_t                       fNVertCont;                  //!<!event vertex number of contributors
  Int_t                       fNVertSPDCont;               //!<!event SPD vertex number of contributors
  BeamType                    fBeamType;                   //!<!event beam type
  AliGenPythiaEventHeader    *fPythiaHeader;               //!<!event Pythia header
  AliGenHerwigEventHeader    *fHerwigHeader;               //!<!event Herwig header
  Float_t                     fPtHard;                     //!<!event \f$ p_{t}\f$-hard
  Int_t                       fPtHardBin;                  //!<!event \f$ p_{t}\f$-hard bin
  Int_t                       fPtHardBinGlobal;            //!<!event \f$ p_{t}\f$-hard bin, detected from filename
  Bool_t                      fPtHardInitialized;          //!<!flag whether the \f$ p_{t}\f$-hard bin was initialized, purely for internal processing
  Int_t                       fNPtHardBins;                ///< Number of \f$ p_{t}\f$-hard bins in the dataset
  TArrayI                     fPtHardBinning;              ///< \f$ p_{t}\f$-hard binning
  Int_t                       fNTrials;                    //!<!event trials
  Float_t                     fXsection;                   //!<!x-section from pythia header
  AliEmcalPythiaInfo         *fPythiaInfo;                 //!<!event parton info

  // Output
  AliEmcalList               *fOutput;                     //!<!output list
  TH1                        *fHistEventCount;             //!<!incoming and selected events
  TH1                        *fHistTrialsAfterSel;         //!<!total number of trials per pt hard bin after selection
  TH1                        *fHistEventsAfterSel;         //!<!total number of events per pt hard bin after selection
  TProfile                   *fHistXsectionAfterSel;       //!<!x section from pythia header
  TH1                        *fHistTrials;                 //!<!trials from pyxsec.root
  TH1                        *fHistEvents;                 //!<!total number of events per pt hard bin
  TProfile                   *fHistXsection;               //!<!x section from pyxsec.root
  TH1                        *fHistPtHard;                 //!<!\f$ p_{t}\f$-hard distribution
  TH2                        *fHistPtHardCorr;             //!<!Correlation between \f$ p_{t}\f$-hard value and bin
  TH2                        *fHistPtHardCorrGlobal;       //!<!Correlation between \f$ p_{t}\f$-hard value and global bin
  TH2                        *fHistPtHardBinCorr;          //!<!Correlation between global and local (per-event) \f$ p_{t}\f$-hard bin
  TH1                        *fHistCentrality;             //!<!event centrality distribution
  TH1                        *fHistZVertex;                //!<!z vertex position
  TH1                        *fHistEventPlane;             //!<!event plane distribution
  TH1                        *fHistEventRejection;         //!<!book keep reasons for rejecting event
  TH1                        *fHistTriggerClasses;         //!<!number of events in each trigger class
  TH1                        *fHistTriggerClassesCorr;     //!<!corrected number of events in each trigger class

 private:
  AliAnalysisTaskEmcal(const AliAnalysisTaskEmcal&);            // not implemented
  AliAnalysisTaskEmcal &operator=(const AliAnalysisTaskEmcal&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcal, 20) // EMCAL base analysis task
  /// \endcond
};

inline Double_t AliAnalysisTaskEmcal::DeltaPhi(Double_t phia, Double_t phib, Double_t rangeMin, Double_t rangeMax) 
{
  Double_t dphi = -999;
  const Double_t tpi = TMath::TwoPi();
  
  if (phia < 0)         phia += tpi;
  else if (phia > tpi) phia -= tpi;
  if (phib < 0)         phib += tpi;
  else if (phib > tpi) phib -= tpi;
  dphi = phib - phia;
  if (dphi < rangeMin)      dphi += tpi;
  else if (dphi > rangeMax) dphi -= tpi;
  
  return dphi;
}

inline void AliAnalysisTaskEmcal::GenerateFixedBinArray(Int_t n, Double_t min, Double_t max, Double_t* array)
{
  Double_t binWidth = (max-min)/n;
  array[0] = min;
  for (Int_t i = 1; i <= n; i++) {
    array[i] = array[i-1]+binWidth;
  }
}

inline Double_t* AliAnalysisTaskEmcal::GenerateFixedBinArray(Int_t n, Double_t min, Double_t max)
{
  Double_t *array = new Double_t[n+1];
  GenerateFixedBinArray(n, min, max, array);
  return array;
}

#endif
