#ifndef ALIANALYSISTASKEMCALLIGHT_H
#define ALIANALYSISTASKEMCALLIGHT_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class TClonesArray;
class TString;
class TList;
class AliEmcalParticle;
class AliMCParticle;
class AliVCluster;
class AliVTrack;
class AliVParticle;
class AliVCaloCells;
class TH1;
class TProfile;
class AliEMCALGeometry;
class AliGenEventHeader;
class AliGenPythiaEventHeader;
class AliVCaloTrigger;
class AliAnalysisUtils;
class AliEMCALTriggerPatchInfo;
class AliAODTrack;


#include <map>
#include <set>
#include <string>
#include <vector>
#include <iostream>

#include "Rtypes.h"

#include "AliEventCuts.h"
#include "AliParticleContainer.h"
#include "AliMCParticleContainer.h"
#include "AliTrackContainer.h"
#include "AliClusterContainer.h"
#include "AliEmcalStringView.h"


#include "AliAnalysisTaskSE.h"
/**
 * @class AliAnalysisTaskEmcalLight
 * @brief Base task in the EMCAL framework (lighter version of AliAnalysisTaskEmcal)
 * @ingroup EMCALCOREFW
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 *
 * This class is the base class for Analysis Tasks using the
 * core EMCAL framework. User tasks can choose to inherit from it
 * as an alternative to AliAnalysisTaskEmcal.
 * In contrast to the normal AliAnalysisTaskSE, the main event
 * loop function to be implemented by the user is called Run.
 * This function is only called in case the event was selected
 * previously.
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
 * of cluster- or particle-containers are not mixed. Containers
 * provide an easy access to content created by other tasks and
 * attached to the event as a TClonesArray.
 *
 * For more information refer to \subpage EMCALAnalysisTask
 */
class AliAnalysisTaskEmcalLight : public AliAnalysisTaskSE {
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
   * @enum EBeamType_t
   * @brief Switch for the beam type
   */
  enum EBeamType_t {
    kNA       = -1,//!< Undefined
    kpp       = 0, //!< Proton-Proton
    kAA       = 1, //!< Nucleus-Nucleus
    kpA       = 2  //!< Proton-Nucleus
  };

  /**
   * @enum ECentralityEstimation_t
   * @brief Switch for the centrality estimation
   */
  enum ECentralityEstimation_t {
    kNoCentrality  = 0, //!< No centrality estimation
    kNewCentrality = 1, //!< New centrality estimation (AliMultSelection, see https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AliMultSelectionCalibStatus for calibration status period-by-period)
    kOldCentrality = 2  //!< Old centrality estimation (AliCentrality, works only on Run-1 PbPb and pPb)
  };

  /**
   * @brief Default constructor.
   */
  AliAnalysisTaskEmcalLight();

  /**
   * Standard constructor. Should be used by the user.
   *
   * Note: This constructor also handles the general histograms. In
   * case the second parameter is true, then general histograms (see
   * UserCreateOutputObjects and FillHistograms) are created and filled
   * by the task, and a container is provided handling the user histograms.
   * @param[in] name Name of the task
   * @param[in] histo If true then general histograms are filled by the task
   */
  AliAnalysisTaskEmcalLight(const char *name, Bool_t histo=kFALSE);

  /**
   * Destructor
   */
  virtual ~AliAnalysisTaskEmcalLight();

  // Containers
  /**
   * Create new particle container and attach it to the task. The name
   * provided to this function must match the name of the array attached
   * to the new container inside the input event.
   * @param[in] branchName Name of the array the container points to
   * @param[in] contName Name of the container points to (optional)
   * @return Pointer to the new particle container
   */
  AliParticleContainer       *AddParticleContainer(EMCAL_STRINGVIEW branchName, EMCAL_STRINGVIEW contName="");
  /**
   * Create new cluster container and attach it to the task. The name
   * provided to this function must match the name of the array attached
   * to the new container inside the input event.
   * @param[in] branchName Name of the array the container points to
   * @param[in] contName Name of the container points to (optional)
   * @return Pointer to the new cluster container
   */
  AliClusterContainer        *AddClusterContainer(EMCAL_STRINGVIEW branchName, EMCAL_STRINGVIEW contName="");
  void                        AdoptParticleContainer(AliParticleContainer* cont)    { fParticleCollArray[cont->GetName()] = cont; }
  void                        AdoptClusterContainer(AliClusterContainer* cont)      { fClusterCollArray[cont->GetName()]  = cont; }

  /**
   * Find particle container attached to this task according to its name
   * @param[in] name Name of the particle container
   * @return Particle container found under the given name
   */
  AliParticleContainer       *GetParticleContainer(EMCAL_STRINGVIEW name)          const;

  /**
   * Find cluster container attached to this task according to its name
   * @param[in] name Name of the cluster container
   * @return Cluster container found under the given name
   */
  AliClusterContainer        *GetClusterContainer(EMCAL_STRINGVIEW name)           const;
  AliMCParticleContainer     *GetMCParticleContainer(EMCAL_STRINGVIEW name)        const { return dynamic_cast<AliMCParticleContainer*>(GetParticleContainer(name)); }
  AliTrackContainer          *GetTrackContainer(EMCAL_STRINGVIEW name)             const { return dynamic_cast<AliTrackContainer*>(GetParticleContainer(name))     ; }
  void                        RemoveParticleContainer(const std::string& name)           { fParticleCollArray.erase(name)                   ; }
  void                        RemoveClusterContainer(const std::string& name)            { fClusterCollArray.erase(name)                       ; }

  // Other input data
  void                        SetCaloCellsName(const char *n)                       { fCaloCellsName     = n                              ; }
  void                        SetCaloTriggerPatchInfoName(const char *n)            { fCaloTriggerPatchInfoName = n                       ; }
  void                        SetCaloTriggersName(const char *n)                    { fCaloTriggersName  = n                              ; }
  void                        SetCentralityEstimator(const char *c)                 { fCentEst           = c                              ; }
  void                        SetIsMonteCarlo(Bool_t i)                             { fIsMonteCarlo      = i                              ; }
  void                        SetIsPythia(Bool_t i);
  void                        SetMCEventHeaderName(const char* name);
  void                        SetForceBeamType(EBeamType_t f)                       { fForceBeamType     = f                              ; }

  // Task configuration
  void                        SetCentralityEstimation(ECentralityEstimation_t b)    { fCentralityEstimation = b                     ; }
  void                        SetMakeGeneralHistograms(Bool_t g)                    { fGeneralHistograms = g                              ; }
  void                        SetNeedEmcalGeom(Bool_t n)                            { fNeedEmcalGeom     = n                              ; }
  void                        SetCentBins(const std::vector<double>& bins)          { fCentBins = std::vector<double>(bins)               ; }
  Int_t                       GetNCentBins()                                  const { return fCentBins.size() > 1 ? fCentBins.size() - 1 : 1; }
  void                        SetSwitchOffLHC15oFaultyBranches(Bool_t b)            { fSwitchOffLHC15oFaultyBranches = b                  ; }

  // Event selection
  void                        SetWarnMissingCentrality(Bool_t doWarn)               { fWarnMissingCentrality = doWarn                     ; }
  void                        SetTriggerSelectionBitMap(UInt_t t)                   { fTriggerSelectionBitMap = t                         ; }
  void                        SetCentRange(Double_t min, Double_t max)              { fMinCent           = min  ; fMaxCent     = max      ; }
  void                        SetVzRange(Double_t min, Double_t max)                { fMinVz             = min  ; fMaxVz       = max      ; }
  void                        SetZvertexDiffValue(Double_t cut)                     { fMaxVzDiff         = cut                            ; }
  void                        SetMinNVertCont(Int_t cut)                            { fMinNVertCont      = cut                            ; }
  void                        SetPtHardRange(Double_t min, Double_t max)            { fMinPtHard         = min  ; fMaxPtHard   = max      ; }
  void                        SetMaxMinimumBiasPtHard(Int_t max)                    { fMaxMinimumBiasPtHard = max                         ; }
  void                        AddAcceptedTriggerClass(const char* trigClass)        { fAcceptedTriggerClasses.insert(trigClass)           ; }
  void                        AddRejectedTriggerClass(const char* trigClass)        { fRejectedTriggerClasses.insert(trigClass)           ; }
  void                        ClearAcceptedTriggerClasses()                         { fAcceptedTriggerClasses.clear()                     ; }
  void                        ClearRejectedTriggerClasses()                         { fRejectedTriggerClasses.clear()                     ; }
  void                        SetMCFilter()                                         { fMCRejectFilter = kTRUE                             ; }
  void                        ResetMCFilter()                                       { fMCRejectFilter = kFALSE                            ; }
  void                        SetJetPtFactor(Float_t f)                             { fPtHardAndJetPtFactor = f                           ; }
  Float_t                     JetPtFactor()                                         { return fPtHardAndJetPtFactor                        ; }
  void                        SetClusterPtFactor(Float_t f)                         { fPtHardAndClusterPtFactor = f                       ; }
  Float_t                     ClusterPtFactor()                                     { return fPtHardAndClusterPtFactor                    ; }
  void                        SetTrackPtFactor(Float_t f)                           { fPtHardAndTrackPtFactor = f                         ; }
  Float_t                     TrackPtFactor()                                       { return fPtHardAndTrackPtFactor                      ; }
  void                        SetEventSelectionAfterRun(Bool_t b)                   { fEventSelectionAfterRun = b                         ; }
  void                        SelectGeneratorName(TString gen)                      { fSelectGeneratorName = gen                          ; }
  void                        SetInhibit(Bool_t s)                                  { fInhibit = s                                        ; }
  void                        SetEventWeightRange(Double_t min, Double_t max)       { fMinimumEventWeight = min; fMaximumEventWeight = max; }
  void                        SetUseAliEmcalList(Bool_t doUse)                      { fUseAliEmcalList = doUse                            ; }
  void                        SetUsePtHardBinScaling(Bool_t b)                      { fUsePtHardBinScaling = b                            ; }

  /**
   * @brief Use internal (old) event selection
   * @param[in] doUse It true use the old internal event selection instead of AliEventCuts
   */
  void                        SetUseBuiltinEventSelection(Bool_t doUse)            { fUseBuiltinEventSelection = doUse                  ; }
  
  AliEventCuts               &GetEventCuts()                                        { return fAliEventCuts; }

  Bool_t IsInhibit() const { return fInhibit; }

 protected:
  void                        SetRejectionReasonLabels(TAxis* axis);

  /**
   * Add object to event
   * @param[in] obj Object to be added
   * @param[in] attempt If true don't handle error
   */
  void                        AddObjectToEvent(TObject *obj, Bool_t attempt = kFALSE);

  /**
   * Read a TClonesArray from event. Attention: Both the
   * name of the array and the name of the object stored inside
   * must match.
   * @param[in] name Name of the array to be read in
   * @param[in] clname Name of the type of the objects stored in the array
   * @return Pointer to the TClonesArray (NULL if not found)
   */
  TClonesArray               *GetArrayFromEvent(const char *name, const char *clname=0);

  /**
   * Get beam type : pp-AA-pA
   * ESDs have it directly, AODs get it from hardcoded run number ranges
   * @return Beam type of the run.
   */
  EBeamType_t                 GetBeamType();

  /**
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
  Bool_t                      PythiaInfoFromFile(const char* currFile, Float_t &fXsec, Float_t &fTrials, Int_t &pthard, Bool_t &useXsecFromHeader);

  /**
   * Determines if a track is inside the EMCal acceptance, using \f$\eta\f$/\f$\phi\f$ at the vertex (no propagation).
   * Includes +/- edges. Useful to determine whether track propagation should be attempted.
   * @param[in] part Particle to check
   * @param[in] edges Size of the edges in \f$\phi\f$ excluded from the EMCAL acceptance
   * @return True if a particle is inside the EMCAL acceptance, false otherwise
   */
  Bool_t                      IsTrackInEmcalAcceptance(AliVParticle* part, Double_t edges=0.9) const;

  /**
   * Filter the mc tails in pt-hard distributions
   * See https://twiki.cern.ch/twiki/bin/view/ALICE/JetMCProductionsCrossSections#How_to_reject_tails_in_the_pT_ha
   * @return kTRUE if it is not a MC outlier
   */
  Bool_t                      CheckMCOutliers();

  // Overloaded AliAnalysisTaskSE methods

  /**
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
   */
  void                        UserCreateOutputObjects();

  /**
   * Event loop, called for each event. The function consists of three
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
   * event, the UserExec performs several basic initilization steps, documented
   * in the functions ExecOnce. Note that this is only done for the first event
   * and only for properties which need the presence of an input event.
   *
   * @param[in] option Not used
   */
  void                        UserExec(Option_t *option);

  /**
   * Notifying the user that the input data file has
   * changed and performing steps needed to be done.
   *
   * This function is of relevance for analysis of
   * Monte-Carlo productions in \f$ p_{t} \f$-hard
   * bins as it reads the pythia cross section and
   * the number of trials from the file pyxsec.root
   * and fills the relevant distributions with
   * the values obtained.
   * @return False if the data tree or the data file
   * doesn't exist, true otherwise
   */
  Bool_t                      UserNotify();

  // Virtual functions, to be overloaded in derived classes

  /**
   * Perform steps needed to initialize the analysis.
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
   * Filling general histrograms. Among the general histograms
   * that are filled only in case of running over MC productions
   * are
   * - \f$ p_{t} \f$-hard bin
   * - Cross section after event selection
   * - Number of trials after event selection
   * - Number of events after event selection
   * In any case the vertex distribution is filled as general
   * histograms. For heavy ion collisions also the centrality
   * distribution and the event plane distribution are filled.
   * @param[in] eventSelected flag that tells the method whether event selection has been performed already (a different set of histograms is filled)
   * @return Always true
   */
  virtual Bool_t              FillGeneralHistograms(Bool_t eventSelected);

  virtual Bool_t              IsEventSelected();

  /**
   * Performing event selection. This contains
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
   * monitoring the rejeciton reason is filled with
   * the bin corresponding to the source of the rejection
   * of the current event.
   *
   * @return True if the event is selected.
   */
  virtual Bool_t              IsEventSelectedInternal();

  /**
   * @brief Perform trigger selection
   * 
   * The default implementation checks if a trigger string is selected or rejected.
   * Users can overwrite the trigger selection method with a custom trigger selection.
   * 
   * @return True if the event is selected as matched trigger, false otherwise
   */
  virtual Bool_t              IsTriggerSelected();

  /**
   * Retrieve objects from event.
   * @return
   */
  virtual Bool_t              RetrieveEventObjects();
  /**
   * This function optionally fills histograms created by the users. Can
   * access data previously handled by the user Run function.
   * @return
   */
  virtual Bool_t              FillHistograms()                  { return kTRUE                 ; }
  /**
   * Run function. This is the core function of the analysis and
   * contains the user code. Therefore users have to implement this
   * function.
   *
   * @return True if event is selected, false otherwise
   */
  virtual Bool_t              Run()                             { return kTRUE                 ; }

  // Static utilities

  /**
   * Calculate \f$\phi\f$ and \f$\eta\f$ difference between a track (t) and a cluster (c). The
   * position of the track is obtained on the EMCAL surface
   * @param[in] t Track to check
   * @param[in] v Cluster to check
   * @param[out] phidiff Distance in \f$\phi\f$ between cluster and track
   * @param[out] etadiff Distance in \f$\eta\f$ between cluster and track
   */
  static void                 GetEtaPhiDiff(const AliVTrack *t, const AliVCluster *v, Double_t &phidiff, Double_t &etadiff);

  /**
   * Get track type encoded from bits 20 and 21.
   * @param[in] t Track to check
   * @return
   */
  static Byte_t               GetTrackType(const AliVTrack *t);

  /**
   * Return track type: 0 = filterBit1, 1 = filterBit2 && ITS, 2 = filterBit2 && !ITS.
   * Returns 3 if filterBit1 and filterBit2 do not test.
   * WARNING: only works with AOD tracks and AOD filter bits must be provided. Otherwise will always return 0.
   * @param aodTrack
   * @param filterBit1
   * @param filterBit2
   * @return
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
   * Generate array with fixed binning within min and max with n bins. The parameter array
   * will contain the bin edges set by this function. Attention, the array needs to be
   * provided from outside with a size of n+1
   * @param[in] n Number of bins
   * @param[in] min Minimum value for the binning
   * @param[in] max Maximum value for the binning
   * @param[out] array Vector where the bins are added
   */
  static std::vector<double>  GenerateFixedBinArray(int n, double min, double max, bool last = true);

  /**
   * Generate array with fixed binning within min and max with n bins. The array containing the bin
   * edges set will be created by this function. Attention, this function does not take care about
   * memory it allocates - the array needs to be deleted outside of this function
   * @param[in] n Number of bins
   * @param[in] min Minimum value for the binning
   * @param[in] max Maximum value for the binning
   * @return Vector containing the bin edges created by this function
   */
  static void                 GenerateFixedBinArray(int n, double min, double max, std::vector<double>& array, bool last = true);

  /**
   * Generate array with logaritmic fixed binning within min and max with n bins. The parameter array
   * will contain the bin edges set by this function. Attention, the array needs to be
   * provided from outside with a size of n+1
   * @param[in] n Number of bins
   * @param[in] min Minimum value for the binning
   * @param[in] max Maximum value for the binning
   * @param[out] array Vector where the bins are added
   */
  static std::vector<double>  GenerateLogFixedBinArray(int n, double min, double max, bool last = true);
  
  /**
   * Generate array with logaritmic fixed binning within min and max with n bins. The array containing the bin
   * edges set will be created by this function. Attention, this function does not take care about
   * memory it allocates - the array needs to be deleted outside of this function
   * @param[in] n Number of bins
   * @param[in] min Minimum value for the binning
   * @param[in] max Maximum value for the binning
   * @return Vector containing the bin edges created by this function
   */
  static void                 GenerateLogFixedBinArray(int n, double min, double max, std::vector<double>& array, bool last = true);

  /**
   * Calculates the fraction of momentum z of part 1 w.r.t. part 2 in the direction of part 2.
   * @param[in] part1 Momentum vector for which the relative fraction is calculated
   * @param[in] part2 Reference momentum vector for the calculation
   * @return Relative fraction of momentum of particle 1 with respect to particle 2
   */
  static Double_t             GetParallelFraction(AliVParticle* part1, AliVParticle* part2);

  /**
   * Calculates the fraction of momentum z of vect 1 w.r.t. part 2 in the direction of part 2.
   * @param[in] vect1 Momentum vector for which the relative fraction is calculated
   * @param[in] part2 Reference momentum vector for the calculation
   * @return Relative fraction of momentum of particle 1 with respect to particle 2
   */
  static Double_t             GetParallelFraction(const TVector3& vect1, AliVParticle* part2);

  /**
   * Determine the beam type based on hard-coded run ranges
   * \param runnumber run number
   * \return enumeration value corresponding to the beam type
   */
  static EBeamType_t          BeamTypeFromRunNumber(Int_t runnumber);

  static Double_t             fgkEMCalDCalPhiDivide;       ///<  phi value used to distinguish between DCal and EMCal

  // Task configuration
  EBeamType_t                 fForceBeamType;              ///< forced beam type
  Bool_t                      fGeneralHistograms;          ///< whether or not it should fill some general histograms
  Bool_t                      fCreateHisto;                ///< whether or not create histograms
  Bool_t                      fNeedEmcalGeom;              ///< whether or not the task needs the emcal geometry
  Bool_t                      fUseBuiltinEventSelection;   ///< use builtin (old) event selection
  std::vector<double>         fCentBins;                   ///< how many centrality bins
  ECentralityEstimation_t     fCentralityEstimation;       ///< Centrality estimation
  AliEventCuts                fAliEventCuts;               ///< Event cut object

  // Input data
  Bool_t                      fIsPythia;                   ///< if it is a PYTHIA production
  Bool_t                      fIsMonteCarlo;               ///< if it is a MC production
  TString                     fMCEventHeaderName;          ///< Looks for MC event properties in a particular MC event type (useful for a MC cocktail production)
  TString                     fCaloCellsName;              ///< name of calo cell collection
  TString                     fCaloTriggersName;           ///< name of calo triggers collection
  TString                     fCaloTriggerPatchInfoName;   ///< trigger patch info array name
  TString                     fCentEst;                    ///< name of the centrality estimator

  std::map<std::string,
  AliParticleContainer*>      fParticleCollArray;          ///< particle/track collection array
  std::map<std::string,
  AliClusterContainer*>       fClusterCollArray;          ///< cluster collection array

  // Event selection
  UInt_t                      fTriggerSelectionBitMap;     ///< trigger selection bit map
  Double_t                    fMinCent;                    ///< min centrality for event selection
  Double_t                    fMaxCent;                    ///< max centrality for event selection
  Double_t                    fMinVz;                      ///< min vertex for event selection
  Double_t                    fMaxVz;                      ///< max vertex for event selection
  Double_t                    fMaxVzDiff;                  ///< upper limit for distance between primary and SPD vertex
  Double_t                    fMinNVertCont;               ///< minumum number of vertex contributors
  Double_t                    fMinPtHard;                  ///< select minimum pt hard (MC)
  Double_t                    fMaxPtHard;                  ///< select maximum pt hard (MC)
  Double_t                    fMaxMinimumBiasPtHard;       ///< maximum pt hard for the minimum bias pt hard bin (MC)
  std::set<std::string>       fAcceptedTriggerClasses;     ///< list of accepted trigger classes
  std::set<std::string>       fRejectedTriggerClasses;     ///< list of accepted trigger classes
  Bool_t                      fMCRejectFilter;             ///< enable the filtering of events by tail rejection
  Float_t                     fPtHardAndJetPtFactor;       ///< Factor between ptHard and jet pT to reject/accept event.
  Float_t                     fPtHardAndClusterPtFactor;   ///< Factor between ptHard and cluster pT to reject/accept event.
  Float_t                     fPtHardAndTrackPtFactor;     ///< Factor between ptHard and track pT to reject/accept event.
  Bool_t                      fSwitchOffLHC15oFaultyBranches; ///< Switch off faulty tree branches in LHC15o AOD trees
  Bool_t                      fEventSelectionAfterRun;     ///< If kTRUE, the event selection is performed after Run() but before FillHistograms()
  Bool_t                      fUseAliEmcalList;            ///< Use AliEmcalList as output object
  Bool_t                      fUsePtHardBinScaling;        ///< Apply pt-hard bin scaling (in case AliEmcalList is used to handle the output)
  TString                     fSelectGeneratorName;        ///< Selects only events produced by a generator that has a name containing a string
  Double_t                    fMinimumEventWeight;         ///< Minimum event weight for the related bookkeping histogram
  Double_t                    fMaximumEventWeight;         ///< Minimum event weight for the related bookkeping histogram

  // Service fields
  Bool_t                      fInhibit;                    //!<!inhibit execution of the task
  Bool_t                      fLocalInitialized;           //!<!whether or not the task has been already initialized
  Bool_t                      fWarnMissingCentrality;      //!<!switch for verbosity in centrality information
  EDataType_t                 fDataType;                   //!<!data type (ESD or AOD)
  AliEMCALGeometry           *fGeom;                       //!<!emcal geometry
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
  ULong_t                     fFiredTriggerBitMap;         //!<!bit map of fired triggers
  std::vector<std::string>    fFiredTriggerClasses;        //!<!trigger classes fired by the current event
  EBeamType_t                 fBeamType;                   //!<!event beam type
  AliGenEventHeader          *fMCHeader;                   //!<!event MC header
  AliGenPythiaEventHeader    *fPythiaHeader;               //!<!event Pythia header
  Bool_t                      fUseXsecFromHeader;          //!<!Switch for using cross section from header (if not found in pythia file)
  Int_t                       fPtHardBin;                  //!<!event pt hard bin
  Double_t                    fPtHard;                     //!<!event pt hard
  Int_t                       fNTrials;                    //!<!event trials
  Float_t                     fXsection;                   //!<!x-section from pythia header
  Float_t                     fEventWeight;                //!<!event weight
  TString                     fGeneratorName;              //!<!name of the MC generator used to produce the current event (only AOD)

  // Output
  TList                      *fOutput;                     //!<!output list

 private:
  std::map<std::string, TH1*> fHistograms;                 //!<!general QA histograms
  TH1* GetGeneralTH1(const char* name, bool warn=false);
  TH2* GetGeneralTH2(const char* name, bool warn=false);
  TProfile* GetGeneralTProfile(const char* name, bool warn=false);

  AliAnalysisTaskEmcalLight(const AliAnalysisTaskEmcalLight&);            // not implemented
  AliAnalysisTaskEmcalLight &operator=(const AliAnalysisTaskEmcalLight&); // not implemented

  ClassDef(AliAnalysisTaskEmcalLight, 5);
};

#endif
