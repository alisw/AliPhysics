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

#include "AliParticleContainer.h"
#include "AliMCParticleContainer.h"
#include "AliTrackContainer.h"
#include "AliClusterContainer.h"

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

  AliAnalysisTaskEmcalLight();
  AliAnalysisTaskEmcalLight(const char *name, Bool_t histo=kFALSE);
  virtual ~AliAnalysisTaskEmcalLight();

  // Containers
  AliParticleContainer       *AddParticleContainer(std::string branchName, std::string contName="");
  AliClusterContainer        *AddClusterContainer(std::string branchName, std::string contName="");
  void                        AdoptParticleContainer(AliParticleContainer* cont)    { fParticleCollArray[cont->GetName()] = cont; }
  void                        AdoptClusterContainer(AliClusterContainer* cont)      { fClusterCollArray[cont->GetName()]  = cont; }
  AliParticleContainer       *GetParticleContainer(std::string name)          const;
  AliClusterContainer        *GetClusterContainer(std::string name)           const;
  AliMCParticleContainer     *GetMCParticleContainer(std::string name)        const { return dynamic_cast<AliMCParticleContainer*>(GetParticleContainer(name)); }
  AliTrackContainer          *GetTrackContainer(std::string name)             const { return dynamic_cast<AliTrackContainer*>(GetParticleContainer(name))     ; }
  void                        RemoveParticleContainer(std::string name)             { fParticleCollArray.erase(name)                      ; }
  void                        RemoveClusterContainer(std::string name)              { fClusterCollArray.erase(name)                       ; }

  // Other input data
  void                        SetCaloCellsName(const char *n)                       { fCaloCellsName     = n                              ; }
  void                        SetCaloTriggerPatchInfoName(const char *n)            { fCaloTriggerPatchInfoName = n                       ; }
  void                        SetCaloTriggersName(const char *n)                    { fCaloTriggersName  = n                              ; }
  void                        SetCentralityEstimator(const char *c)                 { fCentEst           = c                              ; }
  void                        SetIsPythia(Bool_t i)                                 { fIsPythia          = i                              ; }
  void                        SetForceBeamType(EBeamType_t f)                       { fForceBeamType     = f                              ; }

  // Task configuration
  void                        SetCentralityEstimation(ECentralityEstimation_t b)    { fCentralityEstimation = b                     ; }
  void                        SetMakeGeneralHistograms(Bool_t g)                    { fGeneralHistograms = g                              ; }
  void                        SetNeedEmcalGeom(Bool_t n)                            { fNeedEmcalGeom     = n                              ; }
  void                        SetCentBins(const std::vector<double>& bins)          { fCentBins = std::vector<double>(bins)               ; }
  Int_t                       GetNCentBins()                                  const { return fCentBins.size() > 1 ? fCentBins.size() - 1 : 1; }
  void                        SetSwitchOffLHC15oFaultyBranches(Bool_t b)            { fSwitchOffLHC15oFaultyBranches = b                  ; }

  // Event selection
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

  Bool_t IsInhibit() const { return fInhibit; }

 protected:
  void                        SetRejectionReasonLabels(TAxis* axis);
  void                        AddObjectToEvent(TObject *obj, Bool_t attempt = kFALSE);
  TClonesArray               *GetArrayFromEvent(const char *name, const char *clname=0);
  EBeamType_t                 GetBeamType();
  Bool_t                      PythiaInfoFromFile(const char* currFile, Float_t &fXsec, Float_t &fTrials, Int_t &pthard);
  Bool_t                      IsTrackInEmcalAcceptance(AliVParticle* part, Double_t edges=0.9) const;
  Bool_t                      CheckMCOutliers();

  // Overloaded AliAnalysisTaskSE methods
  void                        UserCreateOutputObjects();
  void                        UserExec(Option_t *option);
  Bool_t                      UserNotify();

  // Virtual functions, to be overloaded in derived classes
  virtual void                ExecOnce();
  virtual Bool_t              FillGeneralHistograms(Bool_t eventSelected);
  virtual Bool_t              IsEventSelected();
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
  static void                 GetEtaPhiDiff(const AliVTrack *t, const AliVCluster *v, Double_t &phidiff, Double_t &etadiff);
  static Byte_t               GetTrackType(const AliVTrack *t);
  static Byte_t               GetTrackType(const AliAODTrack *aodTrack, UInt_t filterBit1, UInt_t filterBit2);
  static Double_t             DeltaPhi(Double_t phia, Double_t phib, Double_t rMin = -TMath::Pi()/2, Double_t rMax = 3*TMath::Pi()/2);
  static std::vector<double>  GenerateFixedBinArray(int n, double min, double max, bool last = true);
  static void                 GenerateFixedBinArray(int n, double min, double max, std::vector<double>& array, bool last = true);
  static std::vector<double>  GenerateLogFixedBinArray(int n, double min, double max, bool last = true);
  static void                 GenerateLogFixedBinArray(int n, double min, double max, std::vector<double>& array, bool last = true);
  static Double_t             GetParallelFraction(AliVParticle* part1, AliVParticle* part2);
  static Double_t             GetParallelFraction(const TVector3& vect1, AliVParticle* part2);
  static EBeamType_t          BeamTypeFromRunNumber(Int_t runnumber);

  static Double_t             fgkEMCalDCalPhiDivide;       ///<  phi value used to distinguish between DCal and EMCal

  // Task configuration
  EBeamType_t                 fForceBeamType;              ///< forced beam type
  Bool_t                      fGeneralHistograms;          ///< whether or not it should fill some general histograms
  Bool_t                      fCreateHisto;                ///< whether or not create histograms
  Bool_t                      fNeedEmcalGeom;              ///< whether or not the task needs the emcal geometry
  std::vector<double>         fCentBins;                   ///< how many centrality bins
  ECentralityEstimation_t     fCentralityEstimation;       ///< Centrality estimation

  // Input data
  Bool_t                      fIsPythia;                   ///< if it is a PYTHIA production
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
  TString                     fSelectGeneratorName;        ///< Selects only events produced by a generator that has a name containing a string
  Double_t                    fMinimumEventWeight;         ///< Minimum event weight for the related bookkeping histogram
  Double_t                    fMaximumEventWeight;         ///< Minimum event weight for the related bookkeping histogram

  // Service fields
  Bool_t                      fInhibit;                    //!<!inhibit execution of the task
  Bool_t                      fLocalInitialized;           //!<!whether or not the task has been already initialized
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
  AliGenPythiaEventHeader    *fPythiaHeader;               //!<!event Pythia header
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

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalLight, 4);
  /// \endcond
};

#endif
