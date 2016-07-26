#ifndef ALIANALYSISTASKEMCAL_H
#define ALIANALYSISTASKEMCAL_H

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
class AliEmcalPythiaInfo;

#include "Rtypes.h"

#include "AliParticleContainer.h"
#include "AliMCParticleContainer.h"
#include "AliTrackContainer.h"
#include "AliClusterContainer.h"

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEmcal : public AliAnalysisTaskSE {
 public:

  enum BeamType {
    kNA       = -1,
    kpp       = 0,
    kAA       = 1,
    kpA       = 2
  };

  enum TriggerType {
    kND       = -1,
    kJ1       = 0,
    kJ2       = 1,
    kG1	      = 2,
    kG2       = 3,
    kL0       = 4
  };

  enum TriggerCategory {       // Online trigger categories
    kTriggerLevel0      = 0,   
    kTriggerLevel1Jet   = 1,
    kTriggerLevel1Gamma = 2,
    kTriggerRecalcJet   = 3,   // Recalculated max trigger patch; does not need to be above trigger threshold
    kTriggerRecalcGamma = 4
  };

  enum EMCalTriggerMode_t {
    kNoSpecialTreatment,       // No special treatment for EMCal triggers
    kOverlapWithLowThreshold   // The overlap between low and high threshold trigger is assigned to the lower threshold only
  };

  AliAnalysisTaskEmcal();
  AliAnalysisTaskEmcal(const char *name, Bool_t histo=kFALSE);
  virtual ~AliAnalysisTaskEmcal();

  AliParticleContainer       *AddParticleContainer(const char *n);
  AliTrackContainer          *AddTrackContainer(const char *n);
  AliMCParticleContainer     *AddMCParticleContainer(const char *n);
  AliClusterContainer        *AddClusterContainer(const char *n);
  void                        AdoptParticleContainer(AliParticleContainer* cont)    { fParticleCollArray.Add(cont)                        ; }
  void                        AdoptClusterContainer(AliClusterContainer* cont)      { fClusterCollArray.Add(cont)                         ; }
  AliParticleContainer       *GetParticleContainer(Int_t i=0)         const;
  AliParticleContainer       *GetParticleContainer(const char* name)  const;
  AliClusterContainer        *GetClusterContainer(Int_t i=0)          const;
  AliClusterContainer        *GetClusterContainer(const char* name)   const;
  AliMCParticleContainer     *GetMCParticleContainer(Int_t i=0)               const { return dynamic_cast<AliMCParticleContainer*>(GetParticleContainer(i))   ; }
  AliMCParticleContainer     *GetMCParticleContainer(const char* name)        const { return dynamic_cast<AliMCParticleContainer*>(GetParticleContainer(name)); }
  AliTrackContainer          *GetTrackContainer(Int_t i=0)                    const { return dynamic_cast<AliTrackContainer*>(GetParticleContainer(i))        ; }
  AliTrackContainer          *GetTrackContainer(const char* name)             const { return dynamic_cast<AliTrackContainer*>(GetParticleContainer(name))     ; }
  void                        RemoveParticleContainer(Int_t i=0)                    { fParticleCollArray.RemoveAt(i)                      ; } 
  void                        RemoveClusterContainer(Int_t i=0)                     { fClusterCollArray.RemoveAt(i)                       ; } 
  void                        SetCaloCellsName(const char *n)                       { fCaloCellsName     = n                              ; }
  void                        SetCaloTriggerPatchInfoName(const char *n)            { fCaloTriggerPatchInfoName = n                       ; }
  void                        SetCaloTriggersName(const char *n)                    { fCaloTriggersName  = n                              ; }
  void                        SetCentRange(Double_t min, Double_t max)              { fMinCent           = min  ; fMaxCent = max          ; }
  void                        SetCentralityEstimator(const char *c)                 { fCentEst           = c                              ; }
  void                        SetClusName(const char *n)                            { AddClusterContainer(n)                              ; }
  void                        SetClusPtCut(Double_t cut, Int_t c=0);
  void                        SetClusTimeCut(Double_t min, Double_t max, Int_t c=0);
  void                        SetEventPlaneVsEmcal(Double_t ep)                     { fEventPlaneVsEmcal = ep                             ; }
  void                        SetForceBeamType(BeamType f)                          { fForceBeamType     = f                              ; }
  void                        SetHistoBins(Int_t nbins, Double_t min, Double_t max) { fNbins = nbins; fMinBinPt = min; fMaxBinPt = max    ; }
  void                        SetIsEmbedded(Bool_t i)                               { fIsEmbedded        = i                              ; }
  void                        SetIsPythia(Bool_t i)                                 { fIsPythia          = i                              ; }
  void                        SetMakeGeneralHistograms(Bool_t g)                    { fGeneralHistograms = g                              ; }
  void                        SetMCLabelShift(Int_t s)                              { fMCLabelShift      = s                              ; }
  void                        SetMinMCLabel(Int_t s)                                { fMinMCLabel        = s                              ; }
  void                        SetMinNTrack(Int_t min)                               { fMinNTrack         = min                            ; }
  void                        SetMinPtTrackInEmcal(Double_t min)                    { fMinPtTrackInEmcal = min                            ; }
  virtual void                SetNCentBins(Int_t n)                                 { fNcentBins         = n                              ; }
  void                        SetNeedEmcalGeom(Bool_t n)                            { fNeedEmcalGeom     = n                              ; }
  void                        SetOffTrigger(UInt_t t)                               { fOffTrigger        = t                              ; }
  void                        SetTrackEtaLimits(Double_t min, Double_t max, Int_t c=0);
  void                        SetTrackPhiLimits(Double_t min, Double_t max, Int_t c=0);
  void                        SetTrackPtCut(Double_t cut, Int_t c=0);
  void                        SetTracksName(const char *n)                          { AddParticleContainer(n)                             ; }
  void                        SetTrigClass(const char *n)                           { fTrigClass         = n                              ; } 
  void                        SetTriggerTypeSel(TriggerType t)                      { fTriggerTypeSel    = t                              ; } 
  void                        SetUseAliAnaUtils(Bool_t b, Bool_t bRejPilup = kTRUE) { fUseAliAnaUtils    = b ; fRejectPileup = bRejPilup  ; }
  void                        SetVzRange(Double_t min, Double_t max)                { fMinVz             = min  ; fMaxVz   = max          ; }
  void                        SetUseSPDTrackletVsClusterBG(Bool_t b)                { fTklVsClusSPDCut   = b                              ; }
  void                        SetEMCalTriggerMode(EMCalTriggerMode_t m)             { fEMCalTriggerMode  = m                              ; }
  void                        SetUseNewCentralityEstimation(Bool_t b)               { fUseNewCentralityEstimation = b                     ; }
  void                        SetGeneratePythiaInfoObject(Bool_t b)                 { fGeneratePythiaInfoObject = kTRUE                   ; }
  void                        SetPythiaInfoName(const char *n)                      { fPythiaInfoName    = n                              ; }
  const TString&              GetPythiaInfoName()                             const { return fPythiaInfoName                              ; }
  const AliEmcalPythiaInfo   *GetPythiaInfo()                                 const { return fPythiaInfo                                  ; }

 protected:  
  void                        LoadPythiaInfo(AliVEvent *event);
  void                        SetRejectionReasonLabels(TAxis* axis);
  Bool_t                      AcceptCluster(AliVCluster *clus, Int_t c = 0)      const;
  Bool_t                      AcceptTrack(AliVParticle *track, Int_t c = 0)      const;
  void                        AddObjectToEvent(TObject *obj, Bool_t attempt = kFALSE);
  AliVParticle               *GetAcceptParticleFromArray(Int_t p, Int_t c=0)     const;
  AliVCluster                *GetAcceptClusterFromArray(Int_t cl, Int_t c=0)     const;
  TClonesArray               *GetArrayFromEvent(const char *name, const char *clname=0);
  BeamType                    GetBeamType();
  TClonesArray               *GetParticleArray(Int_t i=0)                        const;
  TClonesArray               *GetClusterArray(Int_t i=0)                         const;
  Int_t                       GetNParticles(Int_t i=0)                           const;
  Int_t                       GetNClusters(Int_t i=0)                            const;
  AliEMCALTriggerPatchInfo   *GetMainTriggerPatch(TriggerCategory triggersel = kTriggerLevel1Jet, Bool_t doOfflinSimple = kFALSE);
  Bool_t		                  HasTriggerType(TriggerType triggersel);
  ULong_t 		                GetTriggerList();
  Bool_t                      PythiaInfoFromFile(const char* currFile, Float_t &fXsec, Float_t &fTrials, Int_t &pthard);
  Bool_t                      IsTrackInEmcalAcceptance(AliVParticle* part, Double_t edges=0.9) const;

  void                        GeneratePythiaInfoObject(AliMCEvent* mcEvent);

  // Overloaded AliAnalysisTaskSE methods
  void                        UserCreateOutputObjects();
  void                        UserExec(Option_t *option);
  Bool_t                      UserNotify();

  // Virtual functions, to be overloaded in derived classes
  virtual void                ExecOnce();
  virtual Bool_t              FillGeneralHistograms();
  virtual Bool_t              IsEventSelected();
  virtual Bool_t              RetrieveEventObjects();
  virtual Bool_t              FillHistograms()                  { return kTRUE                 ; }
  virtual Bool_t              Run()                             { return kTRUE                 ; }
    
  // Static utilities
  static void                 GetEtaPhiDiff(const AliVTrack *t, const AliVCluster *v, Double_t &phidiff, Double_t &etadiff);
  static Byte_t               GetTrackType(const AliVTrack *t);
  static Byte_t               GetTrackType(const AliAODTrack *aodTrack, UInt_t filterBit1, UInt_t filterBit2);
  static Double_t             DeltaPhi(Double_t phia, Double_t phib, Double_t rMin = -TMath::Pi()/2, Double_t rMax = 3*TMath::Pi()/2);
  static Double_t*            GenerateFixedBinArray(Int_t n, Double_t min, Double_t max);
  static void                 GenerateFixedBinArray(Int_t n, Double_t min, Double_t max, Double_t* array);
  static Double_t             GetParallelFraction(AliVParticle* part1, AliVParticle* part2);
  static Double_t             GetParallelFraction(const TVector3& vect1, AliVParticle* part2);

  static Double_t             fgkEMCalDCalPhiDivide;      //  phi value used to distinguish between DCal and EMCal

  // Task configuration
  TString                     fPythiaInfoName;             // name of pythia info object
  BeamType                    fForceBeamType;              // forced beam type
  Bool_t                      fGeneralHistograms;          // whether or not it should fill some general histograms
  Bool_t                      fInitialized;                // whether or not the task has been already initialized
  Bool_t                      fCreateHisto;                // whether or not create histograms
  TString                     fCaloCellsName;              // name of calo cell collection
  TString                     fCaloTriggersName;           // name of calo triggers collection
  TString                     fCaloTriggerPatchInfoName;   // trigger patch info array name
  Double_t                    fMinCent;                    // min centrality for event selection
  Double_t                    fMaxCent;                    // max centrality for event selection
  Double_t                    fMinVz;                      // min vertex for event selection
  Double_t                    fMaxVz;                      // max vertex for event selection
  Double_t                    fTrackPtCut;                 // cut on track pt in event selection
  Int_t                       fMinNTrack;                  // minimum nr of tracks in event with pT>fTrackPtCut
  Bool_t                      fUseAliAnaUtils;             // used for LHC13* data: z-vtx, Ncontributors, z-vtx resolution cuts
  Bool_t                      fRejectPileup;               // Reject pilup using function AliAnalysisUtils::IsPileUpEvent()
  Bool_t                      fTklVsClusSPDCut;            // Apply tracklet-vs-cluster SPD cut to reject background events in pp
  UInt_t                      fOffTrigger;                 // offline trigger for event selection
  TString                     fTrigClass;                  // trigger class name for event selection
  TriggerType                 fTriggerTypeSel;             // trigger type to select based on trigger patches
  Int_t                       fNbins;                      // no. of pt bins
  Double_t                    fMinBinPt;                   // min pt in histograms
  Double_t                    fMaxBinPt;                   // max pt in histograms
  Double_t                    fMinPtTrackInEmcal;          // min pt track in emcal
  Double_t                    fEventPlaneVsEmcal;          // select events which have a certain event plane wrt the emcal
  Double_t                    fMinEventPlane;              // minimum event plane value
  Double_t                    fMaxEventPlane;              // maximum event plane value
  TString                     fCentEst;                    // name of V0 centrality estimator
  Bool_t                      fIsEmbedded;                 // trigger, embedded signal
  Bool_t                      fIsPythia;                   // trigger, if it is a PYTHIA production
  Int_t                       fSelectPtHardBin;            // select one pt hard bin for analysis
  Int_t                       fMinMCLabel;                 // minimum MC label value for the tracks/clusters being considered MC particles
  Int_t                       fMCLabelShift;               // if MC label > fMCLabelShift, MC label -= fMCLabelShift
  Int_t                       fNcentBins;                  // how many centrality bins
  Bool_t                      fNeedEmcalGeom;              // whether or not the task needs the emcal geometry
  TObjArray                   fParticleCollArray;          // particle/track collection array
  TObjArray                   fClusterCollArray;           // cluster collection array
  ULong_t                     fTriggers;                   // list of fired triggers
  EMCalTriggerMode_t          fEMCalTriggerMode;           // EMCal trigger selection mode
  Bool_t                      fUseNewCentralityEstimation; // Use new centrality estimation (for 2015 data)
  Bool_t                      fGeneratePythiaInfoObject;   // Generate Pythia info object

  // Service fields
  AliAnalysisUtils           *fAliAnalysisUtils;           //!vertex selection (optional)
  Bool_t                      fIsEsd;                      //!whether it's an ESD analysis
  AliEMCALGeometry           *fGeom;                       //!emcal geometry
  TClonesArray               *fTracks;                     //!tracks
  TClonesArray               *fCaloClusters;               //!clusters
  AliVCaloCells              *fCaloCells;                  //!cells
  AliVCaloTrigger            *fCaloTriggers;               //!calo triggers
  TClonesArray               *fTriggerPatchInfo;           //!trigger patch info array
  Double_t                    fCent;                       //!event centrality
  Int_t                       fCentBin;                    //!event centrality bin
  Double_t                    fEPV0;                       //!event plane V0
  Double_t                    fEPV0A;                      //!event plane V0A
  Double_t                    fEPV0C;                      //!event plane V0C
  Double_t                    fVertex[3];                  //!event vertex
  Int_t                       fNVertCont;                  //!event vertex number of contributors
  BeamType                    fBeamType;                   //!event beam type
  AliGenPythiaEventHeader    *fPythiaHeader;               //!event Pythia header
  Double_t                    fPtHard;                     //!event pt hard
  Int_t                       fPtHardBin;                  //!event pt hard bin
  Int_t                       fNTrials;                    //!event trials
  Float_t                     fXsection;                   //!x-section from pythia header
  AliEmcalPythiaInfo         *fPythiaInfo;                 //!event parton info

  // Output
  TList                      *fOutput;                     //!output list
  TH1                        *fHistEventCount;             //!incoming and selected events
  TH1                        *fHistTrialsAfterSel;         //!total number of trials per pt hard bin after selection
  TH1                        *fHistEventsAfterSel;         //!total number of events per pt hard bin after selection
  TProfile                   *fHistXsectionAfterSel;       //!x section from pythia header
  TH1                        *fHistTrials;                 //!trials from pyxsec.root
  TH1                        *fHistEvents;                 //!total number of events per pt hard bin
  TProfile                   *fHistXsection;               //!x section from pyxsec.root
  TH1                        *fHistPtHard;                 //!pt hard distribution
  TH1                        *fHistCentrality;             //!event centrality distribution
  TH1                        *fHistZVertex;                //!z vertex position
  TH1                        *fHistEventPlane;             //!event plane distribution
  TH1                        *fHistEventRejection;         //!book keep reasons for rejecting event
  TH1                        *fHistTriggerClasses;         //!number of events in each trigger class

 private:
  AliAnalysisTaskEmcal(const AliAnalysisTaskEmcal&);            // not implemented
  AliAnalysisTaskEmcal &operator=(const AliAnalysisTaskEmcal&); // not implemented

  ClassDef(AliAnalysisTaskEmcal, 14) // EMCAL base analysis task
};

//________________________________________________________________________
inline Double_t AliAnalysisTaskEmcal::DeltaPhi(Double_t phia, Double_t phib, Double_t rangeMin, Double_t rangeMax) 
{
  // Calculate Delta Phi.

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

//________________________________________________________________________
inline void AliAnalysisTaskEmcal::GenerateFixedBinArray(Int_t n, Double_t min, Double_t max, Double_t* array)
{
  Double_t binWidth = (max-min)/n;
  array[0] = min;
  for (Int_t i = 1; i <= n; i++) {
    array[i] = array[i-1]+binWidth;
  }
}

//________________________________________________________________________
inline Double_t* AliAnalysisTaskEmcal::GenerateFixedBinArray(Int_t n, Double_t min, Double_t max)
{
  Double_t *array = new Double_t[n+1];
  GenerateFixedBinArray(n, min, max, array);
  return array;
}

#endif
