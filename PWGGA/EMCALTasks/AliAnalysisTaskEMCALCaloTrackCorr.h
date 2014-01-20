#ifndef ALIANALYSISTASKEMCALCALOTRACKCORR_cxx
#define ALIANALYSISTASKEMCALCALOTRACKCORR_cxx

class TList;
class TH1F;
class TH2F;
class TH1I;
class TString;
class TGeoHMatrix;
class TClonesArray;
class AliEMCALGeometry;
class AliEMCALRecoUtils;
class AliESDtrackCuts;
class AliESDEvent;
class AliMCEvent;
class AliStack;
class AliVCluster;
class AliFiducialCut;
class AliCaloTrackParticle;
class AliCentrality;
class AliEventplane;
class AliAnalysisManager;
class AliInputEventHandler;


#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEMCALCaloTrackCorr : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskEMCALCaloTrackCorr(const char *name = "AliAnalysisTaskEMCALCaloTrackCorr");
  virtual ~AliAnalysisTaskEMCALCaloTrackCorr() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

//  virtual  AliFiducialCut  *GetFiducialCut() { if(!fFidCut)  
//                      fFidCut = new AliFiducialCut(); return  fFidCut ; } 


  enum type {kPtThresholdIC=0, kSumPtInConeIC=1, kPtFracationIC=2, kSumPtFracationIC=3};
  enum particleInCone { kIsolatedNeutralAndCharged=0, kIsolatedOnlyNeutral=1, kIsolatedOnlyCharged=2  };  

  Int_t    GetMinNCells()           const { return fMinNCells  ; }
  Double_t GetMinE()                const { return fMinE       ; }
  Double_t GetMinDistBad()          const { return fMinDistBad ; }  

  Int_t    GetDebug()               const { return fDebug ; }
  void     SetDebug(Int_t deb)            { fDebug = deb  ; }

  Bool_t   IsDataMC()               const { return kMC ; }
  void     SetMC(Bool_t mc)               { kMC = mc   ; }

  TString  GetDataType()            const { return fDataType ; }
  void     SetDataType(TString data)      { fDataType = data ; }

  Int_t    GetHistoPtBins()         const { return fHistoPtBins  ; }
  Float_t  GetHistoPtMin()          const { return fHistoPtMin   ; }
  Float_t  GetHistoPtMax()          const { return fHistoPtMax   ; }

  void    SetHistoPtRangeAndNBins(Float_t min, Float_t max, Int_t n)
          {fHistoPtBins = n; fHistoPtMax = max; fHistoPtMin = min;}

  Int_t    GetHistoPhiBins()        const { return fHistoPhiBins  ; }
  Float_t  GetHistoPhiMin()         const { return fHistoPhiMin   ; }
  Float_t  GetHistoPhiMax()         const { return fHistoPhiMax   ; }

  void    SetHistoPhiRangeAndNBins(Float_t min, Float_t max, Int_t n) 
          {fHistoPhiBins = n; fHistoPhiMax = max; fHistoPhiMin = min;}

  Int_t    GetHistoEtaBins()        const { return fHistoEtaBins  ; }
  Float_t  GetHistoEtaMin()         const { return fHistoEtaMin   ; }
  Float_t  GetHistoEtaMax()         const { return fHistoEtaMax   ; }

  void    SetHistoEtaRangeAndNBins(Float_t min, Float_t max, Int_t n) 
          {fHistoEtaBins = n; fHistoEtaMax = max; fHistoEtaMin = min;}



  Float_t  GetConeR()               const { return fSetConeR          ; }
  Float_t  GetPtThreshold()         const { return fSetPtThreshold    ; }
  Float_t  GetSumPtThreshold()      const { return fSetSumPtThreshold ; }
  Float_t  GetPtFraction()          const { return fSetPtFraction     ; }
  TString  GetICMethod()            const { return fICMethod          ; }
  TString  GetParticleTypeInCone()  const { return fParticlesInCone    ; }
  
  void   SetMinNCells(Int_t n)                { fMinNCells   = n  ; }
  void   SetMinE(Double_t pt)                 { fMinE        = pt ; }
  void   SetMinDistBad(Double_t fn)           { fMinDistBad = fn ; }
  
  void   SetConeR(Float_t r)              { fSetConeR          = r       ; }
  void   SetPtThreshold(Float_t pt)       { fSetPtThreshold    = pt      ; }
  void   SetSumPtThreshold(Float_t s)     { fSetSumPtThreshold = s       ; }
  void   SetPtFraction(Float_t pt)        { fSetPtFraction     = pt      ; }
  void   SetICMethod(Int_t iMethod )      { fICMethod          = iMethod ; }
  void   SetParticleTypeInCone(Int_t i)   { fParticlesInCone    = i       ; }
 
  void   SwitchOnAnaIsolated()            { kDoIsolatedAna = kTRUE  ; }
  void   SwitchOffAnaIsolated()           { kDoIsolatedAna = kFALSE ; }

  void   SwitchOnTrackMultBins()          { kDoTrackMultBins = kTRUE ; }
  void   SwitchOffTrackMultBins()         { kDoTrackMultBins = kTRUE ; }

  void   SetAnaUELeftRightOrNearAway(Bool_t leftright, Bool_t nearaway) {
                                                   kUELeftRight = leftright,
                                                   kUENearAway  = nearaway ; }
  void   SwitchOnTwoTracksCorr()     { kTwoTracksCorr = kTRUE    ; }
  void   SwitchOffTwoTracksCorr()    { kTwoTracksCorr = kFALSE   ; }  

  void   SelectDecayPhotonCorr(Bool_t decayphoton) {kDecayPhotonCorr = decayphoton;}

  void   SetAnaMCTruthOrPrimaryCorr(Bool_t kAnatruth, Bool_t kAnaprimary){
                                             kAnaMCTruthCorr = kAnatruth,
                                             kAnaMCPrimaryCorr = kAnaprimary ; }

  void   SetAnaMCPrimaryParticle(Bool_t kpi0, Bool_t keta, Bool_t kphoton){
                                                    kAnaPi0Prim = kpi0,
                                                    kAnaEtaPrim = keta,
                                                    kAnaPhotonPrim = kphoton ; }
  // Taking the absolute leading as the trigger or not
  Bool_t  DoAbsoluteLeading()         const { return kMakeAbsoluteLeading   ; }
  void    SwitchOnAbsoluteLeading()         { kMakeAbsoluteLeading = kTRUE  ; }
  void    SwitchOffAbsoluteLeading()        { kMakeAbsoluteLeading = kFALSE ; }

  // Taking the near side leading as the trigger or not
  Bool_t  DoNearSideLeading()         const { return kMakeNearSideLeading   ; }
  void    SwitchOnNearSideLeading()         { kMakeNearSideLeading = kTRUE  ; }
  void    SwitchOffNearSideLeading()        { kMakeNearSideLeading = kFALSE ; }
  
  void    SwitchOnInAcceptance()            { kPhotonInAcceptance = kTRUE   ; }
  void    SwitchOffInAcceptance()           { kPhotonInAcceptance = kFALSE  ; }
 
  void    SwitchOnAnaMapping()              { kAnaDecayMapping = kTRUE  ; }
  void    SwitchOffAnaMapping()             { kAnaDecayMapping = kFALSE ; }

//  Bool_t  IsFiducialCutOn()           const { return fCheckFidCut           ; }
//  void    SwitchOnFiducialCut()             { fCheckFidCut = kTRUE;
 //                                  if(!fFidCut)fFidCut = new AliFiducialCut();}
//  void    SwitchOffFiducialCut()            { fCheckFidCut = kFALSE         ; }

  void    SetNTriggPtBins(Int_t nbins)      { fNTriggPtBins = nbins ; }
  Int_t   GetNTriggPtBins()           const { return fNTriggPtBins  ; }
  void    SetTriggerBins(Float_t *ptTriggBins);

  void    SetNAssocPtBins(Int_t mbins)      { fNAssocPtBins = mbins ; }
  Int_t   GetNAssocPtBins()           const { return fNAssocPtBins  ; }
  void    SetAssociatedBins(Float_t *ptAssocBins);

  Float_t GetDeltaPhiMaxCut()         const { return fDeltaPhiMaxCut ; }
  Float_t GetDeltaPhiMinCut()         const { return fDeltaPhiMinCut ; }
  void     SetDeltaPhiCutRange(Float_t phimin, Float_t phimax)
                {fDeltaPhiMaxCut = phimax;  fDeltaPhiMinCut = phimin ; } 
 
  void   SetMixedEventsPool(Int_t n)        { nMixedEvents = n ; }
  void   SwitchOnAnaMixEvent()              { kDoMixEventsAna = kTRUE    ; }
  void   SwitchOffAnaMixEvent()             { kDoMixEventsAna = kFALSE   ; }

  void   SwitchOnFillMesonAOD()             { kDoMesonFill = kTRUE  ; }
  void   SwitchOffFillMesonAOD()            { kDoMesonFill = kFALSE ; }

  void   SwitchOnFillMesonHistos()          { kNeutralMesonHistos = kTRUE  ; }
  void   SwitchOffFillMesonHistos()         { kNeutralMesonHistos = kFALSE ; }

  void   SwitchOnAnaMesonCorr()             { kDoMesonCorrAna = kTRUE  ; }
  void   SwitchOffAnaMesonCorr()            { kDoMesonCorrAna = kFALSE ; }

  void   SwitchOnAnaPhotonCorr()            { kDoPhotonCorrAna = kTRUE  ; }
  void   SwitchOffAnaPhotonCorr()           { kDoPhotonCorrAna = kFALSE ; }  

  void   SwitchOnEventTriggerAtSE()         { kEventTriggerAtSE = kTRUE ; }
  void   SwitchOffEventTriggerAtSE()        { kEventTriggerAtSE = kFALSE; }

  void   SwithchOnPhotonPairTimeCut()       { kPhotonPairTimeCut = kTRUE ;}
  void   SwithchOffPhotonPairTimeCut()      { kPhotonPairTimeCut = kFALSE ;}

  void   SwithchOnPhotonIDCut()             { kDoPhotonIDCut = kTRUE  ; }
  void   SwithchOffPhotonIDCut()            { kDoPhotonIDCut = kFALSE ; }

  TString GetAnaTypeInIsolated()         const { return fAnaTypeInIsolated ; }
  void    SetAnaTypeInIsolated(TString & part) { fAnaTypeInIsolated = part ; }

  Float_t GetDeltaPhiHRSize()           const { return fDeltaPhiHRSize   ; }
  void SetDeltaPhiHRSize(Float_t fHRphi)       { fDeltaPhiHRSize = fHRphi ; }

  Float_t  GetUeDeltaPhiSize()           const { return fUeDeltaPhiSize ; }
  Float_t  GetUeDeltaPhiFix()            const { return fUeDeltaPhiFix  ; }
  void     SetUeDeltaPhiFixAndSize(Float_t uefix, Float_t uesize)
                    { fUeDeltaPhiFix = uefix,  fUeDeltaPhiSize = uesize ; }

  void     SetAssociatedPtBegin(Float_t begin) {fptAssociatedBegin = begin;}

  void     SetLargeCorrTrigger(Float_t ftrigger1, Float_t ftrigger2) 
                 {fptTriggerBegin = ftrigger1, fptTriggerEnd   = ftrigger2; }

  void    SwitchOnAsymmetryCut()             { kDoAsymmetryCut = kTRUE  ; }
  void    SwitchOffAsymmetryCut()            { kDoAsymmetryCut = kFALSE ; }
  void    SetAsymmetryCut(Float_t asycut)    { fAsymmetryCut = asycut ; }

  void    SwitchOnAODHybridTrackSelection()  { kDoSelectHybridTracks = kTRUE ;}
  void    SwitchOffAODHybridTrackSelection() { kDoSelectHybridTracks = kFALSE;}  
  void    SetAnaMesonType(TString fmesontype){ fAnaMesonType = fmesontype ; }

  void    SetCentralityBin(Int_t min, Int_t max)
                       { fCentralityBinMin = min; fCentralityBinMax = max ; }
  void    SetCentralityClass(TString name)   { fCentralityClass   = name  ; }
  
  void    SetEventPlaneMethod(TString m)     { fEventPlaneMethod = m      ; }

  UInt_t  GetEventTriggerMask()        const { return fEventTriggerMaks   ; }
  void    SetEventTriggerMask(UInt_t evtTrig = AliVEvent::kAny)
                                        { fEventTriggerMaks = evtTrig     ; }
  void    SetNCentralityBins(Int_t fcenbins) { fNCentralityBins = fcenbins; }
  void    SetNEventPlaneBins(Int_t fevbins)  { fNEventPlaneBins = fevbins ; }

  void    SetEMCALLambda0Cut(Float_t l0min, Float_t l0max)
                                   { fL0CutMin = l0min, fL0CutMax = l0max ; }
  void    SetClusterTimeCut(Float_t timemin, Float_t timemax){
                             fTimeCutMin = timemin, fTimeCutMax = timemax ; }
  void    SetPhotonPairDeltaTimeCut(Float_t deltatime)
                                         { fPhotonPairTimeCut = deltatime ; }
 
  void    SetTrackMatchedDPhiCut(Float_t dphicut){ fEMCALDPhiCut = dphicut; }
  void    SetTrackMatchedDEtaCut(Float_t detacut){ fEMCALDEtaCut = detacut; }

  void    SetMesonInMassRangeCut(Float_t inmass1, Float_t inmass2){
                        fInvMassMinCut = inmass1, fInvMassMaxCut = inmass2 ; }

  void    SetMesonInMassLeftRangeCut(Float_t leftmin, Float_t leftmax){
                      fLeftBandMinCut = leftmin, fLeftBandMaxCut = leftmax ; }

  void    SetMesonInMassRightRangeCut(Float_t rightmin, Float_t rightmax){
                  fRightBandMinCut = rightmin, fRightBandMaxCut = rightmax ; }
 
  void    SetEMCALGeometryName(TString name) { fEMCALGeomName = name ; }

  void    SetTrackCuts(AliESDtrackCuts * cuts);
 
  void    SetZvertexCut(Float_t fzcut) { fZVertexCut = fzcut ; }

  void    SetTrackFilterMask(ULong_t bit)  { fTrackFilterMask = bit ; }
  
 
private:
  AliAnalysisTaskEMCALCaloTrackCorr(const AliAnalysisTaskEMCALCaloTrackCorr&); // not implemented
  AliAnalysisTaskEMCALCaloTrackCorr& operator=(const AliAnalysisTaskEMCALCaloTrackCorr&); // not implemented

  ////////////Add function
  void    InitParameters();
 
  Bool_t  FillInputEvent();

  Bool_t  SelectPair(AliCaloTrackParticle *mesonCandidate);
  Bool_t  IsolatedPhoton(TClonesArray *fEMCALEventIsolated, 
                         TClonesArray *fCTSEventIsolated,
                         Int_t fIndexPhotonCan,  Double_t ptPhotonCan,
                         Double_t phiPhotonCan, Double_t etaPhotonCan);
  void    FillInputPhoton();
  void    FillInputMeson() ;
  void    FillInputTrack() ;
  Bool_t  MakeChargedCorrelation(Int_t fTrackIndex, Double_t ptTrigg, 
                               Double_t phiTrigg, Double_t etaTrigg);
  void    MakeChargedMixCorrelation(Double_t ptTriggMix,  Double_t phiTriggMix,
                            Double_t etaTriggMix, TList *poolMix);
private:
  enum {kNtriggPtBins=10, kNassocPtBins=10};

  AliAnalysisManager    *fManager;
  AliInputEventHandler  *fInputHandler;

  AliVEvent   *fEvent;
  AliMCEvent    *fMCEvent;
  AliStack      *fStack;
  AliCentrality *fCentrality;
  AliEventplane *fEventPlane;

  AliEMCALRecoUtils *fEMCALRecU;
  AliEMCALGeometry  *fEMCALGeom;
  AliESDtrackCuts   *fESDtrackCuts;
  //  AliFiducialCut    *fFidCut;

  TList    *outputContainer;

  TString  fEMCALGeomName; 
  TString  fCentralityClass;    
  Int_t    fCentralityBinMin;
  Int_t    fCentralityBinMax;
  TString  fEventPlaneMethod;
  UInt_t   fEventTriggerMaks;
  Int_t    fNCentralityBins;
  Int_t    fNEventPlaneBins;
 
  Int_t    fHistoPtBins ; 
  Float_t  fHistoPtMax  ; 
  Float_t  fHistoPtMin  ; 
  Int_t    fHistoPhiBins; 
  Float_t  fHistoPhiMax ;
  Float_t  fHistoPhiMin ;
  Int_t    fHistoEtaBins;
  Float_t  fHistoEtaMax ;
  Float_t  fHistoEtaMin ;
 
  Int_t    fMinNCells;
  Float_t  fMinE;
  Double_t fMinDistBad;
  Float_t  fL0CutMin;
  Float_t  fL0CutMax;
  Float_t  fTimeCutMin;
  Float_t  fTimeCutMax;
  Float_t  fPhotonPairTimeCut;
  Float_t  fEMCALDPhiCut;
  Float_t  fEMCALDEtaCut;
  Float_t  fZVertexCut; 
  Int_t    fDebug;         
  TString  fAnaMesonType;
  Float_t  fAsymmetryCut;
  TString  fDataType;
  ULong_t  fTrackFilterMask;

  Float_t  fInvMassMinCut;
  Float_t  fInvMassMaxCut;
  Float_t  fLeftBandMinCut;
  Float_t  fLeftBandMaxCut;
  Float_t  fRightBandMinCut;
  Float_t  fRightBandMaxCut;

  Bool_t  kMC;
  Bool_t  kNeutralMesonHistos;
  Bool_t  kDoMixEventsAna;
  Bool_t  kDoPhotonCorrAna;
  Bool_t  kDoAsymmetryCut;
  Bool_t  kDoSelectHybridTracks;
  Bool_t  kDoMesonFill;
  Bool_t  kDoMesonCorrAna;
  Bool_t  kDoIsolatedAna;
  Bool_t  kDoTrackMultBins;
  Bool_t  kUELeftRight;
  Bool_t  kUENearAway;
  Bool_t  kDecayPhotonCorr;
  Bool_t  kAnaMCTruthCorr;
  Bool_t  kAnaMCPrimaryCorr;
  Bool_t  kAnaPi0Prim;
  Bool_t  kAnaEtaPrim;
  Bool_t  kAnaPhotonPrim;
  Bool_t  kMakeAbsoluteLeading;
  Bool_t  kMakeNearSideLeading;
  Bool_t  kTwoTracksCorr;
  Bool_t  kPhotonInAcceptance;
  Bool_t  kAnaDecayMapping;
//  Bool_t  fCheckFidCut;
  Bool_t  kEventTriggerAtSE;
  Bool_t  kPhotonPairTimeCut;
  Bool_t  kDoPhotonIDCut;

  TH1F     *fhNEvents;        
  Int_t    fnEvents;
  TH1F     *fhNEventsAnalyized;
  Int_t    fEventAnalyized;
  TClonesArray  *fPhotonEvent;
  TClonesArray  *fPhotonPairEvent;
  TClonesArray  *fCTSEvent;
  Int_t    nPhotonsEMCAL;
  Int_t    nTracksCTS;
  TString  fAnaTypeInIsolated;
  Int_t    nMixedEvents;
  Float_t  fSetConeR;
  Float_t  fSetPtThreshold;
  Float_t  fSetSumPtThreshold;
  Float_t  fSetPtFraction;
  Int_t    fICMethod;
  Int_t    fParticlesInCone;

  Float_t  *fTriggPtArray;
  Int_t    fNTriggPtBins;
  Float_t  fptTriggerBegin;
  Float_t  fptTriggerEnd;
  Float_t  *fAssocPtArray;
  Int_t    fNAssocPtBins;
  Float_t  fptAssociatedBegin;

  Float_t  fDeltaPhiMaxCut;
  Float_t  fDeltaPhiMinCut;
  Float_t  fUeDeltaPhiSize;
  Float_t  fUeDeltaPhiFix;
  Float_t  fDeltaPhiHRSize;

  TH1F * fhPhotonE;
  TH2F * fhPhotonPtPhi;
  TH2F * fhPhotonPtEta;
  TH2F * fhPhotonPhiEta;

  TH1F * fhMesonE;
  TH2F * fhMesonPtPhi;
  TH2F * fhMesonPtEta;
  TH2F * fhMesonPhiEta;

  TH2F * fhAnglePairNoCut;
  TH2F * fhInvMassPairNoCut;
  TH2F * fhAsyNoCut;
  TH2F * fhInvMassPairAsyCut; 
  TH2F * fhAnglePairAsyCut;
  TH2F * fhInvMassPairPhi;
  TH2F * fhInvMassPairEta;
  TH2F * fhInvMassPairAllCut;
  TH2F * fhAnglePairAllCut;
  TH2F * fhAsyAllCut;
  TH2F * fhPi0DecayPhoton1;
  TH2F * fhPi0DecayPhoton1Dphi;
  TH2F * fhDecayPhoton1Pi0Dphi;
  TH2F * fhPi0DecayPhoton2;
  TH2F * fhPi0DecayPhoton2Dphi;
  TH2F * fhDecayPhoton2Pi0Dphi;
  TH2F * fhDecayPhoton1Photon2;
  TH2F * fhDecayPhoton1Photon2Dphi;
  TH2F * fhDecayPhoton2Photon1Dphi;

  TH1F * fhNtracksAll;
  TH1F * fhNtracksEMC7;
  TH1F * fhNtracksAnyINT;
  TH1F * fhNtracksCentral;
  TH1F * fhNtracksSemiCentral;
  TH1F * fhNtracksOtherTirgger;
  TH1F * fhNtracksCorr;
  TH2F * fhTrackPtPhi;
  TH2F * fhTrackPtEta;
  TH2F * fhTrackPhiEta;
  TH2F * fhPtPhiLeading;        //! phi distribution vs pT of leading
  TH2F * fhPtEtaLeading;        //! eta distribution vs pT of leading
  TH2F * fhMixPtPhiLeading;        //! phi distribution vs pT of leading
  TH2F * fhMixPtEtaLeading;        //! eta distribution vs pT of leading
 
  TH2F * fhDPhiTriggPtAssocPt;
  TH2F * fhDEtaTriggPtAssocPt;
  TH2F * fhAssocPtTriggPt;
  TH2F * fhxELogTriggPt;
  TH2F * fhpoutTriggPt;
  TH2F * fhzTTriggPt;
  TH2F * fhxETriggPt;       
  TH2F * fhAssocPtTriggPtHR;
  TH2F * fhxELogTriggPtHR;
  TH2F * fhpoutTriggPtHR;
  TH2F * fhzTTriggPtHR;
  TH2F * fhxETriggPtHR;
  TH2F * fhNUeAssocPtTriggPt;
  TH2F * fhNUepoutTriggPt;
  TH2F * fhNUexETriggPt;
  TH2F * fhNUezTTriggPt;
  TH2F * fhNUexELogTriggPt;
  TH2F * fhNUeDPhiDEta;
  TH2F * fhAUeAssocPtTriggPt;
  TH2F * fhAUepoutTriggPt;
  TH2F * fhAUezTTriggPt; 
  TH2F * fhAUexETriggPt;
  TH2F * fhAUexELogTriggPt;
  TH2F * fhAUeDPhiDEta;

  TH2F * fhMCPtPhiLeading;        //! phi distribution vs pT of leading
  TH2F * fhMCPtEtaLeading;        //! eta distribution vs pT of leading

  TH2F * fhMCAssocPtTriggPt;
  TH2F * fhMCxELogTriggPt;
  TH2F * fhMCpoutTriggPt;
  TH2F * fhMCzTTriggPt;
  TH2F * fhMCxETriggPt;
  TH2F * fhMCAssocPtTriggPtHR;
  TH2F * fhMCxELogTriggPtHR;
  TH2F * fhMCpoutTriggPtHR;
  TH2F * fhMCzTTriggPtHR;
  TH2F * fhMCxETriggPtHR;
  TH2F * fhMCNUeAssocPtTriggPt;
  TH2F * fhMCNUepoutTriggPt;
  TH2F * fhMCNUexETriggPt;
  TH2F * fhMCNUezTTriggPt;
  TH2F * fhMCNUexELogTriggPt;
  TH2F * fhMCAUeAssocPtTriggPt;
  TH2F * fhMCAUepoutTriggPt;
  TH2F * fhMCAUezTTriggPt;
  TH2F * fhMCAUexETriggPt;
  TH2F * fhMCAUexELogTriggPt;

  TH2F * fhDPhiAssocPt15T;
  TH2F * fhDEtaAssocPt15T;
  TH2F * fhMixDPhiAssocPt15T;
  TH2F * fhMixDEtaAssocPt15T;
  TH2F * fhMCDPhiAssocPt15T;
  TH2F * fhMCDEtaAssocPt15T;

  TList* fListMixEvents[10][10][10]; 
 
  TH2F  *fhDPhiTriggPtT[kNassocPtBins];
  TH2F  *fhDEtaTriggPtT[kNassocPtBins];
  TH2F  *fhMixDPhiTriggPtT[kNassocPtBins];
  TH2F  *fhMixDEtaTriggPtT[kNassocPtBins];

  TH2F  *fhDPhiSumPtBin[kNtriggPtBins][kNassocPtBins];
  TH2F  *fhDEtaSumPtBin[kNtriggPtBins][kNassocPtBins];
  TH2F  *fhDPhiDEtaBin[kNtriggPtBins][kNassocPtBins];
  TH2F  *fhMixDPhiDEtaBin[kNtriggPtBins][kNassocPtBins];

  TH2F  *fhDPhiAssocPtA[kNtriggPtBins];
  TH2F  *fhDEtaAssocPtA[kNtriggPtBins];
  TH2F  *fhMixDPhiAssocPtA[kNtriggPtBins];
  TH2F  *fhMixDEtaAssocPtA[kNtriggPtBins];

  TH2F  *fhMCDPhiTriggPtT[kNassocPtBins];
  TH2F  *fhMCDEtaTriggPtT[kNassocPtBins];

  TH2F  *fhMCDPhiSumPtBin[kNtriggPtBins][kNassocPtBins];
  TH2F  *fhMCDEtaSumPtBin[kNtriggPtBins][kNassocPtBins];
  TH2F  *fhMCDPhiDEtaBin[kNtriggPtBins][kNassocPtBins];

  TH2F  *fhMCDPhiAssocPtA[kNtriggPtBins];
  TH2F  *fhMCDEtaAssocPtA[kNtriggPtBins];

  ClassDef(AliAnalysisTaskEMCALCaloTrackCorr, 1); // example of analysis
};

#endif
