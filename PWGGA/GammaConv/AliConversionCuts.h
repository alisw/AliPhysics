#ifndef ALICONVERSIONCUTS_H
#define ALICONVERSIONCUTS_H

// Class handling all kinds of selection cuts for Gamma Conversion analysis
// Authors: Svein Lindal, Daniel Lohner												*

#include "AliAODpidUtil.h"
#include "AliConversionPhotonBase.h"
#include "AliAODConversionMother.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliStack.h"
#include "AliAnalysisCuts.h"
#include "TH1F.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisManager.h"

class AliESDEvent;
class AliAODEvent;
class AliConversionPhotonBase;
class AliKFVertex;
class TH1F;
class TH2F;
class AliPIDResponse;
class AliAnalysisCuts;
class iostream;
class TList;
class AliAnalysisManager;

using namespace std;

class AliConversionCuts : public AliAnalysisCuts {
	
 public: 


  enum cutIds {
	kisHeavyIon,                  
	kCentralityMin,               
	kCentralityMax,               
	kselectV0AND,                 
	kmultiplicityMethod,             
	kremovePileUp,                
	kExtraSignals,	
	kv0FinderType,                
	ketaCut,                                     
	kRCut,                     
	ksinglePtCut,                 
	kclsTPCCut,                   
	kededxSigmaCut,               
	kpidedxSigmaCut,              
	kpiMomdedxSigmaCut,        
	kpiMaxMomdedxSigmaCut,        
	kLowPRejectionSigmaCut,       
	kTOFelectronPID,              
	kQtMaxCut,                    
	kchi2GammaCut,                
	kPsiPair, 
	kdoPhotonAsymmetryCut,
	kCosPAngle,
	kElecShare,
	kToCloseV0s,
	kNCuts
  };

  enum photonCuts {
      kPhotonIn=0,
      kOnFly,
      kNoTracks,
      kTrackCuts,
      kdEdxCuts,
      kConvPointFail,
      kPhotonCuts,
      kPhotonOut
  };


  Bool_t SetCutIds(TString cutString); 
  Int_t fCuts[kNCuts];
  Bool_t SetCut(cutIds cutID, Int_t cut);
  Bool_t UpdateCutString();


  static const char * fgkCutNames[kNCuts];

  Double_t GetCosineOfPointingAngle(const AliConversionPhotonBase * photon, AliVEvent * event) const; 


  Bool_t InitializeCutsFromCutString(const TString analysisCutSelection);
  void SelectCollisionCandidates(UInt_t offlineTriggerMask = AliVEvent::kMB) {fOfflineTriggerMask = offlineTriggerMask;}
  void FillElectonLabelArray(AliAODConversionPhoton* photon, Int_t nV0);
  void SetAcceptedHeader(TList *HeaderList){fHeaderList = HeaderList;}   
  TString *GetFoundHeader(){return fGeneratorNames;}

  Int_t GetEventQuality(){return fEventQuality;}
  Bool_t GetIsFromPileup(){return fRemovePileUp;}
  
  AliConversionCuts(const char *name="V0Cuts", const char * title="V0 Cuts");
  AliConversionCuts(const AliConversionCuts&);
  AliConversionCuts& operator=(const AliConversionCuts&);

  virtual ~AliConversionCuts();                            //virtual destructor

  static AliConversionCuts * GetStandardCuts2010PbPb();
  static AliConversionCuts * GetStandardCuts2010pp();

  virtual Bool_t IsSelected(TObject* /*obj*/){return kTRUE;}
  virtual Bool_t IsSelected(TList* /*list*/) {return kTRUE;}

  TString GetCutNumber();
  
  void GetCentralityRange(Double_t range[2]){range[0]=10*fCentralityMin;range[1]=10*fCentralityMax;};
  
  // Cut Selection
  Bool_t EventIsSelected(AliVEvent *fInputEvent, AliVEvent *fMCEvent);
  Int_t IsEventAcceptedByConversionCut(AliConversionCuts *ReaderCuts, AliVEvent *InputEvent, AliMCEvent *MCEvent, Bool_t isHeavyIon);
  Bool_t PhotonIsSelected(AliConversionPhotonBase * photon, AliVEvent  * event);
  Bool_t PhotonIsSelectedMC(TParticle *particle,AliStack *fMCStack,Bool_t checkForConvertedGamma=kTRUE);
  Bool_t ElectronIsSelectedMC(TParticle *particle,AliStack *fMCStack);
  Bool_t TracksAreSelected(AliVTrack * negTrack, AliVTrack * posTrack);
  Bool_t MesonIsSelected(AliAODConversionMother *pi0,Bool_t IsSignal=kTRUE);
  Bool_t MesonIsSelectedMC(TParticle *fMCMother,AliStack *fMCStack, Bool_t bMCDaughtersInAcceptance=kFALSE);

  void InitAODpidUtil(Int_t type);
  Bool_t InitPIDResponse();
  
  void SetPIDResponse(AliPIDResponse * pidResponse) {fPIDResponse = pidResponse;}
  AliPIDResponse * GetPIDResponse() { return fPIDResponse;}
  
  void PrintCuts();

  void InitCutHistograms(TString name="",Bool_t preCut = kTRUE);
  void SetFillCutHistograms(TString name="",Bool_t preCut = kTRUE){if(!fHistograms){InitCutHistograms(name,preCut);};}
  TList *GetCutHistograms(){return fHistograms;}
  void FillPhotonCutIndex(Int_t photoncut){if(hCutIndex)hCutIndex->Fill(photoncut);}
  void SetEtaShift(Double_t etaShift) {
     fEtaShift = etaShift;
     fLineCutZRSlope = tan(2*atan(exp(-fEtaCut + etaShift)));
     if(fEtaCutMin > -0.1)
        fLineCutZRSlopeMin = tan(2*atan(exp(-fEtaCutMin + etaShift)));
     fDoEtaShift = kTRUE;
  }
  Double_t GetEtaShift() {return fEtaShift;}
  Bool_t GetDoEtaShift(){return fDoEtaShift;}
  void DoEtaShift(Bool_t doEtaShift){fDoEtaShift = doEtaShift;}
  void ForceEtaShift(Int_t forcedShift){
     fForceEtaShift = forcedShift;
     if(forcedShift>0)fDoEtaShift = kTRUE;
  }
  Int_t IsEtaShiftForced() {return fForceEtaShift;}

  static AliVTrack * GetTrack(AliVEvent * event, Int_t label);
  static AliESDtrack *GetESDTrack(AliESDEvent * event, Int_t label);
  
  ///Cut functions
  Bool_t SpecificTrackCuts(AliAODTrack * negTrack, AliAODTrack * posTrack,Int_t &cutIndex);
  Bool_t SpecificTrackCuts(AliESDtrack * negTrack, AliESDtrack * posTrack,Int_t &cutIndex);
  Bool_t AcceptanceCuts(AliConversionPhotonBase *photon);
  Bool_t AcceptanceCut(TParticle *particle, TParticle * ePos,TParticle* eNeg);
  Bool_t dEdxCuts(AliVTrack * track);
  Bool_t ArmenterosQtCut(AliConversionPhotonBase *photon);
  Bool_t AsymmetryCut(AliConversionPhotonBase *photon,AliVEvent *event);
  Bool_t PIDProbabilityCut(AliConversionPhotonBase *photon, AliVEvent * event);
  Bool_t SelectV0Finder(Bool_t onfly){
     if(onfly == fUseOnFlyV0Finder) return kTRUE;
     else return kFALSE;
  }
  Bool_t PhotonCuts(AliConversionPhotonBase *photon,AliVEvent *event);
  Bool_t CorrectedTPCClusterCut(AliConversionPhotonBase *photon, AliVEvent * event);
  Bool_t PsiPairCut(const AliConversionPhotonBase * photon) const;
  Bool_t CosinePAngleCut(const AliConversionPhotonBase * photon, AliVEvent * event) const;
  Bool_t RejectSharedElectronV0s(AliAODConversionPhoton* photon, Int_t nV0, Int_t nV0s);
  Bool_t RejectToCloseV0s(AliAODConversionPhoton* photon, TList *photons, Int_t nV0);
  Int_t IsParticleFromBGEvent(Int_t index, AliStack *MCStack);
  void GetNotRejectedParticles(Int_t rejection, TList *HeaderList, AliMCEvent *MCEvent);

  // Event Cuts
  Bool_t IsCentralitySelected(AliVEvent *fInputEvent, AliVEvent *fMCEvent = NULL);
  Double_t GetCentrality(AliVEvent *event);
  Int_t GetNumberOfContributorsVtx(AliVEvent *event);
  Bool_t VertexZCut(AliVEvent *fInputEvent);
  Bool_t IsTriggerSelected();
  Bool_t HasV0AND(){return fHasV0AND;}
  Bool_t IsSDDFired(){return fIsSDDFired;}
  Int_t IsSpecialTrigger(){return fSpecialTrigger;}

  // Set Individual Cuts
  Bool_t SetRCut(Int_t RCut);
  Bool_t SetV0Finder(Int_t v0FinderType);
  Bool_t SetChi2GammaCut(Int_t chi2GammaCut);
  Bool_t SetTPCdEdxCutPionLine(Int_t pidedxSigmaCut);
  Bool_t SetTPCdEdxCutElectronLine(Int_t ededxSigmaCut);
  Bool_t SetSinglePtCut(Int_t singlePtCut);
  Bool_t SetTPCClusterCut(Int_t clsTPCCut);
  Bool_t SetEtaCut(Int_t etaCut);
  Bool_t SetMinMomPiondEdxCut(Int_t piMinMomdedxSigmaCut);
  Bool_t SetMaxMomPiondEdxCut(Int_t piMaxMomdedxSigmaCut);
  Bool_t SetLowPRejectionCuts(Int_t LowPRejectionSigmaCut);
  Bool_t SetQtMaxCut(Int_t QtMaxCut);
  Bool_t SetTOFElectronPIDCut(Int_t TOFelectronPID);
  Bool_t SetTRDElectronCut(Int_t TRDElectronCut);
  Bool_t SetCentralityMin(Int_t useCentrality);
  Bool_t SetIsHeavyIon(Int_t isHeavyIon);
  Bool_t SetCentralityMax(Int_t centralityBin);
  Bool_t SetPhotonAsymmetryCut(Int_t doPhotonAsymmetryCut);
  Bool_t SetRemovePileUp(Int_t removePileUp);  
  Bool_t SetMultiplicityMethod(Int_t multiplicityMethod);
  Int_t SetSelectSpecialTrigger(Int_t selectSpecialTrigger);
  Bool_t SetCosPAngleCut(Int_t cosCut);
  Bool_t SetPsiPairCut(Int_t psiCut);
  Bool_t SetSharedElectronCut(Int_t sharedElec);
  Bool_t SetToCloseV0sCut(Int_t toClose);
  Bool_t SetRejectExtraSignalsCut(Int_t extraSignal);


  // Request Flags

  Int_t IsHeavyIon(){return fIsHeavyIon;}
  Int_t GetFirstTPCRow(Double_t radius);
  Float_t GetWeightForMeson(TString period, Int_t index, AliStack *MCStack);

  Bool_t UseElecSharingCut(){return fDoSharedElecCut;}
  Bool_t UseToCloseV0sCut(){return fDoToCloseV0sCut;}
  Int_t GetMultiplicityMethod(){return fMultiplicityMethod;}
  Double_t GetEtaCut(){return fEtaCut;}
  Int_t GetSignalRejection(){return fRejectExtraSignals;}
  Int_t GetNAcceptedHeaders(){return fnHeaders; }
  TString * GetAcceptedHeaderNames(){return fGeneratorNames;}
  Int_t * GetAcceptedHeaderStart(){return fNotRejectedStart;}
  Int_t * GetAcceptedHeaderEnd(){return fNotRejectedEnd;}
  TList* GetAcceptedHeader(){return fHeaderList;}
  
  
  protected:
  TList *fHistograms;
  TList *fHeaderList;
  AliPIDResponse *fPIDResponse;


  Int_t fEventQuality; // EventQuality
  //cuts
  Double_t fMaxR; //r cut
  Double_t fMinR; //r cut
  Double_t fEtaCut; //eta cut
  Double_t fEtaCutMin; //eta cut
  Double_t fPtCut; // pt cut
  Double_t fSinglePtCut; // pt cut for electron/positron
  Double_t fMaxZ; //z cut
  Double_t fMinClsTPC; // minimum clusters in the TPC
  Double_t fMinClsTPCToF; // minimum clusters to findable clusters
  Double_t fLineCutZRSlope; //linecut
  Double_t fLineCutZValue; //linecut
  Double_t fLineCutZRSlopeMin; //linecut
  Double_t fLineCutZValueMin; //linecut
  Double_t fChi2CutConversion; //chi2cut
  Double_t fPIDProbabilityCutNegativeParticle;
  Double_t fPIDProbabilityCutPositiveParticle;
  Bool_t   fDodEdxSigmaCut; // flag to use the dEdxCut based on sigmas
  Bool_t   fDoTOFsigmaCut; // flag to use TOF pid cut RRnewTOF
  Double_t fPIDTRDEfficiency; // required electron efficiency for TRD PID
  Bool_t   fDoTRDPID; // flag to use TRD pid
  Double_t fPIDnSigmaAboveElectronLine; // sigma cut
  Double_t fPIDnSigmaBelowElectronLine; // sigma cut
  Double_t fTofPIDnSigmaAboveElectronLine; // sigma cut RRnewTOF
  Double_t fTofPIDnSigmaBelowElectronLine; // sigma cut RRnewTOF 
  Double_t fPIDnSigmaAbovePionLine;     // sigma cut
  Double_t fPIDnSigmaAbovePionLineHighPt;     // sigma cut
  Double_t fPIDMinPnSigmaAbovePionLine; // sigma cut
  Double_t fPIDMaxPnSigmaAbovePionLine; // sigma cut
  Double_t fDoKaonRejectionLowP;   // Kaon rejection at low p
  Double_t fDoProtonRejectionLowP; // Proton rejection at low p
  Double_t fDoPionRejectionLowP;   // Pion rejection at low p
  Double_t fPIDnSigmaAtLowPAroundKaonLine; // sigma cut
  Double_t fPIDnSigmaAtLowPAroundProtonLine; // sigma cut
  Double_t fPIDnSigmaAtLowPAroundPionLine; // sigma cut
  Double_t fPIDMinPKaonRejectionLowP; // Momentum limit to apply kaon rejection
  Double_t fPIDMinPProtonRejectionLowP; // Momentum limit to apply proton rejection
  Double_t fPIDMinPPionRejectionLowP; // Momentum limit to apply proton rejection
  Bool_t   fDoQtGammaSelection; // Select gammas using qtMax
  Bool_t   fDoHighPtQtGammaSelection; // RRnew Select gammas using qtMax for high pT
  Double_t fQtMax; // Maximum Qt from Armenteros to select Gammas
  Double_t fHighPtQtMax; // RRnew Maximum Qt for High pT from Armenteros to select Gammas
  Double_t fPtBorderForQt; // RRnew 
  Double_t fXVertexCut; //vertex cut
  Double_t fYVertexCut; //vertex cut
  Double_t fZVertexCut; // vertexcut
  Double_t fNSigmaMass; //nsigma cut
  Bool_t fUseEtaMinCut; //flag
  Bool_t fUseOnFlyV0Finder; //flag
  Bool_t   fDoPhotonAsymmetryCut; // flag to use the PhotonAsymetryCut
  Double_t fMinPPhotonAsymmetryCut; // Min Momentum for Asymmetry Cut
  Double_t fMinPhotonAsymmetry;  // Asymmetry Cut
  Int_t  fIsHeavyIon;               // flag for heavy ion
  Int_t fDetectorCentrality;	// centrality detecotor V0M or CL1
  Int_t fModCentralityClass; // allows to select smaller centrality classes
  Double_t fMaxVertexZ;    // max z offset of vertex
  Int_t fCentralityMin;  // centrality selection lower bin value
  Int_t fCentralityMax; // centrality selection upper bin value
  Bool_t fUseCorrectedTPCClsInfo; // flag to use corrected tpc cl info
  Bool_t fUseTOFpid; // flag to use tof pid
  Int_t fMultiplicityMethod; // selected multiplicity method
  Int_t fSpecialTrigger; // flag
  Bool_t fRemovePileUp; //flag
  Float_t fOpeningAngle; // min opening angle for meson
  Float_t fPsiPairCut;
  Float_t fCosPAngleCut;
  Bool_t fDoToCloseV0sCut; //
  Int_t fRejectExtraSignals;//
  Double_t fminV0Dist; //
  Bool_t fDoSharedElecCut; //
  UInt_t fOfflineTriggerMask;   //  Task processes collision candidates only
  Bool_t fHasV0AND; // V0AND Offline Trigger
  Bool_t fIsSDDFired; // SDD FIRED to select with SDD events
  TRandom3 fRandom; //
  Int_t fElectronArraySize; // Size of electron array
  Int_t *fElectronLabelArray; //[fElectronArraySize]
  Float_t fConversionPointXArray; // Array with conversion Point x
  Float_t fConversionPointYArray; // Array with conversion Point y
  Float_t fConversionPointZArray; // Array with conversion Point z
  Int_t fnHeaders; // Number of Headers
  Int_t *fNotRejectedStart; //[fnHeaders]
  Int_t *fNotRejectedEnd; //[fnHeaders]
  TString *fGeneratorNames; //[fnHeaders]
  TObjString *fCutString; // cut number used for analysis
  AliAnalysisUtils *fUtils;
  Double_t fEtaShift;
  Bool_t fDoEtaShift;
  Int_t  fForceEtaShift;

  // Histograms
  TH1F *hdEdxCuts;  // bookkeeping for dEdx cuts
  TH2F *hTPCdEdxbefore; // TPC dEdx before cuts
  TH2F *hTPCdEdxafter; // TPC dEdx after cuts
  TH2F *hTPCdEdxSigbefore; // TPC Sigma dEdx before cuts
  TH2F *hTPCdEdxSigafter; // TPC Sigm dEdx after cuts
  TH2F *hTOFbefore; // TOF before cuts
  TH2F *hTOFSigbefore; // TOF Sigma before cuts
  TH2F *hTOFSigafter; // TOF Sigma after cuts
  TH1F *hTrackCuts; // bookkeeping for track cuts
  TH1F *hPhotonCuts; // bookkeeping for photon specific cuts
  TH1F *hInvMassbefore; // e+e- inv mass distribution before cuts
  TH2F *hArmenterosbefore; // armenteros podolanski plot before cuts
  TH1F *hInvMassafter; // e+e- inv mass distribution after cuts
  TH2F *hArmenterosafter;  // armenteros podolanski plot after cuts
  TH1F *hAcceptanceCuts; // bookkeeping for acceptance cuts
  TH1F *hCutIndex; // bookkeeping for cuts
  TH1F *hV0EventCuts; // bookkeeping for event selection cuts
  TH1F *hCentrality; // centrality distribution for selected events
  TH1F *hVertexZ; // vertex z distribution for selected events
  TH1F *hTriggerClass; //fired offline trigger class
  TH1F *hTriggerClassSelected; //selected fired offline trigger class
private:

  ClassDef(AliConversionCuts,3)
};


inline void AliConversionCuts::InitAODpidUtil(Int_t type) {
  if (!fPIDResponse) fPIDResponse = new AliAODpidUtil();
  Double_t alephParameters[5];
  // simulation
  alephParameters[0] = 2.15898e+00/50.;
  alephParameters[1] = 1.75295e+01;
  alephParameters[2] = 3.40030e-09;
  alephParameters[3] = 1.96178e+00;
  alephParameters[4] = 3.91720e+00;
  fPIDResponse->GetTOFResponse().SetTimeResolution(80.);
  
  // data
  if (type==1){
    alephParameters[0] = 0.0283086/0.97;
    alephParameters[1] = 2.63394e+01;
    alephParameters[2] = 5.04114e-11;
    alephParameters[3] = 2.12543e+00;
    alephParameters[4] = 4.88663e+00;
    fPIDResponse->GetTOFResponse().SetTimeResolution(130.);
    fPIDResponse->GetTPCResponse().SetMip(50.);
  }
  
  fPIDResponse->GetTPCResponse().SetBetheBlochParameters(
    alephParameters[0],alephParameters[1],alephParameters[2],
    alephParameters[3],alephParameters[4]);
  
  fPIDResponse->GetTPCResponse().SetSigma(3.79301e-03, 2.21280e+04);
}


#endif
