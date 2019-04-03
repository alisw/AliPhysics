#ifndef ALIDALITZELECTRONCUTS_H
#define ALIDALITZELECTRONCUTS_H

// Class handling all kinds of selection cuts for electrons

// Authors: Svein Lindal, Daniel Lohner												*


#include "AliAODpidUtil.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliMCEvent.h"
#include "AliAnalysisCuts.h"
#include "AliESDtrackCuts.h"
#include "TH1F.h"
#include "AliDalitzData.h"
#include "AliDalitzAODESD.h"
#include "AliDalitzEventMC.h"
#include "AliDalitzAODESDMC.h"

class AliESDEvent;
class AliAODEvent;
class AliConversionPhotonBase;
class AliKFVertex;
class AliKFParticle;
class TH1F;
class TH2F;
class AliPIDResponse;
class AliAnalysisCuts;
class iostream;
class TList;
class AliAnalysisManager;


using namespace std;

class AliDalitzElectronCuts : public AliAnalysisCuts {

 public:


  enum cutIds {
    kMaxChi2TPCConstrainedGlobal,
    kededxSigmaITSCut,
    kededxSigmaTPCCut,
    kpidedxSigmaTPCCut,
    kpiMinMomdedxSigmaTPCCut,
    kpiMaxMomdedxSigmaTPCCut,
    kLowPRejectionSigmaCut,
    kTOFelectronPID,
    kclsITSCut,
    kclsTPCCut,
    ketaCut,
    kPsiPair,
    kRejectSharedElecGamma,
    kMaxChi2PerClusterTPC,
    kMaxChi2PerClusterITS,
    kptCut,
    kDCACut,
    kmassCut,
    kWeights,
    kuseVPhotonMCPSmearing,
    kNCuts
  };


 enum electronCuts {
    kElectronIn=0,
    kNoTracks,
    kTrackCuts,
    kdEdxCuts,
    kElectronOut
 };


  Bool_t SetCutIds(TString cutString);
  Int_t fCuts[kNCuts];
  Bool_t SetCut(cutIds cutID, Int_t cut);
  Bool_t UpdateCutString();
  static const char * fgkCutNames[kNCuts];


  Bool_t InitializeCutsFromCutString(const TString analysisCutSelection);


  AliDalitzElectronCuts(const char *name="ElectronCuts", const char * title="Electron Cuts");
  virtual ~AliDalitzElectronCuts();                            //virtual destructor

  virtual Bool_t IsSelected(TObject* /*obj*/){return kTRUE;}
  virtual Bool_t IsSelected(TList* /*list*/) {return kTRUE;}

  TString GetCutNumber();

    // Cut Selection
  Bool_t ElectronIsSelectedMC(Int_t labelParticle,AliMCEvent *mcEvent);
  Bool_t ElectronIsSelectedMC(Int_t labelParticle,AliMCEvent *mcEvent, AliDalitzEventMC *mcAODESDEvent);
  Bool_t TrackIsSelected(AliDalitzAODESD* lTrack);
  Bool_t ElectronIsSelected(AliESDtrack* lTrack);
  Bool_t ElectronIsSelected(AliDalitzAODESD* lTrack);
  void InitAODpidUtil(Int_t type);
  static AliDalitzElectronCuts * GetStandardCuts2010PbPb();
  static AliDalitzElectronCuts * GetStandardCuts2010pp();
  Bool_t InitPIDResponse();

  void SetPIDResponse(AliPIDResponse * pidResponse) {fPIDResponse = pidResponse;}
  AliPIDResponse * GetPIDResponse() { return fPIDResponse;}

  void PrintCuts();
  Bool_t AcceptedAODESDTrack(AliDalitzAODESD* aliaodtrack);
  void InitCutHistograms(TString name="",Bool_t preCut = kTRUE,TString cutName="");
  void SetFillCutHistograms(TString name="",Bool_t preCut = kTRUE,TString cutName=""){if(!fHistograms){InitCutHistograms(name,preCut,cutName);};}
  TList *GetCutHistograms(){return fHistograms;}

  static AliVTrack * GetTrack(AliVEvent * event, Int_t label);



  ///Cut functions
  Bool_t dEdxCuts(AliVTrack * track);
  Bool_t PIDProbabilityCut(AliConversionPhotonBase *photon, AliVEvent * event);
  Bool_t RejectSharedElecGamma(TList *photons, Int_t indexEle);
  Bool_t IsFromGammaConversion( Double_t psiPair, Double_t deltaPhi );
  Bool_t MassCut(Double_t pi0CandidatePt,Double_t vphotonCandidateMass);

  // Event Cuts

  //Double_t GetPsiPair( const AliESDtrack *trackPos, const AliESDtrack *trackNeg );
  Bool_t SetMaxChi2TPCConstrainedGlobal(Int_t maxChi2);
  Bool_t SetTPCdEdxCutPionLine(Int_t pidedxSigmaCut);
  Bool_t SetTPCdEdxCutElectronLine(Int_t ededxSigmaCut);
  Bool_t SetITSdEdxCutElectronLine(Int_t ededxSigmaCut);
  Bool_t SetMinMomPiondEdxTPCCut(Int_t piMomdedxSigmaCut);
  Bool_t SetMaxMomPiondEdxTPCCut(Int_t piMomdedxSigmaCut);
  Bool_t SetITSClusterCut(Int_t clsITSCut);
  Bool_t SetTPCClusterCut(Int_t clsTPCCut);
  Bool_t SetEtaCut(Int_t etaCut);
  Bool_t SetPtCut(Int_t ptCut);
  Bool_t SetDCACut(Int_t dcaCut);
  Bool_t SetMinMomPiondEdxCut(Int_t piMinMomdedxSigmaCut);
  Bool_t SetMaxMomPiondEdxCut(Int_t piMaxMomdedxSigmaCut);
  Bool_t SetLowPRejectionCuts(Int_t LowPRejectionSigmaCut);
  Bool_t SetTOFElectronPIDCut(Int_t TOFelectronPID);
  Bool_t SetPsiPairCut(Int_t psiCut);
  Bool_t SetRejectSharedElecGamma(Int_t RCut);
  Bool_t SetMaxChi2PerClusterTPC(Int_t maxChi2);
  Bool_t SetMaxChi2PerClusterITS(Int_t maxChi2);
  Bool_t SetMassCut(Int_t massCut);
  Bool_t SetDoWeights(Int_t opc);
  Bool_t SetUseVPhotonMCPmearing(Int_t useMCPSmearing);
  void SetUseCrossedRows(Bool_t opc){fUseCrossedRows = opc;};

  // Request Flags

  Double_t GetEtaCut(){ return  fEtaCut;}
  Double_t GetPsiPairCut(){ return fPsiPairCut; }
  Double_t DoRejectSharedElecGamma(){ return fDoRejectSharedElecGamma;}
  Double_t DoPsiPairCut(){return fDoPsiPairCut;}
  Double_t GetNFindableClustersTPC(AliDalitzAODESD* lTrack);
  Bool_t   DoMassCut(){return  fDoMassCut;}
  Bool_t   DoMassMinCut(){return fDoMassMinCut;}
  Double_t GetMassCutLowPt(){return fMassCutLowPt;}
  Double_t GetMassCutHighPt(){return fMassCutHighPt;}
  Double_t GetPtMinMassCut(){return fMassCutPtMin;}
  Double_t GetMassMinCut(){return fMassMinCut;}
  Bool_t   DoWeights(){return fDoWeights;}
  Bool_t   GetUseVPhotonMCPmearing(){return fUseVPhotonMCPSmearing; }
  Bool_t   GetUseElectronMCSmearing(){ return fUseElectronMCPSmearing;}



  protected:

  TList *fHistograms;
  AliPIDResponse *fPIDResponse;
  AliESDtrackCuts *fesdTrackCuts;

  Double_t fEtaCut; //eta cutç
  Bool_t   fDoEtaCut;
  Double_t fPtMinCut;
  Double_t fPtMaxCut;
  Double_t fRadiusCut; // radius cut
  Double_t fPsiPairCut;
  Double_t fDeltaPhiCutMin;
  Double_t fDeltaPhiCutMax;
  Double_t fMaxDCAVertexz;
  Double_t fMaxDCAVertexxy;
  TString fDCAVertexPt;//Conversion from coordenates to momentum.
  Double_t fMinClsTPC; // minimum clusters in the TPC
  Double_t fMinClsTPCToF; // minimum clusters to findable clusters
  Bool_t   fDodEdxSigmaITSCut; // flag to use the dEdxCut ITS based on sigmas
  Bool_t   fDodEdxSigmaTPCCut; // flag to use the dEdxCut TPC based on sigmas
  Bool_t   fDoTOFsigmaCut; // flag to use TOF pid cut RRnewTOF
  Bool_t   fDoRejectSharedElecGamma; //Reject electrons from the gammas with Radius < RadiusCut
  Bool_t   fDoPsiPairCut; // PsiPair Cut
  Double_t fPIDnSigmaAboveElectronLineITS; // sigma cut
  Double_t fPIDnSigmaBelowElectronLineITS; // sigma cut
  Double_t fPIDnSigmaAboveElectronLineTPC; // sigma cut
  Double_t fPIDnSigmaBelowElectronLineTPC; // sigma cut
  Double_t fPIDnSigmaAbovePionLineTPC;
  Double_t fPIDnSigmaAbovePionLineTPCHighPt;
  Double_t fTofPIDnSigmaAboveElectronLine; // sigma cut RRnewTOF
  Double_t fTofPIDnSigmaBelowElectronLine; // sigma cut RRnewTOF
  Double_t fPIDMinPnSigmaAbovePionLineTPC; // sigma cut
  Double_t fPIDMaxPnSigmaAbovePionLineTPC; // sigma cut
  Double_t fDoKaonRejectionLowP;   // Kaon rejection at low p
  Double_t fDoProtonRejectionLowP; // Proton rejection at low p
  Double_t fDoPionRejectionLowP;   // Pion rejection at low p
  Double_t fPIDnSigmaAtLowPAroundKaonLine; // sigma cut
  Double_t fPIDnSigmaAtLowPAroundProtonLine; // sigma cut
  Double_t fPIDnSigmaAtLowPAroundPionLine; // sigma cut
  Double_t fPIDMinPKaonRejectionLowP; // Momentum limit to apply kaon rejection
  Double_t fPIDMinPProtonRejectionLowP; // Momentum limit to apply proton rejection
  Double_t fPIDMinPPionRejectionLowP; // Momentum limit to apply proton rejection

  Bool_t   fUseCorrectedTPCClsInfo; // flag to use corrected tpc cl info
  Bool_t   fUseCrossedRows;  //UseCrossedRows 2011
  Int_t    fITSCut;
  Bool_t   fUseTOFpid; // flag to use tof pid
  Bool_t   fRequireTOF; //flg to analyze only tracks with TOF signal
  Bool_t   fDoMassCut;
  Bool_t   fDoMassMinCut;
  Double_t fMassCutLowPt;
  Double_t fMassCutHighPt;
  Double_t fMassCutPtMin;
  Double_t fMassMinCut;
  Bool_t   fDoWeights;
  Bool_t   fUseVPhotonMCPSmearing;
  Bool_t   fUseElectronMCPSmearing;


  // Histograms
  TObjString *fCutString; // cut number used for analysis
  TString fCutStringRead;
  TH1F *hCutIndex; // bookkeeping for cuts
  TH1F *hdEdxCuts;  // bookkeeping for dEdx cuts
  TH2F *hITSdEdxbefore; // ITS dEdx before cuts
  TH2F *hITSdEdxafter;
  TH2F *hTPCdEdxbefore; // TPC dEdx before cuts
  TH2F *hTPCdEdxafter; // TPC dEdx after cuts
  TH2F *hTPCdEdxSignalbefore; //TPC dEdx signal before
  TH2F *hTPCdEdxSignalafter; //TPC dEdx signal  after
  TH2F *hTOFbefore; // TOF after cuts
  TH2F *hTOFafter; // TOF after cuts
  TH2F *hTrackDCAxyPtbefore;
  TH2F *hTrackDCAxyPtafter;
  TH2F *hTrackDCAzPtbefore;
  TH2F *hTrackDCAzPtafter;
  TH2F *hTrackNFindClsPtTPCbefore;
  TH2F *hTrackNFindClsPtTPCafter;
  TH1F *hTrackPosEtabeforeDedx;
  TH1F *hTrackNegEtabeforeDedx;
  TH1F *hTrackPosEtaafterDedx;
  TH1F *hTrackNegEtaafterDedx;



private:

  AliDalitzElectronCuts(const AliDalitzElectronCuts&); // not implemented
  AliDalitzElectronCuts& operator=(const AliDalitzElectronCuts&); // not implemented


  ClassDef(AliDalitzElectronCuts,3)
};


inline void AliDalitzElectronCuts::InitAODpidUtil(Int_t type) {
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
