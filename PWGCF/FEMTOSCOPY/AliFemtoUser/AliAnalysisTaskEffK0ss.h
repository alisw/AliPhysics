#ifndef AliAnalysisTaskEffK0SS
#define AliAnalysisTaskEffK0SS



#define MULTBINS 1
#define PARTTYPES 6

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliEventCuts.h"
#include "THnSparse.h"
#include "AliPIDResponse.h"
#include "AliAODpidUtil.h"
class AliAnalysisUtils;

class AliAnalysisTaskEffK0ss :public AliAnalysisTaskSE{
 public:

  enum EventMult {kRefMult=0, kV0M=1, kV0A=2};
  typedef enum EventMult EstEventMult;

  enum PidMethod {kNSigma=0, kNSigmaNoDoubleCounting=1, kExclusivePID=2, kExclusivePIDDiffRejection=3};
  typedef enum PidMethod PidMethod;
 
 AliAnalysisTaskEffK0ss(); /*AliAnalysisTaskSE(name), centrality(0), fHistoList(0),  fMassInvK0sPass(0), fMassInvK0sFail(0),fEtaK0s(0),fPtK0s(0), fCutsK0s(0), fTruePtK0sMC(0), fRecPtK0sMC(0), fDCAtoPrimVtx(0), fIfAliEventCuts(kFALSE), fFB(128), fPidMethod(kExclusivePIDDiffRejection),  fEstEventMult(kV0M), fpidResponse(0), fAODpidUtil(0), fEventCuts(0)


    {

      for(Int_t i = 0; i < MULTBINS*PARTTYPES; i++)  {
	for(Int_t chg=0;chg<2;chg++){
	  fGeneratedMCPrimaries[i][chg] = NULL;
	  fMCPrimariesThatAreReconstructed[i][chg] = NULL;
	  fGeneratedMCPrimaries4D[i][chg] = NULL;
	  fMCPrimariesThatAreReconstructed4D[i][chg] = NULL;
	  fMCPrimariesThatAreReconstructedNoNsigma[i][chg] = NULL;
	  fReconstructedAfterCuts[i][chg] = NULL;
	  fReconstructedNotPrimaries[i][chg] = NULL;
	  fReconstructedPrimaries[i][chg] = NULL;
	  fContamination[i][chg] = NULL;
	}
      }

      for ( Int_t i = 0; i < 11; i++) {
	fHistQA[i] = NULL;
	if(i<3) fHistQA2D[i] = NULL;
	if(i<4) fHistP[i]=NULL;
      }
    }*/

AliAnalysisTaskEffK0ss(TString name,
                       int pidMethod,
                       int filterbit,
                       Double_t nsigmaDaught,
                       Double_t cosPA,
                       Double_t window,
                       Double_t dcaDaughtersMax,
                       Double_t dcaPVMin,
                       Bool_t electronRejection,
                       Double_t nsigmaErej,
		               Double_t minDCApos, 
		               Double_t minDCAneg,
		               Double_t maxctau, 
		               Double_t minV0rad,
		               Bool_t ownDCA, 
		               Float_t dcaxy, 
		               Float_t dcaz
);
  virtual ~AliAnalysisTaskEffK0ss(); // default destructor
  virtual void UserCreateOutputObjects(); // user create output objects
  virtual void UserExec(Option_t *option); // user exec
  //void Terminate(Option_t *option);

  void SetFB(int fb);
  void SetPidMethod(PidMethod method);
  void SetPidMethod(int method);
  int GetPidMethod();
  void SetMultMethod(EstEventMult method);
  void SetAliEventCuts(Bool_t ec);
  void SetIfXiAnalysis(Bool_t xi);
  void SetIfTrackPileUp(Bool_t ifTrackPlp);
  void SetV0PileUpRemoval(Bool_t v0PileUpRemoval);
  void AnalyseCascades(int fcent, AliAODEvent* aodEvent, TClonesArray  *arrayMC);
  void SetNsigmaDaughters(Double_t nsigmaDaught);
  void SetMinCosPointingAngle(Double_t cosPA);
  void SetInvMassWindow(Double_t window);
  void SetMaxDCADaughters(Double_t dcaDaughtersMax);
  void SetMinDCAToPrimVtx(Double_t dcaPVMin);
  void SetElectronRejection(Bool_t electronRejection);
  void SetNsigmaElectronRejection(Double_t nsigmaErej);
  void SetMinDcaPosToPrimVertex(Double_t minDCApos);
  void SetMinDcaNegToPrimVertex(Double_t minDCAneg);
  void SetMaxCTauK0s(Double_t maxctau);
  void SetMinV0Radius(Double_t minV0rad);
  void SetUseDCAcuts(Bool_t ownDCA);
  void SetUseDCAcutsxy(Float_t dcaxy);
  void SetUseDCAcutsz(Float_t dcaz);

  bool IsElectronAM(float nsigmaTPCe, float nsigmaTPCPi, float nsigmaTPCK, float nsigmaTPCP);
  bool IsElectronAM1(float nsigmaTPCe);

 private:
  AliAnalysisTaskEffK0ss(const AliAnalysisTaskEffK0ss &); // copy constructor
  AliAnalysisTaskEffK0ss &operator=(const AliAnalysisTaskEffK0ss &); // operator=
  //AliAODEvent *aodEvent;
  AliCentrality *centrality;
  //AliAODTrack *fTpcTracks;
  //AliAODVertex    *vertex;
  //AliAODVertex    *vtxSPD;
  //AliAODMCParticle *MCtrk;
  //AliESDtrackCuts *fTrackCuts; // ESD track cuts
  TList *fHistoList; // histo list
  //TClonesArray *arrayMC;
  TH1F *fHistEv[4];

  TH1F *fHistQA[11];
  TH1F *fHistP[4];
  TH2F *fHistQA2D[5];
  TH2F *fHistQAPID[5][PARTTYPES][2];
  TH2F *fHistQAPIDFail[5][PARTTYPES][2];
  TH1F* fHistEvCuts[4];
  TH2F *fHistQAK0ss[2];
  TH2F *fOriginK0ss[5][2];
  TH2F *fHistQAXi[2];
  TH2F *fOriginXi[5][2];

  //TObjArray *recoParticleArray;
  TH2F *fGeneratedMCPrimaries[MULTBINS*PARTTYPES][2];
  TH2F *fMCPrimariesThatAreReconstructed[MULTBINS*PARTTYPES][2];
  THnSparseF *fMCPrimariesThatAreReconstructed4D[MULTBINS*PARTTYPES][2];
  THnSparseF *fGeneratedMCPrimaries4D[MULTBINS*PARTTYPES][2];

  TH2F *fMCPrimariesThatAreReconstructedNoNsigma[MULTBINS*PARTTYPES][2];
  TH2F *fReconstructedAfterCuts[MULTBINS*PARTTYPES][2];
  TH2F *fReconstructedNotPrimaries[MULTBINS*PARTTYPES][2];
  TH2F *fReconstructedPrimaries[MULTBINS*PARTTYPES][2];
  TH2F *fContamination[MULTBINS*PARTTYPES][2];
  TH2F *fMisidentification[MULTBINS*PARTTYPES][2];

  TH2F *fPrimVsDCA[MULTBINS*PARTTYPES][2];
  TH2F *fSecWeakVsDCA[MULTBINS*PARTTYPES][2];
  TH2F *fSecMatVsDCA[MULTBINS*PARTTYPES][2];
  TH2F *fFakeVsDCA[MULTBINS*PARTTYPES][2];

  TH2F *fPrimVsCosPointingAngle[MULTBINS*PARTTYPES][2];
  TH2F *fSecWeakVsCosPointingAngle[MULTBINS*PARTTYPES][2];
  TH2F *fSecMatVsCosPointingAngle[MULTBINS*PARTTYPES][2];
  TH2F *fFakeVsCosPointingAngle[MULTBINS*PARTTYPES][2];

  TH2F *fPrimVsDecayRadius[MULTBINS*PARTTYPES][2];
  TH2F *fSecWeakVsDecayRadius[MULTBINS*PARTTYPES][2];
  TH2F *fSecMatVsDecayRadius[MULTBINS*PARTTYPES][2];
  TH2F *fFakeVsDecayRadius[MULTBINS*PARTTYPES][2];

  TH2F *fAllVsDCA[MULTBINS*PARTTYPES][2];
  TH2F *fAllVsCosPointingAngle[MULTBINS*PARTTYPES][2];
  TH2F *fAllVsDecayRadius[MULTBINS*PARTTYPES][2];


  TH1D *fMassInvK0sPass;
  TH1D *fMassInvAntiK0sPass;

  TH1D *fMassInvK0sFail;
  TH1D *fMassInvAntiK0sFail;

  TH1D *fMassInvK0sAfterCuts;
  TH2D *fMassInvK0s;
  TH2D *fMassInvK0sPt;

  TH1D *fEtaK0s;
  TH1D *fPtK0s;
  TH1D *fYAntiK0s;
  TH1D *fPtAntiK0s;


  TH1D *fCutsK0s;
  TH1D *fCutsAntiK0s;

  TH2D *fTruePtK0sMC;
  TH2D *fRecPtK0sMC;
  TH2D *fTruePtAntiK0sMC;
  TH2D *fRecPtAntiK0sMC;


  //********************Xi*******************
  TH1D *fMassInvXimPass;
  TH1D *fMassInvXipPass;

  TH1D *fMassInvXimFail;
  TH1D *fMassInvXipFail;

  TH1D *fYXim;
  TH1D *fPtXim;
  TH1D *fYXip;
  TH1D *fPtXip;

  TH2D *fTruePtXimMC;
  TH2D *fRecPtXimMC;
  TH2D *fTruePtXipMC;
  TH2D *fRecPtXipMC;

  TH1D *fCutsXim;
  TH1D *fCutsXip;
  TObjArray *recoParticleArrayXi;
  TH2F* fPIDKch;      // PID histogram for Kaons
  TH2F* fPIDKeCut;    // PID histogram for Kaons with electron rejection

  //******************************************

  double fDCAtoPrimVtx;
  Bool_t fIfAliEventCuts;
  int    fFB;
  PidMethod    fPidMethod; //PID method
  EstEventMult   fEstEventMult;  // Type of the event multiplicity estimator
  Bool_t fIfXiAnalysis;

  //******************************************

  AliPIDResponse *fpidResponse;
  AliAODpidUtil  *fAODpidUtil;
  AliEventCuts   *fEventCuts;
  ClassDef(AliAnalysisTaskEffK0ss, 0);
  Bool_t fTrackPileUpRemoval;
  Bool_t fV0PileUpRemoval;

  Double_t fNsigmaDaughters;
  Double_t fMinCosPointingAngle;
  Double_t fInvMassWindow;
  Double_t fMaxDcaV0Daughters;
  Double_t fMinDCAToPrimVtx;
  Bool_t   fElectronReject;
  Double_t fNsigmaElectronRejection;

  Double_t fMinDcaPosToPrimVertex;
  Double_t fMinDcaNegToPrimVertex;
  Double_t fMaxCTauK0s;
  Double_t fMinV0Radius;

  Bool_t   fUseDcaCuts;
  Float_t  fDcaXYCut;
  Float_t  fDcaZCut;

  // ---- Other settings mentioned in the constructor ----
  Int_t    fCentMin;
  Int_t    fCentMax;
  Double_t fPVzCut;

};

#endif /* ALIANALYSISTASKPARTICLEEFFICIENCYDY */
