#ifndef ALIANALYSISTASKPARTICLEEFF
#define ALIANALYSISTASKPARTICLEEFF



#define MULTBINS 1
#define PARTTYPES 6

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliEventCuts.h"
#include "THnSparse.h"
#include "AliPIDResponse.h"
#include "AliAODpidUtil.h"
class AliAnalysisUtils;

class AliAnalysisTaskParticleEff :public AliAnalysisTaskSE{
 public:

  enum EventMult {kRefMult=0, kV0M=1, kV0A=2};
  typedef enum EventMult EstEventMult;

  enum PidMethod {kNSigma=0, kNSigmaNoDoubleCounting=1, kExclusivePID=2, kExclusivePIDDiffRejection=3};
  typedef enum PidMethod PidMethod;

 AliAnalysisTaskParticleEff() :  AliAnalysisTaskSE(), centrality(0), fHistoList(0),  fMassInvLambdaPass(0),fMassInvAntiLambdaPass(0), fMassInvLambdaFail(0), fMassInvAntiLambdaFail(0),fEtaLambda(0),fPtLambda(0), fEtaAntiLambda(0),fPtAntiLambda(0), fCutsLambda(0), fCutsAntiLambda(0), fTruePtLambdaMC(0), fRecPtLambdaMC(0), fTruePtAntiLambdaMC(0),fRecPtAntiLambdaMC(0), fMassInvXimPass(0),fMassInvXipPass(0), fMassInvXimFail(0), fMassInvXipFail(0),fEtaXim(0),fPtXim(0), fEtaXip(0),fPtXip(0), fCutsXim(0), fCutsXip(0), recoParticleArrayXi(0), fTruePtXimMC(0), fRecPtXimMC(0), fTruePtXipMC(0),fRecPtXipMC(0), fDCAtoPrimVtx(0), fIfAliEventCuts(kFALSE), fFB(96), fPidMethod(kExclusivePIDDiffRejection), fEstEventMult(kRefMult), fIfXiAnalysis(kFALSE), fpidResponse(0), fAODpidUtil(0), fEventCuts(0), fTrackPileUpRemoval(kFALSE), fV0PileUpRemoval(kFALSE)

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
      }
    }

  AliAnalysisTaskParticleEff(TString name, int pidMethod=-1, int filterbit=96); // default constructor
  virtual ~AliAnalysisTaskParticleEff(); // default destructor
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

 private:
  AliAnalysisTaskParticleEff(const AliAnalysisTaskParticleEff &); // copy constructor
  AliAnalysisTaskParticleEff &operator=(const AliAnalysisTaskParticleEff &); // operator=
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
  TH2F *fHistQA2D[3];
  TH2F *fHistQAPID[5][PARTTYPES][2];
  TH2F *fHistQAPIDFail[5][PARTTYPES][2];
  TH1F* fHistEvCuts[MULTBINS];
  TH2F *fHistQALambdas[2];
  TH2F *fOriginLambdas[5][2];
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


  TH1D *fMassInvLambdaPass;
  TH1D *fMassInvAntiLambdaPass;

  TH1D *fMassInvLambdaFail;
  TH1D *fMassInvAntiLambdaFail;

  TH1D *fEtaLambda;
  TH1D *fPtLambda;
  TH1D *fEtaAntiLambda;
  TH1D *fPtAntiLambda;


  TH1D *fCutsLambda;
  TH1D *fCutsAntiLambda;

  TH2D *fTruePtLambdaMC;
  TH2D *fRecPtLambdaMC;
  TH2D *fTruePtAntiLambdaMC;
  TH2D *fRecPtAntiLambdaMC;


  //********************Xi*******************
  TH1D *fMassInvXimPass;
  TH1D *fMassInvXipPass;

  TH1D *fMassInvXimFail;
  TH1D *fMassInvXipFail;

  TH1D *fEtaXim;
  TH1D *fPtXim;
  TH1D *fEtaXip;
  TH1D *fPtXip;

  TH2D *fTruePtXimMC;
  TH2D *fRecPtXimMC;
  TH2D *fTruePtXipMC;
  TH2D *fRecPtXipMC;

  TH1D *fCutsXim;
  TH1D *fCutsXip;
  TObjArray *recoParticleArrayXi;

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
  ClassDef(AliAnalysisTaskParticleEff, 0);
  Bool_t fTrackPileUpRemoval;
  Bool_t fV0PileUpRemoval;

};

#endif /* ALIANALYSISTASKPARTICLEEFFICIENCY */
