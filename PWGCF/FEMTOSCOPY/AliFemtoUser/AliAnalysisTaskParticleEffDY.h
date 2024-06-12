#ifndef ALIANALYSISTASKPARTICLEEFFDY
#define ALIANALYSISTASKPARTICLEEFFDY



#define MULTBINS 1
#define PARTTYPES 6

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliEventCuts.h"
#include "THnSparse.h"
#include "AliPIDResponse.h"
#include "AliAODpidUtil.h"
class AliAnalysisUtils;

class AliAnalysisTaskParticleEffDY :public AliAnalysisTaskSE{
 public:

  enum EventMult {kRefMult=0, kV0M=1, kV0A=2};
  typedef enum EventMult EstEventMult;

  enum PidMethod {kNSigma=0, kNSigmaNoDoubleCounting=1, kExclusivePID=2, kExclusivePIDDiffRejection=3};
  typedef enum PidMethod PidMethod;

 AliAnalysisTaskParticleEffDY();

  AliAnalysisTaskParticleEffDY(TString name, int pidMethod=3, int filterbit=96); // default constructor
  virtual ~AliAnalysisTaskParticleEffDY(); // default destructor
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
  AliAnalysisTaskParticleEffDY(const AliAnalysisTaskParticleEffDY &); // copy constructor
  AliAnalysisTaskParticleEffDY &operator=(const AliAnalysisTaskParticleEffDY &); // operator=
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
  TH1F *fHistP[4];
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

  TH1D *fYLambda;
  TH1D *fPtLambda;
  TH1D *fYAntiLambda;
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
  ClassDef(AliAnalysisTaskParticleEffDY, 0);
  Bool_t fTrackPileUpRemoval;
  Bool_t fV0PileUpRemoval;
  
  

};
#endif /* ALIANALYSISTASKPARTICLEEFFICIENCYDY */
