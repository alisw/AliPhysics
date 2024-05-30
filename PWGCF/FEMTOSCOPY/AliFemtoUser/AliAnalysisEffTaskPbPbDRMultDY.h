#ifndef ALIANALYSISEFFTASKPBPBDRMULTDY
#define ALIANALYSISEFFTASKPBPBDRMULTDY



#define MULTBINS 5
#define PARTTYPES 6

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliEventCuts.h"
#include "THnSparse.h"
#include "AliPIDResponse.h"
#include "AliAODpidUtil.h"
class AliAnalysisUtils;

class AliAnalysisEffTaskPbPbDRMultDY :public AliAnalysisTaskSE{
public:

  enum EventMult {kRefMult=0, kV0M=1, kV0A=2};
  typedef enum EventMult EstEventMult;

  enum PidMethod {kNSigma=0, kNSigmaNoDoubleCounting=1, kExclusivePID=2, kExclusivePIDDiffRejection=3};
  typedef enum PidMethod PidMethod;

AliAnalysisEffTaskPbPbDRMultDY();
  AliAnalysisEffTaskPbPbDRMultDY(TString name, int pidMethod=3, int filterbit=128); // default constructor
  virtual ~AliAnalysisEffTaskPbPbDRMultDY(); // default destructor
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
  AliAnalysisEffTaskPbPbDRMultDY(const AliAnalysisEffTaskPbPbDRMultDY &); // copy constructor
  AliAnalysisEffTaskPbPbDRMultDY &operator=(const AliAnalysisEffTaskPbPbDRMultDY &); // operator=
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

  //TObjArray *recoParticleArray;
  TH2F *fGeneratedMCPrimaries[MULTBINS*PARTTYPES][2];
  TH2F *fMCPrimariesThatAreReconstructed[MULTBINS*PARTTYPES][2];

  TH2F *fMCPrimariesThatAreReconstructedNoNsigma[MULTBINS*PARTTYPES][2];
  TH2F *fReconstructedAfterCuts[MULTBINS*PARTTYPES][2];
  TH2F *fReconstructedNotPrimaries[MULTBINS*PARTTYPES][2];
  TH2F *fReconstructedPrimaries[MULTBINS*PARTTYPES][2];
  TH2F *fContamination[MULTBINS*PARTTYPES][2];
  TH2F *fMisidentification[MULTBINS*PARTTYPES][2];

  
  //******************************************

  double fDCAtoPrimVtx;
  Bool_t fIfAliEventCuts;
  int    fFB;
  Bool_t fDCAglobalTrack;
  PidMethod    fPidMethod; //PID method
  EstEventMult   fEstEventMult;  // Type of the event multiplicity estimator


  //******************************************

  AliPIDResponse *fpidResponse;
  AliAODpidUtil  *fAODpidUtil;
  AliEventCuts   *fEventCuts;
  ClassDef(AliAnalysisEffTaskPbPbDRMultDY, 0);
  Bool_t fTrackPileUpRemoval;
  Bool_t fV0PileUpRemoval;
  
  

};

#endif /* ALIANALYSISTASKPARTICLEEFFICIENCYDY */

