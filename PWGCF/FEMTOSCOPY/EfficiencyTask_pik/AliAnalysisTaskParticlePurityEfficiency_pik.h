#ifndef ALIANALYSISTASKPARTICLEPURITYEFFICIENCY_PIK
#define ALIANALYSISTASKPARTICLEPURITYEFFICIENCY_PIK



#define CENTRBINS 9
#define PARTTYPES 4

#include "AliAnalysisTaskSE.h"
//#include "/home/przemcio/alice/sw/ubuntu1404_x86-64/AliPhysics/0-1/include/AliAnalysisUtils.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"

class AliAnalysisUtils;

class AliAnalysisTaskParticlePurityEfficiency_pik :public AliAnalysisTaskSE{
 public:
 AliAnalysisTaskParticlePurityEfficiency_pik() : AliAnalysisTaskSE(), centrality(0), fHistoList(0), fHistEv(0), fpidResponse(0), fAODpidUtil(0)
    {
      for(Int_t i = 0; i < CENTRBINS*PARTTYPES; i++) {
	for(Int_t chg = 0; chg < 2; chg++) {
	  fGeneratedMCPrimaries[i][chg] = NULL;
	  fMCPrimariesThatAreReconstructed[i][chg] = NULL;
	  fMCPrimariesThatAreReconstructedNoNsigma[i][chg] = NULL;
	  fReconstructedAfterCuts[i][chg] = NULL;
	  fReconstructedNotPrimaries[i][chg] = NULL;
	  fReconstructedPrimaries[i][chg] = NULL;
	  fContamination[i][chg] = NULL;
	}
      }
    }
  /*AliAnalysisTaskParticleEff() : AliAnalysisTaskSE(),centrality(0), fHistoList(0),  fHistEv(0), fpidResponse(0), fAODpidUtil(0)
    {

      for(Int_t i = 0; i < MULTBINS*PARTTYPES; i++)  {
	for(Int_t chg=0;chg<2;chg++){
	  fGeneratedMCPrimaries[i][chg] = NULL;
	  fMCPrimariesThatAreReconstructed[i][chg] = NULL;
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
      }*/

  AliAnalysisTaskParticlePurityEfficiency_pik(const Char_t *partName); // default constructor
  virtual ~AliAnalysisTaskParticlePurityEfficiency_pik(); // default destructor
  virtual void UserCreateOutputObjects(); // user create output objects
  virtual void UserExec(Option_t *option); // user exec
  //void Terminate(Option_t *option);
  
 private:
  AliAnalysisTaskParticlePurityEfficiency_pik(const AliAnalysisTaskParticlePurityEfficiency_pik &); // copy constructor
  AliAnalysisTaskParticlePurityEfficiency_pik &operator=(const AliAnalysisTaskParticlePurityEfficiency_pik &); // operator=
  AliCentrality *centrality;
  TList *fHistoList;
  TH1F *fHistEv;
 
  TH1F *fHistQA[11];
  TH2F *fHistQA2D[11];
  TH2F *fHistQAPID[5][PARTTYPES][2];
  TH2F *fHistQAPIDFail[5][PARTTYPES][2];
  TH1F* fHistEvCuts[CENTRBINS];
  TH2F *fGeneratedMCPrimaries[CENTRBINS*PARTTYPES][2];
  TH2F *fMCPrimariesThatAreReconstructed[CENTRBINS*PARTTYPES][2];
  TH2F *fMCPrimariesThatAreReconstructedNoNsigma[CENTRBINS*PARTTYPES][2];
  TH2F *fReconstructedAfterCuts[CENTRBINS*PARTTYPES][2];
  TH2F *fReconstructedNotPrimaries[CENTRBINS*PARTTYPES][2];
  TH2F *fReconstructedPrimaries[CENTRBINS*PARTTYPES][2];
  TH2F *fContamination[CENTRBINS*PARTTYPES][2];
  TH3F *fMisidentification[CENTRBINS*PARTTYPES][2];
  TProfile *fProfilePTrueReconstructed[4];
  TProfile *fProfileThetaTrueReconstructed[4];
  TProfile *fProfilePhiTrueReconstructed[4];
  TH2D *fHist2DPTrueReconstructed[4];
  TH2D *fHist2DThetaTrueReconstructed[4];
  TH2D *fHist2DPhiTrueReconstructed[4];
  
  /*TH1F *fPrimVsDCA[MULTBINS*PARTTYPES][2];
  TH1F *fSecWeakVsDCA[MULTBINS*PARTTYPES][2];
  TH1F *fSecMatVsDCA[MULTBINS*PARTTYPES][2];
  TH1F *fFakeVsDCA[MULTBINS*PARTTYPES][2];

  TH1F *fPrimVsCosPointingAngle[MULTBINS*PARTTYPES][2];
  TH1F *fSecWeakVsCosPointingAngle[MULTBINS*PARTTYPES][2];
  TH1F *fSecMatVsCosPointingAngle[MULTBINS*PARTTYPES][2];
  TH1F *fFakeVsCosPointingAngle[MULTBINS*PARTTYPES][2];

  TH1F *fPrimVsDecayRadius[MULTBINS*PARTTYPES][2];
  TH1F *fSecWeakVsDecayRadius[MULTBINS*PARTTYPES][2];
  TH1F *fSecMatVsDecayRadius[MULTBINS*PARTTYPES][2];
  TH1F *fFakeVsDecayRadius[MULTBINS*PARTTYPES][2];


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
  TH2D *fRecPtAntiLambdaMC;*/
  
  
  
  AliPIDResponse *fpidResponse;
  AliAODpidUtil  *fAODpidUtil;
  ClassDef(AliAnalysisTaskParticlePurityEfficiency_pik, 0);

};

#endif /* ALIANALYSISTASKPARTICLEEFFICIENCY */
