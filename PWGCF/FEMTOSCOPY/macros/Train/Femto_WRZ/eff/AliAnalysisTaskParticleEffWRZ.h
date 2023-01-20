#ifndef ALIANALYSISTASKPARTICLEEFFWRZ
#define ALIANALYSISTASKPARTICLEEFFWRZ

#define MULTBINS 1
#define PARTTYPES 5

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "THnSparse.h"
#include "AliPIDResponse.h"
#include "AliAODpidUtil.h"

class AliAnalysisUtils;

class AliAnalysisTaskParticleEffWRZ :public AliAnalysisTaskSE{
 public:
 AliAnalysisTaskParticleEffWRZ() : AliAnalysisTaskSE(),centrality(0), fHistoList(0),  fHistEv(0), fpidResponse(0), fAODpidUtil(0)
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

  AliAnalysisTaskParticleEffWRZ(const Char_t *partName); // default constructor
  virtual ~AliAnalysisTaskParticleEffWRZ(); // default destructor
  virtual void UserCreateOutputObjects(); // user create output objects
  virtual void UserExec(Option_t *option); // user exec
  
 private:
  AliAnalysisTaskParticleEffWRZ(const AliAnalysisTaskParticleEffWRZ &); // copy constructor
  AliAnalysisTaskParticleEffWRZ &operator=(const AliAnalysisTaskParticleEffWRZ &); // operator=
  AliCentrality *centrality;
  TList *fHistoList; // histo list
  TH1F *fHistEv;
 
  TH1F *fHistQA[11];
  TH2F *fHistQA2D[3];
  TH2F *fHistQAPID[5][PARTTYPES][2];
  TH2F *fHistQAPIDFail[5][PARTTYPES][2];
  TH1F* fHistEvCuts[MULTBINS];
  TH2F *fHistQALambdas[2];
  TH2F *fOriginLambdas[5][2];

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

  
  TH2D *fPrim_DCAxy_Pt[MULTBINS*PARTTYPES][2];
  TH2D *fSecMat_DCAxy_Pt[MULTBINS*PARTTYPES][2];
  TH2D *fSecWeak_DCAxy_Pt[MULTBINS*PARTTYPES][2];

  TH2D *fPrim_DCAz_Pt[MULTBINS*PARTTYPES][2];
  TH2D *fSecMat_DCAz_Pt[MULTBINS*PARTTYPES][2];
  TH2D *fSecWeak_DCAz_Pt[MULTBINS*PARTTYPES][2];

  double fDCAtoPrimVtx;
  double get_mass_squared(AliAODTrack *track);
  AliPIDResponse *fpidResponse;
  AliAODpidUtil  *fAODpidUtil;
  ClassDef(AliAnalysisTaskParticleEffWRZ, 0);

};

#endif /* ALIANALYSISTASKPARTICLEEFFICIENCY */
