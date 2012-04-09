//********************************************************************
//********************************************************************
// Author : Zhongbao Yin
//********************************************************************
//********************************************************************
#ifndef AliAnalysisTaskFlowEPCascade_cxx
#define AliAnalysisTaskFlowEPCascade_cxx

#include "AliAnalysisTaskSE.h"
#include "AliEventplane.h"

class AliESDEvent;
//class AliFlowEventCuts;
class AliFlowTrackCuts;
class TH1F;
class TH2F;
class AliPIDResponse;
class TProfile;
class TProfile2D;
class TFile;

class AliAnalysisTaskFlowEPCascade : public AliAnalysisTaskSE {
 public:

  AliAnalysisTaskFlowEPCascade();
    
  AliAnalysisTaskFlowEPCascade(const char *name, Double_t centMin, 
			       Double_t centMax, 
			       Double_t xis[3][2], 
			       Double_t omegas[3][2] );
  virtual ~AliAnalysisTaskFlowEPCascade();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  //  void  SetFlowEventCuts(AliFlowEventCuts* cuts){fFlowEventCuts = cuts;}
  //  void  SetFlowTrackCuts(AliFlowTrackCuts* cuts){fFlowTrackCuts = cuts;}
  void  SetFlowDauTrackCuts(AliFlowTrackCuts* cuts){fCutsDau = cuts;}
  void  SetVertexCut(Double_t vtxCut = 10.){fVtxCut = vtxCut;}

 private:
  // AliFlowEventCuts*     fFlowEventCuts       ;
  //AliFlowTrackCuts*     fFlowTrackCuts       ;
  AliFlowTrackCuts *fCutsDau; // cuts for daughters
  AliPIDResponse * fPIDResponse;

  Double_t fXiBands[3][2]; //
  Double_t fOmegaBands[3][2]; //

  Double_t fMinCent, fMaxCent;
  Double_t fVtxCut;
  // Double_t fEtaCut;
  //Double_t fMinPt;
  //Int_t fNTPCcls;

  TFile * fOADB; //
  Int_t fRun; //! current run checked to load VZERO calibration
  Int_t fICent; //! centrality bin number

  TProfile *fMultV0; //! object containing VZERO calibration info
  Float_t fV0Cpol;          //! loaded by OADB       
  Float_t fV0Apol;          //! loaded by OADB                                
  Float_t fMeanQ[9][2][2];//! and recentering                  
  Float_t fWidthQ[9][2][2]; //! ...                              
  
  TList* fHistList; //!
  TH1I * fhEvent;  //!
  TH1F * fhEPangleVZero; //!
  TH1F * fhEPangleV0A; //!
  TH1F * fhEPangleV0C; //!
  TH1F * fhEPangleTPC; //!

  TH1F *fh1Chi2Xi; //!
  TH1F *fh1DCAXiDaughters; //!
  TH1F *fh1DCABachToPrimVertex; //!
  TH1F *fh1XiCosOfPointingAngle; //!
  TH1F *fh1XiRadius; //!
  
  TH1F *fh1MassLambda; //! mass of lambda as the cascade daughter
  TH1F *fh1V0Chi2; //! for V0 associated to a cascade
  TH1F *fh1V0CosOfPointingAngle; //!
  TH1F *fh1V0Radius; //!
  TH1F *fh1DcaV0DaughtersXi; //!
  TH1F *fh1DcaV0ToPrimVertex; //! DCA of V0 to prim. vtx
  TH1F *fh1DCAPosToPrimVertex; //! V0 positive daughter to prim. vertex
  TH1F *fh1DCANegToPrimVertex; //! V0 neg. daughter to prim. vertex
  
  TH1F *fh1MassXiMinus; //! effective mass under Xi- hyp.
  TH1F *fh1MassXiPlus; //!
  TH1F *fh1MassOmegaMinus; //! effective mass under Omega- hyp.
  TH1F *fh1MassOmegaPlus; //!

  TH1F *fh1MassXi;    //!
  TH1F *fh1MassOmega; //!

  TH1F *fh1XiPt; //! transverse momentum
  TH1F *fh1XiP; //! total momentum
  TH1F *fh1XiBachPt; //!  
  TH1F *fh1XiBachP; //!

  TH1F *fh1ChargeXi;//!
  TH1F *fh1V0toXiCosOfPointingAngle; //! cos of pointing angle btw the V0 mom and the Xi-V0 vtx line

  TH1F *fh1PhiXi; //!
  
  TH2F *fh2Armenteros; //! alpha vs ptArm for casc. candidate
  
  TH2F *fh2MassLambdaVsMassXiMinus; //! Xi- effective mass vs V0 eff. mass
  TH2F *fh2MassXiVsMassOmegaMinus; //!
  TH2F *fh2MassLambdaVsMassXiPlus; //!
  TH2F *fh2MassXiVsMassOmegaPlus; //!

  TH2F *fh2XiRadiusVsMassXiMinus; //! decay radius under Xi- hyp.
  TH2F *fh2XiRadiusVsMassXiPlus; //!
  TH2F *fh2XiRadiusVsMassOmegaMinus; //! decay radius under Omega- hyp.
  TH2F *fh2XiRadiusVsMassOmegaPlus; //!

  TH2F *fh2TPCdEdxOfCascDghters; //! dEdx with the cascade daughters
  
  TH2F *fh2MassVsPtXiMinus; //!
  TH2F *fh2MassVsPtXiPlus; //!
  TH2F *fh2MassVsPtXiAll; //!
  
  TH2F *fh2MassVsPtOmegaMinus; //!
  TH2F *fh2MassVsPtOmegaPlus; //!
  TH2F *fh2MassVsPtOmegaAll; //!

  TH1F * fhXiRapidity; //! 
  TH1F * fhOmegaRapidity; //!
  
  TProfile *fProfXiV2PtV0A[3]; //! three mass bands V0A
  TProfile *fProfOmegaV2PtV0A[3];//!
  TProfile *fProfXiSinePtV0A[3];//!
  TProfile *fProfOmegaSinePtV0A[3];//!
  
  TProfile *fProfXiV2PtV0C[3]; //! three mass bands V0C
  TProfile *fProfOmegaV2PtV0C[3];//!
  TProfile *fProfXiSinePtV0C[3];//!
  TProfile *fProfOmegaSinePtV0C[3];//!

  TProfile *fProfXiV2Pt[3]; //! three mass bands V0A&V0C
  TProfile *fProfOmegaV2Pt[3]; //!
  TProfile *fProfXiSinePt[3]; //!
  TProfile *fProfOmegaSinePt[3]; //!

  TProfile2D *fProf2dXiV2PtV0A[3]; //! three mass bands V0A                 
  TProfile2D *fProf2dOmegaV2PtV0A[3];//! 
  TProfile2D *fProf2dXiV2PtV0C[3]; //! three mass bands V0C
  TProfile2D *fProf2dOmegaV2PtV0C[3];//!
  TProfile2D *fProf2dXiV2Pt[3]; //! three mass bands V0A&V0C
  TProfile2D *fProf2dOmegaV2Pt[3]; //!

  TProfile * fProfResolution;  //!

  TH1F *fh1DistToVtxZAfter; //! After propagation to vertex 
  TH1F *fh1DistToVtxXYAfter; //!                                               
  TH2F *fh2DistToVtxZBeforeVsAfter; //!
  TH2F *fh2DistToVtxXYBeforeVsAfter; //!
  TH2F *fh2PxBeforeVsAfter; //!
  TH2F *fh2PyBeforeVsAfter; //!
  TH2F *fh2PhiPosBeforeVsAfter; //! 
  TH2F *fh2PhiNegBeforeVsAfter; //!
  
  //  void ReadFromESDv0(AliESDEvent *fESD);
  void ReadFromAODv0(AliAODEvent *fAOD);

  //Progate to the primary vertex
  void Propagate(Double_t vv[3], Double_t x[3], Double_t p[3], 
		 Double_t bz, Short_t sign);

  void OpenInfoCalbration(Int_t run );//9 bins: 0-5,5-10,10-20

  AliAnalysisTaskFlowEPCascade(const AliAnalysisTaskFlowEPCascade&); // not implemented
  AliAnalysisTaskFlowEPCascade& operator=(const AliAnalysisTaskFlowEPCascade&); // not implemented
  
  ClassDef(AliAnalysisTaskFlowEPCascade, 1); // example of analysis
};
#endif
