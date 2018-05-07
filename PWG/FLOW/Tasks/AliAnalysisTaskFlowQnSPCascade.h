/////////////////////////////////////////////////////
// AliAnalysisTaskFlowCascade:
// Analysis task to select Xi and Omega candidates for flow analysis.
// Author: 
//////////////////////////////////////////////////////

#ifndef AliAnalysisTaskFlowQnSPCascade_cxx
#define AliAnalysisTaskFlowQnSPCascade_cxx

#include "AliAnalysisTaskSE.h"
#include "AliEventplane.h"
#include "TFile.h"
#include "TTree.h"
#include "AliQnCorrectionsFillEventTask.h"
#include "AliEventCuts.h"

class AliEventCuts;
class AliAnalysis;
class AliQnCorrectionsManager;
class AliQnCorrectionsCutsSet;
class AliQnCorrectionsHistos;
class TList;
class TProfile;
class TGraphErrors;
class AliESDEvent;
class AliFlowTrackCuts;
class TH1F;
class TH2F;
class AliPIDResponse;
class AliQnCorrectionsManager;
class TProfile2D;
class TFile;
class AliAnalysisTaskFlowQnSPCascade : public AliAnalysisTaskSE {
 public:

  AliAnalysisTaskFlowQnSPCascade();
  //AliEventCuts *fEventCuts; /// Event cuts
  AliEventCuts *fEventCuts;//! standard event cuts
  AliAnalysisTaskFlowQnSPCascade(const char *name, Double_t centMin, 
			       Double_t centMax , bool WeightCorrection);
  virtual ~AliAnalysisTaskFlowQnSPCascade();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);


  void  SetOnline(){fOnline = kTRUE;}
  void  SetSepAnalysis(Bool_t value){fSepAnalysis = value;}
  void  SetHarmonicOrder(Int_t value){fHarmonicOrder =value;}
//-----------------------Xi--------------------------------------------------
  void  SetXiEtaMin(Double_t value){fXiPseMin = value;} 
  void  SetXiEtaMax(Double_t value){fXiPseMax = value;}
  void  SetV0RadiusXiMin(Double_t value){fV0RadiusXiMin = value;}
  void  SetV0RadiusXiMax(Double_t value){fV0RadiusXiMax = value;}
  void  SetXiRadiusMin(Double_t value){fXiRadiusMin = value;}
  void  SetXiRadiusMax(Double_t value){fXiRadiusMax = value;}
  void  SetDCAXiDaughtersMax(Double_t value){fdcaXiDaughtersMax = value;}
  void  SetXiCosOfPointingAngleMin(Double_t value){fXiCosOfPointingAngleMin = value;}
  void  SetDCAV0ToPrimaryVtxXiMin(Double_t value){fdcaV0ToPrimaryVtxXiMin = value;}
  void  SetDCABachToPrimaryVtxXiMin(Double_t value){fdcaBachToPrimaryVtxXiMin = value;}
  void  SetLambdaMassWindow(Double_t value){fLambdaMassWind = value;}
  void  SetDCAV0DaughtersXi(Double_t value){fdcaV0DaughtersXi = value;}
  void  SetV0CosOfPointingAngleXiMin(Double_t value){fV0CosOfPointingAngleXiMin = value;}
  void  SetDCAPosToPrimaryVtxXiMin(Double_t value){fdcaPosToPrimaryVtxXiMin = value;}
  void  SetDCANegToPrimaryVtxXiMin(Double_t value){fdcaNegToPrimaryVtxXiMin = value;}
//------------------------V0---------------------------------------------------
//  void  SetV0Eta(Double_t value){tV0Eta = value;}
 // void  SetV0PtMin(Double_t value){tV0Pt = value;} 
  void  SetV0Rapidity(Double_t value){tDecayRapidity = value;}
  void  SetV0DecayRadius(Double_t value){tDecayRad = value;}
  void  SetV0DCADaughtersMax(Double_t value){tV0DCAdaughtersMax = value;}
  void  SetV0CosinePointingAngleMin(Double_t value){tV0CosinePointingAngleMin = value;}
  void  SetV0DCAToPrimVertexMin(Double_t value){tV0DCAToPrimVertexMin = value;}
  void  SetV0LifeTimeMax(Double_t value){tV0LifeTimeMax = value;}
  void  SetV0DecayLength(Double_t value){tDecayLengthV0 = value;}
//-----------------------track-------------------------------------------------  
  void  SetRPTrackFromTPC(Bool_t value){fRPFromTPC = value;}
  void  SetPrimaryTrackEta(Double_t value){fPrimaryTrackEta = value;}
  void  SetTrackEta(Double_t value){fTrackEta = value;}
  void  SetTrackPtMin(Double_t value){fTrackPtMin = value;}
  void  SetTPCNcls(Double_t ncls = 70){fTPCNcls = ncls;}

//-----------------------PID-------------------------------------------------        
  void  SetXiPIDSigma(Double_t value){fXiPIDsigma = value;}
  void  SetV0PIDSigma(Double_t value){fV0PIDsigma = value;}
 
//======================================================================================== 
  void SetUseHybridGlobalTrack(){fUseHybridGlobalTrack = kTRUE;}

  void  SetAnaObjectMultiStrange(){fIsMultiStrange = kTRUE;fIsStrange = kFALSE;}
  void  SetAnaObjectStrange(){fIsStrange = kTRUE;fIsMultiStrange = kFALSE;}
  void  SetExpectedCorrectionPass(const char *pass) 
        { fExpectedCorrectionPass = pass; }
  void  SetAlternativeCorrectionPass(const char *pass) 
        { fAlternativeCorrectionPass = pass; }

 private:
  TString fExpectedCorrectionPass;
  TString fAlternativeCorrectionPass;
  Double_t fMinCent, fMaxCent;
  Bool_t fWeightCorrection;
  TH1F *hOmegaWeight;//!
  TH1F *hXiWeight;//!
//--------------------------param----------------------------------------
  Double_t weightXi[20];
  Double_t weightOmega[20];
  Double_t fPrimaryTrackEta = 0;
  Double_t XiPse = 0;
  Double_t V0RadiusXi = 0;
  Double_t XiRadius = 0;
  Double_t dcaXiDaughters = 0;
  Double_t XiCosOfPointingAngle = 0;
  Double_t dcaV0ToPrimaryVtxXi = 0;
  Double_t dcaBachToPrimaryVtxXi = 0;
  Double_t invMassLambdaAsCascDghter = 0;
  Double_t dcaV0DaughtersXi = 0;
  Double_t V0CosOfPointingAngleXi = 0;
  Double_t dcaPosToPrimaryVtxXi = 0;
  Double_t dcaNegToPrimaryVtxXi = 0;
  Double_t fV0DCAdaughters = 0;
  Double_t fV0CosinePointingAngle = 0;
  Double_t fDecayRad = 0;
  Double_t fV0DCAToPrimVertex = 0;
  Double_t V0LifeTime = 0;               
//-----------------------Xi--------------------------------------------------
  Int_t fHarmonicOrder = 0;
  Double_t fXiPseMin = 0;
  Double_t fXiPseMax = 0;
  Double_t fV0RadiusXiMin = 0;
  Double_t fV0RadiusXiMax = 0;
  Double_t fXiRadiusMin = 0;
  Double_t fXiRadiusMax = 0;
  Double_t fdcaXiDaughtersMax = 0;
  Double_t fXiCosOfPointingAngleMin = 0;
  Double_t fdcaV0ToPrimaryVtxXiMin = 0;
  Double_t fdcaBachToPrimaryVtxXiMin = 0;
  Double_t fLambdaMassWind = 0;
  Double_t fdcaV0DaughtersXi = 0;
  Double_t fV0CosOfPointingAngleXiMin = 0;
  Double_t fdcaPosToPrimaryVtxXiMin = 0;
  Double_t fdcaNegToPrimaryVtxXiMin = 0;
//-----------------V0---------------------------------------------------
  Double_t tV0Eta =0;
  Double_t tV0Pt = 0;
  Double_t tDecayRapidity = 0;
  Double_t tDecayRad = 0;
  Double_t tV0DCAdaughtersMax = 0;
  Double_t tV0CosinePointingAngleMin = 0;
  Double_t tV0DCAToPrimVertexMin = 0;
  Double_t tV0LifeTimeMax = 0;
  Double_t tDecayLengthV0 = 0;

//----------------track-------------------------------------------------  
  Double_t fTPCNcls; // number of TPC clusters   
  Double_t fTrackEta = 0;
  Double_t fTrackPtMin = 0;
  Bool_t   fRPFromTPC;
//-------Calibration----------------
 Bool_t   fRemChV0A;
 TH1D* fMultV0;
 Int_t fNHarm = 0;
 Int_t fRun =0;
/*
   Double_t   fQxnmV0A =0;
    Double_t  fQynmV0A =0;
    Double_t  fQxnsV0A =0;
    Double_t  fQynsV0A =0;
   Double_t   fQxnmV0C =0;
    Double_t  fQynmV0C =0;
    Double_t  fQxnsV0C =0;
    Double_t  fQynsV0C =0;
*/
  TH1D*        fQxnmV0A;            // <Qx2> V0A
  TH1D*        fQynmV0A;            // <Qy2> V0A
 TH1D*        fQxnsV0A;            // sigma Qx2 V0A
 TH1D*        fQynsV0A;            // sigma Qy2 V0A
 TH1D*        fQxnmV0C;            // <Qx2> V0C
 TH1D*        fQynmV0C;            // <Qy2> V0C
 TH1D*        fQxnsV0C;            // sigma Qx2 V0C
 TH1D*        fQynsV0C;            // sigma Qy2 V0C
  

//----------------PID-------------------------------------------------        
  Double_t fXiPIDsigma = 0;
  Double_t fV0PIDsigma = 0;

//================================================================================
  
  AliPIDResponse * fPIDResponse;
  AliQnCorrectionsManager *fFlowQnVectorMgr;
  Bool_t fIsMultiStrange; // 
  Bool_t fIsStrange;     //
  Bool_t fOnline;
  Bool_t fSepAnalysis;
  Bool_t fUseHybridGlobalTrack;
  TList* fHistList; //!


//------------------------------------------Event----------------------------------------------------------//  
  void QAEventInput();
  void QAEventPlane();
  void QATrack();  
//-----------------------------------multi-strange bayon particle------------------------------------------//
  void ReadFromAODCascade(AliAODEvent *fAOD, Double_t psiV0A, Double_t psiV0C, 
                          Double_t psiTPC ,Double_t Qn_Ax ,Double_t Qn_Ay);
  void QACascadeCandidates();
  Bool_t AcceptCascadeCandidate();
  void QACascadeFlow();
  void FlowAnaCascade();
//-------------------------------------strange bayon particle-----------------------------------------------//
  void ReadFromAODv0(AliAODEvent *fAOD, Double_t psiV0A, Double_t psiV0C, Double_t psiTPC);
  void QAV0Candidates();
  void LoadTrack(AliESDtrack *myTrack, Double_t aodChi2NDF=0);
  Bool_t AcceptV0Candidate();
  void FlowAnaV0();
  void SepFlowAnaV0();
  void SepFlowAnaCascade();
//--------------------------------------track selection-----------------------------------------------------//
//  Bool_t IsGoodPrimaryTrack(const AliAODTrack* aodtrack);
  Bool_t IsSelected(AliAODTrack *aodTrk);


//---------------------------------------Propagate in the electromagnetic field-----------------------------//
  void Propagate(Double_t vv[3],Double_t x[3], Double_t p[3], Double_t bz,Short_t sign);
  void OpenInfoCalbration(Int_t run);

//----------------------------------------------------------------------------------------------------------//
  AliAnalysisTaskFlowQnSPCascade(const AliAnalysisTaskFlowQnSPCascade&); // not implemented
  AliAnalysisTaskFlowQnSPCascade& operator=(const AliAnalysisTaskFlowQnSPCascade&); // not implemented
  
  ClassDef(AliAnalysisTaskFlowQnSPCascade, 1); // example of analysis
};
#endif
