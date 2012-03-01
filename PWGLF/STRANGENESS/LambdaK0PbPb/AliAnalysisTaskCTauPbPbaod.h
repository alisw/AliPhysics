#ifndef AliAnalysisTaskCTauPbPbaod_h
#define AliAnalysisTaskCTauPbPbaod_h

#include "AliAnalysisTaskSE.h"

class TH1F;
class TH2F;
class TH3F;
class TList;

class AliAODEvent;
class AliAODv0;

//
//  This is a little task for checking the c*tau of the strange particles 
//

class AliAnalysisTaskCTauPbPbaod : public AliAnalysisTaskSE {

public:

  AliAnalysisTaskCTauPbPbaod(const char *name = "AliAnalysisTaskCTauPbPbaod");
  virtual ~AliAnalysisTaskCTauPbPbaod() {}

  void SetCentrality(Double_t min, Double_t max) {fCMin=min;fCMax=max;} 
  void SetMC(Bool_t isMC=kTRUE) {fIsMC=isMC;} 
  void SetCosPA(Double_t cospa) {fCPA=cospa;} 
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);  


private: 

  AliAnalysisTaskCTauPbPbaod(const AliAnalysisTaskCTauPbPbaod&);           //not implemented
  AliAnalysisTaskCTauPbPbaod& operator=(const AliAnalysisTaskCTauPbPbaod&);//not implemented 

  Bool_t AcceptV0(const AliAODv0 *v0, const AliAODEvent *aod); 

  Bool_t fIsMC;
  Double_t fCMin;       // Min centrality
  Double_t fCMax;       // Max centrality
  Double_t fCPA;        // cos(PA) threshold

  TList       *fOutput; //! The list of histograms

  TH1F *fMult;       //! Track multiplicity
  TH1F *fCosPA;      //! Track multiplicity
  TH2F* fdEdx;       //! dEdx
  TH2F* fdEdxPid;    //! dEdx with PID


  TH2F* fK0sM;       //! Mass for K0s
  TH2F* fK0sSi;      //! Side-band subtracted LvsP  for K0s 
  TH2F* fK0sMC;      //! LvsP for the K0s from the Monte Carlo stack 
  TH2F* fK0sAs;      //! LvsP for the K0s associated with the Monte Carlo 


  TH2F* fLambdaM;    //! Mass for Lambdas
  TH2F* fLambdaSi;   //! Side-band subtrated LvsP for Lambda
  TH2F* fLambdaMC;   //! LvsP for Lambdas from the Monte Carlo stack
  TH2F* fLambdaAs;   //! LvsP for Lambdas associated with the Monte Carlo

  TH1D* fLambdaEff;  //! Efficiency for Lambda  
  TH1D* fLambdaPt;   //! Pt spectrum for Lambda

  TH2F* fLambdaBarM;  //! Mass for anti-Lambdas
  TH2F* fLambdaBarSi; //! Side-band subtrated LvsP for anti-Lambda
  TH2F* fLambdaBarMC; //! LvsP for anti-Lambdas from the Monte Carlo stack
  TH2F* fLambdaBarAs; //! LvsP for anti-Lambdas associated with the Monte Carlo

  TH1D* fLambdaBarEff;  //! Efficiency for anti-Lambda  
  TH1D* fLambdaBarPt;   //! Pt spectrum for anti-Lambda

  TH3F* fLambdaFromXi;//! LvsPvsPxi for Lambdas from Xis associated with MC 
  TH2F* fXiM;         //! Mass for Xis
  TH1F* fXiSiP;       //! Side-band subtracted Pt for reconstructed Xi

  TH3F* fLambdaBarFromXiBar;//! LvsPvsPxi for anti-Lambdas from anti-Xis associated with MC 
  TH2F* fXiBarM;         //! Mass for anti-Xis
  TH1F* fXiBarSiP;       //! Side-band subtracted Pt for reconstructed anti-Xi

  ClassDef(AliAnalysisTaskCTauPbPbaod,4);
};

#endif
