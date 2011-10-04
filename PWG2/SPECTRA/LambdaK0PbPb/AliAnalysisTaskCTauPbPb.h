#ifndef AliAnalysisTaskCTauPbPb_h
#define AliAnalysisTaskCTauPbPb_h

#include "AliAnalysisTaskSE.h"

class TH1F;
class TH2F;
class TH3F;
class TList;

//
//  This is a little task for checking the c*tau of the strange particles 
//

class AliAnalysisTaskCTauPbPb : public AliAnalysisTaskSE {

public:

  AliAnalysisTaskCTauPbPb(const char *name = "AliAnalysisTaskCTauPbPb");
  virtual ~AliAnalysisTaskCTauPbPb() {}

  void SetCentrality(Double_t min, Double_t max) {fCMin=min;fCMax=max;} 
  void SetMC(Bool_t isMC=kTRUE) {fIsMC=isMC;} 
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);  


private: 

  AliAnalysisTaskCTauPbPb(const AliAnalysisTaskCTauPbPb&);           //not implemented
  AliAnalysisTaskCTauPbPb& operator=(const AliAnalysisTaskCTauPbPb&);//not implemented 

  Bool_t fIsMC;
  Double_t fCMin;       // Min centrality
  Double_t fCMax;       // Max centrality

  TList       *fOutput; //! The list of histograms

  TH1F *fMult;       //! Track multiplicity
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

  TH1F* fCPA;   //! cos(PA) side-band subtructed
  TH1F* fDCA;   //! DCA daughters side-band subtructed

  TH1D* fLambdaEff;  //! Efficiency for Lambda  
  TH1D* fLambdaPt;   //! Pt spectrum for Lambda

  TH3F* fLambdaFromXi;//! LvsPvsPxi for Lambdas from Xis associated with MC 
  TH2F* fXiM;         //! Mass for Xis
  TH1F* fXiSiP;       //! Side-band subtracted Pt for reconstructed Xi

  ClassDef(AliAnalysisTaskCTauPbPb,1);
};

#endif
