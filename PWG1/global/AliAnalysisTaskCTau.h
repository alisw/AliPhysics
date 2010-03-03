#ifndef AliAnalysisTaskCTau_h
#define AliAnalysisTaskCTau_h

#include "AliAnalysisTaskSE.h"

class TH1F;
class AliESDEvent;
class TTree;

//
//  This is a little task for checking the c*tau of the strange particles 
//

class AliAnalysisTaskCTau : public AliAnalysisTaskSE {

public:

  AliAnalysisTaskCTau(const char *name = "AliAnalysisTaskCTau");
  virtual ~AliAnalysisTaskCTau() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);  


private: 

  AliAnalysisTaskCTau(const AliAnalysisTaskCTau&);           //not implemented
  AliAnalysisTaskCTau& operator=(const AliAnalysisTaskCTau&);//not implemented 

  TList       *fOutput;   //! The list of histograms

  AliESDEvent    *fESD ; //! ESD event

  TH1F* fK0s;         //! cTau for K0s
  TH1F* fK0sMC;       //! cTau for the K0s in Monte Carlo 

  TH1F* fLambdas;     //! cTau for Lambdas
  TH1F* fLambdasMC;   //! cTau for Lambdas in Monte Carlo

  TH1F* fLambdaBars;  //! cTau for anti-Lambdas
  TH1F* fLambdaBarsMC;//! cTau for anti-Lambdas in Monte Carlo

  TH1F* fXis;         //! cTau for Xis
  TH1F* fXisMC;       //! cTau for Xis in Monte Carlo

  TH1F* fMass;        //! Effective mass for reconstructed V0s
  TH1F* fMassMC;      //! Effective mass for associated V0s

  ClassDef(AliAnalysisTaskCTau,0);
};

#endif
