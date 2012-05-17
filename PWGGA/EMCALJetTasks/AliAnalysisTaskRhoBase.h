#ifndef ALIANALYSISTASKRHOBASE_cxx
#define ALIANALYSISTASKRHOBASE_cxx

// $Id$

class TString;
class TF1;

#include <TParameter.h>

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskRhoBase : public AliAnalysisTaskSE {

 public:
  AliAnalysisTaskRhoBase();
  AliAnalysisTaskRhoBase(const char *name);
  virtual ~AliAnalysisTaskRhoBase() {}
  
  virtual void          UserCreateOutputObjects();
  virtual void          UserExec(Option_t*);
  virtual void          Terminate(Option_t*);

  void                  SetRhoFunction(TF1* rf)                               { fRhoFunction   = rf   ; }
  void                  SetRhoName(const char *name)                          { fRhoName       = name ; }
  
 protected:
  virtual Double_t       GetRhoFactor(Double_t cent);

  TString                fRhoName;                       // name of rho
  TF1                   *fRhoFunction;                   // pre-computed rho as a function of centrality
  Double_t               fCent;                          //!event centrality
  TParameter<Double_t>  *fRho;                           //!per event calculated rho

  AliAnalysisTaskRhoBase(const AliAnalysisTaskRhoBase&);             // not implemented
  AliAnalysisTaskRhoBase& operator=(const AliAnalysisTaskRhoBase&);  // not implemented
  
  ClassDef(AliAnalysisTaskRhoBase, 1); // Rho base task
};
#endif
