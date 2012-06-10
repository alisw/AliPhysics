#ifndef ALIANALYSISTASKRHOBASE_H
#define ALIANALYSISTASKRHOBASE_H

// $Id$

class TString;
class TF1;
class AliRhoParameter;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskRhoBase : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskRhoBase();
  AliAnalysisTaskRhoBase(const char *name);
  virtual ~AliAnalysisTaskRhoBase() {}
  
  void                   UserCreateOutputObjects();
  void                   UserExec(Option_t*);

  const char            *GetRhoName() const                                    { return fRhoName       ; }
  void                   SetRhoFunction(TF1* rf)                               { fRhoFunction   = rf   ; }
  void                   SetRhoName(const char *name)                          { fRhoName       = name ; }

 protected:
  virtual void           DetermineCent();
  virtual void           ExecOnce();
  TString                GetBeamType();
  virtual Double_t       GetRhoFactor(Double_t cent);

  TString                fRhoName;                       // name of rho
  TF1                   *fRhoFunction;                   // pre-computed rho as a function of centrality
  Double_t               fCent;                          //!event centrality
  AliRhoParameter       *fRho;                           //!per event calculated rho
  Bool_t                 fDoCent;                        //!==1 then do centrality
  Bool_t                 fIsInit;                        //!==1 then do init

  AliAnalysisTaskRhoBase(const AliAnalysisTaskRhoBase&);             // not implemented
  AliAnalysisTaskRhoBase& operator=(const AliAnalysisTaskRhoBase&);  // not implemented
  
  ClassDef(AliAnalysisTaskRhoBase, 3); // Rho base task
};
#endif
