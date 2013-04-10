#ifndef ALIANALYSISTASKRHO_H
#define ALIANALYSISTASKRHO_H

// $Id$

#include "AliAnalysisTaskRhoBase.h"

class AliAnalysisTaskRho : public AliAnalysisTaskRhoBase {

 public:
  AliAnalysisTaskRho();
  AliAnalysisTaskRho(const char *name, Bool_t histo=kFALSE);
  virtual ~AliAnalysisTaskRho() {}

  void             SetExcludeLeadJets(UInt_t n)    { fNExclLeadJets = n    ; }

 protected:
  Bool_t           Run();

  UInt_t           fNExclLeadJets;                 // number of leading jets to be excluded from the median calculation

  AliAnalysisTaskRho(const AliAnalysisTaskRho&);             // not implemented
  AliAnalysisTaskRho& operator=(const AliAnalysisTaskRho&);  // not implemented
  
  ClassDef(AliAnalysisTaskRho, 10); // Rho task
};
#endif
