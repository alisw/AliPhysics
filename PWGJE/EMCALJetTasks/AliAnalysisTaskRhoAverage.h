#ifndef ALIANALYSISTASKRHOAVERAGE_H
#define ALIANALYSISTASKRHOAVERAGE_H

// $Id$

#include "AliAnalysisTaskRhoBase.h"

class AliAnalysisTaskRhoAverage : public AliAnalysisTaskRhoBase {

 public:
  AliAnalysisTaskRhoAverage();
  AliAnalysisTaskRhoAverage(const char *name, Bool_t histo=kFALSE);
  virtual ~AliAnalysisTaskRhoAverage() {}

  void             SetRhoType(Int_t t)             { fRhoType       = t    ; }
  void             SetExcludeLeadPart(UInt_t n)    { fNExclLeadPart = n    ; }
  
 protected:
  Bool_t           Run();

  Int_t            fRhoType       ;// rho type: 0 = charged+neutral, 1 = charged, 2 = neutral
  UInt_t           fNExclLeadPart ;// number of leading particles to be excluded from the median calculation

  AliAnalysisTaskRhoAverage(const AliAnalysisTaskRhoAverage&);             // not implemented
  AliAnalysisTaskRhoAverage& operator=(const AliAnalysisTaskRhoAverage&);  // not implemented
  
  ClassDef(AliAnalysisTaskRhoAverage, 3); // Rho task
};
#endif
