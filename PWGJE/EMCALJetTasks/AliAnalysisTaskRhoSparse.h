#ifndef ALIANALYSISTASKRHOSPARSE_H
#define ALIANALYSISTASKRHOSPARSE_H

// $Id: AliAnalysisTaskRho.h 58408 2012-09-03 07:00:58Z loizides $

#include "AliAnalysisTaskRhoBase.h"

class AliAnalysisTaskRhoSparse : public AliAnalysisTaskRhoBase {

 public:
  AliAnalysisTaskRhoSparse();
  AliAnalysisTaskRhoSparse(const char *name, Bool_t histo=kFALSE);
  virtual ~AliAnalysisTaskRhoSparse() {}

  void             UserCreateOutputObjects();
  void             SetExcludeLeadJets(UInt_t n)    { fNExclLeadJets = n    ; }
  void             SetRhoCMS(Bool_t cms)           { fRhoCMS = cms ; }


 protected:
  Bool_t           Run();
  TH2F             *fHistOccCorrvsCent;             //!occupancy correction vs. centrality

  UInt_t           fNExclLeadJets;                 // number of leading jets to be excluded from the median calculation

  Bool_t           fRhoCMS;


  AliAnalysisTaskRhoSparse(const AliAnalysisTaskRhoSparse&);             // not implemented
  AliAnalysisTaskRhoSparse& operator=(const AliAnalysisTaskRhoSparse&);  // not implemented
  
  ClassDef(AliAnalysisTaskRhoSparse, 1); // Rho task
};
#endif
