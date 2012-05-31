#ifndef ALIANALYSISTASKPHICORR_cxx
#define ALIANALYSISTASKPHICORR_cxx
/*
 Simple use case for mixed event analysis
 based on ESD or AOD
 Delta_phi correlation analysis is performed on charged tracks 
 for same and mixed events
 Author: andreas.morsch@cern.ch 
*/


class TH1F;
class TList;
class AliESDEvent;

#include "AliAnalysisTaskME.h"
#include "AliMixedEvent.h"

class AliAnalysisTaskPhiCorr : public AliAnalysisTaskME {
 public:
  AliAnalysisTaskPhiCorr(const char *name = "AliAnalysisTaskPt");
  virtual ~AliAnalysisTaskPhiCorr() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
  TList       *fHists;        // List with histos
  TH1F        *fHistDphiCO;   // Pt spectrum
  TH1F        *fHistDphiUC;   // Pt spectrum
  AliMixedEvent fMixedEvent;  // Mixed event
  
  AliAnalysisTaskPhiCorr(const AliAnalysisTaskPhiCorr&); // not implemented
  AliAnalysisTaskPhiCorr& operator=(const AliAnalysisTaskPhiCorr&); // not implemented
  
  ClassDef(AliAnalysisTaskPhiCorr, 1); // example of analysis
};

#endif
