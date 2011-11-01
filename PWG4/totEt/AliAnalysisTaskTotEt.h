#ifndef ALIANALYSISTASKTOTET_H 
#define ALIANALYSISTASKTOTET_H 
//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Task for analysis
//  - reconstruction and MC output
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________

class AliAnalysisEtReconstructed;
class AliAnalysisEtMonteCarlo;
class AliESDtrackCuts;
class TH2F;
class TList;

#include "AliAnalysisTaskTransverseEnergy.h"

class AliAnalysisTaskTotEt : public AliAnalysisTaskTransverseEnergy {
  
public:
  AliAnalysisTaskTotEt(const char *name = "AliAnalysisTaskTotEt", Bool_t isMc = false);
  virtual ~AliAnalysisTaskTotEt();
  
  //  virtual void   ConnectInputData(Option_t *);
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
    
private:
  
  //Declare private to avoid compilation warning
  AliAnalysisTaskTotEt & operator = (const AliAnalysisTaskTotEt & g) ;//copy assignment
  AliAnalysisTaskTotEt(const AliAnalysisTaskTotEt & g) ; // copy ctor

  AliAnalysisEtReconstructed *fRecAnalysis; // Rec 
  AliAnalysisEtMonteCarlo *fMCAnalysis; // MC
  
  //THnSparseD *fSparseHistRecVsMc; // Hist Rec vs Mc
  //Double_t *fSparseRecVsMc; // Rec vs Mc
  
  ClassDef(AliAnalysisTaskTotEt, 2) 
};

#endif
