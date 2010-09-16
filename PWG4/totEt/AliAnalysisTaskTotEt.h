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
class TH2F;
class TList;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskTotEt : public AliAnalysisTaskSE {
  
public:
  AliAnalysisTaskTotEt(const char *name = "AliAnalysisTaskTotEt");
  virtual ~AliAnalysisTaskTotEt() {}
private:
  //Declare it private to avoid compilation warning
  AliAnalysisTaskTotEt & operator = (const AliAnalysisTaskTotEt & g) ;//cpy assignment
  AliAnalysisTaskTotEt(const AliAnalysisTaskTotEt & g) ; // cpy ctor
  
public:
  
  //  virtual void   ConnectInputData(Option_t *);
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  
private:
  
  TList *fOutputList; //output list
  
  AliAnalysisEtReconstructed *fRecAnalysis; // Rec 
  AliAnalysisEtMonteCarlo *fMCAnalysis; // MC
  
  TH2F *fHistEtRecvsEtMC; // Rec vs MC histo
  
  ClassDef(AliAnalysisTaskTotEt, 1); // example of analysis
};

#endif
