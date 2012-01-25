//_________________________________________________________________________
//  Utility Class for transverse energy studies; charged hadrons
//  Task for analysis
//  - reconstruction and MC output
//
//Created by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
//_________________________________________________________________________
#ifndef ALIANALYSISTASKHADET_H 
#define ALIANALYSISTASKHADET_H 

class AliAnalysisHadEtReconstructed;
class AliAnalysisHadEtMonteCarlo;
class AliESDtrackCuts;
class TH2F;
class TList;

#include "AliAnalysisTaskTransverseEnergy.h"
class AliPWG0Helper;

class AliAnalysisTaskHadEt : public AliAnalysisTaskTransverseEnergy {
public:
  AliAnalysisTaskHadEt(const char *name = "AliAnalysisTaskHadEt", Bool_t isMc = false, TString recoConfigFile = "ConfigHadEtReconstructed.C", TString mcConfigFile = "ConfigHadEtMonteCarlo.C");
    virtual ~AliAnalysisTaskHadEt();

    //  virtual void   ConnectInputData(Option_t *);
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    void IsSim(){fIsSim = kTRUE;}

private:

  //Declare it private to avoid compilation warning
    AliAnalysisTaskHadEt & operator = (const AliAnalysisTaskHadEt & g) ;//cpy assignment
    AliAnalysisTaskHadEt(const AliAnalysisTaskHadEt & g) ; // cpy ctor

    AliAnalysisHadEtReconstructed *fRecAnalysis; // Rec
    AliAnalysisHadEtMonteCarlo *fMCAnalysis; // MC
    Bool_t fIsSim;//Boolean to keep track of whether or not this is running on simulations

    ClassDef(AliAnalysisTaskHadEt, 2); // example of analysis
};

#endif
