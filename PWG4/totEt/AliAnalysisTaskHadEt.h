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

#include "AliAnalysisTaskSE.h"


class AliAnalysisTaskHadEt : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskHadEt(const char *name = "AliAnalysisTaskHadEt");
    virtual ~AliAnalysisTaskHadEt();

    //  virtual void   ConnectInputData(Option_t *);
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);

    AliESDtrackCuts* GetTPCITSTrackCuts(){return (AliESDtrackCuts*) fOutputList->FindObject("fEsdTrackCuts");}
    AliESDtrackCuts* GetTPCOnlyTrackCuts(){return (AliESDtrackCuts*) fOutputList->FindObject("fEsdTrackCutsTPCOnly");}
    AliESDtrackCuts* GetITSTrackCuts(){return (AliESDtrackCuts*) fOutputList->FindObject("fEsdTrackCutsITS");}

private:

  //Declare it private to avoid compilation warning
    AliAnalysisTaskHadEt & operator = (const AliAnalysisTaskHadEt & g) ;//cpy assignment
    AliAnalysisTaskHadEt(const AliAnalysisTaskHadEt & g) ; // cpy ctor

    TList *fOutputList; //output list

    AliAnalysisHadEtReconstructed *fRecAnalysis; // Rec
    AliAnalysisHadEtMonteCarlo *fMCAnalysis; // MC

    TH2F *fHistEtRecvsEtMC; // Rec vs MC histo 
    
    AliESDtrackCuts* fEsdtrackCutsITSTPC; // track cuts ITS&TPC
    AliESDtrackCuts* fEsdtrackCutsTPC; // track cuts TPC
    AliESDtrackCuts* fEsdtrackCutsITS; // track cuts ITS
    
    ClassDef(AliAnalysisTaskHadEt, 1); // example of analysis
};

#endif
