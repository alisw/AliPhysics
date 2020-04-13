
// For: Net Lambda fluctuation analysis via traditional method
// By: Ejiro Umaka Apr 2018
//update Mar 2019

#ifndef AliAnalysisTaskNetLambdaTrad_h
#define AliAnalysisTaskNetLambdaTrad_h


#include "AliAnalysisTaskSE.h"
class TList;
class AliESDEvent;
class AliESDtrack;
class AliAnalysisUtils;
class AliPIDResponse;
class TTree;
class TH1;
class TH2;
class TH3;
class TH3F;
#include "AliEventCuts.h"

class AliAnalysisTaskNetLambdaTrad : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskNetLambdaTrad(const char* name="AliAnalysisTaskNetLambdaTrad");
    virtual ~AliAnalysisTaskNetLambdaTrad ();
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    void SetEventSelection(UInt_t val) {fEvSel = val;}
    
protected:
    AliAnalysisTaskNetLambdaTrad(const  AliAnalysisTaskNetLambdaTrad &task);
    AliAnalysisTaskNetLambdaTrad& operator=(const  AliAnalysisTaskNetLambdaTrad &task);
    
    AliESDEvent* fESD;
    AliPIDResponse* fPIDResponse;
    AliEventCuts fEventCuts;
    TList* fListHist;
    
    TH1D*  fHistEventCounter;
    TH1D*  fHistCentrality;

    
    TH3F*  f3fHistCentVsInvMassLambda1point0;
    TH3F*  f3fHistCentVsInvMassAntiLambda1point0;
    
    TH3F*  f3fHistCentVsInvMassLambda1point0nsigtight;
    TH3F*  f3fHistCentVsInvMassAntiLambda1point0nsigtight;
    
    TH3F*  f3fHistCentVsInvMassLambda1point0negloose;
    TH3F*  f3fHistCentVsInvMassAntiLambda1point0negloose;
    
    TH3F*  f3fHistCentVsInvMassLambda1point0negtight;
    TH3F*  f3fHistCentVsInvMassAntiLambda1point0negtight;
    
    TH3F*  f3fHistCentVsInvMassLambda1point0posloose;
    TH3F*  f3fHistCentVsInvMassAntiLambda1point0posloose;
    
    TH3F*  f3fHistCentVsInvMassLambda1point0postight;
    TH3F*  f3fHistCentVsInvMassAntiLambda1point0postight;

    
    Float_t fCentrality;
    Int_t fTreeVariableLeastNbrCrossedRows;
    Float_t fTreeVariableLeastRatioCrossedRowsOverFindable;
    
    UInt_t fEvSel;
    Int_t  fNptBins;
    
    THnSparse *fPtBinNplusNminusCh;
    THnSparse *fPtBinNplusNminusChnsigtight;
    THnSparse *fPtBinNplusNminusChnegloose;
    THnSparse *fPtBinNplusNminusChnegtight;
    THnSparse *fPtBinNplusNminusChposloose;
    THnSparse *fPtBinNplusNminusChpostight;
    
    Int_t    GetPtBin(Double_t pt);
    ClassDef(AliAnalysisTaskNetLambdaTrad,5);
};


#endif

