
// For: Net Lambda fluctuation analysis via traditional method
// By: Ejiro Umaka Apr 2018

#ifndef AliAnalysisTaskNetLambdaMCTrad_h
#define AliAnalysisTaskNetLambdaMCTrad_h



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

class AliAnalysisTaskNetLambdaMCTrad : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskNetLambdaMCTrad(const char* name="AliAnalysisTaskNetLambdaMCTrad");
    virtual ~AliAnalysisTaskNetLambdaMCTrad(){};
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    
    void SetIsMC(Bool_t val){fIsMC = val;};
    Bool_t GetIsMC(){return fIsMC;};
    void SetEventSelection(UInt_t val) {fEvSel = val;}
    
protected:
    AliAnalysisTaskNetLambdaMCTrad(const  AliAnalysisTaskNetLambdaMCTrad &task);
    AliAnalysisTaskNetLambdaMCTrad& operator=(const  AliAnalysisTaskNetLambdaMCTrad &task);
    
    AliESDEvent* fESD;
    AliPIDResponse* fPIDResponse;
    AliEventCuts fEventCuts;
    TList* fListHist;
    
    TH1D*  fHistEventCounter;
    TH1D*  fHistCentrality;
    
    TH3F*  f3fHistCentVsInvMassLambda1point0;
    TH3F*  f3fHistCentVsInvMassLambda1point0Masscut;
    
    
    
    TH3F*  f3fHistCentVsInvMassAntiLambda1point0;
    TH3F*  f3fHistCentVsInvMassAntiLambda1point0Masscut;
    
    
    
    
    Float_t fCentrality;
    Int_t fTreeVariableLeastNbrCrossedRows;
    Float_t fTreeVariableLeastRatioCrossedRowsOverFindable;
    
    
    UInt_t fEvSel;
    Int_t  fNptBins;
    
    
    THnSparse *fPtBinNplusNminusCh;
    
    
    
    Int_t    GetPtBin(Double_t pt);

    ClassDef(AliAnalysisTaskNetLambdaMCTrad,4);
};


#endif


