
// For: Net Lambda fluctuation analysis via traditional method
// By: Ejiro Umaka Apr 2018
//update Feb 2019

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
    virtual ~AliAnalysisTaskNetLambdaTrad(){};
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
    TTree* fTreeV0;
    
    TH1D*  fHistEventCounter;
    TH1D*  fHistCentrality;
    
    TH2F*  f2fHistRecCentVsPtLambda;
    TH2F*  f2fHistRecCentVsPtAntiLambda;
    TH2F*  f2fHistInvMassVsPtLambda;
    TH2F*  f2fHistInvMassVsPtAntiLambda;
    TH1F*  f1fHistmassctLambda;
    TH1F*  f1fHistmassctAntiLambda;
    TH2F*  f2fHistPtmassctLambda;
    TH2F* f2fHistPtmassctAntiLambda;
    TH2F* f2fHistCentVsInvMassLambda;
    TH2F* f2fHistCentVsInvMassAntiLambda;

    
    Float_t fCentrality;
    
    Int_t fTreeVariableLeastNbrCrossedRows;
    Float_t fTreeVariableLeastRatioCrossedRowsOverFindable;
    
    Float_t fTreeVariableInvMassLambda;
    Float_t fTreeVariableInvMassAntiLambda;
    Float_t fTreeVariableDcaV0Daughters;
    Float_t fTreeVariableDcaV0ToPrimVertex;
    Float_t fTreeVariableDcaPosToPrimVertex;
    Float_t fTreeVariableDcaNegToPrimVertex;
    Float_t fTreeVariableNsigmaPosProton;
    Float_t fTreeVariableNsigmaNegProton;
    Float_t fTreeVariableCentrality;
    
    UInt_t fEvSel;
    Int_t  fNptBins;
    
    THnSparse *fPtBinNplusNminusChALL;
    THnSparse *fPtBinNplusNminusChCut;
    Int_t    GetPtBin(Double_t pt);
    
    ClassDef(AliAnalysisTaskNetLambdaTrad,5);
};


#endif


