// For: Net Lambda fluctuation analysis via traditional method
// By: Ejiro Umaka Apr 2018
//updated 0ct 1

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
    
    void SetIsMC(Bool_t val){fIsMC = val;};
    Bool_t GetIsMC(){return fIsMC;};
    void SetIsAOD(Bool_t val){fIsAOD = val;};
    Bool_t GetIsAOD(){return fIsAOD;};
    void SetEventSelection(UInt_t val) {fEvSel = val;}
    
protected:
    AliAnalysisTaskNetLambdaTrad(const  AliAnalysisTaskNetLambdaTrad &task);
    AliAnalysisTaskNetLambdaTrad& operator=(const  AliAnalysisTaskNetLambdaTrad &task);
    
    AliESDEvent* fESD;
    AliAODEvent* fAOD;
    AliPIDResponse* fPIDResponse;
    AliEventCuts fEventCuts;
    TList* fListHist;
    TTree* fTreeV0;
    
    TH1D*  fHistEventCounter;
    TH1D*  fHistCentrality;
    
    
    TH2F*  f2fHistGenCentVsPtLambda;
    TH2F*  f2fHistGenCentVsPtAntiLambda;
    TH2F*  f2fHistRecCentVsPtLambda;
    TH2F*  f2fHistRecCentVsPtAntiLambda;
    TH2F*  f2fHistInvMassVsPtLambda;
    TH2F*  f2fHistInvMassVsPtLambdaRec;
    TH2F*  f2fHistInvMassVsPtAntiLambda;
    TH2F*  f2fHistInvMassVsPtAntiLambdaRec;
    TH2F*  f2fHistRecPrimariesCentVsPtLambda;
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambda;
    TH2F*  f2fHistmassctLambda;
    TH2F*  f2fHistmassctAntiLambda;
    TH2F*  f2fHistLambdaSecFromWeakDecay;
    TH2F*  f2fHistAntiLambdaSecFromWeakDecay;
    TH2F*  f2fHistLambdaMaterial;
    TH2F*  f2fHistAntiLambdaMaterial;
    TH2F*  f2fHistLambdaMisId;
    TH2F*  f2fHistAntiLambdaMisId;
    TH2F*  f2fHistLRecstat;
    TH2F*  f2fHistARecstat;
    TH2F*  f2fHistLGenstat;
    TH2F*  f2fHistAGenstat;
    
    Float_t fCentrality;
    Int_t fTreeVariablePID;
    Int_t fTreeVariablePIDPositive;
    Int_t fTreeVariablePIDNegative;
    Int_t fTreeVariableLeastNbrCrossedRows;
    Float_t fTreeVariableLeastRatioCrossedRowsOverFindable;
    
    Float_t fTreeVariableInvMassLambda;
    Float_t fTreeVariableInvMassAntiLambda;
    Float_t fTreeVariableDcaV0Daughters;
    Float_t fTreeVariableDcaV0ToPrimVertex;
    Float_t fTreeVariableDcaPosToPrimVertex;
    Float_t fTreeVariableDcaNegToPrimVertex;
    
    Bool_t fIsMC;
    Bool_t fIsAOD;
    UInt_t fEvSel;
    Int_t  fNptBins;
    
    THnSparse *fPtBinNplusNminusCh;
    THnSparse *fPtBinNplusNminusChCut;
    THnSparse *fPtBinNplusNminusChTruth;
    Int_t    GetPtBin(Double_t pt);
    
    
    ClassDef(AliAnalysisTaskNetLambdaTrad,5);
};


#endif
