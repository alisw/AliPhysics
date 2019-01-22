
// For: Net Lambda fluctuation analysis via traditional method
// By: Ejiro Umaka Apr 2018

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
    TH2F*  f2fHistInvMassVsPtAntiLambda;
    TH2F*  f2fHistRecPrimariesCentVsPtLambda;
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambda;
    TH1F*  f1fHistmassctLambda;
    TH1F*  f1fHistmassctAntiLambda;
    TH2F*  f2fHistLambdaSecFromWeakDecay;
    TH2F* f2fHistAntiLambdaSecFromWeakDecay;
    TH2F*  f2fHistLambdaMisId;
    TH2F* f2fHistAntiLambdaMisId;
    TH2F*  f2fHistAntiLambdaRecPt;
    TH2F* f2fHistLambdaRecPt;
    TH2F*  f2fHistInvMassVsPtLambdaRec;
    TH2F* f2fHistInvMassVsPtAntiLambdaRec;
    TH2F*  f2fHistPtmassctLambda;
    TH2F* f2fHistPtmassctAntiLambda;
    TH2F*  f2fHistLambdaMaterial;
    TH2F* f2fHistAntiLambdaMaterial;
    TH2F*  f2fHistLRecstat;
    TH2F* f2fHistARecstat;
    TH2F*  f2fHistLGenstat;
    TH2F* f2fHistAGenstat;
    TH2F*  f2fHistXiPlus;
    TH2F* f2fHistXiMinus;
    TH3F*  f2fHistLambdafromXi;
    TH3F* f2fHistAntiLambdafromXi;
    
    
    
    Float_t fCentrality;
    Int_t fTreeVariablePID;
    Int_t fTreeVariablePIDParent;
    Int_t fTreeVariablePIDPositive;
    Int_t fTreeVariablePIDNegative;
    Int_t fTreeVariablePrimaryStatusMother;

    Int_t fTreeVariableLeastNbrCrossedRows;
    Float_t fTreeVariableLeastRatioCrossedRowsOverFindable;
    
    Float_t fTreeVariableInvMassLambda;
    Float_t fTreeVariableInvMassAntiLambda;
    Float_t fTreeVariablePtMC;
    Float_t fTreeVariablePtMother;
    Float_t fTreeVariableDcaV0Daughters;
    Float_t fTreeVariableDcaV0ToPrimVertex;
    Float_t fTreeVariableDcaPosToPrimVertex;
    Float_t fTreeVariableDcaNegToPrimVertex;
    Float_t fTreeVariableNsigmaPosProton;
    Float_t fTreeVariableNsigmaNegProton;
    Float_t fTreeVariableCentrality;
    
    
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


