
// For: Net Lambda fluctuation analysis via traditional method
// By: Ejiro Umaka Apr 2018
//Update jul 24

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
    TH1F*  fHistPA;
    TH1F*  fHistPAparLambda;
    TH1F*  fHistPAparLambdabar;
    
    TH2F*  f2fHistGenCentVsPtLambda;
    TH2F*  f2fHistGenCentVsPtAntiLambda;
    TH2F*  f2fHistXiPlus;
    TH2F*  f2fHistXiMinus;
    
    TH2F*  f2fHistGenCentVsPtLambdaRap;
    TH2F*  f2fHistGenCentVsPtAntiLambdaRap;
    TH2F*  f2fHistXiPlusRap;
    TH2F*  f2fHistXiMinusRap;
    
    
    
    TH2F*  f2fHistRecPrimariesCentVsPtLambdaFourSigthree;
    
    TH2F*  f2fHistRecPrimariesCentVsPtLambdaFourSigthreeRap;
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreeRap;
    
    
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthree;
    
    TH2F*  f2fHistRecSecCentVsPtLambdaFourSigthree;
    TH2F*  f2fHistRecSecCentVsPtAntiLambdaFourSigthree;
    
    TH2F*  f2fHistRecMatCentVsPtLambdaFourSigthree;
    TH2F*  f2fHistRecMatCentVsPtAntiLambdaFourSigthree;
    
    
    TH3F*  f3fHistLambdafromXiFourSigthree;

    TH3F*  f3fHistAntiLambdafromXiFourSigthree;
    
    TH3F*  f3fHistLambdafromXiFourSigthreeRap;
    TH3F*  f3fHistAntiLambdafromXiFourSigthreeRap;
    
    
    TH3F*  f3fHistCentInvMassVsPtLambdaRecFourSigthreeUntag;
    TH3F*  f3fHistCentInvMassVsPtLambdaRecFourSigthreeUntagCut;
    TH3F*  f3fHistCentInvMassVsPtAntiLambdaRecFourSigthreeUntag;
    TH3F*  f3fHistCentInvMassVsPtAntiLambdaRecFourSigthreeUntagCut;
    
    TH3F*  f3fHistCentInvMassVsPtLambdaRecFourSigthree;
    TH3F*  f3fHistCentInvMassVsPtAntiLambdaRecFourSigthree;
    
    
    
    
    
    Float_t fCentrality;
    Int_t fTreeVariablePID;
    Int_t fTreeVariablePIDParent;
    Int_t fTreeVariablePIDPositive;
    Int_t fTreeVariablePIDNegative;
    Int_t fTreeVariablePrimaryStatusMother;
    
    Int_t fTreeVariableLeastNbrCrossedRows;
    Float_t fTreeVariableLeastRatioCrossedRowsOverFindable;
    
    Bool_t fIsMC;
    UInt_t fEvSel;
    Int_t  fNptBins;
    
    THnSparse *fPtBinNplusNminusChTruth;
    
    THnSparse *fPtBinNplusNminusChRec;
    THnSparse *fPtBinNplusNminusChRecTag;
    
    Int_t    GetPtBin(Double_t pt);
    Double_t MyRapidity(Double_t rE, Double_t rPz) const;
    
    
    ClassDef(AliAnalysisTaskNetLambdaMCTrad,4);
};


#endif



