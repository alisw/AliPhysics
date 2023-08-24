

// For: Net Lambda fluctuation analysis via traditional method
// By: Ejiro Umaka Apr 2018
//Update Mar 2019

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
    virtual ~AliAnalysisTaskNetLambdaMCTrad ();
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
    
    TH2F*  f2fHistGenCentVsPtLambda;
    TH2F*  f2fHistGenCentVsPtAntiLambda;
    TH2F*  f2fHistXiPlus;
    TH2F*  f2fHistXiMinus;

    
    TH2F*  f2fHistRecPrimariesCentVsPtLambdaFourSigthree;
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthree;
    TH3F*  f3fHistLambdafromXiFourSigthree;
    TH3F*  f3fHistAntiLambdafromXiFourSigthree;
    
    TH2F*  f2fHistRecPrimariesCentVsPtLambdaFourSigthreensigtight;
    TH3F*  f3fHistLambdafromXiFourSigthreensigtight;
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreensigtight;
    TH3F*  f3fHistAntiLambdafromXiFourSigthreensigtight;
    
    TH2F*  f2fHistRecPrimariesCentVsPtLambdaFourSigthreenegloose;
    TH3F*  f3fHistLambdafromXiFourSigthreenegloose;
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreenegloose;
    TH3F*  f3fHistAntiLambdafromXiFourSigthreenegloose;

    TH2F*  f2fHistRecPrimariesCentVsPtLambdaFourSigthreenegtight;
    TH3F*  f3fHistLambdafromXiFourSigthreenegtight;
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreenegtight;
    TH3F*  f3fHistAntiLambdafromXiFourSigthreenegtight;
    
    TH2F*  f2fHistRecPrimariesCentVsPtLambdaFourSigthreeposloose;
    TH3F*  f3fHistLambdafromXiFourSigthreeposloose;
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreeposloose;
    TH3F*  f3fHistAntiLambdafromXiFourSigthreeposloose;

    TH2F*  f2fHistRecPrimariesCentVsPtLambdaFourSigthreepostight;
    TH3F*  f3fHistLambdafromXiFourSigthreepostight;
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthreepostight;
    TH3F*  f3fHistAntiLambdafromXiFourSigthreepostight;
    
    
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
    THnSparse *fPtBinNplusNminusChnsigtight;
    THnSparse *fPtBinNplusNminusChnegloose;
    THnSparse *fPtBinNplusNminusChnegtight;
    THnSparse *fPtBinNplusNminusChposloose;
    THnSparse *fPtBinNplusNminusChpostight;
    THnSparse *fPtBinNplusNminusCh;


    
    Int_t    GetPtBin(Double_t pt);
    Double_t MyRapidity(Double_t rE, Double_t rPz) const;
    
    
    ClassDef(AliAnalysisTaskNetLambdaMCTrad,4);
};


#endif





