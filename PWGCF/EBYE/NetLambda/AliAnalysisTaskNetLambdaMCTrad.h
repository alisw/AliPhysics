
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
    
    TH2F*  f2fHistGenCentVsPtLambda;
    TH2F*  f2fHistGenCentVsPtAntiLambda;
    TH2F*  f2fHistXiPlus;
    TH2F*  f2fHistXiMinus;
    TH2F*  f2fHistLRecstat;
    TH2F*  f2fHistARecstat;
    TH2F*  f2fHistLGenstat;
    TH2F*  f2fHistAGenstat;


    
    TH2F*  f2fHistRecPrimariesCentVsPtLambdaFourSigthree;
    TH2F*  f2fHistRecPrimariesCentVsPtLambdaFourloose;
    TH2F*  f2fHistRecPrimariesCentVsPtLambdaFourtight;
    TH2F*  f2fHistRecPrimariesCentVsPtLambdaFournegloose;
    TH2F*  f2fHistRecPrimariesCentVsPtLambdaFournegtight;
    TH2F*  f2fHistRecPrimariesCentVsPtLambdaFourposloose;
    TH2F*  f2fHistRecPrimariesCentVsPtLambdaFourpostight;
    
    
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthree;
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambdaFourloose;
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambdaFourtight;
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambdaFournegloose;
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambdaFournegtight;
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambdaFourposloose;
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambdaFourpostight;


    
    TH2F*  f2fHistRecSecCentVsPtLambdaFourSigthree;
    TH2F*  f2fHistRecSecCentVsPtLambdaFourloose;
    TH2F*  f2fHistRecSecCentVsPtLambdaFourtight;
    TH2F*  f2fHistRecSecCentVsPtAntiLambdaFourSigthree;
    TH2F*  f2fHistRecSecCentVsPtAntiLambdaFourloose;
    TH2F*  f2fHistRecSecCentVsPtAntiLambdaFourtight;
    
    
    TH2F*  f2fHistRecMaterialCentVsPtLambdaFourSigthree;
    TH2F*  f2fHistRecMaterialCentVsPtLambdaFourloose;
    TH2F*  f2fHistRecMaterialCentVsPtLambdaFourtight;
    TH2F*  f2fHistRecMaterialCentVsPtAntiLambdaFourSigthree;
    TH2F*  f2fHistRecMaterialCentVsPtAntiLambdaFourloose;
    TH2F*  f2fHistRecMaterialCentVsPtAntiLambdaFourtight;
    

    TH3F*  f3fHistLambdafromXiFourSigthree;
    TH3F*  f3fHistLambdafromXiFourloose;
    TH3F*  f3fHistLambdafromXiFourtight;
    TH3F*  f3fHistLambdafromXiFournegloose;
    TH3F*  f3fHistLambdafromXiFournegtight;
    TH3F*  f3fHistLambdafromXiFourposloose;
    TH3F*  f3fHistLambdafromXiFourpostight;
    TH3F*  f3fHistAntiLambdafromXiFourSigthree;
    TH3F*  f3fHistAntiLambdafromXiFourloose;
    TH3F*  f3fHistAntiLambdafromXiFourtight;
    TH3F*  f3fHistAntiLambdafromXiFournegloose;
    TH3F*  f3fHistAntiLambdafromXiFournegtight;
    TH3F*  f3fHistAntiLambdafromXiFourposloose;
    TH3F*  f3fHistAntiLambdafromXiFourpostight;


    TH3F*  f3fHistCentInvMassVsPtLambdaRecFourUntagloose;
    TH3F*  f3fHistCentInvMassVsPtLambdaRecFourUntagtight;
    TH3F*  f3fHistCentInvMassVsPtLambdaRecFourSigthreeUntag;
    TH3F*  f3fHistCentInvMassVsPtLambdaRecFourUntagCutloose;
    TH3F*  f3fHistCentInvMassVsPtLambdaRecFourUntagCuttight;
    TH3F*  f3fHistCentInvMassVsPtLambdaRecFourSigthreeUntagCut;
    TH3F*  f3fHistCentInvMassVsPtAntiLambdaRecFourSigthreeUntag;
    TH3F*  f3fHistCentInvMassVsPtAntiLambdaRecFourSigthreeUntagCut;
    TH3F*  f3fHistCentInvMassVsPtAntiLambdaRecFourUntagloose;
    TH3F*  f3fHistCentInvMassVsPtAntiLambdaRecFourUntagCutloose;
    TH3F*  f3fHistCentInvMassVsPtAntiLambdaRecFourUntagtight;
    TH3F*  f3fHistCentInvMassVsPtAntiLambdaRecFourUntagCuttight;

    



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
    
        THnSparse *fPtBinNplusNminusChUNTagFourlooseBKG;
        THnSparse *fPtBinNplusNminusChUNTagFourloose;
        THnSparse *fPtBinNplusNminusChUNTagFourTightBKG;
        THnSparse *fPtBinNplusNminusChUNTagFourTight;
        THnSparse *fPtBinNplusNminusChUNTagFourBKG;
        THnSparse *fPtBinNplusNminusChUNTagFour;
    
    

    Int_t    GetPtBin(Double_t pt);

    ClassDef(AliAnalysisTaskNetLambdaTrad,5);
};


#endif


