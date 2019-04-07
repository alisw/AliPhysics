
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

    TH2F*  f2fHistRecPrimariesCentVsPtLambdaFour;
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambdaFour;
    TH2F*  f2fHistRecPrimariesCentVsPtLambdaThree;
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambdaThree;
    TH2F*  f2fHistRecPrimariesCentVsPtLambdaTwo;
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambdaTwo;
    TH2F*  f2fHistRecPrimariesCentVsPtLambdaOne;
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambdaOne;
    
    TH2F*  f2fHistRecPrimariesCentVsPtLambdaFourSigthree;
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambdaFourSigthree;
    TH2F*  f2fHistRecPrimariesCentVsPtLambdaThreeSigthree;
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambdaThreeSigthree;
    TH2F*  f2fHistRecPrimariesCentVsPtLambdaTwoSigthree;
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambdaTwoSigthree;
    TH2F*  f2fHistRecPrimariesCentVsPtLambdaOneSigthree;
    TH2F*  f2fHistRecPrimariesCentVsPtAntiLambdaOneSigthree;
    
    TH2F*  f2fHistLambdaMisIdFour;
    TH2F*  f2fHistAntiLambdaMisIdFour;
    TH2F*  f2fHistLambdaMisIdThree;
    TH2F*  f2fHistAntiLambdaMisIdThree;
    TH2F*  f2fHistLambdaMisIdTwo;
    TH2F*  f2fHistAntiLambdaMisIdTwo;
    TH2F*  f2fHistLambdaMisIdOne;
    TH2F*  f2fHistAntiLambdaMisIdOne;
    

    TH2F*  f2fHistLambdaMisIdFourSigthree;
    TH2F*  f2fHistAntiLambdaMisIdFourSigthree;
    TH2F*  f2fHistLambdaMisIdThreeSigthree;
    TH2F*  f2fHistAntiLambdaMisIdThreeSigthree;
    TH2F*  f2fHistLambdaMisIdTwoSigthree;
    TH2F*  f2fHistAntiLambdaMisIdTwoSigthree;
    TH2F*  f2fHistLambdaMisIdOneSigthree;
    TH2F*  f2fHistAntiLambdaMisIdOneSigthree;
    
    TH2F*  f2fHistRecSecCentVsPtLambdaFourSigthree;
    TH2F*  f2fHistRecSecCentVsPtAntiLambdaFourSigthree;
    TH2F*  f2fHistRecMaterialCentVsPtLambdaFourSigthree;
    TH2F*  f2fHistRecMaterialCentVsPtAntiLambdaFourSigthree;

    TH3F*  f3fHistLambdafromXiFour;
    TH3F*  f3fHistAntiLambdafromXiFour;
    TH3F*  f3fHistLambdafromXiThree;
    TH3F*  f3fHistAntiLambdafromXiThree;
    TH3F*  f3fHistLambdafromXiTwo;
    TH3F*  f3fHistAntiLambdafromXiTwo;
    TH3F*  f3fHistLambdafromXiOne;
    TH3F*  f3fHistAntiLambdafromXiOne;
    
    TH3F*  f3fHistLambdafromXiFourSigthree;
    TH3F*  f3fHistAntiLambdafromXiFourSigthree;
    TH3F*  f3fHistLambdafromXiThreeSigthree;
    TH3F*  f3fHistAntiLambdafromXiThreeSigthree;
    TH3F*  f3fHistLambdafromXiTwoSigthree;
    TH3F*  f3fHistAntiLambdafromXiTwoSigthree;
    TH3F*  f3fHistLambdafromXiOneSigthree;
    TH3F*  f3fHistAntiLambdafromXiOneSigthree;
    
    
    TH3F*  f3fHistCentInvMassVsPtLambdaRecFour;
    TH3F*  f3fHistCentInvMassVsPtAntiLambdaRecFour;
    TH3F*  f3fHistCentInvMassVsPtLambdaRecThree;
    TH3F*  f3fHistCentInvMassVsPtAntiLambdaRecThree;
    TH3F*  f3fHistCentInvMassVsPtLambdaRecTwo;
    TH3F*  f3fHistCentInvMassVsPtAntiLambdaRecTwo;
    TH3F*  f3fHistCentInvMassVsPtLambdaRecOne;
    TH3F*  f3fHistCentInvMassVsPtAntiLambdaRecOne;
    
    TH3F*  f3fHistCentInvMassVsPtLambdaRecFourSigthree;
    TH3F*  f3fHistCentInvMassVsPtAntiLambdaRecFourSigthree;
    TH3F*  f3fHistCentInvMassVsPtLambdaRecThreeSigthree;
    TH3F*  f3fHistCentInvMassVsPtAntiLambdaRecThreeSigthree;
    TH3F*  f3fHistCentInvMassVsPtLambdaRecTwoSigthree;
    TH3F*  f3fHistCentInvMassVsPtAntiLambdaRecTwoSigthree;
    TH3F*  f3fHistCentInvMassVsPtLambdaRecOneSigthree;
    TH3F*  f3fHistCentInvMassVsPtAntiLambdaRecOneSigthree;
    
    //tagged
    TH3F*  f3fHistCentInvMassVsPtLambdaRecFourSigthreeUntag;
    TH3F*  f3fHistCentInvMassVsPtAntiLambdaRecFourSigthreeUntag;
    TH3F*  f3fHistCentInvMassVsPtLambdaRecFourSigthreeUntagCut;
    TH3F*  f3fHistCentInvMassVsPtAntiLambdaRecFourSigthreeUntagCut;
    
    //VO
    TH3F*  f3fHistCentInvMassVsPtLambdaVOFourSigthree;
    TH3F*  f3fHistCentInvMassVsPtLambdaVOFourSigthreeCut;
    TH3F*  f3fHistCentInvMassVsPtAntiLambdaVOFourSigthree;
    TH3F*  f3fHistCentInvMassVsPtAntiLambdaVOFourSigthreeCut;
    

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
    THnSparse *fPtBinNplusNminusChTagFour;
    THnSparse *fPtBinNplusNminusChTagThree;
    THnSparse *fPtBinNplusNminusChTagTwo;
    THnSparse *fPtBinNplusNminusChTagOne;
    
    THnSparse *fPtBinNplusNminusChTagFourSigThree;
    THnSparse *fPtBinNplusNminusChTagThreeSigThree;
    THnSparse *fPtBinNplusNminusChTagTwoSigThree;
    THnSparse *fPtBinNplusNminusChTagOneSigThree;
    THnSparse *fPtBinNplusNminusChUnTagFour;
    THnSparse *fPtBinNplusNminusChVOFour;

    
    THnSparse *fPtBinNplusNminusChBKGM;
    THnSparse *fPtBinNplusNminusChLF;
    THnSparse *fPtBinNplusNminusChRT;
    
    
    Int_t    GetPtBin(Double_t pt);
    
    
    
    
    ClassDef(AliAnalysisTaskNetLambdaTrad,5);
};


#endif


