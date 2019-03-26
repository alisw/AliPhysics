
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
    
    TH1D*  fHistEventCounter;
    TH1D*  fHistCentrality;
    
    
    
    TH3F*  f3fHistCentVsInvMassLambda1point6;
    TH3F*  f3fHistCentVsInvMassLambda1point0;
    TH3F*  f3fHistCentVsInvMassLambda0point6;
    TH3F*  f3fHistCentVsInvMassLambda0point2;
    
    TH3F*  f3fHistPtmassctLambda1point6;
    TH3F*  f3fHistPtmassctLambda1point0;
    TH3F*  f3fHistPtmassctLambda0point6;
    TH3F*  f3fHistPtmassctLambda0point2;

    TH3F*  f3fHistPtmassctLambdaPosOpoint2;
    TH3F*  f3fHistPtmassctLambdaPosOpoint4;
    TH3F*  f3fHistPtmassctLambdaPosOpoint6;
    TH3F*  f3fHistPtmassctLambdaPosOpoint8;
    
    ///ANTI-LAMBDA
    TH3F*  f3fHistCentVsInvMassAntiLambda1point6;
    TH3F*  f3fHistCentVsInvMassAntiLambda1point0;
    TH3F*  f3fHistCentVsInvMassAntiLambda0point6;
    TH3F*  f3fHistCentVsInvMassAntiLambda0point2;
    
    TH3F*  f3fHistPtmassctAntiLambda1point6;
    TH3F*  f3fHistPtmassctAntiLambda1point0;
    TH3F*  f3fHistPtmassctAntiLambda0point6;
    TH3F*  f3fHistPtmassctAntiLambda0point2;
    
    TH3F*  f3fHistPtmassctAntiLambdaPosOpoint2;
    TH3F*  f3fHistPtmassctAntiLambdaPosOpoint4;
    TH3F*  f3fHistPtmassctAntiLambdaPosOpoint6;
    TH3F*  f3fHistPtmassctAntiLambdaPosOpoint8;

    Float_t fCentrality;
    
    Int_t fTreeVariableLeastNbrCrossedRows;
    Float_t fTreeVariableLeastRatioCrossedRowsOverFindable;
    
  
    UInt_t fEvSel;
    Int_t  fNptBins;

    
    THnSparse *fPtBinNplusNminusChEtaFour;
    THnSparse *fPtBinNplusNminusChEtaThree;
    THnSparse *fPtBinNplusNminusChEtaTwo;
    THnSparse *fPtBinNplusNminusChEtaOne;
    THnSparse *fPtBinNplusNminusChPosEtaFour;
    THnSparse *fPtBinNplusNminusChPosEtaThree;
    THnSparse *fPtBinNplusNminusChPosEtaTwo;
    THnSparse *fPtBinNplusNminusChPosEtaOne;
    THnSparse *fPtBinNplusNminusChBproxy;

    
    Int_t    GetPtBin(Double_t pt);
    
    ClassDef(AliAnalysisTaskNetLambdaTrad,5);
};


#endif


