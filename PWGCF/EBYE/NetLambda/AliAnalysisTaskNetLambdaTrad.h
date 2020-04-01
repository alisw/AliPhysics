
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
    TH1D*  hpVz;
    TH2F*  hlmasseta;
    TH2F*  hamasseta;
 
    
    TH2F*  hetaVO;
    TH2F*  hpetaVO;
    TH2F*  hnetaVO;
    TH2F*  hpxy;
    TH2F*  hprow;
    TH2F*  hnrow;
    TH2F*  hposp;
    TH2F*  hnegp;
    TH2F*  hpospi;
    TH2F*  hnegpi;
    
    TH2F*  lhDCAd;
    TH2F*  ahDCAd;
    TH2F*  lhV0rad;
    TH2F*  ahV0rad;
    
    TH2F*  lhCosp;
    TH2F*  ahCosp;
    TH2F*  lhV0tPV;
    TH2F*  ahV0tPV;
    
    TH2F*  lhntPV;
    TH2F*  ahntPV;
    TH2F*  lhptPV;
    TH2F*  ahptPV;
    
    TH3F*  f3fHistCentVsInvMassLambda1point0;
    TH3F*  f3fHistCentVsInvMassAntiLambda1point0;
    TH3F*  f3fHistCentVsInvMassLambda0point3;
    TH3F*  f3fHistCentVsInvMassAntiLambda0point3;
    
    Float_t fCentrality;
    Int_t fTreeVariableLeastNbrCrossedRows;
    Float_t fTreeVariableLeastRatioCrossedRowsOverFindable;
    
    UInt_t fEvSel;
    Int_t  fNptBins;
    
    THnSparse *fPtBinNplusNminusCh;
    THnSparse *fPtBinNplusNminusCh03;
    Int_t    GetPtBin(Double_t pt);
    ClassDef(AliAnalysisTaskNetLambdaTrad,5);
};


#endif

