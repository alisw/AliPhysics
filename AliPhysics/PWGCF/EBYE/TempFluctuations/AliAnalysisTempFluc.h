#ifndef ALIANALYSISTEMPFLUC_H
#define ALIANALYSISTEMPFLUC_H

class TList;
class TH1F;
class TH2F;
class TH3F;
class AliAODEvent;
class AliAODTrack;
class AliPIDResponse;
class TClonesArray;
class TSpline;
class TString;

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class AliAnalysisTempFluc : public AliAnalysisTaskSE {
 public:
    AliAnalysisTempFluc();
    AliAnalysisTempFluc(const char *name);
    virtual ~AliAnalysisTempFluc();
    
    virtual void     UserCreateOutputObjects();
    virtual void     UserExec(Option_t *option);
    virtual void     Terminate(Option_t *);
    //-----------------------------------
    // NEW Functions to be declared here
    //-----------------------------------
     Bool_t CheckTPC(AliAODTrack *track);
     Bool_t CheckTOF(AliAODTrack * track);
    
    void CalEfficiencyMap();
    void SetAnalysisMode(TString mode) {AnalysisMode=mode;}
    void SetDataType(TString type) {DataType=type;}
    void SetCentralityEstimator(TString Estimator) {fCentralityEstimator=Estimator;}
    void SetCentralityCutL(Double_t CentralityCutL) {fCentralityCutL=CentralityCutL;}
    void SetCentralityCutH(Double_t CentralityCutH) {fCentralityCutH=CentralityCutH;}
    void SetgEtaL(Float_t EtaL) {gEtaL=EtaL;}
    void SetgEtaH(Float_t EtaH) {gEtaH=EtaH;}
    void SetgRapL(Float_t RapL) {gRapL=RapL;}
    void SetgRapH(Float_t RapH) {gRapH=RapH;}
    void SetPtCutL(Double_t PtLow) {fPtL = PtLow;}
    void SetPtCutH(Double_t PtHigh) {fPtH = PtHigh;}
    void SetFilterBit(Int_t Bit) {fBit = Bit;}
    void SetCutNSigmaTPC(Float_t CutNSigmaTPC) {fCutNSigmaTPC = CutNSigmaTPC;}
    void SetCutNSigmaTPCTOF(Float_t CutNSigmaTPCTOF) {fCutNSigmaTPCTOF = CutNSigmaTPCTOF;}
    void CalMeanPtEbyETPC(TH1* h1Pt);
    void CalMeanPtEbyETPCTOF(TH1* h1Pt);
    void CalMeanPtEbyETruth(TH1* h1Pt);
    TSpline *func_slope_expo(const Float_t fitmin, const Float_t fitmax);
    Double_t Getbeta(AliAODTrack *track);
    Double_t CalculateRapidity(AliAODTrack *track , Double_t mass);
    void SetArrayDimension(Int_t dim) {fsizearray = dim;}
    
    void Correct_efficiency(Bool_t correction){fcorrection = correction;}
    void SetCustomBinningptL(TString CustomBinning){fBinningptL=CustomBinning;}
    void SetCustomBinningptH(TString CustomBinning){fBinningptH=CustomBinning;}
    void SetEffcorectionfilePathName(TString efffilename) {fefffilename=efffilename;}

 private:
    Bool_t   AcceptEvent(AliAODEvent *event) const; // accept event
   
    TList           *fOutput;        // Output list
    TH1F            *fHistPtTPC;        // Pt spectrum
    TH1F            *fHistPtTruthdummy;
    TH3F            *fHistTruth;
    TH3F            *fHistRecTPC;
    TH3F            *fHistRecPrimaryTPC;
    TH3F            *fHistRecSecWDTPC;
    TH3F            *fHistRecSecMatTPC;
    TH3F            *fHistRecMisIdTPC;
  
    TH2F            *fHistPM;       //PLusMinus Distribution
    
    TH2F            *fHistdEdx;
    TH2F            *fHistdEdxSigma;
    TH2F            *fNSigmaTPC;
    TH2F            *fNSigmaTPCCut;
    
    TH2F            *fNSigmaTOF;
    TH2F            *fNSigmaTPCTOF;
    TH2F            *fNSigmaTPCTOFCut;
    TH1F            *fHistPtTPCTOF;
    TH3F            *fHistRecTPCTOF;
    TH3F            *fHistRecPrimaryTPCTOF;
    TH3F            *fHistRecSecWDTPCTOF;
    TH3F            *fHistRecSecMatTPCTOF;
    TH3F            *fHistRecMisIdTPCTOF;
    TH2F            *fHistTOF;
    TH1F            *fHistPtTruth;
    
    TH1F            *fEffptTPC;
    TH1F            *fEffptTOF;
    TH1F            *fEffptTPCTOF;
    
    //---------  Any local integer or float  ---------
    Int_t           fNEvt;
    Float_t DCAxy;
    Float_t DCAz;
    
    //-----------------------------------
    // NEW HISTO to be declared here
    //-----------------------------------
    
    TString         fCentralityEstimator;
    Float_t         fCentralityCutL;
    Float_t         fCentralityCutH;
    Float_t         gRapL;
    Float_t         gRapH;
    Float_t         gEtaL;
    Float_t         gEtaH;
    Float_t         fPtL;
    Float_t         fPtH;
    Int_t           fBit;
    Float_t         fCutNSigmaTPC;
    Float_t         fCutNSigmaTPCTOF;
    TClonesArray * mcTruth;
    TString fBinningptL;
    TString fBinningptH;
    TString AnalysisMode;
    TString DataType;
    AliPIDResponse *fPIDResponse;
    AliAODEvent            *fAOD;
    Bool_t fcorrection;
    TString fefffilename;
    
    TH1F          *fQAHist[38];
    TH1F          *MeanpTebyeTPC[10];
    TH1F          *Teff_ebyeTPC[10];
    TH1F          *MeanpTebyeTPCTOF[10];
    TH1F          *Teff_ebyeTPCTOF[10];
    TH1F          *MeanpTebyeTruth[10];
    TH1F          *Teff_ebyeTruth[10];
    TH1F          *No_pion_evt_TPC[10];
    TH1F          *No_pion_evt_TPCTOF[10];
    TH1F          *No_pion_evt_Truth[10];
    Int_t   fsizearray;
    Double_t         mass_pion;
    //-----------------
    AliAnalysisTempFluc(const AliAnalysisTempFluc&); // not implemented
    AliAnalysisTempFluc& operator=(const AliAnalysisTempFluc&); // not implemented
    
    ClassDef(AliAnalysisTempFluc, 1); // example of analysis
};

#endif

