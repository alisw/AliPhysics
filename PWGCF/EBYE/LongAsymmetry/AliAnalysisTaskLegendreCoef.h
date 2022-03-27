////////////////////////////////////////////////////////
// AliAnalysisTaskLegendreCoef:
// Description: Analysis task computes the background 
// and extracts legendre coefficients from eta dist
// Author: Raquel Quishpe (raquel.quishpe@cern.ch)
////////////////////////////////////////////////////////

#ifndef AliAnalysisTaskLegendreCoef_H
#define AliAnalysisTaskLegendreCoef_H
class TH1;
class THn;
class TH1F;
class TH2D;
class TH3D;
class TList;
class TTree;

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#define PI 3.1415927

class AliAnalysisTaskLegendreCoef : public AliAnalysisTaskSE
{
  public:
                    AliAnalysisTaskLegendreCoef();
                    AliAnalysisTaskLegendreCoef(const char *name);
    virtual         ~AliAnalysisTaskLegendreCoef();

    virtual void    UserCreateOutputObjects();
    virtual void    UserExec(Option_t* option);
    virtual void    Terminate(Option_t* option);


    void            SetMCRead(Bool_t flag) { fIsMC = flag; }
    void            SetChi2DoF(Double_t Chi2DoF) { fChi2DoF = Chi2DoF; }
    void            SetNclTPC(Int_t ncl) { fTPCNcls = ncl; }
    void            SetPtLimits(Double_t ptmin, Double_t ptmax) { fPtmin = ptmin; fPtmax=ptmax; }
    void            SetEtaMinLimit(Double_t etalowlimit) { fEtaMin = etalowlimit; }
    void            SetEtaMaxLimit(Double_t etauplimit) { fEtaMax = etauplimit; }
    void            SetFilterBit(Int_t filterbit) { fBit = filterbit; }
    void            SetPileUpRead(Bool_t flag) {fIsPileUpCuts = flag;}
    void            SetBuildBackground(Bool_t flag) {fIsBuildBG = flag; }
    void            SetBuildLegendre(Bool_t flag) {fIsBuildLG = flag; }
    void            GetPosBackground(TH2D* hist) { fPosBackgroundHist = hist; }
    void            GetNegBackground(TH2D* hist) { fNegBackgroundHist = hist; }
    void            GetChargedBackground(TH2D* hist) { fChargedBackgroundHist = hist; }
    void            GetMCPosBackground(TH2D* hist) { fMCPosBackgroundHist = hist; }
    void            GetMCNegBackground(TH2D* hist) { fMCNegBackgroundHist = hist; }
    void            GetMCChargedBackground(TH2D* hist) { fMCChargedBackgroundHist = hist; }
    void            GetNeventsCentHist(TH1D* hist) { fNeventCentHist = hist; }
    void            SetGeneratorName(TString generator) {fGenName = generator; }
    void            SetPileUpLevel(Int_t level) {fPileUpLevel = level; }
    void            SetTPCNCrossedRows(UShort_t crossedrows) { fTPCNCrossedRows = crossedrows; }
    void            SetPVzMinLimit(Float_t pvzmin) {fPVzMin=pvzmin;}
    void            SetPVzMaxLimit(Float_t pvzmax) {fPVzMax=pvzmax;}
    void            SetPVzSign(Int_t sign) {fPVzSign=sign;}//-1 then negative pvz, +1 then positive pvz, 0 then absolute value (to test effect from the TPC membrane)
    void            SetNEtaBins(Int_t Netabins) {fNetabins=Netabins;}//default 16
    void            SetRunOnlyFB(Bool_t flag) {fIsRunFBOnly = flag;}//default 16

  private:
    Double_t GetSingleAnCoef(int order, TH1D *hist); //method to get direct an
    Double_t LegPol(int order, Double_t x);
    void    BuildBackground();
    void    BuildSignal();
    void    BuildCoefficients(TH1D *signal, TH1D *background, Float_t centrality, TString type);
    AliAODEvent*    fAOD;                 //! input event
    TList*          fOutputList;          //! output list

  
    Bool_t fIsMC; //MC flag
    Double_t fChi2DoF; //limit for chi2
    Int_t fTPCNcls; //limit for TPC Ncls
    Double_t fPtmin; //min PT
    Double_t fPtmax; //max PT
    Double_t fEtaMin; //min eta
    Double_t fEtaMax; //max eta
    Int_t fBit; //filter bit - 96 default
    Bool_t fIsPileUpCuts; //pile up cuts flag
    Bool_t fIsBuildBG; //build background flag
    Bool_t fIsBuildLG; //calculates the legendre coefficients
    TH2D* fPosBackgroundHist; //input background histogram for positive tracks (eta,cent)
    TH2D* fNegBackgroundHist; //input background histogram for negative tracks (eta,cent)
    TH2D* fChargedBackgroundHist; //input background histogram for positive and negative tracks (eta,cent)
    TH2D* fMCPosBackgroundHist; //input MC background histogram for positive tracks (eta,cent)
    TH2D* fMCNegBackgroundHist; //input MC background histogram for negative tracks (eta,cent)
    TH2D* fMCChargedBackgroundHist; //input MC background histogram for positive and negative tracks (eta,cent)
    TH1D* fNeventCentHist; //input Nevents vs centrality histogram 
    TString fGenName; //MC generator name
    Int_t fPileUpLevel;
    UShort_t  fTPCNCrossedRows;  
    Float_t fPVzMax; //max PVz
    Float_t fPVzMin; //min PVz
    Int_t fPVzSign; //sign of PVz
    Int_t fNetabins; //number of bins in eta
    Bool_t fIsRunFBOnly; //only filterbit cuts

    AliEventCuts fEventCuts;

    AliAnalysisTaskLegendreCoef(const AliAnalysisTaskLegendreCoef&); // not implemented
    AliAnalysisTaskLegendreCoef& operator=(const AliAnalysisTaskLegendreCoef&); // not implemented

    ClassDef(AliAnalysisTaskLegendreCoef, 2);
};

#endif
