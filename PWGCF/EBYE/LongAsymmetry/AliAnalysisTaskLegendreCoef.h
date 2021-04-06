////////////////////////////////////////////////////////
// AliAnalysisTaskLegendreCoef:
// Description: Analysis task computes the background 
// and extracts legendre coefficients from eta dist
// Author: Raquel Quishpe (raquel.quishpe@cern.ch)
////////////////////////////////////////////////////////

#ifndef AliAnalysisTaskLegendreCoef_H
#define AliAnalysisTaskLegendreCoef_H

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
    void            SetEtaLimit(Double_t etalimit) { fEta = etalimit; }
    void            SetPileUpRead(Bool_t flag) {fIsPileUpCuts = flag;}
    void            SetBuildBackground(Bool_t flag) {fIsBuildBG = flag; }
    void            SetBuildLegendre(Bool_t flag) {fIsBuildLG = flag; }
    void            GetPosBackground(TH2D* hist) { fPosBackgroundHist = hist; }
    void            GetNegBackground(TH2D* hist) { fNegBackgroundHist = hist; }
    void            GetChargedBackground(TH2D* hist) { fChargedBackgroundHist = hist; }
    void            GetMCPosBackground(TH2D* hist) { fMCPosBackgroundHist = hist; }
    void            GetMCNegBackground(TH2D* hist) { fMCNegBackgroundHist = hist; }
    void            GetMCChargedBackground(TH2D* hist) { fMCChargedBackgroundHist = hist; }

  private:
    AliAODEvent*    fAOD;                 //! input event
    TList*          fOutputList;          //! output list
  
    Bool_t fIsMC; //MC flag
    Double_t fChi2DoF; //limit for chi2
    Int_t fTPCNcls; //limit for TPC Ncls
    Double_t fPtmin; //min PT
    Double_t fPtmax; //max PT
    Double_t fEta; //max eta
    Bool_t fIsPileUpCuts; //pile up cuts flag
    Bool_t fIsBuildBG; //build background flag
    Bool_t fIsBuildLG; //calculates the legendre coefficients
    TH2D* fPosBackgroundHist; //input background histogram for positive tracks (eta,cent)
    TH2D* fNegBackgroundHist; //input background histogram for negative tracks (eta,cent)
    TH2D* fChargedBackgroundHist; //input background histogram for positive and negative tracks (eta,cent)
    TH2D* fMCPosBackgroundHist; //input MC background histogram for positive tracks (eta,cent)
    TH2D* fMCNegBackgroundHist; //input MC background histogram for negative tracks (eta,cent)
    TH2D* fMCChargedBackgroundHist; //input MC background histogram for positive and negative tracks (eta,cent)
    AliEventCuts fEventCuts;

    AliAnalysisTaskLegendreCoef(const AliAnalysisTaskLegendreCoef&); // not implemented
    AliAnalysisTaskLegendreCoef& operator=(const AliAnalysisTaskLegendreCoef&); // not implemented

    ClassDef(AliAnalysisTaskLegendreCoef, 2);
};

#endif
