/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// AliAnalysisMultPt:
// Description: Analysis task to get multiplicity
// and pT distributions
// Author: Negin Alizadehvandchali
// (negin.alizadehvandchali@cern.ch)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef AliAnalysisMultPt_H
#define AliAnalysisMultPt_H

class TH1;
class THn;
class TH1F;
class TH2D;
class TH3D;
class TList;
class TTree;

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"

class AliAnalysisMultPt : public AliAnalysisTaskSE
{
    public:
                        AliAnalysisMultPt();
                        AliAnalysisMultPt(const char *name);
        virtual         ~AliAnalysisMultPt();

        virtual void    UserCreateOutputObjects();
        virtual void    UserExec(Option_t* option);
        virtual void    Terminate(Option_t* option);

        
        void            SetMCRead(Bool_t flag) { fIsMC = flag; }
        void            SetPtLimits(Double_t ptmin, Double_t ptmax) { fPtmin = ptmin; fPtmax=ptmax; }
        void            SetEtaMinLimit(Double_t etalowlimit) { fEtaMin = etalowlimit; }
        void            SetEtaMaxLimit(Double_t etauplimit) { fEtaMax = etauplimit; }
        void            SetFilterBit(Int_t filterbit) { fBit = filterbit; }
        void            SetPVzMinLimit(Float_t pvzmin) {fPVzMin=pvzmin;}
        void            SetPVzMaxLimit(Float_t pvzmax) {fPVzMax=pvzmax;}
        void            SetChi2DoF(Double_t Chi2DoF) { fChi2DoF = Chi2DoF; }
        void            SetTPCNCrossedRows(UShort_t crossedrows) { fTPCNCrossedRows = crossedrows; }
        void            SetRunOnlyFB(Bool_t flag) {fIsRunFBOnly = flag;}//default 16
        void            SetNclTPC(Int_t ncl) { fTPCNcls = ncl; }
        void            SetPileUpRead(Bool_t flag) {fIsPileUpCuts = flag;}
        void            SetPileUpLevel(Int_t level) {fPileUpLevel = level; }
        void            SetGeneratorName(TString generator) {fGenName = generator; }
 
    private:
        void            BuildData();
        void            BuildMC();
        AliAODEvent     *fAOD;             // input event
        TList           *fOutputList;      // output list
        TH1D            *MultHist;         // Mult histogram
        TH2D            *MultPtHist;       // pT-Mult histogram
        TH2D            *MultPtHistRec;    // pT-Mult MC histogram
        TH2D            *MultPtHistGen;    // pT-Mult MC histogram
        TH2D            *MultHistRatio;    // Ratio of Gen mult vs rec mult histogram
        TH1D            *fptRatio;         //pt ratio
    
        Bool_t          fIsMC;             // MC flag
        Double_t        fPtmin;            // min PT
        Double_t        fPtmax;            // max PT
        Double_t        fEtaMin;           // min eta
        Double_t        fEtaMax;           // max eta
        Int_t           fBit;              // filter bit
        Float_t         fPVzMax;           //max PVz
        Float_t         fPVzMin;           //min PVz
        Double_t        fChi2DoF;          // limit for chi2
        UShort_t        fTPCNCrossedRows;
        Bool_t          fIsRunFBOnly;      // only filterbit cuts
        Int_t           fTPCNcls;          // limit for TPC Ncls
        Bool_t          fIsPileUpCuts;     // pile up cuts flag
        Int_t           fPileUpLevel;
        TString         fGenName;          // MC generator name
    
    
        AliEventCuts    fEventCuts;

        AliAnalysisMultPt(const AliAnalysisMultPt&); // not implemented
        AliAnalysisMultPt& operator=(const AliAnalysisMultPt&); // not implemented

        ClassDef(AliAnalysisMultPt, 2);
};

#endif
