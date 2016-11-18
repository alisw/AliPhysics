/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#ifndef ALIFLOWANALYSISWITHSIMPLESP_H
#define ALIFLOWANALYSISWITHSIMPLESP_H

class AliFlowTrackSimple;
class AliFlowEventSimple;
class AliFlowCommonHist;
class AliFlowCommonHistResults;

#include  "TList.h"
#include "AliFlowAnalysis.h"

class TH1D;
class TH1F;
class TH2D;
class TProfile;
class TDirectoryFile;

/////////////////////////////////////////////////////////////////////////////
// Description: Simple scalar product method
// author: Redmer Alexander Bertens (rbertens@cern.ch)
//         but cloned from the original Scalar Product task
/////////////////////////////////////////////////////////////////////////////

class AliFlowAnalysisWithSimpleSP : public AliFlowAnalysis {
    public:

        AliFlowAnalysisWithSimpleSP();            //default constructor
        virtual  ~AliFlowAnalysisWithSimpleSP();  //destructor

        void Init();                                       //Define output objects
        void Make(AliFlowEventSimple* anEvent);            //Main routine
        void GetOutputHistograms(TList *outputListHistos); //Copy output objects from TList
        void Finish(Bool_t A = kFALSE);                                     //Fill results
        void WriteHistograms(TDirectoryFile *outputFileName) const; //writes histograms locally (for OnTheFly)


        void SetHarmonic(Int_t iHarmonic)          { fHarmonic = iHarmonic; }
        void SetWeights(Bool_t weights)            { fWeights = weights; }
        void SetScaling(Bool_t scaling)            { fScaling = scaling; }
        void SetApplyCorrectionForNUA(Bool_t iVal) { fApplyCorrectionForNUA = iVal?1:0; }
        void SetNormalizationType(Int_t iVal)      { fNormalizationType = iVal; }
        void SetV0SanityCheck(Bool_t t)            { fV0SanityCheck = t;}
        void SetExternalResolution(Float_t r, Float_t s)      { 
            fExternalResolution = r;
            fExternalResErr = s;}
        Bool_t IsEqualRel(double dX, double dY, double epsilon) {            
            return TMath::Abs(dX - dY) <= epsilon; 
        } 

        void SetDebug(Bool_t bVal)                 { fDebug = bVal; }
        void SetBookOnlyBasicCCH(Bool_t bVal)           { fMinimalBook = bVal; }
        void SetTotalQvector(Int_t iVal)           { fTotalQvector = iVal; }

        void SetUsePhiWeights(Bool_t bVal)        { fUsePhiWeights = bVal; }
        void SetWeightsList(TList* const aWeightsList)  { fWeightsList = (TList*)aWeightsList->Clone(); }

        TList*    GetHistList()      const { return fHistList; }
        TProfile* GetHistProConfig() const { return fHistProConfig; }
        TH1D* GetHistProUQ(Int_t iRFPorPOI, Int_t iPTorETA) const { return fHistProUQ[iRFPorPOI][iPTorETA]; }
        TProfile* GetHistProQaQbNorm() const  { return fHistProQaQbNorm; }
        TProfile* GetHistProNUAq()     const  { return fHistProNUAq; }
        TProfile* GetHistProNUAu(Int_t iRFPorPOI, Int_t iPTorETA, Int_t iIMorRE) const { return fHistProNUAu[iRFPorPOI][iPTorETA][iIMorRE]; }
        TH1D*     GetHistSumOfWeights() const { return fHistSumOfWeights; }
        TProfile* GetHistProUQQaQb( Int_t iRFPorPOI, Int_t iPTorETA ) const { return fHistProUQQaQb[iRFPorPOI][iPTorETA]; }
        TH1D*     GetHistSumOfWeightsu(Int_t iRFPorPOI, Int_t iPTorETA, Int_t iWeight) const { return fHistSumOfWeightsu[iRFPorPOI][iPTorETA][iWeight]; }
        AliFlowCommonHist*        GetCommonHists()    const { return fCommonHists; }
        AliFlowCommonHistResults* GetCommonHistsRes() const { return fCommonHistsRes; }

    private:
        AliFlowAnalysisWithSimpleSP(const AliFlowAnalysisWithSimpleSP& anAnalysis);            //copy constructor
        AliFlowAnalysisWithSimpleSP& operator=(const AliFlowAnalysisWithSimpleSP& anAnalysis); //assignment operator


        Int_t fDebug ;                // flag for analysis: more print statements
        Bool_t fMinimalBook;          // flag to turn off QA and minimize FlowCommonHist
        Int_t fUsePhiWeights;         // use phi weights
        Int_t fApplyCorrectionForNUA; // apply correction for non-uniform acceptance
        Int_t fHarmonic;              // harmonic
        Bool_t fWeights;              // use evnet weights
        Bool_t fScaling;              // use qvec scaling 
        Int_t fNormalizationType;     // 0: EP mode || 1: SP mode
        Bool_t fV0SanityCheck;        // test to suppress NaN/inf values
        Float_t fExternalResolution;  // external resolution (only used when set)
        Float_t fExternalResErr;      // uncertainty of external resolution
        Int_t fTotalQvector;          // 1:Qa 2:Qb 3:QaQb

        TList*     fWeightsList;      // list holding input histograms with phi weights
        TList*     fHistList;         // list to hold all output histograms  
        TProfile*  fHistProConfig;    // configuration values
        TProfile*  fHistProQaQbNorm;  // average of QaQb
        TH1D*      fHistSumOfWeights; // holds sum of Na*Nb and (Na*Nb)^2
        TProfile*  fHistProNUAq;      // NUA related qq
        Float_t    fCentralityWeight; // centrality weight

        //QAHists
        TProfile* fHistProQNorm; // QNorm
        TH1D* fHistProQaQb; // QaQb
        TH1D* fHistProQaQbM; // QaQb/MaMb
        TH2D* fHistMaMb;         // MaMb
        TH2D* fHistQNormQaQbNorm; // QNorm vs QaQbNorm
        TH2D* fHistQaNormMa; // QaNorm Ma
        TH2D* fHistQbNormMb; // QbNorm Mb
        TH1D* fResolution; // Resolution
        TH1D* fHistQaQb; // QaQb
        TH1D* fHistQaQc; // QaQc
        TH1D* fHistQbQc; // QbQc
        TH1D* fHistQaQbCos; // QaQbCos
        TH1I* fHistNumberOfSubtractedDaughters; //how many daughters were subtracted during Q calculation?

        AliFlowCommonHist*        fCommonHists;    // control histograms
        AliFlowCommonHist*        fCommonHistsuQ;  // control histograms
        AliFlowCommonHistResults* fCommonHistsRes; // results histograms

        TH1F*      fPhiWeightsSub[2];           // histogram holding phi weights for subevents
        TProfile*  fHistProUQ[2][2];            // uQ for RP|POI PT|ETA
        TProfile*  fHistProUQQaQb[2][2];        // holds weighted average of <QuQaQb> for RP|POI PT|ETA
        TH1D*      fHistSumOfWeightsu[2][2][3]; // holds sums of 0: Nq', 1: Nq'^2, 2: Nq'*Na*Nb
        TProfile*  fHistProNUAu[2][2][2];          // NUA related qq for RP|POI PT|ETA

        ClassDef(AliFlowAnalysisWithSimpleSP,1)  // class version
};


#endif
