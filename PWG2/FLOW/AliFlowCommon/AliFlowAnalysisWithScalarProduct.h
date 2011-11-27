/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#ifndef ALIFLOWANALYSISWITHSCALARPRODUCT_H
#define ALIFLOWANALYSISWITHSCALARPRODUCT_H

class AliFlowTrackSimple;
class AliFlowEventSimple;
class AliFlowCommonHist;
class AliFlowCommonHistResults;

class TH1D;
class TH2D;
class TProfile;
#include  "TList.h"
class TDirectoryFile;

/////////////////////////////////////////////////////////////////////////////
// Description: Maker to analyze Flow from the Event Plane method.
//              Adaptation based on Scalar Product
// authors: Naomi van del Kolk (kolk@nikhef.nl)
//          Ante Bilandzic (anteb@nikhef.nl)
// mods:    Carlos Perez (cperez@nikhef.nl)
/////////////////////////////////////////////////////////////////////////////
 
class AliFlowAnalysisWithScalarProduct {
 public:
 
   AliFlowAnalysisWithScalarProduct();            //default constructor
   virtual  ~AliFlowAnalysisWithScalarProduct();  //destructor
 
   void Init();                                       //Define output objects
   void Make(AliFlowEventSimple* anEvent);            //Main routine
   void GetOutputHistograms(TList *outputListHistos); //Copy output objects from TList
   void Finish();                                     //Fill results
   void WriteHistograms(TDirectoryFile *outputFileName); //writes histograms locally (for OnTheFly)


   void SetHarmonic(Int_t iHarmonic)          { fHarmonic = iHarmonic; }
   void SetApplyCorrectionForNUA(Bool_t iVal) { fApplyCorrectionForNUA = iVal?1:0; }
   void SetNormalizationType(Int_t iVal)      { fNormalizationType = iVal; }
   void SetDebug(Bool_t bVal)                 { fDebug = bVal; }
   void SetTotalQvector(Int_t iVal)           { fTotalQvector = iVal; }

   void SetUsePhiWeights(Bool_t bVal)        { fUsePhiWeights = bVal; }
   void SetWeightsList(TList* const aWeightsList)  { fWeightsList = (TList*)aWeightsList->Clone(); }
  
   TList*    GetHistList()      { return fHistList; }
   TProfile* GetHistProConfig() { return fHistProConfig; }
   TProfile* GetHistProUQ(Int_t iRFPorPOI, Int_t iPTorETA) { return fHistProUQ[iRFPorPOI][iPTorETA]; }
   TProfile* GetHistProQaQbNorm()   { return fHistProQaQbNorm; }
   TProfile* GetHistProNUAq()   { return fHistProNUAq; }
   TProfile* GetHistProNUAu(Int_t iRFPorPOI, Int_t iPTorETA, Int_t iIMorRE) { return fHistProNUAu[iRFPorPOI][iPTorETA][iIMorRE]; }
   TH1D*     GetHistSumOfWeights() { return fHistSumOfWeights; }
   TProfile* GetHistProUQQaQb( Int_t iRFPorPOI, Int_t iPTorETA ) { return fHistProUQQaQb[iRFPorPOI][iPTorETA]; }
   TH1D*     GetHistSumOfWeightsu(Int_t iRFPorPOI, Int_t iPTorETA, Int_t iWeight) { return fHistSumOfWeightsu[iRFPorPOI][iPTorETA][iWeight]; }
   AliFlowCommonHist*        GetCommonHists()    { return fCommonHists; }
   AliFlowCommonHistResults* GetCommonHistsRes() { return fCommonHistsRes; }
   
 private:
   AliFlowAnalysisWithScalarProduct(const AliFlowAnalysisWithScalarProduct& anAnalysis);            //copy constructor
   AliFlowAnalysisWithScalarProduct& operator=(const AliFlowAnalysisWithScalarProduct& anAnalysis); //assignment operator
   Double_t CalculateStatisticalError( Int_t RFPorPOI, Int_t PTorETA, Int_t bin, Double_t errV );
   Double_t computeResolution( Double_t x );
   Double_t findXi( Double_t res, Double_t prec );

      
   Int_t fDebug ;                // flag for analysis: more print statements
   Int_t fUsePhiWeights;         // use phi weights
   Int_t fApplyCorrectionForNUA; // apply correction for non-uniform acceptance
   Int_t fHarmonic;              // harmonic 
   Int_t fNormalizationType;     // 0: EP mode || 1: SP mode
   Int_t fTotalQvector;          // 1:Qa 2:Qb 3:QaQb

   TList*     fWeightsList;      // list holding input histograms with phi weights
   TList*     fHistList;         // list to hold all output histograms  
   TProfile*  fHistProConfig;    // configuration values
   TProfile*  fHistProQaQbNorm;  // average of QaQb
   TH1D*      fHistSumOfWeights; // holds sum of Na*Nb and (Na*Nb)^2
   TProfile*  fHistProNUAq;      // NUA related qq

   //QAHists
   TProfile* fHistProQNorm;
   TProfile* fHistProQaQb;
   TProfile* fHistProQaQbM;
   TH2D* fHistMaMb;
   TH2D* fHistQNormQaQbNorm;
   TH2D* fHistQaNormMa;
   TH2D* fHistQbNormMb;
   TH1D* fResolution;
   TH1D* fHistQaQb;
   TH1D* fHistQaQbCos;

   AliFlowCommonHist*        fCommonHists;    // control histograms
   AliFlowCommonHist*        fCommonHistsuQ;  // control histograms
   AliFlowCommonHistResults* fCommonHistsRes; // results histograms

   TH1F*      fPhiWeightsSub[2];           // histogram holding phi weights for subevents
   TProfile*  fHistProUQ[2][2];            // uQ for RP|POI PT|ETA
   TProfile*  fHistProUQQaQb[2][2];        // holds weighted average of <QuQaQb> for RP|POI PT|ETA
   TH1D*      fHistSumOfWeightsu[2][2][3]; // holds sums of 0: Nq', 1: Nq'^2, 2: Nq'*Na*Nb
   TProfile*  fHistProNUAu[2][2][2];          // NUA related qq for RP|POI PT|ETA

   ClassDef(AliFlowAnalysisWithScalarProduct,0)  // class version
};
 

#endif
