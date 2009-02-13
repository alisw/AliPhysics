/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#ifndef ALIFLOWANALYSISWITHMCEVENTPLANE_H
#define ALIFLOWANALYSISWITHMCEVENTPLANE_H

#include "TVector2.h"          //called explicitly
#include "AliFlowVector.h"
#include "TString.h"

class AliFlowTrackSimple;
class AliFlowEventSimple;
class AliFlowCommonHist;
class AliFlowCommonHistResults;

class TH1F;
class TH1D;
class TProfile;
class TProfile2D;
class TObjArray;
class TFile;
class TList;
class TComplex;
class Riostream;


// Description: Maker to analyze Flow from the generated MC reaction plane.
//              This class is used to get the real value of the flow 
//              to compare the other methods to when analysing simulated events.

 
class AliFlowAnalysisWithMCEventPlane {

 public:
 
   AliFlowAnalysisWithMCEventPlane();            //default constructor
   virtual  ~AliFlowAnalysisWithMCEventPlane();  //destructor
 
   void      WriteHistograms(TString* outputFileName);
   void      WriteHistograms(TString outputFileName);
   void      Init();                             //defines variables and histograms
   void      Make(AliFlowEventSimple* anEvent, Double_t aRP);   //calculates variables and fills histograms
   void      Finish();                           //saves histograms
   
   void      SetDebug(Bool_t kt)          { this->fDebug = kt ; }
   Bool_t    GetDebug() const             { return this->fDebug ; }

   void      SetEventNumber(Int_t n)      { this->fEventNumber = n; }
   Int_t     GetEventNumber() const       { return this->fEventNumber; }

   // Output 
   TList*    GetHistList() const          { return this->fHistList ; }  
   AliFlowCommonHist* GetCommonHists() const  { return this->fCommonHists; }
   void               SetCommonHists(AliFlowCommonHist* aCommonHist)  
     { this->fCommonHists = aCommonHist; }
   AliFlowCommonHistResults*  GetCommonHistsRes() const  
     { return this->fCommonHistsRes; }
   void      SetCommonHistsRes(AliFlowCommonHistResults* aCommonHistResult) 
     { this->fCommonHistsRes = aCommonHistResult; }
   
   //histograms
   TProfile* GetHistProFlow()             {return this->fHistProFlow; }  
   void      SetHistProFlow(TProfile* aHistProFlow)   
     {this->fHistProFlow = aHistProFlow; }
   TH1F*     GetHistRP()                  {return this->fHistRP; } 
   void      SetHistRP(TH1F* aHistRP)     {this->fHistRP = aHistRP; }
   
   TProfile* GetHistProIntFlow()                          {return this->fHistProIntFlow; } 
   void      SetHistProIntFlow(TProfile* aHistProIntFlow) {this->fHistProIntFlow = aHistProIntFlow; }
   
   TProfile* GetHistProDiffFlowPtRP()                               {return this->fHistProDiffFlowPtRP; } 
   void      SetHistProDiffFlowPtRP(TProfile* aHistProDiffFlowPtRP) {this->fHistProDiffFlowPtRP = aHistProDiffFlowPtRP; } 
   
   TProfile* GetHistProDiffFlowEtaRP()                                {return this->fHistProDiffFlowEtaRP; } 
   void      SetHistProDiffFlowEtaRP(TProfile* aHistProDiffFlowEtaRP) {this->fHistProDiffFlowEtaRP = aHistProDiffFlowEtaRP; } 
   
   TProfile* GetHistProDiffFlowPtPOI()                               {return this->fHistProDiffFlowPtPOI; } 
   void      SetHistProDiffFlowPtPOI(TProfile* aHistProDiffFlowPtPOI) {this->fHistProDiffFlowPtPOI = aHistProDiffFlowPtPOI; } 
   
   TProfile* GetHistProDiffFlowEtaPOI()                                {return this->fHistProDiffFlowEtaPOI; } 
   void      SetHistProDiffFlowEtaPOI(TProfile* aHistProDiffFlowEtaPOI) {this->fHistProDiffFlowEtaPOI = aHistProDiffFlowEtaPOI; }    

 private:
 
   AliFlowAnalysisWithMCEventPlane(const AliFlowAnalysisWithMCEventPlane& aAnalysis);             //copy constructor
   AliFlowAnalysisWithMCEventPlane& operator=(const AliFlowAnalysisWithMCEventPlane& aAnalysis);  //assignment operator 

      
#ifndef __CINT__
   TVector2*    fQsum;              // flow vector sum
   Double_t     fQ2sum;             // flow vector sum squared
#endif /*__CINT__*/

   Int_t        fEventNumber;       // event counter
   Bool_t       fDebug ;            //! flag for lyz analysis: more print statements

   TList*       fHistList;          //list to hold all output histograms  
    
   AliFlowCommonHist* fCommonHists;              //
   AliFlowCommonHistResults* fCommonHistsRes;    //

   TProfile*    fHistProFlow;       //
   TH1F*        fHistRP;            //
   
   TProfile*    fHistProIntFlow;     //profile used to calculate the integrated flow of RP particles
   TProfile*    fHistProDiffFlowPtRP;  //profile used to calculate the differential flow (Pt) of RP particles 
   TProfile*    fHistProDiffFlowEtaRP; //profile used to calculate the differential flow (Eta) of RP particles 
   
   TProfile*    fHistProDiffFlowPtPOI;  //profile used to calculate the differential flow (Pt) of POI particles 
   TProfile*    fHistProDiffFlowEtaPOI; //profile used to calculate the differential flow (Eta) of POI particles    
   
   ClassDef(AliFlowAnalysisWithMCEventPlane,1)  // macro for rootcint
     };
 
     
#endif


