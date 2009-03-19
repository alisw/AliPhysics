/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#ifndef ALIFLOWANALYSISWITHLEEYANGZEROS_H
#define ALIFLOWANALYSISWITHLEEYANGZEROS_H

////////////////////////////////////////////////////////////////////
// Description: Maker to analyze Flow by the LeeYangZeros method
//              One needs to do two runs over the data; 
//              First to calculate the integrated flow 
//              and in the second to calculate the differential flow
// Author: Naomi van der Kolk (kolk@nikhef.nl)
////////////////////////////////////////////////////////////////////


#include "AliFlowVector.h" //needed as include

//class AliFlowTrackSimple;
class AliFlowEventSimple;
class AliFlowLYZHist1; 
class AliFlowLYZHist2;
class AliFlowCommonHist;
class AliFlowCommonHistResults;

class TH1F;
class TH1D;
class TProfile;
class TProfile2D;
class TObjArray;
class TFile;
class TComplex;
class TString;
class TList;
class Riostream;

 
class AliFlowAnalysisWithLeeYangZeros {

 public:
 
   AliFlowAnalysisWithLeeYangZeros();                  //default constructor
   virtual  ~AliFlowAnalysisWithLeeYangZeros();        //destructor
 
   Bool_t    Init();                                   //defines variables and histograms
   Bool_t    Make(AliFlowEventSimple* anEvent);        //calculates variables and fills histograms
   Bool_t    Finish();                                 //saves histograms
   void      WriteHistograms(TString* outputFileName); //writes histograms locally
   void      WriteHistograms(TString outputFileName);  //writes histograms locally
   
   Double_t  GetQtheta(AliFlowVector aQ, Double_t aTheta);
   
   void      SetFirstRun(Bool_t kt)       { this->fFirstRun = kt ; }
   Bool_t    GetFirstRun() const          { return this->fFirstRun ; }
   void      SetUseSum(Bool_t kt)         { this->fUseSum = kt ; }
   Bool_t    GetUseSum() const            { return this->fUseSum ; }
   void      SetDoubleLoop(Bool_t kt)     { this->fDoubleLoop = kt ; }
   Bool_t    GetDoubleLoop() const        { return this->fDoubleLoop ; }
   void      SetDebug(Bool_t kt)          { this->fDebug = kt ; }
   Bool_t    GetDebug() const             { return this->fDebug ; }
   
   void      SetEventNumber(Int_t n)      { this->fEventNumber = n; }
   Int_t     GetEventNumber() const       { return this->fEventNumber; }
   void      SetQ2sum(Double_t d)         { this->fQ2sum = d; }
   Double_t  GetQ2sum() const             { return this->fQ2sum; }

   // Output 
   TList*             GetHistList() const      { return this->fHistList ; }     
   AliFlowCommonHist* GetCommonHists() const   { return this->fCommonHists; }
   void               SetCommonHists(AliFlowCommonHist* const aCommonHist)  
     { this->fCommonHists = aCommonHist; }
   AliFlowCommonHistResults* GetCommonHistsRes() const  
     { return this->fCommonHistsRes; }
   void               SetCommonHistsRes(AliFlowCommonHistResults* const aCommonHistResult) 
     { this->fCommonHistsRes = aCommonHistResult; }
   //AliFlowLYZHist1* GetHist1() const             {return this->fHist1; } 
   void               SetHist1(AliFlowLYZHist1* const aLYZHist1[])  
     {for (Int_t i=0;i<5;i++) {this->fHist1[i] = aLYZHist1[i];} }
   //AliFlowLYZHist2* GetHist2() const             {return this->fHist2; } 
   void               SetHist2RP(AliFlowLYZHist2* const aLYZHist2RP[])  
     {for (Int_t i=0;i<5;i++) {this->fHist2RP[i] = aLYZHist2RP[i];} }
   void               SetHist2POI(AliFlowLYZHist2* const aLYZHist2POI[])  
     {for (Int_t i=0;i<5;i++) {this->fHist2POI[i] = aLYZHist2POI[i];} }

   TProfile*  GetHistProVtheta() const   {return this->fHistProVtheta; } 
   void       SetHistProVtheta(TProfile* const aHistProVtheta)     
     { this->fHistProVtheta = aHistProVtheta; }
   TProfile*  GetHistProVetaRP() const     {return this->fHistProVetaRP; }  
   void       SetHistProVetaRP(TProfile* const aHistProVetaRP)         
     {this->fHistProVetaRP = aHistProVetaRP; }
   TProfile*  GetHistProVetaPOI() const     {return this->fHistProVetaPOI; }  
   void       SetHistProVetaPOI(TProfile* const aHistProVetaPOI)         
     {this->fHistProVetaPOI = aHistProVetaPOI; }
   TProfile*  GetHistProVPtRP() const      {return this->fHistProVPtRP;}
   void       SetHistProVPtRP(TProfile* const aHistProVPtRP)           
     {this->fHistProVPtRP = aHistProVPtRP; }
   TProfile*  GetHistProVPtPOI() const      {return this->fHistProVPtPOI;}
   void       SetHistProVPtPOI(TProfile* const aHistProVPtPOI)           
     {this->fHistProVPtPOI = aHistProVPtPOI; }
   TProfile*  GetHistProR0theta() const  {return this->fHistProR0theta; }
   void       SetHistProR0theta(TProfile* const aHistProR0theta)   
     {this->fHistProR0theta = aHistProR0theta; }
   TProfile*  GetHistProReDenom() const  {return this->fHistProReDenom; } 
   void       SetHistProReDenom(TProfile* const aHistProReDenom)   
     {this->fHistProReDenom = aHistProReDenom; }
   TProfile*  GetHistProImDenom() const  {return this->fHistProImDenom; }
   void       SetHistProImDenom(TProfile* const aHistProImDenom)   
     {this->fHistProImDenom = aHistProImDenom; }
   TProfile*  GetHistProReDtheta() const {return this->fHistProReDtheta; } 
   void       SetHistProReDtheta(TProfile* const aHistProReDtheta) 
     {this->fHistProReDtheta = aHistProReDtheta; }
   TProfile*  GetHistProImDtheta() const {return this->fHistProImDtheta; }
   void       SetHistProImDtheta(TProfile* const aHistProImDtheta) 
     {this->fHistProImDtheta = aHistProImDtheta; }
   TH1F*      GetHistQsumforChi() const  {return this->fHistQsumforChi; }
   void       SetHistQsumforChi(TH1F* const aHistQsumforChi) 
    {this->fHistQsumforChi = aHistQsumforChi; }

   void       SetFirstRunList(TList* const list) { this->fFirstRunList = list; }
   TList*     GetFirstRunList() const            { return this->fFirstRunList; }

 private:

   AliFlowAnalysisWithLeeYangZeros(const AliFlowAnalysisWithLeeYangZeros& aAnalysis);            // copy constructor
   AliFlowAnalysisWithLeeYangZeros& operator=(const AliFlowAnalysisWithLeeYangZeros& aAnalysis); //assignment operator

   Bool_t   MakeControlHistograms(AliFlowEventSimple* anEvent); 
   Bool_t   FillFromFlowEvent(AliFlowEventSimple* anEvent);
   Bool_t   SecondFillFromFlowEvent(AliFlowEventSimple* anEvent);

   TComplex GetGrtheta(AliFlowEventSimple* anEvent, Double_t aR, Double_t aTheta);
   TComplex GetDiffFlow(AliFlowEventSimple* anEvent, Double_t aR, Int_t theta); 
   Double_t GetR0(TH1D* fHistGtheta);
   
#ifndef __CINT__
   
   TVector2*  fQsum;         // flow vector sum              
   Double_t  fQ2sum;         // flow vector sum squared                 
      
#endif /*__CINT__*/

   Int_t        fEventNumber;     // event counter
        
   Bool_t       fFirstRun ;       // flag for lyz analysis: true=first run over data, false=second run 
   Bool_t       fUseSum ;         // flag for lyz analysis: true=use sum gen.function, false=use product gen.function
   Bool_t       fDoubleLoop ;     // flag for studying non flow effects
   Bool_t       fDebug ;          // flag for lyz analysis: more print statements

   TList*       fHistList;        //list to hold all output histograms 
   TList*       fFirstRunList;    //list from first run output
        
   TProfile*    fHistProVtheta;   //hist of V(theta)      
   TProfile*    fHistProVetaRP;   //hist of v(eta) for RP selection
   TProfile*    fHistProVetaPOI;  //hist of v(eta) for POI selection
   TProfile*    fHistProVPtRP;    //hist of v(pt) for RP selection
   TProfile*    fHistProVPtPOI;   //hist of v(pt) for POI selection
   TProfile*    fHistProR0theta;  //hist of r0(theta)    
   TProfile*    fHistProReDenom;  //hist of the real part of the denominator   
   TProfile*    fHistProImDenom;  //hist of the imaginary part of the denominator 
   TProfile*    fHistProReDtheta; //hist of the real part of D^theta   
   TProfile*    fHistProImDtheta; //hist of the imaginary part of D^theta  
   TH1F*        fHistQsumforChi;  //hist of sum of Q-vectors and the sum of the square of Q-vectors
  
    
  //class AliFlowLYZHist1 defines the histograms: fHistProGtheta, fHistProReGtheta, fHistProImGtheta
  AliFlowLYZHist1* fHist1[5];     //array of hist1

  //class AliFlowLYZHist1 defines the histograms: fHistProReNumer, fHistProImNumer, fHistProReNumerPt,
  //fHistProImNumerPt, fHistProReNumer2D, fHistProImNumer2D.
  AliFlowLYZHist2* fHist2RP[5];   //array of hist2
  AliFlowLYZHist2* fHist2POI[5];  //array of hist2

  AliFlowCommonHist*        fCommonHists;     //control histograms
  AliFlowCommonHistResults* fCommonHistsRes;  //final results histograms
 
  ClassDef(AliFlowAnalysisWithLeeYangZeros,1)  // macro for rootcint
    };
 
     
#endif

