/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#ifndef ALIFLOWANALYSISWITHSCALARPRODUCT_H
#define ALIFLOWANALYSISWITHSCALARPRODUCT_H

#include "TString.h"

class AliFlowTrackSimple;
class AliFlowEventSimple;
class AliFlowCommonHist;
class AliFlowCommonHistResults;

class TProfile;
class TList;
class TFile;
class Riostream;


// Description: Maker to analyze Flow from the Scalar Product method.
//               
// author: N. van der Kolk (kolk@nikhef.nl)              

 
class AliFlowAnalysisWithScalarProduct {

 public:
 
   AliFlowAnalysisWithScalarProduct();            //default constructor
 
   virtual  ~AliFlowAnalysisWithScalarProduct();  //destructor
 
   void    Init();                                          //defines variables and histograms
   void    Make(AliFlowEventSimple* anEvent);               //calculates variables and fills histograms
   void    GetOutputHistograms(TList *outputListHistos);    //get pointers to all output histograms (called before Finish()) 
   void    Finish();                                        //saves histograms
   void    WriteHistograms(TString* outputFileName);        //writes histograms locally
   void    WriteHistograms(TString outputFileName);         //writes histograms locally
   void    WriteHistograms(TDirectoryFile *outputFileName); //writes histograms locally
   
   void    SetDebug(Bool_t kt)   { this->fDebug = kt ; }
   Bool_t  GetDebug() const      { return this->fDebug ; }

   //phi weights
   void    SetWeightsList(TList* const aWeightsList)  {this->fWeightsList = (TList*)aWeightsList->Clone();}
   TList*  GetWeightsList() const                     {return this->fWeightsList;}  
   void    SetUsePhiWeights(Bool_t const aPhiW)       {this->fUsePhiWeights = aPhiW;}
   Bool_t  GetUsePhiWeights() const                   {return this->fUsePhiWeights;}

   // Output 
   TList*    GetHistList() const    { return this->fHistList ; }     // Gets output histogram list
  
   TProfile* GetHistProUQetaRP() const {return this->fHistProUQetaRP;}   
   void      SetHistProUQetaRP(TProfile* const aHistProUQetaRP) {this->fHistProUQetaRP = aHistProUQetaRP;}
   TProfile* GetHistProUQetaPOI() const {return this->fHistProUQetaPOI;}
   void      SetHistProUQetaPOI(TProfile* const aHistProUQetaPOI) {this->fHistProUQetaPOI = aHistProUQetaPOI;}
   TProfile* GetHistProUQPtRP() const {return this->fHistProUQPtRP;} 
   void      SetHistProUQPtRP(TProfile* const aHistProUQPtRP) {this->fHistProUQPtRP = aHistProUQPtRP;}
   TProfile* GetHistProUQPtPOI() const {return this->fHistProUQPtPOI;}
   void      SetHistProUQPtPOI(TProfile* const aHistProUQPtPOI) {this->fHistProUQPtPOI = aHistProUQPtPOI;}
   TProfile* GetHistProQaQb() const {return this->fHistProQaQb;}
   void      SetHistProQaQb(TProfile* const aHistProQaQb) {this->fHistProQaQb = aHistProQaQb;}
   TProfile* GetHistProM() const {return this->fHistProM;}
   void      SetHistProM(TProfile* const aHistProM) {this->fHistProM = aHistProM;}
   TH1D*     GetHistM() const {return this->fHistM;}
   void      SetHistM(TH1D* const aHistM) {this->fHistM = aHistM;}

   AliFlowCommonHist* GetCommonHists() const {return this->fCommonHists; }
   void SetCommonHists(AliFlowCommonHist* const someCommonHists) {this->fCommonHists = someCommonHists; }
   AliFlowCommonHistResults* GetCommonHistsRes() const {return this->fCommonHistsRes; }
   void SetCommonHistsRes(AliFlowCommonHistResults* const someCommonHistsRes) {this->fCommonHistsRes = someCommonHistsRes; }
   


 private:
   AliFlowAnalysisWithScalarProduct(const AliFlowAnalysisWithScalarProduct& anAnalysis);            //copy constructor
   AliFlowAnalysisWithScalarProduct& operator=(const AliFlowAnalysisWithScalarProduct& anAnalysis); //assignment operator 
      
   Int_t      fEventNumber;      // event counter
   Bool_t     fDebug ;           // flag for analysis: more print statements

   TList*     fWeightsList;      // list holding input histograms with phi weights
   Bool_t     fUsePhiWeights;    // use phi weights
   TH1F*      fPhiWeights;       // histogram holding phi weights

   TList*     fHistList;         //list to hold all output histograms  
   TProfile*  fHistProUQetaRP;   //uQ(eta) for RP
   TProfile*  fHistProUQetaPOI;  //uQ(eta) for POI
   TProfile*  fHistProUQPtRP;    //uQ(pt) for RP
   TProfile*  fHistProUQPtPOI;   //uQ(pt) for POI
   TProfile*  fHistProQaQb;      //average of QaQb (for event plane resolution)
   TProfile*  fHistProM;         //holds avarage of M-1 and Ma*Mb
   TH1D*      fHistM;            //holds sum of M-1 and Ma*Mb
   AliFlowCommonHist*        fCommonHists;    //control histograms
   AliFlowCommonHistResults* fCommonHistsRes; //results histograms

   ClassDef(AliFlowAnalysisWithScalarProduct,0)  // macro for rootcint
     };
 

#endif
