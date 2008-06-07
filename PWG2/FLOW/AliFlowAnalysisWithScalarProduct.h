/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#ifndef AliFlowAnalysisWithScalarProduct_H
#define AliFlowAnalysisWithScalarProduct_H

#include "TVector2.h"          //called explicitly
#include "TString.h"

#include "AliFlowVector.h"

class AliFlowTrackSimple;
class AliFlowEventSimple;
class AliFlowCommonHist;
//class AliFlowCommonHistResults;

class TProfile;
class TList;
class Riostream;


// Description: Maker to analyze Flow from the Scalar Product method.
//               
// author: N. van der Kolk (kolk@nikhef.nl)              

 
class AliFlowAnalysisWithScalarProduct {

 public:
 
   AliFlowAnalysisWithScalarProduct();   //default constructor
 
   virtual  ~AliFlowAnalysisWithScalarProduct();  //destructor
 
   void    Init();                             //defines variables and histograms
   void    Make(AliFlowEventSimple* anEvent);   //calculates variables and fills histograms
   void    Finish();                           //saves histograms
  
   void      SetDebug(Bool_t kt)            { this->fDebug = kt ; }
   Bool_t    GetDebug() const               { return this->fDebug ; }


   // Output 
   TList*   GetHistList() const           { return this->fHistList ; }     // Gets output histogram list
  
   
 private:
  AliFlowAnalysisWithScalarProduct(const AliFlowAnalysisWithScalarProduct& aAnalysis);
  AliFlowAnalysisWithScalarProduct& operator=(const AliFlowAnalysisWithScalarProduct& aAnalysis); 
   
  AliFlowVector  *fQ;       // flow vector
  TVector2       *fU;       // particle unit vector
   
  Int_t          fEventNumber;  // event counter
         
  Bool_t         fDebug ;            // flag for lyz analysis: more print statements

  TList*         fHistList;      
  TProfile*      fHistProUQ; 

  AliFlowCommonHist* fCommonHists;  

  ClassDef(AliFlowAnalysisWithScalarProduct,0)  // macro for rootcint
    };
 

#endif
