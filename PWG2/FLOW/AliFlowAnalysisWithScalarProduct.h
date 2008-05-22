/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#ifndef AliFlowAnalysisWithScalarProduct_H
#define AliFlowAnalysisWithScalarProduct_H

#include "TVector2.h"          //called explicitly
#include "AliFlowVector.h"
#include "TString.h"

class AliFlowTrackSimple;
class AliFlowEventSimple;
class AliFlowCommonHist;
//class AliFlowCommonHistResults;

class TProfile;
class TFile;
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
  
   void      SetDebug(Bool_t kt)                 { this->fDebug = kt ; }
   Bool_t    GetDebug() const                    { return this->fDebug ; }


   // Output 
   void	    SetHistFileName(TString name) 	{ this->fHistFileName = name ; } // Sets output file name
   TString  GetHistFileName() const		{ return this->fHistFileName ; } // Gets output file name
   TFile*   GetHistFile() const                 { return this->fHistFile ; }     // Gets output file
  
   
 private:
   //cp const
   //ass op
   
   AliFlowVector  fQ;       // flow vector
   TVector2  fU;            // particle unit vector
   
   Int_t     fEventNumber;  // event counter
         
   AliFlowTrackSimple*   fTrack ;     //!
     
   Bool_t       fDebug ;            //! flag for lyz analysis: more print statements

   TString      fHistFileName;      //!
   TFile*       fHistFile;          //!
      
   TProfile*      fHistProUQ;              //!

   AliFlowCommonHist* fCommonHists;              //!
   //AliFlowCommonHistResults* fCommonHistsRes;    //!

   ClassDef(AliFlowAnalysisWithScalarProduct,0)  // macro for rootcint
     };
 
     
#endif
