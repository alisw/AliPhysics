/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#ifndef AliFlowAnalysisWithMCEventPlane_H
#define AliFlowAnalysisWithMCEventPlane_H

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
class TComplex;
class Riostream;


// Description: Maker to analyze Flow from the generated MC reaction plane.
//              This class is used to get the real value of the flow 
//              to compare the other methods to when analysing simulated events.

 
class AliFlowAnalysisWithMCEventPlane {

 public:
 
   AliFlowAnalysisWithMCEventPlane();   //default constructor
 
   virtual  ~AliFlowAnalysisWithMCEventPlane();  //destructor
 
   void    Init();                             //defines variables and histograms
   void    Make(AliFlowEventSimple* anEvent, Double_t fRP);   //calculates variables and fills histograms
   void    Finish();                           //saves histograms
  
   void      SetDebug(Bool_t kt)                 { this->fDebug = kt ; }
   Bool_t    GetDebug() const                    { return this->fDebug ; }


   // Output 
   void	    SetHistFileName(TString name) 	{ this->fHistFileName = name ; } // Sets output file name
   TString  GetHistFileName() const		{ return this->fHistFileName ; } // Gets output file name
   TFile*   GetHistFile() const                 { return this->fHistFile ; }     // Gets output file
  
   
 private:
 
   AliFlowAnalysisWithMCEventPlane(const AliFlowAnalysisWithMCEventPlane& aAnalysis);
   AliFlowAnalysisWithMCEventPlane& operator=(const AliFlowAnalysisWithMCEventPlane& aAnalysis);

      
#ifndef __CINT__
   AliFlowVector  fQ;       // flow vector
   TVector2  fQsum;         // flow vector sum
   Double_t  fQ2sum;        // flow vector sum squared
#endif /*__CINT__*/

   Int_t     fEventNumber;  // event counter
   Int_t     fMult;         // multiplicity
   Int_t     fNbins;        // number of bins
      
   AliFlowTrackSimple*   fTrack ;     //!
     
   Bool_t       fDebug ;            //! flag for lyz analysis: more print statements

   TString      fHistFileName;      //!
   TFile*       fHistFile;          //!
    
   AliFlowCommonHist* fCommonHists;              //!
   AliFlowCommonHistResults* fCommonHistsRes;    //!

   TProfile*  fHistProFlow;         //!
   TH1F*      fHistRP;              //!

   ClassDef(AliFlowAnalysisWithMCEventPlane,0)  // macro for rootcint
     };
 
     
#endif


