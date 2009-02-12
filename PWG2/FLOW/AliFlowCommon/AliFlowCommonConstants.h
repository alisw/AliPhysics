/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#ifndef ALIFLOWCOMMONCONSTANTS_H
#define ALIFLOWCOMMONCONSTANTS_H

#include <TROOT.h>

// AliFlowCommonConstants:
// Description: constants for the Common Histograms in the Flow Analysis
// Author: Naomi van der Kolk (kolk@nikhef.nl)

namespace AliFlowCommonConstants {
 
  // Enumerators
  enum {
    kNbinsMult = 10000,
    kNbinsPt   = 100,   
    kNbinsPhi  = 72,
    kNbinsEta  = 80,
    kNbinsQ    = 500
  };
 
  //getters
  Int_t GetNbinsMult(); 
  Int_t GetNbinsPt();   
  Int_t GetNbinsPhi();  
  Int_t GetNbinsEta();  
  Int_t GetNbinsQ();    

  // Histograms limits
  extern Double_t  fgMultMin ;            
  extern Double_t  fgMultMax ;
  extern Double_t  fgPtMin ;	     
  extern Double_t  fgPtMax ;
  extern Double_t  fgPhiMin ;	     
  extern Double_t  fgPhiMax ;
  extern Double_t  fgEtaMin ;	     
  extern Double_t  fgEtaMax ;	     
  extern Double_t  fgQMin ;	     
  extern Double_t  fgQMax ;
 
  //getters
  Double_t GetMultMin(); 
  Double_t GetMultMax(); 
  Double_t GetPtMin();   
  Double_t GetPtMax();   
  Double_t GetPhiMin();  
  Double_t GetPhiMax();  
  Double_t GetEtaMin();  
  Double_t GetEtaMax();  
  Double_t GetQMin();    
  Double_t GetQMax();    
  


  //ClassDef(AliFlowCommonConstants,0)  // macro for rootcint

}

#endif


