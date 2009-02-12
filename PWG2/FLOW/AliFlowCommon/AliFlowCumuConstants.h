/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#ifndef ALIFLOWCUMUCONSTANTS_H
#define ALIFLOWCUMUCONSTANTS_H

#include <TROOT.h>


// Description: constants for the LYZ flow makers


namespace AliFlowCumuConstants {


 // Enumerators
  enum {
    kQmax        = 11,
    kPmax        = 5,
    kQmax4        = 5,
    kPmax4        = 2,     
    kQmax6        = 7,
    kPmax6        = 3,   
    kQmax8        = 9,
    kPmax8        = 4, 
    kQmax16       = 17,
    kPmax16       = 8,     
             
    kFlow        = 2,  
    kMltpl       = 1
  };
 

 // Histograms limits
  extern Double_t  fgBinWidth;   
  extern Double_t  fgR0;   
  extern Double_t  fgPtMax;
  extern Double_t  fgPtMin;
  
 // Other numerical equations for cumulants 
  extern Bool_t  fgOtherEquations;
}

#endif

