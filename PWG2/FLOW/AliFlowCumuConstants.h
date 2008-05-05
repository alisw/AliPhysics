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
    kQmax        = 17,
    kPmax        = 8,  
    kFlow        = 2,  
    kMltpl       = 1
  };
 

 // Histograms limits
  extern Double_t  fgBinWidth;   
  extern Double_t  fgR0;   
  extern Double_t  fgPtMax ;
  extern Double_t  fgPtMin ;
}

#endif

