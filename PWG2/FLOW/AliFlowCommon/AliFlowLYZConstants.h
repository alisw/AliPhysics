/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliFlowLYZConstants.h 25556 2008-05-02 08:22:51Z snelling $ */

#ifndef ALIFLOWLYZCONSTANTS_H
#define ALIFLOWLYZCONSTANTS_H

#include <TROOT.h>


// Description: constants for the LYZ flow makers


namespace AliFlowLYZConstants {


 // Enumerators
  enum {
    kTheta       = 5,     // number of reference angles theta
    kNbins       = 8000   // number of bins in fHistGtheta (AliFlowLYZHist1)
  };
 

 // Histograms limits
  extern Double_t  fgMin ;   // lower limit for fHistGtheta (AliFlowLYZHist1)
  extern Double_t  fgMax ;   // upper limit for fHistGtheta (AliFlowLYZHist1)
   
 
}

#endif

