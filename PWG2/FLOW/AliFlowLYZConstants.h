/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#ifndef ALIFLOWLYZCONSTANTS_H
#define ALIFLOWLYZCONSTANTS_H

#include <TROOT.h>


// Description: constants for the LYZ flow makers


namespace AliFlowLYZConstants {


 // Enumerators
  enum {
    kTheta       = 5,                // number of reference angles theta
    kNbins       = 500,              // number of bins in fHistGtheta (AliFlowLYZHist1)
    kEtaBins     = 100,		     // number of eta bins in histograms (AliFlowLYZHist2)
    kPtBins      = 100		     // number of pT bins in histograms (AliFlowLYZHist2)
  };
 

 // Histograms limits
  extern Float_t  fgMin ;            // lower limit for fHistGtheta (AliFlowLYZHist1)
  extern Float_t  fgMax ;            // upper limit for fHistGtheta (AliFlowLYZHist1)
  extern Float_t  fgEtaMin ;	     // eta lower limit for histograms
  extern Float_t  fgEtaMax ;	     // eta upper limit for histograms
  extern Float_t  fgPtMin ;	     // pT lower limit for  histograms
  extern Float_t  fgPtMax ;	     // pT upper limit for  histograms
  
 
}

#endif

