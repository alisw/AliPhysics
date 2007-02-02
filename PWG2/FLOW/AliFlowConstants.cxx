//////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//_____________________________________________________________
//
// Description: constants for the flow makers and the flow analysis
//
// Original Authors:                 Art Poskanzer & Raimond Snellings
//

#include "AliFlowConstants.h"

ClassImp(Flow)
  
  Float_t  Flow::fEtaMin    = -2. ;
  Float_t  Flow::fEtaMax    =  2. ;
  Float_t  Flow::fPtMin     =  0. ;
  Float_t  Flow::fPtMax     = 10. ;
  Float_t  Flow::fPtMaxPart =  5. ;
  Float_t  Flow::fPtWgtSaturation =  5. ;
  Float_t  Flow::fEtaMinTpcOnly = -0.9 ;
  Float_t  Flow::fEtaMaxTpcOnly =  0.9 ;

  Float_t  Flow::fEtaMid = 0.5 ;
  Float_t  Flow::fEtaGood = 0.9 ;
  Float_t  Flow::fMaxMult = 2900. ; // Maximum expected multiplicity (now hijing)
  Float_t  Flow::fCentNorm[nCents] = {0.0205,0.044,0.08,0.133,0.203,0.301,0.462,0.596,0.707} ;
  Int_t    Flow::fCent0[nCents]  = {50,100,200,500,1000,2000,5000,8000,10000} ;
  Double_t Flow::fBayesian[nPid] = {1.,1.,1.,1.,1.,1.} ;
 
  Double_t Flow::fMagneticField = 0.4 ;
  Double_t Flow::fCenterOfMassEnergy = 5500. ;
  Short_t  Flow::fBeamMassNumberEast = 208 ;
  Short_t  Flow::fBeamMassNumberWest = 208 ;

  Float_t  Flow::fITSx = 15. ;     	// SDD (?) where dE/dx is measured
  Float_t  Flow::fTPCx = 84.5 ;		// TPC (?) inner wall
  Float_t  Flow::fTRDx = 294.5 ;	// TRD first plate
  Float_t  Flow::fTOFx = 370. ;		// TOF (?)

  Int_t    Flow::fMClabel = 10000 ; 		  // mcLabel < 0 --> not used
  Bool_t   Flow::fDebug = kFALSE ;
  
  Float_t  Flow::fMaxInt = 1000. ;

//////////////////////////////////////////////////////////////////////
