//////////////////////////////////////////////////////////////////////
//
// $Id: AliFlowConstants.cxx 18618 2007-05-16 15:38:22Z snelling $
//
// Author: Emanuele Simili
//
//////////////////////////////////////////////////////////////////////
//_____________________________________________________________
//
// Description: constants for the flow makers and the flow analysis
//  bla bla bla bla bla bla bla bla bla bla bla bla bla bla bla  
//  bla bla bla bla bla bla bla bla bla bla bla bla bla bla ...
//
// Original Authors:                 Art Poskanzer & Raimond Snellings
//

#include "AliFlowConstants.h"

  
  Float_t  AliFlowConstants::fgEtaMin    = -2. ;
  Float_t  AliFlowConstants::fgEtaMax    =  2. ;
  Float_t  AliFlowConstants::fgPtMin     =  0. ;
  Float_t  AliFlowConstants::fgPtMax     = 10. ;
  Float_t  AliFlowConstants::fgPtMaxPart =  5. ;
  Float_t  AliFlowConstants::fgPtWgtSaturation =  5. ;
  Float_t  AliFlowConstants::fgEtaMinTpcOnly = -0.9 ;
  Float_t  AliFlowConstants::fgEtaMaxTpcOnly =  0.9 ;

  Float_t  AliFlowConstants::fgEtaMid = 0.5 ;
  Float_t  AliFlowConstants::fgEtaGood = 0.9 ;
  Float_t  AliFlowConstants::fgMaxMult = 2000. ;   // Maximum expected multiplicity (for fgCentNorm)
  Float_t  AliFlowConstants::fgCentNorm[kCents] = {0.0205,0.044,0.08,0.133,0.203,0.301,0.462,0.596,0.707} ;
  Int_t    AliFlowConstants::fgCent0[kCents]  = {50,100,200,500,1000,2000,5000,8000,10000} ;
  Double_t AliFlowConstants::fgBayesian[kPid] = {1.,1.,1.,1.,1.,1.} ;
 
  Double_t AliFlowConstants::fgMagneticField = 0.4 ;
  Double_t AliFlowConstants::fgCenterOfMassEnergy = 5500. ;
  Short_t  AliFlowConstants::fgBeamMassNumberEast = 208 ;
  Short_t  AliFlowConstants::fgBeamMassNumberWest = 208 ;

  Float_t  AliFlowConstants::fgITSx = 15. ;       // SDD where dE/dx is measured
  Float_t  AliFlowConstants::fgTPCx = 84.5 ;	  // TPC inner wall
  Float_t  AliFlowConstants::fgTRDx = 294.5 ;	  // TRD first plate
  Float_t  AliFlowConstants::fgTOFx = 370. ;	  // TOF 

  Int_t    AliFlowConstants::fgMClabel = 10000 ;  // mcLabel < 0 --> not used
  Bool_t   AliFlowConstants::fgDebug = kFALSE ;
  
  Float_t  AliFlowConstants::fgMaxInt = 1000. ;

//////////////////////////////////////////////////////////////////////
