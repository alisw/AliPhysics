/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//_________________________________________________________________________
// Reconstructed Particle
//*-- Y. Schutz:   SUBATECH 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSRecParticle.h"

ClassImp(AliPHOSRecParticle)


//____________________________________________________________________________
 AliPHOSRecParticle::AliPHOSRecParticle(AliPHOSTrackSegment * ts)
{
  // ctor
 
  fPHOSTrackSegment = new AliPHOSTrackSegment(*ts) ; 
  fE                = ts->GetEnergy() ; 
  TVector3 momdir   = ts->GetMomentumDirection() ;
  fPx               = fE * momdir.X() ; 
  fPy               = fE * momdir.Y() ; 
  fPz               = fE * momdir.Z() ; 

  fType = kUNDEFINED ; // undefined
                            
  if( ts->GetPpsdUp() == 0 ) {     // Neutral

    if( ts->GetPpsdLow() == 0 )    // Neutral  
      fType = kNEUTRAL ;   
    else                           // Gamma
      fType = kGAMMA ;               
  }
  else                            // Charged           
    fType = kCHARGED ;   
  
}

//____________________________________________________________________________
TString AliPHOSRecParticle::Name()
{
  TString  name ; 
  switch (fType) {
  case kGAMMA:
    name = "PHOTON" ;
    break ; 
   case kELECTRON:
     name = "ELECTRON" ;
    break ; 
  case kNEUTRAL:
    name = "NEUTRAL" ;
    break ; 
   case kCHARGEDHADRON:
    name = "CHARGED HADRON" ;
    break ; 
  case kNEUTRALHADRON:
    name = "NEUTRAL HADRON" ; 
    break ; 
  case kNEUTRALEM:
    name = "NEUTRAL EM" ; 
    break ; 
  }
  return name ; 
}

//____________________________________________________________________________
void AliPHOSRecParticle::Print()
{
  cout << "AliPHOSRecParticle > " << "type is  " << Name() << endl 
       << "                     " << "Energy = " << fE << endl 
       << "                     " << "Px     = " << fPx << endl 
       << "                     " << "Py     = " << fPy << endl 
       << "                     " << "Pz     = " << fPz << endl ; 
}
