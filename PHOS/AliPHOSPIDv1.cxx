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

/* $Id$ */

//_________________________________________________________________________
// Implementation version v1 of the PHOS particle identifier 
// Identification is based on information from PPSD and EMC
//                  
//*-- Author: Yves Schutz (SUBATECH)  & Gines Martinez (SUBATECH)


// --- ROOT system ---

// --- Standard library ---

#include <iostream.h>

// --- AliRoot header files ---

#include "AliPHOSPIDv1.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSRecParticle.h"
#include "AliPHOSIndexToObject.h"

ClassImp( AliPHOSPIDv1) 

//____________________________________________________________________________
AliPHOSPIDv1::AliPHOSPIDv1():AliPHOSPID()
{ 
  fCutOnDispersion = 2.0; 
  fCutOnRelativeDistance = 3.0 ;
}

//____________________________________________________________________________
Float_t  AliPHOSPIDv1::GetDistanceInPHOSPlane(AliPHOSEmcRecPoint * emcclu,AliPHOSPpsdRecPoint * PpsdClu, Bool_t &toofar, Option_t *  Axis)
{
  // Calculates the distance between the EMC RecPoint and the PPSD RecPoint
 
  Float_t r = 1.e+10;
  TVector3 vecEmc ;
  TVector3 vecPpsd ;
  
  emcclu->GetLocalPosition(vecEmc) ;
  PpsdClu->GetLocalPosition(vecPpsd)  ; 
  if(emcclu->GetPHOSMod() == PpsdClu->GetPHOSMod())
    { 

      // Correct to difference in CPV and EMC position due to different distance to center.
      // we assume, that particle moves from center
      Float_t dCPV = fGeom->GetIPtoOuterCoverDistance();
      Float_t dEMC = fGeom->GetIPtoCrystalSurface() ;
      dEMC         = dEMC / dCPV ;
      vecPpsd = dEMC * vecPpsd  - vecEmc ; 
      r = vecPpsd.Mag() ;
      if (Axis == "X") r = vecPpsd.X();
      if (Axis == "Y") r = vecPpsd.Y();
      if (Axis == "Z") r = vecPpsd.Z();
      if (Axis == "R") r = vecPpsd.Mag();
      
    } 
  else 
    {
      toofar = kTRUE ;
    }
  return r ;
}

//____________________________________________________________________________
void  AliPHOSPIDv1::MakeParticles(AliPHOSTrackSegment::TrackSegmentsList * trsl, 
				  AliPHOSRecParticle::RecParticlesList * rpl)
{
  // Makes a RecParticle out of a TrackSegment

  TIter next(trsl) ; 
  AliPHOSTrackSegment * tracksegment ; 
  Int_t index = 0 ; 
  AliPHOSRecParticle * rp ; 
  Bool_t tDistance;
  Int_t type ; 
  Int_t showerprofile;  // 0 narrow and 1 wide
  Int_t cpvdetector ;   // 1 hit and 0 no hit
  Int_t pcdetector ;    // 1 hit and 0 no hit

  while ( (tracksegment = (AliPHOSTrackSegment *)next()) ) {
    Int_t module = tracksegment->GetPHOSMod();
    cout << "PHOS module: " << module << endl;
    if ( module <= fGeom->GetNCPVModules()) continue;
    new( (*rpl)[index] ) AliPHOSRecParticle(tracksegment) ;
    rp = (AliPHOSRecParticle *)rpl->At(index) ; 
    AliPHOSEmcRecPoint * recp = tracksegment->GetEmcRecPoint() ;
    AliPHOSPpsdRecPoint * rpcpv = tracksegment->GetPpsdUpRecPoint() ;
    AliPHOSPpsdRecPoint * rppc  = tracksegment->GetPpsdLowRecPoint() ;
 
//     Float_t * lambda = new Float_t[2]; 
//     recp->GetElipsAxis(lambda) ; 

//     // Looking at the lateral development of the shower
//     if ( ( lambda[0] > fLambda1m && lambda[0] < fLambda1M ) && // shower profile cut
// 	 ( lambda[1] > fLambda2m && lambda[1] < fLambda2M ) )	      
//       //    Float_t R ;
//       //R=(lambda[0]-1.386)*(lambda[0]-1.386)+1.707*1.707*(lambda[1]-1.008)*(lambda[1]-1.008) ;
//       //if(R<0.35*0.35)

    Float_t Dispersion;
    Dispersion = recp->GetDispersion();
    if (Dispersion < fCutOnDispersion)
      showerprofile = 0 ;   // NARROW PROFILE   
    else      
      showerprofile = 1 ;// WIDE PROFILE
  

    // Looking at the photon conversion detector
    if( tracksegment->GetPpsdLowRecPoint() == 0 )   
      pcdetector = 0 ;  // No hit
    else{      
      if (GetDistanceInPHOSPlane(recp, rppc, tDistance, "R")  < fCutOnRelativeDistance) 
	pcdetector = 1 ;  // hit
      else
	pcdetector = 0 ;
    }
  
    // Looking at the photon conversion detector
    if( tracksegment->GetPpsdUpRecPoint() == 0 )
      cpvdetector = 0 ;  // No hit
    else{  
      if (GetDistanceInPHOSPlane(recp, rpcpv, tDistance, "R")< fCutOnRelativeDistance) 
	cpvdetector = 1 ;  // Hit
      else
	cpvdetector = 0 ;
    }
     
    type = showerprofile + 2 * pcdetector + 4 * cpvdetector ;
    rp->SetType(type) ; 
    index++ ; 
  }
    
}

//____________________________________________________________________________
void  AliPHOSPIDv1:: Print(const char * opt) 
{
  // Print the parameters used for the particle type identification
  
  cout << "AliPHOSPIDv1 : cuts for the particle idendification based on the shower profile " << endl 
       << fLambda1m << " < value1 < " << fLambda1M << endl 
       << fLambda2m << " < value2 < " << fLambda2M << endl ;  

}

//____________________________________________________________________________
void  AliPHOSPIDv1::SetShowerProfileCuts(Float_t l1m, Float_t l1M, Float_t l2m, Float_t l2M)
{
  // Modifies the parameters used for the particle type identification

  fLambda1m = l1m ; 
  fLambda1M = l1M ; 
  fLambda2m = l2m ; 
  fLambda2M = l2M ; 
}

//____________________________________________________________________________
void  AliPHOSPIDv1::SetRelativeDistanceCut(Float_t CutOnRelativeDistance)
{
  // Modifies the parameters used for the particle type identification

  fCutOnRelativeDistance = CutOnRelativeDistance ; 
 
}

