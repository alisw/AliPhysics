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
//  A Reconstructed Particle in PHOS    
//  To become a general class of AliRoot ?        
//       
//*-- Author: Yves Schutz (SUBATECH)


// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSRecParticle.h"
#include "TPad.h"
#include "AliPHOSIndexToObject.h"

ClassImp(AliPHOSRecParticle)


//____________________________________________________________________________
 AliPHOSRecParticle::AliPHOSRecParticle(AliPHOSTrackSegment * ts)
{
  // ctor
 
  fPHOSTrackSegment = ts->GetIndexInList() ;
  fIndexInList      = -1 ;
  Float_t kenergy   = ts->GetEnergy() ; 
  TVector3 momdir   = ts->GetMomentumDirection() ;
  fPx               = kenergy * momdir.X() ; 
  fPy               = kenergy * momdir.Y() ; 
  fPz               = kenergy * momdir.Z() ; 
  fType             = kUNDEFINED;  
  fE                = kenergy ;    // !!! all particles have mass = 0 
}

//____________________________________________________________________________
 AliPHOSRecParticle::AliPHOSRecParticle(const AliPHOSRecParticle & rp)
{
  // copy ctor

  fPHOSTrackSegment = rp.fPHOSTrackSegment ; 
  fType             = rp.fType ; 
  fIndexInList      = rp.fIndexInList ;

  fPdgCode     = rp.fPdgCode;
  fStatusCode  = rp.fStatusCode;
  fMother[0]   = rp.fMother[0];
  fMother[1]   = rp.fMother[1];
  fDaughter[0] = rp.fDaughter[0];
  fDaughter[1] = rp.fDaughter[1];
  fWeight      = rp.fWeight;
  fCalcMass    = rp.fCalcMass;
  fPx          = rp.fPx;
  fPy          = rp.fPy;
  fPz          = rp.fPz;
  fE           = rp.fE;
  fVx          = rp.fVx;
  fVy          = rp.fVy;
  fVz          = rp.fVz;
  fVt          = rp.fVt;
  fPolarTheta  = rp.fPolarTheta;
  fPolarPhi    = rp.fPolarPhi;
  fParticlePDG = rp.fParticlePDG; 
}

//____________________________________________________________________________
Int_t * AliPHOSRecParticle::GetPrimaries(Int_t & number) 
{
  // Retrieves all the primary particles at the origine of this reconstructed particle

  AliPHOSTrackSegment * ts = GetPHOSTrackSegment() ;

  Int_t emcnumber = 0 ; 
  Int_t * emclist = ts->GetPrimariesEmc(emcnumber) ;
  
  Int_t ppsdlnumber = 0 ;
  Int_t * ppsdllist = ts->GetPrimariesPpsdLow(ppsdlnumber) ;
 
  Int_t ppsdunumber = 0 ; 
  Int_t * ppsdulist = ts->GetPrimariesPpsdUp(ppsdunumber) ;

  number = emcnumber + ppsdlnumber + ppsdunumber ;
  Int_t * list   = new Int_t[number] ;
  
  Int_t index ; 
  for ( index = 0 ; index < emcnumber ; index++)
    list[index] = emclist[index] ;

  Int_t jndex ; 
  for ( jndex = 0 ; jndex < ppsdlnumber ; jndex++) {
    assert(index < number) ;
    list[index] = ppsdllist[jndex] ;
    index++ ; 
  }

  for ( jndex = 0 ; jndex < ppsdunumber ; jndex++) {
    assert(index < number) ;
    list[index] = ppsdulist[jndex] ;
    index++ ; 
  }

  delete emclist ;
  delete ppsdllist ;
  delete ppsdulist ;

  return list ; 
}

//____________________________________________________________________________
AliPHOSTrackSegment * AliPHOSRecParticle::GetPHOSTrackSegment() const 
{
  // Retrieves the PHOS track segment at the origine of this reconstructed particle

  AliPHOSIndexToObject * please = AliPHOSIndexToObject::GetInstance() ;
  return please->GimeTrackSegment( fPHOSTrackSegment ) ;

}



