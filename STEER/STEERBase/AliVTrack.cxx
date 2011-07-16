/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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


//-------------------------------------------------------------------------
//     base class for ESD and AOD tracks
//     Author: A. Dainese
//-------------------------------------------------------------------------

#include <TGeoGlobalMagField.h>

#include "AliMagF.h"
#include "AliVTrack.h"

ClassImp(AliVTrack)

AliVTrack::AliVTrack(const AliVTrack& vTrack) :
  AliVParticle(vTrack) { } // Copy constructor

AliVTrack& AliVTrack::operator=(const AliVTrack& vTrack)
{ if (this!=&vTrack) { 
    AliVParticle::operator=(vTrack); 
  }
  
  return *this; 
}

Double_t AliVTrack::GetBz() const 
{
  // returns Bz component of the magnetic field (kG)
  AliMagF* fld = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!fld) return 0.5*kAlmost0Field;
  double bz;
  if (fld->IsUniform()) bz = fld->SolenoidField();
  else {
    Double_t r[3]; 
    GetXYZ(r); 
    bz = fld->GetBz(r);
  }
  return TMath::Sign(0.5*kAlmost0Field,bz) + bz;
}

void AliVTrack::GetBxByBz(Double_t b[3]) const 
{
  // returns the Bx, By and Bz components of the magnetic field (kG)
  AliMagF* fld = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!fld) {
     b[0] = b[1] = 0.;
     b[2] = 0.5*kAlmost0Field;
     return;
  }

  if (fld->IsUniform()) {
     b[0] = b[1] = 0.;
     b[2] = fld->SolenoidField();
  }  else {
     Double_t r[3]; GetXYZ(r);
     fld->Field(r,b);
  }
  b[2] = (TMath::Sign(0.5*kAlmost0Field,b[2]) + b[2]);
  return;
}

