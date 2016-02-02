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

#include "AliLog.h"
#include "AliITSUGeomTGeo.h"
#include "AliITSUHit.h"

ClassImp(AliITSUHit)

////////////////////////////////////////////////////////////////////////
//
// At the moment the same functionality/data-members as parent AliITShit 
// except the geometry transformation uses AliITSgeomTGeoUp 
//
////////////////////////////////////////////////////////////////////////

//----------------------------------------------------------------------
AliITSUHit::AliITSUHit() : AliITSMFTHit()
{
}

//----------------------------------------------------------------------
AliITSUHit::AliITSUHit(Int_t shunt,Int_t track,Int_t *vol,Float_t edep,
Float_t tof, TLorentzVector &x,TLorentzVector &x0,TLorentzVector &p) 
: AliITSMFTHit( shunt, track, vol, edep, tof, x, x0, p)
{
}

//______________________________________________________________________
AliITSUHit::AliITSUHit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits) 
: AliITSMFTHit(shunt, track, vol, hits)
{
}

//______________________________________________________________________
AliITSUHit::AliITSUHit(const AliITSUHit &h)
: AliITSMFTHit(h)
{
}

//______________________________________________________________________
void AliITSUHit::GetChipID(Int_t &layer,Int_t &stave,Int_t &sstave, Int_t &mod,Int_t &det,
const AliITSUGeomTGeo *geom) const
{
  // Returns the layer stave and detector number lables for this
  // ITS chip. Note: indices start from 0!
  if (!geom) { AliFatal("NULL pointer to the geometry!"); return; }
  geom->GetChipId(fModule,layer,stave,sstave,mod,det);
}  
