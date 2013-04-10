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

#include "AliITSU.h"
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
AliITSUHit::AliITSUHit(Int_t shunt,Int_t track,Int_t *vol,Float_t edep,Float_t tof,
			   TLorentzVector &x,TLorentzVector &x0,TLorentzVector &p) 
: AliITShit(shunt,track,vol,edep,tof,x,x0,p)
{
  // ct-r
}

//______________________________________________________________________
AliITSUHit::AliITSUHit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits) 
  : AliITShit(shunt, track, vol, hits) 
{
  // c-tor
}

//______________________________________________________________________
AliITSUHit::AliITSUHit(const AliITSUHit &h)
: AliITShit(h)
{
  // cp c-tor
}

//______________________________________________________________________
AliITSUHit& AliITSUHit::operator=(const AliITSUHit &h)
{
  // The standard = operator
  if(this == &h) return *this;
  AliITShit::operator=(h);
  return *this;
}

//______________________________________________________________________
void AliITSUHit::GetPositionL(Float_t &x,Float_t &y,Float_t &z,Float_t &tof)
{
  // Returns the position and time of flight of this hit in the local
  // coordinates of this module, and in the units of the Monte Carlo.
  //
  AliITSUGeomTGeo *gm = ((AliITSU*)gAlice->GetDetector("ITS"))->GetITSGeomTGeo();
  if (!gm) AliFatal("NULL pointer to the geometry!");
  double g[3]={fX,fY,fZ},l[3];
  gm->GetMatrixSens(fModule)->MasterToLocal(g,l);
  x = l[0];
  y = l[1];
  z = l[2];
  tof = fTof;
  //
}

//______________________________________________________________________
void AliITSUHit::GetPositionL0(Double_t &x,Double_t &y,Double_t &z,Double_t &tof)
{
  // Returns the initial position and time of flight of this hit 
  // in the local coordinates of this module, and in the units of the 
  AliITSUGeomTGeo *gm = ((AliITSU*)gAlice->GetDetector("ITS"))->GetITSGeomTGeo();
  if (!gm) AliFatal("NULL pointer to the geometry!");
  double g[3]={fx0,fy0,fz0},l[3];  
  gm->GetMatrixSens(fModule)->MasterToLocal(g,l);
  x = l[0];
  y = l[1];
  z = l[2];
  tof = ft0;
}

//______________________________________________________________________
void AliITSUHit::GetDetectorID(Int_t &layer,Int_t &ladder,Int_t &det) const
{
  // Returns the layer ladder and detector number lables for this
  // ITS module. Note: indices start from 0!
  AliITSUGeomTGeo *gm = ((AliITSU*)gAlice->GetDetector("ITS"))->GetITSGeomTGeo();
  if (!gm) AliFatal("NULL pointer to the geometry!");
  gm->GetModuleId(fModule,layer,ladder,det);
}  

//______________________________________________________________________
Int_t AliITSUHit::GetLayer() const
{
  // Returns the layer. Note: indices start from 0!
  AliITSUGeomTGeo *gm = ((AliITSU*)gAlice->GetDetector("ITS"))->GetITSGeomTGeo();
  if (!gm) AliFatal("NULL pointer to the geometry!");
  return gm->GetLayer(fModule);
}  

//______________________________________________________________________
Int_t AliITSUHit::GetLadder() const
{
  // Returns the ladder of TS module. Note: indices start from 0!
  AliITSUGeomTGeo *gm = ((AliITSU*)gAlice->GetDetector("ITS"))->GetITSGeomTGeo();
  if (!gm) AliFatal("NULL pointer to the geometry!");
  return gm->GetLadder(fModule);
}  

//______________________________________________________________________
Int_t AliITSUHit::GetDetector() const
{
  // Returns the detector within the ladder. Note: indices start from 0!
  AliITSUGeomTGeo *gm = ((AliITSU*)gAlice->GetDetector("ITS"))->GetITSGeomTGeo();
  if (!gm) AliFatal("NULL pointer to the geometry!");
  return gm->GetModIdInLadder(fModule);
}  

//______________________________________________________________________
void AliITSUHit::Print(Option_t *option) const 
{
  // print itself
  printf("Mod%4d Tr:%5d DE:%.2e TOF: %.3e| P:%.3f %.3f %.3f |>%.4f %.4f %.4f >%.4f %.4f %.4f\n",
	 fModule,fTrack,fDestep,fTof,fPx,fPy,fPz, fx0,fy0,fz0,fX,fY,fZ);

}
