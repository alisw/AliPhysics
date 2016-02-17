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


/////////////////////////////////////////////////////////////////////////////////
// This class sets the local coordinates via a specific setter. Needed because //
// the AliGeomManager class can not be used for the upgrade code at this stage //
/////////////////////////////////////////////////////////////////////////////////

<<<<<<< HEAD:ITS/UPGRADE/ITSUpgradeBase/AliITSRecPointU.cxx
#include "AliITSRecPointU.h"
#include "AliLog.h"
=======
#include <AliITSRecPointU.h>
#include <AliLog.h>
>>>>>>> feature-itsmft:ITSMFT/ITS/v0/AliITSRecPointU.cxx

ClassImp(AliITSRecPointU)
//_____________________________________________________________
AliITSRecPointU::AliITSRecPointU():
    fXloc(0),
    fZloc(0),
    fModule(0),
    fNTracksIdMC(0)
{
 //
 // Default constructor
 // 
 for(Int_t i=0; i<kMaxLab ; i++) {
    fTrackIdMC[i]=-3;
 }
}
//_____________________________________________________________
AliITSRecPointU::AliITSRecPointU(const AliITSRecPointU& pt):
    fXloc(pt.fXloc),
    fZloc(pt.fZloc),
    fModule(pt.fModule),
    fNTracksIdMC(pt.fNTracksIdMC)
{
  //
  // Copy constructor
  //
  for(Int_t i=0; i<kMaxLab ; i++) {
    fTrackIdMC[i]=pt.fTrackIdMC[i];
  }
}
//______________________________________________________________________
AliITSRecPointU& AliITSRecPointU::operator=(const AliITSRecPointU& source)
{
  //
  // Assignment operator (as in AliITSRecPoint)
  //

  this->~AliITSRecPointU();
  new(this) AliITSRecPointU(source);
  return *this;

}

//______________________________________________________________________________
void AliITSRecPointU::AddTrackID(Int_t tid) { 
  // 
  // Add an MC label (track ID) to the "expanded list"
  //
  if (fNTracksIdMC==kMaxLab) {
    AliWarning("Max. numbers of labels reached!"); 
  } else {
    fTrackIdMC[fNTracksIdMC]=tid; 
    fNTracksIdMC++;
  } 
}

//______________________________________________________________________________
void AliITSRecPointU::Print(Option_t* /*option*/) const
{
  // Print cluster information.
  
  printf("AliITSRecPointU pos=(%.4f, %.4f, %.4f), s_y2=%f, s_z2=%f, s_yz=%f, vol=%hu\n",
         GetX(), GetY(), GetZ(), GetSigmaY2(), GetSigmaZ2(), GetSigmaYZ(), GetVolumeId());
  printf("      MC Track Ids =(");
  if (kMaxLab<=fNTracksIdMC) {
    for (Int_t i=0; i<kMaxLab; i++) { printf("%d,",fTrackIdMC[i]); }
  } else {
    for (Int_t i=0; i<fNTracksIdMC; i++) { printf("%d,",fTrackIdMC[i]); }
  }  
  printf(")\n");
  Float_t g[3];
  if (GetGlobalXYZ(g))
    printf("    global_pos=(%.4f, %.4f, %.4f)\n", g[0], g[1], g[2]);
  
}


