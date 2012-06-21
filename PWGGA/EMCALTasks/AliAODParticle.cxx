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
/* $Id: $ */
 
//_________________________________________________________________________

#include "AliAODParticle.h"
ClassImp(AliAODParticle) 
//===============================================
AliAODParticle::AliAODParticle():TLorentzVector(),
fChargedSign(0),
fL0(0), fL1(0), fModule(0), fBadDist(0), fNCells(0), fClusterTime(0), fClusterID(0), fAODClusterID(0),
kInSSA(0), kInTOF(0),  kInTrackMatched(0),
fPhotonPairDTime(0),  fPhotonPairDModule(0), fPhotonPairAsy(0), fPhotonPairAngle(0), fPhotonPairID0(0), fAODPhotonPairID0(0),fPhotonPairID1(0), fAODPhotonPairID1(0),
kIsLeading(0),   kIsIsolated()
{


} 
//===============================================
AliAODParticle::AliAODParticle(Double_t px,Double_t py,Double_t pz,Double_t energy):
TLorentzVector(px,py,pz,energy),
fChargedSign(0),
fL0(0), fL1(0), fModule(0),fBadDist(0),fNCells(0), fClusterTime(0), fClusterID(0), fAODClusterID(0),
kInSSA(0), kInTOF(0),  kInTrackMatched(0),
fPhotonPairDTime(0),   fPhotonPairDModule(0), fPhotonPairAsy(0), fPhotonPairAngle(0), fPhotonPairID0(0), fAODPhotonPairID0(0), fPhotonPairID1(0),fAODPhotonPairID1(0),
kIsLeading(0),   kIsIsolated(0)
{
  
}
//===============================================
Bool_t AliAODParticle::IsPIDOK(Int_t ipid)const{
  // returns true if photon satisfies given PID criterium
  switch(ipid){
  case 0: return  kTRUE ; //No PID at all
  case 1: return  kInSSA ;   //only shower shape cut
  case 2: return  !kInTrackMatched ;    //Only track matched cut
  case 3: return  kInTOF ;    //Only TOF cut
  case 4: return  kInSSA && !kInTrackMatched ;  //shower shape and CPV
  case 5: return  kInSSA && kInTOF ;  //shower shape and TOF
  case 7: return  !kInTrackMatched && kInTOF;  //CPV and TOF
  case 8: return  kInSSA && !kInTrackMatched && kInTOF; // all 3 cuts
  default: return kFALSE ; //Not known combination
  }
}
