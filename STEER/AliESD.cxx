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

//-----------------------------------------------------------------
//           Implementation of the ESD class
//   This is the class to deal with during the phisical analysis of data
//
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------

#include "AliESD.h"

ClassImp(AliESD)

//_______________________________________________________________________
AliESD::AliESD():
  fEventNumber(0),
  fRunNumber(0),
  fTrigger(0),
  fRecoVersion(0),
  fMagneticField(0),
  fTracks("AliESDtrack",15000),
  fCaloTracks("AliESDCaloTrack",500),
  fMuonTracks("AliESDMuonTrack",30),
  fV0s("AliESDv0",200),
  fCascades("AliESDcascade",20)
{
  Int_t i;
  for (i=0; i<3; i++) fVtx[i]=0.;
  for (i=0; i<6; i++) fCovVtx[i]=0.;
}

void AliESD::SetVertex(const Double_t *vtx, const Double_t *cvtx) {
  //Save the primary vertex position
  Int_t i;
  for (i=0; i<3; i++) fVtx[i]=vtx[i];
  if (cvtx)
  for (i=0; i<6; i++) fCovVtx[i]=cvtx[i];   
}

void AliESD::GetVertex(Double_t *vtx, Double_t *cvtx) const {
  //Get the primary vertex position
  Int_t i;
  for (i=0; i<3; i++) vtx[i]=fVtx[i];
  if (cvtx)
  for (i=0; i<6; i++) cvtx[i]=fCovVtx[i];   
}

