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

//-------------------------------------------------------------------------
//               Implementation of the Cluster class
//
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-------------------------------------------------------------------------

#include "AliCluster.h"

ClassImp(AliCluster)
 
//_____________________________________________________________________________
AliCluster::AliCluster() {
  //default constructor
  fTracks[0]=fTracks[1]=fTracks[2]=-3141593; 
  fY=fZ=fSigmaY2=fSigmaZ2=0.;
}

//_____________________________________________________________________________
AliCluster::AliCluster(Int_t *lab, Float_t *hit) {
  //Creates a simulated cluster
  fTracks[0]  = lab[0];
  fTracks[1]  = lab[1];
  fTracks[2]  = lab[2];
  fY          = hit[0];
  fZ          = hit[1];
  fSigmaY2    = hit[2];
  fSigmaZ2    = hit[3];
}
