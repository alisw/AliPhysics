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

/*
$Log$
Revision 1.1.2.2  2000/06/25 08:38:41  kowal2
Splitted from AliTPCtracking

*/

//-----------------------------------------------------------------
//           Implementation of the TPC cluster class
//
// Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------

#include "AliTPCcluster.h"

ClassImp(AliTPCcluster)
 
//_____________________________________________________________________________
AliTPCcluster::AliTPCcluster() {
  //default constructor
  fTracks[0]=fTracks[1]=fTracks[2]=0; 
  fY=fZ=fQ=fSigmaY2=fSigmaZ2=0.;
}

//_____________________________________________________________________________
AliTPCcluster::AliTPCcluster(Float_t *hits, Int_t *lab)
{
  //
  // Creates a simulated cluster for the TPC
  //
  fTracks[0]  = lab[0];
  fTracks[1]  = lab[1];
  fTracks[2]  = lab[2];
  fY          = hits[0];
  fZ          = hits[1];
  fQ          = hits[2];
  fSigmaY2    = hits[3];
  fSigmaZ2    = hits[4];
}

