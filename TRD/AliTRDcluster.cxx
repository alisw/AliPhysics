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
Revision 1.1.2.1  2000/09/22 14:47:52  cblume
Add the tracking code

*/

#include "AliTRDcluster.h"
#include "AliTRDgeometry.h"
#include "AliTRDrecPoint.h"


ClassImp(AliTRDcluster)
 
//_____________________________________________________________________________
AliTRDcluster::AliTRDcluster() {
  //default constructor

  fDetector = fTimeBin = 0;
  fTracks[0]=fTracks[1]=fTracks[2]=0; 
  fY=fZ=fQ=fSigmaY2=fSigmaZ2=0.;
}

//_____________________________________________________________________________
AliTRDcluster::AliTRDcluster(AliTRDrecPoint *rp)
{
  //
  // constructor from AliTRDrecPoint
  //

  fDetector   = rp->GetDetector();
  fTimeBin    = rp->GetLocalTimeBin();

  fTracks[0]  = rp->GetTrackIndex(0);
  fTracks[1]  = rp->GetTrackIndex(1);
  fTracks[2]  = rp->GetTrackIndex(2);

  fQ          = rp->GetEnergy();

  fY          = rp->GetY();
  fZ          = rp->GetZ();
  fSigmaY2    = rp->GetSigmaY2();
  fSigmaZ2    = rp->GetSigmaZ2();  

  fSigmaY2    = 1;
  fSigmaZ2    = 5;  

}

