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
*/


#include <TMath.h>
#include <TRandom.h>
#include <TVector.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TTUBE.h>

#include "AliITSgeom.h"
#include "AliITS.h"
#include "AliITShit.h"
#include "AliRun.h"


ClassImp(AliITShit)
////////////////////////////////////////////////////////////////////////
// Version: 0
// Written by Rene Brun, Federico Carminati, and Roberto Barbera
//
// Version: 1
// Modified and documented by Bjorn S. Nilsen
// July 11 1999
//
// AliITShit is the hit class for the ITS. Hits are the information
// that comes from a Monte Carlo at each step as a particle mass through
// sensitive detector elements as particles are transported through a
// detector.
//
//Begin_Html
/*
<img src="picts/ITS/AliITShit_Class_Diagram.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This show the relasionships between the ITS hit class and the rest of Aliroot.
</font>
<pre>
*/
//End_Html
//_____________________________________________________________________________
AliITShit::AliITShit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
  AliHit(shunt, track){
  //
  // Create ITS hit
  //     The creator of the AliITShit class. The variables shunt and
  // track are passed to the creator of the AliHit class. See the AliHit
  // class for a full description. the integer array *vol contains, in order,
  // fLayer = vol[0], fDet = vol[1], fLadder = vol[2], fStatus = vol[3].
  // The array *hits contains, in order, fX = hits[0], fY = hits[1], 
  // fZ = hits[2], fPx = hits[3], fPy = hits[4], fPz = hits[5],
  // fDestep = hits[6], and fTof = hits[7].
  //
  fLayer      = vol[0];   // Layer number
  fLadder     = vol[2];   // Ladder number
  fDet        = vol[1];   // Detector number
  fStatus     = vol[3];   // Track status flags
  fX          = hits[0];  // Track X position
  fY          = hits[1];  // Track Y position
  fZ          = hits[2];  // Track Z position
  fPx         = hits[3];  // Track X Momentum
  fPy         = hits[4];  // Track Y Momentum
  fPz         = hits[5];  // Track Z Momentum
  fDestep     = hits[6];  // Track dE/dx for this step
  fTof        = hits[7];  // Track Time of Flight for this step
}

void AliITShit::GetPositionL(Float_t &x,Float_t &y,Float_t &z){
    AliITSgeom *gm = ((AliITS*)gAlice->GetDetector("ITS"))->GetITSgeom();
    Float_t g[3],l[3];

    g[0] = fX;
    g[1] = fY;
    g[2] = fZ;
    gm->GtoL(fLayer,fLadder,fDet,g,l);
    x = l[0];
    y = l[1];
    z = l[2];
    return;
}

void AliITShit::GetPositionL(Float_t &x,Float_t &y,Float_t &z,Float_t &tof){
    AliITSgeom *gm = ((AliITS*)gAlice->GetDetector("ITS"))->GetITSgeom();
    Float_t g[3],l[3];

    g[0] = fX;
    g[1] = fY;
    g[2] = fZ;
    gm->GtoL(fLayer,fLadder,fDet,g,l);
    x = l[0];
    y = l[1];
    z = l[2];
    tof = fTof;
    return;
}

Float_t AliITShit::GetXL(){
    AliITSgeom *gm = ((AliITS*)gAlice->GetDetector("ITS"))->GetITSgeom();
    Float_t g[3],l[3];

    g[0] = fX;
    g[1] = fY;
    g[2] = fZ;
    gm->GtoL(fLayer,fLadder,fDet,g,l);
    return l[0];
}

Float_t AliITShit::GetYL(){
    AliITSgeom *gm = ((AliITS*)gAlice->GetDetector("ITS"))->GetITSgeom();
    Float_t g[3],l[3];

    g[0] = fX;
    g[1] = fY;
    g[2] = fZ;
    gm->GtoL(fLayer,fLadder,fDet,g,l);
    return l[1];
}

Float_t AliITShit::GetZL(){
    AliITSgeom *gm = ((AliITS*)gAlice->GetDetector("ITS"))->GetITSgeom();
    Float_t g[3],l[3];

    g[0] = fX;
    g[1] = fY;
    g[2] = fZ;
    gm->GtoL(fLayer,fLadder,fDet,g,l);
    return l[2];
}

void AliITShit::GetMomentumL(Float_t &px,Float_t &py,Float_t &pz){
    AliITSgeom *gm = ((AliITS*)gAlice->GetDetector("ITS"))->GetITSgeom();
    Float_t g[3],l[3];

    g[0] = fPx;
    g[1] = fPy;
    g[2] = fPz;
    gm->GtoLMomentum(fLayer,fLadder,fDet,g,l);
    px = l[0];
    py = l[1];
    pz = l[2];
    return;
}
