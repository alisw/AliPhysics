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
Revision 1.3.4.2  2000/03/04 23:43:57  nilsen
Fixed up the comments/documentation.

Revision 1.3.4.1  2000/01/12 19:03:32  nilsen
This is the version of the files after the merging done in December 1999.
See the ReadMe110100.txt file for details

Revision 1.3  1999/09/29 09:24:20  fca
Introduction of the Copyright and cvs Log

*/


#include <TMath.h>
#include <TRandom.h>
#include <TVector.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TTUBE.h>

#include "AliRun.h"
#include "AliITSgeom.h"
#include "AliITS.h"
#include "AliITShit.h"


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
////////////////////////////////////////////////////////////////////////
// Inline Member functions:
//
// AliITShit()
//     The default creator of the AliITShit class.
//
// ~AliITShit()
//     The default destructor of the AliITShit class.
//
// int GetTrack()
//     See AliHit for a full description. Returns the track number fTrack
// for this hit.
//
// SetTrack(int track)
//     See AliHit for a full description. Sets the track number fTrack
// for this hit.
//
// Int_t GetTrackStatus()
//     Returns the value of the track status flag fStatus. This flag
// indicates the track status at the time of creating this hit. It is
// made up of the following 8 status bits from highest order to lowest
// order bits
// 0           :  IsTrackAlive():    IsTrackStop():IsTrackDisappeared():
// IsTrackOut():IsTrackExiting():IsTrackEntering():IsTrackInside()     .
// See AliMC for a description of these functions. If the function is
// true then the bit is set to one, otherwise it is zero.
//
// Bool_t StatusInside()
//     Returns kTRUE if the particle producing this hit is still inside
// the present volume. Returns kFalse if this particle will be in another
// volume. {bit IsTrackInside is set or not}
//
// Bool_t StatusEntering()
//     Returns kTRUE if the particle producing this hit is has just enterd
// the present volume. Returns kFalse otherwise.  {bit IsTrackEntering is
// set or not}
//
// Bool_t StatusExiting()
//     Returns kTRUE if the particle producing this hit is will exit
// the present volume. Returns kFalse otherwise. {bit IsTrackExiting is set
// or not}
//
// Bool_t StatusOut()
//     Returns kTRUE if the particle producing this hit is goint exit the
// simulation. Returns kFalse otherwise. {bit IsTrackOut is set or not}
//
// Bool_t StatusDisappeared()
//     Returns kTRUE if the particle producing this hit is going to "disappear"
// for example it has interacted producing some other particles. Returns
//  kFalse otherwise. {bit IsTrackOut is set or not}
//
// Bool_t StatusStop()
//     Returns kTRUE if the particle producing this hit is has dropped below
// its energy cut off producing some other particles. Returns kFalse otherwise.
// {bit IsTrackOut is set or not}
//
// Bool_t StatuAlives()
//     Returns kTRUE if the particle producing this hit is going to continue
// to be transported. Returns kFalse otherwise. {bit IsTrackOut is set or not}
//
// Int_t GetLayer()
//     Returns the layer number, fLayer, for this hit.
//
// Int_t GetLadder()
//     Returns the ladder number, fLadder, for this hit.
//
// Int_t GetDetector()
//     Returns the detector number, fDet, for this hit.
//
// GetDetectorID(Int_t &layer, Int_t &ladder, Int_t &detector)
//     Returns the layer, ladder, and detector numbers, fLayer fLadder fDet,
// in one call.
//
// Float_t GetIonization()
//     Returns the energy lost, fDestep, by the particle creating this hit,
// in the units defined by the Monte Carlo.
//
// GetPositionG(Float_t &x, Float_t &y, Float_t &z)
//     Returns the global position, fX fY fZ, of this hit, in the units
// define by the Monte Carlo.
//
// GetPositionG(Double_t &x, Double_t &y, Double_t &z)
//     Returns the global position, fX fY fZ, of this hit, in the units
// define by the Monte Carlo.
//
// GetPositionG(Float_t &x, Float_t &y, Float_t &z, Float_t &tof)
//     Returns the global position and time of flight, fX fY fZ fTof, of
// this hit, in the units define by the Monte Carlo.
//
// GetPositionG(Double_t &x,Double_t &y,Double_t &z,Double_t &tof)
//     Returns the global position and time of flight, fX fY fZ fTof, of
// this hit, in the units define by the Monte Carlo.
//
// GetPositionL(Double_t &x,Double_t &y,Double_t &z)
//     Returns the local position, fX fY fZ, of this hit in the coordiante
// of this module, in the units define by the Monte Carlo.
//
// GetPositionG(Double_t &x,Double_t &y,Double_t &z,Double_t &tof)
//     Returns the local position and time of flight, fX fY fZ fTof, of
// this hit in the coordinates of this module, in the units define by the
//  Monte Carlo.
//
// Float_t GetXG()
//     Returns the global x position in the units defined by the Monte Carlo.
//
// Float_t GetYG()
//     Returns the global y position in the units defined by the Monte Carlo.
//
// Float_t GetYG()
//     Returns the global z position in the units defined by the Monte Carlo.
//
// Float_t GetTOF()
//     Returns the time of flight, fTof, of this hit, in the units defined
// by the Monte Carlo.
//
// GetMomentumG(Float_t &px, Float_t &py, Float_t &pz)
//     Returns the global momentum, fPx fPy fPz, of the particle that made
// this hit, in the units define by the Monte Carlo.
//
// GetMomentumG(Double_t &px,Double_t &py,Double_t &pz)
//     Returns the global momentum, fPx fPy fPz, of the particle that made
// this hit, in the units define by the Monte Carlo.
//
// GetMomentumL(Double_t &px,Double_t &py,Double_t &pz)
//     Returns the momentum, fPx fPy fPz in coordinate appropreate for this
// specific module, in the units define by the Monte Carlo.
//
// Float_t GetPXG()
//     Returns the global X momentum in the units defined by the Monte Carlo.
//
// Float_t GetPYG()
//     Returns the global Y momentum in the units defined by the Monte Carlo.
//
// Float_t GetPZG()
//     Returns the global Z momentum in the units defined by the Monte Carlo.
//
////////////////////////////////////////////////////////////////////////
//_____________________________________________________________________________
AliITShit::AliITShit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
    AliHit(shunt, track){
////////////////////////////////////////////////////////////////////////
// Create ITS hit
//     The creator of the AliITShit class. The variables shunt and
// track are passed to the creator of the AliHit class. See the AliHit
// class for a full description. the integer array *vol contains, in order,
// fLayer = vol[0], fDet = vol[1], fLadder = vol[2], fStatus = vol[3].
// The array *hits contains, in order, fX = hits[0], fY = hits[1],
// fZ = hits[2], fPx = hits[3], fPy = hits[4], fPz = hits[5],
// fDestep = hits[6], and fTof = hits[7]. In the units of the Monte Carlo
////////////////////////////////////////////////////////////////////////
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
//______________________________________________________________________
void AliITShit::GetPositionL(Float_t &x,Float_t &y,Float_t &z){
////////////////////////////////////////////////////////////////////////
//     Returns the position of this hit in the local coordinates of this
// module, and in the units of the Monte Carlo.
////////////////////////////////////////////////////////////////////////
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
//______________________________________________________________________
void AliITShit::GetPositionL(Float_t &x,Float_t &y,Float_t &z,Float_t &tof){
////////////////////////////////////////////////////////////////////////
//     Returns the position and time of flight of this hit in the local
// coordinates of this module, and in the units of the Monte Carlo.
////////////////////////////////////////////////////////////////////////
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
//______________________________________________________________________
Float_t AliITShit::GetXL(){
////////////////////////////////////////////////////////////////////////
//     Returns the x position of this hit in the local coordinates of this
// module, and in the units of the Monte Carlo.
////////////////////////////////////////////////////////////////////////
    AliITSgeom *gm = ((AliITS*)gAlice->GetDetector("ITS"))->GetITSgeom();
    Float_t g[3],l[3];

    g[0] = fX;
    g[1] = fY;
    g[2] = fZ;
    gm->GtoL(fLayer,fLadder,fDet,g,l);
    return l[0];
}
//______________________________________________________________________
Float_t AliITShit::GetYL(){
////////////////////////////////////////////////////////////////////////
//     Returns the y position of this hit in the local coordinates of this
// module, and in the units of the Monte Carlo.
////////////////////////////////////////////////////////////////////////
    AliITSgeom *gm = ((AliITS*)gAlice->GetDetector("ITS"))->GetITSgeom();
    Float_t g[3],l[3];

    g[0] = fX;
    g[1] = fY;
    g[2] = fZ;
    gm->GtoL(fLayer,fLadder,fDet,g,l);
    return l[1];
}
//______________________________________________________________________
Float_t AliITShit::GetZL(){
////////////////////////////////////////////////////////////////////////
//     Returns the z position of this hit in the local coordinates of this
// module, and in the units of the Monte Carlo.
////////////////////////////////////////////////////////////////////////
    AliITSgeom *gm = ((AliITS*)gAlice->GetDetector("ITS"))->GetITSgeom();
    Float_t g[3],l[3];

    g[0] = fX;
    g[1] = fY;
    g[2] = fZ;
    gm->GtoL(fLayer,fLadder,fDet,g,l);
    return l[2];
}
//______________________________________________________________________
void AliITShit::GetMomentumL(Float_t &px,Float_t &py,Float_t &pz){
////////////////////////////////////////////////////////////////////////
//     Returns the momentum of this hit in the local coordinates of this
// module, and in the units of the Monte Carlo.
////////////////////////////////////////////////////////////////////////
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
//______________________________________________________________________
Float_t AliITShit::GetPXL(){
////////////////////////////////////////////////////////////////////////
//     Returns the X momentum of this hit in the local coordinates of this
// module, and in the units of the Monte Carlo.
////////////////////////////////////////////////////////////////////////
    AliITSgeom *gm = ((AliITS*)gAlice->GetDetector("ITS"))->GetITSgeom();
    Float_t g[3],l[3];

    g[0] = fPx;
    g[1] = fPy;
    g[2] = fPz;
    gm->GtoLMomentum(fLayer,fLadder,fDet,g,l);
    return l[0];
}
//______________________________________________________________________
Float_t AliITShit::GetPYL(){
////////////////////////////////////////////////////////////////////////
//     Returns the Y momentum of this hit in the local coordinates of this
// module, and in the units of the Monte Carlo.
////////////////////////////////////////////////////////////////////////
    AliITSgeom *gm = ((AliITS*)gAlice->GetDetector("ITS"))->GetITSgeom();
    Float_t g[3],l[3];

    g[0] = fPx;
    g[1] = fPy;
    g[2] = fPz;
    gm->GtoLMomentum(fLayer,fLadder,fDet,g,l);
    return l[1];
}
//______________________________________________________________________
Float_t AliITShit::GetPZL(){
////////////////////////////////////////////////////////////////////////
//     Returns the Z momentum of this hit in the local coordinates of this
// module, and in the units of the Monte Carlo.
////////////////////////////////////////////////////////////////////////
    AliITSgeom *gm = ((AliITS*)gAlice->GetDetector("ITS"))->GetITSgeom();
    Float_t g[3],l[3];

    g[0] = fPx;
    g[1] = fPy;
    g[2] = fPz;
    gm->GtoLMomentum(fLayer,fLadder,fDet,g,l);
    return l[2];
}
//___________________________________________________________________________;
Int_t AliITShit::GetModule(){
////////////////////////////////////////////////////////////////////////
//     Returns the module index number of the module where this hit was in.
////////////////////////////////////////////////////////////////////////
    AliITSgeom *gm = ((AliITS*)gAlice->GetDetector("ITS"))->GetITSgeom();

    return gm->GetModuleIndex(fLayer,fLadder,fDet);
}
//______________________________________________________________________
TParticle * AliITShit::GetParticle(){
////////////////////////////////////////////////////////////////////////
//     Returns the pointer to the TParticle for the particle that created
// this hit. From the TParticle all kinds of information about this 
// particle can be found. See the TParticle class.
////////////////////////////////////////////////////////////////////////
    TClonesArray *Parts = gAlice->Particles();
    return ((TParticle*)Parts->UncheckedAt(this->GetTrack()));
}
//___________________________________________________________________________
void AliITShit::Streamer(TBuffer &R__b){
////////////////////////////////////////////////////////////////////////
//     The streamer function for this AliITShit class. This has been set
// up so that future changes to this class can be easly incorporated.
////////////////////////////////////////////////////////////////////////
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); 
      if (R__v==1) {
	  AliHit::Streamer(R__b);
	  R__b >> fStatus;
	  R__b >> fLayer;
	  R__b >> fLadder;
	  R__b >> fDet;
	  R__b >> fPx;
	  R__b >> fPy;
	  R__b >> fPz;
	  R__b >> fDestep;
	  R__b >> fTof;
      }else{ // for futrue changes to this class.
      } // end if R__v==1
   } else {
      R__b.WriteVersion(AliITShit::IsA());
      AliHit::Streamer(R__b);
      R__b << fStatus;
      R__b << fLayer;
      R__b << fLadder;
      R__b << fDet;
      R__b << fPx;
      R__b << fPy;
      R__b << fPz;
      R__b << fDestep;
      R__b << fTof;
   } // end if R__b.IsReading()
}
