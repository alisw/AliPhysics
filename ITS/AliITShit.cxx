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
Revision 1.20  2002/10/22 14:45:42  alibrary
Introducing Riostream.h

Revision 1.19  2002/10/14 14:57:00  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.13.6.3  2002/08/28 15:06:50  alibrary
Updating to v3-09-01

Revision 1.18  2002/08/07 18:37:53  nilsen
Removed endl from print function. should be supplied by user as wanted.

Revision 1.17  2002/06/20 09:10:14  hristov
Data member ft0 initialized

Revision 1.16  2002/06/19 21:12:37  nilsen
Fixed bug with non-zero-ed new data members in constructors. Thanks Jiri
for finding it and pointing it out.

Revision 1.15  2002/06/12 18:59:47  nilsen
Added Starting track location to hit class and related changes to modules.
This is at present still fully backwards compatible since starting hits
are still written to the file. When aliroot v4.0 will be released, this
backwards compatiblity will be broken by removing the enterence hit, and making
the nessesary changes to module at that time.

Revision 1.14  2002/05/19 18:17:03  hristov
Changes needed by ICC/IFC compiler (Intel)

Revision 1.13  2002/03/09 18:35:35  nilsen
Added functions to print out Hit data members.

Revision 1.12  2002/03/08 16:05:05  nilsen
Standeard io streamers added to make debugging et al. easier.

Revision 1.11  2001/01/30 09:23:13  hristov
Streamers removed (R.Brun)

Revision 1.10  2001/01/26 20:01:19  hristov
Major upgrade of AliRoot code

Revision 1.9  2000/10/02 16:32:51  barbera
Automatic streamer used and forward declarations added

Revision 1.3.4.7  2000/10/02 15:54:49  barbera
Automatic streamer used and forward declarations added

Revision 1.8  2000/09/22 12:35:21  nilsen
Traps placed incase it is used without a properly initilized AliITSgeom class.

Revision 1.7  2000/07/10 16:07:18  fca
Release version of ITS code

Revision 1.3.4.2  2000/03/04 23:43:57  nilsen
Fixed up the comments/documentation.

Revision 1.3.4.1  2000/01/12 19:03:32  nilsen
This is the version of the files after the merging done in December 1999.
See the ReadMe110100.txt file for details

Revision 1.3  1999/09/29 09:24:20  fca
Introduction of the Copyright and cvs Log

*/

#include <Riostream.h>

#include <TMath.h>
#include <TRandom.h>
#include <TVector.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TTUBE.h>
#include "TParticle.h"

#include "AliRun.h"
#include "AliITS.h"
#include "AliITSgeom.h"
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
AliITShit::AliITShit():AliHit(){
    // Default Constructor
    // Zero data member just to be safe.
    // Intputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    A default created AliITShit class.

    fStatus = 0; // Track Status
    fLayer  = 0;  // Layer number
    fLadder = 0; // Ladder number
    fDet    = 0;    // Detector number  
    fPx     = 0.0;     // PX of particle at the point of the hit
    fPy     = 0.0;     // PY of particle at the point of the hit
    fPz     = 0.0;     // PZ of particle at the point of the hit
    fDestep = 0.0; // Energy deposited in the current step
    fTof    = 0.0;    // Time of flight at the point of the hit
    fStatus0 = 0; // zero status bit by default.
    fx0     = 0.0;     // Starting point of this step
    fy0     = 0.0;     // Starting point of this step
    fz0     = 0.0;     // Starting point of this step
    ft0     = 0.0;     // Starting point of this step
}
AliITShit::AliITShit(Int_t shunt,Int_t track,Int_t *vol,Float_t edep,
		     Float_t tof,TLorentzVector &x,TLorentzVector &x0,
		     TLorentzVector &p) : AliHit(shunt, track){
////////////////////////////////////////////////////////////////////////
// Create ITS hit
//     The creator of the AliITShit class. The variables shunt and
// track are passed to the creator of the AliHit class. See the AliHit
// class for a full description. In the units of the Monte Carlo
////////////////////////////////////////////////////////////////////////
    // Intputs:
    //    Int_t shunt   See AliHit
    //    Int_t track   Track number, see AliHit
    //    Int_t *vol     Array of integer hit data,
    //                     vol[0] Layer where the hit is, 1-6 typicaly
    //                     vol[1] Ladder where the hit is.
    //                     vol[2] Detector number where the hit is
    //                     vol[3] Set of status bits
    //                     vol[4] Set of status bits at start
    // Outputs:
    //    none.
    // Return:
    //    A default created AliITShit class.

    fLayer      = vol[0];  // Layer number
    fLadder     = vol[2];  // Ladder number
    fDet        = vol[1];  // Detector number
    fStatus     = vol[3];  // Track status flags
    fStatus0    = vol[4];  // Track status flag for start position.
    fX          = x.X();   // Track X global position
    fY          = x.Y();   // Track Y global position
    fZ          = x.Z();   // Track Z global position
    fPx         = p.Px();  // Track X Momentum
    fPy         = p.Py();  // Track Y Momentum
    fPz         = p.Pz();  // Track Z Momentum
    fDestep     = edep;    // Track dE/dx for this step
    fTof        = tof   ;  // Track Time of Flight for this step
    fx0         = x0.X();  // Track X global position
    fy0         = x0.Y();  // Track Y global position
    fz0         = x0.Z();  // Track Z global position
    ft0         = x0.T();     // Starting point of this step
}
//______________________________________________________________________
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
    // Intputs:
    //    Int_t shunt   See AliHit
    //    Int_t track   Track number, see AliHit
    //    Int_t *vol     Array of integer hit data,
    //                     vol[0] Layer where the hit is, 1-6 typicaly
    //                     vol[1] Ladder where the hit is.
    //                     vol[2] Detector number where the hit is
    //                     vol[3] Set of status bits
    //    Float_t *hits   Array of hit information
    //                     hits[0] X global position of this hit
    //                     hits[1] Y global position of this hit
    //                     hits[2] Z global position of this hit
    //                     hits[3] Px global position of this hit
    //                     hits[4] Py global position of this hit
    //                     hits[5] Pz global position of this hit
    //                     hits[6] Energy deposited by this step
    //                     hits[7] Time of flight of this particle at this step
    // Outputs:
    //    none.
    // Return:
    //    A standard created AliITShit class.
  fLayer      = vol[0];   // Layer number
  fLadder     = vol[2];   // Ladder number
  fDet        = vol[1];   // Detector number
  fStatus     = vol[3];   // Track status flags
  fX          = hits[0];  // Track X global position
  fY          = hits[1];  // Track Y global position
  fZ          = hits[2];  // Track Z global position
  fPx         = hits[3];  // Track X Momentum
  fPy         = hits[4];  // Track Y Momentum
  fPz         = hits[5];  // Track Z Momentum
  fDestep     = hits[6];  // Track dE/dx for this step
  fTof        = hits[7];  // Track Time of Flight for this step
  fStatus0 = 0;// Track Status of Starting point
  fx0 = 0.0;     // Starting point of this step
  fy0 = 0.0;     // Starting point of this step
  fz0 = 0.0;     // Starting point of this step
  ft0 = 0.0;     // Starting point of this step
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
    if(gm) {
      gm->GtoL(fLayer,fLadder,fDet,g,l);
      x = l[0];
      y = l[1];
      z = l[2];
    } else {
      Error("AliITShit","NULL pointer to the geometry! return smth else",gm);
      // AliITSv7 - SDD case
      x=fX;
      y=fZ;
      z=fY;
    }
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
    if(gm) {
      gm->GtoL(fLayer,fLadder,fDet,g,l);
      x = l[0];
      y = l[1];
      z = l[2];
    } else {
      Error("AliITShit","NULL pointer to the geometry! return smth else",gm);
      // AliITSv7 - SDD case
      x=fX;
      y=fZ;
      z=fY;
    }
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
    if(gm) {
      gm->GtoL(fLayer,fLadder,fDet,g,l);
      return l[0];
    } else {
      Error("AliITShit","NULL pointer to the geometry! return smth else",gm);
      return fX;
    }
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
    if (gm) {
      gm->GtoL(fLayer,fLadder,fDet,g,l);
      return l[1];
    } else {
      Error("AliITShit","NULL pointer to the geometry! return smth else",gm);
      return fZ;
    }
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
    if(gm) {
      gm->GtoL(fLayer,fLadder,fDet,g,l);
      return l[2];
    } else {
      Error("AliITShit","NULL pointer to the geometry! return smth else",gm);
      return fY;
    }
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
    if (gm) {
      gm->GtoLMomentum(fLayer,fLadder,fDet,g,l);
      px = l[0];
      py = l[1];
      pz = l[2];
    } else {
      Error("AliITShit","NULL pointer to the geometry! return smth else",gm);
      px=fPx;
      py=fPy;
      pz=fPz;
    }
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
    if (gm) {
      gm->GtoLMomentum(fLayer,fLadder,fDet,g,l);
      return l[0];
    } else {
      Error("AliITShit","NULL pointer to the geometry! return smth else",gm);
      return fPx;
    }
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
    if (gm) {
      gm->GtoLMomentum(fLayer,fLadder,fDet,g,l);
      return l[1];
    } else {
      Error("AliITShit","NULL pointer to the geometry! return smth else",gm);
      return fPy;
    }

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
    if (gm) {
      gm->GtoLMomentum(fLayer,fLadder,fDet,g,l);
      return l[2];
    } else {
      Error("AliITShit","NULL pointer to the geometry! return smth else",gm);
      return fPz;
    }

}
//___________________________________________________________________________;
Int_t AliITShit::GetModule(){
////////////////////////////////////////////////////////////////////////
//     Returns the module index number of the module where this hit was in.
////////////////////////////////////////////////////////////////////////
    AliITSgeom *gm = ((AliITS*)gAlice->GetDetector("ITS"))->GetITSgeom();

    if (gm) return gm->GetModuleIndex(fLayer,fLadder,fDet);
    else {
      Error("AliITShit","NULL pointer to the geometry! return smth else",gm);
      return 0;
    }
}
//______________________________________________________________________
TParticle * AliITShit::GetParticle(){
////////////////////////////////////////////////////////////////////////
//     Returns the pointer to the TParticle for the particle that created
// this hit. From the TParticle all kinds of information about this 
// particle can be found. See the TParticle class.
////////////////////////////////////////////////////////////////////////
    return gAlice->Particle(GetTrack());
}  
//----------------------------------------------------------------------
void AliITShit::Print(ostream *os){
////////////////////////////////////////////////////////////////////////
// Standard output format for this class.
////////////////////////////////////////////////////////////////////////
#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif
 
    fmt = os->setf(ios::scientific);  // set scientific floating point output
    *os << fTrack << " " << fX << " " << fY << " " << fZ << " ";
    fmt = os->setf(ios::hex); // set hex for fStatus only.
    *os << fStatus << " ";
    fmt = os->setf(ios::dec); // every thing else decimel.
    *os << fLayer << " " << fLadder << " " << fDet << " ";;
    *os << fPx << " " << fPy << " " << fPz << " ";
    *os << fDestep << " " << fTof;
    *os << " " << fx0 << " " << fy0 << " " << fz0;
//    *os << " " << endl;
    os->flags(fmt); // reset back to old formating.
    return;
}
//----------------------------------------------------------------------
void AliITShit::Read(istream *is){
////////////////////////////////////////////////////////////////////////
// Standard input format for this class.
////////////////////////////////////////////////////////////////////////
 

    *is >> fTrack >> fX >> fY >> fZ;
    *is >> fStatus >> fLayer >> fLadder >> fDet >> fPx >> fPy >> fPz >>
	   fDestep >> fTof;
    *is >> fx0 >> fy0 >> fz0;
    return;
}
//----------------------------------------------------------------------
ostream &operator<<(ostream &os,AliITShit &p){
////////////////////////////////////////////////////////////////////////
// Standard output streaming function.
////////////////////////////////////////////////////////////////////////
 
    p.Print(&os);
    return os;
}
//----------------------------------------------------------------------
istream &operator>>(istream &is,AliITShit &r){
////////////////////////////////////////////////////////////////////////
// Standard input streaming function.
////////////////////////////////////////////////////////////////////////
 
    r.Read(&is);
    return is;
}
//----------------------------------------------------------------------
