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
Revision 1.14  1999/10/22 08:16:49  fca
Correct destructors, thanks to I.Hrivnacova

Revision 1.13  1999/10/06 10:15:19  fca
Correct bug in allocation of layer name and add destructor

Revision 1.12  1999/10/05 08:05:09  fca
Minor corrections for uninitialised variables.

Revision 1.11  1999/09/29 09:24:20  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Inner Traking System version 5                                           //
//  This class contains the base procedures for the Inner Tracking System    //
//                                                                           //
// Authors: R. Barbera, B. S. Nilsen.
// version 5.
// Created September 17 1999.
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <TMath.h>

#include "AliRun.h"
#include "TSystem.h"
#include "AliITShit.h"
#include "AliITS.h"
#include "AliITSv5.h"
#include "AliITSgeom.h"

ClassImp(AliITSv5)
 
//_____________________________________________________________________________
AliITSv5::AliITSv5() {
    //
    // Standard constructor for the ITS
    //
    fId5N = 6;
    fId5Name = new char*[fId5N];
    fId5Name[0] = "ITS1";
    fId5Name[1] = "ITS2";
    fId5Name[2] = "ITS3";
    fId5Name[3] = "ITS4";
    fId5Name[4] = "ITS5";
    fId5Name[5] = "ITS6";
}
//_____________________________________________________________________________
AliITSv5::~AliITSv5() {
    //
    // Standard destructor for the ITS
    //
  delete [] fId5Name;
}
//_____________________________________________________________________________
AliITSv5::AliITSv5(const char *name, const char *title) : AliITS(name, title){
    //
    // Standard constructor for the ITS
    //
    fId5N = 6;
    fId5Name = new char*[fId5N];
    fId5Name[0] = "ITS1";
    fId5Name[1] = "ITS2";
    fId5Name[2] = "ITS3";
    fId5Name[3] = "ITS4";
    fId5Name[4] = "ITS5";
    fId5Name[5] = "ITS6";

    fEuclidMaterial = "$ALICE_ROOT/Euclid/ITSgeometry_5.tme";
    fEuclidGeometry = "$ALICE_ROOT/Euclid/ITSgeometry_5.euc";
}

 
//_____________________________________________________________________________
void AliITSv5::CreateMaterials(){
  //
  // Read materials for the ITS
  //
    char *filtmp;
//
  filtmp = gSystem->ExpandPathName(fEuclidMaterial.Data());
//  FILE *file = fopen(fEuclidMaterial.Data(),"r");
  FILE *file = fopen(filtmp,"r");
  if(file) {
    fclose(file);
//    ReadEuclidMedia(fEuclidMaterial.Data());
    ReadEuclidMedia(filtmp);
  } else {
    Error("CreateMaterials"," THE MEDIA FILE %s DOES NOT EXIST !",
//	  fEuclidMaterial.Data());
	  filtmp);
    exit(1);
  } // end if(file)
}

//_____________________________________________________________________________
void AliITSv5::CreateGeometry(){
//////////////////////////////////////////////////////////////////////
//    This is the geometry used for the ITS Pre-TDR and comes from an 
// Euclid to Geant conversion. The only difference
// is in the details of the ITS supports. The detectors elements, 
// detector numbering, and local and global reference frames are shown in
//  the following figures.
//Begin_Html
/*
<img src="picts/ITS/its1+2_convention_front_5.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This shows the front view of the SPDs.
</font>
<pre>
<img src="picts/ITS/its1+2_convention_side_5.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This shows the perspective view of the SPDs.
</font>
<img src="picts/ITS/its1+2_tree.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This shows the geometry Tree for the SPDs.
</font>
<pre>

<pre>
<img src="picts/ITS/its3+4_convention_front_5.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This shows the front view of the SDDs.
</font>
<pre>
<img src="picts/ITS/its3+4_convention_side_5.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This shows the perspective view of the SDDs.
</font>
<img src="picts/ITS/its3+4_tree.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This shows the geometry Tree for the SDDs.
</font>
<pre>

<pre>
<img src="picts/ITS/its5+6_convention_front_5.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This shows the front view of the SSDs.
</font>
<pre>
<img src="picts/ITS/its5+6_convention_side_5.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This shows the perspective view of the SSDs.
</font>
<pre>
<img src="picts/ITS/its5+6_tree.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This shows the geometry Tree for the SSDs.
</font>
<pre>


<img src="picts/ITS/its_layer1-6_2.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This shows the front view of the whole ITS..
</font>
<pre>

<img src="picts/ITS/its_layer1-6_1.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This shows the perspective view of the whole ITS..
</font>
<pre>

<img src="picts/ITS/its1-6_tree.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This shows the geometry Tree for the whole ITS.
</font>
<pre>
*/
//End_Html
//
//
//      Here are shown the details of the ITS support cones and services.
// First is a GEANT tree showing the organization of all of the volumes
// that make up the ITS supports and services.
//Begin_Html
/*
<img src="picts/ITS/supports_tree.gif">
 */
//End_Html
//     What follows are a number of figures showing what these support
// structures look like.
//Begin_Html
/*

<img src="picts/ITS/supports_3.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This shows the geometry of the supports for the Drift and Strip layers only.
</font>
<pre>

<img src="picts/ITS/supports_2.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This shows the geometry of the supports for the Drift and Strip layers in front cut out.
</font>
<pre>

<img src="picts/ITS/supports_1.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This shows the geometry of the supports for the Drift and Strip layers in a back cut out..
</font>
<pre>

<img src="picts/ITS/suppssd.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This shows the geometry for the Strip layers supports.
</font>
<pre>

<img src="picts/ITS/suppsdd.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This shows the geometry for the Drift layers supports.
</font>
<pre>
 */
//End_Html
////////////////////////////////////////////////////////////////////////
// Read geometry for the ITS
//

    Int_t size;
    char topvol[5];
    char *filtmp;
//
  filtmp = gSystem->ExpandPathName(fEuclidGeometry.Data());
  FILE *file = fopen(filtmp,"r");
  delete [] filtmp;
  if(file) {
    fclose(file);
    printf("Ready to read Euclid geometry file\n");
    ReadEuclid(fEuclidGeometry.Data(),topvol);
    printf("Read in euclid geometries\n");
  } else {
    Error("CreateGeometry"," THE GEOM FILE %s DOES NOT EXIST !",
	  fEuclidGeometry.Data());
    exit(1);
  } // end if(file)
  //
  //---Place the ITS ghost volume ITSV in its mother volume (ALIC) and make it
  //     invisible
  //
  gMC->Gspos("ITSV",1,"ALIC",0,0,0,0,"ONLY");
  //
  //---Outputs the geometry tree in the EUCLID/CAD format
  
    if (fEuclidOut) {
      gMC->WriteEuclid("ITSgeometry", "ITSV", 1, 5);
    } // end if (fEuclidOut)

    filtmp = gSystem->ExpandPathName(fEuclidGeometry.Data());
    size = strlen(filtmp);
    if(size>4){
	filtmp[size-3] = 'd'; // change from .euc to .det
        filtmp[size-2] = 'e';
        filtmp[size-1] = 't';
	file = fopen(filtmp,"r");
	if(file){ // if file exists use it to fill AliITSgeom structure.
	    fclose(file);
	    printf("ready to read .det file %s\n",filtmp);
	    fITSgeom = new AliITSgeom(filtmp);
	}else{
	    fITSgeom = 0;
	    // fill AliITSgeom structure from geant structure just filled above
	}// end if(file)
        delete [] filtmp;
    }// end if(size>4)
    printf("finished with euclid geometrys\n");
}

//_____________________________________________________________________________
void AliITSv5::Init(){
    //
    // Initialise the ITS after it has been created
    //
    Int_t i,j,l;

    fIdN    = fId5N;
    fIdName = new char*[fId5N];
    fIdSens = new Int_t[fId5N];
    for(i=0;i<fId5N;i++) {
	l = strlen(fId5Name[i]);
	fIdName[i] = new char[l+1];
	for(j=0;j<l;j++) fIdName[i][j] = fId5Name[i][j];
	fIdName[i][l] = '\0'; // Null terminate this string.
    } // end for i

    AliITS::Init();
} 
 
//_____________________________________________________________________________
void AliITSv5::StepManager(){
  //
  // Called for every step in the ITS
  //
  Int_t          copy, id;
  Int_t          copy1,copy2;
  Float_t        hits[8];
  Int_t          vol[4];
  TLorentzVector position, momentum;
  TClonesArray   &lhits = *fHits;
  //
  // Track status
  vol[3] = 0;
  if(gMC->IsTrackInside())      vol[3] +=  1;
  if(gMC->IsTrackEntering())    vol[3] +=  2;
  if(gMC->IsTrackExiting())     vol[3] +=  4;
  if(gMC->IsTrackOut())         vol[3] +=  8;
  if(gMC->IsTrackDisappeared()) vol[3] += 16;
  if(gMC->IsTrackStop())        vol[3] += 32;
  if(gMC->IsTrackAlive())       vol[3] += 64;
  //
  // Fill hit structure.
  if(!(gMC->TrackCharge())) return;
  //
  // Only entering charged tracks
  if((id = gMC->CurrentVolID(copy)) == fIdSens[0]) {
      vol[0] = 1;
      id = gMC->CurrentVolOffID(0,copy);
      //detector copy in the ladder = 1<->4  (ITS1)
      vol[1] = copy;
      gMC->CurrentVolOffID(1,copy1);
      //ladder copy in the module   = 1<->2  (I186)
      gMC->CurrentVolOffID(2,copy2);
      //module copy in the layer    = 1<->10 (I132)
      vol[2] = copy1+(copy2-1)*2;//# of ladders in one module  = 2
  } else if(id == fIdSens[1]){
      vol[0] = 2;
      id = gMC->CurrentVolOffID(0,copy);
      //detector copy in the ladder = 1<->4  (ITS2)
      vol[1] = copy;
      gMC->CurrentVolOffID(1,copy1);
      //ladder copy in the module   = 1<->4  (I131)
      gMC->CurrentVolOffID(2,copy2);
      //module copy in the layer    = 1<->10 (I132)
      vol[2] = copy1+(copy2-1)*4;//# of ladders in one module  = 4
  } else if(id == fIdSens[2]){
      vol[0] = 3;
      id = gMC->CurrentVolOffID(1,copy);
      //detector copy in the ladder = 1<->5  (ITS3 is inside I314)
      vol[1] = copy;
      id = gMC->CurrentVolOffID(2,copy);
      //ladder copy in the layer    = 1<->12 (I316)
      vol[2] = copy;
  } else if(id == fIdSens[3]){
      vol[0] = 4;
      id = gMC->CurrentVolOffID(1,copy);
      //detector copy in the ladder = 1<->8  (ITS4 is inside I414)
      vol[1] = copy;
      id = gMC->CurrentVolOffID(2,copy);
      //ladder copy in the layer    = 1<->22 (I417)
      vol[2] = copy;
  }else if(id == fIdSens[4]){
      vol[0] = 5;
      id = gMC->CurrentVolOffID(1,copy);
      //detector copy in the ladder = 1<->23  (ITS5 is inside I562)
      vol[1] = copy;
      id = gMC->CurrentVolOffID(2,copy);
     //ladder copy in the layer    = 1<->34 (I565)
      vol[2] = copy;
  }else if(id == fIdSens[5]){
      vol[0] = 6;
      id = gMC->CurrentVolOffID(1,copy);
      //detector copy in the ladder = 1<->26  (ITS6 is inside I566)
      vol[1] = copy;
      id = gMC->CurrentVolOffID(2,copy);
      //ladder copy in the layer = 1<->38 (I569)
      vol[2] = copy;
  } else {
      return; // not an ITS volume?
  } // end if/else if (gMC->CurentVolID(copy) == fIdSens[i])
//
  gMC->TrackPosition(position);
  gMC->TrackMomentum(momentum);
  hits[0]=position[0];
  hits[1]=position[1];
  hits[2]=position[2];
  hits[3]=momentum[0];
  hits[4]=momentum[1];
  hits[5]=momentum[2];
  hits[6]=gMC->Edep();
  hits[7]=gMC->TrackTime();
  // Fill hit structure with this new hit.
  new(lhits[fNhits++]) AliITShit(fIshunt,gAlice->CurrentTrack(),vol,hits);
  return;
}
//____________________________________________________________________________
void AliITSv5::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliITSv5.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliITS::Streamer(R__b);
      // This information does not need to be read. It is "hard wired"
      // into this class via its creators.
      //R__b >> fId5N;
      //R__b.ReadArray(fId5Name);
   } else {
      R__b.WriteVersion(AliITSv5::IsA());
      AliITS::Streamer(R__b);
      // This information does not need to be saved. It is "hard wired"
      // into this class via its creators.
      //R__b << fId5N;
      //R__b.WriteArray(fId5Name, __COUNTER__);
   }
}

