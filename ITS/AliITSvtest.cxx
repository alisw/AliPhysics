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
Revision 1.1.2.2  2000/03/02 21:53:36  nilsen
to make it compatable with the changes in AliRun/AliModule.

Revision 1.1.2.1  2000/01/12 20:19:03  nilsen
	The changes made with this latest inclusion of code is very large.
Many of the new files were added just in December when P. Cerello added his
SDD simulations to the distrobutions. Also added are some file of P. Skowronski
for SSD cluster finding and ghost RecPoints. None of this "new" code has been
proporly tested. Other code new to this cvs repository is explained in the
ITS Off-line web page. In general the changes are too large to give a resonable
discription of them but probably should be taken as the starting point for
the developement branch (ITS-working).
    B. S. Nilsen

Revision 1.13  1999/10/16 19:49:00  BSN
$Name$
$Author$
$Id$
*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Inner Traking System version Test                                        //
//  This class contains the base procedures for the Inner Tracking System    //
//                                                                           //
// Authors: R. Barbera, B. S. Nilsen.                                        //
// version  Test                                                             //
// Created October 16 1999.                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <TMath.h>

#include "AliRun.h"
#include "TSystem.h"
#include "AliITShit.h"
#include "AliITS.h"
#include "AliITSvtest.h"
#include "AliITSgeom.h"

ClassImp(AliITSvtest)
 
//_____________________________________________________________________________
AliITSvtest::AliITSvtest() {
    //
    // Standard constructor for the ITS
    //
/*
    fIdTestN = 6;
    fIdTestName = new TString[fIdTestN];
    fIdTestName[0] = "ITS1";
    fIdTestName[1] = "ITS2";
    fIdTestName[2] = "ITS3";
    fIdTestName[3] = "ITS4";
    fIdTestName[4] = "ITS5";
    fIdTestName[5] = "ITS6";
*/
}
//____________________________________________________________________________
AliITSvtest::AliITSvtest(const AliITSvtest &source){
////////////////////////////////////////////////////////////////////////
//     Copy Constructor for ITS test version.
////////////////////////////////////////////////////////////////////////
    if(&source == this) return;
    printf("Not allowed to copy AliITSvtest\n");
    return;
}
//_____________________________________________________________________________
AliITSvtest& AliITSvtest::operator=(const AliITSvtest &source){
////////////////////////////////////////////////////////////////////////
//    Assignment operator for the ITS version 1.
////////////////////////////////////////////////////////////////////////
	if(&source == this) return *this;
	printf("Not allowed to copy AliITSvtest\n");
	return *this;
}
//_____________________________________________________________________________
AliITSvtest::~AliITSvtest() {
    //
    // Standard destructor for the ITS
    //
  delete [] fIdTestName;
}
//_____________________________________________________________________________
AliITSvtest::AliITSvtest(const char *fileeuc,const char *filetme,
			 const char *name, const char *title) 
    : AliITS(name, title){
    //
    // Standard constructor for the ITS
    //
    fIdTestN = 6;
    fIdTestName = new TString[fIdTestN];
    fIdTestName[0] = "ITS1";
    fIdTestName[1] = "ITS2";
    fIdTestName[2] = "ITS3";
    fIdTestName[3] = "ITS4";
    fIdTestName[4] = "ITS5";
    fIdTestName[5] = "ITS6";

    fEuclidMaterial = filetme;
    fEuclidGeometry = fileeuc;
//  The .det file for the geometry must have the same name as fileeuc with
//  .euc replaced by .det.
}

 
//_____________________________________________________________________________
void AliITSvtest::CreateMaterials(){
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
//    ReadEuclidMedia(fEuclidMaterial.Data(),this);
    ReadEuclidMedia(filtmp);
  } else {
    Error("CreateMaterials"," THE MEDIA FILE %s DOES NOT EXIST !",
//	  fEuclidMaterial.Data());
	  filtmp);
    exit(1);
  } // end if(file)
}

//_____________________________________________________________________________
void AliITSvtest::CreateGeometry(){
//////////////////////////////////////////////////////////////////////
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
void AliITSvtest::Init(){
    //
    // Initialise the ITS after it has been created
    //
    Int_t i;

    fIdN    = fIdTestN;
    fIdName = new TString[fIdTestN];
    fIdSens = new Int_t[fIdTestN];
    for(i=0;i<fIdTestN;i++) {
      fIdName[i]=fIdTestName[i];
      /*
	l = strlen(fIdTestName[i]);
	fIdName[i] = new char[l+1];
	for(j=0;j<l;j++) fIdName[i][j] = fIdTestName[i][j];
	fIdName[i][l] = '\0'; // Null terminate this string.
      */
    } // end for i

    AliITS::Init();
} 
 
//_____________________________________________________________________________
void AliITSvtest::StepManager(){
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
void AliITSvtest::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliITSvtest.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliITS::Streamer(R__b);
   } else {
      R__b.WriteVersion(AliITSvtest::IsA());
      AliITS::Streamer(R__b);
   }
}

