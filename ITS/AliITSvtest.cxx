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
Revision 1.7  2001/02/09 20:06:26  nilsen
Fixed bug in distructor. Can't distroy fixxed length arrays. Thanks Peter.

Revision 1.6  2001/02/09 00:05:31  nilsen
Added fMajor/MinorVersion variables and made other changes to better make
use of the new code changes in AliITSgeom related classes.

Revision 1.5  2001/01/30 09:23:14  hristov
Streamers removed (R.Brun)

Revision 1.4  2001/01/18 06:25:09  barbera
ITS geometry using test Euclid files

Revision 1.1.2.8  2000/10/05 20:28:18  nilsen
Now using root generated streamer function.

Revision 1.1.2.7  2000/07/31 13:51:22  barbera
Updated from the release

Revision 1.2  2000/07/10 16:07:19  fca
Release version of ITS code

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
#include <TGeometry.h>
#include <TNode.h>
#include <TTUBE.h>
#include <TFile.h>    // only required for Tracking function?
#include <TCanvas.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TBRIK.h>
#include <TSystem.h>

#include "AliMC.h"
#include "AliRun.h"
#include "../TGeant3/TGeant3.h"
#include "AliITSGeant3Geometry.h"
#include "AliITShit.h"
#include "AliITS.h"
#include "AliITSvtest.h"
#include "AliITSgeom.h"
#include "AliITSgeomSPD.h"
#include "AliITSgeomSDD.h"
#include "AliITSgeomSSD.h"

ClassImp(AliITSvtest)
 
//_____________________________________________________________________________
AliITSvtest::AliITSvtest() {
    // Standard constructor for the ITS
    Int_t i;

    fIdN    = 0;
    fIdName = 0;
    fIdSens = 0;
    fEuclidOut    = kFALSE; // Don't write Euclide file
    fGeomDetOut   = kFALSE; // Don't write .det file
    fGeomDetIn    = kTRUE; // Read .det file
    fMajorVersion = IsVersion();
    fMinorVersion = -1;
    for(i=0;i<60;i++) fRead[i] = '\0';
    for(i=0;i<60;i++) fWrite[i] = '\0';
    for(i=0;i<60;i++) fEuclidGeomDet[i] = '\0';
}
//____________________________________________________________________________
AliITSvtest::AliITSvtest(const AliITSvtest &source){
////////////////////////////////////////////////////////////////////////
//     Copy Constructor for ITS test version.
////////////////////////////////////////////////////////////////////////
    if(&source == this) return;
    Warning("Copy Constructor","Not allowed to copy AliITSvtest");
    return;
}
//_____________________________________________________________________________
AliITSvtest& AliITSvtest::operator=(const AliITSvtest &source){
////////////////////////////////////////////////////////////////////////
//    Assignment operator for the ITS version 1.
////////////////////////////////////////////////////////////////////////
	if(&source == this) return *this;
	Warning("= operator","Not allowed to copy AliITSvtest");
	return *this;
}
//_____________________________________________________________________________
AliITSvtest::~AliITSvtest() {
    // Standard destructor for the ITS
}
//_____________________________________________________________________________
AliITSvtest::AliITSvtest(const char *fileeuc,const char *filetme,
			 const char *name, const char *title) 
    : AliITS(name, title){
    //
    // Standard constructor for the ITS
    //
    fIdN    = 6;
    fIdName    = new TString[fIdN];
    fIdName[0] = "ITS1";
    fIdName[1] = "ITS2";
    fIdName[2] = "ITS3";
    fIdName[3] = "ITS4";
    fIdName[4] = "ITS5";
    fIdName[5] = "ITS6";
    fIdSens    = new Int_t[fIdN];
    for (Int_t i=0;i<fIdN;i++) fIdSens[i] = 0;
    fMajorVersion = IsVersion();
    fMinorVersion = 1;
    fEuclidOut    = kFALSE; // Don't write Euclide file
    fGeomDetOut   = kFALSE; // Don't write .det file
    fGeomDetIn    = kTRUE; // Read .det file

    fEuclidMaterial = filetme;
    fEuclidGeometry = fileeuc;
    strncpy(fEuclidGeomDet,"$ALICE_ROOT/ITS/ITSgeometry_PPR.det",60);
    strncpy(fRead,fEuclidGeomDet,60);
    strncpy(fWrite,fEuclidGeomDet,60);
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
    cout <<"finished with euclid geometrys"<< endl;
}
//______________________________________________________________________
void AliITSvtest::InitAliITSgeom(){
//     Based on the geometry tree defined in Geant 3.21, this
// routine initilizes the Class AliITSgeom from the Geant 3.21 ITS geometry
// sturture.
    if(!(dynamic_cast<TGeant3*>(gMC))) {
	Error("InitAliITSgeom",
		"Wrong Monte Carlo. InitAliITSgeom uses TGeant3 calls");
	return;
    } // end if
    cout << "Reading Geometry transformation directly from Geant 3." << endl;
    const Int_t nlayers = 6;
    const Int_t ndeep = 9;
    Int_t itsGeomTreeNames[nlayers][ndeep],lnam[20],lnum[20];
    Int_t nlad[nlayers],ndet[nlayers];
    Double_t t[3],r[10];
    Float_t  par[20],att[20];
    Int_t    npar,natt,idshape,imat,imed;
    AliITSGeant3Geometry *ig = new AliITSGeant3Geometry();
    Int_t mod,lay,lad,det,i,j,k;
    char *names[nlayers][ndeep] = {
     {"ALIC","ITSV","ITSD","IT12","I12B","I10B","I107","I101","ITS1"}, // lay=1
     {"ALIC","ITSV","ITSD","IT12","I12B","I20B","I1D7","I1D1","ITS2"}, // lay=2
     {"ALIC","ITSV","ITSD","IT34","I004","I302","ITS3","    ","    "}, // lay=3
     {"ALIC","ITSV","ITSD","IT34","I005","I402","ITS4","    ","    "}, // lay=4
     {"ALIC","ITSV","ITSD","IT56","I565","I562","ITS5","    ","    "}, // lay=5
     {"ALIC","ITSV","ITSD","IT56","I569","I566","ITS6","    ","    "}};// lay=6
    Int_t itsGeomTreeCopys[nlayers][ndeep] = {{1,1,1,1,10, 2, 4, 1, 1},// lay=1
					      {1,1,1,1,10, 4, 4, 1, 1},// lay=2
					      {1,1,1,1,14, 6, 1, 0, 0},// lay=3
					      {1,1,1,1,22, 8, 1, 0, 0},// lay=4
					      {1,1,1,1,34,22, 1, 0, 0},// lay=5
					      {1,1,1,1,38,25, 1, 0, 0}};//lay=6

    // Sorry, but this is not very pritty code. It should be replaced
    // at some point with a version that can search through the geometry
    // tree its self.
    cout << "Reading Geometry informaton from Geant3 common blocks" << endl;
    for(i=0;i<20;i++) lnam[i] = lnum[i] = 0;
    for(i=0;i<nlayers;i++)for(j=0;j<ndeep;j++) 
	itsGeomTreeNames[i][j] = ig->StringToInt(names[i][j]);
    mod = 0;
    for(i=0;i<nlayers;i++){
	k = 1;
	for(j=0;j<ndeep;j++) if(itsGeomTreeCopys[i][j]!=0)
	    k *= TMath::Abs(itsGeomTreeCopys[i][j]);
	mod += k;
    } // end for i

    if(fITSgeom!=0) delete fITSgeom;
    nlad[0]=20;nlad[1]=40;nlad[2]=14;nlad[3]=22;nlad[4]=34;nlad[5]=38;
    ndet[0]=4;ndet[1]=4;ndet[2]=6;ndet[3]=8;ndet[4]=22;ndet[5]=25;
    fITSgeom = new AliITSgeom(0,6,nlad,ndet,mod);
    mod = -1;
    for(lay=1;lay<=nlayers;lay++){
	for(j=0;j<ndeep;j++) lnam[j] = itsGeomTreeNames[lay-1][j];
	for(j=0;j<ndeep;j++) lnum[j] = itsGeomTreeCopys[lay-1][j];
	switch (lay){
	case 1: case 2: // layers 1 and 2 are a bit special
	    lad = 0;
	    for(j=1;j<=itsGeomTreeCopys[lay-1][4];j++){
		lnum[4] = j;
		for(k=1;k<=itsGeomTreeCopys[lay-1][5];k++){
		    lad++;
		    lnum[5] = k;
		    for(det=1;det<=itsGeomTreeCopys[lay-1][6];det++){
			lnum[6] = det;
			mod++;
			ig->GetGeometry(ndeep,lnam,lnum,t,r,idshape,npar,natt,
					par,att,imat,imed);
			fITSgeom->CreatMatrix(mod,lay,lad,det,kSPD,t,r);
			if(!(fITSgeom->IsShapeDefined((Int_t)kSPD)))
			    if(fMinorVersion==1){
                             fITSgeom->ReSetShape(kSPD,
						  new AliITSgeomSPD425Short());
			    } else if(fMinorVersion==2)
                             fITSgeom->ReSetShape(kSPD,
						  new AliITSgeomSPD425Short());
		    } // end for det
		} // end for k
            } // end for j
	    break;
	case 3: case 4: case 5: case 6: // layers 3-6
	    lnum[6] = 1;
	    for(lad=1;lad<=itsGeomTreeCopys[lay-1][4];lad++){
		lnum[4] = lad;
		for(det=1;det<=itsGeomTreeCopys[lay-1][5];det++){
		    lnum[5] = det;
		    mod++;
		    ig->GetGeometry(7,lnam,lnum,t,r,idshape,npar,natt,
				    par,att,imat,imed);
		    switch (lay){
		    case 3: case 4:
			fITSgeom->CreatMatrix(mod,lay,lad,det,kSDD,t,r);
			if(!(fITSgeom->IsShapeDefined(kSDD))) 
			    fITSgeom->ReSetShape(kSDD,new AliITSgeomSDD256());
			    break;
			case 5:
			    fITSgeom->CreatMatrix(mod,lay,lad,det,kSSD,t,r);
			    if(!(fITSgeom->IsShapeDefined(kSSD))) 
				fITSgeom->ReSetShape(kSSD,new AliITSgeomSSD275and75());
			    break;
			case 6:
			    fITSgeom->CreatMatrix(mod,lay,lad,det,kSSDp,t,r);
			    if(!(fITSgeom->IsShapeDefined(kSSDp))) 
				fITSgeom->ReSetShape(kSSDp,new AliITSgeomSSD75and275());
			    break;
			} // end switch
		} // end for det
	    } // end for lad
	    break;
	} // end switch
    } // end for lay
    return;
}
//_____________________________________________________________________________
void AliITSvtest::Init(){
////////////////////////////////////////////////////////////////////////
//     Initialise the ITS after it has been created.
////////////////////////////////////////////////////////////////////////
    Int_t i;

    cout << endl;
    for(i=0;i<29;i++) cout << "*";cout << " ITSvtest_Init ";
    for(i=0;i<28;i++) cout << "*";cout << endl;
//
    if(fRead[0]=='\0') strncpy(fRead,fEuclidGeomDet,60);
    if(fWrite[0]=='\0') strncpy(fWrite,fEuclidGeomDet,60);
    if(fITSgeom!=0) delete fITSgeom;
    fITSgeom = new AliITSgeom();
    if(fGeomDetIn) fITSgeom->ReadNewFile(fRead);
    if(!fGeomDetIn) this->InitAliITSgeom();
    if(fGeomDetOut) fITSgeom->WriteNewFile(fWrite);
    AliITS::Init();
//
    for(i=0;i<72;i++) cout << "*";
    cout << endl;
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

