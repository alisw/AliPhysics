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
Revision 1.17  2000/05/10 19:57:26  nilsen
Fixed problem with display.C and ITS versions v1 and v3.

Revision 1.16  2000/03/05 00:11:03  nilsen
Fixed merge error.

Revision 1.15  2000/02/23 16:25:21  fca
AliVMC and AliGeant3 classes introduced
ReadEuclid moved from AliRun to AliModule

Revision 1.14.4.3  2000/03/04 23:46:38  nilsen
Fixed up the comments/documentation.

Revision 1.14.4.2  2000/03/02 21:53:02  nilsen
To make it compatable with the changes in AliRun/AliModule.

Revision 1.14.4.1  2000/01/12 19:03:33  nilsen
This is the version of the files after the merging done in December 1999.
See the ReadMe110100.txt file for details

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
//
//  Inner Traking System version 5
//  This class contains the base procedures for the Inner Tracking System
//
// Authors: R. Barbera, B. S. Nilsen.
// version 5.
// Created September 17 1999.
//
///////////////////////////////////////////////////////////////////////////////

// See AliITSv5::StepManager().
#define ALIITSPRINTGEOM 0 // default. don't print out gemetry information
//#define ALIITSPRINTGEOM 1 // print out geometry information

#include <stdio.h>
#include <stdlib.h>
#include <TMath.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TTUBE.h>

#include "AliRun.h"
#include "TSystem.h"
#if ALIITSPRINTGEOM==1
#include "../TGeant3/TGeant3.h"
#endif
#include "AliITShit.h"
#include "AliITS.h"
#include "AliITSv5.h"
#include "AliITSgeom.h"

ClassImp(AliITSv5)
 
//_____________________________________________________________________________
AliITSv5::AliITSv5() {
////////////////////////////////////////////////////////////////////////
//    Standard default constructor for the ITS version 5.
////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////
//    Standard destructor for the ITS version 5.
////////////////////////////////////////////////////////////////////////
  delete [] fId5Name;
}
//_____________________________________________________________________________
AliITSv5::AliITSv5(const char *name, const char *title) : AliITS(name, title){
////////////////////////////////////////////////////////////////////////
//    Standard constructor for the ITS version 5.
////////////////////////////////////////////////////////////////////////
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
void AliITSv5::BuildGeometry(){
  //
  // Build ITS TNODE geometry for event display using detailed geometry.
  // This function builds a simple ITS geometry used by the ROOT macro
  // ITSdisplay.C.

  TNode *Top;
  TNode *nd;
  //const int kColorITS_SPD=kRed;
  //const int kColorITS_SDD=kGreen;
  const int kColorITS_SSD=kBlue;
  //
  Top=gAlice->GetGeometry()->GetNode("alice");
  AliITSgeom  *gm = this->GetITSgeom();

  Int_t       lay,lad,det,i;
  Text_t      name[10];
  Float_t     xg[3];
  Float_t     rt[9];
  Double_t    rtd[9];
  TBRIK       *box;
  TRotMatrix  *rm;
  //TCanvas     *c1 = new TCanvas("c1","ITS");

  for(lay=1;lay<=2;lay++)
   for(lad=1;lad<=gm->GetNladders(lay);lad++)
    for(det=1;det<=gm->GetNdetectors(lay);det++){
          try {
              box  = new TBRIK ("ActiveSPD","Active volume of SPD","SPD SI DET",
		        	    0.64,0.0075,4.19); 
          } catch (...) {
	      cout << "EXCEPTION in box = new TBRIK" << endl;
	      return;
	  }
          gm->GetTrans(lay,lad,det,xg[0],xg[1],xg[2]);
          gm->GetRotMatrix(lay,lad,det,rt);
          //sprintf(name,"ROT%1.1d2.2d2.2d",lay,lad,det);
          for(Int_t i=0;i<9;i++) rtd[i] = rt[i];
          try {
	        rm  = new TRotMatrix(name,name,rtd);
          } catch (...) {
	        cout << "EXCEPTION in   new TRotMatrix" << endl;
                return;
          }
         Top->cd();
	  //sprintf(name,"ND%1.1d2.2d2.2d",lay,lad,det); 
         try {
              nd  = new TNode("SPD"," ",box,xg[0],xg[1],xg[2],rm);
         } catch (...) {
              cout << "EXCEPTION in new TNode" << endl;
              return;
         }
         nd->SetLineColor(kColorITS_SSD);
         fNodes->Add(nd);
    }

  for(lay=3;lay<=3;lay++)
   for(lad=1;lad<=gm->GetNladders(lay);lad++)
    for(det=1;det<=gm->GetNdetectors(lay);det++){
          try {
              box  = new TBRIK ("ActiveSDD","Active volume of SDD","SDD SI DET",
		        	    3.5,0.014,3.763); 
          } catch (...) {
	      cout << "EXCEPTION in box = new TBRIK" << endl;
	      return;
	  }
          gm->GetTrans(lay,lad,det,xg[0],xg[1],xg[2]);
          gm->GetRotMatrix(lay,lad,det,rt);
          //sprintf(name,"ROT%1.1d2.2d2.2d",lay,lad,det);
          for(Int_t i=0;i<9;i++) rtd[i] = rt[i];
          try {
	        rm  = new TRotMatrix(name,name,rtd);
          } catch (...) {
	        cout << "EXCEPTION in   new TRotMatrix" << endl;
                return;
          }
         Top->cd();
	  //sprintf(name,"ND%1.1d2.2d2.2d",lay,lad,det); 
         try {
              nd  = new TNode("SDD"," ",box,xg[0],xg[1],xg[2],rm);
         } catch (...) {
              cout << "EXCEPTION in new TNode" << endl;
              return;
         }
         nd->SetLineColor(kColorITS_SSD);
         fNodes->Add(nd);
    }

  for(lay=4;lay<=4;lay++)
   for(lad=1;lad<=gm->GetNladders(lay);lad++)
    for(det=1;det<=gm->GetNdetectors(lay);det++){
          try {
              box  = new TBRIK ("ActiveSDD","Active volume of SDD","SDD SI DET",
		        	    3.5,0.014,3.763); 
          } catch (...) {
	      cout << "EXCEPTION in box = new TBRIK" << endl;
	      return;
	  }
          gm->GetTrans(lay,lad,det,xg[0],xg[1],xg[2]);
          gm->GetRotMatrix(lay,lad,det,rt);
          //sprintf(name,"ROT%1.1d2.2d2.2d",lay,lad,det);
          for(Int_t i=0;i<9;i++) rtd[i] = rt[i];
          try {
	        rm  = new TRotMatrix(name,name,rtd);
          } catch (...) {
	        cout << "EXCEPTION in   new TRotMatrix" << endl;
                return;
          }
         Top->cd();
	  //sprintf(name,"ND%1.1d2.2d2.2d",lay,lad,det); 
         try {
              nd  = new TNode("SDD"," ",box,xg[0],xg[1],xg[2],rm);
         } catch (...) {
              cout << "EXCEPTION in new TNode" << endl;
              return;
         }
         nd->SetLineColor(kColorITS_SSD);
         fNodes->Add(nd);
    }
 for(lay=5;lay<=5;lay++)
   for(lad=1;lad<=gm->GetNladders(lay);lad++)
    for(det=1;det<=gm->GetNdetectors(lay);det++){
          try {
              box  = new TBRIK ("ActiveSSD","Active volume of SSD","SSD SI DET",
		        	    3.65,0.015,2.0); 
          } catch (...) {
	      cout << "EXCEPTION in box = new TBRIK" << endl;
	      return;
	  }
          gm->GetTrans(lay,lad,det,xg[0],xg[1],xg[2]);
          gm->GetRotMatrix(lay,lad,det,rt);
          //sprintf(name,"ROT%1.1d2.2d2.2d",lay,lad,det);
          for(Int_t i=0;i<9;i++) rtd[i] = rt[i];
          try {
	        rm  = new TRotMatrix(name,name,rtd);
          } catch (...) {
	        cout << "EXCEPTION in   new TRotMatrix" << endl;
                return;
          }
         Top->cd();
	  //sprintf(name,"ND%1.1d2.2d2.2d",lay,lad,det); 
         try {
              nd  = new TNode("SSD"," ",box,xg[0],xg[1],xg[2],rm);
         } catch (...) {
              cout << "EXCEPTION in new TNode" << endl;
              return;
         }
         nd->SetLineColor(kColorITS_SSD);
         fNodes->Add(nd);
    }

 for(lay=6;lay<=6;lay++)
   for(lad=1;lad<=gm->GetNladders(lay);lad++)
    for(det=1;det<=gm->GetNdetectors(lay);det++){
          try {
              box  = new TBRIK ("ActiveSSD","Active volume of SSD","SSD SI DET",
		        	    3.65,0.015,2.0); 
          } catch (...) {
	      cout << "EXCEPTION in box = new TBRIK" << endl;
	      return;
	  }

          gm->GetTrans(lay,lad,det,xg[0],xg[1],xg[2]); 
          gm->GetRotMatrix(lay,lad,det,rt);
          //sprintf(name,"ROT%1.1d2.2d2.2d",lay,lad,det);
          for(i=0;i<9;i++) rtd[i] = rt[i];
          try {
	        rm  = new TRotMatrix(name,name,rtd);
          } catch (...) {
	        cout << "EXCEPTION in   new TRotMatrix" << endl;
                return;
          }
         Top->cd();
          //sprintf(name,"ND%1.1d2.2d2.2d",lay,lad,det); 
         try {
              nd  = new TNode("SSD"," ",box,xg[0],xg[1],xg[2],rm);
         } catch (...) {
              cout << "EXCEPTION in new TNode" << endl;
              return;
         }
         nd->SetLineColor(kColorITS_SSD);
         fNodes->Add(nd);
    }


}
//_____________________________________________________________________________
void AliITSv5::CreateMaterials(){
////////////////////////////////////////////////////////////////////////
//     Read a file containing the materials for the ITS version 5.
////////////////////////////////////////////////////////////////////////
    char *filtmp;

  filtmp = gSystem->ExpandPathName(fEuclidMaterial.Data());

  FILE *file = fopen(filtmp,"r");
  if(file) {
    fclose(file);
    ReadEuclidMedia(filtmp);
  } else {
    Error("CreateMaterials"," THE MEDIA FILE %s DOES NOT EXIST !",filtmp);
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
//
//    Read a file containing the geometry for the ITS version 5.
////////////////////////////////////////////////////////////////////////

    Int_t size;
    char topvol[5];
    char *filtmp;

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
  // Place the ITS ghost volume ITSV in its mother volume (ALIC) and make it
  // invisible
  //
  gMC->Gspos("ITSV",1,"ALIC",0,0,0,0,"ONLY");
  //
  // Outputs the geometry tree in the EUCLID/CAD format if requested to do so
  
    if (fEuclidOut) {
      gMC->WriteEuclid("ITSgeometry", "ITSV", 1, 5);
    } // end if (fEuclidOut)

    // read in the file containing the transformations for the active
    // volumes for the ITS version 5. This is expected to be in a file
    // ending in .det. This geometry is kept in the AliITSgeom class.
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
////////////////////////////////////////////////////////////////////////
//     Initialise the ITS after it has been created.
////////////////////////////////////////////////////////////////////////
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
    fMajorVersion = 5;
    fMinorVersion = 0;
} 
//_____________________________________________________________________________
void AliITSv5::StepManager(){
////////////////////////////////////////////////////////////////////////
//    Called for every step in the ITS, then calles the AliITShit class
// creator with the information to be recoreded about that hit.
//     The value of the macro ALIITSPRINTGEOM if set to 1 will allow the
// printing of information to a file which can be used to create a .det
// file read in by the routine CreateGeometry(). If set to 0 or any other
// value except 1, the default behavior, then no such file is created nor
// it the extra variables and the like used in the printing allocated.
////////////////////////////////////////////////////////////////////////
  Int_t          copy, id;
  Int_t          copy1,copy2;
  Float_t        hits[8];
  Int_t          vol[4];
  TLorentzVector position, momentum;
  TClonesArray   &lhits = *fHits;
#if ALIITSPRINTGEOM==1
  FILE          *fp;
  Int_t         i;
  Float_t       xl[3],xt[3],angl[6];
//  Float_t       par[20],att[20];
  Float_t      mat[9];
  static Bool_t first=kTRUE,printit[6][50][50];
  if(first){ for(copy1=0;copy1<6;copy1++)for(copy2=0;copy2<50;copy2++)
      for(id=0;id<50;id++) printit[copy1][copy2][id] = kTRUE;
  first = kFALSE;
  }
  // end if first
#endif
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
#if ALIITSPRINTGEOM==1
  if(printit[vol[0]][vol[2]][vol[1]]){
      printit[vol[0]][vol[2]][vol[1]] = kFALSE;
      xl[0] = xl[1] = xl[2] = 0.0;
      gMC->Gdtom(xl,xt,1);
      for(i=0;i<9;i++) mat[i] = 0.0;
      mat[0] = mat[4] = mat[8] = 1.0;  // default with identity matrix
      xl[0] = 1.0;
      xl[1] = xl[2] =0.0;
      gMC->Gdtom(xl,&(mat[0]),2);
      xl[1] = 1.0;
      xl[0] = xl[2] =0.0;
      gMC->Gdtom(xl,&(mat[3]),2);
      xl[2] = 1.0;
      xl[1] = xl[0] =0.0;
      gMC->Gdtom(xl,&(mat[6]),2);

      angl[0] = TMath::ACos(mat[2]);
      if(mat[2]==1.0) angl[0] = 0.0;
      angl[1] = TMath::ATan2(mat[1],mat[0]);
      if(angl[1]<0.0) angl[1] += 2.0*TMath::Pi();

      angl[2] = TMath::ACos(mat[5]);
      if(mat[5]==1.0) angl[2] = 0.0;
      angl[3] = TMath::ATan2(mat[4],mat[3]);
      if(angl[3]<0.0) angl[3] += 2.0*TMath::Pi();

      angl[4] = TMath::ACos(mat[8]);
      if(mat[8]==1.0) angl[4] = 0.0;
      angl[5] = TMath::ATan2(mat[7],mat[6]);
      if(angl[5]<0.0) angl[5] += 2.0*TMath::Pi();

      for(i=0;i<6;i++) angl[i] *= 180.0/TMath::Pi(); // degrees
      fp = fopen("ITSgeometry_v5.det","a");
      fprintf(fp,"%2d %2d %2d %9e %9e %9e %9e %9e %9e %9e %9e %9e ",
	      vol[0],vol[2],vol[1], // layer ladder detector
	      xt[0],xt[1],xt[2],    // Translation vector
	      angl[0],angl[1],angl[2],angl[3],angl[4],angl[5] // Geant rotaion
                                                           // angles (degrees)
	      );
      fprintf(fp,"%9e %9e %9e %9e %9e %9e %9e %9e %9e",
	     mat[0],mat[1],mat[2],mat[3],mat[4],mat[5],mat[6],mat[7],mat[8]
	  );  // Adding the rotation matrix.
      fprintf(fp,"\n");
      fclose(fp);
  } // end if printit[layer][ladder][detector]
#endif
  return;
}
//____________________________________________________________________________
void AliITSv5::Streamer(TBuffer &R__b){
////////////////////////////////////////////////////////////////////////
//    A dummy Streamer function for this class AliITSv5. By default it
// only streams the AliITS class as it is required. Since this class
// dosen't contain any "real" data to be saved, it doesn't.
////////////////////////////////////////////////////////////////////////
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); 
      if (R__v==1) {
	  AliITS::Streamer(R__b);
	  // This information does not need to be read. It is "hard wired"
	  // into this class via its creators.
	  //R__b >> fId5N;
	  //R__b.ReadArray(fId5Name);
      }else{
      } // end if R__v==1
   } else {
      R__b.WriteVersion(AliITSv5::IsA());
      AliITS::Streamer(R__b);
      // This information does not need to be saved. It is "hard wired"
      // into this class via its creators.
      //R__b << fId5N;
      //R__b.WriteArray(fId5Name, __COUNTER__);
   } // end if R__b.IsReading()
}

