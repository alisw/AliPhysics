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

/* $Id$ */

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

#include <Riostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeoMatrix.h>

#include "AliRun.h"
#include "AliITSgeom.h"
#include "AliITSgeomSDD.h"
#include "AliITSgeomSPD.h"
#include "AliITSgeomSSD.h"
#include "AliITShit.h"
#include "AliITSvtest.h"
#include "AliRun.h"
#include "AliMC.h"

ClassImp(AliITSvtest)
 
//_____________________________________________________________________________
AliITSvtest::AliITSvtest() {
    // Standard constructor for the ITS
    Int_t i;

    fIdN    = 0;
    fIdName = 0;
    fIdSens = 0;
    SetEUCLID(kFALSE);
    fGeomDetOut   = kFALSE; // Don't write .det file
    fGeomDetIn    = kTRUE; // Read .det file
    fMajorVersion = IsVersion();
    fMinorVersion = -1;
    for(i=0;i<60;i++) fRead[i] = '\0';
    for(i=0;i<60;i++) fWrite[i] = '\0';
    for(i=0;i<60;i++) fEuclidGeomDet[i] = '\0';
}
//____________________________________________________________________________
AliITSvtest::AliITSvtest(const AliITSvtest &source) : AliITS(source){
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
    SetEUCLID(kFALSE);
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
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.
    const Int_t knlayers = 6;
    const Int_t kndeep = 3;
    const AliITSDetector idet[knlayers]={kSPD,kSPD,kSDD,kSDD,kSSD,kSSDp};
    const TString names[2][knlayers] = {
     {"/ALIC_1/ITSV_1/ITSD_1/IT12_1/I12A_%d/I10A_%d/I103_%d/I101_1/ITS1_1", // lay=1
      "/ALIC_1/ITSV_1/ITSD_1/IT12_1/I12A_%d/I20A_%d/I1D3_%d/I1D1_1/ITS2_1", // lay=2
      "/ALIC_1/ITSV_1/ITSD_1/IT34_1/I004_%d/I302_%d/ITS3_%d", // lay=3
      "/ALIC_1/ITSV_1/ITSD_1/IT34_1/I005_%d/I402_%d/ITS4_%d", // lay=4
      "/ALIC_1/ITSV_1/ITSD_1/IT56_1/I565_%d/I562_%d/ITS5_%d", // lay=5
      "/ALIC_1/ITSV_1/ITSD_1/IT56_1/I569_%d/I566_%d/ITS6_%d"},// lay=6
     {"/ALIC_1/ITSV_1/ITSD_1/IT12_1/I12B_%d/I10B_%d/I107_%d/I101_1/ITS1_1", // lay=1
      "/ALIC_1/ITSV_1/ITSD_1/IT12_1/I12B_%d/I20B_%d/I1D7_%d/I1D1_1/ITS2_1", // lay=2
      "/ALIC_1/ITSV_1/ITSD_1/IT34_1/I004_%d/I302_%d/ITS3_%d", // lay=3
      "/ALIC_1/ITSV_1/ITSD_1/IT34_1/I005_%d/I402_%d/ITS4_%d", // lay=4
      "/ALIC_1/ITSV_1/ITSD_1/IT56_1/I565_%d/I562_%d/ITS5_%d", // lay=5
      "/ALIC_1/ITSV_1/ITSD_1/IT56_1/I569_%d/I566_%d/ITS6_%d"}
    };
    const Int_t itsGeomTreeCopys[knlayers][kndeep]= {{10, 2, 4},// lay=1
                                                     {10, 4, 4},// lay=2
                                                     {14, 6, 1},// lay=3
                                                     {22, 8, 1},// lay=4
                                                     {34,22, 1},// lay=5
                                                     {38,25, 1}};//lay=6
    Int_t       nlad[knlayers],ndet[knlayers];
    Int_t       mod,lay,lad=0,det=0,i,j,k,cp0,cp1,cp2;
    TString path,shapeName;
    TGeoHMatrix matrix;
    Double_t trans[3]={3*0.0},rot[10]={9*0.0,1.0};
    TArrayD shapePar;
    TArrayF shapeParF;
    Bool_t shapeDefined[4]={kFALSE,kFALSE,kFALSE,kFALSE};

    AliDebug(1,"Reading Geometry transformation directly from Modler.");
    mod = 0;
    for(i=0;i<knlayers;i++){
        k = 1;
        for(j=0;j<kndeep;j++) if(itsGeomTreeCopys[i][j]!=0)
            k *= TMath::Abs(itsGeomTreeCopys[i][j]);
        mod += k;
    } // end for i

    SetITSgeom(0);
    nlad[0]=20;nlad[1]=40;nlad[2]=14;nlad[3]=22;nlad[4]=34;nlad[5]=38;
    ndet[0]= 4;ndet[1]= 4;ndet[2]= 6;ndet[3]= 8;ndet[4]=22;ndet[5]=25;
    AliITSgeom* geom = new AliITSgeom(0,6,nlad,ndet,mod);
    SetITSgeom(geom);
    mod = 0;
    for(lay=1;lay<=knlayers;lay++){
        for(cp0=1;cp0<=itsGeomTreeCopys[lay-1][0];cp0++){
            for(cp1=1;cp1<=itsGeomTreeCopys[lay-1][1];cp1++){
                for(cp2=1;cp2<=itsGeomTreeCopys[lay-1][2];cp2++){
                    path.Form(names[fMinorVersion-1][lay-1].Data(),
                              cp0,cp1,cp2);
                    switch (lay){
                    case 1:{
                        det = cp2;
                        lad = cp1+2*(cp0-1);
                    }break;
                    case 2:{
                        det = cp2;
                        lad = cp1+4*(cp0-1);
                    } break;
                    case 3: case 4: case 5: case 6:{
                        det = cp1;
                        lad = cp0;
                    } break;
                    } // end switch
                         //AliInfo(Form("path=%s lay=%d lad=%d det=%d",
                         //             path.Data(),lay,lad,det));
                    gMC->GetTransformation(path.Data(),matrix);
                    gMC->GetShape(path.Data(),shapeName,shapePar);
                    shapeParF.Set(shapePar.GetSize());
                    for(i=0;i<shapePar.GetSize();i++) shapeParF[i]=shapePar[i];
                    geom->CreateMatrix(mod,lay,lad,det,idet[lay-1],trans,rot);
                    geom->SetTrans(mod,matrix.GetTranslation());
                    geom->SetRotMatrix(mod,matrix.GetRotationMatrix());
		    geom->GetGeomMatrix(mod)->SetPath(path.Data());
                    switch (lay){
                    case 1: case 2:
			if(!shapeDefined[kSPD]){
                        geom->ReSetShape(kSPD,new AliITSgeomSPD425Short(
                                shapeParF.GetSize(),shapeParF.GetArray()));
			shapeDefined[kSPD] = kTRUE;
                    }break;
                    case 3: case 4:
			if(!shapeDefined[kSDD]){
                        geom->ReSetShape(kSDD,new AliITSgeomSDD256(
                                shapeParF.GetSize(),shapeParF.GetArray()));
			shapeDefined[kSDD] = kTRUE;
                    }break;
                    case 5:
			if(!shapeDefined[kSSD]){
                        geom->ReSetShape(kSSD,new AliITSgeomSSD75and275(
                                shapeParF.GetSize(),shapeParF.GetArray()));
			shapeDefined[kSSD] = kTRUE;
                    }break;
                    case 6:
			if(!shapeDefined[kSSDp]){
                        geom->ReSetShape(kSSDp,new AliITSgeomSSD275and75(
                                shapeParF.GetSize(),shapeParF.GetArray()));
			shapeDefined[kSSDp] = kTRUE;
                    }break;
                    default:{
                    }break;
                    } // end switch
                    mod++;
                } /// end for cp2
            } // end for cp1
        } // end for cp0
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
    if(GetITSgeom()!=0) SetITSgeom(0x0);
    SetITSgeom(new AliITSgeom());
    if(fGeomDetIn) GetITSgeom()->ReadNewFile(fRead);
    if(!fGeomDetIn) this->InitAliITSgeom();
    if(fGeomDetOut) GetITSgeom()->WriteNewFile(fWrite);
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
  new(lhits[fNhits++]) AliITShit(fIshunt,gAlice->GetMCApp()->GetCurrentTrackNumber(),vol,hits);
  return;
}

