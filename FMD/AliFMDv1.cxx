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

 /////////////////////////////////////////////////////////////////////
//                                                                 //
// Forward Multiplicity detector based on Silicon version 0        //
//
//Begin Html       
/*
<img src="gif/AliFMDv0Class.gif">
*/
//End Html
//                                                                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <stdlib.h>

#include <TClonesArray.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGeometry.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TNode.h>
#include <TTUBE.h>
#include <TTree.h>
#include <TVirtualMC.h>

#include "AliFMDdigit.h"
#include "AliFMDhit.h"
#include "AliFMDv0.h"
#include "AliFMDv1.h"
#include "AliMagF.h"
#include "AliRun.h"

ClassImp(AliFMDv1)

//--------------------------------------------------------------------
AliFMDv1::AliFMDv1(const char *name, const char *title):
 AliFMD(name,title)
{
  //
  // Standart constructor for Forward Multiplicity Detector version 0
  //
  fIdSens1=0;
  fIdSens2=0;
  fIdSens3=0;
  fIdSens4=0;
  fIdSens5=0;
//  setBufferSize(128000);
}
//-------------------------------------------------------------------------
void AliFMDv1::CreateGeometry()
{
 //
  // Create the geometry of Forward Multiplicity Detector version 0
  //
  //Detector consists of 6 volumes: 
  // 1st covered pseudorapidity interval from 3.3 to 2.0
  // and placed on 65cm in Z-direction;
  // 2nd - from 2.0 to 1.6 and Z=85 cm;
  // 3d  - the same pseudorapidity interval as the 1st 
  // but on the other side from the interaction point z=-65cm;
  // 4th - simmetricaly with the 2nd : 
  // pseudorapidity from 2.0 to 1.6, Z=-85cm   
  // 5th - from 3.6 to 4.7, Z=-270cm
  // 6th - from 4.5 to 5.5 , Z=-630cm.
  // Each part has 400mkm Si (sensetive area, detector itself),
  // 0.75cm of plastic simulated electronics material,
  // Al support ring 2cm thickness and 1cm width placed on 
  // the outer radius of each Si disk;
  //    
  // begin Html
  /*
   <img src="gif/AliFMDv0.gif">
   */
  //

  
  Int_t *idtmed = fIdtmed->GetArray();
   
  Int_t ifmd;
  Int_t idrotm[999];
  Float_t zfmd,par[3];
  char name[5], nameSi[5], nameSector[5], nameRing[5];

  Float_t rin[6], rout[6],zpos;

  Float_t etain[5]= {3.40, 2.29, 3.68, 2.29, 5.09};
  Float_t etaout[6]={2.01, 1.70, 2.28, 1.70, 3.68};
  //  Float_t z[6]={64., 85., -64., -85., -270., -630};
  Float_t z[6]={62.8, 75.2, -83.4, -75.2, -340.};
  Float_t zDet=0.03;
  Float_t zElectronic=0.1;
  Float_t zSupport=1.;

  Float_t zFMD=1.;
//-------------------------------------------------------------------
 //  FMD 
 //------------------------------------------------------------------
	cout<<" !!!!!!!!!!!New FMD geometry !!!!!!!!!"<<endl;

  AliMatrix(idrotm[901], 90, 0, 90, 90, 180, 0);

  //  gMC->Gsvolu("GSI","TUBE", idtmed[1], par, 0);
  gMC->Gsvolu("GEL ","TUBE", idtmed[4], par, 0);
  gMC->Gsvolu("GSUP","TUBE", idtmed[2], par, 0);

  for (ifmd =0; ifmd < 5; ifmd++){

    sprintf(name,"FMD%d",ifmd+1);
    sprintf(nameSi,"GSI%d",ifmd+1);
    sprintf(nameSector,"GSC%d",ifmd+1);
    sprintf(nameRing,"GRN%d",ifmd+1);

    zfmd=TMath::Abs(z[ifmd]);

    AliFMD::Eta2Radius(etain[ifmd],zfmd,&rin[ifmd]);
    AliFMD::Eta2Radius(etaout[ifmd],zfmd,&rout[ifmd]);
    
    par[0]=rin[ifmd]; // pipe size
    par[1]=rout[ifmd];
    par[2]=zFMD/2;
    gMC->Gsvolu(name,"TUBE", idtmed[3], par, 3);
    gMC->Gsvolu(nameSi,"TUBE", idtmed[1], par, 0);
    
    if (z[ifmd] < 0){  
      gMC->Gspos(name,1,"ALIC",0,0,z[ifmd],0, "ONLY");}
    else { 
      gMC->Gspos(name,1,"ALIC",0,0,z[ifmd],idrotm[901], "ONLY");}
  //Silicon detector
    par[2]=zDet/2;
    zpos=zFMD/2 -par[2];
    gMC->Gsposp(nameSi,ifmd+1,name,0,0,zpos,0, "ONLY",par,3);

    //Granularity
    fSectorsSi1=20;
    fRingsSi1=256;
    fSectorsSi2=40;
    fRingsSi2=128;
   if(ifmd==1||ifmd==3)
      { 
	gMC->Gsdvn(nameSector, nameSi , fSectorsSi2, 2);
	gMC->Gsdvn(nameRing, nameSector, fRingsSi2, 1);
      }
    else
      {
	gMC->Gsdvn(nameSector, nameSi , fSectorsSi1, 2);
	gMC->Gsdvn(nameRing, nameSector , fRingsSi1, 1);
      }

    //Plastic slice for electronics
    par[2]=zElectronic/2;
    zpos=zpos-zDet/2-par[2];
    gMC->Gsposp("GEL ",ifmd+1,name,0,0,zpos,0, "ONLY",par,3);

   //Simple Al support
   par[1]=rout[ifmd];
   par[0]=rout[ifmd]-2;
   par[2]=zSupport/2;
   zpos=zpos-zElectronic/2-par[2];
   //   gMC->Gsposp("GSUP",ifmd+1,name,0,0,zpos,0, "ONLY",par,3);
   

  }  

}    

//------------------------------------------------------------------------
void AliFMDv1::CreateMaterials()
{
 Int_t isxfld   = gAlice->Field()->Integ();
 Float_t sxmgmx = gAlice->Field()->Max();

 // Plastic CH
 Float_t aPlastic[2]={1.01,12.01};
 Float_t zPlastic[2]={1,6};
 Float_t wPlastic[2]={1,1};
 Float_t denPlastic=1.03;
   //
 
 //*** Definition Of avaible FMD materials ***
 AliMaterial(0, "Si chip$", 28.0855,14.,2.33,9.36,999);
 AliMaterial(1, "Al supprt$", 26.980,13.,2.70,8.9,999);
 AliMaterial(2, "FMD Air$", 14.61, 7.3, .001205, 30423.,999); 
 AliMixture( 5, "Plastic$",aPlastic,zPlastic,denPlastic,-2,wPlastic);
 

//**
 AliMedium(1, "Si chip$", 0, 1, isxfld, sxmgmx, 1., .001, 1., .001, .001);
 AliMedium(2, "Al support$", 1, 0, isxfld, sxmgmx, 1., .001, 1., .001, .001);
 AliMedium(3, "FMD air$", 2, 0, isxfld, sxmgmx, 1., .001, 1., .001, .001);
 AliMedium(4, "Plastic$", 5, 0,isxfld, sxmgmx,  10., .01, 1., .003, .003);
 


}
//---------------------------------------------------------------------
void AliFMDv1::DrawDetector()
{
//
// Draw a shaded view of the Forward multiplicity detector version 0
//

//Set ALIC mother transparent
gMC->Gsatt("ALIC","SEEN",0);
//
//Set volumes visible
gMC->Gsatt("FMD1","SEEN",1);
gMC->Gsatt("FMD2","SEEN",1);
gMC->Gsatt("FMD3","SEEN",1);
gMC->Gsatt("FMD4","SEEN",1);
gMC->Gsatt("FMD5","SEEN",1);

//
gMC->Gdopt("hide","on");
gMC->Gdopt("shad","on");
gMC->SetClipBox(".");
gMC->SetClipBox("*",0,1000,-1000,1000,-1000,1000);
gMC->DefaultRange();
gMC->Gdraw("alic",40,30,0,12,9.5,.2,0.2);
gMC->Gdhead(1111,"Forward multiplicity detector");
gMC->Gdopt("hide","off");
}
//-------------------------------------------------------------------
void AliFMDv1::Init()
{
// Initialises version 0 of the Forward Multiplicity Detector
//
AliFMD::Init();
fIdSens1=gMC->VolId("GRN1");
fIdSens2=gMC->VolId("GRN2");
fIdSens3=gMC->VolId("GRN3");
fIdSens4=gMC->VolId("GRN4");
fIdSens5=gMC->VolId("GRN5");
printf("*** FMD version 1 initialized ***\n");
}

//-------------------------------------------------------------------

void AliFMDv1::StepManager()
{
  //
  // Called for every step in the Forward Multiplicity Detector
  //
  Int_t id,copy,copy1,copy2;
  static Float_t hits[9];
  static Int_t vol[3];
  static Float_t de;
  TLorentzVector pos;
  TLorentzVector mom;


  TClonesArray &lhits = *fHits;
  if(!gMC->IsTrackAlive()) return; // particle has disappeared

  Float_t charge = gMC->TrackCharge();
  if(TMath::Abs(charge)<=0.) return; //take only charged particles

  //  printf(" in StepManeger \n");
  id=gMC->CurrentVolID(copy);
  //((TGeant3*)gMC)->Gpcxyz();
  
// Check the sensetive volume
   if(id==fIdSens1||id==fIdSens2||id==fIdSens3||id==fIdSens4||id==fIdSens5)
     {
       if(gMC->IsTrackEntering())
	 {
	   vol[2]=copy;
	   gMC->CurrentVolOffID(1,copy1);
	   vol[1]=copy1;
	   gMC->CurrentVolOffID(2,copy2);
	   vol[0]=copy2;

	   gMC->TrackPosition(pos);
	   hits[0]=pos[0];
	   hits[1]=pos[1];
	   hits[2]=pos[2];

	   gMC->TrackMomentum(mom);
	   hits[3]=mom[0];
	   hits[4]=mom[1];
	   hits[5]=mom[2];

	   Int_t iPart= gMC->TrackPid();
	   Int_t partId=gMC->IdFromPDG(iPart);
	   hits[7]=partId;
	   hits[8]=1e9*gMC->TrackTime();
	   de=0.;
	 }
       if(gMC->IsTrackInside()){
	   de=de+1000.*gMC->Edep();
       }
       
       if(gMC->IsTrackExiting()
	  ||gMC->IsTrackDisappeared()||
	  gMC->IsTrackStop())
	 {
	     hits[6]=de+1000.*gMC->Edep();
      new(lhits[fNhits++]) AliFMDhit(fIshunt,gAlice->GetCurrentTrackNumber(),vol,hits);
	 } // IsTrackExiting()
     }
  }
//--------------------------------------------------------------------------

void AliFMDv1::Response( Float_t Edep)
{
  Float_t I=1.664*0.04*2.33/22400; // = 0.69e-6;
  Float_t chargeOnly=Edep/I;
  //Add noise ~500electrons
  Int_t charge=500;
  if (Edep>0)
     charge=Int_t(gRandom->Gaus(chargeOnly,500));	
 }   






