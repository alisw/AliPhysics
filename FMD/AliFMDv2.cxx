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
#include <TMath.h>
#include <TGeometry.h>
#include <TTUBE.h>
#include <TFile.h>
#include <TTree.h>
#include <TNode.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TDirectory.h>
#include "AliFMDv2.h"
#include "AliFMDv0.h"
#include "AliRun.h"
#include <Riostream.h>
#include "AliMagF.h"
#include "AliFMDhit.h"
#include "AliFMDdigit.h"
#include <stdlib.h>

ClassImp(AliFMDv2)

//--------------------------------------------------------------------
AliFMDv2::AliFMDv2(const char *name, const char *title):
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
void AliFMDv2::CreateGeometry()
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
  Float_t zFMD,par[3],ppcon[5];
    //ppcon1[15],ppcon2[15],ppcon3[15],ppcon4[15],ppcon5[15];
  Float_t z[5]={62.8, 75.2, -83.4, -75.2, -340.};
  Float_t NylonTube[3]={0.2,0.6,0.45};
  Float_t zPCB=0.12; Float_t zHoneyComb=0.5; 
  Float_t zSi=0.03;
 
//-------------------------------------------------------------------
 //  FMD 
 //------------------------------------------------------------------
  //  cout<<" !!!!!!!!!!!New FMD geometry !!!!!!!!!"<<endl;
  

  char nameFMD[5], nameSi[5], nameSector[5], nameRing[5];
  Char_t nameHoney[5], nameHoneyIn[5], nameHoneyOut[5];
  Char_t namePCB[5], nameCopper[5], nameChips[5], nameG10[5];
  Char_t nameLPCB[5], nameLCopper[5], nameLChips[5], nameGL10[5];;
  Float_t rin[5]={4.2,15.4,4.2,15.4,4.2};
  Float_t rout[5]={17.4,28.4,17.4,28.4,17.4};
  Float_t RinHoneyComb[5] ={ 5.15,16.4,  5.15,16.4,  5.15};
  Float_t RoutHoneyComb[5]={20.63,34.92,22.3, 32.02,20.63};
  Float_t zInside;
  Float_t zCooper=0.01; Float_t zChips=0.01;
  Float_t yNylonTube[5]={10,20,10,20,10};


  AliMatrix(idrotm[901], 90, 0, 90, 90, 180, 0);
  
  
  // Nylon tubes
   gMC->Gsvolu("GNYL","TUBE", idtmed[1], NylonTube, 3);  //support nylon tube
   Float_t wideSupport=zSi+3*zPCB+2*NylonTube[2]+zHoneyComb;
   //   cout<<" wideSupport "<<wideSupport<<endl;

  for (ifmd=0; ifmd<5; ifmd++)
    {
      sprintf(nameFMD,"FMD%d",ifmd+1);
      ppcon[0]=0;
      ppcon[1]=360;
      ppcon[2]=4;
      
      //      ppcon[3]=z[ifmd];
      ppcon[3]=-wideSupport/2.;
      ppcon[4]=rin[ifmd];
      ppcon[5]=rout[ifmd];
      
      ppcon[6]=ppcon[3]+zSi+2*zPCB+2*NylonTube[2];
      ppcon[7]=rin[ifmd];
      ppcon[8]=rout[ifmd];
      
      ppcon[9]=ppcon[6];
      ppcon[10]=RinHoneyComb[ifmd];
      ppcon[11]=RoutHoneyComb[ifmd];

      ppcon[12]=ppcon[9]+zHoneyComb+zPCB;
      ppcon[13]=RinHoneyComb[ifmd];
      ppcon[14]=RoutHoneyComb[ifmd];
      zFMD=z[ifmd]+wideSupport/2.;
      //      cout<<" Si "<<ifmd+1<<" "<<zFMD<<" "<<ppcon[12]<<endl;
      gMC->Gsvolu(nameFMD,"PCON",idtmed[0],ppcon,15);
      if (z[ifmd] >0){  
	gMC->Gspos(nameFMD,1,"ALIC",0,0,z[ifmd],0, "ONLY");}
      else { 
	gMC->Gspos(nameFMD,1,"ALIC",0,0,z[ifmd],idrotm[901], "ONLY");}
     //silicon
      sprintf(nameSi,"GSI%d",ifmd+1);
      sprintf(nameSector,"GSC%d",ifmd+1);
      sprintf(nameRing,"GRN%d",ifmd+1);
      //honeycomb support
      sprintf(nameHoney,"GSU%d",ifmd+1);
      gMC->Gsvolu(nameHoney,"TUBE", idtmed[0], par, 0);  //honeycomb 
      sprintf(nameHoneyIn,"GHI%d",ifmd+1);
      gMC->Gsvolu(nameHoneyIn,"TUBE", idtmed[7], par, 0);  //honey comb inside 
      sprintf(nameHoneyOut,"GHO%d",ifmd+1);
      gMC->Gsvolu(nameHoneyOut,"TUBE", idtmed[6], par, 0);  //honey comb skin
      //PCB
      sprintf(namePCB,"GPC%d",ifmd+1);
      gMC->Gsvolu(namePCB,"TUBE", idtmed[0], par, 0); //PCB
      sprintf(nameCopper,"GCO%d",ifmd+1);
      gMC->Gsvolu(nameCopper,"TUBE", idtmed[3], par, 0);  // Cooper
      sprintf(nameChips,"GCH%d",ifmd+1);
      gMC->Gsvolu(nameChips,"TUBE", idtmed[5], par, 0); // Si chips
      sprintf(nameG10,"G10%d",ifmd+1);
      gMC->Gsvolu(nameG10,"TUBE", idtmed[2], par, 0);  //G10 plate
      //last PCB
      sprintf(nameLPCB,"GPL%d",ifmd+1);
      gMC->Gsvolu(nameLPCB,"TUBE", idtmed[0], par, 0); //PCB
      sprintf(nameLCopper,"GCL%d",ifmd+1);
      gMC->Gsvolu(nameLCopper,"TUBE", idtmed[3], par, 0);  // Cooper
      sprintf(nameLChips,"GHL%d",ifmd+1);
      gMC->Gsvolu(nameLChips,"TUBE", idtmed[5], par, 0); // Si chips
      sprintf(nameGL10,"G1L%d",ifmd+1);
      gMC->Gsvolu(nameGL10,"TUBE", idtmed[2], par, 0);  //G10 plate
     
     
      par[0]=rin[ifmd]; // pipe size
      par[1]=rout[ifmd];
      par[2]=zSi/2;
      gMC->Gsvolu(nameSi,"TUBE", idtmed[4], par, 3);
      printf ("rin %f rout %f ZFMD %f\n",par[0],par[1],z[ifmd]);
      zInside=ppcon[3]+par[2];
      //      cout<<" zInside "<<zInside<<endl;
      gMC->Gspos(nameSi,ifmd+1,nameFMD,0,0,zInside,0, "ONLY");
      //PCB 1
      zInside += par[2]+zPCB/2;
      par[2]=zPCB/2;
      //   cout<<"FMD "<<ifmd<<" zInside PCB1 = "<<zInside<<endl;
      gMC->Gsposp(namePCB,1,nameFMD,0,0,zInside,0, "ONLY",par,3);
      zInside += zPCB;
      //  cout<<"FMD "<<ifmd<<" zInside PCB2 = "<<zInside<<endl;
      gMC->Gsposp(namePCB,2,nameFMD,0,0,zInside,0, "ONLY",par,3);
      Float_t NulonTubeBegin=zInside+2*zPCB;
      par[2]=zPCB/2-0.02;
      Float_t zInPCB = -zPCB/2+par[2];
      gMC->Gsposp(nameG10,1,namePCB,0,0,zInPCB,0, "ONLY",par,3);
      zInPCB+=par[2]+zCooper/2 ;
      par[2]=zCooper/2;
      gMC->Gsposp(nameCopper,1,namePCB,0,0,zInPCB,0, "ONLY",par,3);
      zInPCB += zCooper/2 + zChips/2;
      par[2]=zChips/2;
      gMC->Gsposp(nameChips,1,namePCB,0,0,zInPCB,0, "ONLY",par,3);
      //HoneyComb
      zHoneyComb=0.8;   
      par[0] = RinHoneyComb[ifmd];
      par[1] = RoutHoneyComb[ifmd];
      par[2] = zHoneyComb/2;
      //  cout<<"FMD "<<ifmd<<" zInside before HoneyComb = "<<zInside<<" par[2] "<<par[2]<<endl;
      zInside += 2*NylonTube[2]+par[2];
      // cout<<" zInside HoneyComb"<<zInside<<endl;
      gMC->Gsposp(nameHoney,1,nameFMD,0,0,zInside,0, "ONLY",par,3);
      par[2]=0.1/2;
      Float_t zHoney=-zHoneyComb/2+par[2];
      //  cout<<" zHoney shurka 1 "<<zHoney <<endl;
      gMC->Gsposp(nameHoneyOut,1,nameHoney,0,0,zHoney,0,
		  "ONLY",par,3); //shkurki
      zHoney=zHoneyComb/2-par[2];
      //  cout<<" zHoney shurka 2 "<<zHoney <<endl;
      gMC->Gsposp(nameHoneyOut,2,nameHoney,0,0,zHoney,0, "ONLY",par,3);
      par[2]=(zHoneyComb-2.*0.1)/2; //soty vnutri
      gMC->Gsposp(nameHoneyIn,1,nameHoney,0,0,0,0, "ONLY",par,3);
      
      gMC->Gspos("GNYL",1,nameFMD,0,yNylonTube[ifmd],
		 NulonTubeBegin+NylonTube[2]/2.,0, "ONLY");
      gMC->Gspos("GNYL",2,nameFMD,0,-yNylonTube[ifmd],
		 NulonTubeBegin+NylonTube[2]/2.,0, "ONLY");
         
      //last PCB
      par[0]=RoutHoneyComb[ifmd]-9;
      par[1]=RoutHoneyComb[ifmd];
      par[2]=zPCB/2;
      //  cout<<" rin "<<par[0]<<" rout "<<par[1]<<endl;
      zInside += zHoneyComb/2+par[2];
      gMC->Gsposp(nameLPCB,1,nameFMD,0,0,zInside,0, "ONLY",par,3);
      
       par[2]=zPCB/2-0.02;
       zInPCB = -zPCB/2+par[2];
       gMC->Gsposp(nameGL10,1,nameLPCB,0,0,zInPCB,0, "ONLY",par,3);
       zInPCB+=par[2]+zCooper/2 ;
       par[2]=zCooper/2;
       gMC->Gsposp(nameLCopper,1,nameLPCB,0,0,zInPCB,0, "ONLY",par,3);
       zInPCB += zCooper/2 + zChips/2;
       par[2]=zChips/2;
       gMC->Gsposp(nameLChips,1,nameLPCB,0,0,zInPCB,0, "ONLY",par,3);
      
       /*
      if (ifmd==3) {
 	par[0]=0.1;
	par[1]=7.;
	par[2]=0.1;
	gMC->Gsvolu("GST1","BOX ", idtmed[6], par, 3);  //support carbon stick
	gMC->Gspos("GST1",1,nameFMD,0.,49.2-par[1],z[3]-0.7,0,"ONLY");
	gMC->Gspos("GST1",2,nameFMD,0.,-49.2+par[1],z[3]-0.7,0,"ONLY");
	
	Float_t pSupportCone[5];
	pSupportCone[0]=1.5;
	pSupportCone[1]=34.4;
	pSupportCone[2]=34.5;
	pSupportCone[3]=30.;
	pSupportCone[4]=29.9;
	gMC->Gsvolu("GCY1","CONE", idtmed[6],pSupportCone , 5);  //support carbon stick
	gMC->Gspos("GCY1",1,nameFMD,0.,0.,z[3]+1.5,0,"ONLY");
	
      }
     if (ifmd==2) {
 	par[0]=0.1;
	par[1]=1.3;
	par[2]=0.1;
	gMC->Gsvolu("GST3","BOX ", idtmed[6], par, 3);  //support carbon stick
	gMC->Gspos("GST3",1,nameFMD,0.,23.-par[1],z[2]+2.1,0,"ONLY");
	gMC->Gspos("GST3",2,nameFMD,0.,-23.+par[1],z[2]+2.1,0,"ONLY");

	Float_t pSupportCone[5];
	pSupportCone[0]=1.5;
	pSupportCone[1]=30.;
	pSupportCone[2]=29.9;
	pSupportCone[3]=21.;
	pSupportCone[4]=20.9;
	gMC->Gsvolu("GCY2","CONE", idtmed[6],pSupportCone , 5);  //support carbon stick
	gMC->Gspos("GCY2",1,nameFMD,0.,0.,z[2]+1.5,0,"ONLY");
	
     }
       */
       
     //Granularity
    fSectorsSi1=20;
    fRingsSi1=256*3;
    //fRingsSi1=3; // for drawing only
    fSectorsSi2=40;
    fRingsSi2=128*3;
    // fRingsSi2=3; //for  drawing onl
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
    
    }
}    


//------------------------------------------------------------------------


void AliFMDv2::CreateMaterials()
{
 Int_t isxfld   = gAlice->Field()->Integ();
 Float_t sxmgmx = gAlice->Field()->Max();

 // Plastic CH
 Float_t aPlastic[2]={1.01,12.01};
 Float_t zPlastic[2]={1,6};
 Float_t wPlastic[2]={1,1};
 Float_t denPlastic=1.03;
 //     60% SiO2 , 40% G10FR4 
 // PC board
 Float_t apcb[3]  = { 28.0855,15.9994,17.749 };
 Float_t zpcb[3]  = { 14.,8.,8.875 };
 Float_t wpcb[3]  = { .28,.32,.4 };
 Float_t denspcb  = 1.8;
   //
 
 //*** Definition Of avaible FMD materials ***
 AliMaterial(0, "FMD Air$", 14.61, 7.3, .001205, 30423.,999); 
 AliMixture(1, "Plastic$",aPlastic,zPlastic,denPlastic,-2,wPlastic);
 AliMixture(2, "SSD PCB$",   apcb, zpcb, denspcb, 3, wpcb);
 AliMaterial(3, "SSD Copper$", 63.546, 29., 8.96, 1.43, 999.);
 AliMaterial(4, "SSD Si$",      28.0855, 14., 2.33, 9.36, 999.);
 AliMaterial(5, "SSD Si chip$", 28.0855, 14., 2.33, 9.36, 999.);
 AliMaterial(6, "SSD C$",       12.011,   6., 2.265,18.8, 999.);
 AliMaterial(7, "SSD Kapton$", 12.011, 6., 0.01, 31.27, 999.);//honeycomb
  AliMaterial(8, "SSD G10FR4$", 17.749, 8.875, 1.8, 21.822, 999.);
   

//**
 AliMedium(0, "FMD air$", 0, 0, isxfld, sxmgmx, 1., .001, 1., .001, .001);
 AliMedium(1, "Plastic$", 1, 0,isxfld, sxmgmx,  10., .01, 1., .003, .003);
 AliMedium(2, "SSD PCB$", 2, 0, isxfld, sxmgmx, 1., .001, 1., .001, .001);
 AliMedium(3, "SSD Copper$", 3, 0,isxfld, sxmgmx,  10., .01, 1., .003, .003);
 AliMedium(4, "SSD Si$", 4, 1, isxfld, sxmgmx, 1., .001, 1., .001, .001);
 AliMedium(5, "SSD Si chip$", 5, 0,isxfld, sxmgmx,  10., .01, 1., .003, .003);
 AliMedium(6, "SSD C$", 6, 0,isxfld, sxmgmx,  10., .01, 1., .003, .003);
 AliMedium(7, "SSD Kapton$", 7, 0, isxfld, sxmgmx, 1., .001, 1., .001, .001);
 AliMedium(8, "SSD G10FR4$", 8, 0,isxfld, sxmgmx,  10., .01, 1., .003, .003);
 


}

//---------------------------------------------------------------------
void AliFMDv2::DrawDetector()
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
void AliFMDv2::Init()
{
// Initialises version 0 of the Forward Multiplicity Detector
//
AliFMD::Init();
fIdSens1=gMC->VolId("GRN1");
fIdSens2=gMC->VolId("GRN2");
fIdSens3=gMC->VolId("GRN3");
fIdSens4=gMC->VolId("GRN4");
fIdSens5=gMC->VolId("GRN5");
printf("*** FMD version 2 initialized ***\n");
}

//-------------------------------------------------------------------

void AliFMDv2::StepManager()
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
      new(lhits[fNhits++]) AliFMDhit(fIshunt,gAlice->CurrentTrack(),vol,hits);
	 } // IsTrackExiting()
     }
  }

//--------------------------------------------------------------------------

void AliFMDv2::Response( Float_t Edep)
{
  Float_t I=1.664*0.04*2.33/22400; // = 0.69e-6;
  Float_t chargeOnly=Edep/I;
  //Add noise ~500electrons
  Int_t charge=500;
  if (Edep>0)
     charge=Int_t(gRandom->Gaus(chargeOnly,500));	
 }   







