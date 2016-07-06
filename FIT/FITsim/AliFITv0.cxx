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

/* $Id: AliT0v1.cxx 50646 2011-07-18 14:40:16Z alla $ */

/////////////////////////////////////////////////////////////////////
//                                                                 //
// T0 ( T-zero) detector  version 1                             //
//
//Begin Html       
/*
<img src="gif/AliT0v1Class.gif">
*/
//End Html
//                                                                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <stdlib.h>

#include "TGeoCompositeShape.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoVolume.h"
#include "TGeoTube.h"
#include "TGeoBBox.h"
#include "TGeoNode.h"


#include <TGeoGlobalMagField.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TVirtualMC.h>
#include <TString.h>

#include "AliLog.h"
#include "AliMagF.h"
#include "AliRun.h"

#include "AliFITHits.h"
#include "AliFITv0.h"

#include "AliMC.h"
#include "AliCDBLocal.h"
#include "AliCDBStorage.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliTrackReference.h"

ClassImp(AliFITv0)


//--------------------------------------------------------------------
AliFITv0::AliFITv0():  AliFIT(),
		     fIdSens1(0)

{
  //
  // Standart constructor for T0 Detector version 0
}
//--------------------------------------------------------------------
AliFITv0::AliFITv0(const char *name, const char *title):
  AliFIT(name,title),
  fIdSens1(0)
{
  //
  // Standart constructor for T0 Detector version 0
  //
  printf("@@@@@@@@@@@@@@@ AliFITv0::AliFITv0 \n");
  fIshunt = 2; 
}
//_____________________________________________________________________________

AliFITv0::~AliFITv0() 
{
  // desctructor  
}

//-------------------------------------------------------------------------
void AliFITv0::CreateGeometry()
{
  //
  // Create the geometry of FIT Detector version 1 full geometry
  //
  // begin Html
  //

  printf("@@@@@@@@@@@ AliFITv0::CreateGeometry\n");
  Int_t *idtmed = fIdtmed->GetArray();
  /*
  TGeoMedium* kMedAir = gGeoManager->GetMedium("T0_AIR");
  TGeoMedium* kMedMCPGlass = gGeoManager->GetMedium("T0_glass");
  TGeoMedium* kMedOptGlass = gGeoManager->GetMedium("T0_OpticalGlass");
  TGeoMedium* kMedOptGlassCathode = gGeoManager->GetMedium("T0_OpticalGlassCathode");
  */
  Float_t zdetC = 80;
  Float_t zdetA = 373.;
  Int_t idrotm[999];
  Double_t x,y,z;
  Float_t pstart[3] = {6, 20 ,2.6};
  Float_t pinstart[3] = {3,3,2.55};
  Float_t  pmcp[3] = {2.95, 2.95, 1.5}; //MCP
  Float_t ptop[3] = {1.324, 1.324, 1.};//cherenkov radiator
  Float_t preg[3] = {1.324, 1.324, 0.05};//photcathode 
 
  AliMatrix(idrotm[901], 90., 0., 90., 90., 180., 0.);
  
  //-------------------------------------------------------------------
  //  T0 volume 
  //-------------------------------------------------------------------
  
  Float_t x1[20] = {9,  9, 15 ,15 , 9,  
		    3, -3,  3, -3, -9, 
		    -9, -9, -15, -15, -9, 
		    -3, 3, -3, 3, 9}; 
  
  Float_t y1[20] = {3.2, -3.2, 3.2, -3.2, -9.2, 
		    -9, -9, -15, -15, -9.2,
		    -3.2, 3.2, -3.2, 3.2, 9.2,
		    9, 9, 15, 15, 9.2};
  
  
  //mother tube
   TVirtualMC::GetMC()->Gsvolu("0STR","TUBE",idtmed[kAir],pstart,18);
   TVirtualMC::GetMC()->Gspos("0STR",1,"ALIC",0.,0.,-zdetC-pstart[2],idrotm[901],"ONLY");
   TVirtualMC::GetMC()->Gspos("0STR",2,"ALIC",0.,0.,zdetA+pstart[2],0,"ONLY");
   
   //T0 interior
   TVirtualMC::GetMC()->Gsvolu("0INS","BOX",idtmed[kAir],pinstart,3);
    z=-pstart[2]+pinstart[2];
   for (Int_t is=0; is<20; is++) {
     TVirtualMC::GetMC()->Gspos ("0INS", is + 1, "0STR", x1[is], y1[is], z, 0, "ONLY");
     TVirtualMC::GetMC()->Gspos ("0INS", is + 21, "0STR", x1[is], y1[is], z, 0, "ONLY");
     printf(" 0INS is %i x %f y %f z %f \n",is, x1[is],y1[is], z);
   }
   // 
  x=y=0;
  // Entry window (glass)
  TVirtualMC::GetMC()->Gsvolu("0TOP","BOX",idtmed[kAir],ptop,3); //glass
  TVirtualMC::GetMC()->Gsvolu ("0REG", "BOX", idtmed[kSensAir], preg, 3); 
  TVirtualMC::GetMC()->Gsvolu("0MCP","BOX",idtmed[kGlass],pmcp,3); //glass
  Int_t ntops=0;
  Float_t xin=0, yin=0;
  for (Int_t ix=0; ix<4; ix++) {
    xin = - pinstart[0] + 0.35 + (ix+0.5)*2*ptop[0] ;
    for (Int_t iy=0; iy<4; iy++) {
      z = - pinstart[2]+ptop[2];
      yin = - pinstart[1] + 0.35 + (iy+0.5)*2*ptop[1];
      ntops++;
      TVirtualMC::GetMC()->Gspos("0TOP",ntops,"0INS",xin,yin,z,0,"ONLY");
     //   printf(" 0TOP  full x %f y %f z %f \n", xin, yin, z);
      z = -pinstart[2] + 2 * ptop[2] + preg[2];
      TVirtualMC::GetMC()->Gspos ("0REG",ntops, "0INS", xin, yin, z, 0, "ONLY");
      printf(" GEOGEO  %i %i %i %f %f %f %f %f %f", ntops, ix, iy,
	     xin,yin,x1[ntops],y1[ntops],x1[ntops]+xin,y1[ntops]+yin);
    }
  }
  
 // MCP
 //  TGeoVolume* mcp =  gGeoManager->MakeBox("0MCP",kMedMCPGlass, 2.95, 2.95, 1.5);  
   z=-pinstart[2] + 2*ptop[2] + 2*preg[2] + pmcp[2];
   TVirtualMC::GetMC()->Gspos("0MCP",1,"0INS",0,0,z,0,"ONLY");
    
}    
//------------------------------------------------------------------------
void AliFITv0::AddAlignableVolumes() const
{
  //
  // Create entries for alignable volumes associating the symbolic volume
  // name with the corresponding volume path. Needs to be syncronized with
  // eventual changes in the geometry.
  //
  printf("@@@@@@@@@@@@@@@AliFITv0::AddAlignableVolumes()\n");
  TString volPath;
  TString symName, sn;
  TString vpAalign = "/ALIC_1/0STR_1";
  TString vpCalign = "/ALIC_1/0STR_2";
  for (Int_t imod=0; imod<2; imod++)  {
    if (imod==0) {volPath  = vpCalign; symName="/ALIC_1/0STR_1"; }
    if (imod==1) {volPath  = vpAalign; symName="/ALIC_1/0STR_2"; }
    
    AliDebug(2,"--------------------------------------------");
    AliDebug(2,Form("volPath=%s\n",volPath.Data()));
    AliDebug(2,Form("symName=%s\n",symName.Data()));
    AliDebug(2,"--------------------------------------------");
    if(!gGeoManager->SetAlignableEntry(symName.Data(),volPath.Data()))
      AliFatal(Form("Alignable entry %s not created. Volume path %s not valid",
symName.Data(),volPath.Data()));
  }
}   
//------------------------------------------------------------------------
void AliFITv0::CreateMaterials()
{

  printf("@@@@@@@@@@@@AliFITv0::CreateMaterials\n"); 
   Int_t isxfld   = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ();
   Float_t sxmgmx = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max();
   //   Float_t a,z,d,radl,absl,buf[1];
   // Int_t nbuf;
// AIR
                                                                                
   Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
   Float_t zAir[4]={6.,7.,8.,18.};
   Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
   Float_t dAir = 1.20479E-3;
   Float_t dAir1 = 1.20479E-11;
                                                                               
 // PMT glass SiO2
   Float_t aglass[2]={28.0855,15.9994};
   Float_t zglass[2]={14.,8.};
   Float_t wglass[2]={1.,2.};
   Float_t dglass=2.65;

    
//*** Definition Of avaible T0 materials ***
   AliMixture(1, "Vacuum$", aAir, zAir, dAir1,4,wAir);
   AliMixture(2, "Air$", aAir, zAir, dAir,4,wAir);
   AliMixture( 4, "PMT glass   $",aglass,zglass,dglass,-2,wglass);
   
   AliMedium(1, "FIT_Air$", 2, 0, isxfld, sxmgmx, 10., .1, 1., .003, .003);
   AliMedium(22, "FIT_AirSens$", 2, 1, isxfld, sxmgmx, 10., .1, 1., .003, .003);
   AliMedium(3, "FIT_Vacuum$", 1, 0, isxfld, sxmgmx, 10., .01, .1, .003, .003);
   AliMedium(6, "Glass$", 4, 1, isxfld, sxmgmx, 10., .01, .1, .003, .003);
    
   AliDebugClass(1,": ++++++++++++++Medium set++++++++++");
   
   
}


//-------------------------------------------------------------------
void AliFITv0::Init()
{
  // Initialises version 0 of the Forward Multiplicity Detector
  //
  AliFIT::Init();
  fIdSens1=TVirtualMC::GetMC()->VolId("0REG");
  
  AliDebug(1,Form("%s: *** FIT version 0 initialized ***\n",ClassName()));
}

//-------------------------------------------------------------------

void AliFITv0::StepManager()
{
  //
  // Called for every step in the T0 Detector
  //
  Int_t id,copy,copy1;
  static Float_t hits[6];
  static Int_t vol[3];
  TLorentzVector pos;
  TLorentzVector mom;
  
  //   TClonesArray &lhits = *fHits;
  
  if(!fMC->IsTrackAlive()) return; // particle has disappeared
  
  id=fMC->CurrentVolID(copy);
  // Check the sensetive volume
  if(id==fIdSens1 ) { 
    if(fMC->IsTrackEntering()) {
      fMC->CurrentVolOffID(1,copy1);
      vol[1] = copy1;
      vol[0]=copy;
      fMC->TrackPosition(pos);
      hits[0] = pos[0];
      hits[1] = pos[1];
      hits[2] = pos[2];
      if(pos[2]<0) vol[2] = 0;
      else vol[2] = 1 ;
      printf(" volumes pmt %i mcp %i side %i x %f y %f z %f\n",  vol[0], vol[1],  vol[2], hits[0], hits[1], hits[2] );
      
      Float_t etot=fMC->Etot();
      hits[3]=etot;
      Int_t iPart= fMC->TrackPid();
      Int_t partID=fMC->IdFromPDG(iPart);
      hits[4]=partID;
      Float_t ttime=fMC->TrackTime();
      hits[5]=ttime*1e12;
      AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(),vol,hits);
      // Create a track reference at the exit of photocatode
    }
    
    //charge particle 
    if ( TVirtualMC::GetMC()->TrackCharge() )
      AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber(), AliTrackReference::kFIT);
  } //sensitive

}

