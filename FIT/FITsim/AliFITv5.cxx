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
// FIT detector full geometry  version 5                           //
//
//Begin Html       
/*
<img src="gif/AliFITv5Class.gif">
*/
//End Html
//
// V0+ part by Lizardo Valencia Palomo    lvalenci@cern.ch          //
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
#include "TGeoArb8.h"

#include <TGeoGlobalMagField.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TVirtualMC.h>
#include <TString.h>
#include <TParticle.h>

#include "AliLog.h"
#include "AliMagF.h"
#include "AliRun.h"

#include "AliFITHits.h"
#include "AliFITv5.h"

#include "AliMC.h"
#include "AliCDBLocal.h"
#include "AliCDBStorage.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliTrackReference.h"

ClassImp(AliFITv5)


//--------------------------------------------------------------------
AliFITv5::AliFITv5():  AliFIT(),
  fIdSens1(0),
  fIdSens2(0),
  fPMTeff(0x0),
//V0+
  nSectors(8),
  nRings(5),
  fSenseless(-1),
  fV0PlusR0(3.97),//Computed for z = 325 cm from IP
  fV0PlusR1(7.6),//From V0A 
  fV0PlusR2(13.8),//From V0A 
  fV0PlusR3(22.7),//From V0A
  fV0PlusR4(41.3),//From V0A
  fV0PlusR5(72.94),//Computed for z = 325 cm from IP
  fV0PlusR6(72.6),//Needed to compute fV0PlusnMeters
  fV0PlusSciWd(2.5),//From V0A
  fV0PlusFraWd(0.2),//From V0A
  fV0PlusZposition(+325),//Must be changed to specifications from Corrado (meeting 10/11/176)
  fV0PlusnMeters(fV0PlusR6*0.01),//From V0A
  fV0PlusLightYield(93.75),//From V0A
  fV0PlusLightAttenuation(0.05),//From V0A 
  fV0PlusFibToPhot(0.3)//From V0A

{
  //
  // Standart constructor for T0 Detector version 0

  for (Int_t i = 0; i < nSectors; i++)
    {
      for (Int_t j = 0; j < nRings; j++) fIdV0Plus[i][j] = 0;
    }

}
//--------------------------------------------------------------------
AliFITv5::AliFITv5(const char *name, const char *title):
  AliFIT(name,title),
  fIdSens1(0),
  fIdSens2(0),
  fPMTeff(0x0),
//V0+
  nSectors(8),
  nRings(5),
  fSenseless(-1),
  fV0PlusR0(3.97),//Computed for z = 325 cm from IP
  fV0PlusR1(7.6),//From V0A 
  fV0PlusR2(13.8),//From V0A 
  fV0PlusR3(22.7),//From V0A
  fV0PlusR4(41.3),//From V0A
  fV0PlusR5(72.94),//Computed for z = 325 cm from IP
  fV0PlusR6(72.6),//Needed to compute fV0PlusnMeters
  fV0PlusSciWd(2.5),//From V0A
  fV0PlusFraWd(0.2),//From V0A
  fV0PlusZposition(+325),//Must be changed to specifications from Corrado (meeting 10/11/176)
  fV0PlusnMeters(fV0PlusR6*0.01),//From V0A
  fV0PlusLightYield(93.75),//From V0A
  fV0PlusLightAttenuation(0.05),//From V0A 
  fV0PlusFibToPhot(0.3)//From V0A

{
  //
  // Standart constructor for FIT Detector version 0
  //
  fIshunt = 2; 
  SetPMTeff();
}
//_____________________________________________________________________________

AliFITv5::~AliFITv5() 
{
  // desctructor  
}

//-------------------------------------------------------------------------
void AliFITv5::CreateGeometry()
{
  //
  // Create the geometry of FIT Detector version 1 full geometry
  //
  // begin Html
  //

  Int_t *idtmed = fIdtmed->GetArray();
  Float_t zdetC = 85; //center of mother volume
  Float_t zdetA = 333;
  
  Int_t idrotm[999];
  Double_t x,y,z;
  Float_t pstartC[3] = {6., 20 ,5};
  Float_t pstartA[3] = {2.55, 20 ,5};
  Float_t pinstart[3] = {2.95,2.95,4.};
  Float_t  pmcp[3] = {2.949, 2.949, 2.8}; //MCP
  Float_t ptop[3] = {1.324, 1.324, 1.};//cherenkov radiator
  Float_t preg[3] = {1.324, 1.324, 0.05};//photcathode 
 
  Float_t zV0A = 325.;
  Float_t zV0C = 89;
  Float_t pV0Amother[3] = {4.25, 41.25, 0.15};
  Float_t pV0A[3] = {4.3, 41.2, 0.1};

  AliMatrix(idrotm[901], 90., 0., 90., 90., 180., 0.);
  
  //-------------------------------------------------------------------
  //  T0 volume 
  //-------------------------------------------------------------------
  //C side
  /*
  Float_t xc[28] = {-8.95, -3.05, 3.05, 8.95, -14.85,
		    -8.95, -3.05, 3.05, 8.95,  14.85,
		    -14.85, -8.95, 8.95, 14.85, -14.85, 
		    -8.95, 8.95, 14.85, -14.85,-8.95,
		    -3.05, 3.05, 8.95, 14.85, -8.95,
		    -3.05, 3.05, 8.95 };
  
  Float_t yc[28] = {14.85, 14.85,14.85,14.85,8.95,
		    8.95,8.95,8.95,8.95,8.95, 
		    3.05, 3.05,3.05,3.05,-3.05,
		    -3.05,-3.05,-3.05,-8.95,-8.95,
		    -8.95,-8.95,-8.95,-8.95,-14.85, 
		    -14.85,-14.85,-14.85};
  */ 
  //Alpha, Beta, and Gamma angles for the C side Rotations
 
  Double_t ac[28] = { 211.0555505 ,  191.37002946,  348.62997054,  328.9444495 ,
         238.9444495 ,  225.        ,  198.46612083,  341.53387917,
         315.        ,  301.0555505 ,  258.62997054,  251.53387917,
         288.46612083,  281.37002946,  281.37002946,  288.46612083,
         251.53387917,  258.62997054,  301.0555505 ,  315.        ,
         341.53387917,  198.46612083,  225.        ,  238.9444495 ,
         328.9444495 ,  348.62997054,  191.37002946,  211.0555505 };
   Double_t bc[28] = {12.35160924,  10.77300323, -10.77300323, -12.35160924,
         12.35160924,   8.978165  ,   6.68091565,  -6.68091565,
         -8.978165  , -12.35160924,  10.77300323,   6.68091565,
         -6.68091565, -10.77300323,  10.77300323,   6.68091565,
         -6.68091565, -10.77300323,  12.35160924,   8.978165  ,
          6.68091565,  -6.68091565,  -8.978165  , -12.35160924,
         12.35160924,  10.77300323, -10.77300323, -12.35160924};
   Double_t gc[28] = {-211.0555505 , -191.37002946, -348.62997054, -328.9444495 ,
        -238.9444495 , -225.        , -198.46612083, -341.53387917,
        -315.        , -301.0555505 , -258.62997054, -251.53387917,
        -288.46612083, -281.37002946, -281.37002946, -288.46612083,
        -251.53387917, -258.62997054, -301.0555505 , -315.        ,
        -341.53387917, -198.46612083, -225.        , -238.9444495 ,
        -328.9444495 , -348.62997054, -191.37002946, -211.0555505 };
   Double_t xc2[28] = { -9.10385087,  -3.04012128,   3.04012128,   9.10385087,
        -15.11813129,  -9.10385087,  -3.04012128,   3.04012128,
          9.10385087,  15.11813129, -15.11813129,  -9.10385087,
          9.10385087,  15.11813129, -15.11813129,  -9.10385087,
          9.10385087,  15.11813129, -15.11813129,  -9.10385087,
         -3.04012128,   3.04012128,   9.10385087,  15.11813129,
         -9.10385087,  -3.04012128,   3.04012128,   9.10385087};
   Double_t yc2[28] = { 15.11813129,  15.11813129,  15.11813129,  15.11813129,
          9.10385087,   9.10385087,   9.10385087,   9.10385087,
          9.10385087,   9.10385087,   3.04012128,   3.04012128,
          3.04012128,   3.04012128,  -3.04012128,  -3.04012128,
         -3.04012128,  -3.04012128,  -9.10385087,  -9.10385087,
         -9.10385087,  -9.10385087,  -9.10385087,  -9.10385087,
        -15.11813129, -15.11813129, -15.11813129, -15.11813129};
   Double_t zc2[28] = {80.59039648,  81.04597318,  81.04597318,  80.59039648,
         80.59039648,  81.4892005 ,  81.93978009,  81.93978009,
         81.4892005 ,  80.59039648,  81.04597318,  81.93978009,
         81.93978009,  81.04597318,  81.04597318,  81.93978009,
         81.93978009,  81.04597318,  80.59039648,  81.4892005 ,
         81.93978009,  81.93978009,  81.4892005 ,  80.59039648,
         80.59039648,  81.04597318,  81.04597318,  80.59039648};
  //A Side
  /*
  Float_t xa[24] = {-11.8, -5.9, 0, 5.9, 11.8, 
		    -11.8, -5.9, 0, 5.9, 11.8,
		    -11.8, -5.9, 5.9, 11.8, 
		    -11.8, -5.9, 0, 5.9, 11.8,
		    -11.8, -5.9, 0, 5.9, 11.8 };
  
  Float_t ya[24] = { 11.9, 11.9, 11.9, 11.9, 11.9,
		     6.0,   6.0,  6.0, 6.0,  6.0,
		    -0.1, -0.1, 0.1, 0.1,
		    -6.0, -6.0, -6.0, -6.0, -6.0,
		     -11.9, -11.9, -11.9,  -11.9, -11.9 };  
  */  
 Float_t xa[24] = {-11.8, -5.9, 0, 5.9, 11.8, 
		    -11.8, -5.9, 0, 5.9, 11.8,
		    -12.8, -6.9, 6.9, 12.8, 
		    -11.8, -5.9, 0, 5.9, 11.8,
		    -11.8, -5.9, 0, 5.9, 11.8 };
  
  Float_t ya[24] = { 11.9, 11.9, 12.9, 11.9, 11.9,
		     6.0,   6.0,  7.0, 6.0,  6.0,
		    -0.1, -0.1, 0.1, 0.1,
		    -6.0, -6.0, -7.0, -6.0, -6.0,
		     -11.9, -11.9, -12.9,  -11.9, -11.9 }; 

  
 
  
   
  TGeoVolumeAssembly * stlinA = new TGeoVolumeAssembly("0STL");  // A side mother 
  TGeoVolumeAssembly * stlinC = new TGeoVolumeAssembly("0STR");  // C side mother 
 //T0 interior
  TVirtualMC::GetMC()->Gsvolu("0INS","BOX",idtmed[kAir],pinstart,3);
  TGeoVolume *ins = gGeoManager->GetVolume("0INS");
 // 
  TGeoTranslation *tr[52];
  TString nameTr;
  

 //A side Translations
  for (Int_t itr=0; itr<24; itr++) {
    nameTr = Form("0TR%i",itr+1);
    z=-pstartA[2]+pinstart[2];
    tr[itr] = new TGeoTranslation(nameTr.Data(),xa[itr],ya[itr], z );
    printf(" itr %i A %f %f %f \n",itr, xa[itr], ya[itr], z+zdetA);
    tr[itr]->RegisterYourself();
    stlinA->AddNode(ins,itr,tr[itr]);
  }
  
  TGeoRotation *rot[28];
  TString nameRot;
  
  TGeoCombiTrans *com[28];
  TString nameCom;
  
  //C Side Transformations
  for (Int_t itr=24;itr<52; itr++) {
    nameTr = Form("0TR%i",itr+1);
    nameRot = Form("0Rot%i",itr+1);
    //nameCom = Form("0Com%i",itr+1);
    rot[itr-24] = new TGeoRotation(nameRot.Data(),ac[itr-24],bc[itr-24],gc[itr-24]);
    rot[itr-24]->RegisterYourself();
    
    tr[itr] = new TGeoTranslation(nameTr.Data(),xc2[itr-24],yc2[itr-24], (zc2[itr-24]-80.) );
    tr[itr]->Print();
    tr[itr]->RegisterYourself();
      
    //   com[itr-24] = new TGeoCombiTrans(tr[itr],rot[itr-24]);
    com[itr-24] = new TGeoCombiTrans(xc2[itr-24],yc2[itr-24], (zc2[itr-24]-80),rot[itr-24]);
    TGeoHMatrix hm = *com[itr-24];
    TGeoHMatrix *ph = new TGeoHMatrix(hm);
    stlinC->AddNode(ins,itr,ph);
  }
  
  TGeoVolume *alice = gGeoManager->GetVolume("ALIC");
  alice->AddNode(stlinA,1,new TGeoTranslation(0,0, zdetA ) );
  // alice->AddNode(stlinC,1,new TGeoTranslation(0,0, -zdetC ) );
    TGeoRotation * rotC = new TGeoRotation( "rotC",90., 0., 90., 90., 180., 0.);
     alice->AddNode(stlinC,1, new TGeoCombiTrans(0., 0., -zdetC , rotC) );
  
  x=0;
  y=0;
   
  // Entry window (glass)
  TVirtualMC::GetMC()->Gsvolu("0TOP","BOX",idtmed[kOpGlass],ptop,3); //glass
  TGeoVolume *top = gGeoManager->GetVolume("0TOP");
  TVirtualMC::GetMC()->Gsvolu ("0REG", "BOX", idtmed[kOpGlassCathode], preg, 3); 
  TGeoVolume *cat = gGeoManager->GetVolume("0REG");
  TVirtualMC::GetMC()->Gsvolu("0MCP","BOX",idtmed[kGlass],pmcp,3); //glass
  TGeoVolume *mcp = gGeoManager->GetVolume("0MCP");
 
  Int_t ntops=0;
   Float_t xin=0, yin=0;
   for (Int_t ix=0; ix<2; ix++) {
     xin = - pinstart[0] + 0.3 + (ix+0.5)*2*ptop[0] ;
     for (Int_t iy=0; iy<2 ; iy++) {
       z = - pinstart[2]+ptop[2];
       yin = - pinstart[1] + 0.3 + (iy+0.5)*2*ptop[1];
       ntops++;
       ins->AddNode(top, ntops, new TGeoTranslation(xin,yin,z) );
       // printf(" 0TOP  full x %f y %f z %f \n", xin, yin, z);
       z = -pinstart[2] + 2 * ptop[2] + preg[2];
       ins->AddNode(cat, ntops, new TGeoTranslation(xin,yin,z) );

       // printf(" GEOGEO  %i %i %i %f %f %f %f %f %f \n", ntops, ix, iy,
       //  xin,yin,x1[ntops],y1[ntops],x1[ntops]+xin,y1[ntops]+yin);
     }
   }
// MCP
   z=-pinstart[2] + 2*ptop[2] + 2*preg[2] + pmcp[2];
  ins->AddNode(mcp, 1 , new TGeoTranslation(0,0,z) );

  //V0+

  const int kV0PlusColorSci   = 5;
  TGeoMedium *medV0PlusSci = gGeoManager->GetMedium("FIT_V0PlusSci");
  Double_t Pi = TMath::Pi();
  Double_t Sin45    = TMath::Sin(Pi/4.); 
  Double_t Cos45    = TMath::Cos(Pi/4.); 
  Double_t v0PlusPts[16];
  
  //Defining the master volume for V0+
  TGeoVolume *v0Plus = new TGeoVolumeAssembly("V0PLUS");

    /// For boolean sustraction
  for (Int_t i = 0; i < 2; i++) 
    {
      v0PlusPts[0+8*i] = fV0PlusR0-fV0PlusFraWd/2.-fV0PlusFraWd;  v0PlusPts[1+8*i] = -fV0PlusFraWd;
      v0PlusPts[2+8*i] = fV0PlusR0-fV0PlusFraWd/2.-fV0PlusFraWd;  v0PlusPts[3+8*i] = fV0PlusFraWd/2.;
      v0PlusPts[4+8*i] = fV0PlusR5+fV0PlusFraWd/2.+fV0PlusFraWd;  v0PlusPts[5+8*i] = fV0PlusFraWd/2.;
      v0PlusPts[6+8*i] = fV0PlusR5+fV0PlusFraWd/2.+fV0PlusFraWd;  v0PlusPts[7+8*i] = -fV0PlusFraWd;
    }
  new TGeoArb8("sV0PlusCha1",fV0PlusSciWd/1.5,v0PlusPts);
  for (Int_t i = 0; i < 2; i++) 
    {
      v0PlusPts[0+8*i] = fV0PlusR0*Cos45-fV0PlusFraWd;
      v0PlusPts[1+8*i] = (fV0PlusR0-fV0PlusFraWd)*Sin45-fV0PlusFraWd;
      v0PlusPts[2+8*i] = (fV0PlusR0-fV0PlusFraWd/2.)*Cos45-fV0PlusFraWd;
      v0PlusPts[3+8*i] = (fV0PlusR0-fV0PlusFraWd/2.)*Sin45;
      v0PlusPts[4+8*i] = (fV0PlusR5+fV0PlusFraWd/2.)*Cos45+fV0PlusFraWd;
      v0PlusPts[5+8*i] = (fV0PlusR5+fV0PlusFraWd/2.)*Sin45+2.*fV0PlusFraWd;
      v0PlusPts[6+8*i] = (fV0PlusR5+fV0PlusFraWd)*Cos45+fV0PlusFraWd;
      v0PlusPts[7+8*i] = fV0PlusR5*Sin45+fV0PlusFraWd;
    }
  new TGeoArb8("sV0PlusCha2", fV0PlusSciWd/2.+2.*fV0PlusFraWd, v0PlusPts);
  new TGeoCompositeShape("sV0PlusCha","sV0PlusCha1+sV0PlusCha2");

  //Sensitive scintillator
  new TGeoTubeSeg( "sV0PlusR1b", fV0PlusR0+fV0PlusFraWd/2., fV0PlusR1-fV0PlusFraWd/2., fV0PlusSciWd/2., 0, 45);
  new TGeoTubeSeg( "sV0PlusR2b", fV0PlusR1+fV0PlusFraWd/2., fV0PlusR2-fV0PlusFraWd/2., fV0PlusSciWd/2., 0, 45);
  new TGeoTubeSeg( "sV0PlusR3b", fV0PlusR2+fV0PlusFraWd/2., fV0PlusR3-fV0PlusFraWd/2., fV0PlusSciWd/2., 0, 45);
  new TGeoTubeSeg( "sV0PlusR4b", fV0PlusR3+fV0PlusFraWd/2., fV0PlusR4-fV0PlusFraWd/2., fV0PlusSciWd/2., 0, 45);
  new TGeoTubeSeg( "sV0PlusR5b", fV0PlusR4+fV0PlusFraWd/2., fV0PlusR5-fV0PlusFraWd/2., fV0PlusSciWd/2., 0, 45);
  TGeoCompositeShape *sV0PlusR1 = new TGeoCompositeShape("sV0PlusR1","sV0PlusR1b-sV0PlusCha");
  TGeoCompositeShape *sV0PlusR2 = new TGeoCompositeShape("sV0PlusR2","sV0PlusR2b-sV0PlusCha");
  TGeoCompositeShape *sV0PlusR3 = new TGeoCompositeShape("sV0PlusR3","sV0PlusR3b-sV0PlusCha");
  TGeoCompositeShape *sV0PlusR4 = new TGeoCompositeShape("sV0PlusR4","sV0PlusR4b-sV0PlusCha");
  TGeoCompositeShape *sV0PlusR5 = new TGeoCompositeShape("sV0PlusR5","sV0PlusR5b-sV0PlusCha");

  //Definition sector 1
  TGeoVolume *v0Plus1Sec1 = new TGeoVolume("V0Plus1Sec1",sV0PlusR1,medV0PlusSci); 
  TGeoVolume *v0Plus2Sec1 = new TGeoVolume("V0Plus2Sec1",sV0PlusR2,medV0PlusSci);
  TGeoVolume *v0Plus3Sec1 = new TGeoVolume("V0Plus3Sec1",sV0PlusR3,medV0PlusSci);
  TGeoVolume *v0Plus4Sec1 = new TGeoVolume("V0Plus4Sec1",sV0PlusR4,medV0PlusSci);
  TGeoVolume *v0Plus5Sec1 = new TGeoVolume("V0Plus5Sec1",sV0PlusR5,medV0PlusSci); 
  v0Plus1Sec1->SetLineColor(kV0PlusColorSci); 
  v0Plus2Sec1->SetLineColor(kV0PlusColorSci);
  v0Plus3Sec1->SetLineColor(kV0PlusColorSci); 
  v0Plus4Sec1->SetLineColor(kV0PlusColorSci);
  v0Plus5Sec1->SetLineColor(kV0PlusColorSci);
  TGeoVolume *v0PlusSciSec1 = new TGeoVolumeAssembly("V0PlusSciSec1");
  v0PlusSciSec1->AddNode(v0Plus1Sec1,1);
  v0PlusSciSec1->AddNode(v0Plus2Sec1,1);
  v0PlusSciSec1->AddNode(v0Plus3Sec1,1);
  v0PlusSciSec1->AddNode(v0Plus4Sec1,1);
  v0PlusSciSec1->AddNode(v0Plus5Sec1,1);
  TGeoVolume *v0PlusSec1 = new TGeoVolumeAssembly("V0PlusSec1");
  v0PlusSec1->AddNode(v0PlusSciSec1,1);
  
  TGeoRotation *RotSec1 = new TGeoRotation("RotSec1", 90., 0*45., 90., 90.+0*45., 0., 0.);
  v0Plus->AddNode(v0PlusSec1,0+1,RotSec1);  
  
  /// Definition sector 2
  TGeoVolume *v0Plus1Sec2 = new TGeoVolume("V0Plus1Sec2",sV0PlusR1,medV0PlusSci); 
  TGeoVolume *v0Plus2Sec2 = new TGeoVolume("V0Plus2Sec2",sV0PlusR2,medV0PlusSci);
  TGeoVolume *v0Plus3Sec2 = new TGeoVolume("V0Plus3Sec2",sV0PlusR3,medV0PlusSci);
  TGeoVolume *v0Plus4Sec2 = new TGeoVolume("V0Plus4Sec2",sV0PlusR4,medV0PlusSci);
  TGeoVolume *v0Plus5Sec2 = new TGeoVolume("V0Plus5Sec2",sV0PlusR5,medV0PlusSci); 
  v0Plus1Sec2->SetLineColor(kV0PlusColorSci); 
  v0Plus2Sec2->SetLineColor(kV0PlusColorSci);
  v0Plus3Sec2->SetLineColor(kV0PlusColorSci); 
  v0Plus4Sec2->SetLineColor(kV0PlusColorSci);
  v0Plus5Sec2->SetLineColor(kV0PlusColorSci);
  TGeoVolume *v0PlusSciSec2 = new TGeoVolumeAssembly("V0PlusSciSec2");
  v0PlusSciSec2->AddNode(v0Plus1Sec2,1);
  v0PlusSciSec2->AddNode(v0Plus2Sec2,1);
  v0PlusSciSec2->AddNode(v0Plus3Sec2,1);
  v0PlusSciSec2->AddNode(v0Plus4Sec2,1);      
  v0PlusSciSec2->AddNode(v0Plus5Sec2,1);
  TGeoVolume *v0PlusSec2 = new TGeoVolumeAssembly("V0PlusSec2");
  v0PlusSec2->AddNode(v0PlusSciSec2,1);
  
  TGeoRotation *RotSec2 = new TGeoRotation("RotSec2", 90., 1*45., 90., 90.+1*45., 0., 0.);
  v0Plus->AddNode(v0PlusSec2,1+1,RotSec2);  
  
  /// Definition sector 3
  TGeoVolume *v0Plus1Sec3 = new TGeoVolume("V0Plus1Sec3",sV0PlusR1,medV0PlusSci); 
  TGeoVolume *v0Plus2Sec3 = new TGeoVolume("V0Plus2Sec3",sV0PlusR2,medV0PlusSci);
  TGeoVolume *v0Plus3Sec3 = new TGeoVolume("V0Plus3Sec3",sV0PlusR3,medV0PlusSci);
  TGeoVolume *v0Plus4Sec3 = new TGeoVolume("V0Plus4Sec3",sV0PlusR4,medV0PlusSci);
  TGeoVolume *v0Plus5Sec3 = new TGeoVolume("V0Plus5Sec3",sV0PlusR5,medV0PlusSci); 
  v0Plus1Sec3->SetLineColor(kV0PlusColorSci); 
  v0Plus2Sec3->SetLineColor(kV0PlusColorSci);
  v0Plus3Sec3->SetLineColor(kV0PlusColorSci); 
  v0Plus4Sec3->SetLineColor(kV0PlusColorSci);
  v0Plus5Sec3->SetLineColor(kV0PlusColorSci);
  TGeoVolume *v0PlusSciSec3 = new TGeoVolumeAssembly("V0PlusSciSec3");
  v0PlusSciSec3->AddNode(v0Plus1Sec3,1);
  v0PlusSciSec3->AddNode(v0Plus2Sec3,1);
  v0PlusSciSec3->AddNode(v0Plus3Sec3,1);
  v0PlusSciSec3->AddNode(v0Plus4Sec3,1);      
  v0PlusSciSec3->AddNode(v0Plus5Sec3,1);
  TGeoVolume *v0PlusSec3 = new TGeoVolumeAssembly("V0PlusSec3");
  v0PlusSec3->AddNode(v0PlusSciSec3,1);
  
  TGeoRotation *RotSec3 = new TGeoRotation("RotSec3", 90., 2*45., 90., 90.+2*45., 0., 0.);
  v0Plus->AddNode(v0PlusSec3,2+1,RotSec3);
  
  /// Definition sector 4
  TGeoVolume *v0Plus1Sec4 = new TGeoVolume("V0Plus1Sec4",sV0PlusR1,medV0PlusSci); 
  TGeoVolume *v0Plus2Sec4 = new TGeoVolume("V0Plus2Sec4",sV0PlusR2,medV0PlusSci);
  TGeoVolume *v0Plus3Sec4 = new TGeoVolume("V0Plus3Sec4",sV0PlusR3,medV0PlusSci);
  TGeoVolume *v0Plus4Sec4 = new TGeoVolume("V0Plus4Sec4",sV0PlusR4,medV0PlusSci);
  TGeoVolume *v0Plus5Sec4 = new TGeoVolume("V0Plus5Sec4",sV0PlusR5,medV0PlusSci); 
  v0Plus1Sec4->SetLineColor(kV0PlusColorSci); 
  v0Plus2Sec4->SetLineColor(kV0PlusColorSci);
  v0Plus3Sec4->SetLineColor(kV0PlusColorSci); 
  v0Plus4Sec4->SetLineColor(kV0PlusColorSci);
  v0Plus5Sec4->SetLineColor(kV0PlusColorSci);
  TGeoVolume *v0PlusSciSec4 = new TGeoVolumeAssembly("V0PlusSciSec4");
  v0PlusSciSec4->AddNode(v0Plus1Sec4,1);
  v0PlusSciSec4->AddNode(v0Plus2Sec4,1);
  v0PlusSciSec4->AddNode(v0Plus3Sec4,1);
  v0PlusSciSec4->AddNode(v0Plus4Sec4,1);      
  v0PlusSciSec4->AddNode(v0Plus5Sec4,1);
  TGeoVolume *v0PlusSec4 = new TGeoVolumeAssembly("V0PlusSec4");
  v0PlusSec4->AddNode(v0PlusSciSec4,1);
  
  TGeoRotation *RotSec4 = new TGeoRotation("RotSec4", 90., 3*45., 90., 90.+3*45., 0., 0.);
  v0Plus->AddNode(v0PlusSec4,3+1,RotSec4);
  
  /// Definition sector 5
  TGeoVolume *v0Plus1Sec5 = new TGeoVolume("V0Plus1Sec5",sV0PlusR1,medV0PlusSci); 
  TGeoVolume *v0Plus2Sec5 = new TGeoVolume("V0Plus2Sec5",sV0PlusR2,medV0PlusSci);
  TGeoVolume *v0Plus3Sec5 = new TGeoVolume("V0Plus3Sec5",sV0PlusR3,medV0PlusSci);
  TGeoVolume *v0Plus4Sec5 = new TGeoVolume("V0Plus4Sec5",sV0PlusR4,medV0PlusSci);
  TGeoVolume *v0Plus5Sec5 = new TGeoVolume("V0Plus5Sec5",sV0PlusR5,medV0PlusSci); 
  v0Plus1Sec5->SetLineColor(kV0PlusColorSci); 
  v0Plus2Sec5->SetLineColor(kV0PlusColorSci);
  v0Plus3Sec5->SetLineColor(kV0PlusColorSci); 
  v0Plus4Sec5->SetLineColor(kV0PlusColorSci);
  v0Plus5Sec5->SetLineColor(kV0PlusColorSci);
  TGeoVolume *v0PlusSciSec5 = new TGeoVolumeAssembly("V0PlusSciSec5");
  v0PlusSciSec5->AddNode(v0Plus1Sec5,1);
  v0PlusSciSec5->AddNode(v0Plus2Sec5,1);
  v0PlusSciSec5->AddNode(v0Plus3Sec5,1);
  v0PlusSciSec5->AddNode(v0Plus4Sec5,1);      
  v0PlusSciSec5->AddNode(v0Plus5Sec5,1);
  TGeoVolume *v0PlusSec5 = new TGeoVolumeAssembly("V0PlusSec5");
  v0PlusSec5->AddNode(v0PlusSciSec5,1);
 
  TGeoRotation *RotSec5 = new TGeoRotation("RotSec5", 90., 4*45., 90., 90.+4*45., 0., 0.);
  v0Plus->AddNode(v0PlusSec5,4+1,RotSec5);
  
  /// Definition sector 6
  TGeoVolume *v0Plus1Sec6 = new TGeoVolume("V0Plus1Sec6",sV0PlusR1,medV0PlusSci); 
  TGeoVolume *v0Plus2Sec6 = new TGeoVolume("V0Plus2Sec6",sV0PlusR2,medV0PlusSci);
  TGeoVolume *v0Plus3Sec6 = new TGeoVolume("V0Plus3Sec6",sV0PlusR3,medV0PlusSci);
  TGeoVolume *v0Plus4Sec6 = new TGeoVolume("V0Plus4Sec6",sV0PlusR4,medV0PlusSci);
  TGeoVolume *v0Plus5Sec6 = new TGeoVolume("V0Plus5Sec6",sV0PlusR5,medV0PlusSci); 
  v0Plus1Sec6->SetLineColor(kV0PlusColorSci); 
  v0Plus2Sec6->SetLineColor(kV0PlusColorSci);
  v0Plus3Sec6->SetLineColor(kV0PlusColorSci); 
  v0Plus4Sec6->SetLineColor(kV0PlusColorSci);
  v0Plus5Sec6->SetLineColor(kV0PlusColorSci);
  TGeoVolume *v0PlusSciSec6 = new TGeoVolumeAssembly("V0PlusSciSec6");
  v0PlusSciSec6->AddNode(v0Plus1Sec6,1);
  v0PlusSciSec6->AddNode(v0Plus2Sec6,1);
  v0PlusSciSec6->AddNode(v0Plus3Sec6,1);
  v0PlusSciSec6->AddNode(v0Plus4Sec6,1);      
  v0PlusSciSec6->AddNode(v0Plus5Sec6,1);
  TGeoVolume *v0PlusSec6 = new TGeoVolumeAssembly("V0PlusSec6");
  v0PlusSec6->AddNode(v0PlusSciSec6,1);
  
  TGeoRotation *RotSec6 = new TGeoRotation("RotSec6", 90., 5*45., 90., 90.+5*45., 0., 0.);
  v0Plus->AddNode(v0PlusSec6,5+1,RotSec6);
  
  /// Definition sector 7
  TGeoVolume *v0Plus1Sec7 = new TGeoVolume("V0Plus1Sec7",sV0PlusR1,medV0PlusSci); 
  TGeoVolume *v0Plus2Sec7 = new TGeoVolume("V0Plus2Sec7",sV0PlusR2,medV0PlusSci);
  TGeoVolume *v0Plus3Sec7 = new TGeoVolume("V0Plus3Sec7",sV0PlusR3,medV0PlusSci);
  TGeoVolume *v0Plus4Sec7 = new TGeoVolume("V0Plus4Sec7",sV0PlusR4,medV0PlusSci);
  TGeoVolume *v0Plus5Sec7 = new TGeoVolume("V0Plus5Sec7",sV0PlusR5,medV0PlusSci); 
  v0Plus1Sec7->SetLineColor(kV0PlusColorSci); 
  v0Plus2Sec7->SetLineColor(kV0PlusColorSci);
  v0Plus3Sec7->SetLineColor(kV0PlusColorSci); 
  v0Plus4Sec7->SetLineColor(kV0PlusColorSci);
  v0Plus5Sec7->SetLineColor(kV0PlusColorSci);
  TGeoVolume *v0PlusSciSec7 = new TGeoVolumeAssembly("V0PlusSciSec7");
  v0PlusSciSec7->AddNode(v0Plus1Sec7,1);
  v0PlusSciSec7->AddNode(v0Plus2Sec7,1);
  v0PlusSciSec7->AddNode(v0Plus3Sec7,1);
  v0PlusSciSec7->AddNode(v0Plus4Sec7,1);      
  v0PlusSciSec7->AddNode(v0Plus5Sec7,1);
  TGeoVolume *v0PlusSec7 = new TGeoVolumeAssembly("V0PlusSec7");
  v0PlusSec7->AddNode(v0PlusSciSec7,1);
 
  TGeoRotation *RotSec7 = new TGeoRotation("RotSec7", 90., 6*45., 90., 90.+6*45., 0., 0.);
  v0Plus->AddNode(v0PlusSec7,6+1,RotSec7);
  
  /// Definition sector 8
  TGeoVolume *v0Plus1Sec8 = new TGeoVolume("V0Plus1Sec8",sV0PlusR1,medV0PlusSci); 
  TGeoVolume *v0Plus2Sec8 = new TGeoVolume("V0Plus2Sec8",sV0PlusR2,medV0PlusSci);
  TGeoVolume *v0Plus3Sec8 = new TGeoVolume("V0Plus3Sec8",sV0PlusR3,medV0PlusSci);
  TGeoVolume *v0Plus4Sec8 = new TGeoVolume("V0Plus4Sec8",sV0PlusR4,medV0PlusSci);
  TGeoVolume *v0Plus5Sec8 = new TGeoVolume("V0Plus5Sec8",sV0PlusR5,medV0PlusSci); 
  v0Plus1Sec8->SetLineColor(kV0PlusColorSci); 
  v0Plus2Sec8->SetLineColor(kV0PlusColorSci);
  v0Plus3Sec8->SetLineColor(kV0PlusColorSci); 
  v0Plus4Sec8->SetLineColor(kV0PlusColorSci);
  v0Plus5Sec8->SetLineColor(kV0PlusColorSci);
  TGeoVolume *v0PlusSciSec8 = new TGeoVolumeAssembly("V0PlusSciSec8");
  v0PlusSciSec8->AddNode(v0Plus1Sec8,1);
  v0PlusSciSec8->AddNode(v0Plus2Sec8,1);
  v0PlusSciSec8->AddNode(v0Plus3Sec8,1);
  v0PlusSciSec8->AddNode(v0Plus4Sec8,1);      
  v0PlusSciSec8->AddNode(v0Plus5Sec8,1);
  TGeoVolume *v0PlusSec8 = new TGeoVolumeAssembly("V0PlusSec8");
  v0PlusSec8->AddNode(v0PlusSciSec8,1);
  
  TGeoRotation *RotSec8 = new TGeoRotation("RotSec8", 90., 7*45., 90., 90.+7*45., 0., 0.);
  v0Plus->AddNode(v0PlusSec8,7+1,RotSec8);

  alice->AddNode(v0Plus,1,new TGeoTranslation(0, 0, fV0PlusZposition));
 
}    
//------------------------------------------------------------------------
void AliFITv5::AddAlignableVolumes() const
{
  //
  // Create entries for alignable volumes associating the symbolic volume
  // name with the corresponding volume path. Needs to be syncronized with
  // eventual changes in the geometry.
  //
  TString volPath;
  TString symName, sn;
  TString vpAalign = "/ALIC_1/0STL_1";
  TString vpCalign = "/ALIC_1/0STR_1";
  for (Int_t imod=0; imod<2; imod++)  {
    if (imod==0) {volPath  = vpCalign; symName="/ALIC_1/0STL"; }
    if (imod==1) {volPath  = vpAalign; symName="/ALIC_1/0STR"; }
    
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
void AliFITv5::CreateMaterials()
{

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
 // Radiator  glass SiO2
   Float_t aglass[2]={28.0855,15.9994};
   Float_t zglass[2]={14.,8.};
   Float_t wglass[2]={1.,2.};
   Float_t dglass=2.65;
 // MCP glass SiO2
   Float_t dglass_mcp=1.3;
//*** Definition Of avaible T0 materials ***
   AliMixture(1, "Vacuum$", aAir, zAir, dAir1,4,wAir);
   AliMixture(2, "Air$", aAir, zAir, dAir,4,wAir);
   AliMixture( 4, "MCP glass   $",aglass,zglass,dglass_mcp,-2,wglass);
   AliMixture( 24, "Radiator Optical glass$",aglass,zglass,dglass,-2,wglass);
   
   AliMedium(1, "Air$", 2, 0, isxfld, sxmgmx, 10., .1, 1., .003, .003);
   AliMedium(3, "Vacuum$", 1, 0, isxfld, sxmgmx, 10., .01, .1, .003, .003);
   AliMedium(6, "Glass$", 4, 0, isxfld, sxmgmx, 10., .01, .1, .003, .003);
  
   AliMedium(16, "OpticalGlass$", 24, 1, isxfld, sxmgmx, 10., .01, .1, .003, .003);
    AliMedium(19, "OpticalGlassCathode$", 24, 1, isxfld, sxmgmx, 10., .01, .1, .003, .003);
   AliMedium(22, "SensAir$", 2, 1, isxfld, sxmgmx, 10., .1, 1., .003, .003);

   //V0+

  // Parameters for simulation scope
  Int_t     FieldType       = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ();     // Field type 
  Double_t  MaxField        = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max();       // Field max.
  Double_t  MaxBending      = 10;    // Max Angle
  Double_t  MaxStepSize     = 0.01;  // Max step size 
  Double_t  MaxEnergyLoss   = 1;     // Max Delta E
  Double_t  Precision       = 0.003; // Precision
  Double_t  MinStepSize     = 0.003; // Minimum step size 

  Int_t    Id;
  Double_t A, Z, RadLength, AbsLength;
  Float_t Density, as[4], zs[4], ws[4];

// Parameters for V0Plusscintilator: BC404
   as[0] = 1.00794;	as[1] = 12.011;
   zs[0] = 1.;		zs[1] = 6.;
   ws[0] = 5.21;	ws[1] = 4.74;
   Density      = 1.032;
   Id           = 5;
   AliMixture(Id, "V0PlusSci", as, zs, Density, -2, ws);
   AliMedium(Id,  "V0PlusSci", Id, 1, FieldType, MaxField, MaxBending, MaxStepSize, MaxEnergyLoss, Precision, MinStepSize);
   
   AliDebugClass(1,": ++++++++++++++Medium set++++++++++");
   
   
}

//-------------------------------------------------------------------
void AliFITv5::DefineOpticalProperties()
{


// Optical properties definition.
   Int_t *idtmed = fIdtmed->GetArray();
// Definition Cherenkov parameters
   int i;
   const Int_t kNbins=31;
   
   Float_t rindexSiO2[kNbins],  efficAll[kNbins], rindexAir[kNbins], absorAir[kNbins],rindexCathodeNext[kNbins], absorbCathodeNext[kNbins];
   Double_t efficMet[kNbins], aReflMet[kNbins];
           
   // quartz 20mm
   Float_t aAbsSiO2[kNbins]={29.0, 28.6, 28.3, 27.7, 27.3, 26.7, 26.4, 
			     25.9, 25.3, 24.9, 24.5, 23.7, 
			     23.2, 22.8, 22.4, 21.8, 21.3,
			     22.8, 22.1, 21.7, 21.2, 20.5, 
			     19.9, 19.3, 18.7, 18.0, 17.1, 
			     16.3, 15.3, 14.3, 14.3   };
          
   Float_t aPckov[kNbins]  ={3.87, 3.94, 4.02, 4.11, 4.19, 4.29, 4.38,
			     4.48, 4.58, 4.69, 4.81, 4.93, 
			     5.05, 5.19, 5.33, 5.48, 5.63,
			     5.8,  5.97, 6.16, 6.36, 6.57, 
			     6.8,  7.04, 7.3,  7.58, 7.89, 
			     8.22, 8.57, 8.97, 9.39 };  
  Double_t dPckov[kNbins]  ={3.87, 3.94, 4.02, 4.11, 4.19, 4.29, 4.38,
                             4.48, 4.58, 4.69, 4.81, 4.93,
                             5.05, 5.19, 5.33, 5.48, 5.63,
                             5.8,  5.97, 6.16, 6.36, 6.57,
                             6.8,  7.04, 7.3,  7.58, 7.89,
                             8.22, 8.57, 8.97, 9.39 };

   /*     
   Float_t effCathode[kNbins]={0.11, 0.13, 0.15, 0.16, 0.18, 0.19, 0.20,
			      0.21, 0.22, 0.23, 0.24, 0.26, 
			      0.27, 0.29, 0.30, 0.29, 0.29, 
			      0.28, 0.28, 0.27, 0.26, 0.25, 
			      0.25, 0.23, 0.20, 0.19, 0.17,
			      0.17, 0.17, 0.2, 0.23};
   */     
   //  Float_t aAbsSiO2[kNbins]; //quartz 30mm
 for(i=0;i<kNbins;i++)

    {
      aPckov[i]=aPckov[i]*1e-9;//Photons energy bins 4 eV - 8.5 eV step 0.1 eV
      dPckov[i]=dPckov[i]*1e-9;//Photons energy bins 4 eV - 8.5 eV step 0.1 eV 
      //      rindexAir[i]=0.0001;
      rindexAir[i] = 1.;
      rindexSiO2[i]=1.458; //refractive index for qwarts
      rindexCathodeNext[i]=0;
      efficAll[i]=1.;
      efficMet[i]=0.;
      aReflMet[i]=1.;
      //      aAbsSiO2[i]=28.5; //quartz 30mm
      absorAir[i]=0.3;      
      absorbCathodeNext[i]=0;

    }
  
  TVirtualMC::GetMC()->SetCerenkov (idtmed[kOpGlass], kNbins, aPckov, aAbsSiO2, efficAll, rindexSiO2 );
   // TVirtualMC::GetMC()->SetCerenkov (idtmed[kOpGlassCathode], kNbins, aPckov, aAbsSiO2, effCathode, rindexSiO2 );
   TVirtualMC::GetMC()->SetCerenkov (idtmed[kOpGlassCathode], kNbins, aPckov, aAbsSiO2,efficAll , rindexSiO2 );
   //  TVirtualMC::GetMC()->SetCerenkov (idtmed[kOpAir], kNbins, aPckov,absorAir , efficAll,rindexAir );
   //   TVirtualMC::GetMC()->SetCerenkov (idtmed[kOpAirNext], kNbins, aPckov,absorbCathodeNext , efficAll, rindexCathodeNext);

   //Define a boarder for radiator optical properties
   TVirtualMC::GetMC()->DefineOpSurface("surfRd", kUnified /*kGlisur*/,kDielectric_metal,kPolished, 0.);
   TVirtualMC::GetMC()->SetMaterialProperty("surfRd", "EFFICIENCY", kNbins, dPckov, efficMet);
   TVirtualMC::GetMC()->SetMaterialProperty("surfRd", "REFLECTIVITY", kNbins, dPckov, aReflMet);


}

//-------------------------------------------------------------------
void AliFITv5::Init()
{
// Initialises version 0 of the Forward Multiplicity Detector
//
  AliFIT::Init();
  fIdSens1=TVirtualMC::GetMC()->VolId("0REG");

  //Defining the sensitive volumes
  for (Int_t Sec = 0; Sec < nSectors; Sec++)
    {
      for (Int_t Ring = 0; Ring < nRings; Ring++)
	{
	  fIdV0Plus[Sec][Ring] = TVirtualMC::GetMC()->VolId(Form("V0Plus%dSec%d",Ring+1,Sec+1));
	}
    }

   AliDebug(1,Form("%s: *** FIT version 1 initialized ***\n",ClassName()));
}

//-------------------------------------------------------------------

void AliFITv5::StepManager()
{
  //
  // Called for every step in the T0 Detector
  //
  Int_t id,copy,copy1;
  static Float_t hits[13];
  static Int_t vol[3];
  TLorentzVector pos;
  TLorentzVector mom;
  //   TClonesArray &lhits = *fHits;
  
  if(!TVirtualMC::GetMC()->IsTrackAlive()) return; // particle has disappeared
  
  id=TVirtualMC::GetMC()->CurrentVolID(copy);  
  // Check the sensetive volume
  if(id==fIdSens1 ) { 
    if(TVirtualMC::GetMC()->IsTrackEntering()) {
      TVirtualMC::GetMC()->CurrentVolOffID(1,copy1);
      vol[1] = copy1;
      vol[0]=copy;
      TVirtualMC::GetMC()->TrackPosition(pos);
      TVirtualMC::GetMC()->TrackMomentum(mom);
      Float_t Pt  = TMath::Sqrt( mom.Px() * mom.Px() + mom.Py() * mom.Py() );
      TParticle *Particle = gAlice->GetMCApp()->Particle(gAlice->GetMCApp()->GetCurrentTrackNumber());
      hits[0] = pos[0];
      hits[1] = pos[1];
      hits[2] = pos[2];
      if(pos[2]<0) vol[2] = 0;
      else vol[2] = 1;
      Float_t etot=TVirtualMC::GetMC()->Etot();
      hits[3]=etot;
      Int_t iPart= TVirtualMC::GetMC()->TrackPid();
      Int_t partID=TVirtualMC::GetMC()->IdFromPDG(iPart);
      hits[4]=partID;
      Float_t ttime=TVirtualMC::GetMC()->TrackTime();
      hits[5]=ttime*1e12;
      hits[6] = TVirtualMC::GetMC()->TrackCharge();
      hits[7] = mom.Px();
      hits[8] = mom.Py();
      hits[9] = mom.Pz();
      hits[10] = fSenseless;//Energy loss is sensless for T0+
      hits[11] = fSenseless;//Track length is sensless for T0+
      hits[12] = fSenseless;//Photon production for V0+
      if (TVirtualMC::GetMC()->TrackPid() == 50000050)   // If particles is photon then ...
	{
	  if(RegisterPhotoE(hits[3])) {
	    fIshunt = 2;
	    AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(),vol,hits);
	    // Create a track reference at the exit of photocatode
	  }	    
	}
      //charge particle HITS
      if ( TVirtualMC::GetMC()->TrackCharge() ) {
	fIshunt = 0;
	AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(),vol,hits);
	//	printf("T0 :::volumes pmt %i mcp %i vol %i x %f y %f z %f particle %f all \n",  vol[0], vol[1],  vol[2], hits[0], hits[1], hits[2], hits[4]);
      }
            
      //charge particle TrackReference
      if ( TVirtualMC::GetMC()->TrackCharge() )
	AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber(), AliTrackReference::kFIT);
      
    }// trck entering		
  } //sensitive




  //V0+

  if ( !gMC->TrackCharge() || !gMC->IsTrackAlive() ) return;//Only interested in charged and alive tracks

  //Check if there is a hit in any of the V0+ sensitive volumes defined in Init
  Bool_t IsId = kFALSE;
  for (Int_t i = 0; i < nSectors; i++)
    {
      for (Int_t j = 0; j < nRings; j++)
	{
	  if ( id == fIdV0Plus[i][j] ) 
	    {
	      IsId = kTRUE;
	      break;
	    }
	}
    }

  if ( IsId == kTRUE ) //If yes, then perform the following
    {

      //Defining V0+ ring numbers using the sensitive volumes
      Int_t RingNumber;      

      if ( (id == fIdV0Plus[0][0]) || (id == fIdV0Plus[1][0]) || (id == fIdV0Plus[2][0]) || (id == fIdV0Plus[3][0]) || (id == fIdV0Plus[4][0]) || (id == fIdV0Plus[5][0]) || (id == fIdV0Plus[6][0]) || (id == fIdV0Plus[7][0]) ) RingNumber = 1;

      else if ( (id == fIdV0Plus[0][1]) || (id == fIdV0Plus[1][1]) || (id == fIdV0Plus[2][1]) || (id == fIdV0Plus[3][1]) || (id == fIdV0Plus[4][1]) || (id == fIdV0Plus[5][1]) || (id == fIdV0Plus[6][1]) || (id == fIdV0Plus[7][1]) ) RingNumber = 2;

      else if ( (id == fIdV0Plus[0][2]) || (id == fIdV0Plus[1][2]) || (id == fIdV0Plus[2][2]) || (id == fIdV0Plus[3][2]) || (id == fIdV0Plus[4][2]) || (id == fIdV0Plus[5][2]) || (id == fIdV0Plus[6][2]) || (id == fIdV0Plus[7][2]) ) RingNumber = 3;

      else if ( (id == fIdV0Plus[0][3]) || (id == fIdV0Plus[1][3]) || (id == fIdV0Plus[2][3]) || (id == fIdV0Plus[3][3]) || (id == fIdV0Plus[4][3]) || (id == fIdV0Plus[5][3]) || (id == fIdV0Plus[6][3]) || (id == fIdV0Plus[7][3]) ) RingNumber = 4;

      else if ( (id == fIdV0Plus[0][4]) || (id == fIdV0Plus[1][4]) || (id == fIdV0Plus[2][4]) || (id == fIdV0Plus[3][4]) || (id == fIdV0Plus[4][4]) || (id == fIdV0Plus[5][4]) || (id == fIdV0Plus[6][4]) || (id == fIdV0Plus[7][4]) ) RingNumber = 5;

      else RingNumber = 0;

      if ( RingNumber )
	{

	  if(TVirtualMC::GetMC()->IsTrackEntering()) 
	    {
	      TVirtualMC::GetMC()->TrackPosition(pos);
	      TVirtualMC::GetMC()->TrackMomentum(mom);
	      Float_t Pt  = TMath::Sqrt( mom.Px() * mom.Px() + mom.Py() * mom.Py() );
	      TParticle *Particle = gAlice->GetMCApp()->Particle(gAlice->GetMCApp()->GetCurrentTrackNumber());
	      hits[0] = pos[0];
	      hits[1] = pos[1];
	      hits[2] = pos[2];
	      Float_t etot=TVirtualMC::GetMC()->Etot();
	      hits[3]=etot;
	      Int_t iPart= TVirtualMC::GetMC()->TrackPid();
	      Int_t partID=TVirtualMC::GetMC()->IdFromPDG(iPart);
	      hits[4]=partID;
	      Float_t ttime=TVirtualMC::GetMC()->TrackTime();
	      hits[5]=ttime*1e12;
	      hits[6] = TVirtualMC::GetMC()->TrackCharge();
	      hits[7] = mom.Px();
	      hits[8] = mom.Py();
	      hits[9] = mom.Pz();
	    }//Track entering
	  if( gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared())
	    {
	      Float_t Eloss, Tlength;    
	      Float_t EnergyDep = TVirtualMC::GetMC()->Edep();
	      Float_t Step   = TVirtualMC::GetMC()->TrackStep();//Energy loss in the current step
	      Int_t nPhotonsInStep = 0;
	      Int_t nPhotons = 0;
	      nPhotonsInStep = Int_t(EnergyDep / (fV0PlusLightYield*1e-9) );
	      nPhotons = nPhotonsInStep - Int_t((Float_t(nPhotonsInStep) * fV0PlusLightAttenuation * fV0PlusnMeters));
	      nPhotons = nPhotons - Int_t( Float_t(nPhotons) * fV0PlusFibToPhot);
	      Eloss += EnergyDep;
	      Tlength += Step;
	      hits[10] = Eloss;
	      hits[11] = Tlength;
	      hits[12] = nPhotons;
	      vol[0] = GetCellId(vol);
	      vol[1] = RingNumber;
	      vol[2] = 2;
	      fIshunt = 0;
	      AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(), vol, hits);
	      AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber(), AliTrackReference::kFIT);
	      Tlength = 0.0;
	      Eloss = 0.0;
	    }//Exiting, stopped or disappeared track
	}//Ring number
    }
  
}


//------------------------------------------------------------------------
Bool_t AliFITv5::RegisterPhotoE(Double_t energy)
{
  
  
  //  Float_t hc=197.326960*1.e6; //mev*nm
  Double_t hc=1.973*1.e-6; //gev*nm
  Float_t lambda=hc/energy;
  Float_t eff = fPMTeff->Eval(lambda);
  Double_t  p = gRandom->Rndm();
  
  if (p > eff)
    return kFALSE;
  
  return kTRUE;
}

//----------------------------------------------------------------------------

void AliFITv5::SetPMTeff()
{
  Float_t lambda[50];
  Float_t eff[50 ] = {0,        0,       0.23619,  0.202909, 0.177913, 
		    0.175667, 0.17856, 0.190769, 0.206667, 0.230286,
		    0.252276, 0.256267,0.26,     0.27125,  0.281818,
		    0.288118, 0.294057,0.296222, 0.301622, 0.290421, 
		    0.276615, 0.2666,  0.248,    0.23619,  0.227814, 
		    0.219818, 0.206667,0.194087, 0.184681, 0.167917, 
		    0.154367, 0.1364,  0.109412, 0.0834615,0.0725283, 
		    0.0642963,0.05861, 0.0465,   0.0413333,0.032069, 
		    0.0252203,0.02066, 0.016262, 0.012,    0.00590476,
		    0.003875, 0.00190, 0,        0,        0          } ;
  for (Int_t i=0; i<50; i++) lambda[i]=200+10*i; 

  fPMTeff = new TGraph(50,lambda,eff);
 
}

//-------------------------------------------------------------------------
Int_t AliFITv5::GetCellId(Int_t *vol)
{

  fCellId = 0;
  
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus1Sec1")) fCellId = 1;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus1Sec2")) fCellId = 2;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus1Sec3")) fCellId = 3;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus1Sec4")) fCellId = 4;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus1Sec5")) fCellId = 5;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus1Sec6")) fCellId = 6;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus1Sec7")) fCellId = 7;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus1Sec8")) fCellId = 8;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus2Sec1")) fCellId = 9;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus2Sec2")) fCellId = 10;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus2Sec3")) fCellId = 11;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus2Sec4")) fCellId = 12;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus2Sec5")) fCellId = 13;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus2Sec6")) fCellId = 14;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus2Sec7")) fCellId = 15;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus2Sec8")) fCellId = 16;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus3Sec1")) fCellId = 17;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus3Sec2")) fCellId = 18;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus3Sec3")) fCellId = 19;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus3Sec4")) fCellId = 20;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus3Sec5")) fCellId = 21;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus3Sec6")) fCellId = 22;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus3Sec7")) fCellId = 23;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus3Sec8")) fCellId = 24;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus4Sec1")) fCellId = 25;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus4Sec2")) fCellId = 26;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus4Sec3")) fCellId = 27;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus4Sec4")) fCellId = 28;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus4Sec5")) fCellId = 29;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus4Sec6")) fCellId = 30;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus4Sec7")) fCellId = 31;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus4Sec8")) fCellId = 32;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus5Sec1")) fCellId = 33;	
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus5Sec2")) fCellId = 34;	
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus5Sec3")) fCellId = 35;	
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus5Sec4")) fCellId = 36;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus5Sec5")) fCellId = 37;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus5Sec6")) fCellId = 38;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus5Sec7")) fCellId = 39;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus5Sec8")) fCellId = 40;
 
  return fCellId;

}
