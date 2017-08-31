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
<img src="gif/AliFITv6Class.gif">
*/
//End Html
// Alla.Maevskaya@cern.ch
////T0+ optical propreties from Maciej and Noa
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
#include "AliFITv7.h"

#include "AliMC.h"
#include "AliCDBLocal.h"
#include "AliCDBStorage.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliTrackReference.h"

ClassImp(AliFITv7)

using std::cout;
using std::endl;

//--------------------------------------------------------------------
AliFITv7::AliFITv7():  AliFIT(),
  fIdSens1(0),
  fIdSens2(0),
  fPMTeff(0x0),
//V0+
  nSectors(16),
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
AliFITv7::AliFITv7(const char *name, const char *title):
  AliFIT(name,title),
  fIdSens1(0),
  fIdSens2(0),
  fPMTeff(0x0),
//V0+
  nSectors(16),
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
  //  SetPMTeff();
}
//_____________________________________________________________________________

AliFITv7::~AliFITv7() 
{
  // desctructor  
}

//-------------------------------------------------------------------------
void AliFITv7::CreateGeometry()
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
  //  Float_t pstartC[3] = {6., 20 ,5};
  //  Float_t pstartA[3] = {2.55, 20 ,5};
  Float_t pstartC[3] = {20., 20 ,5};
  Float_t pstartA[3] = {20, 20 ,5};
  Float_t pinstart[3] = {2.95,2.95,4.34};
  Float_t pmcp[3] = {2.949, 2.949, 2.8}; //MCP

  AliMatrix(idrotm[901], 90., 0., 90., 90., 180., 0.);
  
  //-------------------------------------------------------------------
  //-------------------------------------------------------------------
  //C side Concave Geometry

  Double_t crad = 82.;		//define concave c-side radius here
    
  Double_t dP = 3.31735114408; 	// Work in Progress side length

  //uniform angle between detector faces==
  Double_t btta = 2*TMath::ATan(dP/crad);
 
  //get noncompensated translation data
  Double_t grdin[6] = {-3, -2, -1, 1, 2, 3};  
  Double_t gridpoints[6];
  for(Int_t i = 0; i < 6; i++){
    gridpoints[i] = crad*TMath::Sin((1 - 1/(2*TMath::Abs(grdin[i])))*grdin[i]*btta);
  } 

  std::vector<Double_t> xi,yi;
    
  for(Int_t j = 5; j >= 0; j--){
      for(Int_t i = 0; i < 6; i++){
          if(!(((j == 5 || j == 0) && (i == 5 || i == 0)) || 
               ((j == 3 || j == 2) && (i == 3 || i == 2))))
          {
              xi.push_back(gridpoints[i]);
              yi.push_back(gridpoints[j]);
          }
      }
      
  }
  
  Double_t zi[28];
  for(Int_t i = 0; i < 28; i++) {
    zi[i] = TMath::Sqrt(TMath::Power(crad, 2) - TMath::Power(xi[i], 2) - TMath::Power(yi[i], 2));
  }

  //get rotation data
  Double_t ac[28], bc[28], gc[28];
  for(Int_t i = 0; i < 28; i++) {
    ac[i] = TMath::ATan(yi[i]/xi[i]) - TMath::Pi()/2 + 2*TMath::Pi();
    if(xi[i] < 0){
      bc[i] = TMath::ACos(zi[i]/crad);
    }
    else {
      bc[i] = -1 * TMath::ACos(zi[i]/crad);
    }
  }
  Double_t xc2[28], yc2[28], zc2[28];
    
    
  //compensation based on node position within individual detector geometries
  //determine compensated radius
  Double_t rcomp = crad + pstartC[2] / 2.0; //
  for(Int_t i = 0; i < 28; i++) {
    //Get compensated translation data
    xc2[i] = rcomp*TMath::Cos(ac[i] + TMath::Pi()/2)*TMath::Sin(-1*bc[i]);
    yc2[i] = rcomp*TMath::Sin(ac[i] + TMath::Pi()/2)*TMath::Sin(-1*bc[i]);
    zc2[i] = rcomp*TMath::Cos(bc[i]);

    //Convert angles to degrees
    ac[i]*=180/TMath::Pi();
    bc[i]*=180/TMath::Pi();
    gc[i] = -1 * ac[i];
  }
  
  //A Side
 
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
 //FIT interior
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
  
  SetOneMCP(ins);
  
  SetVZEROGeo(alice);
  
}
//_________________________________________
void AliFITv7::SetOneMCP(TGeoVolume *ins)
{
  cout<<"AliFITv7::SetOneMCP "<<ins<<endl;
  Double_t x,y,z;
  Double_t crad = 82.;		//define concave c-side radius here
  Double_t dP = 3.31735114408; 	// Work in Progress side length
  
  Int_t *idtmed = fIdtmed->GetArray();
  Float_t pinstart[3] = {2.95,2.95,4.34};
  Float_t pmcp[3] = {2.949, 2.949, 2.8}; //MCP
  Float_t ptop[3] = {1.324, 1.324, 1.};//cherenkov radiator
  Float_t ptopref[3] = {1.3241, 1.3241, 1.};//cherenkov radiator wrapped with reflection 
  Float_t preg[3] = {1.324, 1.324, 0.005};//photcathode 
  Double_t prfv[3]= {0.0002,1.323, 1.};//vertical refracting layer bettwen radiators and bettwen radiator and not optical Air
  Double_t prfh[3]= {1.323,0.0002, 1.};//horizontal refracting layer bettwen radiators a 
  Double_t pal[3]= {2.648,2.648, 0.25};  // 5mm Al top on th eeach radiator
  // Entry window (glass)
  TVirtualMC::GetMC()->Gsvolu("0TOP","BOX",idtmed[kOpGlass],ptop,3); //glass radiator
  TGeoVolume *top = gGeoManager->GetVolume("0TOP");
  TVirtualMC::GetMC()->Gsvolu("0TRE","BOX",idtmed[kAir],ptopref,3); //air: wrapped  radiator
  TGeoVolume *topref = gGeoManager->GetVolume("0TRE");
  TVirtualMC::GetMC()->Gsvolu ("0REG", "BOX", idtmed[kOpGlassCathode], preg, 3); 
  TGeoVolume *cat = gGeoManager->GetVolume("0REG");
  TVirtualMC::GetMC()->Gsvolu("0MCP","BOX",idtmed[kGlass],pmcp,3); //glass
  TGeoVolume *mcp = gGeoManager->GetVolume("0MCP");
  TVirtualMC::GetMC()->Gsvolu("0RFV","BOX",idtmed[kOpAir],prfv,3); //optical Air vertical
  TGeoVolume *rfv = gGeoManager->GetVolume("0RFV");
  TVirtualMC::GetMC()->Gsvolu("0RFH","BOX",idtmed[kOpAir],prfh,3); //optical Air horizontal
  TGeoVolume *rfh = gGeoManager->GetVolume("0RFH");
  
   TVirtualMC::GetMC()->Gsvolu("0PAL","BOX",idtmed[kAl],pal,3); // 5mmAl top on the radiator
  TGeoVolume *altop = gGeoManager->GetVolume("0PAL");
  
  Double_t thet = TMath::ATan(dP/crad);
  Double_t rat = TMath::Tan(thet)/2.0;
    
  //Al housing definition  
  Double_t mgon[16];
    
  mgon[0]  = -45;
  mgon[1]  = 360.0;
  mgon[2]  = 4;
  mgon[3]  = 4;
  
  z = -pinstart[2] + 2*pal[2];
  mgon[4]  = z;
  mgon[5]  = 2*ptop[0] + preg[2];
  mgon[6]  = dP+rat*z*4/3;
  
  z = -pinstart[2] + 2*pal[2] + 2*ptopref[2];
  mgon[7]  = z;
  mgon[8]  = mgon[5];
  mgon[9]  = dP+z*rat;
  mgon[10] = z;
  mgon[11] = pmcp[0] + preg[2];
  mgon[12] = mgon[9];
  
  z = -pinstart[2] + 2*pal[2] + 2*ptopref[2] + 2*preg[2] + 2*pmcp[2];
  mgon[13] = z;
  mgon[14] = mgon[11];
  mgon[15] = dP+z*rat*pmcp[2]*9/10;  
  
  TVirtualMC::GetMC()->Gsvolu("0SUP","PGON",idtmed[kAl], mgon, 16); //Al Housing for Support Structure
  TGeoVolume *alsup = gGeoManager->GetVolume("0SUP");


  //wrapped radiator +  reflectiong layers 
  Int_t ntops=0, nrfvs=0, nrfhs=0;
  Float_t xin=0, yin=0, xinv=0, yinv=0,xinh=0,yinh=0;
  x=y=z=0;
  topref->AddNode(top, 1, new TGeoTranslation(0,0,0) );
  xinv = -ptop[0] - prfv[0];
  topref->AddNode(rfv, 1, new TGeoTranslation(xinv,0,0) );
  printf(" GEOGEO  refv %f ,  0,0 \n",xinv);
  xinv = ptop[0] + prfv[0];
  topref->AddNode(rfv, 2, new TGeoTranslation(xinv,0,0) );
  printf(" GEOGEO  refv %f ,  0,0 \n",xinv);
  yinv = -ptop[1] - prfh[1];
  topref->AddNode(rfh, 1, new TGeoTranslation(0,yinv,0) );
  printf(" GEOGEO  refh  ,  0, %f, 0 \n",yinv);
  yinv = ptop[1] + prfh[1];
  topref->AddNode(rfh, 2, new TGeoTranslation(0,yinv,0) );
  
  //container for radiator, cathod 
  for (Int_t ix=0; ix<2; ix++) {
    xin = - pinstart[0] + 0.3 + (ix+0.5)*2*ptopref[0] ;
    for (Int_t iy=0; iy<2 ; iy++) {
      z = - pinstart[2] + 2*pal[2] + ptopref[2];
      yin = - pinstart[1] + 0.3 + (iy+0.5)*2*ptopref[1];
      ntops++;
      ins->AddNode(topref, ntops, new TGeoTranslation(xin,yin,z) );
      printf(" 0TOP  full %i x %f y %f z %f \n", ntops, xin, yin, z);
      z = -pinstart[2]   + 2*pal[2] + 2 * ptopref[2] + preg[2];
      ins->AddNode(cat, ntops, new TGeoTranslation(xin, yin, z) );
      // cat->Print();
      printf(" GEOGEO  CATHOD x=%f , y= %f z= %f num  %i\n", xin, yin, z, ntops);
    }
  }
  //Al top
  z=-pinstart[2] + pal[2];
  ins->AddNode(altop, 1 , new TGeoTranslation(0,0,z) );
  
  // MCP
  z=-pinstart[2] + 2*pal[2] + 2*ptopref[2] + 2*preg[2] + pmcp[2];
  //   z=-pinstart[2] + 2*ptopref[2] + preg[2];
  ins->AddNode(mcp, 1 , new TGeoTranslation(0,0,z) );
  
  // Al Housing for Support Structure
  ins->AddNode(alsup,1);
}
  
//_________________________________________
void AliFITv7::SetVZEROGeo(TGeoVolume *alice)
{ 
  //V0+
  cout<<" AliFITv7::SetVZEROGeo "<<alice<<endl;
  const int kV0PlusColorSci   = 5;
  TGeoMedium *medV0PlusSci = gGeoManager->GetMedium("FIT_V0PlusSci");
  Double_t Pi = TMath::Pi();
  Double_t Sin22p5    = TMath::Sin(Pi/8.); 
  Double_t Cos22p5    = TMath::Cos(Pi/8.); 
  Double_t v0PlusPts[16];
  
    //Defining the master volume for V0Plus
    TGeoVolume *v0Plus = new TGeoVolumeAssembly("V0LE");
        
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
	v0PlusPts[0+8*i] = fV0PlusR0*Cos22p5-fV0PlusFraWd;
	v0PlusPts[1+8*i] = (fV0PlusR0-fV0PlusFraWd)*Sin22p5-fV0PlusFraWd;
	v0PlusPts[2+8*i] = (fV0PlusR0-fV0PlusFraWd/2.)*Cos22p5-fV0PlusFraWd;
	v0PlusPts[3+8*i] = (fV0PlusR0-fV0PlusFraWd/2.)*Sin22p5;
	v0PlusPts[4+8*i] = (fV0PlusR5+fV0PlusFraWd/2.)*Cos22p5+fV0PlusFraWd;
	v0PlusPts[5+8*i] = (fV0PlusR5+fV0PlusFraWd/2.)*Sin22p5+2.*fV0PlusFraWd;
	v0PlusPts[6+8*i] = (fV0PlusR5+fV0PlusFraWd)*Cos22p5+fV0PlusFraWd;
	v0PlusPts[7+8*i] = fV0PlusR5*Sin22p5+fV0PlusFraWd;
      }
    new TGeoArb8("sV0PlusCha2", fV0PlusSciWd/2.+2.*fV0PlusFraWd, v0PlusPts);
    new TGeoCompositeShape("sV0PlusCha","sV0PlusCha1+sV0PlusCha2");

    //Sensitive scintillator
    new TGeoTubeSeg( "sV0PlusR1b", fV0PlusR0+fV0PlusFraWd/2., fV0PlusR1-fV0PlusFraWd/2., fV0PlusSciWd/2., 0, 22.5);
    new TGeoTubeSeg( "sV0PlusR2b", fV0PlusR1+fV0PlusFraWd/2., fV0PlusR2-fV0PlusFraWd/2., fV0PlusSciWd/2., 0, 22.5);
    new TGeoTubeSeg( "sV0PlusR3b", fV0PlusR2+fV0PlusFraWd/2., fV0PlusR3-fV0PlusFraWd/2., fV0PlusSciWd/2., 0, 22.5);
    new TGeoTubeSeg( "sV0PlusR4b", fV0PlusR3+fV0PlusFraWd/2., fV0PlusR4-fV0PlusFraWd/2., fV0PlusSciWd/2., 0, 22.5);
    new TGeoTubeSeg( "sV0PlusR5b", fV0PlusR4+fV0PlusFraWd/2., fV0PlusR5-fV0PlusFraWd/2., fV0PlusSciWd/2., 0, 22.5);
    TGeoCompositeShape *sV0PlusR1 = new TGeoCompositeShape("sV0PlusR1","sV0PlusR1b-sV0PlusCha");
    TGeoCompositeShape *sV0PlusR2 = new TGeoCompositeShape("sV0PlusR2","sV0PlusR2b-sV0PlusCha");
    TGeoCompositeShape *sV0PlusR3 = new TGeoCompositeShape("sV0PlusR3","sV0PlusR3b-sV0PlusCha");
    TGeoCompositeShape *sV0PlusR4 = new TGeoCompositeShape("sV0PlusR4","sV0PlusR4b-sV0PlusCha");
    TGeoCompositeShape *sV0PlusR5 = new TGeoCompositeShape("sV0PlusR5","sV0PlusR5b-sV0PlusCha");
    

    /// Definition sector 1
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
 
    TGeoRotation *RotSec1 = new TGeoRotation("RotSec1", 90., 0*22.5, 90., 90.+0*22.5, 0., 0.);
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
 
    TGeoRotation *RotSec2 = new TGeoRotation("RotSec2", 90., 1*22.5, 90., 90.+1*22.5, 0., 0.);
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
 
    TGeoRotation *RotSec3 = new TGeoRotation("RotSec3", 90., 2*22.5, 90., 90.+2*22.5, 0., 0.);
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
 
    TGeoRotation *RotSec4 = new TGeoRotation("RotSec4", 90., 3*22.5, 90., 90.+3*22.5, 0., 0.);
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
 
    TGeoRotation *RotSec5 = new TGeoRotation("RotSec5", 90., 4*22.5, 90., 90.+4*22.5, 0., 0.);
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
 
    TGeoRotation *RotSec6 = new TGeoRotation("RotSec6", 90., 5*22.5, 90., 90.+5*22.5, 0., 0.);
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
 
    TGeoRotation *RotSec7 = new TGeoRotation("RotSec7", 90., 6*22.5, 90., 90.+6*22.5, 0., 0.);
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
 
    TGeoRotation *RotSec8 = new TGeoRotation("RotSec8", 90., 7*22.5, 90., 90.+7*22.5, 0., 0.);
    v0Plus->AddNode(v0PlusSec8,7+1,RotSec8);

    /// Definition sector 9
    TGeoVolume *v0Plus1Sec9 = new TGeoVolume("V0Plus1Sec9",sV0PlusR1,medV0PlusSci); 
    TGeoVolume *v0Plus2Sec9 = new TGeoVolume("V0Plus2Sec9",sV0PlusR2,medV0PlusSci);
    TGeoVolume *v0Plus3Sec9 = new TGeoVolume("V0Plus3Sec9",sV0PlusR3,medV0PlusSci);
    TGeoVolume *v0Plus4Sec9 = new TGeoVolume("V0Plus4Sec9",sV0PlusR4,medV0PlusSci);
    TGeoVolume *v0Plus5Sec9 = new TGeoVolume("V0Plus5Sec9",sV0PlusR5,medV0PlusSci); 
    v0Plus1Sec9->SetLineColor(kV0PlusColorSci); 
    v0Plus2Sec9->SetLineColor(kV0PlusColorSci);
    v0Plus3Sec9->SetLineColor(kV0PlusColorSci); 
    v0Plus4Sec9->SetLineColor(kV0PlusColorSci);
    v0Plus5Sec9->SetLineColor(kV0PlusColorSci);
    TGeoVolume *v0PlusSciSec9 = new TGeoVolumeAssembly("V0PlusSciSec9");
    v0PlusSciSec9->AddNode(v0Plus1Sec9,1);
    v0PlusSciSec9->AddNode(v0Plus2Sec9,1);
    v0PlusSciSec9->AddNode(v0Plus3Sec9,1);
    v0PlusSciSec9->AddNode(v0Plus4Sec9,1);      
    v0PlusSciSec9->AddNode(v0Plus5Sec9,1);
    TGeoVolume *v0PlusSec9 = new TGeoVolumeAssembly("V0PlusSec9");
    v0PlusSec9->AddNode(v0PlusSciSec9,1);
 
    TGeoRotation *RotSec9 = new TGeoRotation("RotSec9", 90., 8*22.5, 90., 90.+8*22.5, 0., 0.);
    v0Plus->AddNode(v0PlusSec9,8+1,RotSec9);

    /// Definition sector 10
    TGeoVolume *v0Plus1Sec10 = new TGeoVolume("V0Plus1Sec10",sV0PlusR1,medV0PlusSci); 
    TGeoVolume *v0Plus2Sec10 = new TGeoVolume("V0Plus2Sec10",sV0PlusR2,medV0PlusSci);
    TGeoVolume *v0Plus3Sec10 = new TGeoVolume("V0Plus3Sec10",sV0PlusR3,medV0PlusSci);
    TGeoVolume *v0Plus4Sec10 = new TGeoVolume("V0Plus4Sec10",sV0PlusR4,medV0PlusSci);
    TGeoVolume *v0Plus5Sec10 = new TGeoVolume("V0Plus5Sec10",sV0PlusR5,medV0PlusSci); 
    v0Plus1Sec10->SetLineColor(kV0PlusColorSci); 
    v0Plus2Sec10->SetLineColor(kV0PlusColorSci);
    v0Plus3Sec10->SetLineColor(kV0PlusColorSci); 
    v0Plus4Sec10->SetLineColor(kV0PlusColorSci);
    v0Plus5Sec10->SetLineColor(kV0PlusColorSci);
    TGeoVolume *v0PlusSciSec10 = new TGeoVolumeAssembly("V0PlusSciSec10");
    v0PlusSciSec10->AddNode(v0Plus1Sec10,1);
    v0PlusSciSec10->AddNode(v0Plus2Sec10,1);
    v0PlusSciSec10->AddNode(v0Plus3Sec10,1);
    v0PlusSciSec10->AddNode(v0Plus4Sec10,1);      
    v0PlusSciSec10->AddNode(v0Plus5Sec10,1);
    TGeoVolume *v0PlusSec10 = new TGeoVolumeAssembly("V0PlusSec10");
    v0PlusSec10->AddNode(v0PlusSciSec10,1);
 
    TGeoRotation *RotSec10 = new TGeoRotation("RotSec10", 90., 9*22.5, 90., 90.+9*22.5, 0., 0.);
    v0Plus->AddNode(v0PlusSec10,9+1,RotSec10);

    /// Definition sector 11
    TGeoVolume *v0Plus1Sec11 = new TGeoVolume("V0Plus1Sec11",sV0PlusR1,medV0PlusSci); 
    TGeoVolume *v0Plus2Sec11 = new TGeoVolume("V0Plus2Sec11",sV0PlusR2,medV0PlusSci);
    TGeoVolume *v0Plus3Sec11 = new TGeoVolume("V0Plus3Sec11",sV0PlusR3,medV0PlusSci);
    TGeoVolume *v0Plus4Sec11 = new TGeoVolume("V0Plus4Sec11",sV0PlusR4,medV0PlusSci);
    TGeoVolume *v0Plus5Sec11 = new TGeoVolume("V0Plus5Sec11",sV0PlusR5,medV0PlusSci); 
    v0Plus1Sec11->SetLineColor(kV0PlusColorSci); 
    v0Plus2Sec11->SetLineColor(kV0PlusColorSci);
    v0Plus3Sec11->SetLineColor(kV0PlusColorSci); 
    v0Plus4Sec11->SetLineColor(kV0PlusColorSci);
    v0Plus5Sec11->SetLineColor(kV0PlusColorSci);
    TGeoVolume *v0PlusSciSec11 = new TGeoVolumeAssembly("V0PlusSciSec11");
    v0PlusSciSec11->AddNode(v0Plus1Sec11,1);
    v0PlusSciSec11->AddNode(v0Plus2Sec11,1);
    v0PlusSciSec11->AddNode(v0Plus3Sec11,1);
    v0PlusSciSec11->AddNode(v0Plus4Sec11,1);      
    v0PlusSciSec11->AddNode(v0Plus5Sec11,1);
    TGeoVolume *v0PlusSec11 = new TGeoVolumeAssembly("V0PlusSec11");
    v0PlusSec11->AddNode(v0PlusSciSec11,1);
 
    TGeoRotation *RotSec11 = new TGeoRotation("RotSec11", 90., 10*22.5, 90., 90.+10*22.5, 0., 0.);
    v0Plus->AddNode(v0PlusSec11,10+1,RotSec11);

    /// Definition sector 12
    TGeoVolume *v0Plus1Sec12 = new TGeoVolume("V0Plus1Sec12",sV0PlusR1,medV0PlusSci); 
    TGeoVolume *v0Plus2Sec12 = new TGeoVolume("V0Plus2Sec12",sV0PlusR2,medV0PlusSci);
    TGeoVolume *v0Plus3Sec12 = new TGeoVolume("V0Plus3Sec12",sV0PlusR3,medV0PlusSci);
    TGeoVolume *v0Plus4Sec12 = new TGeoVolume("V0Plus4Sec12",sV0PlusR4,medV0PlusSci);
    TGeoVolume *v0Plus5Sec12 = new TGeoVolume("V0Plus5Sec12",sV0PlusR5,medV0PlusSci); 
    v0Plus1Sec12->SetLineColor(kV0PlusColorSci); 
    v0Plus2Sec12->SetLineColor(kV0PlusColorSci);
    v0Plus3Sec12->SetLineColor(kV0PlusColorSci); 
    v0Plus4Sec12->SetLineColor(kV0PlusColorSci);
    v0Plus5Sec12->SetLineColor(kV0PlusColorSci);
    TGeoVolume *v0PlusSciSec12 = new TGeoVolumeAssembly("V0PlusSciSec12");
    v0PlusSciSec12->AddNode(v0Plus1Sec12,1);
    v0PlusSciSec12->AddNode(v0Plus2Sec12,1);
    v0PlusSciSec12->AddNode(v0Plus3Sec12,1);
    v0PlusSciSec12->AddNode(v0Plus4Sec12,1);      
    v0PlusSciSec12->AddNode(v0Plus5Sec12,1);
    TGeoVolume *v0PlusSec12 = new TGeoVolumeAssembly("V0PlusSec12");
    v0PlusSec12->AddNode(v0PlusSciSec12,1);
 
    TGeoRotation *RotSec12 = new TGeoRotation("RotSec12", 90., 11*22.5, 90., 90.+11*22.5, 0., 0.);
    v0Plus->AddNode(v0PlusSec12,11+1,RotSec12);

    /// Definition sector 13
    TGeoVolume *v0Plus1Sec13 = new TGeoVolume("V0Plus1Sec13",sV0PlusR1,medV0PlusSci); 
    TGeoVolume *v0Plus2Sec13 = new TGeoVolume("V0Plus2Sec13",sV0PlusR2,medV0PlusSci);
    TGeoVolume *v0Plus3Sec13 = new TGeoVolume("V0Plus3Sec13",sV0PlusR3,medV0PlusSci);
    TGeoVolume *v0Plus4Sec13 = new TGeoVolume("V0Plus4Sec13",sV0PlusR4,medV0PlusSci);
    TGeoVolume *v0Plus5Sec13 = new TGeoVolume("V0Plus5Sec13",sV0PlusR5,medV0PlusSci); 
    v0Plus1Sec13->SetLineColor(kV0PlusColorSci); 
    v0Plus2Sec13->SetLineColor(kV0PlusColorSci);
    v0Plus3Sec13->SetLineColor(kV0PlusColorSci); 
    v0Plus4Sec13->SetLineColor(kV0PlusColorSci);
    v0Plus5Sec13->SetLineColor(kV0PlusColorSci);
    TGeoVolume *v0PlusSciSec13 = new TGeoVolumeAssembly("V0PlusSciSec13");
    v0PlusSciSec13->AddNode(v0Plus1Sec13,1);
    v0PlusSciSec13->AddNode(v0Plus2Sec13,1);
    v0PlusSciSec13->AddNode(v0Plus3Sec13,1);
    v0PlusSciSec13->AddNode(v0Plus4Sec13,1);      
    v0PlusSciSec13->AddNode(v0Plus5Sec13,1);
    TGeoVolume *v0PlusSec13 = new TGeoVolumeAssembly("V0PlusSec13");
    v0PlusSec13->AddNode(v0PlusSciSec13,1);
 
    TGeoRotation *RotSec13 = new TGeoRotation("RotSec13", 90., 12*22.5, 90., 90.+12*22.5, 0., 0.);
    v0Plus->AddNode(v0PlusSec13,12+1,RotSec13);

    /// Definition sector 14
    TGeoVolume *v0Plus1Sec14 = new TGeoVolume("V0Plus1Sec14",sV0PlusR1,medV0PlusSci); 
    TGeoVolume *v0Plus2Sec14 = new TGeoVolume("V0Plus2Sec14",sV0PlusR2,medV0PlusSci);
    TGeoVolume *v0Plus3Sec14 = new TGeoVolume("V0Plus3Sec14",sV0PlusR3,medV0PlusSci);
    TGeoVolume *v0Plus4Sec14 = new TGeoVolume("V0Plus4Sec14",sV0PlusR4,medV0PlusSci);
    TGeoVolume *v0Plus5Sec14 = new TGeoVolume("V0Plus5Sec14",sV0PlusR5,medV0PlusSci); 
    v0Plus1Sec14->SetLineColor(kV0PlusColorSci); 
    v0Plus2Sec14->SetLineColor(kV0PlusColorSci);
    v0Plus3Sec14->SetLineColor(kV0PlusColorSci); 
    v0Plus4Sec14->SetLineColor(kV0PlusColorSci);
    v0Plus5Sec14->SetLineColor(kV0PlusColorSci);
    TGeoVolume *v0PlusSciSec14 = new TGeoVolumeAssembly("V0PlusSciSec14");
    v0PlusSciSec14->AddNode(v0Plus1Sec14,1);
    v0PlusSciSec14->AddNode(v0Plus2Sec14,1);
    v0PlusSciSec14->AddNode(v0Plus3Sec14,1);
    v0PlusSciSec14->AddNode(v0Plus4Sec14,1);      
    v0PlusSciSec14->AddNode(v0Plus5Sec14,1);
    TGeoVolume *v0PlusSec14 = new TGeoVolumeAssembly("V0PlusSec14");
    v0PlusSec14->AddNode(v0PlusSciSec14,1);
 
    TGeoRotation *RotSec14 = new TGeoRotation("RotSec14", 90., 13*22.5, 90., 90.+13*22.5, 0., 0.);
    v0Plus->AddNode(v0PlusSec14,13+1,RotSec14);

    /// Definition sector 15
    TGeoVolume *v0Plus1Sec15 = new TGeoVolume("V0Plus1Sec15",sV0PlusR1,medV0PlusSci); 
    TGeoVolume *v0Plus2Sec15 = new TGeoVolume("V0Plus2Sec15",sV0PlusR2,medV0PlusSci);
    TGeoVolume *v0Plus3Sec15 = new TGeoVolume("V0Plus3Sec15",sV0PlusR3,medV0PlusSci);
    TGeoVolume *v0Plus4Sec15 = new TGeoVolume("V0Plus4Sec15",sV0PlusR4,medV0PlusSci);
    TGeoVolume *v0Plus5Sec15 = new TGeoVolume("V0Plus5Sec15",sV0PlusR5,medV0PlusSci); 
    v0Plus1Sec15->SetLineColor(kV0PlusColorSci); 
    v0Plus2Sec15->SetLineColor(kV0PlusColorSci);
    v0Plus3Sec15->SetLineColor(kV0PlusColorSci); 
    v0Plus4Sec15->SetLineColor(kV0PlusColorSci);
    v0Plus5Sec15->SetLineColor(kV0PlusColorSci);
    TGeoVolume *v0PlusSciSec15 = new TGeoVolumeAssembly("V0PlusSciSec15");
    v0PlusSciSec15->AddNode(v0Plus1Sec15,1);
    v0PlusSciSec15->AddNode(v0Plus2Sec15,1);
    v0PlusSciSec15->AddNode(v0Plus3Sec15,1);
    v0PlusSciSec15->AddNode(v0Plus4Sec15,1);      
    v0PlusSciSec15->AddNode(v0Plus5Sec15,1);
    TGeoVolume *v0PlusSec15 = new TGeoVolumeAssembly("V0PlusSec15");
    v0PlusSec15->AddNode(v0PlusSciSec15,1);
 
    TGeoRotation *RotSec15 = new TGeoRotation("RotSec15", 90., 14*22.5, 90., 90.+14*22.5, 0., 0.);
    v0Plus->AddNode(v0PlusSec15,14+1,RotSec15);

    /// Definition sector 16
    TGeoVolume *v0Plus1Sec16 = new TGeoVolume("V0Plus1Sec16",sV0PlusR1,medV0PlusSci); 
    TGeoVolume *v0Plus2Sec16 = new TGeoVolume("V0Plus2Sec16",sV0PlusR2,medV0PlusSci);
    TGeoVolume *v0Plus3Sec16 = new TGeoVolume("V0Plus3Sec16",sV0PlusR3,medV0PlusSci);
    TGeoVolume *v0Plus4Sec16 = new TGeoVolume("V0Plus4Sec16",sV0PlusR4,medV0PlusSci);
    TGeoVolume *v0Plus5Sec16 = new TGeoVolume("V0Plus5Sec16",sV0PlusR5,medV0PlusSci); 
    v0Plus1Sec16->SetLineColor(kV0PlusColorSci); 
    v0Plus2Sec16->SetLineColor(kV0PlusColorSci);
    v0Plus3Sec16->SetLineColor(kV0PlusColorSci); 
    v0Plus4Sec16->SetLineColor(kV0PlusColorSci);
    v0Plus5Sec16->SetLineColor(kV0PlusColorSci);
    TGeoVolume *v0PlusSciSec16 = new TGeoVolumeAssembly("V0PlusSciSec16");
    v0PlusSciSec16->AddNode(v0Plus1Sec16,1);
    v0PlusSciSec16->AddNode(v0Plus2Sec16,1);
    v0PlusSciSec16->AddNode(v0Plus3Sec16,1);
    v0PlusSciSec16->AddNode(v0Plus4Sec16,1);      
    v0PlusSciSec16->AddNode(v0Plus5Sec16,1);
    TGeoVolume *v0PlusSec16 = new TGeoVolumeAssembly("V0PlusSec16");
    v0PlusSec16->AddNode(v0PlusSciSec16,1);
 
    TGeoRotation *RotSec16 = new TGeoRotation("RotSec16", 90., 15*22.5, 90., 90.+15*22.5, 0., 0.);
    v0Plus->AddNode(v0PlusSec16,15+1,RotSec16);

    alice->AddNode(v0Plus,1,new TGeoTranslation(0, 0, fV0PlusZposition));
 
}    
//------------------------------------------------------------------------
void AliFITv7::AddAlignableVolumes() const
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
void AliFITv7::CreateMaterials()
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
  AliMaterial(11, "Aliminium$", 26.98, 13.0, 2.7, 8.9,999); 
 
  AliMedium(1, "Air$", 2, 0, isxfld, sxmgmx, 10., .1, 1., .003, .003);
  AliMedium(3, "Vacuum$", 1, 0, isxfld, sxmgmx, 10., .01, .1, .003, .003);
  AliMedium(6, "Glass$", 4, 0, isxfld, sxmgmx, 10., .01, .1, .003, .003);
  
  AliMedium(7, "OpAir$", 2, 0, isxfld, sxmgmx, 10., .1, 1., .003, .003);
  
  AliMedium(15, "Aluminium$", 11, 0, isxfld, sxmgmx, 10., .01, 1., .003, .003);  
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
void AliFITv7::DefineOpticalProperties()
{

// Path of the optical properties input file
  TString optPropPath = "$(ALICE_ROOT)/FIT/sim/quartzOptProperties.txt";
  optPropPath = gSystem->ExpandPathName(optPropPath.Data()); // Expand $(ALICE_ROOT) into real system path

// Optical properties definition.
  Int_t *idtmed = fIdtmed->GetArray();

// Prepare pointers for arrays read from the input file
  Float_t *aPckov=NULL;
  Double_t *dPckov=NULL;
  Float_t *aAbsSiO2=NULL;
  Float_t *rindexSiO2=NULL;
  Float_t *qeff = NULL;
  Int_t kNbins=0;
  ReadOptProperties(optPropPath.Data(), &aPckov, &dPckov, &aAbsSiO2, &rindexSiO2, &qeff, kNbins);
  // set QE
   fPMTeff = new TGraph(kNbins,aPckov,qeff);

// Prepare pointers for arrays with constant and hardcoded values (independent on wavelength)
  Float_t *efficAll=NULL;
  Float_t *rindexAir=NULL;
  Float_t *absorAir=NULL;
  Float_t *rindexCathodeNext=NULL;
  Float_t *absorbCathodeNext=NULL;
  Double_t *efficMet=NULL;
  Double_t *aReflMet=NULL;
  FillOtherOptProperties(&efficAll, &rindexAir, &absorAir, &rindexCathodeNext,
    &absorbCathodeNext, &efficMet, &aReflMet, kNbins);
  
  TVirtualMC::GetMC()->SetCerenkov (idtmed[kOpGlass],        kNbins, aPckov, aAbsSiO2, efficAll, rindexSiO2);
  // TVirtualMC::GetMC()->SetCerenkov (idtmed[kOpGlassCathode], kNbins, aPckov, aAbsSiO2, effCathode, rindexSiO2);
  TVirtualMC::GetMC()->SetCerenkov (idtmed[kOpGlassCathode], kNbins, aPckov, aAbsSiO2, efficAll, rindexSiO2);
  // TVirtualMC::GetMC()->SetCerenkov (idtmed[kOpAir],       kNbins, aPckov, absorAir ,efficAll, rindexAir);
  // TVirtualMC::GetMC()->SetCerenkov (idtmed[kOpAirNext],   kNbins, aPckov, absorbCathodeNext, efficAll, rindexCathodeNext);

//Define a border for radiator optical properties
// TODO: Maciek: The following 3 lines just generate warnings and do nothing else - could be deleted
  TVirtualMC::GetMC()->DefineOpSurface("surfRd", kUnified /*kGlisur*/,kDielectric_metal,kPolished, 0.);
  TVirtualMC::GetMC()->SetMaterialProperty("surfRd", "EFFICIENCY", kNbins, dPckov, efficMet);
  TVirtualMC::GetMC()->SetMaterialProperty("surfRd", "REFLECTIVITY", kNbins, dPckov, aReflMet);

  DeleteOptPropertiesArr(&aPckov, &dPckov, &aAbsSiO2, &rindexSiO2, &efficAll, &rindexAir, &absorAir, &rindexCathodeNext, &absorbCathodeNext, &efficMet, &aReflMet);
}


//-------------------------------------------------------------------
void AliFITv7::Init()
{
// Initialises version 0 of the Forward Multiplicity Detector
//
  AliFIT::Init();
    fIdSens1=TVirtualMC::GetMC()->VolId("0REG");
  // fIdSens1=TVirtualMC::GetMC()->VolId("0TOP");

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

void AliFITv7::StepManager()
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
  // printf("T0 :::volumes %i %s \n", id, TVirtualMC::GetMC()->CurrentVolName() ); 
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
      //      printf("T0 :::volumes pmt %i mcp %i vol %i x %f y %f z %f particle %f all \n",  vol[0], vol[1],  vol[2], hits[0], hits[1], hits[2], hits[4]);
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

      if ( (id == fIdV0Plus[0][0]) || (id == fIdV0Plus[1][0]) || (id == fIdV0Plus[2][0]) || (id == fIdV0Plus[3][0]) || (id == fIdV0Plus[4][0]) || (id == fIdV0Plus[5][0]) || (id == fIdV0Plus[6][0]) || (id == fIdV0Plus[7][0]) || (id == fIdV0Plus[8][0]) || (id == fIdV0Plus[9][0]) || (id == fIdV0Plus[10][0]) || (id == fIdV0Plus[11][0]) || (id == fIdV0Plus[12][0]) || (id == fIdV0Plus[13][0]) || (id == fIdV0Plus[14][0]) || (id == fIdV0Plus[15][0]) ) RingNumber = 1;

      else if ( (id == fIdV0Plus[0][1]) || (id == fIdV0Plus[1][1]) || (id == fIdV0Plus[2][1]) || (id == fIdV0Plus[3][1]) || (id == fIdV0Plus[4][1]) || (id == fIdV0Plus[5][1]) || (id == fIdV0Plus[6][1]) || (id == fIdV0Plus[7][1]) || (id == fIdV0Plus[8][1]) || (id == fIdV0Plus[9][1]) || (id == fIdV0Plus[10][1]) || (id == fIdV0Plus[11][1]) || (id == fIdV0Plus[12][1]) || (id == fIdV0Plus[13][1]) || (id == fIdV0Plus[14][1]) || (id == fIdV0Plus[15][1]) ) RingNumber = 2;

      else if ( (id == fIdV0Plus[0][2]) || (id == fIdV0Plus[1][2]) || (id == fIdV0Plus[2][2]) || (id == fIdV0Plus[3][2]) || (id == fIdV0Plus[4][2]) || (id == fIdV0Plus[5][2]) || (id == fIdV0Plus[6][2]) || (id == fIdV0Plus[7][2]) || (id == fIdV0Plus[8][2]) || (id == fIdV0Plus[9][2]) || (id == fIdV0Plus[10][2]) || (id == fIdV0Plus[11][2]) || (id == fIdV0Plus[12][2]) || (id == fIdV0Plus[13][2]) || (id == fIdV0Plus[14][2]) || (id == fIdV0Plus[15][2]) ) RingNumber = 3;

      else if ( (id == fIdV0Plus[0][3]) || (id == fIdV0Plus[1][3]) || (id == fIdV0Plus[2][3]) || (id == fIdV0Plus[3][3]) || (id == fIdV0Plus[4][3]) || (id == fIdV0Plus[5][3]) || (id == fIdV0Plus[6][3]) || (id == fIdV0Plus[7][3]) || (id == fIdV0Plus[8][3]) || (id == fIdV0Plus[9][3]) || (id == fIdV0Plus[10][3]) || (id == fIdV0Plus[11][3]) || (id == fIdV0Plus[12][3]) || (id == fIdV0Plus[13][3]) || (id == fIdV0Plus[14][3]) || (id == fIdV0Plus[15][3]) ) RingNumber = 4;

      else if ( (id == fIdV0Plus[0][4]) || (id == fIdV0Plus[1][4]) || (id == fIdV0Plus[2][4]) || (id == fIdV0Plus[3][4]) || (id == fIdV0Plus[4][4]) || (id == fIdV0Plus[5][4]) || (id == fIdV0Plus[6][4]) || (id == fIdV0Plus[7][4]) || (id == fIdV0Plus[8][4]) || (id == fIdV0Plus[9][4]) || (id == fIdV0Plus[10][4]) || (id == fIdV0Plus[11][4]) || (id == fIdV0Plus[12][4]) || (id == fIdV0Plus[13][4]) || (id == fIdV0Plus[14][4]) || (id == fIdV0Plus[15][4]) ) RingNumber = 5;

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
Bool_t AliFITv7::RegisterPhotoE(Double_t energy)
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

//-------------------------------------------------------------------------
Int_t AliFITv7::GetCellId(Int_t *vol)
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
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus1Sec9")) fCellId = 9;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus1Sec10")) fCellId = 10;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus1Sec11")) fCellId = 11;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus1Sec12")) fCellId = 12;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus1Sec13")) fCellId = 13;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus1Sec14")) fCellId = 14;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus1Sec15")) fCellId = 15;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus1Sec16")) fCellId = 16;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus2Sec1")) fCellId = 17;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus2Sec2")) fCellId = 18;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus2Sec3")) fCellId = 19;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus2Sec4")) fCellId = 20;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus2Sec5")) fCellId = 21;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus2Sec6")) fCellId = 22;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus2Sec7")) fCellId = 23;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus2Sec8")) fCellId = 24;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus2Sec9")) fCellId = 25;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus2Sec10")) fCellId = 26;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus2Sec11")) fCellId = 27;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus2Sec12")) fCellId = 28;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus2Sec13")) fCellId = 29;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus2Sec14")) fCellId = 30;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus2Sec15")) fCellId = 31;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus2Sec16")) fCellId = 32;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus3Sec1")) fCellId = 33;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus3Sec2")) fCellId = 34;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus3Sec3")) fCellId = 35;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus3Sec4")) fCellId = 36;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus3Sec5")) fCellId = 37;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus3Sec6")) fCellId = 38;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus3Sec7")) fCellId = 39;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus3Sec8")) fCellId = 40;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus3Sec9")) fCellId = 41;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus3Sec10")) fCellId = 42;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus3Sec11")) fCellId = 43;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus3Sec12")) fCellId = 44;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus3Sec13")) fCellId = 45;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus3Sec14")) fCellId = 46;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus3Sec15")) fCellId = 47;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus3Sec16")) fCellId = 48;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus4Sec1")) fCellId = 49;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus4Sec2")) fCellId = 50;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus4Sec3")) fCellId = 51;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus4Sec4")) fCellId = 52;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus4Sec5")) fCellId = 53;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus4Sec6")) fCellId = 54;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus4Sec7")) fCellId = 55;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus4Sec8")) fCellId = 56;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus4Sec9")) fCellId = 57;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus4Sec10")) fCellId = 58;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus4Sec11")) fCellId = 59;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus4Sec12")) fCellId = 60;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus4Sec13")) fCellId = 61;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus4Sec14")) fCellId = 62;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus4Sec15")) fCellId = 63;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus4Sec16")) fCellId = 64;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus5Sec1")) fCellId = 65;	
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus5Sec2")) fCellId = 66;	
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus5Sec3")) fCellId = 67;	
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus5Sec4")) fCellId = 68;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus5Sec5")) fCellId = 69;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus5Sec6")) fCellId = 70;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus5Sec7")) fCellId = 71;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus5Sec8")) fCellId = 72;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus5Sec9")) fCellId = 73;	
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus5Sec10")) fCellId = 74;	
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus5Sec11")) fCellId = 75;	
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus5Sec12")) fCellId = 76;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus5Sec13")) fCellId = 77;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus5Sec14")) fCellId = 78;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus5Sec15")) fCellId = 79;
  if (gMC->CurrentVolID(vol[2]) == gMC->VolId("V0Plus5Sec16")) fCellId = 80;
 
  return fCellId;

}


//-----------------------------------------------------------------

Int_t AliFITv7::ReadOptProperties(const std::string filePath, Float_t **e, Double_t **de,
				  Float_t **abs, Float_t **n, Float_t **qe, Int_t &kNbins) const{
  std::ifstream infile;
  infile.open(filePath.c_str());

  // Check if file is opened correctly
  if(infile.fail()==true){
          AliFatal(Form("Error opening ascii file: %s", filePath.c_str()));
    return -1;
  }

  std::string comment; // dummy, used just to read 4 first lines and move the cursor to the 5th, otherwise unused
  if(!getline(infile,comment)){ // first comment line
          AliFatal(Form("Error opening ascii file (it is probably a folder!): %s", filePath.c_str()));
    return -2;
  }
  getline(infile,comment); // 2nd comment line

  // Get number of elements required for the array
  infile >> kNbins;
  if(kNbins<0 || kNbins>1e4){
          AliFatal(Form("Input arraySize out of range 0..1e4: %i. Check input file: %s", kNbins, filePath.c_str()));
    return -4;
  }

  // Allocate memory required for arrays
  *e = new Float_t[kNbins];
  *de = new Double_t[kNbins];
  *abs = new Float_t[kNbins];
  *n = new Float_t[kNbins];
  *qe = new Float_t[kNbins];

  getline(infile,comment); // finish 3rd line after the nEntries are read
  getline(infile,comment); // 4th comment line

  // read the main body of the file (table of values: energy, absorption length and refractive index)
  Int_t iLine=0;
  std::string sLine;
  getline(infile, sLine);
  while(!infile.eof()){
    if(iLine >= kNbins){
            AliFatal(Form("Line number: %i reaches range of declared arraySize: %i. Check input file: %s", iLine, kNbins, filePath.c_str()));
      return -5;
    }
    std::stringstream ssLine(sLine);
    ssLine >> (*de)[iLine];
    (*de)[iLine] *= 1e-9; // Convert eV -> GeV immediately
    (*e)[iLine] = static_cast<Float_t> ((*de)[iLine]); // same value, different precision
    ssLine >> (*abs)[iLine];
    ssLine >> (*n)[iLine];
    ssLine >> (*qe)[iLine];
    if(!(ssLine.good() || ssLine.eof())){ // check if there were problems with numbers conversion
            AliFatal(Form("Error while reading line %i: %s", iLine, ssLine.str().c_str()));
      return -6;
    }
    getline(infile, sLine);
    iLine++;
  }
  if(iLine != kNbins){
          AliFatal(Form("Total number of lines %i is different than declared %i. Check input file: %s", iLine, kNbins, filePath.c_str()));
    return -7;
  }
 
  AliInfo(Form("Optical properties taken from the file: %s. Number of lines read: %i",filePath.c_str(),iLine));
  return 0;
}




void AliFITv7::FillOtherOptProperties(Float_t **efficAll, Float_t **rindexAir, Float_t **absorAir,
    Float_t **rindexCathodeNext, Float_t **absorbCathodeNext,
    Double_t **efficMet, Double_t **aReflMet, const Int_t kNbins) const{
  // Allocate memory for these arrays according to the required size
  *efficAll = new Float_t[kNbins];
  *rindexAir = new Float_t[kNbins];
  *absorAir = new Float_t[kNbins];
  *rindexCathodeNext = new Float_t[kNbins];
  *absorbCathodeNext = new Float_t[kNbins];
  *efficMet = new Double_t[kNbins];
  *aReflMet = new Double_t[kNbins];

  // Set constant values to the arrays
  for(Int_t i=0; i<kNbins; i++)
  {
    (*efficAll)[i]=1.;
    (*rindexAir)[i] = 1.;
    (*absorAir)[i]=0.3;      
    (*rindexCathodeNext)[i]=0;
    (*absorbCathodeNext)[i]=0;
    (*efficMet)[i]=0.;
    (*aReflMet)[i]=1.;
  }
}

void AliFITv7::DeleteOptPropertiesArr(Float_t **e, Double_t **de, Float_t **abs,
    Float_t **n, Float_t **efficAll, Float_t **rindexAir, Float_t **absorAir,
    Float_t **rindexCathodeNext, Float_t **absorbCathodeNext,
    Double_t **efficMet, Double_t **aReflMet) const{
  delete [] (*e);
  delete [] (*de);
  delete [] (*abs);
  delete [] (*n);
  delete [] (*efficAll);
  delete [] (*rindexAir);
  delete [] (*absorAir);
  delete [] (*rindexCathodeNext);
  delete [] (*absorbCathodeNext);
  delete [] (*efficMet);
  delete [] (*aReflMet);
}


