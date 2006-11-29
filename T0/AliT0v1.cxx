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

#include <TGeometry.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TNode.h>
#include <TTUBE.h>
#include <TVirtualMC.h>
#include "TGeoManager.h"
#include "TGeoPgon.h"
#include "TGeoTube.h"
#include "TGeoCompositeShape.h"
#include <TGraph.h>

#include "AliLog.h"
#include "AliMagF.h"
#include "AliRun.h"

#include "AliT0hit.h"
#include "AliT0v1.h"
#include "AliT0Parameters.h"

#include "AliMC.h"
#include "AliCDBLocal.h"
#include "AliCDBStorage.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"

ClassImp(AliT0v1)

//--------------------------------------------------------------------
AliT0v1::AliT0v1(const char *name, const char *title):
  AliT0(name,title),
  fIdSens1(0)
{
  //
  // Standart constructor for T0 Detector version 0
  //
  fIshunt = 2;
  AliT0Parameters* param = AliT0Parameters::Instance();
  for (Int_t i=0; i<24; i++){
    TGraph* grEff = param ->GetPMTeff(i);
    fEffPMT.AddAtAndExpand(grEff,i);
  }
}
//_____________________________________________________________________________

AliT0v1::~AliT0v1() 
{
  // desctructor  
}

//-------------------------------------------------------------------------
void AliT0v1::CreateGeometry()
{
  //
  // Create the geometry of T0 Detector version 0
  //
  // begin Html
  /*
   <img src="gif/AliT0v0.gif">
  */
  //


   Int_t *idtmed = fIdtmed->GetArray();
   /*
   AliT0Parameters* param = AliT0Parameters::Instance();
   param->Init();
  Float_t zdetC = param->GetZposition(0);
  Float_t zdetA = param->GetZposition(1);
   */
   Float_t zdetC = 69.7;
   Float_t zdetA = 373.;

  Int_t is;
  Int_t idrotm[999];
  Float_t x,y,z;
  //C T0 mother volume
  Float_t pstartR[18]={0., 360., 5., 
		       -6.8, 4.25, 10., //-76.5+0.00+69.7
		       -0.75 , 4.5, 10.,  // -76.5+6.05+69.7
		       1.25   , 4.5, 10.,  //-76.5+8.05+69.7
		       1.25 , 5.1, 10., //-76.5+8.05+69.7 
		       6.8 , 5.1, 10.};  //-62.9+0.00+69.7
  
  Float_t pstart[3]={4., 12.5,6.8};
  Float_t pinstart[3]={0.,1.6,6.5};
  Float_t ppmt[3]={0.,1.5,3.5};
  Float_t ptop[3]={0.,1.,1.0};
  Float_t preg[3]={0., 1.0, 0.005}; //photcathode dobavil bogdanov
  Float_t ptopout[3]={1.,1.01, 1.};  //metal air around raiator for photos reflection 
  Float_t pbot[3]={0.6,1.2,0.1};
  Float_t pglass[3]={1.2,1.3,2.};
  Float_t pcer[3]={0.9,1.1,1.35};
  Float_t psteel[3]={0.9,1.1,0.15};
  Float_t psupport1[3] = {4.51,4.6,4.0};//C kozhuh vnutri
  Float_t psupport2[3] = {9.4,9.5,4.0};// snaruzhi  C
  Float_t psupport3[3] = {4.51,9.5,0.05};//kryshki  C
  Float_t psupport5[3] = {1.44,1.5,6.5}; // stakanchik dlai feu  C
  Float_t psupport6[3] = {0,1.4,0.04}; //kryshechka stakanchika  Al
  Float_t psupport7[3] = {1.5,1.6,0.4}; //kolechko snaruzhu stakanchika Al
   
  // Mother Volume katushka dlia krepezha vokrug truby k Absorbru
    AliMatrix(idrotm[901], 90., 0., 90., 90., 180., 0.);
  Float_t ppcon[70]; 
    ppcon[0]  =   0;
    ppcon[1]  = 360;
    ppcon[2]  =  13;
//  1: 
    ppcon[3]  =  7.2;
    ppcon[4]  =   4.4;
    ppcon[5]  =   4.5;
//  2
    ppcon[6]  = 8.2;
    ppcon[7]  = 4.4;
    ppcon[8]  = 4.5;
//  3
    ppcon[9]  = 8.2;
    ppcon[10] = 4.4;
    ppcon[11] = 5.1;

//  4

    ppcon[12] = 8.3; 
    ppcon[13] = 4.4;
    ppcon[14] = 5.1;
//  5

    ppcon[15] = 8.3;
    ppcon[16] = 4.9;
    ppcon[17] = 5.1;
    
//  6
    ppcon[18] = 15.; 
    ppcon[19] = 4.9;
    ppcon[20] = 5.1;
    
//  7
    ppcon[21] = 15.; 
    ppcon[22] = 3.15  ;
    ppcon[23] = 5.1;

/// 8
    ppcon[24] = 15.01; 
    ppcon[25] = 3.15;
    ppcon[26] = 5.1;

/// 9
    ppcon[27] = 15.01; 
    ppcon[28] = 3.15;
    ppcon[29] = 3.25;

    
/// 9
     ppcon[30] = 19.71;
     ppcon[31] = 3.15;
     ppcon[32] = 3.25; 
//  10
    ppcon[33] = 19.8; 
    ppcon[34] = 3.15;
    ppcon[35] = 3.25;
    //  11
    ppcon[36] = 19.8;
    ppcon[37] = 3.15;
    ppcon[38] = 7.6;
        
   //  14

    ppcon[39] =  19.91;
    ppcon[40] =  3.15;
    ppcon[41] =  7.6;


    gMC->Gsvolu("0SUP", "PCON", idtmed[kAir], ppcon,42);
    z = -zdetC;
    gMC->Gspos("0SUP",1,"ALIC",0.,0.,z,idrotm[901],"ONLY");
   
//-------------------------------------------------------------------
 //  T0 volume 
 //-------------------------------------------------------------------
  
    
    gMC->Gsvolu("0STR","PCON",idtmed[kAir],pstartR,18);
    gMC->Gsvolu("0STL","TUBE",idtmed[kAir],pstart,3);
    gMC->Gspos("0STR",1,"ALIC",0.,0.,-zdetC+pstartR[3],idrotm[901],"ONLY");
    // gMC->Gspos("0STL",1,"ALIC",0.,0.,zdetA+pstart[2],0,"ONLY");

//T0 interior
   gMC->Gsvolu("0INS","TUBE",idtmed[kAir],pinstart,3);
   gMC->Gsvolu("0PMT","TUBE",idtmed[kAir],ppmt,3);     
          
    gMC->Gsvolu("0SU1","TUBE",idtmed[kC],psupport1,3);//C kozhuh vnutri
    gMC->Gsvolu("0SU2","TUBE",idtmed[kC],psupport2,3);// snaruzhi  C
    gMC->Gsvolu("0SU3","TUBE",idtmed[kC],psupport3,3);//kryshka perednaiai  C
    gMC->Gsvolu("0SU4","TUBE",idtmed[kC],psupport3,3);//kryshka zadnaiai  C
    //    gMC->Gsvolu("0SU5","TUBE",idtmed[kAir],psupport4,3);// dyrki dlia feu v zadnej kryshke Air
    gMC->Gsvolu("0SU6","TUBE",idtmed[kC],psupport5,3);// stakanchik dlai feu  C
    gMC->Gsvolu("0SU7","TUBE",idtmed[kAl],psupport6,3);//kryshechka stakanchika  Al
    gMC->Gsvolu("0SU8","TUBE",idtmed[kAl],psupport7,3);//kolechko snaruzhu stakanchika Al
   
 

   //non-absorber side support  and T0A !!!!!!!!
    
    TGeoPgon * supPgon = new TGeoPgon("supPgon",0.,360.,360,4);
    supPgon->DefineSection(0, 0, 4.1, 5.5);
    supPgon->DefineSection(1, 10.5 , 4.1, 5.5);
    supPgon->DefineSection(2, 10.5 , 4.1, 4.9);
    supPgon->DefineSection(3, 12.5 , 4.1, 4.9);
    TGeoTranslation *trPgon = new TGeoTranslation("trPgon",0,0,0);
    trPgon->RegisterYourself();
 
    TGeoVolumeAssembly * stlin = new TGeoVolumeAssembly("OSTL");//empty segment
    TGeoVolume *ins = gGeoManager->GetVolume("0INS");

    Double_t phimin = TMath::ACos((16-4.8)/16.) * (180 / TMath::Pi()) ;
    //  TGeoTubeSeg *hole = new TGeoTubeSeg("hole", 0, 1.6, 0, -phimin+90, phimin+90);
    TGeoTube*hole = new TGeoTube("hole", 0, 1.61, 6.5);
    TGeoTranslation *tr [12];
    Float_t angle  = 2 * TMath::Pi() / 12;
    Char_t nameTr[4];
    for (Int_t itr=0; itr<12; itr++) {
      sprintf (nameTr,"tr%i",itr+1);
      x = 6.5 * TMath::Sin(itr * angle);
      y = 6.5 * TMath::Cos(itr * angle);
      tr[itr] = new TGeoTranslation(nameTr,x,y,6.5);
      tr[itr]->RegisterYourself();
       stlin->AddNode(ins,itr+13,tr[itr]);
    }
    TGeoCompositeShape *supsh = new TGeoCompositeShape("supsh","supPgon:trPgon-(hole:tr1+hole:tr2+hole:tr3+hole:tr4+hole:tr5+hole:tr6+hole:tr7+hole:tr8+hole:tr9+hole:tr10+hole:tr11+hole:tr12)");
  
    TGeoMedium *medal = gGeoManager->GetMedium("T0_Aluminium$");
    TGeoVolume *supA = new TGeoVolume("0SUA",supsh,medal);
    stlin->AddNode(supA,1,new TGeoTranslation(0,0,0));
 
    TGeoVolume *alice = gGeoManager->GetVolume("ALIC");
    alice->AddNode(stlin,1,new TGeoTranslation(0,0, zdetA ));
            
    // T0A finished

    //T0C 
// first ring: 12 units of Scintillator+PMT+divider
  Float_t  theta  = (180 / TMath::Pi()) * TMath::ATan(6.5 / zdetC);
  Float_t phi[3];
   for (is=0; is<12; is++)
      {  

	x = 6.5 * TMath::Sin(is * angle);
	y = 6.5 * TMath::Cos(is * angle);
	
	phi[0] = -30 * is;
	phi[1] = 90 - is * 30;
	phi[2] = 90 - is * 30;
	for (Int_t j = 0; j < 3; j++)
	  if (phi[j] < 0)  phi[j] += 360;
	
	AliMatrix (idrotm[902 + is], 90.,         phi[0],
		   90. + theta, phi[1],
		   theta,       phi[2]);  
	z=-pstart[2]+pinstart[2]+0.2;
	gMC->Gspos ("0INS", is + 1, "0STR", x, y, z, idrotm[902 + is], "ONLY");
  }	
     
   x=0;
   y=0;
    z=-pinstart[2]+ppmt[2]+psupport6[2]*2;
   gMC->Gspos("0PMT",1,"0INS",x,y,z,0,"ONLY");
   // PMT
   
   // Entry window (glass)
   gMC->Gsvolu("0TOP","TUBE",idtmed[kOpGlass],ptop,3); //glass
   z=-ppmt[2]+ptop[2];
   gMC->Gspos("0TOP",1,"0PMT",0,0,z,0,"ONLY");
   //"metal" air volume to simulate reclection  
   gMC->Gsvolu("0TOO","TUBE",idtmed[kOpAir],ptopout,3); 
   gMC->Gspos("0TOO",1,"0PMT",0,0,z,0,"ONLY");

   //Fotokatod
   gMC->Gsvolu ("0REG", "TUBE", idtmed[kOpGlassCathode], preg, 3); 
   z = -ppmt[2] + 2 * ptop[2] + preg[2]; 
   gMC->Gspos ("0REG", 1, "0PMT", 0, 0, z, 0, "ONLY"); 
  // Bottom glass
   gMC->Gsvolu("0BOT","TUBE",idtmed[kGlass],pbot,3);
   z=ppmt[2]-pbot[2];
   gMC->Gspos("0BOT",1,"0PMT",0,0,z,0,"ONLY");
   // Side cylinder glass
   gMC->Gsvolu("0OUT","TUBE",idtmed[kGlass],pglass,3);
   z=ppmt[2]-pglass[2];
   gMC->Gspos("0OUT",1,"0PMT",0,0,z,0,"ONLY");
   //PMT electrodes support structure
   gMC->Gsvolu("0CER","TUBE",idtmed[kCer],pcer,3);
   gMC->Gsvolu("0STE","TUBE",idtmed[kSteel],psteel,3);
   z=-ppmt[2]+2*ptop[2]+0.3 + pcer[2];
   gMC->Gspos("0CER",1,"0PMT",0,0,z,0,"ONLY");
   z +=psteel[2]+pcer[2];
   gMC->Gspos("0STE",1,"0PMT",0,0,z,0,"ONLY");
    
   //Support absorber (C) side
   z=-pstart[2]+psupport1[2]+0.1;
   gMC->Gspos("0SU1",1,"0STR",0,0,z,0,"ONLY"); //C kozhuh snaruzhi
   gMC->Gspos("0SU2",1,"0STR",0,0,z,0,"ONLY"); //C kozhuh vnutri
   z=-pstart[2]+psupport3[2]+0.1;
   gMC->Gspos("0SU3",1,"0STR",0,0,z,0,"ONLY"); //peredniaia kryshka
   z=-pstart[2]+2.*psupport1[2]-0.05;//+0.1;
   gMC->Gspos("0SU4",1,"0STR",0,0,z,0,"MANY"); //zadnaiai kryshka
   gMC->Gspos("0SU6",1,"0INS",0,0,0,0,"ONLY");//C stakanchik dlia feu 
   z=-pinstart[2]+psupport6[2];
   gMC->Gspos("0SU7",1,"0INS",0,0,z,0,"ONLY"); //Al kryshechka 
   
   z=pinstart[2]-psupport7[2];
   gMC->Gspos("0SU8",1,"0INS",0,0,z,0,"ONLY"); //Al kolechko
   
   
   Float_t par[3];
   par[0]=4.4;
   par[1]=4.5;
   par[2]=0.5;
   gMC->Gsvolu("0SC0","TUBE",idtmed[kC],par,3);
   z=ppcon[3]+par[2];
   gMC->Gspos("0SC0",1,"0SUP",0,0,z,0,"ONLY"); 
   z += par[2];
   par[0]=4.4;
   par[1]=5.1;
   par[2]=0.05;
   gMC->Gsvolu("0SC1","TUBE",idtmed[kC],par,3);
   z += par[2];
   gMC->Gspos("0SC1",1,"0SUP",0,0,z,0,"ONLY"); 
   z=z+par[2];
   par[0]=4.9;
   par[1]=5.0;
   par[2]=6.7/2;
    gMC->Gsvolu("0SC2","TUBE",idtmed[kC],par,3);
    z += par[2];
    gMC->Gspos("0SC2",1,"0SUP",0,0,z,0,"ONLY"); 
    z += par[2];
   
    par[0]=3.15;
    par[1]=4.9;
    par[2]=0.1/2;
    gMC->Gsvolu("0SA1","TUBE",idtmed[kAl],par,3);
    
    z += par[2];
    cout<<" z 0SA1 "<<z<<endl;
    gMC->Gspos("0SA1",1,"0SUP",0,0,z,0,"ONLY"); 
    z=z+par[2];
    par[0]=3.15;
    par[1]=3.16;
    par[2]=4.7/2;
    cout<<" z 0SA2 "<<z<<endl;
    gMC->Gsvolu("0SA2","TUBE",idtmed[kAl],par,3);
    z += par[2];
    gMC->Gspos("0SA2",1,"0SUP",0,0,z,0,"ONLY"); 
    z=z+par[2];
    par[0]=3.16; // eta chast' prikruchena k absorberu
    par[1]=7.5;
    par[2]=0.1;
    cout<<" z 0SA23"<<z<<endl;
    gMC->Gsvolu("0SA3","TUBE",idtmed[kAl],par,3);
    z += par[2];
    gMC->Gspos("0SA3",1,"0SUP",0,0,z,0,"ONLY"); 
    par[0]=3.16; // gvozdi eta chast' prikruchena k absorberu
    par[1]=7.5;
    par[2]=0.01;
    cout<<" z 0SN2 "<<z<<endl;
    gMC->Gsvolu("0SN2","TUBE",idtmed[kSteel],par,3);
    gMC->Gspos("0SN2",1,"0SUP",0,0,z,0,"ONLY"); 
 
   

 
}    
//------------------------------------------------------------------------
void AliT0v1::AddAlignableVolumes() const
{
  //
  // Create entries for alignable volumes associating the symbolic volume
  // name with the corresponding volume path. Needs to be syncronized with
  // eventual changes in the geometry.
  //

  TString volPath;
  TString symName, sn;

  TString vpA  = "/ALIC_1/0STL_1/0INS_";
  TString vpC  = "/ALIC_1/0STR_1/0INS_";
  TString vpInside  = "/0PMT_1/0TOP_1";


  for (Int_t imod=0; imod<24; imod++)
    {
      if (imod < 12) 
	{volPath  = vpC; sn="T0/C/PMT";}
      else  
	{volPath  = vpA; sn="T0/A/PMT";}
      volPath += imod+1;
      volPath += vpInside;
      
      symName = sn;
      symName += imod+1;
      
      AliDebug(2,"--------------------------------------------");
      AliDebug(2,Form("Alignable object %d", imod));
      AliDebug(2,Form("volPath=%s\n",volPath.Data()));
      AliDebug(2,Form("symName=%s\n",symName.Data()));
      AliDebug(2,"--------------------------------------------");
      gGeoManager->SetAlignableEntry(symName.Data(),volPath.Data());
      // if(!gGeoManager->SetAlignableEntry(symName.Data(),volPath.Data()))
      //	AliFatal(Form("Alignable entry %s not created. Volume path %s not valid", 
      //	      symName.Data(),volPath.Data()));
      
    }
}   
//------------------------------------------------------------------------
void AliT0v1::CreateMaterials()
{
   Int_t isxfld   = gAlice->Field()->Integ();
   Float_t sxmgmx = gAlice->Field()->Max();
   //   Float_t a,z,d,radl,absl,buf[1];
   // Int_t nbuf;
   Int_t *idtmed = fIdtmed->GetArray();

// AIR
                                                                                
   Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
   Float_t zAir[4]={6.,7.,8.,18.};
   Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
   Float_t dAir = 1.20479E-3;
                                                                                
// Scintillator CH
   Float_t ascin[2]={1.01,12.01};
   Float_t zscin[2]={1,6};
   Float_t wscin[2]={1,1};
   Float_t denscin=1.03;
//Lucite C(CH3)CO2CH3
   Float_t alucite[3]={1.01,12.01,15.999};
   Float_t zlucite[3]={1,6,8};
   Float_t wlucite[3]={8,5,2};
   Float_t denlucite=1.16;
 // PMT glass SiO2
   Float_t aglass[2]={28.0855,15.9994};
   Float_t zglass[2]={14.,8.};
   Float_t wglass[2]={1.,2.};
   Float_t dglass=2.65;
// Ceramic   97.2% Al2O3 , 2.8% SiO2
   // Float_t wcer[2]={0.972,0.028};
   Float_t aCeramic[2]  = { 26.981539,15.9994 };
   Float_t zCeramic[2]  = { 13.,8. };
   Float_t wCeramic[2]  = { 2.,3. };
   Float_t denscer  = 3.6;

// Brass 80% Cu, 20% Zn
   Float_t abrass[2] = {63.546,65.39};
   Float_t zbrass[2] = {29,30};
   Float_t wbrass[2] = {0.8,0.2};
   Float_t denbrass=8.96;

//Ribber C6H12S
   Float_t aribber[3] = {12.,1.,32.};
   Float_t zribber[3] = {6.,1.,16.};
   Float_t wribber[3] = {6.,12.,1.};
   Float_t denribber=0.8;
// Support inside 
   Float_t asupport[2] = {12.,1.};
   Float_t zsupport[2] = {6.,1.};
   Float_t wsupport[2] = {1.,1.};
   Float_t densupport=0.1;
    
//*** Definition Of avaible T0 materials ***
   AliMaterial(0, "Steel$", 55.850,26.,7.87,1.76,999);
   AliMaterial(1, "Vacuum$", 1.e-16,1.e-16,1.e-16,1.e16,999);
   AliMixture(2, "Air$", aAir, zAir, dAir,4,wAir);
   AliMaterial(10, "CarbonPlastic$", 12.01, 6.0, 2.26, 18.8,999); 
   AliMaterial(11, "Aliminium$", 26.98, 13.0, 2.7, 8.9,999); 

   AliMixture( 3, "Ceramic  $",aCeramic, zCeramic, denscer, -2, wCeramic);
   AliMixture( 4, "PMT glass   $",aglass,zglass,dglass,-2,wglass);
   /* 
  char namate[21]="";
  
   gMC->Gfmate((*fIdmate)[3], namate, a, z, d, radl, absl, buf, nbuf);
   acer[0]=a;  

   zcer[0]=z;
   gMC->Gfmate((*fIdmate)[4], namate, a, z, d, radl, absl, buf, nbuf);
   acer[1]=a;
   zcer[1]=z;
   AliMixture( 11, "Ceramic    $", acer, zcer, denscer, 2, wcer);
   */ 
    AliMixture( 5, "Scintillator$",ascin,zscin,denscin,-2,wscin);
   AliMixture( 6, "Brass    $", abrass, zbrass, denbrass, 2, wbrass);
   AliMixture( 7, "Ribber $",aribber,zribber,denribber,-3,wribber);
   AliMixture( 8, "Lucite$",alucite,zlucite,denlucite,-3,wlucite);
   AliMixture( 9, "Penoplast$",asupport,zsupport,densupport,-2,wsupport);
   AliMixture( 21, "PMT Optical glass   $",aglass,zglass,dglass,-2,wglass);
   AliMixture( 24, "PMT Optical glass cathode $",aglass,zglass,dglass,-2,wglass);
    AliMixture(22, "T0 Opt Air$", aAir, zAir, dAir,4,wAir);
   AliMixture(23, "T0 Opt Air Next$", aAir, zAir, dAir,4,wAir);
   
   AliMedium(1, "T0 Air$", 2, 0, isxfld, sxmgmx, 10., .1, 1., .003, .003);
   AliMedium(2, "Scintillator$", 5, 1, isxfld, sxmgmx, 10., .01, 1., .003, .003);
   AliMedium(3, "Vacuum$", 1, 0, isxfld, sxmgmx, 10., .01, .1, .003, .003);
   AliMedium(4, "Ceramic$", 3, 0, isxfld, sxmgmx, 10., .01, .1, .003, .003);
   AliMedium(6, "Glass$", 4, 1, isxfld, sxmgmx, 10., .01, .1, .003, .003);
   AliMedium(8, "Steel$", 0, 0, isxfld, sxmgmx, 1., .001, 1., .001, .001);
   AliMedium(9, "Ribber  $", 7, 0, isxfld, sxmgmx, 10., .01, .1, .003, .003);
   AliMedium(11, "Brass  $", 6, 0, isxfld, sxmgmx, 10., .01, .1, .003, .003);
   AliMedium(12, "Lucite$", 8, 1, isxfld, sxmgmx, 10., .01, 1., .003, .003);  
   AliMedium(13, "CarbonPlastic$", 10, 0, isxfld, sxmgmx, 10., .01, 1., .003, .003);  
   AliMedium(14, "PenoPlast$", 9, 0, isxfld, sxmgmx, 10., .01, 1., .003, .003);  
   AliMedium(15, "Aluminium$", 11, 0, isxfld, sxmgmx, 10., .01, 1., .003, .003);  
   AliMedium(16, "OpticalGlass$", 21, 1, isxfld, sxmgmx, 10., .01, .1, .003, .003);
   AliMedium(17, "T0 OpAir$", 22, 0, isxfld, sxmgmx, 10., .1, 1., .003, .003);
   AliMedium(18, "T0 OpAirNext$", 23, 0, isxfld, sxmgmx, 10., .1, 1., .003, .003);
    AliMedium(19, "OpticalGlassCathode$", 24, 1, isxfld, sxmgmx, 10., .01, .1, .003, .003);


// Definition Cherenkov parameters
   int i;
   const Int_t kNbins=31;
   
   Float_t rindexSiO2[kNbins],  efficAll[kNbins], rindexAir[kNbins], absorAir[kNbins],rindexCathodeNext[kNbins], absorbCathodeNext[kNbins];
            
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
      rindexAir[i]=0;
      rindexSiO2[i]=1.458; //refractive index for qwarts
      rindexCathodeNext[i]=0;
      efficAll[i]=1.;

      //      aAbsSiO2[i]=28.5; //quartz 30mm
      absorAir[i]=0.3;      
      absorbCathodeNext[i]=0;

    }
  
  
   gMC->SetCerenkov (idtmed[kOpGlass], kNbins, aPckov, aAbsSiO2, efficAll, rindexSiO2 );
   // gMC->SetCerenkov (idtmed[kOpGlassCathode], kNbins, aPckov, aAbsSiO2, effCathode, rindexSiO2 );
   gMC->SetCerenkov (idtmed[kOpGlassCathode], kNbins, aPckov, aAbsSiO2,efficAll , rindexSiO2 );
  gMC->SetCerenkov (idtmed[kOpAir], kNbins, aPckov,absorAir , efficAll,rindexAir );
   gMC->SetCerenkov (idtmed[kOpAirNext], kNbins, aPckov,absorbCathodeNext , efficAll, rindexCathodeNext);
   
   AliDebugClass(1,": ++++++++++++++Medium set++++++++++");
   
   
}
//---------------------------------------------------------------------
void AliT0v1::DrawDetector()
{
//
// Draw a shaded view of the Forward multiplicity detector version 0
//
  
  //Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN",0);
  //
  //Set volumes visible
  //  gMC->Gsatt("0STR","SEEN",0);
  //  gMC->Gsatt("0INS","SEEN",0);
  // gMC->Gsatt("0PMT","SEEN",1);
  // gMC->Gsatt("0DIV","SEEN",1);
  //
  gMC->Gdopt("hide","off");
  gMC->Gdopt("shad","on");
  gMC->SetClipBox(".");
  gMC->SetClipBox("*",0,1000,-1000,1000,-1000,1000);
  gMC->DefaultRange();
  gMC->Gdraw("alic",90,0,0,-35,9.5,.6,0.6);
  //gMC->Gdraw("alic",0,0,0,10,9.5,.8,0.8); //other side view
  gMC->Gdhead(1111,"T-Zero detector");
  gMC->Gdopt("hide","off");
}

//-------------------------------------------------------------------
void AliT0v1::Init()
{
// Initialises version 0 of the Forward Multiplicity Detector
//
//Int_t *idtmed  = gAlice->Idtmed();
  AliT0::Init();
  fIdSens1=gMC->VolId("0REG");

   AliDebug(1,Form("%s: *** T0 version 1 initialized ***\n",ClassName()));

}

//-------------------------------------------------------------------

void AliT0v1::StepManager()
{
  //
  // Called for every step in the T0 Detector
  //
  Int_t id,copy,copy1;
  static Float_t hits[6];
  static Int_t vol[2];
  TLorentzVector pos;
  
  //   TClonesArray &lhits = *fHits;
  
  if(!gMC->IsTrackAlive()) return; // particle has disappeared
  
  
    // If particles is photon then ...
  
 if (gMC->TrackPid() == 50000050)
  {
      id=gMC->CurrentVolID(copy);
      
      
      // Check the sensetive volume
      if(id==fIdSens1 ) { 
	if(gMC->IsTrackEntering()) {
	    gMC->CurrentVolOffID(2,copy);
	    vol[1]=copy;
	    gMC->CurrentVolOffID(3,copy1);
	    vol[0]=copy1;
	    gMC->TrackPosition(pos);
	    hits[0] = pos[0];
	    hits[1] = pos[1];
	    hits[2] = pos[2];
	    if(pos[2]<0) vol[0]=2;
	    if(pos[2]>=0) vol[0]=1; 
	    
	    Float_t etot=gMC->Etot();
	    hits[3]=etot;
	    Int_t iPart= gMC->TrackPid();
	    Int_t partID=gMC->IdFromPDG(iPart);
	    hits[4]=partID;
	    Float_t ttime=gMC->TrackTime();
	    hits[5]=ttime*1e12;
	    if(RegisterPhotoE(vol[1]-1,hits[3])) 
	      AddHit(fIshunt,vol,hits);
	    
	  }
	  //  cout<< gAlice->GetMCApp()->GetCurrentTrackNumber()<<" hit added "<<endl;
	
	/*		
		printf("track(%i) alive(%i) disap(%i) enter(%i) exit(%i) inside(%i) out(%i) stop(%i) new(%i) \n",
		       gAlice->GetMCApp()->GetCurrentTrackNumber(),
	       gMC->IsTrackAlive(),
	       gMC->IsTrackDisappeared(),
	       gMC->IsTrackEntering(),
	       gMC->IsTrackExiting(),
	       gMC->IsTrackInside(),
	       gMC->IsTrackOut(),
	       gMC->IsTrackStop(),
	       gMC->IsNewTrack());
	*/		
      } //sensitive
  } //photon
}



//------------------------------------------------------------------------
Bool_t AliT0v1::RegisterPhotoE(Int_t ipmt,Double_t energy)
{

  
  //  Float_t hc=197.326960*1.e6; //mev*nm
  Double_t hc=1.973*1.e-6; //gev*nm
  Float_t lambda=hc/energy;
  Float_t eff = ((TGraph*) fEffPMT.At(ipmt))->Eval(lambda);
  Double_t  p = gRandom->Rndm();

  if (p > eff)
    return kFALSE;
  
  return kTRUE;
}

//----------------------------------------------------------------------------






