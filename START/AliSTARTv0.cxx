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
// START ( T-zero) detector  version 0                        //
//
//Begin Html       
/*
<img src="gif/AliSTARTv0Class.gif">
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

#include "AliMagF.h"
#include "AliRun.h"
#include "AliSTARThit.h"
#include "AliSTARTv0.h"

ClassImp(AliSTARTv0)

//--------------------------------------------------------------------
AliSTARTv0::AliSTARTv0(const char *name, const char *title):
 AliSTART(name,title)
{
  //
  // Standart constructor for START Detector version 0
  //
  fIdSens1=0;
//  setBufferSize(128000);
}
//-------------------------------------------------------------------------
void AliSTARTv0::CreateGeometry()
{
  //
  // Create the geometry of START Detector version 0
  // Full geometry with support structure according START prototype
  //
  // begin Html
  /*
   <img src="gif/AliSTARTv0.gif">
  */
  //



  Int_t *idtmed = fIdtmed->GetArray();
  
  Int_t is;
  Int_t idrotm[999];
  Float_t x,y,z;

  Float_t pstart[3]={4.5,10.,4.0};
  Float_t pinstart[3]={0.,1.6,6.5};
  Float_t ppmt[3]={0.,1.3,3.5};
  Float_t pdivider[3]={0.,1.2,1.75};
  Float_t pdiv2[3]={0.,1.2,1.25};
  Float_t pdiv1[3]={0.6,1.2,0.5};
  Float_t ptop[3]={0.,1.3,1.5};
  Float_t pbot[3]={0.6,1.2,0.1};
  Float_t pglass[3]={1.2,1.3,2.};
  Float_t pcer[3]={0.9,1.1,0.09};
  Float_t psteel[3]={0.9,1.1,0.01};
  Float_t ppins[3]={0.6,1.2,0.014};
  Float_t phole[3]={0.6,1.2,0.015};
  Float_t pknob[3]={0.5,0.6,0.4};
  Float_t pknob_vac[3]={0.,0.5,0.4};
  Float_t pknob_bot[3]={0.,0.6,0.05};
  Float_t pribber[3] = {0.,1.2,2.413/2.};
  Float_t presist[3] = {0.,1.2,0.087/2.};
  Float_t psupport1[3] = {4.5,4.6,4.0};//C kozhuh vnutri
  Float_t psupport2[3] = {9.4,9.5,4.0};// snaruzhi  C
  Float_t psupport3[3] = {4.5,9.5,0.05};//kryshki  C
  Float_t psupport4[3] = {0,1.6,0.05};// dyrki dlia feu v zadnej kryshke Air
  Float_t psupport5[3] = {1.44,1.5,6.5}; // stakanchik dlai feu  C
  Float_t psupport6[3] = {0,1.5,0.05}; //kryshechka stakanchika  Al
  Float_t psupport7[3] = {1.5,1.6,0.6}; //kolechko snaruzhu stakanchika Al
  // Mother Volume katushka dlia krepezha vokrug truby k Absorbru
  AliMatrix(idrotm[901], 90., 0., 90., 90., 180., 0.);
  Float_t ppcon[70]; 
    ppcon[0]  =   0;
    ppcon[1]  = 360;
    ppcon[2]  =  15;
//  1: 
    ppcon[3]  =  14.1/2;
    ppcon[4]  =   4.4;
    ppcon[5]  =   4.5;
//  2
    ppcon[6]  = ppcon[3]+1.;
    ppcon[7]  =   4.4;
    ppcon[8]  =   4.5;
//  3
    ppcon[9]  = ppcon[6];
    ppcon[10] =   4.4;
    ppcon[11] =   5.3;

//  4

    ppcon[12]  = ppcon[9]+0.1;
    ppcon[13] =   4.4;
    ppcon[14] =   5.3;
//  5

    ppcon[15]  = ppcon[12];
    ppcon[16] =   5.1;
    ppcon[17] =   5.3;
    
//  6
    ppcon[18]  = ppcon[12]+5.8;
    ppcon[19] =   5.1;
    ppcon[20] =   5.3;
    
//  7
    ppcon[21]  = ppcon[18];
    ppcon[22] =   5.1;
    ppcon[23] =   5.9;
    
//  8
    ppcon[24]  = ppcon[21]+1.5;
    ppcon[25] =   5.1;
    ppcon[26] =   5.9;
//  9
    ppcon[27]  = ppcon[24]+0.5;
    ppcon[28] =   5.1;
    ppcon[29] =   5.9;

//  10
    ppcon[30]  = ppcon[24]+0.5;
    ppcon[31] =   3.15;
    ppcon[32] =   5.9;

//  11
    ppcon[33]  = ppcon[30];
    ppcon[34] =   3.15;
    ppcon[35] =   5.9;
//  12
    ppcon[36]  = ppcon[30];
    ppcon[37] =   3.15;
    ppcon[38] =   3.25;
    
//  13
    ppcon[39]  = ppcon[36]+4.5;
    ppcon[40] =   3.15;
    ppcon[41] =   3.25;
//  14
    ppcon[42]  = ppcon[36]+4.5;
    ppcon[43] =   3.15;
    ppcon[44] =   7.6;
//  15
    ppcon[45]  = ppcon[42]+0.9;
    ppcon[46] =   3.15;
    ppcon[47] =   7.6;
    
    gMC->Gsvolu("0SUP", "PCON", idtmed[kAir], ppcon,48);
    z=69.7; //+14.1/2;
    gMC->Gspos("0SUP",1, "ALIC", 0,0,z,idrotm[901],"ONLY");

  Float_t zdetRight=69.7,zdetLeft=350;
 //-------------------------------------------------------------------
 //  START volume 
 //-------------------------------------------------------------------
  
    
    /*<<<<<<< AliSTARTv0.cxx
    gMC->Gsvolu("0STR","TUBE",idtmed[kAir],pstart,3);
    gMC->Gsvolu("0STL","TUBE",idtmed[kAir],pstart,3);
    gMC->Gspos("0STR",1,"ALIC",0.,0.,zdetRight+pstart[2],0,"ONLY");
    gMC->Gspos("0STL",1,"ALIC",0.,0.,-zdetLeft-pstart[2],idrotm[901],"ONLY");
    =======*/
    gMC->Gsvolu("0STR","TUBE",idtmed[1],pstart,3);
    gMC->Gsvolu("0STL","TUBE",idtmed[1],pstart,3);
    gMC->Gspos("0STR",1,"ALIC",0.,0.,-zdetRight-pstart[2],idrotm[901],"ONLY");
    gMC->Gspos("0STL",1,"ALIC",0.,0.,zdetLeft+pstart[2],0,"ONLY");
    //>>>>>>> 1.19

//START interior
    gMC->Gsvolu("0INS","TUBE",idtmed[kAir],pinstart,3);
    gMC->Gsvolu("0PMT","TUBE",idtmed[kVac],ppmt,3);     
    gMC->Gsvolu("0DIV","TUBE",idtmed[kVac],pdivider,3);     
    gMC->Gsvolu("0SU1","TUBE",idtmed[kC],psupport1,3);//C kozhuh vnutri
    gMC->Gsvolu("0SU2","TUBE",idtmed[kC],psupport2,3);// snaruzhi  C
    gMC->Gsvolu("0SU3","TUBE",idtmed[kC],psupport3,3);//kryshka perednaiai  C
    gMC->Gsvolu("0SU4","TUBE",idtmed[kC],psupport3,3);//kryshka zadnaiai  C
    gMC->Gsvolu("0SU5","TUBE",idtmed[kAir],psupport4,3);// dyrki dlia feu v zadnej kryshke Air
    cout<<" 0Su6 >> "<<psupport5[0]<<" "<<psupport5[1]<<" "<<psupport5[2]<<endl;
    gMC->Gsvolu("0SU6","TUBE",idtmed[kC],psupport5,3);// stakanchik dlai feu  C
    gMC->Gsvolu("0SU7","TUBE",idtmed[kAl],psupport6,3);//kryshechka stakanchika  Al
    gMC->Gsvolu("0SU8","TUBE",idtmed[kAl],psupport7,3);//kolechko snaruzhu stakanchika Al
       
// first ring: 12 units of Scintillator+PMT+divider
  Float_t  theta  = (180 / TMath::Pi()) * TMath::ATan(6.5 / zdetRight);
  Float_t angel  = 2 * TMath::Pi() / 12;
  Float_t phi[3];
   for (is=0; is<12; is++)
      {  

	x = 6.5 * TMath::Sin(is * angel);
	y = 6.5 * TMath::Cos(is * angel);
	
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
     gMC->Gspos ("0INS", is + 13, "0STL", x, y, z, 0, "ONLY");
      }	

   
      x=0;
    y=0;
    z=-pinstart[2]+ppmt[2]+2.*psupport6[2];
    gMC->Gspos("0PMT",1,"0INS",x,y,z,0,"ONLY");
    z=z+pdivider[2]+ppmt[2];
    gMC->Gspos("0DIV",1,"0INS",x,y,z,0,"ONLY");
      
// PMT
      
    // Entry window (glass)
    gMC->Gsvolu("0TOP","TUBE",idtmed[kGlass],ptop,3); //glass
     //   gMC->Gsvolu("0TOP","TUBE",idtmed[12],ptop,3); //lucite
     z=-ppmt[2]+ptop[2];
    gMC->Gspos("0TOP",1,"0PMT",0,0,z,0,"ONLY");
    //     printf("Z PTOP %f -ppmt[2] %f ptop[2] %f\n",z,-ppmt[2],ptop[2]);
    
    // Bottom glass
    gMC->Gsvolu("0BOT","TUBE",idtmed[kGlass],pbot,3);
    z=ppmt[2]-pbot[2];
    if(fDebug) printf("%s: Z bottom %f\n",ClassName(),z);
    gMC->Gspos("0BOT",1,"0PMT",0,0,z,0,"ONLY");
    // Side cylinder glass
    gMC->Gsvolu("0OUT","TUBE",idtmed[kGlass],pglass,3);
    z=ppmt[2]-pglass[2];
    //      printf("Z glass %f\n",z);
    gMC->Gspos("0OUT",1,"0PMT",0,0,z,0,"ONLY");
    //PMT electrodes support structure
    gMC->Gsvolu("0CER","TUBE",idtmed[kCer],pcer,3);
    gMC->Gsvolu("0STE","TUBE",idtmed[kSteel],psteel,3);
    z=-ppmt[2]+2*ptop[2]+0.3;;
    //      printf("Z Cer 1 %f\n",z);
    for (is=1; is<=15; is++)
      {
	z=z+psteel[2]+pcer[2];
	gMC->Gspos("0CER",is,"0PMT",0,0,z,0,"ONLY");
	z=z+psteel[2]+pcer[2];
	gMC->Gspos("0STE",is,"0PMT",0,0,z,0,"ONLY");
      }
    
    // Divider
    // Knob at the bottom of PMT baloon
    
    gMC->Gsvolu("0NB","TUBE",idtmed[6],pknob,3);
    z=-pdivider[2]+pknob[2];
    //      printf("zknob %f\n",z);
    gMC->Gspos("0NB",1,"0DIV",0,0,z,0,"ONLY");
    gMC->Gsvolu("0KB","TUBE",idtmed[kGlass],pknob_bot,3);
    z=-pdivider[2]+2*pknob[2]+pknob_bot[2];
    //      printf(knobbot %f\n",z);
    gMC->Gspos("0KB",1,"0DIV ",0,0,z,0,"ONLY");
    gMC->Gsvolu("0VAC","TUBE",idtmed[kVac],pknob_vac,3);
    z=-pdivider[2]+pknob_vac[2];
    //      printf("knobvac %f\n",z);
    gMC->Gspos("0VAC",1,"0DIV",0,0,z,0,"ONLY");
    //Steel pins + pin holes
    gMC->Gsvolu("0PIN","TUBE",idtmed[kSteel],ppins,3);
    z=-pdivider[2]+ppins[2];
    gMC->Gspos("0PIN",1,"0DIV",0,0,z,0,"ONLY");
    gMC->Gsvolu("0HOL","TUBE",idtmed[kBrass],phole,3);
    z=-pdivider[2]+2*ppins[2]+phole[2];
    gMC->Gspos("0HOL",1,"0DIV",0,0,z,0,"ONLY");
    
    //Socket
    gMC->Gsvolu("0V1","TUBE",idtmed[kCer],pdiv1,3);
    z=-pdivider[2]+pdiv1[2];
    gMC->Gspos("0V1",1,"0DIV",0,0,z,0,"ONLY");
    //Resistors
    gMC->Gsvolu("0V2","TUBE",idtmed[kAir],pdiv2,3);
    z=pdivider[2]-pdiv2[2];
    gMC->Gspos("0V2",1,"0DIV",0,0,z,0,"ONLY");
    gMC->Gsvolu("0RS","TUBE",idtmed[kCer],presist,3);
    z=-pdiv2[2]+presist[2];
    gMC->Gspos("0RS",1,"0V2",0,0,z,0,"ONLY");
    gMC->Gsvolu("0RB","TUBE",idtmed[kRibber],pribber,3);
    z=pdiv2[2]-pribber[2];
    gMC->Gspos("0RB",1,"0V2",0,0,z,0,"ONLY");


    //Support  left side

    z=-pstart[2]+psupport1[2];
    gMC->Gspos("0SU1",2,"0STL",0,0,z,0,"ONLY"); //C kozhuh snaruzhi
    gMC->Gspos("0SU2",2,"0STL",0,0,z,0,"ONLY"); //C kozhuh vnutri
    z=-pstart[2]+psupport3[2];
    gMC->Gspos("0SU3",2,"0STL",0,0,z,0,"ONLY"); //peredniaia kryshka
    z=-pstart[2]+2.*psupport1[2];
    gMC->Gspos("0SU4",2,"0STL",0,0,z,0,"ONLY"); //zadnaiai kryshka
    //Dyrki dlia feu v zadnej kryshke
    z=0;
    
    for (is=1; is<=12; is++)
      {  
	x=6.5*TMath::Sin(is*angel);
	y=6.5 *TMath::Cos(is*angel);
	gMC->Gspos("0SU5",is+12,"0SU4",x,y,z,0,"ONLY");

      }	

    //Support  right side

    z=-pstart[2]+psupport1[2];
    gMC->Gspos("0SU1",1,"0STR",0,0,z,0,"ONLY"); //C kozhuh snaruzhi
    gMC->Gspos("0SU2",1,"0STR",0,0,z,0,"ONLY"); //C kozhuh vnutri
    z=-pstart[2]+psupport3[2];
    gMC->Gspos("0SU3",1,"0STR",0,0,z,0,"ONLY"); //peredniaia kryshka
    z=-pstart[2]+2.*psupport1[2];
    gMC->Gspos("0SU4",1,"0STR",0,0,z,0,"ONLY"); //zadnaiai kryshka
    //Dyrki dlia feu v zadnej kryshke
    z=0;
    
    for (is=1; is<=12; is++)
      {  
	x=(6.5+7.8*TMath::Tan(5.4*3.1415/180)) *TMath::Sin(is*angel);
	y=(6.5+7.8*TMath::Tan(5.4*3.1415/180)) *TMath::Cos(is*angel);
	gMC->Gspos("0SU5",is,"0SU4",x,y,z,0,"ONLY");

      }	

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
    par[0]=5.1;
    par[1]=5.2;
    par[2]=5.8/2;
    gMC->Gsvolu("0SC2","TUBE",idtmed[kC],par,3);
    z += par[2];
    gMC->Gspos("0SC2",1,"0SUP",0,0,z,0,"ONLY"); 
    z += par[2];
    Float_t parC[5];
    parC[0]=0.25;
    parC[1]=5.1;
    parC[2]=5.2;
    parC[3]=5.5;
    parC[4]=5.6;
    gMC->Gsvolu("0SC3","CONE",idtmed[kC],parC,5);
    z += parC[0];
    gMC->Gspos("0SC3",1,"0SUP",0,0,z,0,"ONLY"); 
    z += parC[0];
    par[0]=5.5;
    par[1]=5.6;
    par[2]=1.2/2;
    gMC->Gsvolu("0SC4","TUBE",idtmed[kC],par,3);
    z += par[2];
    gMC->Gspos("0SC4",1,"0SUP",0,0,z,0,"ONLY"); 
    par[0]=5.1;
    par[1]=5.5;
    par[2]=1.2/2;
    gMC->Gsvolu("0SA0","TUBE",idtmed[kAl],par,3);
    gMC->Gspos("0SA0",1,"0SUP",0,0,z,0,"ONLY"); 
    //gvozdi dlia skruchivaniia Al i C parts
    par[0]=5.75; 
    par[1]=5.78;
    gMC->Gsvolu("0SN1","TUBE",idtmed[kSteel],par,3);
    gMC->Gspos("0SN1",1,"0SUP",0,0,z,0,"ONLY"); 
    z += par[2];
    par[0]=3.15;
    par[1]=5.5;
    par[2]=0.1;
    gMC->Gsvolu("0SA1","TUBE",idtmed[kAl],par,3);
    z += par[2];
    gMC->Gspos("0SA1",1,"0SUP",0,0,z,0,"ONLY"); 
    z=z+par[2];
    par[0]=3.15;
    par[1]=3.16;
    par[2]=4.7/2;
    gMC->Gsvolu("0SA2","TUBE",idtmed[kAl],par,3);
    z += par[2];
    gMC->Gspos("0SA2",1,"0SUP",0,0,z,0,"ONLY"); 
    z=z+par[2];
    par[0]=3.16; // eta chast' prikruchena k absorberu
    par[1]=7.5;
    par[2]=0.3;
    gMC->Gsvolu("0SA3","TUBE",idtmed[kAl],par,3);
    z += par[2];
    gMC->Gspos("0SA3",1,"0SUP",0,0,z,0,"ONLY"); 
    par[0]=3.16; // gvozdi eta chast' prikruchena k absorberu
    par[1]=7.5;
    par[2]=0.01;
    gMC->Gsvolu("0SN2","TUBE",idtmed[kSteel],par,3);
    gMC->Gspos("0SN2",1,"0SUP",0,0,z,0,"ONLY"); 
       
 
}    
//------------------------------------------------------------------------
void AliSTARTv0::CreateMaterials()
{
   Int_t isxfld   = gAlice->Field()->Integ();
   Float_t sxmgmx = gAlice->Field()->Max();
   Float_t a,z,d,radl,absl,buf[1];
   Int_t nbuf;

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
   Float_t acer[2],zcer[2],wcer[2]={0.972,0.028};
   Float_t aal2o3[2]  = { 26.981539,15.9994 };
   Float_t zal2o3[2]  = { 13.,8. };
   Float_t wal2o3[2]  = { 2.,3. };
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
// AIR
   Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
   Float_t zAir[4]={6.,7.,8.,18.};
   Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
   Float_t dAir = 1.20479E-3;
     
//*** Definition Of avaible START materials ***
   AliMaterial(0, "START Steel$", 55.850,26.,7.87,1.76,999);
   AliMaterial(1, "START Vacuum$", 1.e-16,1.e-16,1.e-16,1.e16,999);
   AliMixture(2, "START Air$", aAir, zAir, dAir,4,wAir);

   AliMaterial(10, "CarbonPlastic$", 12.01, 6.0, 2.26, 18.8,999); 
   AliMaterial(11, "Aliminium$", 26.98, 13.0, 2.7, 8.9,999); 

   AliMixture( 3, "Al2O3   $", aal2o3, zal2o3, denscer, -2, wal2o3);
   AliMixture( 4, "PMT glass   $",aglass,zglass,dglass,-2,wglass);
   char namate[21]="";
   gMC->Gfmate((*fIdmate)[3], namate, a, z, d, radl, absl, buf, nbuf);
   acer[0]=a;
   zcer[0]=z;
   gMC->Gfmate((*fIdmate)[4], namate, a, z, d, radl, absl, buf, nbuf);
   acer[1]=a;
   zcer[1]=z;
   
   AliMixture( 9, "Ceramic    $", acer, zcer, denscer, 2, wcer);
   AliMixture( 5, "Scintillator$",ascin,zscin,denscin,-2,wscin);
   AliMixture( 6, "Brass    $", abrass, zbrass, denbrass, 2, wbrass);
   
   AliMixture( 7, "Ribber $",aribber,zribber,denribber,-3,wribber);
   AliMixture( 8, "Lucite$",alucite,zlucite,denlucite,-3,wlucite);
   AliMixture( 9, "Penoplast$",asupport,zsupport,densupport,-2,wsupport);
 
   
   AliMedium(1, "START Air$", 2, 0, isxfld, sxmgmx, 10., .1, 1., .003, .003);
   AliMedium(2, "Scintillator$", 5, 1, isxfld, sxmgmx, 10., .01, 1., .003, .003);
   AliMedium(3, "Vacuum$", 1, 0, isxfld, sxmgmx, 10., .01, .1, .003, .003);
   AliMedium(4, "Ceramic$", 9, 0, isxfld, sxmgmx, 10., .01, .1, .003, .003);
   AliMedium(6, "Glass$", 4, 1, isxfld, sxmgmx, 10., .01, .1, .003, .003);
   AliMedium(8, "Steel$", 0, 0, isxfld, sxmgmx, 1., .001, 1., .001, .001);
   AliMedium(9, "Ribber  $", 7, 0, isxfld, sxmgmx, 10., .01, .1, .003, .003);
   AliMedium(11, "Brass  $", 6, 0, isxfld, sxmgmx, 10., .01, .1, .003, .003);
   AliMedium(12, "Lucite$", 8, 1, isxfld, sxmgmx, 10., .01, 1., .003, .003);  
   AliMedium(13, "CarbonPlastic$", 10, 0, isxfld, sxmgmx, 10., .01, 1., .003, .003);  
   AliMedium(14, "PenoPlast$", 9, 0, isxfld, sxmgmx, 10., .01, 1., .003, .003);  
   AliMedium(15, "Aluminium$", 11, 0, isxfld, sxmgmx, 10., .01, 1., .003, .003);  

   if(fDebug) cout<<ClassName()<<": ++++++++++++++Medium set++++++++++"<<endl;

}
//---------------------------------------------------------------------
void AliSTARTv0::DrawModule()
{
//
// Draw a shaded view of the Forward multiplicity detector version 0
//
  
  //Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN",0);
  //
  //Set volumes visible
  gMC->Gsatt("0STA","SEEN",0);
  gMC->Gsatt("0PMT","SEEN",1);
  gMC->Gsatt("0DIV","SEEN",1);
  //
  gMC->Gdopt("hide","on");
  gMC->Gdopt("shad","on");
  gMC->SetClipBox(".");
  gMC->SetClipBox("*",0,1000,-1000,1000,-1000,1000);
  gMC->DefaultRange();
  gMC->Gdraw("alic",40,30,0,12,9.5,.7,0.7);
  gMC->Gdhead(1111,"T-Zero detector");
  gMC->Gdopt("hide","off");
}

//-------------------------------------------------------------------
void AliSTARTv0::Init()
{
// Initialises version 0 of the Forward Multiplicity Detector
//
//Int_t *idtmed  = gAlice->Idtmed();
  AliSTART::Init();
  fIdSens1=gMC->VolId("0TOP");
  printf("*** START version 0 initialized ***\n");
 
}

//-------------------------------------------------------------------

void AliSTARTv0::StepManager()
{
  //
  // Called for every step in the START Detector
  // See AliSTARTv1

}
//---------------------------------------------------------------------










