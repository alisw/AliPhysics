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

// Authors: David Guez, Ivana Hrivnacova, Marion MacCormick; IPN Orsay
//
// Class AliMUONv2
// ---------------
// Inherits from AliMUONv1 but with a more detailed
// geometrical description of station 1 

#include <algorithm>
#include <string>

#include <TVector2.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TArrayI.h>
#include <Riostream.h>
#include <TSystem.h>

#include "AliMpFiles.h"
#include "AliMpReader.h"
#include "AliMpSector.h"
#include "AliMpRow.h"
#include "AliMpVRowSegment.h"
#include "AliMpMotifMap.h"
#include "AliMpMotifPosition.h"

#include "AliMUONv2.h"
#include "AliMUONConstants.h"
#include "AliMUONHit.h"
#include "AliRun.h"
#include "AliMagF.h"
#include "AliConst.h" 
#include "AliMUONChamber.h"

ClassImp(AliMUONv2)

// Thickness Constants
const GReal_t AliMUONv2::fgkHzPadPlane=0.0148/2.;     //Pad plane
const GReal_t AliMUONv2::fgkHzFoam = 2.083/2.;        //Foam of mechanicalplane
const GReal_t AliMUONv2::fgkHzFR4 = 0.0031/2.;        //FR4 of mechanical plane
const GReal_t AliMUONv2::fgkHzSnPb = 0.0091/2.;       //Pad/Kapton connection (66 pt)
const GReal_t AliMUONv2::fgkHzKapton = 0.0122/2.;     //Kapton
const GReal_t AliMUONv2::fgkHzBergPlastic = 0.3062/2.;//Berg connector
const GReal_t AliMUONv2::fgkHzBergCopper = 0.1882/2.; //Berg connector
const GReal_t AliMUONv2::fgkHzDaughter = 0.0156/2.;   //Daughter board
const GReal_t AliMUONv2::fgkHzGas = 0.2/2.;           //Gas thickness

// Quadrant Mother volume - TUBS1 - Middle layer of model
const GReal_t AliMUONv2::fgkMotherIR1 = 18.3;
const GReal_t AliMUONv2::fgkMotherOR1 = 105.673;   
const GReal_t AliMUONv2::fgkMotherThick1 = 6.5/2;  
const GReal_t AliMUONv2::fgkMotherPhiL1 = 0.; 
const GReal_t AliMUONv2::fgkMotherPhiU1 = 90.;

// Quadrant Mother volume - TUBS2 - near and far layers of model
const GReal_t AliMUONv2::fgkMotherIR2 = 20.7;   
const GReal_t AliMUONv2::fgkMotherOR2 = 100.073;   
const GReal_t AliMUONv2::fgkMotherThick2 = 3.0/2; 
const GReal_t AliMUONv2::fgkMotherPhiL2 = 0.; 
const GReal_t AliMUONv2::fgkMotherPhiU2 = 90.;

// Sensitive copper pads, foam layer, PCB and electronics model parameters
const GReal_t AliMUONv2::fgkHxHole=1.5/2.;
const GReal_t AliMUONv2::fgkHyHole=6./2.;
const GReal_t AliMUONv2::fgkHxBergPlastic=0.74/2.;
const GReal_t AliMUONv2::fgkHyBergPlastic=5.09/2.;
const GReal_t AliMUONv2::fgkHxBergCopper=0.25/2.;
const GReal_t AliMUONv2::fgkHyBergCopper=3.6/2.;
const GReal_t AliMUONv2::fgkHxKapton=0.8/2.;
const GReal_t AliMUONv2::fgkHyKapton=5.7/2.;
const GReal_t AliMUONv2::fgkHxDaughter=2.3/2.;
const GReal_t AliMUONv2::fgkHyDaughter=6.3/2.;
const GReal_t AliMUONv2::fgkOffsetX=1.46;
const GReal_t AliMUONv2::fgkOffsetY=0.71;
const GReal_t AliMUONv2::fgkDeltaFilleEtamX=1.46;
const GReal_t AliMUONv2::fgkDeltaFilleEtamY=0.051;

const GReal_t AliMUONv2::fgkDeltaQuadLHC=2.6;  // LHC Origin wrt Quadrant Origin
const GReal_t AliMUONv2::fgkFrameOffset=5.0;  

const char* AliMUONv2::fgkHoleName="MCHL";      
const char* AliMUONv2::fgkDaughterName="MCDB";  
const char  AliMUONv2::fgkFoamLayerSuffix='F';  // prefix for automatic volume naming
const char* AliMUONv2::fgkQuadrantMLayerName="SQM";
const char* AliMUONv2::fgkQuadrantNLayerName="SQN";
const char* AliMUONv2::fgkQuadrantFLayerName="SQF";

//______________________________________________________________________________
AliMUONv2::AliMUONv2()
  : AliMUONv1()
{
// Default Constructor
// --
   fChamberV2[0] = 0;
   fChamberV2[1] = 0;
   
   // keep secondaries
   SetIshunt(0);
 
   // set path to mapping data files
   if (! gSystem->Getenv("MINSTALL")) {    
     TString dirPath = gSystem->Getenv("ALICE_ROOT");
     dirPath += "/MUON/mapping"; 
     AliMpFiles::Instance()->SetTopPath(dirPath);
     gSystem->Setenv("MINSTALL", dirPath.Data());
     //cout << "AliMpFiles top path set to " << dirPath << endl;	  
   }
   //else
   //  cout << gSystem->Getenv("MINSTALL") << endl;	  
}
 
//______________________________________________________________________________
AliMUONv2::AliMUONv2(const char *name, const char *title)
  : AliMUONv1(name,title)
{
   fChamberV2[0] = 0;
   fChamberV2[1] = 0;
   
   // keep secondaries
   SetIshunt(0);
   
   // set path to mapping data files
   if (! gSystem->Getenv("MINSTALL")) {    
     TString dirPath = gSystem->Getenv("ALICE_ROOT");
     dirPath += "/MUON/mapping"; 
     AliMpFiles::Instance()->SetTopPath(dirPath);
     gSystem->Setenv("MINSTALL", dirPath.Data());
     //cout << "AliMpFiles top path set to " << dirPath << endl;	  
   }
   //else
   //  cout << gSystem->Getenv("MINSTALL") << endl;	    	  
}
 
//______________________________________________________________________________
AliMUONv2::AliMUONv2(const AliMUONv2& rMUON):AliMUONv1(rMUON)
{
// Dummy copy constructor
}

//______________________________________________________________________________
AliMUONv2::~AliMUONv2()
{
// Destructor
}

//
//  Private methods
//

//______________________________________________________________________________
void AliMUONv2::CreateHole()
{
// Create all the elements found inside a foam hole
// --
  Int_t* idtmed = fIdtmed->GetArray()-1099;
  Int_t idAir  = idtmed[1100];    // medium 1
  Int_t idCopper  = idtmed[1109]; // medium 10 = copper

  GReal_t par[3];
  GReal_t posX,posY,posZ;
  
  par[0] = fgkHxHole;
  par[1] = fgkHyHole;
  par[2] = fgkHzFoam;
  gMC->Gsvolu(fgkHoleName,"BOX",idAir,par,3);

  par[0] = fgkHxKapton;
  par[1] = fgkHyKapton;
  par[2] = fgkHzSnPb;
  gMC->Gsvolu("SNPB", "BOX", idCopper, par, 3);
  posX = 0.;
  posY = 0.;
  posZ = -fgkHzFoam+fgkHzSnPb;
  gMC->Gspos("SNPB",1,fgkHoleName, posX, posY, posZ, 0,"ONLY");

  par[0] = fgkHxHole;
  par[1] = fgkHyBergPlastic;
  par[2] = fgkHzKapton;
  gMC->Gsvolu("KAPT", "BOX", idCopper, par, 3);
  posX = 0.;
  posY = 0.;
  posZ = 0.;
  gMC->Gspos("KAPT",1,fgkHoleName, posX, posY, posZ, 0,"ONLY");
}

//______________________________________________________________________________
void AliMUONv2::CreateDaughterBoard()
{
// Create all the elements in a daughter board
// --
  Int_t* idtmed = fIdtmed->GetArray()-1099;
  Int_t idAir  = idtmed[1100]; // medium 1
  Int_t idCopper  = idtmed[1109]; // medium 10 = copper
  Int_t idPlastic  =idtmed[1116]; // medium 17 = Plastic

  GReal_t par[3];
  GReal_t posX,posY,posZ;

  par[0]=fgkHxDaughter;
  par[1]=fgkHyDaughter;
  par[2]=TotalHzDaughter();
  gMC->Gsvolu(fgkDaughterName,"BOX",idAir,par,3);
  
  par[0]=fgkHxBergPlastic;
  par[1]=fgkHyBergPlastic;
  par[2]=fgkHzBergPlastic;
  gMC->Gsvolu("BRGP","BOX",idPlastic,par,3);
  posX=0.;
  posY=0.;
  posZ = -TotalHzDaughter() + fgkHzBergPlastic;
  gMC->Gspos("BRGP",1,fgkDaughterName,posX,posY,posZ,0,"ONLY");

  par[0]=fgkHxBergCopper;
  par[1]=fgkHyBergCopper;
  par[2]=fgkHzBergCopper;
  gMC->Gsvolu("BRGC","BOX",idCopper,par,3);
  posX=0.;
  posY=0.;
  posZ=0.;
  gMC->Gspos("BRGC",1,"BRGP",posX,posY,posZ,0,"ONLY");

  par[0]=fgkHxDaughter;
  par[1]=fgkHyDaughter;
  par[2]=fgkHzDaughter;
  gMC->Gsvolu("DGHT","BOX",idCopper,par,3);
  posX=0.;
  posY=0.;
  posZ = -TotalHzDaughter() + 2.*fgkHzBergPlastic + fgkHzDaughter;
  gMC->Gspos("DGHT",1,fgkDaughterName,posX,posY,posZ,0,"ONLY");
}

//______________________________________________________________________________
void AliMUONv2::CreateInnerLayers()
{
// Create the layer of sensitive volumes with gas
// and the copper layer.
// --

// Gas Medium
  Int_t* idtmed = fIdtmed->GetArray()-1099; 
  Int_t idArCO2  = idtmed[1108];  // medium 9 (ArCO2 80%) 
  Int_t idCopper  = idtmed[1109]; // medium 10 = copper

  Float_t par[11];

//Make gas volume - composed of 11 trapezoids
// section 1 of 11
    par[0] = fgkHzGas;
    par[1] = 0.;
    par[2] = 0.;
    par[3] = 71.33/2.;
    par[4] = 9.76/2.;
    par[5] = 48.77/2.;
    par[6] = 15.3;
    par[7] = 71.33/2.;
    par[8] = 9.76/2.;
    par[9] = 48.77/2.;
    par[10] = 15.3;        

  gMC->Gsvolu("SA1G", "TRAP", idArCO2, par, 11);
  gMC->Gsvolu("SA2G", "TRAP", idArCO2, par, 11);
  
  par[0] = fgkHzPadPlane;
  gMC->Gsvolu("SA1C", "TRAP", idCopper,par, 11);

// section 2 of 11  
    par[0] = fgkHzGas;
    par[1] = 0.;
    par[2] = 0.;
    par[3] = 79.68/2.;
    par[4] = 10.4/2.;
    par[5] = 57.0/2.;
    par[6] = 0.;  
    par[7] = 79.68/2.; 
    par[8] = 10.4/2.;
    par[9] = 57.0/2.;
    par[10] = 0.;  
  gMC->Gsvolu("SB1G", "TRAP", idArCO2, par, 11);  
  gMC->Gsvolu("SB2G", "TRAP", idArCO2, par, 11);

  par[0] = fgkHzPadPlane;
  gMC->Gsvolu("SB1C", "TRAP", idCopper,par, 11);

// section 3 of 11
    par[0] = fgkHzGas;
    par[1] = 0.;
    par[2] = 0.;
    par[3] = 71.33/2.;
    par[4] = 48.77/2.;
    par[5] = 9.73/2.;
    par[6] = -15.3;
    par[7] = 71.33/2.;
    par[8] = 48.77/2.;
    par[9] = 9.73/2.;
    par[10] = -15.3;   
 
  gMC->Gsvolu("SC1G", "TRAP", idArCO2, par, 11);  
  gMC->Gsvolu("SC2G", "TRAP", idArCO2, par, 11);

  par[0] = fgkHzPadPlane;
  gMC->Gsvolu("SC1C", "TRAP", idCopper,par, 11);

// section 4 of 11
    par[0] = fgkHzGas;
    par[1] = 0.;
    par[2] = 0.;
    par[3] = 6.00/2.;
    par[4] = 0.;
    par[5] = 1.56/2.;
    par[6] = 7.41; 
    par[7] = 6.00/2.; 
    par[8] = 0.;
    par[9] = 1.56/2.;
    par[10] = 7.41;    
  gMC->Gsvolu("SD1G", "TRAP", idArCO2, par, 11);  
  gMC->Gsvolu("SD2G", "TRAP", idArCO2, par, 11);

  par[0] = fgkHzPadPlane;
  gMC->Gsvolu("SD1C", "TRAP", idCopper,par, 11);

// section 5 of 11  
    par[0] = fgkHzGas;
    par[1] = 0.;
    par[2] = 0.;
    par[3] = 1.516/2.;
    par[4] = 0.;
    par[5] = 0.829/2.;
    par[6] = 15.3;
    par[7] = 1.516/2.;
    par[8] = 0.;
    par[9] = 0.829/2.;
    par[10] = 15.3;   
  gMC->Gsvolu("SE1G", "TRAP", idArCO2, par, 11);  
  gMC->Gsvolu("SE2G", "TRAP", idArCO2, par, 11);

  par[0] = fgkHzPadPlane;
  gMC->Gsvolu("SE1C", "TRAP", idCopper,par, 11);

// section 6 of 11
    par[0] = fgkHzGas;
    par[1] = 0.;
    par[2] = 0.;
    par[3] = 3.92/2.;
    par[4] = 0.;
    par[5] = 0.562/2.;
    par[6] = -4.1;
    par[7] = 3.92/2.;
    par[8] = 0.;
    par[9] = 0.562/2.;
    par[10] = -4.1;   
  gMC->Gsvolu("SF1G", "TRAP", idArCO2, par, 11);  
  gMC->Gsvolu("SF2G", "TRAP", idArCO2, par, 11);
    
  par[0] = fgkHzPadPlane;
  gMC->Gsvolu("SF1C", "TRAP", idCopper,par, 11);

// section 7 of 11
    par[0] = fgkHzGas;
    par[1] = 0.;
    par[2] = 0.;
    par[3] = 0.941/2.;
    par[4] = 0.562/2.;
    par[5] = 0.;
    par[6] = -16.6; 
    par[7] = 0.941/2.;
    par[8] = 0.562/2.;
    par[9] = 0.;
    par[10] =-16.6;    
  gMC->Gsvolu("SG1G", "TRAP", idArCO2, par, 11);  
  gMC->Gsvolu("SG2G", "TRAP", idArCO2, par, 11);

  par[0] = fgkHzPadPlane;
  gMC->Gsvolu("SG1C", "TRAP", idCopper,par, 11);

// section 8 of 11
    par[0] = fgkHzGas;
    par[1] = 0.;
    par[2] = 0.;
    par[3] = 3.94/2.;
    par[4] = 0.57/2.;
    par[5] = 0.;
    par[6] = 4.14; 
    par[7] = 3.94/2.; 
    par[8] = 0.57/2.;
    par[9] = 0.;
    par[10] = 4.14;    
  gMC->Gsvolu("SH1G", "TRAP", idArCO2, par, 11);  
  gMC->Gsvolu("SH2G", "TRAP", idArCO2, par, 11);

  par[0] = fgkHzPadPlane;
  gMC->Gsvolu("SH1C", "TRAP", idCopper,par, 11);

// section 9 of 11  
    par[0] = fgkHzGas;
    par[1] = 0.;
    par[2] = 0.;
    par[3] = 0.95/2.;
    par[4] = 0.;
    par[5] = 0.57/2;
    par[6] = 16.7;
    par[7] = 0.95/2.;
    par[8] = 0.;
    par[9] = 0.57/2;
    par[10] = 16.7;   
  gMC->Gsvolu("SI1G", "TRAP", idArCO2, par, 11);  
  gMC->Gsvolu("SI2G", "TRAP", idArCO2, par, 11);

  par[0] = fgkHzPadPlane;
  gMC->Gsvolu("SI1C", "TRAP", idCopper,par, 11);

// section 10 of 11
    par[0] = fgkHzGas;
    par[1] = 0.;
    par[2] = 0.;
    par[3] = 1.49/2.;
    par[4] = 0.;
    par[5] = 0.817/2.;
    par[6] = -15.4;
    par[7] = 1.49/2.;
    par[8] = 0.;
    par[9] = 0.817/2.;
    par[10] = -15.4;   
  gMC->Gsvolu("SJ1G", "TRAP", idArCO2, par, 11);  
  gMC->Gsvolu("SJ2G", "TRAP", idArCO2, par, 11);
    
  par[0] = fgkHzPadPlane;
  gMC->Gsvolu("SJ1C", "TRAP", idCopper,par, 11);

// section 11 of 11
    par[0] = fgkHzGas;
    par[1] = 0.;
    par[2] = 0.;
    par[3] = 5.93/2.;
    par[4] = 0.;
    par[5] = 1.49/2.;
    par[6] = -7.16; 
    par[7] = 5.93/2.;
    par[8] = 0.;
    par[9] = 1.49/2.;
    par[10] = -7.16;    
  gMC->Gsvolu("SK1G", "TRAP", idArCO2, par, 11);  
  gMC->Gsvolu("SK2G", "TRAP", idArCO2, par, 11);

  par[0] = fgkHzPadPlane;
  gMC->Gsvolu("SK1C", "TRAP", idCopper,par, 11);
}

//______________________________________________________________________________
void AliMUONv2::CreateQuadrant(Int_t chamber)
{
// create the quadrant (bending and non-bending planes)
// for the given chamber
// --

  CreateFrame(chamber);

  TSpecialMap specialMap;
  specialMap[1001] = AliMUONSt1SpecialMotif(TVector2( 0.1, 0.84), 90.);
  specialMap[1002] = AliMUONSt1SpecialMotif(TVector2( 0.5, 0.36));
  specialMap[1003] = AliMUONSt1SpecialMotif(TVector2(1.01, 0.36));
  AliMpReader reader1(kStation1, kBendingPlane);
  AliMpSector* sector1 = reader1.BuildSector();

  Bool_t reflectZ = true;
  TVector3 where = TVector3(2.5+0.1+0.56+0.001, 2.5+0.1+0.001, 0.);
  PlaceSector(sector1, specialMap, where, reflectZ, chamber);
  
  specialMap.clear();
  specialMap[4001] = AliMUONSt1SpecialMotif(TVector2(1.01,0.59),90.);
  specialMap[4002] = AliMUONSt1SpecialMotif(TVector2(1.96, 0.17));
  specialMap[4003] = AliMUONSt1SpecialMotif(TVector2(1.61,-1.18));
  specialMap[4004] = AliMUONSt1SpecialMotif(TVector2(0.2 ,-0.08));
  specialMap[4005] = AliMUONSt1SpecialMotif(TVector2(0.2 , 0.25));
  specialMap[4006] = AliMUONSt1SpecialMotif(TVector2(0.28, 0.21));
  AliMpReader reader2(kStation1, kNonBendingPlane);
  AliMpSector* sector2 = reader2.BuildSector();

  reflectZ = false;
  where = TVector3(where.X()+0.63/2.,where.Y()+0.42/2., 0.); //add a half pad shift
  PlaceSector(sector2, specialMap, where, reflectZ, chamber);
}

//______________________________________________________________________________
void AliMUONv2::CreateFoamBox(const char* name,const  TVector2& dimensions)
{
// create all the elements in the copper plane
// --

  Int_t* idtmed = fIdtmed->GetArray()-1099;
  Int_t idAir  = idtmed[1100]; // medium 1
  Int_t idFoam = idtmed[1115]; // medium 16 = Foam
  Int_t idFR4  = idtmed[1114]; // medium 15 = FR4

  // mother volume
  GReal_t par[3];
  par[0] = dimensions.X();
  par[1] = dimensions.Y();
  par[2] = TotalHzPlane();
  gMC->Gsvolu(name,"BOX",idAir,par,3);
  
  // foam layer
  GReal_t posX,posY,posZ;
  char eName[5];
  strcpy(eName,name);
  eName[3]=fgkFoamLayerSuffix;
  par[0] = dimensions.X();
  par[1] = dimensions.Y();
  par[2] = fgkHzFoam;
  gMC->Gsvolu(eName,"BOX",idFoam,par,3);
  posX=0.;
  posY=0.;
  posZ = -TotalHzPlane() + fgkHzFoam;
  gMC->Gspos(eName,1,name,posX,posY,posZ,0,"ONLY");

  // mechanical plane FR4 layer
  eName[3]='R';
  par[0] = dimensions.X();
  par[1] = dimensions.Y();
  par[2] = fgkHzFR4;
  gMC->Gsvolu(eName,"BOX",idFR4,par,3);
  posX=0.;
  posY=0.;
  posZ = -TotalHzPlane()+ 2.*fgkHzFoam + fgkHzFR4;
  gMC->Gspos(eName,1,name,posX,posY,posZ,0,"ONLY");
}

//______________________________________________________________________________
void AliMUONv2::CreatePlaneSegment(const char* name,const  TVector2& dimensions,
      	      	      	      	   Int_t nofHoles)
{
// Create a segment of a plane (this includes a foam layer, 
// holes in the foam to feed the kaptons through, kapton connectors
// and the mother board.)
// --
  
  CreateFoamBox(name,dimensions);

  char eName[5];
  strcpy(eName,name);
  eName[3]=fgkFoamLayerSuffix;
  
  for (Int_t holeNum=0;holeNum<nofHoles;holeNum++) {
    GReal_t posX = ((2.*holeNum+1.)/nofHoles-1.)*dimensions.X();
    GReal_t posY = 0.;
    GReal_t posZ = 0.;
  
    gMC->Gspos(fgkHoleName,holeNum+1,eName,posX,posY,posZ,0,"ONLY");
  }
}

//______________________________________________________________________________
void AliMUONv2::CreateFrame(Int_t chamber)
{
// Create the non-sensitive elements of the frame for the  <chamber>
//
// 
// Model and notation:
//
// The Quadrant volume name starts with SQ
// The volume segments are numbered 00 to XX.
//
//                              OutTopFrame
//                               (SQ02-16) 
//                              ------------  
//             OutEdgeFrame   /              |
//             (SQ17-24)     /               |  InVFrame (SQ00-01) 
//                          /                |
//                          |                |   
//               OutVFrame  |            _- - 
//               (SQ25-39)  |           |   InArcFrame (SQ42-45)
//                          |           |
//                          -------------
//                        InHFrame (SQ40-41)
//                          
//
// 06 February 2003 - Overlapping volumes resolved.
// One quarter chamber is comprised of three TUBS volumes: SQMx, SQNx, and SQFx,
// where SQMx is the Quadrant Middle layer for chamber <x> ( posZ in [-3.25,3.25]),
// SQNx is the Quadrant Near side layer for chamber <x> ( posZ in [-6.25,3-.25) ), and
// SQFx is the Quadrant Far side layer for chamber <x> ( posZ in (3.25,6.25] ).
//---

const Float_t fgkNearFarLHC=2.4;    // Near and Far TUBS Origin wrt LHC Origin

  // tracking medias
  Int_t* idtmed = fIdtmed->GetArray()-1099;
  
  Int_t idAir  = idtmed[1100];       // medium 1
  Int_t idFrameEpoxy = idtmed[1115]; // medium 16 = Frame Epoxy ME730
  Int_t idInox = idtmed[1116];       // medium 17 Stainless Steel (18%Cr,9%Ni,Fe)
  Int_t idFR4 = idtmed[1110];        // medium 11 FR4
  Int_t idCopper = idtmed[1109];     // medium 10 Copper
  Int_t idAlu = idtmed[1103];        // medium 4 Aluminium
  
  
// Rotation Matrices  
      Int_t rot1, rot2, rot3;    
      
//   Rotation matrices  
     AliMatrix(rot1,  90.,  90., 90., 180.,  0., 0.); // +90 deg in x-y plane
     AliMatrix(rot2,  90.,  45., 90., 135.,  0., 0.); // +45 deg in x-y plane 
     AliMatrix(rot3,  90.,  45., 90., 315.,180., 0.); // +45 deg in x-y + rotation 180° around y

//   Translation matrices ... NOT USED  
//     AliMatrix(trans1, 90.,   0., 90.,  90.,   0., 0.); // X-> X; Y -> Y; Z -> Z
//     AliMatrix(trans2, 90., 180., 90.,  90., 180., 0.); // X->-X; Y -> Y; Z ->-Z
//     AliMatrix(trans3, 90., 180., 90., 270.,   0., 0.); // X->-X; Y ->-Y; Z -> Z
//     AliMatrix(trans4, 90.,   0., 90., 270., 180., 0.); // X-> X; Y ->-Y; Z ->-Z
//  
      // ___________________Volume thicknesses________________________

  const Float_t hzFrameThickness = 1.59/2.;     //equivalent thickness
  const Float_t hzOuterFrameEpoxy = 1.19/2.;    //equivalent thickness
  const Float_t hzOuterFrameInox = 0.1/2.;      //equivalent thickness
  const Float_t hzFoam = 2.083/2.;              //evaluated elsewhere
  
// Pertaining to the top outer area 
  const Float_t hzTopAnodeSteel1 = 0.185/2.;    //equivalent thickness
  const Float_t hzTopAnodeSteel2 = 0.51/2.;     //equivalent thickness  
  const Float_t hzAnodeFR4 = 0.08/2.;           //equivalent thickness
  const Float_t hzTopEarthFaceCu = 0.364/2.;    //equivalent thickness
  const Float_t hzTopEarthProfileCu = 1.1/2.;   //equivalent thickness
  const Float_t hzTopPositionerSteel = 1.45/2.; //should really be 2.125/2.; 
  const Float_t hzTopGasSupportAl = 0.85/2.;    //equivalent thickness
  
// Pertaining to the vertical outer area  
  const Float_t hzVerticalCradleAl = 0.8/2.;     //equivalent thickness
  const Float_t hzLateralSightAl = 0.975/2.;     //equivalent thickness
  const Float_t hzLateralPosnInoxFace = 2.125/2.;//equivalent thickness
  const Float_t hzLatPosInoxProfM = 6.4/2.;      //equivalent thickness
  const Float_t hzLatPosInoxProfNF = 1.45/2.;    //equivalent thickness
  const Float_t hzLateralPosnAl = 0.5/2.;        //equivalent thickness
  const Float_t hzVertEarthFaceCu = 0.367/2.;    //equivalent thickness
  const Float_t hzVertBarSteel = 0.198/2.;       //equivalent thickness
  const Float_t hzVertEarthProfCu = 1.1/2.;      //equivalent thickness

      //_______________Parameter definitions in sequence _________

// InVFrame parameters
  const Float_t hxInVFrame  = 1.85/2.;
  const Float_t hyInVFrame  = 73.95/2.;
  const Float_t hzInVFrame  = hzFrameThickness;

//Flat 7.5mm vertical section
  const Float_t hxV1mm  = 0.75/2.;
  const Float_t hyV1mm  = 1.85/2.;
  const Float_t hzV1mm  = hzFrameThickness;

// OuterTopFrame Structure 
//
// FRAME
// The frame is composed of a cuboid and two trapezoids 
// (TopFrameAnode, TopFrameAnodeA, TopFrameAnodeB). 
// Each shape is composed of two layers (Epoxy and Inox) and 
// takes the frame's inner anode circuitry into account in the material budget.
//
// ANODE
// The overhanging anode part is composed froma cuboid and two trapezoids 
// (TopAnode, TopAnode1, and TopAnode2). These surfaces neglect implanted
// resistors, but accounts for the major Cu, Pb/Sn, and FR4 material
// contributions.  
// The stainless steel anode supports have been included.
//
// EARTHING (TopEarthFace, TopEarthProfile)
// Al GAS SUPPORT (TopGasSupport)
//  
// ALIGNMENT (TopPositioner) - Alignment system, three sights per quarter 
// chamber. This sight is forseen for the alignment of the horizontal level 
// (parallel to the OY axis of LHC). Its position will be evaluated relative 
// to a system of sights places on the cradles;
//
//---
  
//TopFrameAnode parameters - cuboid, 2 layers
  const Float_t hxTFA = 34.1433/2.;
  const Float_t hyTFA = 7.75/2.;
  const Float_t hzTFAE = hzOuterFrameEpoxy;     // layer 1 thickness
  const Float_t hzTFAI = hzOuterFrameInox;      // layer 3 thickness
  
// TopFrameAnodeA parameters - trapezoid, 2 layers
  const Float_t hzFAAE = hzOuterFrameEpoxy;     // layer 1 thickness
  const Float_t hzFAAI = hzOuterFrameInox;      // layer 3 thickness
  const Float_t tetFAA = 0.;
  const Float_t phiFAA = 0.;
  const Float_t h1FAA = 8.7/2.;
  const Float_t bl1FAA = 4.35/2.;
  const Float_t tl1FAA =  7.75/2.;
  const Float_t alp1FAA = 11.06; 
  const Float_t h2FAA = 8.7/2.;
  const Float_t bl2FAA = 4.35/2.;
  const Float_t tl2FAA = 7.75/2.;
  const Float_t alp2FAA = 11.06;  
  
// TopFrameAnodeB parameters - trapezoid, 2 layers
  const Float_t hzFABE = hzOuterFrameEpoxy;     // layer 1 thickness
  const Float_t hzFABI = hzOuterFrameInox;      // layer 3 thickness
  const Float_t tetFAB = 0.;
  const Float_t phiFAB = 0.;
  const Float_t h1FAB = 8.70/2.;
  const Float_t bl1FAB = 0.;
  const Float_t tl1FAB = 4.35/2.;
  const Float_t alp1FAB = 14.03; 
  const Float_t h2FAB = 8.70/2.;
  const Float_t bl2FAB = 0.;
  const Float_t tl2FAB = 4.35/2.;
  const Float_t alp2FAB = 14.03;  
  
// TopAnode parameters - cuboid (part 1 of 3 parts)
  const Float_t hxTA1 = 16.2/2.;
  const Float_t hyTA1 = 3.5/2.;
  const Float_t hzTA11 = hzTopAnodeSteel1;   // layer 1
  const Float_t hzTA12 = hzAnodeFR4;         // layer 2 

// TopAnode parameters - trapezoid 1 (part 2 of 3 parts)
  const Float_t hzTA21 = hzTopAnodeSteel2;   // layer 1 
  const Float_t hzTA22 = hzAnodeFR4;         // layer 2 
  const Float_t tetTA2 = 0.;
  const Float_t phiTA2= 0.;
  const Float_t h1TA2 = 7.268/2.;
  const Float_t bl1TA2 = 2.03/2.;
  const Float_t tl1TA2 = 3.5/2.;
  const Float_t alp1TA2 = 5.78; 
  const Float_t h2TA2 = 7.268/2.;
  const Float_t bl2TA2 = 2.03/2.;
  const Float_t tl2TA2 = 3.5/2.;
  const Float_t alp2TA2 = 5.78;  

// TopAnode parameters - trapezoid 2 (part 3 of 3 parts)
  const Float_t hzTA3 = hzAnodeFR4;       // layer 1 
  const Float_t tetTA3 = 0.;
  const Float_t phiTA3 = 0.;
  const Float_t h1TA3 = 7.268/2.;
  const Float_t bl1TA3 = 0.;
  const Float_t tl1TA3 = 2.03/2.;
  const Float_t alp1TA3 = 7.95; 
  const Float_t h2TA3 = 7.268/2.;
  const Float_t bl2TA3 = 0.;
  const Float_t tl2TA3 = 2.03/2.;
  const Float_t alp2TA3 = 7.95;  
  
// TopEarthFace parameters - single trapezoid
  const Float_t hzTEF = hzTopEarthFaceCu;
  const Float_t tetTEF = 0.;
  const Float_t phiTEF = 0.;
  const Float_t h1TEF = 1.200/2.;
  const Float_t bl1TEF = 21.323/2.;
  const Float_t tl1TEF = 17.963/2.;
  const Float_t alp1TEF = -54.46; 
  const Float_t h2TEF = 1.200/2.;
  const Float_t bl2TEF = 21.323/2.;
  const Float_t tl2TEF = 17.963/2.;
  const Float_t alp2TEF = -54.46;

// TopEarthProfile parameters - single trapezoid
  const Float_t hzTEP = hzTopEarthProfileCu;
  const Float_t tetTEP = 0.;
  const Float_t phiTEP = 0.;
  const Float_t h1TEP = 0.40/2.;
  const Float_t bl1TEP = 31.766/2.;
  const Float_t tl1TEP = 30.535/2.;
  const Float_t alp1TEP = -56.98; 
  const Float_t h2TEP = 0.40/2.;
  const Float_t bl2TEP = 31.766/2.;
  const Float_t tl2TEP = 30.535/2.;
  const Float_t alp2TEP = -56.98;

// TopPositioner parameters - single Stainless Steel trapezoid 
  const Float_t hzTP = hzTopPositionerSteel;
  const Float_t tetTP = 0.;
  const Float_t phiTP = 0.;
  const Float_t h1TP = 3.00/2.;
  const Float_t bl1TP = 7.023/2.;
  const Float_t tl1TP = 7.314/2.;
  const Float_t alp1TP = 2.78; 
  const Float_t h2TP = 3.00/2.;
  const Float_t bl2TP = 7.023/2.;
  const Float_t tl2TP = 7.314/2.;
  const Float_t alp2TP = 2.78;

// TopGasSupport parameters - single cuboid 
  const Float_t hxTGS  = 8.50/2.;
  const Float_t hyTGS  = 3.00/2.;
  const Float_t hzTGS  = hzTopGasSupportAl;
    
// OutEdgeFrame parameters - 4 trapezoidal sections, 2 layers of material
//
//---

// Trapezoid 1
  const Float_t hzOETFE = hzOuterFrameEpoxy;    // layer 1 
  const Float_t hzOETFI = hzOuterFrameInox;     // layer 3
   
  const Float_t tetOETF = 0.;            // common to all 4 trapezoids
  const Float_t phiOETF = 0.;            // common to all 4 trapezoids

  const Float_t h1OETF = 7.196/2.;       // common to all 4 trapezoids
  const Float_t h2OETF = 7.196/2.;       // common to all 4 trapezoids   
  
  const Float_t bl1OETF1 = 3.75/2; 
  const Float_t tl1OETF1 = 3.996/2.;
  const Float_t alp1OETF1 = 0.98;

  const Float_t bl2OETF1 = 3.75/2;
  const Float_t tl2OETF1 = 3.996/2.;
  const Float_t alp2OETF1 = 0.98;
  
// Trapezoid 2
  const Float_t bl1OETF2 = 3.01/2.;
  const Float_t tl1OETF2 = 3.75/2;
  const Float_t alp1OETF2 = 2.94;
      
  const Float_t bl2OETF2 = 3.01/2.;
  const Float_t tl2OETF2 = 3.75/2;
  const Float_t alp2OETF2 = 2.94; 
 
// Trapezoid 3
  const Float_t bl1OETF3 = 1.767/2.;
  const Float_t tl1OETF3 = 3.01/2.;
  const Float_t alp1OETF3 = 4.94;
      
  const Float_t bl2OETF3 = 1.767/2.;
  const Float_t tl2OETF3 = 3.01/2.; 
  const Float_t alp2OETF3 = 4.94; 
  
// Trapezoid 4
  const Float_t bl1OETF4 = 0.;
  const Float_t tl1OETF4 = 1.77/2.;
  const Float_t alp1OETF4 = 7.01;
      
  const Float_t bl2OETF4 = 0.;
  const Float_t tl2OETF4 = 1.77/2.;
  const Float_t alp2OETF4 =  7.01;   
  
// Frame Structure (OutVFrame):
//
// OutVFrame and corner (OutVFrame cuboid, OutVFrame trapezoid)
// EARTHING (VertEarthFaceCu,VertEarthSteel,VertEarthProfCu),
// DETECTOR POSITIONNING (SuppLateralPositionner, LateralPositionner),
// CRADLE (VertCradle), and
// ALIGNMENT (LateralSightSupport, LateralSight) 
//
//---

// OutVFrame parameters - cuboid
  const Float_t hxOutVFrame = 1.85/2.;
  const Float_t hyOutVFrame = 46.23/2.;
  const Float_t hzOutVFrame = hzFrameThickness;

// OutVFrame corner parameters - trapezoid
  const Float_t hzOCTF = hzFrameThickness;
  const Float_t tetOCTF = 0.;
  const Float_t phiOCTF = 0.;
  const Float_t h1OCTF = 1.85/2.;
  const Float_t bl1OCTF = 0.;
  const Float_t tl1OCTF = 3.66/2.;
  const Float_t alp1OCTF = 44.67; 
  const Float_t h2OCTF = 1.85/2.;
  const Float_t bl2OCTF = 0.;
  const Float_t tl2OCTF = 3.66/2.;
  const Float_t alp2OCTF = 44.67;  
  
// VertEarthFaceCu parameters - single trapezoid
  const Float_t hzVFC = hzVertEarthFaceCu;
  const Float_t tetVFC = 0.;
  const Float_t phiVFC = 0.;
  const Float_t h1VFC = 1.200/2.;
  const Float_t bl1VFC = 46.11/2.;
  const Float_t tl1VFC = 48.236/2.;
  const Float_t alp1VFC = 41.54; 
  const Float_t h2VFC = 1.200/2.;
  const Float_t bl2VFC = 46.11/2.;
  const Float_t tl2VFC = 48.236/2.;
  const Float_t alp2VFC = 41.54;
    
// VertEarthSteel parameters - single trapezoid
  const Float_t hzVES = hzVertBarSteel;
  const Float_t tetVES = 0.;
  const Float_t phiVES = 0.;
  const Float_t h1VES = 1.200/2.;
  const Float_t bl1VES = 30.486/2.;
  const Float_t tl1VES = 32.777/2.;
  const Float_t alp1VES = 43.67; 
  const Float_t h2VES = 1.200/2.;
  const Float_t bl2VES = 30.486/2.;
  const Float_t tl2VES = 32.777/2.;
  const Float_t alp2VES = 43.67;

// VertEarthProfCu parameters - single trapezoid
  const Float_t hzVPC = hzVertEarthProfCu;
  const Float_t tetVPC = 0.;
  const Float_t phiVPC = 0.;
  const Float_t h1VPC = 0.400/2.;
  const Float_t bl1VPC = 29.287/2.;
  const Float_t tl1VPC = 30.091/2.;
  const Float_t alp1VPC = 45.14; 
  const Float_t h2VPC = 0.400/2.;
  const Float_t bl2VPC = 29.287/2.;
  const Float_t tl2VPC = 30.091/2.;
  const Float_t alp2VPC = 45.14;

// SuppLateralPositionner - single cuboid
  const Float_t hxSLP  = 2.80/2.;
  const Float_t hySLP  = 5.00/2.;
  const Float_t hzSLP  = hzLateralPosnAl;
  
// LateralPositionner - squared off U bend, face view
  const Float_t hxLPF  = 5.2/2.;
  const Float_t hyLPF  = 3.0/2.;
  const Float_t hzLPF  = hzLateralPosnInoxFace;
  
// LateralPositionner - squared off U bend, profile view
  const Float_t hxLPP  = 0.425/2.;
  const Float_t hyLPP  = 3.0/2.;
  const Float_t hzLPP  = hzLatPosInoxProfM;  // middle layer
  const Float_t hzLPNF  = hzLatPosInoxProfNF; // near and far layers
           
// VertCradle, 3 layers (copies), each composed of 4 trapezoids
// VertCradleA
  const Float_t hzVC1 = hzVerticalCradleAl;
  const Float_t tetVC1 = 0.;
  const Float_t phiVC1 = 0.;
  const Float_t h1VC1 = 10.25/2.;
  const Float_t bl1VC1 = 3.70/2.;
  const Float_t tl1VC1 = 0.;
  const Float_t alp1VC1 = -10.23; 
  const Float_t h2VC1 = 10.25/2.;
  const Float_t bl2VC1 = 3.70/2.;
  const Float_t tl2VC1 = 0.;
  const Float_t alp2VC1 = -10.23;
        
// VertCradleB
  const Float_t hzVC2 = hzVerticalCradleAl;
  const Float_t tetVC2 = 0.;
  const Float_t phiVC2 = 0.;
  const Float_t h1VC2 = 10.25/2.;
  const Float_t bl1VC2 = 6.266/2.;
  const Float_t tl1VC2 = 3.70/2.;
  const Float_t alp1VC2 = -7.13; 
  const Float_t h2VC2 = 10.25/2.;
  const Float_t bl2VC2 = 6.266/2.;
  const Float_t tl2VC2 = 3.70/2.;
  const Float_t alp2VC2 = -7.13;
  
// VertCradleC
  const Float_t hzVC3 = hzVerticalCradleAl;
  const Float_t tetVC3 = 0.;
  const Float_t phiVC3 = 0.;
  const Float_t h1VC3 = 10.25/2.;
  const Float_t bl1VC3 = 7.75/2.;
  const Float_t tl1VC3 = 6.266/2.;
  const Float_t alp1VC3 = -4.14; 
  const Float_t h2VC3 = 10.25/2.;
  const Float_t bl2VC3 = 7.75/2.;
  const Float_t tl2VC3 = 6.266/2.;
  const Float_t alp2VC3 = -4.14;

// VertCradleD
  const Float_t hzVC4 = hzVerticalCradleAl;
  const Float_t tetVC4 = 0.;
  const Float_t phiVC4 = 0.;
  const Float_t h1VC4 = 10.27/2.;
  const Float_t bl1VC4 = 8.273/2.;
  const Float_t tl1VC4 = 7.75/2.;
  const Float_t alp1VC4 = -1.46; 
  const Float_t h2VC4 = 10.27/2.;
  const Float_t bl2VC4 = 8.273/2.;
  const Float_t tl2VC4 = 7.75/2.;
  const Float_t alp2VC4 = -1.46;
  
// LateralSightSupport - single trapezoid
  const Float_t hzVSS = hzLateralSightAl;
  const Float_t tetVSS = 0.;
  const Float_t phiVSS = 0.;
  const Float_t h1VSS = 5.00/2.;
  const Float_t bl1VSS = 7.747/2;
  const Float_t tl1VSS = 7.188/2.;
  const Float_t alp1VSS = -3.20; 
  const Float_t h2VSS = 5.00/2.;
  const Float_t bl2VSS = 7.747/2.;
  const Float_t tl2VSS = 7.188/2.;
  const Float_t alp2VSS = -3.20;  
  
// LateralSight (reference point) - 3 per quadrant, only 1 programmed for now
  const Float_t VSInRad  = 0.6;
  const Float_t VSOutRad  = 1.3;
  const Float_t VSLen  = hzFrameThickness; 
  
//---

// InHFrame parameters
  const Float_t hxInHFrame  = 75.8/2.;
  const Float_t hyInHFrame  = 1.85/2.;
  const Float_t hzInHFrame  = hzFrameThickness;
 
//Flat 7.5mm horizontal section
  const Float_t hxH1mm  = 1.85/2.;
  const Float_t hyH1mm  = 0.75/2.;
  const Float_t hzH1mm  = hzFrameThickness;

//---

// InArcFrame parameters
  const Float_t IAF  = 15.70;
  const Float_t OAF  = 17.55;
  const Float_t hzAF  = hzFrameThickness;
  const Float_t AFphi1  = 0.0;
  const Float_t AFphi2  = 90.0;

//---

// ScrewsInFrame parameters HEAD
  const Float_t SCRUHMI  = 0.;
  const Float_t SCRUHMA  = 0.690/2.;
  const Float_t SCRUHLE  = 0.4/2.;
// ScrewsInFrame parameters MIDDLE
  const Float_t SCRUMMI  = 0.;
  const Float_t SCRUMMA  = 0.39/2.;
  const Float_t SCRUMLE  = hzFrameThickness;
// ScrewsInFrame parameters NUT
  const Float_t SCRUNMI  = 0.;
  const Float_t SCRUNMA  = 0.78/2.;
  const Float_t SCRUNLE  = 0.8/2.;   
  
       // ___________________Make volumes________________________

 Float_t par[11];
 Float_t posX,posY,posZ;

// Quadrant volume TUBS1, positioned at the end
  par[0] = fgkMotherIR1;
  par[1] = fgkMotherOR1; 
  par[2] = fgkMotherThick1;  
  par[3] = fgkMotherPhiL1; 
  par[4] = fgkMotherPhiU1;
  gMC->Gsvolu(QuadrantMLayerName(chamber),"TUBS",idAir,par,5);

// Quadrant volume TUBS2, positioned at the end
  par[0] = fgkMotherIR2;
  par[1] = fgkMotherOR2; 
  par[2] = fgkMotherThick2;  
  par[3] = fgkMotherPhiL2; 
  par[4] = fgkMotherPhiU2;

  gMC->Gsvolu(QuadrantNLayerName(chamber),"TUBS",idAir,par,5); 
  gMC->Gsvolu(QuadrantFLayerName(chamber),"TUBS",idAir,par,5); 

   if (chamber==1) {   
    // InVFrame  
    par[0] = hxInVFrame;
    par[1] = hyInVFrame;
    par[2] = hzInVFrame;
    gMC->Gsvolu("SQ00","BOX",idFrameEpoxy,par,3);

    //Flat 1mm vertical section
    par[0] = hxV1mm;
    par[1] = hyV1mm;
    par[2] = hzV1mm;
    gMC->Gsvolu("SQ01","BOX",idFrameEpoxy,par,3); 
 
// OutTopFrame 
//
// - 3 components (a cuboid and 2 trapezes) and 2 layers (Epoxy/Inox)
//
//---

    // TopFrameAnode - layer 1 of 2 
    par[0] = hxTFA;
    par[1] = hyTFA;
    par[2] = hzTFAE;
    gMC->Gsvolu("SQ02","BOX",idFrameEpoxy,par,3);
    
    // TopFrameAnode - layer 2 of 2 
    par[2] = hzTFAI;
    gMC->Gsvolu("SQ03","BOX",idInox,par,3);
            
    // TopFrameAnodeA - layer 1 of 2  
    par[0] = hzFAAE;
    par[1] = tetFAA;
    par[2] = phiFAA;
    par[3] = h1FAA;
    par[4] = bl1FAA;
    par[5] = tl1FAA;
    par[6] = alp1FAA;
    par[7] = h2FAA;
    par[8] = bl2FAA;
    par[9] = tl2FAA;
    par[10] = alp2FAA;    
    gMC->Gsvolu("SQ04","TRAP",idFrameEpoxy,par,11);    

    // TopFrameAnodeA - layer 2 of 2
    par[0] = hzFAAI;    
    gMC->Gsvolu("SQ05","TRAP",idInox,par,11); 
      
    // TopFrameAnodeB - layer 1 of 2
    par[0] = hzFABE;
    par[1] = tetFAB;
    par[2] = phiFAB;
    par[3] = h1FAB;
    par[4] = bl1FAB;
    par[5] = tl1FAB;
    par[6] = alp1FAB;
    par[7] = h2FAB;
    par[8] = bl2FAB;
    par[9] = tl2FAB;
    par[10] = alp2FAB;
    gMC->Gsvolu("SQ06","TRAP",idFrameEpoxy,par,11);     

    // OutTopTrapFrameB - layer 2 of 2
    par[0] = hzFABI;   
    gMC->Gsvolu("SQ07","TRAP",idInox,par,11);

    // TopAnode1 -  layer 1 of 2
    par[0] = hxTA1;
    par[1] = hyTA1;
    par[2] = hzTA11;    
    gMC->Gsvolu("SQ08","BOX",idInox,par,3); 
    
    // TopAnode1 -  layer 2 of 2
    par[2] = hzTA12;    
    gMC->Gsvolu("SQ09","BOX",idFR4,par,11); 

    // TopAnode2 -  layer 1 of 2
    par[0] = hzTA21;
    par[1] = tetTA2;
    par[2] = phiTA2;
    par[3] = h1TA2;
    par[4] = bl1TA2;
    par[5] = tl1TA2;
    par[6] = alp1TA2;
    par[7] = h2TA2;
    par[8] = bl2TA2;
    par[9] = tl2TA2;
    par[10] = alp2TA2;    
    gMC->Gsvolu("SQ10","TRAP",idInox,par,11); 
 
    // TopAnode2 -  layer 2 of 2
    par[0] = hzTA22;    
    gMC->Gsvolu("SQ11","TRAP",idFR4,par,11);   

    // TopAnode3 -  layer 1 of 1 
    par[0] = hzTA3;
    par[1] = tetTA3;
    par[2] = phiTA3;
    par[3] = h1TA3;
    par[4] = bl1TA3;
    par[5] = tl1TA3;
    par[6] = alp1TA3;
    par[7] = h2TA3;
    par[8] = bl2TA3;
    par[9] = tl2TA3;
    par[10] = alp2TA3;    
    gMC->Gsvolu("SQ12","TRAP",idFR4,par,11); 

    // TopEarthFace 
    par[0] = hzTEF;
    par[1] = tetTEF;
    par[2] = phiTEF;
    par[3] = h1TEF;
    par[4] = bl1TEF;
    par[5] = tl1TEF;
    par[6] = alp1TEF;
    par[7] = h2TEF;
    par[8] = bl2TEF;
    par[9] = tl2TEF;
    par[10] = alp2TEF;    
    gMC->Gsvolu("SQ13","TRAP",idCopper,par,11);   

    // TopEarthProfile 
    par[0] = hzTEP;
    par[1] = tetTEP;
    par[2] = phiTEP;
    par[3] = h1TEP;
    par[4] = bl1TEP;
    par[5] = tl1TEP;
    par[6] = alp1TEP;
    par[7] = h2TEP;
    par[8] = bl2TEP;
    par[9] = tl2TEP;
    par[10] = alp2TEP;
    gMC->Gsvolu("SQ14","TRAP",idCopper,par,11);       

    // TopGasSupport  
    par[0] = hxTGS;
    par[1] = hyTGS;
    par[2] = hzTGS;
    gMC->Gsvolu("SQ15","BOX",idAlu,par,3);

    // TopPositioner parameters - single Stainless Steel trapezoid 
    par[0] = hzTP;
    par[1] = tetTP; 
    par[2] = phiTP;
    par[3] = h1TP;
    par[4] = bl1TP; 
    par[5] = tl1TP; 
    par[6] = alp1TP;
    par[7] = h2TP;
    par[8] = bl2TP; 
    par[9] = tl2TP; 
    par[10] = alp2TP;     
    gMC->Gsvolu("SQ16","TRAP",idInox,par,11);       

//
// OutEdgeTrapFrame Epoxy = (4 trapezes)*2 copies*2 layers (Epoxy/Inox)
//
//---
    // Trapezoid 1 - 2 layers
    par[1] = tetOETF;
    par[2] = phiOETF;
    par[3] = h1OETF;
    par[4] = bl1OETF1;
    par[5] = tl1OETF1;
    par[6] = alp1OETF1;
    par[7] = h2OETF;
    par[8] = bl2OETF1;
    par[9] = tl2OETF1;
    par[10] = alp2OETF1; 
           
    par[0] = hzOETFE;             
    gMC->Gsvolu("SQ17","TRAP",idFrameEpoxy,par,11); 
    par[0] = hzOETFI;
    gMC->Gsvolu("SQ18","TRAP",idInox,par,11);
    
    // Trapezoid 2 - 2 layers
    par[4] = bl1OETF2;
    par[5] = tl1OETF2;
    par[6] = alp1OETF2;

    par[8] = bl2OETF2;
    par[9] = tl2OETF2;
    par[10] = alp2OETF2; 
    
    par[0] = hzOETFE;    
    gMC->Gsvolu("SQ19","TRAP",idFrameEpoxy,par,11);    
    par[0] = hzOETFI;    
    gMC->Gsvolu("SQ20","TRAP",idInox,par,11);     
    
    // Trapezoid 3 - 2 layers
    par[4] = bl1OETF3;
    par[5] = tl1OETF3;
    par[6] = alp1OETF3;

    par[8] = bl2OETF3;
    par[9] = tl2OETF3;
    par[10] = alp2OETF3; 
 
    par[0] = hzOETFE;    
    gMC->Gsvolu("SQ21","TRAP",idFrameEpoxy,par,11);   
    par[0] = hzOETFI;    
    gMC->Gsvolu("SQ22","TRAP",idInox,par,11);     
    
    // Trapezoid 4 - 2 layers

    par[4] = bl1OETF4;
    par[5] = tl1OETF4;
    par[6] = alp1OETF4;

    par[8] = bl2OETF4;
    par[9] = tl2OETF4;
    par[10] = alp2OETF4;  
   
    par[0] = hzOETFE;    
    gMC->Gsvolu("SQ23","TRAP",idFrameEpoxy,par,11);    
    par[0] = hzOETFI;    
    gMC->Gsvolu("SQ24","TRAP",idInox,par,11);     
             
//---
    // OutVFrame    
    par[0] = hxOutVFrame;
    par[1] = hyOutVFrame;
    par[2] = hzOutVFrame;
    gMC->Gsvolu("SQ25","BOX",idFrameEpoxy,par,3);
    
    // OutVFrame corner  
    par[0] = hzOCTF;
    par[1] = tetOCTF;
    par[2] = phiOCTF;
    par[3] = h1OCTF;
    par[4] = bl1OCTF;
    par[5] = tl1OCTF;
    par[6] = alp1OCTF;
    par[7] = h2OCTF;
    par[8] = bl2OCTF;
    par[9] = tl2OCTF;
    par[10] = alp2OCTF;    
    gMC->Gsvolu("SQ26","TRAP",idFrameEpoxy,par,11);
 
    // EarthFaceCu trapezoid
    par[0] = hzVFC;
    par[1] = tetVFC;
    par[2] = phiVFC;
    par[3] = h1VFC;
    par[4] = bl1VFC;
    par[5] = tl1VFC;
    par[6] = alp1VFC;
    par[7] = h2VFC;
    par[8] = bl2VFC;
    par[9] = tl2VFC;
    par[10] = alp2VFC;   
    gMC->Gsvolu("SQ27","TRAP",idCopper,par,11);     

    // VertEarthSteel trapezoid
    par[0] = hzVES;
    par[1] = tetVES;
    par[2] = phiVES;
    par[3] = h1VES;
    par[4] = bl1VES;
    par[5] = tl1VES;
    par[6] = alp1VES;
    par[7] = h2VES;
    par[8] = bl2VES;
    par[9] = tl2VES;
    par[10] = alp2VES;    
    gMC->Gsvolu("SQ28","TRAP",idInox,par,11); 

    // VertEarthProfCu trapezoid       
    par[0] = hzVPC;
    par[1] = tetVPC;
    par[2] = phiVPC;
    par[3] = h1VPC;
    par[4] = bl1VPC;
    par[5] = tl1VPC;
    par[6] = alp1VPC;
    par[7] = h2VPC;
    par[8] = bl2VPC;
    par[9] = tl2VPC;
    par[10] = alp2VPC;
    gMC->Gsvolu("SQ29","TRAP",idCopper,par,11);

    // SuppLateralPositionner cuboid    
    par[0] = hxSLP;
    par[1] = hySLP;
    par[2] = hzSLP;
    gMC->Gsvolu("SQ30","BOX",idAlu,par,3);

    // LateralPositionerFace
    par[0] = hxLPF;
    par[1] = hyLPF;
    par[2] = hzLPF;
    gMC->Gsvolu("SQ31","BOX",idInox,par,3);

    // LateralPositionerProfile
    par[0] = hxLPP;
    par[1] = hyLPP;
    par[2] = hzLPP;
    gMC->Gsvolu("SQ32","BOX",idInox,par,3); // middle layer
    
    par[0] = hxLPP;
    par[1] = hyLPP;
    par[2] = hzLPNF;
    gMC->Gsvolu("SQ33","BOX",idInox,par,3); // near and far layers

    // VertCradleA - 1st trapezoid
    par[0] = hzVC1;
    par[1] = tetVC1;
    par[2] = phiVC1;
    par[3] = h1VC1;
    par[4] = bl1VC1;
    par[5] = tl1VC1;
    par[6] = alp1VC1;
    par[7] = h2VC1;
    par[8] = bl2VC1;
    par[9] = tl2VC1;
    par[10] = alp2VC1;
    gMC->Gsvolu("SQ34","TRAP",idAlu,par,11); 
    
    // VertCradleB - 2nd trapezoid
    par[0] = hzVC2;
    par[1] = tetVC2;
    par[2] = phiVC2;
    par[3] = h1VC2;
    par[4] = bl1VC2;
    par[5] = tl1VC2;
    par[6] = alp1VC2;
    par[7] = h2VC2;
    par[8] = bl2VC2;
    par[9] = tl2VC2;
    par[10] = alp2VC2;
    gMC->Gsvolu("SQ35","TRAP",idAlu,par,11);  
       
    // VertCradleC - 3rd trapezoid
    par[0] = hzVC3;
    par[1] = tetVC3;
    par[2] = phiVC3;
    par[3] = h1VC3;
    par[4] = bl1VC3;
    par[5] = tl1VC3;
    par[6] = alp1VC3;
    par[7] = h2VC3;
    par[8] = bl2VC3;
    par[9] = tl2VC3;
    par[10] = alp2VC3;    
    gMC->Gsvolu("SQ36","TRAP",idAlu,par,11);  

    // VertCradleD - 4th trapezoid
    par[0] = hzVC4;
    par[1] = tetVC4;
    par[2] = phiVC4;
    par[3] = h1VC4;
    par[4] = bl1VC4;
    par[5] = tl1VC4;
    par[6] = alp1VC4;
    par[7] = h2VC4;
    par[8] = bl2VC4;
    par[9] = tl2VC4;
    par[10] = alp2VC4;    
    gMC->Gsvolu("SQ37","TRAP",idAlu,par,11);  
          
    // LateralSightSupport trapezoid
    par[0] = hzVSS;
    par[1] = tetVSS;
    par[2] = phiVSS;
    par[3] = h1VSS;
    par[4] = bl1VSS;
    par[5] = tl1VSS;
    par[6] = alp1VSS;
    par[7] = h2VSS;
    par[8] = bl2VSS;
    par[9] = tl2VSS;
    par[10] = alp2VSS;
    gMC->Gsvolu("SQ38","TRAP",idAlu,par,11);

    // LateralSight
    par[0] = VSInRad;
    par[1] = VSOutRad;
    par[2] = VSLen;       
    gMC->Gsvolu("SQ39","TUBE",idFrameEpoxy,par,3);   

//---
    // InHFrame
    par[0] = hxInHFrame;
    par[1] = hyInHFrame;
    par[2] = hzInHFrame;
    gMC->Gsvolu("SQ40","BOX",idFrameEpoxy,par,3);

    //Flat 7.5mm horizontal section
    par[0] = hxH1mm;
    par[1] = hyH1mm;
    par[2] = hzH1mm;
    gMC->Gsvolu("SQ41","BOX",idFrameEpoxy,par,3);

    // InArcFrame 
    par[0] = IAF;
    par[1] = OAF; 
    par[2] = hzAF;  
    par[3] = AFphi1; 
    par[4] = AFphi2;

    gMC->Gsvolu("SQ42","TUBS",idFrameEpoxy,par,5);

//---
    // ScrewsInFrame - 3 sections in order to avoid overlapping volumes
    // Screw Head, in air
    par[0] = SCRUHMI;
    par[1] = SCRUHMA; 
    par[2] = SCRUHLE;  

    gMC->Gsvolu("SQ43","TUBE",idInox,par,3);
    
    // Middle part, in the Epoxy
    par[0] = SCRUMMI;
    par[1] = SCRUMMA;
    par[2] = SCRUMLE;
    gMC->Gsvolu("SQ44","TUBE",idInox,par,3);
    
    // Screw nut, in air
    par[0] = SCRUNMI;
    par[1] = SCRUNMA;
    par[2] = SCRUNLE;   
    gMC->Gsvolu("SQ45","TUBE",idInox,par,3);     
   }
              
// __________________Place volumes in the quadrant ____________ 
        
    // InVFrame  
    posX = hxInVFrame;
    posY = 2.0*hyInHFrame+2.*hyH1mm+IAF+hyInVFrame;        
    posZ = 0.;
    gMC->Gspos("SQ00",1,QuadrantMLayerName(chamber),posX, posY, posZ, 0, "ONLY"); 

    //Flat 7.5mm vertical section
    posX = 2.0*hxInVFrame+hxV1mm;
    posY = 2.0*hyInHFrame+2.*hyH1mm+IAF+hyV1mm;
    posZ = 0.;
    gMC->Gspos("SQ01",1,QuadrantMLayerName(chamber),posX, posY, posZ,0, "ONLY"); 
    
    // TopFrameAnode place 2 layers of TopFrameAnode cuboids  
    posX = hxTFA;
    posY = 2.*hyInHFrame+2.*hyH1mm+IAF+2.*hyInVFrame+hyTFA;   
    posZ = hzOuterFrameInox;
    gMC->Gspos("SQ02",1,QuadrantMLayerName(chamber),posX, posY, posZ,0,"ONLY"); 
    posZ = posZ+hzOuterFrameInox;
    gMC->Gspos("SQ03",1,QuadrantMLayerName(chamber),posX, posY, posZ,0,"ONLY");
    
    // place 2 layers of TopFrameAnodeA trapezoids 
    posX = 35.8932+fgkDeltaQuadLHC;
    posY = 92.6745+fgkDeltaQuadLHC;
    posZ = hzOuterFrameInox; 
    gMC->Gspos("SQ04",1,QuadrantMLayerName(chamber),posX, posY, posZ, rot1,"ONLY");
    posZ = posZ+hzOuterFrameInox;
    gMC->Gspos("SQ05",1,QuadrantMLayerName(chamber),posX, posY, posZ, rot1,"ONLY");
    
    // place 2 layers of TopFrameAnodeB trapezoids 
    posX = 44.593+fgkDeltaQuadLHC;
    posY = 90.737+fgkDeltaQuadLHC;
    posZ = hzOuterFrameInox; 
    gMC->Gspos("SQ06",1,QuadrantMLayerName(chamber),posX, posY, posZ, rot1,"ONLY");
    posZ = posZ+hzOuterFrameInox;
    gMC->Gspos("SQ07",1,QuadrantMLayerName(chamber),posX, posY, posZ, rot1,"ONLY");    

    // TopAnode1 place 2 layers  
    posX = 6.8+fgkDeltaQuadLHC;
    posY = 99.85+fgkDeltaQuadLHC;
    posZ = -1.*hzAnodeFR4;
    gMC->Gspos("SQ08",1,QuadrantMLayerName(chamber),posX, posY, posZ, 0,"ONLY");  
    posZ = posZ+hzTopAnodeSteel1;
    gMC->Gspos("SQ09",1,QuadrantMLayerName(chamber),posX, posY, posZ, 0,"ONLY");    
         
    // TopAnode2 place 2 layers
    posX = 18.534+fgkDeltaQuadLHC;
    posY = 99.482+fgkDeltaQuadLHC; 
    posZ = -1.*hzAnodeFR4;    
    gMC->Gspos("SQ10",1,QuadrantMLayerName(chamber),posX, posY, posZ, rot1,"ONLY");
    posZ = posZ+hzTopAnodeSteel2;    
    gMC->Gspos("SQ11",1,QuadrantMLayerName(chamber),posX, posY, posZ, rot1,"ONLY");       
    
    // TopAnode3 place 1 layer
    posX = 25.80+fgkDeltaQuadLHC;
    posY = 98.61+fgkDeltaQuadLHC;
    posZ = 0.;    
    gMC->Gspos("SQ12",1,QuadrantMLayerName(chamber),posX, posY, posZ, rot1,"ONLY");  
          
    // TopEarthFace - 2 copies
    posX = 23.122+fgkDeltaQuadLHC;
    posY = 96.90+fgkDeltaQuadLHC;
    posZ = hzOuterFrameEpoxy+hzOuterFrameInox+hzTopEarthFaceCu;
    gMC->Gspos("SQ13",1,QuadrantMLayerName(chamber),posX, posY, posZ, 0,"ONLY");
    posZ = -1.*posZ;
    gMC->Gspos("SQ13",2,QuadrantMLayerName(chamber),posX, posY, posZ, 0,"ONLY");

    // TopEarthProfile 
    posX = 14.475+fgkDeltaQuadLHC;
    posY = 97.900+fgkDeltaQuadLHC; 
    posZ = hzTopEarthProfileCu;
    gMC->Gspos("SQ14",1,QuadrantMLayerName(chamber),posX, posY, posZ, 0,"ONLY");
    posZ = -1.0*posZ;
    gMC->Gspos("SQ14",2,QuadrantMLayerName(chamber),posX, posY, posZ, 0,"ONLY");

    // TopGasSupport - 2 copies                            
    posX = 4.9500+fgkDeltaQuadLHC;
    posY = 96.200+fgkDeltaQuadLHC;
    posZ = hzOuterFrameEpoxy+hzOuterFrameInox+hzTopGasSupportAl;
    gMC->Gspos("SQ15",1,QuadrantMLayerName(chamber),posX, posY, posZ, 0,"ONLY");
    posZ = -1.*posZ;
    gMC->Gspos("SQ15",2,QuadrantMLayerName(chamber),posX, posY, posZ, 0,"ONLY");
    
    // TopPositioner parameters - single Stainless Steel trapezoid - 2 copies
    posX = 7.60+fgkDeltaQuadLHC;
    posY = 98.98+fgkDeltaQuadLHC;   
    posZ = hzOuterFrameEpoxy+hzOuterFrameInox+2.*hzTopGasSupportAl+hzTopPositionerSteel;
    gMC->Gspos("SQ16",1,QuadrantMLayerName(chamber),posX, posY, posZ, rot1,"ONLY");
    posZ = -1.*posZ;
    gMC->Gspos("SQ16",2,QuadrantMLayerName(chamber),posX, posY, posZ, rot1,"ONLY"); 

    // OutEdgeFrame 
    Float_t XCenter[8]; 
    Float_t YCenter[8];
    
    XCenter[0] = 73.201 + fgkDeltaQuadLHC;
    XCenter[1] = 78.124 + fgkDeltaQuadLHC; 
    XCenter[2] = 82.862 + fgkDeltaQuadLHC;
    XCenter[3] = 87.418 + fgkDeltaQuadLHC; 
    
    YCenter[0] = 68.122 + fgkDeltaQuadLHC;
    YCenter[1] = 62.860 + fgkDeltaQuadLHC;   
    YCenter[2] = 57.420 + fgkDeltaQuadLHC;
    YCenter[3] = 51.800 + fgkDeltaQuadLHC; 
      
    XCenter[4] = 68.122 + fgkDeltaQuadLHC;
    XCenter[5] = 62.860 + fgkDeltaQuadLHC; 
    XCenter[6] = 57.420 + fgkDeltaQuadLHC;
    XCenter[7] = 51.800 + fgkDeltaQuadLHC; 
    
    YCenter[4] = 73.210 + fgkDeltaQuadLHC;
    YCenter[5] = 78.124 + fgkDeltaQuadLHC; 
    YCenter[6] = 82.862 + fgkDeltaQuadLHC;
    YCenter[7] = 87.418 + fgkDeltaQuadLHC; 
      
    posZ = -1.0*hzOuterFrameInox;     
    gMC->Gspos("SQ17",1,QuadrantMLayerName(chamber), XCenter[0], YCenter[0], posZ, rot2,"ONLY");
    gMC->Gspos("SQ17",2,QuadrantMLayerName(chamber), XCenter[4], YCenter[4], posZ, rot3,"ONLY");

    gMC->Gspos("SQ19",1,QuadrantMLayerName(chamber), XCenter[1], YCenter[1], posZ, rot2,"ONLY");   
    gMC->Gspos("SQ19",2,QuadrantMLayerName(chamber), XCenter[5], YCenter[5], posZ, rot3,"ONLY");

    gMC->Gspos("SQ21",1,QuadrantMLayerName(chamber), XCenter[2], YCenter[2], posZ, rot2,"ONLY");
    gMC->Gspos("SQ21",2,QuadrantMLayerName(chamber), XCenter[6], YCenter[6], posZ, rot3,"ONLY");
    
    gMC->Gspos("SQ23",1,QuadrantMLayerName(chamber), XCenter[3], YCenter[3], posZ, rot2,"ONLY");
    gMC->Gspos("SQ23",2,QuadrantMLayerName(chamber), XCenter[7], YCenter[7], posZ, rot3,"ONLY");
     
    posZ = posZ+hzOuterFrameEpoxy;
   
    gMC->Gspos("SQ18",1,QuadrantMLayerName(chamber), XCenter[0], YCenter[0], posZ, rot2,"ONLY");
    gMC->Gspos("SQ18",2,QuadrantMLayerName(chamber), XCenter[4], YCenter[4], posZ, rot3,"ONLY");
    
    gMC->Gspos("SQ20",1,QuadrantMLayerName(chamber), XCenter[1], YCenter[1], posZ, rot2,"ONLY");   
    gMC->Gspos("SQ20",2,QuadrantMLayerName(chamber), XCenter[5], YCenter[5], posZ, rot3,"ONLY");

    gMC->Gspos("SQ22",1,QuadrantMLayerName(chamber), XCenter[2], YCenter[2], posZ, rot2,"ONLY");
    gMC->Gspos("SQ22",2,QuadrantMLayerName(chamber), XCenter[6], YCenter[6], posZ, rot3,"ONLY");
       
    gMC->Gspos("SQ24",1,QuadrantMLayerName(chamber), XCenter[3], YCenter[3], posZ, rot2,"ONLY");
    gMC->Gspos("SQ24",2,QuadrantMLayerName(chamber), XCenter[7], YCenter[7], posZ, rot3,"ONLY");  

//---    
        
// OutVFrame
    posX = 2.*hxInVFrame+IAF+2.*hxInHFrame-hxOutVFrame+2.*hxV1mm;
    posY = 2.*hyInHFrame+hyOutVFrame;    
    posZ = 0.;              
    gMC->Gspos("SQ25",1,QuadrantMLayerName(chamber),posX, posY, posZ, 0, "ONLY"); 

    const Float_t TOPY = posY+hyOutVFrame;
    const Float_t OUTX = posX;

// OutVFrame corner
    posX = OUTX;
    posY = TOPY+((bl1OCTF+tl1OCTF)/2.);
    posZ = 0.;     
    gMC->Gspos("SQ26",1,QuadrantMLayerName(chamber),posX, posY, posZ, rot1,"ONLY"); 

// VertEarthFaceCu - 2 copies
    posX = 89.4000+fgkDeltaQuadLHC;
    posY = 25.79+fgkDeltaQuadLHC;    
    posZ = hzFrameThickness+2.0*hzFoam+hzVertEarthFaceCu;              
    gMC->Gspos("SQ27",1,QuadrantMLayerName(chamber),posX, posY, posZ, rot1, "ONLY"); 
    posZ = -1.0*posZ; 
    gMC->Gspos("SQ27",2,QuadrantMLayerName(chamber),posX, posY, posZ, rot1, "ONLY"); 
    
// VertEarthSteel - 2 copies
    posX = 91.00+fgkDeltaQuadLHC;
    posY = 30.616+fgkDeltaQuadLHC;    
    posZ = hzFrameThickness+2.0*hzFoam+hzVertBarSteel;              
    gMC->Gspos("SQ28",1,QuadrantMLayerName(chamber),posX, posY, posZ, rot1, "ONLY"); 
    posZ = -1.0*posZ;              
    gMC->Gspos("SQ28",2,QuadrantMLayerName(chamber),posX, posY, posZ, rot1, "ONLY");
 
// VertEarthProfCu - 2 copies
    posX = 92.000+fgkDeltaQuadLHC;
    posY = 29.64+fgkDeltaQuadLHC;    
    posZ = hzFrameThickness;              
    gMC->Gspos("SQ29",1,QuadrantMLayerName(chamber),posX, posY, posZ, rot1, "ONLY"); 
    posZ = -1.0*posZ;    
    gMC->Gspos("SQ29",2,QuadrantMLayerName(chamber),posX, posY, posZ, rot1, "ONLY"); 

// SuppLateralPositionner - 2 copies 
    posX = 90.2-fgkNearFarLHC;
    posY = 5.00-fgkNearFarLHC;    
    posZ = hzLateralPosnAl-fgkMotherThick2;             
    gMC->Gspos("SQ30",1,QuadrantFLayerName(chamber),posX, posY, posZ, 0, "ONLY"); 
    posZ = -1.0*posZ;            
    gMC->Gspos("SQ30",2,QuadrantNLayerName(chamber),posX, posY, posZ, 0, "ONLY"); 

// LateralPositionner - 2 copies - Face view
    posX = 92.175-fgkNearFarLHC-2.*hxLPP;
    posY = 5.00-fgkNearFarLHC;   
    posZ =2.0*hzLateralPosnAl+hzLateralPosnInoxFace-fgkMotherThick2;              
    gMC->Gspos("SQ31",1,QuadrantFLayerName(chamber),posX, posY, posZ, 0, "ONLY"); 
    posZ = -1.0*posZ;             
    gMC->Gspos("SQ31",2,QuadrantNLayerName(chamber),posX, posY, posZ, 0, "ONLY"); 

// LateralPositionner -  Profile view   
    posX = 92.175+fgkDeltaQuadLHC+hxLPF-hxLPP;
    posY = 5.00+fgkDeltaQuadLHC;    
    posZ = 0.;              
    gMC->Gspos("SQ32",1,QuadrantMLayerName(chamber),posX, posY, posZ, 0, "ONLY"); // middle layer

    posX = 92.175-fgkNearFarLHC+hxLPF-hxLPP; 
    posY = 5.0000-fgkNearFarLHC;    
    posZ = fgkMotherThick2-hzLPNF;              
    gMC->Gspos("SQ33",1,QuadrantNLayerName(chamber),posX, posY, posZ, 0, "ONLY"); // near layer
    posZ = -1.*posZ;
    gMC->Gspos("SQ33",2,QuadrantFLayerName(chamber),posX, posY, posZ, 0, "ONLY"); // far layer
      
// VertCradleA  1st Trapezoid - 3 copies
    posX = 95.73+fgkDeltaQuadLHC;
    posY = 33.26+fgkDeltaQuadLHC; 
    posZ = 0.;              
    gMC->Gspos("SQ34",2,QuadrantMLayerName(chamber),posX, posY, posZ, 0, "ONLY");  

    posX = 95.73-fgkNearFarLHC;
    posY = 33.26-fgkNearFarLHC;
    posZ = 2.0*hzLateralSightAl+hzVerticalCradleAl-fgkMotherThick2;               
    gMC->Gspos("SQ34",1,QuadrantNLayerName(chamber),posX, posY, posZ, 0, "ONLY");
    posZ = -1.0*posZ;              
    gMC->Gspos("SQ34",3,QuadrantFLayerName(chamber),posX, posY, posZ, 0, "ONLY");

// VertCradleB  2nd Trapezoid - 3 copies
    posX = 97.29+fgkDeltaQuadLHC;
    posY = 23.02+fgkDeltaQuadLHC;    
    posZ = 0.;              
    gMC->Gspos("SQ35",2,QuadrantMLayerName(chamber),posX, posY, posZ, 0, "ONLY");

    posX = 97.29-fgkNearFarLHC;
    posY = 23.02-fgkNearFarLHC;   
    posZ = 2.0*hzLateralSightAl+hzVerticalCradleAl-fgkMotherThick2;          
    gMC->Gspos("SQ35",1,QuadrantNLayerName(chamber),posX, posY, posZ, 0, "ONLY");    
    posZ = -1.0*posZ;          
    gMC->Gspos("SQ35",3,QuadrantFLayerName(chamber),posX, posY, posZ, 0, "ONLY");

// OutVertCradleC  3rd Trapeze - 3 copies
    posX = 98.31+fgkDeltaQuadLHC;
    posY = 12.77+fgkDeltaQuadLHC;  
    posZ = 0.;              
    gMC->Gspos("SQ36",2,QuadrantMLayerName(chamber),posX, posY, posZ, 0, "ONLY");

    posX = 98.31-fgkNearFarLHC;
    posY = 12.77-fgkNearFarLHC;        

    posZ = 2.0*hzLateralSightAl+hzVerticalCradleAl-fgkMotherThick2;         
    gMC->Gspos("SQ36",1,QuadrantNLayerName(chamber),posX, posY, posZ, 0, "ONLY");       
    posZ = -1.0*posZ;
    gMC->Gspos("SQ36",3,QuadrantFLayerName(chamber),posX, posY, posZ, 0, "ONLY");  

// OutVertCradleD  4th Trapeze - 3 copies
    posX = 98.81+fgkDeltaQuadLHC;
    posY = 2.52+fgkDeltaQuadLHC;    
    posZ = 0.;              
    gMC->Gspos("SQ37",2,QuadrantMLayerName(chamber),posX, posY, posZ, 0, "ONLY");
   
    posZ = fgkMotherThick1-hzVerticalCradleAl;                
    gMC->Gspos("SQ37",1,QuadrantMLayerName(chamber),posX, posY, posZ, 0, "ONLY");
    posZ = -1.0*posZ;          
    gMC->Gspos("SQ37",3,QuadrantMLayerName(chamber),posX, posY, posZ, 0, "ONLY");          
             
// LateralSightSupport - 2 copies
    posX = 98.53-fgkNearFarLHC;
    posY = 10.00-fgkNearFarLHC;    
    posZ = hzLateralSightAl-fgkMotherThick2;
    gMC->Gspos("SQ38",1,QuadrantNLayerName(chamber),posX, posY, posZ, 0, "ONLY"); 
    posZ = -1.0*posZ;             
    gMC->Gspos("SQ38",2,QuadrantFLayerName(chamber),posX, posY, posZ, 0, "ONLY"); 
    
// Mire placement
    posX = 92.84+fgkDeltaQuadLHC;  
    posY = 8.13+fgkDeltaQuadLHC;
    posZ = 0.;
    gMC->Gspos("SQ39",1,QuadrantMLayerName(chamber),posX, posY, posZ, 0,"ONLY");    

//---

// InHFrame
    posX = 2.0*hxInVFrame+2.*hxV1mm+IAF+hxInHFrame;
    posY = hyInHFrame;
    posZ = 0.;       
    gMC->Gspos("SQ40",1,QuadrantMLayerName(chamber),posX, posY, posZ, 0, "ONLY"); 
 
// Flat 7.5mm horizontal section
    posX = 2.0*hxInVFrame+2.*hxV1mm+IAF+hxH1mm;
    posY = 2.0*hyInHFrame+hyH1mm;
    posZ = 0.;
    gMC->Gspos("SQ41",1,QuadrantMLayerName(chamber),posX, posY, posZ,0, "ONLY"); 
        
// InArcFrame 
    posX = 2.0*hxInVFrame+2.*hxV1mm;
    posY = 2.0*hyInHFrame+2.*hyH1mm;
    posZ = 0.;    
    gMC->Gspos("SQ42",1,QuadrantMLayerName(chamber),posX, posY, posZ,0, "ONLY"); 

// ScrewsInFrame - in sensitive volume

     Float_t scruX[64];
     Float_t scruY[64]; 
         
// Screws on IHEpoxyFrame

     const Int_t NumberOfScrewsIH = 14;    // no. of screws on the IHEpoxyFrame
     const Float_t offX = 5.;              // inter-screw distance 

     // first screw coordinates 
     scruX[0] = 21.07;                  
     scruY[0] = -2.23; 
     // other screw coordinates      
     for (Int_t i = 1;i<NumberOfScrewsIH;i++){   
     scruX[i] = scruX[i-1]+offX; 
     scruY[i] = scruY[0];
     }    
     // Position the volumes on the frames
     for (Int_t i = 0;i<NumberOfScrewsIH;i++){
     posX = fgkDeltaQuadLHC + scruX[i];
     posY = fgkDeltaQuadLHC + scruY[i];
     posZ = 0.;   
     gMC->Gspos("SQ43",i+1,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ-hzInHFrame-SCRUHLE, 0, "ONLY");      
     gMC->Gspos("SQ44",i+1,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ, 0, "ONLY");
     gMC->Gspos("SQ45",i+1,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ+hzInHFrame+SCRUNLE, 0, "ONLY"); 
     }
     // special screw coordinates
     scruX[63] = 16.3;  
     scruY[63] = -2.23; 
     posX = fgkDeltaQuadLHC + scruX[63];
     posY = fgkDeltaQuadLHC + scruY[63];
     posZ = 0.;            
     gMC->Gspos("SQ43",64,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ-hzInHFrame-SCRUHLE, 0, "ONLY");
     gMC->Gspos("SQ44",64,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ, 0, "ONLY"); 
     gMC->Gspos("SQ45",64,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ+hzInHFrame+SCRUNLE, 0, "ONLY");  
     
// Screws on the IVEpoxyFrame
  
    const Int_t NumberOfScrewsIV = 15;     // no. of screws on the IVEpoxyFrame
    const Float_t offY = 5.;               // inter-screw distance 
    Int_t FirstScrew = 58;
    Int_t LastScrew = 44;
 
    // first (special) screw coordinates
    scruX[FirstScrew-1] = -2.23; 
    scruY[FirstScrew-1] = 16.3; 
    // second (repetitive) screw coordinates
    scruX[FirstScrew-2] = -2.23; 
    scruY[FirstScrew-2] = 21.07;     
    // other screw coordinates      
    for (Int_t i = FirstScrew-3;i>LastScrew-2;i--){   
    scruX[i] = scruX[FirstScrew-2];
    scruY[i] = scruY[i+1]+offY;
    }
    
    for (Int_t i = 0;i<NumberOfScrewsIV;i++){
    posX = fgkDeltaQuadLHC + scruX[i+LastScrew-1];
    posY = fgkDeltaQuadLHC + scruY[i+LastScrew-1];
    posZ = 0.;       
    gMC->Gspos("SQ43",i+LastScrew,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ-hzInHFrame-SCRUHLE, 0, "ONLY");     
    gMC->Gspos("SQ44",i+LastScrew,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ, 0, "ONLY"); 
    gMC->Gspos("SQ45",i+LastScrew,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ+hzInHFrame+SCRUNLE, 0, "ONLY");
    }    
    
// Screws on the OVEpoxyFrame
  
    const Int_t NumberOfScrewsOV = 10;     // no. of screws on the OVEpoxyFrame

    FirstScrew = 15;
    LastScrew = 25;
 
    // first (repetitive) screw coordinates
    scruX[FirstScrew-1] = 90.9; 
    scruY[FirstScrew-1] = -2.23;  // true value
 
    // other screw coordinates      
    for (Int_t i = FirstScrew; i<LastScrew; i++ ){   
    scruX[i] = scruX[FirstScrew-1];
    scruY[i] = scruY[i-1]+offY;
    }
    for (Int_t i = 0;i<NumberOfScrewsOV;i++){
    posX = fgkDeltaQuadLHC + scruX[i+FirstScrew-1];
    posY = fgkDeltaQuadLHC + scruY[i+FirstScrew-1];
    posZ = 0.;   
    gMC->Gspos("SQ43",i+FirstScrew,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ-hzInHFrame-SCRUHLE, 0, "ONLY");     
    gMC->Gspos("SQ44",i+FirstScrew,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ, 0, "ONLY"); 
    gMC->Gspos("SQ45",i+FirstScrew,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ+hzInHFrame+SCRUNLE, 0, "ONLY"); 
    }
      
// Inner Arc of Frame, screw positions and numbers-1
   scruX[62] = 16.009; scruY[62]  = 1.401;
   scruX[61] = 14.564; scruY[61]  = 6.791;
   scruX[60] = 11.363; scruY[60]  = 11.363;
   scruX[59] = 6.791 ; scruY[59]  = 14.564;
   scruX[58] = 1.401 ; scruY[58]  = 16.009;
    
    for (Int_t i = 0;i<5;i++){
    posX = fgkDeltaQuadLHC + scruX[i+58];
    posY = fgkDeltaQuadLHC + scruY[i+58];
    posZ = 0.;   
    gMC->Gspos("SQ43",i+58+1,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ-hzInHFrame-SCRUHLE, 0, "ONLY");    
    gMC->Gspos("SQ44",i+58+1,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ, 0, "ONLY");
    gMC->Gspos("SQ45",i+58+1,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ+hzInHFrame+SCRUNLE, 0, "ONLY");
    }
}

//______________________________________________________________________________
void AliMUONv2::PlaceInnerLayers(Int_t chamber)
{
// Place the gas and copper layers for the specified chamber.
// --

// Rotation Matrices 
  Int_t rot1, rot2, rot3, rot4;   

  AliMatrix(rot1,  90., 315., 90.,  45., 0., 0.); // -45 deg
  AliMatrix(rot2,  90.,  90., 90., 180., 0., 0.); //  90 deg
  AliMatrix(rot3,  90., 270., 90.,   0., 0., 0.); // -90 deg 
  AliMatrix(rot4,  90.,  45., 90., 135., 0., 0.); //  deg 

  GReal_t x;
  GReal_t y;
  GReal_t zg = 0.;
  GReal_t zc = fgkHzGas + fgkHzPadPlane;
  Int_t dpos = (chamber-1)*2;
  TString name;
  
  x = 14.53 + fgkDeltaQuadLHC;
  y = 53.34 + fgkDeltaQuadLHC;
  name = GasVolumeName("SAG", chamber);
  gMC->Gspos(name,1,QuadrantMLayerName(chamber),x,y,zg,0,"ONLY");
  gMC->Gspos("SA1C", 1+dpos, QuadrantMLayerName(chamber),x,y, zc,0,"ONLY");
  gMC->Gspos("SA1C", 2+dpos, QuadrantMLayerName(chamber),x,y,-zc,0,"ONLY");

  x = 40.67 + fgkDeltaQuadLHC;
  y = 40.66 + fgkDeltaQuadLHC;    
  name = GasVolumeName("SBG", chamber);
  gMC->Gspos(name,1,QuadrantMLayerName(chamber),x,y,zg,rot1,"ONLY"); 
  gMC->Gspos("SB1C", 1+dpos ,QuadrantMLayerName(chamber),x,y, zc,rot1,"ONLY");
  gMC->Gspos("SB1C", 2+dpos, QuadrantMLayerName(chamber),x,y,-zc,rot1,"ONLY");

  x = 53.34 + fgkDeltaQuadLHC;
  y = 14.52 + fgkDeltaQuadLHC; 
  name = GasVolumeName("SCG", chamber);
  gMC->Gspos(name,1,QuadrantMLayerName(chamber),x,y,zg,rot2,"ONLY");
  gMC->Gspos("SC1C", 1+dpos ,QuadrantMLayerName(chamber),x,y, zc,rot2,"ONLY");
  gMC->Gspos("SC1C", 2+dpos ,QuadrantMLayerName(chamber),x,y,-zc,rot2,"ONLY");

  x = 5.83 + fgkDeltaQuadLHC;
  y = 17.29 + fgkDeltaQuadLHC;
  name = GasVolumeName("SDG", chamber);
  gMC->Gspos(name,1,QuadrantMLayerName(chamber),x,y,zg,rot3,"ONLY");
  gMC->Gspos("SD1C", 1+dpos ,QuadrantMLayerName(chamber),x,y, zc,rot3,"ONLY");
  gMC->Gspos("SD1C", 2+dpos ,QuadrantMLayerName(chamber),x,y,-zc,rot3,"ONLY");

  x = 9.04 + fgkDeltaQuadLHC;
  y = 16.91 + fgkDeltaQuadLHC; 
  name = GasVolumeName("SEG", chamber);
  gMC->Gspos(name,1,QuadrantMLayerName(chamber),x,y,zg,0,"ONLY");
  gMC->Gspos("SE1C", 1+dpos ,QuadrantMLayerName(chamber),x,y, zc,0,"ONLY");
  gMC->Gspos("SE1C", 2+dpos ,QuadrantMLayerName(chamber),x,y,-zc,0,"ONLY");

  x = 10.12 + fgkDeltaQuadLHC;
  y = 14.67 + fgkDeltaQuadLHC;  
  name = GasVolumeName("SFG", chamber);
  gMC->Gspos(name,1,QuadrantMLayerName(chamber),x,y,zg,rot4,"ONLY");   
  gMC->Gspos("SF1C", 1+dpos ,QuadrantMLayerName(chamber),x,y, zc,rot4,"ONLY");
  gMC->Gspos("SF1C", 2+dpos ,QuadrantMLayerName(chamber),x,y,-zc,rot4,"ONLY");

  x = 8.2042 + fgkDeltaQuadLHC;
  y = 16.19 + fgkDeltaQuadLHC;
  name = GasVolumeName("SGG", chamber);
  gMC->Gspos(name,1,QuadrantMLayerName(chamber),x,y,zg,rot4,"ONLY");
  gMC->Gspos("SG1C", 1+dpos ,QuadrantMLayerName(chamber),x,y, zc,rot4,"ONLY");
  gMC->Gspos("SG1C", 2+dpos ,QuadrantMLayerName(chamber),x,y,-zc,rot4,"ONLY");

  x = 14.68 + fgkDeltaQuadLHC;
  y = 10.10 + fgkDeltaQuadLHC;
  name = GasVolumeName("SHG", chamber);
  gMC->Gspos(name,1,QuadrantMLayerName(chamber),x,y,zg,rot4,"ONLY");
  gMC->Gspos("SH1C", 1+dpos ,QuadrantMLayerName(chamber),x,y, zc,rot4,"ONLY");
  gMC->Gspos("SH1C", 2+dpos ,QuadrantMLayerName(chamber),x,y,-zc,rot4,"ONLY");

  x = 16.21 + fgkDeltaQuadLHC;
  y = 8.17 + fgkDeltaQuadLHC;
  name = GasVolumeName("SIG", chamber);
  gMC->Gspos(name,1,QuadrantMLayerName(chamber),x,y,zg,rot4,"ONLY");
  gMC->Gspos("SI1C", 1+dpos ,QuadrantMLayerName(chamber),x,y, zc,rot4,"ONLY");
  gMC->Gspos("SI1C", 2+dpos ,QuadrantMLayerName(chamber),x,y,-zc,rot4,"ONLY");

  x = 16.92 + fgkDeltaQuadLHC;
  y = 9.02 + fgkDeltaQuadLHC;
  name = GasVolumeName("SJG", chamber);
  gMC->Gspos(name,1,QuadrantMLayerName(chamber),x,y,zg,rot3,"ONLY");
  gMC->Gspos("SJ1C", 1+dpos ,QuadrantMLayerName(chamber),x,y, zc,rot3,"ONLY");
  gMC->Gspos("SJ1C", 2+dpos ,QuadrantMLayerName(chamber),x,y,-zc,rot3,"ONLY");

  x =  17.30 + fgkDeltaQuadLHC;
  y =  5.85 + fgkDeltaQuadLHC;
  name = GasVolumeName("SKG", chamber);
  gMC->Gspos(name,1,QuadrantMLayerName(chamber),x,y,zg,0,"ONLY");
  gMC->Gspos("SK1C", 1+dpos ,QuadrantMLayerName(chamber),x,y, zc,0,"ONLY");
  gMC->Gspos("SK1C", 2+dpos ,QuadrantMLayerName(chamber),x,y,-zc,0,"ONLY");
}

//______________________________________________________________________________
void AliMUONv2::PlaceSector(AliMpSector* sector,TSpecialMap specialMap, 
                            const TVector3& where, Bool_t reflectZ, Int_t chamber)
{
// Place all the segments in the mother volume, at the position defined
// by the sector's data.
// --

  static Int_t segNum=1;
  Int_t sgn;
  Int_t reflZ;
  Int_t rotMat;

  if (!reflectZ) {
    sgn= 1;
    reflZ=0;                                     // no reflection along z... nothing
    AliMatrix(rotMat,  90.,90.,90,180.,0.,0.);   // 90° rotation around z, NO reflection along z
  } else  {
    sgn=-1;
    AliMatrix(reflZ,  90.,0.,90,90.,180.,0.);    // reflection along z
    AliMatrix(rotMat,  90.,90.,90,180.,180.,0.); // 90° rotation around z AND reflection along z
  }
  
  GReal_t posX,posY,posZ;
  
  IntVector already_done;
  for (Int_t irow=0;irow<sector->GetNofRows();irow++){ // for each row
    AliMpRow* row = sector->GetRow(irow);


    for (Int_t iseg=0;iseg<row->GetNofRowSegments();iseg++){ // for each row segment
      AliMpVRowSegment* seg = row->GetRowSegment(iseg);
      char segName[5];
      
      TSpecialMap::iterator iter 
        = specialMap.find(seg->GetMotifPositionId(0));

      if ( iter == specialMap.end()){ //if this is a normal segment (ie. not part of <specialMap>)
      
        // create the cathode part
        sprintf(segName,"%.3dM", segNum);
        CreatePlaneSegment(segName, seg->Dimensions()/10., seg->GetNofMotifs());
  
        posX = where.X() + seg->Position().X()/10.;
        posY = where.Y() + seg->Position().Y()/10.;
        posZ = where.Z() + sgn * (TotalHzPlane() + fgkHzGas + 2.*fgkHzPadPlane);
        gMC->Gspos(segName, 1, QuadrantMLayerName(chamber), posX, posY, posZ, reflZ, "ONLY");

        // and place all the daughter boards of this segment
        for (Int_t motifNum=0;motifNum<seg->GetNofMotifs();motifNum++) {
          Int_t motifPosId = seg->GetMotifPositionId(motifNum);
          AliMpMotifPosition* motifPos = 
            sector->GetMotifMap()->FindMotifPosition(motifPosId);
  
          posX = where.X() + motifPos->Position().X()/10.+fgkOffsetX;
          posY = where.Y() + motifPos->Position().Y()/10.+fgkOffsetY;
	  posZ = where.Z() + sgn * (fgkMotherThick1 - TotalHzDaughter()); 
          gMC->Gspos(fgkDaughterName, motifPosId, QuadrantMLayerName(chamber), posX, posY, posZ, reflZ, "ONLY");
        }  
        segNum++;
	
      } else { 

        // if this is a special segment	
        for (Int_t motifNum=0;motifNum<seg->GetNofMotifs();motifNum++) {// for each motif

          Int_t motifPosId = seg->GetMotifPositionId(motifNum);
          
          if (find(already_done.begin(),already_done.end(),motifPosId)
              != already_done.end()) continue; // don't treat the same motif twice
          
          AliMUONSt1SpecialMotif spMot = specialMap[motifPosId];
          AliMpMotifPosition* motifPos = sector->GetMotifMap()->FindMotifPosition(motifPosId);

          // place the hole for the motif, wrt the requested rotation angle
          Int_t rot = ( spMot.GetRotAngle()<0.1 ) ? reflZ:rotMat;

          posX = where.X() + motifPos->Position().X()/10.+spMot.GetDelta().X();
          posY = where.Y() + motifPos->Position().Y()/10.+spMot.GetDelta().Y();
          posZ = where.Z() + sgn * (TotalHzPlane() + fgkHzGas + 2.*fgkHzPadPlane);
          gMC->Gspos(fgkHoleName, motifPosId, QuadrantMLayerName(chamber), posX, posY, posZ, rot, "ONLY");

          // then place the daughter board for the motif, wrt the requested rotation angle
          posX = posX+fgkDeltaFilleEtamX;
          posY = posY+fgkDeltaFilleEtamY;
	  posZ = where.Z() + sgn * (fgkMotherThick1 - TotalHzDaughter()); 
          gMC->Gspos(fgkDaughterName, motifPosId, QuadrantMLayerName(chamber), posX, posY, posZ, rot, "ONLY");

          already_done.push_back(motifPosId);// mark this motif as done
	}		
      }// end of special motif case
    }
  }
} 

//______________________________________________________________________________
TString AliMUONv2::GasVolumeName(const TString& name, Int_t chamber) const
{
// Inserts the chamber number into the name.
// ---

  TString newString(name);
 
  TString number(""); 
  number += chamber;

  newString.Insert(2, number);
  
  return newString;
}

//______________________________________________________________________________
Bool_t AliMUONv2::IsInChamber(Int_t ich, Int_t volGid) const
{
// True if volume <volGid> is part of the sensitive 
// volumes of chamber <ich> 
// ---
  for (Int_t i = 0; i < fChamberV2[ich]->GetSize(); i++) {
      if (fChamberV2[ich]->At(i) == volGid) return kTRUE;
  }
  return kFALSE;
}

//
// protected methods
//

//______________________________________________________________________________
Int_t  AliMUONv2::GetChamberId(Int_t volId) const
{
// Check if the volume with specified  volId is a sensitive volume (gas) 
// of some chamber and returns the chamber number;
// if not sensitive volume - return 0.
// ---

  for (Int_t i = 1; i <=2; i++) 
     if (IsInChamber(i-1,volId)) return i;
  
  for (Int_t i = 3; i <= AliMUONConstants::NCh(); i++)
    if (volId==((AliMUONChamber*)(*fChambers)[i-1])->GetGid()) return i;

  return 0;
}

//
// public methods
//

//______________________________________________________________________________
void AliMUONv2::CreateMaterials()
{
// --- Define the various mixtures for GEANT ---
  
  //     Ar-CO2 gas (80%+20%)
  Float_t ag1[2]   = { 39.95,44.01};
  Float_t zg1[2]   = { 18.,22.};
  Float_t dg1      = .001821;
  Float_t wg1[2]   = { .8,0.2};
  // use wg1 weighting factors (6th arg > 0)
  AliMixture(22, "ArCO2 80%$", ag1, zg1, dg1, 2, wg1);  
  
  // Ar-buthane-freon gas -- trigger chambers 
  Float_t atr1[4]  = { 39.95,12.01,1.01,19. };
  Float_t ztr1[4]  = { 18.,6.,1.,9. };
  Float_t wtr1[4]  = { .56,.1262857,.2857143,.028 };
  Float_t dtr1     = .002599;
  AliMixture(23, "Ar-freon $", atr1, ztr1, dtr1, 4, wtr1);

  // Rohacell 51  - imide methacrylique
  Float_t aRohacell51[4] = {12.01,1.01,16.00,14.01}; 
  Float_t zRohacell51[4] = {6.,1.,8.,7.}; 
  Float_t dRohacell51 = 0.052;
  Float_t wRohacell51[4] = {9.,13.,2.,1.};  
  // use relative A (molecular) values (6th arg < 0)
  AliMixture(32, "FOAM$",aRohacell51,zRohacell51,dRohacell51,-4,wRohacell51);  
   
  Float_t aSnPb[2] = {118.69,207.19};
  Float_t zSnPb[2] = {50,82};
  Float_t dSnPb = 8.926;
  Float_t wSnPb[2] = {0.6, 0.4} ;
  // use wSnPb weighting factors (6th arg > 0)
  AliMixture(35, "SnPb$", aSnPb,zSnPb,dSnPb,2,wSnPb);

  // plastic definition from K5, Freiburg (found on web)
  Float_t aPlastic[2]={1.01,12.01};
  Float_t zPlastic[2]={1,6};
  Float_t denPlastic=1.107;
  Float_t wPlastic[2]={1,1};
  // use relative A (molecular) values (6th arg < 0)...no other info...
  AliMixture( 33, "Plastic$",aPlastic,zPlastic,denPlastic,-2,wPlastic);
 
  // from CERN note NUFACT Note023, Oct.2000 
  // Inox/Stainless Steel (18%Cr, 9%Ni)
  Float_t aInox[3] = {55.847,51.9961,58.6934};  
  Float_t zInox[3] = {26.,24.,28.};
  Float_t denInox = 7.930;
  Float_t wInox[3] = {0.73,0.18,0.09}; 
  // use wInox weighting factors (6th arg > 0) 
  AliMixture(37, "StainlessSteel$",aInox,zInox,denInox,3,wInox);   

  // bakelite 
  Float_t abak[3] = {12.01 , 1.01 , 16.};
  Float_t zbak[3] = {6.     , 1.   , 8.};
  Float_t wbak[3] = {6.     , 6.   , 1.}; 
  Float_t dbak = 1.4;
  AliMixture(19, "Bakelite$", abak, zbak, dbak, -3, wbak);
  
  // Ar-Isobutane gas (80%+20%) 
  Float_t ag[3]    = { 39.95,12.01,1.01 };
  Float_t zg[3]    = { 18.,6.,1. };
  Float_t wg[3]    = { .8,.057,.143 };
  Float_t dg       = .0019596;
  AliMixture(20, "ArC4H10 GAS$", ag, zg, dg, 3, wg);

  //     Ar-Isobutane-Forane-SF6 gas (49%+7%+40%+4%) -- trigger 
  Float_t atrig[5] = { 39.95,12.01,1.01,19.,32.066 };
  Float_t ztrig[5] = { 18.,6.,1.,9.,16. };
  Float_t wtrig[5] = { .49,1.08,1.5,1.84,0.04 };
  Float_t dtrig    = .0031463;
  AliMixture(21, "TRIG GAS$", atrig, ztrig, dtrig, -5, wtrig);
      
// --- Define the various AliMaterials for GEANT ---
  // from PDG and "The Particle Detector BriefBook", Bock and Vasilescu, P.18  
  AliMaterial( 9, "Aluminium$", 26.98, 13., 2.7, -8.9, 26.1);
  AliMaterial(10, "Aluminium$", 26.98, 13., 2.7, -8.9, 26.1);
  AliMaterial(15, "air$", 14.61, 7.3, .001205, -30423.24, 67500);
  AliMaterial(30, "Copper$", 63.546,29.,8.96,-1.43,9.6);
  AliMaterial(31, "FR4$", 17.749, 8.875, 1.7, -19.4, 999.);    // from DPG
  AliMaterial(34, "Kapton$", 12.01,6,1.42,-28.6,999);          // from DPG
 // Density of FrameEpoxy only from manufacturer's specifications
 // Frame composite epoxy , X0 in g/cm**2 (guestimation!)
  AliMaterial(36, "FrameEpoxy",12.24,6.0,1.85,-19.14,999);// use 16.75cm
  
// --- Define the tracking medias (AliMediums) for GEANT ---
  GReal_t epsil  = .001;       // Tracking precision,
  GReal_t stemax = -1.;        // Maximum displacement for multiple scat
  GReal_t tmaxfd = -20.;       // Maximum angle due to field deflection
  GReal_t deemax = -.3;        // Maximum fractional energy loss, DLS
  GReal_t stmin  = -.8;
  GReal_t maxStepAlu = 0.001;  // from AliMUON.cxx
  GReal_t maxDestepAlu = -1.;  // from AliMUON.cxx
  GReal_t maxStepGas=0.01;     // from AliMUON.cxx

  Int_t iSXFLD   = gAlice->Field()->Integ();
  Float_t sXMGMX = gAlice->Field()->Max();

  AliMedium(1, "AIR_CH_US$", 15, 1, iSXFLD, sXMGMX, tmaxfd,
	    stemax, deemax, epsil, stmin);
  AliMedium(4, "ALU_CH_US$", 9, 0, iSXFLD, sXMGMX, tmaxfd,
	    maxStepAlu, maxDestepAlu, epsil, stmin);
  AliMedium(5, "ALU_CH_US$", 10, 0, iSXFLD, sXMGMX, tmaxfd,
	    maxStepAlu,maxDestepAlu, epsil, stmin);  
  AliMedium(6, "AR_CH_US          ", 20, 1, iSXFLD, sXMGMX, 
              tmaxfd, fMaxStepGas,fMaxDestepGas, epsil, stmin);
              
  //    Ar-Isobuthane-Forane-SF6 gas 
  AliMedium(7, "GAS_CH_TRIGGER    ", 21, 1, iSXFLD, sXMGMX, 
               tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(8, "BAKE_CH_TRIGGER   ", 19, 0, iSXFLD, sXMGMX, 
               tmaxfd, fMaxStepAlu, fMaxDestepAlu, epsil, stmin);

  AliMedium(9, "ArCO2 80%$", 22, 1, iSXFLD, sXMGMX, tmaxfd, maxStepGas,
	      maxDestepAlu, epsil, stmin);
  AliMedium(10, "COPPER_CH$", 30, 0, iSXFLD, sXMGMX, tmaxfd,
	      maxStepAlu, maxDestepAlu, epsil, stmin);
  AliMedium(11, "PCB_COPPER        ", 31, 0, iSXFLD, sXMGMX, tmaxfd, 
                fMaxStepAlu, fMaxDestepAlu, epsil, stmin);
  AliMedium(12, "VETRONITE         ", 32, 0, iSXFLD, sXMGMX, tmaxfd, 
                fMaxStepAlu, fMaxDestepAlu, epsil, stmin);
  AliMedium(13, "CARBON            ", 33, 0, iSXFLD, sXMGMX, tmaxfd, 
                fMaxStepAlu, fMaxDestepAlu, epsil, stmin);
  AliMedium(14, "Rohacell          ", 34, 0, iSXFLD, sXMGMX, tmaxfd, 
	      fMaxStepAlu, fMaxDestepAlu, epsil, stmin);
  AliMedium(15, "FR4_CH$",  31, 0,iSXFLD, sXMGMX, 10., .01,.1, .003, .003);
  AliMedium(16, "FOAM_CH$", 32, 0,
	      iSXFLD, sXMGMX, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;
  AliMedium(17, "Plastic$", 33, 0,iSXFLD, sXMGMX,  10., .01, 1., .003, .003);
  AliMedium(18, "Kapton$", 34, 0,iSXFLD, sXMGMX,  10., .01, 1., .003, .003);
  AliMedium(19, "SnPb$", 35, 0,iSXFLD, sXMGMX,  10., .01,  1., .003, .003);
  AliMedium(20, "FrameCH$", 36, 1,iSXFLD, sXMGMX,  10., .001, 0.001, .001, .001);
  AliMedium(21, "InoxBolts$", 37,1,iSXFLD, sXMGMX,  10., .01, 1., .003, .003);

// store the parameters
  Float_t a, z, dens, absl;
  char matName[5];
  AliGetMaterial(30,matName,a,z,dens,fRadlCopper,absl);
  AliGetMaterial(31,matName,a,z,dens,fRadlFR4,absl);
  AliGetMaterial(32,matName,a,z,dens,fRadlFoam,absl);
}

//______________________________________________________________________________
void AliMUONv2::CreateGeometry()
{
// Create the GEANT geometry for the dimuon arm.
// Use the parent's method for stations 2, 3, 4 and 5.
// Use the detailed code for the first station.
// --
  cout << "AliMUONv2::CreateGeometry()" << endl;
  cout << "_________________________________________" << endl;

  // Create basic volumes
  // 
  CreateHole();
  CreateDaughterBoard();
  CreateInnerLayers();
  
  // Create reflexion matrices
  //
  Int_t reflXZ, reflYZ, reflXY;
  AliMatrix(reflXZ,  90.,  180., 90., 90., 180., 0.);
  AliMatrix(reflYZ,  90., 0., 90.,-90., 180., 0.);
  AliMatrix(reflXY,  90., 180., 90., 270., 0., 0.);

  // Define transformations for each quadrant
  // 
  //     II. |  I.
  //   _____ | ____
  //         |
  //    III. |  IV.
  // 
  Int_t rotm[4];
  rotm[0]=0;       // quadrant I
  rotm[1]=reflXZ;  // quadrant II
  rotm[2]=reflXY;  // quadrant III
  rotm[3]=reflYZ;  // quadrant IV
  
  TVector3 scale[4];  
  scale[0] = TVector3( 1,  1,  1);  // quadrant I
  scale[1] = TVector3(-1,  1, -1);  // quadrant II
  scale[2] = TVector3(-1, -1,  1);  // quadrant III
  scale[3] = TVector3( 1, -1, -1);  // quadrant IV
  
  // Shift in Z of the middle layer
  Double_t deltaZ = 6.5/2.;         

  // Position of quadrant I wrt to the chamber position
  TVector3 pos0(-fgkDeltaQuadLHC, -fgkDeltaQuadLHC, deltaZ);

  // Shift for near/far layers
  GReal_t  shiftXY = fgkFrameOffset;
  GReal_t  shiftZ  = fgkMotherThick1+fgkMotherThick2;

  // Build two chambers
  //
  for (Int_t ich=1; ich<3; ich++) {

    // Create quadrant volume
    CreateQuadrant(ich);

    // Place gas volumes
    PlaceInnerLayers(ich);
    
    // Place the quadrant
    for (Int_t i=0; i<4; i++) {

      // Middle layer
      GReal_t posx = pos0.X() * scale[i].X();
      GReal_t posy = pos0.Y() * scale[i].Y();
      GReal_t posz = pos0.Z() * scale[i].Z() + AliMUONConstants::DefaultChamberZ(ich-1);
      gMC->Gspos(QuadrantMLayerName(ich), i+1, "ALIC", posx, posy, posz, rotm[i], "ONLY");

      // Near/far layers
      Real_t  posx2 = posx + shiftXY * scale[i].X();
      Real_t  posy2 = posy + shiftXY * scale[i].Y();
      Real_t  posz2 = posz - scale[i].Z()*shiftZ;
      gMC->Gspos(QuadrantNLayerName(ich), i+1, "ALIC", posx2, posy2, posz2, rotm[i],"ONLY");
    
      posz2 = posz + scale[i].Z()*shiftZ;
      gMC->Gspos(QuadrantFLayerName(ich), i+1, "ALIC", posx2, posy2, posz2, rotm[i],"ONLY");
   }
 }     

  static Int_t stations[5]={0,1,1,1,1};
  fStations=stations;
  AliMUONv1::CreateGeometry();
}

//______________________________________________________________________________
void AliMUONv2::Init() 
{
   // Initialize Station 1 Tracking Chambers
 
   //
   // Set the chamber (sensitive region) GEANT identifier
    fChamberV2[0] = new TArrayI(11);             // Chamber 1 sensitive volume Id array
    fChamberV2[1] = new TArrayI(11);             // Chamber 2 sensitive volume Id array

    AddChamberGid(0,gMC->VolId("SA1G"),0);
    AddChamberGid(0,gMC->VolId("SB1G"),1);
    AddChamberGid(0,gMC->VolId("SC1G"),2);
    AddChamberGid(0,gMC->VolId("SD1G"),3);
    AddChamberGid(0,gMC->VolId("SE1G"),4);
    AddChamberGid(0,gMC->VolId("SF1G"),5);
    AddChamberGid(0,gMC->VolId("SG1G"),6);
    AddChamberGid(0,gMC->VolId("SH1G"),7);
    AddChamberGid(0,gMC->VolId("SI1G"),8);
    AddChamberGid(0,gMC->VolId("SJ1G"),9);
    AddChamberGid(0,gMC->VolId("SK1G"),10);

    
    AddChamberGid(1,gMC->VolId("SA2G"),0);
    AddChamberGid(1,gMC->VolId("SB2G"),1);
    AddChamberGid(1,gMC->VolId("SC2G"),2);
    AddChamberGid(1,gMC->VolId("SD2G"),3);
    AddChamberGid(1,gMC->VolId("SE2G"),4);
    AddChamberGid(1,gMC->VolId("SF2G"),5);
    AddChamberGid(1,gMC->VolId("SG2G"),6);
    AddChamberGid(1,gMC->VolId("SH2G"),7);
    AddChamberGid(1,gMC->VolId("SI2G"),8);
    AddChamberGid(1,gMC->VolId("SJ2G"),9);
    AddChamberGid(1,gMC->VolId("SK2G"),10);

  Int_t i;
// now do the other stations as in AliMUONv1
   for (i=0; i<AliMUONConstants::NCh(); i++) {
       ( (AliMUONChamber*) (*fChambers)[i])->Init();
   }
   
   //
   // Set the chamber (sensitive region) GEANT identifier
   ((AliMUONChamber*)(*fChambers)[0])->SetGid(-1);  // joker
   ((AliMUONChamber*)(*fChambers)[1])->SetGid(-1);  // joker

   ((AliMUONChamber*)(*fChambers)[2])->SetGid(gMC->VolId("S03G"));
   ((AliMUONChamber*)(*fChambers)[3])->SetGid(gMC->VolId("S04G"));

   ((AliMUONChamber*)(*fChambers)[4])->SetGid(gMC->VolId("S05G"));
   ((AliMUONChamber*)(*fChambers)[5])->SetGid(gMC->VolId("S06G"));

   ((AliMUONChamber*)(*fChambers)[6])->SetGid(gMC->VolId("S07G"));
   ((AliMUONChamber*)(*fChambers)[7])->SetGid(gMC->VolId("S08G"));

   ((AliMUONChamber*)(*fChambers)[8])->SetGid(gMC->VolId("S09G"));
   ((AliMUONChamber*)(*fChambers)[9])->SetGid(gMC->VolId("S10G"));

   ((AliMUONChamber*)(*fChambers)[10])->SetGid(gMC->VolId("SG1A"));
   ((AliMUONChamber*)(*fChambers)[11])->SetGid(gMC->VolId("SG2A"));
   ((AliMUONChamber*)(*fChambers)[12])->SetGid(gMC->VolId("SG3A"));
   ((AliMUONChamber*)(*fChambers)[13])->SetGid(gMC->VolId("SG4A"));

}

