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
*/

// Authors: David Guez, Ivana Hrivnacova, Marion MacCormick; IPN Orsay
//
// Class AliMUONv2
// ---------------
// Inherits from AliMUONv1 but with a more detailed
// geometrical description of station 1 

#include <algorithm>
#include <vector>

#include <string.h>

#include <TVector2.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>


#include "AliMUONv2.h"
#include "AliMUONConstants.h"
#include "AliMUONHit.h"
#include "AliRun.h"
#include "AliMagF.h"
#include "AliConst.h" 

#include <MReader.h>
#include <MSector.h>
#include <MRow.h>
#include <MVRowSegment.h>
#include <MVMotif.h>
#include <MMotifMap.h>
#include <MMotifPosition.h>
#include <MMotifType.h>
#include <MIntPair.h>


#include <Riostream.h>

ClassImp(AliMUONv2)



const GReal_t AliMUONv2::fgkHzPadPlane=0.0148/2.; //Pad plane
const GReal_t AliMUONv2::fgkHzFoam = 2.083/2.; // Foam of mechanicalplane
const GReal_t AliMUONv2::fgkHzFR4 = 0.0031/2.; // FR4 of mechanical plane
const GReal_t AliMUONv2::fgkHzSnPb = 0.0091/2.; //Pad/Kapton connection (66 pt)
const GReal_t AliMUONv2::fgkHzKapton = 0.0122/2.; //Kapton
const GReal_t AliMUONv2::fgkHzBergPlastic = 0.3062/2.; //Berg connector
const GReal_t AliMUONv2::fgkHzBergCopper = 0.1882/2.; //Berg connector
const GReal_t AliMUONv2::fgkHzDaughter = 0.0156/2.; //Daughter board
const GReal_t AliMUONv2::fgkHzGas = 0.2/2.; //Gas
const GReal_t AliMUONv2::fgkHxQuadrant = 94.69/2.; //Box surrounding quadrant
const GReal_t AliMUONv2::fgkHyQuadrant = 100.31/2.; //Box surrounding quadrant
const GReal_t AliMUONv2::fgkMotherIR = 18.3;
const GReal_t AliMUONv2::fgkMotherOR = 104.974;   
const GReal_t AliMUONv2::fgkMotherThick = 6.5/2; //6.5cm between two quadrants 
const GReal_t AliMUONv2::fgkMotherPhiL = 0.; 
const GReal_t AliMUONv2::fgkMotherPhiU = 90.;


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
const GReal_t AliMUONv2::fgkDeltaFilleEtamX=0.;
const GReal_t AliMUONv2::fgkDeltaFilleEtamY=0.5;

const char* AliMUONv2::fgkHoleName="MCHL";
const char* AliMUONv2::fgkDaughterName="MCDB";
const char  AliMUONv2::fgkFoamLayerSuffix='F';
const char* AliMUONv2::fgkQuadrantName="QUA";

//___________________________________________
AliMUONv2::AliMUONv2()
  : AliMUONv1(),
    fIdSens(0)
{
// Default Constructor
//

   // keep secondaries
   SetIshunt(0);
}
 
//___________________________________________
AliMUONv2::AliMUONv2(const char *name, const char *title)
  : AliMUONv1(name,title),
    fIdSens(0)
{
   // keep secondaries
   SetIshunt(0);

   // create hits array
   fHits     = new TClonesArray("AliMUONHit",1000);
   gAlice->AddHitList(fHits);
}
 
//___________________________________________
AliMUONv2::AliMUONv2(const AliMUONv2& rMUON)
{
// Dummy copy constructor
}

//___________________________________________
AliMUONv2::~AliMUONv2()
{
// Destructor
    if(fDebug) printf("%s: Calling AliMUONv2 destructor !!!\n",ClassName());
    
    if (fHits) {
      fHits->Delete();
      delete fHits;
    }
}

//___________________________________________
void AliMUONv2::CreateGeometry()
{
// Create the GEANT geometry for the dimuon arm.
// Use the parent's method for stations 2, 3, 4 and 5.
// Use the detailed code for the first station.
// --
  cout << "AliMUONv2::CreateGeometry()" << endl;
  cout << "_________________________________________" << endl;

  // tracking medias
  Int_t* idtmed = fIdtmed->GetArray()-1099;
  Int_t idAir  = idtmed[1100];       // medium 1



  //create reflexion matrix
  Int_t reflXZ,reflYZ,reflXY;
  AliMatrix(reflXZ,  90.,  180., 90., 90., 180., 0.);
  AliMatrix(reflYZ,  90., 0., 90.,-90., 180., 0.);
  AliMatrix(reflXY,  90., 180., 90., 270., 0., 0.);

  CreateHole();
  CreateDaughterBoard();

  //build the quadrant
  CreateQuadrant(1);
  CreateQuadrant(2);
  GReal_t par[5];
  par[0] = fgkMotherIR;
  par[1] = fgkMotherOR; 
  par[2] = fgkHzGas;  
  par[3] = fgkMotherPhiL; 
  par[4] = fgkMotherPhiU;

  gMC->Gsvolu("S01G","TUBS",idAir,par,5);
  gMC->Gspos("S01G",1,QuadrantName(1),0,0,0,0,"ONLY");
  gMC->Gsvolu("S02G","TUBS",idAir,par,5);
  gMC->Gspos("S02G",1,QuadrantName(2),0,0,0,0,"ONLY");

  // place the four copies of it
  
  //parameters
  TVector3 pos[4];
  Int_t rotm[4];
  
  Double_t deltaZ = 6.5/2.; 
  pos[0]=TVector3(-2.6,-2.6,AliMUONConstants::DefaultChamberZ(0)+deltaZ);
  rotm[0]=0;
  pos[1]=TVector3(2.6,-2.6,AliMUONConstants::DefaultChamberZ(0)-deltaZ);
  rotm[1]=reflXZ;
  pos[2]=TVector3(2.6,2.6,AliMUONConstants::DefaultChamberZ(0)+deltaZ);
  rotm[2]=reflXY;
  pos[3]=TVector3(-2.6,2.6,AliMUONConstants::DefaultChamberZ(0)-deltaZ);
  rotm[3]=reflYZ;
  
  //placing
  GReal_t posX,posY,posZ;
  for (Int_t i=0;i<4;i++){
    posX=pos[i].X();
    posY=pos[i].Y();
    // the 1st chamber
    posZ=pos[i].Z();
    gMC->Gspos(QuadrantName(1),i+1,"ALIC",posX,posY,posZ,rotm[i],"ONLY");
    // the 2nd chamber
    posZ=pos[i].Z()+AliMUONConstants::DefaultChamberZ(1)
         -AliMUONConstants::DefaultChamberZ(0);
    gMC->Gspos(QuadrantName(2),i+5,"ALIC",posX,posY,posZ,rotm[i],"ONLY");
  }
  static Int_t stations[5]={0,1,1,1,1};
  fStations=stations;
  AliMUONv1::CreateGeometry();
}
//___________________________________________
void AliMUONv2::CreateQuadrant(Int_t chamber)
{
// create the quadrant (bending and non-bending planes)
// for the given chamber
// --
  CreateFrame(chamber);

  TSpecialMap specialMap;  

  TVector3 where;
  specialMap[1001] = AliMUONSt1SpecialMotif(TVector2(0.1 ,0.84),90.);
  specialMap[1002] = AliMUONSt1SpecialMotif(TVector2(0.5 ,0.76));
  specialMap[1003] = AliMUONSt1SpecialMotif(TVector2(1.01,0.76));
  MReader reader1(kBendingPlane);
  MSector* sector1 = reader1.BuildSector();
  where=TVector3(0.185+2.6,-0.52+2.6,totalHz()+fgkHzGas);
  PlaceSector(sector1,specialMap,where,chamber);
  

  specialMap.clear();
  specialMap[4001] = AliMUONSt1SpecialMotif(TVector2(1.01,0.59),90.);
  specialMap[4002] = AliMUONSt1SpecialMotif(TVector2(1.96, 0.17));
  specialMap[4003] = AliMUONSt1SpecialMotif(TVector2(1.61,-1.18));
  specialMap[4004] = AliMUONSt1SpecialMotif(TVector2(0.2 ,-0.08));
  specialMap[4005] = AliMUONSt1SpecialMotif(TVector2(0.2 , 0.25));
  specialMap[4006] = AliMUONSt1SpecialMotif(TVector2(0.28, 0.21));

  MReader reader2(kNonBendingPlane);
  MSector* sector2 = reader2.BuildSector();
  where=TVector3(-0.13+2.6,-0.1775+2.6,-totalHz()-fgkHzGas);
  PlaceSector(sector2,specialMap,where,chamber);
}

//___________________________________________
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

  //	 Ar-buthane-freon gas -- trigger chambers 
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
  
  //     Ar-CO2 gas 
  Float_t agas[3]  = { 39.95,12.01,16. };
  Float_t zgas[3]  = { 18.,6.,8. };
  Float_t wgas[3]  = { .74,.086684,.173316 };
  Float_t dgas     = .0018327;
  AliMixture(24, "ArCO2 GAS$", agas, zgas, dgas, 3, wgas); 
   
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
  AliMaterial(36, "FrameEpoxy",12.24,6.0,1.85,-19.14,999);
  
// --- Define the tracking medias (AliMediums) for GEANT ---
  GReal_t epsil  = .001; // Tracking precision,
  GReal_t stemax = -1.;  // Maximum displacement for multiple scat
  GReal_t tmaxfd = -20.; // Maximum angle due to field deflection
  GReal_t deemax = -.3;  // Maximum fractional energy loss, DLS
  GReal_t stmin  = -.8;
  GReal_t maxStepAlu = 0.001;  // from AliMUON.cxx
  GReal_t maxDestepAlu = -1.; // from AliMUON.cxx
  GReal_t maxStepGas=0.01;  // from AliMUON.cxx

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

  AliMedium(9, "ARG_CO2$", 22, 1, iSXFLD, sXMGMX, tmaxfd, maxStepGas,
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
}

void AliMUONv2::CreateFrame(Int_t chamber)
{
// Create the non-sensitive elements of the frame for the  <chamber>


// Matrices
      Int_t idrotm[1199];   
// To Be Checked....
//   Rotation matrices in the x-y plane  
//   phi =   -45 deg
     AliMatrix(idrotm[1101],  90.,   315., 90.,  45., 0., 0.);
//   phi =  -90 deg
     AliMatrix(idrotm[1102],  90.,  90., 90., 180., 0., 0.);
//   theta =  +90 deg
     AliMatrix(idrotm[1103],  180.,  0., 90., 90.,90., 0.);
      
//   phi =   +45 deg
     AliMatrix(idrotm[1104],  90.,   45., 90.,  135., 0., 0.);
//   phi =  +45 deg + rotation 180° around Y 
     AliMatrix(idrotm[1105],  90.,  45., 90., 315., 180., 0.);

//   Translation matrices in the x-y plane  
//   X -> X ; Y -> Y; Z -> Z
     AliMatrix(idrotm[1110],  90.,   0., 90.,  90., 0., 0.);
//   X->-X; Y -> Y; Z -> -Z
     AliMatrix(idrotm[1111],  90.,  180., 90., 90., 180., 0.);
//   X->-X; Y ->-Y; Z -> Z
     AliMatrix(idrotm[1112],  90., 180., 90., 270., 0., 0.);
//   X->X; Y ->-Y; Z -> -Z
     AliMatrix(idrotm[1113],  90., 0., 90., 270., 180., 0.);
//   
 
  // tracking medias
  Int_t* idtmed = fIdtmed->GetArray()-1099;

  Int_t idAir  = idtmed[1100];       // medium 1
  Int_t idSolder = idtmed[1118];     // medium 19
  Int_t idFrameEpoxy = idtmed[1119]; // medium 20 = Frame Epoxy ME730
  Int_t idInox = idtmed[1120];       // medium 21 Stainless Steel (18%Cr,9%Ni,Fe)

//________________________________________________________________
// 
// Original model:
//
// Epoxy Frame segments
//  
//                         OutTopTrapFrame
//                              ------------ |
//      OutEdgeTrapFrame      /              |
//                           /               |  InVFrame
//       OutCornerTrapFrame /                |
//                          |              --|   
//             OutVFrame    |            _-  InArcFrame
//                          |___________|
//                             InHFrame
//                          
//
// DATE: 27 NOVEMBER 2002 - MODEL UPDATED. THE MOTHER VOLUME IS A TUBS.
//__________________________________________________________________

  const Float_t hzFrameThickness = 1.186/2.;   //equivalent thickness
  const Float_t hzOuterFrameEpoxy = 1.23/2.;   //equivalent thickness
  const Float_t hzOuterFrameSolder = 0.032/2.; //equivalent thickness
  const Float_t hzOuterFrameInox = 0.035/2.;   //equivalent thickness
      
  // InHFrame parameters
  const Float_t hxInHFrame  = 75.8/2.;
  const Float_t hyInHFrame  = 2.5/2.;
  const Float_t hzInHFrame  = hzFrameThickness;
 
  // InVFrame parameters
  const Float_t hxInVFrame  = 2.5/2.;
  const Float_t hyInVFrame  = 73.3/2.;
  const Float_t hzInVFrame  = hzFrameThickness;

  //Flat 1mm vertical section
  const Float_t hxV1mm  = 0.1/2.;
  const Float_t hyV1mm  = 2.5/2.;
  const Float_t hzV1mm  = hzFrameThickness;
  
  //Flat 1mm horizontal section
  const Float_t hxH1mm  = 2.5/2.;
  const Float_t hyH1mm  = 0.1/2.;
  const Float_t hzH1mm  = hzFrameThickness;
   
  // InArcFrame parameters
  const Float_t IAF  = 15.70;
  const Float_t OAF  = 18.20;
  const Float_t hzAF  = hzFrameThickness;
  const Float_t AFphi1  = 0.0;
  const Float_t AFphi2  = 90.0;

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

  // OutVFrame parameters
  const Float_t hxOutVFrame  = 2.5/2.;
  const Float_t hyOutVFrame  = 43.927/2.;
  const Float_t hzOutVFrame  = hzFrameThickness;
  
/////////////////////////////////////////////////////// 

// OutTopCuboidFrame Epoxy parameters
  const Float_t hxOutHFrame  = 31.9304/2.;
  const Float_t hyOutHFrame  = 8.4/2.;
  const Float_t hzOutHFrame  = hzOuterFrameEpoxy;
  
// OutTopTrapFrameA Epoxy parameters
  const Float_t hzOTTFA = hzOuterFrameEpoxy;
  const Float_t tetOTTFA = 0.;
  const Float_t phiOTTFA = 0.;
  const Float_t h1OTTFA = 9.6716/2.;
  const Float_t bl1OTTFA = 4.7787/2.;
  const Float_t tl1OTTFA =  8.4/2.;
  const Float_t alp1OTTFA = 10.6; 
  const Float_t h2OTTFA = 9.6716/2.;
  const Float_t bl2OTTFA = 4.7787/2.;
  const Float_t tl2OTTFA = 8.4/2.;
  const Float_t alp2OTTFA = 10.6;  
  
// OutTopTrapFrameB Epoxy parameters
  const Float_t hzOTTFB = hzOuterFrameEpoxy;
  const Float_t tetOTTFB = 0.;
  const Float_t phiOTTFB = 0.;
  const Float_t h1OTTFB = 9.6716/2.;
  const Float_t bl1OTTFB = 0.;
  const Float_t tl1OTTFB =  4.7787/2.;
  const Float_t alp1OTTFB = 13.88; 
  const Float_t h2OTTFB = 9.6716/2.;
  const Float_t bl2OTTFB = 0.;
  const Float_t tl2OTTFB = 4.7787/2.;
  const Float_t alp2OTTFB = 13.88;

  // OutTopTrapFrame Solder parameters
  const Float_t hzOTTFS = hzOuterFrameSolder;
  
  // OutTopTrapFrame Inox parameters
  const Float_t hzOTTFI = hzOuterFrameInox;
 
/////////////////////////////////////////////////////////  

  // OutEdgeTrapFrame parameters 
  // TRAPEZE 1
  const Float_t hzOETF = hzOuterFrameEpoxy;   
  const Float_t tetOETF = 0.;
  const Float_t phiOETF = 0.;
  const Float_t h1OETF = 7.129/2.; 
  const Float_t bl1OETF1 = 3.705/2; 
  const Float_t tl1OETF1 = 3.9471/2.;
  const Float_t alp1OETF1 = 0.97;
  const Float_t h2OETF = 7.129/2.;    
  const Float_t bl2OETF1 = 3.705/2;
  const Float_t tl2OETF1 = 3.9471/2.;
  const Float_t alp2OETF1 = 0.97;
  
  // TRAPEZE 2
  const Float_t bl1OETF2 = 2.9744/2.;
  const Float_t tl1OETF2 =  3.705/2;
  const Float_t alp1OETF2 =  2.93;
      
  const Float_t bl2OETF2 =  2.9744/2.;
  const Float_t tl2OETF2 =  3.705/2;
  const Float_t alp2OETF2 =  2.93; 
   
  // TRAPEZE 3
  const Float_t bl1OETF3 = 1.7455/2.;
  const Float_t tl1OETF3 =  2.9744/2.;
  const Float_t alp1OETF3 =  4.93;
      
  const Float_t bl2OETF3 = 1.7455/2.;
  const Float_t tl2OETF3 = 2.9744/2.; 
  const Float_t alp2OETF3 =  4.93; 
  
  // TRAPEZE 4
  const Float_t bl1OETF4 = 0.;
  const Float_t tl1OETF4 =  1.7455/2.;
  const Float_t alp1OETF4 =  6.98;
      
  const Float_t bl2OETF4 = 0.;
  const Float_t tl2OETF4 = 1.7455/2.;
  const Float_t alp2OETF4 =  6.98;   
 
 ///////////////////////////////////////////////////
    
  // OutCornerTrapFrame parameters
  const Float_t hzOCTF = hzFrameThickness;
  const Float_t tetOCTF = 0.;
  const Float_t phiOCTF = 0.;
  const Float_t h1OCTF = 2.5/2.;
  const Float_t bl1OCTF = 0.;
  const Float_t tl1OCTF =  4.7469/2.;
  const Float_t alp1OCTF = 43.51; 
  const Float_t h2OCTF = 2.5/2.;
  const Float_t bl2OCTF = 0.;
  const Float_t tl2OCTF = 4.7469/2.;
  const Float_t alp2OCTF = 43.51;  

  // Reference point - MIRE (3 per quadrant, only 1 programmed)
  const Float_t MIREInRad  = 0.6;
  const Float_t MIREOutRad  = 1.3;
  const Float_t MIRELen  = hzFrameThickness;
       
  Float_t par[11];
  Float_t posX,posY,posZ;
       
// ___________________Make volumes________________________

// Quadrant volume
    par[0] = fgkMotherIR;
    par[1] = fgkMotherOR; 
    par[2] = fgkMotherThick;  
    par[3] = fgkMotherPhiL; 
    par[4] = fgkMotherPhiU;

    //Quadrant volume.....positionned at the end
    gMC->Gsvolu(QuadrantName(chamber),"TUBS",idAir,par,5);
            
// _______________________________________________________
// InVFrame  
   if (chamber==1) {   
      // simple epoxy layer...must be inside sensitive surface for tracking
      par[0] = hxInVFrame;
      par[1] = hyInVFrame;
      par[2] = hzInVFrame;
      
      gMC->Gsvolu("IVEF","BOX",idFrameEpoxy,par,3);

  // InHFrame
      par[0] = hxInHFrame;
      par[1] = hyInHFrame;
      par[2] = hzInHFrame;

      gMC->Gsvolu("IHEF","BOX",idFrameEpoxy,par,3);

  //Flat 1mm vertical section
      par[0] = hxV1mm;
      par[1] = hyV1mm;
      par[2] = hzV1mm;

      gMC->Gsvolu("FVMM","BOX",idFrameEpoxy,par,3);

  //Flat 1mm vertical section
      par[0] = hxH1mm;
      par[1] = hyH1mm;
      par[2] = hzH1mm;

      gMC->Gsvolu("FHMM","BOX",idFrameEpoxy,par,3);

  // OutVFrame
      par[0] = hxOutVFrame;
      par[1] = hyOutVFrame;
      par[2] = hzOutVFrame;

      gMC->Gsvolu("OVEF","BOX",idFrameEpoxy,par,3);

  // InArcFrame 
      par[0] = IAF;
      par[1] = OAF; 
      par[2] = hzAF;  
      par[3] = AFphi1; 
      par[4] = AFphi2 ;

      gMC->Gsvolu("IAF1","TUBS",idFrameEpoxy,par,5);

  // ScrewsInFrame - 3 sections in order to avoid overlapping volumes
  // Screw Head, in air
      par[0] = SCRUHMI;
      par[1] = SCRUHMA; 
      par[2] = SCRUHLE;  

      gMC->Gsvolu("SCRH","TUBE",idInox,par,3);

  // Middle part, in the Epoxy
      par[0] = SCRUMMI;
      par[1] = SCRUMMA;
      par[2] = SCRUMLE;
      
      gMC->Gsvolu("SCRM","TUBE",idInox,par,3);

  // Screw nut, in air
      par[0] = SCRUNMI;
      par[1] = SCRUNMA;
      par[2] = SCRUNLE; 
        
      gMC->Gsvolu("SCRN","TUBE",idInox,par,3);   

  ///////////////////////////////////////////

  // OutTopTrapFrame Epoxy 
  // - 3 components (cuboid and 2 trapezes) and 3 layers (Epoxy/Inox/Solder)

  // OutTopCuboidFrame Epoxy 
      par[0] = hxOutHFrame;
      par[1] = hyOutHFrame;
      par[2] = hzOutHFrame;
      
      gMC->Gsvolu("OHEA","BOX",idFrameEpoxy,par,3);

  // OutTopTrapFrameA Epoxy parameters    
      par[0] = hzOTTFA;
      par[1] = tetOTTFA;
      par[2] = phiOTTFA;
      par[3] = h1OTTFA;
      par[4] = bl1OTTFA;
      par[5] = tl1OTTFA;
      par[6] = alp1OTTFA;
      par[7] = h2OTTFA;
      par[8] = bl2OTTFA;
      par[9] = tl2OTTFA;
      par[10] = alp2OTTFA;

      gMC->Gsvolu("OHEB","TRAP",idFrameEpoxy,par,11);    

  // OutTopTrapFrameB Epoxy parameters    
      par[0] = hzOTTFB;
      par[1] = tetOTTFB;
      par[2] = phiOTTFB;
      par[3] = h1OTTFB;
      par[4] = bl1OTTFB;
      par[5] = tl1OTTFB;
      par[6] = alp1OTTFB;
      par[7] = h2OTTFB;
      par[8] = bl2OTTFB;
      par[9] = tl2OTTFB;
      par[10] = alp2OTTFB;

      gMC->Gsvolu("OHEC","TRAP",idFrameEpoxy,par,11);    

  // OutTopCuboidFrame Solder 
      par[0] = hxOutHFrame;
      par[1] = hyOutHFrame;
      par[2] = hzOTTFS;
      gMC->Gsvolu("OHSA","BOX",idSolder,par,3);

  // OutTopTrapFrameA Solder 
      par[0] = hzOTTFS;
      par[1] = tetOTTFA;
      par[2] = phiOTTFA;
      par[3] = h1OTTFA;
      par[4] = bl1OTTFA;
      par[5] = tl1OTTFA;
      par[6] = alp1OTTFA;
      par[7] = h2OTTFA;
      par[8] = bl2OTTFA;
      par[9] = tl2OTTFA;
      par[10] = alp2OTTFA;

      gMC->Gsvolu("OHSB","TRAP",idSolder,par,11); 


  // OutTopTrapFrameB Solder 
      par[0] = hzOTTFI;
      par[1] = tetOTTFB;
      par[2] = phiOTTFB;
      par[3] = h1OTTFB;
      par[4] = bl1OTTFB;
      par[5] = tl1OTTFB;
      par[6] = alp1OTTFB;
      par[7] = h2OTTFB;
      par[8] = bl2OTTFB;
      par[9] = tl2OTTFB;
      par[10] = alp2OTTFB;

      gMC->Gsvolu("OHSC","TRAP",idSolder,par,11);  

  // OutTopCuboidFrame Solder 
      par[0] = hxOutHFrame;
      par[1] = hyOutHFrame;
      par[2] = hzOTTFI;
      gMC->Gsvolu("OHIA","BOX",idInox,par,3);

  // OutTopTrapFrameA Inox 
      par[0] = hzOTTFI;
      par[1] = tetOTTFA;
      par[2] = phiOTTFA;
      par[3] = h1OTTFA;
      par[4] = bl1OTTFA;
      par[5] = tl1OTTFA;
      par[6] = alp1OTTFA;
      par[7] = h2OTTFA;
      par[8] = bl2OTTFA;
      par[9] = tl2OTTFA;
      par[10] = alp2OTTFA;

      gMC->Gsvolu("OHIB","TRAP",idInox,par,11); 

  // OutTopTrapFrameB Inox 
      par[0] = hzOTTFI;
      par[1] = tetOTTFB;
      par[2] = phiOTTFB;
      par[3] = h1OTTFB;
      par[4] = bl1OTTFB;
      par[5] = tl1OTTFB;
      par[6] = alp1OTTFB;
      par[7] = h2OTTFB;
      par[8] = bl2OTTFB;
      par[9] = tl2OTTFB;
      par[10] = alp2OTTFB;

      gMC->Gsvolu("OHIC","TRAP",idInox,par,11);

  /////////////////////////////////

  // OutEdgeTrapFrame Epoxy = (4 trapezes)*2 copies*3 layers (Epoxy/Inox/Solder)
  // TRAPEZE 1
      par[0] = hzOETF;
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

      gMC->Gsvolu("EDE1","TRAP",idFrameEpoxy,par,11);  

  // TRAPEZE 2
      par[4] = bl1OETF2;
      par[5] = tl1OETF2;
      par[6] = alp1OETF2;

      par[8] = bl2OETF2;
      par[9] = tl2OETF2;
      par[10] = alp2OETF2; 
     gMC->Gsvolu("EDE2","TRAP",idFrameEpoxy,par,11);  

  // TRAPEZE 3
      par[4] = bl1OETF3;
      par[5] = tl1OETF3;
      par[6] = alp1OETF3;

      par[8] = bl2OETF3;
      par[9] = tl2OETF3;
      par[10] = alp2OETF3; 
     gMC->Gsvolu("EDE3","TRAP",idFrameEpoxy,par,11);  

  // TRAPEZE 4
      par[4] = bl1OETF4;
      par[5] = tl1OETF4;
      par[6] = alp1OETF4;

      par[8] = bl2OETF4;
      par[9] = tl2OETF4;
      par[10] = alp2OETF4; 
     gMC->Gsvolu("EDE4","TRAP",idFrameEpoxy,par,11);          

  ////////////////////////////////

    // TRAPEZE 1
      par[0] = hzOuterFrameInox;
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

      gMC->Gsvolu("EDI1","TRAP",idInox,par,11);  

  // TRAPEZE 2
      par[4] = bl1OETF2;
      par[5] = tl1OETF2;
      par[6] = alp1OETF2;

      par[8] = bl2OETF2;
      par[9] = tl2OETF2;
      par[10] = alp2OETF2; 
     gMC->Gsvolu("EDI2","TRAP",idInox,par,11);  

  // TRAPEZE 3
      par[4] = bl1OETF3;
      par[5] = tl1OETF3;
      par[6] = alp1OETF3;

      par[8] = bl2OETF3;
      par[9] = tl2OETF3;
      par[10] = alp2OETF3; 
     gMC->Gsvolu("EDI3","TRAP",idInox,par,11);  

  // TRAPEZE 4
      par[4] = bl1OETF4;
      par[5] = tl1OETF4;
      par[6] = alp1OETF4;

      par[8] = bl2OETF4;
      par[9] = tl2OETF4;
      par[10] = alp2OETF4; 
     gMC->Gsvolu("EDI4","TRAP",idInox,par,11);          


  ////////////////////////////////

    // TRAPEZE 1
      par[0] = hzOuterFrameSolder;
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

      gMC->Gsvolu("EDS1","TRAP",idSolder,par,11);

  // TRAPEZE 2
      par[4] = bl1OETF2;
      par[5] = tl1OETF2;
      par[6] = alp1OETF2;

      par[8] = bl2OETF2;
      par[9] = tl2OETF2;
      par[10] = alp2OETF2; 
     gMC->Gsvolu("EDS2","TRAP",idSolder,par,11);  

  // TRAPEZE 3
      par[4] = bl1OETF3;
      par[5] = tl1OETF3;
      par[6] = alp1OETF3;

      par[8] = bl2OETF3;
      par[9] = tl2OETF3;
      par[10] = alp2OETF3; 
     gMC->Gsvolu("EDS3","TRAP",idSolder,par,11);  

  // TRAPEZE 4
      par[4] = bl1OETF4;
      par[5] = tl1OETF4;
      par[6] = alp1OETF4;

      par[8] = bl2OETF4;
      par[9] = tl2OETF4;
      par[10] = alp2OETF4; 
     gMC->Gsvolu("EDS4","TRAP",idSolder,par,11);          

  //////////////////////////////////

  // OutCornerTrapFrame  
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

      gMC->Gsvolu("TCOR","TRAP",idFrameEpoxy,par,11);

  // MIRE
      par[0] = MIREInRad;
      par[1] = MIREOutRad;
      par[2] = MIRELen;    

      gMC->Gsvolu("MIRE","TUBE",idFrameEpoxy,par,3);   
  }
             
// __________________Place volumes in the quadrant ____________ 

    const Float_t BeamOX = 2.6;
    const Float_t BeamOY = 2.6;  

 // Coordinates of the frame corner wrt the beam axis defined at (0,0)
    Float_t QuadOX = 0.;
    Float_t QuadOY = 0.;
    
// InVFrame  
    posX = hxInVFrame;
    posY = 2.0*hyInHFrame+2.*hyH1mm+IAF+hyInVFrame;        
    posZ = 0.;   
    gMC->Gspos("IVEF",1,QuadrantName(chamber),posX, posY, posZ, 0, "ONLY"); 

// InHFrame
    posX = 2.0*hxInVFrame+2.*hxV1mm+IAF+hxInHFrame;
    posY = hyInHFrame;
    posZ = 0.;       
    gMC->Gspos("IHEF",1,QuadrantName(chamber),posX, posY, posZ, 0, "ONLY"); 
    
// OutVFrame
    posX = 2.*hxInVFrame+IAF+2.*hxInHFrame-hxOutVFrame+2.*hxV1mm;
    posY = 2.*hyInHFrame+hyOutVFrame;    
    posZ = 0.;              
    gMC->Gspos("OVEF",1,QuadrantName(chamber),posX, posY, posZ, 0, "ONLY"); 
//    cout << " Outer vertical frame at " << posX << " " << posY << " "
//           << posZ << " and half length " << hyOutVFrame << endl;           
    const Float_t TOPY = posY+hyOutVFrame;
    const Float_t OUTX = posX;

//Flat 1mm vertical section
    posX = 2.0*hxInVFrame+hxV1mm;
    posY = 2.0*hyInHFrame+2.*hyH1mm+IAF+hyV1mm;
    posZ = 0.;
    gMC->Gspos("FVMM",1,QuadrantName(chamber),posX, posY, posZ,0, "ONLY"); 
 
// Flat 1mm horizontal section
    posX = 2.0*hxInVFrame+2.*hxV1mm+IAF+hxH1mm;
    posY = 2.0*hyInHFrame+hyH1mm;
    posZ = 0.;
    gMC->Gspos("FHMM",1,QuadrantName(chamber),posX, posY, posZ,0, "ONLY"); 
        
// InArcFrame 
    posX = 2.0*hxInVFrame+2.*hxV1mm;
    posY = 2.0*hyInHFrame+2.*hyH1mm;
    posZ = 0.;    
    gMC->Gspos("IAF1",1,QuadrantName(chamber),posX, posY, posZ,0, "ONLY"); 

    
// ScrewsInFrame  

// Only place screws that are inside the sensitive volume.

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
     posX = BeamOX + scruX[i];
     posY = BeamOY + scruY[i];
     posZ = 0.;   
     gMC->Gspos("SCRH",i+1,QuadrantName(chamber),posX-QuadOX+0.1, posY-QuadOY+0.1, posZ-hzInHFrame-SCRUHLE, 0, "ONLY");      
     gMC->Gspos("SCRM",i+1,QuadrantName(chamber),posX-QuadOX+0.1, posY-QuadOY+0.1, posZ, 0, "ONLY");
     gMC->Gspos("SCRN",i+1,QuadrantName(chamber),posX-QuadOX+0.1, posY-QuadOY+0.1, posZ+hzInHFrame+SCRUNLE, 0, "ONLY"); 
     }
     // special screw coordinates
     scruX[63] = 16.3;  
     scruY[63] = -2.23; 
     posX = BeamOX + scruX[63];
     posY = BeamOY + scruY[63];
     posZ = 0.;            
     gMC->Gspos("SCRH",64,QuadrantName(chamber),posX-QuadOX+0.1, posY-QuadOY+0.1, posZ-hzInHFrame-SCRUHLE, 0, "ONLY");
     gMC->Gspos("SCRM",64,QuadrantName(chamber),posX-QuadOX+0.1, posY-QuadOY+0.1, posZ, 0, "ONLY"); 
     gMC->Gspos("SCRN",64,QuadrantName(chamber),posX-QuadOX+0.1, posY-QuadOY+0.1, posZ+hzInHFrame+SCRUNLE, 0, "ONLY");  
     
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
    posX = BeamOX + scruX[i+LastScrew-1];
    posY = BeamOY + scruY[i+LastScrew-1];
    posZ = 0.;       
    gMC->Gspos("SCRH",i+LastScrew,QuadrantName(chamber),posX-QuadOX+0.1, posY-QuadOY+0.1, posZ-hzInHFrame-SCRUHLE, 0, "ONLY");     
    gMC->Gspos("SCRM",i+LastScrew,QuadrantName(chamber),posX-QuadOX+0.1, posY-QuadOY+0.1, posZ, 0, "ONLY"); 
    gMC->Gspos("SCRN",i+LastScrew,QuadrantName(chamber),posX-QuadOX+0.1, posY-QuadOY+0.1, posZ+hzInHFrame+SCRUNLE, 0, "ONLY");
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
    posX = BeamOX + scruX[i+FirstScrew-1];
    posY = BeamOY + scruY[i+FirstScrew-1];
    posZ = 0.;   
    gMC->Gspos("SCRH",i+FirstScrew,QuadrantName(chamber),posX-QuadOX+0.1, posY-QuadOY+0.1, posZ-hzInHFrame-SCRUHLE, 0, "ONLY");     
    gMC->Gspos("SCRM",i+FirstScrew,QuadrantName(chamber),posX-QuadOX+0.1, posY-QuadOY+0.1, posZ, 0, "ONLY"); 
    gMC->Gspos("SCRN",i+FirstScrew,QuadrantName(chamber),posX-QuadOX+0.1, posY-QuadOY+0.1, posZ+hzInHFrame+SCRUNLE, 0, "ONLY"); 
    }
      
// Inner Arc of Frame, screw positions and numbers-1
   scruX[62] = 16.009; scruY[62]  = 1.401;
   scruX[61] = 14.564; scruY[61]  = 6.791;
   scruX[60] = 11.363; scruY[60]  = 11.363;
   scruX[59] = 6.791 ; scruY[59]  = 14.564;
   scruX[58] = 1.401 ; scruY[58]  = 16.009;
    
    for (Int_t i = 0;i<5;i++){
    posX = BeamOX + scruX[i+58];
    posY = BeamOY + scruY[i+58];
    posZ = 0.;   
    gMC->Gspos("SCRH",i+58+1,QuadrantName(chamber),posX-QuadOX+0.1, posY-QuadOY+0.1, posZ-hzInHFrame-SCRUHLE, 0, "ONLY");    
    gMC->Gspos("SCRM",i+58+1,QuadrantName(chamber),posX-QuadOX+0.1, posY-QuadOY+0.1, posZ, 0, "ONLY");
    gMC->Gspos("SCRN",i+58+1,QuadrantName(chamber),posX-QuadOX+0.1, posY-QuadOY+0.1, posZ+hzInHFrame+SCRUNLE, 0, "ONLY");
    }  

// OutTopTrapFrame   
    posX = hxOutHFrame;
    posY = 2.*hyInHFrame+IAF+2.*hyInVFrame+hyOutHFrame+2.*hyH1mm;
    posZ = 0.;

// place 3 layers of cuboids    
    posZ = posZ-(hzOuterFrameSolder+hzOuterFrameInox);
    gMC->Gspos("OHEA",1,QuadrantName(chamber),posX, posY, posZ,0,"ONLY"); 
    posZ = posZ+hzOuterFrameInox+hzOuterFrameSolder;
    gMC->Gspos("OHSA",1,QuadrantName(chamber),posX, posY, posZ,0,"ONLY");
    posZ = posZ+hzOuterFrameSolder+hzOuterFrameInox;
    gMC->Gspos("OHIA",1,QuadrantName(chamber),posX, posY, posZ,0,"ONLY");
    
// place 3 layers of trapezoid A
    posX = 34.1663+2.6;
    posY = 92.2946+2.6;
    posZ = 0.; 
    
    posZ = posZ-(hzOuterFrameSolder+hzOuterFrameInox); 
    gMC->Gspos("OHEB",1,QuadrantName(chamber),posX, posY, posZ, idrotm[1102],"ONLY"); 
    posZ = posZ+hzOuterFrameInox+hzOuterFrameSolder;
    gMC->Gspos("OHSB",1,QuadrantName(chamber),posX, posY, posZ, idrotm[1102],"ONLY");
    posZ = posZ+hzOuterFrameSolder+hzOuterFrameInox;
    gMC->Gspos("OHIB",1,QuadrantName(chamber),posX, posY, posZ, idrotm[1102],"ONLY");
    
// place 3 layers of trapezoid B
    posX = 43.8379+2.6;
    posY = 90.1946+2.6;
    posZ = 0.; 
    
    posZ = posZ-(hzOuterFrameSolder+hzOuterFrameInox); 
    gMC->Gspos("OHEC",1,QuadrantName(chamber),posX, posY, posZ, idrotm[1102],"ONLY"); 
    posZ = posZ+hzOuterFrameInox+hzOuterFrameSolder;
    gMC->Gspos("OHSC",1,QuadrantName(chamber),posX, posY, posZ, idrotm[1102],"ONLY");
    posZ = posZ+hzOuterFrameSolder+hzOuterFrameInox;
    gMC->Gspos("OHIC",1,QuadrantName(chamber),posX, posY, posZ, idrotm[1102],"ONLY");    
              
///////////////////////////////////////////////////             
        
// OutEdgeTrapFrame  
 
    const Float_t refY = 70.13685;
    const Float_t refX = 71.43685;
    
    posX = refX;
    posY = refY;
    posZ = 0.; 
    
    Float_t XCenter[8]; 
    Float_t YCenter[8];
    
    XCenter[0] = 72.7099 + 2.6;
    XCenter[1] = 77.5787 + 2.6; 
    XCenter[2] = 82.2732 + 2.6;
    XCenter[3] = 86.7882 + 2.6; 
    
    YCenter[0] = 67.6691 + 2.6;
    YCenter[1] = 62.4564 + 2.6; 
    YCenter[2] = 57.0693 + 2.6;
    YCenter[3] = 51.5027 + 2.6; 
      
    XCenter[4] = 67.6691 + 2.6;
    XCenter[5] = 62.4564 + 2.6; 
    XCenter[6] = 57.0693 + 2.6;
    XCenter[7] = 51.5027 + 2.6; 
    
    YCenter[4] = 72.7099 + 2.6;
    YCenter[5] = 77.5787 + 2.6; 
    YCenter[6] = 82.2732 + 2.6;
    YCenter[7] = 86.7882 + 2.6; 
      
    posZ = posZ-(hzOuterFrameSolder+hzOuterFrameInox);     
    gMC->Gspos("EDE1",1,QuadrantName(chamber), XCenter[0], YCenter[0], posZ, idrotm[1104],"ONLY");
    gMC->Gspos("EDE1",2,QuadrantName(chamber), XCenter[4], YCenter[4],posZ, idrotm[1105],"ONLY");


    gMC->Gspos("EDE2",1,QuadrantName(chamber), XCenter[1], YCenter[1], posZ, idrotm[1104],"ONLY");   
    gMC->Gspos("EDE2",2,QuadrantName(chamber), XCenter[5], YCenter[5], posZ, idrotm[1105],"ONLY");

    gMC->Gspos("EDE3",1,QuadrantName(chamber), XCenter[2], YCenter[2], posZ, idrotm[1104],"ONLY");
    gMC->Gspos("EDE3",2,QuadrantName(chamber), XCenter[6], YCenter[6], posZ, idrotm[1105],"ONLY");
    
    
    gMC->Gspos("EDE4",1,QuadrantName(chamber), XCenter[3], YCenter[3], posZ, idrotm[1104],"ONLY");
    gMC->Gspos("EDE4",2,QuadrantName(chamber), XCenter[7], YCenter[7], posZ, idrotm[1105],"ONLY");          


     
    posZ = posZ+hzOuterFrameEpoxy+hzOuterFrameInox;
     
    gMC->Gspos("EDI1",1,QuadrantName(chamber), XCenter[0], YCenter[0], posZ, idrotm[1104],"ONLY");
    gMC->Gspos("EDI1",2,QuadrantName(chamber), XCenter[4], YCenter[4], posZ, idrotm[1105],"ONLY");
    
    gMC->Gspos("EDI2",1,QuadrantName(chamber), XCenter[1], YCenter[1], posZ, idrotm[1104],"ONLY");   
    gMC->Gspos("EDI2",2,QuadrantName(chamber), XCenter[5], YCenter[5], posZ, idrotm[1105],"ONLY");

    gMC->Gspos("EDI3",1,QuadrantName(chamber), XCenter[2], YCenter[2], posZ, idrotm[1104],"ONLY");
    gMC->Gspos("EDI3",2,QuadrantName(chamber), XCenter[6], YCenter[6], posZ, idrotm[1105],"ONLY");
       
    gMC->Gspos("EDI4",1,QuadrantName(chamber), XCenter[3], YCenter[3], posZ, idrotm[1104],"ONLY");
    gMC->Gspos("EDI4",2,QuadrantName(chamber), XCenter[7], YCenter[7], posZ, idrotm[1105],"ONLY");  
           

    posZ = posZ+hzOuterFrameInox+hzOuterFrameSolder;
       
    gMC->Gspos("EDS1",1,QuadrantName(chamber), XCenter[0], YCenter[0], posZ, idrotm[1104],"ONLY");
    gMC->Gspos("EDS1",2,QuadrantName(chamber), XCenter[4], YCenter[4], posZ, idrotm[1105],"ONLY");
    
    gMC->Gspos("EDS2",1,QuadrantName(chamber), XCenter[1], YCenter[1], posZ, idrotm[1104],"ONLY");   
    gMC->Gspos("EDS2",2,QuadrantName(chamber), XCenter[5], YCenter[5], posZ, idrotm[1105],"ONLY");

    gMC->Gspos("EDS3",1,QuadrantName(chamber), XCenter[2], YCenter[2], posZ, idrotm[1104],"ONLY");
    gMC->Gspos("EDS3",2,QuadrantName(chamber), XCenter[6], YCenter[6], posZ, idrotm[1105],"ONLY");
       
    gMC->Gspos("EDS4",1,QuadrantName(chamber), XCenter[3], YCenter[3], posZ, idrotm[1104],"ONLY");
    gMC->Gspos("EDS4",2,QuadrantName(chamber), XCenter[7], YCenter[7], posZ, idrotm[1105],"ONLY");  
    
  
// OutCornerTrapFrame  
    posX = OUTX;
    posY = TOPY+((bl1OCTF+tl1OCTF)/2.);
    posZ = 0.;     
    gMC->Gspos("TCOR",1,QuadrantName(chamber),posX, posY, posZ, idrotm[1102],"ONLY"); 
    
// Mire placement
    posX = OUTX+hxOutVFrame+1.3;
    posY = 8.13+2.6;
    posZ = 0.;
   gMC->Gspos("MIRE",1,QuadrantName(chamber),posX, posY, posZ, 0,"ONLY");
    
}

//___________________________________________
void AliMUONv2::CreateHole()
{
// Create all the element inside a foam hole
// --
  Int_t* idtmed = fIdtmed->GetArray()-1099;
  Int_t idAir  = idtmed[1100]; // medium 1
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

//___________________________________________
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
  par[2]=totalHzDaughter();
  gMC->Gsvolu(fgkDaughterName,"BOX",idAir,par,3);
  
  par[0]=fgkHxBergPlastic;
  par[1]=fgkHyBergPlastic;
  par[2]=fgkHzBergPlastic;
  gMC->Gsvolu("BRGP","BOX",idPlastic,par,3);
  posX=0.;
  posY=0.;
  posZ = -totalHzDaughter() + fgkHzBergPlastic;
  gMC->Gspos("BRGP",1,fgkDaughterName,posX,posY,posZ,0,"ONLY");

  par[0]=fgkHxBergCopper;
  par[1]=fgkHyBergCopper;
  par[2]=fgkHzBergCopper;
  gMC->Gsvolu("BRGC","BOX",idCopper,par,3);
  posX=0.;
  posY=0.;
  posZ=0.;
  gMC->Gspos("BRGC",1,"BRGC",posX,posY,posZ,0,"ONLY");

  par[0]=fgkHxDaughter;
  par[1]=fgkHyDaughter;
  par[2]=fgkHzDaughter;
  gMC->Gsvolu("DGHT","BOX",idCopper,par,3);
  posX=0.;
  posY=0.;
  posZ = -totalHzDaughter() + 2.*fgkHzBergPlastic + fgkHzDaughter;
  gMC->Gspos("DGHT",1,fgkDaughterName,posX,posY,posZ,0,"ONLY");
}

//___________________________________________
void AliMUONv2::CreatePlaneBox(const char* name,const  TVector2& dimensions)
{
// create all the elements in the copper plane
// --
  Int_t* idtmed = fIdtmed->GetArray()-1099;
  Int_t idAir  = idtmed[1100]; // medium 1
  Int_t idCopper  = idtmed[1109]; // medium 10 = copper
  Int_t idFoam = idtmed[1115]; // medium 16 = Foam
  Int_t idFR4 =idtmed[1114]; // medium 15 = FR4

  GReal_t par[3];
  GReal_t posX,posY,posZ;

  // mother volume
  par[0] = dimensions.X();
  par[1] = dimensions.Y();
  par[2] = totalHzPlane();
  gMC->Gsvolu(name,"BOX",idAir,par,3);
  
  // pad plane
  char* planeName = strdup(name);
  planeName[3]='P';
  par[0] = dimensions.X();
  par[1] = dimensions.Y();
  par[2] = fgkHzPadPlane;
  gMC->Gsvolu(planeName,"BOX",idCopper,par,3);
  posX=0.;
  posY=0.;
  posZ = -totalHzPlane()+fgkHzPadPlane;
  gMC->Gspos(planeName,1,name,posX,posY,posZ,0,"ONLY");
  
  //foam layer
  char* foamName = strdup(name);
  foamName[3]=fgkFoamLayerSuffix;
  par[0] = dimensions.X();
  par[1] = dimensions.Y();
  par[2] = fgkHzFoam;
  gMC->Gsvolu(foamName,"BOX",idFoam,par,3);
  posX=0.;
  posY=0.;
  posZ = -totalHzPlane()+2.*fgkHzPadPlane+fgkHzFoam;
  gMC->Gspos(foamName,1,name,posX,posY,posZ,0,"ONLY");

  // mechanical plane FR4 layer
  char* fr4Name = strdup(name);
  fr4Name[3]='R';
  par[0] = dimensions.X();
  par[1] = dimensions.Y();
  par[2] = fgkHzFR4;
  gMC->Gsvolu(fr4Name,"BOX",idFR4,par,3);
  posX=0.;
  posY=0.;
  posZ = -totalHzPlane()+2.*fgkHzPadPlane+2.*fgkHzFoam+fgkHzFR4;
  gMC->Gspos(fr4Name,1,name,posX,posY,posZ,0,"ONLY");
}

//___________________________________________
void AliMUONv2::CreatePlaneSegment(const char* name,const  TVector2& dimensions
      	      	      	      	  ,Int_t nofHoles)
{
// Create a segment of a plane (this includes a copper layer, foam, hole, 
// and kapton as well as the mother board.)
// --
  static Int_t holeNum=1;
    
  GReal_t posX,posY,posZ;
  
  CreatePlaneBox(name,dimensions);
  char holeName[5];
  strcpy(holeName,name);
  holeName[3]=fgkFoamLayerSuffix;
  // <dname> is a motif on the pad plane
  char* dname = strdup(name);
  dname[3]='D';
  gMC->Gsdvn(dname,holeName,nofHoles,1);
  
  posX=0.;
  posY=0.;
  posZ= fgkHzPadPlane;
  gMC->Gspos(fgkHoleName,holeNum++,dname,posX,posY,posZ,0,"ONLY");
}

//___________________________________________
void AliMUONv2::CreateDaughterSegment(const char* name,const  TVector2& dimensions,
                              Int_t nofHoles)
{
// Create a segment of a daughter board layer
// --
  static Int_t holeNum=1;
  Int_t* idtmed = fIdtmed->GetArray()-1099;
  Int_t idAir  = idtmed[1100]; // medium 1

  GReal_t par[3];
  GReal_t posX,posY,posZ;
  
  par[0] = dimensions.X();
  par[1] = dimensions.Y();
  par[2] = totalHzDaughter();
  gMC->Gsvolu(name,"BOX",idAir,par,3);
  
  // <dname> is a motif on pad plane
  char* dname = strdup(name);
  dname[3]='D';
  gMC->Gsdvn(dname,name,nofHoles,1);
  
  posX=0.;
  posY=0.;
  posZ=0.;
  gMC->Gspos(fgkDaughterName,holeNum++,dname,posX,posY,posZ,0,"ONLY");
}

//___________________________________________
void AliMUONv2::PlaceSector(MSector* sector,TSpecialMap specialMap
                           ,const TVector3& where,Int_t chamber)
{
// Place all the segments in the mother volume, in the position defined
// by the sector's data.
// --
  static Int_t segNum=1;
  
  GReal_t posX,posY,posZ;
  
  vector<int> already_done;
  Int_t rotNum;
  AliMatrix(rotNum,  90.,90.,90,180.,0.,0.);
  
  for (Int_t irow=0;irow<sector->GetNofRows();irow++){
    MRow* row = sector->GetRow(irow);
    for (Int_t iseg=0;iseg<row->GetNofRowSegments();iseg++){
      MVRowSegment* seg = row->GetRowSegment(iseg);
      char* segName;
      
      TSpecialMap::iterator iter 
        = specialMap.find(seg->GetMotifPositionId(0));
      if ( iter == specialMap.end()){
        segName = strdup(Form("%.3dM",segNum*2));
        CreatePlaneSegment(segName,seg->Dimensions()/10.,seg->GetNofMotifs());     
        posX = where.X()+seg->Position().X()/10.-fgkOffsetX/2.;
        posY = where.Y()+seg->Position().Y()/10.-fgkOffsetY/2.;
        posZ = where.Z()-totalHz()+totalHzPlane();
        gMC->Gspos(segName,1,QuadrantName(chamber),posX,posY,posZ,0,"ONLY");

        segName = strdup(Form("%.3dM",segNum*2+1));
        CreateDaughterSegment(segName,seg->Dimensions()/10.,seg->GetNofMotifs());     
        posX = where.X()+seg->Position().X()/10.+fgkOffsetX/2.;
        posY = where.Y()+seg->Position().Y()/10.-+fgkOffsetY/2.;
        posZ = where.Z()-totalHz()+2.*totalHzPlane()+totalHzDaughter();
        gMC->Gspos(segName,1,QuadrantName(chamber),posX,posY,posZ,0,"ONLY");
      	segNum++;
      } else {
      	for (Int_t motifNum=0;motifNum<seg->GetNofMotifs();motifNum++) {
	  Int_t motifPosId = seg->GetMotifPositionId(motifNum);
	  
	  if (find(already_done.begin(),already_done.end(),motifPosId)
	      != already_done.end()) continue;
	  
            AliMUONSt1SpecialMotif spMot = specialMap[motifPosId];

	  MMotifPosition* motifPos = 
	    sector->GetMotifMap()->FindMotifPosition(motifPosId);
	  MMotifType* mType = motifPos->GetMotif()->GetMotifType();
	  for (Int_t line=0;line<mType->GetNofPadsY();line++) {
	    Int_t iMin=-1,iMax=-1,col;
	    for (col=0;col<mType->GetNofPadsX();col++) {
	      if (mType->HasPad(MIntPair(col,line))) {
		iMin=col;
		break;
	      }
	    }
	    for (col=mType->GetNofPadsX()-1;col>=0;col--) {
	      if (mType->HasPad(MIntPair(col,line))) {
		iMax=col;
		break;
	      }
	    }
	    if ( (iMin>=0) && (iMax>=0) ) {
	      TVector2 dim =  
	       motifPos->GetMotif()->GetPadDimensions(MIntPair(iMin,line))/10.;
	      TVector2 lineDim(dim.X()*(iMax-iMin+1),dim.Y());
	      char* boxName = strdup(Form("%.3dM",segNum*2));
	      CreatePlaneBox(boxName,lineDim);
      	      TVector2 posLLline = TVector2(dim.X()*2.*iMin,dim.Y()*2.*line)
	      	      	      	    +motifPos->Position()/10.-motifPos->Dimensions()/10.;
      	      TVector2 centerLine = posLLline + lineDim;
	      posX = where.X()+centerLine.X()-fgkOffsetX/2.;
	      posY = where.Y()+centerLine.Y()-fgkOffsetY/2.;
	      posZ = where.Z()-totalHz()+totalHzPlane();
	      gMC->Gspos(boxName,1,QuadrantName(chamber),posX,posY,posZ,0,"ONLY");
	    }
	    
	    segNum++;
	  }
	  
	  Int_t rot = ( spMot.GetRotAngle()<0.1 ) ? 0:rotNum;

	  posX = where.X()+motifPos->Position().X()/10.-fgkOffsetX/2.+spMot.GetDelta().X();
	  posY = where.Y()+motifPos->Position().Y()/10.-fgkOffsetY/2.+spMot.GetDelta().Y();
      	  posZ = where.Z()-fgkHzPadPlane;
      	  gMC->Gspos(fgkHoleName,motifPosId,QuadrantName(chamber),posX,posY,posZ,rot,"ONLY");
	  
	  posX = where.X()+motifPos->Position().X()/10.+fgkOffsetX/2.-fgkDeltaFilleEtamX+spMot.GetDelta().X();
	  posY = where.Y()+motifPos->Position().Y()/10.-+fgkOffsetY/2.-fgkDeltaFilleEtamY+spMot.GetDelta().Y();
      	  posZ = where.Z()-totalHz()+2.*totalHzPlane()+totalHzDaughter();
      	  gMC->Gspos(fgkDaughterName,motifPosId,QuadrantName(chamber),posX,posY,posZ,rot,"ONLY");

      	  already_done.push_back(motifPosId);		  
	}
      }
    }
  } 
}
