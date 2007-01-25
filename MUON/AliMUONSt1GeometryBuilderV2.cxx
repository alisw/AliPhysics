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

// $Id$
//
// Authors: David Guez, Ivana Hrivnacova, Marion MacCormick; IPN Orsay
//
// Class AliMUONSt1GeometryBuilderV2
// ---------------------------------
// MUON Station1 detailed geometry construction class.
// (Originally defined in AliMUONv2.cxx - now removed.)
// Included in AliRoot 2004/01/23

#include "AliMUONSt1GeometryBuilderV2.h"
#include "AliMUONSt1SpecialMotif.h"
#include "AliMUON.h"
#include "AliMUONConstants.h"
#include "AliMUONGeometryModule.h"
#include "AliMUONGeometryEnvelopeStore.h"

#include "AliMpContainers.h"
#include "AliMpConstants.h"
#include "AliMpFiles.h"
#include "AliMpSectorReader.h"
#include "AliMpSector.h"
#include "AliMpRow.h"
#include "AliMpVRowSegment.h"
#include "AliMpMotifMap.h"
#include "AliMpMotifPosition.h"
#include "AliMpPlaneType.h"

#include "AliRun.h"
#include "AliMagF.h"
#include "AliLog.h"

#include <TVector2.h>
#include <TVector3.h>
#include <TGeoMatrix.h>
#include <TClonesArray.h>
#include <Riostream.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoTube.h>
#include <TGeoCompositeShape.h>

#ifdef WITH_STL
  #include <vector>
#endif

#ifdef WITH_ROOT
  #include "TArrayI.h"
#endif

/// \cond CLASSIMP
ClassImp(AliMUONSt1GeometryBuilderV2)
/// \endcond

// Thickness Constants
const GReal_t AliMUONSt1GeometryBuilderV2::fgkHzPadPlane=0.0148/2.;     //Pad plane
const GReal_t AliMUONSt1GeometryBuilderV2::fgkHzFoam = 2.503/2.;        //Foam of mechanicalplane
const GReal_t AliMUONSt1GeometryBuilderV2::fgkHzFR4 = 0.062/2.;         //FR4 of mechanical plane
const GReal_t AliMUONSt1GeometryBuilderV2::fgkHzSnPb = 0.0091/2.;       //Pad/Kapton connection (66 pt)
const GReal_t AliMUONSt1GeometryBuilderV2::fgkHzKapton = 0.0122/2.;     //Kapton
const GReal_t AliMUONSt1GeometryBuilderV2::fgkHzBergPlastic = 0.3062/2.;//Berg connector
const GReal_t AliMUONSt1GeometryBuilderV2::fgkHzBergCopper = 0.1882/2.; //Berg connector
const GReal_t AliMUONSt1GeometryBuilderV2::fgkHzDaughter = 0.0156/2.;   //Daughter board
const GReal_t AliMUONSt1GeometryBuilderV2::fgkHzGas = 0.42/2.;          //Gas thickness

// Quadrant Mother volume - TUBS1 - Middle layer of model
const GReal_t AliMUONSt1GeometryBuilderV2::fgkMotherIR1 = 18.3;
const GReal_t AliMUONSt1GeometryBuilderV2::fgkMotherOR1 = 105.673;   
const GReal_t AliMUONSt1GeometryBuilderV2::fgkMotherThick1 = 6.5/2;  
const GReal_t AliMUONSt1GeometryBuilderV2::fgkMotherPhiL1 = 0.; 
const GReal_t AliMUONSt1GeometryBuilderV2::fgkMotherPhiU1 = 90.;

// Quadrant Mother volume - TUBS2 - near and far layers of model
const GReal_t AliMUONSt1GeometryBuilderV2::fgkMotherIR2 = 20.7;   
const GReal_t AliMUONSt1GeometryBuilderV2::fgkMotherOR2 = 100.073;   
const GReal_t AliMUONSt1GeometryBuilderV2::fgkMotherThick2 = 3.0/2; 
const GReal_t AliMUONSt1GeometryBuilderV2::fgkMotherPhiL2 = 0.; 
const GReal_t AliMUONSt1GeometryBuilderV2::fgkMotherPhiU2 = 90.;

// Sensitive copper pads, foam layer, PCB and electronics model parameters
const GReal_t AliMUONSt1GeometryBuilderV2::fgkHxHole=1.5/2.;
const GReal_t AliMUONSt1GeometryBuilderV2::fgkHyHole=6./2.;
const GReal_t AliMUONSt1GeometryBuilderV2::fgkHxBergPlastic=0.74/2.;
const GReal_t AliMUONSt1GeometryBuilderV2::fgkHyBergPlastic=5.09/2.;
const GReal_t AliMUONSt1GeometryBuilderV2::fgkHxBergCopper=0.25/2.;
const GReal_t AliMUONSt1GeometryBuilderV2::fgkHyBergCopper=3.6/2.;
const GReal_t AliMUONSt1GeometryBuilderV2::fgkHxKapton=0.8/2.;
const GReal_t AliMUONSt1GeometryBuilderV2::fgkHyKapton=5.7/2.;
const GReal_t AliMUONSt1GeometryBuilderV2::fgkHxDaughter=2.3/2.;
const GReal_t AliMUONSt1GeometryBuilderV2::fgkHyDaughter=6.3/2.;
const GReal_t AliMUONSt1GeometryBuilderV2::fgkOffsetX=1.46;
const GReal_t AliMUONSt1GeometryBuilderV2::fgkOffsetY=0.71;
const GReal_t AliMUONSt1GeometryBuilderV2::fgkDeltaFilleEtamX=1.46;
const GReal_t AliMUONSt1GeometryBuilderV2::fgkDeltaFilleEtamY=0.051;

const GReal_t AliMUONSt1GeometryBuilderV2::fgkDeltaQuadLHC=2.6;  // LHC Origin wrt Quadrant Origin
const GReal_t AliMUONSt1GeometryBuilderV2::fgkFrameOffset=5.2;
              // Fix (1) of overlap SQN* layers with SQM* ones (was 5.0)
	      
// Pad planes offsets
const GReal_t AliMUONSt1GeometryBuilderV2::fgkPadXOffsetBP =  0.50 - 0.63/2; // = 0.185
const GReal_t AliMUONSt1GeometryBuilderV2::fgkPadYOffsetBP = -0.31 - 0.42/2; // =-0.52

const char* AliMUONSt1GeometryBuilderV2::fgkHoleName="SCHL";      
const char* AliMUONSt1GeometryBuilderV2::fgkDaughterName="SCDB";  
const char* AliMUONSt1GeometryBuilderV2::fgkQuadrantEnvelopeName="SE";
const char* AliMUONSt1GeometryBuilderV2::fgkQuadrantMLayerName="SQM";
const char* AliMUONSt1GeometryBuilderV2::fgkQuadrantNLayerName="SQN";
const char* AliMUONSt1GeometryBuilderV2::fgkQuadrantFLayerName="SQF";
const Int_t AliMUONSt1GeometryBuilderV2::fgkFoamBoxNameOffset=200; 
const Int_t AliMUONSt1GeometryBuilderV2::fgkFR4BoxNameOffset=400; 
const Int_t AliMUONSt1GeometryBuilderV2::fgkDaughterCopyNoOffset=1000;

//______________________________________________________________________________
AliMUONSt1GeometryBuilderV2::AliMUONSt1GeometryBuilderV2(AliMUON* muon)
  : AliMUONVGeometryBuilder(0, 2),
    fMUON(muon)
{
/// Standard constructor

   // set path to mapping data files
   if (! gSystem->Getenv("MINSTALL")) {    
     TString dirPath = gSystem->Getenv("ALICE_ROOT");
     dirPath += "/MUON/mapping"; 
     AliMpFiles::SetTopPath(dirPath);
     gSystem->Setenv("MINSTALL", dirPath.Data());
     //cout << "AliMpFiles top path set to " << dirPath << endl;	  
   }
   //else
   //  cout << gSystem->Getenv("MINSTALL") << endl;	    	  
}
 
//______________________________________________________________________________
AliMUONSt1GeometryBuilderV2::AliMUONSt1GeometryBuilderV2()
  : AliMUONVGeometryBuilder(),
    fMUON(0)
{
/// Default Constructor
}

//______________________________________________________________________________
AliMUONSt1GeometryBuilderV2::~AliMUONSt1GeometryBuilderV2()
{
/// Destructor
}


//
//  Private methods
//

//______________________________________________________________________________
TString 
AliMUONSt1GeometryBuilderV2::QuadrantEnvelopeName(Int_t chamber, Int_t quadrant) const
{ 
/// Generate unique envelope name from chamber Id and quadrant number

  return Form("%s%d", Form("%s%d",fgkQuadrantEnvelopeName,chamber), quadrant); 
}

//______________________________________________________________________________
void AliMUONSt1GeometryBuilderV2::CreateHole()
{
/// Create all the elements found inside a foam hole

  Int_t* idtmed = fMUON->GetIdtmed()->GetArray()-1099;
  Int_t idAir  = idtmed[1100];      // medium 1
  //Int_t idCopper  = idtmed[1109]; // medium 10 = copper 
  Int_t idCopper  = idtmed[1121]; // medium 22 = copper 

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
  gMC->Gsvolu("SKPT", "BOX", idCopper, par, 3);
  posX = 0.;
  posY = 0.;
  posZ = 0.;
  gMC->Gspos("SKPT",1,fgkHoleName, posX, posY, posZ, 0,"ONLY");
}

//______________________________________________________________________________
void AliMUONSt1GeometryBuilderV2::CreateDaughterBoard()
{
/// Create all the elements in a daughter board

  Int_t* idtmed = fMUON->GetIdtmed()->GetArray()-1099;
  Int_t idAir  = idtmed[1100]; // medium 1
  //Int_t idCopper  = idtmed[1109]; // medium 10 = copper
  //Int_t idPlastic  =idtmed[1116]; // medium 17 = Plastic
  Int_t idCopper  = idtmed[1121]; // medium 22 = copper
  Int_t idPlastic  =idtmed[1127]; // medium 28 = Plastic

  GReal_t par[3];
  GReal_t posX,posY,posZ;

  par[0]=fgkHxDaughter;
  par[1]=fgkHyDaughter;
  par[2]=TotalHzDaughter();
  gMC->Gsvolu(fgkDaughterName,"BOX",idAir,par,3);
  
  par[0]=fgkHxBergPlastic;
  par[1]=fgkHyBergPlastic;
  par[2]=fgkHzBergPlastic;
  gMC->Gsvolu("SBGP","BOX",idPlastic,par,3);
  posX=0.;
  posY=0.;
  posZ = -TotalHzDaughter() + fgkHzBergPlastic;
  gMC->Gspos("SBGP",1,fgkDaughterName,posX,posY,posZ,0,"ONLY");

  par[0]=fgkHxBergCopper;
  par[1]=fgkHyBergCopper;
  par[2]=fgkHzBergCopper;
  gMC->Gsvolu("SBGC","BOX",idCopper,par,3);
  posX=0.;
  posY=0.;
  posZ=0.;
  gMC->Gspos("SBGC",1,"SBGP",posX,posY,posZ,0,"ONLY");

  par[0]=fgkHxDaughter;
  par[1]=fgkHyDaughter;
  par[2]=fgkHzDaughter;
  gMC->Gsvolu("SDGH","BOX",idCopper,par,3);
  posX=0.;
  posY=0.;
  posZ = -TotalHzDaughter() + 2.*fgkHzBergPlastic + fgkHzDaughter;
  gMC->Gspos("SDGH",1,fgkDaughterName,posX,posY,posZ,0,"ONLY");
}

//______________________________________________________________________________
void AliMUONSt1GeometryBuilderV2::CreateInnerLayers()
{
/// Create the layer of sensitive volumes with gas
/// and the copper layer.

// Gas Medium
  Int_t* idtmed = fMUON->GetIdtmed()->GetArray()-1099; 
  //Int_t idArCO2  = idtmed[1108];  // medium 9 (ArCO2 80%) 
  //Int_t idCopper  = idtmed[1109]; // medium 10 = copper
  Int_t idArCO2   = idtmed[1124]; // medium 25 (ArCO2 80%) 
  Int_t idCopper  = idtmed[1121]; // medium 22 = copper

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
void AliMUONSt1GeometryBuilderV2::CreateQuadrant(Int_t chamber)
{
/// Create the quadrant (bending and non-bending planes)
/// for the given chamber

  CreateFrame(chamber);

#ifdef WITH_STL
  SpecialMap specialMap;
  specialMap[76] = AliMUONSt1SpecialMotif(TVector2( 0.1, 0.84), 90.);
  specialMap[75] = AliMUONSt1SpecialMotif(TVector2( 0.5, 0.36));
  specialMap[47] = AliMUONSt1SpecialMotif(TVector2(1.01, 0.36));
#endif
  
#ifdef WITH_ROOT
  SpecialMap specialMap;
  specialMap.Add(76, (Long_t) new AliMUONSt1SpecialMotif(TVector2( 0.1, 0.84), 90.));
  specialMap.Add(75, (Long_t) new AliMUONSt1SpecialMotif(TVector2( 0.5, 0.36)));
  specialMap.Add(47, (Long_t) new AliMUONSt1SpecialMotif(TVector2(1.01, 0.36)));
#endif

  AliMpSectorReader reader1(AliMp::kStation1, AliMp::kBendingPlane);
  AliMpSector* sector1 = reader1.BuildSector();

  //Bool_t reflectZ = true;
  Bool_t reflectZ = false;
  //TVector3 where = TVector3(2.5+0.1+0.56+0.001, 2.5+0.1+0.001, 0.);
  TVector3 where = TVector3(fgkDeltaQuadLHC + fgkPadXOffsetBP, 
                            fgkDeltaQuadLHC + fgkPadYOffsetBP, 0.);
  PlaceSector(sector1, specialMap, where, reflectZ, chamber);
  
#ifdef WITH_STL
  specialMap.clear();
  specialMap[76] = AliMUONSt1SpecialMotif(TVector2(1.01,0.59),90.);
  specialMap[75] = AliMUONSt1SpecialMotif(TVector2(1.96, 0.17));
  specialMap[47] = AliMUONSt1SpecialMotif(TVector2(2.18,-0.98));
  specialMap[20] = AliMUONSt1SpecialMotif(TVector2(0.2 ,-0.08));
  specialMap[46] = AliMUONSt1SpecialMotif(TVector2(0.2 , 0.25));
  specialMap[74] = AliMUONSt1SpecialMotif(TVector2(0.28, 0.21));
      // Fix (7) - overlap of SQ42 with MCHL (after moving the whole sector
      // in the true position)   
      // Was: specialMap[47] = AliMUONSt1SpecialMotif(TVector2(1.61,-1.18));
#endif

#ifdef WITH_ROOT
  Int_t nb = AliMpConstants::ManuMask(AliMp::kNonBendingPlane);
  specialMap.Delete();
  specialMap.Add(76 | nb,(Long_t) new AliMUONSt1SpecialMotif(TVector2(1.01,0.59),90.));
  specialMap.Add(75 | nb,(Long_t) new AliMUONSt1SpecialMotif(TVector2(1.96, 0.17)));
  specialMap.Add(47 | nb,(Long_t) new AliMUONSt1SpecialMotif(TVector2(2.18,-0.98)));
  specialMap.Add(20 | nb,(Long_t) new AliMUONSt1SpecialMotif(TVector2(0.2 ,-0.08)));
  specialMap.Add(46 | nb,(Long_t) new AliMUONSt1SpecialMotif(TVector2(0.2 , 0.25)));
  specialMap.Add(74 | nb,(Long_t) new AliMUONSt1SpecialMotif(TVector2(0.28, 0.21)));  
      // Fix (7) - overlap of SQ42 with MCHL (after moving the whole sector
      // in the true position)   
      // Was: specialMap.Add(47,(Long_t) new AliMUONSt1SpecialMotif(TVector2(1.61,-1.18)));
#endif

  AliMpSectorReader reader2(AliMp::kStation1, AliMp::kNonBendingPlane);
  AliMpSector* sector2 = reader2.BuildSector();
  
  //reflectZ = false;
  reflectZ = true;
  TVector2 offset = sector2->Position();
  where = TVector3(where.X()+offset.X(), where.Y()+offset.Y(), 0.); 
      // Add the half-pad shift of the non-bending plane wrt bending plane
      // (The shift is defined in the mapping as sector offset)
      // Fix (4) - was TVector3(where.X()+0.63/2, ... - now it is -0.63/2
  PlaceSector(sector2, specialMap, where, reflectZ, chamber);

#ifdef WITH_ROOT
  specialMap.Delete();
#endif
}

//______________________________________________________________________________
void AliMUONSt1GeometryBuilderV2::CreateFoamBox(
                                        Int_t segNumber,
                                        const  TVector2& dimensions)
{
/// Create all the elements in the copper plane

  Int_t* idtmed = fMUON->GetIdtmed()->GetArray()-1099;
  Int_t idAir  = idtmed[1100]; // medium 1
  //Int_t idFoam = idtmed[1115]; // medium 16 = Foam
  //Int_t idFR4  = idtmed[1114]; // medium 15 = FR4
  Int_t idFoam = idtmed[1125]; // medium 26 = Foam
  Int_t idFR4  = idtmed[1122]; // medium 23 = FR4

  // mother volume
  GReal_t par[3];
  par[0] = dimensions.X();
  par[1] = dimensions.Y();
  par[2] = TotalHzPlane();
  gMC->Gsvolu(PlaneSegmentName(segNumber).Data(),"BOX",idAir,par,3);
  
  // foam layer
  par[0] = dimensions.X();
  par[1] = dimensions.Y();
  par[2] = fgkHzFoam;
  gMC->Gsvolu(FoamBoxName(segNumber).Data(),"BOX",idFoam,par,3);
  GReal_t posX,posY,posZ;
  posX=0.;
  posY=0.;
  posZ = -TotalHzPlane() + fgkHzFoam;
  gMC->Gspos(FoamBoxName(segNumber).Data(),1, 
             PlaneSegmentName(segNumber).Data(),posX,posY,posZ,0,"ONLY");

  // mechanical plane FR4 layer
  par[0] = dimensions.X();
  par[1] = dimensions.Y();
  par[2] = fgkHzFR4;
  gMC->Gsvolu(FR4BoxName(segNumber).Data(),"BOX",idFR4,par,3);
  posX=0.;
  posY=0.;
  posZ = -TotalHzPlane()+ 2.*fgkHzFoam + fgkHzFR4;
  gMC->Gspos(FR4BoxName(segNumber).Data(),1,
             PlaneSegmentName(segNumber).Data(),posX,posY,posZ,0,"ONLY");
}

//______________________________________________________________________________
void AliMUONSt1GeometryBuilderV2::CreatePlaneSegment(Int_t segNumber,
                                    const  TVector2& dimensions,
      	      	      	      	    Int_t nofHoles)
{
/// Create a segment of a plane (this includes a foam layer, 
/// holes in the foam to feed the kaptons through, kapton connectors
/// and the mother board.)
  
  CreateFoamBox(segNumber,dimensions);

  for (Int_t holeNum=0;holeNum<nofHoles;holeNum++) {
    GReal_t posX = ((2.*holeNum+1.)/nofHoles-1.)*dimensions.X();
    GReal_t posY = 0.;
    GReal_t posZ = 0.;
  
    gMC->Gspos(fgkHoleName,holeNum+1,
               FoamBoxName(segNumber).Data(),posX,posY,posZ,0,"ONLY");
  }
}

//______________________________________________________________________________
void AliMUONSt1GeometryBuilderV2::CreateFrame(Int_t chamber)
{
/// Create the non-sensitive elements of the frame for the \a chamber
///
/// Model and notation:                                                     \n
///                                                                         \n
/// The Quadrant volume name starts with SQ                                 \n
/// The volume segments are numbered 00 to XX                               \n
///                                                                         \n
///                              OutTopFrame                                \n
///                               (SQ02-16)                                 \n 
///                              ------------                               \n
///             OutEdgeFrame   /              |                             \n
///             (SQ17-24)     /               |  InVFrame (SQ00-01)         \n 
///                          /                |                             \n
///                          |                |                             \n 
///               OutVFrame  |            _- -                              \n
///               (SQ25-39)  |           |   InArcFrame (SQ42-45)           \n
///                          |           |                                  \n 
///                          -------------                                  \n 
///                        InHFrame (SQ40-41)                               \n 
///                                                                         \n                         
///                                                                         \n
/// 06 February 2003 - Overlapping volumes resolved.                        \n
/// One quarter chamber is comprised of three TUBS volumes: SQMx, SQNx, and SQFx,
/// where SQMx is the Quadrant Middle layer for chamber \a chamber ( posZ in [-3.25,3.25]),
/// SQNx is the Quadrant Near side layer for chamber \a chamber ( posZ in [-6.25,3-.25) ), and
/// SQFx is the Quadrant Far side layer for chamber \a chamber ( posZ in (3.25,6.25] ).

  const Float_t kNearFarLHC=2.4;    // Near and Far TUBS Origin wrt LHC Origin

  // tracking medias
  Int_t* idtmed = fMUON->GetIdtmed()->GetArray()-1099;
  
  Int_t idAir  = idtmed[1100];       // medium 1
  //Int_t idFrameEpoxy = idtmed[1115]; // medium 16 = Frame Epoxy ME730
  //Int_t idInox = idtmed[1116];       // medium 17 Stainless Steel (18%Cr,9%Ni,Fe)
  //Int_t idFR4 = idtmed[1110];        // medium 11 FR4
  //Int_t idCopper = idtmed[1109];     // medium 10 Copper
  //Int_t idAlu = idtmed[1103];        // medium 4 Aluminium
  Int_t idFrameEpoxy = idtmed[1123]; // medium 24 = Frame Epoxy ME730  // was 20 not 16
  Int_t idInox = idtmed[1128];       // medium 29 Stainless Steel (18%Cr,9%Ni,Fe) // was 21 not 17
  Int_t idFR4 = idtmed[1122];        // medium 23 FR4  // was 15 not 11
  Int_t idCopper = idtmed[1121];     // medium 22 Copper
  Int_t idAlu = idtmed[1120];        // medium 21 Aluminium
  
  
// Rotation Matrices  
      Int_t rot1, rot2, rot3;    
      
//   Rotation matrices  
     fMUON->AliMatrix(rot1,  90.,  90., 90., 180.,  0., 0.); // +90 deg in x-y plane
     fMUON->AliMatrix(rot2,  90.,  45., 90., 135.,  0., 0.); // +45 deg in x-y plane 
     fMUON->AliMatrix(rot3,  90.,  45., 90., 315.,180., 0.); // +45 deg in x-y + rotation 180° around y

//   Translation matrices ... NOT USED  
//     fMUON->AliMatrix(trans1, 90.,   0., 90.,  90.,   0., 0.); // X-> X; Y -> Y; Z -> Z
//     fMUON->AliMatrix(trans2, 90., 180., 90.,  90., 180., 0.); // X->-X; Y -> Y; Z ->-Z
//     fMUON->AliMatrix(trans3, 90., 180., 90., 270.,   0., 0.); // X->-X; Y ->-Y; Z -> Z
//     fMUON->AliMatrix(trans4, 90.,   0., 90., 270., 180., 0.); // X-> X; Y ->-Y; Z ->-Z
//  
      // ___________________Volume thicknesses________________________

  const Float_t kHzFrameThickness = 1.59/2.;     //equivalent thickness
  const Float_t kHzOuterFrameEpoxy = 1.19/2.;    //equivalent thickness
  const Float_t kHzOuterFrameInox = 0.1/2.;      //equivalent thickness
  const Float_t kHzFoam = 2.083/2.;              //evaluated elsewhere
                                                 // CHECK with fgkHzFoam
  
// Pertaining to the top outer area 
  const Float_t kHzTopAnodeSteel1 = 0.185/2.;    //equivalent thickness
  const Float_t kHzTopAnodeSteel2 = 0.51/2.;     //equivalent thickness  
  const Float_t kHzAnodeFR4 = 0.08/2.;           //equivalent thickness
  const Float_t kHzTopEarthFaceCu = 0.364/2.;    //equivalent thickness
  const Float_t kHzTopEarthProfileCu = 1.1/2.;   //equivalent thickness
  const Float_t kHzTopPositionerSteel = 1.45/2.; //should really be 2.125/2.; 
  const Float_t kHzTopGasSupportAl = 0.85/2.;    //equivalent thickness
  
// Pertaining to the vertical outer area  
  const Float_t kHzVerticalCradleAl = 0.8/2.;     //equivalent thickness
  const Float_t kHzLateralSightAl = 0.975/2.;     //equivalent thickness
  const Float_t kHzLateralPosnInoxFace = 2.125/2.;//equivalent thickness
  const Float_t kHzLatPosInoxProfM = 6.4/2.;      //equivalent thickness
  const Float_t kHzLatPosInoxProfNF = 1.45/2.;    //equivalent thickness
  const Float_t kHzLateralPosnAl = 0.5/2.;        //equivalent thickness
  const Float_t kHzVertEarthFaceCu = 0.367/2.;    //equivalent thickness
  const Float_t kHzVertBarSteel = 0.198/2.;       //equivalent thickness
  const Float_t kHzVertEarthProfCu = 1.1/2.;      //equivalent thickness

      //_______________Parameter definitions in sequence _________

// InVFrame parameters
  const Float_t kHxInVFrame  = 1.85/2.;
  const Float_t kHyInVFrame  = 73.95/2.;
  const Float_t kHzInVFrame  = kHzFrameThickness;

//Flat 7.5mm vertical section
  const Float_t kHxV1mm  = 0.75/2.;
  const Float_t kHyV1mm  = 1.85/2.;
  const Float_t kHzV1mm  = kHzFrameThickness;

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
  const Float_t kHxTFA = 34.1433/2.;
  const Float_t kHyTFA = 7.75/2.;
  const Float_t kHzTFAE = kHzOuterFrameEpoxy;     // layer 1 thickness
  const Float_t kHzTFAI = kHzOuterFrameInox;      // layer 3 thickness
  
// TopFrameAnodeA parameters - trapezoid, 2 layers
  const Float_t kHzFAAE = kHzOuterFrameEpoxy;     // layer 1 thickness
  const Float_t kHzFAAI = kHzOuterFrameInox;      // layer 3 thickness
  const Float_t kTetFAA = 0.;
  const Float_t kPhiFAA = 0.;
  const Float_t kH1FAA = 8.7/2.;
  const Float_t kBl1FAA = 4.35/2.;
  const Float_t kTl1FAA =  7.75/2.;
  const Float_t kAlp1FAA = 11.06; 
  const Float_t kH2FAA = 8.7/2.;
  const Float_t kBl2FAA = 4.35/2.;
  const Float_t kTl2FAA = 7.75/2.;
  const Float_t kAlp2FAA = 11.06;  
  
// TopFrameAnodeB parameters - trapezoid, 2 layers
  const Float_t kHzFABE = kHzOuterFrameEpoxy;     // layer 1 thickness
  const Float_t kHzFABI = kHzOuterFrameInox;      // layer 3 thickness
  const Float_t kTetFAB = 0.;
  const Float_t kPhiFAB = 0.;
  const Float_t kH1FAB = 8.70/2.;
  const Float_t kBl1FAB = 0.;
  const Float_t kTl1FAB = 4.35/2.;
  const Float_t kAlp1FAB = 14.03; 
  const Float_t kH2FAB = 8.70/2.;
  const Float_t kBl2FAB = 0.;
  const Float_t kTl2FAB = 4.35/2.;
  const Float_t kAlp2FAB = 14.03;  
  
// TopAnode parameters - cuboid (part 1 of 3 parts)
  const Float_t kHxTA1 = 16.2/2.;
  const Float_t kHyTA1 = 3.5/2.;
  const Float_t kHzTA11 = kHzTopAnodeSteel1;   // layer 1
  const Float_t kHzTA12 = kHzAnodeFR4;         // layer 2 

// TopAnode parameters - trapezoid 1 (part 2 of 3 parts)
  const Float_t kHzTA21 = kHzTopAnodeSteel2;   // layer 1 
  const Float_t kHzTA22 = kHzAnodeFR4;         // layer 2 
  const Float_t kTetTA2 = 0.;
  const Float_t kPhiTA2= 0.;
  const Float_t kH1TA2 = 7.268/2.;
  const Float_t kBl1TA2 = 2.03/2.;
  const Float_t kTl1TA2 = 3.5/2.;
  const Float_t kAlp1TA2 = 5.78; 
  const Float_t kH2TA2 = 7.268/2.;
  const Float_t kBl2TA2 = 2.03/2.;
  const Float_t kTl2TA2 = 3.5/2.;
  const Float_t kAlp2TA2 = 5.78;  

// TopAnode parameters - trapezoid 2 (part 3 of 3 parts)
  const Float_t kHzTA3 = kHzAnodeFR4;       // layer 1 
  const Float_t kTetTA3 = 0.;
  const Float_t kPhiTA3 = 0.;
  const Float_t kH1TA3 = 7.268/2.;
  const Float_t kBl1TA3 = 0.;
  const Float_t kTl1TA3 = 2.03/2.;
  const Float_t kAlp1TA3 = 7.95; 
  const Float_t kH2TA3 = 7.268/2.;
  const Float_t kBl2TA3 = 0.;
  const Float_t kTl2TA3 = 2.03/2.;
  const Float_t kAlp2TA3 = 7.95;  
  
// TopEarthFace parameters - single trapezoid
  const Float_t kHzTEF = kHzTopEarthFaceCu;
  const Float_t kTetTEF = 0.;
  const Float_t kPhiTEF = 0.;
  const Float_t kH1TEF = 1.200/2.;
  const Float_t kBl1TEF = 21.323/2.;
  const Float_t kTl1TEF = 17.963/2.;
  const Float_t kAlp1TEF = -54.46; 
  const Float_t kH2TEF = 1.200/2.;
  const Float_t kBl2TEF = 21.323/2.;
  const Float_t kTl2TEF = 17.963/2.;
  const Float_t kAlp2TEF = -54.46;

// TopEarthProfile parameters - single trapezoid
  const Float_t kHzTEP = kHzTopEarthProfileCu;
  const Float_t kTetTEP = 0.;
  const Float_t kPhiTEP = 0.;
  const Float_t kH1TEP = 0.40/2.;
  const Float_t kBl1TEP = 31.766/2.;
  const Float_t kTl1TEP = 30.535/2.;
  const Float_t kAlp1TEP = -56.98; 
  const Float_t kH2TEP = 0.40/2.;
  const Float_t kBl2TEP = 31.766/2.;
  const Float_t kTl2TEP = 30.535/2.;
  const Float_t kAlp2TEP = -56.98;

// TopPositioner parameters - single Stainless Steel trapezoid 
  const Float_t kHzTP = kHzTopPositionerSteel;
  const Float_t kTetTP = 0.;
  const Float_t kPhiTP = 0.;
  const Float_t kH1TP = 3.00/2.;
  const Float_t kBl1TP = 7.023/2.;
  const Float_t kTl1TP = 7.314/2.;
  const Float_t kAlp1TP = 2.78; 
  const Float_t kH2TP = 3.00/2.;
  const Float_t kBl2TP = 7.023/2.;
  const Float_t kTl2TP = 7.314/2.;
  const Float_t kAlp2TP = 2.78;

// TopGasSupport parameters - single cuboid 
  const Float_t kHxTGS  = 8.50/2.;
  const Float_t kHyTGS  = 3.00/2.;
  const Float_t kHzTGS  = kHzTopGasSupportAl;
    
// OutEdgeFrame parameters - 4 trapezoidal sections, 2 layers of material
//
//---

// Trapezoid 1
  const Float_t kHzOETFE = kHzOuterFrameEpoxy;    // layer 1 
  const Float_t kHzOETFI = kHzOuterFrameInox;     // layer 3
   
  const Float_t kTetOETF = 0.;            // common to all 4 trapezoids
  const Float_t kPhiOETF = 0.;            // common to all 4 trapezoids

  const Float_t kH1OETF = 7.196/2.;       // common to all 4 trapezoids
  const Float_t kH2OETF = 7.196/2.;       // common to all 4 trapezoids   
  
  const Float_t kBl1OETF1 = 3.75/2; 
  const Float_t kTl1OETF1 = 3.996/2.;
  const Float_t kAlp1OETF1 = 0.98;

  const Float_t kBl2OETF1 = 3.75/2;
  const Float_t kTl2OETF1 = 3.996/2.;
  const Float_t kAlp2OETF1 = 0.98;
  
// Trapezoid 2
  const Float_t kBl1OETF2 = 3.01/2.;
  const Float_t kTl1OETF2 = 3.75/2;
  const Float_t kAlp1OETF2 = 2.94;
      
  const Float_t kBl2OETF2 = 3.01/2.;
  const Float_t kTl2OETF2 = 3.75/2;
  const Float_t kAlp2OETF2 = 2.94; 
 
// Trapezoid 3
  //const Float_t kBl1OETF3 = 1.767/2.;
  //const Float_t kTl1OETF3 = 3.01/2.;
  const Float_t kBl1OETF3 = 1.117/2.;
  const Float_t kTl1OETF3 = 2.36/2.;
  const Float_t kAlp1OETF3 = 4.94;
        // Fix (5) - overlap of SQ21 with 041M and 125M
      
  //const Float_t kBl2OETF3 = 1.767/2.;
  //const Float_t kTl2OETF3 = 3.01/2.; 
  const Float_t kBl2OETF3 = 1.117/2.;
  const Float_t kTl2OETF3 = 2.36/2.;
  const Float_t kAlp2OETF3 = 4.94; 
        // Fix (5) - overlap of SQ21 with 041M and 125M
  
// Trapezoid 4
  const Float_t kBl1OETF4 = 0.;
  const Float_t kTl1OETF4 = 1.77/2.;
  const Float_t kAlp1OETF4 = 7.01;
      
  const Float_t kBl2OETF4 = 0.;
  const Float_t kTl2OETF4 = 1.77/2.;
  const Float_t kAlp2OETF4 =  7.01;   
  
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
  const Float_t kHxOutVFrame = 1.85/2.;
  const Float_t kHyOutVFrame = 46.23/2.;
  const Float_t kHzOutVFrame = kHzFrameThickness;

// OutVFrame corner parameters - trapezoid
  const Float_t kHzOCTF = kHzFrameThickness;
  const Float_t kTetOCTF = 0.;
  const Float_t kPhiOCTF = 0.;
  const Float_t kH1OCTF = 1.85/2.;
  const Float_t kBl1OCTF = 0.;
  const Float_t kTl1OCTF = 3.66/2.;
  const Float_t kAlp1OCTF = 44.67; 
  const Float_t kH2OCTF = 1.85/2.;
  const Float_t kBl2OCTF = 0.;
  const Float_t kTl2OCTF = 3.66/2.;
  const Float_t kAlp2OCTF = 44.67;  
  
// VertEarthFaceCu parameters - single trapezoid
  const Float_t kHzVFC = kHzVertEarthFaceCu;
  const Float_t kTetVFC = 0.;
  const Float_t kPhiVFC = 0.;
  const Float_t kH1VFC = 1.200/2.;
  const Float_t kBl1VFC = 46.11/2.;
  const Float_t kTl1VFC = 48.236/2.;
  const Float_t kAlp1VFC = 41.54; 
  const Float_t kH2VFC = 1.200/2.;
  const Float_t kBl2VFC = 46.11/2.;
  const Float_t kTl2VFC = 48.236/2.;
  const Float_t kAlp2VFC = 41.54;
    
// VertEarthSteel parameters - single trapezoid
  const Float_t kHzVES = kHzVertBarSteel;
  const Float_t kTetVES = 0.;
  const Float_t kPhiVES = 0.;
  const Float_t kH1VES = 1.200/2.;
  const Float_t kBl1VES = 30.486/2.;
  const Float_t kTl1VES = 32.777/2.;
  const Float_t kAlp1VES = 43.67; 
  const Float_t kH2VES = 1.200/2.;
  const Float_t kBl2VES = 30.486/2.;
  const Float_t kTl2VES = 32.777/2.;
  const Float_t kAlp2VES = 43.67;

// VertEarthProfCu parameters - single trapezoid
  const Float_t kHzVPC = kHzVertEarthProfCu;
  const Float_t kTetVPC = 0.;
  const Float_t kPhiVPC = 0.;
  const Float_t kH1VPC = 0.400/2.;
  const Float_t kBl1VPC = 29.287/2.;
  const Float_t kTl1VPC = 30.091/2.;
  const Float_t kAlp1VPC = 45.14; 
  const Float_t kH2VPC = 0.400/2.;
  const Float_t kBl2VPC = 29.287/2.;
  const Float_t kTl2VPC = 30.091/2.;
  const Float_t kAlp2VPC = 45.14;

// SuppLateralPositionner - single cuboid
  const Float_t kHxSLP  = 2.80/2.;
  const Float_t kHySLP  = 5.00/2.;
  const Float_t kHzSLP  = kHzLateralPosnAl;
  
// LateralPositionner - squared off U bend, face view
  const Float_t kHxLPF  = 5.2/2.;
  const Float_t kHyLPF  = 3.0/2.;
  const Float_t kHzLPF  = kHzLateralPosnInoxFace;
  
// LateralPositionner - squared off U bend, profile view
  const Float_t kHxLPP  = 0.425/2.;
  const Float_t kHyLPP  = 3.0/2.;
  const Float_t kHzLPP  = kHzLatPosInoxProfM;  // middle layer
  const Float_t kHzLPNF  = kHzLatPosInoxProfNF; // near and far layers
           
// VertCradle, 3 layers (copies), each composed of 4 trapezoids
// VertCradleA
  const Float_t kHzVC1 = kHzVerticalCradleAl;
  const Float_t kTetVC1 = 0.;
  const Float_t kPhiVC1 = 0.;
  const Float_t kH1VC1 = 10.25/2.;
  const Float_t kBl1VC1 = 3.70/2.;
  const Float_t kTl1VC1 = 0.;
  const Float_t kAlp1VC1 = -10.23; 
  const Float_t kH2VC1 = 10.25/2.;
  const Float_t kBl2VC1 = 3.70/2.;
  const Float_t kTl2VC1 = 0.;
  const Float_t kAlp2VC1 = -10.23;
        
// VertCradleB
  const Float_t kHzVC2 = kHzVerticalCradleAl;
  const Float_t kTetVC2 = 0.;
  const Float_t kPhiVC2 = 0.;
  const Float_t kH1VC2 = 10.25/2.;
  const Float_t kBl1VC2 = 6.266/2.;
  const Float_t kTl1VC2 = 3.70/2.;
  const Float_t kAlp1VC2 = -7.13; 
  const Float_t kH2VC2 = 10.25/2.;
  const Float_t kBl2VC2 = 6.266/2.;
  const Float_t kTl2VC2 = 3.70/2.;
  const Float_t kAlp2VC2 = -7.13;
  
// VertCradleC
  const Float_t kHzVC3 = kHzVerticalCradleAl;
  const Float_t kTetVC3 = 0.;
  const Float_t kPhiVC3 = 0.;
  const Float_t kH1VC3 = 10.25/2.;
  const Float_t kBl1VC3 = 7.75/2.;
  const Float_t kTl1VC3 = 6.266/2.;
  const Float_t kAlp1VC3 = -4.14; 
  const Float_t kH2VC3 = 10.25/2.;
  const Float_t kBl2VC3 = 7.75/2.;
  const Float_t kTl2VC3 = 6.266/2.;
  const Float_t kAlp2VC3 = -4.14;

// VertCradleD
  const Float_t kHzVC4 = kHzVerticalCradleAl;
  const Float_t kTetVC4 = 0.;
  const Float_t kPhiVC4 = 0.;
  const Float_t kH1VC4 = 10.27/2.;
  const Float_t kBl1VC4 = 8.273/2.;
  const Float_t kTl1VC4 = 7.75/2.;
  const Float_t kAlp1VC4 = -1.46; 
  const Float_t kH2VC4 = 10.27/2.;
  const Float_t kBl2VC4 = 8.273/2.;
  const Float_t kTl2VC4 = 7.75/2.;
  const Float_t kAlp2VC4 = -1.46;
  
// LateralSightSupport - single trapezoid
  const Float_t kHzVSS = kHzLateralSightAl;
  const Float_t kTetVSS = 0.;
  const Float_t kPhiVSS = 0.;
  const Float_t kH1VSS = 5.00/2.;
  const Float_t kBl1VSS = 7.747/2;
  const Float_t kTl1VSS = 7.188/2.;
  const Float_t kAlp1VSS = -3.20; 
  const Float_t kH2VSS = 5.00/2.;
  const Float_t kBl2VSS = 7.747/2.;
  const Float_t kTl2VSS = 7.188/2.;
  const Float_t kAlp2VSS = -3.20;  
  
// LateralSight (reference point) - 3 per quadrant, only 1 programmed for now
  const Float_t kVSInRad  = 0.6;
  const Float_t kVSOutRad  = 1.3;
  const Float_t kVSLen  = kHzFrameThickness; 
  
//---

// InHFrame parameters
  const Float_t kHxInHFrame  = 75.8/2.;
  const Float_t kHyInHFrame  = 1.85/2.;
  const Float_t kHzInHFrame  = kHzFrameThickness;
 
//Flat 7.5mm horizontal section
  const Float_t kHxH1mm  = 1.85/2.;
  const Float_t kHyH1mm  = 0.75/2.;
  const Float_t kHzH1mm  = kHzFrameThickness;

//---

// InArcFrame parameters
  const Float_t kIAF  = 15.70;
  const Float_t kOAF  = 17.55;
  const Float_t kHzAF  = kHzFrameThickness;
  const Float_t kAFphi1  = 0.0;
  const Float_t kAFphi2  = 90.0;

//---

// ScrewsInFrame parameters HEAD
  const Float_t kSCRUHMI  = 0.;
  const Float_t kSCRUHMA  = 0.690/2.;
  const Float_t kSCRUHLE  = 0.4/2.;
// ScrewsInFrame parameters MIDDLE
  const Float_t kSCRUMMI  = 0.;
  const Float_t kSCRUMMA  = 0.39/2.;
  const Float_t kSCRUMLE  = kHzFrameThickness;
// ScrewsInFrame parameters NUT
  const Float_t kSCRUNMI  = 0.;
  const Float_t kSCRUNMA  = 0.78/2.;
  const Float_t kSCRUNLE  = 0.8/2.;   
  
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

// Replace the volume shape with a composite shape
// with substracted overlap with beam shield (YMOT)

  if ( gMC->IsRootGeometrySupported() &&
       TString(gMC->ClassName()) != "TGeant4") { 

    // Get shape
    TGeoVolume* mlayer 
      = gGeoManager->FindVolumeFast(QuadrantMLayerName(chamber));
    if ( !mlayer ) {
      AliErrorStream() 
         << "Quadrant volume " << QuadrantMLayerName(chamber) << " not found" 
	 << endl;
    }
    else {
      TGeoShape* quadrant = mlayer->GetShape();
      quadrant->SetName("quadrant");	 

      // Beam shield recess
      par[0] = 0;
      par[1] = 15.4; 
      par[2] = fgkMotherThick1;  
      new TGeoTube("shield_tube", par[0], par[1], par[2]);
  
      // Displacement
      posX = 2.6;
      posY = 2.6;
      posZ = 0;
      TGeoTranslation* displacement 
        = new TGeoTranslation("TR", posX, posY, posZ);
      displacement->RegisterYourself();

      // Composite shape
      TGeoShape* composite
      = new TGeoCompositeShape("composite", "quadrant-shield_tube:TR"); 
      
      // Reset shape to volume      
      mlayer->SetShape(composite);
    }
  }    

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
    par[0] = kHxInVFrame;
    par[1] = kHyInVFrame;
    par[2] = kHzInVFrame;
    gMC->Gsvolu("SQ00","BOX",idFrameEpoxy,par,3);

    //Flat 1mm vertical section
    par[0] = kHxV1mm;
    par[1] = kHyV1mm;
    par[2] = kHzV1mm;
    gMC->Gsvolu("SQ01","BOX",idFrameEpoxy,par,3); 
 
// OutTopFrame 
//
// - 3 components (a cuboid and 2 trapezes) and 2 layers (Epoxy/Inox)
//
//---

    // TopFrameAnode - layer 1 of 2 
    par[0] = kHxTFA;
    par[1] = kHyTFA;
    par[2] = kHzTFAE;
    gMC->Gsvolu("SQ02","BOX",idFrameEpoxy,par,3);
    
    // TopFrameAnode - layer 2 of 2 
    par[2] = kHzTFAI;
    gMC->Gsvolu("SQ03","BOX",idInox,par,3);
            
    // TopFrameAnodeA - layer 1 of 2  
    par[0] = kHzFAAE;
    par[1] = kTetFAA;
    par[2] = kPhiFAA;
    par[3] = kH1FAA;
    par[4] = kBl1FAA;
    par[5] = kTl1FAA;
    par[6] = kAlp1FAA;
    par[7] = kH2FAA;
    par[8] = kBl2FAA;
    par[9] = kTl2FAA;
    par[10] = kAlp2FAA;    
    gMC->Gsvolu("SQ04","TRAP",idFrameEpoxy,par,11);    

    // TopFrameAnodeA - layer 2 of 2
    par[0] = kHzFAAI;    
    gMC->Gsvolu("SQ05","TRAP",idInox,par,11); 
      
    // TopFrameAnodeB - layer 1 of 2
    par[0] = kHzFABE;
    par[1] = kTetFAB;
    par[2] = kPhiFAB;
    par[3] = kH1FAB;
    par[4] = kBl1FAB;
    par[5] = kTl1FAB;
    par[6] = kAlp1FAB;
    par[7] = kH2FAB;
    par[8] = kBl2FAB;
    par[9] = kTl2FAB;
    par[10] = kAlp2FAB;
    gMC->Gsvolu("SQ06","TRAP",idFrameEpoxy,par,11);     

    // OutTopTrapFrameB - layer 2 of 2
    par[0] = kHzFABI;   
    gMC->Gsvolu("SQ07","TRAP",idInox,par,11);

    // TopAnode1 -  layer 1 of 2
    par[0] = kHxTA1;
    par[1] = kHyTA1;
    par[2] = kHzTA11;    
    gMC->Gsvolu("SQ08","BOX",idInox,par,3); 
    
    // TopAnode1 -  layer 2 of 2
    par[2] = kHzTA12;    
    gMC->Gsvolu("SQ09","BOX",idFR4,par,11); 

    // TopAnode2 -  layer 1 of 2
    par[0] = kHzTA21;
    par[1] = kTetTA2;
    par[2] = kPhiTA2;
    par[3] = kH1TA2;
    par[4] = kBl1TA2;
    par[5] = kTl1TA2;
    par[6] = kAlp1TA2;
    par[7] = kH2TA2;
    par[8] = kBl2TA2;
    par[9] = kTl2TA2;
    par[10] = kAlp2TA2;    
    gMC->Gsvolu("SQ10","TRAP",idInox,par,11); 
 
    // TopAnode2 -  layer 2 of 2
    par[0] = kHzTA22;    
    gMC->Gsvolu("SQ11","TRAP",idFR4,par,11);   

    // TopAnode3 -  layer 1 of 1 
    par[0] = kHzTA3;
    par[1] = kTetTA3;
    par[2] = kPhiTA3;
    par[3] = kH1TA3;
    par[4] = kBl1TA3;
    par[5] = kTl1TA3;
    par[6] = kAlp1TA3;
    par[7] = kH2TA3;
    par[8] = kBl2TA3;
    par[9] = kTl2TA3;
    par[10] = kAlp2TA3;    
    gMC->Gsvolu("SQ12","TRAP",idFR4,par,11); 

    // TopEarthFace 
    par[0] = kHzTEF;
    par[1] = kTetTEF;
    par[2] = kPhiTEF;
    par[3] = kH1TEF;
    par[4] = kBl1TEF;
    par[5] = kTl1TEF;
    par[6] = kAlp1TEF;
    par[7] = kH2TEF;
    par[8] = kBl2TEF;
    par[9] = kTl2TEF;
    par[10] = kAlp2TEF;    
    gMC->Gsvolu("SQ13","TRAP",idCopper,par,11);   

    // TopEarthProfile 
    par[0] = kHzTEP;
    par[1] = kTetTEP;
    par[2] = kPhiTEP;
    par[3] = kH1TEP;
    par[4] = kBl1TEP;
    par[5] = kTl1TEP;
    par[6] = kAlp1TEP;
    par[7] = kH2TEP;
    par[8] = kBl2TEP;
    par[9] = kTl2TEP;
    par[10] = kAlp2TEP;
    gMC->Gsvolu("SQ14","TRAP",idCopper,par,11);       

    // TopGasSupport  
    par[0] = kHxTGS;
    par[1] = kHyTGS;
    par[2] = kHzTGS;
    gMC->Gsvolu("SQ15","BOX",idAlu,par,3);

    // TopPositioner parameters - single Stainless Steel trapezoid 
    par[0] = kHzTP;
    par[1] = kTetTP; 
    par[2] = kPhiTP;
    par[3] = kH1TP;
    par[4] = kBl1TP; 
    par[5] = kTl1TP; 
    par[6] = kAlp1TP;
    par[7] = kH2TP;
    par[8] = kBl2TP; 
    par[9] = kTl2TP; 
    par[10] = kAlp2TP;     
    gMC->Gsvolu("SQ16","TRAP",idInox,par,11);       

//
// OutEdgeTrapFrame Epoxy = (4 trapezes)*2 copies*2 layers (Epoxy/Inox)
//
//---
    // Trapezoid 1 - 2 layers
    par[1] = kTetOETF;
    par[2] = kPhiOETF;
    par[3] = kH1OETF;
    par[4] = kBl1OETF1;
    par[5] = kTl1OETF1;
    par[6] = kAlp1OETF1;
    par[7] = kH2OETF;
    par[8] = kBl2OETF1;
    par[9] = kTl2OETF1;
    par[10] = kAlp2OETF1; 
           
    par[0] = kHzOETFE;             
    gMC->Gsvolu("SQ17","TRAP",idFrameEpoxy,par,11); 
    par[0] = kHzOETFI;
    gMC->Gsvolu("SQ18","TRAP",idInox,par,11);
    
    // Trapezoid 2 - 2 layers
    par[4] = kBl1OETF2;
    par[5] = kTl1OETF2;
    par[6] = kAlp1OETF2;

    par[8] = kBl2OETF2;
    par[9] = kTl2OETF2;
    par[10] = kAlp2OETF2; 
    
    par[0] = kHzOETFE;    
    gMC->Gsvolu("SQ19","TRAP",idFrameEpoxy,par,11);    
    par[0] = kHzOETFI;    
    gMC->Gsvolu("SQ20","TRAP",idInox,par,11);     
    
    // Trapezoid 3 - 2 layers
    par[4] = kBl1OETF3;
    par[5] = kTl1OETF3;
    par[6] = kAlp1OETF3;

    par[8] = kBl2OETF3;
    par[9] = kTl2OETF3;
    par[10] = kAlp2OETF3; 
 
    par[0] = kHzOETFE;    
    gMC->Gsvolu("SQ21","TRAP",idFrameEpoxy,par,11);   
    par[0] = kHzOETFI;    
    gMC->Gsvolu("SQ22","TRAP",idInox,par,11);     
    
    // Trapezoid 4 - 2 layers

    par[4] = kBl1OETF4;
    par[5] = kTl1OETF4;
    par[6] = kAlp1OETF4;

    par[8] = kBl2OETF4;
    par[9] = kTl2OETF4;
    par[10] = kAlp2OETF4;  
   
    par[0] = kHzOETFE;    
    gMC->Gsvolu("SQ23","TRAP",idFrameEpoxy,par,11);    
    par[0] = kHzOETFI;    
    gMC->Gsvolu("SQ24","TRAP",idInox,par,11);     
             
//---
    // OutVFrame    
    par[0] = kHxOutVFrame;
    par[1] = kHyOutVFrame;
    par[2] = kHzOutVFrame;
    gMC->Gsvolu("SQ25","BOX",idFrameEpoxy,par,3);
        
    // OutVFrame corner  
    par[0] = kHzOCTF;
    par[1] = kTetOCTF;
    par[2] = kPhiOCTF;
    par[3] = kH1OCTF;
    par[4] = kBl1OCTF;
    par[5] = kTl1OCTF;
    par[6] = kAlp1OCTF;
    par[7] = kH2OCTF;
    par[8] = kBl2OCTF;
    par[9] = kTl2OCTF;
    par[10] = kAlp2OCTF;    
    gMC->Gsvolu("SQ26","TRAP",idFrameEpoxy,par,11);
 
    // EarthFaceCu trapezoid
    par[0] = kHzVFC;
    par[1] = kTetVFC;
    par[2] = kPhiVFC;
    par[3] = kH1VFC;
    par[4] = kBl1VFC;
    par[5] = kTl1VFC;
    par[6] = kAlp1VFC;
    par[7] = kH2VFC;
    par[8] = kBl2VFC;
    par[9] = kTl2VFC;
    par[10] = kAlp2VFC;   
    gMC->Gsvolu("SQ27","TRAP",idCopper,par,11);     

    // VertEarthSteel trapezoid
    par[0] = kHzVES;
    par[1] = kTetVES;
    par[2] = kPhiVES;
    par[3] = kH1VES;
    par[4] = kBl1VES;
    par[5] = kTl1VES;
    par[6] = kAlp1VES;
    par[7] = kH2VES;
    par[8] = kBl2VES;
    par[9] = kTl2VES;
    par[10] = kAlp2VES;    
    gMC->Gsvolu("SQ28","TRAP",idInox,par,11); 

    // VertEarthProfCu trapezoid       
    par[0] = kHzVPC;
    par[1] = kTetVPC;
    par[2] = kPhiVPC;
    par[3] = kH1VPC;
    par[4] = kBl1VPC;
    par[5] = kTl1VPC;
    par[6] = kAlp1VPC;
    par[7] = kH2VPC;
    par[8] = kBl2VPC;
    par[9] = kTl2VPC;
    par[10] = kAlp2VPC;
    gMC->Gsvolu("SQ29","TRAP",idCopper,par,11);

    // SuppLateralPositionner cuboid    
    par[0] = kHxSLP;
    par[1] = kHySLP;
    par[2] = kHzSLP;
    gMC->Gsvolu("SQ30","BOX",idAlu,par,3);

    // LateralPositionerFace
    par[0] = kHxLPF;
    par[1] = kHyLPF;
    par[2] = kHzLPF;
    gMC->Gsvolu("SQ31","BOX",idInox,par,3);

    // LateralPositionerProfile
    par[0] = kHxLPP;
    par[1] = kHyLPP;
    par[2] = kHzLPP;
    gMC->Gsvolu("SQ32","BOX",idInox,par,3); // middle layer
    
    par[0] = kHxLPP;
    par[1] = kHyLPP;
    par[2] = kHzLPNF;
    gMC->Gsvolu("SQ33","BOX",idInox,par,3); // near and far layers

    // VertCradleA - 1st trapezoid
    par[0] = kHzVC1;
    par[1] = kTetVC1;
    par[2] = kPhiVC1;
    par[3] = kH1VC1;
    par[4] = kBl1VC1;
    par[5] = kTl1VC1;
    par[6] = kAlp1VC1;
    par[7] = kH2VC1;
    par[8] = kBl2VC1;
    par[9] = kTl2VC1;
    par[10] = kAlp2VC1;
    gMC->Gsvolu("SQ34","TRAP",idAlu,par,11); 
    
    // VertCradleB - 2nd trapezoid
    par[0] = kHzVC2;
    par[1] = kTetVC2;
    par[2] = kPhiVC2;
    par[3] = kH1VC2;
    par[4] = kBl1VC2;
    par[5] = kTl1VC2;
    par[6] = kAlp1VC2;
    par[7] = kH2VC2;
    par[8] = kBl2VC2;
    par[9] = kTl2VC2;
    par[10] = kAlp2VC2;
    gMC->Gsvolu("SQ35","TRAP",idAlu,par,11);  
       
    // VertCradleC - 3rd trapezoid
    par[0] = kHzVC3;
    par[1] = kTetVC3;
    par[2] = kPhiVC3;
    par[3] = kH1VC3;
    par[4] = kBl1VC3;
    par[5] = kTl1VC3;
    par[6] = kAlp1VC3;
    par[7] = kH2VC3;
    par[8] = kBl2VC3;
    par[9] = kTl2VC3;
    par[10] = kAlp2VC3;    
    gMC->Gsvolu("SQ36","TRAP",idAlu,par,11);  

    // VertCradleD - 4th trapezoid
    par[0] = kHzVC4;
    par[1] = kTetVC4;
    par[2] = kPhiVC4;
    par[3] = kH1VC4;
    par[4] = kBl1VC4;
    par[5] = kTl1VC4;
    par[6] = kAlp1VC4;
    par[7] = kH2VC4;
    par[8] = kBl2VC4;
    par[9] = kTl2VC4;
    par[10] = kAlp2VC4;    
    gMC->Gsvolu("SQ37","TRAP",idAlu,par,11);  
          
    // LateralSightSupport trapezoid
    par[0] = kHzVSS;
    par[1] = kTetVSS;
    par[2] = kPhiVSS;
    par[3] = kH1VSS;
    par[4] = kBl1VSS;
    par[5] = kTl1VSS;
    par[6] = kAlp1VSS;
    par[7] = kH2VSS;
    par[8] = kBl2VSS;
    par[9] = kTl2VSS;
    par[10] = kAlp2VSS;
    gMC->Gsvolu("SQ38","TRAP",idAlu,par,11);

    // LateralSight
    par[0] = kVSInRad;
    par[1] = kVSOutRad;
    par[2] = kVSLen;       
    gMC->Gsvolu("SQ39","TUBE",idFrameEpoxy,par,3);   

//---
    // InHFrame
    par[0] = kHxInHFrame;
    par[1] = kHyInHFrame;
    par[2] = kHzInHFrame;
    gMC->Gsvolu("SQ40","BOX",idFrameEpoxy,par,3);

    //Flat 7.5mm horizontal section
    par[0] = kHxH1mm;
    par[1] = kHyH1mm;
    par[2] = kHzH1mm;
    gMC->Gsvolu("SQ41","BOX",idFrameEpoxy,par,3);

    // InArcFrame 
    par[0] = kIAF;
    par[1] = kOAF; 
    par[2] = kHzAF;  
    par[3] = kAFphi1; 
    par[4] = kAFphi2;

    gMC->Gsvolu("SQ42","TUBS",idFrameEpoxy,par,5);

//---
    // ScrewsInFrame - 3 sections in order to avoid overlapping volumes
    // Screw Head, in air
    par[0] = kSCRUHMI;
    par[1] = kSCRUHMA; 
    par[2] = kSCRUHLE;  

    gMC->Gsvolu("SQ43","TUBE",idInox,par,3);
    
    // Middle part, in the Epoxy
    par[0] = kSCRUMMI;
    par[1] = kSCRUMMA;
    par[2] = kSCRUMLE;
    gMC->Gsvolu("SQ44","TUBE",idInox,par,3);
    
    // Screw nut, in air
    par[0] = kSCRUNMI;
    par[1] = kSCRUNMA;
    par[2] = kSCRUNLE;   
    gMC->Gsvolu("SQ45","TUBE",idInox,par,3);     
   }
              
// __________________Place volumes in the quadrant ____________ 
        
    // InVFrame  
    posX = kHxInVFrame;
    posY = 2.0*kHyInHFrame+2.*kHyH1mm+kIAF+kHyInVFrame;        
    posZ = 0.;
    gMC->Gspos("SQ00",1,QuadrantMLayerName(chamber),posX, posY, posZ, 0, "ONLY"); 

// keep memory of the mid position. Used for placing screws
    const GReal_t kMidVposX = posX;
    const GReal_t kMidVposY = posY;
    const GReal_t kMidVposZ = posZ;

    //Flat 7.5mm vertical section
    posX = 2.0*kHxInVFrame+kHxV1mm;
    posY = 2.0*kHyInHFrame+2.*kHyH1mm+kIAF+kHyV1mm;
    posZ = 0.;
    gMC->Gspos("SQ01",1,QuadrantMLayerName(chamber),posX, posY, posZ,0, "ONLY"); 
    
    // TopFrameAnode place 2 layers of TopFrameAnode cuboids  
    posX = kHxTFA;
    posY = 2.*kHyInHFrame+2.*kHyH1mm+kIAF+2.*kHyInVFrame+kHyTFA;   
    posZ = kHzOuterFrameInox;
    gMC->Gspos("SQ02",1,QuadrantMLayerName(chamber),posX, posY, posZ,0,"ONLY"); 
    posZ = posZ+kHzOuterFrameInox;
    gMC->Gspos("SQ03",1,QuadrantMLayerName(chamber),posX, posY, posZ,0,"ONLY");
    
    // place 2 layers of TopFrameAnodeA trapezoids 
    posX = 35.8932+fgkDeltaQuadLHC;
    posY = 92.6745+fgkDeltaQuadLHC;
    posZ = kHzOuterFrameInox; 
    gMC->Gspos("SQ04",1,QuadrantMLayerName(chamber),posX, posY, posZ, rot1,"ONLY");
    posZ = posZ+kHzOuterFrameInox;
    gMC->Gspos("SQ05",1,QuadrantMLayerName(chamber),posX, posY, posZ, rot1,"ONLY");
    
    // place 2 layers of TopFrameAnodeB trapezoids 
    posX = 44.593+fgkDeltaQuadLHC;
    posY = 90.737+fgkDeltaQuadLHC;
    posZ = kHzOuterFrameInox; 
    gMC->Gspos("SQ06",1,QuadrantMLayerName(chamber),posX, posY, posZ, rot1,"ONLY");
    posZ = posZ+kHzOuterFrameInox;
    gMC->Gspos("SQ07",1,QuadrantMLayerName(chamber),posX, posY, posZ, rot1,"ONLY");    

    // TopAnode1 place 2 layers  
    posX = 6.8+fgkDeltaQuadLHC;
    posY = 99.85+fgkDeltaQuadLHC;
    posZ = -1.*kHzAnodeFR4;
    gMC->Gspos("SQ08",1,QuadrantMLayerName(chamber),posX, posY, posZ, 0,"ONLY");  
    posZ = posZ+kHzTopAnodeSteel1;
    gMC->Gspos("SQ09",1,QuadrantMLayerName(chamber),posX, posY, posZ, 0,"ONLY");    
         
    // TopAnode2 place 2 layers
    posX = 18.534+fgkDeltaQuadLHC;
    posY = 99.482+fgkDeltaQuadLHC; 
    posZ = -1.*kHzAnodeFR4;    
    gMC->Gspos("SQ10",1,QuadrantMLayerName(chamber),posX, posY, posZ, rot1,"ONLY");
    posZ = posZ+kHzTopAnodeSteel2;    
    gMC->Gspos("SQ11",1,QuadrantMLayerName(chamber),posX, posY, posZ, rot1,"ONLY");       
    
    // TopAnode3 place 1 layer
    posX = 25.80+fgkDeltaQuadLHC;
    posY = 98.61+fgkDeltaQuadLHC;
    posZ = 0.;    
    gMC->Gspos("SQ12",1,QuadrantMLayerName(chamber),posX, posY, posZ, rot1,"ONLY");  
          
    // TopEarthFace - 2 copies
    posX = 23.122+fgkDeltaQuadLHC;
    posY = 96.90+fgkDeltaQuadLHC;
    posZ = kHzOuterFrameEpoxy+kHzOuterFrameInox+kHzTopEarthFaceCu;
    gMC->Gspos("SQ13",1,QuadrantMLayerName(chamber),posX, posY, posZ, 0,"ONLY");
    posZ = -1.*posZ;
    gMC->Gspos("SQ13",2,QuadrantMLayerName(chamber),posX, posY, posZ, 0,"ONLY");

    // TopEarthProfile 
    posX = 14.475+fgkDeltaQuadLHC;
    posY = 97.900+fgkDeltaQuadLHC; 
    posZ = kHzTopEarthProfileCu;
    gMC->Gspos("SQ14",1,QuadrantMLayerName(chamber),posX, posY, posZ, 0,"ONLY");
    posZ = -1.0*posZ;
    gMC->Gspos("SQ14",2,QuadrantMLayerName(chamber),posX, posY, posZ, 0,"ONLY");

    // TopGasSupport - 2 copies                            
    posX = 4.9500+fgkDeltaQuadLHC;
    posY = 96.200+fgkDeltaQuadLHC;
    posZ = kHzOuterFrameEpoxy+kHzOuterFrameInox+kHzTopGasSupportAl;
    gMC->Gspos("SQ15",1,QuadrantMLayerName(chamber),posX, posY, posZ, 0,"ONLY");
    posZ = -1.*posZ;
    gMC->Gspos("SQ15",2,QuadrantMLayerName(chamber),posX, posY, posZ, 0,"ONLY");
    
    // TopPositioner parameters - single Stainless Steel trapezoid - 2 copies
    posX = 7.60+fgkDeltaQuadLHC;
    posY = 98.98+fgkDeltaQuadLHC;   
    posZ = kHzOuterFrameEpoxy+kHzOuterFrameInox+2.*kHzTopGasSupportAl+kHzTopPositionerSteel;
    gMC->Gspos("SQ16",1,QuadrantMLayerName(chamber),posX, posY, posZ, rot1,"ONLY");
    posZ = -1.*posZ;
    gMC->Gspos("SQ16",2,QuadrantMLayerName(chamber),posX, posY, posZ, rot1,"ONLY"); 

    // OutEdgeFrame 
    Float_t xCenter[8]; 
    Float_t yCenter[8];
    
    xCenter[0] = 73.201 + fgkDeltaQuadLHC;
    xCenter[1] = 78.124 + fgkDeltaQuadLHC; 
    //xCenter[2] = 82.862 + fgkDeltaQuadLHC;
    xCenter[2] = 83.102 + fgkDeltaQuadLHC;
    xCenter[3] = 87.418 + fgkDeltaQuadLHC; 
        // Fix (5) - overlap of SQ21 with 041M and 125M
    
    yCenter[0] = 68.122 + fgkDeltaQuadLHC;
    yCenter[1] = 62.860 + fgkDeltaQuadLHC;   
    //yCenter[2] = 57.420 + fgkDeltaQuadLHC;
    yCenter[2] = 57.660 + fgkDeltaQuadLHC;
    yCenter[3] = 51.800 + fgkDeltaQuadLHC; 
        // Fix (5) - overlap of SQ21 with 041M and 125M
      
    xCenter[4] = 68.122 + fgkDeltaQuadLHC;
    xCenter[5] = 62.860 + fgkDeltaQuadLHC; 
    xCenter[6] = 57.420 + fgkDeltaQuadLHC;
    xCenter[7] = 51.800 + fgkDeltaQuadLHC; 
    
    yCenter[4] = 73.210 + fgkDeltaQuadLHC;
    yCenter[5] = 78.124 + fgkDeltaQuadLHC; 
    yCenter[6] = 82.862 + fgkDeltaQuadLHC;
    yCenter[7] = 87.418 + fgkDeltaQuadLHC; 
      
    posZ = -1.0*kHzOuterFrameInox;     
    gMC->Gspos("SQ17",1,QuadrantMLayerName(chamber), xCenter[0], yCenter[0], posZ, rot2,"ONLY");
    gMC->Gspos("SQ17",2,QuadrantMLayerName(chamber), xCenter[4], yCenter[4], posZ, rot3,"ONLY");

    gMC->Gspos("SQ19",1,QuadrantMLayerName(chamber), xCenter[1], yCenter[1], posZ, rot2,"ONLY");   
    gMC->Gspos("SQ19",2,QuadrantMLayerName(chamber), xCenter[5], yCenter[5], posZ, rot3,"ONLY");

    gMC->Gspos("SQ21",1,QuadrantMLayerName(chamber), xCenter[2], yCenter[2], posZ, rot2,"ONLY");
    gMC->Gspos("SQ21",2,QuadrantMLayerName(chamber), xCenter[6], yCenter[6], posZ, rot3,"ONLY");
    
    gMC->Gspos("SQ23",1,QuadrantMLayerName(chamber), xCenter[3], yCenter[3], posZ, rot2,"ONLY");
    gMC->Gspos("SQ23",2,QuadrantMLayerName(chamber), xCenter[7], yCenter[7], posZ, rot3,"ONLY");
     
    posZ = posZ+kHzOuterFrameEpoxy;
   
    gMC->Gspos("SQ18",1,QuadrantMLayerName(chamber), xCenter[0], yCenter[0], posZ, rot2,"ONLY");
    gMC->Gspos("SQ18",2,QuadrantMLayerName(chamber), xCenter[4], yCenter[4], posZ, rot3,"ONLY");
    
    gMC->Gspos("SQ20",1,QuadrantMLayerName(chamber), xCenter[1], yCenter[1], posZ, rot2,"ONLY");   
    gMC->Gspos("SQ20",2,QuadrantMLayerName(chamber), xCenter[5], yCenter[5], posZ, rot3,"ONLY");

    gMC->Gspos("SQ22",1,QuadrantMLayerName(chamber), xCenter[2], yCenter[2], posZ, rot2,"ONLY");
    gMC->Gspos("SQ22",2,QuadrantMLayerName(chamber), xCenter[6], yCenter[6], posZ, rot3,"ONLY");
       
    gMC->Gspos("SQ24",1,QuadrantMLayerName(chamber), xCenter[3], yCenter[3], posZ, rot2,"ONLY");
    gMC->Gspos("SQ24",2,QuadrantMLayerName(chamber), xCenter[7], yCenter[7], posZ, rot3,"ONLY");  

//---    
        
// OutVFrame
    posX = 2.*kHxInVFrame+kIAF+2.*kHxInHFrame-kHxOutVFrame+2.*kHxV1mm;
    posY = 2.*kHyInHFrame+kHyOutVFrame;    
    posZ = 0.;              
    gMC->Gspos("SQ25",1,QuadrantMLayerName(chamber),posX, posY, posZ, 0, "ONLY"); 

 // keep memory of the mid position. Used for placing screws
    const GReal_t kMidOVposX = posX;
    const GReal_t kMidOVposY = posY;
    const GReal_t kMidOVposZ = posZ;

    const Float_t kTOPY = posY+kHyOutVFrame;
    const Float_t kOUTX = posX;

// OutVFrame corner
    posX = kOUTX;
    posY = kTOPY+((kBl1OCTF+kTl1OCTF)/2.);
    posZ = 0.;     
    gMC->Gspos("SQ26",1,QuadrantMLayerName(chamber),posX, posY, posZ, rot1,"ONLY"); 

// VertEarthFaceCu - 2 copies
    posX = 89.4000+fgkDeltaQuadLHC;
    posY = 25.79+fgkDeltaQuadLHC;    
    posZ = kHzFrameThickness+2.0*kHzFoam+kHzVertEarthFaceCu;              
    gMC->Gspos("SQ27",1,QuadrantMLayerName(chamber),posX, posY, posZ, rot1, "ONLY"); 
    posZ = -1.0*posZ; 
    gMC->Gspos("SQ27",2,QuadrantMLayerName(chamber),posX, posY, posZ, rot1, "ONLY"); 
    
// VertEarthSteel - 2 copies
    posX = 91.00+fgkDeltaQuadLHC;
    posY = 30.616+fgkDeltaQuadLHC;    
    posZ = kHzFrameThickness+2.0*kHzFoam+kHzVertBarSteel;              
    gMC->Gspos("SQ28",1,QuadrantMLayerName(chamber),posX, posY, posZ, rot1, "ONLY"); 
    posZ = -1.0*posZ;              
    gMC->Gspos("SQ28",2,QuadrantMLayerName(chamber),posX, posY, posZ, rot1, "ONLY");
 
// VertEarthProfCu - 2 copies
    posX = 92.000+fgkDeltaQuadLHC;
    posY = 29.64+fgkDeltaQuadLHC;    
    posZ = kHzFrameThickness;              
    gMC->Gspos("SQ29",1,QuadrantMLayerName(chamber),posX, posY, posZ, rot1, "ONLY"); 
    posZ = -1.0*posZ;    
    gMC->Gspos("SQ29",2,QuadrantMLayerName(chamber),posX, posY, posZ, rot1, "ONLY"); 

// SuppLateralPositionner - 2 copies 
    posX = 90.2-kNearFarLHC;
    posY = 5.00-kNearFarLHC;    
    posZ = kHzLateralPosnAl-fgkMotherThick2;             
    gMC->Gspos("SQ30",1,QuadrantFLayerName(chamber),posX, posY, posZ, 0, "ONLY"); 
    posZ = -1.0*posZ;            
    gMC->Gspos("SQ30",2,QuadrantNLayerName(chamber),posX, posY, posZ, 0, "ONLY"); 

// LateralPositionner - 2 copies - Face view
    posX = 92.175-kNearFarLHC-2.*kHxLPP;
    posY = 5.00-kNearFarLHC;   
    posZ =2.0*kHzLateralPosnAl+kHzLateralPosnInoxFace-fgkMotherThick2;              
    gMC->Gspos("SQ31",1,QuadrantFLayerName(chamber),posX, posY, posZ, 0, "ONLY"); 
    posZ = -1.0*posZ;             
    gMC->Gspos("SQ31",2,QuadrantNLayerName(chamber),posX, posY, posZ, 0, "ONLY"); 

// LateralPositionner -  Profile view   
    posX = 92.175+fgkDeltaQuadLHC+kHxLPF-kHxLPP;
    posY = 5.00+fgkDeltaQuadLHC;    
    posZ = 0.;              
    gMC->Gspos("SQ32",1,QuadrantMLayerName(chamber),posX, posY, posZ, 0, "ONLY"); // middle layer

    posX = 92.175-kNearFarLHC+kHxLPF-kHxLPP; 
    posY = 5.0000-kNearFarLHC;    
    posZ = fgkMotherThick2-kHzLPNF;              
    gMC->Gspos("SQ33",1,QuadrantNLayerName(chamber),posX, posY, posZ, 0, "ONLY"); // near layer
    posZ = -1.*posZ;
    gMC->Gspos("SQ33",2,QuadrantFLayerName(chamber),posX, posY, posZ, 0, "ONLY"); // far layer
      
// VertCradleA  1st Trapezoid - 3 copies
    posX = 95.73+fgkDeltaQuadLHC;
    posY = 33.26+fgkDeltaQuadLHC; 
    posZ = 0.;              
    gMC->Gspos("SQ34",2,QuadrantMLayerName(chamber),posX, posY, posZ, 0, "ONLY");  

    posX = 95.73-kNearFarLHC;
    posY = 33.26-kNearFarLHC;
    posZ = 2.0*kHzLateralSightAl+kHzVerticalCradleAl-fgkMotherThick2;               
    gMC->Gspos("SQ34",1,QuadrantNLayerName(chamber),posX, posY, posZ, 0, "ONLY");
    posZ = -1.0*posZ;              
    gMC->Gspos("SQ34",3,QuadrantFLayerName(chamber),posX, posY, posZ, 0, "ONLY");

// VertCradleB  2nd Trapezoid - 3 copies
    posX = 97.29+fgkDeltaQuadLHC;
    posY = 23.02+fgkDeltaQuadLHC;    
    posZ = 0.;              
    gMC->Gspos("SQ35",2,QuadrantMLayerName(chamber),posX, posY, posZ, 0, "ONLY");

    posX = 97.29-kNearFarLHC;
    posY = 23.02-kNearFarLHC;   
    posZ = 2.0*kHzLateralSightAl+kHzVerticalCradleAl-fgkMotherThick2;          
    gMC->Gspos("SQ35",1,QuadrantNLayerName(chamber),posX, posY, posZ, 0, "ONLY");    
    posZ = -1.0*posZ;          
    gMC->Gspos("SQ35",3,QuadrantFLayerName(chamber),posX, posY, posZ, 0, "ONLY");

// OutVertCradleC  3rd Trapeze - 3 copies
    posX = 98.31+fgkDeltaQuadLHC;
    posY = 12.77+fgkDeltaQuadLHC;  
    posZ = 0.;              
    gMC->Gspos("SQ36",2,QuadrantMLayerName(chamber),posX, posY, posZ, 0, "ONLY");

    posX = 98.05-kNearFarLHC;
    posY = 12.77-kNearFarLHC;        
    posZ = 2.0*kHzLateralSightAl+kHzVerticalCradleAl-fgkMotherThick2;         
           // Fix (2) of extrusion SQ36 from SQN1, SQN2, SQF1, SQF2 
	   // (was posX = 98.31 ...)
    gMC->Gspos("SQ36",1,QuadrantNLayerName(chamber),posX, posY, posZ, 0, "ONLY");       
    posZ = -1.0*posZ;
    gMC->Gspos("SQ36",3,QuadrantFLayerName(chamber),posX, posY, posZ, 0, "ONLY");  

// OutVertCradleD  4th Trapeze - 3 copies
    posX = 98.81+fgkDeltaQuadLHC;
    posY = 2.52+fgkDeltaQuadLHC;    
    posZ = 0.;              
    gMC->Gspos("SQ37",2,QuadrantMLayerName(chamber),posX, posY, posZ, 0, "ONLY");
   
    posZ = fgkMotherThick1-kHzVerticalCradleAl;                
    gMC->Gspos("SQ37",1,QuadrantMLayerName(chamber),posX, posY, posZ, 0, "ONLY");
    posZ = -1.0*posZ;          
    gMC->Gspos("SQ37",3,QuadrantMLayerName(chamber),posX, posY, posZ, 0, "ONLY");          
             
// LateralSightSupport - 2 copies
    posX = 98.33-kNearFarLHC;
    posY = 10.00-kNearFarLHC;    
    posZ = kHzLateralSightAl-fgkMotherThick2;
           // Fix (3) of extrusion SQ38 from SQN1, SQN2, SQF1, SQF2 
           // (was posX = 98.53 ...)
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
    posX = 2.0*kHxInVFrame+2.*kHxV1mm+kIAF+kHxInHFrame;
    posY = kHyInHFrame;
    posZ = 0.;       
    gMC->Gspos("SQ40",1,QuadrantMLayerName(chamber),posX, posY, posZ, 0, "ONLY"); 
 
 // keep memory of the mid position. Used for placing screws
    const GReal_t kMidHposX = posX;
    const GReal_t kMidHposY = posY;
    const GReal_t kMidHposZ = posZ;

// Flat 7.5mm horizontal section
    posX = 2.0*kHxInVFrame+2.*kHxV1mm+kIAF+kHxH1mm;
    posY = 2.0*kHyInHFrame+kHyH1mm;
    posZ = 0.;
    gMC->Gspos("SQ41",1,QuadrantMLayerName(chamber),posX, posY, posZ,0, "ONLY"); 
        
// InArcFrame 
    posX = 2.0*kHxInVFrame+2.*kHxV1mm;
    posY = 2.0*kHyInHFrame+2.*kHyH1mm;
    posZ = 0.;    
    gMC->Gspos("SQ42",1,QuadrantMLayerName(chamber),posX, posY, posZ,0, "ONLY"); 

// keep memory of the mid position. Used for placing screws
    const GReal_t kMidArcposX = posX;
    const GReal_t kMidArcposY = posY;
    const GReal_t kMidArcposZ = posZ;

// ScrewsInFrame - in sensitive volume

     Float_t scruX[64];
     Float_t scruY[64]; 
         
// Screws on IHEpoxyFrame

     const Int_t kNumberOfScrewsIH = 14;    // no. of screws on the IHEpoxyFrame
     const Float_t kOffX = 5.;              // inter-screw distance 

     // first screw coordinates 
     scruX[0] = 21.07;                  
     scruY[0] = -2.23; 
     // other screw coordinates      
     for (Int_t i = 1;i<kNumberOfScrewsIH;i++){   
     scruX[i] = scruX[i-1]+kOffX; 
     scruY[i] = scruY[0];
     }    
     // Position the volumes on the frames
     for (Int_t i = 0;i<kNumberOfScrewsIH;i++){
     posX = fgkDeltaQuadLHC + scruX[i];
     posY = fgkDeltaQuadLHC + scruY[i];
     posZ = 0.;   
     gMC->Gspos("SQ43",i+1,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ-kHzInHFrame-kSCRUHLE, 0, "ONLY");      
     if (chamber==1)
       gMC->Gspos("SQ44",i+1,"SQ40",posX+0.1-kMidHposX, posY+0.1-kMidHposY, posZ-kMidHposZ, 0, "ONLY");
     gMC->Gspos("SQ45",i+1,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ+kHzInHFrame+kSCRUNLE, 0, "ONLY"); 
     }
     // special screw coordinates
     scruX[63] = 16.3;  
     scruY[63] = -2.23; 
     posX = fgkDeltaQuadLHC + scruX[63];
     posY = fgkDeltaQuadLHC + scruY[63];
     posZ = 0.;            
     gMC->Gspos("SQ43",64,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ-kHzInHFrame-kSCRUHLE, 0, "ONLY");
     if (chamber==1)
       gMC->Gspos("SQ44",64,"SQ40",posX+0.1-kMidHposX, posY+0.1-kMidHposY, posZ-kMidHposZ, 0, "ONLY"); 
     gMC->Gspos("SQ45",64,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ+kHzInHFrame+kSCRUNLE, 0, "ONLY");  
     
// Screws on the IVEpoxyFrame
  
    const Int_t kNumberOfScrewsIV = 15;     // no. of screws on the IVEpoxyFrame
    const Float_t kOffY = 5.;               // inter-screw distance 
    Int_t firstScrew = 58;
    Int_t lastScrew = 44;
 
    // first (special) screw coordinates
    scruX[firstScrew-1] = -2.23; 
    scruY[firstScrew-1] = 16.3; 
    // second (repetitive) screw coordinates
    scruX[firstScrew-2] = -2.23; 
    scruY[firstScrew-2] = 21.07;     
    // other screw coordinates      
    for (Int_t i = firstScrew-3;i>lastScrew-2;i--){   
    scruX[i] = scruX[firstScrew-2];
    scruY[i] = scruY[i+1]+kOffY;
    }
    
    for (Int_t i = 0;i<kNumberOfScrewsIV;i++){
    posX = fgkDeltaQuadLHC + scruX[i+lastScrew-1];
    posY = fgkDeltaQuadLHC + scruY[i+lastScrew-1];
    posZ = 0.;       
    gMC->Gspos("SQ43",i+lastScrew,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ-kHzInHFrame-kSCRUHLE, 0, "ONLY");     
    if (chamber==1)
      gMC->Gspos("SQ44",i+lastScrew,"SQ00",posX+0.1-kMidVposX, posY+0.1-kMidVposY, posZ-kMidVposZ, 0, "ONLY"); 
    gMC->Gspos("SQ45",i+lastScrew,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ+kHzInHFrame+kSCRUNLE, 0, "ONLY");
    }    
    
// Screws on the OVEpoxyFrame
  
    const Int_t kNumberOfScrewsOV = 10;     // no. of screws on the OVEpoxyFrame

    firstScrew = 15;
    lastScrew = 25;
 
    // first (repetitive) screw coordinates
    // notes: 1st screw should be placed in volume 40 (InnerHorizFrame)
    scruX[firstScrew-1] = 90.9; 
    scruY[firstScrew-1] = -2.23;  // true value
 
    // other screw coordinates      
    for (Int_t i = firstScrew; i<lastScrew; i++ ){   
    scruX[i] = scruX[firstScrew-1];
    scruY[i] = scruY[i-1]+kOffY;
    }
    for (Int_t i = 1;i<kNumberOfScrewsOV;i++){
    posX = fgkDeltaQuadLHC + scruX[i+firstScrew-1];
    posY = fgkDeltaQuadLHC + scruY[i+firstScrew-1];
    posZ = 0.;   
    gMC->Gspos("SQ43",i+firstScrew,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ-kHzInHFrame-kSCRUHLE, 0, "ONLY");     
    // ??
    if (chamber==1)
      gMC->Gspos("SQ44",i+firstScrew,"SQ25",posX+0.1-kMidOVposX, posY+0.1-kMidOVposY, posZ-kMidOVposZ, 0, "ONLY"); 
    gMC->Gspos("SQ45",i+firstScrew,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ+kHzInHFrame+kSCRUNLE, 0, "ONLY"); 
    }
    // special case for 1st screw, inside the horizontal frame (volume 40)
    posX = fgkDeltaQuadLHC + scruX[firstScrew-1];
    posY = fgkDeltaQuadLHC + scruY[firstScrew-1];
    posZ = 0.;   
    if (chamber==1)
      gMC->Gspos("SQ44",firstScrew,"SQ40",posX+0.1-kMidHposX, posY+0.1-kMidHposY, posZ-kMidHposZ, 0, "ONLY"); 
          
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
    gMC->Gspos("SQ43",i+58+1,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ-kHzInHFrame-kSCRUHLE, 0, "ONLY");    
    if (chamber==1)
      gMC->Gspos("SQ44",i+58+1,"SQ42",posX+0.1-kMidArcposX, posY+0.1-kMidArcposY, posZ-kMidArcposZ, 0, "ONLY");
    gMC->Gspos("SQ45",i+58+1,QuadrantMLayerName(chamber),posX+0.1, posY+0.1, posZ+kHzInHFrame+kSCRUNLE, 0, "ONLY");
    }
}

//______________________________________________________________________________
void AliMUONSt1GeometryBuilderV2::PlaceInnerLayers(Int_t chamber)
{
/// Place the gas and copper layers for the specified chamber.

// Rotation Matrices 
  Int_t rot1, rot2, rot3, rot4;   

  fMUON->AliMatrix(rot1,  90., 315., 90.,  45., 0., 0.); // -45 deg
  fMUON->AliMatrix(rot2,  90.,  90., 90., 180., 0., 0.); //  90 deg
  fMUON->AliMatrix(rot3,  90., 270., 90.,   0., 0., 0.); // -90 deg 
  fMUON->AliMatrix(rot4,  90.,  45., 90., 135., 0., 0.); //  deg 

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
void AliMUONSt1GeometryBuilderV2::PlaceSector(AliMpSector* sector,SpecialMap specialMap, 
                            const TVector3& where, Bool_t reflectZ, Int_t chamber)
{
/// Place all the segments in the mother volume, at the position defined
/// by the sector's data.

/// \cond SKIP

  static Int_t segNum=1;
  Int_t sgn;
  Int_t reflZ;
  Int_t rotMat;

  if (!reflectZ) {
    sgn= 1;
    reflZ=0;                                     // no reflection along z... nothing
    fMUON->AliMatrix(rotMat,  90.,90.,90,180.,0.,0.);   // 90° rotation around z, NO reflection along z
  } else  {
    sgn=-1;
    fMUON->AliMatrix(reflZ,  90.,0.,90,90.,180.,0.);    // reflection along z
    fMUON->AliMatrix(rotMat,  90.,90.,90,180.,180.,0.); // 90° rotation around z AND reflection along z
  }
  
  GReal_t posX,posY,posZ;
  
#ifdef WITH_STL  
  vector<Int_t> alreadyDone;
#endif

#ifdef WITH_ROOT  
  TArrayI alreadyDone(20);
  Int_t nofAlreadyDone = 0;
#endif  

  for (Int_t irow=0;irow<sector->GetNofRows();irow++){ // for each row
    AliMpRow* row = sector->GetRow(irow);


    for (Int_t iseg=0;iseg<row->GetNofRowSegments();iseg++){ // for each row segment
      AliMpVRowSegment* seg = row->GetRowSegment(iseg);
      
#ifdef WITH_STL 
      SpecialMap::iterator iter 
        = specialMap.find(seg->GetMotifPositionId(0));

      if ( iter == specialMap.end()){ //if this is a normal segment (ie. not part of <specialMap>)
#endif  
      
#ifdef WITH_ROOT  
      Long_t value = specialMap.GetValue(seg->GetMotifPositionId(0));

      if ( value == 0 ){ //if this is a normal segment (ie. not part of <specialMap>)
#endif  
      
        // create the cathode part
        CreatePlaneSegment(segNum, seg->Dimensions(), seg->GetNofMotifs());
  
        posX = where.X() + seg->Position().X();
        posY = where.Y() + seg->Position().Y();
        posZ = where.Z() + sgn * (TotalHzPlane() + fgkHzGas + 2.*fgkHzPadPlane);
        gMC->Gspos(PlaneSegmentName(segNum).Data(), 1, 
	           QuadrantMLayerName(chamber), posX, posY, posZ, reflZ, "ONLY");

        // and place all the daughter boards of this segment
        for (Int_t motifNum=0;motifNum<seg->GetNofMotifs();motifNum++) {

	  // Copy number
          Int_t motifPosId = seg->GetMotifPositionId(motifNum);
          AliMpMotifPosition* motifPos = 
            sector->GetMotifMap()->FindMotifPosition(motifPosId);
	  Int_t copyNo = motifPosId;
	  if ( sector->GetDirection() == AliMp::kX) copyNo += fgkDaughterCopyNoOffset;
  
          // Position
          posX = where.X() + motifPos->Position().X() + fgkOffsetX;
          posY = where.Y() + motifPos->Position().Y() + fgkOffsetY;
	  posZ = where.Z() + sgn * (fgkMotherThick1 - TotalHzDaughter()); 

          gMC->Gspos(fgkDaughterName, copyNo, QuadrantMLayerName(chamber), posX, posY, posZ, reflZ, "ONLY");
        }  
        segNum++;
	
      } else { 

        // if this is a special segment	
        for (Int_t motifNum=0;motifNum<seg->GetNofMotifs();motifNum++) {// for each motif

          Int_t motifPosId = seg->GetMotifPositionId(motifNum);
          
#ifdef WITH_STL
          if (find(alreadyDone.begin(),alreadyDone.end(),motifPosId)
              != alreadyDone.end()) continue; // don't treat the same motif twice

          AliMUONSt1SpecialMotif spMot = specialMap[motifPosId];
#endif
#ifdef WITH_ROOT
          Bool_t isDone = false;
	  Int_t i=0;
	  while (i<nofAlreadyDone && !isDone) {
	    if (alreadyDone.At(i) == motifPosId) isDone=true;
	    i++;
	  }  
	  if (isDone) continue; // don't treat the same motif twice

          AliMUONSt1SpecialMotif spMot = *((AliMUONSt1SpecialMotif*)specialMap.GetValue(motifPosId));
#endif
          // check
	  // cout << chamber << " processing special motif: " << motifPosId << endl;  

          AliMpMotifPosition* motifPos = sector->GetMotifMap()->FindMotifPosition(motifPosId);

          // Copy number
	  Int_t copyNo = motifPosId;
	  if ( sector->GetDirection() == AliMp::kX) copyNo += fgkDaughterCopyNoOffset;

          // place the hole for the motif, wrt the requested rotation angle
          Int_t rot = ( spMot.GetRotAngle()<0.1 ) ? reflZ:rotMat;

          posX = where.X() + motifPos->Position().X() + spMot.GetDelta().X();
          posY = where.Y() + motifPos->Position().Y() + spMot.GetDelta().Y();
          posZ = where.Z() + sgn * (TotalHzPlane() + fgkHzGas + 2.*fgkHzPadPlane);
          gMC->Gspos(fgkHoleName, copyNo, QuadrantMLayerName(chamber), posX, posY, posZ, rot, "ONLY");

          // then place the daughter board for the motif, wrt the requested rotation angle
          posX = posX+fgkDeltaFilleEtamX;
          posY = posY+fgkDeltaFilleEtamY;
	  posZ = where.Z() + sgn * (fgkMotherThick1 - TotalHzDaughter()); 
          gMC->Gspos(fgkDaughterName, copyNo, QuadrantMLayerName(chamber), posX, posY, posZ, rot, "ONLY");

#ifdef WITH_STL
          alreadyDone.push_back(motifPosId);// mark this motif as done
#endif
#ifdef WITH_ROOT
          if (nofAlreadyDone == alreadyDone.GetSize()) 
	     alreadyDone.Set(2*nofAlreadyDone); 
          alreadyDone.AddAt(motifPosId, nofAlreadyDone++);       	  
#endif
          // check
	  // cout << chamber << " processed motifPosId: " << motifPosId << endl;
	}		
      }// end of special motif case
    }
  }
/// \endcond
} 

//______________________________________________________________________________
TString AliMUONSt1GeometryBuilderV2::GasVolumeName(const TString& name, Int_t chamber) const
{
/// Insert the chamber number into the name.

  TString newString(name);
 
  TString number(""); 
  number += chamber;

  newString.Insert(2, number);
  
  return newString;
}

//
// public methods
//

//______________________________________________________________________________
void AliMUONSt1GeometryBuilderV2::CreateMaterials()
{
/// Define materials specific to station 1

// Materials and medias defined in MUONv1:
//
//  AliMaterial( 9, "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
//  AliMaterial(10, "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
//  AliMaterial(15, "AIR$      ", 14.61, 7.3, .001205, 30423.24, 67500);
//  AliMixture( 19, "Bakelite$", abak, zbak, dbak, -3, wbak);
//  AliMixture( 20, "ArC4H10 GAS$", ag, zg, dg, 3, wg);
//  AliMixture( 21, "TRIG GAS$", atrig, ztrig, dtrig, -5, wtrig);
//  AliMixture( 22, "ArCO2 80%$", ag1, zg1, dg1, 3, wg1);
//  AliMixture( 23, "Ar-freon $", atr1, ztr1, dtr1, 4, wtr1);
//  AliMixture( 24, "ArCO2 GAS$", agas, zgas, dgas, 3, wgas);
//  AliMaterial(31, "COPPER$",   63.54,    29.,   8.96,  1.4, 0.);
//  AliMixture( 32, "Vetronite$",aglass, zglass, dglass,    5, wglass);
//  AliMaterial(33, "Carbon$",   12.01,     6.,  2.265, 18.8, 49.9);
//  AliMixture( 34, "Rohacell$", arohac, zrohac, drohac,   -4, wrohac); 

//  AliMedium( 1, "AIR_CH_US         ",  15, 1, iSXFLD, ...
//  AliMedium( 4, "ALU_CH_US          ",  9, 0, iSXFLD, ... 
//  AliMedium( 5, "ALU_CH_US          ", 10, 0, iSXFLD, ... 
//  AliMedium( 6, "AR_CH_US          ",  20, 1, iSXFLD, ... 
//  AliMedium( 7, "GAS_CH_TRIGGER    ",  21, 1, iSXFLD, ... 
//  AliMedium( 8, "BAKE_CH_TRIGGER   ",  19, 0, iSXFLD, ... 
//  AliMedium( 9, "ARG_CO2   ",          22, 1, iSXFLD, ... 
//  AliMedium(11, "PCB_COPPER        ",  31, 0, iSXFLD, ... 
//  AliMedium(12, "VETRONITE         ",  32, 0, iSXFLD, ... 
//  AliMedium(13, "CARBON            ",  33, 0, iSXFLD, ... 
//  AliMedium(14, "Rohacell          ",  34, 0, iSXFLD, ... 

  //
  // --- Define materials for GEANT ---
  //

  fMUON->AliMaterial(41, "Aluminium II$", 26.98, 13., 2.7, -8.9, 26.1);
       // was id: 9
       // from PDG and "The Particle Detector BriefBook", Bock and Vasilescu, P.18  
        // ??? same but the last but one argument < 0 
  //
  // --- Define mixtures for GEANT ---
  //

  //     Ar-CO2 gas II (80%+20%)
  Float_t ag1[2]   = { 39.95,  44.01};
  Float_t zg1[2]   = { 18., 22.};
  Float_t wg1[2]   = { .8, 0.2};
  Float_t dg1      = .001821;
  fMUON->AliMixture(45, "ArCO2 II 80%$", ag1, zg1, dg1, 2, wg1);  
            // was id: 22
            // use wg1 weighting factors (6th arg > 0)

  // Rohacell 51  II - imide methacrylique
  Float_t aRohacell51[4] = { 12.01, 1.01, 16.00, 14.01}; 
  Float_t zRohacell51[4] = { 6., 1., 8., 7.}; 
  Float_t wRohacell51[4] = { 9., 13., 2., 1.};  
  Float_t dRohacell51 = 0.052;
  fMUON->AliMixture(46, "FOAM$",aRohacell51,zRohacell51,dRohacell51,-4,wRohacell51);  
            // was id: 32
            // use relative A (molecular) values (6th arg < 0)
   
  Float_t aSnPb[2] = { 118.69, 207.19};
  Float_t zSnPb[2] = { 50, 82};
  Float_t wSnPb[2] = { 0.6, 0.4} ;
  Float_t dSnPb = 8.926;
  fMUON->AliMixture(47, "SnPb$", aSnPb,zSnPb,dSnPb,2,wSnPb);
            // was id: 35
            // use wSnPb weighting factors (6th arg > 0)

  // plastic definition from K5, Freiburg (found on web)
  Float_t aPlastic[2]={ 1.01, 12.01};
  Float_t zPlastic[2]={ 1, 6};
  Float_t wPlastic[2]={ 1, 1};
  Float_t denPlastic=1.107;
  fMUON->AliMixture(48, "Plastic$",aPlastic,zPlastic,denPlastic,-2,wPlastic);
            // was id: 33
            // use relative A (molecular) values (6th arg < 0)...no other info...
 
  // Not used, to be removed
  //
       // was id: 34

  // Inox/Stainless Steel (18%Cr, 9%Ni)
  Float_t aInox[3] = {55.847, 51.9961, 58.6934};  
  Float_t zInox[3] = {26., 24., 28.};
  Float_t wInox[3] = {0.73, 0.18, 0.09}; 
  Float_t denInox = 7.930;
  fMUON->AliMixture(50, "StainlessSteel$",aInox,zInox,denInox,3,wInox);   
            // was id: 37
            // use wInox weighting factors (6th arg > 0) 
            // from CERN note NUFACT Note023, Oct.2000 
  //
  // End - Not used, to be removed

  //
  // --- Define the tracking medias for GEANT ---
  // 

  GReal_t epsil  = .001;       // Tracking precision,
  //GReal_t stemax = -1.;        // Maximum displacement for multiple scat
  GReal_t tmaxfd = -20.;       // Maximum angle due to field deflection
  //GReal_t deemax = -.3;        // Maximum fractional energy loss, DLS
  GReal_t stmin  = -.8;
  GReal_t maxStepAlu   = fMUON->GetMaxStepAlu();
  GReal_t maxDestepAlu = fMUON->GetMaxDestepAlu();
  GReal_t maxStepGas   = fMUON->GetMaxStepGas();
  Int_t iSXFLD   = gAlice->Field()->PrecInteg();
  Float_t sXMGMX = gAlice->Field()->Max();

  fMUON->AliMedium(21, "ALU_II$",    41, 0, iSXFLD, sXMGMX, 
                   tmaxfd, maxStepAlu, maxDestepAlu, epsil, stmin);

		   // was med: 15  mat: 31 
  fMUON->AliMedium(24, "FrameCH$",   44, 1, iSXFLD, sXMGMX, 
                   10.0, 0.001, 0.001, 0.001, 0.001);
		   // was med: 20  mat: 36
  fMUON->AliMedium(25, "ARG_CO2_II", 45, 1, iSXFLD, sXMGMX,
                   tmaxfd, maxStepGas, maxDestepAlu, epsil, stmin);
		   // was med: 9   mat: 22
  fMUON->AliMedium(26, "FOAM_CH$",   46, 0, iSXFLD, sXMGMX,
                   10.0,  0.1, 0.1, 0.1, 0.1, 0, 0) ;
		   // was med: 16  mat: 32
  fMUON->AliMedium(27, "SnPb$",      47, 0, iSXFLD, sXMGMX,  
                   10.0, 0.01, 1.0, 0.003, 0.003);
		   // was med: 19  mat: 35
  fMUON->AliMedium(28, "Plastic$",   48, 0, iSXFLD, sXMGMX,
                   10.0, 0.01, 1.0, 0.003, 0.003);
		   // was med: 17  mat: 33

  // Not used, to be romoved
  //

  fMUON->AliMedium(30, "InoxBolts$", 50, 1, iSXFLD, sXMGMX, 
                   10.0, 0.01, 1.0, 0.003, 0.003);
		   // was med: 21  mat: 37
  //
  // End - Not used, to be removed
}

//______________________________________________________________________________
void AliMUONSt1GeometryBuilderV2::CreateGeometry()
{
/// Create the detailed GEANT geometry for the dimuon arm station1

  AliDebug(1,"Called");

  // Define chamber volumes as virtual
  // 

  // Create basic volumes
  // 
  CreateHole();
  CreateDaughterBoard();
  CreateInnerLayers();
  
  // Create reflexion matrices
  //
/*
  Int_t reflXZ, reflYZ, reflXY;
  fMUON->AliMatrix(reflXZ,  90.,  180., 90., 90., 180., 0.);
  fMUON->AliMatrix(reflYZ,  90., 0., 90.,-90., 180., 0.);
  fMUON->AliMatrix(reflXY,  90., 180., 90., 270., 0., 0.);
*/
  // Define transformations for each quadrant
  // In old coordinate system:        In new coordinate system:
  // 
  // 
  //     II. |  I.                   I. |  II. 
  //         |                    (101) | (100)
  //   _____ | ____               _____ | ____                         
  //         |                          |
  //    III. |  IV.                 IV. | III.
  //                              (102) | (103) 
  // 
/*
  Int_t rotm[4];
  rotm[0]=0;       // quadrant I
  rotm[1]=reflXZ;  // quadrant II
  rotm[2]=reflXY;  // quadrant III
  rotm[3]=reflYZ;  // quadrant IV
*/
  TGeoRotation rotm[4]; 
  rotm[0] = TGeoRotation("identity");
  rotm[1] = TGeoRotation("reflXZ", 90.,  180., 90., 90., 180., 0.);
  rotm[2] = TGeoRotation("reflXY", 90., 180., 90., 270., 0., 0.);
  rotm[3] = TGeoRotation("reflYZ", 90., 0., 90.,-90., 180., 0.);
  
  TVector3 scale[4];  
  scale[0] = TVector3( 1,  1,  1);  // quadrant I
  scale[1] = TVector3(-1,  1, -1);  // quadrant II
  scale[2] = TVector3(-1, -1,  1);  // quadrant III
  scale[3] = TVector3( 1, -1, -1);  // quadrant IV
  
  Int_t  detElemId[4];  
  detElemId[0] =  1;  // quadrant I
  detElemId[1] =  0;  // quadrant II
  detElemId[2] =  3;  // quadrant III
  detElemId[3] =  2;  // quadrant IV
  
  // Shift in Z of the middle layer
  Double_t deltaZ = 7.5/2.;         

  // Position of quadrant I wrt to the chamber position
  // TVector3 pos0(-fgkDeltaQuadLHC, -fgkDeltaQuadLHC, deltaZ);

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

      // DE envelope
      GReal_t posx0, posy0, posz0;
      posx0 = fgkPadXOffsetBP * scale[i].X();
      posy0 = fgkPadYOffsetBP * scale[i].Y();;
      posz0 = deltaZ * scale[i].Z();
      GetEnvelopes(ich-1)
        ->AddEnvelope(QuadrantEnvelopeName(ich,i), detElemId[i] + ich*100, true,
	              TGeoTranslation(posx0, posy0, posz0), rotm[i]);

      // Middle layer
      GReal_t posx, posy, posz;
      posx = -fgkDeltaQuadLHC - fgkPadXOffsetBP;
      posy = -fgkDeltaQuadLHC - fgkPadYOffsetBP;
      posz = 0.;
      GetEnvelopes(ich-1)
        ->AddEnvelopeConstituent(QuadrantMLayerName(ich), QuadrantEnvelopeName(ich,i),
	             i+1, TGeoTranslation(posx, posy, posz));

      // Near/far layers
      GReal_t  posx2 = posx + shiftXY;;
      GReal_t  posy2 = posy + shiftXY;;
      GReal_t  posz2 = posz - shiftZ;;
      //gMC->Gspos(QuadrantNLayerName(ich), i+1, "ALIC", posx2, posy2, posz2, rotm[i],"ONLY");
      GetEnvelopes(ich-1)
        ->AddEnvelopeConstituent(QuadrantNLayerName(ich), QuadrantEnvelopeName(ich,i),
	             i+1, TGeoTranslation(posx2, posy2, posz2)); 
    
      posz2 = posz + shiftZ;      
      //gMC->Gspos(QuadrantFLayerName(ich), i+1, "ALIC", posx2, posy2, posz2, rotm[i],"ONLY");
      GetEnvelopes(ich-1)
        ->AddEnvelopeConstituent(QuadrantFLayerName(ich), QuadrantEnvelopeName(ich,i), 
	             i+1, TGeoTranslation(posx2, posy2, posz2)); 
   }
 }     
}

//______________________________________________________________________________
void AliMUONSt1GeometryBuilderV2::SetTransformations() 
{
/// Define the transformations for the station2 chambers.

  if (gAlice->GetModule("SHIL")) {
    SetMotherVolume(0, "YOUT1");
    SetMotherVolume(1, "YOUT1");
  }  

  SetVolume(0, "SC01", true);
  SetVolume(1, "SC02", true);

  Double_t zpos1 = - AliMUONConstants::DefaultChamberZ(0); 
  SetTranslation(0, TGeoTranslation(0., 0., zpos1));

  Double_t zpos2 = - AliMUONConstants::DefaultChamberZ(1); 
  SetTranslation(1, TGeoTranslation(0., 0., zpos2));
}

//______________________________________________________________________________
void AliMUONSt1GeometryBuilderV2::SetSensitiveVolumes()
{
/// Define the sensitive volumes for station2 chambers.

  GetGeometry(0)->SetSensitiveVolume("SA1G");
  GetGeometry(0)->SetSensitiveVolume("SB1G");
  GetGeometry(0)->SetSensitiveVolume("SC1G");
  GetGeometry(0)->SetSensitiveVolume("SD1G");
  GetGeometry(0)->SetSensitiveVolume("SE1G");
  GetGeometry(0)->SetSensitiveVolume("SF1G");
  GetGeometry(0)->SetSensitiveVolume("SG1G");
  GetGeometry(0)->SetSensitiveVolume("SH1G");
  GetGeometry(0)->SetSensitiveVolume("SI1G");
  GetGeometry(0)->SetSensitiveVolume("SJ1G");
  GetGeometry(0)->SetSensitiveVolume("SK1G");
    
  GetGeometry(1)->SetSensitiveVolume("SA2G");
  GetGeometry(1)->SetSensitiveVolume("SB2G");
  GetGeometry(1)->SetSensitiveVolume("SC2G");
  GetGeometry(1)->SetSensitiveVolume("SD2G");
  GetGeometry(1)->SetSensitiveVolume("SE2G");
  GetGeometry(1)->SetSensitiveVolume("SF2G");
  GetGeometry(1)->SetSensitiveVolume("SG2G");
  GetGeometry(1)->SetSensitiveVolume("SH2G");
  GetGeometry(1)->SetSensitiveVolume("SI2G");
  GetGeometry(1)->SetSensitiveVolume("SJ2G");
  GetGeometry(1)->SetSensitiveVolume("SK2G");
}

