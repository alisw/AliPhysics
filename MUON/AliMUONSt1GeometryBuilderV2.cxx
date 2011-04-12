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
//-----------------------------------------------------------------------------
// Class AliMUONSt1GeometryBuilderV2
// ---------------------------------
// MUON Station1 detailed geometry construction class.
// (Originally defined in AliMUONv2.cxx - now removed.)
// Included in AliRoot 2004/01/23
// Authors: David Guez, Ivana Hrivnacova, Marion MacCormick; IPN Orsay
//-----------------------------------------------------------------------------

#include "AliMUONSt1GeometryBuilderV2.h"
#include "AliMUONSt1SpecialMotif.h"
#include "AliMUON.h"
#include "AliMUONConstants.h"
#include "AliMUONGeometryModule.h"
#include "AliMUONGeometryEnvelopeStore.h"

#include "AliMpSegmentation.h"
#include "AliMpDEManager.h"
#include "AliMpConstants.h"
#include "AliMpCDB.h"
#include "AliMpSector.h"
#include "AliMpRow.h"
#include "AliMpVRowSegment.h"
#include "AliMpMotifMap.h"
#include "AliMpMotifPosition.h"
#include "AliMpPlaneType.h"

#include "AliRun.h"
#include "AliMagF.h"
#include "AliLog.h"

#include <Riostream.h>
#include <TClonesArray.h>
#include <TGeoCompositeShape.h>
#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoTube.h>
#include <TGeoVolume.h>
#include <TGeoXtru.h>
#include <TSystem.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TVirtualMC.h>
#include <TArrayI.h>

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
const GReal_t AliMUONSt1GeometryBuilderV2::fgkDeltaFilleEtamX=1.00;
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
const char* AliMUONSt1GeometryBuilderV2::fgkQuadrantMFLayerName="SQMF";
const Int_t AliMUONSt1GeometryBuilderV2::fgkFoamBoxNameOffset=200; 
const Int_t AliMUONSt1GeometryBuilderV2::fgkFR4BoxNameOffset=400; 
const Int_t AliMUONSt1GeometryBuilderV2::fgkDaughterCopyNoOffset=1000;

//______________________________________________________________________________
AliMUONSt1GeometryBuilderV2::AliMUONSt1GeometryBuilderV2(AliMUON* muon)
  : AliMUONVGeometryBuilder(0, 2),
    fMUON(muon)
{
/// Standard constructor
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
/// The shape of the sensitive area is defined as an extruded
/// solid substracted with tube (to get inner circular shape). 

  TGeoMedium* kMedArCO2  = gGeoManager->GetMedium("MUON_ARG_CO2");
  TGeoMedium* kMedCopper = gGeoManager->GetMedium("MUON_COPPER_II");

  Double_t rmin = 0.0;
  Double_t rmax = fgkMotherIR1;
  Double_t hz   = fgkHzPadPlane + fgkHzGas;
  new TGeoTube("cutTube",rmin, rmax, hz); 

  Double_t maxXY = 89.0; 
  Double_t xy1   = 77.33;
  Double_t xy2   = 48.77;
  Double_t dxy1  = maxXY - xy1;
    
  Int_t nz = 2;
  Int_t nv = 6;
  Double_t vx[6] = {  0.0,   0.0,   xy2, maxXY, maxXY, dxy1 };
  Double_t vy[6] = { dxy1, maxXY, maxXY,   xy2,   0.0,  0.0 };

  TGeoXtru* xtruS1 = new TGeoXtru(nz);
  xtruS1->SetName("xtruS1");
  xtruS1->DefinePolygon(nv, vx, vy);
  xtruS1->DefineSection(0, -fgkHzGas,  0.0, 0.0, 1.0); 
  xtruS1->DefineSection(1,  fgkHzGas,  0.0, 0.0, 1.0); 
  TGeoCompositeShape* layerS1 = new TGeoCompositeShape("layerS1", "xtruS1-cutTube");
  new TGeoVolume("SA1G", layerS1, kMedArCO2 );
  
  TGeoXtru* xtruS2 = new TGeoXtru(nz);
  xtruS2->SetName("xtruS2");
  xtruS2->DefinePolygon(nv, vx, vy);
  xtruS2->DefineSection(0, -fgkHzGas,  0.0, 0.0, 1.0); 
  xtruS2->DefineSection(1,  fgkHzGas,  0.0, 0.0, 1.0); 
  TGeoCompositeShape* layerS2 = new TGeoCompositeShape("layerS2", "xtruS2-cutTube");
  new TGeoVolume("SA2G", layerS2, kMedArCO2 );

  TGeoXtru* xtruS3 = new TGeoXtru(nz);
  xtruS3->SetName("xtruS3");
  xtruS3->DefinePolygon(nv, vx, vy);
  xtruS3->DefineSection(0, -fgkHzPadPlane,  0.0, 0.0, 1.0); 
  xtruS3->DefineSection(1,  fgkHzPadPlane,  0.0, 0.0, 1.0); 
  TGeoCompositeShape* layerS3 = new TGeoCompositeShape("layerS3", "xtruS3-cutTube");
  new TGeoVolume("SA1C", layerS3, kMedCopper );
}  
  

//______________________________________________________________________________
void AliMUONSt1GeometryBuilderV2::CreateSpacer0()
{
/// The spacer volumes are defined according to the input prepared by Nicole Willis
/// without any modifications
///                                                                       <pre>
/// No.    Type  Material Center (mm)            Dimensions (mm) (half lengths)
///  5     BOX   EPOXY    408.2  430.4 522.41    5.75  1.5   25.5
///  5P    BOX   EPOXY    408.2  445.4 522.41    5.75  1.5   25.5
///  6     BOX   EPOXY    408.2  437.9 519.76    5.75  15.0   1.0
///  6P    BOX   EPOXY    408.2  437.9 525.06    5.75  15.0   1.0
///  7     CYL   INOX     408.2  437.9 522.41    r=3.0  hz=20.63
///                                                                      </pre>

  // tracking medias
  Int_t* idtmed = fMUON->GetIdtmed()->GetArray()-1099;
  Int_t idFrameEpoxy = idtmed[1123]; // medium 24 = Frame Epoxy ME730  // was 20 not 16
  Int_t idInox = idtmed[1128];       // medium 29 Stainless Steel (18%Cr,9%Ni,Fe) // was 21 not 17

  GReal_t par[3];
  par[0] = 0.575;
  par[1] = 0.150;
  par[2] = 2.550;
  gMC->Gsvolu("Spacer05","BOX",idFrameEpoxy,par,3);

  par[0] = 0.575;
  par[1] = 1.500;
  par[2] = 0.100;
  gMC->Gsvolu("Spacer06","BOX",idFrameEpoxy,par,3);

  par[0] = 0.000;
  par[1] = 0.300;
  par[2] = 2.063;
  gMC->Gsvolu("Spacer07","TUBE",idInox,par,3);
}  


//______________________________________________________________________________
void AliMUONSt1GeometryBuilderV2::CreateSpacer()
{
/// The spacer volumes are defined according to the input prepared by Nicole Willis
/// with modifications needed to fit into existing geometry.
///                                                                       <pre>
/// No.    Type  Material Center (mm)            Dimensions (mm) (half lengths)
///  5     BOX   EPOXY    408.2  430.4 522.41    5.75  1.5   25.5
///  5P    BOX   EPOXY    408.2  445.4 522.41    5.75  1.5   25.5
///  6     BOX   EPOXY    408.2  437.9 519.76    5.75  15.0   1.0
///  6P    BOX   EPOXY    408.2  437.9 525.06    5.75  15.0   1.0
///  7     CYL   INOX     408.2  437.9 522.41    r=3.0  hz=20.63
///                                                                      </pre>
/// To fit in existing volumes the volumes 5 and 7 are represented by 2 volumes
/// with half size in z (5A, &A); the dimensions of the volume 5A were also modified
/// to avoid overlaps (x made smaller, y larger to abotain the identical volume)   

  // tracking medias
  Int_t* idtmed = fMUON->GetIdtmed()->GetArray()-1099;
  Int_t idFrameEpoxy = idtmed[1123]; // medium 24 = Frame Epoxy ME730  // was 20 not 16
  Int_t idInox = idtmed[1128];       // medium 29 Stainless Steel (18%Cr,9%Ni,Fe) // was 21 not 17

  //GReal_t par[3];
  //par[0] = 0.575;
  //par[1] = 0.150;
  //par[2] = 2.550;
  //gMC->Gsvolu("Spacer5","BOX",idFrameEpoxy,par,3);

  GReal_t par[3];
  par[0] = 0.510;
  par[1] = 0.170;
  par[2] = 1.1515;
  gMC->Gsvolu("Spacer5A","BOX",idFrameEpoxy,par,3);

  par[0] = 0.510;
  par[1] = 1.500;
  par[2] = 0.100;
  gMC->Gsvolu("Spacer6","BOX",idFrameEpoxy,par,3);

  //par[0] = 0.000;
  //par[1] = 0.300;
  //par[2] = 2.063;
  //gMC->Gsvolu("Spacer7","TUBE",idInox,par,3);

  par[0] = 0.000;
  par[1] = 0.300;
  par[2] = 1.0315;
  gMC->Gsvolu("Spacer7A","TUBE",idInox,par,3);
}  

//______________________________________________________________________________
void AliMUONSt1GeometryBuilderV2::CreateQuadrant(Int_t chamber)
{
/// Create the quadrant (bending and non-bending planes)
/// for the given chamber

  // CreateQuadrantLayersAsVolumes(chamber);
  CreateQuadrantLayersAsAssemblies(chamber);

  CreateFrame(chamber);
  
  TExMap specialMap;
  specialMap.Add(76, (Long_t) new AliMUONSt1SpecialMotif(TVector2( 0.1, 0.72), 90.));
  specialMap.Add(75, (Long_t) new AliMUONSt1SpecialMotif(TVector2( 0.7, 0.36)));
  specialMap.Add(47, (Long_t) new AliMUONSt1SpecialMotif(TVector2(1.01, 0.36)));

  // Load mapping from OCDB
  if ( ! AliMpSegmentation::Instance() ) {
    AliFatal("Mapping has to be loaded first !");
  }
       
  const AliMpSector* kSector1 
    = AliMpSegmentation::Instance()->GetSector(100, AliMpDEManager::GetCathod(100, AliMp::kBendingPlane));
  if ( ! kSector1 ) {
    AliFatal("Could not access sector segmentation !");
  }

  //Bool_t reflectZ = true;
  Bool_t reflectZ = false;
  //TVector3 where = TVector3(2.5+0.1+0.56+0.001, 2.5+0.1+0.001, 0.);
  TVector3 where = TVector3(fgkDeltaQuadLHC + fgkPadXOffsetBP, 
                            fgkDeltaQuadLHC + fgkPadYOffsetBP, 0.);
  PlaceSector(kSector1, specialMap, where, reflectZ, chamber);
  
  Int_t nb = AliMpConstants::ManuMask(AliMp::kNonBendingPlane);
  TExMapIter it(&specialMap);
#if ROOT_SVN_REVISION >= 29598
  Long64_t key;
  Long64_t value;
#else
  Long_t key;
  Long_t value;
#endif  
  
  while ( it.Next(key,value) == kTRUE ) { 
    delete reinterpret_cast<AliMUONSt1SpecialMotif*>(value);
  }
  specialMap.Delete();
  specialMap.Add(76 | nb,(Long_t) new AliMUONSt1SpecialMotif(TVector2(1.01,0.51),90.));
  specialMap.Add(75 | nb,(Long_t) new AliMUONSt1SpecialMotif(TVector2(2.20,-0.08)));
  specialMap.Add(47 | nb,(Long_t) new AliMUONSt1SpecialMotif(TVector2(2.40,-1.11)));
  specialMap.Add(20 | nb,(Long_t) new AliMUONSt1SpecialMotif(TVector2(0.2 ,-0.08)));
  specialMap.Add(46 | nb,(Long_t) new AliMUONSt1SpecialMotif(TVector2(0.92 , 0.17)));
  specialMap.Add(74 | nb,(Long_t) new AliMUONSt1SpecialMotif(TVector2(0.405, -0.10)));  
      // Fix (7) - overlap of SQ42 with MCHL (after moving the whole sector
      // in the true position)   

  const AliMpSector* kSector2 
    = AliMpSegmentation::Instance()
          ->GetSector(100, AliMpDEManager::GetCathod(100, AliMp::kNonBendingPlane));
  if ( ! kSector2 ) {
    AliFatal("Could not access sector !");
  }

  //reflectZ = false;
  reflectZ = true;
  TVector2 offset = TVector2(kSector2->GetPositionX(), kSector2->GetPositionY());
  where = TVector3(where.X()+offset.X(), where.Y()+offset.Y(), 0.); 
      // Add the half-pad shift of the non-bending plane wrt bending plane
      // (The shift is defined in the mapping as sector offset)
      // Fix (4) - was TVector3(where.X()+0.63/2, ... - now it is -0.63/2
  PlaceSector(kSector2, specialMap, where, reflectZ, chamber);

  it.Reset();
  while ( it.Next(key,value) == kTRUE ) {
    delete reinterpret_cast<AliMUONSt1SpecialMotif*>(value);
  }
  specialMap.Delete();
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
  
  // Place spacer in the concrete plane segments:
  // S225 (in S025), S267 (in S067) in chamber1 and S309 (in S109). S351(in S151) 
  // in chamber2
  // The segments were found as those which caused overlaps when we placed
  // the spacer in global coordinates via PlaceSpacer0 
  //
  //    <posXYZ   X_Y_Z=" 12.6000;   0.75000;   0.0000"> <volume name="Spacer5A"/>
  //    <posXYZ   X_Y_Z=" 12.6000;  -0.75000;   0.0000"> <volume name="Spacer5A"/>
  //    <posXYZ   X_Y_Z=" 12.6000;   0.0000;    1.1515"> <volume name="Spacer6"/>
  //    <posXYZ   X_Y_Z=" 12.6000;   0.0000;    0.0000"> <volume name="Spacer7A"/>

  if ( FoamBoxName(segNumber) == "S225" || 
       FoamBoxName(segNumber) == "S267" ||
       FoamBoxName(segNumber) == "S309" ||
       FoamBoxName(segNumber) == "S351"    )  
  {
    GReal_t posX =  12.6;
    GReal_t posY =  0.75;
    GReal_t posZ = -0.1;
    if ( FoamBoxName(segNumber) == "S267" || 
         FoamBoxName(segNumber) == "S351" ) posY += fgkPadYOffsetBP;    
    gMC->Gspos("Spacer5A", 1, FoamBoxName(segNumber).Data(), posX, posY, posZ,0, "ONLY");

    posY = -0.75;
    if ( FoamBoxName(segNumber) == "S267" || 
         FoamBoxName(segNumber) == "S351" ) posY += fgkPadYOffsetBP;    
    gMC->Gspos("Spacer5A", 2, FoamBoxName(segNumber).Data(), posX, posY, posZ,0, "ONLY");

    posY = 0.0;
    posZ = 1.1515;
    if ( FoamBoxName(segNumber) == "S267" || 
         FoamBoxName(segNumber) == "S351" ) posY += fgkPadYOffsetBP;    
    gMC->Gspos("Spacer6",  1, FoamBoxName(segNumber).Data(), posX, posY, posZ,0, "ONLY");

    posY = 0.0;
    posZ = 0.0;
    if ( FoamBoxName(segNumber) == "S267" || 
         FoamBoxName(segNumber) == "S351" ) posY += fgkPadYOffsetBP;    
    gMC->Gspos("Spacer7A", 1, FoamBoxName(segNumber).Data(), posX, posY, posZ,0, "ONLY");
  }  

  for (Int_t holeNum=0;holeNum<nofHoles;holeNum++) {
    GReal_t posX = ((2.*holeNum+1.)/nofHoles-1.)*dimensions.X();
    GReal_t posY = 0.;
    GReal_t posZ = 0.;
 
    gMC->Gspos(fgkHoleName,holeNum+1,
               FoamBoxName(segNumber).Data(),posX,posY,posZ,0,"ONLY");
  }
}

//______________________________________________________________________________
void AliMUONSt1GeometryBuilderV2::CreateQuadrantLayersAsVolumes(Int_t chamber)
{
/// Create the three main layers as real volumes.
/// Not used anymore.

  // tracking medias
  Int_t* idtmed = fMUON->GetIdtmed()->GetArray()-1099;
  Int_t idAir  = idtmed[1100];       // medium 1

  Float_t par[11];
  Float_t posX,posY,posZ;

// Quadrant volume TUBS1, positioned at the end
  par[0] = fgkMotherIR1;
  par[1] = fgkMotherOR1; 
  par[2] = fgkMotherThick1;  
  par[3] = fgkMotherPhiL1; 
  par[4] = fgkMotherPhiU1;
  gMC->Gsvolu(QuadrantMLayerName(chamber),"TUBS",idAir,par,5);
  // gMC->Gsvolu(QuadrantMFLayerName(chamber),"TUBS",idAir,par,5);

// Replace the volume shape with a composite shape
// with substracted overlap with beam shield (YMOT)

  if ( gMC->IsRootGeometrySupported() ) { 

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

    TGeoVolume* malayer 
      = gGeoManager->FindVolumeFast(QuadrantMFLayerName(chamber));
    if ( !malayer ) {
      AliErrorStream() 
         << "Quadrant volume " << QuadrantMFLayerName(chamber) << " not found" 
	 << endl;
    }
    else {
      TGeoShape* quadrant = malayer->GetShape();
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
      malayer->SetShape(composite);
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
}  

//______________________________________________________________________________
void AliMUONSt1GeometryBuilderV2::CreateQuadrantLayersAsAssemblies(Int_t chamber)
{
/// Create the three main layers as assemblies

  gGeoManager->MakeVolumeAssembly(QuadrantMLayerName(chamber).Data()); 
  gGeoManager->MakeVolumeAssembly(QuadrantMFLayerName(chamber).Data()); 
  gGeoManager->MakeVolumeAssembly(QuadrantNLayerName(chamber).Data()); 
  gGeoManager->MakeVolumeAssembly(QuadrantFLayerName(chamber).Data()); 
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

  // TString quadrantMLayerName = QuadrantMLayerName(chamber);

  TString quadrantMLayerName = QuadrantMFLayerName(chamber);
  TString quadrantNLayerName = QuadrantNLayerName(chamber);
  TString quadrantFLayerName = QuadrantFLayerName(chamber);

  const Float_t kNearFarLHC=2.4;    // Near and Far TUBS Origin wrt LHC Origin

  // tracking medias
  Int_t* idtmed = fMUON->GetIdtmed()->GetArray()-1099;
  
  //Int_t idAir  = idtmed[1100];       // medium 1
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
  
  
    TGeoMedium* kMedEpoxy = gGeoManager->GetMedium("MUON_FrameCH$");
    TGeoMedium* kMedInox  = gGeoManager->GetMedium("MUON_Kapton");
    TGeoMedium* kMedAlu   = gGeoManager->GetMedium("MUON_ALU_II$");


// Rotation Matrices  
      Int_t rot1, rot2, rot3, rot4;    
      
//   Rotation matrices  
     fMUON->AliMatrix(rot1,  90.,  90., 90., 180.,  0., 0.); // +90 deg in x-y plane
     fMUON->AliMatrix(rot2,  90.,  45., 90., 135.,  0., 0.); // +45 deg in x-y plane 
     fMUON->AliMatrix(rot3,  90.,  45., 90., 315.,180., 0.); // +45 deg in x-y + rotation 180° around y
     fMUON->AliMatrix(rot4,  90., 315., 90.,  45.,  0., 0.); // -45 deg in x-y plane 

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
  
// TopFrameAnode parameters - 2 trapezoids, 2 layers
// (redefined with TGeoXtru shape)
  const Float_t kH1FAA = 8.7/2.;
  const Float_t kTl1FAB = 4.35/2.;
  const Float_t kTl1FAA = 7.75/2.;

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
// (redefined with TGeoXtru shape)
//
  const Float_t kH1OETF = 7.196/2.;       // common to all 4 trapezoids
  const Float_t kTl1OETF1 = 3.996/2.;     // Trapezoid 1
  const Float_t kTl1OETF2 = 3.75/2;       // Trapezoid 2
  const Float_t kTl1OETF3 = 3.01/2.;      // Trapezoid 3
  const Float_t kTl1OETF4 = 1.77/2.;      // Trapezoid 4


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
// (redefined with TGeoXtru shape)
//
  const Float_t kH1VC1 = 10.25/2.;  // all cradles
  const Float_t kBl1VC1 = 3.70/2.;  // VertCradleA
  const Float_t kBl1VC2 = 6.266/2.; // VertCradleB
  const Float_t kBl1VC3 = 7.75/2.;  // VertCradleC

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
            

    // Common declarations for TGeoXtru parameters
    Double_t dx, dx0, dx1, dx2, dx3; 
    Double_t dy, dy1, dy2, dy3, dy4;
    Double_t vx[16];
    Double_t vy[16];
    Int_t nz;
    Int_t nv;

    // SQ04to06 and SQ05to07

    dx  =  2.*kH1FAA; 
    dy1 =  2.*kTl1FAA;
    dy2 =  2.*kTl1FAB;
    
    nz =  2;
    nv = 5;
    vx[0]  =   0.0;  vy[0]  =  0.0;
    vx[1]  =   0.0;  vy[1]  =  dy1;
    vx[2]  =    dx;  vy[2]  =  dy2;
    vx[3]  =  2*dx;  vy[3]  =  0.0;
    vx[4]  =    dx;  vy[4]  =  0.0;

    // Shift center in the middle
    for ( Int_t i=0; i<nv; i++ ) { 
      vx[i] -= dx;
      vy[i] -= 0.5*dy1;
    }  
  
    TGeoXtru* xtruS5 = new TGeoXtru(nz);
    xtruS5->DefinePolygon(nv, vx, vy);
    xtruS5->DefineSection(0, -kHzOuterFrameEpoxy,  0.0, 0.0, 1.0); 
    xtruS5->DefineSection(1,  kHzOuterFrameEpoxy,  0.0, 0.0, 1.0); 
    new TGeoVolume("SQ04toSQ06", xtruS5, kMedEpoxy);

    TGeoXtru* xtruS6 = new TGeoXtru(nz);
    xtruS6->DefinePolygon(nv, vx, vy);
    xtruS6->DefineSection(0, -kHzOuterFrameInox,  0.0, 0.0, 1.0); 
    xtruS6->DefineSection(1,  kHzOuterFrameInox,  0.0, 0.0, 1.0); 
    new TGeoVolume("SQ05toSQ07", xtruS6, kMedInox);


    // TopAnode1 -  layer 1 of 2
    par[0] = kHxTA1;
    par[1] = kHyTA1;
    par[2] = kHzTA11;    
    gMC->Gsvolu("SQ08","BOX",idInox,par,3); 
    
    // TopAnode1 -  layer 2 of 2
    par[2] = kHzTA12;    
    gMC->Gsvolu("SQ09","BOX",idFR4,par,3); 

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
// (redefined with TGeoXtru shape )
//---

    dx  = 2.*kH1OETF;
    dy1 = 2.*kTl1OETF4;
    dy2 = 2.*kTl1OETF3; 
    dy3 = 2.*kTl1OETF2;
    dy4 = 2.*kTl1OETF1;
    
    nz =  2;
    nv = 16;
    vx[0]  = -4*dx;  vy[0]  = 0.0;
    vx[1]  = -3*dx;  vy[1]  = dy1;
    vx[2]  = -2*dx;  vy[2]  = dy2;
    vx[3]  = -1*dx;  vy[3]  = dy3;
    vx[4]  =   0.0;  vy[4]  = dy4;
    vx[5]  =    dx;  vy[5]  = dy3; 
    vx[6]  =  2*dx;  vy[6]  = dy2;
    vx[7]  =  3*dx;  vy[7]  = dy1;
    vx[8]  =  4*dx;  vy[8]  = 0.0;
    vx[9]  =  3*dx;  vy[9]  = 0.0;
    vx[10] =  2*dx;  vy[10] = 0.0;
    vx[11] =    dx;  vy[11] = 0.0;
    vx[12] =   0.0;  vy[12] = 0.0;
    vx[13] = -1*dx;  vy[13] = 0.0;
    vx[14] = -2*dx;  vy[14] = 0.0;
    vx[15] = -3*dx;  vy[15] = 0.0;

    // Shift center in the middle
    for ( Int_t i=0; i<nv; i++ ) vy[i] += dy4/2.0;
  
    TGeoXtru* xtruS1 = new TGeoXtru(nz);
    xtruS1->DefinePolygon(nv, vx, vy);
    xtruS1->DefineSection(0, -kHzOuterFrameEpoxy,  0.0, 0.0, 1.0); 
    xtruS1->DefineSection(1,  kHzOuterFrameEpoxy,  0.0, 0.0, 1.0); 
    new TGeoVolume("SQ17to23", xtruS1, kMedEpoxy );

    TGeoXtru* xtruS2 = new TGeoXtru(nz);
    xtruS2->DefinePolygon(nv, vx, vy);
    xtruS2->DefineSection(0, -kHzOuterFrameInox,  0.0, 0.0, 1.0); 
    xtruS2->DefineSection(1,  kHzOuterFrameInox,  0.0, 0.0, 1.0); 
    new TGeoVolume("SQ18to24", xtruS2, kMedInox );

//
// OutEdgeTrapFrame Epoxy = (4 trapezes)*2 copies*2 layers (Epoxy/Inox)
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

    dy  = 2.*kH1VC1;
    dx0 = 2.*kBl1VC4;
    dx1 = 2.*kBl1VC3;
    dx2 = 2.*kBl1VC2; 
    dx3 = 2.*kBl1VC1;   
    
    // VertCradle
    // (Trapezoids SQ34 to SQ36 or SQ37 redefined with TGeoXtru shape)

    nz =  2;
    nv = 7;
    vx[0]  =   0.0;  vy[0]  =  0.0;
    vx[1]  =   0.0;  vy[1]  =   dy;
    vx[2]  =   0.0;  vy[2]  = 2*dy;
    vx[3]  =   0.0;  vy[3]  = 3*dy;
    vx[4]  =   dx3;  vy[4]  = 2*dy;
    vx[5]  =   dx2;  vy[5]  =   dy; 
    vx[6]  =   dx1;  vy[6]  =  0.0;

    // Shift center in the middle
    for ( Int_t i=0; i<nv; i++ ) { 
      vx[i] -= dx1/2.0;
      vy[i] -= 1.5*dy;
    }  
  
    TGeoXtru* xtruS3 = new TGeoXtru(nz);
    xtruS3->DefinePolygon(nv, vx, vy);
    xtruS3->DefineSection(0, -kHzVerticalCradleAl,  0.0, 0.0, 1.0); 
    xtruS3->DefineSection(1,  kHzVerticalCradleAl,  0.0, 0.0, 1.0); 
    new TGeoVolume("SQ34to36", xtruS3, kMedAlu);

    // Trapezoids SQ34 to SQ37;
    // (keeping the same coordinate system as for SQ34to36)

    nz =  2;
    nv = 9;
    vx[0]  =   0.0;  vy[0]  =-1.0*dy;
    vx[1]  =   0.0;  vy[1]  =  0.0;
    vx[2]  =   0.0;  vy[2]  =   dy;
    vx[3]  =   0.0;  vy[3]  = 2*dy;
    vx[4]  =   0.0;  vy[4]  = 3*dy;
    vx[5]  =   dx3;  vy[5]  = 2*dy;
    vx[6]  =   dx2;  vy[6]  =   dy; 
    vx[7]  =   dx1;  vy[7]  =  0.0;
    vx[8]  =   dx0;  vy[8]  =-1.0*dy;

    // Shift center in the middle (of SQ34to36!!)
    for ( Int_t i=0; i<nv; i++ ) { 
      vx[i] -= dx1/2.0;
      vy[i] -= 1.5*dy;
    }  
  
    TGeoXtru* xtruS4 = new TGeoXtru(nz);
    xtruS4->DefinePolygon(nv, vx, vy);
    xtruS4->DefineSection(0, -kHzVerticalCradleAl,  0.0, 0.0, 1.0); 
    xtruS4->DefineSection(1,  kHzVerticalCradleAl,  0.0, 0.0, 1.0); 
    new TGeoVolume("SQ34to37", xtruS4, kMedAlu);

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
    gMC->Gspos("SQ00",1,quadrantMLayerName,posX, posY, posZ, 0, "ONLY"); 

// keep memory of the mid position. Used for placing screws
    const GReal_t kMidVposX = posX;
    const GReal_t kMidVposY = posY;
    const GReal_t kMidVposZ = posZ;

    //Flat 7.5mm vertical section
    posX = 2.0*kHxInVFrame+kHxV1mm;
    posY = 2.0*kHyInHFrame+2.*kHyH1mm+kIAF+kHyV1mm;
    posZ = 0.;
    gMC->Gspos("SQ01",1,quadrantMLayerName,posX, posY, posZ,0, "ONLY"); 
    
    // TopFrameAnode place 2 layers of TopFrameAnode cuboids  
    posX = kHxTFA;
    posY = 2.*kHyInHFrame+2.*kHyH1mm+kIAF+2.*kHyInVFrame+kHyTFA;   
    posZ = -kHzOuterFrameInox;
    gMC->Gspos("SQ02",1,quadrantMLayerName,posX, posY, posZ,0,"ONLY"); 
    posZ = kHzOuterFrameEpoxy;
    gMC->Gspos("SQ03",1,quadrantMLayerName,posX, posY, posZ,0,"ONLY");
    
    // TopFrameAnode - place 2 layers of 2 trapezoids 
    // (SQ04 - SQ07)
    posX += kHxTFA + 2.*kH1FAA;
    posZ = -kHzOuterFrameInox; 
    gMC->Gspos("SQ04toSQ06",1,quadrantMLayerName,posX, posY, posZ, 0,"ONLY");
    posZ = kHzOuterFrameEpoxy;
    gMC->Gspos("SQ05toSQ07",1,quadrantMLayerName,posX, posY, posZ, 0,"ONLY");

    // TopAnode1 place 2 layers  
    posX = 6.8+fgkDeltaQuadLHC;
    posY = 99.85+fgkDeltaQuadLHC;
    posZ = -1.*kHzAnodeFR4;
    gMC->Gspos("SQ08",1,quadrantMLayerName,posX, posY, posZ, 0,"ONLY");  
    posZ = kHzTopAnodeSteel1;
    gMC->Gspos("SQ09",1,quadrantMLayerName,posX, posY, posZ, 0,"ONLY");    
         
    // TopAnode2 place 2 layers
    posX = 18.534+fgkDeltaQuadLHC;
    posY = 99.482+fgkDeltaQuadLHC; 
    posZ = -1.*kHzAnodeFR4;    
    // shift up to solve overlap with SQ14
    posY += 0.1;
    gMC->Gspos("SQ10",1,quadrantMLayerName,posX, posY, posZ, rot1,"ONLY");
    posZ = kHzTopAnodeSteel2;    
    gMC->Gspos("SQ11",1,quadrantMLayerName,posX, posY, posZ, rot1,"ONLY");       
    
    // TopAnode3 place 1 layer
    posX = 25.804+fgkDeltaQuadLHC;
    posY = 98.61+fgkDeltaQuadLHC;
    posZ = 0.;    
    gMC->Gspos("SQ12",1,quadrantMLayerName,posX, posY, posZ, rot1,"ONLY");  
          
    // TopEarthFace - 2 copies
    posX = 23.122+fgkDeltaQuadLHC;
    posY = 96.90+fgkDeltaQuadLHC;
    posZ = kHzOuterFrameEpoxy+kHzOuterFrameInox+kHzTopEarthFaceCu;
    gMC->Gspos("SQ13",1,quadrantMLayerName,posX, posY, posZ, 0,"ONLY");
    posZ = -1.*posZ;
    gMC->Gspos("SQ13",2,quadrantMLayerName,posX, posY, posZ, 0,"ONLY");

    // TopEarthProfile 
    posX = 14.475+fgkDeltaQuadLHC;
    posY = 97.900+fgkDeltaQuadLHC; 
    posZ = kHzTopEarthProfileCu;
    gMC->Gspos("SQ14",1,quadrantMLayerName,posX, posY, posZ, 0,"ONLY");
    posZ = -1.0*posZ;
    gMC->Gspos("SQ14",2,quadrantMLayerName,posX, posY, posZ, 0,"ONLY");

    // TopGasSupport - 2 copies                            
    posX = 4.9500+fgkDeltaQuadLHC;
    posY = 96.200+fgkDeltaQuadLHC;
    posZ = kHzOuterFrameEpoxy+kHzOuterFrameInox+kHzTopGasSupportAl;
    gMC->Gspos("SQ15",1,quadrantMLayerName,posX, posY, posZ, 0,"ONLY");
    posZ = -1.*posZ;
    gMC->Gspos("SQ15",2,quadrantMLayerName,posX, posY, posZ, 0,"ONLY");
    
    // TopPositioner parameters - single Stainless Steel trapezoid - 2 copies
    posX = 7.60+fgkDeltaQuadLHC;
    posY = 98.98+fgkDeltaQuadLHC;   
    posZ = kHzOuterFrameEpoxy+kHzOuterFrameInox+2.*kHzTopGasSupportAl+kHzTopPositionerSteel;
    gMC->Gspos("SQ16",1,quadrantMLayerName,posX, posY, posZ, rot1,"ONLY");
    posZ = -1.*posZ;
    gMC->Gspos("SQ16",2,quadrantMLayerName,posX, posY, posZ, rot1,"ONLY"); 

    // OutEdgeFrame 

    posZ = -1.0*kHzOuterFrameInox;     
    //Double_t xCenterAll = 70.6615;
    Double_t xCenterAll = 70.500;
    Double_t yCenterAll = 70.350;
    gMC->Gspos("SQ17to23",1,quadrantMLayerName, xCenterAll, yCenterAll, posZ, rot4,"ONLY");
     
    posZ = kHzOuterFrameEpoxy;
    gMC->Gspos("SQ18to24",1,quadrantMLayerName, xCenterAll, yCenterAll, posZ, rot4,"ONLY");
    
//---    
        
// OutVFrame
    posX = 2.*kHxInVFrame+kIAF+2.*kHxInHFrame-kHxOutVFrame+2.*kHxV1mm;
    posY = 2.*kHyInHFrame+kHyOutVFrame;    
    posZ = 0.;              
    gMC->Gspos("SQ25",1,quadrantMLayerName,posX, posY, posZ, 0, "ONLY"); 

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
    // shift to solve overlap with SQ17to23 and SQ18to24
    posX += 0.02;
    gMC->Gspos("SQ26",1,quadrantMLayerName,posX, posY, posZ, rot1,"ONLY"); 

// VertEarthFaceCu - 2 copies
    posX = 89.4000+fgkDeltaQuadLHC;
    posY = 25.79+fgkDeltaQuadLHC;    
    posZ = kHzFrameThickness+2.0*kHzFoam+kHzVertEarthFaceCu;              
    gMC->Gspos("SQ27",1,quadrantMLayerName,posX, posY, posZ, rot1, "ONLY"); 
    posZ = -1.0*posZ; 
    gMC->Gspos("SQ27",2,quadrantMLayerName,posX, posY, posZ, rot1, "ONLY"); 
    
// VertEarthSteel - 2 copies
    posX = 91.00+fgkDeltaQuadLHC;
    posY = 30.616+fgkDeltaQuadLHC;    
    posZ = kHzFrameThickness+2.0*kHzFoam+kHzVertBarSteel;              
    gMC->Gspos("SQ28",1,quadrantMLayerName,posX, posY, posZ, rot1, "ONLY"); 
    posZ = -1.0*posZ;              
    gMC->Gspos("SQ28",2,quadrantMLayerName,posX, posY, posZ, rot1, "ONLY");
 
// VertEarthProfCu - 2 copies
    posX = 92.000+fgkDeltaQuadLHC;
    posY = 29.64+fgkDeltaQuadLHC;    
    posZ = kHzFrameThickness;              
    gMC->Gspos("SQ29",1,quadrantMLayerName,posX, posY, posZ, rot1, "ONLY"); 
    posZ = -1.0*posZ;    
    gMC->Gspos("SQ29",2,quadrantMLayerName,posX, posY, posZ, rot1, "ONLY"); 

// SuppLateralPositionner - 2 copies 
    posX = 90.2-kNearFarLHC;
    posY = 5.00-kNearFarLHC;    
    posZ = kHzLateralPosnAl-fgkMotherThick2;             
    gMC->Gspos("SQ30",1,quadrantFLayerName,posX, posY, posZ, 0, "ONLY"); 
    posZ = -1.0*posZ;            
    gMC->Gspos("SQ30",2,quadrantNLayerName,posX, posY, posZ, 0, "ONLY"); 

// LateralPositionner - 2 copies - Face view
    posX = 92.175-kNearFarLHC-2.*kHxLPP;
    posY = 5.00-kNearFarLHC;   
    posZ =2.0*kHzLateralPosnAl+kHzLateralPosnInoxFace-fgkMotherThick2;              
    gMC->Gspos("SQ31",1,quadrantFLayerName,posX, posY, posZ, 0, "ONLY"); 
    posZ = -1.0*posZ;             
    gMC->Gspos("SQ31",2,quadrantNLayerName,posX, posY, posZ, 0, "ONLY"); 

// LateralPositionner -  Profile view   
    posX = 92.175+fgkDeltaQuadLHC+kHxLPF-kHxLPP;
    posY = 5.00+fgkDeltaQuadLHC;    
    posZ = 0.;              
    gMC->Gspos("SQ32",1,quadrantMLayerName,posX, posY, posZ, 0, "ONLY"); // middle layer

    posX = 92.175-kNearFarLHC+kHxLPF-kHxLPP; 
    posY = 5.0000-kNearFarLHC;    
    posZ = fgkMotherThick2-kHzLPNF;              
    gMC->Gspos("SQ33",1,quadrantNLayerName,posX, posY, posZ, 0, "ONLY"); // near layer
    posZ = -1.*posZ;
    gMC->Gspos("SQ33",2,quadrantFLayerName,posX, posY, posZ, 0, "ONLY"); // far layer
      

// VertCradle - 3 (or 4 ) trapezoids redefined with TGeoXtru shape

    posX = 97.29+fgkDeltaQuadLHC;
    posY = 23.02+fgkDeltaQuadLHC;    
    posZ = 0.;          
    posX += 1.39311;
    gMC->Gspos("SQ34to37",2,quadrantMLayerName,posX, posY, posZ, 0, "ONLY");  

    posX = 97.29-kNearFarLHC;
    posY = 23.02-kNearFarLHC;   
    posZ = 2.0*kHzLateralSightAl+kHzVerticalCradleAl-fgkMotherThick2;          
    posX += 1.39311;
    gMC->Gspos("SQ34to36",1,quadrantNLayerName,posX, posY, posZ, 0, "ONLY");

    posZ = -1.0*posZ;              
    gMC->Gspos("SQ34to36",3,quadrantFLayerName,posX, posY, posZ, 0, "ONLY");


// OutVertCradleD  4th Trapeze - 3 copies

    posX = 98.81+fgkDeltaQuadLHC;
    posY = 2.52+fgkDeltaQuadLHC;    
    posZ = fgkMotherThick1-kHzVerticalCradleAl;                
    gMC->Gspos("SQ37",1,quadrantMLayerName,posX, posY, posZ, 0, "ONLY");
    posZ = -1.0*posZ;          
    gMC->Gspos("SQ37",3,quadrantMLayerName,posX, posY, posZ, 0, "ONLY");          
             
// LateralSightSupport - 2 copies
    posX = 98.33-kNearFarLHC;
    posY = 10.00-kNearFarLHC;    
    posZ = kHzLateralSightAl-fgkMotherThick2;
           // Fix (3) of extrusion SQ38 from SQN1, SQN2, SQF1, SQF2 
           // (was posX = 98.53 ...)
    gMC->Gspos("SQ38",1,quadrantNLayerName,posX, posY, posZ, 0, "ONLY"); 
    posZ = -1.0*posZ;             
    gMC->Gspos("SQ38",2,quadrantFLayerName,posX, posY, posZ, 0, "ONLY"); 
    
// Mire placement
    posX = 92.84+fgkDeltaQuadLHC;  
    posY = 8.13+fgkDeltaQuadLHC;
    posZ = 0.;
    gMC->Gspos("SQ39",1,quadrantMLayerName,posX, posY, posZ, 0,"ONLY");    

//---

// InHFrame
    posX = 2.0*kHxInVFrame+2.*kHxV1mm+kIAF+kHxInHFrame;
    posY = kHyInHFrame;
    posZ = 0.;       
    gMC->Gspos("SQ40",1,quadrantMLayerName,posX, posY, posZ, 0, "ONLY"); 
 
 // keep memory of the mid position. Used for placing screws
    const GReal_t kMidHposX = posX;
    const GReal_t kMidHposY = posY;
    const GReal_t kMidHposZ = posZ;

// Flat 7.5mm horizontal section
    posX = 2.0*kHxInVFrame+2.*kHxV1mm+kIAF+kHxH1mm;
    posY = 2.0*kHyInHFrame+kHyH1mm;
    posZ = 0.;
    gMC->Gspos("SQ41",1,quadrantMLayerName,posX, posY, posZ,0, "ONLY"); 
        
// InArcFrame 
    posX = 2.0*kHxInVFrame+2.*kHxV1mm;
    posY = 2.0*kHyInHFrame+2.*kHyH1mm;
    posZ = 0.;    
    gMC->Gspos("SQ42",1,quadrantMLayerName,posX, posY, posZ,0, "ONLY"); 

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
     gMC->Gspos("SQ43",i+1,quadrantMLayerName,posX+0.1, posY+0.1, posZ-kHzInHFrame-kSCRUHLE, 0, "ONLY");      
     if (chamber==1)
       gMC->Gspos("SQ44",i+1,"SQ40",posX+0.1-kMidHposX, posY+0.1-kMidHposY, posZ-kMidHposZ, 0, "ONLY");
     gMC->Gspos("SQ45",i+1,quadrantMLayerName,posX+0.1, posY+0.1, posZ+kHzInHFrame+kSCRUNLE, 0, "ONLY"); 
     }
     // special screw coordinates
     scruX[63] = 16.3;  
     scruY[63] = -2.23; 
     posX = fgkDeltaQuadLHC + scruX[63];
     posY = fgkDeltaQuadLHC + scruY[63];
     posZ = 0.;            
     gMC->Gspos("SQ43",64,quadrantMLayerName,posX+0.1, posY+0.1, posZ-kHzInHFrame-kSCRUHLE, 0, "ONLY");
     if (chamber==1)
       gMC->Gspos("SQ44",64,"SQ40",posX+0.1-kMidHposX, posY+0.1-kMidHposY, posZ-kMidHposZ, 0, "ONLY"); 
     gMC->Gspos("SQ45",64,quadrantMLayerName,posX+0.1, posY+0.1, posZ+kHzInHFrame+kSCRUNLE, 0, "ONLY");  
     
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
    gMC->Gspos("SQ43",i+lastScrew,quadrantMLayerName,posX+0.1, posY+0.1, posZ-kHzInHFrame-kSCRUHLE, 0, "ONLY");     
    if (chamber==1)
      gMC->Gspos("SQ44",i+lastScrew,"SQ00",posX+0.1-kMidVposX, posY+0.1-kMidVposY, posZ-kMidVposZ, 0, "ONLY"); 
    gMC->Gspos("SQ45",i+lastScrew,quadrantMLayerName,posX+0.1, posY+0.1, posZ+kHzInHFrame+kSCRUNLE, 0, "ONLY");
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
    gMC->Gspos("SQ43",i+firstScrew,quadrantMLayerName,posX+0.1, posY+0.1, posZ-kHzInHFrame-kSCRUHLE, 0, "ONLY");     
    // ??
    if (chamber==1)
      gMC->Gspos("SQ44",i+firstScrew,"SQ25",posX+0.1-kMidOVposX, posY+0.1-kMidOVposY, posZ-kMidOVposZ, 0, "ONLY"); 
    gMC->Gspos("SQ45",i+firstScrew,quadrantMLayerName,posX+0.1, posY+0.1, posZ+kHzInHFrame+kSCRUNLE, 0, "ONLY"); 
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
    gMC->Gspos("SQ43",i+58+1,quadrantMLayerName,posX+0.1, posY+0.1, posZ-kHzInHFrame-kSCRUHLE, 0, "ONLY");    
    if (chamber==1)
      gMC->Gspos("SQ44",i+58+1,"SQ42",posX+0.1-kMidArcposX, posY+0.1-kMidArcposY, posZ-kMidArcposZ, 0, "ONLY");
    gMC->Gspos("SQ45",i+58+1,quadrantMLayerName,posX+0.1, posY+0.1, posZ+kHzInHFrame+kSCRUNLE, 0, "ONLY");
    }
}
//______________________________________________________________________________
void AliMUONSt1GeometryBuilderV2::PlaceInnerLayers(Int_t chamber)
{
/// Place the gas and copper layers for the specified chamber.

   GReal_t x  = fgkDeltaQuadLHC;
   GReal_t y  = fgkDeltaQuadLHC;
   GReal_t zg = 0.0;
   GReal_t zc = fgkHzGas + fgkHzPadPlane;
   Int_t dpos = (chamber-1)*2;
 
   TString name = GasVolumeName("SAG", chamber);
   gMC->Gspos(name,1,QuadrantMLayerName(chamber),x,y,zg,0,"ONLY");
   gMC->Gspos("SA1C", 1+dpos, QuadrantMLayerName(chamber),x,y, zc,0,"ONLY");
   gMC->Gspos("SA1C", 2+dpos, QuadrantMLayerName(chamber),x,y,-zc,0,"ONLY");
}

//______________________________________________________________________________
void AliMUONSt1GeometryBuilderV2::PlaceSpacer0(Int_t chamber)
{
/// Place the spacer defined in global positions
/// !! This method should be used only to find out the right mother volume
/// for the spacer if geometry is changed and the plane segment volumes
/// will change their numbering

  // Global position of mother volume for the QuadrantMLayer
  // SQM1: (-2.6, -2.6, -522.41)
  // SQM2: (-2.6, -2.6, -541.49)
  GReal_t mx =  2.6;
  GReal_t my = -2.6;
  GReal_t mz =  522.41;
  
  GReal_t x, y, z;
  x = 40.82  - mx;
  y = 43.04  - my;
  z = 522.41 - mz;
  AliDebugStream(2) << "spacer05 pos1: " << x << ", " << y << ", " << z << endl;
  gMC->Gspos("Spacer05", 1, QuadrantMLayerName(chamber), x, y, z, 0, "ONLY");

  y = 44.54  - my;
  AliDebugStream(2) << "spacer05 pos2: " << x << ", " << y << ", " << z << endl;
  gMC->Gspos("Spacer05", 2, QuadrantMLayerName(chamber), x, y, z, 0, "ONLY");

  x = 40.82  - mx;
  y = 43.79  - my;
  z = 519.76 - mz;
  AliDebugStream(2) << "spacer06 pos1: " << x << ", " << y << ", " << z << endl;
  gMC->Gspos("Spacer06", 1, QuadrantMLayerName(chamber), x, y, z, 0, "ONLY");

  z = 525.06 - mz;
  AliDebugStream(2) << "spacer06 pos2: " << x << ", " << y << ", " << z << endl;
  gMC->Gspos("Spacer06", 2, QuadrantMLayerName(chamber), x, y, z, 0, "ONLY");

  x = 40.82  - mx;
  y = 43.79  - my;
  z = 522.41 - mz;
  AliDebugStream(2) << "spacer07 pos1: " << x << ", " << y << ", " << z << endl;
  gMC->Gspos("Spacer07", 1, QuadrantMLayerName(chamber), x, y, z, 0, "ONLY");
}

//______________________________________________________________________________
void AliMUONSt1GeometryBuilderV2::PlaceSector(const AliMpSector* sector,
                            TExMap specialMap, 
                            const TVector3& where, Bool_t reflectZ, Int_t chamber)
{
/// Place all the segments in the mother volume, at the position defined
/// by the sector's data.                                                      \n
/// The lines with comments COMMENT OUT BEGIN/END indicates blocks
/// which can be commented out in order to reduce the number of volumes
/// in a sector to the plane segments corresponding to regular motifs only.

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
  
  TArrayI alreadyDone(20);
  Int_t nofAlreadyDone = 0;

  for (Int_t irow=0;irow<sector->GetNofRows();irow++){ // for each row
    AliMpRow* row = sector->GetRow(irow);


    for (Int_t iseg=0;iseg<row->GetNofRowSegments();iseg++){ // for each row segment
      AliMpVRowSegment* seg = row->GetRowSegment(iseg);
      
      Long_t value = specialMap.GetValue(seg->GetMotifPositionId(0));

      if ( value == 0 ){ //if this is a normal segment (ie. not part of <specialMap>)
      
        // create the cathode part
        CreatePlaneSegment(segNum, TVector2(seg->GetDimensionX(),seg->GetDimensionY()), 
                           seg->GetNofMotifs());
  
        posX = where.X() + seg->GetPositionX();
        posY = where.Y() + seg->GetPositionY();
        posZ = where.Z() + sgn * (TotalHzPlane() + fgkHzGas + 2.*fgkHzPadPlane);
        gMC->Gspos(PlaneSegmentName(segNum).Data(), 1, 
	           QuadrantMLayerName(chamber), posX, posY, posZ, reflZ, "ONLY");

        // and place all the daughter boards of this segment

// COMMENT OUT BEGIN
        for (Int_t motifNum=0;motifNum<seg->GetNofMotifs();motifNum++) {

	  // Copy number
          Int_t motifPosId = seg->GetMotifPositionId(motifNum);
          AliMpMotifPosition* motifPos = 
            sector->GetMotifMap()->FindMotifPosition(motifPosId);
	  Int_t copyNo = motifPosId;
	  if ( sector->GetDirection() == AliMp::kX) copyNo += fgkDaughterCopyNoOffset;
  
          // Position
          posX = where.X() + motifPos->GetPositionX() + fgkOffsetX;
          posY = where.Y() + motifPos->GetPositionY() + fgkOffsetY;
	  posZ = where.Z() + sgn * (fgkMotherThick1 - TotalHzDaughter()); 
          gMC->Gspos(fgkDaughterName, copyNo, QuadrantMLayerName(chamber), posX, posY, posZ, reflZ, "ONLY");
        }  
// COMMENT OUT END

        segNum++;
	
      } else { 

// COMMENT OUT BEGIN
        // if this is a special segment	
        for (Int_t motifNum=0;motifNum<seg->GetNofMotifs();motifNum++) {// for each motif

          Int_t motifPosId = seg->GetMotifPositionId(motifNum);
          
          Bool_t isDone = false;
	  Int_t i=0;
	  while (i<nofAlreadyDone && !isDone) {
	    if (alreadyDone.At(i) == motifPosId) isDone=true;
	    i++;
	  }  
	  if (isDone) continue; // don't treat the same motif twice

          AliMUONSt1SpecialMotif spMot = *((AliMUONSt1SpecialMotif*)specialMap.GetValue(motifPosId));
	  AliDebugStream(2) << chamber << " processing special motif: " << motifPosId << endl;  

          AliMpMotifPosition* motifPos = sector->GetMotifMap()->FindMotifPosition(motifPosId);

          // Copy number
	  Int_t copyNo = motifPosId;
	  if ( sector->GetDirection() == AliMp::kX) copyNo += fgkDaughterCopyNoOffset;

          // place the hole for the motif, wrt the requested rotation angle
          Int_t rot = ( spMot.GetRotAngle()<0.1 ) ? reflZ:rotMat;

          posX = where.X() + motifPos->GetPositionX() + spMot.GetDelta().X();
          posY = where.Y() + motifPos->GetPositionY() + spMot.GetDelta().Y();
          posZ = where.Z() + sgn * (TotalHzPlane() + fgkHzGas + 2.*fgkHzPadPlane);
          // Shift the hole for special motif 46 to avoid debording into S047
          if ( copyNo == 2070 ) {
            posX -= 0.1;
            posY -= 0.1;
          }  
          gMC->Gspos(fgkHoleName, copyNo, QuadrantMLayerName(chamber), posX, posY, posZ, rot, "ONLY");

          // then place the daughter board for the motif, wrt the requested rotation angle
          posX = posX+fgkDeltaFilleEtamX;
          posY = posY+fgkDeltaFilleEtamY;
          // Do not shift the daughter board
          if ( copyNo == 2070 ) {
            posX += 0.1;
            posY += 0.1;
          }  
	  posZ = where.Z() + sgn * (fgkMotherThick1 - TotalHzDaughter()); 
          gMC->Gspos(fgkDaughterName, copyNo, QuadrantMLayerName(chamber), posX, posY, posZ, rot, "ONLY");

          if (nofAlreadyDone == alreadyDone.GetSize()) 
	     alreadyDone.Set(2*nofAlreadyDone); 
          alreadyDone.AddAt(motifPosId, nofAlreadyDone++);       	  

	  AliDebugStream(2) << chamber << " processed motifPosId: " << motifPosId << endl;
	}		
// COMMENT OUT END
 
      }// end of special motif case
    }
  }
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
//  AliMedium(24, "FrameCH$          ",  44, 1, iSXFLD, ...

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

 //  //     Ar-CO2 gas II (80%+20%)   
//   Float_t ag1[2]   = { 39.95,  44.01};
//   Float_t zg1[2]   = { 18., 22.};
//   Float_t wg1[2]   = { .8, 0.2};
//   Float_t dg1      = .001821;
//   fMUON->AliMixture(45, "ArCO2 II 80%$", ag1, zg1, dg1, 2, wg1);  
//             // was id: 22
//             // use wg1 weighting factors (6th arg > 0)

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
  // GReal_t maxStepGas   = fMUON->GetMaxStepGas();
  Int_t iSXFLD   = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->PrecInteg();
  Float_t sXMGMX = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max();

  fMUON->AliMedium(21, "ALU_II$",    41, 0, iSXFLD, sXMGMX, 
                   tmaxfd, maxStepAlu, maxDestepAlu, epsil, stmin);

		   // was med: 20  mat: 36
 //  fMUON->AliMedium(25, "ARG_CO2_II", 45, 1, iSXFLD, sXMGMX,
//                    tmaxfd, maxStepGas, maxDestepAlu, epsil, stmin);
// 		   // was med: 9   mat: 22
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
  // CreateSpacer0();
  CreateSpacer();
  
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
  scale[0] = TVector3( 1,  1, -1);  // quadrant I
  scale[1] = TVector3(-1,  1,  1);  // quadrant II
  scale[2] = TVector3(-1, -1, -1);  // quadrant III
  scale[3] = TVector3( 1, -1,  1);  // quadrant IV
  
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
  //for (Int_t ich=1; ich<2; ich++) {

    // Create quadrant volume
    CreateQuadrant(ich);

    // Place gas volumes
    PlaceInnerLayers(ich);
    
    // Place the quadrant
    for (Int_t i=0; i<4; i++) {
    //for (Int_t i=1; i<2; i++) {
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
      GetEnvelopes(ich-1)
        ->AddEnvelopeConstituent(QuadrantMFLayerName(ich), QuadrantEnvelopeName(ich,i),
	             i+5, TGeoTranslation(posx, posy, posz));

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

      // Place spacer in global coordinates in the first non rotated quadrant
      // if ( detElemId[i] == 0 ) PlaceSpacer0(ich);
               // !! This placement should be used only to find out the right mother volume
               // for the spacer if geometry is changed and the plane segment volumes
               // will change their numbering
               // The call to the method CreateSpacer0(); above haa to be uncommented, too
   }
 }     
}

//______________________________________________________________________________
void AliMUONSt1GeometryBuilderV2::SetVolumes() 
{
/// Define the volumes for the station2 chambers.

  if (gAlice->GetModule("SHIL")) {
    SetMotherVolume(0, "YOUT1");
    SetMotherVolume(1, "YOUT1");
  }  

  SetVolume(0, "SC01", true);
  SetVolume(1, "SC02", true);
}

//______________________________________________________________________________
void AliMUONSt1GeometryBuilderV2::SetTransformations() 
{
/// Define the transformations for the station2 chambers.

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
  GetGeometry(1)->SetSensitiveVolume("SA2G");
}

