
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


//-------------------------------------------------------------------------
//  Beam pipe class for ALICE MFT upgrade
//  This version uses TGeo
//  Authors:
//  F. Manso 
//  A. Morsch
//  R. Tieulent
//-------------------------------------------------------------------------


#include <Riostream.h>

#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoVolume.h>
#include <TGeoTorus.h>
#include <TGeoTube.h>
#include <TGeoCone.h>
#include <TGeoPcon.h>
#include <TGeoBBox.h>
#include <TGeoXtru.h>
#include <TGeoCompositeShape.h>
#include <TGeoGlobalMagField.h>

#include "AliConst.h"
#include "AliMagF.h"
#include "AliPIPEv4.h"
#include "AliRun.h"
#include "AliLog.h"

ClassImp(AliPIPEv4)

//_____________________________________________________________________________
AliPIPEv4::AliPIPEv4(): AliPIPE()              

{
  // Constructor
}

//_____________________________________________________________________________
AliPIPEv4::AliPIPEv4(const char *name, const char *title)
  : AliPIPE(name,title)
{
  // Constructor
}


 
//___________________________________________
void AliPIPEv4::CreateGeometry()
{
  AliDebug(1,"Create PIPEv4 geometry");
  //
  //  Class describing the beam pipe geometry
  //
  Float_t z, zsh, z0;
  //
  // Rotation Matrices
  //
  const Float_t  kDegRad = TMath::Pi() / 180.;
  // Rotation by 180 deg
  TGeoRotation* rot180        = new TGeoRotation("rot180", 90., 180.,  90.,  90., 180.,   0.);
  TGeoRotation* rotyz         = new TGeoRotation("rotyz",  90., 180.,   0., 180.,  90.,  90.);
  TGeoRotation* rotxz         = new TGeoRotation("rotxz",   0.,   0.,  90.,  90.,  90., 180.);
  //TGeoRotation* rot045        = new TGeoRotation("rot045", 90.,  45.,  90., 135.,   0.,   0.);
  //TGeoRotation* rot135        = new TGeoRotation("rot135", 90. ,135.,  90., 225.,   0.,   0.);
  //TGeoRotation* rot225        = new TGeoRotation("rot225", 90. ,225.,  90., 315.,   0.,   0.);
  //TGeoRotation* rot315        = new TGeoRotation("rot315", 90. ,315.,  90.,  45.,   0.,   0.);    
  //
  // Media
  //const TGeoMedium* kMedSi      =  gGeoManager->GetMedium("PIPE_SILICON"); //FM
  const TGeoMedium* kMedAir     =  gGeoManager->GetMedium("PIPE_AIR");
  const TGeoMedium* kMedAirHigh =  gGeoManager->GetMedium("PIPE_AIR_HIGH");
  const TGeoMedium* kMedVac     =  gGeoManager->GetMedium("PIPE_VACUUM");    
  const TGeoMedium* kMedInsu    =  gGeoManager->GetMedium("PIPE_INS_C0");    
  const TGeoMedium* kMedSteel   =  gGeoManager->GetMedium("PIPE_INOX");        
  const TGeoMedium* kMedBe      =  gGeoManager->GetMedium("PIPE_BE");       
  const TGeoMedium* kMedCu      =  gGeoManager->GetMedium("PIPE_CU");        
  //const TGeoMedium* kMedKapton  =  gGeoManager->GetMedium("PIPE_KAPTON");        
  //const TGeoMedium* kMedAco     =  gGeoManager->GetMedium("PIPE_ANTICORODAL");        
  //const TGeoMedium* kMedNEG     =  gGeoManager->GetMedium("PIPE_NEG COATING"); 
  //const TGeoMedium* kMedAlu     =  gGeoManager->GetMedium("PIPE_ALU");    // fm       
  const TGeoMedium* kMedAlu2219 =  gGeoManager->GetMedium("PIPE_AA2219");   // fm     
  const TGeoMedium* kMedRohacell =  gGeoManager->GetMedium("PIPE_ROHACELL");
  const TGeoMedium* kMedPolyimide =  gGeoManager->GetMedium("PIPE_POLYIMIDE");
  const TGeoMedium* kMedCarbonFiber =  gGeoManager->GetMedium("PIPE_M55J6K");

  // Top volume
  TGeoVolume* top    = gGeoManager->GetVolume("ALIC");
  //
  //
  ////////////////////////////////////////////////////////////////////////////////     
  //                                                                            //
  //                                  The Central Vacuum system                 // 
  //                                                                            //
  ////////////////////////////////////////////////////////////////////////////////
  //
  //
  //  The ALICE central beam-pipe according to drawing         LHCVC2C_0001 
  //  Drawings of sub-elements:
  //  
  //  Pos 7 - Minimised Flange:                                LHCVFX_P0025
  //  Pos 6 - Standard Flange:                                 STDVFUHV0009
  //  Pos 8 - Bellow:                                          LHCVBX__0001
  //
  //  Absolute z-coordinates -82.0 - 400.0 cm 
  //  Total length:                                          482.0 cm
  //  It consists of 3 main parts:
  //  CP/2 The flange on the non-absorber side:               36.5 cm  
  //  CP/1 The central Be pipe:                              405.0 cm 
  //  CP/3 The double-bellow and flange on the absorber side: 40.5 cm 

  //
  /*
  //  Starting position in z
  const Float_t kCPz0      = -400.0;
  //  Length of the CP/1 section
  const Float_t kCP1Length =  405.0;    
  //  Length of the CP/2 section    
  const Float_t kCP2Length =   36.5;
  //  Length of the CP/3 section    
  const Float_t kCP3Length =   40.5;
  //  Position of the CP/2 section    
  //    const Float_t kCP2pos    = kCPz0 + kCP2Length / 2.;
  //  Position of the CP/3 section
  const Float_t kCP3pos    = kCPz0 + kCP2Length + kCP1Length + kCP3Length/2.;
  */


  //////////////////// NEW BEAM PIPE GEOMETRY FOR MuonForwardTracker ,
  // Authors: F. Manso, R. Tieulent 
  // Drawings from C. Gargiulo :
  // \\cern.ch\dfs\Workspaces\c\cgargiul\EXPERIMENT\ALICE\ALICE_MECHANICS\ALICE_DATA_PACKAGE\IN\DETECTORS\ITS_UPGRADE\1-DESIGN\3D_cad_model\R14_20140311_ALI\
  //
  //------------------- Pipe version 4.7 March 2014 -----------------------------

  TGeoVolumeAssembly * beamPipeCsideSection = new TGeoVolumeAssembly("BeamPipeCsideSection");
  
  Float_t fBeryliumSectionOuterRadius = 1.9;
  Float_t fBeryliumSectionZmax        =  44.4;
  Float_t fBeryliumSectionZmin        = -44.4;
  Float_t fBeryliumSectionThickness   = 0.08;

  Float_t fBellowSectionOuterRadius   = 2.15;
  Float_t fCSideBPSOuterRadius        = 2.22;
  Float_t fCSideBPSWallThickness      = 0.15;
  Float_t fBellowSectionZmax          = -55.35;
  Float_t fBellowOuterRadius          = 2.8;
  Float_t fFirstConeAngle             = 15. * TMath::DegToRad();
  Float_t fChangeThicknessAngle       = 45. * TMath::DegToRad();
  //  Float_t fCSideBPSLength             = 3.53;
  Float_t fCSideBPSLength             = 3.53+1.52;
  Float_t fDzFirstCone = (fCSideBPSOuterRadius - fBeryliumSectionOuterRadius) / TMath::Tan(fFirstConeAngle);
  //  Float_t fReduceThicknessPartAfterBPSLength = 1.52;
  Float_t fReduceThicknessPartAfterBPSLength = 0.;
  Float_t fThinPartBeforeBellowLength = 1.025;

  
  Float_t fDistanceBetweenBellows = 2.5;
  
  
  Float_t fAdaptConeZmax = -77.43;
  Float_t fAdaptConeZmin = -80.6;
  Float_t fAdaptConeRmax = 3.0;
  Float_t fFlangeRmax = 4.3;
  Float_t fFlangeLength = 1.4;
  
  Float_t fBellowPlieRadius = 0.17;  // radius of bellow plies
  Float_t fBellowPlieThickness = 0.03;  // Thickness of bellow plies 300 microns
  Int_t fNBellowConvolutions = 7;
  
  

  Float_t fZ1 = fBeryliumSectionZmin; // z of Be - Al jonction on the C-side
  Float_t fZ2 = fBellowSectionZmax +fDzFirstCone ; // z of end of small diameter part (beginning of first cone before the bellow
  Float_t fZ3 = fBellowSectionZmax +(fCSideBPSOuterRadius - fBellowSectionOuterRadius) / TMath::Tan(fFirstConeAngle); // z of End of first cone part with 0.8mm thickness
  Float_t fZ4 = fBellowSectionZmax; // z of End of first Cone
  Float_t fZ5 = fBellowSectionZmax - fCSideBPSLength; // z of End of Beam Pipe support section
  Float_t fZ6 = fBellowSectionZmax - fCSideBPSLength - (fCSideBPSOuterRadius-fBellowSectionOuterRadius) / TMath::Tan(fChangeThicknessAngle); // z of End of Beam Pipe support section after reduction of thickness
  Float_t fZ7 = fZ6 - fReduceThicknessPartAfterBPSLength ; // Z of end of 800 microns section after Beam Pipe Support
  Float_t fZ8 = fZ7 - (fBeryliumSectionThickness-fBellowPlieThickness) / TMath::Tan(fChangeThicknessAngle);
  Float_t fZ9 = fZ7 - fThinPartBeforeBellowLength; // Z of the start of first bellow
  Float_t fFirstBellowZmax = fZ9;
  
  //---------------- Be pipe around the IP ----------
  TGeoPcon* berylliumTube = new TGeoPcon(0., 360., 2);
  berylliumTube->DefineSection(0,fBeryliumSectionZmax,fBeryliumSectionOuterRadius-fBeryliumSectionThickness,fBeryliumSectionOuterRadius);
  berylliumTube->DefineSection(1,fBeryliumSectionZmin,fBeryliumSectionOuterRadius-fBeryliumSectionThickness,fBeryliumSectionOuterRadius);
  TGeoVolume* voberylliumTube = new TGeoVolume("berylliumTube",berylliumTube,kMedBe);
  voberylliumTube->SetLineColor(kRed);
  beamPipeCsideSection->AddNode(voberylliumTube,1,new TGeoTranslation(0., 0., 0.));

  TGeoPcon* berylliumTubeVacuum = new TGeoPcon(0., 360., 2);
  berylliumTubeVacuum->DefineSection(0,fBeryliumSectionZmax, 0.,fBeryliumSectionOuterRadius-fBeryliumSectionThickness);
  berylliumTubeVacuum->DefineSection(1,fBeryliumSectionZmin, 0.,fBeryliumSectionOuterRadius-fBeryliumSectionThickness);
  TGeoVolume* voberylliumTubeVacuum = new TGeoVolume("berylliumTubeVacuum",berylliumTubeVacuum,kMedVac);
  voberylliumTubeVacuum->SetVisibility(0);voberylliumTubeVacuum->SetLineColor(kGreen);
  beamPipeCsideSection->AddNode(voberylliumTubeVacuum,1,new TGeoTranslation(0., 0., 0.));
  //-------------------------------------------------


  //----------------  Al tube ------------------
  TGeoPcon* aluBeforeBellows = new TGeoPcon(0., 360., 9);
  aluBeforeBellows->DefineSection(0,fZ1, fBeryliumSectionOuterRadius-fBeryliumSectionThickness,fBeryliumSectionOuterRadius);
  aluBeforeBellows->DefineSection(1,fZ2,fBeryliumSectionOuterRadius-fBeryliumSectionThickness,fBeryliumSectionOuterRadius);
  aluBeforeBellows->DefineSection(2,fZ3,fBellowSectionOuterRadius-fBeryliumSectionThickness,fBellowSectionOuterRadius);
  aluBeforeBellows->DefineSection(3,fZ4,fCSideBPSOuterRadius-fCSideBPSWallThickness,fCSideBPSOuterRadius);
  aluBeforeBellows->DefineSection(4,fZ5,fCSideBPSOuterRadius-fCSideBPSWallThickness,fCSideBPSOuterRadius);
  aluBeforeBellows->DefineSection(5,fZ6,fBellowSectionOuterRadius-fBeryliumSectionThickness,fBellowSectionOuterRadius);
  aluBeforeBellows->DefineSection(6,fZ7,fBellowSectionOuterRadius-fBeryliumSectionThickness,fBellowSectionOuterRadius);
  aluBeforeBellows->DefineSection(7,fZ8,fBellowSectionOuterRadius-fBeryliumSectionThickness,fBellowSectionOuterRadius-fBeryliumSectionThickness+fBellowPlieThickness);
  aluBeforeBellows->DefineSection(8,fZ9,fBellowSectionOuterRadius-fBeryliumSectionThickness,fBellowSectionOuterRadius-fBeryliumSectionThickness+fBellowPlieThickness);
  TGeoVolume* voaluBeforeBellows = new TGeoVolume("aluBeforeBellows",aluBeforeBellows,kMedAlu2219);
  voaluBeforeBellows->SetLineColor(kBlue);
  beamPipeCsideSection->AddNode(voaluBeforeBellows,1,new TGeoTranslation(0., 0., 0.));

  TGeoPcon* aluBeforeBellowsVacuum = new TGeoPcon(0., 360., 7);
  aluBeforeBellowsVacuum->DefineSection(0,fZ1,0.,fBeryliumSectionOuterRadius-fBeryliumSectionThickness);
  aluBeforeBellowsVacuum->DefineSection(1,fZ2,0.,fBeryliumSectionOuterRadius-fBeryliumSectionThickness);
  aluBeforeBellowsVacuum->DefineSection(2,fZ3,0.,fBellowSectionOuterRadius-fBeryliumSectionThickness);
  aluBeforeBellowsVacuum->DefineSection(3,fZ4,0.,fCSideBPSOuterRadius-fCSideBPSWallThickness);
  aluBeforeBellowsVacuum->DefineSection(4,fZ5,0.,fCSideBPSOuterRadius-fCSideBPSWallThickness);
  aluBeforeBellowsVacuum->DefineSection(5,fZ6,0.,fBellowSectionOuterRadius-fBeryliumSectionThickness);
  aluBeforeBellowsVacuum->DefineSection(6,fZ9,0.,fBellowSectionOuterRadius-fBeryliumSectionThickness);
  TGeoVolume* voaluBeforeBellowsVacuum = new TGeoVolume("aluBeforeBellowsVacuum",aluBeforeBellowsVacuum,kMedVac);
  voaluBeforeBellowsVacuum->SetVisibility(0);voaluBeforeBellowsVacuum->SetLineColor(kGreen);
  beamPipeCsideSection->AddNode(voaluBeforeBellowsVacuum,1,new TGeoTranslation(0., 0., 0.));
  //-------------------------------------------------

  
  Float_t fBellowLength = fNBellowConvolutions * (4.*fBellowPlieRadius - 2. *fBellowPlieThickness);
  // ------------------ First Bellow  --------------------
  TGeoVolume* vobellows1 = MakeBellowCside("bellows1", fNBellowConvolutions, fBellowSectionOuterRadius-fBeryliumSectionThickness, fBellowOuterRadius, fBellowPlieRadius ,fBellowPlieThickness);
  beamPipeCsideSection->AddNode(vobellows1, 1, new TGeoTranslation(0., 0., fFirstBellowZmax-fBellowLength/2. - 2.*fBellowPlieRadius));
  //------------------------------------------------------

  Float_t fZ10 = fFirstBellowZmax - fBellowLength; // End of First bellow
  Float_t fZ12 = fZ10 - fThinPartBeforeBellowLength;
  Float_t fZ11 = fZ12 +  (fBeryliumSectionThickness-fBellowPlieThickness) / TMath::Tan(fChangeThicknessAngle); // End of 300 microns thickness part after first bellow
  Float_t fZ13 = fZ12  - fDistanceBetweenBellows;
  Float_t fZ14 = fZ13 -(fBeryliumSectionThickness-fBellowPlieThickness) / TMath::Tan(fChangeThicknessAngle);
  Float_t fZ15 = fZ14 -fThinPartBeforeBellowLength;
  Float_t fSecondBellowZmax = fZ15;

  
  //---------- Al tube between the bellows ----------
  TGeoPcon* tube4 = new TGeoPcon(0., 360., 6);
  tube4->DefineSection(0,fZ10, fBellowSectionOuterRadius-fBeryliumSectionThickness,fBellowSectionOuterRadius-fBeryliumSectionThickness+fBellowPlieThickness);
  tube4->DefineSection(1,fZ11,fBellowSectionOuterRadius-fBeryliumSectionThickness,fBellowSectionOuterRadius-fBeryliumSectionThickness+fBellowPlieThickness);
  tube4->DefineSection(2,fZ12,fBellowSectionOuterRadius-fBeryliumSectionThickness,fBellowSectionOuterRadius);
  tube4->DefineSection(3,fZ13,fBellowSectionOuterRadius-fBeryliumSectionThickness,fBellowSectionOuterRadius);
  tube4->DefineSection(4,fZ14,fBellowSectionOuterRadius-fBeryliumSectionThickness,fBellowSectionOuterRadius-fBeryliumSectionThickness+fBellowPlieThickness);
  tube4->DefineSection(5,fZ15,fBellowSectionOuterRadius-fBeryliumSectionThickness,fBellowSectionOuterRadius-fBeryliumSectionThickness+fBellowPlieThickness);
  TGeoVolume* votube4 = new TGeoVolume("votube4",tube4,kMedAlu2219);
  votube4->SetLineColor(kBlue);
  beamPipeCsideSection->AddNode(votube4,1,new TGeoTranslation(0., 0., 0.));
  
  TGeoPcon* tube4Vacuum = new TGeoPcon(0., 360., 2);
  tube4Vacuum->DefineSection(0,fZ10,0., fBellowSectionOuterRadius-fBeryliumSectionThickness);
  tube4Vacuum->DefineSection(1,fZ15,0.,fBellowSectionOuterRadius-fBeryliumSectionThickness);
  TGeoVolume* votube4Vacuum = new TGeoVolume("tube4Vacuum",tube4Vacuum,kMedVac);
  votube4Vacuum->SetVisibility(0);

  beamPipeCsideSection->AddNode(votube4Vacuum,1,new TGeoTranslation(0., 0., 0.));
  
  
  // ------------------ Second Bellow --------------------
  TGeoVolume* vobellows2 = MakeBellowCside("bellows2", fNBellowConvolutions, fBellowSectionOuterRadius-fBeryliumSectionThickness, fBellowOuterRadius, fBellowPlieRadius ,fBellowPlieThickness);
  beamPipeCsideSection->AddNode(vobellows2, 1, new TGeoTranslation(0., 0., fSecondBellowZmax-fBellowLength/2. - 2.*fBellowPlieRadius));
  // -----------------------------------------------------
 
  Float_t fZ16 = fSecondBellowZmax - fBellowLength; // End of Second bellow
  Float_t fZ18 = fZ16 - fThinPartBeforeBellowLength;
  Float_t fZ17 = fZ18 +  (fBeryliumSectionThickness-fBellowPlieThickness) / TMath::Tan(fChangeThicknessAngle); // End of 300 microns thickness part after first bellow
  Float_t fZ19 = fAdaptConeZmax; // Start of the Adpation Cone
  Float_t fZ20 = fAdaptConeZmin; // End of the Adpation Cone
  Float_t fZ21 = fAdaptConeZmin - fFlangeLength; // End of the Flange
  
  
  //----------- 15 deg Conical adaptator + flange ----------
  TGeoPcon* adaptator = new TGeoPcon(0., 360., 7);
  adaptator->DefineSection(0,fZ16, fBellowSectionOuterRadius-fBeryliumSectionThickness,fBellowSectionOuterRadius-fBeryliumSectionThickness+fBellowPlieThickness);
  adaptator->DefineSection(1,fZ17, fBellowSectionOuterRadius-fBeryliumSectionThickness,fBellowSectionOuterRadius-fBeryliumSectionThickness+fBellowPlieThickness);
  adaptator->DefineSection(2,fZ18, fBellowSectionOuterRadius-fBeryliumSectionThickness,fBellowSectionOuterRadius);
  adaptator->DefineSection(3,fZ19, fBellowSectionOuterRadius-fBeryliumSectionThickness,fBellowSectionOuterRadius);
  adaptator->DefineSection(4,fZ20, fAdaptConeRmax-fBeryliumSectionThickness,fAdaptConeRmax);
  adaptator->DefineSection(5,fZ20, fAdaptConeRmax-fBeryliumSectionThickness,fFlangeRmax);
  adaptator->DefineSection(6,fZ21, fAdaptConeRmax-fBeryliumSectionThickness,fFlangeRmax);
  TGeoVolume* voadaptator = new TGeoVolume("voadaptator",adaptator,kMedAlu2219);
  voadaptator->SetLineColor(kBlue);
  beamPipeCsideSection->AddNode(voadaptator,1,new TGeoTranslation(0., 0., 0.));

  TGeoPcon* adaptatorvide = new TGeoPcon(0., 360., 4);
  adaptatorvide->DefineSection(0,fZ16, 0., fBellowSectionOuterRadius-fBeryliumSectionThickness);
  adaptatorvide->DefineSection(1,fZ19, 0., fBellowSectionOuterRadius-fBeryliumSectionThickness);
  adaptatorvide->DefineSection(2,fZ20, 0., fAdaptConeRmax-fBeryliumSectionThickness);
  adaptatorvide->DefineSection(3,fZ21, 0., fAdaptConeRmax-fBeryliumSectionThickness);
  TGeoVolume* voadaptatorvide = new TGeoVolume("voadaptatorvide",adaptatorvide,kMedVac);
  voadaptatorvide->SetVisibility(0);
//  voadaptatorvide->SetLineColor(kGreen);
  beamPipeCsideSection->AddNode(voadaptatorvide,1,new TGeoTranslation(0., 0., 0.));
  //------------------------------------------------------

  top->AddNode(beamPipeCsideSection,1);
 
  ///////////////////////////////////
  //    Beam Pipe support          //
  ///////////////////////////////////

  
  // Beam Pipe Support
  TGeoVolume *beamPipeSupport = new TGeoVolumeAssembly("BeamPipeSupport");
  Float_t beamPipesupportZpos = fZ5;
  
  // Dimensions :
  
  Float_t supportXdim= 20.67;
  Float_t beamPipeRingZdim = 4.0;
  Float_t vespelRmax = 2.3;
  Float_t vespelRmin = 2.22;
  Float_t beampipeCarbonCollarRmin = 2.4;
  Float_t beampipeCarbonCollarRmax = 2.7;
  
  Float_t fixationCarbonCollarRmin = 1.5;
  Float_t fixationCarbonCollarRmax = 1.7;
  Float_t fixationCarbonCollarDZ = 2.5;
  
  
  Float_t skinThickness = 0.1;
  Float_t skinXdim = 14.25;
  Float_t skinYdim = 1.;
  Float_t skinZdim = fixationCarbonCollarDZ;
  Float_t carbonEarsXdim = 1.01;
  Float_t carbonEarsYdim = 0.2;
  Float_t carbonEarsZdim = fixationCarbonCollarDZ;
  
  // Support Bar
  TGeoVolumeAssembly *supportBar = new TGeoVolumeAssembly("BPS_SupportBar");
  
  TGeoBBox * carbonSkinBPS = new TGeoBBox(skinXdim/2.,skinYdim/2.,skinZdim/2.);
  carbonSkinBPS->SetName("carbonSkinBPS");
  
  TGeoBBox * foambarBPS = new TGeoBBox("foambarBPS",skinXdim/2.-skinThickness,skinYdim/2.-skinThickness,skinZdim/2.-skinThickness/2.);
  TGeoBBox * carbonEarsBPS = new TGeoBBox(carbonEarsXdim/2.,carbonEarsYdim/2.,carbonEarsZdim/2.);
  carbonEarsBPS->SetName("carbonEarsBPS");
  
  TGeoTranslation * transBP1 = new TGeoTranslation("transBP1",(skinXdim+carbonEarsXdim)/2.,0.,0.);
  transBP1->RegisterYourself();
  TGeoTranslation * transBP2 = new TGeoTranslation("transBP2",-(skinXdim+carbonEarsXdim)/2.,0.,0.);
  transBP2->RegisterYourself();
  TGeoCompositeShape *supportBarCarbon = new TGeoCompositeShape("BPS_supportBarCarbon", "(carbonSkinBPS-foambarBPS)+carbonEarsBPS:transBP1+carbonEarsBPS:transBP2");
  
  TGeoVolume *supportBarCarbonVol = new TGeoVolume("BPS_supportBarCarbon",supportBarCarbon,kMedCarbonFiber);
  supportBarCarbonVol->SetLineColor(kGray+3);
  
  supportBar->AddNode(supportBarCarbonVol, 1, new TGeoTranslation(skinXdim/2.+carbonEarsXdim+beampipeCarbonCollarRmax,0,0));
  supportBar->AddNode(supportBarCarbonVol, 2, new TGeoTranslation(-(skinXdim/2.+carbonEarsXdim+beampipeCarbonCollarRmax),0,0));
  
  TGeoVolume *foamVol = new TGeoVolume("supportBarFoam",foambarBPS,kMedRohacell);
  foamVol->SetLineColor(kGray);
  supportBar->AddNode(foamVol, 1, new TGeoTranslation(skinXdim/2.+carbonEarsXdim+beampipeCarbonCollarRmax,0,0));
  supportBar->AddNode(foamVol, 2, new TGeoTranslation(-(skinXdim/2.+carbonEarsXdim+beampipeCarbonCollarRmax),0,0));
  
  beamPipeSupport->AddNode(supportBar,1);
  
  
  // Fixation to wings
  
  TGeoVolumeAssembly *fixationToWings = new TGeoVolumeAssembly("BPS_fixationToWings");
  
  Float_t delatX = 0.1;
  
  TGeoTubeSeg * fixationTube = new TGeoTubeSeg(fixationCarbonCollarRmin,fixationCarbonCollarRmax,fixationCarbonCollarDZ/2.,-90.,90.);
  fixationTube->SetName("fixationTube");
  TGeoBBox * fixationToBar = new TGeoBBox(carbonEarsXdim/2.+delatX,carbonEarsYdim/2.,carbonEarsZdim/2.);
  fixationToBar->SetName("fixationToBar");
  
  TGeoTranslation * transBP3 = new TGeoTranslation("transBP3",fixationCarbonCollarRmax+carbonEarsXdim/2.-delatX,carbonEarsYdim,0.);
  transBP3->RegisterYourself();
  TGeoTranslation * transBP4 = new TGeoTranslation("transBP4",fixationCarbonCollarRmax+carbonEarsXdim/2.-delatX,-carbonEarsYdim,0.);
  transBP4->RegisterYourself();
  TGeoCompositeShape *fixationToWing = new TGeoCompositeShape("fixationToWing", "fixationTube+fixationToBar:transBP3+fixationToBar:transBP4");
  
  TGeoVolume *fixationToWingVol = new TGeoVolume("fixationToWing",fixationToWing,kMedCarbonFiber);
  fixationToWingVol->SetLineColor(kGray+2);
  
  
  fixationToWings->AddNode(fixationToWingVol,1, new TGeoTranslation(-supportXdim,0,0));
  fixationToWings->AddNode(fixationToWingVol,2, new TGeoCombiTrans(+supportXdim,0,0,new TGeoRotation("rot",0.,0.,180.)));
  
  
  beamPipeSupport->AddNode(fixationToWings,1);
  
  
  // Fixation to pipe
  
  TGeoVolumeAssembly *fixationToPipe = new TGeoVolumeAssembly("fixationToPipe");
  
  TGeoTubeSeg * pipeSupportTubeCarbon = new TGeoTubeSeg(beampipeCarbonCollarRmin,beampipeCarbonCollarRmax,fixationCarbonCollarDZ/2.,0.,180.);
  pipeSupportTubeCarbon->SetName("pipeSupportTubeCarbon");
  
  TGeoBBox * fixationTubeToBar = new TGeoBBox(carbonEarsXdim/2.+delatX,carbonEarsYdim/2.,carbonEarsZdim/2.);
  fixationTubeToBar->SetName("fixationTubeToBar");
  TGeoBBox * hole = new TGeoBBox((beampipeCarbonCollarRmax-vespelRmin)/2.,carbonEarsYdim/2.,carbonEarsZdim/2.+1e-3);
  hole->SetName("hole");
  
  TGeoTranslation * transBP5 = new TGeoTranslation("transBP5",beampipeCarbonCollarRmax+carbonEarsXdim/2.-delatX,carbonEarsYdim,0.);
  transBP5->RegisterYourself();
  TGeoTranslation * transBP6 = new TGeoTranslation("transBP6",-(beampipeCarbonCollarRmax+carbonEarsXdim/2.-delatX),carbonEarsYdim,0.);
  transBP6->RegisterYourself();
  TGeoTranslation * transBP7 = new TGeoTranslation("transBP7",(beampipeCarbonCollarRmax+vespelRmin)/2.,0.,0.);
  transBP7->RegisterYourself();
  TGeoTranslation * transBP8 = new TGeoTranslation("transBP8",-((beampipeCarbonCollarRmax+vespelRmin)/2.),0.,0.);
  transBP8->RegisterYourself();
  TGeoCompositeShape *halfFixationToPipe = new TGeoCompositeShape("halfFixationToPipe", "(pipeSupportTubeCarbon-hole:transBP7-hole:transBP8)+fixationTubeToBar:transBP5+fixationTubeToBar:transBP6");
  
  TGeoVolume *halfFixationToPipeVol = new TGeoVolume("halfFixationToPipe",halfFixationToPipe,kMedCarbonFiber);
  halfFixationToPipeVol->SetLineColor(kRed+2);
  
  fixationToPipe->AddNode(halfFixationToPipeVol,1);
  fixationToPipe->AddNode(halfFixationToPipeVol,2, new TGeoCombiTrans(0,0,0,new TGeoRotation("rot",0.,0.,180.)));
  
  beamPipeSupport->AddNode(fixationToPipe,1);
  
  
  // Beam Pipe Ring
  
  TGeoVolumeAssembly *beamPipeRing = new TGeoVolumeAssembly("beamPipeRing");
  
  TGeoTube * beamPipeRingCarbon = new TGeoTube(vespelRmax,beampipeCarbonCollarRmin,beamPipeRingZdim/2.);
  TGeoVolume *beamPipeRingCarbonVol = new TGeoVolume("beamPipeRingCarbon",beamPipeRingCarbon,kMedCarbonFiber);
  beamPipeRingCarbonVol->SetLineColor(kGreen+2);
  beamPipeRing->AddNode(beamPipeRingCarbonVol,1, new TGeoTranslation(0.,0,(beamPipeRingZdim-fixationCarbonCollarDZ)/2.));
  
  TGeoTube * beamPipeRingVespel = new TGeoTube(vespelRmin,vespelRmax,beamPipeRingZdim/2.);
  TGeoVolume *beamPipeRingVespelVol = new TGeoVolume("beamPipeRingVespel",beamPipeRingVespel,kMedPolyimide);
  beamPipeRingVespelVol->SetLineColor(kGreen+4);
  beamPipeRing->AddNode(beamPipeRingVespelVol,1, new TGeoTranslation(0.,0,(beamPipeRingZdim-fixationCarbonCollarDZ)/2.));
  
  beamPipeSupport->AddNode(beamPipeRing,1);
  beamPipeSupport->SetVisibility(0);
  
  top->AddNode(beamPipeSupport,1,new TGeoTranslation(0.,0,beamPipesupportZpos+fixationCarbonCollarDZ/2.));
  
  
  
  ///////////// END NEW BEAM PIPE GEOMETRY fOR MFT ////////////////////
  
  ////////////////////////////////////////////////////////////////////////////////
  // Side A section after Beryllium
  // !!!!!! THIS PART NEED TO BE WORKED OUT !!!!!
  ////////////////////////////////////////////////////////////////////////////////

  //----------------  Al tube ------------------
  Float_t fSmallRadiusZmax =fBeryliumSectionZmax + 20.43;
  Float_t fLargeRadiusZmin =fSmallRadiusZmax + 2.61;
  Float_t fAdaptConeSideAZmin =fLargeRadiusZmin + 200.; // THIS PART NEED TO BE WORKED OUT
  Float_t fAluSideARmax = 2.5;

  TGeoPcon* aluSideA = new TGeoPcon(0., 360., 4);
  aluSideA->DefineSection(0,fBeryliumSectionZmax, fBeryliumSectionOuterRadius-fBeryliumSectionThickness,fBeryliumSectionOuterRadius);
  aluSideA->DefineSection(1,fSmallRadiusZmax, fBeryliumSectionOuterRadius-fBeryliumSectionThickness,fBeryliumSectionOuterRadius);
  aluSideA->DefineSection(2,fLargeRadiusZmin, fAluSideARmax-fBeryliumSectionThickness,fAluSideARmax);
  aluSideA->DefineSection(3,fAdaptConeSideAZmin, fAluSideARmax-fBeryliumSectionThickness,fAluSideARmax);
  TGeoVolume* voaluSideA = new TGeoVolume("aluSideA",aluSideA,kMedAlu2219);
  voaluSideA->SetLineColor(kBlue);
  top->AddNode(voaluSideA,1,new TGeoTranslation(0., 0., 0.));
  
  TGeoPcon* aluSideAVac = new TGeoPcon(0., 360., 4);
  aluSideAVac->DefineSection(0,fBeryliumSectionZmax, 0., fBeryliumSectionOuterRadius-fBeryliumSectionThickness);
  aluSideAVac->DefineSection(1,fSmallRadiusZmax, 0., fBeryliumSectionOuterRadius-fBeryliumSectionThickness);
  aluSideAVac->DefineSection(2,fLargeRadiusZmin, 0., fAluSideARmax-fBeryliumSectionThickness);
  aluSideAVac->DefineSection(3,fAdaptConeSideAZmin, 0., fAluSideARmax-fBeryliumSectionThickness);
  TGeoVolume* voaluSideAVac = new TGeoVolume("aluSideAVac",aluSideAVac,kMedVac);
  voaluSideAVac->SetLineColor(kGreen);
  voaluSideAVac->SetVisibility(0);
  top->AddNode(voaluSideAVac,1,new TGeoTranslation(0., 0., 0.));
  //-------------------------------------------------

  
  ////////////////////////////////////////////////////////////////////////////////     
  //                                                                            //
  //                                  RB24/1                                    // 
  //                                                                            //
  ////////////////////////////////////////////////////////////////////////////////
  //
  //
  // Drawing LHCVC2U_0001
  // Copper Tube RB24/1      393.5 cm 
  // Warm module VMACA        18.0 cm
  // Annular Ion Pump         35.0 cm
  // Valve                     7.5 cm
  // Warm module VMABC        28.0 cm
  // ================================
  //                         462.0 cm
  //

    
  // Copper Tube RB24/1
  const Float_t  kRB24CuTubeL   = 393.5;
  const Float_t  kRB24CuTubeRi  = 8.0/2.;
  const Float_t  kRB24CuTubeRo  = 8.4/2.;
  const Float_t  kRB24CuTubeFRo = 7.6;
  const Float_t  kRB24CuTubeFL  = 1.86;

  TGeoVolume* voRB24CuTubeM = new TGeoVolume("voRB24CuTubeM", 
					     new TGeoTube(0., kRB24CuTubeRo, kRB24CuTubeL/2.), kMedVac);
  voRB24CuTubeM->SetVisibility(0);
  TGeoVolume* voRB24CuTube  = new TGeoVolume("voRB24CuTube", 
					     new TGeoTube(kRB24CuTubeRi, kRB24CuTubeRo, kRB24CuTubeL/2.), kMedCu);
  voRB24CuTubeM->AddNode(voRB24CuTube, 1, gGeoIdentity);
  // Air outside tube with higher transport cuts
  TGeoVolume* voRB24CuTubeA  = new TGeoVolume("voRB24CuTubeA", 
					      new TGeoTube(25., 100., kRB24CuTubeL/2.), kMedAirHigh);
  voRB24CuTubeA->SetVisibility(0);
  // Simplified DN 100 Flange
  TGeoVolume* voRB24CuTubeF = new TGeoVolume("voRB24CuTubeF", 
					     new TGeoTube(kRB24CuTubeRo, kRB24CuTubeFRo, kRB24CuTubeFL/2.), kMedSteel);

  // Warm Module Type VMACA
  // LHCVMACA_0002
  // 
  // Pos 1 Warm Bellows DN100       LHCVBU__0012
  // Pos 2 RF Contact   D80         LHCVSR__0005
  // Pos 3 Trans. Tube Flange       LHCVSR__0065
  // [Pos 4 Hex. Countersunk Screw   Bossard BN4719]
  // [Pos 5 Tension spring           LHCVSR__0011]
  //
  //
  //
  // Pos1    Warm Bellows DN100
  // Pos1.1  Bellows                  LHCVBU__0006
  //
  //
  // Connection Tubes    
  // Connection tube inner r
  const Float_t kRB24B1ConTubeRin        = 10.0/2.;
  // Connection tube outer r
  const Float_t kRB24B1ConTubeRou        = 10.3/2.;
  // Connection tube length
  const Float_t kRB24B1ConTubeL          =  2.5;
  // 
  const Float_t kRB24B1CompL             = 16.00;    // Length of the compensator
  const Float_t kRB24B1BellowRi          = 10.25/2.; // Bellow inner radius        
  const Float_t kRB24B1BellowRo          = 11.40/2.; // Bellow outer radius        
  const Int_t   kRB24B1NumberOfPlies     = 27;       // Number of plies            
  const Float_t kRB24B1BellowUndL        = 11.00;    // Length of undulated region 
  const Float_t kRB24B1PlieThickness     =  0.015;   // Plie thickness             

  const Float_t kRB24B1PlieRadius = 
    (kRB24B1BellowUndL + (2. *  kRB24B1NumberOfPlies - 2.) * kRB24B1PlieThickness) / (4. * kRB24B1NumberOfPlies);

  const Float_t kRB24B1ProtTubeThickness = 0.02;     // Thickness of the protection tube
  const Float_t kRB24B1ProtTubeLength    = 4.2;      // Length of the protection tube

  const Float_t kRB24B1RFlangeL          = 1.86;     // Length of the flanges
  const Float_t kRB24B1RFlangeLO         = 0.26;     // Flange overlap
  const Float_t kRB24B1RFlangeRO         = 11.18/2;  // Inner radius at Flange overlap    
  const Float_t kRB24B1RFlangeRou        = 15.20/2.; // Outer radius of flange
  const Float_t kRB24B1RFlangeRecess     = 0.98;     // Flange recess
  const Float_t kRB24B1L                 = kRB24B1CompL +  2. * (kRB24B1RFlangeL - kRB24B1RFlangeRecess);
    
  ///      
  //
  // Bellow mother volume
  TGeoPcon* shRB24B1BellowM = new TGeoPcon(0., 360., 14);
  // Connection Tube and Flange
  z = 0.;
  shRB24B1BellowM->DefineSection( 0, z, 0.,               kRB24B1RFlangeRou);
  z += kRB24B1RFlangeLO;
  shRB24B1BellowM->DefineSection( 1, z, 0.,               kRB24B1RFlangeRou);
  shRB24B1BellowM->DefineSection( 2, z, 0.,               kRB24B1RFlangeRou);    
  z = kRB24B1RFlangeL;
  shRB24B1BellowM->DefineSection( 3, z, 0.,               kRB24B1RFlangeRou);    
  shRB24B1BellowM->DefineSection( 4, z, 0.,               kRB24B1ConTubeRou);
  z = kRB24B1ConTubeL +  kRB24B1RFlangeL - kRB24B1RFlangeRecess;
  shRB24B1BellowM->DefineSection( 5, z, 0.,               kRB24B1ConTubeRou);
  // Plie
  shRB24B1BellowM->DefineSection( 6, z, 0.,               kRB24B1BellowRo + kRB24B1ProtTubeThickness);
  z += kRB24B1BellowUndL;
  shRB24B1BellowM->DefineSection( 7, z, 0.,               kRB24B1BellowRo + kRB24B1ProtTubeThickness);
  shRB24B1BellowM->DefineSection( 8, z, 0.,               kRB24B1ConTubeRou);
  // Connection Tube and Flange
  z = kRB24B1L - shRB24B1BellowM->GetZ(3);
  shRB24B1BellowM->DefineSection( 9, z, 0.,               kRB24B1ConTubeRou);
  shRB24B1BellowM->DefineSection(10, z, 0.,               kRB24B1RFlangeRou);
  z = kRB24B1L - shRB24B1BellowM->GetZ(1);
  shRB24B1BellowM->DefineSection(11, z, 0.,               kRB24B1RFlangeRou);
  shRB24B1BellowM->DefineSection(12, z, 0.,               kRB24B1RFlangeRou);
  z = kRB24B1L - shRB24B1BellowM->GetZ(0);
  shRB24B1BellowM->DefineSection(13, z, 0.,               kRB24B1RFlangeRou);

  TGeoVolume* voRB24B1BellowM = new TGeoVolume("RB24B1BellowM", shRB24B1BellowM, kMedVac);
  voRB24B1BellowM->SetVisibility(0);
  //
  // Bellow Section    
  TGeoVolume* voRB24B1Bellow 
    = MakeBellow("RB24B1", kRB24B1NumberOfPlies, kRB24B1BellowRi, kRB24B1BellowRo, 
		 kRB24B1BellowUndL, kRB24B1PlieRadius ,kRB24B1PlieThickness);
  voRB24B1Bellow->SetVisibility(0);
    
  //
  // End Parts (connection tube)
  TGeoVolume* voRB24B1CT = new TGeoVolume("RB24B1CT", new TGeoTube(kRB24B1ConTubeRin, kRB24B1ConTubeRou,  kRB24B1ConTubeL/2.), kMedSteel); 
  //
  // Protection Tube      
  TGeoVolume* voRB24B1PT = new TGeoVolume("RB24B1PT", new TGeoTube(kRB24B1BellowRo, kRB24B1BellowRo + kRB24B1ProtTubeThickness,  
								   kRB24B1ProtTubeLength / 2.), kMedSteel);
    
  z = kRB24B1ConTubeL/2. +  (kRB24B1RFlangeL - kRB24B1RFlangeRecess);
    
  voRB24B1BellowM->AddNode(voRB24B1CT, 1, new TGeoTranslation(0., 0., z));
  z += (kRB24B1ConTubeL/2.+ kRB24B1BellowUndL/2.);
  voRB24B1BellowM->AddNode(voRB24B1Bellow, 1, new TGeoTranslation(0., 0., z));
  z += (kRB24B1BellowUndL/2. + kRB24B1ConTubeL/2);
  voRB24B1BellowM->AddNode(voRB24B1CT, 2, new TGeoTranslation(0., 0., z));
  z =  kRB24B1ConTubeL +  kRB24B1ProtTubeLength / 2. + 1. + kRB24B1RFlangeLO;
  voRB24B1BellowM->AddNode(voRB24B1PT, 1, new TGeoTranslation(0., 0., z));
  z +=  kRB24B1ProtTubeLength + 0.6;
  voRB24B1BellowM->AddNode(voRB24B1PT, 2, new TGeoTranslation(0., 0., z));

                 

  // Pos 1/2 Rotatable Flange         LHCVBU__0013
  // Pos 1/3 Flange DN100/103         LHCVBU__0018
  // The two flanges can be represented by the same volume
  // Outer Radius (including the outer movable ring).
  // The inner ring has a diameter of 12.04 cm

  
  TGeoPcon* shRB24B1RFlange = new TGeoPcon(0., 360., 10);
  z = 0.;
  shRB24B1RFlange->DefineSection(0, z, 10.30/2., kRB24B1RFlangeRou);
  z += 0.55;  // 5.5 mm added for outer ring
  z += 0.43;
  shRB24B1RFlange->DefineSection(1, z, 10.30/2., kRB24B1RFlangeRou);
  shRB24B1RFlange->DefineSection(2, z, 10.06/2., kRB24B1RFlangeRou);    
  z += 0.15;
  shRB24B1RFlange->DefineSection(3, z, 10.06/2., kRB24B1RFlangeRou);    
  // In reality this part is rounded
  shRB24B1RFlange->DefineSection(4, z, 10.91/2., kRB24B1RFlangeRou);    
  z += 0.15;
  shRB24B1RFlange->DefineSection(5, z, 10.91/2., kRB24B1RFlangeRou);    
  shRB24B1RFlange->DefineSection(6, z, 10.06/2., kRB24B1RFlangeRou);    
  z += 0.32;
  shRB24B1RFlange->DefineSection(7, z, 10.06/2., kRB24B1RFlangeRou);    
  shRB24B1RFlange->DefineSection(8, z, kRB24B1RFlangeRO, kRB24B1RFlangeRou);    
  z += kRB24B1RFlangeLO;
  shRB24B1RFlange->DefineSection(9, z, kRB24B1RFlangeRO, kRB24B1RFlangeRou);    
    
  TGeoVolume* voRB24B1RFlange = new TGeoVolume("RB24B1RFlange", shRB24B1RFlange, kMedSteel);

    
  z = kRB24B1L - kRB24B1RFlangeL;
  voRB24B1BellowM->AddNode(voRB24B1RFlange, 1, new TGeoTranslation(0., 0., z));
  z = kRB24B1RFlangeL;
  voRB24B1BellowM->AddNode(voRB24B1RFlange, 2, new TGeoCombiTrans(0., 0., z, rot180));
  //
  // Pos 2 RF Contact   D80         LHCVSR__0005
  //
  // Pos 2.1 RF Contact Flange      LHCVSR__0003
  //
  TGeoPcon* shRB24B1RCTFlange = new TGeoPcon(0., 360., 6);
  const Float_t kRB24B1RCTFlangeRin  = 8.06/2. + 0.05;  // Inner radius
  const Float_t kRB24B1RCTFlangeL    = 1.45;            // Length
    
  z = 0.;
  shRB24B1RCTFlange->DefineSection(0, z, kRB24B1RCTFlangeRin,  8.20/2.);
  z += 0.15;
  shRB24B1RCTFlange->DefineSection(1, z, kRB24B1RCTFlangeRin,  8.20/2.);
  shRB24B1RCTFlange->DefineSection(2, z, kRB24B1RCTFlangeRin,  8.60/2.);
  z += 1.05;
  shRB24B1RCTFlange->DefineSection(3, z, kRB24B1RCTFlangeRin,  8.60/2.);
  shRB24B1RCTFlange->DefineSection(4, z, kRB24B1RCTFlangeRin, 11.16/2.);
  z += 0.25;
  shRB24B1RCTFlange->DefineSection(5, z, kRB24B1RCTFlangeRin, 11.16/2.);
  TGeoVolume* voRB24B1RCTFlange = new TGeoVolume("RB24B1RCTFlange", shRB24B1RCTFlange, kMedCu);
  z = kRB24B1L - kRB24B1RCTFlangeL;
    
  voRB24B1BellowM->AddNode(voRB24B1RCTFlange, 1, new TGeoTranslation(0., 0., z));
  //
  // Pos 2.2 RF-Contact        LHCVSR__0004
  //
  TGeoPcon* shRB24B1RCT = new TGeoPcon(0., 360., 3);
  const Float_t kRB24B1RCTRin  = 8.00/2.;        // Inner radius
  const Float_t kRB24B1RCTCRin = 8.99/2.;        // Max. inner radius conical section
  const Float_t kRB24B1RCTL    = 11.78;          // Length
  const Float_t kRB24B1RCTSL   = 10.48;          // Length of straight section
  const Float_t kRB24B1RCTd    =  0.03;          // Thickness
    
  z = 0;
  shRB24B1RCT->DefineSection(0, z,  kRB24B1RCTCRin,  kRB24B1RCTCRin + kRB24B1RCTd);
  z =  kRB24B1RCTL -  kRB24B1RCTSL;
  // In the (VSR0004) this section is straight in (LHCVC2U_0001) it is conical ????
  shRB24B1RCT->DefineSection(1, z,  kRB24B1RCTRin + 0.35,  kRB24B1RCTRin + 0.35 + kRB24B1RCTd);
  z = kRB24B1RCTL - 0.03;
  shRB24B1RCT->DefineSection(2, z,  kRB24B1RCTRin,  kRB24B1RCTRin + kRB24B1RCTd);

  TGeoVolume* voRB24B1RCT = new TGeoVolume("RB24B1RCT", shRB24B1RCT, kMedCu);
  z = kRB24B1L - kRB24B1RCTL - 0.45;
  voRB24B1BellowM->AddNode(voRB24B1RCT, 1, new TGeoTranslation(0., 0., z));    

  //
  // Pos 3 Trans. Tube Flange       LHCVSR__0065
  //
  // Pos 3.1 Transition Tube D53    LHCVSR__0064
  // Pos 3.2 Transition Flange      LHCVSR__0060
  // Pos 3.3 Transition Tube        LHCVSR__0058
  TGeoPcon* shRB24B1TTF = new TGeoPcon(0., 360., 7);
  // Flange
  z = 0.;
  shRB24B1TTF->DefineSection(0, z,  6.30/2., 11.16/2.);
  z += 0.25;
  shRB24B1TTF->DefineSection(1, z,  6.30/2., 11.16/2.);
  shRB24B1TTF->DefineSection(2, z,  6.30/2.,  9.3/2.);
  z += 0.55;
  shRB24B1TTF->DefineSection(3, z,  6.30/2.,  9.3/2.);
  // Tube
  shRB24B1TTF->DefineSection(4, z,  6.30/2.,  6.7/2.);
  z += 5.80;
  shRB24B1TTF->DefineSection(5, z,  6.30/2.,  6.7/2.);
  // Transition Tube
  z += 3.75;
  shRB24B1TTF->DefineSection(6, z,  8.05/2.,  8.45/2.);
  TGeoVolume* voRB24B1TTF = new TGeoVolume("RB24B1TTF", shRB24B1TTF, kMedSteel);
  z =  0.;
  voRB24B1BellowM->AddNode(voRB24B1TTF, 1, new TGeoTranslation(0., 0., z));    

  // Annular Ion Pump        
  // LHCVC2U_0003
  //
  // Pos  1 Rotable Flange         LHCVFX__0031
  // Pos  2 RF Screen Tube         LHCVC2U_0005
  // Pos  3 Shell                  LHCVC2U_0007
  // Pos  4 Extruded Shell         LHCVC2U_0006
  // Pos  5 Feedthrough Tube       LHCVC2U_0004
  // Pos  6 Tubulated Flange       STDVFUHV0021
  // Pos  7 Fixed Flange           LHCVFX__0032
  // Pos  8 Pumping Elements

  //
  // Pos 1 Rotable Flange          LHCVFX__0031
  // pos 7 Fixed Flange            LHCVFX__0032
  //
  //  Mother volume
  const Float_t kRB24AIpML = 35.;
    
  TGeoVolume* voRB24AIpM = new TGeoVolume("voRB24AIpM", new TGeoTube(0., 10., kRB24AIpML/2.), kMedAir);
  voRB24AIpM->SetVisibility(0);
    
  //
  // Length 35 cm
  // Flange 2 x 1.98 =   3.96
  // Tube            =  32.84
  //==========================
  //                    36.80
  // Overlap 2 * 0.90 =  1.80
                        
  const Float_t kRB24IpRFD1     =  0.68;    // Length of section 1
  const Float_t kRB24IpRFD2     =  0.30;    // Length of section 2						     
  const Float_t kRB24IpRFD3     =  0.10;    // Length of section 3						           
  const Float_t kRB24IpRFD4     =  0.35;    // Length of section 4						           
  const Float_t kRB24IpRFD5     =  0.55;    // Length of section 5						           
    
  const Float_t kRB24IpRFRo     = 15.20/2.; // Flange outer radius 
  const Float_t kRB24IpRFRi1    =  6.30/2.; // Flange inner radius section 1
  const Float_t kRB24IpRFRi2    =  6.00/2.; // Flange inner radius section 2
  const Float_t kRB24IpRFRi3    =  5.84/2.; // Flange inner radius section 3    
  const Float_t kRB24IpRFRi4    =  6.00/2.; // Flange inner radius section 1
  const Float_t kRB24IpRFRi5    = 10.50/2.; // Flange inner radius section 2

  TGeoPcon* shRB24IpRF = new TGeoPcon(0., 360., 9);
  z0 = 0.;
  shRB24IpRF->DefineSection(0, z0, kRB24IpRFRi1, kRB24IpRFRo);
  z0 += kRB24IpRFD1;
  shRB24IpRF->DefineSection(1, z0, kRB24IpRFRi2, kRB24IpRFRo);
  z0 += kRB24IpRFD2;
  shRB24IpRF->DefineSection(2, z0, kRB24IpRFRi2, kRB24IpRFRo);
  shRB24IpRF->DefineSection(3, z0, kRB24IpRFRi3, kRB24IpRFRo);
  z0 += kRB24IpRFD3;
  shRB24IpRF->DefineSection(4, z0, kRB24IpRFRi3, kRB24IpRFRo);
  shRB24IpRF->DefineSection(5, z0, kRB24IpRFRi4, kRB24IpRFRo);
  z0 += kRB24IpRFD4;
  shRB24IpRF->DefineSection(6, z0, kRB24IpRFRi4, kRB24IpRFRo);
  shRB24IpRF->DefineSection(7, z0, kRB24IpRFRi5, kRB24IpRFRo);
  z0 += kRB24IpRFD5;
  shRB24IpRF->DefineSection(8, z0, kRB24IpRFRi5, kRB24IpRFRo);

  TGeoVolume* voRB24IpRF = new TGeoVolume("RB24IpRF", shRB24IpRF, kMedSteel);
    
  //
  // Pos  2 RF Screen Tube         LHCVC2U_0005
  //

  //
  // Tube
  Float_t kRB24IpSTTL  = 32.84;            // Total length of the tube
  Float_t kRB24IpSTTRi =  5.80/2.;         // Inner Radius
  Float_t kRB24IpSTTRo =  6.00/2.;         // Outer Radius
  TGeoVolume* voRB24IpSTT = new TGeoVolume("RB24IpSTT", new TGeoTube(kRB24IpSTTRi, kRB24IpSTTRo, kRB24IpSTTL/2.), kMedSteel);
  // Screen
  Float_t kRB24IpSTCL  =  0.4;             // Lenth of the crochet detail
  // Length of the screen 
  Float_t kRB24IpSTSL  =  9.00 - 2. * kRB24IpSTCL; 
  // Rel. position of the screen 
  Float_t kRB24IpSTSZ  =  7.00 + kRB24IpSTCL; 
  TGeoVolume* voRB24IpSTS = new TGeoVolume("RB24IpSTS", new TGeoTube(kRB24IpSTTRi, kRB24IpSTTRo, kRB24IpSTSL/2.), kMedSteel);
  // Vacuum
  TGeoVolume* voRB24IpSTV = new TGeoVolume("RB24IpSTV", new TGeoTube(0., kRB24IpSTTRi, kRB24AIpML/2.), kMedVac);
  //
  voRB24IpSTT->AddNode(voRB24IpSTS, 1, new TGeoTranslation(0., 0., kRB24IpSTSZ -  kRB24IpSTTL/2. +  kRB24IpSTSL/2.));
    
  // Crochets
  // Inner radius
  Float_t kRB24IpSTCRi  = kRB24IpSTTRo + 0.25;
  // Outer radius
  Float_t kRB24IpSTCRo  = kRB24IpSTTRo + 0.35;
  // Length of 1stsection
  Float_t kRB24IpSTCL1  = 0.15;
  // Length of 2nd section
  Float_t kRB24IpSTCL2  = 0.15;
  // Length of 3rd section
  Float_t kRB24IpSTCL3  = 0.10;
  // Rel. position of 1st Crochet


  TGeoPcon* shRB24IpSTC = new TGeoPcon(0., 360., 5);
  z0 = 0;
  shRB24IpSTC->DefineSection(0, z0, kRB24IpSTCRi, kRB24IpSTCRo);
  z0 += kRB24IpSTCL1;
  shRB24IpSTC->DefineSection(1, z0, kRB24IpSTCRi, kRB24IpSTCRo);
  shRB24IpSTC->DefineSection(2, z0, kRB24IpSTTRo, kRB24IpSTCRo);
  z0 += kRB24IpSTCL2;
  shRB24IpSTC->DefineSection(3, z0, kRB24IpSTTRo, kRB24IpSTCRo);
  z0 += kRB24IpSTCL3;
  shRB24IpSTC->DefineSection(4, z0, kRB24IpSTTRo, kRB24IpSTTRo + 0.001);
  TGeoVolume* voRB24IpSTC = new TGeoVolume("RB24IpSTC", shRB24IpSTC, kMedSteel);

  // Pos  3 Shell                  LHCVC2U_0007
  // Pos  4 Extruded Shell         LHCVC2U_0006
  Float_t kRB24IpShellL     =  4.45;    // Length of the Shell
  Float_t kRB24IpShellD     =  0.10;    // Wall thickness of the shell
  Float_t kRB24IpShellCTRi  =  6.70/2.; // Inner radius of the connection tube
  Float_t kRB24IpShellCTL   =  1.56;    // Length of the connection tube
  Float_t kRB24IpShellCARi  = 17.80/2.; // Inner radius of the cavity
  Float_t kRB24IpShellCCRo  = 18.20/2.; // Inner radius at the centre

  TGeoPcon* shRB24IpShell = new TGeoPcon(0., 360., 7);
  z0 = 0;
  shRB24IpShell->DefineSection(0, z0, kRB24IpShellCTRi, kRB24IpShellCTRi + kRB24IpShellD);
  z0 +=  kRB24IpShellCTL;
  shRB24IpShell->DefineSection(1, z0, kRB24IpShellCTRi, kRB24IpShellCTRi + kRB24IpShellD);
  shRB24IpShell->DefineSection(2, z0, kRB24IpShellCTRi, kRB24IpShellCARi + kRB24IpShellD);
  z0 += kRB24IpShellD;
  shRB24IpShell->DefineSection(3, z0, kRB24IpShellCARi, kRB24IpShellCARi + kRB24IpShellD);
  z0 = kRB24IpShellL - kRB24IpShellD;
  shRB24IpShell->DefineSection(4, z0, kRB24IpShellCARi, kRB24IpShellCARi + kRB24IpShellD);
  shRB24IpShell->DefineSection(5, z0, kRB24IpShellCARi, kRB24IpShellCCRo);
  z0 = kRB24IpShellL;
  shRB24IpShell->DefineSection(6, z0, kRB24IpShellCARi, kRB24IpShellCCRo);
  TGeoVolume* voRB24IpShell = new TGeoVolume("RB24IpShell", shRB24IpShell, kMedSteel);
    
  TGeoPcon* shRB24IpShellM   = MakeMotherFromTemplate(shRB24IpShell, 0, 6, kRB24IpShellCTRi , 13);
    
    
  for (Int_t i = 0; i < 6; i++) {
    z = 2. * kRB24IpShellL  - shRB24IpShellM->GetZ(5-i);
    Float_t rmin = shRB24IpShellM->GetRmin(5-i);
    Float_t rmax = shRB24IpShellM->GetRmax(5-i);
    shRB24IpShellM->DefineSection(7+i, z, rmin, rmax);
  }
    
  TGeoVolume* voRB24IpShellM = new TGeoVolume("RB24IpShellM", shRB24IpShellM, kMedVac);
  voRB24IpShellM->SetVisibility(0);
  voRB24IpShellM->AddNode(voRB24IpShell, 1, gGeoIdentity);
  voRB24IpShellM->AddNode(voRB24IpShell, 2, new TGeoCombiTrans(0., 0., 2. * kRB24IpShellL, rot180));
  //
  // Pos  8 Pumping Elements
  //
  //  Anode array
  TGeoVolume* voRB24IpPE = new TGeoVolume("voRB24IpPE", new TGeoTube(0.9, 1., 2.54/2.), kMedSteel);
  Float_t kRB24IpPEAR = 5.5;
    
  for (Int_t i = 0; i < 15; i++) {
    Float_t phi = Float_t(i) * 24.;
    Float_t x   =  kRB24IpPEAR * TMath::Cos(kDegRad * phi);
    Float_t y   =  kRB24IpPEAR * TMath::Sin(kDegRad * phi);
    voRB24IpShellM->AddNode(voRB24IpPE, i+1, new TGeoTranslation(x, y, kRB24IpShellL));
  }
    
    
  //
  //  Cathodes
  //
  // Here we could add some Ti strips

  // Postioning of elements
  voRB24AIpM->AddNode(voRB24IpRF,     1, new TGeoTranslation(0., 0., -kRB24AIpML/2.));
  voRB24AIpM->AddNode(voRB24IpRF,     2, new TGeoCombiTrans (0., 0., +kRB24AIpML/2., rot180));
  voRB24AIpM->AddNode(voRB24IpSTT,    1, new TGeoTranslation(0., 0., 0.));
  voRB24AIpM->AddNode(voRB24IpSTV,    1, new TGeoTranslation(0., 0., 0.));
  voRB24AIpM->AddNode(voRB24IpShellM, 1, new TGeoTranslation(0., 0., -kRB24AIpML/2. +  8.13));
  voRB24AIpM->AddNode(voRB24IpSTC,    1, new TGeoTranslation(0., 0., 8.13 - kRB24AIpML/2.));
  voRB24AIpM->AddNode(voRB24IpSTC,    2, new TGeoCombiTrans (0., 0., 8.14 + 8.9 - kRB24AIpML/2., rot180));
    
  //
  // Valve
  // VAC Series 47 DN 63 with manual actuator
  //
  const Float_t kRB24ValveWz = 7.5;
  const Float_t kRB24ValveDN = 10.0/2.;
  //
  //  Body containing the valve plate
  //
  const Float_t kRB24ValveBoWx =  15.6;
  const Float_t kRB24ValveBoWy = (21.5 + 23.1 - 5.);
  const Float_t kRB24ValveBoWz =  4.6;
  const Float_t kRB24ValveBoD  =  0.5;

  TGeoVolume* voRB24ValveBoM =
    new TGeoVolume("RB24ValveBoM", 
		   new TGeoBBox( kRB24ValveBoWx/2.,  kRB24ValveBoWy/2., kRB24ValveBoWz/2.), kMedAir);
  voRB24ValveBoM->SetVisibility(0);
  TGeoVolume* voRB24ValveBo =
    new TGeoVolume("RB24ValveBo", 
		   new TGeoBBox( kRB24ValveBoWx/2.,  kRB24ValveBoWy/2., kRB24ValveBoWz/2.), kMedSteel);
  voRB24ValveBoM->AddNode(voRB24ValveBo, 1, gGeoIdentity);
  //
  // Inner volume
  //
  TGeoVolume* voRB24ValveBoI = new TGeoVolume("RB24ValveBoI", 
					      new TGeoBBox( kRB24ValveBoWx/2. -  kRB24ValveBoD,  
							    kRB24ValveBoWy/2. -  kRB24ValveBoD/2., 
							    kRB24ValveBoWz/2. -  kRB24ValveBoD), 
					      kMedVac);
  voRB24ValveBo->AddNode(voRB24ValveBoI, 1, new TGeoTranslation(0., kRB24ValveBoD/2., 0.));
  //
  // Opening and Flanges
  const Float_t  kRB24ValveFlRo = 18./2.;
  const Float_t  kRB24ValveFlD  = 1.45;    
  TGeoVolume* voRB24ValveBoA = new TGeoVolume("RB24ValveBoA", 
					      new TGeoTube(0., kRB24ValveDN/2., kRB24ValveBoD/2.), kMedVac);
  voRB24ValveBo->AddNode(voRB24ValveBoA, 1, new TGeoTranslation(0., - kRB24ValveBoWy/2. + 21.5, -kRB24ValveBoWz/2. +  kRB24ValveBoD/2.));
  voRB24ValveBo->AddNode(voRB24ValveBoA, 2, new TGeoTranslation(0., - kRB24ValveBoWy/2. + 21.5, +kRB24ValveBoWz/2. -  kRB24ValveBoD/2.));
 
  TGeoVolume* voRB24ValveFl  = new TGeoVolume("RB24ValveFl",  new TGeoTube(kRB24ValveDN/2.,  kRB24ValveFlRo, kRB24ValveFlD/2.), kMedSteel);
  TGeoVolume* voRB24ValveFlI = new TGeoVolume("RB24ValveFlI", new TGeoTube(0.,               kRB24ValveFlRo, kRB24ValveFlD/2.), kMedVac);
  voRB24ValveFlI->AddNode(voRB24ValveFl, 1, gGeoIdentity);
    
  //
  // Actuator Flange
  const Float_t kRB24ValveAFlWx =  18.9;
  const Float_t kRB24ValveAFlWy =   5.0;
  const Float_t kRB24ValveAFlWz =   7.7;
  TGeoVolume* voRB24ValveAFl = new TGeoVolume("RB24ValveAFl", new TGeoBBox(kRB24ValveAFlWx/2., kRB24ValveAFlWy/2., kRB24ValveAFlWz/2.), kMedSteel);
  //
  // Actuator Tube
  const Float_t kRB24ValveATRo = 9.7/2.;
  const Float_t kRB24ValveATH  = 16.6;
  TGeoVolume* voRB24ValveAT = new TGeoVolume("RB24ValveAT", new TGeoTube(kRB24ValveATRo -  2. * kRB24ValveBoD,kRB24ValveATRo,  kRB24ValveATH/2.), 
					     kMedSteel);
  //
  // Manual Actuator (my best guess)
  TGeoVolume* voRB24ValveMA1 = new TGeoVolume("RB24ValveMA1", new TGeoCone(2.5/2., 0., 0.5, 4.5, 5.), kMedSteel);
  TGeoVolume* voRB24ValveMA2 = new TGeoVolume("RB24ValveMA2", new TGeoTorus(5., 0., 1.25), kMedSteel);
  TGeoVolume* voRB24ValveMA3 = new TGeoVolume("RB24ValveMA3", new TGeoTube (0., 1.25, 2.5), kMedSteel);
    

  //
  // Position all volumes
  Float_t y0;
  TGeoVolumeAssembly*  voRB24ValveMo = new TGeoVolumeAssembly("RB24ValveMo");
  voRB24ValveMo->AddNode(voRB24ValveFl,  1, new TGeoTranslation(0., 0., - 7.5/2. + kRB24ValveFlD/2.));
  voRB24ValveMo->AddNode(voRB24ValveFl,  2, new TGeoTranslation(0., 0., + 7.5/2. - kRB24ValveFlD/2.));
  y0 = -21.5;
  voRB24ValveMo->AddNode(voRB24ValveBoM, 1, new TGeoTranslation(0., y0 + kRB24ValveBoWy/2.,   0.));
  y0 +=  kRB24ValveBoWy;
  voRB24ValveMo->AddNode(voRB24ValveAFl, 1, new TGeoTranslation(0., y0 +  kRB24ValveAFlWy/2., 0.));
  y0 +=  kRB24ValveAFlWy;
  voRB24ValveMo->AddNode(voRB24ValveAT,  1, new TGeoCombiTrans(0.,  y0 + kRB24ValveATH/2.,    0., rotyz));
  y0 += kRB24ValveATH;
  voRB24ValveMo->AddNode(voRB24ValveMA1, 1, new TGeoCombiTrans(0.,  y0 + 2.5/2.,    0., rotyz));
  y0 += 2.5;
  voRB24ValveMo->AddNode(voRB24ValveMA2, 1, new TGeoCombiTrans(0.,  y0 + 2.5/2.,    0., rotyz));
  y0 += 2.5;
  voRB24ValveMo->AddNode(voRB24ValveMA3, 1, new TGeoCombiTrans(5./TMath::Sqrt(2.),  y0 + 5.0/2., 5./TMath::Sqrt(2.), rotyz));
  //
  // Warm Module Type VMABC
  // LHCVMABC_0002
  // 
  //
  //
  // Flange                  1.00
  // Central Piece          11.50
  // Bellow                 14.50
  // End Flange              1.00
  //===================================
  // Total                  28.00 
  //                        
  // Pos 1 Warm Bellows DN100       LHCVBU__0016
  // Pos 2 Trans. Tube Flange       LHCVSR__0062
  // Pos 3 RF Contact   D63         LHCVSR__0057
  // [Pos 4 Hex. Countersunk Screw   Bossard BN4719]
  // [Pos 5 Tension spring           LHCVSR__00239]
  //

  // Pos 1 Warm Bellows DN100                   LHCVBU__0016
  // Pos 1.1 Right Body 2 Ports with Support    LHCVBU__0014
  //
  // Tube 1
  const Float_t kRB24VMABCRBT1Ri = 10.0/2.;
  const Float_t kRB24VMABCRBT1Ro = 10.3/2.;
  const Float_t kRB24VMABCRBT1L  = 11.5;   
  const Float_t kRB24VMABCRBT1L2 = 8.;
  const Float_t kRB24VMABCL      = 28.;
    
  TGeoTube* shRB24VMABCRBT1 = new TGeoTube(kRB24VMABCRBT1Ri, kRB24VMABCRBT1Ro, kRB24VMABCRBT1L/2.);
  shRB24VMABCRBT1->SetName("RB24VMABCRBT1");
  TGeoTube* shRB24VMABCRBT1o = new TGeoTube(0., kRB24VMABCRBT1Ro,  kRB24VMABCRBT1L/2.);
  shRB24VMABCRBT1o->SetName("RB24VMABCRBT1o");
  TGeoTube* shRB24VMABCRBT1o2 = new TGeoTube(0., kRB24VMABCRBT1Ro + 0.3, kRB24VMABCRBT1L/2.);
  shRB24VMABCRBT1o2->SetName("RB24VMABCRBT1o2");
  // Lower inforcement 
  TGeoVolume*  voRB24VMABCRBT12  = new TGeoVolume("RB24VMABCRBT12", 
						  new TGeoTubeSeg(kRB24VMABCRBT1Ro, kRB24VMABCRBT1Ro + 0.3, kRB24VMABCRBT1L2/2., 220., 320.)
						  , kMedSteel);
  //
  // Tube 2
  const Float_t kRB24VMABCRBT2Ri =   6.0/2.;
  const Float_t kRB24VMABCRBT2Ro =   6.3/2.;
  const Float_t kRB24VMABCRBF2Ro =  11.4/2.;
  const Float_t kRB24VMABCRBT2L  =   5.95 + 2.; // 2. cm added for welding    
  const Float_t kRB24VMABCRBF2L  =   1.75;
  TGeoTube* shRB24VMABCRBT2 = new TGeoTube(kRB24VMABCRBT2Ri, kRB24VMABCRBT2Ro,  kRB24VMABCRBT2L/2.);
  shRB24VMABCRBT2->SetName("RB24VMABCRBT2");
  TGeoTube* shRB24VMABCRBT2i = new TGeoTube(0., kRB24VMABCRBT2Ri, kRB24VMABCRBT2L/2. + 2.);
  shRB24VMABCRBT2i->SetName("RB24VMABCRBT2i");
  TGeoCombiTrans* tRBT2 = new TGeoCombiTrans(-11.5 + kRB24VMABCRBT2L/2., 0., 7.2 - kRB24VMABCRBT1L/2.  , rotxz);
  tRBT2->SetName("tRBT2");
  tRBT2->RegisterYourself();
  TGeoCompositeShape* shRB24VMABCRBT2c =  new TGeoCompositeShape("shRB24VMABCRBT2c","RB24VMABCRBT2:tRBT2-RB24VMABCRBT1o");
  TGeoVolume* voRB24VMABCRBT2 = new TGeoVolume("shRB24VMABCRBT2", shRB24VMABCRBT2c, kMedSteel);
  // Flange
  // Pos 1.4 Flange DN63                        LHCVBU__0008
  TGeoVolume* voRB24VMABCRBF2 = new TGeoVolume("RB24VMABCRBF2", 
					       new TGeoTube(kRB24VMABCRBT2Ro, kRB24VMABCRBF2Ro, kRB24VMABCRBF2L/2.), kMedSteel);
  // DN63 Blank Flange (my best guess)
  TGeoVolume* voRB24VMABCRBF2B = new TGeoVolume("RB24VMABCRBF2B", 
						new TGeoTube(0., kRB24VMABCRBF2Ro, kRB24VMABCRBF2L/2.), kMedSteel);
  //
  // Tube 3
  const Float_t kRB24VMABCRBT3Ri =  3.5/2.;
  const Float_t kRB24VMABCRBT3Ro =  3.8/2.;
  const Float_t kRB24VMABCRBF3Ro =  7.0/2.;
  const Float_t kRB24VMABCRBT3L  =  4.95 + 2.; // 2. cm added for welding    
  const Float_t kRB24VMABCRBF3L  =  1.27;
  TGeoTube* shRB24VMABCRBT3 = new TGeoTube(kRB24VMABCRBT3Ri, kRB24VMABCRBT3Ro,  kRB24VMABCRBT3L/2);
  shRB24VMABCRBT3->SetName("RB24VMABCRBT3");
  TGeoTube* shRB24VMABCRBT3i = new TGeoTube(0., kRB24VMABCRBT3Ri, kRB24VMABCRBT3L/2. + 2.);
  shRB24VMABCRBT3i->SetName("RB24VMABCRBT3i");
  TGeoCombiTrans* tRBT3 = new TGeoCombiTrans(0., 10.5 - kRB24VMABCRBT3L/2., 7.2 - kRB24VMABCRBT1L/2.  , rotyz);
  tRBT3->SetName("tRBT3");
  tRBT3->RegisterYourself();
  TGeoCompositeShape* shRB24VMABCRBT3c =  new TGeoCompositeShape("shRB24VMABCRBT3c","RB24VMABCRBT3:tRBT3-RB24VMABCRBT1o");
  TGeoVolume* voRB24VMABCRBT3 = new TGeoVolume("shRB24VMABCRBT3", shRB24VMABCRBT3c, kMedSteel);
  // Flange
  // Pos 1.4 Flange DN35                        LHCVBU__0007
  TGeoVolume* voRB24VMABCRBF3 = new TGeoVolume("RB24VMABCRBF3", 
					       new TGeoTube(kRB24VMABCRBT3Ro, kRB24VMABCRBF3Ro, kRB24VMABCRBF3L/2.), kMedSteel);
  //
  // Tube 4
  const Float_t kRB24VMABCRBT4Ri =  6.0/2.;
  const Float_t kRB24VMABCRBT4Ro =  6.4/2.;
  const Float_t kRB24VMABCRBT4L  =  6.6;    
  TGeoTube* shRB24VMABCRBT4 = new TGeoTube(kRB24VMABCRBT4Ri, kRB24VMABCRBT4Ro,  kRB24VMABCRBT4L/2.);
  shRB24VMABCRBT4->SetName("RB24VMABCRBT4");
  TGeoCombiTrans* tRBT4 = new TGeoCombiTrans(0.,-11.+kRB24VMABCRBT4L/2., 7.2 - kRB24VMABCRBT1L/2.  , rotyz);
  tRBT4->SetName("tRBT4");
  tRBT4->RegisterYourself();
  TGeoCompositeShape* shRB24VMABCRBT4c =  new TGeoCompositeShape("shRB24VMABCRBT4c","RB24VMABCRBT4:tRBT4-RB24VMABCRBT1o2");
  TGeoVolume* voRB24VMABCRBT4 = new TGeoVolume("shRB24VMABCRBT4", shRB24VMABCRBT4c, kMedSteel);
  TGeoCompositeShape* shRB24VMABCRB = new TGeoCompositeShape("shRB24VMABCRB", "RB24VMABCRBT1-(RB24VMABCRBT2i:tRBT2+RB24VMABCRBT3i:tRBT3)");
  TGeoVolume* voRB24VMABCRBI = new TGeoVolume("RB24VMABCRBI", shRB24VMABCRB, kMedSteel);
  //
  // Plate
  const Float_t kRB24VMABCRBBx = 16.0;
  const Float_t kRB24VMABCRBBy =  1.5;
  const Float_t kRB24VMABCRBBz = 15.0;
    
  // Relative position of tubes
  const Float_t  kRB24VMABCTz =   7.2;
  // Relative position of plate
  const Float_t  kRB24VMABCPz =   3.6;
  const Float_t  kRB24VMABCPy = -12.5;
    
  TGeoVolume* voRB24VMABCRBP = new TGeoVolume("RB24VMABCRBP", new TGeoBBox(kRB24VMABCRBBx/2., kRB24VMABCRBBy/2., kRB24VMABCRBBz/2.), kMedSteel);
  //
  // Pirani Gauge (my best guess)
  //
  TGeoPcon* shRB24VMABCPirani = new TGeoPcon(0., 360., 15);
  // DN35/16 Coupling
  z = 0;
  shRB24VMABCPirani->DefineSection( 0, z,  0.8 , kRB24VMABCRBF3Ro);
  z += kRB24VMABCRBF3L; // 1.3
  shRB24VMABCPirani->DefineSection( 1, z,  0.8 , kRB24VMABCRBF3Ro);
  shRB24VMABCPirani->DefineSection( 2, z,  0.8 , 1.0);
  // Pipe
  z += 2.8;
  shRB24VMABCPirani->DefineSection( 3, z,  0.8 , 1.0);
  // Flange
  shRB24VMABCPirani->DefineSection( 4, z,  0.8 , 1.75);
  z += 1.6;
  shRB24VMABCPirani->DefineSection( 5, z,  0.8 , 1.75);
  shRB24VMABCPirani->DefineSection( 6, z,  0.8 , 1.0);
  z += 5.2;
  shRB24VMABCPirani->DefineSection( 7, z,  0.8 , 1.0);
  shRB24VMABCPirani->DefineSection( 8, z,  0.8 , 2.5);    
  z += 2.0;
  shRB24VMABCPirani->DefineSection( 9, z,  0.80, 2.50);    
  shRB24VMABCPirani->DefineSection(10, z,  1.55, 1.75);    
  z += 5.7;
  shRB24VMABCPirani->DefineSection(11, z,  1.55, 1.75);    
  shRB24VMABCPirani->DefineSection(11, z,  0.00, 1.75);    
  z += 0.2;
  shRB24VMABCPirani->DefineSection(12, z,  0.00, 1.75);    
  shRB24VMABCPirani->DefineSection(13, z,  0.00, 0.75);    
  z += 0.5;
  shRB24VMABCPirani->DefineSection(14, z,  0.00, 0.75);  
  TGeoVolume* voRB24VMABCPirani = new TGeoVolume("RB24VMABCPirani", shRB24VMABCPirani, kMedSteel);
  //
  //
  // 
    
    
  //
  // Positioning of elements
  TGeoVolumeAssembly* voRB24VMABCRB = new TGeoVolumeAssembly("RB24VMABCRB");
  //
  voRB24VMABCRB->AddNode(voRB24VMABCRBI,   1, gGeoIdentity);
  // Plate
  voRB24VMABCRB->AddNode(voRB24VMABCRBP,   1, new TGeoTranslation(0., kRB24VMABCPy +  kRB24VMABCRBBy /2., 
								  kRB24VMABCRBBz/2. - kRB24VMABCRBT1L/2. +  kRB24VMABCPz));
  // Tube 2
  voRB24VMABCRB->AddNode(voRB24VMABCRBT2,  1, gGeoIdentity);
  // Flange Tube 2
  voRB24VMABCRB->AddNode(voRB24VMABCRBF2,  1, new TGeoCombiTrans(kRB24VMABCPy + kRB24VMABCRBF2L/2., 0.,  kRB24VMABCTz - kRB24VMABCRBT1L/2., rotxz));
  // Blank Flange Tube 2
  voRB24VMABCRB->AddNode(voRB24VMABCRBF2B, 1, new TGeoCombiTrans(kRB24VMABCPy- kRB24VMABCRBF2L/2., 0.,  kRB24VMABCTz - kRB24VMABCRBT1L/2., rotxz));    
  // Tube 3
  voRB24VMABCRB->AddNode(voRB24VMABCRBT3,  1, gGeoIdentity);
  // Flange Tube 3
  voRB24VMABCRB->AddNode(voRB24VMABCRBF3,  1, new TGeoCombiTrans(0.,   11.2 - kRB24VMABCRBF3L/2.,  kRB24VMABCTz - kRB24VMABCRBT1L/2., rotyz));
  // Pirani Gauge
  voRB24VMABCRB->AddNode(voRB24VMABCPirani, 1, new  TGeoCombiTrans(0., 11.2,  kRB24VMABCTz - kRB24VMABCRBT1L/2., rotyz));
  // Tube 4
  voRB24VMABCRB->AddNode(voRB24VMABCRBT4,  1, gGeoIdentity);
  // Inforcement 
  voRB24VMABCRB->AddNode(voRB24VMABCRBT12, 1, new TGeoTranslation(0., 0., kRB24VMABCRBT1L2/2. - kRB24VMABCRBT1L/2. + 2.8));
    

  // Pos 1.3 Bellows with end part              LHCVBU__0002
  //
  // Connection Tube    
  // Connection tube inner r
  const Float_t kRB24VMABBEConTubeRin        = 10.0/2.;
  // Connection tube outer r
  const Float_t kRB24VMABBEConTubeRou        = 10.3/2.;
  // Connection tube length
  const Float_t kRB24VMABBEConTubeL1         =  0.9;
  const Float_t kRB24VMABBEConTubeL2         =  2.6;
  //  const Float_t RB24VMABBEBellowL            =  kRB24VMABBEConTubeL1 + kRB24VMABBEConTubeL2 + kRB24B1BellowUndL;
    
  // Mother volume
  TGeoPcon* shRB24VMABBEBellowM = new TGeoPcon(0., 360., 6);
  // Connection Tube and Flange
  z = 0.;
  shRB24VMABBEBellowM->DefineSection( 0, z, kRB24VMABBEConTubeRin,  kRB24VMABBEConTubeRou);
  z += kRB24VMABBEConTubeL1;
  shRB24VMABBEBellowM->DefineSection( 1, z, kRB24VMABBEConTubeRin, kRB24VMABBEConTubeRou);
  shRB24VMABBEBellowM->DefineSection( 2, z, kRB24B1BellowRi,       kRB24B1BellowRo + kRB24B1ProtTubeThickness);
  z += kRB24B1BellowUndL;
  shRB24VMABBEBellowM->DefineSection( 3, z, kRB24B1BellowRi,       kRB24B1BellowRo + kRB24B1ProtTubeThickness);
  shRB24VMABBEBellowM->DefineSection( 4, z, kRB24VMABBEConTubeRin,  kRB24VMABBEConTubeRou);
  z += kRB24VMABBEConTubeL2;
  shRB24VMABBEBellowM->DefineSection( 5, z, kRB24VMABBEConTubeRin,  kRB24VMABBEConTubeRou);
  TGeoVolume* voRB24VMABBEBellowM = new TGeoVolume("RB24VMABBEBellowM", shRB24VMABBEBellowM, kMedVac);
  voRB24VMABBEBellowM->SetVisibility(0);
    
  //  Connection tube left
  TGeoVolume* voRB24VMABBECT1 = new TGeoVolume("RB24VMABBECT1", 
					       new TGeoTube(kRB24VMABBEConTubeRin, kRB24VMABBEConTubeRou,kRB24VMABBEConTubeL1/2.),
					       kMedSteel);
  //  Connection tube right
  TGeoVolume* voRB24VMABBECT2 = new TGeoVolume("RB24VMABBECT2", 
					       new TGeoTube(kRB24VMABBEConTubeRin, kRB24VMABBEConTubeRou,kRB24VMABBEConTubeL2/2.),
					       kMedSteel);
  z = kRB24VMABBEConTubeL1/2.;
  voRB24VMABBEBellowM->AddNode(voRB24VMABBECT1, 1, new TGeoTranslation(0., 0., z));
  z += kRB24VMABBEConTubeL1/2.;
  z += kRB24B1BellowUndL/2.;
  voRB24VMABBEBellowM->AddNode(voRB24B1Bellow, 2, new TGeoTranslation(0., 0., z));
  z += kRB24B1BellowUndL/2.;
  z += kRB24VMABBEConTubeL2/2.;
  voRB24VMABBEBellowM->AddNode(voRB24VMABBECT2, 1, new TGeoTranslation(0., 0., z));
  z += kRB24VMABBEConTubeL2/2.;

  voRB24VMABCRB->AddNode(voRB24VMABBEBellowM, 1, new TGeoTranslation(0., 0., kRB24VMABCRBT1L/2.));

  // Pos 1.2 Rotable flange                     LHCVBU__0013[*]
  // Front
  voRB24VMABCRB->AddNode(voRB24B1RFlange,  3, new TGeoCombiTrans(0., 0., - kRB24VMABCRBT1L/2. + 0.86, rot180));
  // End
  z =  kRB24VMABCRBT1L/2. + kRB24B1BellowUndL +kRB24VMABBEConTubeL1 +  kRB24VMABBEConTubeL2;
  voRB24VMABCRB->AddNode(voRB24B1RFlange,  4, new TGeoTranslation(0., 0., z - 0.86));


  // Pos 2    Trans. Tube Flange       LHCVSR__0062
  // Pos 2.1  Transition Tube          LHCVSR__0063
  // Pos 2.2  Transition Flange        LHCVSR__0060
  //
  // Transition Tube with Flange
  TGeoPcon* shRB24VMABCTT = new TGeoPcon(0., 360., 7);
  z = 0.;
  shRB24VMABCTT->DefineSection(0, z, 6.3/2., 11.16/2.);
  z += 0.25;
  shRB24VMABCTT->DefineSection(1, z, 6.3/2., 11.16/2.);
  shRB24VMABCTT->DefineSection(2, z, 6.3/2.,  9.30/2.);
  z += 0.25;
  shRB24VMABCTT->DefineSection(3, z, 6.3/2.,  9.30/2.);
  shRB24VMABCTT->DefineSection(4, z, 6.3/2.,  6.70/2.);
  z += (20.35 - 0.63);
  shRB24VMABCTT->DefineSection(5, z, 6.3/2.,  6.7/2.);
  z += 0.63;
  shRB24VMABCTT->DefineSection(6, z, 6.3/2.,  6.7/2.);
  TGeoVolume* voRB24VMABCTT = new TGeoVolume("RB24VMABCTT", shRB24VMABCTT, kMedSteel);
  voRB24VMABCRB->AddNode(voRB24VMABCTT, 1, new TGeoTranslation(0., 0., - kRB24VMABCRBT1L/2.-1.));

  // Pos 3   RF Contact   D63         LHCVSR__0057
  // Pos 3.1 RF Contact Flange        LHCVSR__0017
  //
  TGeoPcon* shRB24VMABCCTFlange = new TGeoPcon(0., 360., 6);
  const Float_t kRB24VMABCCTFlangeRin  = 6.36/2.;  // Inner radius
  const Float_t kRB24VMABCCTFlangeL    = 1.30;     // Length
    
  z = 0.;
  shRB24VMABCCTFlange->DefineSection(0, z, kRB24VMABCCTFlangeRin,  6.5/2.);
  z += 0.15;
  shRB24VMABCCTFlange->DefineSection(1, z, kRB24VMABCCTFlangeRin,  6.5/2.);
  shRB24VMABCCTFlange->DefineSection(2, z, kRB24VMABCCTFlangeRin,  6.9/2.);
  z += 0.9;
  shRB24VMABCCTFlange->DefineSection(3, z, kRB24VMABCCTFlangeRin,  6.9/2.);
  shRB24VMABCCTFlange->DefineSection(4, z, kRB24VMABCCTFlangeRin, 11.16/2.);
  z += 0.25;
  shRB24VMABCCTFlange->DefineSection(5, z, kRB24VMABCCTFlangeRin, 11.16/2.);
  TGeoVolume* voRB24VMABCCTFlange = new TGeoVolume("RB24VMABCCTFlange", shRB24VMABCCTFlange, kMedCu);
  //
  // Pos 3.2 RF-Contact        LHCVSR__0056
  //
  TGeoPcon* shRB24VMABCCT = new TGeoPcon(0., 360., 4);
  const Float_t kRB24VMABCCTRin  = 6.30/2.;        // Inner radius
  const Float_t kRB24VMABCCTCRin = 7.29/2.;        // Max. inner radius conical section
  const Float_t kRB24VMABCCTL    = 11.88;          // Length
  const Float_t kRB24VMABCCTSL   = 10.48;          // Length of straight section
  const Float_t kRB24VMABCCTd    =  0.03;          // Thickness
  z = 0;
  shRB24VMABCCT->DefineSection(0, z,  kRB24VMABCCTCRin,  kRB24VMABCCTCRin + kRB24VMABCCTd);
  z =  kRB24VMABCCTL -  kRB24VMABCCTSL;
  shRB24VMABCCT->DefineSection(1, z,  kRB24VMABCCTRin + 0.35,  kRB24VMABCCTRin + 0.35 + kRB24VMABCCTd);
  z = kRB24VMABCCTL  -  kRB24VMABCCTFlangeL;
  shRB24VMABCCT->DefineSection(2, z,  kRB24VMABCCTRin,  kRB24VMABCCTRin + kRB24VMABCCTd);
  z = kRB24VMABCCTL;
  shRB24VMABCCT->DefineSection(3, z,  kRB24VMABCCTRin,  kRB24VMABCCTRin + kRB24VMABCCTd);

  TGeoVolume* voRB24VMABCCT = new TGeoVolume("RB24VMABCCT", shRB24VMABCCT, kMedCu);
    
  TGeoVolumeAssembly* voRB24VMABRFCT = new TGeoVolumeAssembly("RB24VMABRFCT");
  voRB24VMABRFCT->AddNode(voRB24VMABCCT,        1, gGeoIdentity);
  voRB24VMABRFCT->AddNode( voRB24VMABCCTFlange, 1, new TGeoTranslation(0., 0.,  kRB24VMABCCTL - kRB24VMABCCTFlangeL));

  z =  kRB24VMABCRBT1L/2. + kRB24B1BellowUndL + kRB24VMABBEConTubeL1 +  kRB24VMABBEConTubeL2 - kRB24VMABCCTL + 1.;    
  voRB24VMABCRB->AddNode(voRB24VMABRFCT, 1, new TGeoTranslation(0., 0., z));


  //
  // Assembling RB24/1
  //    
  TGeoVolumeAssembly* voRB24 = new TGeoVolumeAssembly("RB24");
  // Cu Tube with two simplified flanges
  voRB24->AddNode(voRB24CuTubeM, 1, gGeoIdentity);
  voRB24->AddNode(voRB24CuTubeA, 1, gGeoIdentity);
  z = - kRB24CuTubeL/2 + kRB24CuTubeFL/2.;
  voRB24->AddNode(voRB24CuTubeF, 1, new TGeoTranslation(0., 0., z));
  z = + kRB24CuTubeL/2 - kRB24CuTubeFL/2.;
  voRB24->AddNode(voRB24CuTubeF, 2, new TGeoTranslation(0., 0., z));
  // VMABC close to compensator magnet
  z = - kRB24CuTubeL/2. -  (kRB24VMABCL - kRB24VMABCRBT1L/2) + 1.;
    
  voRB24->AddNode(voRB24VMABCRB, 2, new TGeoTranslation(0., 0., z));
  // Bellow
  z =  kRB24CuTubeL/2;
  voRB24->AddNode(voRB24B1BellowM, 1, new TGeoTranslation(0., 0., z));
  z +=  (kRB24B1L +  kRB24AIpML/2.);
  // Annular ion pump
  voRB24->AddNode(voRB24AIpM, 1, new TGeoTranslation(0., 0., z));
  z +=  (kRB24AIpML/2. +  kRB24ValveWz/2.);
  // Valve
  voRB24->AddNode(voRB24ValveMo, 1, new TGeoTranslation(0., 0., z));
  z += (kRB24ValveWz/2.+ kRB24VMABCRBT1L/2. + 1.);
  // VMABC close to forward detectors
  voRB24->AddNode(voRB24VMABCRB, 3, new TGeoTranslation(0., 0., z));
  //
  //   RB24/2
  //     
  // Copper Tube RB24/2
  const Float_t  kRB242CuTubeL  = 330.0;
    
  TGeoVolume* voRB242CuTubeM = new TGeoVolume("voRB242CuTubeM", 
					      new TGeoTube(0., kRB24CuTubeRo, kRB242CuTubeL/2.), kMedVac);
  voRB24CuTubeM->SetVisibility(0);
  TGeoVolume* voRB242CuTube = new TGeoVolume("voRB242CuTube", 
					     new TGeoTube(kRB24CuTubeRi, kRB24CuTubeRo, kRB242CuTubeL/2.), kMedCu);
  voRB242CuTubeM->AddNode(voRB242CuTube, 1, gGeoIdentity);
    

  TGeoVolumeAssembly* voRB242 = new TGeoVolumeAssembly("RB242");
  voRB242->AddNode(voRB242CuTube, 1, gGeoIdentity);
  z = - kRB242CuTubeL/2 + kRB24CuTubeFL/2.;
  voRB242->AddNode(voRB24CuTubeF, 3, new TGeoTranslation(0., 0., z));
  z = + kRB242CuTubeL/2 - kRB24CuTubeFL/2.;
  voRB242->AddNode(voRB24CuTubeF, 4, new TGeoTranslation(0., 0., z));
  z = - kRB24CuTubeL/2 - kRB24VMABCL - kRB242CuTubeL/2.;
  voRB24->AddNode(voRB242, 1, new TGeoTranslation(0., 0., z));
  //
  //   RB24/3
  //     
  // Copper Tube RB24/3
  const Float_t  kRB243CuTubeL  = 303.35;
    
  TGeoVolume* voRB243CuTubeM = new TGeoVolume("voRB243CuTubeM", 
					      new TGeoTube(0., kRB24CuTubeRo, kRB243CuTubeL/2.), kMedVac);
  voRB24CuTubeM->SetVisibility(0);
  TGeoVolume* voRB243CuTube = new TGeoVolume("voRB243CuTube", 
					     new TGeoTube(kRB24CuTubeRi, kRB24CuTubeRo, kRB243CuTubeL/2.), kMedCu);
  voRB243CuTubeM->AddNode(voRB243CuTube, 1, gGeoIdentity);
    

  TGeoVolumeAssembly* voRB243  = new TGeoVolumeAssembly("RB243");
  TGeoVolumeAssembly* voRB243A = new TGeoVolumeAssembly("RB243A");
    
  voRB243A->AddNode(voRB243CuTube, 1, gGeoIdentity);
  z = - kRB243CuTubeL/2 + kRB24CuTubeFL/2.;
  voRB243A->AddNode(voRB24CuTubeF, 5, new TGeoTranslation(0., 0., z));
  z = + kRB243CuTubeL/2 - kRB24CuTubeFL/2.;
  voRB243A->AddNode(voRB24CuTubeF,    6, new TGeoTranslation(0., 0., z));
  z = + kRB243CuTubeL/2;
  voRB243A->AddNode(voRB24B1BellowM,  2, new TGeoTranslation(0., 0., z));    

  z = - kRB243CuTubeL/2.  - kRB24B1L;
  voRB243->AddNode(voRB243A, 1, new TGeoTranslation(0., 0., z));    
  z = - (1.5 * kRB243CuTubeL + 2. * kRB24B1L);
  voRB243->AddNode(voRB243A, 2, new TGeoTranslation(0., 0., z));    

  z = - 2. * (kRB243CuTubeL + kRB24B1L) - (kRB24VMABCL - kRB24VMABCRBT1L/2) + 1.;
  voRB243->AddNode(voRB24VMABCRB, 3, new TGeoTranslation(0., 0., z));    
    
  z = - kRB24CuTubeL/2 - kRB24VMABCL - kRB242CuTubeL;
  voRB24->AddNode(voRB243, 1, new TGeoTranslation(0., 0., z));


  //
  //
  top->AddNode(voRB24, 1, new TGeoCombiTrans(0., 0., kRB24CuTubeL/2 + 88.5 + 400., rot180));


  // 
  ////////////////////////////////////////////////////////////////////////////////     
  //                                                                            //
  //                                  The Absorber Vacuum system                // 
  //                                                                            //
  ////////////////////////////////////////////////////////////////////////////////
  //
  //    Rotable Flange starts at:            82.00 cm from IP      
  //    Length of rotable flange section:    10.68 cm             
  //    Weld                                  0.08 cm                  
  //    Length of straight section          207.21 cm
  //    =======================================================================
  //                                        299.97 cm  [0.03 cm missing ?]
  //    Length of opening cone              252.09 cm
  //    Weld                                  0.15 cm                
  //    Length of compensator                30.54 cm
  //    Weld                                  0.15 cm                
  //    Length of fixed flange  2.13 - 0.97   1.16 cm
  //    ======================================================================= 
  //                                        584.06 cm [584.80 installed] [0.74 cm missing]
  //    RB26/3
  //    Length of split flange  2.13 - 1.2    0.93 cm
  //    Weld                                  0.15 cm                
  //    Length of fixed point section        16.07 cm               
  //    Weld                                  0.15 cm                
  //    Length of opening cone              629.20 cm
  //    Weld                                  0.30 cm                
  //    Kength of the compensator            41.70 cm
  //    Weld                                  0.30 cm                
  //    Length of fixed flange  2.99 - 1.72   1.27 cm
  // =================================================
  //    Length of RB26/3                    690.07 cm [689.20 installed] [0.87 cm too much] 
  //
  //    RB26/4-5
  //    Length of split flange  2.13 - 1.2    0.93 cm
  //    Weld                                  0.15 cm                
  //    Length of fixed point section        16.07 cm               
  //    Weld                                  0.15 cm                
  //    Length of opening cone              629.20 cm
  //    Weld                                  0.30 cm                
  //    Length of closing cone
  //    Weld
  //    Lenth of straight section 
  //    Kength of the compensator            41.70 cm
  //    Weld                                  0.30 cm                
  //    Length of fixed flange  2.99 - 1.72   1.27 cm
  // =================================================
  //    Length of RB26/3                    690.07 cm [689.20 installed] [0.87 cm too much] 
      
  ///////////////////////////////////////////
  //                                       //
  //    RB26/1-2                           //  
  //    Drawing LHCV2a_0050 [as installed] //
  //    Drawing LHCV2a_0008                //
  //    Drawing LHCV2a_0001                //
  ///////////////////////////////////////////
  //    Pos1 Vacuum Tubes   LHCVC2A__0010
  //    Pos2 Compensator    LHCVC2A__0064
  //    Pos3 Rotable Flange LHCVFX___0016
  //    Pos4 Fixed Flange   LHCVFX___0006
  //    Pos5 Bellow Tooling LHCVFX___0003
  //
  //             
  //
  ///////////////////////////////////
  //    RB26/1-2 Vacuum Tubes      //
  //    Drawing  LHCVC2a_0010      //
  ///////////////////////////////////
  const Float_t kRB26s12TubeL = 459.45; // 0.15 cm added for welding       
  //
  // Add 1 cm on outer diameter for insulation
  //
  TGeoPcon* shRB26s12Tube = new TGeoPcon(0., 360., 5);
  // Section 1: straight section
  shRB26s12Tube->DefineSection(0,   0.00,         5.84/2.,  6.00/2.);
  shRB26s12Tube->DefineSection(1, 207.21,         5.84/2.,  6.00/2.);      
  // Section 2: 0.72 deg opening cone
  shRB26s12Tube->DefineSection(2, 207.21,         5.84/2.,  6.14/2.);      
  shRB26s12Tube->DefineSection(3, 452.30,        12.00/2., 12.30/2.);      
  shRB26s12Tube->DefineSection(4, kRB26s12TubeL, 12.00/2., 12.30/2.); 
  TGeoVolume* voRB26s12Tube  = new TGeoVolume("RB26s12Tube", shRB26s12Tube, kMedSteel);
  // Add the insulation layer    
  TGeoVolume* voRB26s12TubeIns = new TGeoVolume("RB26s12TubeIns", MakeInsulationFromTemplate(shRB26s12Tube), kMedInsu); 
  voRB26s12Tube->AddNode(voRB26s12TubeIns, 1, gGeoIdentity);

 
  TGeoVolume* voRB26s12TubeM  = new TGeoVolume("RB26s12TubeM", MakeMotherFromTemplate(shRB26s12Tube), kMedVac);
  voRB26s12TubeM->AddNode(voRB26s12Tube, 1, gGeoIdentity);
      

      
  ///////////////////////////////////
  //    RB26/2   Axial Compensator //
  //    Drawing  LHCVC2a_0064      //
  ///////////////////////////////////
  const Float_t kRB26s2CompL             = 30.65;    // Length of the compensator
  const Float_t kRB26s2BellowRo          = 14.38/2.; // Bellow outer radius        [Pos 1]
  const Float_t kRB26s2BellowRi          = 12.12/2.; // Bellow inner radius        [Pos 1] 
  const Int_t   kRB26s2NumberOfPlies     = 14;       // Number of plies            [Pos 1] 
  const Float_t kRB26s2BellowUndL        = 10.00;    // Length of undulated region [Pos 1]  [+10 mm installed including pretension ?] 
  const Float_t kRB26s2PlieThickness     =  0.025;   // Plie thickness             [Pos 1]
  const Float_t kRB26s2ConnectionPlieR   =  0.21;    // Connection plie radius     [Pos 1] 
  //  Plie radius
  const Float_t kRB26s2PlieR = 
    (kRB26s2BellowUndL - 4. *  kRB26s2ConnectionPlieR + 2. * kRB26s2PlieThickness + 
     (2. *  kRB26s2NumberOfPlies - 2.) * kRB26s2PlieThickness) / (4. * kRB26s2NumberOfPlies - 2.);
  const Float_t kRB26s2CompTubeInnerR    = 12.00/2.;  // Connection tubes inner radius     [Pos 2 + 3]
  const Float_t kRB26s2CompTubeOuterR    = 12.30/2.;  // Connection tubes outer radius     [Pos 2 + 3]
  const Float_t kRB26s2WeldingTubeLeftL  =  9.00/2.;  // Left connection tube half length  [Pos 2]
  const Float_t kRB26s2WeldingTubeRightL = 11.65/2.;  // Right connection tube half length [Pos 3]  [+ 0.15 cm for welding]
  const Float_t kRB26s2RingOuterR        = 18.10/2.;  // Ring inner radius                 [Pos 4]
  const Float_t kRB26s2RingL             =  0.40/2.;  // Ring half length                  [Pos 4]
  const Float_t kRB26s2RingZ             =  6.50   ;  // Ring z-position                   [Pos 4]
  const Float_t kRB26s2ProtOuterR        = 18.20/2.;  // Protection tube outer radius      [Pos 5]
  const Float_t kRB26s2ProtL             = 15.00/2.;  // Protection tube half length       [Pos 5]
  const Float_t kRB26s2ProtZ             =  6.70   ;  // Protection tube z-position        [Pos 5]
   
      
  // Mother volume
  //
  TGeoPcon* shRB26s2Compensator  = new TGeoPcon(0., 360., 6);
  shRB26s2Compensator->DefineSection( 0,   0.0, 0., kRB26s2CompTubeOuterR);
  shRB26s2Compensator->DefineSection( 1,   kRB26s2RingZ, 0., kRB26s2CompTubeOuterR);      
  shRB26s2Compensator->DefineSection( 2,   kRB26s2RingZ, 0., kRB26s2ProtOuterR);      
  shRB26s2Compensator->DefineSection( 3,   kRB26s2ProtZ + 2. * kRB26s2ProtL, 0., kRB26s2ProtOuterR);            
  shRB26s2Compensator->DefineSection( 4,   kRB26s2ProtZ + 2. * kRB26s2ProtL, 0., kRB26s2CompTubeOuterR);
  shRB26s2Compensator->DefineSection( 5,   kRB26s2CompL                    , 0., kRB26s2CompTubeOuterR);            
  TGeoVolume* voRB26s2Compensator  = new TGeoVolume("RB26s2Compensator", shRB26s2Compensator, kMedVac);
            
  //
  // [Pos 1] Bellow
  //      
  //
  TGeoVolume* voRB26s2Bellow = new TGeoVolume("RB26s2Bellow", new TGeoTube(kRB26s2BellowRi, kRB26s2BellowRo, kRB26s2BellowUndL/2.), kMedVac);
  //      
  //  Upper part of the undulation
  //
  TGeoTorus* shRB26s2PlieTorusU  =  new TGeoTorus(kRB26s2BellowRo - kRB26s2PlieR, kRB26s2PlieR - kRB26s2PlieThickness, kRB26s2PlieR);
  shRB26s2PlieTorusU->SetName("RB26s2TorusU");
  TGeoTube*  shRB26s2PlieTubeU   =  new TGeoTube (kRB26s2BellowRo - kRB26s2PlieR, kRB26s2BellowRo, kRB26s2PlieR);
  shRB26s2PlieTubeU->SetName("RB26s2TubeU");
  TGeoCompositeShape*  shRB26s2UpperPlie = new TGeoCompositeShape("RB26s2UpperPlie", "RB26s2TorusU*RB26s2TubeU");
 
  TGeoVolume* voRB26s2WiggleU = new TGeoVolume("RB26s2UpperPlie", shRB26s2UpperPlie, kMedSteel);
  //
  // Lower part of the undulation
  TGeoTorus* shRB26s2PlieTorusL =  new TGeoTorus(kRB26s2BellowRi + kRB26s2PlieR, kRB26s2PlieR - kRB26s2PlieThickness, kRB26s2PlieR);
  shRB26s2PlieTorusL->SetName("RB26s2TorusL");
  TGeoTube*  shRB26s2PlieTubeL   =  new TGeoTube (kRB26s2BellowRi, kRB26s2BellowRi + kRB26s2PlieR, kRB26s2PlieR);
  shRB26s2PlieTubeL->SetName("RB26s2TubeL");
  TGeoCompositeShape*  shRB26s2LowerPlie = new TGeoCompositeShape("RB26s2LowerPlie", "RB26s2TorusL*RB26s2TubeL");
      
  TGeoVolume* voRB26s2WiggleL = new TGeoVolume("RB26s2LowerPlie", shRB26s2LowerPlie, kMedSteel); 

  //
  // Connection between upper and lower part of undulation
  TGeoVolume* voRB26s2WiggleC1 = new TGeoVolume("RB26s2PlieConn1",  
						new TGeoTube(kRB26s2BellowRi + kRB26s2PlieR, 
							     kRB26s2BellowRo - kRB26s2PlieR, kRB26s2PlieThickness / 2.), kMedSteel);
  //
  // One wiggle
  TGeoVolumeAssembly* voRB26s2Wiggle = new TGeoVolumeAssembly("RB26s2Wiggle");
  z0 =  -  kRB26s2PlieThickness / 2.;
  voRB26s2Wiggle->AddNode(voRB26s2WiggleC1,  1 , new TGeoTranslation(0., 0., z0));
  z0 += kRB26s2PlieR -  kRB26s2PlieThickness / 2.;
  voRB26s2Wiggle->AddNode(voRB26s2WiggleU,   1 , new TGeoTranslation(0., 0., z0));
  z0 += kRB26s2PlieR -  kRB26s2PlieThickness / 2.;
  voRB26s2Wiggle->AddNode(voRB26s2WiggleC1,  2 , new TGeoTranslation(0., 0., z0));
  z0 += kRB26s2PlieR -  kRB26s2PlieThickness;
  voRB26s2Wiggle->AddNode(voRB26s2WiggleL ,  1 , new TGeoTranslation(0., 0., z0));
  // Positioning of the volumes
  z0   = - kRB26s2BellowUndL/2.+ kRB26s2ConnectionPlieR;
  voRB26s2Bellow->AddNode(voRB26s2WiggleL, 1, new TGeoTranslation(0., 0., z0));
  z0  +=  kRB26s2ConnectionPlieR;
  zsh  = 4. *  kRB26s2PlieR -  2. * kRB26s2PlieThickness;
  for (Int_t iw = 0; iw < kRB26s2NumberOfPlies; iw++) {
    Float_t zpos =  z0 + iw * zsh;	
    voRB26s2Bellow->AddNode(voRB26s2Wiggle,  iw + 1, new TGeoTranslation(0., 0., zpos -  kRB26s2PlieThickness));	
  }

  voRB26s2Compensator->AddNode(voRB26s2Bellow, 1,  new TGeoTranslation(0., 0., 2. * kRB26s2WeldingTubeLeftL + kRB26s2BellowUndL/2.));
      
  //
  // [Pos 2] Left Welding Tube
  //      
  TGeoTube* shRB26s2CompLeftTube = new TGeoTube(kRB26s2CompTubeInnerR, kRB26s2CompTubeOuterR, kRB26s2WeldingTubeLeftL);
  TGeoVolume* voRB26s2CompLeftTube = new TGeoVolume("RB26s2CompLeftTube", shRB26s2CompLeftTube, kMedSteel);
  voRB26s2Compensator->AddNode(voRB26s2CompLeftTube, 1,  new TGeoTranslation(0., 0., kRB26s2WeldingTubeLeftL));
  //
  // [Pos 3] Right Welding Tube
  //      
  TGeoTube* shRB26s2CompRightTube = new TGeoTube(kRB26s2CompTubeInnerR, kRB26s2CompTubeOuterR, kRB26s2WeldingTubeRightL);
  TGeoVolume* voRB26s2CompRightTube = new TGeoVolume("RB26s2CompRightTube", shRB26s2CompRightTube, kMedSteel);
  voRB26s2Compensator->AddNode(voRB26s2CompRightTube,  1, new TGeoTranslation(0., 0.,  kRB26s2CompL - kRB26s2WeldingTubeRightL));
  //
  // [Pos 4] Ring
  //      
  TGeoTube* shRB26s2CompRing = new TGeoTube(kRB26s2CompTubeOuterR, kRB26s2RingOuterR, kRB26s2RingL);
  TGeoVolume* voRB26s2CompRing = new TGeoVolume("RB26s2CompRing", shRB26s2CompRing, kMedSteel);
  voRB26s2Compensator->AddNode(voRB26s2CompRing,  1, new TGeoTranslation(0., 0., kRB26s2RingZ + kRB26s2RingL));

  //
  // [Pos 5] Outer Protecting Tube
  //      
  TGeoTube* shRB26s2CompProtTube = new TGeoTube(kRB26s2RingOuterR, kRB26s2ProtOuterR, kRB26s2ProtL);
  TGeoVolume* voRB26s2CompProtTube = new TGeoVolume("RB26s2CompProtTube", shRB26s2CompProtTube, kMedSteel);
  voRB26s2Compensator->AddNode(voRB26s2CompProtTube, 1,  new TGeoTranslation(0., 0., kRB26s2ProtZ + kRB26s2ProtL));
      
  ///////////////////////////////////
  //    Rotable Flange             //
  //    Drawing  LHCVFX_0016       //
  /////////////////////////////////// 
  const Float_t kRB26s1RFlangeTubeRi    = 5.84/2.;  // Tube inner radius
  const Float_t kRB26s1RFlangeTubeRo    = 6.00/2.;  // Tube outer radius

  // Pos 1 Clamp Ring          LHCVFX__0015
  const Float_t kRB26s1RFlangeCrL       = 1.40     ; // Lenth of the clamp ring
  const Float_t kRB26s1RFlangeCrRi1     = 6.72/2.; // Ring inner radius section 1
  const Float_t kRB26s1RFlangeCrRi2     = 6.06/2.; // Ring inner radius section 2
  const Float_t kRB26s1RFlangeCrRo      = 8.60/2.;// Ring outer radius 
  const Float_t kRB26s1RFlangeCrD       = 0.800    ; // Width section 1
      
  TGeoPcon* shRB26s1RFlangeCr = new TGeoPcon(0., 360., 4);
  z0 = 0.;
  shRB26s1RFlangeCr->DefineSection(0, z0, kRB26s1RFlangeCrRi1, kRB26s1RFlangeCrRo);
  z0 += kRB26s1RFlangeCrD;
  shRB26s1RFlangeCr->DefineSection(1, z0, kRB26s1RFlangeCrRi1, kRB26s1RFlangeCrRo);
  shRB26s1RFlangeCr->DefineSection(2, z0, kRB26s1RFlangeCrRi2, kRB26s1RFlangeCrRo);      
  z0 = kRB26s1RFlangeCrL;
  shRB26s1RFlangeCr->DefineSection(3, z0, kRB26s1RFlangeCrRi2, kRB26s1RFlangeCrRo);
  TGeoVolume* voRB26s1RFlangeCr =  
    new TGeoVolume("RB26s1RFlangeCr", shRB26s1RFlangeCr, kMedSteel);

  // Pos 2 Insert              LHCVFX__0015
  const Float_t kRB26s1RFlangeIsL       = 4.88     ; // Lenth of the insert
  const Float_t kRB26s1RFlangeIsR       = 6.70/2.  ; // Ring radius
  const Float_t kRB26s1RFlangeIsD       = 0.80     ; // Ring Width

  TGeoPcon* shRB26s1RFlangeIs = new TGeoPcon(0., 360., 4);
  z0 = 0.;
  shRB26s1RFlangeIs->DefineSection(0, z0, kRB26s1RFlangeTubeRi, kRB26s1RFlangeIsR);
  z0 += kRB26s1RFlangeIsD;
  shRB26s1RFlangeIs->DefineSection(1, z0, kRB26s1RFlangeTubeRi, kRB26s1RFlangeIsR);
  shRB26s1RFlangeIs->DefineSection(2, z0, kRB26s1RFlangeTubeRi, kRB26s1RFlangeTubeRo);      
  z0 = kRB26s1RFlangeIsL;
  shRB26s1RFlangeIs->DefineSection(3, z0, kRB26s1RFlangeTubeRi, kRB26s1RFlangeTubeRo);
  TGeoVolume* voRB26s1RFlangeIs =  
    new TGeoVolume("RB26s1RFlangeIs", shRB26s1RFlangeIs, kMedSteel);
  // 4.88 + 3.7 = 8.58 (8.7 to avoid overlap)
  // Pos 3 Fixed Point Section LHCVC2A_0021
  const Float_t kRB26s1RFlangeFpL       = 5.88     ; // Length of the fixed point section (0.08 cm added for welding)
  const Float_t kRB26s1RFlangeFpZ       = 3.82     ; // Position of the ring
  const Float_t kRB26s1RFlangeFpD       = 0.59     ; // Width of the ring
  const Float_t kRB26s1RFlangeFpR       = 7.00/2.  ; // Radius of the ring
      
  TGeoPcon* shRB26s1RFlangeFp = new TGeoPcon(0., 360., 6);
  z0 = 0.;
  shRB26s1RFlangeFp->DefineSection(0, z0, kRB26s1RFlangeTubeRi, kRB26s1RFlangeTubeRo);
  z0 += kRB26s1RFlangeFpZ;
  shRB26s1RFlangeFp->DefineSection(1, z0, kRB26s1RFlangeTubeRi, kRB26s1RFlangeTubeRo);      
  shRB26s1RFlangeFp->DefineSection(2, z0, kRB26s1RFlangeTubeRi, kRB26s1RFlangeFpR);      	  
  z0 += kRB26s1RFlangeFpD;
  shRB26s1RFlangeFp->DefineSection(3, z0, kRB26s1RFlangeTubeRi, kRB26s1RFlangeFpR);      	  
  shRB26s1RFlangeFp->DefineSection(4, z0, kRB26s1RFlangeTubeRi, kRB26s1RFlangeTubeRo);
  z0 = kRB26s1RFlangeFpL;
  shRB26s1RFlangeFp->DefineSection(5, z0, kRB26s1RFlangeTubeRi, kRB26s1RFlangeTubeRo);
  TGeoVolume* voRB26s1RFlangeFp = new TGeoVolume("RB26s1RFlangeFp", shRB26s1RFlangeFp, kMedSteel);
      	     
  // Put everything in a mother volume
  TGeoPcon* shRB26s1RFlange = new TGeoPcon(0., 360., 8);
  z0 =  0.;
  shRB26s1RFlange->DefineSection(0, z0, 0., kRB26s1RFlangeCrRo);
  z0 += kRB26s1RFlangeCrL;
  shRB26s1RFlange->DefineSection(1, z0, 0., kRB26s1RFlangeCrRo);
  shRB26s1RFlange->DefineSection(2, z0, 0., kRB26s1RFlangeTubeRo);
  z0 = kRB26s1RFlangeIsL + kRB26s1RFlangeFpZ;
  shRB26s1RFlange->DefineSection(3, z0, 0., kRB26s1RFlangeTubeRo);      
  shRB26s1RFlange->DefineSection(4, z0, 0., kRB26s1RFlangeFpR);
  z0 += kRB26s1RFlangeFpD;
  shRB26s1RFlange->DefineSection(5, z0, 0., kRB26s1RFlangeFpR);      	  
  shRB26s1RFlange->DefineSection(6, z0, 0., kRB26s1RFlangeTubeRo);
  z0 = kRB26s1RFlangeIsL + kRB26s1RFlangeFpL;
  shRB26s1RFlange->DefineSection(7, z0, 0., kRB26s1RFlangeTubeRo);
  TGeoVolume* voRB26s1RFlange = new TGeoVolume("RB26s1RFlange", shRB26s1RFlange, kMedVac);

  voRB26s1RFlange->AddNode(voRB26s1RFlangeIs, 1, gGeoIdentity);
  voRB26s1RFlange->AddNode(voRB26s1RFlangeCr, 1, gGeoIdentity);
  voRB26s1RFlange->AddNode(voRB26s1RFlangeFp, 1, new TGeoTranslation(0., 0., kRB26s1RFlangeIsL));
      
  ///////////////////////////////////
  //    Fixed Flange               //
  //    Drawing  LHCVFX_0006       //
  /////////////////////////////////// 
  const Float_t kRB26s2FFlangeL      =  2.13;    // Length of the flange
  const Float_t kRB26s2FFlangeD1     =  0.97;    // Length of section 1
  const Float_t kRB26s2FFlangeD2     =  0.29;    // Length of section 2						     
  const Float_t kRB26s2FFlangeD3     =  0.87;    // Length of section 3						           
  const Float_t kRB26s2FFlangeRo     = 17.15/2.; // Flange outer radius 
  const Float_t kRB26s2FFlangeRi1    = 12.30/2.; // Flange inner radius section 1
  const Float_t kRB26s2FFlangeRi2    = 12.00/2.; // Flange inner radius section 2
  const Float_t kRB26s2FFlangeRi3    = 12.30/2.; // Flange inner radius section 3
  z0 = 0;
  TGeoPcon* shRB26s2FFlange = new TGeoPcon(0., 360., 6);
  z0 = 0.;
  shRB26s2FFlange->DefineSection(0, z0, kRB26s2FFlangeRi1, kRB26s2FFlangeRo);
  z0 += kRB26s2FFlangeD1;
  shRB26s2FFlange->DefineSection(1, z0, kRB26s2FFlangeRi1, kRB26s2FFlangeRo);
  shRB26s2FFlange->DefineSection(2, z0, kRB26s2FFlangeRi2, kRB26s2FFlangeRo);
  z0 += kRB26s2FFlangeD2;
  shRB26s2FFlange->DefineSection(3, z0, kRB26s2FFlangeRi2, kRB26s2FFlangeRo);
  shRB26s2FFlange->DefineSection(4, z0, kRB26s2FFlangeRi3, kRB26s2FFlangeRo);
  z0 += kRB26s2FFlangeD3;
  shRB26s2FFlange->DefineSection(5, z0, kRB26s2FFlangeRi3, kRB26s2FFlangeRo);
  TGeoVolume* voRB26s2FFlange = new TGeoVolume("RB26s2FFlange", shRB26s2FFlange, kMedSteel);

  TGeoVolume* voRB26s2FFlangeM = new TGeoVolume("RB26s2FFlangeM", MakeMotherFromTemplate(shRB26s2FFlange, 2, 5), kMedVac);
  voRB26s2FFlangeM->AddNode(voRB26s2FFlange, 1, gGeoIdentity);
      
      

  ////////////////////////////////////////
  //                                    //
  //    RB26/3                          //  
  //    Drawing LHCV2a_0048             //
  //    Drawing LHCV2a_0002             //
  ////////////////////////////////////////    
  //
  //    Pos 1 Vacuum Tubes      LHCVC2A__0003
  //    Pos 2 Fixed Point       LHCVFX___0005
  //    Pos 3 Split Flange      LHCVFX___0007
  //    Pos 4 Fixed Flange      LHCVFX___0004
  //    Pos 5 Axial Compensator LHCVC2A__0065
  //
  //
  //
  //
  ///////////////////////////////////
  //    Vacuum Tube                //
  //    Drawing  LHCVC2A_0003      //
  /////////////////////////////////// 
  const Float_t kRB26s3TubeL  = 629.35 + 0.3; // 0.3 cm added for welding
  const Float_t kRB26s3TubeR1 =  12./2.;
  const Float_t kRB26s3TubeR2 =  kRB26s3TubeR1 + 215.8 * TMath::Tan(0.829 / 180. * TMath::Pi());
      
      
  TGeoPcon* shRB26s3Tube = new TGeoPcon(0., 360., 7);
  // Section 1: straight section
  shRB26s3Tube->DefineSection(0,   0.00, kRB26s3TubeR1, kRB26s3TubeR1 + 0.15);
  shRB26s3Tube->DefineSection(1,   2.00, kRB26s3TubeR1, kRB26s3TubeR1 + 0.15);      
  // Section 2: 0.829 deg opening cone
  shRB26s3Tube->DefineSection(2,   2.00, kRB26s3TubeR1, kRB26s3TubeR1 + 0.20);
      
  shRB26s3Tube->DefineSection(3, 217.80, kRB26s3TubeR2, kRB26s3TubeR2 + 0.20);
  shRB26s3Tube->DefineSection(4, 217.80, kRB26s3TubeR2, kRB26s3TubeR2 + 0.30);      

  shRB26s3Tube->DefineSection(5, 622.20,       30.00/2., 30.60/2.);      
  shRB26s3Tube->DefineSection(6, kRB26s3TubeL, 30.00/2., 30.60/2.); 

  TGeoVolume* voRB26s3Tube = new TGeoVolume("RB26s3Tube", shRB26s3Tube, kMedSteel);
  //    Add the insulation layer
  TGeoVolume* voRB26s3TubeIns = new TGeoVolume("RB26s3TubeIns", MakeInsulationFromTemplate(shRB26s3Tube), kMedInsu); 
  voRB26s3Tube->AddNode(voRB26s3TubeIns, 1, gGeoIdentity);

  TGeoVolume* voRB26s3TubeM  = new TGeoVolume("RB26s3TubeM", MakeMotherFromTemplate(shRB26s3Tube), kMedVac);
  voRB26s3TubeM->AddNode(voRB26s3Tube, 1, gGeoIdentity);

      

  ///////////////////////////////////
  //    Fixed Point                //
  //    Drawing  LHCVFX_0005       //
  /////////////////////////////////// 
  const Float_t kRB26s3FixedPointL       = 16.37     ; // Length of the fixed point section (0.3 cm added for welding)
  const Float_t kRB26s3FixedPointZ       =  9.72     ; // Position of the ring (0.15 cm added for welding)
  const Float_t kRB26s3FixedPointD       =  0.595    ; // Width of the ring
  const Float_t kRB26s3FixedPointR       = 13.30/2.  ; // Radius of the ring
  const Float_t kRB26s3FixedPointRi      = 12.00/2.  ; // Inner radius of the tube
  const Float_t kRB26s3FixedPointRo1     = 12.30/2.  ; // Outer radius of the tube (in)
  const Float_t kRB26s3FixedPointRo2     = 12.40/2.  ; // Outer radius of the tube (out)
  const Float_t kRB26s3FixedPointDs      =  1.5      ; // Width of straight section behind ring
  const Float_t kRB26s3FixedPointDc      =  3.15     ; // Width of conical  section behind ring (0.15 cm added for welding)      
      
  TGeoPcon* shRB26s3FixedPoint = new TGeoPcon(0., 360., 8);
  z0 = 0.;
  shRB26s3FixedPoint->DefineSection(0, z0, kRB26s3FixedPointRi, kRB26s3FixedPointRo1);
  z0 += kRB26s3FixedPointZ;
  shRB26s3FixedPoint->DefineSection(1, z0, kRB26s3FixedPointRi, kRB26s3FixedPointRo1);      
  shRB26s3FixedPoint->DefineSection(2, z0, kRB26s3FixedPointRi, kRB26s3FixedPointR);      	  
  z0 += kRB26s3FixedPointD;
  shRB26s3FixedPoint->DefineSection(3, z0, kRB26s3FixedPointRi, kRB26s3FixedPointR);      	  
  shRB26s3FixedPoint->DefineSection(4, z0, kRB26s3FixedPointRi, kRB26s3FixedPointRo1);
  z0 += kRB26s3FixedPointDs;
  shRB26s3FixedPoint->DefineSection(5, z0, kRB26s3FixedPointRi, kRB26s3FixedPointRo1);
  z0 += kRB26s3FixedPointDc;
  shRB26s3FixedPoint->DefineSection(6, z0, kRB26s3FixedPointRi, kRB26s3FixedPointRo2);
  z0 = kRB26s3FixedPointL;
  shRB26s3FixedPoint->DefineSection(7, z0, kRB26s3FixedPointRi, kRB26s3FixedPointRo2);
  TGeoVolume* voRB26s3FixedPoint = new TGeoVolume("RB26s3FixedPoint", shRB26s3FixedPoint, kMedSteel);

  TGeoVolume* voRB26s3FixedPointM = new TGeoVolume("RB26s3FixedPointM", MakeMotherFromTemplate(shRB26s3FixedPoint), kMedVac);
  voRB26s3FixedPointM->AddNode(voRB26s3FixedPoint, 1, gGeoIdentity);
      
  ///////////////////////////////////
  //    Split Flange               //
  //    Drawing  LHCVFX_0005       //
  /////////////////////////////////// 
  const Float_t kRB26s3SFlangeL      =  2.13;        // Length of the flange
  const Float_t kRB26s3SFlangeD1     =  0.57;        // Length of section 1
  const Float_t kRB26s3SFlangeD2     =  0.36;        // Length of section 2						     
  const Float_t kRB26s3SFlangeD3     =  0.50 + 0.70; // Length of section 3						           
  const Float_t kRB26s3SFlangeRo     = 17.15/2.;     // Flange outer radius 
  const Float_t kRB26s3SFlangeRi1    = 12.30/2.;     // Flange inner radius section 1
  const Float_t kRB26s3SFlangeRi2    = 12.00/2.;     // Flange inner radius section 2
  const Float_t kRB26s3SFlangeRi3    = 12.30/2.;     // Flange inner radius section 3
  z0 = 0;
  TGeoPcon* shRB26s3SFlange = new TGeoPcon(0., 360., 6);
  z0 = 0.;
  shRB26s3SFlange->DefineSection(0, z0, kRB26s3SFlangeRi1, kRB26s3SFlangeRo);
  z0 += kRB26s3SFlangeD1;
  shRB26s3SFlange->DefineSection(1, z0, kRB26s3SFlangeRi1, kRB26s3SFlangeRo);
  shRB26s3SFlange->DefineSection(2, z0, kRB26s3SFlangeRi2, kRB26s3SFlangeRo);
  z0 += kRB26s3SFlangeD2;
  shRB26s3SFlange->DefineSection(3, z0, kRB26s3SFlangeRi2, kRB26s3SFlangeRo);
  shRB26s3SFlange->DefineSection(4, z0, kRB26s3SFlangeRi3, kRB26s3SFlangeRo);
  z0 += kRB26s3SFlangeD3;
  shRB26s3SFlange->DefineSection(5, z0, kRB26s3SFlangeRi3, kRB26s3SFlangeRo);
  TGeoVolume* voRB26s3SFlange = new TGeoVolume("RB26s3SFlange", shRB26s3SFlange, kMedSteel);

  TGeoVolume* voRB26s3SFlangeM = new TGeoVolume("RB26s3SFlangeM", MakeMotherFromTemplate(shRB26s3SFlange, 0, 3), kMedVac);
  voRB26s3SFlangeM->AddNode(voRB26s3SFlange, 1, gGeoIdentity);
        
  ///////////////////////////////////
  //    RB26/3   Fixed Flange      //
  //    Drawing  LHCVFX___0004     //
  /////////////////////////////////// 
  const Float_t kRB26s3FFlangeL      =  2.99;    // Length of the flange
  const Float_t kRB26s3FFlangeD1     =  1.72;    // Length of section 1
  const Float_t kRB26s3FFlangeD2     =  0.30;    // Length of section 2						     
  const Float_t kRB26s3FFlangeD3     =  0.97;    // Length of section 3						           
  const Float_t kRB26s3FFlangeRo     = 36.20/2.; // Flange outer radius 
  const Float_t kRB26s3FFlangeRi1    = 30.60/2.; // Flange inner radius section 1
  const Float_t kRB26s3FFlangeRi2    = 30.00/2.; // Flange inner radius section 2
  const Float_t kRB26s3FFlangeRi3    = 30.60/2.; // Flange inner radius section 3
  z0 = 0;
  TGeoPcon* shRB26s3FFlange = new TGeoPcon(0., 360., 6);
  z0 = 0.;
  shRB26s3FFlange->DefineSection(0, z0, kRB26s3FFlangeRi1, kRB26s3FFlangeRo);
  z0 += kRB26s3FFlangeD1;
  shRB26s3FFlange->DefineSection(1, z0, kRB26s3FFlangeRi1, kRB26s3FFlangeRo);
  shRB26s3FFlange->DefineSection(2, z0, kRB26s3FFlangeRi2, kRB26s3FFlangeRo);
  z0 += kRB26s3FFlangeD2;
  shRB26s3FFlange->DefineSection(3, z0, kRB26s3FFlangeRi2, kRB26s3FFlangeRo);
  shRB26s3FFlange->DefineSection(4, z0, kRB26s3FFlangeRi3, kRB26s3FFlangeRo);
  z0 += kRB26s3FFlangeD3;
  shRB26s3FFlange->DefineSection(5, z0, kRB26s3FFlangeRi3, kRB26s3FFlangeRo);
  TGeoVolume* voRB26s3FFlange = new TGeoVolume("RB26s3FFlange", shRB26s3FFlange, kMedSteel);
      
  TGeoVolume* voRB26s3FFlangeM = new TGeoVolume("RB26s3FFlangeM", MakeMotherFromTemplate(shRB26s3FFlange, 2, 5), kMedVac);
  voRB26s3FFlangeM->AddNode(voRB26s3FFlange, 1, gGeoIdentity);
            


  ///////////////////////////////////
  //    RB26/3   Axial Compensator //
  //    Drawing  LHCVC2a_0065      //
  /////////////////////////////////// 
  const Float_t kRB26s3CompL              = 42.0;     // Length of the compensator (0.3 cm added for welding)
  const Float_t kRB26s3BellowRo           = 34.00/2.; // Bellow outer radius        [Pos 1]
  const Float_t kRB26s3BellowRi           = 30.10/2.; // Bellow inner radius        [Pos 1] 
  const Int_t   kRB26s3NumberOfPlies      = 13;       // Number of plies            [Pos 1] 
  const Float_t kRB26s3BellowUndL         = 17.70;    // Length of undulated region [Pos 1] 
  const Float_t kRB26s3PlieThickness      =  0.06;    // Plie thickness             [Pos 1]
  const Float_t kRB26s3ConnectionPlieR    =  0.21;    // Connection plie radius     [Pos 1] 
  //  Plie radius
  const Float_t kRB26s3PlieR = 
    (kRB26s3BellowUndL - 4. *  kRB26s3ConnectionPlieR + 2. * kRB26s3PlieThickness + 
     (2. *  kRB26s3NumberOfPlies - 2.) * kRB26s3PlieThickness) / (4. * kRB26s3NumberOfPlies - 2.);

  //
  // The welding tubes have 3 sections with different radii and 2 transition regions.
  // Section 1: connection to the outside
  // Section 2: commection to the bellow
  // Section 3: between 1 and 2
  const Float_t kRB26s3CompTubeInnerR1    = 30.0/2.;  // Outer Connection tubes inner radius     [Pos 4 + 3]
  const Float_t kRB26s3CompTubeOuterR1    = 30.6/2.;  // Outer Connection tubes outer radius     [Pos 4 + 3]
  const Float_t kRB26s3CompTubeInnerR2    = 29.4/2.;  // Connection tubes inner radius           [Pos 4 + 3]
  const Float_t kRB26s3CompTubeOuterR2    = 30.0/2.;  // Connection tubes outer radius           [Pos 4 + 3]
  const Float_t kRB26s3CompTubeInnerR3    = 30.6/2.;  // Connection tubes inner radius at bellow [Pos 4 + 3]
  const Float_t kRB26s3CompTubeOuterR3    = 32.2/2.;  // Connection tubes outer radius at bellow [Pos 4 + 3]
 
  const Float_t kRB26s3WeldingTubeLeftL1  =  2.0;     // Left connection tube length             [Pos 4]
  const Float_t kRB26s3WeldingTubeLeftL2  =  3.4;     // Left connection tube length             [Pos 4]
  const Float_t kRB26s3WeldingTubeLeftL   =  7.0;     // Left connection tube total length       [Pos 4]
  const Float_t kRB26s3WeldingTubeRightL1 =  2.3;     // Right connection tube length            [Pos 3] (0.3 cm added for welding)
  const Float_t kRB26s3WeldingTubeRightL2 = 13.4;     // Right connection tube length            [Pos 3]

  const Float_t kRB26s3WeldingTubeT1      =  0.6;     // Length of first r-transition            [Pos 4 + 3]
  const Float_t kRB26s3WeldingTubeT2      =  1.0;     // Length of 2nd   r-transition            [Pos 4 + 3]       

      
      
  const Float_t kRB26s3RingOuterR         = 36.1/2.;  // Ring inner radius                       [Pos 4]
  const Float_t kRB26s3RingL              =  0.8/2.;  // Ring half length                        [Pos 4]
  const Float_t kRB26s3RingZ              =  3.7   ;  // Ring z-position                         [Pos 4]
  const Float_t kRB26s3ProtOuterR         = 36.2/2.;  // Protection tube outer radius            [Pos 2]
  const Float_t kRB26s3ProtL              = 27.0/2.;  // Protection tube half length             [Pos 2]
  const Float_t kRB26s3ProtZ              =  4.0   ;  // Protection tube z-position              [Pos 2]
   
      
  // Mother volume
  //
  TGeoPcon* shRB26s3Compensator  = new TGeoPcon(0., 360., 6);
  shRB26s3Compensator->DefineSection( 0,   0.0, 0., kRB26s3CompTubeOuterR1);
  shRB26s3Compensator->DefineSection( 1,   kRB26s3RingZ, 0., kRB26s3CompTubeOuterR1);      
  shRB26s3Compensator->DefineSection( 2,   kRB26s3RingZ, 0., kRB26s3ProtOuterR);      
  shRB26s3Compensator->DefineSection( 3,   kRB26s3ProtZ + 2. * kRB26s3ProtL, 0., kRB26s3ProtOuterR);            
  shRB26s3Compensator->DefineSection( 4,   kRB26s3ProtZ + 2. * kRB26s3ProtL, 0., kRB26s3CompTubeOuterR1);
  shRB26s3Compensator->DefineSection( 5,   kRB26s3CompL                    , 0., kRB26s3CompTubeOuterR1);            
  TGeoVolume* voRB26s3Compensator  =  
    new TGeoVolume("RB26s3Compensator", shRB26s3Compensator, kMedVac);
            
  //
  // [Pos 1] Bellow
  //      
  //
  TGeoVolume* voRB26s3Bellow = new TGeoVolume("RB26s3Bellow", 
					      new TGeoTube(kRB26s3BellowRi, kRB26s3BellowRo, kRB26s3BellowUndL/2.), kMedVac);
  //      
  //  Upper part of the undulation
  //
  TGeoTorus* shRB26s3PlieTorusU  =  new TGeoTorus(kRB26s3BellowRo - kRB26s3PlieR, kRB26s3PlieR - kRB26s3PlieThickness, kRB26s3PlieR);
  shRB26s3PlieTorusU->SetName("RB26s3TorusU");
  TGeoTube*  shRB26s3PlieTubeU   =  new TGeoTube (kRB26s3BellowRo - kRB26s3PlieR, kRB26s3BellowRo, kRB26s3PlieR);
  shRB26s3PlieTubeU->SetName("RB26s3TubeU");
  TGeoCompositeShape*  shRB26s3UpperPlie = new TGeoCompositeShape("RB26s3UpperPlie", "RB26s3TorusU*RB26s3TubeU");
 
  TGeoVolume* voRB26s3WiggleU = new TGeoVolume("RB26s3UpperPlie", shRB26s3UpperPlie, kMedSteel);
  //
  // Lower part of the undulation
  TGeoTorus* shRB26s3PlieTorusL =  new TGeoTorus(kRB26s3BellowRi + kRB26s3PlieR, kRB26s3PlieR - kRB26s3PlieThickness, kRB26s3PlieR);
  shRB26s3PlieTorusL->SetName("RB26s3TorusL");
  TGeoTube*  shRB26s3PlieTubeL   =  new TGeoTube (kRB26s3BellowRi, kRB26s3BellowRi + kRB26s3PlieR, kRB26s3PlieR);
  shRB26s3PlieTubeL->SetName("RB26s3TubeL");
  TGeoCompositeShape*  shRB26s3LowerPlie = new TGeoCompositeShape("RB26s3LowerPlie", "RB26s3TorusL*RB26s3TubeL");
      
  TGeoVolume* voRB26s3WiggleL = new TGeoVolume("RB26s3LowerPlie", shRB26s3LowerPlie, kMedSteel); 

  //
  // Connection between upper and lower part of undulation
  TGeoVolume* voRB26s3WiggleC1 = new TGeoVolume("RB26s3PlieConn1",  
						new TGeoTube(kRB26s3BellowRi + kRB26s3PlieR, 
							     kRB26s3BellowRo - kRB26s3PlieR, kRB26s3PlieThickness / 2.), kMedSteel);
  //
  // One wiggle
  TGeoVolumeAssembly* voRB26s3Wiggle = new TGeoVolumeAssembly("RB26s3Wiggle");
  z0 =  -  kRB26s3PlieThickness / 2.;
  voRB26s3Wiggle->AddNode(voRB26s3WiggleC1,  1 , new TGeoTranslation(0., 0., z0));
  z0 += kRB26s3PlieR -  kRB26s3PlieThickness / 2.;
  voRB26s3Wiggle->AddNode(voRB26s3WiggleU,   1 , new TGeoTranslation(0., 0., z0));
  z0 += kRB26s3PlieR -  kRB26s3PlieThickness / 2.;
  voRB26s3Wiggle->AddNode(voRB26s3WiggleC1,  2 , new TGeoTranslation(0., 0., z0));
  z0 += kRB26s3PlieR -  kRB26s3PlieThickness;
  voRB26s3Wiggle->AddNode(voRB26s3WiggleL,  1 , new TGeoTranslation(0., 0., z0));
  // Positioning of the volumes
  z0   = - kRB26s3BellowUndL/2.+ kRB26s3ConnectionPlieR;
  voRB26s3Bellow->AddNode(voRB26s3WiggleL, 1, new TGeoTranslation(0., 0., z0));
  z0  +=  kRB26s3ConnectionPlieR;
  zsh  = 4. *  kRB26s3PlieR -  2. * kRB26s3PlieThickness;
  for (Int_t iw = 0; iw < kRB26s3NumberOfPlies; iw++) {
    Float_t zpos =  z0 + iw * zsh;	
    voRB26s3Bellow->AddNode(voRB26s3Wiggle,  iw + 1, new TGeoTranslation(0., 0., zpos -  kRB26s3PlieThickness));	
  }

  voRB26s3Compensator->AddNode(voRB26s3Bellow, 1,  new TGeoTranslation(0., 0., kRB26s3WeldingTubeLeftL + kRB26s3BellowUndL/2.));


  //
  // [Pos 2] Outer Protecting Tube
  //      
  TGeoTube* shRB26s3CompProtTube = new TGeoTube(kRB26s3RingOuterR, kRB26s3ProtOuterR, kRB26s3ProtL);
  TGeoVolume* voRB26s3CompProtTube =  
    new TGeoVolume("RB26s3CompProtTube", shRB26s3CompProtTube, kMedSteel);
  voRB26s3Compensator->AddNode(voRB26s3CompProtTube, 1,  new TGeoTranslation(0., 0., kRB26s3ProtZ + kRB26s3ProtL));
      

  //
  // [Pos 3] Right Welding Tube
  //      
  TGeoPcon* shRB26s3CompRightTube = new TGeoPcon(0., 360., 5);
  z0 = 0.;
  shRB26s3CompRightTube->DefineSection(0, z0,  kRB26s3CompTubeInnerR3, kRB26s3CompTubeOuterR3);
  z0 += kRB26s3WeldingTubeT2;
  shRB26s3CompRightTube->DefineSection(1, z0,  kRB26s3CompTubeInnerR2, kRB26s3CompTubeOuterR2);
  z0 += kRB26s3WeldingTubeRightL2;
  shRB26s3CompRightTube->DefineSection(2, z0,  kRB26s3CompTubeInnerR2, kRB26s3CompTubeOuterR2);
  z0 += kRB26s3WeldingTubeT1;
  shRB26s3CompRightTube->DefineSection(3, z0,  kRB26s3CompTubeInnerR1, kRB26s3CompTubeOuterR1);
  z0 += kRB26s3WeldingTubeRightL1;
  shRB26s3CompRightTube->DefineSection(4, z0,  kRB26s3CompTubeInnerR1, kRB26s3CompTubeOuterR1);
      
  TGeoVolume* voRB26s3CompRightTube =  
    new TGeoVolume("RB26s3CompRightTube", shRB26s3CompRightTube, kMedSteel);
  voRB26s3Compensator->AddNode(voRB26s3CompRightTube,  1, new TGeoTranslation(0., 0.,  kRB26s3CompL - z0));

  //
  // [Pos 4] Left Welding Tube
  //      
  TGeoPcon* shRB26s3CompLeftTube = new TGeoPcon(0., 360., 5);
  z0 = 0.;
  shRB26s3CompLeftTube->DefineSection(0, z0,  kRB26s3CompTubeInnerR1, kRB26s3CompTubeOuterR1);
  z0 += kRB26s3WeldingTubeLeftL1;
  shRB26s3CompLeftTube->DefineSection(1, z0,  kRB26s3CompTubeInnerR1, kRB26s3CompTubeOuterR1);
  z0 += kRB26s3WeldingTubeT1;
  shRB26s3CompLeftTube->DefineSection(2, z0,  kRB26s3CompTubeInnerR2, kRB26s3CompTubeOuterR2);
  z0 += kRB26s3WeldingTubeLeftL2;
  shRB26s3CompLeftTube->DefineSection(3, z0,  kRB26s3CompTubeInnerR2, kRB26s3CompTubeOuterR2);
  z0 += kRB26s3WeldingTubeT2;
  shRB26s3CompLeftTube->DefineSection(4, z0,  kRB26s3CompTubeInnerR3, kRB26s3CompTubeOuterR3);

  TGeoVolume* voRB26s3CompLeftTube =  
    new TGeoVolume("RB26s3CompLeftTube", shRB26s3CompLeftTube, kMedSteel);
  voRB26s3Compensator->AddNode(voRB26s3CompLeftTube, 1,  gGeoIdentity);
  //
  // [Pos 5] Ring
  //      
  TGeoTube* shRB26s3CompRing = new TGeoTube(kRB26s3CompTubeOuterR2, kRB26s3RingOuterR, kRB26s3RingL);
  TGeoVolume* voRB26s3CompRing =  
    new TGeoVolume("RB26s3CompRing", shRB26s3CompRing, kMedSteel);
  voRB26s3Compensator->AddNode(voRB26s3CompRing,  1, new TGeoTranslation(0., 0., kRB26s3RingZ + kRB26s3RingL));



  ///////////////////////////////////////////
  //                                       //
  //    RB26/4-5                           //  
  //    Drawing LHCV2a_0012 [as installed] //
  ////////////////////////////////////////////
  //    Pos1 Vacuum Tubes        LHCVC2A__0014
  //    Pos2 Compensator         LHCVC2A__0066
  //    Pos3 Fixed Point Section LHCVC2A__0016
  //    Pos4 Split Flange        LHCVFX___0005
  //    Pos5 RotableFlange       LHCVFX___0009
  ////////////////////////////////////////////

  ///////////////////////////////////
  //    RB26/4-5 Vacuum Tubes      //
  //    Drawing  LHCVC2a_0014      //
  /////////////////////////////////// 
  const Float_t kRB26s45TubeL = 593.12 + 0.3; // 0.3 cm added for welding
      
  TGeoPcon* shRB26s45Tube = new TGeoPcon(0., 360., 11);
  // Section 1: straight section
  shRB26s45Tube->DefineSection( 0,   0.00, 30.00/2., 30.60/2.);
  shRB26s45Tube->DefineSection( 1,   1.20, 30.00/2., 30.60/2.);
  shRB26s45Tube->DefineSection( 2,   1.20, 30.00/2., 30.80/2.);
  shRB26s45Tube->DefineSection( 3,  25.10, 30.00/2., 30.80/2.);      
  // Section 2: 0.932 deg opening cone
  shRB26s45Tube->DefineSection( 4, 486.10, 45.00/2., 45.80/2.);      
  // Section 3: straight section 4 mm 
  shRB26s45Tube->DefineSection( 5, 512.10, 45.00/2., 45.80/2.);
  // Section 4: straight section 3 mm
  shRB26s45Tube->DefineSection( 6, 512.10, 45.00/2., 45.60/2.);
  shRB26s45Tube->DefineSection( 7, 527.70, 45.00/2., 45.60/2.);
  // Section 4: closing cone 
  shRB26s45Tube->DefineSection( 8, 591.30, 10.00/2., 10.60/2.);      
  shRB26s45Tube->DefineSection( 9, 591.89, 10.00/2., 10.30/2.);      

  shRB26s45Tube->DefineSection(10, kRB26s45TubeL, 10.00/2., 10.30/2.);      
  TGeoVolume* voRB26s45Tube  =  
    new TGeoVolume("RB26s45Tube", shRB26s45Tube, kMedSteel);

  TGeoVolume* voRB26s45TubeM  = new TGeoVolume("RB26s45TubeM", MakeMotherFromTemplate(shRB26s45Tube), kMedVac);
  voRB26s45TubeM->AddNode(voRB26s45Tube, 1, gGeoIdentity);
            
      

  ///////////////////////////////////
  //    RB26/5   Axial Compensator //
  //    Drawing  LHCVC2a_0066      //
  /////////////////////////////////// 
  const Float_t kRB26s5CompL             = 27.60;    // Length of the compensator (0.30 cm added for welding)
  const Float_t kRB26s5BellowRo          = 12.48/2.; // Bellow outer radius        [Pos 1]
  const Float_t kRB26s5BellowRi          = 10.32/2.; // Bellow inner radius        [Pos 1] 
  const Int_t   kRB26s5NumberOfPlies     = 15;       // Number of plies            [Pos 1] 
  const Float_t kRB26s5BellowUndL        = 10.50;    // Length of undulated region [Pos 1] 
  const Float_t kRB26s5PlieThickness     =  0.025;   // Plie thickness             [Pos 1]
  const Float_t kRB26s5ConnectionPlieR   =  0.21;    // Connection plie radius     [Pos 1] 
  const Float_t kRB26s5ConnectionR       = 11.2/2.;  // Bellow connection radius   [Pos 1] 
  //  Plie radius
  const Float_t kRB26s5PlieR = 
    (kRB26s5BellowUndL - 4. *  kRB26s5ConnectionPlieR + 2. * kRB26s5PlieThickness + 
     (2. *  kRB26s5NumberOfPlies - 2.) * kRB26s5PlieThickness) / (4. * kRB26s5NumberOfPlies - 2.);
  const Float_t kRB26s5CompTubeInnerR    = 10.00/2.;  // Connection tubes inner radius     [Pos 2 + 3]
  const Float_t kRB26s5CompTubeOuterR    = 10.30/2.;  // Connection tubes outer radius     [Pos 2 + 3]
  const Float_t kRB26s5WeldingTubeLeftL  =  3.70/2.;  // Left connection tube half length  [Pos 2]
  const Float_t kRB26s5WeldingTubeRightL = 13.40/2.;  // Right connection tube half length [Pos 3]   (0.3 cm added for welding)
  const Float_t kRB26s5RingInnerR        = 11.2/2.;   // Ring inner radius                 [Pos 4]
  const Float_t kRB26s5RingOuterR        = 16.0/2.;   // Ring inner radius                 [Pos 4]
  const Float_t kRB26s5RingL             =  0.4/2.;   // Ring half length                  [Pos 4]
  const Float_t kRB26s5RingZ             = 14.97;     // Ring z-position                   [Pos 4]
  const Float_t kRB26s5ProtOuterR        = 16.2/2.;   // Protection tube outer radius      [Pos 5]
  const Float_t kRB26s5ProtL             = 13.0/2.;   // Protection tube half length       [Pos 5]
  const Float_t kRB26s5ProtZ             =  2.17;     // Protection tube z-position        [Pos 5]
  const Float_t kRB26s5DetailZR          = 11.3/2.;   // Detail Z max radius
      
      
  // Mother volume
  //
  TGeoPcon* shRB26s5Compensator  = new TGeoPcon(0., 360., 8);
  shRB26s5Compensator->DefineSection( 0,   0.0,                                                  0., kRB26s5CompTubeOuterR);
  shRB26s5Compensator->DefineSection( 1,   kRB26s5ProtZ,                                         0., kRB26s5CompTubeOuterR);      
  shRB26s5Compensator->DefineSection( 2,   kRB26s5ProtZ,                                         0., kRB26s5ProtOuterR);
  shRB26s5Compensator->DefineSection( 3,   kRB26s5ProtZ + 2. * kRB26s5ProtL + 2. * kRB26s5RingL, 0., kRB26s5ProtOuterR);      
  shRB26s5Compensator->DefineSection( 4,   kRB26s5ProtZ + 2. * kRB26s5ProtL + 2. * kRB26s5RingL, 0., kRB26s5DetailZR);
  shRB26s5Compensator->DefineSection( 5,   kRB26s5CompL - 8.,                                    0., kRB26s5DetailZR);
  shRB26s5Compensator->DefineSection( 6,   kRB26s5CompL - 8.,                                    0., kRB26s5CompTubeOuterR);            
  shRB26s5Compensator->DefineSection( 7,   kRB26s5CompL,                                         0., kRB26s5CompTubeOuterR);            
  TGeoVolume* voRB26s5Compensator  = new TGeoVolume("RB26s5Compensator", shRB26s5Compensator, kMedVac);
            
  //
  // [Pos 1] Bellow
  //      
  //
  TGeoVolume* voRB26s5Bellow = new TGeoVolume("RB26s5Bellow", 
					      new TGeoTube(kRB26s5BellowRi, kRB26s5BellowRo, kRB26s5BellowUndL/2.), kMedVac);
  //      
  //  Upper part of the undulation
  //
  TGeoTorus* shRB26s5PlieTorusU  =  new TGeoTorus(kRB26s5BellowRo - kRB26s5PlieR, kRB26s5PlieR - kRB26s5PlieThickness, kRB26s5PlieR);
  shRB26s5PlieTorusU->SetName("RB26s5TorusU");
  TGeoTube*  shRB26s5PlieTubeU   =  new TGeoTube (kRB26s5BellowRo - kRB26s5PlieR, kRB26s5BellowRo, kRB26s5PlieR);
  shRB26s5PlieTubeU->SetName("RB26s5TubeU");
  TGeoCompositeShape*  shRB26s5UpperPlie = new TGeoCompositeShape("RB26s5UpperPlie", "RB26s5TorusU*RB26s5TubeU");
 
  TGeoVolume* voRB26s5WiggleU = new TGeoVolume("RB26s5UpperPlie", shRB26s5UpperPlie, kMedSteel);
  //
  // Lower part of the undulation
  TGeoTorus* shRB26s5PlieTorusL =  new TGeoTorus(kRB26s5BellowRi + kRB26s5PlieR, kRB26s5PlieR - kRB26s5PlieThickness, kRB26s5PlieR);
  shRB26s5PlieTorusL->SetName("RB26s5TorusL");
  TGeoTube*  shRB26s5PlieTubeL   =  new TGeoTube (kRB26s5BellowRi, kRB26s5BellowRi + kRB26s5PlieR, kRB26s5PlieR);
  shRB26s5PlieTubeL->SetName("RB26s5TubeL");
  TGeoCompositeShape*  shRB26s5LowerPlie = new TGeoCompositeShape("RB26s5LowerPlie", "RB26s5TorusL*RB26s5TubeL");
      
  TGeoVolume* voRB26s5WiggleL = new TGeoVolume("RB26s5LowerPlie", shRB26s5LowerPlie, kMedSteel); 

  //
  // Connection between upper and lower part of undulation
  TGeoVolume* voRB26s5WiggleC1 = new TGeoVolume("RB26s5PlieConn1",  
						new TGeoTube(kRB26s5BellowRi + kRB26s5PlieR, 
							     kRB26s5BellowRo - kRB26s5PlieR, kRB26s5PlieThickness / 2.), kMedSteel);
  //
  // One wiggle
  TGeoVolumeAssembly* voRB26s5Wiggle = new TGeoVolumeAssembly("RB26s5Wiggle");
  z0 =  -  kRB26s5PlieThickness / 2.;
  voRB26s5Wiggle->AddNode(voRB26s5WiggleC1,  1 , new TGeoTranslation(0., 0., z0));
  z0 += kRB26s5PlieR -  kRB26s5PlieThickness / 2.;
  voRB26s5Wiggle->AddNode(voRB26s5WiggleU,   1 , new TGeoTranslation(0., 0., z0));
  z0 += kRB26s5PlieR -  kRB26s5PlieThickness / 2.;
  voRB26s5Wiggle->AddNode(voRB26s5WiggleC1,  2 , new TGeoTranslation(0., 0., z0));
  z0 += kRB26s5PlieR -  kRB26s5PlieThickness;
  voRB26s5Wiggle->AddNode(voRB26s5WiggleL ,  1 , new TGeoTranslation(0., 0., z0));
  // Positioning of the volumes
  z0   = - kRB26s5BellowUndL/2.+ kRB26s5ConnectionPlieR;
  voRB26s5Bellow->AddNode(voRB26s5WiggleL, 1, new TGeoTranslation(0., 0., z0));
  z0  +=  kRB26s5ConnectionPlieR;
  zsh  = 4. *  kRB26s5PlieR -  2. * kRB26s5PlieThickness;
  for (Int_t iw = 0; iw < kRB26s5NumberOfPlies; iw++) {
    Float_t zpos =  z0 + iw * zsh;	
    voRB26s5Bellow->AddNode(voRB26s5Wiggle,  iw + 1, new TGeoTranslation(0., 0., zpos -  kRB26s5PlieThickness));	
  }

  voRB26s5Compensator->AddNode(voRB26s5Bellow, 1,  new TGeoTranslation(0., 0., 2. * kRB26s5WeldingTubeLeftL + kRB26s5BellowUndL/2.));
      
  //
  // [Pos 2] Left Welding Tube
  //      
  TGeoPcon* shRB26s5CompLeftTube = new TGeoPcon(0., 360., 3);
  z0 = 0;
  shRB26s5CompLeftTube->DefineSection(0, z0, kRB26s5CompTubeInnerR, kRB26s5CompTubeOuterR);
  z0 += 2 * kRB26s5WeldingTubeLeftL - ( kRB26s5ConnectionR - kRB26s5CompTubeOuterR);
  shRB26s5CompLeftTube->DefineSection(1, z0, kRB26s5CompTubeInnerR, kRB26s5CompTubeOuterR);
  z0 += ( kRB26s5ConnectionR - kRB26s5CompTubeOuterR);
  shRB26s5CompLeftTube->DefineSection(2, z0, kRB26s5ConnectionR - 0.15, kRB26s5ConnectionR);
  TGeoVolume* voRB26s5CompLeftTube = new TGeoVolume("RB26s5CompLeftTube", shRB26s5CompLeftTube, kMedSteel);
  voRB26s5Compensator->AddNode(voRB26s5CompLeftTube, 1,  gGeoIdentity);
  //
  // [Pos 3] Right Welding Tube
  //      
  TGeoPcon* shRB26s5CompRightTube = new TGeoPcon(0., 360., 11);
  // Detail Z
  shRB26s5CompRightTube->DefineSection( 0, 0.  , kRB26s5CompTubeInnerR + 0.22, 11.2/2.);
  shRB26s5CompRightTube->DefineSection( 1, 0.05, kRB26s5CompTubeInnerR + 0.18, 11.2/2.);
  shRB26s5CompRightTube->DefineSection( 2, 0.22, kRB26s5CompTubeInnerR       , 11.2/2. - 0.22);
  shRB26s5CompRightTube->DefineSection( 3, 0.44, kRB26s5CompTubeInnerR       , 11.2/2.);
  shRB26s5CompRightTube->DefineSection( 4, 1.70, kRB26s5CompTubeInnerR       , 11.2/2.);
  shRB26s5CompRightTube->DefineSection( 5, 2.10, kRB26s5CompTubeInnerR       , kRB26s5CompTubeOuterR);
  shRB26s5CompRightTube->DefineSection( 6, 2.80, kRB26s5CompTubeInnerR       , kRB26s5CompTubeOuterR);
  shRB26s5CompRightTube->DefineSection( 7, 2.80, kRB26s5CompTubeInnerR       , 11.3/2.);
  shRB26s5CompRightTube->DefineSection( 8, 3.40, kRB26s5CompTubeInnerR       , 11.3/2.);
  // Normal pipe
  shRB26s5CompRightTube->DefineSection( 9, 3.50, kRB26s5CompTubeInnerR       , kRB26s5CompTubeOuterR);
  shRB26s5CompRightTube->DefineSection(10, 2. * kRB26s5WeldingTubeRightL, kRB26s5CompTubeInnerR, kRB26s5CompTubeOuterR);
      
  TGeoVolume* voRB26s5CompRightTube =  
    new TGeoVolume("RB26s5CompRightTube", shRB26s5CompRightTube, kMedSteel);
  voRB26s5Compensator->AddNode(voRB26s5CompRightTube,  1, 
			       new TGeoTranslation(0., 0.,  kRB26s5CompL - 2. * kRB26s5WeldingTubeRightL));
  //
  // [Pos 4] Ring
  //      
  TGeoTube* shRB26s5CompRing = new TGeoTube(kRB26s5RingInnerR, kRB26s5RingOuterR, kRB26s5RingL);
  TGeoVolume* voRB26s5CompRing =  
    new TGeoVolume("RB26s5CompRing", shRB26s5CompRing, kMedSteel);
  voRB26s5Compensator->AddNode(voRB26s5CompRing,  1, new TGeoTranslation(0., 0., kRB26s5RingZ + kRB26s5RingL));

  //
  // [Pos 5] Outer Protecting Tube
  //      
  TGeoTube* shRB26s5CompProtTube = new TGeoTube(kRB26s5RingOuterR, kRB26s5ProtOuterR, kRB26s5ProtL);
  TGeoVolume* voRB26s5CompProtTube =  
    new TGeoVolume("RB26s5CompProtTube", shRB26s5CompProtTube, kMedSteel);
  voRB26s5Compensator->AddNode(voRB26s5CompProtTube, 1,  new TGeoTranslation(0., 0., kRB26s5ProtZ + kRB26s5ProtL));

  ///////////////////////////////////////
  //    RB26/4   Fixed Point Section   //
  //    Drawing  LHCVC2a_0016          //
  /////////////////////////////////////// 
  const Float_t kRB26s4TubeRi            =  30.30/2. ; // Tube inner radius  (0.3 cm added for welding)
  const Float_t kRB26s4TubeRo            =  30.60/2. ; // Tube outer radius      
  const Float_t kRB26s4FixedPointL       =  12.63    ; // Length of the fixed point section
  const Float_t kRB26s4FixedPointZ       =  10.53    ; // Position of the ring (0.15 added for welding)
  const Float_t kRB26s4FixedPointD       =   0.595   ; // Width of the ring
  const Float_t kRB26s4FixedPointR       =  31.60/2. ; // Radius of the ring
      
  TGeoPcon* shRB26s4FixedPoint = new TGeoPcon(0., 360., 6);
  z0 = 0.;
  shRB26s4FixedPoint->DefineSection(0, z0, kRB26s4TubeRi, kRB26s4TubeRo);
  z0 += kRB26s4FixedPointZ;
  shRB26s4FixedPoint->DefineSection(1, z0, kRB26s4TubeRi, kRB26s4TubeRo);      
  shRB26s4FixedPoint->DefineSection(2, z0, kRB26s4TubeRi, kRB26s4FixedPointR);      	  
  z0 += kRB26s4FixedPointD;
  shRB26s4FixedPoint->DefineSection(3, z0, kRB26s4TubeRi, kRB26s4FixedPointR);      	  
  shRB26s4FixedPoint->DefineSection(4, z0, kRB26s4TubeRi, kRB26s4TubeRo);
  z0 = kRB26s4FixedPointL;
  shRB26s4FixedPoint->DefineSection(5, z0, kRB26s4TubeRi, kRB26s4TubeRo);
  TGeoVolume* voRB26s4FixedPoint = new TGeoVolume("RB26s4FixedPoint", shRB26s4FixedPoint, kMedSteel);
      
  TGeoVolume* voRB26s4FixedPointM = new TGeoVolume("RB26s4FixedPointM", MakeMotherFromTemplate(shRB26s4FixedPoint), kMedVac);
  voRB26s4FixedPointM->AddNode(voRB26s4FixedPoint, 1, gGeoIdentity);
            

  ///////////////////////////////////////
  //    RB26/4   Split Flange          //
  //    Drawing  LHCVFX__0005          //
  /////////////////////////////////////// 
  const Float_t kRB26s4SFlangeL      =  2.99;        // Length of the flange
  const Float_t kRB26s4SFlangeD1     =  0.85;        // Length of section 1
  const Float_t kRB26s4SFlangeD2     =  0.36;        // Length of section 2						     
  const Float_t kRB26s4SFlangeD3     =  0.73 + 1.05; // Length of section 3						           
  const Float_t kRB26s4SFlangeRo     = 36.20/2.;     // Flange outer radius 
  const Float_t kRB26s4SFlangeRi1    = 30.60/2.;     // Flange inner radius section 1
  const Float_t kRB26s4SFlangeRi2    = 30.00/2.;     // Flange inner radius section 2
  const Float_t kRB26s4SFlangeRi3    = 30.60/2.;     // Flange inner radius section 3
  z0 = 0;
  TGeoPcon* shRB26s4SFlange = new TGeoPcon(0., 360., 6);
  z0 = 0.;
  shRB26s4SFlange->DefineSection(0, z0, kRB26s4SFlangeRi1, kRB26s4SFlangeRo);
  z0 += kRB26s4SFlangeD1;
  shRB26s4SFlange->DefineSection(1, z0, kRB26s4SFlangeRi1, kRB26s4SFlangeRo);
  shRB26s4SFlange->DefineSection(2, z0, kRB26s4SFlangeRi2, kRB26s4SFlangeRo);
  z0 += kRB26s4SFlangeD2;
  shRB26s4SFlange->DefineSection(3, z0, kRB26s4SFlangeRi2, kRB26s4SFlangeRo);
  shRB26s4SFlange->DefineSection(4, z0, kRB26s4SFlangeRi3, kRB26s4SFlangeRo);
  z0 += kRB26s4SFlangeD3;
  shRB26s4SFlange->DefineSection(5, z0, kRB26s4SFlangeRi3, kRB26s4SFlangeRo);
  TGeoVolume* voRB26s4SFlange = new TGeoVolume("RB26s4SFlange", shRB26s4SFlange, kMedSteel);

  TGeoVolume* voRB26s4SFlangeM = new TGeoVolume("RB26s4SFlangeM", MakeMotherFromTemplate(shRB26s4SFlange, 0, 3), kMedVac);
  voRB26s4SFlangeM->AddNode(voRB26s4SFlange, 1, gGeoIdentity);
      
  ///////////////////////////////////////
  //    RB26/5   Rotable Flange        //
  //    Drawing  LHCVFX__0009          //
  /////////////////////////////////////// 
  const Float_t kRB26s5RFlangeL      =  1.86;    // Length of the flange
  const Float_t kRB26s5RFlangeD1     =  0.61;    // Length of section 1
  const Float_t kRB26s5RFlangeD2     =  0.15;    // Length of section 2						     
  const Float_t kRB26s5RFlangeD3     =  0.60;    // Length of section 3						           
  const Float_t kRB26s5RFlangeD4     =  0.50;    // Length of section 4						           
  const Float_t kRB26s5RFlangeRo     = 15.20/2.; // Flange outer radius 
  const Float_t kRB26s5RFlangeRi1    = 10.30/2.; // Flange inner radius section 1
  const Float_t kRB26s5RFlangeRi2    = 10.00/2.; // Flange inner radius section 2
  const Float_t kRB26s5RFlangeRi3    = 10.30/2.; // Flange inner radius section 3
  const Float_t kRB26s5RFlangeRi4    = 10.50/2.; // Flange inner radius section 4

  z0 = 0;
  TGeoPcon* shRB26s5RFlange = new TGeoPcon(0., 360., 8);
  z0 = 0.;
  shRB26s5RFlange->DefineSection(0, z0, kRB26s5RFlangeRi4, kRB26s5RFlangeRo);
  z0 += kRB26s5RFlangeD4;
  shRB26s5RFlange->DefineSection(1, z0, kRB26s5RFlangeRi4, kRB26s5RFlangeRo);
  shRB26s5RFlange->DefineSection(2, z0, kRB26s5RFlangeRi3, kRB26s5RFlangeRo);
  z0 += kRB26s5RFlangeD3;
  shRB26s5RFlange->DefineSection(3, z0, kRB26s5RFlangeRi3, kRB26s5RFlangeRo);
  shRB26s5RFlange->DefineSection(4, z0, kRB26s5RFlangeRi2, kRB26s5RFlangeRo);
  z0 += kRB26s5RFlangeD2;
  shRB26s5RFlange->DefineSection(5, z0, kRB26s5RFlangeRi2, kRB26s5RFlangeRo);
  shRB26s5RFlange->DefineSection(6, z0, kRB26s5RFlangeRi1, kRB26s5RFlangeRo);
  z0 += kRB26s5RFlangeD1;
  shRB26s5RFlange->DefineSection(7, z0, kRB26s5RFlangeRi1, kRB26s5RFlangeRo);
  TGeoVolume* voRB26s5RFlange = new TGeoVolume("RB26s5RFlange", shRB26s5RFlange, kMedSteel);

  TGeoVolume* voRB26s5RFlangeM = new TGeoVolume("RB26s5RFlangeM", MakeMotherFromTemplate(shRB26s5RFlange, 4, 7), kMedVac);
  voRB26s5RFlangeM->AddNode(voRB26s5RFlange, 1, gGeoIdentity);

  //      
  // Assemble RB26/1-2
  //
  TGeoVolumeAssembly* asRB26s12 = new TGeoVolumeAssembly("RB26s12"); 
  z0 = 0.;
  asRB26s12->AddNode(voRB26s1RFlange,       1, gGeoIdentity);
  z0 += kRB26s1RFlangeIsL + kRB26s1RFlangeFpL;
  asRB26s12->AddNode(voRB26s12TubeM,         1, new TGeoTranslation(0., 0., z0));
  z0 += kRB26s12TubeL;
  asRB26s12->AddNode(voRB26s2Compensator,   1, new TGeoTranslation(0., 0., z0));
  z0 += kRB26s2CompL;
  z0 -= kRB26s2FFlangeD1;
  asRB26s12->AddNode(voRB26s2FFlangeM,       1, new TGeoTranslation(0., 0., z0));
  z0 += kRB26s2FFlangeL;
  const Float_t kRB26s12L = z0;

  //
  // Assemble RB26/3
  //
  TGeoVolumeAssembly* asRB26s3 = new TGeoVolumeAssembly("RB26s3"); 
  z0 = 0.;
  asRB26s3->AddNode(voRB26s3SFlangeM,      1, gGeoIdentity);
  z0 +=  kRB26s3SFlangeL;
  z0 -=  kRB26s3SFlangeD3;
  asRB26s3->AddNode(voRB26s3FixedPointM,   1, new TGeoTranslation(0., 0., z0));
  z0 += kRB26s3FixedPointL;
  asRB26s3->AddNode(voRB26s3TubeM,         1, new TGeoTranslation(0., 0., z0));
  z0 += kRB26s3TubeL;
  asRB26s3->AddNode(voRB26s3Compensator,   1, new TGeoTranslation(0., 0., z0));
  z0 += kRB26s3CompL;
  z0 -= kRB26s3FFlangeD1;
  asRB26s3->AddNode(voRB26s3FFlangeM,      1, new TGeoTranslation(0., 0., z0));
  z0 += kRB26s3FFlangeL;
  const Float_t kRB26s3L = z0;
      

  //
  // Assemble RB26/4-5
  //
  TGeoVolumeAssembly* asRB26s45 = new TGeoVolumeAssembly("RB26s45"); 
  z0 = 0.;
  asRB26s45->AddNode(voRB26s4SFlangeM,       1, gGeoIdentity);
  z0 +=  kRB26s4SFlangeL;
  z0 -=  kRB26s4SFlangeD3;
  asRB26s45->AddNode(voRB26s4FixedPointM,    1, new TGeoTranslation(0., 0., z0));
  z0 += kRB26s4FixedPointL;
  asRB26s45->AddNode(voRB26s45TubeM,         1, new TGeoTranslation(0., 0., z0));
  z0 += kRB26s45TubeL;
  asRB26s45->AddNode(voRB26s5Compensator,    1, new TGeoTranslation(0., 0., z0));
  z0 += kRB26s5CompL;
  z0 -= kRB26s5RFlangeD3;
  z0 -= kRB26s5RFlangeD4;
  asRB26s45->AddNode(voRB26s5RFlangeM,       1, new TGeoTranslation(0., 0., z0));
  z0 += kRB26s5RFlangeL;
  const Float_t kRB26s45L = z0;
      
  //
  // Assemble RB26
  //
  TGeoVolumeAssembly* asRB26Pipe = new TGeoVolumeAssembly("RB26Pipe"); 
  z0 = 0.;
  asRB26Pipe->AddNode(asRB26s12,       1, new TGeoTranslation(0., 0., z0));
  z0 +=  kRB26s12L;
  asRB26Pipe->AddNode(asRB26s3,        1, new TGeoTranslation(0., 0., z0));
  z0 +=  kRB26s3L;
  asRB26Pipe->AddNode(asRB26s45,       1, new TGeoTranslation(0., 0., z0));
  z0 +=  kRB26s45L;
  top->AddNode(asRB26Pipe, 1, new TGeoCombiTrans(0., 0., -82., rot180));
}



//___________________________________________
void AliPIPEv4::CreateMaterials()
{
  //
  // Define materials for beam pipe
  //

  AliDebugClass(1,"Create PIPEv4 materials");
  Int_t   isxfld = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ();
  Float_t sxmgmx = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max();
  // Steel (Inox)  
  Float_t asteel[4] = { 55.847,51.9961,58.6934,28.0855 };
  Float_t zsteel[4] = { 26.,24.,28.,14. };
  Float_t wsteel[4] = { .715,.18,.1,.005 };
  // AlBe - alloy 
  Float_t aAlBe[2] = { 26.98, 9.01};             // al=2.702 be=1.8477      
  Float_t zAlBe[2] = { 13.00, 4.00};
  Float_t wAlBe[2] = { 0.4, 0.6};
  //
  // Polyamid
  Float_t aPA[4] = {16., 14., 12.,  1.};
  Float_t zPA[4] = { 8.,  7.,  6.,  1.};
  Float_t wPA[4] = { 1.,  1.,  6., 11.};
  //
  // Polyimide film
  Float_t aPI[4] = {16., 14., 12.,  1.};
  Float_t zPI[4] = { 8.,  7.,  6.,  1.};
  Float_t wPI[4] = { 5.,  2.,  22., 10.};
  // Rohacell
  Float_t aRohacell[4] = {16., 14., 12.,  1.};
  Float_t zRohacell[4] = { 8.,  7.,  6.,  1.};
  Float_t wRohacell[4] = { 2.,  1.,  9., 13.};
  //
  // Air 
  //
  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir = 1.20479E-3;
  Float_t dAir1 = 1.20479E-10;
  //
  // Insulation powder
  //                    Si         O       Ti     Al
  Float_t ains[4] ={28.0855, 15.9994, 47.867,  26.982};
  Float_t zins[4] ={14.,      8.    , 22.   ,  13.   };
  Float_t wins[4] ={ 0.3019,  0.4887,  0.1914,  0.018};
  //
  //
  // Anticorodal
  //
  // Al Si7 Mg 0.6
  //
  Float_t aaco[3] ={26.982, 28.0855, 24.035};
  Float_t zaco[3] ={13.,    14.    , 12.   };
  Float_t waco[3] ={ 0.924,  0.07,  0.006};
  // Kapton
  //
  Float_t aKapton[4]={1.00794,12.0107, 14.010,15.9994};
  Float_t zKapton[4]={1.,6.,7.,8.};
  Float_t wKapton[4]={0.026362,0.69113,0.07327,0.209235};
  Float_t dKapton = 1.42;
  // NEG coating
  //                  Ti     V      Zr
  Float_t aNEG[4] = {47.87, 50.94, 91.24};
  Float_t zNEG[4] = {22.00, 23.00, 40.00};
  Float_t wNEG[4] = {1./3., 1./3., 1./3.};  
  Float_t dNEG = 5.6; // ?

  //---------------------------------
  // Aluminium AA 5083 for MFT: Al Manganese(Mn) Magnesium(Mg) Chrome(Cr)
  Float_t aALU5083[4]={26.982, 54.938, 24.305, 51.996};  // Mg pas meme a que la ligne Anticorodal!
  Float_t zALU5083[4] ={13., 25., 12., 24.};
  Float_t wALU5083[4] ={0.947, 0.007, 0.044, 0.0015};
  // Aluminium AA 2219 for MFT: Al Cu Mn Ti V Zr
  Float_t aALU2219[6]={26.982, 63.546, 54.938, 47.867, 50.941, 91.224};
  Float_t zALU2219[6] ={13., 29., 25., 22., 23., 40.};
  Float_t wALU2219[6] ={0.93, 0.063, 0.003, 0.0006, 0.001, 0.0018};
  //---------------------------------

  //
  // Silicon for ITS UPGRADE
  AliMaterial(2,  "SILICON$",28.09 , 14.00 , 2.33 , 9.36 , 45.); 

  //
  //     Berillium 
  AliMaterial(5, "BERILLIUM$", 9.01, 4., 1.848, 35.3, 36.7);
  //
  //     Carbon 
  AliMaterial(6,  "CARBON$   ", 12.01, 6., 2.265, 18.8, 49.9);
  //
  //     Aluminum 
  AliMaterial(9,  "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
  //
  //     Copper 
  AliMaterial(10, "COPPER", 63.55, 29, 8.96, 1.43, 85.6/8.96);
  //
  //     Air 
  AliMixture(15, "AIR$      ", aAir, zAir, dAir, 4, wAir);
  AliMixture(35, "AIR_HIGH$ ", aAir, zAir, dAir, 4, wAir);
  //
  //     Vacuum 
  AliMixture(16, "VACUUM$ ", aAir, zAir, dAir1, 4, wAir);
  //
  //     stainless Steel 
  AliMixture(19, "STAINLESS STEEL$", asteel, zsteel, 7.88, 4, wsteel);
  //
  //     reduced density steel to approximate pump getter material
  AliMixture(20, "GETTER$", asteel, zsteel, 1.00, 4, wsteel);
  //     Al-Be alloy
  //     
  AliMixture(21, "AlBe$", aAlBe, zAlBe, 2.07, 2, wAlBe); 
  //     Polyamid
  //   
  AliMixture(22, "PA$", aPA, zPA, 1.14, -4, wPA);
  //
  //     Kapton
  AliMixture(23, "KAPTON", aKapton, zKapton, dKapton, 4, wKapton);
  // Anticorodal 
  AliMixture(24, "ANTICORODAL", aaco, zaco, 2.66, 3, waco);
  
  //
  //     Insulation powder 
  AliMixture(14, "INSULATION0$", ains, zins, 0.41, 4, wins);
  AliMixture(34, "INSULATION1$", ains, zins, 0.41, 4, wins);
  AliMixture(54, "INSULATION2$", ains, zins, 0.41, 4, wins);

  //    NEG
  AliMixture(25, "NEG COATING", aNEG, zNEG, dNEG, -3, wNEG);
  
  //---------------------------------
  //  Aluminium AA5083 for MFT
  AliMixture(63, "ALUMINIUM5083$",aALU5083,zALU5083, 2.66 ,4,wALU5083); // from aubertduval.fr
  // Aluminium AA2219 for MFT
  AliMixture(64, "ALUMINIUM2219$",aALU2219,zALU2219, 2.84 ,6,wALU2219); // from aubertduval.fr
  //---------------------------------
  //     Polyimide Film
  //
  AliMixture(65, "PI$", aPI, zPI, 1.42, -4, wPI);
  //---------------------------------
  //     Carbon Fiber M55J
  AliMaterial(66,"M55J6K$",12.0107,6,1.92,999,999);
  // Rohacell  C9 H13 N1 O2  0.03 g/cm^3
  AliMixture(67,"Rohacell$", aRohacell, zRohacell, 0.03, -4, wRohacell);

  // ****************
  //     Defines tracking media parameters. 
  //
  Float_t epsil  = .001;    // Tracking precision, 
  Float_t stemax = -0.01;   // Maximum displacement for multiple scat 
  Float_t tmaxfd = -20.;    // Maximum angle due to field deflection 
  Float_t deemax = -.3;     // Maximum fractional energy loss, DLS 
  Float_t stmin  = -.8;
  // *************** 
  //
  // Silicon for ITS UPGRADE
  AliMedium(2,  "SILICON",   2, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);


  //    Beryllium 
  
  AliMedium(5, "BE",       5, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);

  //    Carbon 
  AliMedium(6, "C",        6, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Aluminum 
  AliMedium(9, "ALU",      9, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //    Copper 
  AliMedium(10, "CU",      10, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Air 
  AliMedium(15, "AIR",     15, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(35, "AIR_HIGH",35, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Vacuum 
  AliMedium(16, "VACUUM", 16, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Steel 
  AliMedium(19, "INOX",   19, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //    Getter 
  AliMedium(20, "GETTER", 20, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //   AlBe - Aloy 
  AliMedium(21, "AlBe"  , 21, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //   Polyamid
  AliMedium(22, "PA"  ,   22, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //   Antocorodal
  AliMedium(24, "ANTICORODAL",   24, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //    Insulation Powder 
  AliMedium(14, "INS_C0          ", 14, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(34, "INS_C1          ", 34, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(54, "INS_C2          ", 54, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //
  //   KAPTON
  AliMedium(23, "KAPTON", 23, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);

  //
  //   NEG
  AliMedium(25, "NEG COATING", 25, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);

  //----------------- for the MFT ----------------------
  AliMedium(63,"AA5083", 63, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(64,"AA2219", 64, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //----------------------------------------------------
  AliMedium(65,"POLYIMIDE", 65, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  //---------------------------------
  //     Carbon Fiber M55J
  AliMedium(66,  "M55J6K",66,0,isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
  // Rohacell
  AliMedium(67,"ROHACELL",67,0,isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);


}


TGeoPcon* AliPIPEv4::MakeMotherFromTemplate(TGeoPcon* shape, Int_t imin, Int_t imax, Float_t r0, Int_t nz)
{
  //
  //  Create a mother shape from a template setting some min radii to 0
  //
  Int_t nz0 = shape->GetNz();
  // if nz > -1 the number of planes is given by nz
  if (nz != -1) nz0 = nz;
  TGeoPcon* mother = new TGeoPcon(0., 360., nz0);

  if (imin == -1 || imax == -1) {
    imin = 0;
    imax = shape->GetNz();
  } else if (imax >= nz0) {
    imax = nz0 - 1;
    printf("Warning: imax reset to nz-1 %5d %5d %5d %5d\n", imin, imax, nz, nz0);
  }
    

    
  for (Int_t i = 0;  i < shape->GetNz(); i++) {
    Double_t rmin = shape->GetRmin(i);
    if ((i >= imin) && (i <= imax) ) rmin = r0;
    Double_t rmax = shape->GetRmax(i);
    Double_t    z = shape->GetZ(i);
    mother->DefineSection(i, z, rmin, rmax);
  }
  return mother;
    
}

TGeoPcon* AliPIPEv4::MakeInsulationFromTemplate(TGeoPcon* shape)
{
  //
  //  Create an beam pipe insulation layer shape from a template
  //
  Int_t nz = shape->GetNz();
  TGeoPcon* insu = new TGeoPcon(0., 360., nz);
    
  for (Int_t i = 0;  i < nz; i++) {
    Double_t    z = shape->GetZ(i);
    Double_t rmin = shape->GetRmin(i);
    Double_t rmax = shape->GetRmax(i);
    rmax += 0.5;
    shape->DefineSection(i, z, rmin, rmax);	
    rmin  = rmax - 0.5;
    insu->DefineSection(i, z, rmin, rmax);	
  }
  return insu;
    
}


TGeoVolume* AliPIPEv4::MakeBellow(const char* ext, Int_t nc, Float_t rMin, Float_t rMax, Float_t dU, Float_t rPlie, Float_t dPlie)
{
  // nc     Number of convolution
  // rMin   Inner radius of the bellow
  // rMax   Outer radius of the bellow
  // dU     Undulation length
  // rPlie  Plie radius
  // dPlie  Plie thickness
  const TGeoMedium* kMedVac    =  gGeoManager->GetMedium("PIPE_VACUUM");    
  //const TGeoMedium* kMedSteel  =  gGeoManager->GetMedium("PIPE_INOX");   
  const TGeoMedium* kMedAlu5083 =  gGeoManager->GetMedium("PIPE_AA5083"); //fm       

  char name[64], nameA[64], nameB[64], bools[64];
  snprintf(name, 64, "%sBellowUS", ext);
  TGeoVolume* voBellow = new TGeoVolume(name, new TGeoTube(rMin, rMax, dU/2.), kMedVac);
  //      
  //  Upper part of the undulation
  //
  TGeoTorus* shPlieTorusU  =  new TGeoTorus(rMax - rPlie, rPlie - dPlie, rPlie);
  snprintf(nameA, 64, "%sTorusU", ext);
  shPlieTorusU->SetName(nameA);
  TGeoTube*  shPlieTubeU   =  new TGeoTube (rMax - rPlie, rMax, rPlie);
  snprintf(nameB, 64, "%sTubeU", ext);
  shPlieTubeU->SetName(nameB);
  snprintf(name, 64, "%sUpperPlie", ext);
  snprintf(bools, 64, "%s*%s", nameA, nameB);
  TGeoCompositeShape*  shUpperPlie = new TGeoCompositeShape(name, bools);
    
  TGeoVolume* voWiggleU = new TGeoVolume(name, shUpperPlie, kMedAlu5083);
  voWiggleU->SetLineColor(kOrange); // fm
  
  // Lower part of the undulation
  TGeoTorus* shPlieTorusL =  new TGeoTorus(rMin + rPlie, rPlie - dPlie, rPlie);
  snprintf(nameA, 64, "%sTorusL", ext);
  shPlieTorusL->SetName(nameA);
  TGeoTube*  shPlieTubeL  =  new TGeoTube (rMin, rMin + rPlie, rPlie);
  snprintf(nameB, 64, "%sTubeL", ext);
  shPlieTubeL->SetName(nameB);
  snprintf(name, 64, "%sLowerPlie", ext);
  snprintf(bools, 64, "%s*%s", nameA, nameB);
  TGeoCompositeShape*  shLowerPlie = new TGeoCompositeShape(name, bools);
    
  TGeoVolume* voWiggleL = new TGeoVolume(name, shLowerPlie, kMedAlu5083); 
  voWiggleL->SetLineColor(kOrange); // fm
  
  // Connection between upper and lower part of undulation
  snprintf(name, 64, "%sPlieConn1", ext);
  TGeoVolume* voWiggleC1 = new TGeoVolume(name, new TGeoTube(rMin + rPlie, rMax - rPlie, dPlie/2.), kMedAlu5083);
  voWiggleC1->SetLineColor(kOrange); // fm
  
  // One wiggle
  Float_t dz = rPlie -  dPlie / 2.;
  Float_t z0 = -  dPlie / 2.;
  snprintf(name, 64, "%sWiggle", ext);
  TGeoVolumeAssembly* asWiggle = new TGeoVolumeAssembly(name);

  asWiggle->AddNode(voWiggleC1,  1 , new TGeoTranslation(0., 0., z0));
  z0 += dz;
  asWiggle->AddNode(voWiggleU,   1 , new TGeoTranslation(0., 0., z0));
  z0 += dz;
  asWiggle->AddNode(voWiggleC1,  2 , new TGeoTranslation(0., 0., z0));
  z0 += dz;
  asWiggle->AddNode(voWiggleL ,  1 , new TGeoTranslation(0., 0., z0));
  // Positioning of the volumes
  z0   = - dU / 2.+ rPlie;
  ////////////voBellow->AddNode(voWiggleL, 2, new TGeoTranslation(0., 0., z0));   removing the first 1/2 plie, fm
  z0  +=  rPlie;
  Float_t zsh  = 4. *  rPlie -  2. * dPlie;
  for (Int_t iw = 0; iw < nc; iw++) {
    Float_t zpos =  z0 + iw * zsh;	
    voBellow->AddNode(asWiggle,  iw + 1, new TGeoTranslation(0., 0., zpos - dPlie));

  }
  return voBellow;
}

TGeoVolume* AliPIPEv4::MakeBellowCside(const char* ext, Int_t nc, Float_t rMin, Float_t rMax, Float_t rPlie, Float_t dPlie)
{
  // nc     Number of convolution
  // rMin   Inner radius of the bellow
  // rMax   Outer radius of the bellow
  // dU     Undulation length
  // rPlie  Plie radius
  // dPlie  Plie thickness
  const TGeoMedium* kMedVac    =  gGeoManager->GetMedium("PIPE_VACUUM");
  //const TGeoMedium* kMedSteel  =  gGeoManager->GetMedium("PIPE_INOX");
  const TGeoMedium* kMedAlu5083 =  gGeoManager->GetMedium("PIPE_AA5083"); //fm
  
  Float_t dU = nc * (4.*rPlie - 2. *dPlie);
  
  char name[64], nameA[64], nameB[64], bools[64];
  snprintf(name, 64, "%sBellowUS", ext);
//  TGeoVolume* voBellow = new TGeoVolume(name, new TGeoTube(rMin, rMax, dU/2.), kMedVac);
  TGeoVolumeAssembly *voBellow = new TGeoVolumeAssembly(name);
  //
  //  Upper part of the undulation
  //
  
  TGeoTorus* shPlieTorusU  =  new TGeoTorus(rMax - rPlie, rPlie - dPlie, rPlie);
  snprintf(nameA, 64, "%sTorusU", ext);
  shPlieTorusU->SetName(nameA);
  TGeoTube*  shPlieTubeU   =  new TGeoTube (rMax - rPlie, rMax, rPlie);
  snprintf(nameB, 64, "%sTubeU", ext);
  shPlieTubeU->SetName(nameB);
  snprintf(name, 64, "%sUpperPlie", ext);
  snprintf(bools, 64, "%s*%s", nameA, nameB);
  TGeoCompositeShape*  shUpperPlie = new TGeoCompositeShape(name, bools);
  
  TGeoVolume* voWiggleU = new TGeoVolume(name, shUpperPlie, kMedAlu5083);
  voWiggleU->SetLineColor(kOrange); // fm
  
  // First Lower part of the ondulation
  TGeoTorus* shPlieTorusL =  new TGeoTorus(rMin + rPlie, rPlie - dPlie, rPlie);
  snprintf(nameA, 64, "%sTorusL", ext);
  shPlieTorusL->SetName(nameA);
  TGeoTranslation *t1 = new TGeoTranslation("t1",0,0,-rPlie/2.);
  t1->RegisterYourself();

  TGeoTube*  shPlieTubeL  =  new TGeoTube (rMin, rMin + rPlie, rPlie/2.);
  snprintf(nameB, 64, "%sTubeL", ext);
  shPlieTubeL->SetName(nameB);
  snprintf(name, 64, "%sLowerPlie", ext);
  snprintf(bools, 64, "%s*%s:t1", nameA, nameB);
  TGeoCompositeShape*  shLowerPlie1 = new TGeoCompositeShape(name, bools);
  
  TGeoVolume* voWiggleL1 = new TGeoVolume(name, shLowerPlie1, kMedAlu5083);
  voWiggleL1->SetLineColor(kOrange); // fm
  
  // Second Lower part of the undulation
  TGeoTranslation *t2 = new TGeoTranslation("t2",0,0,rPlie/2.);
  t2->RegisterYourself();
  
  snprintf(bools, 64, "%s*%s:t2", nameA, nameB);
  TGeoCompositeShape*  shLowerPlie2 = new TGeoCompositeShape(name, bools);
  
  TGeoVolume* voWiggleL2 = new TGeoVolume(name, shLowerPlie2, kMedAlu5083);
  voWiggleL2->SetLineColor(kOrange); // fm
  
  // Connection between upper and lower part of undulation
  snprintf(name, 64, "%sPlieConn1", ext);
  TGeoVolume* voWiggleC1 = new TGeoVolume(name, new TGeoTube(rMin + rPlie, rMax - rPlie, dPlie/2.), kMedAlu5083);
  voWiggleC1->SetLineColor(kOrange); // fm

  //
  // Vacuum Part
  //
  
  //--Upper part of the ondulation
  
  TGeoTorus* vacPlieTorusU  =  new TGeoTorus(rMax - rPlie, 0., rPlie- dPlie);
  snprintf(nameA, 64, "%svacTorusU", ext);
  vacPlieTorusU->SetName(nameA);
  TGeoTube*  vacPlieTubeU   =  new TGeoTube (0., rMax- rPlie, rPlie-dPlie);
  snprintf(nameB, 64, "%svacTubeU", ext);
  vacPlieTubeU->SetName(nameB);
  snprintf(name, 64, "%svacUpperPlie", ext);
  snprintf(bools, 64, "%s+%s", nameA, nameB);
  TGeoCompositeShape*  vacUpperPlie = new TGeoCompositeShape(name, bools);
  
  TGeoVolume* voVacWiggleU = new TGeoVolume(name, vacUpperPlie, kMedVac);
  voVacWiggleU->SetVisibility(0);

  
  // First Lower part of the undulation
  TGeoTorus* vacPlieTorusL =  new TGeoTorus(rMin + rPlie, 0., rPlie);
  snprintf(nameA, 64, "%svacTorusL", ext);
  vacPlieTorusL->SetName(nameA);
  
  TGeoTube*  vacPlieTubeL  =  new TGeoTube (0., rMin + rPlie, rPlie/2.);
  snprintf(nameB, 64, "%svacTubeL", ext);
  vacPlieTubeL->SetName(nameB);
  snprintf(name, 64, "%svacLowerPlie", ext);
  snprintf(bools, 64, "%s:t1-%s", nameB, nameA);
  TGeoCompositeShape*  vacLowerPlie1 = new TGeoCompositeShape(name, bools);
  
  TGeoVolume* voVacWiggleL1 = new TGeoVolume(name, vacLowerPlie1, kMedVac);
  voVacWiggleL1->SetVisibility(0);
  
  
  // Second Lower part of the undulation
  
  snprintf(bools, 64, "%s:t2-%s", nameB, nameA);
  TGeoCompositeShape*  vacLowerPlie2 = new TGeoCompositeShape(name, bools);
  
  TGeoVolume* voVacWiggleL2 = new TGeoVolume(name, vacLowerPlie2, kMedVac);
  voVacWiggleL2->SetVisibility(0);

  
  // One wiggle
  Float_t dz = rPlie -  dPlie / 2.;
  Float_t z0 = 2.*rPlie;
  snprintf(name, 64, "%sWiggle", ext);
  TGeoVolumeAssembly* asWiggle = new TGeoVolumeAssembly(name);
  
  asWiggle->AddNode(voWiggleL1 ,  1 , new TGeoTranslation(0., 0., z0));
  asWiggle->AddNode(voVacWiggleL1 ,  1 , new TGeoTranslation(0., 0., z0));
  z0 -= dz;
  asWiggle->AddNode(voWiggleC1,  1 , new TGeoTranslation(0., 0., z0));
  z0 -= dz;
  asWiggle->AddNode(voWiggleU,   1 , new TGeoTranslation(0., 0., z0));
  asWiggle->AddNode(voVacWiggleU,   1 , new TGeoTranslation(0., 0., z0));
  z0 -= dz;
  asWiggle->AddNode(voWiggleC1,  2 , new TGeoTranslation(0., 0., z0));
  z0 -= dz;
  asWiggle->AddNode(voWiggleL2 ,  1 , new TGeoTranslation(0., 0., z0));
  asWiggle->AddNode(voVacWiggleL2 ,  1 , new TGeoTranslation(0., 0., z0));


  
  
  // Positioning of the volumes
  z0   = + dU / 2.;
  Float_t zsh  = 4. *  dz;
  //for (Int_t iw = 0; iw < 1; iw++) {
     for (Int_t iw = 0; iw < nc; iw++) {
    Float_t zpos =  z0 - iw * zsh;
    voBellow->AddNode(asWiggle,  iw + 1, new TGeoTranslation(0., 0., zpos));
    
  }
  return voBellow;
}
//TGeoVolume* AliPIPEv4::MakeBellowCside(const char* ext, Int_t nc, Float_t rMin, Float_t rMax, Float_t dU, Float_t rPlie, Float_t dPlie)
//{
//  // nc     Number of convolution
//  // rMin   Inner radius of the bellow
//  // rMax   Outer radius of the bellow
//  // dU     Undulation length
//  // rPlie  Plie radius
//  // dPlie  Plie thickness
//  const TGeoMedium* kMedVac    =  gGeoManager->GetMedium("PIPE_VACUUM");
//  //const TGeoMedium* kMedSteel  =  gGeoManager->GetMedium("PIPE_INOX");
//  const TGeoMedium* kMedAlu5083 =  gGeoManager->GetMedium("PIPE_AA5083"); //fm
//  
//  char name[64], nameA[64], nameB[64], bools[64];
//  snprintf(name, 64, "%sBellowUS", ext);
//  TGeoVolume* voBellow = new TGeoVolume(name, new TGeoTube(rMin, rMax, dU/2.), kMedVac);
//  //
//  //  Upper part of the undulation
//  //
//  
//  TGeoTorus* shPlieTorusU  =  new TGeoTorus(rMax - rPlie, rPlie - dPlie, rPlie);
//  snprintf(nameA, 64, "%sTorusU", ext);
//  shPlieTorusU->SetName(nameA);
//  TGeoTube*  shPlieTubeU   =  new TGeoTube (rMax - rPlie, rMax, rPlie);
//  snprintf(nameB, 64, "%sTubeU", ext);
//  shPlieTubeU->SetName(nameB);
//  snprintf(name, 64, "%sUpperPlie", ext);
//  snprintf(bools, 64, "%s*%s", nameA, nameB);
//  TGeoCompositeShape*  shUpperPlie = new TGeoCompositeShape(name, bools);
//  
//  TGeoVolume* voWiggleU = new TGeoVolume(name, shUpperPlie, kMedAlu5083);
//  voWiggleU->SetLineColor(kOrange); // fm
//  
//  // First Lower part of the undulation
//  TGeoTorus* shPlieTorusL =  new TGeoTorus(rMin + rPlie, rPlie - dPlie, rPlie);
//  snprintf(nameA, 64, "%sTorusL", ext);
//  shPlieTorusL->SetName(nameA);
//  TGeoTranslation *t1 = new TGeoTranslation("t1",0,0,-rPlie/2.);
//  t1->RegisterYourself();
//  
//  TGeoTube*  shPlieTubeL  =  new TGeoTube (rMin, rMin + rPlie, rPlie/2.);
//  snprintf(nameB, 64, "%sTubeL", ext);
//  shPlieTubeL->SetName(nameB);
//  snprintf(name, 64, "%sLowerPlie", ext);
//  snprintf(bools, 64, "%s*%s:t1", nameA, nameB);
//  TGeoCompositeShape*  shLowerPlie1 = new TGeoCompositeShape(name, bools);
//  
//  TGeoVolume* voWiggleL1 = new TGeoVolume(name, shLowerPlie1, kMedAlu5083);
//  voWiggleL1->SetLineColor(kOrange); // fm
//  
//  // Second Lower part of the undulation
//  TGeoTranslation *t2 = new TGeoTranslation("t2",0,0,rPlie/2.);
//  t2->RegisterYourself();
//  
//  snprintf(bools, 64, "%s*%s:t2", nameA, nameB);
//  TGeoCompositeShape*  shLowerPlie2 = new TGeoCompositeShape(name, bools);
//  
//  TGeoVolume* voWiggleL2 = new TGeoVolume(name, shLowerPlie2, kMedAlu5083);
//  voWiggleL2->SetLineColor(kOrange); // fm
//  
//  // Connection between upper and lower part of undulation
//  snprintf(name, 64, "%sPlieConn1", ext);
//  TGeoVolume* voWiggleC1 = new TGeoVolume(name, new TGeoTube(rMin + rPlie, rMax - rPlie, dPlie/2.), kMedAlu5083);
//  voWiggleC1->SetLineColor(kOrange); // fm
//  
//  // One wiggle
//  Float_t dz = rPlie -  dPlie / 2.;
//  Float_t z0 = 2.*rPlie;
//  snprintf(name, 64, "%sWiggle", ext);
//  TGeoVolumeAssembly* asWiggle = new TGeoVolumeAssembly(name);
//  
//  asWiggle->AddNode(voWiggleL1 ,  1 , new TGeoTranslation(0., 0., z0));
//  //  z0 -= dz;
//  //  asWiggle->AddNode(voWiggleC1,  1 , new TGeoTranslation(0., 0., z0));
//  //  z0 -= dz;
//  //  asWiggle->AddNode(voWiggleU,   1 , new TGeoTranslation(0., 0., z0));
//  //  z0 -= dz;
//  //  asWiggle->AddNode(voWiggleC1,  2 , new TGeoTranslation(0., 0., z0));
//  //  z0 -= dz;
//  //  asWiggle->AddNode(voWiggleL2 ,  1 , new TGeoTranslation(0., 0., z0));
//  // Positioning of the volumes
//  z0   = + dU / 2.;
//  Float_t zsh  = 4. *  dz;
//  for (Int_t iw = 0; iw < 1; iw++) {
//    // for (Int_t iw = 0; iw < nc; iw++) {
//    Float_t zpos =  z0 - iw * zsh;
//    voBellow->AddNode(asWiggle,  iw + 1, new TGeoTranslation(0., 0., zpos));
//    
//  }
//  return voBellow;
//}


