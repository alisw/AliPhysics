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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Muon ABSOrber                                                            //
//  This class contains the description of the muon absorber geometry        //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliABSOClass.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>The responsible person for this module is
<a href="mailto:andreas.morsch@cern.ch">Andreas Morsch</a>.
</font>
<pre>
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TVirtualMC.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoMatrix.h>
#include <TGeoCompositeShape.h>
#include <TGeoBBox.h>
#include <TGeoXtru.h>
#include <TGeoTube.h>
#include <TGeoPgon.h>
#include <TGeoArb8.h>
#include <TGeoMedium.h>

#include "AliABSOv0.h"
#include "AliConst.h"
#include "AliRun.h"
#include "AliLog.h"

ClassImp(AliABSOv0)
 
//_____________________________________________________________________________
AliABSOv0::AliABSOv0()
{
  //
  // Default constructor
  //
}
 
//_____________________________________________________________________________
AliABSOv0::AliABSOv0(const char *name, const char *title)
       : AliABSO(name,title)
{
  //
  // Standard constructor
  //
}
 
//_____________________________________________________________________________
void AliABSOv0::CreateGeometry()
{
    //
    // Creation of the geometry of the muon absorber
    //
    //Begin_Html
    /*
      <img src="picts/AliABSOv0Tree.gif">
    */
    //End_Html
    //Begin_Html
    /*
      <img src="picts/AliABSOv0.gif">
    */
    //End_Html
    
    //
    //

    enum {kC=1605, kAl=1608, kFe=1609, kCu=1610, kW=1611, kPb=1612,
	  kNiCuW=1620, kVacuum=1615, kAir=1614, kConcrete=1616,
	  kPolyCH2=1617, kSteel=1618, kInsulation=1613, kPolyCc=1619};	  
    
    Int_t *idtmed = fIdtmed->GetArray()-1599;
    
    Float_t par[24], cpar[5], cpar0[5], pcpar[12], tpar[3], tpar0[3]; 
    Float_t dz;
#include "ABSOSHILConst.h"
#include "ABSOConst.h"
//
// Structure of Tracking Region
//
  Float_t dzFe = 25.;

// 3 < theta < 9
    fNLayers[0] = 5; 
    fMLayers[0][0]  = kAir;              fZLayers[0][0] = kZAbsStart;
    fMLayers[0][1]  = kC;                fZLayers[0][1] = kZAbsCc;             
    fMLayers[0][2]  = kConcrete;         fZLayers[0][2] = kZRear - kDRear - dzFe;
    fMLayers[0][3]  = kSteel;            fZLayers[0][3] = kZRear - kDRear;
    fMLayers[0][4]  = kSteel;            fZLayers[0][4] = kZRear;
// 2 < theta < 3
    fNLayers[1] = 6; 

    fMLayers[1][0] = kAir          ;      fZLayers[1][0] = fZLayers[0][0] - 10.;
    fMLayers[1][1] = kAl           ;      fZLayers[1][1] = fZLayers[0][0];
    fMLayers[1][2] = fMLayers[0][1];      fZLayers[1][2] = fZLayers[0][1];
    fMLayers[1][3] = fMLayers[0][2];      fZLayers[1][3] = fZLayers[0][2];
    fMLayers[1][4] = fMLayers[0][3];      fZLayers[1][4] = fZLayers[0][3];
    fMLayers[1][5] = kNiCuW;              fZLayers[1][5] = fZLayers[0][4];
//    

    Float_t dTube = 0.1;                     // tube thickness
    Float_t dInsu = 0.5;                     // insulation thickness
    Float_t dEnve = 0.1;                     // protective envelope thickness


// Mother volume and outer shielding: Pb
  par[0]  = 0.;
  par[1]  = 360.;
  par[2]  = 7.;
    
  par[21] = (kZRear - kZAbsStart) / 2.;
  par[22] = kRAbs;
  par[23] = kZAbsStart * TMath::Tan(kTheta1);

  par[18] = par[21] - (kZNose - kZAbsStart);
  par[19] = kRAbs;
  par[20] = kZNose * TMath::Tan(kTheta1);

  par[15] = par[21] - (kZConeTPC - kZAbsStart);
  par[16] = kRAbs;
  par[17] = par[20] - (par[15] - par[18]) * TMath::Tan(kTheta2);

  par[12] = par[21]  - (kZOpen - kZAbsStart);
  par[13] = kRAbs;
  par[14] = par[17] - (par[12] - par[15]) * TMath::Tan(kAccMax);

  par[9]  = par[21]  - (kZRear - kDRear - kZAbsStart);
  par[10] = kRAbs   - (par[9] - par[12]) * TMath::Tan(kThetaOpen1) ;
  par[11] = par[14] - (par[9] - par[12]) * TMath::Tan(kAccMax);

  par[6]  = par[21]  - (kZRear - kDRear - kZAbsStart);
  par[7]  = (kZRear - kDRear) * TMath::Tan(kAccMin);
  par[8]  = par[14] - (par[6] - par[12]) * TMath::Tan(kAccMax);

  par[3] =  - par[21];
  par[4] = kZRear  * TMath::Tan(kAccMin);
  par[5] = par[8] - (par[3] - par[6]) * TMath::Tan(kAccMax);
  gMC->Gsvolu("ABSS", "PCON", idtmed[kPb+40], par, 24);

  for (Int_t i = 22; i > 7; i -= 3) par[i]  = 0;

  gMC->Gsvolu("ABSM", "PCON", idtmed[kVacuum+40], par, 24);
  gMC->Gspos("ABSS", 1, "ABSM", 0., 0., 0., 0, "ONLY");

//
// Steel envelope
//
  par[4] = par[5]  - kDSteel;
  par[7] = par[8]  - kDSteel;
  par[10]= par[11] - kDSteel;  
  par[13]= par[14] - kDSteel;  
  par[16]= par[17] - kDSteel;  
  par[19]= par[20] - kDSteel;  
  par[22]= par[23] - kDSteel;  

  gMC->Gsvolu("ABST", "PCON", idtmed[kSteel], par, 24);
  gMC->Gspos("ABST", 1, "ABSS", 0., 0., 0., 0, "ONLY");
//
// Polyethylene shield
// 
  cpar[0] = (kZRear - kZConeTPC) / 2.;
  cpar[1] = kZRear * TMath::Tan(kAccMax);
  cpar[2] = cpar[1] + kDPoly;
  cpar[3] = kZConeTPC * TMath::Tan(kAccMax);
  cpar[4] = cpar[3] + kDPoly;

  gMC->Gsvolu("APOL", "CONE", idtmed[kPolyCH2+40], cpar, 5);
  dz = - (kZRear - kZAbsStart) / 2. + cpar[0];
  gMC->Gspos("APOL", 1, "ABSS", 0., 0., dz, 0, "ONLY");

//
// Tungsten nose to protect TPC
// 
  cpar[0] = (kZNose - kZAbsStart) / 2.;
  cpar[1] = kZNose * TMath::Tan(kAccMax);
  cpar[2] = kZNose * TMath::Tan(kTheta1) - kDSteel;
  cpar[3] = kZAbsStart * TMath::Tan(kAccMax);
  cpar[4] = kZAbsStart * TMath::Tan(kTheta1) - kDSteel;

  gMC->Gsvolu("ANOS", "CONE", idtmed[kNiCuW], cpar, 5);
  //
  dz =  (kZRear - kZAbsStart) / 2. - cpar[0];
  gMC->Gspos("ANOS", 1, "ABSS", 0., 0., dz, 0, "ONLY");
  //
  // Tungsten inner shield
  //
  Float_t zW = kZTwoDeg + .1;
  Float_t dZ = zW + (kZRear - kDRear - zW) / 2.;
  //
  pcpar[0]  = 0.;
  pcpar[1]  = 360.;
  pcpar[2]  = 3.;
  pcpar[9]  = - (zW - dZ);
  pcpar[10] = kRAbs;
  pcpar[11] = zW * TMath::Tan(kAccMin);
  pcpar[6]  = - (kZOpen - dZ);
  pcpar[7]  = kRAbs;
  pcpar[8]  = kZOpen * TMath::Tan(kAccMin);
  pcpar[3]  = - (kZRear - kDRear - dZ);
  pcpar[4]  = kRAbs + (kZRear - kDRear - kZOpen) * TMath::Tan(kThetaOpen1);
  pcpar[5]  = (kZRear - kDRear) * TMath::Tan(kAccMin);
  
  gMC->Gsvolu("AWIN", "PCON", idtmed[kNiCuW+40], pcpar, 12);
  dz = -(zW + kZRear - kDRear) / 2 + (kZAbsStart + kZRear) / 2.;
  gMC->Gspos("AWIN", 1, "ABSS", 0., 0., dz, 0, "ONLY");
//
// First part replaced by Carbon  
//
  cpar[0] = (200.-zW)/2.;

  cpar[1] = kRAbs;
  cpar[2] = 200. * TMath::Tan(kAccMin);
  cpar[3] = kRAbs;
  cpar[4] = pcpar[11];

  gMC->Gsvolu("ACNO", "CONE", idtmed[kC], cpar, 5);
  dz = - (zW - dZ+cpar[0]);
  gMC->Gspos("ACNO", 1, "AWIN", 0., 0., dz, 0, "ONLY");

/*  
  Float_t zWW = 383.5;
  cpar[0] = (kZRear-kDRear-zWW)/2.;
  cpar[1] = kRAbs + (zWW-kZOpen) *  TMath::Tan(kThetaOpen1);
  cpar[2] =  zWW * TMath::Tan(kAccMin);
  cpar[3] = pcpar[10];
  cpar[4] = pcpar[11];
  gMC->Gsvolu("AWNO", "CONE", idtmed[kCu+40], cpar, 5);
  dz = zWW-dZ+cpar[0];
  
  gMC->Gspos("AWNO", 1, "AWIN", 0., 0., dz, 0, "ONLY");
*/
  //
  //     Inner tracking region
  //
  //
  //
  pcpar[0]  = 0.;
  pcpar[1]  = 360.;
  pcpar[2]  = 3.;
  pcpar[9]  = (kZRear - kZAbsStart) / 2.;
  pcpar[10] = kRAbs;
  pcpar[11] = kZAbsStart * TMath::Tan(kAccMax);
  pcpar[6]  = pcpar[9] - (kZTwoDeg - kZAbsStart);
  pcpar[7]  = kRAbs;
  pcpar[8]  = kZTwoDeg * TMath::Tan(kAccMax);
  pcpar[3]  = - pcpar[9];
  pcpar[4]  = kZRear * TMath::Tan(kAccMin);
  pcpar[5]  = kZRear * TMath::Tan(kAccMax);
  gMC->Gsvolu("AITR", "PCON", idtmed[fMLayers[0][4]], pcpar, 12);
  //
  // special Pb medium for last 5 cm of Pb
  Float_t zr = kZRear - 2. - 0.001;
  cpar[0] = 1.0;
  cpar[3] = zr * TMath::Tan(kThetaR);
  cpar[4] = zr * TMath::Tan(kAccMax);
  cpar[1] = cpar[3] + TMath::Tan(kThetaR) * 2;
  cpar[2] = cpar[4] + TMath::Tan(kAccMax) * 2;
  
  gMC->Gsvolu("ARPB", "CONE", idtmed[fMLayers[0][4]], cpar, 5);
  dz= - (kZRear - kZAbsStart) / 2. + cpar[0] - 0.001;
  gMC->Gspos("ARPB", 1, "AITR", 0., 0., dz, 0, "ONLY");
  //
  //     concrete cone: concrete 
  //
  pcpar[3]  = pcpar[9] - (kZRear - kDRear - kZAbsStart);
  pcpar[4] = (kZRear-kDRear) * TMath::Tan(kAccMin);
  pcpar[5] = (kZRear-kDRear) * TMath::Tan(kAccMax);
  gMC->Gsvolu("ACON", "PCON", idtmed[fMLayers[0][2]+40], pcpar, 12);
  gMC->Gspos("ACON", 1, "AITR", 0., 0., 0., 0, "ONLY");
//
//    Fe Cone 
//
  zr = kZRear - kDRear - dzFe;

  cpar[0] = dzFe/2.;
  cpar[3] = zr * TMath::Tan(kAccMin);
  cpar[4] = zr * TMath::Tan(kAccMax);
  cpar[1] = cpar[3] + TMath::Tan(kAccMin) * dzFe;
  cpar[2] = cpar[4] + TMath::Tan(kAccMax) * dzFe;

  gMC->Gsvolu("ACFE", "CONE",idtmed[fMLayers[0][3]], cpar, 5);

  dz = - (kZRear - kZAbsStart) / 2. + kDRear + dzFe / 2.;

  gMC->Gspos("ACFE", 1, "ACON", 0., 0., dz, 0, "ONLY");

  
  //
  //
  //     carbon cone: carbon
  //
  pcpar[3]   = pcpar[9] - (kZAbsCc - kZAbsStart);
  pcpar[4]   = kZAbsCc * TMath::Tan(kAccMin);
  pcpar[5]   = kZAbsCc * TMath::Tan(kAccMax);
  gMC->Gsvolu("ACAR", "PCON", idtmed[fMLayers[0][1]+40], pcpar, 12);
  gMC->Gspos("ACAR", 1, "ACON", 0., 0., 0., 0, "ONLY");
 //
 //     carbon cone outer region
 //
  cpar[0]  = 10.;
  cpar[3]  = kRAbs;
  cpar[4]  = kZAbsStart * TMath::Tan(kAccMax);
  cpar[1]  = kRAbs;
  cpar[2]  = cpar[4] + 2. * cpar[0] * TMath::Tan(kAccMax);

  gMC->Gsvolu("ACAO", "CONE", idtmed[fMLayers[0][1]], cpar, 5);
  dz= (kZRear-kZAbsStart) / 2. - cpar[0];
  gMC->Gspos("ACAO", 1, "ACAR", 0., 0., dz, 0, "ONLY");
  //
  //     inner W shield
  Float_t epsi  = 0.;
  Float_t repsi = 1.;
  
  zr = kZRear - (kDRear - epsi);
  cpar[0] = (kDRear - epsi) / 2.;
  cpar[3] = zr * TMath::Tan(kAccMin);
  cpar[4] = zr * TMath::Tan(kThetaR * repsi);
  cpar[1] = cpar[3] + TMath::Tan(kAccMin) * (kDRear - epsi);
  cpar[2] = cpar[4] + TMath::Tan(kThetaR * repsi) * (kDRear - epsi);

  gMC->Gsvolu("ARW0", "CONE", idtmed[fMLayers[1][5]+40], cpar, 5);
  dz= - (kZRear - kZAbsStart) / 2. + cpar[0];
  gMC->Gspos("ARW0", 1, "AITR", 0., 0., dz, 0, "ONLY");
  //
  // special W medium for last 5 cm of W
  zr = kZRear - 5;
  cpar[0] = 2.5;
  cpar[3] = zr * TMath::Tan(kAccMin);
  cpar[4] = zr * TMath::Tan(kThetaR * repsi);
  cpar[1] = cpar[3] + TMath::Tan(kAccMin) * 5.;
  cpar[2] = cpar[4] + TMath::Tan(kThetaR*repsi) * 5.;

  gMC->Gsvolu("ARW1", "CONE", idtmed[fMLayers[1][5]+20], cpar, 5);
  dz = - (kDRear-epsi) / 2. + cpar[0];
  gMC->Gspos("ARW1", 1, "ARW0", 0., 0., dz, 0, "ONLY");
  //
  // Cu
  Float_t drMin = TMath::Tan(kThetaR) * 5;
  Float_t drMax = TMath::Tan(kAccMax) * 5;
  gMC->Gsvolu("ARPE", "CONE", idtmed[fMLayers[0][4]], cpar, 0);
  cpar[0] = 2.5;

  for (Int_t i = 0; i < 3; i++) {
      zr = kZRear - kDRear + 5 + i * 10.;
      cpar[3] = zr * TMath::Tan(kThetaR);
      cpar[4] = zr * TMath::Tan(kAccMax);
      cpar[1] = cpar[3] + drMin;
      cpar[2] = cpar[4] + drMax;
      dz = - (kZRear - kZAbsStart) / 2. + cpar[0] + 5. + (2 - i)*10;
      gMC->Gsposp("ARPE", i+1, "AITR", 0., 0., dz, 0, "ONLY",cpar,5);
  }

  gMC->Gspos("AITR", 1, "ABSS", 0., 0., 0., 0, "ONLY");	
  dz = - (kZRear - kZAbsStart) / 2. - kZAbsStart;
  gMC->Gspos("ABSM", 1, "ALIC", 0., 0., dz, 0, "ONLY");	
//
//
// vacuum system
//
// pipe and heating jackets
//
//
// cylindrical piece
  tpar0[2] = (kZOpen-kZAbsStart)/2;
  tpar0[0] = kRVacu;
  tpar0[1] = kRVacu + dTube + dInsu + dEnve;
  gMC->Gsvolu("AV11", "TUBE", idtmed[kSteel+40], tpar0, 3);
//
// insulation

  tpar[2] = tpar0[2];
  tpar[0] = kRVacu  + dTube;
  tpar[1] = tpar[0] + dInsu;
  gMC->Gsvolu("AI11", "TUBE", idtmed[kInsulation+40], tpar, 3);
  gMC->Gspos("AI11", 1, "AV11", 0., 0., 0., 0, "ONLY"); 
//
  dz = (kZRear - kZAbsStart) / 2. - tpar0[2];
  gMC->Gspos("AV11", 1, "ABSM", 0., 0., dz, 0, "ONLY"); 
//
// conical piece

  cpar0[0] = (kZRear - kDRear - kZOpen) / 2.;
  cpar0[3] = kRVacu - 0.05;
  cpar0[4] = kRVacu + dTube + dInsu + dEnve;
  Float_t dR = 2. * cpar0[0] * TMath::Tan(kThetaOpen1);
  cpar0[1]=cpar0[3] + dR;
  cpar0[2]=cpar0[4] + dR;
  gMC->Gsvolu("AV21", "CONE", idtmed[kSteel+40], cpar0, 5);
  dTube += 0.05;

//
// insulation
  cpar[0] = cpar0[0];
  cpar[1] = cpar0[1] + dTube;
  cpar[2] = cpar0[1] + dTube + dInsu;
  cpar[3] = cpar0[3] + dTube;
  cpar[4] = cpar0[3] + dTube + dInsu;

  gMC->Gsvolu("AI21", "CONE", idtmed[kInsulation+40], cpar, 5);
  gMC->Gspos("AI21", 1, "AV21", 0., 0., 0., 0, "ONLY"); 
  
  dz = - (kZRear - kZAbsStart) / 2. + cpar0[0] + kDRear;
  gMC->Gspos("AV21", 1, "ABSM", 0., 0., dz, 0, "ONLY"); 
////////////////////////////////////////////////////
//                                                // 
//    Front Absorber Support Structure FASS       // 
//                                                //
//    Drawing ALIP2A__0035                        //
//    Drawing ALIP2A__0089                        //
//    Drawing ALIP2A__0090                        //
//    Drawing ALIP2A__0109                        //
////////////////////////////////////////////////////
  TGeoTranslation* vec0 = new TGeoTranslation(0., 0., 0.);
  
  TGeoVolumeAssembly* voFass   = new TGeoVolumeAssembly("Fass");
  const Float_t kDegRad        = TMath::Pi() / 180.;
  const TGeoMedium* kMedSteel  = gGeoManager->GetMedium("ABSO_ST_C0");
  const TGeoMedium* kMedAlu    = gGeoManager->GetMedium("ABSO_ALU_C0");
  
  const Float_t kFassUBFlangeH = 380.;
  const Float_t kFassUBFlangeW =  77.;
  
  const Float_t kFassUMFlangeH = 380.;
  const Float_t kFassUMFlangeB = 246.;
  const Float_t kFassUMFlangeT =  10.;
  const Float_t kFassUMFalpha  = - TMath::ATan((kFassUMFlangeB-kFassUMFlangeT)/ kFassUMFlangeH / 2.) / kDegRad;      
// Upper back   flange
// B1
// 380 x 77
  TGeoVolume* voFassUBFlange = new TGeoVolume("FassUBFlange", new TGeoBBox(kFassUBFlangeW/2., 
									   kFassUBFlangeH/2., 3./2.), kMedSteel);
  voFass->AddNode(voFassUBFlange, 1, new TGeoTranslation(+1.5 + kFassUBFlangeW/2., 
							 180. + kFassUBFlangeH/2.,
							 kFassUMFlangeB - 1.5));
  voFass->AddNode(voFassUBFlange, 2, new TGeoTranslation(-1.5 - kFassUBFlangeW/2., 
							 180. + kFassUBFlangeH/2.,
							 kFassUMFlangeB - 1.5));
  
  
// Lower back   flange
// Upper median flange
//    Drawing ALIP2A__0090                        //
//    Drawing ALIP2A__0089                        //
//    A2
  
  TGeoVolume* voFassUMFlange = new TGeoVolume("FassUMFlange", 
					      new TGeoTrap(kFassUMFlangeH/2., kFassUMFalpha,  
							   0., 1.5, 
							   kFassUMFlangeB/2., kFassUMFlangeB/2.,
							   0., 1.5, 
							   kFassUMFlangeT/2., kFassUMFlangeT/2.,
							   0.), kMedSteel);
  
  TGeoRotation* rotFass1 = new TGeoRotation("rotFass1", 180., 0., 90., 0., 90., 90.);
  voFass->AddNode(voFassUMFlange,1 , 
		  new TGeoCombiTrans(0., 180. + kFassUMFlangeH/2., -(kFassUMFlangeB+kFassUMFlangeT)/4. + kFassUMFlangeB, 
				     rotFass1));
  
  
// Lower median flange
//    Drawing ALIP2A__0090                        //
//    Drawing ALIP2A__0089                        //
//    A1
  const Float_t kFassLMFlangeH = 242.;
  const Float_t kFassLMFlangeB = 246.;
  const Float_t kFassLMFlangeT =  43.;
  const Float_t kFassLMFalpha  = - TMath::ATan((kFassLMFlangeB-kFassLMFlangeT)/ kFassLMFlangeH / 2.) / kDegRad;
  TGeoVolume* voFassLMFlange = new TGeoVolume("FassLMFlange", 
					      new TGeoTrap(kFassLMFlangeH/2., kFassLMFalpha,  
							   0., 1.5, 
							   kFassLMFlangeB/2., kFassLMFlangeB/2.,
							   0., 1.5, 
							   kFassLMFlangeT/2., kFassLMFlangeT/2.,
							   0.), kMedSteel);
  TGeoRotation* rotFass2 = new TGeoRotation("rotFass2", 180., 0., 90., 0., 90., 270.);
  voFass->AddNode(voFassLMFlange, 1, 
		  new TGeoCombiTrans(0., -180. - kFassLMFlangeH/2., -(kFassLMFlangeB+kFassLMFlangeT)/4. + kFassLMFlangeB, 
				     rotFass2));
  
// Stiffeners
// Support Plate
//
// Central cone
  TGeoPgon* shFassCone = new TGeoPgon(22.5, 360., 8, 4);
  shFassCone->DefineSection(0,   0.,   0., 180.);
  shFassCone->DefineSection(1,   3.,   0., 180.);
  shFassCone->DefineSection(2,   3., 177., 180.);
  shFassCone->DefineSection(3, 246., 177., 180.);
  shFassCone->SetName("FassCone");
  
  TGeoBBox* shFassWindow = new TGeoBBox( 190., 53., 28.);
  shFassWindow->SetName("FassWindow");
  TGeoTranslation* tFassWindow = new TGeoTranslation("tFassWindow", 0., 0., 78.);
  tFassWindow->RegisterYourself();
  
  TGeoTube* shFassApperture = new TGeoTube(0., 104., 3.01);
  shFassApperture->SetName("FassApperture");
  
  TGeoCompositeShape* shFassCentral = 
      new TGeoCompositeShape("shFassCentral", "FassCone-(FassWindow:tFassWindow+FassApperture)");
  
  TGeoVolume* voFassCentral = new TGeoVolume("FassCentral", shFassCentral, kMedSteel);
  voFass->AddNode(voFassCentral, 1, vec0);
  
//
// Aluminum ring
//
  TGeoVolume* voFassAlRing = new TGeoVolume("FassAlRing", new TGeoTube(100., 180., 10.), kMedAlu);
  voFass->AddNode(voFassAlRing, 1, new TGeoTranslation(0., 0., -11.));
  TGeoRotation* rotxz = new TGeoRotation("rotxz",  90., 0., 90., 90., 180., 0.);
  gGeoManager->GetVolume("ALIC")->AddNode(voFass, 1, new TGeoCombiTrans(0., 0., -388.45 - 90., rotxz));
}

//_____________________________________________________________________________

void AliABSOv0::Init()
{
  //
  // Initialisation of the muon absorber after it has been built
  Int_t i;
  //
  if(AliLog::GetGlobalDebugLevel()>0) {
    printf("\n%s: ",ClassName());
    for(i=0;i<35;i++) printf("*");
    printf(" ABSOv0_INIT ");
    for(i=0;i<35;i++) printf("*");
    printf("\n%s: ",ClassName());
    //
    for(i=0;i<80;i++) printf("*");
    printf("\n");
  }
}
 









