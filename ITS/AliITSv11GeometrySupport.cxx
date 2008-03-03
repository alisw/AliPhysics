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

// This class Defines the Geometry for the ITS services and support cones
// outside of the ceneteral volume (except for the Ceneteral support 
// cylinders. Other classes define the rest of the ITS. Specificaly the ITS
// The SSD support cone,SSD Support centeral cylinder, SDD support cone,
// The SDD cupport centeral cylinder, the SPD Thermal Sheald, The supports
// and cable trays on both the RB26 (muon dump) and RB24 sides, and all of
// the cabling from the ladders/stave ends out past the TPC. 

/* $Id$ */
// General Root includes
#include <TMath.h>
// Root Geometry includes
//#include <AliLog.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoPcon.h>
#include <TGeoCone.h>
#include <TGeoTube.h> // contaings TGeoTubeSeg
#include <TGeoArb8.h>
#include <TGeoXtru.h>
#include <TGeoCompositeShape.h>
#include <TGeoMatrix.h>
#include "AliITSv11GeometrySupport.h"

ClassImp(AliITSv11GeometrySupport)

#define SQ(A) (A)*(A)

//______________________________________________________________________
void AliITSv11GeometrySupport::SPDCone(TGeoVolume *moth,TGeoManager *mgr)
{
//
// Creates the SPD thermal shield as a volume assembly
// and adds it to the mother volume
// (this is actually a merge of the previous SPDThermalSheald method
// of AliITSv11GeometrySupport.cxx,v 1.9 2007/06/06 and the
// CreateSPDThermalShield method of AliITSv11Hybrid)
//
// Input:
//         moth : the TGeoVolume owing the volume structure
//         mgr  : the GeoManager (default gGeoManager)
// Output:
//
// Created:         ???          ???
// Updated:      11 Dec 2007  Mario Sitta
//
// Technical data are taken from:  ALICE-Thermal Screen "Cone transition"
// (thermal-screen1_a3.ps), "Cylinder" (thermal-screen2_a3.ps), "Half
// assembly" (thermal-screen3_a3.ps), "Flange" (thermal-screen4_a3.ps)


  // Dimensions of the Central shield
  const Double_t kHalfLengthCentral  = 400.*fgkmm;
  const Double_t kThicknessCentral   = 0.4*fgkmm;
  const Double_t kInnerRadiusCentral = 8.1475*fgkcm;
  const Double_t kOuterRadiusCentral = 9.9255*fgkcm;
  const Double_t kInnerACentral = 3.1674*fgkcm;
  const Double_t kInnerBCentral = 2.023 *fgkcm;
  const Double_t kOuterACentral = 2.4374*fgkcm;
  const Double_t kOuterBCentral = 3.8162*fgkcm;
  // Dimensions of the EndCap shield
  const Double_t kHalfLengthEndCap  = 25.*fgkmm;
  const Double_t kThicknessEndCap   = 2.0*fgkmm;
  const Double_t kInnerRadiusEndCap = 8.0775*fgkcm;
  const Double_t kOuterRadiusEndCap = 9.9955*fgkcm;
  const Double_t kInnerAEndCap = 3.1453*fgkcm;
  const Double_t kInnerBEndCap = 2.0009*fgkcm;
  const Double_t kOuterAEndCap = 2.4596*fgkcm;
  const Double_t kOuterBEndCap = 3.8384*fgkcm;
  // Dimensions of the Cone shield
  const Double_t kHalfLengthCone  = 145.*fgkmm;
  const Double_t kThicknessCone   = 0.3*fgkmm;
  const Double_t kInnerRadialCone = 37.3*fgkcm;
  const Double_t kOuterRadialCone = 39.0*fgkcm;
  const Double_t kInnerACone = 14.2344*fgkcm;
  //  const Double_t kInnerBCone =  9.0915*fgkcm;
  const Double_t kOuterACone =  9.5058*fgkcm;
  //  const Double_t kOuterBCone = 14.8831*fgkcm;
  // Dimensions of the Flange's Ring and Wing
  const Double_t kHalfLengthRing  = 7.5*fgkmm;
  const Double_t kThicknessRing   = 0.3*fgkmm;
  const Double_t kInnerRadiusRing = 37.3*fgkcm;
  const Double_t kOuterRadiusRing = 42.0*fgkcm;
  const Double_t kOuterRadiusWing = 49.25*fgkcm;
  const Double_t kWideWing  = 6.0*fgkcm;
  const Double_t kThetaWing = 45.0;
  // Common data
  const Double_t kTheta = 36.0*TMath::DegToRad();
  const Double_t kThicknessOmega = 0.3*fgkmm;

  // Local variables
  Double_t x, y;
  Double_t xshld[24], yshld[24];
  Double_t xair[24] , yair[24];
  Double_t xomega[48], yomega[48];
  //  Double_t *xyarb8;

  // The entire shield is made up of two half central shields
  // symmetric with respect to the XZ plane, four half end cap
  // shields, again symmetric with respect to the XZ plane, and four
  // half cones, symmetric with respect to the XZ plane too.

  TGeoVolumeAssembly *vM = new TGeoVolumeAssembly("ITSspdThermalShield");

  // The central half shield: a half tube of carbon fiber,
  // a similar but proportionally smaller half tube of air inside it,
  // and a Omega-shaped carbon fiber insert inside the air.
  // They are all XTru shapes

  TGeoXtru *centralshape = new TGeoXtru(2);

  CreateSPDThermalShape(kInnerACentral,kInnerBCentral,kInnerRadiusCentral,
			kOuterACentral,kOuterBCentral,kOuterRadiusCentral,
			kTheta,xshld,yshld);

  centralshape->DefinePolygon(24,xshld,yshld);
  centralshape->DefineSection(0,-kHalfLengthCentral);
  centralshape->DefineSection(1, kHalfLengthCentral);

  // Now rescale to get the air volume dimensions
    InsidePoint(xshld[23], yshld[23],
		xshld[ 0], yshld[ 0],
		xshld[ 1], yshld[ 1], kThicknessCentral,
		xair[0], yair[0]);
  for (Int_t i=1; i<23; i++) {
    InsidePoint(xshld[i-1], yshld[i-1],
		xshld[ i ], yshld[ i ],
		xshld[i+1], yshld[i+1], kThicknessCentral,
		xair[i], yair[i]);
  }
    InsidePoint(xshld[22], yshld[22],
		xshld[23], yshld[23],
		xshld[ 0], yshld[ 0], kThicknessCentral,
		xair[23], yair[23]);

  // Create the air shape
  TGeoXtru *centralairshape = new TGeoXtru(2);

  centralairshape->DefinePolygon(24,xair,yair);
  centralairshape->DefineSection(0,-kHalfLengthCentral);
  centralairshape->DefineSection(1, kHalfLengthCentral);

  // Create the Omega insert
  TGeoXtru *centralomegashape = new TGeoXtru(2);

  CreateSPDOmegaShape(xair,yair,kTheta,kThicknessOmega,xomega,yomega);

  centralomegashape->DefinePolygon(48,xomega,yomega);
  centralomegashape->DefineSection(0,-kHalfLengthCentral);
  centralomegashape->DefineSection(1, kHalfLengthCentral);

  // The end cap half shield: a half tube of carbon fiber,
  // a similar but proportionally smaller half tube of air inside it,
  // and a Omega-shaped carbon fiber insert inside the air.
  // They are all XTru shapes

  TGeoXtru *endcapshape = new TGeoXtru(2);

  CreateSPDThermalShape(kInnerAEndCap,kInnerBEndCap,kInnerRadiusEndCap,
			kOuterAEndCap,kOuterBEndCap,kOuterRadiusEndCap,
			kTheta,xshld,yshld);

  endcapshape->DefinePolygon(24,xshld,yshld);
  endcapshape->DefineSection(0,-kHalfLengthEndCap);
  endcapshape->DefineSection(1, kHalfLengthEndCap);

  // Now rescale to get the air volume dimensions
    InsidePoint(xshld[23], yshld[23],
		xshld[ 0], yshld[ 0],
		xshld[ 1], yshld[ 1], kThicknessEndCap,
		xair[0], yair[0]);
  for (Int_t i=1; i<23; i++) {
    InsidePoint(xshld[i-1], yshld[i-1],
		xshld[ i ], yshld[ i ],
		xshld[i+1], yshld[i+1], kThicknessEndCap,
		xair[i], yair[i]);
  }
    InsidePoint(xshld[22], yshld[22],
		xshld[23], yshld[23],
		xshld[ 0], yshld[ 0], kThicknessEndCap,
		xair[23], yair[23]);

  // Create the air shape
  TGeoXtru *endcapairshape = new TGeoXtru(2);

  endcapairshape->DefinePolygon(24,xair,yair);
  endcapairshape->DefineSection(0,-kHalfLengthEndCap);
  endcapairshape->DefineSection(1, kHalfLengthEndCap);

  // Create the Omega insert
  TGeoXtru *endcapomegashape = new TGeoXtru(2);

  CreateSPDOmegaShape(xair,yair,kTheta,kThicknessOmega,xomega,yomega);

  endcapomegashape->DefinePolygon(48,xomega,yomega);
  endcapomegashape->DefineSection(0,-kHalfLengthEndCap);
  endcapomegashape->DefineSection(1, kHalfLengthEndCap);

  // The cone half shield is more complex since there is no basic
  // TGeo shape to describe it correctly. So it is made of a series
  // of TGeoArb8 shapes filled with air, which all together make up the
  // the cone AND its internal insert. Part of the following code is
  // adapted from SPDThermalSheald method.

  // Filled portions
  TGeoArb8 *sC1 = new TGeoArb8(kHalfLengthCone);
  TGeoArb8 *sC2 = new TGeoArb8(kHalfLengthCone);

  CreateSPDThermalShape(kInnerACentral,kInnerBCentral,kInnerRadiusCentral,
			kOuterACentral,kOuterBCentral,kOuterRadiusCentral,
			kTheta,xshld,yshld);

  sC1->SetVertex(0,xshld[12],yshld[12]);
  sC1->SetVertex(1,xshld[11],yshld[11]);
  sC1->SetVertex(2,xshld[ 0],yshld[ 0]);
  sC1->SetVertex(3,xshld[23],yshld[23]);

  sC2->SetVertex(0,xshld[11],yshld[11]);
  sC2->SetVertex(1,xshld[10],yshld[10]);
  sC2->SetVertex(2,xshld[ 1],yshld[ 1]);
  sC2->SetVertex(3,xshld[ 0],yshld[ 0]);

  // Drawings give only the radius, convert it to the apothegm
  Double_t kInnerRadiusCone = TMath::Sqrt(kInnerRadialCone*kInnerRadialCone
					  - 0.25*kInnerACone*kInnerACone);
  Double_t kOuterRadiusCone = TMath::Sqrt(kOuterRadialCone*kOuterRadialCone
					  - 0.25*kOuterACone*kOuterACone);

  Double_t xco[4], yco[4], xci[4], yci[4];

  for (Int_t i=0; i<2; i++) {
    Double_t th = i*kTheta*TMath::RadToDeg();
    xco[2*i  ] = kOuterRadiusCone*SinD(th) - 0.5*kOuterACone*CosD(th);
    yco[2*i  ] = kOuterRadiusCone*CosD(th) + 0.5*kOuterACone*SinD(th);
    xci[2*i  ] = kInnerRadiusCone*SinD(th) - 0.5*kInnerACone*CosD(th);
    yci[2*i  ] = kInnerRadiusCone*CosD(th) + 0.5*kInnerACone*SinD(th);
    xco[2*i+1] = kOuterRadiusCone*SinD(th) + 0.5*kOuterACone*CosD(th);
    yco[2*i+1] = kOuterRadiusCone*CosD(th) - 0.5*kOuterACone*SinD(th);
    xci[2*i+1] = kInnerRadiusCone*SinD(th) + 0.5*kInnerACone*CosD(th);
    yci[2*i+1] = kInnerRadiusCone*CosD(th) - 0.5*kInnerACone*SinD(th);
  }

  sC1->SetVertex(4,xco[0],yco[0]);
  sC1->SetVertex(5,xco[1],yco[1]);
  sC1->SetVertex(6,xci[1],yci[1]);
  sC1->SetVertex(7,xci[0],yci[0]);

  sC2->SetVertex(4,xco[1],yco[1]);
  sC2->SetVertex(5,xco[2],yco[2]);
  sC2->SetVertex(6,xci[2],yci[2]);
  sC2->SetVertex(7,xci[1],yci[1]);

  // Air holes
  TGeoArb8 *sCh1 = new TGeoArb8(kHalfLengthCone);
  TGeoArb8 *sCh2 = new TGeoArb8(kHalfLengthCone);

  for(Int_t i=0; i<4; i++){
    InsidePoint(sC1->GetVertices()[((i+3)%4)*2+0],
		sC1->GetVertices()[((i+3)%4)*2+1],
		sC1->GetVertices()[i*2+0],
		sC1->GetVertices()[i*2+1],
		sC1->GetVertices()[((i+1)%4)*2+0],
		sC1->GetVertices()[((i+1)%4)*2+1],-kThicknessCone,x,y);
    sCh1->SetVertex(i,x,y);

    InsidePoint(sC1->GetVertices()[((i+3)%4 +4)*2+0],
		sC1->GetVertices()[((i+3)%4 +4)*2+1],
		sC1->GetVertices()[(i+4)*2+0],
		sC1->GetVertices()[(i+4)*2+1],
		sC1->GetVertices()[((i+1)%4 +4)*2+0],
		sC1->GetVertices()[((i+1)%4 +4)*2+1],-kThicknessCone,x,y);
    sCh1->SetVertex(i+4,x,y);

    InsidePoint(sC2->GetVertices()[((i+3)%4)*2+0],
		sC2->GetVertices()[((i+3)%4)*2+1],
		sC2->GetVertices()[i*2+0],
		sC2->GetVertices()[i*2+1],
		sC2->GetVertices()[((i+1)%4)*2+0],
		sC2->GetVertices()[((i+1)%4)*2+1],-kThicknessCone,x,y);
    sCh2->SetVertex(i,x,y);

    InsidePoint(sC2->GetVertices()[((i+3)%4 +4)*2+0],
		sC2->GetVertices()[((i+3)%4 +4)*2+1],
		sC2->GetVertices()[(i+4)*2+0],
		sC2->GetVertices()[(i+4)*2+1],
		sC2->GetVertices()[((i+1)%4 +4)*2+0],
		sC2->GetVertices()[((i+1)%4 +4)*2+1],-kThicknessCone,x,y);
    sCh2->SetVertex(i+4,x,y);
  }

  // Finally the carbon fiber Ring with its Wings and their
  // stesalite inserts. They are Tube and TubeSeg shapes

  TGeoTube *ringshape = new TGeoTube(kInnerRadiusRing,kOuterRadiusRing,
				     kHalfLengthRing);

  TGeoTube *ringinsertshape = new TGeoTube(kInnerRadiusRing+kThicknessRing,
					   kOuterRadiusRing-kThicknessRing,
					   kHalfLengthRing-kThicknessRing);

  Double_t angleWideWing, angleWideWingThickness;
  angleWideWing = (kWideWing/kOuterRadiusWing)*TMath::RadToDeg();
  angleWideWingThickness = (kThicknessRing/kOuterRadiusWing)*TMath::RadToDeg();

  TGeoTubeSeg *wingshape = new TGeoTubeSeg(kOuterRadiusRing,kOuterRadiusWing,
					   kHalfLengthRing, 0, angleWideWing);

  TGeoTubeSeg *winginsertshape = new TGeoTubeSeg(kOuterRadiusRing,
             kOuterRadiusWing-kThicknessRing, kHalfLengthRing-kThicknessRing,
             angleWideWingThickness, angleWideWing-angleWideWingThickness);


  // We have the shapes: now create the real volumes

  TGeoMedium *medSPDcf  = mgr->GetMedium("ITS_SPD shield$");
  TGeoMedium *medSPDair = mgr->GetMedium("ITS_SPD AIR$");
  TGeoMedium *medSPDste = mgr->GetMedium("ITS_G10FR4$"); // stesalite

  TGeoVolume *centralshield = new TGeoVolume("SPDcentralshield",
					     centralshape,medSPDcf);
  centralshield->SetVisibility(kTRUE);
  centralshield->SetLineColor(7);
  centralshield->SetLineWidth(1);

  TGeoVolume *centralairshield = new TGeoVolume("SPDcentralairshield",
						centralairshape,medSPDair);
  centralairshield->SetVisibility(kTRUE);
  centralairshield->SetLineColor(5); // Yellow
  centralairshield->SetLineWidth(1);
  centralairshield->SetFillColor(centralairshield->GetLineColor());
  centralairshield->SetFillStyle(4090); // 90% transparent

  TGeoVolume *centralomega = new TGeoVolume("SPDcentralomega",
					     centralomegashape,medSPDcf);
  centralomega->SetVisibility(kTRUE);
  centralomega->SetLineColor(7);
  centralomega->SetLineWidth(1);

  centralairshield->AddNode(centralomega,1,0);
  centralshield->AddNode(centralairshield,1,0);

  TGeoVolume *endcapshield = new TGeoVolume("SPDendcapshield",
					     endcapshape,medSPDcf);
  endcapshield->SetVisibility(kTRUE);
  endcapshield->SetLineColor(7);
  endcapshield->SetLineWidth(1);

  TGeoVolume *endcapairshield = new TGeoVolume("SPDendcapairshield",
						endcapairshape,medSPDair);
  endcapairshield->SetVisibility(kTRUE);
  endcapairshield->SetLineColor(5); // Yellow
  endcapairshield->SetLineWidth(1);
  endcapairshield->SetFillColor(endcapairshield->GetLineColor());
  endcapairshield->SetFillStyle(4090); // 90% transparent

  TGeoVolume *endcapomega = new TGeoVolume("SPDendcapomega",
					   endcapomegashape,medSPDcf);
  endcapomega->SetVisibility(kTRUE);
  endcapomega->SetLineColor(7);
  endcapomega->SetLineWidth(1);

  endcapairshield->AddNode(endcapomega,1,0);
  endcapshield->AddNode(endcapairshield,1,0);

  TGeoVolume *vC1 = new TGeoVolume("SPDconeshieldV1",sC1,medSPDcf);
  vC1->SetVisibility(kTRUE);
  vC1->SetLineColor(7);
  vC1->SetLineWidth(1);

  TGeoVolume *vCh1 = new TGeoVolume("SPDconeshieldH1",sCh1,medSPDair);

  vCh1->SetVisibility(kTRUE);
  vCh1->SetLineColor(5); // Yellow
  vCh1->SetLineWidth(1);
  vCh1->SetFillColor(vCh1->GetLineColor());
  vCh1->SetFillStyle(4090); // 90% transparent

  vC1->AddNode(vCh1,1,0);

  TGeoVolume *vC2 = new TGeoVolume("SPDconeshieldV2",sC2,medSPDcf);

  vC2->SetVisibility(kTRUE);
  vC2->SetLineColor(7);
  vC2->SetLineWidth(1);

  TGeoVolume *vCh2 = new TGeoVolume("SPDconeshieldH2",sCh2,medSPDair);

  vCh2->SetVisibility(kTRUE);
  vCh2->SetLineColor(5); // Yellow
  vCh2->SetLineWidth(1);
  vCh2->SetFillColor(vCh2->GetLineColor());
  vCh2->SetFillStyle(4090); // 90% transparent

  vC2->AddNode(vCh2,1,0);

  TGeoVolume *ring = new TGeoVolume("SPDshieldring",ringshape,medSPDcf);
  ring->SetVisibility(kTRUE);
  ring->SetLineColor(7);
  ring->SetLineWidth(1);

  TGeoVolume *ringinsert = new TGeoVolume("SPDshieldringinsert",
					  ringinsertshape,medSPDste);
  ringinsert->SetVisibility(kTRUE);
  ringinsert->SetLineColor(3); // Green
//  ringinsert->SetLineWidth(1);
  ringinsert->SetFillColor(ringinsert->GetLineColor());
  ringinsert->SetFillStyle(4010); // 10% transparent

  ring->AddNode(ringinsert,1,0);

  TGeoVolume *wing = new TGeoVolume("SPDshieldringwing",wingshape,medSPDcf);
  wing->SetVisibility(kTRUE);
  wing->SetLineColor(7);
  wing->SetLineWidth(1);

  TGeoVolume *winginsert = new TGeoVolume("SPDshieldringinsert",
					  winginsertshape,medSPDste);
  winginsert->SetVisibility(kTRUE);
  winginsert->SetLineColor(3); // Green
//  winginsert->SetLineWidth(1);
  winginsert->SetFillColor(winginsert->GetLineColor());
  winginsert->SetFillStyle(4010); // 10% transparent

  wing->AddNode(winginsert,1,0);


  // Add all volumes in the assembly
  vM->AddNode(centralshield,1,0);
  vM->AddNode(centralshield,2,new TGeoRotation("",180,0,0));

  vM->AddNode(endcapshield,1,
	      new TGeoTranslation(0,0, kHalfLengthCentral+kHalfLengthEndCap));
  vM->AddNode(endcapshield,2,
	      new TGeoTranslation(0,0,-kHalfLengthCentral-kHalfLengthEndCap));
  vM->AddNode(endcapshield,3,new TGeoCombiTrans(
              0, 0, kHalfLengthCentral+kHalfLengthEndCap,
	      new TGeoRotation("",180,0,0)     ) );
  vM->AddNode(endcapshield,4,new TGeoCombiTrans(
              0, 0,-kHalfLengthCentral-kHalfLengthEndCap,
	      new TGeoRotation("",180,0,0)     ) );

  for (Int_t i=0; i<10; i++) {
    Double_t thetaC12 = kTheta*TMath::RadToDeg();
    vM->AddNode(vC1,2*i+1, new TGeoCombiTrans(
               0, 0,  kHalfLengthCentral+2*kHalfLengthEndCap+kHalfLengthCone,
	       new TGeoRotation("",0,  0,i*thetaC12)   ) );
    vM->AddNode(vC1,2*i+2, new TGeoCombiTrans(
               0, 0, -kHalfLengthCentral-2*kHalfLengthEndCap-kHalfLengthCone,
	       new TGeoRotation("",0,180,i*thetaC12)   ) );
    vM->AddNode(vC2,2*i+1, new TGeoCombiTrans(
               0, 0,  kHalfLengthCentral+2*kHalfLengthEndCap+kHalfLengthCone,
	       new TGeoRotation("",0,  0,i*thetaC12)   ) );
    vM->AddNode(vC2,2*i+2, new TGeoCombiTrans(
               0, 0, -kHalfLengthCentral-2*kHalfLengthEndCap-kHalfLengthCone,
	       new TGeoRotation("",0,180,i*thetaC12)   ) );
  }

  vM->AddNode(ring,1,new TGeoTranslation(0, 0,
	      kHalfLengthCentral+2*kHalfLengthEndCap+2*kHalfLengthCone
             +kHalfLengthRing));
  vM->AddNode(ring,2,new TGeoTranslation(0, 0,
	     -kHalfLengthCentral-2*kHalfLengthEndCap-2*kHalfLengthCone
             -kHalfLengthRing));

  for (Int_t i=0; i<4; i++) {
    Double_t thetaW = kThetaWing*(2*i+1);
    vM->AddNode(wing,2*i+1,new TGeoCombiTrans(0, 0,
	      kHalfLengthCentral+2*kHalfLengthEndCap+2*kHalfLengthCone
             +kHalfLengthRing, new TGeoRotation("",thetaW,0,0)  ));
    vM->AddNode(wing,2*i+2,new TGeoCombiTrans(0, 0,
	     -kHalfLengthCentral-2*kHalfLengthEndCap-2*kHalfLengthCone
             -kHalfLengthRing, new TGeoRotation("",thetaW,0,0)  ));
  }

  // Some debugging if requested
  if(GetDebug(1)){
    vM->PrintNodes();
    vM->InspectShape();
  }

  // Finally put the entire shield in the mother volume
  moth->AddNode(vM,1,0);

  return;
}

//______________________________________________________________________
void AliITSv11GeometrySupport::CreateSPDThermalShape(
     Double_t ina, Double_t inb, Double_t inr,
     Double_t oua, Double_t oub, Double_t our,
     Double_t   t, Double_t *x , Double_t *y )
{
//
// Creates the proper sequence of X and Y coordinates to determine
// the base XTru polygon for the SPD thermal shapes
//
// Input:
//        ina, inb : inner shape sides
//        inr      : inner radius
//        oua, oub : outer shape sides
//        our      : outer radius
//        t        : theta angle
//
// Output:
//        x, y : coordinate vectors [24]
//
// Created:      14 Nov 2007  Mario Sitta
// Updated:      11 Dec 2007  Mario Sitta
//
  Double_t xlocal[6],ylocal[6];

  //Create the first inner quadrant (X > 0)
  FillSPDXtruShape(ina,inb,inr,t,xlocal,ylocal);
  for (Int_t i=0; i<6; i++) {
    x[i] = xlocal[i];
    y[i] = ylocal[i];
  }

  // Then reflex on the second quadrant (X < 0)
  for (Int_t i=0; i<6; i++) {
    x[23-i] = -x[i];
    y[23-i] =  y[i];
  }

  // Now create the first outer quadrant (X > 0)
  FillSPDXtruShape(oua,oub,our,t,xlocal,ylocal);
  for (Int_t i=0; i<6; i++) {
    x[11-i] = xlocal[i];
    y[11-i] = ylocal[i];
  }

  // Finally reflex on the second quadrant (X < 0)
  for (Int_t i=0; i<6; i++) {
    x[12+i] = -x[11-i];
    y[12+i] =  y[11-i];
  }

  return;
}

//______________________________________________________________________
void AliITSv11GeometrySupport::CreateSPDOmegaShape(
                             Double_t *xin, Double_t *yin, Double_t  t,
			     Double_t    d, Double_t   *x, Double_t *y)
{
//
// Creates the proper sequence of X and Y coordinates to determine
// the SPD Omega XTru polygon
//
// Input:
//        xin, yin : coordinates of the air volume
//        d        : Omega shape thickness
//        t        : theta angle
//
// Output:
//        x, y     : coordinate vectors [48]
//
// Created:      17 Nov 2007  Mario Sitta
// Updated:      11 Dec 2007  Mario Sitta
//
  Double_t xlocal[6],ylocal[6];

  // First determine various parameters
  Double_t ina = TMath::Sqrt( (xin[23]-xin[0])*(xin[23]-xin[0]) +
			      (yin[23]-yin[0])*(yin[23]-yin[0]) );
  Double_t inb = TMath::Sqrt( (xin[ 1]-xin[0])*(xin[ 1]-xin[0]) +
			      (yin[ 1]-yin[0])*(yin[ 1]-yin[0]) );
  Double_t inr = yin[0];
  Double_t oua = TMath::Sqrt( (xin[12]-xin[11])*(xin[12]-xin[11]) +
			      (yin[12]-yin[11])*(yin[12]-yin[11]) );
  Double_t oub = TMath::Sqrt( (xin[10]-xin[11])*(xin[10]-xin[11]) +
			      (yin[10]-yin[11])*(yin[10]-yin[11]) );
  Double_t our = yin[11];

  //Create the first inner pseudo-quadrant
  FillSPDXtruShape(ina,inb,inr,t,xlocal,ylocal);
  x[ 1] = xlocal[0];
  y[ 1] = ylocal[0];

  x[ 2] = xlocal[1];
  y[ 2] = ylocal[1];

  x[ 5] = xlocal[2];
  y[ 5] = ylocal[2];

  x[ 6] = xlocal[3];
  y[ 6] = ylocal[3];

  x[ 9] = xlocal[4];
  y[ 9] = ylocal[4];

  x[10] = xlocal[5];
  y[10] = ylocal[5];

  //Create the first outer pseudo-quadrant
  FillSPDXtruShape(oua,oub,our,t,xlocal,ylocal);
  x[23] = xlocal[0];
  y[23] = ylocal[0];

  x[20] = xlocal[1];
  y[20] = ylocal[1];

  x[19] = xlocal[2];
  y[19] = ylocal[2];

  x[16] = xlocal[3];
  y[16] = ylocal[3];

  x[15] = xlocal[4];
  y[15] = ylocal[4];

  x[11] = xlocal[5];
  y[11] = ylocal[5];

  //Create the second inner pseudo-quadrant
  FillSPDXtruShape(ina+2*d,inb-2*d,inr+d,t,xlocal,ylocal);
  x[22] = xlocal[0];
  y[22] = ylocal[0];

  x[21] = xlocal[1];
  y[21] = ylocal[1];

  x[18] = xlocal[2];
  y[18] = ylocal[2];

  x[17] = xlocal[3];
  y[17] = ylocal[3];

  x[14] = xlocal[4];
  y[14] = ylocal[4];

  x[13] = xlocal[5];
  y[13] = ylocal[5];

  //Create the second outer pseudo-quadrant
  FillSPDXtruShape(oua-2*d,oub+2*d,our-d,t,xlocal,ylocal);
  x[ 0] = xlocal[0];
  y[ 0] = ylocal[0];

  x[ 3] = xlocal[1];
  y[ 3] = ylocal[1];

  x[ 4] = xlocal[2];
  y[ 4] = ylocal[2];

  x[ 7] = xlocal[3];
  y[ 7] = ylocal[3];

  x[ 8] = xlocal[4];
  y[ 8] = ylocal[4];

  x[12] = xlocal[5];
  y[12] = ylocal[5];

  // These need to be fixed explicitly
  y[10] = yin[5];
  y[11] = yin[6];
  x[12] = x[11];
  y[12] = y[11] + d;
  x[13] = x[10] + d;
  y[13] = y[12];

  // Finally reflex on the negative side
  for (Int_t i=0; i<24; i++) {
    x[24+i] = -x[23-i];
    y[24+i] =  y[23-i];
  }

  // Wow ! We've finished
  return;
}

//______________________________________________________________________
void AliITSv11GeometrySupport::FillSPDXtruShape(Double_t a, Double_t b,
						Double_t r, Double_t t,
						Double_t *x, Double_t *y)
{
//
// Creates the partial sequence of X and Y coordinates to determine
// the lateral part of the SPD thermal shield
//
// Input:
//        a, b : shape sides
//        r    : radius
//        t    : theta angle
//
// Output:
//        x, y : coordinate vectors [6]
//
// Created:      14 Nov 2007  Mario Sitta
//
  x[0] = a/2;
  y[0] = r;

  x[1] = x[0] + b * TMath::Cos(t/2);
  y[1] = y[0] - b * TMath::Sin(t/2);

  x[2] = x[1] + a * TMath::Cos(t);
  y[2] = y[1] - a * TMath::Sin(t);

  x[3] = x[2] + b * TMath::Cos(3*t/2);
  y[3] = y[2] - b * TMath::Sin(3*t/2);

  x[4] = x[3] + a * TMath::Cos(2*t);
  y[4] = y[3] - a * TMath::Sin(2*t);

  x[5] = x[4];
  y[5] = 0.;

  return;
}

//______________________________________________________________________
void AliITSv11GeometrySupport::SDDCone(TGeoVolume *moth,TGeoManager *mgr)
{
//
// Creates the SDD support cone and cylinder geometry as a
// volume assembly and adds it to the mother volume
// (part of this code is taken or anyway inspired to SDDCone method
// of AliITSv11GeometrySupport.cxx,v 1.9 2007/06/06)
//
// Input:
//         moth : the TGeoVolume owing the volume structure
//         mgr  : the GeoManager (default gGeoManager)
// Output:
//
// Created:         ???       Bjorn S. Nilsen
// Updated:      18 Feb 2008  Mario Sitta
//
// Technical data are taken from:  "Supporto Generale Settore SDD"
// (technical drawings ALR-0816/1-B), "Supporto Globale Settore SDD"
// (technical drawings ALR-0816/2A, ALR-0816/2B, ALR-0816/2C, ALR-0816/2D), 
// private communication with B. Giraudo

  // Dimensions of the Central cylinder and flanges
  const Double_t kCylinderHalfLength = (790.0/2)*fgkmm;
  const Double_t kCylinderInnerR     = (210.0/2)*fgkmm;
  const Double_t kCylinderOuterR     = (231.0/2)*fgkmm;
  const Double_t kFlangeHalfLength   = ( 15.0/2)*fgkmm;
  const Double_t kFlangeInnerR       = (210.5/2)*fgkmm;
  const Double_t kFlangeOuterR       = (230.5/2)*fgkmm;
  const Double_t kInsertoHalfLength  =
                                     kCylinderHalfLength - 2*kFlangeHalfLength;
//  const Double_t kCFThickness        = kFlangeInnerR - kCylinderInnerR;
  const Double_t kBoltDiameter       =       6.0*fgkmm; // M6 screw
  const Double_t kBoltDepth          =       6.0*fgkmm; // In the flange
  const Double_t kBoltRadius         = (220.0/2)*fgkmm; // Radius in flange
  const Double_t kThetaBolt          =      30.0*fgkDegree;
  const Int_t    kNBolts             = (Int_t)(360.0/kThetaBolt);
  // Dimensions of the Cone
  const Double_t kConeROutMin        = (540.0/2)*fgkmm;
  const Double_t kConeROutMax        = (560.0/2)*fgkmm;
  const Double_t kConeRCurv          =      15.0*fgkmm; // Radius of curvature
  const Double_t kConeRinMin         = (210.0/2)*fgkmm;
  const Double_t kConeRinMax         = (216.0/2)*fgkmm;
  const Double_t kConeRinCylinder    = (231.0/2)*fgkmm;
  const Double_t kConeZCylinder      =     186.0*fgkmm;
  const Double_t kConeZOuterMilled   =      23.0*fgkmm;
  const Double_t kConeDZin           =      15.0*fgkmm; // ???
  const Double_t kConeThickness      =      10.5*fgkmm; // Rohacell + Carb.Fib.
  const Double_t kConeTheta          =      45.0*fgkDegree; // SDD cone angle
  const Double_t kSinConeTheta       =
                                     TMath::Sin(kConeTheta*TMath::DegToRad());
  const Double_t kCosConeTheta       =
                                     TMath::Cos(kConeTheta*TMath::DegToRad());
  const Double_t kTanConeTheta       =
                                     TMath::Tan(kConeTheta*TMath::DegToRad());
  // Dimensions of the Cone Inserts
  const Double_t kConeCFThickness       = 1.5*fgkmm; // Carbon fiber thickness
  // Dimensions of the Cone Holes
  const Double_t kHole1RMin          = (450.0/2)*fgkmm;
//  const Double_t kHole1RMax          = (528.0/2)*fgkmm;
  const Double_t kHole1RMax          = (527.4/2)*fgkmm; // ??? No overlaps !
  const Double_t kHole2RMin          = (280.0/2)*fgkmm;
  const Double_t kHole2RMax          = (375.0/2)*fgkmm;
  const Double_t kHole1Phi           =      25.0*fgkDegree;
  const Double_t kHole2Phi           =      50.0*fgkDegree;
  const Double_t kHole3RMin          =     205.0*fgkmm;
  const Double_t kHole3DeltaR        =        15*fgkmm;
  const Double_t kHole3Width         =        30*fgkmm;
  const Int_t    kNHole3             =         6      ;
  const Double_t kHole4RMin          =     116.0*fgkmm;
  const Double_t kHole4DeltaR        =        15*fgkmm;
  //  const Double_t kHole4Width         =        30*fgkmm;
  // const Int_t    kNHole4             =         3      ;

  // Local variables
  Double_t x, y, z, t, dza, rmin, rmax;


  // The master volume which holds everything
  TGeoVolumeAssembly *vM = new TGeoVolumeAssembly("ITSsddCone");

  // Recover the needed materials
  TGeoMedium *medSDDcf  = mgr->GetMedium("ITS_SDD C (M55J)$");
  TGeoMedium *medSDDair = mgr->GetMedium("ITS_SDD AIR$");
  TGeoMedium *medSDDste = mgr->GetMedium("ITS_G10FR4$"); // stesalite
  TGeoMedium *medSDDroh = mgr->GetMedium("ITS_ROHACELL$");
  TGeoMedium *medSDDss  = mgr->GetMedium("ITS_INOX$");

  // First define the geometrical shapes

  // Central cylinder with its internal foam and the lateral flanges:
  // a carbon fiber Tube which contains a rohacell Tube and two
  // stesalite Tube's
  TGeoTube *cylindershape = new TGeoTube(kCylinderInnerR,kCylinderOuterR,
					 kCylinderHalfLength);

  TGeoTube *insertoshape = new TGeoTube(kFlangeInnerR,kFlangeOuterR,
					kInsertoHalfLength);

  TGeoTube *flangeshape = new TGeoTube(kFlangeInnerR,kFlangeOuterR,
				       kFlangeHalfLength);

  // The flange bolt: it is a Tube
  TGeoTube *boltshape = new TGeoTube(0.0, 0.5*kBoltDiameter, 0.5*kBoltDepth);

  // Debug if requested
  if (GetDebug(1)) {
    cylindershape->InspectShape();
    insertoshape->InspectShape();
    flangeshape->InspectShape();
    boltshape->InspectShape();
  }


  // We have the shapes: now create the real volumes

  TGeoVolume *cfcylinder = new TGeoVolume("SDDCarbonFiberCylinder",
					  cylindershape,medSDDcf);
  cfcylinder->SetVisibility(kTRUE);
  cfcylinder->SetLineColor(4); // Blue
  cfcylinder->SetLineWidth(1);
  cfcylinder->SetFillColor(cfcylinder->GetLineColor());
  cfcylinder->SetFillStyle(4000); // 0% transparent

  TGeoVolume *foamcylinder = new TGeoVolume("SDDFoamCylinder",
					    insertoshape,medSDDroh);
  foamcylinder->SetVisibility(kTRUE);
  foamcylinder->SetLineColor(3); // Green
  foamcylinder->SetLineWidth(1);
  foamcylinder->SetFillColor(foamcylinder->GetLineColor());
  foamcylinder->SetFillStyle(4050); // 50% transparent

  TGeoVolume *flangecylinder = new TGeoVolume("SDDFlangeCylinder",
					      flangeshape,medSDDste);
  flangecylinder->SetVisibility(kTRUE);
  flangecylinder->SetLineColor(2); // Red
  flangecylinder->SetLineWidth(1);
  flangecylinder->SetFillColor(flangecylinder->GetLineColor());
  flangecylinder->SetFillStyle(4050); // 50% transparent

  TGeoVolume *bolt = new TGeoVolume("SDDFlangeBolt",boltshape,medSDDss);
  bolt->SetVisibility(kTRUE);
  bolt->SetLineColor(1);  // Black
  bolt->SetLineWidth(1);
  bolt->SetFillColor(bolt->GetLineColor());
  bolt->SetFillStyle(4050); // 50% transparent

  // Mount up the cylinder
  for(Int_t i=0; i<kNBolts; i++){
    t = kThetaBolt*i;
    x = kBoltRadius*TMath::Cos(t);
    y = kBoltRadius*TMath::Sin(t);
    z = kFlangeHalfLength-kBoltDepth;
    flangecylinder->AddNode(bolt, i+1, new TGeoTranslation("",x,y,z));
  }

  cfcylinder->AddNode(foamcylinder,1,0);
  cfcylinder->AddNode(flangecylinder,1,
	      new TGeoTranslation(0, 0, kInsertoHalfLength+kFlangeHalfLength));
  cfcylinder->AddNode(flangecylinder,2,new TGeoCombiTrans(
              0, 0, -kInsertoHalfLength-kFlangeHalfLength,
	      new TGeoRotation("",0,180,0)     ) );


  // SDD Support Cone with its internal inserts: a carbon fiber Pcon
  // with holes which contains a stesalite Pcon which on turn contains a
  // rohacell Pcon

  dza = kConeThickness/kSinConeTheta-(kConeROutMax-kConeROutMin)/kTanConeTheta;

  TGeoPcon *coneshape = new TGeoPcon(0.0, 360.0, 12);

  coneshape->Z(0)     = 0.0;
  coneshape->Rmin(0)  = kConeROutMin;
  coneshape->Rmax(0)  = kConeROutMax;

  coneshape->Z(1)     = kConeZOuterMilled - dza;
  coneshape->Rmin(1)  = coneshape->GetRmin(0);
  coneshape->Rmax(1)  = coneshape->GetRmax(0);

  coneshape->Z(2)     = kConeZOuterMilled;
  coneshape->Rmax(2)  = coneshape->GetRmax(0);

  RadiusOfCurvature(kConeRCurv,0.,coneshape->GetZ(1),
		    coneshape->GetRmin(1),kConeTheta,z,rmin);
  coneshape->Z(3)     = z;
  coneshape->Rmin(3)  = rmin;

  coneshape->Rmin(2)  = RminFrom2Points(coneshape,3,1,coneshape->GetZ(2));

  RadiusOfCurvature(kConeRCurv,0.,coneshape->GetZ(2),
		    coneshape->GetRmax(2),kConeTheta,z,rmax);
  coneshape->Z(4)     = z;
  coneshape->Rmax(4)  = rmax;
  coneshape->Rmin(4)  = RminFromZpCone(coneshape,3,kConeTheta,
				       coneshape->GetZ(4),0.0);

  coneshape->Rmax(3)  = RmaxFrom2Points(coneshape,4,2,coneshape->GetZ(3));

  coneshape->Rmin(7)  = kConeRinMin;

  coneshape->Rmin(8)  = kConeRinMin;

  RadiusOfCurvature(kConeRCurv,90.0,0.0,kConeRinMax,90.0-kConeTheta,z,rmax);
  coneshape->Rmax(8)  = rmax;
  coneshape->Z(8)     = ZFromRmaxpCone(coneshape,4,kConeTheta,
				       coneshape->GetRmax(8));

  coneshape->Z(9)     = kConeZCylinder;
  coneshape->Rmin(9)  = kConeRinMin;

  coneshape->Z(10)    = coneshape->GetZ(9);
  coneshape->Rmin(10) = kConeRinCylinder;

  coneshape->Rmin(11) = kConeRinCylinder;
  coneshape->Rmax(11) = coneshape->GetRmin(11);

  rmin                = coneshape->GetRmin(8);
  RadiusOfCurvature(kConeRCurv,90.0-kConeTheta,
		    coneshape->GetZ(8),coneshape->GetRmax(8),90.0,z,rmax);
  rmax                = kConeRinMax;
  coneshape->Z(11)    = z + (coneshape->GetZ(8)-z)*
                   (coneshape->GetRmax(11)-rmax)/(coneshape->GetRmax(8)-rmax);

  coneshape->Rmax(9)  = RmaxFrom2Points(coneshape,11,8,coneshape->GetZ(9));

  coneshape->Rmax(10) = coneshape->GetRmax(9);

  coneshape->Z(6)     = z - kConeDZin;
  coneshape->Z(7)     = coneshape->GetZ(6);

  coneshape->Rmax(6)  = RmaxFromZpCone(coneshape,4,kConeTheta,
				       coneshape->GetZ(6));

  coneshape->Rmax(7)  = coneshape->GetRmax(6);

  RadiusOfCurvature(kConeRCurv,90.,
		    coneshape->GetZ(6),0.0,90.0-kConeTheta,z,rmin);
  coneshape->Z(5)     = z;
  coneshape->Rmin(5)  = RminFromZpCone(coneshape,3,kConeTheta,z);
  coneshape->Rmax(5)  = RmaxFromZpCone(coneshape,4,kConeTheta,z);

  RadiusOfCurvature(kConeRCurv,90.-kConeTheta,
		    0.0,coneshape->Rmin(5),90.0,z,rmin);
  coneshape->Rmin(6)  = rmin;

  // SDD Cone Insert: another Pcon
  Double_t x0, y0, x1, y1, x2, y2;
  TGeoPcon *coneinsertshape = new TGeoPcon(0.0, 360.0, 9);

  coneinsertshape->Z(0)    = coneshape->GetZ(0) + kConeCFThickness;
  coneinsertshape->Rmin(0) = coneshape->GetRmin(0) + kConeCFThickness;
  coneinsertshape->Rmax(0) = coneshape->GetRmax(0) - kConeCFThickness;

  x0 = coneshape->GetZ(0); y0 = coneshape->GetRmin(0);
  x1 = coneshape->GetZ(1); y1 = coneshape->GetRmin(1);
  x2 = coneshape->GetZ(2); y2 = coneshape->GetRmin(2);
  InsidePoint(x0, y0, x1, y1, x2, y2,  kConeCFThickness, z, rmin);
  coneinsertshape->Z(1)    = z;
  coneinsertshape->Rmin(1) = rmin;
  coneinsertshape->Rmax(1) = coneinsertshape->GetRmax(0);

  x0 = coneshape->GetZ(1); y0 = coneshape->GetRmax(1);
  x1 = coneshape->GetZ(2); y1 = coneshape->GetRmax(2);
  x2 = coneshape->GetZ(3); y2 = coneshape->GetRmax(3);
  InsidePoint(x0, y0, x1, y1, x2, y2, -kConeCFThickness, z, rmax);
  coneinsertshape->Z(2)    = z;
  coneinsertshape->Rmax(2) = rmax;

  x0 = coneshape->GetZ(2); y0 = coneshape->GetRmin(2);
  x1 = coneshape->GetZ(3); y1 = coneshape->GetRmin(3);
  x2 = coneshape->GetZ(4); y2 = coneshape->GetRmin(4);
  InsidePoint(x0, y0, x1, y1, x2, y2,  kConeCFThickness, z, rmin);
  coneinsertshape->Z(3)    = z;
  coneinsertshape->Rmin(3) = rmin;

  x0 = coneinsertshape->GetZ(1); y0 = coneinsertshape->GetRmin(1);
  x1 = coneinsertshape->GetZ(3); y1 = coneinsertshape->GetRmin(3);
  coneinsertshape->Rmin(2) = Yfrom2Points(x0, y0, x1, y1,
					  coneinsertshape->Z(2));

  x0 = coneshape->GetZ(3); y0 = coneshape->GetRmax(3);
  x1 = coneshape->GetZ(4); y1 = coneshape->GetRmax(4);
  x2 = coneshape->GetZ(5); y2 = coneshape->GetRmax(5);
  InsidePoint(x0, y0, x1, y1, x2, y2, -kConeCFThickness, z, rmax);
  coneinsertshape->Z(4)    = z;
  coneinsertshape->Rmax(4) = rmax;

  x0 = coneinsertshape->GetZ(2); y0 = coneinsertshape->GetRmax(2);
  x1 = coneinsertshape->GetZ(4); y1 = coneinsertshape->GetRmax(4);
  coneinsertshape->Rmax(3) = Yfrom2Points(x0, y0, x1, y1,
					  coneinsertshape->Z(3));

  x0 = coneshape->GetZ(4); y0 = coneshape->GetRmin(4);
  x1 = coneshape->GetZ(5); y1 = coneshape->GetRmin(5);
  x2 = coneshape->GetZ(6); y2 = coneshape->GetRmin(6);
  InsidePoint(x0, y0, x1, y1, x2, y2,  kConeCFThickness, z, rmin);
  coneinsertshape->Z(5)    = z;
  coneinsertshape->Rmin(5) = rmin;
  coneinsertshape->Rmax(5) = coneinsertshape->GetRmax(4) -
          kTanConeTheta*(coneinsertshape->GetZ(5) - coneinsertshape->GetZ(4));

  x0 = coneinsertshape->GetZ(3); y0 = coneinsertshape->GetRmin(3);
  x1 = coneinsertshape->GetZ(5); y1 = coneinsertshape->GetRmin(5);
  coneinsertshape->Rmin(4) = Yfrom2Points(x0, y0, x1, y1,
					  coneinsertshape->Z(4));

  x0 = coneshape->GetZ(5); y0 = coneshape->GetRmin(5);
  x1 = coneshape->GetZ(6); y1 = coneshape->GetRmin(6);
  x2 = coneshape->GetZ(7); y2 = coneshape->GetRmin(7);
  InsidePoint(x0, y0, x1, y1, x2, y2,  kConeCFThickness, z, rmin);
  coneinsertshape->Z(6)    = z;
  coneinsertshape->Rmin(6) = rmin;
  coneinsertshape->Rmax(6) = coneinsertshape->GetRmax(4) -
          kTanConeTheta*(coneinsertshape->GetZ(6) - coneinsertshape->GetZ(4));

  coneinsertshape->Z(7)    = coneinsertshape->GetZ(6);
  coneinsertshape->Rmin(7) = coneshape->GetRmin(7) + kConeCFThickness;
  coneinsertshape->Rmax(7) = coneinsertshape->GetRmax(6);

  coneinsertshape->Z(8)    = coneshape->GetZ(9) - kConeCFThickness;
  coneinsertshape->Rmin(8) = coneinsertshape->GetRmin(7);
  coneinsertshape->Rmax(8) = coneinsertshape->GetRmax(4) -
          kTanConeTheta*(coneinsertshape->GetZ(8) - coneinsertshape->GetZ(4));

  // SDD Cone Foam: another Pcon
  TGeoPcon *conefoamshape = new TGeoPcon(0.0, 360.0, 4);

  RadiusOfCurvature(kConeRCurv+kConeCFThickness,0.0,coneinsertshape->GetZ(1),
		    coneinsertshape->GetRmin(1),kConeTheta,z,rmin);

  conefoamshape->Z(0)    = z;
  conefoamshape->Rmin(0) = rmin;
  conefoamshape->Rmax(0) = conefoamshape->GetRmin(0);

  conefoamshape->Z(1)    = conefoamshape->GetZ(0)+
                         (kConeThickness-2.0*kConeCFThickness)/kSinConeTheta;
  conefoamshape->Rmin(1) = RminFromZpCone(coneinsertshape,3,kConeTheta,
					  conefoamshape->GetZ(1));
  conefoamshape->Rmax(1) = RmaxFromZpCone(coneinsertshape,4,kConeTheta,
					  conefoamshape->GetZ(1));

  conefoamshape->Z(2)    = coneshape->GetZ(5)-kConeCFThickness;
  conefoamshape->Rmin(2) = RminFromZpCone(coneinsertshape,3,kConeTheta,
					  conefoamshape->GetZ(2));
  conefoamshape->Rmax(2) = RmaxFromZpCone(coneinsertshape,4,kConeTheta,
					  conefoamshape->GetZ(2));

  conefoamshape->Z(3)    = coneinsertshape->GetZ(5)+
                         (kConeThickness-2.0*kConeCFThickness)*kCosConeTheta;
  conefoamshape->Rmax(3) = RmaxFromZpCone(coneinsertshape,4,kConeTheta,
					  conefoamshape->GetZ(3));
  conefoamshape->Rmin(3) = conefoamshape->GetRmax(3);

  // SDD Cone Holes: Pcon's
  TGeoPcon *hole1shape = new TGeoPcon(-kHole1Phi/2., kHole1Phi, 4);

  hole1shape->Rmin(0) = kHole1RMax;
  hole1shape->Rmax(0) = hole1shape->GetRmin(0);
  hole1shape->Z(0)    = ZFromRminpCone(coneshape,3,kConeTheta,
				       hole1shape->GetRmin(0));

  hole1shape->Rmax(1) = hole1shape->GetRmax(0);
  hole1shape->Z(1)    = ZFromRmaxpCone(coneshape,4,kConeTheta,
				       hole1shape->GetRmax(1));
  hole1shape->Rmin(1) = RminFromZpCone(coneshape,3,kConeTheta,
				       hole1shape->GetZ(1));

  hole1shape->Rmin(2) = kHole1RMin;
  hole1shape->Z(2)    = ZFromRminpCone(coneshape,3,kConeTheta,
				       hole1shape->GetRmin(2));
  hole1shape->Rmax(2) = RmaxFromZpCone(coneshape,4,kConeTheta,
				       hole1shape->GetZ(2));

  hole1shape->Rmin(3) = hole1shape->GetRmin(2);
  hole1shape->Rmax(3) = hole1shape->GetRmin(3);
  hole1shape->Z(3)    = ZFromRmaxpCone(coneshape,4,kConeTheta,
				       hole1shape->GetRmax(3));

  TGeoPcon *hole2shape = new TGeoPcon(-kHole2Phi/2., kHole2Phi, 4);

  hole2shape->Rmin(0) = kHole2RMax;
  hole2shape->Rmax(0) = hole2shape->GetRmin(0);
  hole2shape->Z(0)    = ZFromRminpCone(coneshape,3,kConeTheta,
				       hole2shape->GetRmin(0));

  hole2shape->Rmax(1) = hole2shape->GetRmax(0);
  hole2shape->Z(1)    = ZFromRmaxpCone(coneshape,4,kConeTheta,
				       hole2shape->GetRmax(1));
  hole2shape->Rmin(1) = RminFromZpCone(coneshape,3,kConeTheta,
				       hole2shape->GetZ(1));

  hole2shape->Rmin(2) = kHole2RMin;
  hole2shape->Z(2)    = ZFromRminpCone(coneshape,3,kConeTheta,
				       hole2shape->GetRmin(2));
  hole2shape->Rmax(2) = RmaxFromZpCone(coneshape,4,kConeTheta,
				       hole2shape->GetZ(2));

  hole2shape->Rmin(3) = hole2shape->GetRmin(2);
  hole2shape->Rmax(3) = hole2shape->GetRmin(3);
  hole2shape->Z(3)    = ZFromRmaxpCone(coneshape,4,kConeTheta,
				       hole2shape->GetRmax(3));

  Double_t holePhi;
  holePhi = (kHole3Width/kHole3RMin)*TMath::RadToDeg();

  TGeoPcon *hole3shape = new TGeoPcon(-holePhi/2., holePhi, 4);

  hole3shape->Rmin(0) = kHole3RMin + kHole3DeltaR;
  hole3shape->Rmax(0) = hole3shape->GetRmin(0);
  hole3shape->Z(0)    = ZFromRminpCone(coneshape,3,kConeTheta,
				       hole3shape->GetRmin(0));

  hole3shape->Rmax(1) = hole3shape->GetRmax(0);
  hole3shape->Z(1)    = ZFromRmaxpCone(coneshape,4,kConeTheta,
				       hole3shape->GetRmax(1));
  hole3shape->Rmin(1) = RminFromZpCone(coneshape,3,kConeTheta,
				       hole3shape->GetZ(1));

  hole3shape->Rmin(2) = kHole3RMin;
  hole3shape->Z(2)    = ZFromRminpCone(coneshape,3,kConeTheta,
				       hole3shape->GetRmin(2));
  hole3shape->Rmax(2) = RmaxFromZpCone(coneshape,4,kConeTheta,
				       hole3shape->GetZ(2));

  hole3shape->Rmin(3) = hole3shape->GetRmin(2);
  hole3shape->Rmax(3) = hole3shape->GetRmin(3);
  hole3shape->Z(3)    = ZFromRmaxpCone(coneshape,4,kConeTheta,
				       hole3shape->GetRmax(3));

  TGeoPcon *hole4shape = new TGeoPcon(-holePhi/2., holePhi, 4);

  hole4shape->Rmin(0) = kHole4RMin + kHole4DeltaR;
  hole4shape->Rmax(0) = hole4shape->GetRmin(0);
  hole4shape->Z(0)    = ZFromRminpCone(coneshape,3,kConeTheta,
				       hole4shape->GetRmin(0));

  hole4shape->Rmax(1) = hole4shape->GetRmax(0);
  hole4shape->Z(1)    = ZFromRmaxpCone(coneshape,4,kConeTheta,
				       hole4shape->GetRmax(1));
  hole4shape->Rmin(1) = RminFromZpCone(coneshape,3,kConeTheta,
				       hole4shape->GetZ(1));

  hole4shape->Rmin(2) = kHole4RMin;
  hole4shape->Z(2)    = ZFromRminpCone(coneshape,3,kConeTheta,
				       hole4shape->GetRmin(2));
  hole4shape->Rmax(2) = RmaxFromZpCone(coneshape,4,kConeTheta,
				       hole4shape->GetZ(2));

  hole4shape->Rmin(3) = hole4shape->GetRmin(2);
  hole4shape->Rmax(3) = hole4shape->GetRmin(3);
  hole4shape->Z(3)    = ZFromRmaxpCone(coneshape,4,kConeTheta,
				       hole4shape->GetRmax(3));

  // Debug if requested
  if (GetDebug(1)) {
    coneshape->InspectShape();
    coneinsertshape->InspectShape();
    conefoamshape->InspectShape();
    hole1shape->InspectShape();
    hole2shape->InspectShape();
  }


  // We have the shapes: now create the real volumes

  TGeoVolume *cfcone = new TGeoVolume("SDDCarbonFiberCone",
				      coneshape,medSDDcf);
  cfcone->SetVisibility(kTRUE);
  cfcone->SetLineColor(4); // Blue
  cfcone->SetLineWidth(1);
  cfcone->SetFillColor(cfcone->GetLineColor());
  cfcone->SetFillStyle(4000); // 0% transparent

  TGeoVolume *cfconeinsert = new TGeoVolume("SDDCarbonFiberConeInsert",
					    coneinsertshape,medSDDste);
  cfconeinsert->SetVisibility(kTRUE);
  cfconeinsert->SetLineColor(2); // Red
  cfconeinsert->SetLineWidth(1);
  cfconeinsert->SetFillColor(cfconeinsert->GetLineColor());
  cfconeinsert->SetFillStyle(4050); // 50% transparent

  TGeoVolume *cfconefoam = new TGeoVolume("SDDCarbonFiberConeFoam",
					  conefoamshape,medSDDroh);
  cfconefoam->SetVisibility(kTRUE);
  cfconefoam->SetLineColor(7); // Light blue
  cfconefoam->SetLineWidth(1);
  cfconefoam->SetFillColor(cfconefoam->GetLineColor());
  cfconefoam->SetFillStyle(4050); // 50% transparent

  TGeoVolume *hole1 = new TGeoVolume("SDDCableHole1",
				     hole1shape,medSDDair);
  hole1->SetVisibility(kTRUE);
  hole1->SetLineColor(5); // Yellow
  hole1->SetLineWidth(1);
  hole1->SetFillColor(hole1->GetLineColor());
  hole1->SetFillStyle(4090); // 90% transparent

  TGeoVolume *hole2 = new TGeoVolume("SDDCableHole2",
				     hole2shape,medSDDair);
  hole2->SetVisibility(kTRUE);
  hole2->SetLineColor(5); // Yellow
  hole2->SetLineWidth(1);
  hole2->SetFillColor(hole2->GetLineColor());
  hole2->SetFillStyle(4090); // 90% transparent

  TGeoVolume *hole3 = new TGeoVolume("SDDCableHole3",
				     hole3shape,medSDDair);
  hole3->SetVisibility(kTRUE);
  hole3->SetLineColor(5); // Yellow
  hole3->SetLineWidth(1);
  hole3->SetFillColor(hole3->GetLineColor());
  hole3->SetFillStyle(4090); // 90% transparent

  TGeoVolume *hole4 = new TGeoVolume("SDDCableHole4",
				     hole4shape,medSDDair);
  hole4->SetVisibility(kTRUE);
  hole4->SetLineColor(5); // Yellow
  hole4->SetLineWidth(1);
  hole4->SetFillColor(hole4->GetLineColor());
  hole4->SetFillStyle(4090); // 90% transparent

  // Mount up a cone
  cfconeinsert->AddNode(cfconefoam,1,0);

  cfcone->AddNode(cfconeinsert,1,0);

  for (Int_t i=0; i<12; i++) {
    Double_t phiH = i*30.0;
    cfcone->AddNode(hole1, i+1, new TGeoRotation("", 0, 0, phiH));
  }

  for (Int_t i=0; i<6; i++) {
    Double_t phiH = i*60.0;
    cfcone->AddNode(hole2, i+1, new TGeoRotation("", 0, 0, phiH));
  }

  for (Int_t i=0; i<kNHole3; i++) {
    Double_t phiH0 = 360./(Double_t)kNHole3;
    Double_t phiH  = i*phiH0 + 0.5*phiH0;
    cfcone->AddNode(hole3, i+1, new TGeoRotation("", phiH, 0, 0));
  }
/*
  for (Int_t i=0; i<kNHole4; i++) {
    Double_t phiH0 = 360./(Double_t)kNHole4;
    Double_t phiH  = i*phiH0 + 0.25*phiH0;
    cfcone->AddNode(hole4, i+1, new TGeoRotation("", phiH, 0, 0));
  }
*/
  // Add all volumes in the assembly
  vM->AddNode(cfcylinder,1,0);

  z = coneshape->Z(9);
  vM->AddNode(cfcone,1,new TGeoTranslation(0, 0, -z - kCylinderHalfLength));
  vM->AddNode(cfcone,2,new TGeoCombiTrans (0, 0,  z + kCylinderHalfLength,
		       new TGeoRotation("", 0, 180, 0)                   ));

  // Some debugging if requested
  if(GetDebug(1)){
    vM->PrintNodes();
    vM->InspectShape();
  }

  // Finally put the entire shield in the mother volume
  moth->AddNode(vM,1,0);

  return;
}

//______________________________________________________________________
void AliITSv11GeometrySupport::SSDCone(TGeoVolume *moth,TGeoManager *mgr){
    // Define the detail SSD support cone geometry.
    // Inputs:
    //   TGeoVolume  *moth  The mother volume to place this object.
    //   TGeoManager *mgr   A pointer to the Geo-Manager default gGeoManager
    // Outputs:
    //  none.
    // Return:
    //  none.
    //
    Int_t i,j;
    Double_t t,t0,dt,x,y,z,vl[3],vg[3],x0,y0,rmin,rmax;
    TGeoMedium *medSSDcf  = 0; // SSD support cone Carbon Fiber materal number.
    TGeoMedium *medSSDfs  = 0; // SSD support cone inserto stesalite 4411w.
    TGeoMedium *medSSDfo  = 0; // SSD support cone foam, Rohacell 50A.
    TGeoMedium *medSSDss  = 0; // SSD support cone screw material,Stainless
    TGeoMedium *medSSDair = 0; // SSD support cone Air
    TGeoMedium *medSSDal  = 0; // SSD support cone SDD mounting bracket Al
    medSSDcf = mgr->GetMedium("ITSssdCarbonFiber");
    medSSDfs = mgr->GetMedium("ITSssdStaselite4411w");
    medSSDfo = mgr->GetMedium("ITSssdRohacell50A");
    medSSDss = mgr->GetMedium("ITSssdStainlessSteal");
    medSSDair= mgr->GetMedium("ITSssdAir");
    medSSDal = mgr->GetMedium("ITSssdAl");
    //
    // SSD Central cylinder/Thermal Sheald.
    const Double_t kcylZlength     = 1140.0*fgkmm; //
    const Double_t kcylZFoamlength = 1020.0*fgkmm; //
    const Double_t kcylROuter      = 0.5*595.0*fgkmm; //
    const Double_t kcylRInner      = 0.5*560.5*fgkmm; //
    const Double_t kcylCthick      = 0.64*fgkmm; //
    const Double_t kcylFoamThick   = 5.0*fgkmm; //
    const Double_t kcylRholes      = 0.5*570.0*fgkmm;
    const Double_t kcylZM6         = 6.0*fgkmm; //
    const Double_t kcylRM6         = 0.5*6.0*fgkmm;
    const Double_t kcylPhi0M6      = 4.5*fgkDegree;
    const Int_t    kcylNM6         = 40;
    const Double_t kcylZPin        = 10.0*fgkmm;
    const Double_t kcylRPin        = 0.5*4.0*fgkmm;
    const Double_t kcylPhi0Pin     = (90.0+4.5)*fgkDegree;
    const Int_t    kcylNPin        = 2;
    //
    TGeoPcon *sCA,*sCB;
    TGeoTube *sCC,*sCD,*sCE;
    //
    //Begin_Html
    /*
      <img src="picts/ITS/file_name.gif">
      <P>
      <FONT FACE'"TIMES">
      ITS SSD centreal support and thermal sheal cylinder.
      </FONT>
      </P>
     */
    //End_Html
    //
    sCC = new TGeoTube("ITS SSD Thermal Centeral Rohacell CylinderCC",
                       kcylROuter-kcylCthick-kcylFoamThick,
                       kcylROuter-kcylCthick,0.5*kcylZFoamlength);
    sCA = new TGeoPcon("ITS SSD Thermal Centeral Carbon Fiber CylinderCA",
                       0.0,360.0,6);
    sCB = new TGeoPcon("ITS SSD Thermal Centeral Stesalite CylinderCB",
                       0.0,360.0,6);
    sCA->Z(0)    = -0.5*kcylZlength;
    sCA->Rmin(0) = kcylRInner;
    sCA->Rmax(0) = kcylROuter;
    sCA->Z(1)    = sCA->GetZ(0) + kcylZM6;
    sCA->Rmin(1) = sCA->GetRmin(0);
    sCA->Rmax(1) = sCA->GetRmax(0);
    sCA->Z(2)    = -0.5*kcylZFoamlength;
    sCA->Rmin(2) = kcylROuter - 2.0*kcylCthick-kcylFoamThick;
    sCA->Rmax(2) = sCA->GetRmax(0);
    sCA->Z(3)    = -sCA->GetZ(2);
    sCA->Rmin(3) = sCA->GetRmin(2);
    sCA->Rmax(3) = sCA->GetRmax(2);
    sCA->Z(4)    = -sCA->GetZ(1);
    sCA->Rmin(4) = sCA->GetRmin(1);
    sCA->Rmax(4) = sCA->GetRmax(1);
    sCA->Z(5)    = -sCA->GetZ(0);
    sCA->Rmin(5) = sCA->GetRmin(0);
    sCA->Rmax(5) = sCA->GetRmax(0);
    //
    sCB->Z(0)    = sCA->GetZ(0);
    sCB->Rmin(0) = sCA->GetRmin(0) + kcylCthick;
    sCB->Rmax(0) = sCA->GetRmax(0) - kcylCthick;
    sCB->Z(1)    = sCA->GetZ(1);
    sCB->Rmin(1) = sCA->GetRmin(1) + kcylCthick;
    sCB->Rmax(1) = sCA->GetRmax(1) - kcylCthick;
    sCB->Z(2)    = sCA->GetZ(2);
    sCB->Rmin(2) = sCA->GetRmin(2) + kcylCthick;
    sCB->Rmax(2) = sCA->GetRmax(2) - kcylCthick;
    sCB->Z(3)    = sCA->GetZ(3);
    sCB->Rmin(3) = sCA->GetRmin(3) + kcylCthick;
    sCB->Rmax(3) = sCA->GetRmax(3) - kcylCthick;
    sCB->Z(4)    = sCA->GetZ(4);
    sCB->Rmin(4) = sCA->GetRmin(4) + kcylCthick;
    sCB->Rmax(4) = sCA->GetRmax(4) - kcylCthick;
    sCB->Z(5)    = sCA->GetZ(5);
    sCB->Rmin(5) = sCA->GetRmin(5) + kcylCthick;
    sCB->Rmax(5) = sCA->GetRmax(5) - kcylCthick;
    //
    sCD = new TGeoTube("ITS SSD Thermal Centeral Cylinder M6 screwCD",
                      0.0,kcylRM6,0.5*kcylZM6);
    sCE = new TGeoTube("ITS SSD Thermal Centeral Cylinder PinCE",
                      0.0,kcylRPin,0.5*kcylZPin);
    //
    if(GetDebug(1)){
        sCA->InspectShape();
        sCB->InspectShape();
        sCC->InspectShape();
        sCD->InspectShape();
        sCE->InspectShape();
    } // end if GetDegut()
    TGeoVolume *vCA,*vCB,*vCC,*vCD,*vCE;
    vCA = new TGeoVolume("ITSssdCentCylCA",sCA,medSSDcf);
    vCA->SetVisibility(kTRUE);
    vCA->SetLineColor(4); // blue
    vCA->SetLineWidth(1);
    vCA->SetFillColor(vCA->GetLineColor());
    vCA->SetFillStyle(4000); // 0% transparent
    vCB = new TGeoVolume("ITSssdCentCylCB",sCB,medSSDfs);
    vCB->SetVisibility(kTRUE);
    vCB->SetLineColor(2); // red
    vCB->SetLineWidth(1);
    vCB->SetFillColor(vCB->GetLineColor());
    vCB->SetFillStyle(4050); // 50% transparent
    vCC = new TGeoVolume("ITSssdCentCylCC",sCC,medSSDfo);
    vCC->SetVisibility(kTRUE);
    vCC->SetLineColor(3); // green
    vCC->SetLineWidth(1);
    vCC->SetFillColor(vCC->GetLineColor());
    vCC->SetFillStyle(4050); // 50% transparent
    vCD = new TGeoVolume("ITSssdCentCylCD",sCD,medSSDss);
    vCD->SetVisibility(kTRUE);
    vCD->SetLineColor(1); // black
    vCD->SetLineWidth(1);
    vCD->SetFillColor(vCD->GetLineColor());
    vCD->SetFillStyle(4000); // 0% transparent
    vCE = new TGeoVolume("ITSssdCentCylCE",sCE,medSSDss);
    vCE->SetVisibility(kTRUE);
    vCE->SetLineColor(1); // black
    vCE->SetLineWidth(1);
    vCE->SetFillColor(vCE->GetLineColor());
    vCE->SetFillStyle(4000); // 0% transparent
    // Insert Bolt and Pins in both the Cone and Cylinder at the same time.
    vCB->AddNode(vCC,1,0);
    vCA->AddNode(vCB,1,0);
    moth->AddNode(vCA,1,0);
    if(GetDebug(1)){
        vCA->PrintNodes();
        vCB->PrintNodes();
        vCC->PrintNodes();
        vCD->PrintNodes();
        vCE->PrintNodes();
    } // end if
    //
    // SSD Cone
    // Data from Drawings ALR 0743/2E "Supporto Globale Settore SSD" and 
    // ALR 0743/2A "Supporto Generale Settore SSD".
    //
    const Double_t kconThick            = 13.0*fgkmm; // Thickness of Cone.
    const Double_t kconCthick           = 0.75*fgkmm; // Car. finber thickness
    const Double_t kconRCurv0           = 10.0*fgkmm; // Radius of curvature.
    const Double_t kconRCurv1           = 25.0*fgkmm; // Radius of curvature.
    const Double_t kconT                = 39.0*fgkDegree; // angle of SSD cone.
    const Double_t kconZOuterRing       = 47.0*fgkmm;
    const Double_t kconZOuterRingMill   = kconZOuterRing-5.0*fgkmm;
    const Double_t kconZToCylinder      = 170.0*fgkmm;
    const Double_t kconZLengthMill      = 171.5*fgkmm;
    const Double_t kconZLength          = 176.5*fgkmm-
                                          (kconZOuterRing-kconZOuterRingMill);
    //const Double_t kconZInnerRing       = 161.5*fgkmm-
    //                                     (kconZOuterRing-kconZOuterRingMill);
    const Double_t kconZOuterRingInside = 30.25*fgkmm-
                                          (kconZOuterRing-kconZOuterRingMill);
    const Double_t kconZDisplacement    = kconZToCylinder + 0.5*kcylZlength;
    const Double_t kconROuterMax        = 0.5*985.0*fgkmm;
    const Double_t kconROuterMin        = 0.5*945.0*fgkmm;
    const Double_t kconRCylOuterMill    = 0.5*597.0*fgkmm;
    const Double_t kconRInnerMin        = 0.5*562.0*fgkmm;
    //const Double_t kconRCentCurv0       = 0.5*927.0*fgkmm;
    const Double_t kconRCentCurv1       = 0.5*593.0*fgkmm;
    const Double_t kconRCentCurv2       = 0.5*578.0*fgkmm;
    // Foam core.
    const Double_t kconRohacellL0       = 112.3*fgkmm;
    const Double_t kconRohacellL1       = 58.4*fgkmm;
    // Screws and pins in outer SSD cone ring
    const Double_t kconROutHoles        = 0.5*965.0*fgkmm;
    const Double_t kconRScrewM5by12     = 0.5*5.0*fgkmm;
    const Double_t kconLScrewM5by12     = 0.5*12.0*fgkmm;
    const Int_t    kconNScrewM5by12     = 2;
    const Double_t kconRPinO6           = 0.5*6.0*fgkmm;
    const Double_t kconLPinO6           = 0.5*10.0*fgkmm;
    const Int_t    kconNPinO6           = 3;
    const Int_t    kconNRailScrews      = 4;
    const Int_t    kconNRailPins        = 2;
    const Int_t    kconNmounts          = 4;
    const Double_t kconMountPhi0        = 9.0*fgkDegree; // degrees
    //
    const Double_t kconCableHoleROut    = 0.5*920.0*fgkmm;
    const Double_t kconCableHoleRinner  = 0.5*800.0*fgkmm;
    const Double_t kconCableHoleWidth   = 200.0*fgkmm;
    const Double_t kconCableHoleAngle   = 42.0*fgkDegree;
    //const Double_t kconCableHolePhi0    = 90.0/4.0*fgkDegree;
    //const Int_t    kconNCableHoles      = 8;
    const Double_t kconCoolHoleWidth    = 40.0*fgkmm;
    const Double_t kconCoolHoleHight    = 30.0*fgkmm;
    const Double_t kconCoolHoleRmin     = 350.0*fgkmm;
    //const Double_t kconCoolHolephi0     = 90.0/4.0*fgkDegree;
    //const Int_t    kconNCoolHoles       = 8;
    const Double_t kconMountHoleWidth   = 20.0*fgkmm;
    const Double_t kconMountHoleHight   = 20.0*fgkmm;
    const Double_t kconMountHoleRmin    = 317.5*fgkmm;
    //const Double_t kconMountHolephi0    = 0.0*fgkDegree;
    //const Int_t    kconNMountHoles      = 6;
    // SSD cone Wings with holes.
    const Double_t kconWingRmax         = 527.5*fgkmm;
    const Double_t kconWingWidth        = 70.0*fgkmm;
    const Double_t kconWingThick        = 10.0*fgkmm;
    const Double_t kconWingPhi0         = 45.0*fgkDegree;
    //const Int_t    kconNWings           = 4;
    // SSD-SDD Thermal/Mechanical cylinder mounts
    const Double_t kconRM6Head          = 0.5*8.0*fgkmm;
    const Double_t kconZM6Head          = 8.5*fgkmm;
    //
    // SSD-SDD Mounting bracket
    const Double_t ksupPRmin            = 0.5*539.0*fgkmm;// see SDD RoutMin
    const Double_t ksupPRmax            = 0.5*585.0*fgkmm;
    const Double_t ksupPZ               = 4.0*fgkmm;
    const Double_t ksupPPhi1            = (-0.5*70.*fgkmm/ksupPRmax)*fgkRadian;
    const Double_t ksupPPhi2            = -ksupPPhi1;
    //
    const Double_t kSinkconTc           = SinD(kconT);
    const Double_t kCoskconTc           = CosD(kconT);
    //
    TGeoPcon *sA0,*sB0,*sC0,*sF0,*sQ;
    TGeoConeSeg *sAh1,*sBh1;
    TGeoArb8 *sAh2,*sBh2;
    TGeoBBox *sAh3,*sBh3,*sAh4,*sBh4;
    TGeoConeSeg *sG,*sH;
    TGeoTubeSeg *sT;
    TGeoTube *sD,*sE,*sR,*sS;
    TGeoCompositeShape *sA,*sB,*sC,*sF;
    //
    // Lets start with the upper left outer carbon fiber surface.
    // Between za[2],rmaxa[2] and za[4],rmaxa[4] there is a curved section
    // given by rmaxa = rmaxa[2]-r*Sind(t) for 0<=t<=kconT and 
    // za = za[2] + r*Cosd(t) for 0<=t<=kconT. Simularly between za[1],rmina[1
    // and za[3],rmina[3] there is a curve section given by
    // rmina = rmina[1]-r*Sind(t) for 0<=t<=kconT and za = za[1]+r&Sind(t)
    // for t<=0<=kconT. These curves have been replaced by straight lines
    // between the equivelent points for simplicity.
    // Poly-cone Volume sA0. Top part of SSD cone Carbon Fiber.
    sA0 = new TGeoPcon("ITSssdSuportConeCarbonFiberSurfaceA0",0.0,360.0,15);
    sA0->Z(0)    = 0.0;
    sA0->Rmin(0) = kconROuterMin;
    sA0->Rmax(0) = kconROuterMax;
    sA0->Z(1)    = kconZOuterRingInside-kconRCurv0;
    sA0->Rmin(1) = sA0->GetRmin(0);
    sA0->Rmax(1) = sA0->GetRmax(0);
    sA0->Z(2)    = kconZOuterRingInside;
    sA0->Rmin(2) = sA0->GetRmin(1)-kconRCurv0;
    sA0->Rmax(2) = sA0->GetRmax(0);
    sA0->Z(3)    = sA0->GetZ(2);
    sA0->Rmin(3) = -1000; // See Below
    sA0->Rmax(3) = sA0->GetRmax(0);
    sA0->Z(4)    = kconZOuterRingMill-kconRCurv0;
    sA0->Rmin(4) = -1000; // See Below
    sA0->Rmax(4) = sA0->GetRmax(0);
    sA0->Z(5)    = kconZOuterRingMill;
    sA0->Rmin(5) = -1000; // See Below
    sA0->Rmax(5) = sA0->GetRmax(4) - kconRCurv0;
    sA0->Z(6)    = sA0->GetZ(5);
    sA0->Rmin(6) = -1000; // See Below
    sA0->Rmax(6) = -1000; // See Below
    sA0->Z(7)    = sA0->GetZ(6)+kconRCurv0*(1.-kCoskconTc);
    sA0->Rmin(7) = -1000; // See Below
    sA0->Rmax(7) = -1000; // See Below
    sA0->Z(8)    = -1000; // See Below
    sA0->Rmin(8) = kconRCentCurv2+kconRCurv1*kSinkconTc; // See Below
    sA0->Rmax(8) = -1000; // See Below
    sA0->Z(9)    = -1000; // See Below
    sA0->Rmin(9) = kconRCentCurv2;
    sA0->Rmax(9) = -1000; // See Below
    sA0->Z(10)   = -1000; // See Below
    sA0->Rmin(10)= kconRInnerMin;
    sA0->Rmax(10)= -1000; // See Below
    sA0->Z(11)   = kconZLengthMill-kconRCurv0*(1.0-kCoskconTc);
    sA0->Rmin(11)= sA0->GetRmin(10);
    sA0->Rmax(11)= kconRCentCurv1+kconRCurv0*kSinkconTc;
    sA0->Z(12)   = kconZToCylinder;
    sA0->Rmin(12)= sA0->GetRmin(10);
    sA0->Rmax(12)= -1000; // See Below
    sA0->Z(13)   = sA0->GetZ(12);
    sA0->Rmin(13)= kconRCylOuterMill;
    sA0->Rmax(13)= -1000; // See Below
    z            = kconZLengthMill;
    rmin         = kconRCentCurv1;
    rmax         = rmin;
    sA0->Z(14)   = -1000; // See Below
    sA0->Rmin(14)= sA0->GetRmin(13);
    sA0->Rmax(14)= sA0->GetRmin(14);
    // Compute values undefined above
    sA0->Z(14)   = Xfrom2Points(sA0->GetZ(11),sA0->GetRmax(11),z,rmax,
                               sA0->GetRmax(14));
    sA0->Z(8)    = ZFromRmaxpCone(sA0,11,90.-kconT,sA0->GetRmin(8),-kconThick);
    sA0->Rmax(8) = RmaxFromZpCone(sA0,11,90.-kconT,sA0->GetZ(8),0.0);
    sA0->Z(9)    = sA0->GetZ(8)+kconRCurv1*(1.-kCoskconTc);
    sA0->Z(10)   = sA0->GetZ(9);
    sA0->Rmin(3) = RminFromZpCone(sA0,8,90.-kconT,sA0->GetZ(3),0.0);
    sA0->Rmin(4) = RminFromZpCone(sA0,3,90.-kconT,sA0->GetZ(4),0.0);
    sA0->Rmin(5) = RminFromZpCone(sA0,3,90.-kconT,sA0->GetZ(5),0.0);
    sA0->Rmin(7) = RminFromZpCone(sA0,3,90.-kconT,sA0->GetZ(7),0.0);
    sA0->Rmax(7) = RmaxFromZpCone(sA0,11,90.-kconT,sA0->GetZ(7),0.0);
    sA0->Rmin(6) = sA0->GetRmin(5);
    sA0->Rmax(6) = RmaxFromZpCone(sA0,11,90.-kconT,sA0->GetZ(7),0.0);
    sA0->Rmax(9) = RmaxFromZpCone(sA0,11,90.-kconT,sA0->GetZ(9),0.0);
    sA0->Rmax(10)= sA0->GetRmax(9);
    t = TanD(270.+kconT);
    sA0->Rmax(12)= RmaxFrom2Points(sA0,11,14,sA0->GetZ(12));
    sA0->Rmax(13)= sA0->GetRmax(12);
    //
    // Poly-cone Volume B. Stesalite inside volume sA0.
    // Now lets define the Inserto Stesalite 4411w material volume.
    // Poly-cone Volume sA0. Top part of SSD cone Carbon Fiber.
    sB0 = new TGeoPcon("ITSssdSuportConeStaseliteB0",0.0,360.0,15);
    //
    sB0->Z(0)    = sA0->GetZ(0);
    sB0->Rmin(0) = sA0->GetRmin(0) + kconCthick;
    sB0->Rmax(0) = sA0->GetRmax(0) - kconCthick;
    //printf("A0#%d ",1);
    InsidePoint(sA0,0,1,2,kconCthick,sB0,1,kFALSE); // Rmin
    sB0->Rmax(1) = sB0->Rmax(0);
    //printf("A0#%d ",2);
    InsidePoint(sA0,1,2,3,kconCthick,sB0,2,kFALSE); // Rmin
    sB0->Rmax(2) = sB0->Rmax(0);
    //printf("A0#%d ",3);
    InsidePoint(sA0,2,3,9,kconCthick,sB0,3,kFALSE);
    sB0->Rmax(3) = sB0->Rmax(0);
    //printf("A0#%d ",4);
    InsidePoint(sA0,0,4,5,kconCthick,sB0,4,kTRUE); // Rmax
    sB0->Rmin(4) = -1000.; // see Bellow
    //printf("A0#%d ",5);
    InsidePoint(sA0,4,5,6,kconCthick,sB0,5,kTRUE); // Rmax
    sB0->Rmin(5) = -1000.; // see Bellow
    //printf("A0#%d ",6);
    InsidePoint(sA0,5,6,7,kconCthick,sB0,6,kTRUE); // Rmax
    sB0->Rmin(6) = -1000.; // see Bellow
    //printf("A0#%d ",7);
    InsidePoint(sA0,6,7,11,kconCthick,sB0,7,kTRUE); // Rmax
    sB0->Rmin(7) = -1000.; // see Bellow
    //printf("A0#%d ",8);
    InsidePoint(sA0,3,8,9,kconCthick,sB0,8,kFALSE); // Rmin
    sB0->Rmax(8) = -1000.; // see Bellow
    //printf("A0#%d ",9);
    InsidePoint(sA0,8,9,10,kconCthick,sB0,9,kFALSE); // Rmin
    sB0->Rmax(9) = -1000.; // see Bellow
    sB0->Z(10)   = sA0->GetZ(10) + kconCthick;
    sB0->Rmin(10)= sA0->GetRmin(10);
    sB0->Rmax(10)= -1000.; // see Bellow
    //printf("A0#%d ",11);
    InsidePoint(sA0,7,11,14,kconCthick,sB0,11,kTRUE); // Rmax
    sB0->Rmin(11)= sA0->GetRmin(10);
    sB0->Z(12)    = sA0->GetZ(12);
    sB0->Rmin(12)= sA0->GetRmin(12);
    sB0->Rmax(12)= -1000.; // see Bellow
    sB0->Z(13)   = sA0->GetZ(13);
    sB0->Rmin(13)= sA0->GetRmin(13);
    sB0->Rmax(13)= -1000.; // see Bellow
    sB0->Z(14)   = sA0->GetZ(14) - kconCthick;
    sB0->Rmin(14)= sA0->GetRmin(14);
    sB0->Rmax(14)= sB0->Rmin(14); // Close?
    sB0->Rmin(4) = RminFrom2Points(sB0,3,8,sB0->GetZ(4));
    sB0->Rmin(5) = RminFrom2Points(sB0,3,8,sB0->GetZ(5));
    sB0->Rmin(6) = sB0->GetRmin(5);
    sB0->Rmin(7) = RminFrom2Points(sB0,3,8,sB0->GetZ(7));
    sB0->Rmax(8) = RmaxFrom2Points(sB0,7,11,sB0->GetZ(8));
    sB0->Rmax(9) = RmaxFrom2Points(sB0,7,11,sB0->GetZ(9));
    sB0->Rmax(10)= sB0->GetRmax(9);
    sB0->Rmax(12)= RmaxFrom2Points(sB0,11,14,sB0->GetZ(12));
    sB0->Rmax(13)= RmaxFrom2Points(sB0,11,14,sB0->GetZ(13));
    //
    // Poly-cone Volume sC0. Foam inside volume sA0.
    // Now lets define the Rohacell foam material volume.
    sC0 = new TGeoPcon("ITSssdSuportConeRohacellC0",0.0,360.0,4);
    sC0->Z(1)    = sB0->GetZ(7);
    sC0->Rmax(1) = sB0->GetRmax(7);
    sC0->Rmin(1) = RminFrom2Points(sB0,3,8,sC0->GetZ(1));
    sC0->Rmin(0) = sC0->GetRmax(1);
    sC0->Rmax(0) = sC0->GetRmin(0);
    sC0->Z(0)    = Zfrom2MinPoints(sB0,3,8,sC0->Rmin(0));
    t = kconThick-2.0*kconCthick;
    sC0->Rmax(3) = sC0->GetRmax(0)-kCoskconTc*TMath::Sqrt(
                             kconRohacellL0*kconRohacellL0-t*t)+t*kSinkconTc;
    sC0->Rmin(3) = sC0->GetRmax(3);
    sC0->Z(3)    = ZFromRmaxpCone(sB0,11,90.-kconT,sC0->GetRmax(3),0.0);;
    sC0->Rmin(2) = sC0->GetRmin(3);
    sC0->Z(2)    = ZFromRminpCone(sB0,3,90.-kconT,sC0->GetRmin(2),0.0);
    sC0->Rmax(2) = RmaxFromZpCone(sB0,11,90.0-kconT,sC0->GetZ(2),0.0);
    //
    // Poly-cone Volume sF0.  Second Foam inside volume sA0.
    // Now lets define the Rohacell foam material volume.
    sF0 = new TGeoPcon("ITSssdSuportConeRohacellCF0",0.0,360.0,4);
    sF0->Z(2)    = sB0->GetZ(8);
    sF0->Rmin(2) = sB0->GetRmin(8);
    sF0->Rmax(2) = sB0->GetRmax(8);
    sF0->Z(0)    = sF0->GetZ(2)-kconRohacellL1*kSinkconTc;
    sF0->Rmin(0) = sF0->GetRmin(2)+kconRohacellL1*kCoskconTc;
    sF0->Rmax(0) = sF0->GetRmin(0);
    sF0->Z(1)    = ZFromRmaxpCone(sB0,11,90.-kconT,sF0->GetRmax(0),0.0);;
    sF0->Rmax(1) = sF0->GetRmax(0);
    sF0->Rmin(1) = RminFrom2Points(sB0,3,8,sF0->GetZ(1));
    sF0->Rmax(3) = sF0->GetRmin(2)+(kconThick-2.0*kconCthick)*kCoskconTc;
    sF0->Rmin(3) = sF0->GetRmax(3);
    sF0->Z(3)    = ZFromRmaxpCone(sB0,11,90.-kconT,sF0->GetRmax(3),0.0);
    // Holes for Cables to pass Through is created by the intersection
    // between a cone segment and an Arb8, One for the volume sA0 and a
    // larger one for the volumes sB0 and sC0, so that the surface is covered
    // in carbon figer (volume sA0).
    sAh1 = new TGeoConeSeg("ITSssdCableHoleAh1",
                           0.5*kconZLength,kconCableHoleRinner,
                           kconCableHoleROut,kconCableHoleRinner,
                           kconCableHoleROut,
                           90.-(0.5*kconCableHoleWidth/
                                kconCableHoleROut)*fgkRadian,
                           90.+(0.5*kconCableHoleWidth/
                                kconCableHoleROut)*fgkRadian);
    sBh1 = new TGeoConeSeg("ITSssdCableHoleBh1",0.5*kconZLength,
                           kconCableHoleRinner-kconCthick,
                           kconCableHoleROut+kconCthick,
                           kconCableHoleRinner-kconCthick,
                           kconCableHoleROut+kconCthick,
                           90.-(((0.5*kconCableHoleWidth+kconCthick)/
                                 (kconCableHoleROut+kconCthick)))*fgkRadian,
                           90.+(((0.5*kconCableHoleWidth+kconCthick)/
                                 (kconCableHoleROut+kconCthick)))*fgkRadian);
    x0 = sAh1->GetRmax1()*CosD(sAh1->GetPhi2());
    y0 = sAh1->GetRmax1()*SinD(sAh1->GetPhi2());
    sAh2 = new TGeoArb8("ITSssdCableHoleAh2",0.5*kconZLength);
    y  = sAh1->GetRmax1();
    x  = x0+(y-y0)/TanD(90.0+kconCableHoleAngle);
    sAh2->SetVertex(0,x,y);
    y  = sAh1->GetRmin1()*SinD(sAh1->GetPhi2());
    x  = x0+(y-y0)/TanD(90.0+kconCableHoleAngle);
    sAh2->SetVertex(3,x,y);
    x0 = sAh1->GetRmax1()*CosD(sAh1->GetPhi1());
    y0 = sAh1->GetRmax1()*SinD(sAh1->GetPhi1());
    y  = sAh1->GetRmax1();
    x  = x0+(y-y0)/TanD(90.0-kconCableHoleAngle);
    sAh2->SetVertex(1,x,y);
    y  = sAh1->GetRmin1()*SinD(sAh1->GetPhi1());
    x  = x0+(y-y0)/TanD(90.0-kconCableHoleAngle);
    sAh2->SetVertex(2,x,y);
    //
    x0 = sBh1->GetRmax1()*CosD(sBh1->GetPhi2());
    y0 = sBh1->GetRmax1()*SinD(sBh1->GetPhi2());
    sBh2 = new TGeoArb8("ITSssdCableHoleBh2",0.5*kconZLength);
    y  = sBh1->GetRmax1();
    x  = x0+(y-y0)/TanD(90.0+kconCableHoleAngle);
    sBh2->SetVertex(0,x,y);
    y  = sBh1->GetRmin1()*SinD(sBh1->GetPhi2());
    x  = x0+(y-y0)/TanD(90.0+kconCableHoleAngle);
    sBh2->SetVertex(3,x,y);
    x0 = sBh1->GetRmax1()*CosD(sBh1->GetPhi1());
    y0 = sBh1->GetRmax1()*SinD(sBh1->GetPhi1());
    y  = sBh1->GetRmax1();
    x  = x0+(y-y0)/TanD(90.0-kconCableHoleAngle);
    sBh2->SetVertex(1,x,y);
    y  = sBh1->GetRmin1()*SinD(sBh1->GetPhi1());
    x  = x0+(y-y0)/TanD(90.0-kconCableHoleAngle);
    sBh2->SetVertex(2,x,y);
    for(i=0;i<4;i++){ // define points at +dz
        sAh2->SetVertex(i+4,(sAh2->GetVertices())[2*i],
                           (sAh2->GetVertices())[1+2*i]);
        sBh2->SetVertex(i+4,(sBh2->GetVertices())[2*i],
                           (sBh2->GetVertices())[1+2*i]);
    } // end for i
    sAh3 = new TGeoBBox("ITSssdCoolingHoleAh3",0.5*kconCoolHoleWidth,
                        0.5*kconCoolHoleHight,kconZLength);
    sBh3 = new TGeoBBox("ITSssdCoolingHoleBh3",
                        0.5*kconCoolHoleWidth+kconCthick,
                        0.5*kconCoolHoleHight+kconCthick,kconZLength);
    sAh4 = new TGeoBBox("ITSssdMountingPostHoleAh4",0.5*kconMountHoleWidth,
                        0.5*kconMountHoleHight,0.5*kconZLength);
    z = sF0->GetZ(0)-sF0->GetZ(sF0->GetNz()-1);
    if(z<0.0) z = -z;
    sBh4 = new TGeoBBox("ITSssdMountingPostHoleBh4",
                        0.5*kconMountHoleWidth+kconCthick,
                        0.5*kconMountHoleHight+kconCthick,0.5*z);
    // SSD Cone Wings
    sG = new TGeoConeSeg("ITSssdWingCarbonFiberSurfaceG",
                         0.5*kconWingThick,kconROuterMax-kconCthick,
                         kconWingRmax,kconROuterMax-kconCthick,kconWingRmax,
                      kconWingPhi0-(0.5*kconWingWidth/kconWingRmax)*fgkRadian,
                      kconWingPhi0+(0.5*kconWingWidth/kconWingRmax)*fgkRadian);
    sH = new TGeoConeSeg("ITSssdWingStaseliteH",
                         0.5*kconWingThick-kconCthick,kconROuterMax-kconCthick,
                         kconWingRmax-kconCthick,
                         kconROuterMax-kconCthick,
                         kconWingRmax-kconCthick,
                         kconWingPhi0-((0.5*kconWingWidth-kconCthick)/
                                       (kconWingRmax-kconCthick))*fgkRadian,
                         kconWingPhi0+((0.5*kconWingWidth-kconCthick)/
                                       (kconWingRmax-kconCthick))*fgkRadian);
    // SDD support plate, SSD side.
    //Poly-cone Volume sT.
    sT = new TGeoTubeSeg("ITSssdsddMountingBracketT",ksupPRmin,ksupPRmax,
                         0.5*ksupPZ,ksupPPhi1,ksupPPhi2);
    //
    TGeoRotation *rotZ225 =new TGeoRotation("ITSssdConeZ225", 0.0,0.0, 22.5);
    rotZ225->RegisterYourself();
    TGeoRotation *rotZ675 =new TGeoRotation("ITSssdConeZ675", 0.0,0.0, 67.5);
    rotZ675->RegisterYourself();
    TGeoRotation *rotZ90  =new TGeoRotation("ITSssdConeZ90",  0.0,0.0, 90.0);
    rotZ90->RegisterYourself();
    TGeoRotation *rotZ1125=new TGeoRotation("ITSssdConeZ1125",0.0,0.0,112.5);
    rotZ1125->RegisterYourself();
    TGeoRotation *rotZ1575=new TGeoRotation("ITSssdConeZ1575",0.0,0.0,157.5);
    rotZ1575->RegisterYourself();
    TGeoRotation *rotZ180 =new TGeoRotation("ITSssdConeZ180", 0.0,0.0,180.0);
    rotZ180->RegisterYourself();
    TGeoRotation *rotZ2025=new TGeoRotation("ITSssdConeZ2025",0.0,0.0,202.5);
    rotZ2025->RegisterYourself();
    TGeoRotation *rotZ2475=new TGeoRotation("ITSssdConeZ2475",0.0,0.0,247.5);
    rotZ2475->RegisterYourself();
    TGeoRotation *rotZ270 =new TGeoRotation("ITSssdConeZ270", 0.0,0.0,270.0);
    rotZ270->RegisterYourself();
    TGeoRotation *rotZ2925=new TGeoRotation("ITSssdConeZ2925",0.0,0.0,292.5);
    rotZ2925->RegisterYourself();
    TGeoRotation *rotZ3375=new TGeoRotation("ITSssdConeZ3375",0.0,0.0,337.5);
    rotZ3375->RegisterYourself();
    //
    vl[0] = 0.0;vl[1] = kconCoolHoleRmin+0.5*kconCoolHoleHight;vl[2] = 0.0;
    rotZ225->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA225  = new TGeoCombiTrans("ITSssdConeTZ225",vg[0],
                                                     vg[1],vg[2],rotZ225);
    rotranA225->RegisterYourself();
    rotZ675->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA675  = new TGeoCombiTrans("ITSssdConeTZ675", vg[0],
                                                     vg[1],vg[2],rotZ675);
    rotranA675->RegisterYourself();
    rotZ1125->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA1125 = new TGeoCombiTrans("ITSssdConeTZ1125",vg[0],
                                                     vg[1],vg[2],rotZ1125);
    rotranA1125->RegisterYourself();
    rotZ1575->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA1575 = new TGeoCombiTrans("ITSssdConeTZ1575",vg[0],
                                                     vg[1],vg[2],rotZ1575);
    rotranA1575->RegisterYourself();
    rotZ2025->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA2025 = new TGeoCombiTrans("ITSssdConeTZ2025",vg[0],
                                                     vg[1],vg[2],rotZ2025);
    rotranA2025->RegisterYourself();
    rotZ2475->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA2475 = new TGeoCombiTrans("ITSssdConeTZ2475",vg[0],
                                                     vg[1],vg[2],rotZ2475);
    rotranA2475->RegisterYourself();
    rotZ2925->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA2925 = new TGeoCombiTrans("ITSssdConeTZ2925",vg[0],
                                                     vg[1],vg[2],rotZ2925);
    rotranA2925->RegisterYourself();
    rotZ3375->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA3375 = new TGeoCombiTrans("ITSssdConeTZ3375",vg[0],
                                                     vg[1],vg[2],rotZ3375);
    rotranA3375->RegisterYourself();
    TGeoRotation *rotZ30  = new TGeoRotation("ITSssdConeZ30", 0.0,0.0, 30.0);
    TGeoRotation *rotZ60  = new TGeoRotation("ITSssdConeZ60", 0.0,0.0, 60.0);
    //TGeoRotation *rotZ120 = new TGeoRotation("ITSssdConeZ120",0.0,0.0,120.0);
    TGeoRotation *rotZ150 = new TGeoRotation("ITSssdConeZ150",0.0,0.0,150.0);
    TGeoRotation *rotZ210 = new TGeoRotation("ITSssdConeZ210",0.0,0.0,210.0);
    //TGeoRotation *rotZ240 = new TGeoRotation("ITSssdConeZ240",0.0,0.0,240.0);
    TGeoRotation *rotZ300 = new TGeoRotation("ITSssdConeZ300",0.0,0.0,300.0);
    TGeoRotation *rotZ330 = new TGeoRotation("ITSssdConeZ330",0.0,0.0,330.0);
    vl[0] = kconMountHoleRmin+0.5*kconMountHoleHight; vl[1] = 0.0; vl[2] = 0.0;
    for(i=0;i<sF0->GetNz();i++) vl[2] += sF0->GetZ(i);
    vl[2] /= (Double_t)(sF0->GetNz());
    rotZ30->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA30 = new TGeoCombiTrans("ITSssdConeTZ30",vg[0],
                                                      vg[1],vg[2],rotZ30);
    rotranA30->RegisterYourself();
    rotZ90->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA90  = new TGeoCombiTrans("ITSssdConeTZ90", vg[0],
                                                     vg[1],vg[2],rotZ90);
    rotranA90->RegisterYourself();
    rotZ150->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA150 = new TGeoCombiTrans("ITSssdConeTZ150",vg[0],
                                                     vg[1],vg[2],rotZ150);
    rotranA150->RegisterYourself();
    rotZ210->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA210 = new TGeoCombiTrans("ITSssdConeTZ210",vg[0],
                                                     vg[1],vg[2],rotZ210);
    rotranA210->RegisterYourself();
    rotZ270->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA270 = new TGeoCombiTrans("ITSssdConeTZ270",vg[0],
                                                     vg[1],vg[2],rotZ270);
    rotranA270->RegisterYourself();
    rotZ330->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranA330 = new TGeoCombiTrans("ITSssdConeTZ330",vg[0],
                                                     vg[1],vg[2],rotZ330);
    rotranA330->RegisterYourself();
    vl[0] = 0.0; vl[1] = 0.0; vl[2] = sA0->GetZ(10)+sT->GetDz();
    rotZ60->LocalToMaster(vl,vg);
    TGeoCombiTrans *rotranBrTZ60  = new TGeoCombiTrans("ITSssdConeBrTZ60",
                                                  vg[0],vg[1],vg[2],rotZ60);
    rotranBrTZ60->RegisterYourself();
    TGeoCombiTrans *rotranBrTZ180 = new TGeoCombiTrans("ITSssdConeBrTZ180",
                                                  vg[0],vg[1],vg[2],rotZ180);
    rotranBrTZ180->RegisterYourself();
    TGeoCombiTrans *rotranBrTZ300 = new TGeoCombiTrans("ITSssdConeBrTZ300",
                                                  vg[0],vg[1],vg[2],rotZ300);
    rotranBrTZ300->RegisterYourself();
    if(GetDebug(1)){
        rotZ225->Print();
        rotZ675->Print();
        rotZ90->Print();
        rotZ1125->Print();
        rotZ1575->Print();
        rotZ180->Print();
        rotZ2025->Print();
        rotZ2475->Print();
        rotZ270->Print();
        rotZ2925->Print();
        rotZ3375->Print();
        rotranA225->Print();
        rotranA675->Print();
        rotranA1125->Print();
        rotranA1575->Print();
        rotranA2025->Print();
        rotranA2475->Print();
        rotranA2925->Print();
        rotranA3375->Print();
        rotZ60->Print();
        rotZ300->Print();
        rotranA30->Print();
        rotranA90->Print();
        rotranA150->Print();
        rotranA210->Print();
        rotranA270->Print();
        rotranA330->Print();
        rotranBrTZ60->Print();
        rotranBrTZ180->Print();
        rotranBrTZ300->Print();
    } // end if GetDebug(1)
    sA = new TGeoCompositeShape("ITSssdSuportConeCarbonFiberSurfaceA",
        "(((((((((((((((((((((((((((("
        "ITSssdSuportConeCarbonFiberSurfaceA0 +"
        "ITSssdWingCarbonFiberSurfaceG) +"
        "ITSssdWingCarbonFiberSurfaceG:ITSssdConeZ90) +"
        "ITSssdWingCarbonFiberSurfaceG:ITSssdConeZ180) +"
        "ITSssdWingCarbonFiberSurfaceG:ITSssdConeZ270) -"
        "(ITSssdCableHoleAh1:ITSssdConeZ225*ITSssdCableHoleAh2:ITSssdConeZ225)) -"
        "(ITSssdCableHoleAh1:ITSssdConeZ675*ITSssdCableHoleAh2:ITSssdConeZ675)) -"
        "(ITSssdCableHoleAh1:ITSssdConeZ1125*ITSssdCableHoleAh2:ITSssdConeZ1125)) -"
        "(ITSssdCableHoleAh1:ITSssdConeZ1575*ITSssdCableHoleAh2:ITSssdConeZ1575)) -"
        "(ITSssdCableHoleAh1:ITSssdConeZ2025*ITSssdCableHoleAh2:ITSssdConeZ2025)) -"
        "(ITSssdCableHoleAh1:ITSssdConeZ2475*ITSssdCableHoleAh2:ITSssdConeZ2475)) -"
        "(ITSssdCableHoleAh1:ITSssdConeZ2925*ITSssdCableHoleAh2:ITSssdConeZ2925)) -"
        "(ITSssdCableHoleAh1:ITSssdConeZ3375*ITSssdCableHoleAh2:ITSssdConeZ3375)) -"
        "ITSssdCoolingHoleAh3:ITSssdConeTZ225) -"
        "ITSssdCoolingHoleAh3:ITSssdConeTZ675) -"
        "ITSssdCoolingHoleAh3:ITSssdConeTZ1125) -"
        "ITSssdCoolingHoleAh3:ITSssdConeTZ1575) -"
        "ITSssdCoolingHoleAh3:ITSssdConeTZ2025) -"
        "ITSssdCoolingHoleAh3:ITSssdConeTZ2475) -"
        "ITSssdCoolingHoleAh3:ITSssdConeTZ2925) -"
        "ITSssdCoolingHoleAh3:ITSssdConeTZ3375) -"
        "ITSssdMountingPostHoleAh4:ITSssdConeTZ30) -"
        "ITSssdMountingPostHoleAh4:ITSssdConeTZ90) -"
        "ITSssdMountingPostHoleAh4:ITSssdConeTZ150) -"
        "ITSssdMountingPostHoleAh4:ITSssdConeTZ210) -"
        "ITSssdMountingPostHoleAh4:ITSssdConeTZ270) -"
        "ITSssdMountingPostHoleAh4:ITSssdConeTZ330) -"
        "ITSssdsddMountingBracketT:ITSssdConeBrTZ60) -"
        "ITSssdsddMountingBracketT:ITSssdConeBrTZ180) -"
        "ITSssdsddMountingBracketT:ITSssdConeBrTZ300"
        );
    sB = new TGeoCompositeShape("ITSssdSuportConeStaseliteB",
        "(((((((((((((((((((((((((((("
        "ITSssdSuportConeStaseliteB0 +"
        "ITSssdWingStaseliteH) +"
        "ITSssdWingStaseliteH:ITSssdConeZ90) +"
        "ITSssdWingStaseliteH:ITSssdConeZ180) +"
        "ITSssdWingStaseliteH:ITSssdConeZ270) -"
        "(ITSssdCableHoleBh1:ITSssdConeZ225*ITSssdCableHoleBh2:ITSssdConeZ225)) -"
        "(ITSssdCableHoleBh1:ITSssdConeZ675*ITSssdCableHoleBh2:ITSssdConeZ675)) -"
        "(ITSssdCableHoleBh1:ITSssdConeZ1125*ITSssdCableHoleBh2:ITSssdConeZ1125)) -"
        "(ITSssdCableHoleBh1:ITSssdConeZ1575*ITSssdCableHoleBh2:ITSssdConeZ1575)) -"
        "(ITSssdCableHoleBh1:ITSssdConeZ2025*ITSssdCableHoleBh2:ITSssdConeZ2025)) -"
        "(ITSssdCableHoleBh1:ITSssdConeZ2475*ITSssdCableHoleBh2:ITSssdConeZ2475)) -"
        "(ITSssdCableHoleBh1:ITSssdConeZ2925*ITSssdCableHoleBh2:ITSssdConeZ2925)) -"
        "(ITSssdCableHoleBh1:ITSssdConeZ3375*ITSssdCableHoleBh2:ITSssdConeZ3375)) -"
        "ITSssdCoolingHoleBh3:ITSssdConeTZ225) -"
        "ITSssdCoolingHoleBh3:ITSssdConeTZ675) -"
        "ITSssdCoolingHoleBh3:ITSssdConeTZ1125) -"
        "ITSssdCoolingHoleBh3:ITSssdConeTZ1575) -"
        "ITSssdCoolingHoleBh3:ITSssdConeTZ2025) -"
        "ITSssdCoolingHoleBh3:ITSssdConeTZ2475) -"
        "ITSssdCoolingHoleBh3:ITSssdConeTZ2925) -"
        "ITSssdCoolingHoleBh3:ITSssdConeTZ3375) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ30) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ90) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ150) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ210) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ270) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ330) -"
        "ITSssdsddMountingBracketT:ITSssdConeBrTZ60) -"
        "ITSssdsddMountingBracketT:ITSssdConeBrTZ180) -"
        "ITSssdsddMountingBracketT:ITSssdConeBrTZ300"
        );
    sC = new TGeoCompositeShape("ITSssdSuportConeRohacellC",
      "((((((("
      "ITSssdSuportConeRohacellC0 -"
      "ITSssdCableHoleBh1:ITSssdConeZ225*ITSssdCableHoleBh2:ITSssdConeZ225) -"
      "ITSssdCableHoleBh1:ITSssdConeZ675*ITSssdCableHoleBh2:ITSssdConeZ675) -"
      "ITSssdCableHoleBh1:ITSssdConeZ1125*ITSssdCableHoleBh2:ITSssdConeZ1125) -"
      "ITSssdCableHoleBh1:ITSssdConeZ1575*ITSssdCableHoleBh2:ITSssdConeZ1575) -"
      "ITSssdCableHoleBh1:ITSssdConeZ2025*ITSssdCableHoleBh2:ITSssdConeZ2025) -"
      "ITSssdCableHoleBh1:ITSssdConeZ2475*ITSssdCableHoleBh2:ITSssdConeZ2475) -"
      "ITSssdCableHoleBh1:ITSssdConeZ2925*ITSssdCableHoleBh2:ITSssdConeZ2925) -"
      "ITSssdCableHoleBh1:ITSssdConeZ3375*ITSssdCableHoleBh2:ITSssdConeZ3375 "
        );
    sF = new TGeoCompositeShape("ITSssdSuportConeRohacellCF",
        "((((("
        "ITSssdSuportConeRohacellCF0 -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ30) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ90) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ150) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ210) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ270) -"
        "ITSssdMountingPostHoleBh4:ITSssdConeTZ330"
        );
    //
    // In volume SCB, th Inserto Stesalite 4411w material volume, there
    // are a number of Stainless steel screw and pin studs which will be
    // filled with screws/studs.
    sD = new TGeoTube("ITS Screw+stud used to mount things to the SSD "
                      "support cone",
                      0.0,kconRScrewM5by12,kconLScrewM5by12);
    sE = new TGeoTube("ITS pin used to mount things to the "
                      "SSD support cone",0.0,kconRPinO6,kconLPinO6);
    // Bolt heads holding the SSD-SDD tube to the SSD cone.
    // Bolt -- PolyCone
    //Poly-cone Volume sQ.
    sQ = new TGeoPcon("ITS SSD Thermal sheald M6 screw headQ",0.0,360.0,4);
    sQ->Z(0)    = sA0->GetZ(12);
    sQ->Rmin(0) = 0.0;
    sQ->Rmax(0) = kcylRM6;
    sQ->Z(1)    = sA0->GetZ(10) + kconZM6Head;
    sQ->Rmin(1) = 0.0;
    sQ->Rmax(1) = kcylRM6;
    sQ->Z(2)    = sQ->GetZ(1);
    sQ->Rmin(2) = 0.0;
    sQ->Rmax(2) = kconRM6Head;
    sQ->Z(3)    = sA0->GetZ(10)+ksupPZ;
    sQ->Rmin(3) = 0.0;
    sQ->Rmax(3) = kconRM6Head;
    // air infront of bolt (stasolit Volume K) -- Tube
    sR = new TGeoTube("ITS Air in front of bolt (in stasolit)R",
                      sQ->GetRmin(3),sQ->GetRmax(3),0.5*(ksupPZ-kconCthick));
    // air infront of bolt (carbon fiber volume I) -- Tube
    sS = new TGeoTube("ITS Air in front of Stainless Steal Screw end, M6S",
                      sQ->GetRmin(3),sQ->GetRmax(3),0.5*kconCthick);
    //
    if(GetDebug(1)){
        sA0->InspectShape();
        sB0->InspectShape();
        sC0->InspectShape();
        sF0->InspectShape();
        sQ->InspectShape();
        sAh1->InspectShape();
        sBh1->InspectShape();
        sAh2->InspectShape();
        sBh2->InspectShape();
        sAh3->InspectShape();
        sBh3->InspectShape();
        sAh4->InspectShape();
        sBh4->InspectShape();
        sG->InspectShape();
        sH->InspectShape();
        sT->InspectShape();
        sD->InspectShape();
        sE->InspectShape();
        sR->InspectShape();
        sS->InspectShape();
        sA->InspectShape();
        sB->InspectShape();
        sC->InspectShape();
        sF->InspectShape();
    } // end if GetDebug(1)
    TGeoVolume *vA,*vB,*vC,*vD,*vE,*vF,*vQ,*vR,*vS,*vT;
    //
    vA = new TGeoVolume("ITSssdConeA",sA,medSSDcf); // Carbon Fiber
    vA->SetVisibility(kTRUE);
    vA->SetLineColor(4); // blue
    vA->SetLineWidth(1);
    vA->SetFillColor(vA->GetLineColor());
    vA->SetFillStyle(4050); // 50% transparent
    vB = new TGeoVolume("ITSssdConeB",sB,medSSDfs); // Staselite
    vB->SetVisibility(kTRUE);
    vB->SetLineColor(2); // red
    vB->SetLineWidth(1);
    vB->SetFillColor(vB->GetLineColor());
    vB->SetFillStyle(4050); // 50% transparent
    vC = new TGeoVolume("ITSssdConeC",sC,medSSDfo); // Rohacell
    vC->SetVisibility(kTRUE);
    vC->SetLineColor(3); // green
    vC->SetLineWidth(1);
    vC->SetFillColor(vC->GetLineColor());
    vC->SetFillStyle(4050); // 50% transparent
    vF = new TGeoVolume("ITSssdConeF",sF,medSSDfo); // Rohacell;
    vF->SetVisibility(kTRUE);
    vF->SetLineColor(3); // green
    vF->SetLineWidth(1);
    vF->SetFillColor(vF->GetLineColor());
    vF->SetFillStyle(4050); // 50% transparent
    vD = new TGeoVolume("ITSssdConeD",sD,medSSDss);
    vD->SetVisibility(kTRUE);
    vD->SetLineColor(1); // black
    vD->SetLineWidth(1);
    vD->SetFillColor(vD->GetLineColor());
    vD->SetFillStyle(4000); // 0% transparent
    vE = new TGeoVolume("ITSssdConeE",sE,medSSDss);
    vE->SetVisibility(kTRUE);
    vE->SetLineColor(1); // black
    vE->SetLineWidth(1);
    vE->SetFillColor(vE->GetLineColor());
    vE->SetFillStyle(4000); // 0% transparent
    vQ = new TGeoVolume("ITSssdConeQ",sQ,medSSDss);
    vQ->SetVisibility(kTRUE);
    vQ->SetLineColor(1); // black
    vQ->SetLineWidth(1);
    vQ->SetFillColor(vQ->GetLineColor());
    vQ->SetFillStyle(4000); // 0% transparent
    vR = new TGeoVolume("ITSssdConeR",sR,medSSDair);
    vR->SetVisibility(kTRUE);
    vR->SetLineColor(5); // yellow
    vR->SetLineWidth(1);
    vR->SetFillColor(vR->GetLineColor());
    vR->SetFillStyle(4090); // 90% transparent
    vS = new TGeoVolume("ITSssdConeS",sS,medSSDair);
    vS->SetVisibility(kTRUE);
    vS->SetLineColor(5); // yellow
    vS->SetLineWidth(1);
    vS->SetFillColor(vS->GetLineColor());
    vS->SetFillStyle(4090); // 90% transparent
    vT = new TGeoVolume("ITSssdsddMountingBracket",sT,medSSDal);
    vT->SetVisibility(kTRUE);
    vT->SetLineColor(5); // yellow
    vT->SetLineWidth(1);
    vT->SetFillColor(vT->GetLineColor());
    vT->SetFillStyle(4000); // 0% transparent
    //
    TGeoCombiTrans *rotran;
    TGeoTranslation *tran;
    tran = new TGeoTranslation("ITSssdConeTrans",0.0,0.0,-kconZDisplacement);
    TGeoRotation *rotY180 = new TGeoRotation("",0.0,180.0,0.0);
    TGeoCombiTrans *flip  = new TGeoCombiTrans("ITSssdConeFlip",
                                           0.0,0.0,kconZDisplacement,rotY180);
    //delete rotY180;// rot not explicity used in AddNode functions.
    //
    //
    //
    //
    vA->AddNode(vB,1,0);
    vB->AddNode(vC,1,0);
    vB->AddNode(vF,1,0);
    moth->AddNode(vA,1,tran); // RB24 side
    moth->AddNode(vA,2,flip); // RB26 side (Absorber)
    //
    //
    //
    // Insert Bolt and Pins in both the Cone and Cylinder at the same time.
    Int_t nCopyCDv=0,nCopyCEv=0,nCopyQv=0,nCopyvR=0,nCopySv=0,nCopyTv=0;
    Int_t nCopyvD=0,nCopyvE=0;
    z = sCB->GetZ(0)+sCD->GetDz(); // sCB->GetZ(0)<0!
    dt = (360.0/((Double_t)kcylNPin));
    for(i=0;i<kcylNPin;i++){
        t = ((Double_t)i)*dt;
        x = kcylRholes*CosD(t+kcylPhi0Pin);
        y = kcylRholes*SinD(t+kcylPhi0Pin);
        tran = new TGeoTranslation("",x,y,z);
        vCB->AddNode(vCD,++nCopyCDv,tran);
        tran = new TGeoTranslation("",x,y,-z);
        vCB->AddNode(vCD,++nCopyCDv,tran);
    } // end for i
    dt = (360.0/((Double_t)kcylNM6));
    for(i=0;i<kcylNM6;i++){
        t = ((Double_t)i)*dt;
        x = kcylRholes*CosD(t+kcylPhi0M6);
        y = kcylRholes*SinD(t+kcylPhi0M6);
        z = sCB->GetZ(0)+sCE->GetDz(); // sCB->GetZ()<0!
        tran = new TGeoTranslation("",x,y,z);
        vCB->AddNode(vCE,++nCopyCEv,tran);
        tran = new TGeoTranslation("",x,y,-z);
        vCB->AddNode(vCE,++nCopyCEv,tran);
        tran = new TGeoTranslation("",x,y,0.0);
        vB->AddNode(vQ,++nCopyQv,tran);
        if(!((t<rotranBrTZ60->GetRotation()->GetPhiRotation()+sT->GetPhi2()&&
             t>rotranBrTZ60->GetRotation()->GetPhiRotation()-sT->GetPhi1())||
            (t<rotranBrTZ180->GetRotation()->GetPhiRotation()+sT->GetPhi2()&&
             t>rotranBrTZ180->GetRotation()->GetPhiRotation()-sT->GetPhi1())||
            (t<rotranBrTZ300->GetRotation()->GetPhiRotation()+sT->GetPhi2()&&
             t>rotranBrTZ300->GetRotation()->GetPhiRotation()-sT->GetPhi1()))){
            // If not at an angle where the bracket sT is located.
            tran = new TGeoTranslation("",x,y,sB0->GetZ(10)-sR->GetDz());
            vB->AddNode(vR,++nCopyvR,tran);
            tran = new TGeoTranslation("",x,y,sA0->GetZ(10)-sS->GetDz());
            vA->AddNode(vS,++nCopySv,tran);
        } // end if
    } // end for i
    // Add the mounting brackets to the RB24 side only.
    vl[0] = 0.0;
    vl[1] = 0.0;
    vl[2] = sA0->GetZ(10)+kconZDisplacement-sT->GetDz();
    rotZ60->LocalToMaster(vl,vg);
    rotran = new TGeoCombiTrans("",vg[0],vg[1],vg[2],rotZ60);
    moth->AddNode(vT,++nCopyTv,rotran);
    rotZ180->LocalToMaster(vl,vg);
    rotran = new TGeoCombiTrans("",vg[0],vg[1],vg[2],rotZ180);
    moth->AddNode(vT,++nCopyTv,rotran);
    rotZ300->LocalToMaster(vl,vg);
    rotran = new TGeoCombiTrans("",vg[0],vg[1],vg[2],rotZ300);
    moth->AddNode(vT,++nCopyTv,rotran);
    //
    Double_t da[] = {-3.5,-1.5,1.5,3.5};
    for(i=0;i<2;i++){ // Mounting for ITS-TPC bracket or ITS-Rails
        t0 = 180.*((Double_t)i);
        for(j=-kconNScrewM5by12/2;j<=kconNScrewM5by12/2;j++)if(j!=0){
                    //screws per ITS-TPC brkt
            t = t0 + 5.0*((Double_t)j);
            tran = new TGeoTranslation("",kconROutHoles*CosD(t),
                                          kconROutHoles*SinD(t),
                                          sB0->GetZ(0)+sE->GetDz());
            vB->AddNode(vE,++nCopyvE,tran);
        } // end or j
        for(j=-kconNPinO6/2;j<=kconNPinO6/2;j++){ // pins per ITS-TPC bracket
            t = t0 + 3.0*((Double_t)j);
            tran = new TGeoTranslation("",kconROutHoles*CosD(t),
                                          kconROutHoles*SinD(t),
                                          sB0->GetZ(0)+sD->GetDz());
            vB->AddNode(vD,++nCopyvD,tran);
        } // end or j
        t0 = (-5.5+191.*((Double_t)i));
        for(j=0;j<kconNRailScrews;j++){ // screws per ITS-rail bracket
            t = t0+da[j];
            tran = new TGeoTranslation("",kconROutHoles*CosD(t),
                                          kconROutHoles*SinD(t),
                                          sB0->GetZ(0)+sE->GetDz());
            vB->AddNode(vE,++nCopyvE,tran);
        } // end or j
        t0 = (95.5+191.*((Double_t)i));
        for(j=-kconNRailPins/2;j<=kconNRailPins/2;j++)if(j!=0){ 
             // pins per ITS-rail bracket
            t = t0+(5.5*((Double_t)j));
            tran = new TGeoTranslation("",kconROutHoles*CosD(t),
                                          kconROutHoles*SinD(t),
                                          sB0->GetZ(0)+sD->GetDz());
            vB->AddNode(vD,++nCopyvD,tran);
        } // end or j
    } // end for i
    for(i=0;i<kconNmounts;i++){ 
                // mounting points for SPD-cone+Beam-pipe support
        t0 = (45.0+((Double_t)i)*360./((Double_t)kconNmounts));
        for(j=-1;j<=1;j++)if(j!=0){ // 2 screws per bracket
            t = t0+((Double_t)j)*0.5*kconMountPhi0;
            tran = new TGeoTranslation("",kconROutHoles*CosD(t),
                                          kconROutHoles*SinD(t),
                                          sB0->GetZ(0)+sD->GetDz());
            vB->AddNode(vD,++nCopyvD,tran);
        } // end for j
        for(j=0;j<1;j++){ // 1 pin per bracket
            t = t0;
            tran = new TGeoTranslation("",kconROutHoles*CosD(t),
                                          kconROutHoles*SinD(t),
                                          sB0->GetZ(0)+sD->GetDz());
            vB->AddNode(vE,++nCopyvE,tran);
        } // end for j
    } // end for i
    if(GetDebug(1)){
        vA->PrintNodes();
        vB->PrintNodes();
        vC->PrintNodes();
        vD->PrintNodes();
        vE->PrintNodes();
        vF->PrintNodes();
        vQ->PrintNodes();
        vR->PrintNodes();
        vS->PrintNodes();
        vT->PrintNodes();
    } // end if
}

//______________________________________________________________________
void AliITSv11GeometrySupport::ServicesCableSupport(TGeoVolume *moth,
                                                    TGeoManager *mgr){
    // Define the detail ITS cable support trays on both the RB24 and 
    // RB26 sides..
    // Inputs:
    //   TGeoVolume  *moth  The mother volume to place this object.
    //   TGeoManager *mgr   A pointer to the Geo-Manager default gGeoManager
    // Outputs:
    //  none.
    // Return:
    //  none.
    // Based on the Drawings SSup_201A.jpg unless otherwise stated, 
    // Volumes A..., 
    TGeoMedium *medSUPcf    = 0; // SUP support cone Carbon Fiber materal nbr.
    TGeoMedium *medSUPfs    = 0; // SUP support cone inserto stesalite 4411w.
    TGeoMedium *medSUPfo    = 0; // SUP support cone foam, Rohacell 50A.
    TGeoMedium *medSUPss    = 0; // SUP support cone screw material,Stainless
    TGeoMedium *medSUPair   = 0; // SUP support cone Air
    TGeoMedium *medSUPal    = 0; // SUP support cone SDD mounting bracket Al
    TGeoMedium *medSUPwater = 0; // SUP support cone Water
    medSUPcf    = mgr->GetMedium("ITSssdCarbonFiber");
    medSUPfs    = mgr->GetMedium("ITSssdStaselite4411w");
    medSUPfo    = mgr->GetMedium("ITSssdRohacell50A");
    medSUPss    = mgr->GetMedium("ITSssdStainlessSteal");
    medSUPair   = mgr->GetMedium("ITSssdAir");
    medSUPal    = mgr->GetMedium("ITSssdAl");
    medSUPwater = mgr->GetMedium("ITSssdWater");
    //
    Int_t i,j,iRmin;
    Double_t x,y,z,t,t0,dt,di,r,l,local[3],master[3];
    Char_t name[100];
    Double_t r1,r2,m;
    // RB 24, Open Side.
    const Double_t kfrm24Z0           = 900*fgkmm;//SSup_203A.jpg
    const Double_t kfrm24Thss         = 5.0*fgkmm;
    const Double_t kfrm24Rss          = 444.5*fgkmm-kfrm24Thss; //SSup_204A.jpg
    const Double_t kfrm24Width        = 10.0*fgkmm;
    const Double_t kfrm24Hight        = 10.0*fgkmm;
    const Double_t kfrm24Phi0         = 15.2*fgkDegree; // SSup_602A.jpg
    const Double_t kfrm24Phi1         = (90.0-7.6)*fgkDegree; // SSup_802A.jpg
    const Double_t kfrm24ZssSection   = (415.0-10.0)*fgkmm;
    const Int_t    kfrm24NZsections   = 4;
    const Int_t    kfrm24NPhiSections = 4;
    const Int_t    kfrm24NPhi         = 4;
    // These numbers are guessed at.
    const Double_t kfrm24ZfracAngle   =  0.55; // frational z length to brack
    const Double_t kfrm24Angle        =  10.0*fgkDegree; // Guessed at
    //
    TGeoTubeSeg *sA24[kfrm24NZsections+1];
    TGeoArb8    *sB24[kfrm24NZsections+1];
    Double_t zA24[kfrm24NZsections+1];
    l = 4.*kfrm24ZssSection+5*kfrm24Width;
    j = iRmin = 0;
    for(i=0;i<kfrm24NZsections+1;i++){
        sprintf(name,"ITS sup Cable tray support frame radial section A24[%d]",
                i);
        r1 = kfrm24Rss;
        if(i==0) zA24[i] = kfrm24Width;
        else zA24[i] = zA24[i-1] + kfrm24ZssSection + kfrm24Width;
        if(zA24[i]>l*kfrm24ZfracAngle){ // break, radii get larger
            r1 = kfrm24Rss + (zA24[i]-kfrm24ZfracAngle*l)*SinD(kfrm24Angle);
        } // end if
        r2 = r1+kfrm24Thss;
        sA24[i] = new TGeoTubeSeg(name,r1,r2,0.5*kfrm24Width,kfrm24Phi0,
                                  kfrm24Phi1);
        if(i>0)if(sA24[i-1]->GetRmin()==sA24[i]->GetRmin()) j = iRmin = i;
    } // end for i
    for(i=0;i<kfrm24NZsections;i++){
        sprintf(name,"ITS sup Cable tray support frame Z section B24[%d]",i);
        sB24[i] = new TGeoArb8(name,0.5*kfrm24ZssSection);
        sB24[i]->SetVertex(0,sA24[i]->GetRmin(),0.5*kfrm24Hight);
        sB24[i]->SetVertex(1,sA24[i]->GetRmax(),0.5*kfrm24Hight);
        sB24[i]->SetVertex(2,sA24[i]->GetRmin(),-0.5*kfrm24Hight);
        sB24[i]->SetVertex(3,sA24[i]->GetRmax(),-0.5*kfrm24Hight);
        sB24[i]->SetVertex(4,sA24[i+1]->GetRmin(),0.5*kfrm24Hight);
        sB24[i]->SetVertex(5,sA24[i+1]->GetRmax(),0.5*kfrm24Hight);
        sB24[i]->SetVertex(6,sA24[i+1]->GetRmin(),-0.5*kfrm24Hight);
        sB24[i]->SetVertex(7,sA24[i+1]->GetRmax(),-0.5*kfrm24Hight);
    } // end for i
    if(GetDebug(1)){
        for(i=0;i<kfrm24NZsections+1;i++) sA24[i]->InspectShape();
        for(i=0;i<kfrm24NZsections;i++)   sB24[i]->InspectShape();
    } // end if GetDebug(1)
    TGeoVolume *vA24[kfrm24NZsections+1],*vB24[kfrm24NZsections];
    TGeoVolumeAssembly *vM24;
    TGeoTranslation *tran;
    TGeoRotation    *rot,*rot1;
    TGeoCombiTrans  *tranrot;
    //
    for(i=0;i<kfrm24NZsections+1;i++){
        vA24[i] = 0;
        sprintf(name,"ITSsupFrameA24[%d]",i);
        vA24[i] = new TGeoVolume(name,sA24[i],medSUPss);
        vA24[i]->SetVisibility(kTRUE);
        vA24[i]->SetLineColor(1); // black
        vA24[i]->SetLineWidth(1);
        vA24[i]->SetFillColor(vA24[i]->GetLineColor());
        vA24[i]->SetFillStyle(4000); // 0% transparent
    } // end for i
    for(i=0;i<kfrm24NZsections;i++){
        vB24[i] = 0;
        sprintf(name,"ITSsupFrameB24[%d]",i);
        vB24[i] = new TGeoVolume(name,sB24[i],medSUPss);
        vB24[i]->SetVisibility(kTRUE);
        vB24[i]->SetLineColor(1); // black
        vB24[i]->SetLineWidth(1);
        vB24[i]->SetFillColor(vB24[i]->GetLineColor());
        vB24[i]->SetFillStyle(4000); // 0% transparent
    } // end for i
    vM24 = new TGeoVolumeAssembly("ITSsupFrameM24");
    //vM24->SetVisibility(kTRUE);
    //vM24->SetLineColor(7); // light blue
    //vM24->SetLineWidth(1);
    //vM24->SetFillColor(vM24->GetLineColor());
    //vM24->SetFillStyle(4090); // 90% transparent
    //
    Int_t ncopyB24[kfrm24NPhiSections];
    t0 = kfrm24Phi0;
    dt = (kfrm24Phi1-kfrm24Phi0)/((Double_t)kfrm24NPhiSections);
    for(i=0;i<=kfrm24NZsections;i++){
        z = zA24[i];
        tran = new TGeoTranslation("",0.0,0.0,z);
        vM24->AddNode(vA24[i],1,tran);
       if(i<kfrm24NZsections){
           ncopyB24[i] = 1;
           for(j=0;j<=kfrm24NPhiSections;j++){
               t = t0 + ((Double_t)j)*dt;
               rot = new TGeoRotation("",0.0,0.0,t);
               tranrot = new TGeoCombiTrans("",0.0,0.0,z+sB24[i]->GetDz(),rot);
               //delete rot;// rot not explicity used in AddNode functions.
               vM24->AddNode(vB24[i],ncopyB24[i]++,tranrot);
           } // end for j
       } // end if
    } // end for i
    tran = new TGeoTranslation("",0.0,0.0,kfrm24Z0);
    moth->AddNode(vM24,1,tran);
    for(i=1;i<kfrm24NPhi;i++){
        di = (Double_t) i;
        rot = new TGeoRotation("",0.0,0.0,90.0*di);
        tranrot = new TGeoCombiTrans("",0.0,0.0,kfrm24Z0,rot);
        //delete rot;// rot not explicity used in AddNode functions.
        moth->AddNode(vM24,i+1,tranrot);
    } // end for i
    if(GetDebug(1)){
        for(i=0;i<kfrm24NZsections+1;i++) vA24[i]->PrintNodes();
        for(i=0;i<kfrm24NZsections;i++) vB24[i]->PrintNodes();
        vM24->PrintNodes();
    } // end if
    //==================================================================
    // RB24 Cable Tray
    const Double_t kct24WidthBottom   = 44.0*fgkmm; // Serv-C_208.jpg
    const Double_t kct24WidthTop      = 46.0*fgkmm; // Serv-C_208.jpg
    const Double_t kct24Hight         = 51.0*fgkmm; // Serv-C_208.jpg
    const Double_t kct24AlThick       = 1.0*fgkmm; // Serv-C_208.jpg
    const Double_t kct24CapWidth      = 46.0*fgkmm; // Serv-C_208.jpg
    const Double_t kct24CapEar        = 5.0*fgkmm; // Guess
    const Double_t kct24Rmin          = 455.0*fgkmm; // Serv-C_203.jpg
    const Double_t kct24CoolSectionH  = 470.0*fgkmm-kct24Rmin;// Serv-C_203.jpg
    const Double_t kct24CoolCableDivEar = 2.0*fgkmm; // Guess
    const Int_t kct24Ntrays           = 48; // Serv-C_205.jpg
    //const Int_t kct24Ntubes           = 3; // Serv-C_208.jpg
    // Patch Pannels for RB 24 side
    const Double_t kft24PPHightSPDFMD = 72.0*fgkmm; // Serv-C_SPD/FMD.jpg
    const Double_t kft24PPHightSDDSSD = 104.0*fgkmm; // Serv-C_SDD/SSD.jpg
    const Double_t kft24PPlength      = 350.0*fgkmm;//Serv-C_SPD/SDD/SSD/FMD_1.jpg
    const Double_t kft24Theta         = 2.0*TMath::ATan2(kct24WidthBottom,
                                                 2.0*kct24Rmin)*fgkRadian; //
    const Int_t    kft24NPatchPannels = 20; //
    //
    Double_t xp[12],yp[12];
    TGeoPcon *sMT24;
    TGeoXtru *sT24,*sTs24,*sTl24,*sTt24,*sU24,*sVl24,*sVs24,*sW24;
    TGeoXtru *s3PP24,*s2PP24,*sV3PP24,*sV2PP24;
    // Outer Tray Full
    sT24 = new TGeoXtru(3);
    sT24->SetName("ITS sup Full Cable Tray for RB24 Side T24");
    xp[0]  = -0.5*kct24WidthBottom;
    yp[0]  = sA24[0]->GetRmax();
    yp[1]  = yp[0] + kct24Hight-kct24CapEar;
    xp[1]  = Xfrom2Points(xp[0],yp[0],-0.5*kct24WidthTop+kct24AlThick,
                          yp[0]+kct24Hight,yp[1]);
    yp[2]  = yp[1];
    xp[2]  = xp[1]-kct24AlThick;
    xp[3]  = -0.5*kct24CapWidth;
    yp[3]  = yp[0] + kct24Hight;
    xp[4]  = -xp[3];
    yp[4]  =  yp[3];
    xp[5]  = -xp[2];
    yp[5]  =  yp[2];
    xp[6]  = -xp[1];
    yp[6]  =  yp[1];
    xp[7]  = -xp[0];
    yp[7]  =  yp[0];
    sT24->DefinePolygon(8,xp,yp);
    sT24->DefineSection(0,zA24[0]-kfrm24Width,0.0,0.0,1.0);
    sT24->DefineSection(1,zA24[iRmin],0.0,0.0,1.0);
    sT24->DefineSection(2,zA24[kfrm24NZsections]+kfrm24Width,0.0,
                      sA24[kfrm24NZsections]->GetRmax()-sA24[0]->GetRmin());
    // RB 24 full tray no divider (for ALG and T0-V0 cables?)
    sW24 = new TGeoXtru(3);
    sW24->SetName("ITS sup Cable Tray No Divider for RB24 Side W24");
    xp[0] = sT24->GetX(0) + kct24AlThick;
    yp[0] = sT24->GetY(0) + kct24AlThick;
    yp[1] = sT24->GetY(3) - kct24AlThick;
    xp[1] = Xfrom2Points(sT24->GetX(0),sT24->GetY(0),sT24->GetX(1),
                         sT24->GetY(1),yp[1]) + kct24AlThick;
    xp[2] = -xp[1];
    yp[2] =  yp[1];
    xp[3] = -xp[0];
    yp[3] =  yp[0];
    sW24->DefinePolygon(4,xp,yp);
    for(i=0;i<sT24->GetNz();i++){
        sW24->DefineSection(i,sT24->GetZ(i),sT24->GetXOffset(i),
                            sT24->GetYOffset(i),sT24->GetScale(i));
    } // end for i
    // Outer Tray Short
    sTs24 = new TGeoXtru(3);
    sTs24->SetName("ITS sup Short Cable Tray for RB24 Side Ts24");
    yp[0]  = sT24->GetY(0) + kct24CoolSectionH;
    xp[0]  = Xfrom2Points(sT24->GetX(0),sT24->GetY(0),sT24->GetX(1),
                         sT24->GetY(1),yp[0]);
    for(i=1;i<7;i++){
        xp[i]  = sT24->GetX(i);
        yp[i]  = sT24->GetY(i);
    } // end for i
    xp[7]  = -xp[0];
    yp[7]  =  yp[0];
    sTs24->DefinePolygon(8,xp,yp);
    sTs24->DefineSection(0,zA24[0] -kfrm24Width+kft24PPlength);
    sTs24->DefineSection(1,zA24[iRmin]);
    sTs24->DefineSection(2,zA24[kfrm24NZsections]+kfrm24Width,
                         sT24->GetXOffset(2),
                         sT24->GetYOffset(2),sT24->GetScale(2));
    // Outer Tray Long
    sTl24 = new TGeoXtru(3);
    sTl24->SetName("ITS sup Long Cable Tray for RB24 Side Tl24");
    for(i=0;i<8;i++){
    xp[i]  = sTs24->GetX(i);
    yp[i]  = sTs24->GetY(i);
    } // End for i
    sTl24->DefinePolygon(8,xp,yp);
    sTl24->DefineSection(0,zA24[0]-kfrm24Width,0.0,0.0,1.0);
    sTl24->DefineSection(1,zA24[iRmin],0.0,0.0,1.0);
    sTl24->DefineSection(2,zA24[kfrm24NZsections]+kfrm24Width,0.0,
                     sA24[kfrm24NZsections]->GetRmax()-sA24[0]->GetRmin(),1.0);
    // Outer Tray for air Tubes
    sTt24 = new TGeoXtru(3);
    sTt24->SetName("ITS sup Long Air Tube Tray for RB24 Side Tt24");
    xp[0]  = sT24->GetX(0);
    yp[0]  = sT24->GetY(0);
    xp[1]  = sTl24->GetX(0);
    yp[1]  = sTl24->GetY(0);
    xp[2]  = -xp[1];
    yp[2]  =  yp[1];
    xp[3]  = -xp[0];
    yp[3]  =  yp[0];
    sTt24->DefinePolygon(4,xp,yp);
    sTt24->DefineSection(0,zA24[0]-kfrm24Width,0.0,0.0,1.0);
    sTt24->DefineSection(1,zA24[iRmin],0.0,0.0,1.0);
    sTt24->DefineSection(2,zA24[kfrm24NZsections]+kfrm24Width,0.0,
                         sA24[kfrm24NZsections]->GetRmax()-sA24[0]->GetRmin());
    // Inner opening for cooling (lower) {inside sTt24}
    sU24 = new TGeoXtru(3);
    sU24->SetName("ITS sup Cable Tray Cooling tube space RB24 Side U24");
    xp[0] = sTt24->GetX(0) + kct24AlThick;
    yp[0] = sTt24->GetY(0) + kct24AlThick;
    xp[1] = sTt24->GetX(1) + kct24AlThick;
    yp[1] = sTt24->GetY(1) - kct24AlThick;
    xp[2] = -xp[1];
    yp[2] =  yp[1];
    xp[3] = -xp[0];
    yp[3] =  yp[0];
    sU24->DefinePolygon(4,xp,yp);
    for(i=0;i<sTt24->GetNz();i++){
        sU24->DefineSection(i,sTt24->GetZ(i),sTt24->GetXOffset(i),
                            sTt24->GetYOffset(i),sTt24->GetScale(i));
    } // end for i
    // Inner opening for cables (upper) {inside sTl24}
    sVl24 = new TGeoXtru(3);
    sVl24->SetName("ITS sup Cable Tray Cable space RB24 Side Vl24");
    xp[0] = sTl24->GetX(0)+2.0*kct24AlThick;
    yp[0] = sTl24->GetY(0);
    yp[1] = yp[0] + kct24CoolCableDivEar;
    xp[1] = Xfrom2Points(sTl24->GetX(0),sTl24->GetY(0),
                         sTl24->GetX(1),sTl24->GetY(1),yp[1])+2.0*kct24AlThick;
    yp[2] = yp[1];
    xp[2] = xp[1] - kct24AlThick;
    yp[3] = sTl24->GetY(3) - kct24AlThick;
    xp[3] = Xfrom2Points(sTl24->GetX(0),sTl24->GetY(0),sTl24->GetX(1),
                         sTl24->GetY(1),yp[3]) + kct24AlThick;
    xp[4] = -xp[3];
    yp[4] =  yp[3];
    xp[5] = -xp[2];
    yp[5] =  yp[2];
    xp[6] = -xp[1];
    yp[6] =  yp[1];
    xp[7] = -xp[0];
    yp[7] =  yp[0];
    sVl24->DefinePolygon(8,xp,yp);
    for(i=0;i<sTl24->GetNz();i++){
        sVl24->DefineSection(i,sTl24->GetZ(i),sTl24->GetXOffset(i),
                            sTl24->GetYOffset(i),sTl24->GetScale(i));
    } // end for i
    // Inner opening for cables (upper) {inside sTs24}
    sVs24 = new TGeoXtru(3);
    sVs24->SetName("ITS sup Cable Tray Cable space RB24 Side Vs24");
    sVs24->DefinePolygon(8,xp,yp);
    for(i=0;i<8;i++){
    xp[i]  = sVl24->GetX(i);
    yp[i]  = sVl24->GetY(i);
    } // end for i
    for(i=0;i<sTl24->GetNz();i++){
        sVs24->DefineSection(i,sTs24->GetZ(i),sTs24->GetXOffset(i),
                            sTs24->GetYOffset(i),sTs24->GetScale(i));
    } // end for i
    //------------------------------------------------------------------
    // Patch Pannels on RB 24 Side
    rot  = new TGeoRotation("",0.0,0.0,-kft24Theta); // Gets Used later as well
    rot1 = new TGeoRotation("",0.0,0.0,kft24Theta);  // Gets Used later as well
    s3PP24 = new TGeoXtru(2);
    s3PP24->SetName("ITS sup 3 bay pach pannel RB24 side 3PP24");
    yp[5]  = sT24->GetY(7) + kct24CoolSectionH;
    xp[5]  = Xfrom2Points(sT24->GetX(7),sT24->GetY(7),sT24->GetX(6),
                          sT24->GetY(6),yp[6]);
    yp[6]  = sT24->GetY(0) + kct24CoolSectionH;
    xp[6]  =  Xfrom2Points(sT24->GetX(0),sT24->GetY(0),sT24->GetX(1),
                          sT24->GetY(1),yp[9]);
    local[0] = xp[6]; local[1] = yp[6]; local[2] = 0.0;
    rot1->LocalToMaster(local,master);
    xp[0]  = master[0];
    yp[0]  = master[1];
    local[0] = xp[6]; local[1] = yp[6] + kft24PPHightSDDSSD; local[2] = 0.0;
    rot1->LocalToMaster(local,master);
    xp[1]  = master[0];
    yp[1]  = master[1];
    xp[2]  = -xp[1];
    yp[2]  =  yp[1];
    xp[3]  = -xp[0];
    yp[3]  =  yp[0];
    local[0] = xp[6]; local[1] = yp[6]; local[2] = 0.0;
    rot1->MasterToLocal(local,master);
    xp[4]  = master[0];
    yp[4]  = master[1];
    local[0] = xp[5]; local[1] = yp[5]; local[2] = 0.0;
    rot1->LocalToMaster(local,master);
    xp[7]  = master[0];
    yp[7]  = master[1];
    s3PP24->DefinePolygon(8,xp,yp);
    s3PP24->DefineSection(0,0.0);
    s3PP24->DefineSection(1,kft24PPlength);
    //
    s2PP24 = new TGeoXtru(2);
    s2PP24->SetName("ITS sup 2 bay pach pannel RB24 side 2PP24");
    local[1] = sTl24->GetY(3); local[2] = 0.0;
    local[0] = Xfrom2Points(sTl24->GetX(0),sTl24->GetY(0),
                            sTl24->GetX(1),sTl24->GetY(1),local[1]);
    rot1->LocalToMaster(local,master);
    xp[0]  = master[0];
    yp[0]  = master[1];
    local[1] = sTl24->GetY(3) + kft24PPHightSPDFMD; local[2] = 0.0;
    local[0] = Xfrom2Points(sTl24->GetX(0),sTl24->GetY(0),
                            sTl24->GetX(1),sTl24->GetY(1),local[1]);
    rot1->LocalToMaster(local,master);
    xp[1]  = master[0];
    yp[1]  = master[1];
    yp[2]  = sTl24->GetY(4) + kft24PPHightSPDFMD;
    xp[2]  = Xfrom2Points(sTl24->GetX(6),sTl24->GetY(6),
                          sTl24->GetX(7),sTl24->GetY(7),yp[2]);
    yp[3]  = sTl24->GetY(7);
    xp[3]  = Xfrom2Points(sTl24->GetX(6),sTl24->GetY(6),
                          sTl24->GetX(7),sTl24->GetY(7),yp[3]);
    xp[4]  = sTl24->GetX(3);
    yp[4]  = sTl24->GetY(3);
    local[0] = sTl24->GetX(4);local[1] = sTl24->GetY(4); local[2] = 0.0;
    rot1->LocalToMaster(local,master);
    xp[5]  = master[0];
    yp[5]  = master[1];
    s2PP24->DefinePolygon(6,xp,yp);
    s2PP24->DefineSection(0,0.0);
    s2PP24->DefineSection(1,kft24PPlength);
    //
    sV3PP24 = new TGeoXtru(2);
    sV3PP24->SetName("ITS sup Patch Pannel 3 Bay inside Rb24 side V3PP24");
    xp[0] = s3PP24->GetX(0) + kct24AlThick;
    yp[0] = s3PP24->GetY(0) + kct24AlThick;
    local[1] = s3PP24->GetY(6) + kft24PPHightSDDSSD - kct24AlThick;local[2]=0.;
    local[0] = Xfrom2Points(sTl24->GetX(0),sTl24->GetY(0),
                           sTl24->GetX(1),sTl24->GetY(1),local[1]);
    rot1->LocalToMaster(local,master);
    xp[1] = master[0];
    yp[1] = master[1];
    xp[2] = -xp[1];
    yp[2] =  yp[1];
    xp[3] = -xp[0];
    yp[3] =  yp[0];
    xp[4] = s3PP24->GetX(4);
    yp[4] = s3PP24->GetY(4);
    xp[5] = s3PP24->GetX(5);
    yp[5] = s3PP24->GetY(5);
    xp[6] = s3PP24->GetX(6);
    yp[6] = s3PP24->GetY(6);
    xp[7] = s3PP24->GetX(7);
    yp[7] = s3PP24->GetY(7);
    sV3PP24->DefinePolygon(8,xp,yp);
    sV3PP24->DefineSection(0,s3PP24->GetZ(0),s3PP24->GetXOffset(0),
                           s3PP24->GetYOffset(0),s3PP24->GetScale(0));
    sV3PP24->DefineSection(1,s3PP24->GetZ(1),s3PP24->GetXOffset(1),
                           s3PP24->GetYOffset(1),s3PP24->GetScale(1));
    //
    sV2PP24 = new TGeoXtru(2);
    sV2PP24->SetName("ITS sup Patch Pannel 2 Bay inside Rb24 side V2PP24");
    xp[0] = s2PP24->GetX(0) + kct24AlThick;
    yp[0] = s2PP24->GetY(0) + kct24AlThick;
    local[1] = sTl24->GetY(3) + kft24PPHightSPDFMD - kct24AlThick;local[2]=0.;
    local[0] = Xfrom2Points(sTl24->GetX(0),sTl24->GetY(0),
                           sTl24->GetX(1),sTl24->GetY(1),local[1]);
    rot1->LocalToMaster(local,master);
    xp[1] = master[0];
    yp[1] = master[1];
    yp[2] = sTl24->GetY(4) + kft24PPHightSPDFMD - kct24AlThick;
    xp[2] = Xfrom2Points(sTl24->GetX(6),sTl24->GetY(6),
                           sTl24->GetX(7),sTl24->GetY(7),yp[2]);
    yp[3] = sTl24->GetY(4);
    xp[3] = Xfrom2Points(sTl24->GetX(6),sTl24->GetY(6),
                           sTl24->GetX(7),sTl24->GetY(7),yp[3]);;
    xp[4] = s2PP24->GetX(4);
    yp[4] = s2PP24->GetY(4);
    xp[5] = s2PP24->GetX(5);
    yp[5] = s2PP24->GetY(5);
    sV2PP24->DefinePolygon(6,xp,yp);
    sV2PP24->DefineSection(0,s2PP24->GetZ(0),s2PP24->GetXOffset(0),
                           s2PP24->GetYOffset(0),s2PP24->GetScale(0));
    sV2PP24->DefineSection(1,s2PP24->GetZ(1),s2PP24->GetXOffset(1),
                           s2PP24->GetYOffset(1),s2PP24->GetScale(1));
    // RB 24 Tray Mother Volume
    sMT24 = new TGeoPcon("ITS sup Cable Tray Mother Volume RB24 MT24",
                         0.0,360.0,5);
    sMT24->Z(0)    = 0.0;
    sMT24->Rmin(0) = sA24[0]->GetRmax();
    sMT24->Rmax(0) = TMath::Max(TMath::Hypot(s3PP24->GetX(1),s3PP24->GetY(1)),
                                TMath::Hypot(s2PP24->GetX(1),s2PP24->GetY(1)));

    sMT24->Z(1)    = sMT24->GetZ(0) + kft24PPlength;
    sMT24->Rmin(1) = sMT24->GetRmin(0);
    sMT24->Rmax(1) = sMT24->GetRmax(0);
    sMT24->Z(2)    = sMT24->GetZ(1);
    sMT24->Rmin(2) = sMT24->GetRmin(0);
    sMT24->Rmax(2) = sMT24->GetRmax(0) - kft24PPHightSPDFMD;

    sMT24->Z(3)    = sMT24->GetZ(0) + zA24[iRmin] - zA24[0] -kfrm24Width;
    sMT24->Rmin(3) = sA24[iRmin]->GetRmin();
    sMT24->Rmax(3) = TMath::Hypot(sT24->GetX(3),sT24->GetY(3));
    sMT24->Z(4)    = sMT24->GetZ(0) + zA24[kfrm24NZsections] + kfrm24Width  - 
        zA24[0] -kfrm24Width;
    sMT24->Rmin(4) = sA24[kfrm24NZsections]->GetRmax();
    sMT24->Rmax(4) = TMath::Hypot(sT24->GetX(3)+sT24->GetXOffset(2),
                                  sT24->GetY(3)+sT24->GetYOffset(2));
    //
    if(GetDebug(1)){
        sT24->InspectShape();
        sW24->InspectShape();
        sTl24->InspectShape();
        sTs24->InspectShape();
        sTt24->InspectShape();
        sU24->InspectShape();
        sVl24->InspectShape();
        sVs24->InspectShape();
        s3PP24->InspectShape();
        s2PP24->InspectShape();
        sV3PP24->InspectShape();
        sV2PP24->InspectShape();
        sMT24->InspectShape();
    } // end if GetDebug(1)
    //
    TGeoVolume *vC24[kct24Ntrays],*vT24[kct24Ntrays],*vPP24[kft24NPatchPannels];
    TGeoVolume *vWTV024,*vW24,*vU24,*vUFMD24,*vVl24,*vVlFMD24,*vVs24;
    TGeoVolume *vV3PP24,*vV2PP24,*vV2PPFMD24;
    TGeoVolumeAssembly *vMT24;
    vMT24 = new TGeoVolumeAssembly("ITSsupCableTrayMotherMT24");
    //vMT24->SetVisibility(kTRUE);
    //vMT24->SetLineColor(8); // white
    //vMT24->SetLineWidth(1);
    //vMT24->SetFillColor(vMT24->GetLineColor());
    //vMT24->SetFillStyle(4100); // 100% transparent
    //
    vU24 = new TGeoVolume("ITSsupCableTrayLowerU24",sU24,medSUPair);
    vU24->SetVisibility(kTRUE);
    vU24->SetLineColor(7); // light blue
    vU24->SetLineWidth(1);
    vU24->SetFillColor(vU24->GetLineColor());
    vU24->SetFillStyle(4090); // 90% transparent
    vUFMD24 = new TGeoVolume("FMDsupCableTrayLowerU24",sU24,medSUPair);
    vUFMD24->SetVisibility(kTRUE);
    vUFMD24->SetLineColor(7); // light blue
    vUFMD24->SetLineWidth(1);
    vUFMD24->SetFillColor(vUFMD24->GetLineColor());
    vUFMD24->SetFillStyle(4090); // 90% transparent
    vVl24 = new TGeoVolume("ITSsupCableTrayUpperV24",sVl24,medSUPair);
    vVl24->SetVisibility(kTRUE);
    vVl24->SetLineColor(7); // light blue
    vVl24->SetLineWidth(1);
    vVl24->SetFillColor(vVl24->GetLineColor());
    vVl24->SetFillStyle(4090); // 90% transparent
    vVlFMD24 = new TGeoVolume("FMDsupCableTrayUpperVl24",sVl24,medSUPair);
    vVlFMD24->SetVisibility(kTRUE);
    vVlFMD24->SetLineColor(7); // light blue
    vVlFMD24->SetLineWidth(1);
    vVlFMD24->SetFillColor(vVlFMD24->GetLineColor());
    vVlFMD24->SetFillStyle(4090); // 90% transparent
    vVs24 = new TGeoVolume("ITSsupCableTrayUpperVs24",sVs24,medSUPair);
    vVs24->SetVisibility(kTRUE);
    vVs24->SetLineColor(7); // light blue
    vVs24->SetLineWidth(1);
    vVs24->SetFillColor(vVs24->GetLineColor());
    vVs24->SetFillStyle(4090); // 90% transparent
    vW24 = new TGeoVolume("ITSsupCableTrayUpperW24",sW24,medSUPair);
    vW24->SetVisibility(kTRUE);
    vW24->SetLineColor(7); // light blue
    vW24->SetLineWidth(1);
    vW24->SetFillColor(vW24->GetLineColor());
    vW24->SetFillStyle(4090); // 90% transparent
    //
    vWTV024 = new TGeoVolume("V0supCableTrayUpperWTV024",sW24,medSUPair);
    vWTV024->SetVisibility(kTRUE);
    vWTV024->SetLineColor(7); // light blue
    vWTV024->SetLineWidth(1);
    vWTV024->SetFillColor(vWTV024->GetLineColor());
    vWTV024->SetFillStyle(4090); // 90% transparent
    //
    vV3PP24 = new TGeoVolume("ITSsup3BayPachPannelInsideV3PP24",sV3PP24,medSUPair);
    vV3PP24->SetVisibility(kTRUE);
    vV3PP24->SetLineColor(8); // white
    vV3PP24->SetLineWidth(1);
    vV3PP24->SetFillColor(vV3PP24->GetLineColor());
    vV3PP24->SetFillStyle(4100); // 100% transparent
    vV2PP24 = new TGeoVolume("ITSsup2BayPachPannelInsideV2PP24",sV2PP24,medSUPair);
    vV2PP24->SetVisibility(kTRUE);
    vV2PP24->SetLineColor(8); // white
    vV2PP24->SetLineWidth(1);
    vV2PP24->SetFillColor(vV2PP24->GetLineColor());
    vV2PP24->SetFillStyle(4100); // 100% transparent
    vV2PPFMD24 = new TGeoVolume("FMDsup2BayPachPannelInsideV2PP24",sV2PP24,medSUPair);
    vV2PPFMD24->SetVisibility(kTRUE);
    vV2PPFMD24->SetLineColor(8); // white
    vV2PPFMD24->SetLineWidth(1);
    vV2PPFMD24->SetFillColor(vV2PPFMD24->GetLineColor());
    vV2PPFMD24->SetFillStyle(4100); // 100% transparent
    //
    //delete rot;
    //delete rot1;
    //
    Double_t tha[kct24Ntrays],thb[kft24NPatchPannels];
    for(i=0;i<kct24Ntrays/4;i++) {
        if(i==0) tha[0] = 17.0+0.5*kft24Theta;
        else tha[i] = tha[i-1] + kft24Theta;
        tha[i+  kct24Ntrays/4] =  90.0 + tha[i];
        tha[i+  kct24Ntrays/2] = 180.0 + tha[i];
        tha[i+3*kct24Ntrays/4] = 270.0 + tha[i];
    } // end for i
    if(GetDebug(1)) for(i=0;i<kct24Ntrays;i++) Info("ServicesCableSupport",
                                                  "tha[%d]=%f",i,tha[i]);
    Char_t *airName[kct24Ntrays]={"FMD0","SDD0","SSD0","SSD1","SPD0","SPD1",
                                  "TV00","SDD1","SDD2","SPD2","SPD3","ALG0",
                                  "SPD4","SPD5","SSD2","SSD3","SPD6","SPD7",
                                  "TV01","SDD3","SDD4","SPD8","SPD9","ALG1",
                                  "FMD1","SDD5","SSD4","SSD5","SPDA","SPDB",
                                  "TV02","SDD6","SDD7","SPDC","SPDD","ALG2",
                                  "SPDE","SPDF","SSD6","SSD7","SPDG","SPDH",
                                  "TV03","SDD8","SDD9","SPDI","SPDJ","ALG3"};
    Char_t *trayName[kct24Ntrays]={"FMD0","SSD0","SSD1","SSD2","SSD3","SPD0",
                                   "TV00","SDD0","SDD1","SDD2","SPD1","ALG0",
                                   "SPD2","SSD4","SSD5","SSD6","SSD7","SPD3",
                                   "TV01","SDD3","SDD4","SDD5","SPD4","ALG1",
                                   "FMD1","SSD8","SSD9","SSDA","SSDB","SPD5",
                                   "TV02","SDD6","SDD7","SDD8","SPD6","ALG2",
                                   "SPD7","SSDC","SSDD","SSDE","SSDF","SPD8",
                                   "TV03","SDD9","SDDA","SDDB","SPD9","ALG3"};
    //
    //Int_t ncopyW24=1,ncopyU24=1,ncopyV24=1;
    j = 0;
    for(i=0;i<kct24Ntrays;i++){
        if(strncmp(trayName[i],"FMD",3)==0){
            sprintf(name,"FMDsupCableTrayT24[%s]",trayName[i]);
            vT24[i] = new TGeoVolume(name,sTl24,medSUPal);
            vT24[i]->AddNode(vVlFMD24,1,0);
        }else if(strncmp(trayName[i],"TV0",3)==0){
            sprintf(name,"V0supCableTrayT24[%s]",trayName[i]);
            vT24[i] = new TGeoVolume(name,sT24,medSUPal);
            vT24[i]->AddNode(vWTV024,1,0);
        }else if(strncmp(trayName[i],"ALG",3)==0){ // ITS Alignment Channel
            sprintf(name,"ITSsupCableTrayT24[%s]",trayName[i]);
            vT24[i] = new TGeoVolume(name,sT24,medSUPal);
            vT24[i]->AddNode(vW24,1,0);
        }else  if(strncmp(trayName[i],"SPD",3)==0){ /*ITS SPD*/
            sprintf(name,"ITSsupCableTrayT24[%s]",trayName[i]);
            vT24[i] = new TGeoVolume(name,sTl24,medSUPal);
            vT24[i]->AddNode(vVl24,1,0);
        }else { /*ITS*/
            sprintf(name,"ITSsupCableTrayT24[%s]",trayName[i]);
            vT24[i] = new TGeoVolume(name,sTs24,medSUPal); /// replace solid
            vT24[i]->AddNode(vVs24,1,0);
        } // end if
        vT24[i]->SetVisibility(kTRUE);
        vT24[i]->SetLineColor(6); // purple
        vT24[i]->SetLineWidth(1);
        vT24[i]->SetFillColor(vT24[i]->GetLineColor());
        vT24[i]->SetFillStyle(4000); // 0% transparent
        rot = new TGeoRotation("",0.0,0.0,tha[i]-90.0);
        if(GetDebug(1)) rot->Print();
        vMT24->AddNode(vT24[i],1,rot);
        //
        if(strncmp(trayName[i],"FMD",3)==0){
            sprintf(name,"FMDsupAirTubeTrayT24[%s]",airName[i]);
            vC24[j] = new TGeoVolume(name,sTt24,medSUPair);
            vC24[j]->AddNode(vUFMD24,1,0);
        }else if(strncmp(trayName[i],"TV0",3)==0){
            continue;
        }else if(strncmp(trayName[i],"ALG",3)==0){
            continue;
        }else{ /*ITS*/
            sprintf(name,"ITSsupAirTubTrayT24[%s]",airName[i]);
            vC24[j] = new TGeoVolume(name,sTt24,medSUPair);
            vC24[j]->AddNode(vU24,1,0);
        } // end if
        vC24[j]->SetVisibility(kTRUE);
        vC24[j]->SetLineColor(6); // purple
        vC24[j]->SetLineWidth(1);
        vC24[j]->SetFillColor(vC24[j]->GetLineColor());
        vC24[j]->SetFillStyle(4000); // 0% transparent
        vMT24->AddNode(vC24[j++],1,rot);
    } // end for i
    for(i=0;i<kft24NPatchPannels/4;i++) {
        if(i==0) thb[0] = 17.0+0.5*kft24Theta;
        else{
            if(i%2) thb[i] = thb[i-1] + 3.0*kft24Theta;
            else thb[i] = thb[i-1] + 2.0*kft24Theta;
        } // end if-else
        thb[i+  kft24NPatchPannels/4] =  90.0 + thb[i];
        thb[i+  kft24NPatchPannels/2] = 180.0 + thb[i];
        thb[i+3*kft24NPatchPannels/4] = 270.0 + thb[i];
    } // end for i
    Char_t *pachName[kft24NPatchPannels]={"FMD0","SSD0","SPD0","SDD0","SPD1",
                                          "SPD2","SSD1","SPD3","SDD1","SPD4",
                                          "FMD1","SSD2","SPD5","SDD2","SPD6",
                                          "SPD7","SSD3","SPD8","SDD3","SPD9"};
    for(i=0;i<kft24NPatchPannels;i++){
        if(strncmp(pachName[i],"FMD",3)==0){
            sprintf(name,"FMDsupPatchPannelPP24[%s]",pachName[i]);
            vPP24[i] = new TGeoVolume(name,s2PP24,medSUPal);
            vPP24[i]->AddNode(vV2PPFMD24,1,0);
        }else if(strncmp(pachName[i],"SPD",3)==0){ /*ITS SPD*/
            sprintf(name,"ITSsupPathcPannelPP24[%s]",pachName[i]);
            vPP24[i] = new TGeoVolume(name,s2PP24,medSUPal);
            vPP24[i]->AddNode(vV2PP24,1,0);
        }else { /*ITS*/
            sprintf(name,"ITSsupPathcPannelPP24[%s]",pachName[i]);
            vPP24[i] = new TGeoVolume(name,s3PP24,medSUPal); /// replace solid
            vPP24[i]->AddNode(vV3PP24,1,0);
        } // end if
        vPP24[i]->SetVisibility(kTRUE);
        vPP24[i]->SetLineColor(6); // purple
        vPP24[i]->SetLineWidth(1);
        vPP24[i]->SetFillColor(vPP24[i]->GetLineColor());
        vPP24[i]->SetFillStyle(4000); // 0% transparent
        rot = new TGeoRotation("",0.0,0.0,thb[i]-90.0);
        if(GetDebug(1)) rot->Print();
        vMT24->AddNode(vPP24[i],1,rot);
    } // end for i
    tran = new TGeoTranslation("",0.0,0.0,kfrm24Z0);
    moth->AddNode(vMT24,1,tran);
    if(GetDebug(1)){
        for(i=0;i<kct24Ntrays;i++) vT24[i]->PrintNodes();
        for(i=0;i<kct24Ntrays-8;i++) vC24[i]->PrintNodes();
        vU24->PrintNodes();
        vUFMD24->PrintNodes();
        vVl24->PrintNodes();
        vVlFMD24->PrintNodes();
        vVs24->PrintNodes();
        vW24->PrintNodes();
        vWTV024->PrintNodes();
        vMT24->PrintNodes();
    } // end if
    //==================================================================
    //
    // RB 26, Muon Absober side
    const Double_t kfrm26Z0           = -900*fgkmm;//SSup_203A.jpg
    const Double_t kfrm26Thss         = 5.0*fgkmm;
    const Double_t kfrm26R0ss         = 444.5*fgkmm-kfrm26Thss; //SSup_204A.jpg
    const Double_t kfrm26R1ss         = 601.6*fgkmm-kfrm26Thss; //SSup_208A.jpg
    const Double_t kfrm26Width        = 10.0*fgkmm;
    //const Double_t kfrm26Hight       = 10.0*fgkmm;
    const Double_t kfrm26Phi0         = 15.2*fgkDegree; // SSup_602A.jpg
    const Double_t kfrm26Phi1         = (90.0-7.6)*fgkDegree; // SSup_802A.jpg
    const Double_t kfrm26ZssSection   = (415.0-10.0)*fgkmm;
    const Int_t    kfrm26NZsections   = 4;
    const Int_t    kfrm26NPhiSections = 4;
    const Int_t    kfrm26NPhi         = 4;
    TGeoConeSeg *sA26[kfrm26NZsections+1];//,*sM26;//Cylinderial support structure
    TGeoArb8     *sB26; // Cylinderial support structure
    /*
    sM26 = new TGeoConeSeg("ITS sup Cable tray support frame mother volume "
                          "M26",0.5*(4.*kfrm26ZssSection+5*kfrm26Width),
                          kfrm26R1ss,kfrm26R1ss+kfrm26Thss,
                          kfrm26R0ss,kfrm26R0ss+kfrm26Thss,
                          kfrm26Phi0,kfrm26Phi1);
    */
    m = -((kfrm26R1ss-kfrm26R0ss)/
         (((Double_t)kfrm26NZsections)*(kfrm26ZssSection+kfrm26Width)));
    for(i=0;i<kfrm26NZsections+1;i++){
        di = ((Double_t) i)*(kfrm26ZssSection+kfrm26Width);
        sprintf(name,
                "ITS sup Cable tray support frame radial section A26[%d]",i);
        r1 = kfrm26R1ss+m*di;
        r2 = kfrm26R1ss+m*(di+kfrm26Width);
        sA26[i] = new TGeoConeSeg(name,0.5*kfrm26Width,r2,r2+kfrm26Thss,
                                 r1,r1+kfrm26Thss,kfrm26Phi0,kfrm26Phi1);
    } // end for i
    sB26 = new TGeoArb8("ITS sup Cable tray support frame Z section B26",
                       0.5*kfrm26ZssSection);
    r = 0.25*(sA26[0]->GetRmax1()+sA26[0]->GetRmin1()+
              sA26[1]->GetRmax2()+sA26[1]->GetRmin2());
    sB26->SetVertex(0,sA26[0]->GetRmax2()-r,+0.5*kfrm26Width);
    sB26->SetVertex(1,sA26[0]->GetRmax2()-r,-0.5*kfrm26Width);
    sB26->SetVertex(2,sA26[0]->GetRmin2()-r,-0.5*kfrm26Width);
    sB26->SetVertex(3,sA26[0]->GetRmin2()-r,+0.5*kfrm26Width);
    sB26->SetVertex(4,sA26[1]->GetRmax1()-r,+0.5*kfrm26Width);
    sB26->SetVertex(5,sA26[1]->GetRmax1()-r,-0.5*kfrm26Width);
    sB26->SetVertex(6,sA26[1]->GetRmin1()-r,-0.5*kfrm26Width);
    sB26->SetVertex(7,sA26[1]->GetRmin1()-r,+0.5*kfrm26Width);
    if(GetDebug(1)){
        for(i=0;i<kfrm26NZsections+1;i++) sA26[i]->InspectShape();
        //sM26->InspectShape();
        sB26->InspectShape();
    } // end if GetDebug(1)
    //
    TGeoVolume *vA26[kfrm26NZsections+1],*vB26;
    TGeoVolumeAssembly *vM26;
    //
    for(i=0;i<kfrm26NZsections+1;i++){
        sprintf(name,"ITSsupFrameA26[%d]",i);
        vA26[i] = new TGeoVolume(name,sA26[i],medSUPss);
        vA26[i]->SetVisibility(kTRUE);
        vA26[i]->SetLineColor(1); // black
        vA26[i]->SetLineWidth(1);
        vA26[i]->SetFillColor(vA26[i]->GetLineColor());
        vA26[i]->SetFillStyle(4000); // 0% transparent
    } // end for i
    vB26 = new TGeoVolume("ITSsupFrameB26",sB26,medSUPss);
    vB26->SetVisibility(kTRUE);
    vB26->SetLineColor(1); // black
    vB26->SetLineWidth(1);
    vB26->SetFillColor(vB26->GetLineColor());
    vB26->SetFillStyle(4000); // 0% transparent
    vM26 = new TGeoVolumeAssembly("ITSsupFrameM26");
    //vM26 = new TGeoVolume("ITSsupFrameM26",sM26,medSUPair);
    //vM26->SetVisibility(kTRUE);
    //vM26->SetLineColor(7); // light blue
    //vM26->SetLineWidth(1);
    //vM26->SetFillColor(vM26->GetLineColor());
    //vM26->SetFillStyle(4090); // 90% transparent
    //
    Int_t ncopyB26=1;
    t0 = kfrm26Phi0;
    dt = (kfrm26Phi1-kfrm26Phi0)/((Double_t)kfrm26NPhiSections);
    for(i=0;i<=kfrm26NZsections;i++){
        di = ((Double_t) i)*(kfrm26ZssSection+kfrm26Width);
        z = 0.5*(4.*kfrm26ZssSection+5*kfrm26Width);
        z = -z+sA26[i]->GetDz() + di;
        tran = new TGeoTranslation("",0.0,0.0,z);
        vM26->AddNode(vA26[i],1,tran);
        z = z+sB26->GetDz();
        if(i<kfrm26NZsections)for(j=0;j<=kfrm26NPhiSections;j++){
            r = 0.25*(sA26[i]->GetRmax1()+sA26[i]->GetRmin1()+
                      sA26[i+1]->GetRmax2()+sA26[i+1]->GetRmin2());
            t = t0 + ((Double_t)j)*dt;
            rot = new TGeoRotation("",0.0,0.0,t);
            y = r*SinD(t);
            x = r*CosD(t);
            tranrot = new TGeoCombiTrans("",x,y,z,rot);
            //delete rot; // rot not explicity used in AddNode functions.
            vM26->AddNode(vB26,ncopyB26++,tranrot);
        } // end for j
    } // end for i
    tran = new TGeoTranslation("",0.0,0.0,kfrm26Z0-0.5*(4.*kfrm26ZssSection+5*kfrm26Width));
    moth->AddNode(vM26,1,tran);
    for(i=1;i<kfrm26NPhi;i++){
        rot = new TGeoRotation("",0.0,0.0,90.0*((Double_t)i));
        tranrot = new TGeoCombiTrans(*tran,*rot);
        //delete rot; // rot not explicity used in AddNode functions.
        moth->AddNode(vM26,i+1,tranrot);
    } // end for i
    if(GetDebug(1)){
        for(i=0;i<kfrm26NZsections+1;i++) vA26[i]->PrintNodes();
        vB26->PrintNodes();
        vM26->PrintNodes();
    } // end if
}
