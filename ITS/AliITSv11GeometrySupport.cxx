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
void AliITSv11GeometrySupport::SPDCone(TGeoVolume *moth,const TGeoManager *mgr)
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
  const Double_t kWideWing      = 6.0*fgkcm;
  const Double_t kThetaWing     = 45.0;
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

  CreateSPDOmegaShape(xair,yair,kThicknessOmega,xomega,yomega);

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

  CreateSPDOmegaShape(xair,yair,kThicknessOmega,xomega,yomega);

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
    Double_t thetaW = kThetaWing*(2*i+1) - angleWideWing/2.;
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
     Double_t   t, Double_t *x , Double_t *y ) const
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
                   const Double_t *xin, const Double_t *yin, Double_t  d,
		   Double_t   *x, Double_t *y)
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
// Updated:      20 Feb 2009  Mario Sitta       New algorithm (the old one
//                                              gives erroneous vertexes)
//

  // This vector contains the index of those points which coincide
  // with the corresponding points in the air shape
  Int_t indexAir2Omega[12] = {1, 2, 5, 6, 9, 10, 11, 15, 16, 19, 20, 23};

  // First fill those vertexes corresponding to
  // the edges aligned to the air shape edges
  for (Int_t j=0; j<12; j++) {
    x[*(indexAir2Omega+j)] = xin[j];
    y[*(indexAir2Omega+j)] = yin[j];
  }

  // Now get the coordinates of the first inner point
  PointFromParallelLines(x[23],y[23],x[1],y[1],d,x[0],y[0]);

  // Knowing this, the second internal point can be determined
  InsidePoint(x[0],y[0],x[1],y[1],x[2],y[2],d,x[22],y[22]);

  // The third point is now computable
  ReflectPoint(x[1],y[1],x[2],y[2],x[22],y[22],x[21],y[21]);

  // Repeat this logic
  InsidePoint(x[21],y[21],x[20],y[20],x[19],y[19],-d,x[3],y[3]);

  ReflectPoint(x[20],y[20],x[19],y[19],x[3],y[3],x[4],y[4]);

  InsidePoint(x[4],y[4],x[5],y[5],x[6],y[6],d,x[18],y[18]);

  ReflectPoint(x[5],y[5],x[6],y[6],x[18],y[18],x[17],y[17]);

  InsidePoint(x[17],y[17],x[16],y[16],x[15],y[15],-d,x[7],y[7]);

  ReflectPoint(x[16],y[16],x[15],y[15],x[7],y[7],x[8],y[8]);

  InsidePoint(x[8],y[8],x[9],y[9],x[10],y[10],d,x[14],y[14]);

  // These need to be fixed explicitly
  x[12] = x[11];
  y[12] = y[11] + d;
  x[13] = x[10] + d;
  y[13] = y[12];

  // Finally reflect on the negative side
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
						Double_t *x, Double_t *y) const
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
void AliITSv11GeometrySupport::PointFromParallelLines(Double_t x1, Double_t y1,
			      Double_t x2, Double_t y2, Double_t d,
			      Double_t &x, Double_t &y) const
{
//
// Determines the X and Y of the first internal point of the Omega shape
// (i.e. the coordinates of a point given two parallel lines passing by
// two points and placed at a known distance)
//
// Input:
//        x1, y1 : first point
//        x2, y2 : second point
//        d      : distance between the two lines
//
// Output:
//        x, y   : coordinate of the point
//
// Created:      22 Feb 2009  Mario Sitta
//
//Begin_Html
/*
<img src="ITS/doc/PointFromParallelLines.gif">
*/
//End_Html

  // The slope of the paralles lines at a distance d
  Double_t m; 

  // The parameters of the solving equation
  // a x^2 - 2 b x + c = 0
  Double_t a = (x1 - x2)*(x1 - x2) - d*d;
  Double_t b = (x1 - x2)*(y1 - y2);
  Double_t c = (y1 - y2)*(y1 - y2) - d*d;

  // (delta4 is Delta/4 because we use the reduced formula)
  Double_t delta4 = b*b - a*c;

  // Compute the slope of the two parallel lines
  // (one of the two possible slopes, the one with the smaller
  // absolute value is needed)
  if (delta4 < 0) { // Should never happen with our data, but just to be sure
    x = -1;         // x is expected positive, so this flags an error
    return;
  } else
    m = (b + TMath::Sqrt(delta4))/a;  // b is negative with our data

  // Finally compute the coordinates of the point
  x = x2 + (y1 - y2 - d)/m;
  y = y1 - d;

  // Done
  return;
}

//______________________________________________________________________
void AliITSv11GeometrySupport::ReflectPoint(Double_t x1, Double_t y1,
					    Double_t x2, Double_t y2,
					    Double_t x3, Double_t y3,
					    Double_t &x, Double_t &y) const
{
//
// Given two points (x1,y1) and (x2,y2), determines the point (x,y)
// lying on the line parallel to the line passing by these points,
// at a distance d and passing by the point (x3,y3), which is symmetric to
// the third point with respect to the axis of the segment delimited by
// the two first points.
//
// Input:
//        x1, y1 : first point
//        x2, y2 : second point
//        x3, y3 : third point
//        d      : distance between the two lines
//
// Output:
//        x, y   : coordinate of the reflected point
//
// Created:      22 Feb 2009  Mario Sitta
//
//Begin_Html
/*
<img src="ITS/doc/ReflectPoint.gif">
*/
//End_Html

  // The slope of the line passing by the first two points
  Double_t k = (y2 - y1)/(x2 - x1);

  // The middle point of the segment 1-2
  Double_t xK = (x1 + x2)/2.;
  Double_t yK = (y1 + y2)/2.;

  // The intercept between the axis of the segment 1-2 and the line
  // passing by 3 and parallel to the line passing by 1-2
  Double_t xH = (k*k*x3 + k*(yK - y3) + xK)/(k*k + 1);
  Double_t yH = k*(xH - x3) + y3;

  // The point symmetric to 3 with respect to H
  x = 2*xH - x3;
  y = 2*yH - y3;

  // Done
  return;
}

//______________________________________________________________________
void AliITSv11GeometrySupport::SDDCone(TGeoVolume *moth,const TGeoManager *mgr)
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
// Updated:      25 Jul 2008  Mario Sitta   SDDCarbonFiberCone simpler
// Updated:      10 Jun 2010  Mario Sitta   Cables across cone holes added
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
  const Double_t kConeRCurv          =      10.0*fgkmm; // Radius of curvature
  const Double_t kConeRinMin         = (210.0/2)*fgkmm;
//  const Double_t kConeRinMax         = (216.0/2)*fgkmm;
  const Double_t kConeRinCylinder    = (231.0/2)*fgkmm;
  const Double_t kConeZCylinder      =     192.0*fgkmm;
  const Double_t kConeZOuterMilled   =      23.0*fgkmm;
  const Double_t kConeDZin           =      15.0*fgkmm; // ???
  const Double_t kConeThickness      =      10.0*fgkmm; // Rohacell + Carb.Fib.
  const Double_t kConeTheta          =      45.0*fgkDegree; // SDD cone angle
  const Double_t kSinConeTheta       =
                                     TMath::Sin(kConeTheta*TMath::DegToRad());
  const Double_t kCosConeTheta       =
                                     TMath::Cos(kConeTheta*TMath::DegToRad());
  const Double_t kTanConeTheta       =
                                     TMath::Tan(kConeTheta*TMath::DegToRad());
  // Dimensions of the Cone Inserts
  const Double_t kConeCFThickness    =       1.5*fgkmm;//Carbon fiber thickness
  // Dimensions of the Cone Holes
  const Double_t kHole1RMin          = (450.0/2)*fgkmm;
  const Double_t kHole1RMax          = (530.0/2)*fgkmm;
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
  const Double_t kHole4Width         =        30*fgkmm;
  //  const Int_t    kNHole4             =         3      ;
  // Fraction of materials in holes
  const Double_t kHolePlasticFrac    =       0.55846;
  const Double_t kHoleCuFrac         =       0.06319;
  const Double_t kHoleGlassFrac      =       0.02652;

  // Local variables
  Double_t x, y, z, t, dza, rmin, rmax;


  // Recover the needed materials
  TGeoMedium *medSDDcf    = mgr->GetMedium("ITS_SDD C (M55J)$");
  TGeoMedium *medSDDair   = mgr->GetMedium("ITS_SDD AIR$");
  TGeoMedium *medSDDste   = mgr->GetMedium("ITS_G10FR4$"); // stesalite
  TGeoMedium *medSDDroh   = mgr->GetMedium("ITS_ROHACELL$");
  TGeoMedium *medSDDss    = mgr->GetMedium("ITS_INOX$");
  TGeoMedium *medSDDplast = mgr->GetMedium("ITS_SDDKAPTON (POLYCH2)$");
  TGeoMedium *medSDDCu    = mgr->GetMedium("ITS_COPPER$");
  TGeoMedium *medSDDglass = mgr->GetMedium("ITS_SDD OPTICFIB$");

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
    x = kBoltRadius*CosD(t);
    y = kBoltRadius*SinD(t);
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

  TGeoPcon *coneshape = new TGeoPcon(0.0, 360.0, 10);

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

  coneshape->Z(6)     = kConeZCylinder - kConeDZin;

  RadiusOfCurvature(kConeRCurv,90.0,coneshape->GetZ(6),0.0,
		    90.0-kConeTheta,z,rmin);
  coneshape->Z(5)     = z;
  coneshape->Rmin(5)  = RminFromZpCone(coneshape,3,kConeTheta,z);
  coneshape->Rmax(5)  = RmaxFromZpCone(coneshape,4,kConeTheta,z);

  RadiusOfCurvature(kConeRCurv,90.-kConeTheta,
		    0.0,coneshape->Rmin(5),90.0,z,rmin);
  coneshape->Rmin(6)  = rmin;
  coneshape->Rmax(6)  = RmaxFromZpCone(coneshape,4,kConeTheta,
				       coneshape->GetZ(6));

  coneshape->Z(7)     = coneshape->GetZ(6);
  coneshape->Rmin(7)  = kConeRinMin;
  coneshape->Rmax(7)  = coneshape->GetRmax(6);

  coneshape->Rmin(8)  = kConeRinMin;

  RadiusOfCurvature(kConeRCurv,90.0,kConeZCylinder,kConeRinCylinder,
		    90.0-kConeTheta,z,rmax);
  coneshape->Z(8)     = z;
  coneshape->Rmax(8)  = rmax;

  coneshape->Z(9)     = kConeZCylinder;
  coneshape->Rmin(9)  = kConeRinMin;
  coneshape->Rmax(9)  = kConeRinCylinder;


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
  // A single hole volume gives an overlap with coneinsert, so
  // three contiguous volumes are created: one to be put in the cone foam
  // and two in the cone carbon fiber envelope
  TGeoPcon *hole1shape = new TGeoPcon(-kHole1Phi/2., kHole1Phi, 4);

  hole1shape->Rmin(0) = kHole1RMax;
  hole1shape->Rmax(0) = hole1shape->GetRmin(0);
  hole1shape->Z(0)    = ZFromRminpCone(conefoamshape,0,kConeTheta,
				       hole1shape->GetRmin(0));

  hole1shape->Rmax(1) = hole1shape->GetRmax(0);
  hole1shape->Z(1)    = ZFromRmaxpCone(conefoamshape,3,kConeTheta,
				       hole1shape->GetRmax(1));
  hole1shape->Rmin(1) = RminFromZpCone(conefoamshape,1,kConeTheta,
				       hole1shape->GetZ(1));

  hole1shape->Rmin(2) = kHole1RMin;
  hole1shape->Z(2)    = ZFromRminpCone(conefoamshape,1,kConeTheta,
				       hole1shape->GetRmin(2));
  hole1shape->Rmax(2) = RmaxFromZpCone(conefoamshape,3,kConeTheta,
				       hole1shape->GetZ(2));

  hole1shape->Rmin(3) = hole1shape->GetRmin(2);
  hole1shape->Rmax(3) = hole1shape->GetRmin(3);
  hole1shape->Z(3)    = ZFromRmaxpCone(conefoamshape,3,kConeTheta,
				       hole1shape->GetRmax(3));

  TGeoPcon *hole11shape = new TGeoPcon(-kHole1Phi/2., kHole1Phi, 4);

  hole11shape->Rmin(0) = kHole1RMax;
  hole11shape->Rmax(0) = hole11shape->GetRmin(0);
  hole11shape->Z(0)    = ZFromRminpCone(coneshape,3,kConeTheta,
					hole11shape->GetRmin(0));

  hole11shape->Rmax(1) = hole11shape->GetRmax(0);
  hole11shape->Z(1)    = ZFromRminpCone(coneinsertshape,3,kConeTheta,
					hole11shape->GetRmax(1));
  hole11shape->Rmin(1) = RminFromZpCone(coneshape,3,kConeTheta,
					hole11shape->GetZ(1));

  hole11shape->Rmin(2) = kHole1RMin;
  hole11shape->Z(2)    = ZFromRminpCone(coneshape,3,kConeTheta,
					hole11shape->GetRmin(2));
  hole11shape->Rmax(2) = RminFromZpCone(coneinsertshape,3,kConeTheta,
					hole11shape->GetZ(2));

  hole11shape->Rmin(3) = hole11shape->GetRmin(2);
  hole11shape->Rmax(3) = hole11shape->GetRmin(3);
  hole11shape->Z(3)    = ZFromRminpCone(coneinsertshape,3,kConeTheta,
					hole11shape->GetRmax(3));

  TGeoPcon *hole12shape = new TGeoPcon(-kHole1Phi/2., kHole1Phi, 4);

  hole12shape->Rmin(0) = kHole1RMax;
  hole12shape->Rmax(0) = hole12shape->GetRmin(0);
  hole12shape->Z(0)    = ZFromRmaxpCone(coneinsertshape,4,kConeTheta,
					hole12shape->GetRmin(0));

  hole12shape->Rmax(1) = hole12shape->GetRmax(0);
  hole12shape->Z(1)    = ZFromRmaxpCone(coneshape,4,kConeTheta,
					hole12shape->GetRmax(1));
  hole12shape->Rmin(1) = RmaxFromZpCone(coneinsertshape,4,kConeTheta,
					hole12shape->GetZ(1));

  hole12shape->Rmin(2) = kHole1RMin;
  hole12shape->Z(2)    = ZFromRmaxpCone(coneinsertshape,4,kConeTheta,
					hole12shape->GetRmin(2));
  hole12shape->Rmax(2) = RmaxFromZpCone(coneshape,4,kConeTheta,
					hole12shape->GetZ(2));

  hole12shape->Rmin(3) = hole12shape->GetRmin(2);
  hole12shape->Rmax(3) = hole12shape->GetRmin(3);
  hole12shape->Z(3)    = ZFromRmaxpCone(coneshape,4,kConeTheta,
					hole12shape->GetRmax(3));

  //
  TGeoPcon *hole2shape = new TGeoPcon(-kHole2Phi/2., kHole2Phi, 4);

  hole2shape->Rmin(0) = kHole2RMax;
  hole2shape->Rmax(0) = hole2shape->GetRmin(0);
  hole2shape->Z(0)    = ZFromRminpCone(conefoamshape,0,kConeTheta,
				       hole2shape->GetRmin(0));

  hole2shape->Rmax(1) = hole2shape->GetRmax(0);
  hole2shape->Z(1)    = ZFromRmaxpCone(conefoamshape,3,kConeTheta,
				       hole2shape->GetRmax(1));
  hole2shape->Rmin(1) = RminFromZpCone(conefoamshape,1,kConeTheta,
				       hole2shape->GetZ(1));

  hole2shape->Rmin(2) = kHole2RMin;
  hole2shape->Z(2)    = ZFromRminpCone(conefoamshape,1,kConeTheta,
				       hole2shape->GetRmin(2));
  hole2shape->Rmax(2) = RmaxFromZpCone(conefoamshape,3,kConeTheta,
				       hole2shape->GetZ(2));

  hole2shape->Rmin(3) = hole2shape->GetRmin(2);
  hole2shape->Rmax(3) = hole2shape->GetRmin(3);
  hole2shape->Z(3)    = ZFromRmaxpCone(conefoamshape,3,kConeTheta,
				       hole2shape->GetRmax(3));

  TGeoPcon *hole21shape = new TGeoPcon(-kHole2Phi/2., kHole2Phi, 4);

  hole21shape->Rmin(0) = kHole2RMax;
  hole21shape->Rmax(0) = hole21shape->GetRmin(0);
  hole21shape->Z(0)    = ZFromRminpCone(coneshape,3,kConeTheta,
					hole21shape->GetRmin(0));

  hole21shape->Rmax(1) = hole21shape->GetRmax(0);
  hole21shape->Z(1)    = ZFromRminpCone(coneinsertshape,3,kConeTheta,
					hole21shape->GetRmax(1));
  hole21shape->Rmin(1) = RminFromZpCone(coneshape,3,kConeTheta,
					hole21shape->GetZ(1));

  hole21shape->Rmin(2) = kHole2RMin;
  hole21shape->Z(2)    = ZFromRminpCone(coneshape,3,kConeTheta,
					hole21shape->GetRmin(2));
  hole21shape->Rmax(2) = RminFromZpCone(coneinsertshape,3,kConeTheta,
					hole21shape->GetZ(2));

  hole21shape->Rmin(3) = hole21shape->GetRmin(2);
  hole21shape->Rmax(3) = hole21shape->GetRmin(3);
  hole21shape->Z(3)    = ZFromRminpCone(coneinsertshape,3,kConeTheta,
					hole21shape->GetRmax(3));

  TGeoPcon *hole22shape = new TGeoPcon(-kHole2Phi/2., kHole2Phi, 4);

  hole22shape->Rmin(0) = kHole2RMax;
  hole22shape->Rmax(0) = hole22shape->GetRmin(0);
  hole22shape->Z(0)    = ZFromRmaxpCone(coneinsertshape,4,kConeTheta,
					hole22shape->GetRmin(0));

  hole22shape->Rmax(1) = hole22shape->GetRmax(0);
  hole22shape->Z(1)    = ZFromRmaxpCone(coneshape,4,kConeTheta,
					hole22shape->GetRmax(1));
  hole22shape->Rmin(1) = RmaxFromZpCone(coneinsertshape,4,kConeTheta,
					hole22shape->GetZ(1));

  hole22shape->Rmin(2) = kHole2RMin;
  hole22shape->Z(2)    = ZFromRmaxpCone(coneinsertshape,4,kConeTheta,
					hole22shape->GetRmin(2));
  hole22shape->Rmax(2) = RmaxFromZpCone(coneshape,4,kConeTheta,
					hole22shape->GetZ(2));

  hole22shape->Rmin(3) = hole22shape->GetRmin(2);
  hole22shape->Rmax(3) = hole22shape->GetRmin(3);
  hole22shape->Z(3)    = ZFromRmaxpCone(coneshape,4,kConeTheta,
					hole22shape->GetRmax(3));

  //
  Double_t holePhi;
  holePhi = (kHole3Width/kHole3RMin)*TMath::RadToDeg();

  TGeoPcon *hole3shape = new TGeoPcon(-holePhi/2., holePhi, 4);

  hole3shape->Rmin(0) = kHole3RMin + kHole3DeltaR;
  hole3shape->Rmax(0) = hole3shape->GetRmin(0);
  hole3shape->Z(0)    = ZFromRminpCone(conefoamshape,0,kConeTheta,
				       hole3shape->GetRmin(0));

  hole3shape->Rmax(1) = hole3shape->GetRmax(0);
  hole3shape->Z(1)    = ZFromRmaxpCone(conefoamshape,3,kConeTheta,
				       hole3shape->GetRmax(1));
  hole3shape->Rmin(1) = RminFromZpCone(conefoamshape,1,kConeTheta,
				       hole3shape->GetZ(1));

  hole3shape->Rmin(2) = kHole3RMin;
  hole3shape->Z(2)    = ZFromRminpCone(conefoamshape,1,kConeTheta,
				       hole3shape->GetRmin(2));
  hole3shape->Rmax(2) = RmaxFromZpCone(conefoamshape,3,kConeTheta,
				       hole3shape->GetZ(2));

  hole3shape->Rmin(3) = hole3shape->GetRmin(2);
  hole3shape->Rmax(3) = hole3shape->GetRmin(3);
  hole3shape->Z(3)    = ZFromRmaxpCone(conefoamshape,3,kConeTheta,
				       hole3shape->GetRmax(3));

  TGeoPcon *hole31shape = new TGeoPcon(-holePhi/2., holePhi, 4);

  hole31shape->Rmin(0) = kHole3RMin + kHole3DeltaR;
  hole31shape->Rmax(0) = hole31shape->GetRmin(0);
  hole31shape->Z(0)    = ZFromRminpCone(coneshape,3,kConeTheta,
					hole31shape->GetRmin(0));

  hole31shape->Rmax(1) = hole31shape->GetRmax(0);
  hole31shape->Z(1)    = ZFromRminpCone(coneinsertshape,3,kConeTheta,
					hole31shape->GetRmax(1));
  hole31shape->Rmin(1) = RminFromZpCone(coneshape,3,kConeTheta,
					hole31shape->GetZ(1));

  hole31shape->Rmin(2) = kHole3RMin;
  hole31shape->Z(2)    = ZFromRminpCone(coneshape,3,kConeTheta,
					hole31shape->GetRmin(2));
  hole31shape->Rmax(2) = RminFromZpCone(coneinsertshape,3,kConeTheta,
					hole31shape->GetZ(2));

  hole31shape->Rmin(3) = hole31shape->GetRmin(2);
  hole31shape->Rmax(3) = hole31shape->GetRmin(3);
  hole31shape->Z(3)    = ZFromRminpCone(coneinsertshape,3,kConeTheta,
					hole31shape->GetRmax(3));

  TGeoPcon *hole32shape = new TGeoPcon(-holePhi/2., holePhi, 4);

  hole32shape->Rmin(0) = kHole3RMin + kHole3DeltaR;
  hole32shape->Rmax(0) = hole32shape->GetRmin(0);
  hole32shape->Z(0)    = ZFromRmaxpCone(coneinsertshape,4,kConeTheta,
					hole32shape->GetRmin(0));

  hole32shape->Rmax(1) = hole32shape->GetRmax(0);
  hole32shape->Z(1)    = ZFromRmaxpCone(coneshape,4,kConeTheta,
					hole32shape->GetRmax(1));
  hole32shape->Rmin(1) = RmaxFromZpCone(coneinsertshape,4,kConeTheta,
					hole32shape->GetZ(1));

  hole32shape->Rmin(2) = kHole3RMin;
  hole32shape->Z(2)    = ZFromRmaxpCone(coneinsertshape,4,kConeTheta,
					hole32shape->GetRmin(2));
  hole32shape->Rmax(2) = RmaxFromZpCone(coneshape,4,kConeTheta,
					hole32shape->GetZ(2));

  hole32shape->Rmin(3) = hole32shape->GetRmin(2);
  hole32shape->Rmax(3) = hole32shape->GetRmin(3);
  hole32shape->Z(3)    = ZFromRmaxpCone(coneshape,4,kConeTheta,
					hole32shape->GetRmax(3));

  //
  holePhi = (kHole4Width/kHole4RMin)*TMath::RadToDeg();

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

  // Cables to be put inside the holes: Pcon's
  // (fractions are manually computed from AliITSv11GeometrySDD::SDDCables
  TGeoPcon *hole1plastshape = new TGeoPcon(-kHole1Phi/2., kHole1Phi, 4);

  hole1plastshape->Rmin(0) = hole1shape->GetRmin(0);
  hole1plastshape->Rmax(0) = hole1shape->GetRmax(0);
  hole1plastshape->Z(0)    = hole1shape->GetZ(0);

  hole1plastshape->Rmin(1) = hole1shape->GetRmin(1);
  hole1plastshape->Rmax(1) = hole1shape->GetRmax(1);
  hole1plastshape->Z(1)    = hole1shape->GetZ(1);

  dza = hole1plastshape->GetRmax(0) - (kHole1RMax-kHole1RMin)*kHolePlasticFrac;

  hole1plastshape->Rmin(2) = dza;
  hole1plastshape->Z(2)    = ZFromRminpCone(conefoamshape,1,kConeTheta,
					    hole1plastshape->GetRmin(2));
  hole1plastshape->Rmax(2) = RmaxFromZpCone(conefoamshape,3,kConeTheta,
					    hole1plastshape->GetZ(2));

  hole1plastshape->Rmin(3) = hole1plastshape->GetRmin(2);
  hole1plastshape->Rmax(3) = hole1plastshape->GetRmin(3);
  hole1plastshape->Z(3)    = ZFromRmaxpCone(conefoamshape,3,kConeTheta,
					    hole1plastshape->GetRmax(3));

  TGeoPcon *hole1Cushape = new TGeoPcon(-kHole1Phi/2., kHole1Phi, 4);

  hole1Cushape->Rmin(0) = hole1plastshape->GetRmin(2);
  hole1Cushape->Rmax(0) = hole1Cushape->GetRmin(0);
  hole1Cushape->Z(0)    = hole1plastshape->GetZ(2);

  dza = hole1Cushape->GetRmax(0) - (kHole1RMax-kHole1RMin)*kHoleCuFrac;

  hole1Cushape->Rmin(1) = dza;
  hole1Cushape->Rmax(1) = hole1Cushape->GetRmax(0);
  hole1Cushape->Z(1)    = ZFromRminpCone(conefoamshape,1,kConeTheta,
					 hole1Cushape->GetRmin(1));

  hole1Cushape->Rmax(2) = hole1Cushape->GetRmax(0);
  hole1Cushape->Rmin(2) = hole1Cushape->GetRmin(1);
  hole1Cushape->Z(2)    = hole1plastshape->GetZ(3);

  hole1Cushape->Rmin(3) = hole1Cushape->GetRmin(1);
  hole1Cushape->Rmax(3) = hole1Cushape->GetRmin(3);
  hole1Cushape->Z(3)    = ZFromRmaxpCone(conefoamshape,3,kConeTheta,
					 hole1Cushape->GetRmax(3));

  TGeoPcon *hole1glassshape = new TGeoPcon(-kHole1Phi/2., kHole1Phi, 4);

  hole1glassshape->Rmin(0) = hole1Cushape->GetRmin(1);
  hole1glassshape->Rmax(0) = hole1glassshape->GetRmin(0);
  hole1glassshape->Z(0)    = hole1Cushape->GetZ(1);

  dza = hole1glassshape->GetRmax(0) - (kHole1RMax-kHole1RMin)*kHoleGlassFrac;

  hole1glassshape->Rmin(1) = dza;
  hole1glassshape->Rmax(1) = hole1glassshape->GetRmax(0);
  hole1glassshape->Z(1)    = ZFromRminpCone(conefoamshape,1,kConeTheta,
					    hole1glassshape->GetRmin(1));

  hole1glassshape->Rmax(2) = hole1glassshape->GetRmax(0);
  hole1glassshape->Rmin(2) = hole1glassshape->GetRmin(1);
  hole1glassshape->Z(2)    = hole1Cushape->GetZ(3);

  hole1glassshape->Rmin(3) = hole1glassshape->GetRmin(1);
  hole1glassshape->Rmax(3) = hole1glassshape->GetRmin(3);
  hole1glassshape->Z(3)    = ZFromRmaxpCone(conefoamshape,3,kConeTheta,
					    hole1glassshape->GetRmax(3));
  //
  TGeoPcon *hole2plastshape = new TGeoPcon(-kHole2Phi/2., kHole2Phi, 4);

  hole2plastshape->Rmin(0) = hole2shape->GetRmin(0);
  hole2plastshape->Rmax(0) = hole2shape->GetRmax(0);
  hole2plastshape->Z(0)    = hole2shape->GetZ(0);

  hole2plastshape->Rmin(1) = hole2shape->GetRmin(1);
  hole2plastshape->Rmax(1) = hole2shape->GetRmax(1);
  hole2plastshape->Z(1)    = hole2shape->GetZ(1);

  dza = hole2plastshape->GetRmax(0) - (kHole2RMax-kHole2RMin)*kHolePlasticFrac;

  hole2plastshape->Rmin(2) = dza;
  hole2plastshape->Z(2)    = ZFromRminpCone(conefoamshape,1,kConeTheta,
					    hole2plastshape->GetRmin(2));
  hole2plastshape->Rmax(2) = RmaxFromZpCone(conefoamshape,3,kConeTheta,
					    hole2plastshape->GetZ(2));

  hole2plastshape->Rmin(3) = hole2plastshape->GetRmin(2);
  hole2plastshape->Rmax(3) = hole2plastshape->GetRmin(3);
  hole2plastshape->Z(3)    = ZFromRmaxpCone(conefoamshape,3,kConeTheta,
					    hole2plastshape->GetRmax(3));

  TGeoPcon *hole2Cushape = new TGeoPcon(-kHole2Phi/2., kHole2Phi, 4);

  hole2Cushape->Rmin(0) = hole2plastshape->GetRmin(2);
  hole2Cushape->Rmax(0) = hole2Cushape->GetRmin(0);
  hole2Cushape->Z(0)    = hole2plastshape->GetZ(2);

  dza = hole2Cushape->GetRmax(0) - (kHole2RMax-kHole2RMin)*kHoleCuFrac;

  hole2Cushape->Rmin(1) = dza;
  hole2Cushape->Rmax(1) = hole2Cushape->GetRmax(0);
  hole2Cushape->Z(1)    = ZFromRminpCone(conefoamshape,1,kConeTheta,
					 hole2Cushape->GetRmin(1));

  hole2Cushape->Rmax(2) = hole2Cushape->GetRmax(0);
  hole2Cushape->Rmin(2) = hole2Cushape->GetRmin(1);
  hole2Cushape->Z(2)    = hole2plastshape->GetZ(3);

  hole2Cushape->Rmin(3) = hole2Cushape->GetRmin(1);
  hole2Cushape->Rmax(3) = hole2Cushape->GetRmin(3);
  hole2Cushape->Z(3)    = ZFromRmaxpCone(conefoamshape,3,kConeTheta,
					 hole2Cushape->GetRmax(3));

  TGeoPcon *hole2glassshape = new TGeoPcon(-kHole2Phi/2., kHole2Phi, 4);

  hole2glassshape->Rmin(0) = hole2Cushape->GetRmin(1);
  hole2glassshape->Rmax(0) = hole2glassshape->GetRmin(0);
  hole2glassshape->Z(0)    = hole2Cushape->GetZ(1);

  dza = hole2glassshape->GetRmax(0) - (kHole2RMax-kHole2RMin)*kHoleGlassFrac;

  hole2glassshape->Rmin(1) = dza;
  hole2glassshape->Rmax(1) = hole2glassshape->GetRmax(0);
  hole2glassshape->Z(1)    = ZFromRminpCone(conefoamshape,1,kConeTheta,
					    hole2glassshape->GetRmin(1));

  hole2glassshape->Rmax(2) = hole2glassshape->GetRmax(0);
  hole2glassshape->Rmin(2) = hole2glassshape->GetRmin(1);
  hole2glassshape->Z(2)    = hole2Cushape->GetZ(3);

  hole2glassshape->Rmin(3) = hole2glassshape->GetRmin(1);
  hole2glassshape->Rmax(3) = hole2glassshape->GetRmin(3);
  hole2glassshape->Z(3)    = ZFromRmaxpCone(conefoamshape,3,kConeTheta,
					    hole2glassshape->GetRmax(3));


  // Debug if requested
  if (GetDebug(1)) {
    coneshape->InspectShape();
    coneinsertshape->InspectShape();
    conefoamshape->InspectShape();
    hole1shape->InspectShape();
    hole2shape->InspectShape();
    hole3shape->InspectShape();
    hole4shape->InspectShape();
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

  TGeoVolume *hole11 = new TGeoVolume("SDDCableHole11",
				      hole11shape,medSDDair);
  hole11->SetVisibility(kTRUE);
  hole11->SetLineColor(5); // Yellow
  hole11->SetLineWidth(1);
  hole11->SetFillColor(hole11->GetLineColor());
  hole11->SetFillStyle(4090); // 90% transparent

  TGeoVolume *hole12 = new TGeoVolume("SDDCableHole12",
				      hole12shape,medSDDair);
  hole12->SetVisibility(kTRUE);
  hole12->SetLineColor(5); // Yellow
  hole12->SetLineWidth(1);
  hole12->SetFillColor(hole12->GetLineColor());
  hole12->SetFillStyle(4090); // 90% transparent

  TGeoVolume *hole1plast = new TGeoVolume("SDDCableHole1Plast",
					  hole1plastshape,medSDDplast);
  hole1plast->SetVisibility(kTRUE);
  hole1plast->SetLineColor(kBlue);
  hole1plast->SetLineWidth(1);
  hole1plast->SetFillColor(hole1plast->GetLineColor());
  hole1plast->SetFillStyle(4090); // 90% transparent

  TGeoVolume *hole1Cu = new TGeoVolume("SDDCableHole1Cu",
				       hole1Cushape,medSDDCu);
  hole1Cu->SetVisibility(kTRUE);
  hole1Cu->SetLineColor(kRed);
  hole1Cu->SetLineWidth(1);
  hole1Cu->SetFillColor(hole1Cu->GetLineColor());
  hole1Cu->SetFillStyle(4090); // 90% transparent

  TGeoVolume *hole1glass = new TGeoVolume("SDDCableHole1glass",
					  hole1glassshape,medSDDglass);
  hole1glass->SetVisibility(kTRUE);
  hole1glass->SetLineColor(kGreen);
  hole1glass->SetLineWidth(1);
  hole1glass->SetFillColor(hole1glass->GetLineColor());
  hole1glass->SetFillStyle(4090); // 90% transparent

  TGeoVolume *hole2 = new TGeoVolume("SDDCableHole2",
				     hole2shape,medSDDair);
  hole2->SetVisibility(kTRUE);
  hole2->SetLineColor(5); // Yellow
  hole2->SetLineWidth(1);
  hole2->SetFillColor(hole2->GetLineColor());
  hole2->SetFillStyle(4090); // 90% transparent

  TGeoVolume *hole21 = new TGeoVolume("SDDCableHole21",
				      hole21shape,medSDDair);
  hole21->SetVisibility(kTRUE);
  hole21->SetLineColor(5); // Yellow
  hole21->SetLineWidth(1);
  hole21->SetFillColor(hole21->GetLineColor());
  hole21->SetFillStyle(4090); // 90% transparent

  TGeoVolume *hole22 = new TGeoVolume("SDDCableHole22",
				      hole22shape,medSDDair);
  hole22->SetVisibility(kTRUE);
  hole22->SetLineColor(5); // Yellow
  hole22->SetLineWidth(1);
  hole22->SetFillColor(hole22->GetLineColor());
  hole22->SetFillStyle(4090); // 90% transparent

  TGeoVolume *hole2plast = new TGeoVolume("SDDCableHole2Plast",
					  hole2plastshape,medSDDplast);
  hole2plast->SetVisibility(kTRUE);
  hole2plast->SetLineColor(kBlue);
  hole2plast->SetLineWidth(1);
  hole2plast->SetFillColor(hole2plast->GetLineColor());
  hole2plast->SetFillStyle(4090); // 90% transparent

  TGeoVolume *hole2Cu = new TGeoVolume("SDDCableHole2Cu",
				       hole2Cushape,medSDDCu);
  hole2Cu->SetVisibility(kTRUE);
  hole2Cu->SetLineColor(kRed);
  hole2Cu->SetLineWidth(1);
  hole2Cu->SetFillColor(hole2Cu->GetLineColor());
  hole2Cu->SetFillStyle(4090); // 90% transparent

  TGeoVolume *hole2glass = new TGeoVolume("SDDCableHole2glass",
					  hole2glassshape,medSDDglass);
  hole2glass->SetVisibility(kTRUE);
  hole2glass->SetLineColor(kGreen);
  hole2glass->SetLineWidth(1);
  hole2glass->SetFillColor(hole2glass->GetLineColor());
  hole2glass->SetFillStyle(4090); // 90% transparent

  TGeoVolume *hole3 = new TGeoVolume("SDDCableHole3",
				     hole3shape,medSDDair);
  hole3->SetVisibility(kTRUE);
  hole3->SetLineColor(5); // Yellow
  hole3->SetLineWidth(1);
  hole3->SetFillColor(hole3->GetLineColor());
  hole3->SetFillStyle(4090); // 90% transparent

  TGeoVolume *hole31 = new TGeoVolume("SDDCableHole31",
				      hole31shape,medSDDair);
  hole31->SetVisibility(kTRUE);
  hole31->SetLineColor(5); // Yellow
  hole31->SetLineWidth(1);
  hole31->SetFillColor(hole31->GetLineColor());
  hole31->SetFillStyle(4090); // 90% transparent

  TGeoVolume *hole32 = new TGeoVolume("SDDCableHole32",
				      hole32shape,medSDDair);
  hole32->SetVisibility(kTRUE);
  hole32->SetLineColor(5); // Yellow
  hole32->SetLineWidth(1);
  hole32->SetFillColor(hole32->GetLineColor());
  hole32->SetFillStyle(4090); // 90% transparent

  TGeoVolume *hole4 = new TGeoVolume("SDDCableHole4",
				     hole4shape,medSDDair);
  hole4->SetVisibility(kTRUE);
  hole4->SetLineColor(5); // Yellow
  hole4->SetLineWidth(1);
  hole4->SetFillColor(hole4->GetLineColor());
  hole4->SetFillStyle(4090); // 90% transparent

  // Mount up a cone
  cfconeinsert->AddNode(cfconefoam,1,0);

  hole1->AddNode(hole1plast, 1, 0);
  hole1->AddNode(hole1Cu, 1, 0);
  hole1->AddNode(hole1glass, 1, 0);

  hole2->AddNode(hole2plast, 1, 0);
  hole2->AddNode(hole2Cu, 1, 0);
  hole2->AddNode(hole2glass, 1, 0);

  for (Int_t i=0; i<12; i++) {
    Double_t phiH = i*30.0;
    cfconefoam->AddNode(hole1 , i+1, new TGeoRotation("", 0, 0, phiH));
        cfcone->AddNode(hole11, i+1, new TGeoRotation("", 0, 0, phiH));
        cfcone->AddNode(hole12, i+1, new TGeoRotation("", 0, 0, phiH));
  }

  for (Int_t i=0; i<6; i++) {
    Double_t phiH = i*60.0;
    cfconefoam->AddNode(hole2 , i+1, new TGeoRotation("", 0, 0, phiH));
        cfcone->AddNode(hole21, i+1, new TGeoRotation("", 0, 0, phiH));
        cfcone->AddNode(hole22, i+1, new TGeoRotation("", 0, 0, phiH));
  }

  for (Int_t i=0; i<kNHole3; i++) {
    Double_t phiH0 = 360./(Double_t)kNHole3;
    Double_t phiH  = i*phiH0 + 0.5*phiH0;
    cfconefoam->AddNode(hole3 , i+1, new TGeoRotation("", phiH, 0, 0));
        cfcone->AddNode(hole31, i+1, new TGeoRotation("", phiH, 0, 0));
        cfcone->AddNode(hole32, i+1, new TGeoRotation("", phiH, 0, 0));
  }

  cfcone->AddNode(cfconeinsert,1,0);

/*
  for (Int_t i=0; i<kNHole4; i++) {
    Double_t phiH0 = 360./(Double_t)kNHole4;
    Double_t phiH  = i*phiH0 + 0.25*phiH0;
    cfcone->AddNode(hole4, i+1, new TGeoRotation("", phiH, 0, 0));
  }
*/
  // Finally put everything in the mother volume
  moth->AddNode(cfcylinder,1,0);

  z = coneshape->Z(9);
  moth->AddNode(cfcone,1,new TGeoTranslation(0, 0, -z - kCylinderHalfLength));
  moth->AddNode(cfcone,2,new TGeoCombiTrans (0, 0,  z + kCylinderHalfLength,
			 new TGeoRotation("", 0, 180, 0)                   ));


  return;
}

//______________________________________________________________________
void AliITSv11GeometrySupport::SSDCone(TGeoVolume *moth,const TGeoManager *mgr)
{
//
// Creates the SSD support cone and cylinder geometry. as a
// volume assembly and adds it to the mother volume
// (part of this code is taken or anyway inspired to SSDCone method
// of AliITSv11GeometrySupport.cxx,v 1.9 2007/06/06)
//
// Input:
//         moth : the TGeoVolume owing the volume structure
//         mgr  : the GeoManager (default gGeoManager)
// Output:
//
// Created:         ???       Bjorn S. Nilsen
// Updated:      08 Mar 2008  Mario Sitta
//
// Technical data are taken from:  "ITS Supporto Generale" (technical
// drawings ALR3-0743/1, ALR3-0743/1A and ALR3-0743/1B), "Supporto Generale
// Settore SSD" (technical drawings ALR3-0743/2A and ALR3-0743/2E), private
// communication with B. Giraudo
//
// Updated:      11 Apr 2008  Mario Sitta
// Measures from drawings give overlaps with SPD thermal shield wings,
// so the terminal part of the SSD cone was reduced
//
// Updated:      30 Mar 2010  Mario Sitta
// Following M. van Leeuwen's suggestion on material budget, the thickness
// of the carbon fiber cylinder was increased from 0.6 to 0.625mm

  // Dimensions of the Central cylinder and flanges
  const Double_t kCylinderHalfLength   = (1144.0/2) *fgkmm;
  const Double_t kCylinderOuterRadius  = ( 595.0/2) *fgkmm;
  const Double_t kCylinderThickness    =       0.625*fgkmm;
  const Double_t kFoamHalfLength       = (1020.0/2) *fgkmm;
  const Double_t kFoamThickness        =        5.0 *fgkmm;
  const Double_t kFlangeHalfLength     =
                                      (kCylinderHalfLength-kFoamHalfLength)/2.;
  const Double_t kFlangeInnerRadius    = ( 563.0/2) *fgkmm;
  // Dimensions of the Cone
  const Double_t kConeROuterMin        = ( 957.0/2) *fgkmm;
  const Double_t kConeROuterMax        = ( 997.0/2) *fgkmm;
  const Double_t kConeRInnerMin        = ( 564.0/2) *fgkmm;
  const Double_t kConeRCurv1           =       10.0 *fgkmm;
  const Double_t kConeRCurv2           =       25.0 *fgkmm;
  const Double_t kConeCent1RCurv2      = ( 578.0/2) *fgkmm;
  const Double_t kConeCent2RCurv2      = ( 592.0/2) *fgkmm;
//  const Double_t kConeZOuterRing       =       47.0 *fgkmm;
//  const Double_t kConeZOuterRingInside =       30.25*fgkmm;
//  const Double_t kConeZInnerRing       =      161.5 *fgkmm;
//  const Double_t kConeZLength          =      176.5 *fgkmm;
  const Double_t kConeZOuterRing       =       38.5 *fgkmm;
  const Double_t kConeZOuterRingInside =       22.2 *fgkmm;
  const Double_t kConeZInnerRing       =      153.0 *fgkmm;
  const Double_t kConeZLength          =      168.0 *fgkmm;
  const Double_t kConeZPosition        = kConeZLength + kCylinderHalfLength;
  const Double_t kConeThickness        =       13.0 *fgkmm; // Cone thickness
  const Double_t kConeTheta            =       39.1 *fgkDegree; // Cone angle
  const Double_t kSinConeTheta         =
                                      TMath::Sin(kConeTheta*TMath::DegToRad());
  const Double_t kCosConeTheta         =
                                      TMath::Cos(kConeTheta*TMath::DegToRad());
  // Dimensions of the Foam cores
  const Double_t kConeFoam1Length      =      112.3 *fgkmm;
  const Double_t kConeFoam2Length      =       58.4 *fgkmm;
  // Dimensions of the Cone Holes
  const Double_t kCoolingHoleWidth     =       40.0 *fgkmm;
  const Double_t kCoolingHoleHight     =       30.0 *fgkmm;
  const Double_t kCoolingHoleRmin      =      350.0 *fgkmm;
  const Double_t kCoolingHolePhi       =       45.0 *fgkDegree;
  const Double_t kMountingHoleWidth    =       20.0 *fgkmm;
  const Double_t kMountingHoleHight    =       20.0 *fgkmm;
  const Double_t kMountingHoleRmin     =      317.5 *fgkmm;
  const Double_t kMountingHolePhi      =       60.0 *fgkDegree;
  const Double_t kCableHoleRin         = ( 800.0/2) *fgkmm;
  const Double_t kCableHoleRout        = ( 920.0/2) *fgkmm;
  const Double_t kCableHoleWidth       =      200.0 *fgkmm;
//  const Double_t kCableHoleAngle       =       42.0 *fgkDegree;
  // Dimensions of the Cone Wings
  const Double_t kWingRmax             =      527.5 *fgkmm;
  const Double_t kWingWidth            =       70.0 *fgkmm;
  const Double_t kWingHalfThick        = (  10.0/2) *fgkmm;
  const Double_t kThetaWing            =       45.0 *fgkDegree;
  // Dimensions of the SSD-SDD Mounting Brackets
  const Double_t kBracketRmin          = ( 541.0/2) *fgkmm;// See SDD ROutMin
  const Double_t kBracketRmax          = ( 585.0/2) *fgkmm;
  const Double_t kBracketHalfLength    = (   4.0/2) *fgkmm;
  const Double_t kBracketPhi           = (70.*fgkmm/kBracketRmax)*fgkRadian;
  // Common data
  const Double_t kCFThickness          =        0.75*fgkmm; //Carb. fib. thick.


  // Local variables
  Double_t rmin1, rmin2, rmax, z;

  //
  //Begin_Html
  /*
    <img src="picts/ITS/file_name.gif">
    <P>
    <FONT FACE'"TIMES">
    ITS SSD central support and thermal shield cylinder.
    </FONT>
    </P>
  */
  //End_Html
  //

  // Central cylinder with its internal foam and the lateral flanges:
  // a carbon fiber Pcon which contains a rohacell Tube and two
  // stesalite Cone's
  TGeoPcon *externalcylshape = new TGeoPcon(0,360,4);

  rmax  = kCylinderOuterRadius;
  rmin1 = kFlangeInnerRadius - kCylinderThickness;
  rmin2 = rmax - 2*kCylinderThickness - kFoamThickness;
  externalcylshape->DefineSection(0,-kCylinderHalfLength,rmin1,rmax);
  externalcylshape->DefineSection(1,-kFoamHalfLength    ,rmin2,rmax);
  externalcylshape->DefineSection(2, kFoamHalfLength    ,rmin2,rmax);
  externalcylshape->DefineSection(3, kCylinderHalfLength,rmin1,rmax);

  rmax  = kCylinderOuterRadius - kCylinderThickness;
  rmin1 = rmax - kFoamThickness;
  TGeoTube *foamshape = new TGeoTube(rmin1,rmax,kFoamHalfLength);

  rmax  = kCylinderOuterRadius - kCylinderThickness;
  rmin1 = rmax - kFoamThickness;
  rmin2 = kFlangeInnerRadius;
  TGeoCone *flangeshape = new TGeoCone(kFlangeHalfLength,
				       rmin1,rmax,rmin2,rmax);


  // We have the shapes: now create the real volumes

  TGeoMedium *medSSDcf  = mgr->GetMedium("ITS_SSD C (M55J)$");
  TGeoMedium *medSSDair = mgr->GetMedium("ITS_SSD AIR$");
  TGeoMedium *medSSDste = mgr->GetMedium("ITS_G10FR4$"); // stesalite
  TGeoMedium *medSSDroh = mgr->GetMedium("ITS_ROHACELL$");
  TGeoMedium *medSSDal  = mgr->GetMedium("ITS_ALUMINUM$");

  TGeoVolume *cfcylinder = new TGeoVolume("SSDexternalcylinder",
					   externalcylshape,medSSDcf);
  cfcylinder->SetVisibility(kTRUE);
  cfcylinder->SetLineColor(4); // blue
  cfcylinder->SetLineWidth(1);
  cfcylinder->SetFillColor(cfcylinder->GetLineColor());
  cfcylinder->SetFillStyle(4000); // 0% transparent

  TGeoVolume *foamcylinder = new TGeoVolume("SSDfoamcylinder",
					    foamshape,medSSDroh);
  foamcylinder->SetVisibility(kTRUE);
  foamcylinder->SetLineColor(3); // green
  foamcylinder->SetLineWidth(1);
  foamcylinder->SetFillColor(foamcylinder->GetLineColor());
  foamcylinder->SetFillStyle(4050); // 50% transparent

  TGeoVolume *flangecylinder = new TGeoVolume("SSDflangecylinder",
					      flangeshape,medSSDste);
  flangecylinder->SetVisibility(kTRUE);
  flangecylinder->SetLineColor(2); // red
  flangecylinder->SetLineWidth(1);
  flangecylinder->SetFillColor(flangecylinder->GetLineColor());
  flangecylinder->SetFillStyle(4050); // 50% transparent

  // Mount up the cylinder
  cfcylinder->AddNode(foamcylinder,1,0);
  cfcylinder->AddNode(flangecylinder,1,
	      new TGeoTranslation(0, 0, kFoamHalfLength+kFlangeHalfLength));
  cfcylinder->AddNode(flangecylinder,2,new TGeoCombiTrans(
              0, 0, -kFoamHalfLength-kFlangeHalfLength,
	      new TGeoRotation("",0,180,0)     ) );


  // The whole Cone as an assembly
  TGeoVolumeAssembly *vC = new TGeoVolumeAssembly("ITSssdCone");


  // SSD Support Cone with its internal inserts: a carbon fiber Pcon
  // with holes which contains a stesalite Pcon which on turn contains a
  // rohacell Pcon
  TGeoPcon *coneshape = new TGeoPcon(0.0, 360.0, 12);

  coneshape->Z(0)     = 0.0;
  coneshape->Rmin(0)  = kConeROuterMin;
  coneshape->Rmax(0)  = kConeROuterMax;

  coneshape->Z(1)     = kConeZOuterRingInside - kConeRCurv1;
  coneshape->Rmin(1)  = coneshape->GetRmin(0);
  coneshape->Rmax(1)  = coneshape->GetRmax(0);

  coneshape->Z(2)     = kConeZOuterRingInside;
  coneshape->Rmin(2)  = coneshape->GetRmin(1) - kConeRCurv1;
  coneshape->Rmax(2)  = coneshape->GetRmax(0);

  coneshape->Z(3)     = coneshape->GetZ(2);
  coneshape->Rmax(3)  = coneshape->GetRmax(0);

  coneshape->Z(4)     = kConeZOuterRing - kConeRCurv1;
  coneshape->Rmax(4)  = coneshape->GetRmax(0);

  coneshape->Z(5)     = kConeZOuterRing;
  coneshape->Rmax(5)  = coneshape->GetRmax(4) - kConeRCurv1;

  coneshape->Z(6)     = coneshape->GetZ(5);

  RadiusOfCurvature(kConeRCurv2,90.0,kConeZInnerRing,kConeCent1RCurv2,
		    90.0-kConeTheta,z,rmin1);
  coneshape->Z(7)     = z;
  coneshape->Rmin(7)  = rmin1;

  coneshape->Rmin(3)  = RminFromZpCone(coneshape,7,90.-kConeTheta,
				       coneshape->GetZ(3));

  coneshape->Rmin(4)  = RminFrom2Points(coneshape,3,7,coneshape->GetZ(4));

  coneshape->Rmin(5)  = RminFrom2Points(coneshape,3,7,coneshape->GetZ(5));

  coneshape->Rmin(6) = coneshape->GetRmin(5);

  coneshape->Z(8)     = kConeZInnerRing;
  coneshape->Rmin(8)  = kConeCent1RCurv2;

  coneshape->Z(9)     = coneshape->GetZ(8);
  coneshape->Rmin(9)  = kConeRInnerMin;

  RadiusOfCurvature(kConeRCurv2,90.0,kConeZLength,kConeCent2RCurv2,
		    90.0-kConeTheta,z,rmax);

  coneshape->Z(10)    = z;
  coneshape->Rmin(10) = coneshape->GetRmin(9);
  coneshape->Rmax(10) = rmax;

  coneshape->Rmax(6)  = RmaxFromZpCone(coneshape,10,90.-kConeTheta,
				       coneshape->GetZ(6));

  coneshape->Rmax(7)  = RmaxFrom2Points(coneshape,6,10,coneshape->GetZ(7));

  coneshape->Rmax(8)  = RmaxFrom2Points(coneshape,6,10,coneshape->GetZ(8));

  coneshape->Rmax(9)  = coneshape->GetRmax(8);

  coneshape->Z(11)    = kConeZLength;
  coneshape->Rmin(11) = coneshape->GetRmin(10);
  coneshape->Rmax(11) = kConeCent2RCurv2;

  // SSD Cone Insert: another Pcon
  Double_t x0, y0, x1, y1, x2, y2;
  TGeoPcon *coneinsertshape = new TGeoPcon(0.0,360.0,12);

  coneinsertshape->Z(0)     = coneshape->GetZ(0) + kCFThickness;
  coneinsertshape->Rmin(0)  = coneshape->GetRmin(0) + kCFThickness;
  coneinsertshape->Rmax(0)  = coneshape->GetRmax(0) - kCFThickness;

  x0 = coneshape->GetZ(0); y0 = coneshape->GetRmin(0);
  x1 = coneshape->GetZ(1); y1 = coneshape->GetRmin(1);
  x2 = coneshape->GetZ(2); y2 = coneshape->GetRmin(2);
  InsidePoint(x0, y0, x1, y1, x2, y2,  kCFThickness, z, rmin1);
  coneinsertshape->Z(1)     = z;
  coneinsertshape->Rmin(1)  = rmin1;
  coneinsertshape->Rmax(1)  = coneinsertshape->GetRmax(0);

  x0 = coneshape->GetZ(1); y0 = coneshape->GetRmin(1);
  x1 = coneshape->GetZ(2); y1 = coneshape->GetRmin(2);
  x2 = coneshape->GetZ(3); y2 = coneshape->GetRmin(3);
  InsidePoint(x0, y0, x1, y1, x2, y2,  kCFThickness, z, rmin1);
  coneinsertshape->Z(2)     = z;
  coneinsertshape->Rmin(2)  = rmin1;
  coneinsertshape->Rmax(2)  = coneinsertshape->GetRmax(1);

  x0 = coneshape->GetZ(2); y0 = coneshape->GetRmin(2);
  x1 = coneshape->GetZ(3); y1 = coneshape->GetRmin(3);
  x2 = coneshape->GetZ(4); y2 = coneshape->GetRmin(4);
  InsidePoint(x0, y0, x1, y1, x2, y2,  kCFThickness, z, rmin1);
  coneinsertshape->Z(3)     = z;
  coneinsertshape->Rmin(3)  = rmin1;
  coneinsertshape->Rmax(3)  = coneinsertshape->GetRmax(2);

  x0 = coneshape->GetZ(3); y0 = coneshape->GetRmax(3);
  x1 = coneshape->GetZ(4); y1 = coneshape->GetRmax(4);
  x2 = coneshape->GetZ(5); y2 = coneshape->GetRmax(5);
  InsidePoint(x0, y0, x1, y1, x2, y2, -kCFThickness, z, rmax);
  coneinsertshape->Z(4)     = z;
  coneinsertshape->Rmax(4)  = rmax;

  x0 = coneshape->GetZ(4); y0 = coneshape->GetRmax(4);
  x1 = coneshape->GetZ(5); y1 = coneshape->GetRmax(5);
  x2 = coneshape->GetZ(6); y2 = coneshape->GetRmax(6);
  InsidePoint(x0, y0, x1, y1, x2, y2, -kCFThickness, z, rmax);
  coneinsertshape->Z(5)     = z;
  coneinsertshape->Rmax(5)  = rmax;

  x0 = coneshape->GetZ(5); y0 = coneshape->GetRmax(5);
  x1 = coneshape->GetZ(6); y1 = coneshape->GetRmax(6);
  x2 = coneshape->GetZ(7); y2 = coneshape->GetRmax(7);
  InsidePoint(x0, y0, x1, y1, x2, y2, -kCFThickness, z, rmax);
  coneinsertshape->Z(6)     = z;
  coneinsertshape->Rmax(6)  = rmax;

  x0 = coneshape->GetZ(6); y0 = coneshape->GetRmin(6);
  x1 = coneshape->GetZ(7); y1 = coneshape->GetRmin(7);
  x2 = coneshape->GetZ(8); y2 = coneshape->GetRmin(8);
  InsidePoint(x0, y0, x1, y1, x2, y2,  kCFThickness, z, rmin1);
  coneinsertshape->Z(7)     = z;
  coneinsertshape->Rmin(7)  = rmin1;

  coneinsertshape->Rmin(4)  = RminFrom2Points(coneinsertshape,3,7,
					      coneinsertshape->GetZ(4));

  coneinsertshape->Rmin(5)  = RminFrom2Points(coneinsertshape,3,7,
					      coneinsertshape->GetZ(5));

  coneinsertshape->Rmin(6)  = coneinsertshape->GetRmin(5);

  x0 = coneshape->GetZ(7); y0 = coneshape->GetRmin(7);
  x1 = coneshape->GetZ(8); y1 = coneshape->GetRmin(8);
  x2 = coneshape->GetZ(9); y2 = coneshape->GetRmin(9);
  InsidePoint(x0, y0, x1, y1, x2, y2,  kCFThickness, z, rmin1);
  coneinsertshape->Z(8)     = z;
  coneinsertshape->Rmin(8)  = rmin1;

  x0 = coneshape->GetZ( 8); y0 = coneshape->GetRmin( 8);
  x1 = coneshape->GetZ( 9); y1 = coneshape->GetRmin( 9);
  x2 = coneshape->GetZ(10); y2 = coneshape->GetRmin(10);
  InsidePoint(x0, y0, x1, y1, x2, y2,  kCFThickness, z, rmin1);
  coneinsertshape->Z(9)     = z;
  coneinsertshape->Rmin(9)  = rmin1;

  x0 = coneshape->GetZ( 9); y0 = coneshape->GetRmax( 9);
  x1 = coneshape->GetZ(10); y1 = coneshape->GetRmax(10);
  x2 = coneshape->GetZ(11); y2 = coneshape->GetRmax(11);
  InsidePoint(x0, y0, x1, y1, x2, y2, -kCFThickness, z, rmax);
  coneinsertshape->Z(10)    = z;
  coneinsertshape->Rmax(10) = rmax;
  coneinsertshape->Rmin(10) = coneinsertshape->GetRmin(9);

  coneinsertshape->Rmax(7)  = RmaxFrom2Points(coneinsertshape,6,10,
					      coneinsertshape->GetZ(7));

  coneinsertshape->Rmax(8)  = RmaxFrom2Points(coneinsertshape,6,10,
					      coneinsertshape->GetZ(8));

  coneinsertshape->Rmax(9)  = coneinsertshape->GetRmax(8);

  x0 = coneshape->GetZ(10); y0 = coneshape->GetRmax(10);
  x1 = coneshape->GetZ(11); y1 = coneshape->GetRmax(11);
  x2 = coneshape->GetZ(11); y2 = coneshape->GetRmin(11);
  InsidePoint(x0, y0, x1, y1, x2, y2, -kCFThickness, z, rmax);
  coneinsertshape->Z(11)    = z;
  coneinsertshape->Rmax(11) = rmax;
  coneinsertshape->Rmin(11) = coneinsertshape->GetRmin(10);

  // SSD Cone Foams: two other Pcon's
  TGeoPcon *conefoam1shape = new TGeoPcon(0.0, 360.0, 4);

  conefoam1shape->Z(0)    = coneinsertshape->GetZ(3);
  conefoam1shape->Rmin(0) = coneinsertshape->GetRmin(3);
  conefoam1shape->Rmax(0) = conefoam1shape->GetRmin(0);

  conefoam1shape->Rmax(1) = conefoam1shape->GetRmax(0);
  conefoam1shape->Z(1)    = ZFromRmaxpCone(coneinsertshape,7,90.-kConeTheta,
					   conefoam1shape->GetRmax(1));
  conefoam1shape->Rmin(1) = RminFromZpCone(coneinsertshape,3,90.-kConeTheta,
					   conefoam1shape->GetZ(1));

  Double_t t = kConeThickness - 2*kCFThickness;
  conefoam1shape->Rmin(2) = conefoam1shape->GetRmax(0) -
                           (kConeFoam1Length*kCosConeTheta - t*kSinConeTheta);
  conefoam1shape->Z(2)    = ZFromRminpCone(coneinsertshape,3,90.-kConeTheta,
					   conefoam1shape->GetRmin(2));
  conefoam1shape->Rmax(2) = RmaxFromZpCone(coneinsertshape,7,90.-kConeTheta,
					   conefoam1shape->GetZ(2));

  conefoam1shape->Rmin(3) = conefoam1shape->GetRmin(2);
  conefoam1shape->Rmax(3) = conefoam1shape->GetRmin(3);
  conefoam1shape->Z(3)    = ZFromRmaxpCone(coneinsertshape,7,90.-kConeTheta,
					   conefoam1shape->GetRmax(3));

  TGeoPcon *conefoam2shape = new TGeoPcon(0.0, 360.0, 4);

  conefoam2shape->Z(3)    = coneinsertshape->GetZ(10);
  conefoam2shape->Rmin(3) = coneinsertshape->GetRmax(10);
  conefoam2shape->Rmax(3) = conefoam2shape->GetRmin(3);

  conefoam2shape->Rmin(2) = conefoam2shape->GetRmin(3);
  conefoam2shape->Z(2)    = ZFromRminpCone(coneinsertshape,3,90.-kConeTheta,
					   conefoam2shape->GetRmin(2));
  conefoam2shape->Rmax(2) = RmaxFromZpCone(coneinsertshape,7,90.-kConeTheta,
					   conefoam2shape->GetZ(2));

  conefoam2shape->Rmin(0) = conefoam2shape->GetRmax(2) +
                           (kConeFoam2Length*kCosConeTheta - t*kSinConeTheta);
  conefoam2shape->Rmax(0) = conefoam2shape->GetRmin(0);
  conefoam2shape->Z(0)    = ZFromRminpCone(coneinsertshape,3,90.-kConeTheta,
					   conefoam2shape->GetRmin(0));

  conefoam2shape->Rmax(1) = conefoam2shape->GetRmax(0);
  conefoam2shape->Z(1)    = ZFromRmaxpCone(coneinsertshape,7,90.-kConeTheta,
					   conefoam2shape->GetRmax(1));
  conefoam2shape->Rmin(1) = RminFromZpCone(coneinsertshape,3,90.-kConeTheta,
					   conefoam2shape->GetZ(1));

  // SSD Cone Holes: Pcon's
  // A single hole volume gives an overlap with coneinsert, so
  // three contiguous volumes are created: one to be put in coneinsert
  // and two in the cone carbon fiber envelope
  Double_t holePhi;
  holePhi = (kCoolingHoleWidth/kCoolingHoleRmin)*TMath::RadToDeg();

  TGeoPcon *coolingholeshape = new TGeoPcon(-holePhi/2., holePhi, 4);

  coolingholeshape->Rmin(0) = kCoolingHoleRmin + kCoolingHoleHight;
  coolingholeshape->Rmax(0) = coolingholeshape->GetRmin(0);
  coolingholeshape->Z(0)    = ZFromRminpCone(coneinsertshape,3,90.-kConeTheta,
					     coolingholeshape->GetRmin(0));

  coolingholeshape->Rmax(1) = coolingholeshape->GetRmax(0);
  coolingholeshape->Z(1)    = ZFromRmaxpCone(coneinsertshape,7,90.-kConeTheta,
					     coolingholeshape->GetRmax(1));
  coolingholeshape->Rmin(1) = RminFromZpCone(coneinsertshape,3,90.-kConeTheta,
					     coolingholeshape->GetZ(1));

  coolingholeshape->Rmin(2) = kCoolingHoleRmin;
  coolingholeshape->Z(2)    = ZFromRminpCone(coneinsertshape,3,90.-kConeTheta,
					     coolingholeshape->GetRmin(2));
  coolingholeshape->Rmax(2) = RmaxFromZpCone(coneinsertshape,7,90.-kConeTheta,
					     coolingholeshape->GetZ(2));

  coolingholeshape->Rmin(3) = coolingholeshape->GetRmin(2);
  coolingholeshape->Rmax(3) = coolingholeshape->GetRmin(3);
  coolingholeshape->Z(3)    = ZFromRmaxpCone(coneinsertshape,7,90.-kConeTheta,
					     coolingholeshape->GetRmax(3));

  TGeoPcon *coolinghole2shape = new TGeoPcon(-holePhi/2., holePhi, 4);

  coolinghole2shape->Rmin(0) = kCoolingHoleRmin + kCoolingHoleHight;
  coolinghole2shape->Rmax(0) = coolinghole2shape->GetRmin(0);
  coolinghole2shape->Z(0)    = ZFromRminpCone(coneshape,3,90.-kConeTheta,
					      coolinghole2shape->GetRmin(0));

  coolinghole2shape->Rmax(1) = coolinghole2shape->GetRmax(0);
  coolinghole2shape->Z(1)    = coolingholeshape->GetZ(0);
  coolinghole2shape->Rmin(1) = RminFromZpCone(coneshape,3,90.-kConeTheta,
					      coolinghole2shape->GetZ(1));

  coolinghole2shape->Rmin(2) = kCoolingHoleRmin;
  coolinghole2shape->Z(2)    = ZFromRminpCone(coneshape,3,90.-kConeTheta,
					      coolinghole2shape->GetRmin(2));
  coolinghole2shape->Rmax(2) = RminFromZpCone(coneinsertshape,3,90.-kConeTheta,
					      coolinghole2shape->GetZ(2));

  coolinghole2shape->Rmin(3) = coolinghole2shape->GetRmin(2);
  coolinghole2shape->Rmax(3) = coolinghole2shape->GetRmin(3);
  coolinghole2shape->Z(3)    = coolingholeshape->GetZ(2);

  TGeoPcon *coolinghole3shape = new TGeoPcon(-holePhi/2., holePhi, 4);

  coolinghole3shape->Rmin(0) = kCoolingHoleRmin + kCoolingHoleHight;
  coolinghole3shape->Rmax(0) = coolinghole3shape->GetRmin(0);
  coolinghole3shape->Z(0)    = coolingholeshape->GetZ(1);

  coolinghole3shape->Rmax(1) = coolinghole3shape->GetRmax(0);
  coolinghole3shape->Z(1)    = ZFromRmaxpCone(coneshape,7,90.-kConeTheta,
					      coolinghole3shape->GetRmax(1));
  coolinghole3shape->Rmin(1) = RmaxFromZpCone(coneinsertshape,7,90.-kConeTheta,
					      coolinghole3shape->GetZ(1));

  coolinghole3shape->Rmin(2) = kCoolingHoleRmin;
  coolinghole3shape->Z(2)    = coolingholeshape->GetZ(3);
  coolinghole3shape->Rmax(2) = RmaxFromZpCone(coneshape,7,90.-kConeTheta,
					      coolinghole3shape->GetZ(2));

  coolinghole3shape->Rmin(3) = coolinghole3shape->GetRmin(2);
  coolinghole3shape->Rmax(3) = coolinghole3shape->GetRmin(3);
  coolinghole3shape->Z(3)    = ZFromRmaxpCone(coneshape,7,90.-kConeTheta,
					      coolinghole3shape->GetRmax(3));

  //
  holePhi = (kMountingHoleWidth/kMountingHoleRmin)*TMath::RadToDeg();

  TGeoPcon *mountingholeshape = new TGeoPcon(-holePhi/2., holePhi, 4);

  mountingholeshape->Rmin(0) = kMountingHoleRmin + kMountingHoleHight;
  mountingholeshape->Rmax(0) = mountingholeshape->GetRmin(0);
  mountingholeshape->Z(0)    = ZFromRminpCone(coneinsertshape,3,90.-kConeTheta,
					      mountingholeshape->GetRmin(0));

  mountingholeshape->Rmin(1) = kMountingHoleRmin;
  mountingholeshape->Rmax(1) = mountingholeshape->GetRmax(0);
  mountingholeshape->Z(1)    = ZFromRminpCone(coneinsertshape,3,90.-kConeTheta,
					      mountingholeshape->GetRmin(1));

  mountingholeshape->Rmin(2) = mountingholeshape->GetRmin(1);
  mountingholeshape->Rmax(2) = mountingholeshape->GetRmax(1);
  mountingholeshape->Z(2)    = ZFromRmaxpCone(coneinsertshape,7,90.-kConeTheta,
					      mountingholeshape->GetRmax(2));

  mountingholeshape->Rmin(3) = mountingholeshape->GetRmin(2);
  mountingholeshape->Rmax(3) = mountingholeshape->GetRmin(3);
  mountingholeshape->Z(3)    = ZFromRmaxpCone(coneinsertshape,7,90.-kConeTheta,
					      mountingholeshape->GetRmax(3));

  TGeoPcon *mountinghole2shape = new TGeoPcon(-holePhi/2., holePhi, 4);

  mountinghole2shape->Rmin(0) = kMountingHoleRmin + kMountingHoleHight;
  mountinghole2shape->Rmax(0) = mountingholeshape->GetRmin(0);
  mountinghole2shape->Z(0)    = ZFromRminpCone(coneshape,3,90.-kConeTheta,
					       mountinghole2shape->GetRmin(0));

  mountinghole2shape->Rmax(1) = mountinghole2shape->GetRmax(0);
  mountinghole2shape->Z(1)    = mountingholeshape->Z(0);
  mountinghole2shape->Rmin(1) = RminFromZpCone(coneshape,3,90.-kConeTheta,
					       mountinghole2shape->GetZ(1));

  mountinghole2shape->Rmin(2) = kMountingHoleRmin;
  mountinghole2shape->Z(2)    = ZFromRminpCone(coneshape,3,90.-kConeTheta,
					       mountinghole2shape->GetRmin(2));
  mountinghole2shape->Rmax(2) = RminFromZpCone(coneinsertshape,3,90.-kConeTheta,
					       mountinghole2shape->GetZ(2));

  mountinghole2shape->Rmin(3) = mountinghole2shape->Rmin(2);
  mountinghole2shape->Rmax(3) = mountinghole2shape->Rmin(3);
  mountinghole2shape->Z(3)    = mountingholeshape->Z(1);

  TGeoPcon *mountinghole3shape = new TGeoPcon(-holePhi/2., holePhi, 4);

  mountinghole3shape->Rmin(0) = kMountingHoleRmin + kMountingHoleHight;
  mountinghole3shape->Rmax(0) = mountingholeshape->GetRmin(0);
  mountinghole3shape->Z(0)    = mountingholeshape->GetZ(2);

  mountinghole3shape->Rmax(1) = mountinghole3shape->GetRmax(0);
  mountinghole3shape->Z(1)    = ZFromRmaxpCone(coneshape,7,90.-kConeTheta,
					       mountinghole3shape->GetRmax(1));
  mountinghole3shape->Rmin(1) = RmaxFromZpCone(coneinsertshape,7,90.-kConeTheta,
					       mountinghole3shape->GetZ(1));

  mountinghole3shape->Rmin(2) = kMountingHoleRmin;
  mountinghole3shape->Z(2)    = mountingholeshape->Z(3);
  mountinghole3shape->Rmax(2) = RmaxFromZpCone(coneshape,7,90.-kConeTheta,
					       mountinghole3shape->GetZ(2));

  mountinghole3shape->Rmin(3) = mountinghole3shape->Rmin(2);
  mountinghole3shape->Rmax(3) = mountinghole3shape->Rmin(3);
  mountinghole3shape->Z(3)    = ZFromRmaxpCone(coneshape,7,90.-kConeTheta,
					       mountinghole3shape->GetRmax(3));

  // The Cable Hole is even more complicated, a Composite Shape
  // is unavoidable here (gosh!)
  TGeoPcon *coneshapecopy = new TGeoPcon("conecopy",0.0, 360.0, 12);

  for (Int_t i=0; i<12; i++) {
    coneshapecopy->Rmin(i) = coneshape->GetRmin(i);
    coneshapecopy->Rmax(i) = coneshape->GetRmax(i);
    coneshapecopy->Z(i)    = coneshape->GetZ(i);
  }

  holePhi = (kCableHoleWidth/kCableHoleRout)*TMath::RadToDeg();
  TGeoConeSeg *chCS = new TGeoConeSeg("chCS", 0.5*kConeZLength,
				      kCableHoleRin, kCableHoleRout,
				      kCableHoleRin, kCableHoleRout,
				      -0.5*holePhi, 0.5*holePhi);

  TGeoCompositeShape *cableholeshape = new TGeoCompositeShape(
					   "SSDCableHoleShape",
					   "conecopy*chCS");

  if(GetDebug(1)){
    chCS->InspectShape();
    cableholeshape->InspectShape();
  }

  // SSD Cone Wings: Tube and TubeSeg shapes
  Double_t angleWideWing, angleWideWingThickness;
  angleWideWing = (kWingWidth/kWingRmax)*TMath::RadToDeg();
  angleWideWingThickness = (kCFThickness/kWingRmax)*TMath::RadToDeg();

  TGeoTubeSeg *wingshape = new TGeoTubeSeg(kConeROuterMax, kWingRmax,
					   kWingHalfThick,
					   0, angleWideWing);

  TGeoTubeSeg *winginsertshape = new TGeoTubeSeg(kConeROuterMax,
				 kWingRmax-kCFThickness,
				 kWingHalfThick-kCFThickness,
				 angleWideWingThickness,
				 angleWideWing-angleWideWingThickness);

  // SDD support plate, SSD side (Mounting Bracket): a TubeSeg
  TGeoTubeSeg *bracketshape = new TGeoTubeSeg(kBracketRmin, kBracketRmax,
                            kBracketHalfLength, -kBracketPhi/2, kBracketPhi/2);


  // We have the shapes: now create the real volumes

  TGeoVolume *cfcone = new TGeoVolume("SSDCarbonFiberCone",
				      coneshape,medSSDcf);
  cfcone->SetVisibility(kTRUE);
  cfcone->SetLineColor(4); // Blue
  cfcone->SetLineWidth(1);
  cfcone->SetFillColor(cfcone->GetLineColor());
  cfcone->SetFillStyle(4000); // 0% transparent

  TGeoVolume *cfconeinsert = new TGeoVolume("SSDCarbonFiberConeInsert",
					    coneinsertshape,medSSDste);
  cfconeinsert->SetVisibility(kTRUE);
  cfconeinsert->SetLineColor(2); // Red
  cfconeinsert->SetLineWidth(1);
  cfconeinsert->SetFillColor(cfconeinsert->GetLineColor());
  cfconeinsert->SetFillStyle(4050); // 50% transparent

  TGeoVolume *cfconefoam1 = new TGeoVolume("SSDCarbonFiberConeFoam1",
					    conefoam1shape,medSSDroh);
  cfconefoam1->SetVisibility(kTRUE);
  cfconefoam1->SetLineColor(3); // Green
  cfconefoam1->SetLineWidth(1);
  cfconefoam1->SetFillColor(cfconefoam1->GetLineColor());
  cfconefoam1->SetFillStyle(4050); // 50% transparent

  TGeoVolume *cfconefoam2 = new TGeoVolume("SSDCarbonFiberConeFoam2",
					    conefoam2shape,medSSDroh);
  cfconefoam2->SetVisibility(kTRUE);
  cfconefoam2->SetLineColor(3); // Green
  cfconefoam2->SetLineWidth(1);
  cfconefoam2->SetFillColor(cfconefoam2->GetLineColor());
  cfconefoam2->SetFillStyle(4050); // 50% transparent

  TGeoVolume *coolinghole = new TGeoVolume("SSDCoolingHole",
					   coolingholeshape,medSSDair);
  coolinghole->SetVisibility(kTRUE);
  coolinghole->SetLineColor(5); // Yellow
  coolinghole->SetLineWidth(1);
  coolinghole->SetFillColor(coolinghole->GetLineColor());
  coolinghole->SetFillStyle(4090); // 90% transparent

  TGeoVolume *coolinghole2 = new TGeoVolume("SSDCoolingHole2",
					    coolinghole2shape,medSSDair);
  coolinghole2->SetVisibility(kTRUE);
  coolinghole2->SetLineColor(5); // Yellow
  coolinghole2->SetLineWidth(1);
  coolinghole2->SetFillColor(coolinghole2->GetLineColor());
  coolinghole2->SetFillStyle(4090); // 90% transparent

  TGeoVolume *coolinghole3 = new TGeoVolume("SSDCoolingHole3",
					    coolinghole3shape,medSSDair);
  coolinghole3->SetVisibility(kTRUE);
  coolinghole3->SetLineColor(5); // Yellow
  coolinghole3->SetLineWidth(1);
  coolinghole3->SetFillColor(coolinghole3->GetLineColor());
  coolinghole3->SetFillStyle(4090); // 90% transparent

  TGeoVolume *mountinghole = new TGeoVolume("SSDMountingHole",
					    mountingholeshape,medSSDair);
  mountinghole->SetVisibility(kTRUE);
  mountinghole->SetLineColor(5); // Yellow
  mountinghole->SetLineWidth(1);
  mountinghole->SetFillColor(mountinghole->GetLineColor());
  mountinghole->SetFillStyle(4090); // 90% transparent

  TGeoVolume *mountinghole2 = new TGeoVolume("SSDMountingHole2",
					     mountinghole2shape,medSSDair);
  mountinghole2->SetVisibility(kTRUE);
  mountinghole2->SetLineColor(5); // Yellow
  mountinghole2->SetLineWidth(1);
  mountinghole2->SetFillColor(mountinghole2->GetLineColor());
  mountinghole2->SetFillStyle(4090); // 90% transparent

  TGeoVolume *mountinghole3 = new TGeoVolume("SSDMountingHole3",
					     mountinghole3shape,medSSDair);
  mountinghole3->SetVisibility(kTRUE);
  mountinghole3->SetLineColor(5); // Yellow
  mountinghole3->SetLineWidth(1);
  mountinghole3->SetFillColor(mountinghole3->GetLineColor());
  mountinghole3->SetFillStyle(4090); // 90% transparent

  TGeoVolume *wing = new TGeoVolume("SSDWing",wingshape,medSSDcf);
  wing->SetVisibility(kTRUE);
  wing->SetLineColor(4); // Blue
  wing->SetLineWidth(1);
  wing->SetFillColor(wing->GetLineColor());
  wing->SetFillStyle(4000); // 0% transparent

  TGeoVolume *cablehole = new TGeoVolume("SSDCableHole",
					 cableholeshape,medSSDair);
  cablehole->SetVisibility(kTRUE);
  cablehole->SetLineColor(5); // Yellow
  cablehole->SetLineWidth(1);
  cablehole->SetFillColor(cablehole->GetLineColor());
  cablehole->SetFillStyle(4090); // 90% transparent

  TGeoVolume *winginsert = new TGeoVolume("SSDWingInsert",
					  winginsertshape,medSSDste);
  winginsert->SetVisibility(kTRUE);
  winginsert->SetLineColor(2); // Red
  winginsert->SetLineWidth(1);
  winginsert->SetFillColor(winginsert->GetLineColor());
  winginsert->SetFillStyle(4050); // 50% transparent

  TGeoVolume *bracket = new TGeoVolume("SSDMountingBracket",
				       bracketshape,medSSDal);
  bracket->SetVisibility(kTRUE);
  bracket->SetLineColor(6); // Purple
  bracket->SetLineWidth(1);
  bracket->SetFillColor(bracket->GetLineColor());
  bracket->SetFillStyle(4000); // 0% transparent

  // Mount up a cone
  for (Int_t i=0; i<(Int_t)(360./kMountingHolePhi); i++) {
    Double_t phiH = i*kMountingHolePhi + 0.5*kMountingHolePhi;
    cfconefoam2->AddNode(mountinghole,i+1, new TGeoRotation("", phiH, 0, 0));
  }

  for (Int_t i=0; i<(Int_t)(360./kCoolingHolePhi); i++) {
    Double_t phiH = i*kCoolingHolePhi + 0.5*kCoolingHolePhi;
    cfconeinsert->AddNodeOverlap(coolinghole,i+1, new TGeoRotation("", phiH, 0, 0));
  }

  cfconeinsert->AddNode(cfconefoam1,1,0);
  cfconeinsert->AddNode(cfconefoam2,1,0);

  cfcone->AddNode(cfconeinsert,1,0);

  for (Int_t i=0; i<(Int_t)(360./kCoolingHolePhi); i++) {
    Double_t phiH = i*kCoolingHolePhi + 0.5*kCoolingHolePhi;
    cfcone->AddNode(coolinghole2,i+1, new TGeoRotation("", phiH, 0, 0));
    cfcone->AddNode(coolinghole3,i+1, new TGeoRotation("", phiH, 0, 0));
    cfcone->AddNodeOverlap(cablehole,i+1, new TGeoRotation("", phiH, 0, 0));
  }

  for (Int_t i=0; i<(Int_t)(360./kMountingHolePhi); i++) {
    Double_t phiH = i*kMountingHolePhi + 0.5*kMountingHolePhi;
    cfcone->AddNode(mountinghole2,i+1, new TGeoRotation("", phiH, 0, 0));
    cfcone->AddNode(mountinghole3,i+1, new TGeoRotation("", phiH, 0, 0));
  }

  wing->AddNode(winginsert,1,0);

  // Add all volumes in the Cone assembly
  vC->AddNode(cfcone,1,new TGeoTranslation(0,0,-kConeZPosition));

  for (Int_t i=0; i<4; i++) {
    Double_t thetaW = kThetaWing + 90.*i + angleWideWing/2.;
    vC->AddNode(wing, i+1, new TGeoCombiTrans(0, 0, -kConeZPosition+kWingHalfThick,
			   new TGeoRotation("",thetaW,180,0)));
  }

  Double_t zBracket = kConeZPosition - coneshape->GetZ(9) +
                      2*bracketshape->GetDz();
  for (Int_t i=0; i<3; i++) {
    Double_t thetaB = 60 + 120.*i;
    vC->AddNode(bracket, i+1, new TGeoCombiTrans(0, 0, -zBracket,
			      new TGeoRotation("",thetaB,0,0)));
  }

  // Finally put everything in the mother volume
  moth->AddNode(cfcylinder,1,0);

  moth->AddNode(vC, 1, 0 );
  moth->AddNode(vC, 2, new TGeoRotation("",180, 180, 0) );

  // Some debugging if requested
  if(GetDebug(1)){
    vC->PrintNodes();
    vC->InspectShape();
  }

  return;
}

//______________________________________________________________________
void AliITSv11GeometrySupport::ServicesCableSupport(TGeoVolume *moth,
                                                    TGeoManager *mgr){
//
// Creates the cable trays which are outside the ITS support cones
// but still inside the TPC
// This is now a stearing routine, the actual work is done by three
// specialized methods to avoid a really huge unique method
//
// Input:
//         moth : the TGeoVolume owing the volume structure
//         mgr  : the GeoManager (default gGeoManager)
// Output:
//
// Created:      15 Nov 2009  Mario Sitta
//

  TraySupportsSideA(moth, mgr);

  ServicesCableSupportSPD(moth, mgr);
  ServicesCableSupportSDD(moth, mgr);
  ServicesCableSupportSSD(moth, mgr);

  return;
}

//______________________________________________________________________
void AliITSv11GeometrySupport::TraySupportsSideA(TGeoVolume *moth,
					   const TGeoManager *mgr){
//
// Creates the structure supporting the ITS cable trays on Side A
//
// Input:
//         moth : the TGeoVolume owing the volume structure
//         mgr  : the GeoManager (default gGeoManager)
// Output:
//
// Created:      14 Dec 2009  Mario Sitta
// Updated:      26 Feb 2010  Mario Sitta
//
// Technical data are taken from AutoCAD drawings, L.Simonetti technical
// drawings and other (oral) information given by F.Tosello
//

  // Dimensions and positions of the A-Side Cable Tray Support Ring
  // (0872/G/A/01)
  const Double_t kSuppRingYTrans      =  110.00 *fgkmm;
  const Double_t kSuppRingZTrans      =(1011.00+435.00) *fgkmm;
  const Double_t kSuppForwYTrans      =  185.00 *fgkmm;

  const Double_t kExtSuppRingSpace1   =   33.00 *fgkmm;
  const Double_t kExtSuppRingSpace2   =   45.00 *fgkmm;
  const Double_t kExtSuppRingSpcAbov  =   30.00 *fgkmm;
  const Double_t kExtSuppRingBase     =  491.50 *fgkmm;
  const Double_t kExtSuppRingInward   =   35.00 *fgkmm;
  const Double_t kExtSuppRingRmax     =  540.00 *fgkmm;
  const Double_t kExtSuppRingRint1    =  465.00 *fgkmm;
  const Double_t kExtSuppRingRint2    =  467.00 *fgkmm;
  const Double_t kExtSuppRingInnerHi  =  450.00 *fgkmm;
  const Double_t kExtSuppRingInWide   =  100.00 *fgkmm;
  const Double_t kExtSuppRingR7       =    7.00 *fgkmm;
  const Double_t kExtSuppRingR5       =    5.00 *fgkmm;
  const Double_t kExtSuppRingThick    =   20.00 *fgkmm;

  const Double_t kExtSuppRingSpcAng   =   10.50 *TMath::DegToRad();
  const Double_t kExtSuppRingPartPhi  =   15.00 *TMath::DegToRad();
  const Double_t kExtSuppRingIntAng   =    7.00 *TMath::DegToRad();
  const Double_t kExtSuppRingBaseAng  =   75.00 *TMath::DegToRad();
  const Double_t kExtSuppRingR7Ang    =  100.00 *TMath::DegToRad(); // Guessed

  const Int_t    kExtSuppRingNPtsArc  =   10; // N.points to approximate arc

  const Double_t kIntSuppRingThick1   =   15.00 *fgkmm;
  const Double_t kIntSuppRingThick2   =   13.00 *fgkmm;
  const Double_t kIntSuppRingInward   =   24.00 *fgkmm;
  const Double_t kIntSuppRingThick    =   20.00 *fgkmm;

  const Double_t kSuppCylHeight       =  340.00 *fgkmm;
  const Double_t kSuppCylRint         =  475.00 *fgkmm;
  const Double_t kSuppCylRext         =  478.00 *fgkmm;
  const Double_t kSuppCylDispl        =  137.70 *fgkmm;

  const Double_t kSuppSpacerHeight    =   30.00 *fgkmm;
  const Double_t kSuppSpacerThick     =   10.00 *fgkmm;

  const Double_t kSuppSpacerAngle     =   15.00;  // Degrees

  const Double_t kSuppForwRingRint1   =  500.00 *fgkmm;
  const Double_t kSuppForwRingRint2   =  540.00 *fgkmm;
  const Double_t kSuppForwRingRext    =  560.00 *fgkmm;
  const Double_t kSuppForwRingThikAll =   50.00 *fgkmm;
  const Double_t kSuppForwRingThikInt =   20.00 *fgkmm;

  // (0872/G/B/01)
  const Double_t kSuppForwConeRmin    =  558.00 *fgkmm;
  const Double_t kSuppForwConeRmax    =  681.00 *fgkmm;
  const Double_t kSuppForwConeLen1    =  318.00 *fgkmm;
  const Double_t kSuppForwConeLen2    =  662.00 *fgkmm;
  const Double_t kSuppForwConeThick   =    3.00 *fgkmm;

  const Double_t kSuppBackRingPlacTop =   90.00 *fgkmm;
  const Double_t kSuppBackRingPlacSid =   50.00 *fgkmm;
  const Double_t kSuppBackRingHeight  =  760.00 *fgkmm;
  const Double_t kSuppBackRingRext    =  760.00 *fgkmm;
  const Double_t kSuppBackRingRint    =  685.00 *fgkmm;
//  const Double_t kSuppBackRingRint2   =  675.00 *fgkmm;
  const Double_t kSuppBackRingR10     =   10.00 *fgkmm;
  const Double_t kSuppBackRingBase    =  739.00 *fgkmm;
  const Double_t kSuppBackRingThikAll =   50.00 *fgkmm;
  const Double_t kSuppBackRingThick1  =   20.00 *fgkmm;
  const Double_t kSuppBackRingThick2  =   20.00 *fgkmm;

//  const Double_t kSuppBackRingPlacAng =   10.00 *TMath::DegToRad();
  const Double_t kSuppBackRingPlacAng =   10.25 *TMath::DegToRad();//Fix ovlp.
  const Double_t kSuppBackRing2ndAng1 =   78.40 *TMath::DegToRad();
  const Double_t kSuppBackRing2ndAng2 =   45.00 *TMath::DegToRad();

  const Int_t    kSuppBackRingNPtsArc =   10; // N.points to approximate arc

  // (0872/G/C/01)
  const Double_t kRearSuppZTransGlob  =(1011.00+9315.00-6040.00) *fgkmm;
  const Double_t kBackRodZTrans       = 2420.00 *fgkmm;

  const Double_t kBackRodLength       = 1160.00 *fgkmm;
  const Double_t kBackRodThickLen     =   20.00 *fgkmm;
  const Double_t kBackRodDiameter     =   20.00 *fgkmm;

  const Double_t kSuppRearRingRint    =  360.00 *fgkmm;
  const Double_t kSuppRearRingRext1   =  410.00 *fgkmm;
  const Double_t kSuppRearRingRext2   =  414.00 *fgkmm;
  const Double_t kSuppRearRingHeight  =  397.00 *fgkmm;
  const Double_t kSuppRearRingTopWide =  111.87 *fgkmm;
  const Double_t kSuppRearRingBase    =  451.50 *fgkmm;
  const Double_t kSuppRearRingBaseHi  =   58.00 *fgkmm;
  const Double_t kSuppRearRingSideHi  =   52.00 *fgkmm;
  const Double_t kSuppRearRingInside  =   40.00 *fgkmm;
  const Double_t kSuppRearRingInsideHi=   12.00 *fgkmm;
  const Double_t kSuppRearRingThick   =   20.00 *fgkmm;
  const Double_t kSuppRearRingXRodHole=  441.50 *fgkmm;
  const Double_t kSuppRearRingYRodHole=   42.00 *fgkmm;

  const Double_t kSuppRearRing1stAng  =   22.00 *TMath::DegToRad();
  const Double_t kSuppRearRingStepAng =   15.00 *TMath::DegToRad();

  const Int_t    kSuppRearRingNPtsArc =   10; // N.points to approximate arc


  // Local variables
  Double_t xprof[2*(15+kExtSuppRingNPtsArc)],yprof[2*(15+kExtSuppRingNPtsArc)];
  Double_t slp1, slp2, phi, xm, ym;
  Double_t xloc, yloc, zloc, rmin, rmax, deltaR;
  Int_t npoints;


  // The whole support as an assembly
  TGeoVolumeAssembly *trayASuppStruct = new TGeoVolumeAssembly("ITSsuppSideAStructure");
  

  // First create all needed shapes

  // The External Ring (part of 0872/G/A/01): a really complex Xtru
  TGeoXtru *extSuppRing = new TGeoXtru(2);

  // First the upper notch...
  xprof[ 0] = kExtSuppRingSpace1;
  yprof[ 0] = kExtSuppRingInnerHi + kExtSuppRingSpcAbov;

  slp1 = TMath::Tan(TMath::Pi()/2 - kExtSuppRingSpcAng);
  IntersectCircle(slp1, xprof[0], yprof[0], kExtSuppRingRmax, 0., 0.,
		  xprof[5], yprof[5], xm, ym); // Ignore dummy xm,ym

  xprof[ 4] = xprof[5];
  yprof[ 4] = yprof[5] - kExtSuppRingR5/TMath::Tan(kExtSuppRingSpcAng);
  xprof[ 3] = xprof[4] - kExtSuppRingR5*(1 - TMath::Cos(TMath::Pi()/6));
  yprof[ 3] = yprof[4] - kExtSuppRingR5*(    TMath::Sin(TMath::Pi()/6));
  xprof[ 2] = xprof[4] - kExtSuppRingR5*(1 - TMath::Cos(TMath::Pi()/3));
  yprof[ 2] = yprof[4] - kExtSuppRingR5*(    TMath::Sin(TMath::Pi()/3));
  xprof[ 1] = xprof[4] - kExtSuppRingR5;
  yprof[ 1] = yprof[4] - kExtSuppRingR5;

  Int_t indx = 5+kExtSuppRingNPtsArc;
  // ...then the external arc, approximated with segments,...
  xprof[indx] = kExtSuppRingBase;
  yprof[indx] = TMath::Sqrt(kExtSuppRingRmax*kExtSuppRingRmax -
			    kExtSuppRingBase*kExtSuppRingBase);
  Double_t alphamin = TMath::ASin(kExtSuppRingSpace2/kExtSuppRingRmax);
  Double_t alphamax = TMath::Pi()/2 -
		    TMath::ASin(yprof[5+kExtSuppRingNPtsArc]/kExtSuppRingRmax);

  for (Int_t jp = 1; jp < kExtSuppRingNPtsArc; jp++) {
    Double_t alpha = jp*(alphamax-alphamin)/kExtSuppRingNPtsArc;
    xprof[5+jp] = kExtSuppRingRmax*TMath::Sin(alpha);
    yprof[5+jp] = kExtSuppRingRmax*TMath::Cos(alpha);
  }
  // ...and finally the interior profile
  xprof[indx+1] = kExtSuppRingBase;
  yprof[indx+1] = kSuppRingYTrans;
  xprof[indx+2] = xprof[indx+1] - kExtSuppRingInward;
  yprof[indx+2] = yprof[indx+1];

  phi  = TMath::Pi()/2 - 4*kExtSuppRingPartPhi - kExtSuppRingIntAng;
  slp1 = TMath::Tan(TMath::Pi() - kExtSuppRingBaseAng);
  slp2 = TMath::Tan(TMath::Pi()/2 + phi);
  xm   = kExtSuppRingRint2*TMath::Cos(phi);
  ym   = kExtSuppRingRint2*TMath::Sin(phi);
  IntersectLines(slp1, xprof[indx+2], yprof[indx+2], slp2, xm, ym,
		 xprof[indx+3], yprof[indx+3]);

  slp1 = slp2;
  phi += kExtSuppRingPartPhi;
  slp2 = TMath::Tan(TMath::Pi()/2 + phi);
  xm   = kExtSuppRingRint1*TMath::Cos(phi);
  ym   = kExtSuppRingRint1*TMath::Sin(phi);
  IntersectLines(slp1, xprof[indx+3], yprof[indx+3], slp2, xm, ym,
		 xprof[indx+4], yprof[indx+4]);
  
  slp1 = slp2;
  phi += kExtSuppRingPartPhi;
  slp2 = TMath::Tan(TMath::Pi()/2 + phi);
  xm   = kExtSuppRingRint2*TMath::Cos(phi);
  ym   = kExtSuppRingRint2*TMath::Sin(phi);
  IntersectLines(slp1, xprof[indx+4], yprof[indx+4], slp2, xm, ym,
		 xprof[indx+5], yprof[indx+5]);
  
  slp1 = slp2;
  phi += kExtSuppRingPartPhi;
  slp2 = TMath::Tan(TMath::Pi()/2 + phi);
  xm   = kExtSuppRingRint1*TMath::Cos(phi);
  ym   = kExtSuppRingRint1*TMath::Sin(phi);
  IntersectLines(slp1, xprof[indx+5], yprof[indx+5], slp2, xm, ym,
		 xprof[indx+6], yprof[indx+6]);
  
  xprof[indx+9] = kExtSuppRingInWide;
  yprof[indx+9] = kExtSuppRingInnerHi;
  xprof[indx+8] = xprof[indx+9] +
		  (1 - TMath::Cos(kExtSuppRingR7Ang/2))*kExtSuppRingR7;
  yprof[indx+8] = yprof[indx+9] +
		  (    TMath::Sin(kExtSuppRingR7Ang/2))*kExtSuppRingR7;
  xprof[indx+7] = xprof[indx+9] +
		  (1 + TMath::Cos(kExtSuppRingR7Ang  ))*kExtSuppRingR7;
  yprof[indx+7] = yprof[indx+9] +
		  (    TMath::Sin(kExtSuppRingR7Ang  ))*kExtSuppRingR7;
  // Gosh, we did the right side! now reflex on the left side
  npoints = (sizeof(xprof)/sizeof(Double_t))/2;
  for (Int_t jp = 0; jp < npoints; jp++) {
    xprof[npoints+jp] = -xprof[npoints-1-jp];
    yprof[npoints+jp] =  yprof[npoints-1-jp];
  }
  // wow! now the actual Xtru
  extSuppRing->DefinePolygon(2*npoints, xprof, yprof);
  extSuppRing->DefineSection(0,0);
  extSuppRing->DefineSection(1,kExtSuppRingThick);

  // The Internal Ring (part of 0872/G/A/01): another complex Xtru
  TGeoXtru *intSuppRing = new TGeoXtru(2);

  // First the external profile...
  npoints = 0;

  slp1 = 0;
  phi  = TMath::Pi()/2 - kExtSuppRingPartPhi - kExtSuppRingIntAng;
  slp2 = TMath::Tan(TMath::Pi()/2 + phi);
  xm   = (kExtSuppRingRint1+kIntSuppRingThick1)*TMath::Cos(phi);
  ym   = (kExtSuppRingRint1+kIntSuppRingThick1)*TMath::Sin(phi);
  IntersectLines(slp1,  0, kExtSuppRingInnerHi+kExtSuppRingSpcAbov,
		 slp2, xm, ym,
		 xprof[npoints], yprof[npoints]);
  npoints++;

  slp1 = slp2;
  phi -= kExtSuppRingPartPhi;
  slp2 = TMath::Tan(TMath::Pi()/2 + phi);
  xm   = (kExtSuppRingRint2+kIntSuppRingThick2)*TMath::Cos(phi);
  ym   = (kExtSuppRingRint2+kIntSuppRingThick2)*TMath::Sin(phi);
  IntersectLines(slp1, xprof[npoints-1], yprof[npoints-1],
		 slp2, xm, ym,
		 xprof[npoints], yprof[npoints]);
  npoints++;

  slp1 = slp2;
  phi -= kExtSuppRingPartPhi;
  slp2 = TMath::Tan(TMath::Pi()/2 + phi);
  xm   = (kExtSuppRingRint1+kIntSuppRingThick1)*TMath::Cos(phi);
  ym   = (kExtSuppRingRint1+kIntSuppRingThick1)*TMath::Sin(phi);
  IntersectLines(slp1, xprof[npoints-1], yprof[npoints-1],
		 slp2, xm, ym,
		 xprof[npoints], yprof[npoints]);
  npoints++;

  slp1 = slp2;
  phi -= kExtSuppRingPartPhi;
  slp2 = TMath::Tan(TMath::Pi()/2 + phi);
  xm   = (kExtSuppRingRint2+kIntSuppRingThick2)*TMath::Cos(phi);
  ym   = (kExtSuppRingRint2+kIntSuppRingThick2)*TMath::Sin(phi);
  IntersectLines(slp1, xprof[npoints-1], yprof[npoints-1],
		 slp2, xm, ym,
		 xprof[npoints], yprof[npoints]);
  npoints++;

  xprof[npoints] = kExtSuppRingBase-kIntSuppRingInward;
  yprof[npoints] = Yfrom2Points(xprof[npoints-1], yprof[npoints-1], xm, ym,
				xprof[npoints]);
  npoints++;

  xprof[npoints] = xprof[npoints-1];
  yprof[npoints] = kSuppRingYTrans;
  npoints++;
  // ...and then the interior profile, which is identical to extSuppRing one
  for (Int_t jp=0; jp < 8; jp++) {
    xprof[npoints] = extSuppRing->GetX(17+jp);
    yprof[npoints] = extSuppRing->GetY(17+jp);
    npoints++;
  }
  // We did the right side! now reflex on the left side
  for (Int_t jp = 0; jp < npoints; jp++) {
    xprof[npoints+jp] = -xprof[npoints-1-jp];
    yprof[npoints+jp] =  yprof[npoints-1-jp];
  }
  // And now the actual Xtru
  intSuppRing->DefinePolygon(2*npoints, xprof, yprof);
  intSuppRing->DefineSection(0,0);
  intSuppRing->DefineSection(1,kIntSuppRingThick);

  // The intermediate cylinder (0872/G/A/03): a TubeSeg
  alphamin = TMath::ASin(kSuppCylDispl/kSuppCylRint)*TMath::RadToDeg();
  alphamax = 180 - alphamin;
  TGeoTubeSeg *interCylind = new TGeoTubeSeg(kSuppCylRint, kSuppCylRext,
				     kSuppCylHeight/2, alphamin, alphamax);

  // The spacer (0872/G/A/03): a simple Xtru
  TGeoXtru *suppSpacer = new TGeoXtru(2);

  xprof[0] = kSuppSpacerHeight;
  yprof[0] = kSuppSpacerThick;
  xprof[1] = xprof[0];
  yprof[1] = 0;
  xprof[2] = 0;
  yprof[2] = 0;
  xprof[3] = kSuppSpacerThick*SinD(kSuppSpacerAngle);
  yprof[3] = yprof[0];

  suppSpacer->DefinePolygon(4, xprof, yprof);
  suppSpacer->DefineSection(0,-kSuppCylHeight/2);
  suppSpacer->DefineSection(1, kSuppCylHeight/2);

  // The forward ring (0872/G/B/02): a Pcon (slight oversimplification)
  Double_t rmean = (kSuppForwRingRint1+kSuppForwRingRext)/2;
  alphamin = TMath::ASin(kSuppForwYTrans/rmean)*TMath::RadToDeg();
  alphamax = 180 - alphamin;

  TGeoPcon *forwardRing = new TGeoPcon(alphamin,alphamax-alphamin,4);

  forwardRing->DefineSection(0,0,
			     kSuppForwRingRint1,kSuppForwRingRext);
  forwardRing->DefineSection(1,kSuppForwRingThikInt,
			     kSuppForwRingRint1,kSuppForwRingRext);
  forwardRing->DefineSection(2,kSuppForwRingThikInt,
			     kSuppForwRingRint2,kSuppForwRingRext);
  forwardRing->DefineSection(3,kSuppForwRingThikAll,
			     kSuppForwRingRint2,kSuppForwRingRext);

  // The forward cone (0872/G/B/03): a TGeoPcon
  TGeoPcon *forwardCone = new TGeoPcon(alphamin,alphamax-alphamin,3);

  forwardCone->DefineSection(0,0,
			     kSuppForwConeRmin-kSuppForwConeThick,
			     kSuppForwConeRmin);
  forwardCone->DefineSection(1,kSuppForwConeLen1,
			     kSuppForwConeRmin-kSuppForwConeThick,
			     kSuppForwConeRmin);
  forwardCone->DefineSection(2,kSuppForwConeLen1+kSuppForwConeLen2,
			     kSuppForwConeRmax-kSuppForwConeThick,
			     kSuppForwConeRmax);

  // The first part of the Back Ring (part of 0872/G/B/01): a complex Xtru
  TGeoXtru *firstSuppBackRing = new TGeoXtru(2);

  // First the external profile... (the arc is approximated with segments)
  npoints = 0;

  xprof[npoints] = kSuppBackRingPlacTop;
  yprof[npoints] = kSuppBackRingHeight;
  npoints++;

  alphamax = TMath::Pi()/2 - TMath::ASin(kSuppBackRingPlacTop/kSuppBackRingRext);
  alphamin = TMath::ASin((kSuppForwYTrans+kSuppBackRingPlacSid)/kSuppBackRingRext);

  xprof[npoints] = xprof[npoints-1];
  yprof[npoints] = kSuppBackRingRext*TMath::Sin(alphamax);
  npoints++;

  for (Int_t jp = 1; jp <= kSuppBackRingNPtsArc; jp++) {
    Double_t alpha = alphamax - jp*(alphamax-alphamin)/kSuppBackRingNPtsArc;
    xprof[npoints] = kSuppBackRingRext*TMath::Cos(alpha);
    yprof[npoints] = kSuppBackRingRext*TMath::Sin(alpha);
    npoints++;
  }

  xprof[npoints] = kSuppBackRingBase -
		   kSuppBackRingPlacSid*TMath::Tan(kSuppBackRingPlacAng);
  yprof[npoints] = yprof[npoints-1];
  npoints++;

  xprof[npoints] = kSuppBackRingBase;
  yprof[npoints] = kSuppForwYTrans;
  npoints++;
  // ...then the internal profile (the arc is approximated with segments)
  alphamin = TMath::ASin(kSuppForwYTrans/kSuppBackRingRint);
  alphamax = TMath::Pi()/2;

  for (Int_t jp = 0; jp < kSuppBackRingNPtsArc; jp++) {
    Double_t alpha = alphamin + jp*(alphamax-alphamin)/kSuppBackRingNPtsArc;
    xprof[npoints] = kSuppBackRingRint*TMath::Cos(alpha);
    yprof[npoints] = kSuppBackRingRint*TMath::Sin(alpha);
    npoints++;
  }

  xprof[npoints] = 0;
  yprof[npoints] = kSuppBackRingRint;
  npoints++;
  // We did the right side! now reflex on the left side (except last point)
  for (Int_t jp = 0; jp < npoints-1; jp++) {
    xprof[npoints+jp] = -xprof[npoints-jp-2];
    yprof[npoints+jp] =  yprof[npoints-jp-2];
  }
  // And now the actual Xtru
  firstSuppBackRing->DefinePolygon(2*npoints-1, xprof, yprof);
  firstSuppBackRing->DefineSection(0,0);
  firstSuppBackRing->DefineSection(1,kSuppBackRingThick1);

  // The second part of the Back Ring (part of 0872/G/B/01): a Pcon
  // (slight oversimplification)
  alphamin = TMath::ASin(kSuppForwYTrans/kSuppBackRingRint)*TMath::RadToDeg();
  alphamax = 180 - alphamin;

  TGeoPcon *secondSuppBackRing = new TGeoPcon(alphamin,alphamax-alphamin,6);

  deltaR = kSuppBackRingThick2/TMath::Sin(kSuppBackRing2ndAng1);
  rmin = kSuppBackRingRint - kSuppBackRingThick1/TMath::Tan(kSuppBackRing2ndAng1);
  rmax = rmin + deltaR + kSuppBackRingR10*TMath::Sin(kSuppBackRing2ndAng1);
  secondSuppBackRing->DefineSection(0, 0, rmin, rmax);

  zloc = kSuppBackRingR10*(1 - TMath::Cos(kSuppBackRing2ndAng1/3));
  rmax -= kSuppBackRingR10*TMath::Sin(kSuppBackRing2ndAng1/3);
  rmin = secondSuppBackRing->GetRmin(0) - zloc/TMath::Tan(kSuppBackRing2ndAng1);
  secondSuppBackRing->DefineSection(1, zloc, rmin, rmax);

  zloc = kSuppBackRingR10*(1 - TMath::Cos(kSuppBackRing2ndAng1*2/3));
  rmax = secondSuppBackRing->GetRmax(0) - kSuppBackRingR10*TMath::Sin(kSuppBackRing2ndAng1*2/3);
  rmin = secondSuppBackRing->GetRmin(0) - zloc/TMath::Tan(kSuppBackRing2ndAng1);
  secondSuppBackRing->DefineSection(2, zloc, rmin, rmax);

  zloc = kSuppBackRingR10*(1 - TMath::Cos(kSuppBackRing2ndAng1));
  rmax = secondSuppBackRing->GetRmax(0) - kSuppBackRingR10*TMath::Sin(kSuppBackRing2ndAng1);
  rmin = secondSuppBackRing->GetRmin(0) - zloc/TMath::Tan(kSuppBackRing2ndAng1);
  secondSuppBackRing->DefineSection(3, zloc, rmin, rmax);

  slp1 = TMath::Tan(kSuppBackRing2ndAng2);
  slp2 = TMath::Tan(TMath::Pi()/2 + kSuppBackRing2ndAng1);
  IntersectLines(-slp1,kSuppBackRingThikAll,deltaR/2,
		  slp2,kSuppBackRingThikAll,deltaR,
		  xm, ym);

  zloc = xm - kSuppBackRingThick1;
  rmin = secondSuppBackRing->GetRmin(0) - zloc/TMath::Tan(kSuppBackRing2ndAng1);
  rmax = rmin + deltaR;
  secondSuppBackRing->DefineSection(4, zloc, rmin, rmax);

  zloc = kSuppBackRingThikAll - kSuppBackRingThick1;
  rmin = secondSuppBackRing->GetRmin(0) - zloc/TMath::Tan(kSuppBackRing2ndAng1);
  rmax = rmin + deltaR/2;
  secondSuppBackRing->DefineSection(5, zloc, rmin, rmax);

  // The supporting rod: a Tube
  TGeoTube *suppRod = new TGeoTube(0, kBackRodDiameter/2,
				   (kBackRodLength - kBackRodThickLen)/2);

  // The Back Ring (0872/G/C/01): another complex Xtru
  TGeoXtru *suppRearRing = new TGeoXtru(2);

  // First the external profile...
  npoints = 0;

  xprof[npoints] = kSuppRearRingTopWide;
  yprof[npoints] = kSuppRearRingHeight;
  npoints++;

  phi = kSuppRearRing1stAng;
  slp1 = TMath::Tan(TMath::Pi() - phi);
  phi += kSuppRearRingStepAng;
  slp2 = TMath::Tan(TMath::Pi() - phi);
  xm = kSuppRearRingRext2*TMath::Sin(phi);
  ym = kSuppRearRingRext2*TMath::Cos(phi);
  IntersectLines(slp1, kSuppRearRingTopWide, kSuppRearRingHeight,
		 slp2, xm, ym,
		 xprof[npoints], yprof[npoints]);
  npoints++;

  slp1 = slp2;
  phi += kSuppRearRingStepAng;
  slp2 = TMath::Tan(TMath::Pi() - phi);
  xm = kSuppRearRingRext1*TMath::Sin(phi);
  ym = kSuppRearRingRext1*TMath::Cos(phi);
  IntersectLines(slp1, xprof[npoints-1], yprof[npoints-1],
		 slp2, xm, ym,
		 xprof[npoints], yprof[npoints]);
  npoints++;

  slp1 = slp2;
  phi += kSuppRearRingStepAng;
  slp2 = TMath::Tan(TMath::Pi() - phi);
  xm = kSuppRearRingRext2*TMath::Sin(phi);
  ym = kSuppRearRingRext2*TMath::Cos(phi);
  IntersectLines(slp1, xprof[npoints-1], yprof[npoints-1],
		 slp2, xm, ym,
		 xprof[npoints], yprof[npoints]);
  npoints++;

  slp1 = slp2;
  slp2 = 0;
  xm = kSuppRearRingBase;
  ym = kSuppRearRingBaseHi + kSuppRearRingSideHi;
  IntersectLines(slp1, xprof[npoints-1], yprof[npoints-1],
		 slp2, xm, ym,
		 xprof[npoints], yprof[npoints]);
  npoints++;

  xprof[npoints] = kSuppRearRingBase;
  yprof[npoints] = kSuppRearRingBaseHi + kSuppRearRingSideHi;
  npoints++;
  xprof[npoints] = xprof[npoints - 1];
  yprof[npoints] = kSuppRearRingBaseHi;
  npoints++;
  xprof[npoints] = xprof[npoints - 1] - kSuppRearRingInside;
  yprof[npoints] = yprof[npoints - 1];
  npoints++;
  xprof[npoints] = xprof[npoints - 1];
  yprof[npoints] = yprof[npoints - 1] + kSuppRearRingInsideHi;
  npoints++;
  // ...then the internal arc, approximated with segments,...
  xprof[npoints] = kSuppRearRingRint;
  yprof[npoints] = yprof[npoints - 1];

  alphamin = TMath::ASin(kSuppRearRingBaseHi/kSuppRearRingRint);
  alphamax = TMath::Pi()/2;

  for (Int_t jp = 1; jp < kSuppRearRingNPtsArc; jp++) {
    Double_t alpha = alphamin + jp*(alphamax-alphamin)/kSuppRearRingNPtsArc;
    xprof[npoints+jp] = kSuppRearRingRint*TMath::Cos(alpha);
    yprof[npoints+jp] = kSuppRearRingRint*TMath::Sin(alpha);
  }

  xprof[npoints+kSuppRearRingNPtsArc] = 0;
  yprof[npoints+kSuppRearRingNPtsArc] = kSuppRearRingRint;
  // We did the right side! now reflex on the left side
  Int_t nTotalPoints = npoints+kSuppRearRingNPtsArc;
  for (Int_t jp = 0; jp < nTotalPoints; jp++) {
    xprof[nTotalPoints+1+jp] = -xprof[nTotalPoints-1-jp];
    yprof[nTotalPoints+1+jp] =  yprof[nTotalPoints-1-jp];
  }

  // And now the actual Xtru
  suppRearRing->DefinePolygon(2*nTotalPoints+1, xprof, yprof);
  suppRearRing->DefineSection(0,0);
  suppRearRing->DefineSection(1,kSuppRearRingThick);


  // We have all shapes: now create the real volumes
  TGeoMedium *medAl = mgr->GetMedium("ITS_ANTICORODAL$");

  TGeoVolume *sideAExtSuppRing = new TGeoVolume("ITSsuppSideAExtSuppRing",
						 extSuppRing, medAl);

  sideAExtSuppRing->SetVisibility(kTRUE);
  sideAExtSuppRing->SetLineColor(kMagenta+1);
  sideAExtSuppRing->SetLineWidth(1);
  sideAExtSuppRing->SetFillColor(sideAExtSuppRing->GetLineColor());
  sideAExtSuppRing->SetFillStyle(4000); // 0% transparent

  TGeoVolume *sideAIntSuppRing = new TGeoVolume("ITSsuppSideAIntSuppRing",
						 intSuppRing, medAl);

  sideAIntSuppRing->SetVisibility(kTRUE);
  sideAIntSuppRing->SetLineColor(kMagenta+1);
  sideAIntSuppRing->SetLineWidth(1);
  sideAIntSuppRing->SetFillColor(sideAIntSuppRing->GetLineColor());
  sideAIntSuppRing->SetFillStyle(4000); // 0% transparent

  TGeoVolume *sideASuppCyl = new TGeoVolume("ITSsuppSideASuppCyl",
					    interCylind, medAl);

  sideASuppCyl->SetVisibility(kTRUE);
  sideASuppCyl->SetLineColor(kMagenta+1);
  sideASuppCyl->SetLineWidth(1);
  sideASuppCyl->SetFillColor(sideASuppCyl->GetLineColor());
  sideASuppCyl->SetFillStyle(4000); // 0% transparent

  TGeoVolume *sideASuppSpacer = new TGeoVolume("ITSsuppSideASuppSpacer",
					       suppSpacer, medAl);

  sideASuppSpacer->SetVisibility(kTRUE);
  sideASuppSpacer->SetLineColor(kMagenta+1);
  sideASuppSpacer->SetLineWidth(1);
  sideASuppSpacer->SetFillColor(sideASuppSpacer->GetLineColor());
  sideASuppSpacer->SetFillStyle(4000); // 0% transparent

  TGeoVolume *sideASuppForwRing = new TGeoVolume("ITSsuppSideASuppForwRing",
						 forwardRing, medAl);

  sideASuppForwRing->SetVisibility(kTRUE);
  sideASuppForwRing->SetLineColor(kMagenta+1);
  sideASuppForwRing->SetLineWidth(1);
  sideASuppForwRing->SetFillColor(sideASuppForwRing->GetLineColor());
  sideASuppForwRing->SetFillStyle(4000); // 0% transparent

  TGeoVolume *sideASuppForwCone = new TGeoVolume("ITSsuppSideASuppForwCone",
						 forwardCone, medAl);

  sideASuppForwCone->SetVisibility(kTRUE);
  sideASuppForwCone->SetLineColor(kMagenta+1);
  sideASuppForwCone->SetLineWidth(1);
  sideASuppForwCone->SetFillColor(sideASuppForwCone->GetLineColor());
  sideASuppForwCone->SetFillStyle(4000); // 0% transparent

  TGeoVolume *sideAFirstSuppBackRing = new TGeoVolume("ITSsuppSideAFirstSuppBackRing",
						     firstSuppBackRing, medAl);

  sideAFirstSuppBackRing->SetVisibility(kTRUE);
  sideAFirstSuppBackRing->SetLineColor(kMagenta+1);
  sideAFirstSuppBackRing->SetLineWidth(1);
  sideAFirstSuppBackRing->SetFillColor(sideAFirstSuppBackRing->GetLineColor());
  sideAFirstSuppBackRing->SetFillStyle(4000); // 0% transparent

  TGeoVolume *sideASecondSuppBackRing = new TGeoVolume("ITSsuppSideASecondSuppBackRing",
						       secondSuppBackRing, medAl);

  sideASecondSuppBackRing->SetVisibility(kTRUE);
  sideASecondSuppBackRing->SetLineColor(kMagenta+1);
  sideASecondSuppBackRing->SetLineWidth(1);
  sideASecondSuppBackRing->SetFillColor(sideASecondSuppBackRing->GetLineColor());
  sideASecondSuppBackRing->SetFillStyle(4000); // 0% transparent

  TGeoVolume *sideASuppRod = new TGeoVolume("ITSsuppSideASuppRod",
					    suppRod, medAl);

  sideASuppRod->SetVisibility(kTRUE);
  sideASuppRod->SetLineColor(kMagenta+1);
  sideASuppRod->SetLineWidth(1);
  sideASuppRod->SetFillColor(sideASuppRod->GetLineColor());
  sideASuppRod->SetFillStyle(4000); // 0% transparent

  TGeoVolume *sideASuppRearRing = new TGeoVolume("ITSsuppSideASuppRearRing",
						 suppRearRing, medAl);

  sideASuppRearRing->SetVisibility(kTRUE);
  sideASuppRearRing->SetLineColor(kMagenta+1);
  sideASuppRearRing->SetLineWidth(1);
  sideASuppRearRing->SetFillColor(sideASuppRearRing->GetLineColor());
  sideASuppRearRing->SetFillStyle(4000); // 0% transparent


  // Now build up the support structure
  zloc = kSuppRingZTrans;
  trayASuppStruct->AddNode(sideAExtSuppRing, 1,
			   new TGeoTranslation(0, 0, zloc) );
  trayASuppStruct->AddNode(sideAExtSuppRing, 2,
			   new TGeoCombiTrans( 0, 0, zloc,
					       new TGeoRotation("",180,0,0)));

  zloc += kExtSuppRingThick;
  trayASuppStruct->AddNode(sideAIntSuppRing, 1,
			   new TGeoTranslation(0, 0, zloc) );
  trayASuppStruct->AddNode(sideAIntSuppRing, 2,
			   new TGeoCombiTrans( 0, 0, zloc,
					       new TGeoRotation("",180,0,0)));

  xloc = kExtSuppRingBase - kIntSuppRingInward;
  yloc = kSuppRingYTrans;
  zloc += (kIntSuppRingThick + kSuppCylHeight/2);
  trayASuppStruct->AddNode(sideASuppCyl, 1,
			   new TGeoTranslation(0, 0, zloc) );
  trayASuppStruct->AddNode(sideASuppCyl, 2,
			   new TGeoCombiTrans( 0, 0, zloc,
					       new TGeoRotation("",180,0,0)));
  trayASuppStruct->AddNode(sideASuppSpacer, 1,
			   new TGeoCombiTrans( xloc, yloc, zloc,
			   new TGeoRotation("",90+kSuppSpacerAngle,0,0)));
  trayASuppStruct->AddNode(sideASuppSpacer, 2,
			   new TGeoCombiTrans(-xloc, yloc, zloc,
			   new TGeoRotation("",0,180,kSuppSpacerAngle-90)));
  trayASuppStruct->AddNode(sideASuppSpacer, 3,
			   new TGeoCombiTrans( xloc,-yloc, zloc,
			   new TGeoRotation("",180,180,kSuppSpacerAngle-90)));
  trayASuppStruct->AddNode(sideASuppSpacer, 4,
			   new TGeoCombiTrans(-xloc,-yloc, zloc,
			   new TGeoRotation("",270+kSuppSpacerAngle,0,0)));


  zloc += kSuppCylHeight/2;
  trayASuppStruct->AddNode(sideAIntSuppRing, 3,
			   new TGeoTranslation(0, 0, zloc) );
  trayASuppStruct->AddNode(sideAIntSuppRing, 4,
			   new TGeoCombiTrans( 0, 0, zloc,
					       new TGeoRotation("",180,0,0)));

  zloc += kIntSuppRingThick;
  trayASuppStruct->AddNode(sideAExtSuppRing, 3,
			   new TGeoTranslation(0, 0, zloc) );
  trayASuppStruct->AddNode(sideAExtSuppRing, 4,
			   new TGeoCombiTrans( 0, 0, zloc,
					       new TGeoRotation("",180,0,0)));

  zloc += kExtSuppRingThick;
  trayASuppStruct->AddNode(sideASuppForwRing, 1,
			   new TGeoTranslation(0, 0, zloc) );
  trayASuppStruct->AddNode(sideASuppForwRing, 2,
			   new TGeoCombiTrans( 0, 0, zloc,
					       new TGeoRotation("",180,0,0)));

  zloc += kSuppForwRingThikAll;
  trayASuppStruct->AddNode(sideASuppForwCone, 1,
			   new TGeoTranslation(0, 0, zloc) );
  trayASuppStruct->AddNode(sideASuppForwCone, 2,
			   new TGeoCombiTrans( 0, 0, zloc,
					       new TGeoRotation("",180,0,0)));

  zloc += (kSuppForwConeLen1+kSuppForwConeLen2);
  trayASuppStruct->AddNode(sideAFirstSuppBackRing, 1,
			   new TGeoTranslation(0, 0, zloc) );
  trayASuppStruct->AddNode(sideAFirstSuppBackRing, 2,
			   new TGeoCombiTrans( 0, 0, zloc,
					       new TGeoRotation("",180,0,0)));

  zloc += kSuppBackRingThick1;
  trayASuppStruct->AddNode(sideASecondSuppBackRing, 1,
			   new TGeoTranslation(0, 0, zloc) );
  trayASuppStruct->AddNode(sideASecondSuppBackRing, 2,
			   new TGeoCombiTrans( 0, 0, zloc,
					       new TGeoRotation("",180,0,0)));

  xloc = kSuppRearRingXRodHole;
  yloc = kSuppRearRingBaseHi + kSuppRearRingYRodHole;
  zloc = kRearSuppZTransGlob - kBackRodZTrans + suppRod->GetDz();
  trayASuppStruct->AddNode(sideASuppRod, 1,
			   new TGeoTranslation( xloc, yloc, zloc) );
  trayASuppStruct->AddNode(sideASuppRod, 2,
			   new TGeoTranslation(-xloc, yloc, zloc) );
  trayASuppStruct->AddNode(sideASuppRod, 3,
			   new TGeoTranslation( xloc,-yloc, zloc) );
  trayASuppStruct->AddNode(sideASuppRod, 4,
			   new TGeoTranslation(-xloc,-yloc, zloc) );

  zloc += suppRod->GetDz();
  trayASuppStruct->AddNode(sideASuppRearRing, 1,
			   new TGeoTranslation( 0, 0, zloc) );
  trayASuppStruct->AddNode(sideASuppRearRing, 2,
			   new TGeoCombiTrans( 0, 0, zloc,
					       new TGeoRotation("",180,0,0)));


  // Finally put everything in the mother volume
  moth->AddNode(trayASuppStruct,1,0);

  return;
}

//______________________________________________________________________
void AliITSv11GeometrySupport::ServicesCableSupportSPD(TGeoVolume *moth,
						       TGeoManager *mgr){
//
// Creates the all SPD cable trays which are outside the ITS support cones
// but still inside the TPC
// In order to avoid a huge monolithic routine, this method actually
// calls inner methods to create and assemble the various (macro)pieces
//
// Input:
//         moth : the TGeoVolume owing the volume structure
//         mgr  : the GeoManager (default gGeoManager)
// Output:
//
// Created:         ???       Bjorn S. Nilsen
// Updated:      15 Nov 2009  Mario Sitta
//
// Technical data are taken from AutoCAD drawings and other (oral)
// information given by F.Tosello
//

  SPDCableTraysSideA(moth, mgr);
  SPDCableTraysSideC(moth, mgr);

}

//______________________________________________________________________
void AliITSv11GeometrySupport::ServicesCableSupportSDD(TGeoVolume *moth,
						       TGeoManager *mgr){
//
// Creates the all SDD cable trays which are outside the ITS support cones
// but still inside the TPC
// In order to avoid a huge monolithic routine, this method actually
// calls inner methods to create and assemble the various (macro)pieces
//
// Input:
//         moth : the TGeoVolume owing the volume structure
//         mgr  : the GeoManager (default gGeoManager)
// Output:
//
// Created:      14 Dec 2009  Mario Sitta
//

  SDDCableTraysSideA(moth, mgr);
  SDDCableTraysSideC(moth, mgr);

  return;
}

//______________________________________________________________________
void AliITSv11GeometrySupport::ServicesCableSupportSSD(TGeoVolume *moth,
						       TGeoManager *mgr){
//
// Creates the SSD cable trays which are outside the ITS support cones
// but still inside the TPC
// In order to avoid a huge monolithic routine, this method actually
// calls inner methods to create and assemble the various (macro)pieces
//
// Input:
//         moth : the TGeoVolume owing the volume structure
//         mgr  : the GeoManager (default gGeoManager)
// Output:
//
// Created:      15 Nov 2009  Mario Sitta
//

  SSDCableTraysSideA(moth, mgr);
  SSDCableTraysSideC(moth, mgr);

  return;
}

//______________________________________________________________________
void AliITSv11GeometrySupport::SPDCableTraysSideA(TGeoVolume *moth,
					    const TGeoManager *mgr){
//
// Creates the SPD cable trays which are outside the ITS support cones
// but still inside the TPC on Side A
// (part of this code is taken or anyway inspired to ServicesCableSupport
// method of AliITSv11GeometrySupport.cxx,v 1.9 2007/06/06)
//
// Input:
//         moth : the TGeoVolume owing the volume structure
//         mgr  : the GeoManager (default gGeoManager)
// Output:
//
// Created:      15 Feb 2010  Mario Sitta
// Updated:      10 Jun 2010  Mario Sitta  Freon inside cooling pipes
// Updated:      08 Sep 2010  Mario Sitta
// Updated:      14 Sep 2010  Mario Sitta  Cables prolonged till cone
//
// Technical data are taken from AutoCAD drawings, L.Simonetti technical
// drawings and other (oral) information given by F.Tosello and D.Elia
// (small differences with blueprints - e.g. -0.07mm in R1Trans and
// R2Trans - fix small overlaps; they are then compensated in positioning
// the Rear Tray to avoid its own overlaps with the rear supporting ring)
// Optical fibers and voltage cables are approximated with mean materials
// and square cross sections, but preserving the total material budget.
//

  // Overall position and rotation of the A-Side Cable Trays
  // (parts of 0872/G/D)
  const Double_t kTrayAR1Trans           =  396.93 *fgkmm;
  const Double_t kTrayAR2Trans           =  413.93 *fgkmm;
  const Double_t kTrayAZTrans            = 1011.00 *fgkmm;
  const Double_t kTrayAZRot              = (180-169.5);// Degrees
  const Double_t kTrayAFirstRotAng       =   22.00;    // Degrees
  const Double_t kTrayASecondRotAng      =   15.00;    // Degrees

  const Double_t kForwardTrayWide        =   94.00 *fgkmm;//!!!TO BE CHECKED!!!
  const Double_t kForwardTrayFirstHigh   =   83.00 *fgkmm;//!!!TO BE CHECKED!!!
  const Double_t kForwardTraySecondHigh  =   52.70 *fgkmm;//!!!TO BE CHECKED!!!
  const Double_t kForwardTrayTotalLen    =  853.00 *fgkmm;
  const Double_t kForwardTrayFirstLen    =  435.00 *fgkmm;
  const Double_t kForwardTrayWingWide    =   16.00 *fgkmm;//!!!TO BE CHECKED!!!
  const Double_t kForwardTrayInterSpace  =   18.00 *fgkmm;//!!!TO BE CHECKED!!!
  const Double_t kForwardTrayThick       =    2.00 *fgkmm;

  const Int_t    kForwardSideNpoints     =    6;

  const Double_t kExternalTrayLen        = 1200.00 *fgkmm;
  const Double_t kExternalTrayWide       = kForwardTrayWide;
  const Double_t kExternalTrayHigh       = kForwardTraySecondHigh;
  const Double_t kExternalTrayThick      = kForwardTrayThick;

  const Double_t kCoolingTubeRmin        =    2.00 *fgkmm;
  const Double_t kCoolingTubeRmax        =    3.00 *fgkmm;

  const Double_t kOpticalFibersSect      =    8.696*fgkmm;//!!!ESTIMATED!!!
  const Double_t kLowVoltageCableSectCu  =    7.675*fgkmm;// Computed
  const Double_t kLowVoltageCableHighPUR =    1.000*fgkmm;// Computed
  const Double_t kHiVoltageCableSectCu   =    1.535*fgkmm;// Computed
  const Double_t kHiVoltageCableHighPUR  =    0.500*fgkmm;// Computed
  const Double_t kCoaxCableSectCu        =    6.024*fgkmm;// Computed
  const Double_t kCoaxCableHighMeg       =    5.695*fgkmm;// Computed

  const Double_t kTrayCCablesRot         =   75.000*fgkDegree;// Computed
  const Double_t kTrayCCablesZLenOut     =  227.000*fgkmm;// Computed


  // Local variables
  Double_t xprof[kForwardSideNpoints], yprof[kForwardSideNpoints];
  Double_t xloc, yloc, zloc, alpharot;


  // The two tray components as assemblies
  TGeoVolumeAssembly *cableTrayAForw =
    new TGeoVolumeAssembly("ITSsupportSPDTrayAForwRear");
  TGeoVolumeAssembly *cableTrayAExt =
    new TGeoVolumeAssembly("ITSsupportSPDTrayAExt");
  

  // First create all needed shapes

  // The lower face of the forward tray: a BBox
  TGeoBBox *forwTrayLowerFace = new TGeoBBox(kForwardTrayWide/2,
					     kForwardTrayThick/2,
					     kForwardTrayTotalLen/2);

  // The side face of the forward tray: a Xtru
  TGeoXtru *forwTraySideFace = new TGeoXtru(2);
  forwTraySideFace->SetName("ITSsuppSPDForwTraySide");

  xprof[0] = 0;
  yprof[0] = kForwardTrayThick;
  xprof[1] = kForwardTrayTotalLen;
  yprof[1] = yprof[0];
  xprof[2] = xprof[1];
  yprof[2] = kForwardTraySecondHigh - kForwardTrayThick;
  xprof[3] = kForwardTrayFirstLen;
  yprof[3] = yprof[2];
  xprof[4] = xprof[3];
  yprof[4] = kForwardTrayFirstHigh - kForwardTrayThick;
  xprof[5] = xprof[0];
  yprof[5] = yprof[4];

  forwTraySideFace->DefinePolygon(6, xprof, yprof);
  forwTraySideFace->DefineSection(0, 0);
  forwTraySideFace->DefineSection(1, kForwardTrayThick);

  // The covers of the forward tray: two BBox's
  TGeoBBox *forwTrayShortCover = new TGeoBBox(kForwardTrayWide/2,
					      kForwardTrayThick/2,
					      kForwardTrayFirstLen/2);

  TGeoBBox *forwTrayLongCover = new TGeoBBox(kForwardTrayWide/2,
					     kForwardTrayThick/2,
			     (kForwardTrayTotalLen - kForwardTrayFirstLen)/2);

  // Each small wing of the forward tray: a BBox
  TGeoBBox *forwTrayWing = new TGeoBBox(kForwardTrayWingWide/2,
			     (kForwardTrayFirstHigh-kForwardTraySecondHigh)/2,
					kForwardTrayThick/2);

  // The internal plane of the forward tray: a BBox
  TGeoBBox *forwTrayPlane = new TGeoBBox(kForwardTrayWide/2-kForwardTrayThick,
					 kForwardTrayThick/2,
					 kForwardTrayTotalLen/2);

  // The internal wall of the forward tray: a BBox
  TGeoBBox *forwTrayWall = new TGeoBBox(kForwardTrayThick/2,
				 (kForwardTrayInterSpace-kForwardTrayThick)/2,
					kForwardTrayTotalLen/2);

  // Each horizontal face of the external tray: a BBox
  TGeoBBox *extTrayHorFace = new TGeoBBox(kExternalTrayWide/2-kExternalTrayThick,
					  kExternalTrayThick/2,
					  kExternalTrayLen/2);

  // Each vertical face of the external tray: a BBox
  TGeoBBox *extTrayVerFace = new TGeoBBox(kExternalTrayThick/2,
					  kExternalTrayHigh/2,
					  kExternalTrayLen/2);

  // The internal wall of the external tray: a BBox
  TGeoBBox *extTrayWall = new TGeoBBox(kExternalTrayThick/2,
				 (kForwardTrayInterSpace-kExternalTrayThick)/2,
				       kExternalTrayLen/2);

  // The cooling tube inside the forward tray: a Tube
  Double_t zelong = (kForwardTraySecondHigh - 2*kForwardTrayThick
		- 2*forwTrayWall->GetDY() - kCoolingTubeRmax)*SinD(kTrayAZRot);
  Double_t zlen = (zelong + kForwardTrayTotalLen)/2;
  TGeoTube *coolTubeForw = new TGeoTube(0, kCoolingTubeRmax, zlen);

  // The freon inside the forward tray tubes: a Tube
  TGeoTube *freonTubeForw = new TGeoTube(0, kCoolingTubeRmin, zlen);

  // The cooling tube inside the external tray: a Ctub
  TGeoCtub *coolTubeExt = new TGeoCtub(0, kCoolingTubeRmax,
				       kExternalTrayLen/2, 0, 360,
				       0, SinD(kTrayAZRot),-CosD(kTrayAZRot),
				       0,                0,               1);

  // The freon inside the forward tray tubes: a Tube
  TGeoCtub *freonTubeExt = new TGeoCtub(0, kCoolingTubeRmin,
					kExternalTrayLen/2, 0, 360,
					0, SinD(kTrayAZRot),-CosD(kTrayAZRot),
					0,                0,               1);

  // The optical fibers inside the forward tray: a Xtru
  TGeoXtru *optFibsForw = new TGeoXtru(2);

  xprof[0] = -kTrayCCablesZLenOut;
  yprof[0] = xprof[0]/TanD(kTrayCCablesRot);
  xprof[1] = 0;
  yprof[1] = 0;
  xprof[2] = kForwardTrayTotalLen;
  yprof[2] = yprof[1];
  xprof[3] = xprof[2];
  yprof[3] = yprof[2] + kOpticalFibersSect;
  xprof[4] = xprof[1];
  yprof[4] = yprof[3];
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + kOpticalFibersSect;

  optFibsForw->DefinePolygon(6, xprof, yprof);
  optFibsForw->DefineSection(0,-kOpticalFibersSect/2);
  optFibsForw->DefineSection(1, kOpticalFibersSect/2);

  // The optical fibers inside the external tray: a Xtru
  TGeoXtru *optFibsExt = new TGeoXtru(2);
  optFibsExt->SetName("ITSsuppSPDExtTrayOptFibs");

  yprof[0] = -kExternalTrayHigh + 2*kExternalTrayThick
	   + 2*forwTrayWall->GetDY();
  xprof[0] = yprof[0]*TanD(kTrayAZRot);
  xprof[1] = kExternalTrayLen;
  yprof[1] = yprof[0];
  xprof[2] = xprof[1];
  yprof[2] = yprof[1] + kOpticalFibersSect;
  yprof[3] = yprof[2];
  xprof[3] = yprof[2]*TanD(kTrayAZRot);

  optFibsExt->DefinePolygon(4, xprof, yprof);
  optFibsExt->DefineSection(0, 0);
  optFibsExt->DefineSection(1, kOpticalFibersSect);

  // The Low Voltage cables inside the forward tray: two Xtru
  TGeoXtru *lowCablesForwCu = new TGeoXtru(2);

  xprof[0] = -kTrayCCablesZLenOut;
  yprof[0] = xprof[0]/TanD(kTrayCCablesRot);
  xprof[1] = 0;
  yprof[1] = 0;
  xprof[2] = kForwardTrayTotalLen;
  yprof[2] = yprof[1];
  xprof[3] = xprof[2];
  yprof[3] = yprof[2] + kLowVoltageCableSectCu/2;
  xprof[4] = xprof[1];
  yprof[4] = yprof[3];
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + kLowVoltageCableSectCu/2;

  lowCablesForwCu->DefinePolygon(6, xprof, yprof);
  lowCablesForwCu->DefineSection(0,-kLowVoltageCableSectCu);
  lowCablesForwCu->DefineSection(1, kLowVoltageCableSectCu);

  TGeoXtru *lowCablesForwPUR = new TGeoXtru(2);

  xprof[0] = lowCablesForwCu->GetX(5);
  yprof[0] = lowCablesForwCu->GetY(5);
  xprof[1] = lowCablesForwCu->GetX(4);
  yprof[1] = lowCablesForwCu->GetY(4);
  xprof[2] = lowCablesForwCu->GetX(3);
  yprof[2] = lowCablesForwCu->GetY(3);
  xprof[3] = xprof[2];
  yprof[3] = yprof[2] + kLowVoltageCableHighPUR/2;
  xprof[4] = xprof[1];
  yprof[4] = yprof[3];
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + kLowVoltageCableHighPUR/2;

  lowCablesForwPUR->DefinePolygon(6, xprof, yprof);
  lowCablesForwPUR->DefineSection(0,-kLowVoltageCableSectCu);
  lowCablesForwPUR->DefineSection(1, kLowVoltageCableSectCu);

  // The Low Voltage inside the external tray: two Xtru
  TGeoXtru *lowCablesExtCu = new TGeoXtru(2);
  lowCablesExtCu->SetName("ITSsuppSPDExtTrayLowVoltageCu");

  yprof[0] = -kExternalTrayHigh + 2*kExternalTrayThick
	   + 2*forwTrayWall->GetDY();
  xprof[0] = yprof[0]*TanD(kTrayAZRot);
  xprof[1] = kExternalTrayLen;
  yprof[1] = yprof[0];
  xprof[2] = xprof[1];
  yprof[2] = yprof[1] + kLowVoltageCableSectCu/2;
  yprof[3] = yprof[2];
  xprof[3] = yprof[2]*TanD(kTrayAZRot);

  lowCablesExtCu->DefinePolygon(4, xprof, yprof);
  lowCablesExtCu->DefineSection(0, 0);
  lowCablesExtCu->DefineSection(1, kLowVoltageCableSectCu*2);

  TGeoXtru *lowCablesExtPUR = new TGeoXtru(2);
  lowCablesExtPUR->SetName("ITSsuppSPDExtTrayLowVoltagePUR");

  xprof[0] = lowCablesExtCu->GetX(3);
  yprof[0] = lowCablesExtCu->GetY(3);
  xprof[1] = lowCablesExtCu->GetX(2);
  yprof[1] = lowCablesExtCu->GetY(2);
  xprof[2] = xprof[1];
  yprof[2] = yprof[1] + kLowVoltageCableHighPUR/2;
  yprof[3] = yprof[2];
  xprof[3] = yprof[2]*TanD(kTrayAZRot);

  lowCablesExtPUR->DefinePolygon(4, xprof, yprof);
  lowCablesExtPUR->DefineSection(0, 0);
  lowCablesExtPUR->DefineSection(1, kLowVoltageCableSectCu*2);

  // The High Voltage cables inside the forward tray: two Xtru
  TGeoXtru *hiCablesForwCu = new TGeoXtru(2);

  xprof[0] = -kTrayCCablesZLenOut;
  yprof[0] = xprof[0]/TanD(kTrayCCablesRot);
  xprof[1] = 0;
  yprof[1] = 0;
  xprof[2] = kForwardTrayTotalLen;
  yprof[2] = yprof[1];
  xprof[3] = xprof[2];
  yprof[3] = yprof[2] + kHiVoltageCableSectCu/2;
  xprof[4] = xprof[1];
  yprof[4] = yprof[3];
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + kHiVoltageCableSectCu/2;

  hiCablesForwCu->DefinePolygon(6, xprof, yprof);
  hiCablesForwCu->DefineSection(0,-kHiVoltageCableSectCu);
  hiCablesForwCu->DefineSection(1, kHiVoltageCableSectCu);

  TGeoXtru *hiCablesForwPUR = new TGeoXtru(2);

  xprof[0] = hiCablesForwCu->GetX(5);
  yprof[0] = hiCablesForwCu->GetY(5);
  xprof[1] = hiCablesForwCu->GetX(4);
  yprof[1] = hiCablesForwCu->GetY(4);
  xprof[2] = hiCablesForwCu->GetX(3);
  yprof[2] = hiCablesForwCu->GetY(3);
  xprof[3] = xprof[2];
  yprof[3] = yprof[2] + kHiVoltageCableHighPUR/2;
  xprof[4] = xprof[1];
  yprof[4] = yprof[3];
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + kHiVoltageCableHighPUR/2;

  hiCablesForwPUR->DefinePolygon(6, xprof, yprof);
  hiCablesForwPUR->DefineSection(0,-kHiVoltageCableSectCu);
  hiCablesForwPUR->DefineSection(1, kHiVoltageCableSectCu);

  // The High Voltage inside the external tray: two Xtru
  TGeoXtru *hiCablesExtCu = new TGeoXtru(2);
  hiCablesExtCu->SetName("ITSsuppSPDExtTrayHiVoltageCu");

  yprof[0] = -kExternalTrayHigh + 2*kExternalTrayThick
	   + 2*forwTrayWall->GetDY();
  xprof[0] = yprof[0]*TanD(kTrayAZRot);
  xprof[1] = kExternalTrayLen;
  yprof[1] = yprof[0];
  xprof[2] = xprof[1];
  yprof[2] = yprof[1] + kHiVoltageCableSectCu/2;
  yprof[3] = yprof[2];
  xprof[3] = yprof[2]*TanD(kTrayAZRot);

  hiCablesExtCu->DefinePolygon(4, xprof, yprof);
  hiCablesExtCu->DefineSection(0, 0);
  hiCablesExtCu->DefineSection(1, kHiVoltageCableSectCu*2);

  TGeoXtru *hiCablesExtPUR = new TGeoXtru(2);
  hiCablesExtPUR->SetName("ITSsuppSPDExtTrayHiVoltagePUR");

  xprof[0] = hiCablesExtCu->GetX(3);
  yprof[0] = hiCablesExtCu->GetY(3);
  xprof[1] = hiCablesExtCu->GetX(2);
  yprof[1] = hiCablesExtCu->GetY(2);
  xprof[2] = xprof[1];
  yprof[2] = yprof[1] + kHiVoltageCableHighPUR/2;
  yprof[3] = yprof[2];
  xprof[3] = yprof[2]*TanD(kTrayAZRot);

  hiCablesExtPUR->DefinePolygon(4, xprof, yprof);
  hiCablesExtPUR->DefineSection(0, 0);
  hiCablesExtPUR->DefineSection(1, kHiVoltageCableSectCu*2);

  // The Coaxial cables inside the forward tray: two Xtru
  TGeoXtru *coaxCablesForwCu = new TGeoXtru(2);
  coaxCablesForwCu->SetName("ITSsuppSPDForwTrayCoaxCu");

  xprof[0] = -kTrayCCablesZLenOut;
  yprof[0] = xprof[0]/TanD(kTrayCCablesRot);
  xprof[1] = 0;
  yprof[1] = 0;
  xprof[2] = kForwardTrayTotalLen;
  yprof[2] = yprof[1];
  xprof[3] = xprof[2];
  yprof[3] = yprof[2] + kCoaxCableSectCu/2;
  xprof[4] = xprof[1];
  yprof[4] = yprof[3];
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + kCoaxCableSectCu/2;

  coaxCablesForwCu->DefinePolygon(6, xprof, yprof);
  coaxCablesForwCu->DefineSection(0,-kCoaxCableSectCu);
  coaxCablesForwCu->DefineSection(1, kCoaxCableSectCu);

  TGeoXtru *coaxCablesForwMeg = new TGeoXtru(2);
  coaxCablesForwMeg->SetName("ITSsuppSPDForwTrayCoaxMeg");

  xprof[0] = coaxCablesForwCu->GetX(5);
  yprof[0] = coaxCablesForwCu->GetY(5);
  xprof[1] = coaxCablesForwCu->GetX(4);
  yprof[1] = coaxCablesForwCu->GetY(4);
  xprof[2] = coaxCablesForwCu->GetX(3);
  yprof[2] = coaxCablesForwCu->GetY(3);
  xprof[3] = xprof[2];
  yprof[3] = yprof[2] + kCoaxCableHighMeg/2;
  xprof[4] = xprof[1];
  yprof[4] = yprof[3];
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + kCoaxCableHighMeg/2;

  coaxCablesForwMeg->DefinePolygon(6, xprof, yprof);
  coaxCablesForwMeg->DefineSection(0,-kCoaxCableSectCu);
  coaxCablesForwMeg->DefineSection(1, kCoaxCableSectCu);

  // The Coaxial inside the external tray: two Xtru
  TGeoXtru *coaxCablesExtCu = new TGeoXtru(2);
  coaxCablesExtCu->SetName("ITSsuppSPDExtTrayCoaxCu");

  yprof[0] = -kExternalTrayHigh + 2*kExternalTrayThick
	   + 2*forwTrayWall->GetDY();
  xprof[0] = yprof[0]*TanD(kTrayAZRot);
  xprof[1] = kExternalTrayLen;
  yprof[1] = yprof[0];
  xprof[2] = xprof[1];
  yprof[2] = yprof[1] + kCoaxCableSectCu/2;
  yprof[3] = yprof[2];
  xprof[3] = yprof[2]*TanD(kTrayAZRot);

  coaxCablesExtCu->DefinePolygon(4, xprof, yprof);
  coaxCablesExtCu->DefineSection(0, 0);
  coaxCablesExtCu->DefineSection(1, kCoaxCableSectCu*2);

  TGeoXtru *coaxCablesExtMeg = new TGeoXtru(2);
  coaxCablesExtMeg->SetName("ITSsuppSPDExtTrayCoaxMeg");

  xprof[0] = coaxCablesExtCu->GetX(3);
  yprof[0] = coaxCablesExtCu->GetY(3);
  xprof[1] = coaxCablesExtCu->GetX(2);
  yprof[1] = coaxCablesExtCu->GetY(2);
  xprof[2] = xprof[1];
  yprof[2] = yprof[1] + kCoaxCableHighMeg/2;
  yprof[3] = yprof[2];
  xprof[3] = yprof[2]*TanD(kTrayAZRot);

  coaxCablesExtMeg->DefinePolygon(4, xprof, yprof);
  coaxCablesExtMeg->DefineSection(0, 0);
  coaxCablesExtMeg->DefineSection(1, kCoaxCableSectCu*2);


  // We have all shapes: now create the real volumes
  TGeoMedium *medAl    = mgr->GetMedium("ITS_ALUMINUM$");
  TGeoMedium *medIn    = mgr->GetMedium("ITS_INOX$");
  TGeoMedium *medFreon = mgr->GetMedium("ITS_GASEOUS FREON$");
  TGeoMedium *medFibs  = mgr->GetMedium("ITS_SDD OPTICFIB$");//!TO BE CHECKED!
  TGeoMedium *medCu    = mgr->GetMedium("ITS_COPPER$");
  TGeoMedium *medPUR   = mgr->GetMedium("ITS_POLYURETHANE$");
  TGeoMedium *medMeg   = mgr->GetMedium("ITS_MEGOLON$");

  TGeoVolume *forwTrayABase = new TGeoVolume("ITSsuppSPDSideAForwTrayABase",
					    forwTrayLowerFace, medAl);

  forwTrayABase->SetVisibility(kTRUE);
  forwTrayABase->SetLineColor(6); // Purple
  forwTrayABase->SetLineWidth(1);
  forwTrayABase->SetFillColor(forwTrayABase->GetLineColor());
  forwTrayABase->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwTrayASide = new TGeoVolume("ITSsuppSPDSideAForwTrayASide",
					    forwTraySideFace, medAl);

  forwTrayASide->SetVisibility(kTRUE);
  forwTrayASide->SetLineColor(6); // Purple
  forwTrayASide->SetLineWidth(1);
  forwTrayASide->SetFillColor(forwTrayASide->GetLineColor());
  forwTrayASide->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwTrayACoverShort = new TGeoVolume("ITSsuppSPDSideAForwTrayASC",
						  forwTrayShortCover, medAl);

  forwTrayACoverShort->SetVisibility(kTRUE);
  forwTrayACoverShort->SetLineColor(6); // Purple
  forwTrayACoverShort->SetLineWidth(1);
  forwTrayACoverShort->SetFillColor(forwTrayACoverShort->GetLineColor());
  forwTrayACoverShort->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwTrayACoverLong = new TGeoVolume("ITSsuppSPDSideAForwTrayALC",
						 forwTrayLongCover, medAl);

  forwTrayACoverLong->SetVisibility(kTRUE);
  forwTrayACoverLong->SetLineColor(6); // Purple
  forwTrayACoverLong->SetLineWidth(1);
  forwTrayACoverLong->SetFillColor(forwTrayACoverLong->GetLineColor());
  forwTrayACoverLong->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwTrayAWing = new TGeoVolume("ITSsuppSPDSideAForwTrayAWing",
					     forwTrayWing, medAl);

  forwTrayAWing->SetVisibility(kTRUE);
  forwTrayAWing->SetLineColor(6); // Purple
  forwTrayAWing->SetLineWidth(1);
  forwTrayAWing->SetFillColor(forwTrayAWing->GetLineColor());
  forwTrayAWing->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwTrayAPlane = new TGeoVolume("ITSsuppSPDSideAForwTrayAPlane",
					      forwTrayPlane, medAl);

  forwTrayAPlane->SetVisibility(kTRUE);
  forwTrayAPlane->SetLineColor(6); // Purple
  forwTrayAPlane->SetLineWidth(1);
  forwTrayAPlane->SetFillColor(forwTrayAPlane->GetLineColor());
  forwTrayAPlane->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwTrayAWall = new TGeoVolume("ITSsuppSPDSideAForwTrayAWall",
					     forwTrayWall, medAl);

  forwTrayAWall->SetVisibility(kTRUE);
  forwTrayAWall->SetLineColor(6); // Purple
  forwTrayAWall->SetLineWidth(1);
  forwTrayAWall->SetFillColor(forwTrayAWall->GetLineColor());
  forwTrayAWall->SetFillStyle(4000); // 0% transparent

  TGeoVolume *extTrayAHorFace = new TGeoVolume("ITSsuppSPDSideAExtTrayHorFace",
					       extTrayHorFace, medAl);

  extTrayAHorFace->SetVisibility(kTRUE);
  extTrayAHorFace->SetLineColor(6); // Purple
  extTrayAHorFace->SetLineWidth(1);
  extTrayAHorFace->SetFillColor(extTrayAHorFace->GetLineColor());
  extTrayAHorFace->SetFillStyle(4000); // 0% transparent

  TGeoVolume *extTrayAVerFace = new TGeoVolume("ITSsuppSPDSideAExtTrayVerFace",
					       extTrayVerFace, medAl);

  extTrayAVerFace->SetVisibility(kTRUE);
  extTrayAVerFace->SetLineColor(6); // Purple
  extTrayAVerFace->SetLineWidth(1);
  extTrayAVerFace->SetFillColor(extTrayAVerFace->GetLineColor());
  extTrayAVerFace->SetFillStyle(4000); // 0% transparent

  TGeoVolume *extTrayAWall = new TGeoVolume("ITSsuppSPDSideAExtTrayWall",
					    extTrayWall, medAl);

  extTrayAWall->SetVisibility(kTRUE);
  extTrayAWall->SetLineColor(6); // Purple
  extTrayAWall->SetLineWidth(1);
  extTrayAWall->SetFillColor(extTrayAWall->GetLineColor());
  extTrayAWall->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwCoolTube = new TGeoVolume("ITSsuppSPDSideAForwTrayCoolTube",
					    coolTubeForw, medIn);

  forwCoolTube->SetVisibility(kTRUE);
  forwCoolTube->SetLineColor(kGray); // as in GeometrySPD
  forwCoolTube->SetLineWidth(1);
  forwCoolTube->SetFillColor(forwCoolTube->GetLineColor());
  forwCoolTube->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwCoolFreon = new TGeoVolume("ITSsuppSPDSideAForwTrayFreon",
					     freonTubeForw, medFreon);

  forwCoolFreon->SetVisibility(kTRUE);
  forwCoolFreon->SetLineColor(kBlue); // Blue
  forwCoolFreon->SetLineWidth(1);
  forwCoolFreon->SetFillColor(forwCoolFreon->GetLineColor());
  forwCoolFreon->SetFillStyle(4000); // 0% transparent

  TGeoVolume *extCoolTube = new TGeoVolume("ITSsuppSPDSideAExtTrayCoolTube",
					   coolTubeExt, medIn);

  extCoolTube->SetVisibility(kTRUE);
  extCoolTube->SetLineColor(kGray); // as in GeometrySPD
  extCoolTube->SetLineWidth(1);
  extCoolTube->SetFillColor(extCoolTube->GetLineColor());
  extCoolTube->SetFillStyle(4000); // 0% transparent

  TGeoVolume *extCoolFreon = new TGeoVolume("ITSsuppSPDSideAExtTrayFreon",
					    freonTubeExt, medFreon);

  extCoolFreon->SetVisibility(kTRUE);
  extCoolFreon->SetLineColor(kBlue); // Blue
  extCoolFreon->SetLineWidth(1);
  extCoolFreon->SetFillColor(extCoolFreon->GetLineColor());
  extCoolFreon->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwOptFibs = new TGeoVolume("ITSsuppSPDSideAForwTrayOptFibs",
					   optFibsForw, medFibs);

  forwOptFibs->SetVisibility(kTRUE);
  forwOptFibs->SetLineColor(kOrange); // Orange
  forwOptFibs->SetLineWidth(1);
  forwOptFibs->SetFillColor(forwOptFibs->GetLineColor());
  forwOptFibs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *extOptFibs = new TGeoVolume("ITSsuppSPDSideAExtTrayOptFibs",
					  optFibsExt, medFibs);

  extOptFibs->SetVisibility(kTRUE);
  extOptFibs->SetLineColor(kOrange); // Orange
  extOptFibs->SetLineWidth(1);
  extOptFibs->SetFillColor(extOptFibs->GetLineColor());
  extOptFibs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwLowCabsCu = new TGeoVolume("ITSsuppSPDSideAForwLowCabsCu",
					     lowCablesForwCu, medCu);

  forwLowCabsCu->SetVisibility(kTRUE);
  forwLowCabsCu->SetLineColor(kRed); // Red
  forwLowCabsCu->SetLineWidth(1);
  forwLowCabsCu->SetFillColor(forwLowCabsCu->GetLineColor());
  forwLowCabsCu->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwLowCabsPUR = new TGeoVolume("ITSsuppSPDSideAForwLowCabsPUR",
					      lowCablesForwPUR, medPUR);

  forwLowCabsPUR->SetVisibility(kTRUE);
  forwLowCabsPUR->SetLineColor(kBlack); // Black
  forwLowCabsPUR->SetLineWidth(1);
  forwLowCabsPUR->SetFillColor(forwLowCabsPUR->GetLineColor());
  forwLowCabsPUR->SetFillStyle(4000); // 0% transparent

  TGeoVolume *extLowCabsCu = new TGeoVolume("ITSsuppSPDSideAExtLowCabsCu",
					    lowCablesExtCu, medCu);

  extLowCabsCu->SetVisibility(kTRUE);
  extLowCabsCu->SetLineColor(kRed); // Red
  extLowCabsCu->SetLineWidth(1);
  extLowCabsCu->SetFillColor(extLowCabsCu->GetLineColor());
  extLowCabsCu->SetFillStyle(4000); // 0% transparent

  TGeoVolume *extLowCabsPUR = new TGeoVolume("ITSsuppSPDSideAExtLowCabsPUR",
					     lowCablesExtPUR, medPUR);

  extLowCabsPUR->SetVisibility(kTRUE);
  extLowCabsPUR->SetLineColor(kBlack); // Black
  extLowCabsPUR->SetLineWidth(1);
  extLowCabsPUR->SetFillColor(extLowCabsPUR->GetLineColor());
  extLowCabsPUR->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwHiCabsCu = new TGeoVolume("ITSsuppSPDSideAForwTrayHiCabsCu",
					    hiCablesForwCu, medCu);

  forwHiCabsCu->SetVisibility(kTRUE);
  forwHiCabsCu->SetLineColor(kRed); // Red
  forwHiCabsCu->SetLineWidth(1);
  forwHiCabsCu->SetFillColor(forwHiCabsCu->GetLineColor());
  forwHiCabsCu->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwHiCabsPUR = new TGeoVolume("ITSsuppSPDSideAForwTrayHiCabsPUR",
					     hiCablesForwPUR, medPUR);

  forwHiCabsPUR->SetVisibility(kTRUE);
  forwHiCabsPUR->SetLineColor(kBlack); // Black
  forwHiCabsPUR->SetLineWidth(1);
  forwHiCabsPUR->SetFillColor(forwHiCabsPUR->GetLineColor());
  forwHiCabsPUR->SetFillStyle(4000); // 0% transparent

  TGeoVolume *extHiCabsCu = new TGeoVolume("ITSsuppSPDSideAExtTrayHiCabsCu",
					   hiCablesExtCu, medCu);

  extHiCabsCu->SetVisibility(kTRUE);
  extHiCabsCu->SetLineColor(kRed); // Red
  extHiCabsCu->SetLineWidth(1);
  extHiCabsCu->SetFillColor(extHiCabsCu->GetLineColor());
  extHiCabsCu->SetFillStyle(4000); // 0% transparent

  TGeoVolume *extHiCabsPUR = new TGeoVolume("ITSsuppSPDSideAExtTrayHiCabsPUR",
					    hiCablesExtPUR, medPUR);

  extHiCabsPUR->SetVisibility(kTRUE);
  extHiCabsPUR->SetLineColor(kBlack); // Black
  extHiCabsPUR->SetLineWidth(1);
  extHiCabsPUR->SetFillColor(extHiCabsPUR->GetLineColor());
  extHiCabsPUR->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwCoaxCu = new TGeoVolume("ITSsuppSPDSideAForwTrayCoaxCu",
					  coaxCablesForwCu, medCu);

  forwCoaxCu->SetVisibility(kTRUE);
  forwCoaxCu->SetLineColor(kRed); // Red
  forwCoaxCu->SetLineWidth(1);
  forwCoaxCu->SetFillColor(forwCoaxCu->GetLineColor());
  forwCoaxCu->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwCoaxMeg = new TGeoVolume("ITSsuppSPDSideAForwTrayCoaxMeg",
					   coaxCablesForwMeg, medMeg);

  forwCoaxMeg->SetVisibility(kTRUE);
  forwCoaxMeg->SetLineColor(kBlack); // Black
  forwCoaxMeg->SetLineWidth(1);
  forwCoaxMeg->SetFillColor(forwCoaxMeg->GetLineColor());
  forwCoaxMeg->SetFillStyle(4000); // 0% transparent

  TGeoVolume *extCoaxCu = new TGeoVolume("ITSsuppSPDSideAExtTrayCoaxCu",
					 coaxCablesExtCu, medCu);

  extCoaxCu->SetVisibility(kTRUE);
  extCoaxCu->SetLineColor(kRed); // Red
  extCoaxCu->SetLineWidth(1);
  extCoaxCu->SetFillColor(extCoaxCu->GetLineColor());
  extCoaxCu->SetFillStyle(4000); // 0% transparent

  TGeoVolume *extCoaxMeg = new TGeoVolume("ITSsuppSPDSideAExtTrayCoaxMeg",
					  coaxCablesExtMeg, medMeg);

  extCoaxMeg->SetVisibility(kTRUE);
  extCoaxMeg->SetLineColor(kBlack); // Black
  extCoaxMeg->SetLineWidth(1);
  extCoaxMeg->SetFillColor(extCoaxMeg->GetLineColor());
  extCoaxMeg->SetFillStyle(4000); // 0% transparent


  // Now build up the trays
  yloc = forwTrayLowerFace->GetDY();
  zloc = forwTrayLowerFace->GetDZ();
  cableTrayAForw->AddNode(forwTrayABase, 1,
		      new TGeoTranslation(0, yloc, zloc));

  xloc = kForwardTrayWide/2;
  cableTrayAForw->AddNode(forwTrayASide, 1,
		      new TGeoCombiTrans( xloc, 0, 0,
					 new TGeoRotation("",90,-90,-90)));
  cableTrayAForw->AddNode(forwTrayASide, 2,
		      new TGeoCombiTrans(-xloc+kForwardTrayThick, 0, 0,
					 new TGeoRotation("",90,-90,-90)));

  yloc = kForwardTrayFirstHigh - forwTrayShortCover->GetDY();
  zloc = forwTrayShortCover->GetDZ();
  cableTrayAForw->AddNode(forwTrayACoverShort, 1,
		      new TGeoTranslation(0, yloc, zloc));

  yloc = kForwardTraySecondHigh - forwTrayLongCover->GetDY();
  zloc = kForwardTrayFirstLen + forwTrayLongCover->GetDZ();
  cableTrayAForw->AddNode(forwTrayACoverLong, 1,
		      new TGeoTranslation(0, yloc, zloc));

  xloc = kForwardTrayWide/2 - kForwardTrayThick - forwTrayWing->GetDX();
  yloc = kForwardTrayFirstHigh - kForwardTrayThick - forwTrayWing->GetDY();
  zloc = kForwardTrayFirstLen - forwTrayWing->GetDZ();
  cableTrayAForw->AddNode(forwTrayAWing, 1,
		      new TGeoTranslation( xloc, yloc, zloc));
  cableTrayAForw->AddNode(forwTrayAWing, 2,
		      new TGeoTranslation(-xloc, yloc, zloc));

  yloc = kForwardTrayThick + kForwardTrayInterSpace - forwTrayPlane->GetDY();
  zloc = forwTrayPlane->GetDZ();
  cableTrayAForw->AddNode(forwTrayAPlane, 1,
		      new TGeoTranslation(0, yloc, zloc));

  yloc = kForwardTrayThick + forwTrayWall->GetDY();
  zloc = forwTrayWall->GetDZ();
  cableTrayAForw->AddNode(forwTrayAWall, 1,
		      new TGeoTranslation(0, yloc, zloc));

  forwCoolTube->AddNode(forwCoolFreon, 1, 0);

  yloc = 2*kForwardTrayThick + 2*forwTrayWall->GetDY()
       + coolTubeForw->GetRmax();
  zloc = coolTubeForw->GetDz();
  cableTrayAForw->AddNode(forwCoolTube, 1,
		      new TGeoTranslation(0, yloc, zloc));

  xloc = optFibsForw->GetZ(1) + coolTubeForw->GetRmax();
  yloc = 2*kForwardTrayThick + 2*forwTrayWall->GetDY();
  cableTrayAForw->AddNode(forwOptFibs, 1,
		      new TGeoCombiTrans( xloc, yloc, 0,
					 new TGeoRotation("",-90.,90.,90.)));

  xloc = 2*optFibsForw->GetZ(1) + lowCablesForwCu->GetZ(1) +
	 coolTubeForw->GetRmax();
  yloc = 2*kForwardTrayThick + 2*forwTrayWall->GetDY();
  cableTrayAForw->AddNode(forwLowCabsCu, 1,
		      new TGeoCombiTrans( xloc, yloc, 0,
					 new TGeoRotation("",-90.,90.,90.)));
  cableTrayAForw->AddNode(forwLowCabsPUR, 1,
		      new TGeoCombiTrans( xloc, yloc, 0,
					 new TGeoRotation("",-90.,90.,90.)));

  xloc = 2*optFibsForw->GetZ(1) + 2*lowCablesForwCu->GetZ(1) +
	 hiCablesForwCu->GetZ(1) + coolTubeForw->GetRmax();
  yloc = 2*kForwardTrayThick + 2*forwTrayWall->GetDY();
  cableTrayAForw->AddNode(forwHiCabsCu, 1,
		      new TGeoCombiTrans( xloc, yloc, 0,
					 new TGeoRotation("",-90.,90.,90.)));
  cableTrayAForw->AddNode(forwHiCabsPUR, 1,
		      new TGeoCombiTrans( xloc, yloc, 0,
					 new TGeoRotation("",-90.,90.,90.)));

  xloc = coaxCablesForwCu->GetZ(1) + coolTubeForw->GetRmax();
  yloc = 2*kForwardTrayThick + 2*forwTrayWall->GetDY();
  cableTrayAForw->AddNode(forwCoaxCu, 1,
		      new TGeoCombiTrans(-xloc, yloc, 0,
					 new TGeoRotation("",-90.,90.,90.)));
  cableTrayAForw->AddNode(forwCoaxMeg, 1,
		      new TGeoCombiTrans(-xloc, yloc, 0,
					 new TGeoRotation("",-90.,90.,90.)));

  // To simplify following placement in MARS, origin is on top
  yloc = -kExternalTrayHigh + kExternalTrayThick/2;
  zloc = kExternalTrayLen/2;
  cableTrayAExt->AddNode(extTrayAHorFace, 1,
		      new TGeoTranslation( 0, yloc, zloc));

  xloc = kExternalTrayWide/2 - kExternalTrayThick/2;
  yloc = -kExternalTrayHigh/2;
  cableTrayAExt->AddNode(extTrayAVerFace, 1,
		      new TGeoTranslation( xloc, yloc, zloc));
  cableTrayAExt->AddNode(extTrayAVerFace, 2,
		      new TGeoTranslation(-xloc, yloc, zloc));

  yloc = -kExternalTrayThick/2;
  cableTrayAExt->AddNode(extTrayAHorFace, 2,
		      new TGeoTranslation( 0, yloc, zloc));

  yloc = -kExternalTrayHigh
       + kExternalTrayThick + kForwardTrayInterSpace - kExternalTrayThick/2;
  cableTrayAExt->AddNode(extTrayAHorFace, 3,
		      new TGeoTranslation( 0, yloc, zloc));

  yloc = -kExternalTrayHigh + kExternalTrayThick + extTrayWall->GetDY();
  cableTrayAExt->AddNode(extTrayAWall, 1,
		      new TGeoTranslation( 0, yloc, zloc));

  extCoolTube->AddNode(extCoolFreon, 1, 0);

  yloc = -kExternalTrayHigh + 2*kExternalTrayThick + 2*extTrayWall->GetDY()
       + coolTubeExt->GetRmax();
  zloc = coolTubeExt->GetDz();
  cableTrayAExt->AddNode(extCoolTube, 1,
		      new TGeoTranslation(0, yloc, zloc));

  xloc = optFibsExt->GetZ(1) + coolTubeExt->GetRmax();
  cableTrayAExt->AddNode(extOptFibs, 1,
		      new TGeoCombiTrans( xloc, 0, 0,
					 new TGeoRotation("",90,-90,-90)));

  xloc = coolTubeExt->GetRmax();
  cableTrayAExt->AddNode(extLowCabsCu, 1,
		      new TGeoCombiTrans(-xloc, 0, 0,
					 new TGeoRotation("",90,-90,-90)));
  cableTrayAExt->AddNode(extLowCabsPUR, 1,
		      new TGeoCombiTrans(-xloc, 0, 0,
					 new TGeoRotation("",90,-90,-90)));

  xloc = lowCablesExtCu->GetZ(1) + coolTubeExt->GetRmax();
  cableTrayAExt->AddNode(extHiCabsCu, 1,
		      new TGeoCombiTrans(-xloc, 0, 0,
					 new TGeoRotation("",90,-90,-90)));
  cableTrayAExt->AddNode(extHiCabsPUR, 1,
		      new TGeoCombiTrans(-xloc, 0, 0,
					 new TGeoRotation("",90,-90,-90)));

  xloc = coaxCablesExtCu->GetZ(1) + optFibsExt->GetZ(1) +
	 coolTubeExt->GetRmax();
  cableTrayAExt->AddNode(extCoaxCu, 1,
		      new TGeoCombiTrans( xloc, 0, 0,
					 new TGeoRotation("",90,-90,-90)));
  cableTrayAExt->AddNode(extCoaxMeg, 1,
		      new TGeoCombiTrans( xloc, 0, 0,
					 new TGeoRotation("",90,-90,-90)));


  // Finally put everything in the mother volume
  Double_t rExtTray = kTrayAR2Trans + kExternalTrayHigh;

  moth->AddNode(cableTrayAForw,1,
		new TGeoTranslation( 0, kTrayAR1Trans, kTrayAZTrans));
  moth->AddNode(cableTrayAForw,2,
		new TGeoCombiTrans(  0,-kTrayAR1Trans, kTrayAZTrans,
				    new TGeoRotation("",180, 0, 0)));

  yloc = kTrayAR1Trans + kExternalTrayHigh;
  zloc = kTrayAZTrans + kForwardTrayTotalLen;
  moth->AddNode(cableTrayAExt,1,
		new TGeoCombiTrans( 0, yloc, zloc,
				    new TGeoRotation("",  0,-kTrayAZRot, 0)));
  moth->AddNode(cableTrayAExt,2,
		new TGeoCombiTrans( 0,-yloc, zloc,
				    new TGeoRotation("",180,-kTrayAZRot, 0)));

  alpharot = kTrayAFirstRotAng + kTrayASecondRotAng;
  xloc = kTrayAR2Trans*SinD(alpharot);
  yloc = kTrayAR2Trans*CosD(alpharot);
  moth->AddNode(cableTrayAForw,3,
			    new TGeoCombiTrans( xloc, yloc, kTrayAZTrans,
			    new TGeoRotation("",-alpharot,0,0)   )   );
  xloc = rExtTray*SinD(alpharot);
  yloc = rExtTray*CosD(alpharot);
  moth->AddNode(cableTrayAExt,3,
			    new TGeoCombiTrans( xloc, yloc, zloc,
			    new TGeoRotation("",-alpharot,-kTrayAZRot,0)  )  );

  alpharot += 180;
  xloc = kTrayAR2Trans*SinD(alpharot);
  yloc = kTrayAR2Trans*CosD(alpharot);
  moth->AddNode(cableTrayAForw,4,
			    new TGeoCombiTrans( xloc, yloc, kTrayAZTrans,
			    new TGeoRotation("",-alpharot,0,0)   )   );
  xloc = rExtTray*SinD(alpharot);
  yloc = rExtTray*CosD(alpharot);
  moth->AddNode(cableTrayAExt,4,
			    new TGeoCombiTrans( xloc, yloc, zloc,
			    new TGeoRotation("",-alpharot,-kTrayAZRot,0)  )  );

  alpharot = - kTrayAFirstRotAng - kTrayASecondRotAng;
  xloc = kTrayAR2Trans*SinD(alpharot);
  yloc = kTrayAR2Trans*CosD(alpharot);
  moth->AddNode(cableTrayAForw,5,
			    new TGeoCombiTrans( xloc, yloc, kTrayAZTrans,
			    new TGeoRotation("",-alpharot,0,0)   )   );
  xloc = rExtTray*SinD(alpharot);
  yloc = rExtTray*CosD(alpharot);
  moth->AddNode(cableTrayAExt,5,
			    new TGeoCombiTrans( xloc, yloc, zloc,
			    new TGeoRotation("",-alpharot,-kTrayAZRot,0)  )  );

  alpharot += 180;
  xloc = kTrayAR2Trans*SinD(alpharot);
  yloc = kTrayAR2Trans*CosD(alpharot);
  moth->AddNode(cableTrayAForw,6,
			    new TGeoCombiTrans( xloc, yloc, kTrayAZTrans,
			    new TGeoRotation("",-alpharot,0,0)   )   );
  xloc = rExtTray*SinD(alpharot);
  yloc = rExtTray*CosD(alpharot);
  moth->AddNode(cableTrayAExt,6,
			    new TGeoCombiTrans( xloc, yloc, zloc,
			    new TGeoRotation("",-alpharot,-kTrayAZRot,0)  )  );

  alpharot = kTrayAFirstRotAng + 3*kTrayASecondRotAng;
  xloc = kTrayAR2Trans*SinD(alpharot);
  yloc = kTrayAR2Trans*CosD(alpharot);
  moth->AddNode(cableTrayAForw,7,
			    new TGeoCombiTrans( xloc, yloc, kTrayAZTrans,
			    new TGeoRotation("",-alpharot,0,0)   )   );
  xloc = rExtTray*SinD(alpharot);
  yloc = rExtTray*CosD(alpharot);
  moth->AddNode(cableTrayAExt,7,
			    new TGeoCombiTrans( xloc, yloc, zloc,
			    new TGeoRotation("",-alpharot,-kTrayAZRot,0)  )  );

  alpharot += 180;
  xloc = kTrayAR2Trans*SinD(alpharot);
  yloc = kTrayAR2Trans*CosD(alpharot);
  moth->AddNode(cableTrayAForw,8,
			    new TGeoCombiTrans( xloc, yloc, kTrayAZTrans,
			    new TGeoRotation("",-alpharot,0,0)   )   );
  xloc = rExtTray*SinD(alpharot);
  yloc = rExtTray*CosD(alpharot);
  moth->AddNode(cableTrayAExt,8,
			    new TGeoCombiTrans( xloc, yloc, zloc,
			    new TGeoRotation("",-alpharot,-kTrayAZRot,0)  )  );

  alpharot = - kTrayAFirstRotAng - 3*kTrayASecondRotAng;
  xloc = kTrayAR2Trans*SinD(alpharot);
  yloc = kTrayAR2Trans*CosD(alpharot);
  moth->AddNode(cableTrayAForw,9,
			    new TGeoCombiTrans( xloc, yloc, kTrayAZTrans,
			    new TGeoRotation("",-alpharot,0,0)   )   );
  xloc = rExtTray*SinD(alpharot);
  yloc = rExtTray*CosD(alpharot);
  moth->AddNode(cableTrayAExt,9,
			    new TGeoCombiTrans( xloc, yloc, zloc,
			    new TGeoRotation("",-alpharot,-kTrayAZRot,0)  )  );

  alpharot += 180;
  xloc = kTrayAR2Trans*SinD(alpharot);
  yloc = kTrayAR2Trans*CosD(alpharot);
  moth->AddNode(cableTrayAForw,10,
			    new TGeoCombiTrans( xloc, yloc, kTrayAZTrans,
			    new TGeoRotation("",-alpharot,0,0)   )   );
  xloc = rExtTray*SinD(alpharot);
  yloc = rExtTray*CosD(alpharot);
  moth->AddNode(cableTrayAExt,10,
			    new TGeoCombiTrans( xloc, yloc, zloc,
			    new TGeoRotation("",-alpharot,-kTrayAZRot,0)  )  );


  return;
}

//______________________________________________________________________
void AliITSv11GeometrySupport::SPDCableTraysSideC(TGeoVolume *moth,
					    const TGeoManager *mgr){
//
// Creates the SPD cable trays which are outside the ITS support cones
// but still inside the TPC on Side C
// (part of this code is taken or anyway inspired to ServicesCableSupport
// method of AliITSv11GeometrySupport.cxx,v 1.9 2007/06/06)
//
// Input:
//         moth : the TGeoVolume owing the volume structure
//         mgr  : the GeoManager (default gGeoManager)
// Output:
//
// Return:
//
// Created:         ???       Bjorn S. Nilsen
// Updated:      22 Apr 2010  Mario Sitta
// Updated:      10 Jun 2010  Mario Sitta  Freon inside cooling pipes
// Updated:      08 Sep 2010  Mario Sitta
// Updated:      14 Sep 2010  Mario Sitta  Cables prolonged till cone
//
// Technical data are taken from AutoCAD drawings and other (oral)
// information given by D.Elia
// Optical fibers and voltage cables are approximated with mean materials
// and square cross sections, but preserving the total material budget.
//

  // Dimensions and positions of the C-Side Cable Tray elements
  const Int_t    kNumTraysSideC       =   10;

  const Double_t kTrayCCablesOutRot   =   75.000 *fgkDegree;// Computed
  const Double_t kTrayCCablesZLenOut  =  245.000 *fgkmm;// Computed

  const Double_t kTrayCHalfWide       =    6.350 *fgkcm;
  const Double_t kTrayCLength1        =  172.800 *fgkcm;
  const Double_t kTrayCLength2        =  189.300 *fgkcm;
  const Double_t kTrayCFirstLen       =  435.000 *fgkmm;
  const Double_t kTrayCFirstHigh      =   83.000 *fgkmm;//!!!TO BE CHECKED!!!
  const Double_t kTrayCSecondHigh     =   52.700 *fgkmm;//!!!TO BE CHECKED!!!
  const Double_t kTrayCThick          =    0.200 *fgkcm;
  const Double_t kTrayCInterSpace     =   18.000 *fgkmm;//!!!TO BE CHECKED!!!
  const Double_t kTrayCFoldAngle      =    5.000 *fgkDegree;

  const Double_t kCoolingTubeRmin     =    2.000 *fgkmm;
  const Double_t kCoolingTubeRmax     =    3.000 *fgkmm;
  const Double_t kOpticalFibersSect   =    8.696 *fgkmm;//!!!ESTIMATED!!!
  const Double_t kLowVoltCableSectCu  =    7.675 *fgkmm;// Computed
  const Double_t kLowVoltCableHighPUR =    1.000 *fgkmm;// Computed
  const Double_t kHiVoltCableSectCu   =    1.535 *fgkmm;// Computed
  const Double_t kHiVoltCableHighPUR  =    0.500 *fgkmm;// Computed
  const Double_t kCoaxCableSectCu     =    6.024 *fgkmm;// Computed
  const Double_t kCoaxCableHighMeg    =    5.695 *fgkmm;// Computed

  // Overall position and rotation of the C-Side Cable Trays
  const Double_t kTraySideCRPos       =   45.300 *fgkcm;
  const Double_t kTraySideCZPos       = -102.400 *fgkcm;
  const Double_t kTraySideCAlphaRot[kNumTraysSideC/2]  =
    {    0.0,      41.0,     -41.0,      76.0,      -76.0};
  // From position of the other trays


  // Local variables
  Double_t xprof[8], yprof[8];
  Double_t xloc, yloc, zloc, delta, alpharot;


  // The single C-Side Cable tray as an assembly
  TGeoVolumeAssembly *cableTrayC = new TGeoVolumeAssembly("ITSsupportSPDTrayC");

  // First create all needed shapes

  // The Cable Tray lower face: a Xtru
  TGeoXtru *sideCHorFace = new TGeoXtru(2);

  xprof[0] = 0.;
  yprof[0] = 0.;
  xprof[1] = kTrayCLength1;
  yprof[1] = 0.;
  xprof[2] = xprof[1] + kTrayCLength2*CosD(kTrayCFoldAngle);
  yprof[2] = yprof[1] + kTrayCLength2*SinD(kTrayCFoldAngle);
  xprof[3] = xprof[2] - kTrayCThick*SinD(kTrayCFoldAngle);
  yprof[3] = yprof[2] + kTrayCThick*CosD(kTrayCFoldAngle);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      kTrayCThick , xprof[4], yprof[4]);
  xprof[5] = 0.;
  yprof[5] = kTrayCThick;

  delta = kTrayCHalfWide - kTrayCThick;

  sideCHorFace->DefinePolygon(6, xprof, yprof);
  sideCHorFace->DefineSection(0,-delta);
  sideCHorFace->DefineSection(1, delta);

  // The Cable Tray middle face: a Xtru
  // (somehow duplicate of HorFace, but in this way avoid an overlap with Wall)
  TGeoXtru *sideCMidFace = new TGeoXtru(2);

  xprof[0] = 0.;
  yprof[0] = kTrayCInterSpace + kTrayCThick;
  xprof[1] = kTrayCLength1;
  yprof[1] = yprof[0];
  xprof[2] = xprof[1] + kTrayCLength2*CosD(kTrayCFoldAngle);
  yprof[2] = yprof[1] + kTrayCLength2*SinD(kTrayCFoldAngle);
  xprof[3] = xprof[2] - kTrayCThick*SinD(kTrayCFoldAngle);
  yprof[3] = yprof[2] + kTrayCThick*CosD(kTrayCFoldAngle);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      kTrayCThick , xprof[4], yprof[4]);
  xprof[5] = 0.;
  yprof[5] = yprof[0] + kTrayCThick;

  delta = kTrayCHalfWide - kTrayCThick;

  sideCMidFace->DefinePolygon(6, xprof, yprof);
  sideCMidFace->DefineSection(0,-delta);
  sideCMidFace->DefineSection(1, delta);

  // The Cable Tray lower face: a Xtru
  TGeoXtru *sideCSideFace = new TGeoXtru(2);

  xprof[0] = 0.;
  yprof[0] = 0.;
  xprof[1] = kTrayCLength1;
  yprof[1] = 0.;
  xprof[2] = xprof[1] + kTrayCLength2*CosD(kTrayCFoldAngle);
  yprof[2] = yprof[1] + kTrayCLength2*SinD(kTrayCFoldAngle);
  xprof[3] = xprof[2] - kTrayCSecondHigh*SinD(kTrayCFoldAngle);
  yprof[3] = yprof[2] + kTrayCSecondHigh*CosD(kTrayCFoldAngle);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      kTrayCSecondHigh , xprof[4], yprof[4]);
  xprof[5] = kTrayCFirstLen;
  yprof[5] = kTrayCSecondHigh;
  xprof[6] = xprof[5];
  yprof[6] = kTrayCFirstHigh;
  xprof[7] = xprof[0];
  yprof[7] = yprof[6];

  sideCSideFace->DefinePolygon(8, xprof, yprof);
  sideCSideFace->DefineSection(0, 0);
  sideCSideFace->DefineSection(1, kTrayCThick);

  // The short cover: a BBox
  TGeoBBox *sideCShortCover = new TGeoBBox(kTrayCFirstLen/2,
					   kTrayCThick/2,
					   kTrayCHalfWide-kTrayCThick);

  // The long cover: a Xtru
  TGeoXtru *sideCLongCover = new TGeoXtru(2);

  xprof[5] = sideCSideFace->GetX(5);
  yprof[5] = sideCSideFace->GetY(5);
  xprof[4] = sideCSideFace->GetX(4);
  yprof[4] = sideCSideFace->GetY(4);
  xprof[3] = sideCSideFace->GetX(3);
  yprof[3] = sideCSideFace->GetY(3);
  xprof[2] = xprof[3] + kTrayCThick*SinD(kTrayCFoldAngle);
  yprof[2] = yprof[3] - kTrayCThick*CosD(kTrayCFoldAngle);
  InsidePoint(xprof[5], yprof[5], xprof[4], yprof[4], xprof[3], yprof[3],
	     -kTrayCThick , xprof[1], yprof[1]);
  xprof[0] = xprof[5];
  yprof[0] = yprof[5] - kTrayCThick;

  delta = kTrayCHalfWide - kTrayCThick;

  sideCLongCover->DefinePolygon(6, xprof, yprof);
  sideCLongCover->DefineSection(0,-delta);
  sideCLongCover->DefineSection(1, delta);

  // The internal wall: a Xtru
  TGeoXtru *intWall = new TGeoXtru(2);

  xprof[0] = sideCHorFace->GetX(5);
  yprof[0] = sideCHorFace->GetY(5);
  xprof[1] = sideCHorFace->GetX(4);
  yprof[1] = sideCHorFace->GetY(4);
  xprof[2] = sideCHorFace->GetX(3);
  yprof[2] = sideCHorFace->GetY(3);
  xprof[3] = sideCMidFace->GetX(2);
  yprof[3] = sideCMidFace->GetY(2);
  xprof[4] = sideCMidFace->GetX(1);
  yprof[4] = sideCMidFace->GetY(1);
  xprof[5] = sideCMidFace->GetX(0);
  yprof[5] = sideCMidFace->GetY(0);

  intWall->DefinePolygon(6, xprof, yprof);
  intWall->DefineSection(0,-kTrayCThick/2);
  intWall->DefineSection(1, kTrayCThick/2);

  // The horizontal part of the cooling tube inside the tray: a Tube
  delta = sideCMidFace->GetX(4) - sideCMidFace->GetX(5);
  TGeoTube *horTube = new TGeoTube(0, kCoolingTubeRmax, delta/2);

  // The freon inside the horizontal part of the cooling tube: a Tube
  TGeoTube *horFreon = new TGeoTube(0, kCoolingTubeRmin, delta/2);

  // The inclined part of the cooling tube inside the tray: a Ctub
  Double_t x3, y3, x4, y4;
  x3 = sideCMidFace->GetX(3);
  y3 = sideCMidFace->GetY(3);
  x4 = sideCMidFace->GetX(4);
  y4 = sideCMidFace->GetY(4);
  delta = TMath::Sqrt( (x4 - x3 + kCoolingTubeRmax*SinD(kTrayCFoldAngle))*
		       (x4 - x3 + kCoolingTubeRmax*SinD(kTrayCFoldAngle))    +
       (y4 + kCoolingTubeRmax - y3 - kCoolingTubeRmax*SinD(kTrayCFoldAngle))*
       (y4 + kCoolingTubeRmax - y3 - kCoolingTubeRmax*SinD(kTrayCFoldAngle)) );

  TGeoCtub *incTube = new TGeoCtub(0, kCoolingTubeRmax, delta/2, 0, 360,
			       0, SinD(kTrayCFoldAngle),-CosD(kTrayCFoldAngle),
			       0,                     0,                    1);

  // The freon inside the inclined part of the cooling tube: a Ctub
  TGeoCtub *incFreon = new TGeoCtub(0, kCoolingTubeRmin, delta/2, 0, 360,
			       0, SinD(kTrayCFoldAngle),-CosD(kTrayCFoldAngle),
			       0,                     0,                    1);

  // The part of the cooling tube outside the tray: a Ctub
  TGeoCtub *outTube = new TGeoCtub(0, kCoolingTubeRmax,
			0.5*kTrayCCablesZLenOut/SinD(kTrayCCablesOutRot),
			0, 360,
			0,                        0,                      -1,
			0,-SinD(kTrayCCablesOutRot), CosD(kTrayCCablesOutRot));

  // The freon inside the part of the cooling tube outside the tray: a Ctub
  TGeoCtub *outFreon = new TGeoCtub(0, kCoolingTubeRmin,
			outTube->GetDz(),
			0, 360,
			0,                        0,                      -1,
			0,-SinD(kTrayCCablesOutRot), CosD(kTrayCCablesOutRot));

  // The optical fibers inside the tray: a Xtru
  TGeoXtru *optFibs = new TGeoXtru(2);

  xprof[0] = -kTrayCCablesZLenOut;
  yprof[0] = xprof[0]/TanD(kTrayCCablesOutRot);
  xprof[1] = sideCMidFace->GetX(5);
  yprof[1] = sideCMidFace->GetY(5);
  xprof[2] = sideCMidFace->GetX(4);
  yprof[2] = sideCMidFace->GetY(4);
  xprof[3] = sideCMidFace->GetX(3);
  yprof[3] = sideCMidFace->GetY(3);
  xprof[4] = xprof[3] - kOpticalFibersSect*SinD(kTrayCFoldAngle);
  yprof[4] = yprof[3] + kOpticalFibersSect*CosD(kTrayCFoldAngle);
  InsidePoint(xprof[1], yprof[1], xprof[2], yprof[2], xprof[3], yprof[3],
	      kOpticalFibersSect , xprof[5], yprof[5]);
  xprof[6] = 0.;
  yprof[6] = yprof[1] + kOpticalFibersSect;
  xprof[7] = xprof[0];
  yprof[7] = yprof[0] + kOpticalFibersSect;

  optFibs->DefinePolygon(8, xprof, yprof);
  optFibs->DefineSection(0, 0);
  optFibs->DefineSection(1, kOpticalFibersSect);

  // The low voltage cables inside the tray: two Xtru
  TGeoXtru *lowCablesCu = new TGeoXtru(2);

  xprof[0] = -kTrayCCablesZLenOut;
  yprof[0] = xprof[0]/TanD(kTrayCCablesOutRot);
  xprof[1] = sideCMidFace->GetX(5);
  yprof[1] = sideCMidFace->GetY(5);
  xprof[2] = sideCMidFace->GetX(4);
  yprof[2] = sideCMidFace->GetY(4);
  xprof[3] = sideCMidFace->GetX(3);
  yprof[3] = sideCMidFace->GetY(3);
  xprof[4] = xprof[3] - kLowVoltCableSectCu*SinD(kTrayCFoldAngle);
  yprof[4] = yprof[3] + kLowVoltCableSectCu*CosD(kTrayCFoldAngle);
  InsidePoint(xprof[1], yprof[1], xprof[2], yprof[2], xprof[3], yprof[3],
	      kLowVoltCableSectCu , xprof[5], yprof[5]);
  xprof[6] = 0.;
  yprof[6] = yprof[1] + kLowVoltCableSectCu;
  xprof[7] = xprof[0];
  yprof[7] = yprof[0] + kLowVoltCableSectCu;

  lowCablesCu->DefinePolygon(8, xprof, yprof);
  lowCablesCu->DefineSection(0, 0);
  lowCablesCu->DefineSection(1, kLowVoltCableSectCu);

  TGeoXtru *lowCablesPUR = new TGeoXtru(2);

  xprof[0] = lowCablesCu->GetX(7);
  yprof[0] = lowCablesCu->GetY(7);
  xprof[1] = lowCablesCu->GetX(6);
  yprof[1] = lowCablesCu->GetY(6);
  xprof[2] = lowCablesCu->GetX(5);
  yprof[2] = lowCablesCu->GetY(5);
  xprof[3] = lowCablesCu->GetX(4);
  yprof[3] = lowCablesCu->GetY(4);
  xprof[4] = xprof[3] - kLowVoltCableHighPUR*SinD(kTrayCFoldAngle);
  yprof[4] = yprof[3] + kLowVoltCableHighPUR*CosD(kTrayCFoldAngle);
  InsidePoint(xprof[1], yprof[1], xprof[2], yprof[2], xprof[3], yprof[3],
	      kLowVoltCableHighPUR , xprof[5], yprof[5]);
  xprof[6] = 0.;
  yprof[6] = yprof[1] + kLowVoltCableHighPUR;
  xprof[7] = xprof[0];
  yprof[7] = yprof[0] + kLowVoltCableHighPUR;

  lowCablesPUR->DefinePolygon(8, xprof, yprof);
  lowCablesPUR->DefineSection(0, 0);
  lowCablesPUR->DefineSection(1, kLowVoltCableSectCu);

  // The high voltage cables inside the tray: two Xtru
  TGeoXtru *hiCablesCu = new TGeoXtru(2);

  xprof[0] = -kTrayCCablesZLenOut;
  yprof[0] = xprof[0]/TanD(kTrayCCablesOutRot);
  xprof[1] = sideCMidFace->GetX(5);
  yprof[1] = sideCMidFace->GetY(5);
  xprof[2] = sideCMidFace->GetX(4);
  yprof[2] = sideCMidFace->GetY(4);
  xprof[3] = sideCMidFace->GetX(3);
  yprof[3] = sideCMidFace->GetY(3);
  xprof[4] = xprof[3] - kHiVoltCableSectCu*SinD(kTrayCFoldAngle);
  yprof[4] = yprof[3] + kHiVoltCableSectCu*CosD(kTrayCFoldAngle);
  InsidePoint(xprof[1], yprof[1], xprof[2], yprof[2], xprof[3], yprof[3],
	      kHiVoltCableSectCu , xprof[5], yprof[5]);
  xprof[6] = 0.;
  yprof[6] = yprof[1] + kHiVoltCableSectCu;
  xprof[7] = xprof[0];
  yprof[7] = yprof[0] + kHiVoltCableSectCu;

  hiCablesCu->DefinePolygon(8, xprof, yprof);
  hiCablesCu->DefineSection(0, 0);
  hiCablesCu->DefineSection(1, kHiVoltCableSectCu);

  TGeoXtru *hiCablesPUR = new TGeoXtru(2);

  xprof[0] = hiCablesCu->GetX(7);
  yprof[0] = hiCablesCu->GetY(7);
  xprof[1] = hiCablesCu->GetX(6);
  yprof[1] = hiCablesCu->GetY(6);
  xprof[2] = hiCablesCu->GetX(5);
  yprof[2] = hiCablesCu->GetY(5);
  xprof[3] = hiCablesCu->GetX(4);
  yprof[3] = hiCablesCu->GetY(4);
  xprof[4] = xprof[3] - kHiVoltCableHighPUR*SinD(kTrayCFoldAngle);
  yprof[4] = yprof[3] + kHiVoltCableHighPUR*CosD(kTrayCFoldAngle);
  InsidePoint(xprof[1], yprof[1], xprof[2], yprof[2], xprof[3], yprof[3],
	      kHiVoltCableHighPUR , xprof[5], yprof[5]);
  xprof[6] = 0.;
  yprof[6] = yprof[1] + kHiVoltCableHighPUR;
  xprof[7] = xprof[0];
  yprof[7] = yprof[0] + kHiVoltCableHighPUR;

  hiCablesPUR->DefinePolygon(8, xprof, yprof);
  hiCablesPUR->DefineSection(0, 0);
  hiCablesPUR->DefineSection(1, kHiVoltCableSectCu);

  // The coaxial cables inside the tray: two Xtru
  TGeoXtru *coaxCablesCu = new TGeoXtru(2);

  xprof[0] = -kTrayCCablesZLenOut;
  yprof[0] = xprof[0]/TanD(kTrayCCablesOutRot);
  xprof[1] = sideCMidFace->GetX(5);
  yprof[1] = sideCMidFace->GetY(5);
  xprof[2] = sideCMidFace->GetX(4);
  yprof[2] = sideCMidFace->GetY(4);
  xprof[3] = sideCMidFace->GetX(3);
  yprof[3] = sideCMidFace->GetY(3);
  xprof[4] = xprof[3] - kCoaxCableSectCu*SinD(kTrayCFoldAngle);
  yprof[4] = yprof[3] + kCoaxCableSectCu*CosD(kTrayCFoldAngle);
  InsidePoint(xprof[1], yprof[1], xprof[2], yprof[2], xprof[3], yprof[3],
	      kCoaxCableSectCu , xprof[5], yprof[5]);
  xprof[6] = 0.;
  yprof[6] = yprof[1] + kCoaxCableSectCu;
  xprof[7] = xprof[0];
  yprof[7] = yprof[0] + kCoaxCableSectCu;

  coaxCablesCu->DefinePolygon(8, xprof, yprof);
  coaxCablesCu->DefineSection(0, 0);
  coaxCablesCu->DefineSection(1, kCoaxCableSectCu);

  TGeoXtru *coaxCablesMeg = new TGeoXtru(2);

  xprof[0] = coaxCablesCu->GetX(7);
  yprof[0] = coaxCablesCu->GetY(7);
  xprof[1] = coaxCablesCu->GetX(6);
  yprof[1] = coaxCablesCu->GetY(6);
  xprof[2] = coaxCablesCu->GetX(5);
  yprof[2] = coaxCablesCu->GetY(5);
  xprof[3] = coaxCablesCu->GetX(4);
  yprof[3] = coaxCablesCu->GetY(4);
  xprof[4] = xprof[3] - kCoaxCableHighMeg*SinD(kTrayCFoldAngle);
  yprof[4] = yprof[3] + kCoaxCableHighMeg*CosD(kTrayCFoldAngle);
  InsidePoint(xprof[1], yprof[1], xprof[2], yprof[2], xprof[3], yprof[3],
	      kCoaxCableHighMeg , xprof[5], yprof[5]);
  xprof[6] = 0.;
  yprof[6] = yprof[1] + kCoaxCableHighMeg;
  xprof[7] = xprof[0];
  yprof[7] = yprof[0] + kCoaxCableHighMeg;

  coaxCablesMeg->DefinePolygon(8, xprof, yprof);
  coaxCablesMeg->DefineSection(0, 0);
  coaxCablesMeg->DefineSection(1, kCoaxCableSectCu);


  // We have all shapes: now create the real volumes
  TGeoMedium *medAl   = mgr->GetMedium("ITS_ALUMINUM$");
  TGeoMedium *medIn   = mgr->GetMedium("ITS_INOX$");
  TGeoMedium *medFr   = mgr->GetMedium("ITS_Freon$");
  TGeoMedium *medFibs = mgr->GetMedium("ITS_SDD OPTICFIB$");//!!TO BE CHECKED!!
  TGeoMedium *medCu   = mgr->GetMedium("ITS_COPPER$");
  TGeoMedium *medPUR  = mgr->GetMedium("ITS_POLYURETHANE$");
  TGeoMedium *medMeg  = mgr->GetMedium("ITS_MEGOLON$");

  TGeoVolume *traySideCHorFace  = new TGeoVolume("ITSsuppSPDTraySideCHor",
						 sideCHorFace, medAl);

  traySideCHorFace->SetVisibility(kTRUE);
  traySideCHorFace->SetLineColor(6); // Purple
  traySideCHorFace->SetLineWidth(1);
  traySideCHorFace->SetFillColor(traySideCHorFace->GetLineColor());
  traySideCHorFace->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCMidFace  = new TGeoVolume("ITSsuppSPDTraySideCMid",
						 sideCMidFace, medAl);

  traySideCMidFace->SetVisibility(kTRUE);
  traySideCMidFace->SetLineColor(6); // Purple
  traySideCMidFace->SetLineWidth(1);
  traySideCMidFace->SetFillColor(traySideCMidFace->GetLineColor());
  traySideCMidFace->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCSideFace  = new TGeoVolume("ITSsuppSPDTraySideCSide",
						  sideCSideFace, medAl);

  traySideCSideFace->SetVisibility(kTRUE);
  traySideCSideFace->SetLineColor(6); // Purple
  traySideCSideFace->SetLineWidth(1);
  traySideCSideFace->SetFillColor(traySideCSideFace->GetLineColor());
  traySideCSideFace->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCShortCover  = new TGeoVolume("ITSsuppSPDTraySideCShCov",
						    sideCShortCover, medAl);

  traySideCShortCover->SetVisibility(kTRUE);
  traySideCShortCover->SetLineColor(6); // Purple
  traySideCShortCover->SetLineWidth(1);
  traySideCShortCover->SetFillColor(traySideCShortCover->GetLineColor());
  traySideCShortCover->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCLongCover  = new TGeoVolume("ITSsuppSPDTraySideCLnCov",
						   sideCLongCover, medAl);

  traySideCLongCover->SetVisibility(kTRUE);
  traySideCLongCover->SetLineColor(6); // Purple
  traySideCLongCover->SetLineWidth(1);
  traySideCLongCover->SetFillColor(traySideCLongCover->GetLineColor());
  traySideCLongCover->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCIntWall  = new TGeoVolume("ITSsuppSPDTraySideCWall",
						 intWall, medAl);

  traySideCIntWall->SetVisibility(kTRUE);
  traySideCIntWall->SetLineColor(6); // Purple
  traySideCIntWall->SetLineWidth(1);
  traySideCIntWall->SetFillColor(traySideCIntWall->GetLineColor());
  traySideCIntWall->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCHorTube = new TGeoVolume("ITSsuppSPDTraySideCHorTube",
						horTube, medIn);

  traySideCHorTube->SetVisibility(kTRUE);
  traySideCHorTube->SetLineColor(kGray); // as in GeometrySPD
  traySideCHorTube->SetLineWidth(1);
  traySideCHorTube->SetFillColor(traySideCHorTube->GetLineColor());
  traySideCHorTube->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCHorFreon = new TGeoVolume("ITSsuppSPDTraySideCHorFreon",
						 horFreon, medFr);

  traySideCHorFreon->SetVisibility(kTRUE);
  traySideCHorFreon->SetLineColor(kBlue); // Blue
  traySideCHorFreon->SetLineWidth(1);
  traySideCHorFreon->SetFillColor(traySideCHorFreon->GetLineColor());
  traySideCHorFreon->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCIncTube = new TGeoVolume("ITSsuppSPDTraySideCIncTube",
						incTube, medIn);

  traySideCIncTube->SetVisibility(kTRUE);
  traySideCIncTube->SetLineColor(kGray); // as in GeometrySPD
  traySideCIncTube->SetLineWidth(1);
  traySideCIncTube->SetFillColor(traySideCIncTube->GetLineColor());
  traySideCIncTube->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCIncFreon = new TGeoVolume("ITSsuppSPDTraySideCIncFreon",
						 incFreon, medFr);

  traySideCIncFreon->SetVisibility(kTRUE);
  traySideCIncFreon->SetLineColor(kBlue); // Blue
  traySideCIncFreon->SetLineWidth(1);
  traySideCIncFreon->SetFillColor(traySideCIncFreon->GetLineColor());
  traySideCIncFreon->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCOutTube = new TGeoVolume("ITSsuppSPDTraySideCOutTube",
						outTube, medIn);

  traySideCOutTube->SetVisibility(kTRUE);
  traySideCOutTube->SetLineColor(kGray); // as in GeometrySPD
  traySideCOutTube->SetLineWidth(1);
  traySideCOutTube->SetFillColor(traySideCOutTube->GetLineColor());
  traySideCOutTube->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCOutFreon = new TGeoVolume("ITSsuppSPDTraySideCOutFreon",
						 outFreon, medFr);

  traySideCOutFreon->SetVisibility(kTRUE);
  traySideCOutFreon->SetLineColor(kBlue); // Blue
  traySideCOutFreon->SetLineWidth(1);
  traySideCOutFreon->SetFillColor(traySideCOutFreon->GetLineColor());
  traySideCOutFreon->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCOptFibs = new TGeoVolume("ITSsuppSPDTraySideCOptFibs",
						optFibs, medFibs);

  traySideCOptFibs->SetVisibility(kTRUE);
  traySideCOptFibs->SetLineColor(kOrange); // Orange
  traySideCOptFibs->SetLineWidth(1);
  traySideCOptFibs->SetFillColor(traySideCOptFibs->GetLineColor());
  traySideCOptFibs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCLowCabsCu = new TGeoVolume("ITSsuppSPDTraySideCLVCu",
						  lowCablesCu, medCu);

  traySideCLowCabsCu->SetVisibility(kTRUE);
  traySideCLowCabsCu->SetLineColor(kRed); // Red
  traySideCLowCabsCu->SetLineWidth(1);
  traySideCLowCabsCu->SetFillColor(traySideCLowCabsCu->GetLineColor());
  traySideCLowCabsCu->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCLowCabsPUR = new TGeoVolume("ITSsuppSPDTraySideCLVPUR",
						   lowCablesPUR, medPUR);

  traySideCLowCabsPUR->SetVisibility(kTRUE);
  traySideCLowCabsPUR->SetLineColor(kBlack); // Black
  traySideCLowCabsPUR->SetLineWidth(1);
  traySideCLowCabsPUR->SetFillColor(traySideCLowCabsPUR->GetLineColor());
  traySideCLowCabsPUR->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCHiCabsCu = new TGeoVolume("ITSsuppSPDTraySideCHVCu",
						 hiCablesCu, medCu);

  traySideCHiCabsCu->SetVisibility(kTRUE);
  traySideCHiCabsCu->SetLineColor(kRed); // Red
  traySideCHiCabsCu->SetLineWidth(1);
  traySideCHiCabsCu->SetFillColor(traySideCHiCabsCu->GetLineColor());
  traySideCHiCabsCu->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCHiCabsPUR = new TGeoVolume("ITSsuppSPDTraySideCHVPUR",
						  hiCablesPUR, medPUR);

  traySideCHiCabsPUR->SetVisibility(kTRUE);
  traySideCHiCabsPUR->SetLineColor(kBlack); // Black
  traySideCHiCabsPUR->SetLineWidth(1);
  traySideCHiCabsPUR->SetFillColor(traySideCHiCabsPUR->GetLineColor());
  traySideCHiCabsPUR->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCCoaxCu = new TGeoVolume("ITSsuppSPDTraySideCCoaxCu",
					       coaxCablesCu, medCu);

  traySideCCoaxCu->SetVisibility(kTRUE);
  traySideCCoaxCu->SetLineColor(kRed); // Red
  traySideCCoaxCu->SetLineWidth(1);
  traySideCCoaxCu->SetFillColor(traySideCCoaxCu->GetLineColor());
  traySideCCoaxCu->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCCoaxMeg = new TGeoVolume("ITSsuppSPDTraySideCCoaxMeg",
						coaxCablesMeg, medMeg);

  traySideCCoaxMeg->SetVisibility(kTRUE);
  traySideCCoaxMeg->SetLineColor(kBlack); // Black
  traySideCCoaxMeg->SetLineWidth(1);
  traySideCCoaxMeg->SetFillColor(traySideCCoaxMeg->GetLineColor());
  traySideCCoaxMeg->SetFillStyle(4000); // 0% transparent


  // Now build up the trays
  cableTrayC->AddNode(traySideCHorFace,1,0);

  cableTrayC->AddNode(traySideCMidFace,1,0);

  zloc = kTrayCHalfWide - kTrayCThick;
  cableTrayC->AddNode(traySideCSideFace, 1,
		      new TGeoTranslation( 0, 0, zloc));
  zloc = -kTrayCHalfWide;
  cableTrayC->AddNode(traySideCSideFace, 2,
		      new TGeoTranslation( 0, 0, zloc));

  xloc = sideCShortCover->GetDX();
  yloc = kTrayCFirstHigh - sideCShortCover->GetDY();
  cableTrayC->AddNode(traySideCShortCover, 1,
		      new TGeoTranslation( xloc, yloc, 0));

  cableTrayC->AddNode(traySideCLongCover,1,0);

  cableTrayC->AddNode(traySideCIntWall,1,0);

  traySideCHorTube->AddNode(traySideCHorFreon, 1, 0);
  traySideCIncTube->AddNode(traySideCIncFreon, 1, 0);
  traySideCOutTube->AddNode(traySideCOutFreon, 1, 0);

  xloc = horTube->GetDz();
  yloc = sideCMidFace->GetY(5) + horTube->GetRmax();
  cableTrayC->AddNode(traySideCHorTube, 1,
		      new TGeoCombiTrans( xloc, yloc, 0,
		      new TGeoRotation("",-90.,-90.,90.)));

  xloc = sideCMidFace->GetX(4) + (incTube->GetDz())*CosD(kTrayCFoldAngle);
  yloc = sideCMidFace->GetY(4) +  incTube->GetRmax() +
	    (incTube->GetDz())*SinD(kTrayCFoldAngle)+0.005;//Avoid small ovrlp
  cableTrayC->AddNode(traySideCIncTube, 1,
		      new TGeoCombiTrans( xloc, yloc, 0,
		      new TGeoRotation("",-90.+kTrayCFoldAngle,-90.,90.)));

  xloc = -kTrayCCablesZLenOut/2 - outTube->GetRmax();
  yloc = xloc/TanD(kTrayCCablesOutRot) + sideCMidFace->GetY(4) -
	 2*outTube->GetRmax();
  cableTrayC->AddNode(traySideCOutTube, 1,
		      new TGeoCombiTrans( xloc, yloc, 0,
		      new TGeoRotation("",-70.,-90.,90.)));

  zloc = horTube->GetRmax();
  cableTrayC->AddNode(traySideCOptFibs, 1,
		      new TGeoTranslation( 0, 0, zloc));

  zloc = kLowVoltCableSectCu + horTube->GetRmax();
  cableTrayC->AddNode(traySideCLowCabsCu, 1,
		      new TGeoTranslation( 0, 0,-zloc));
  cableTrayC->AddNode(traySideCLowCabsPUR, 1,
		      new TGeoTranslation( 0, 0,-zloc));

  zloc = kHiVoltCableSectCu + kLowVoltCableSectCu + horTube->GetRmax();
  cableTrayC->AddNode(traySideCHiCabsCu, 1,
		      new TGeoTranslation( 0, 0,-zloc));
  cableTrayC->AddNode(traySideCHiCabsPUR, 1,
		      new TGeoTranslation( 0, 0,-zloc));

  zloc = kOpticalFibersSect + kCoaxCableSectCu + horTube->GetRmax();
  cableTrayC->AddNode(traySideCCoaxCu, 1,
		      new TGeoTranslation( 0, 0, zloc));
  cableTrayC->AddNode(traySideCCoaxMeg, 1,
		      new TGeoTranslation( 0, 0, zloc));


  // Finally put everything in the mother volume
  for (Int_t jt = 0; jt < kNumTraysSideC/2; jt++) {
    alpharot = kTraySideCAlphaRot[jt];

    xloc = kTraySideCRPos*SinD(alpharot);
    yloc = kTraySideCRPos*CosD(alpharot);
    moth->AddNode(cableTrayC,2*jt+1,
		new TGeoCombiTrans(-xloc, yloc, kTraySideCZPos,
		new TGeoRotation("",-90.+alpharot,-90.,90.+kTrayCFoldAngle)));
    alpharot += 180;
    xloc = kTraySideCRPos*SinD(alpharot);
    yloc = kTraySideCRPos*CosD(alpharot);
    moth->AddNode(cableTrayC,2*jt+2,
		new TGeoCombiTrans(-xloc, yloc, kTraySideCZPos,
		new TGeoRotation("",-90.+alpharot,-90.,90.+kTrayCFoldAngle)));
  }


  return;
}

//______________________________________________________________________
void AliITSv11GeometrySupport::SDDCableTraysSideA(TGeoVolume *moth,
						  TGeoManager *mgr){
//
// Creates the SDD cable trays which are outside the ITS support cones
// but still inside the TPC on Side A
// (part of this code is taken or anyway inspired to ServicesCableSupport
// method of AliITSv11GeometrySupport.cxx,v 1.9 2007/06/06)
//
// Input:
//         moth : the TGeoVolume owing the volume structure
//         mgr  : the GeoManager (default gGeoManager)
// Output:
//
// Created:         ???       Bjorn S. Nilsen
// Updated:       5 Jan 2010  Mario Sitta
// Updated:      26 Feb 2010  Mario Sitta
// Updated:      06 Sep 2010  Mario Sitta
//
// Technical data are taken from AutoCAD drawings, L.Simonetti technical
// drawings and other (oral) information given by F.Tosello
//

  // Overall position and rotation of the A-Side Cable Trays
  // (parts of 0872/G/D)
  const Double_t kTrayARTrans            =  408.35 *fgkmm;
  const Double_t kTrayAZTrans            = 1011.00 *fgkmm;
  const Double_t kTrayAZToSupportRing    =  435.00 *fgkmm;
  const Double_t kExternTrayYTrans       =   96.00 *fgkmm; // Computed
  const Double_t kExternTrayZTrans       =  823.00 *fgkmm;
  const Double_t kExternCoverYTrans      =    2.00 *fgkmm;
  const Double_t kTrayAZRot              = (180-169.5);// Degrees
  const Double_t kTrayAFirstRotAng       =   22.00;    // Degrees
  const Double_t kTrayASecondRotAng      =   15.00;    // Degrees

  const Double_t kForwardTrayThick       =    2.00 *fgkmm;
  const Double_t kForwardTrayTailHeight  =  100.00 *fgkmm; // Computed
  const Double_t kForwardTrayTotalHeight =  170.00 *fgkmm; // Computed
  const Double_t kForwardTrayUpperLength =  405.00 *fgkmm; // Computed
  const Double_t kForwardCoverLength     =  380.00 *fgkmm;
  const Double_t kForwardCoverWide       =  133.00 *fgkmm;
  const Double_t kForwardCoverHeight     =   10.00 *fgkmm;
  const Double_t kForwardCoverThick      =    1.00 *fgkmm;

  const Double_t kExternTrayTotalLen     = 1200.00 *fgkmm;
  const Double_t kExternTrayTotalHeight  =   52.00 *fgkmm;
  const Double_t kExternCoverLen         = kExternTrayTotalLen;
  const Double_t kExternCoverThick       =    5.00 *fgkmm;
  const Double_t kExternCoverSideThick   =    3.00 *fgkmm;

  const Int_t    kForwardTrayNpoints     =    8;

  // Dimensions and positions of the Cable Tray elements
  const Double_t kSideACoolManifWide     =    8.23 *fgkcm;
  const Double_t kSideACoolManifHigh     =    8.06 *fgkcm;
  const Double_t kSideACoolManifLen      =    3.90 *fgkcm;
  const Double_t kSideACoolManifPOMFrac  =    0.0054;
  const Double_t kSideACoolManifSteelFrac=    0.8850;
  const Double_t kSideACoolManifWaterFrac=    0.0913;
  const Double_t kSideACoolManifAlFrac   =    0.0183;

  const Double_t kSideACoolTubesWide     =    9.07 *fgkcm;
  const Double_t kSideACoolTubesHigh     =    1.88 *fgkcm;
  const Double_t kSideACoolTubesTrans    =    0.88 *fgkcm;
  const Double_t kSideACoolTubesPURFrac  =    0.5897;
  const Double_t kSideACoolTubesWaterFrac=    0.4101;
  const Double_t kSideACoolTubesAirFrac  =    0.0002;

  const Double_t kSideAOptConnWide       =    0.90    *fgkcm;
  const Double_t kSideAOptConnLen        =    1.37    *fgkcm;
  const Double_t kSideAOptConnPBTFrac    =    0.5010;
  const Double_t kSideAOptConnSteelFrac  =    0.1784;
  const Double_t kSideAOptConnAlFrac     =    0.3206;

  const Double_t kSideAOptFibsWide       =    0.71    *fgkcm;
  const Double_t kSideAOptFibsHigh       =    3.20    *fgkcm;

  const Double_t kSideAInputCablesWide   =   12.50    *fgkcm;
  const Double_t kSideAInputCablesHigh   =    1.24    *fgkcm;
  const Double_t kSideAInputCablesLen    =   25.20    *fgkcm;
  const Double_t kSideAInputCablesYTrans =    1.15    *fgkcm;
  const Double_t kSideAInputCablesCu     =    0.7404;
  const Double_t kSideAInputCablesPlast  =    0.1269;
  const Double_t kSideAInputCablesAl     =    0.0057;
  const Double_t kSideAInputCablesKapton =    0.0172;
  const Double_t kSideAInputCablesPOLYAX =    0.1098;

  const Double_t kSideAOutputCablesWide  =    8.30    *fgkcm;
  const Double_t kSideAOutputCablesHigh  =    1.56    *fgkcm;
  const Double_t kSideAOutputCablesCu    =    0.6783;
  const Double_t kSideAOutputCablesPlast =    0.1605;
  const Double_t kSideAOutputCablesAl    =    0.0078;
  const Double_t kSideAOutputCablesKapton=    0.0232;
  const Double_t kSideAOutputCablesPOLYAX=    0.1302;

  const Double_t kSideAPCBBoardsWide     =   12.50    *fgkcm;
  const Double_t kSideAPCBBoardsHigh     =    6.32    *fgkcm;
  const Double_t kSideAPCBBoardsLen      =   24.00    *fgkcm;
  const Double_t kSideAPCBBoardsYTrans   =    0.75    *fgkcm;
  const Double_t kSideAPCBBoardsCu       =    0.3864;
  const Double_t kSideAPCBBoardsEpoxy    =    0.1486;
  const Double_t kSideAPCBBoardsPlast    =    0.0578;
  const Double_t kSideAPCBBoardsSteel    =    0.1521;
  const Double_t kSideAPCBBoardsPPS      =    0.2551;


  // Local variables
  Double_t xprof[kForwardTrayNpoints], yprof[kForwardTrayNpoints];
  Double_t xloc, yloc, zloc, alpharot, height;


  // The whole tray as an assembly
  TGeoVolumeAssembly *cableTrayA = new TGeoVolumeAssembly("ITSsupportSDDTrayA");
  

  // First create all needed shapes

  // The forward tray is very complex and deserves a dedicated method
  CreateSDDForwardTraySideA(cableTrayA,mgr);

  // The forward cover: a Xtru
  TGeoXtru *forwardCover = new TGeoXtru(2);
  forwardCover->SetName("ITSsuppSDDForwCover");

  xprof[0] = kForwardCoverWide/2;
  yprof[0] = kForwardCoverHeight;
  xprof[1] = xprof[0];
  yprof[1] = 0;
  xprof[2] = xprof[1] - kForwardCoverThick;
  yprof[2] = yprof[1];
  xprof[3] = xprof[2];
  yprof[3] = yprof[0] - kForwardCoverThick;

  // We did the right side, now reflex on the left side
  for (Int_t jp = 0; jp < 4; jp++) {
    xprof[4+jp] = -xprof[3-jp];
    yprof[4+jp] =  yprof[3-jp];
  }

  forwardCover->DefinePolygon(8, xprof, yprof);
  forwardCover->DefineSection(0, 0);
  forwardCover->DefineSection(1, kForwardCoverLength);

  // The external tray (as 0872/G/D/03): a Xtru
  TGeoXtru *externalTray = CreateSDDSSDTraysSideA(kExternTrayTotalLen,
						  kExternTrayTotalHeight);

  // The external covers: a Composite Shape
  TGeoCompositeShape *externCover = CreateTrayAExternalCover(kExternCoverLen);

  // Now the volumes inside it
  // The cooling manifold: four boxes
  TGeoBBox *coolManifPOM = new TGeoBBox(kSideACoolManifWide/2,
		 kSideACoolManifPOMFrac*kSideACoolManifHigh/2,
					kSideACoolManifLen/2);

  TGeoBBox *coolManifSteel = new TGeoBBox(kSideACoolManifWide/2,
		 kSideACoolManifSteelFrac*kSideACoolManifHigh/2,
					  kSideACoolManifLen/2);

  TGeoBBox *coolManifWater = new TGeoBBox(kSideACoolManifWide/2,
		 kSideACoolManifWaterFrac*kSideACoolManifHigh/2,
					  kSideACoolManifLen/2);

  TGeoBBox *coolManifAl = new TGeoBBox(kSideACoolManifWide/2,
		 kSideACoolManifAlFrac*kSideACoolManifHigh/2,
				       kSideACoolManifLen/2);

  // The cooling tubes: three Xtru's
  TGeoXtru *coolTubesPUR = new TGeoXtru(2);

  height = kSideACoolTubesHigh*kSideACoolTubesPURFrac;

  xprof[0] = kSideACoolManifLen;
  yprof[0] = kForwardTrayThick + kSideACoolTubesTrans;
  xprof[2] = kExternTrayZTrans + kForwardTrayTotalHeight*SinD(kTrayAZRot) +
	     kExternTrayTotalLen*CosD(kTrayAZRot) - xprof[0]/2;
  yprof[2] = kForwardTrayTotalHeight*(1 - CosD(kTrayAZRot)) +
	     kExternTrayYTrans - kExternTrayTotalHeight*CosD(kTrayAZRot) +
	     kExternTrayTotalLen*SinD(kTrayAZRot) + yprof[0];
  IntersectLines(              0 , xprof[0], yprof[0],
		 TanD(kTrayAZRot), xprof[2], yprof[2],
				   xprof[1], yprof[1]);
  xprof[3] = xprof[2] - height*SinD(kTrayAZRot);
  yprof[3] = yprof[2] + height*CosD(kTrayAZRot);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      height, xprof[4], yprof[4]);
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + height;

  coolTubesPUR->DefinePolygon(6, xprof, yprof);
  coolTubesPUR->DefineSection(0,-kSideACoolTubesWide/2);
  coolTubesPUR->DefineSection(1, kSideACoolTubesWide/2);

  TGeoXtru *coolTubesWater = new TGeoXtru(2);

  height = kSideACoolTubesHigh*kSideACoolTubesWaterFrac;

  xprof[0] = coolTubesPUR->GetX(5);
  yprof[0] = coolTubesPUR->GetY(5);
  xprof[1] = coolTubesPUR->GetX(4);
  yprof[1] = coolTubesPUR->GetY(4);
  xprof[2] = coolTubesPUR->GetX(3);
  yprof[2] = coolTubesPUR->GetY(3);
  xprof[3] = xprof[2] - height*SinD(kTrayAZRot);
  yprof[3] = yprof[2] + height*CosD(kTrayAZRot);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      height, xprof[4], yprof[4]);
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + height;

  coolTubesWater->DefinePolygon(6, xprof, yprof);
  coolTubesWater->DefineSection(0,-kSideACoolTubesWide/2);
  coolTubesWater->DefineSection(1, kSideACoolTubesWide/2);

  TGeoXtru *coolTubesAir = new TGeoXtru(2);

  height = kSideACoolTubesHigh*kSideACoolTubesAirFrac;

  xprof[0] = coolTubesWater->GetX(5);
  yprof[0] = coolTubesWater->GetY(5);
  xprof[1] = coolTubesWater->GetX(4);
  yprof[1] = coolTubesWater->GetY(4);
  xprof[2] = coolTubesWater->GetX(3);
  yprof[2] = coolTubesWater->GetY(3);
  xprof[3] = xprof[2] - height*SinD(kTrayAZRot);
  yprof[3] = yprof[2] + height*CosD(kTrayAZRot);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      height, xprof[4], yprof[4]);
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + height;

  coolTubesAir->DefinePolygon(6, xprof, yprof);
  coolTubesAir->DefineSection(0,-kSideACoolTubesWide/2);
  coolTubesAir->DefineSection(1, kSideACoolTubesWide/2);

  // The optical fiber connectors: three boxes
  TGeoBBox *optConnPBT = new TGeoBBox(kSideAOptConnWide/2,
		 kSideAOptConnPBTFrac*kSideACoolManifHigh/2,
				      kSideAOptConnLen/2);

  TGeoBBox *optConnSteel = new TGeoBBox(kSideAOptConnWide/2,
		 kSideAOptConnSteelFrac*kSideACoolManifHigh/2,
					kSideAOptConnLen/2);

  TGeoBBox *optConnAl = new TGeoBBox(kSideAOptConnWide/2,
		 kSideAOptConnAlFrac*kSideACoolManifHigh/2,
				     kSideAOptConnLen/2);

  // The optical fibers: a Xtru
  TGeoXtru *opticalFibs = new TGeoXtru(2);

  xprof[0] = kSideAOptConnLen;
  yprof[0] = coolTubesPUR->GetY(0);
  xprof[1] = coolTubesPUR->GetX(1);
  yprof[1] = coolTubesPUR->GetY(1);
  xprof[2] = coolTubesPUR->GetX(2);
  yprof[2] = coolTubesPUR->GetY(2);
  xprof[3] = xprof[2] - kSideAOptFibsHigh*SinD(kTrayAZRot);
  yprof[3] = yprof[2] + kSideAOptFibsHigh*CosD(kTrayAZRot);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      kSideAOptFibsHigh, xprof[4], yprof[4]);
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + kSideAOptFibsHigh;

  opticalFibs->DefinePolygon(6, xprof, yprof);
  opticalFibs->DefineSection(0,-kSideAOptFibsWide/2);
  opticalFibs->DefineSection(1, kSideAOptFibsWide/2);

  // The input cables: five boxes
  TGeoBBox *inputCabsCu = new TGeoBBox(kSideAInputCablesWide/2,
		   kSideAInputCablesCu*kSideAInputCablesHigh/2,
				       kSideAInputCablesLen/2);

  TGeoBBox *inputCabsPlast = new TGeoBBox(kSideAInputCablesWide/2,
		   kSideAInputCablesPlast*kSideAInputCablesHigh/2,
					  kSideAInputCablesLen/2);

  TGeoBBox *inputCabsAl = new TGeoBBox(kSideAInputCablesWide/2,
		   kSideAInputCablesAl*kSideAInputCablesHigh/2,
				       kSideAInputCablesLen/2);

  TGeoBBox *inputCabsKapton = new TGeoBBox(kSideAInputCablesWide/2,
		   kSideAInputCablesKapton*kSideAInputCablesHigh/2,
					   kSideAInputCablesLen/2);

  TGeoBBox *inputCabsPOLYAX = new TGeoBBox(kSideAInputCablesWide/2,
		   kSideAInputCablesPOLYAX*kSideAInputCablesHigh/2,
					   kSideAInputCablesLen/2);

  // The output cables: five Xtru
  TGeoXtru *outputCabsCu = new TGeoXtru(2);

  height = kSideAOutputCablesCu*kSideAOutputCablesHigh;

  xprof[0] = kSideAInputCablesLen/2 + kSideAPCBBoardsLen/2;
  yprof[0] = coolTubesAir->GetY(5);
  xprof[1] = coolTubesAir->GetX(4);
  yprof[1] = coolTubesAir->GetY(4);
  xprof[2] = coolTubesAir->GetX(3);
  yprof[2] = coolTubesAir->GetY(3);
  xprof[3] = xprof[2] - height*SinD(kTrayAZRot);
  yprof[3] = yprof[2] + height*CosD(kTrayAZRot);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      height, xprof[4], yprof[4]);
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + height;

  outputCabsCu->DefinePolygon(6, xprof, yprof);
  outputCabsCu->DefineSection(0,-kSideAOutputCablesWide/2);
  outputCabsCu->DefineSection(1, kSideAOutputCablesWide/2);

  TGeoXtru *outputCabsPlast = new TGeoXtru(2);

  height = kSideAOutputCablesPlast*kSideAOutputCablesHigh;

  xprof[0] = outputCabsCu->GetX(5);
  yprof[0] = outputCabsCu->GetY(5);
  xprof[1] = outputCabsCu->GetX(4);
  yprof[1] = outputCabsCu->GetY(4);
  xprof[2] = outputCabsCu->GetX(3);
  yprof[2] = outputCabsCu->GetY(3);
  xprof[3] = xprof[2] - height*SinD(kTrayAZRot);
  yprof[3] = yprof[2] + height*CosD(kTrayAZRot);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      height, xprof[4], yprof[4]);
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + height;

  outputCabsPlast->DefinePolygon(6, xprof, yprof);
  outputCabsPlast->DefineSection(0,-kSideAOutputCablesWide/2);
  outputCabsPlast->DefineSection(1, kSideAOutputCablesWide/2);

  TGeoXtru *outputCabsAl = new TGeoXtru(2);

  height = kSideAOutputCablesAl*kSideAOutputCablesHigh;

  xprof[0] = outputCabsPlast->GetX(5);
  yprof[0] = outputCabsPlast->GetY(5);
  xprof[1] = outputCabsPlast->GetX(4);
  yprof[1] = outputCabsPlast->GetY(4);
  xprof[2] = outputCabsPlast->GetX(3);
  yprof[2] = outputCabsPlast->GetY(3);
  xprof[3] = xprof[2] - height*SinD(kTrayAZRot);
  yprof[3] = yprof[2] + height*CosD(kTrayAZRot);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      height, xprof[4], yprof[4]);
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + height;

  outputCabsAl->DefinePolygon(6, xprof, yprof);
  outputCabsAl->DefineSection(0,-kSideAOutputCablesWide/2);
  outputCabsAl->DefineSection(1, kSideAOutputCablesWide/2);

  TGeoXtru *outputCabsKapton = new TGeoXtru(2);

  height = kSideAOutputCablesKapton*kSideAOutputCablesHigh;

  xprof[0] = outputCabsAl->GetX(5);
  yprof[0] = outputCabsAl->GetY(5);
  xprof[1] = outputCabsAl->GetX(4);
  yprof[1] = outputCabsAl->GetY(4);
  xprof[2] = outputCabsAl->GetX(3);
  yprof[2] = outputCabsAl->GetY(3);
  xprof[3] = xprof[2] - height*SinD(kTrayAZRot);
  yprof[3] = yprof[2] + height*CosD(kTrayAZRot);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      height, xprof[4], yprof[4]);
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + height;

  outputCabsKapton->DefinePolygon(6, xprof, yprof);
  outputCabsKapton->DefineSection(0,-kSideAOutputCablesWide/2);
  outputCabsKapton->DefineSection(1, kSideAOutputCablesWide/2);

  TGeoXtru *outputCabsPOLYAX = new TGeoXtru(2);

  height = kSideAOutputCablesPOLYAX*kSideAOutputCablesHigh;

  xprof[0] = outputCabsKapton->GetX(5);
  yprof[0] = outputCabsKapton->GetY(5);
  xprof[1] = outputCabsKapton->GetX(4);
  yprof[1] = outputCabsKapton->GetY(4);
  xprof[2] = outputCabsKapton->GetX(3);
  yprof[2] = outputCabsKapton->GetY(3);
  xprof[3] = xprof[2] - height*SinD(kTrayAZRot);
  yprof[3] = yprof[2] + height*CosD(kTrayAZRot);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      height, xprof[4], yprof[4]);
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + height;

  outputCabsPOLYAX->DefinePolygon(6, xprof, yprof);
  outputCabsPOLYAX->DefineSection(0,-kSideAOutputCablesWide/2);
  outputCabsPOLYAX->DefineSection(1, kSideAOutputCablesWide/2);

  // The PCB boards: five boxes
  TGeoBBox *pcbBoardsCu = new TGeoBBox(kSideAPCBBoardsWide/2,
		     kSideAPCBBoardsCu*kSideAPCBBoardsHigh/2,
				       kSideAPCBBoardsLen/2);

  TGeoBBox *pcbBoardsEpoxy = new TGeoBBox(kSideAPCBBoardsWide/2,
		     kSideAPCBBoardsEpoxy*kSideAPCBBoardsHigh/2,
					  kSideAPCBBoardsLen/2);

  TGeoBBox *pcbBoardsPlast = new TGeoBBox(kSideAPCBBoardsWide/2,
		     kSideAPCBBoardsPlast*kSideAPCBBoardsHigh/2,
					  kSideAPCBBoardsLen/2);

  TGeoBBox *pcbBoardsSteel = new TGeoBBox(kSideAPCBBoardsWide/2,
		     kSideAPCBBoardsSteel*kSideAPCBBoardsHigh/2,
					  kSideAPCBBoardsLen/2);

  TGeoBBox *pcbBoardsPPS = new TGeoBBox(kSideAPCBBoardsWide/2,
		     kSideAPCBBoardsPPS*kSideAPCBBoardsHigh/2,
					kSideAPCBBoardsLen/2);


  // We have all shapes: now create the real volumes
  TGeoMedium *medAl     = mgr->GetMedium("ITS_ALUMINUM$");
  TGeoMedium *medAntic  = mgr->GetMedium("ITS_ANTICORODAL$");
  TGeoMedium *medPOM    = mgr->GetMedium("ITS_POLYOXYMETHYLENE$");
  TGeoMedium *medSteel  = mgr->GetMedium("ITS_INOX$");
  TGeoMedium *medWater  = mgr->GetMedium("ITS_WATER$");
  TGeoMedium *medPUR    = mgr->GetMedium("ITS_POLYURETHANE$");
  TGeoMedium *medAir    = mgr->GetMedium("ITS_AIR$");
  TGeoMedium *medPBT    = mgr->GetMedium("ITS_PBT$");
  TGeoMedium *medOptFib = mgr->GetMedium("ITS_SDD OPTICFIB$");
  TGeoMedium *medCu     = mgr->GetMedium("ITS_COPPER$");
  TGeoMedium *medKapton = mgr->GetMedium("ITS_SDDKAPTON (POLYCH2)$");
  TGeoMedium *medPOLYAX = mgr->GetMedium("ITS_POLYAX$");
  TGeoMedium *medPPS    = mgr->GetMedium("ITS_PPS$");
  TGeoMedium *medEpoxy  = mgr->GetMedium("ITS_EPOXY$");

  TGeoVolume *forwardTrayCover = new TGeoVolume("ITSsuppSDDSideAForwTrayCover",
						forwardCover, medAl);

  forwardTrayCover->SetVisibility(kTRUE);
  forwardTrayCover->SetLineColor(kMagenta+1); // Purple
  forwardTrayCover->SetLineWidth(1);
  forwardTrayCover->SetFillColor(forwardTrayCover->GetLineColor());
  forwardTrayCover->SetFillStyle(4000); // 0% transparent

  TGeoVolume *externalTraySDD = new TGeoVolume("ITSsuppSDDSideAExternalTray",
					       externalTray, medAl);

  externalTraySDD->SetVisibility(kTRUE);
  externalTraySDD->SetLineColor(6); // Purple
  externalTraySDD->SetLineWidth(1);
  externalTraySDD->SetFillColor(externalTraySDD->GetLineColor());
  externalTraySDD->SetFillStyle(4000); // 0% transparent

  TGeoVolume *externTrayCover = new TGeoVolume("ITSsuppSDDSideAExtTrayCover",
					       externCover, medAntic);

  externTrayCover->SetVisibility(kTRUE);
  externTrayCover->SetLineColor(kMagenta+1); // Purple
  externTrayCover->SetLineWidth(1);
  externTrayCover->SetFillColor(externTrayCover->GetLineColor());
  externTrayCover->SetFillStyle(4000); // 0% transparent

  TGeoVolume *pomCoolManif = new TGeoVolume("ITSsuppSDDSideACoolManifPOM",
					    coolManifPOM, medPOM);

  pomCoolManif->SetVisibility(kTRUE);
  pomCoolManif->SetLineColor(kRed); // Red
  pomCoolManif->SetLineWidth(1);
  pomCoolManif->SetFillColor(pomCoolManif->GetLineColor());
  pomCoolManif->SetFillStyle(4000); // 0% transparent

  TGeoVolume *steelCoolManif = new TGeoVolume("ITSsuppSDDSideACoolManifSteel",
					      coolManifSteel, medSteel);

  steelCoolManif->SetVisibility(kTRUE);
  steelCoolManif->SetLineColor(kBlue); // Blue
  steelCoolManif->SetLineWidth(1);
  steelCoolManif->SetFillColor(steelCoolManif->GetLineColor());
  steelCoolManif->SetFillStyle(4000); // 0% transparent

  TGeoVolume *waterCoolManif = new TGeoVolume("ITSsuppSDDSideACoolManifWater",
					      coolManifWater, medWater);

  waterCoolManif->SetVisibility(kTRUE);
  waterCoolManif->SetLineColor(33); // Light Blue
  waterCoolManif->SetLineWidth(1);
  waterCoolManif->SetFillColor(waterCoolManif->GetLineColor());
  waterCoolManif->SetFillStyle(4000); // 0% transparent

  TGeoVolume *alCoolManif = new TGeoVolume("ITSsuppSDDSideACoolManifAl",
					   coolManifAl, medAl);

  alCoolManif->SetVisibility(kTRUE);
  alCoolManif->SetLineColor(6); // Purple
  alCoolManif->SetLineWidth(1);
  alCoolManif->SetFillColor(alCoolManif->GetLineColor());
  alCoolManif->SetFillStyle(4000); // 0% transparent

  TGeoVolume *purCoolTubes = new TGeoVolume("ITSsuppSDDSideACoolTubesPUR",
					    coolTubesPUR, medPUR);

  purCoolTubes->SetVisibility(kTRUE);
  purCoolTubes->SetLineColor(kRed); // Red
  purCoolTubes->SetLineWidth(1);
  purCoolTubes->SetFillColor(purCoolTubes->GetLineColor());
  purCoolTubes->SetFillStyle(4000); // 0% transparent

  TGeoVolume *waterCoolTubes = new TGeoVolume("ITSsuppSDDSideACoolTubesWater",
					      coolTubesWater, medWater);

  waterCoolTubes->SetVisibility(kTRUE);
  waterCoolTubes->SetLineColor(33); // Light Blue
  waterCoolTubes->SetLineWidth(1);
  waterCoolTubes->SetFillColor(waterCoolTubes->GetLineColor());
  waterCoolTubes->SetFillStyle(4000); // 0% transparent

  TGeoVolume *airCoolTubes = new TGeoVolume("ITSsuppSDDSideACoolTubesAir",
					    coolTubesAir, medAir);

  airCoolTubes->SetVisibility(kTRUE);
  airCoolTubes->SetLineColor(41);
  airCoolTubes->SetLineWidth(1);
  airCoolTubes->SetFillColor(airCoolTubes->GetLineColor());
  airCoolTubes->SetFillStyle(4000); // 0% transparent

  TGeoVolume *pbtOptConn = new TGeoVolume("ITSsuppSDDSideAOptConnPBT",
					  optConnPBT, medPBT);

  pbtOptConn->SetVisibility(kTRUE);
  pbtOptConn->SetLineColor(kRed); // Red
  pbtOptConn->SetLineWidth(1);
  pbtOptConn->SetFillColor(pbtOptConn->GetLineColor());
  pbtOptConn->SetFillStyle(4000); // 0% transparent

  TGeoVolume *steelOptConn = new TGeoVolume("ITSsuppSDDSideAOptConnSteel",
					    optConnSteel, medSteel);

  steelOptConn->SetVisibility(kTRUE);
  steelOptConn->SetLineColor(kBlue); // Blue
  steelOptConn->SetLineWidth(1);
  steelOptConn->SetFillColor(steelOptConn->GetLineColor());
  steelOptConn->SetFillStyle(4000); // 0% transparent

  TGeoVolume *alOptConn = new TGeoVolume("ITSsuppSDDSideAOptConnAl",
					 optConnAl, medAl);

  alOptConn->SetVisibility(kTRUE);
  alOptConn->SetLineColor(6); // Purple
  alOptConn->SetLineWidth(1);
  alOptConn->SetFillColor(alOptConn->GetLineColor());
  alOptConn->SetFillStyle(4000); // 0% transparent

  TGeoVolume *optFibs = new TGeoVolume("ITSsuppSDDSideAOptFibs",
				       opticalFibs, medOptFib);

  optFibs->SetVisibility(kTRUE);
  optFibs->SetLineColor(kOrange+2); // Orange
  optFibs->SetLineWidth(1);
  optFibs->SetFillColor(optFibs->GetLineColor());
  optFibs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *cuInputCabs = new TGeoVolume("ITSsuppSDDSideAInputCabsCu",
					   inputCabsCu, medCu);

  cuInputCabs->SetVisibility(kTRUE);
  cuInputCabs->SetLineColor(kBlack); // Black
  cuInputCabs->SetLineWidth(1);
  cuInputCabs->SetFillColor(cuInputCabs->GetLineColor());
  cuInputCabs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *plastInputCabs = new TGeoVolume("ITSsuppSDDSideAInputCabsPlast",
					      inputCabsPlast, medPUR);

  plastInputCabs->SetVisibility(kTRUE);
  plastInputCabs->SetLineColor(kRed); // Red
  plastInputCabs->SetLineWidth(1);
  plastInputCabs->SetFillColor(plastInputCabs->GetLineColor());
  plastInputCabs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *alInputCabs = new TGeoVolume("ITSsuppSDDSideAInputCabsAl",
					   inputCabsAl, medAl);

  alInputCabs->SetVisibility(kTRUE);
  alInputCabs->SetLineColor(6); // Purple
  alInputCabs->SetLineWidth(1);
  alInputCabs->SetFillColor(alInputCabs->GetLineColor());
  alInputCabs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *kaptonInputCabs = new TGeoVolume("ITSsuppSDDSideAInputCabsKapton",
					       inputCabsKapton, medKapton);

  kaptonInputCabs->SetVisibility(kTRUE);
  kaptonInputCabs->SetLineColor(14); // 
  kaptonInputCabs->SetLineWidth(1);
  kaptonInputCabs->SetFillColor(kaptonInputCabs->GetLineColor());
  kaptonInputCabs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *polyaxInputCabs = new TGeoVolume("ITSsuppSDDSideAInputCabsPOLYAX",
					       inputCabsPOLYAX, medPOLYAX);

  polyaxInputCabs->SetVisibility(kTRUE);
  polyaxInputCabs->SetLineColor(34); // 
  polyaxInputCabs->SetLineWidth(1);
  polyaxInputCabs->SetFillColor(polyaxInputCabs->GetLineColor());
  polyaxInputCabs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *cuOutputCabs = new TGeoVolume("ITSsuppSDDSideAOutputCabsCu",
					    outputCabsCu, medCu);

  cuOutputCabs->SetVisibility(kTRUE);
  cuOutputCabs->SetLineColor(kBlack); // Black
  cuOutputCabs->SetLineWidth(1);
  cuOutputCabs->SetFillColor(cuOutputCabs->GetLineColor());
  cuOutputCabs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *plastOutputCabs = new TGeoVolume("ITSsuppSDDSideAOutputCabsPlast",
					       outputCabsPlast, medPUR);

  plastOutputCabs->SetVisibility(kTRUE);
  plastOutputCabs->SetLineColor(kRed); // Red
  plastOutputCabs->SetLineWidth(1);
  plastOutputCabs->SetFillColor(plastOutputCabs->GetLineColor());
  plastOutputCabs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *alOutputCabs = new TGeoVolume("ITSsuppSDDSideAOutputCabsAl",
					    outputCabsAl, medAl);

  alOutputCabs->SetVisibility(kTRUE);
  alOutputCabs->SetLineColor(6); // Purple
  alOutputCabs->SetLineWidth(1);
  alOutputCabs->SetFillColor(alOutputCabs->GetLineColor());
  alOutputCabs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *kaptonOutputCabs = new TGeoVolume("ITSsuppSDDSideAOutputCabsKapton",
						outputCabsKapton, medKapton);

  kaptonOutputCabs->SetVisibility(kTRUE);
  kaptonOutputCabs->SetLineColor(14); // 
  kaptonOutputCabs->SetLineWidth(1);
  kaptonOutputCabs->SetFillColor(kaptonOutputCabs->GetLineColor());
  kaptonOutputCabs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *polyaxOutputCabs = new TGeoVolume("ITSsuppSDDSideAOutputCabsPOLYAX",
						outputCabsPOLYAX, medPOLYAX);

  polyaxOutputCabs->SetVisibility(kTRUE);
  polyaxOutputCabs->SetLineColor(34); // 
  polyaxOutputCabs->SetLineWidth(1);
  polyaxOutputCabs->SetFillColor(polyaxOutputCabs->GetLineColor());
  polyaxOutputCabs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *cuPCBBoards = new TGeoVolume("ITSsuppSDDSideAPCBBoardsCu",
					   pcbBoardsCu, medCu);

  cuPCBBoards->SetVisibility(kTRUE);
  cuPCBBoards->SetLineColor(kBlack); // Black
  cuPCBBoards->SetLineWidth(1);
  cuPCBBoards->SetFillColor(cuPCBBoards->GetLineColor());
  cuPCBBoards->SetFillStyle(4000); // 0% transparent

  TGeoVolume *epoxyPCBBoards = new TGeoVolume("ITSsuppSDDSideAPCBBoardsEpoxy",
					      pcbBoardsEpoxy, medEpoxy);

  epoxyPCBBoards->SetVisibility(kTRUE);
  epoxyPCBBoards->SetLineColor(22); //
  epoxyPCBBoards->SetLineWidth(1);
  epoxyPCBBoards->SetFillColor(epoxyPCBBoards->GetLineColor());
  epoxyPCBBoards->SetFillStyle(4000); // 0% transparent

  TGeoVolume *plastPCBBoards = new TGeoVolume("ITSsuppSDDSideAPCBBoardsPlast",
					      pcbBoardsPlast, medPUR);

  plastPCBBoards->SetVisibility(kTRUE);
  plastPCBBoards->SetLineColor(kRed); // Red
  plastPCBBoards->SetLineWidth(1);
  plastPCBBoards->SetFillColor(plastPCBBoards->GetLineColor());
  plastPCBBoards->SetFillStyle(4000); // 0% transparent

  TGeoVolume *steelPCBBoards = new TGeoVolume("ITSsuppSDDSideAPCBBoardsSteel",
					      pcbBoardsSteel, medSteel);

  steelPCBBoards->SetVisibility(kTRUE);
  steelPCBBoards->SetLineColor(kBlue); // Blue
  steelPCBBoards->SetLineWidth(1);
  steelPCBBoards->SetFillColor(steelPCBBoards->GetLineColor());
  steelPCBBoards->SetFillStyle(4000); // 0% transparent

  TGeoVolume *ppsPCBBoards = new TGeoVolume("ITSsuppSDDSideAPCBBoardsPPS",
					    pcbBoardsPPS, medPPS);

  ppsPCBBoards->SetVisibility(kTRUE);
  ppsPCBBoards->SetLineColor(kGreen); // Green
  ppsPCBBoards->SetLineWidth(1);
  ppsPCBBoards->SetFillColor(ppsPCBBoards->GetLineColor());
  ppsPCBBoards->SetFillStyle(4000); // 0% transparent


  // Now build up the tray
  yloc = kForwardTrayTotalHeight - forwardCover->GetY(3);
  zloc = kForwardTrayUpperLength - kForwardCoverLength;
  cableTrayA->AddNode(forwardTrayCover, 1,
		      new TGeoTranslation( 0, yloc, zloc) );

  Double_t totalhi = kExternTrayTotalHeight + kExternCoverThick
		   - kExternCoverYTrans;

  yloc = totalhi*(1 - CosD(kTrayAZRot)) + kExternTrayYTrans -
	 kExternTrayTotalHeight*CosD(kTrayAZRot);
  zloc = kExternTrayZTrans + totalhi*SinD(kTrayAZRot);
  cableTrayA->AddNode(externalTraySDD, 1,
		      new TGeoCombiTrans( 0, yloc, zloc,
		      new TGeoRotation("", 0,-kTrayAZRot, 0)        ) );

  yloc = kExternCoverThick*(1 - CosD(kTrayAZRot)) + kExternTrayYTrans -
	 kExternCoverYTrans*CosD(kTrayAZRot)/2-0.01;
  zloc = kExternTrayZTrans + kExternCoverThick*SinD(kTrayAZRot);
  cableTrayA->AddNode(externTrayCover,1,
		      new TGeoCombiTrans( 0, yloc, zloc,
		      new TGeoRotation("", 0,-kTrayAZRot, 0)        ) );

  yloc = kForwardTrayThick + coolManifPOM->GetDY();
  zloc = coolManifPOM->GetDZ();
  cableTrayA->AddNode(pomCoolManif, 1,
		      new TGeoTranslation( 0, yloc, zloc) );

  yloc += coolManifPOM->GetDY() + coolManifSteel->GetDY();
  cableTrayA->AddNode(steelCoolManif, 1,
		      new TGeoTranslation( 0, yloc, zloc) );

  yloc += coolManifSteel->GetDY() + coolManifWater->GetDY();
  cableTrayA->AddNode(waterCoolManif, 1,
		      new TGeoTranslation( 0, yloc, zloc) );

  yloc += coolManifWater->GetDY() + coolManifAl->GetDY();
  cableTrayA->AddNode(alCoolManif, 1,
		      new TGeoTranslation( 0, yloc, zloc) );

  cableTrayA->AddNode(purCoolTubes,1,
		      new TGeoCombiTrans( 0, 0, 0,
		      new TGeoRotation("",-90, 90, 90)        ) );
  cableTrayA->AddNode(waterCoolTubes,1,
		      new TGeoCombiTrans( 0, 0, 0,
		      new TGeoRotation("",-90, 90, 90)        ) );
  cableTrayA->AddNode(airCoolTubes,1,
		      new TGeoCombiTrans( 0, 0, 0,
		      new TGeoRotation("",-90, 90, 90)        ) );

  xloc = coolManifPOM->GetDX() + optConnPBT->GetDX();
  yloc = kForwardTrayThick + optConnPBT->GetDY();
  zloc = optConnPBT->GetDZ();
  cableTrayA->AddNode(pbtOptConn, 1,
		      new TGeoTranslation( xloc, yloc, zloc) );
  cableTrayA->AddNode(pbtOptConn, 2,
		      new TGeoTranslation(-xloc, yloc, zloc) );

  yloc += optConnPBT->GetDY() + optConnSteel->GetDY();
  cableTrayA->AddNode(steelOptConn, 1,
		      new TGeoTranslation( xloc, yloc, zloc) );
  cableTrayA->AddNode(steelOptConn, 2,
		      new TGeoTranslation(-xloc, yloc, zloc) );

  yloc += optConnSteel->GetDY() + optConnAl->GetDY();
  cableTrayA->AddNode(alOptConn, 1,
		      new TGeoTranslation( xloc, yloc, zloc) );
  cableTrayA->AddNode(alOptConn, 2,
		      new TGeoTranslation(-xloc, yloc, zloc) );


  xloc = kSideACoolTubesWide/2 + kSideAOptFibsWide/2;
  cableTrayA->AddNode(optFibs,1,
		      new TGeoCombiTrans( xloc, 0, 0,
		      new TGeoRotation("",-90, 90, 90)        ) );
  cableTrayA->AddNode(optFibs,2,
		      new TGeoCombiTrans(-xloc, 0, 0,
		      new TGeoRotation("",-90, 90, 90)        ) );

  yloc = kForwardTrayTotalHeight - forwardCover->GetY(3) -
	 kSideAInputCablesYTrans - inputCabsPOLYAX->GetDY();
  zloc = inputCabsPOLYAX->GetDZ();
  cableTrayA->AddNode(polyaxInputCabs, 1,
		      new TGeoTranslation( 0, yloc, zloc) );

  yloc -= (inputCabsPOLYAX->GetDY() + inputCabsKapton->GetDY());
  cableTrayA->AddNode(kaptonInputCabs, 1,
		      new TGeoTranslation( 0, yloc, zloc) );

  yloc -= (inputCabsKapton->GetDY() + inputCabsAl->GetDY());
  cableTrayA->AddNode(alInputCabs, 1,
		      new TGeoTranslation( 0, yloc, zloc) );

  yloc -= (inputCabsAl->GetDY() + inputCabsPlast->GetDY());
  cableTrayA->AddNode(plastInputCabs, 1,
		      new TGeoTranslation( 0, yloc, zloc) );

  yloc -= (inputCabsPlast->GetDY() + inputCabsCu->GetDY());
  cableTrayA->AddNode(cuInputCabs, 1,
		      new TGeoTranslation( 0, yloc, zloc) );

  yloc -= (inputCabsCu->GetDY()+pcbBoardsPPS->GetDY()+kSideAPCBBoardsYTrans);
  zloc += pcbBoardsPPS->GetDZ();
  cableTrayA->AddNode(ppsPCBBoards, 1,
		      new TGeoTranslation( 0, yloc, zloc) );

  yloc -= (pcbBoardsPPS->GetDY()+pcbBoardsSteel->GetDY());
  cableTrayA->AddNode(steelPCBBoards, 1,
		      new TGeoTranslation( 0, yloc, zloc) );

  yloc -= (pcbBoardsSteel->GetDY()+pcbBoardsPlast->GetDY());
  cableTrayA->AddNode(plastPCBBoards, 1,
		      new TGeoTranslation( 0, yloc, zloc) );

  yloc -= (pcbBoardsPlast->GetDY()+pcbBoardsEpoxy->GetDY());
  cableTrayA->AddNode(epoxyPCBBoards, 1,
		      new TGeoTranslation( 0, yloc, zloc) );

  yloc -= (pcbBoardsEpoxy->GetDY()+pcbBoardsCu->GetDY());
  cableTrayA->AddNode(cuPCBBoards, 1,
		      new TGeoTranslation( 0, yloc, zloc) );

  cableTrayA->AddNode(cuOutputCabs,1,
		      new TGeoCombiTrans( 0, 0, 0,
		      new TGeoRotation("",-90, 90, 90)        ) );
  cableTrayA->AddNode(plastOutputCabs,1,
		      new TGeoCombiTrans( 0, 0, 0,
		      new TGeoRotation("",-90, 90, 90)        ) );
  cableTrayA->AddNode(alOutputCabs,1,
		      new TGeoCombiTrans( 0, 0, 0,
		      new TGeoRotation("",-90, 90, 90)        ) );
  cableTrayA->AddNode(kaptonOutputCabs,1,
		      new TGeoCombiTrans( 0, 0, 0,
		      new TGeoRotation("",-90, 90, 90)        ) );
  cableTrayA->AddNode(polyaxOutputCabs,1,
		      new TGeoCombiTrans( 0, 0, 0,
		      new TGeoRotation("",-90, 90, 90)        ) );


  // Finally put everything in the mother volume
  Double_t rforw = kTrayARTrans + kExternTrayTotalHeight +
		   kExternCoverSideThick -
		   kForwardTrayTailHeight;

  alpharot = -kTrayAFirstRotAng;
  xloc = rforw*SinD(alpharot);
  yloc = rforw*CosD(alpharot);
  zloc = kTrayAZTrans + kTrayAZToSupportRing - kForwardTrayUpperLength;

  moth->AddNode(cableTrayA,1,
			    new TGeoCombiTrans( xloc, yloc, zloc,
			    new TGeoRotation("",-alpharot,0,0)   )   );

  alpharot += 180;
  xloc = rforw*SinD(alpharot);
  yloc = rforw*CosD(alpharot);
  moth->AddNode(cableTrayA,2,
			    new TGeoCombiTrans( xloc, yloc, zloc,
			    new TGeoRotation("",-alpharot,0,0)   )   );

  alpharot = kTrayAFirstRotAng + 2*kTrayASecondRotAng;
  xloc = rforw*SinD(alpharot);
  yloc = rforw*CosD(alpharot);
  moth->AddNode(cableTrayA,3,
			    new TGeoCombiTrans( xloc, yloc, zloc,
			    new TGeoRotation("",-alpharot,0,0)   )   );

  alpharot += 180;
  xloc = rforw*SinD(alpharot);
  yloc = rforw*CosD(alpharot);
  moth->AddNode(cableTrayA,4,
			    new TGeoCombiTrans( xloc, yloc, zloc,
			    new TGeoRotation("",-alpharot,0,0)   )   );


  return;
}

//______________________________________________________________________
void AliITSv11GeometrySupport::SDDCableTraysSideC(TGeoVolume *moth,
					    const TGeoManager *mgr){
//
// Creates the SDD cable trays which are outside the ITS support cones
// but still inside the TPC on Side C
// (part of this code is taken or anyway inspired to ServicesCableSupport
// method of AliITSv11GeometrySupport.cxx,v 1.9 2007/06/06)
//
// Input:
//         moth : the TGeoVolume owing the volume structure
//         mgr  : the GeoManager (default gGeoManager)
// Output:
//
// Created:         ???       Bjorn S. Nilsen
// Updated:      17 Apr 2010  Mario Sitta
//
// Technical data are taken from AutoCAD drawings and other (oral)
// information given by F.Tosello
//

  // Dimensions and positions of the C-Side Cable Tray
  // (Change accordingly to CreateSDDSSDTraysSideC !)
  const Int_t    kNumTraySideC           =    4;

  const Double_t kSideCHalfThick         =    0.100   *fgkcm;
  const Double_t kSideCLength1           =  172.800   *fgkcm;
  const Double_t kSideCLength2           =  189.300   *fgkcm;
  const Double_t kBarCoolRmax            =    0.4     *fgkcm;
  const Double_t kXShiftBarCool          =   13.00    *fgkcm;

  const Double_t kSideCFoldAngle         =    5.00 *fgkDegree;

  // Dimensions and positions of the Cable Tray elements
  const Double_t kSideCCoolManifHalfX    =    4.25    *fgkcm;
  const Double_t kSideCCoolManifHalfY    =    4.03    *fgkcm;
  const Double_t kSideCCoolManifHalfZ    =    2.17    *fgkcm;
  const Double_t kSideCCoolManifPOMFrac  =    0.0051;
  const Double_t kSideCCoolManifSteelFrac=    0.8502;
  const Double_t kSideCCoolManifWaterFrac=    0.0868;
  const Double_t kSideCCoolManifAlFrac   =    0.0579;

  const Double_t kSideCCoolTubesHigh     =    1.88    *fgkcm;
  const Double_t kSideCCoolTubesTrans    =    0.85    *fgkcm;
  const Double_t kSideCCoolTubesPURFrac  =    0.5884;
  const Double_t kSideCCoolTubesWaterFrac=    0.4114;
  const Double_t kSideCCoolTubesAirFrac  =    0.0002;

  const Double_t kSideCOptConnHalfX      =    0.90    *fgkcm;
  const Double_t kSideCOptConnHalfZ      =    1.37    *fgkcm;
  const Double_t kSideCOptConnPBTFrac    =    0.6798;
  const Double_t kSideCOptConnSteelFrac  =    0.2421;
  const Double_t kSideCOptConnAlFrac     =    0.0781;

  const Double_t kSideCOptFibsWide       =    0.71    *fgkcm;
  const Double_t kSideCOptFibsHigh       =    3.20    *fgkcm;
  const Double_t kSideCOptFibsTrans      =    0.20    *fgkcm;

  const Double_t kSideCInputCablesLen    =   31.45    *fgkcm;
  const Double_t kSideCInputCablesWide   =   12.50    *fgkcm;
  const Double_t kSideCInputCablesHigh   =    0.95    *fgkcm;
  const Double_t kSideCInputCablesTrans  =    1.15    *fgkcm;
  const Double_t kSideCInputCablesCu     =    0.7405;
  const Double_t kSideCInputCablesPlast  =    0.1268;
  const Double_t kSideCInputCablesAl     =    0.0057;
  const Double_t kSideCInputCablesKapton =    0.0172;
  const Double_t kSideCInputCablesPOLYAX =    0.1098;

  const Double_t kSideCOutputCablesX0    =   27.40    *fgkcm;
  const Double_t kSideCOutputCablesWide  =    8.30    *fgkcm;
  const Double_t kSideCOutputCablesHigh  =    1.18    *fgkcm;
  const Double_t kSideCOutputCablesCu    =    0.6775;
  const Double_t kSideCOutputCablesPlast =    0.1613;
  const Double_t kSideCOutputCablesAl    =    0.0078;
  const Double_t kSideCOutputCablesKapton=    0.0234;
  const Double_t kSideCOutputCablesPOLYAX=    0.1300;

  const Double_t kSideCPCBBoardsHalfX    =    6.30    *fgkcm;
  const Double_t kSideCPCBBoardsHalfY    =    2.00    *fgkcm;
  const Double_t kSideCPCBBoardsHalfZ    =   21.93    *fgkcm;
  const Double_t kSideCPCBBoardsCu       =    0.3864;
  const Double_t kSideCPCBBoardsEpoxy    =    0.1491;
  const Double_t kSideCPCBBoardsPlast    =    0.0579;
  const Double_t kSideCPCBBoardsSteel    =    0.1517;
  const Double_t kSideCPCBBoardsPPS      =    0.2549;

  // Overall position and rotation of the C-Side Cable Trays
  const Double_t kTraySideCRPos          =   45.30    *fgkcm;
  const Double_t kTraySideCZPos          = -102.40    *fgkcm;
  const Double_t kTraySideCAlphaRot[kNumTraySideC]  = {    -23.0,      59.0,
    /* from SSD tray position */		       180.-23.0, 180.+59.0};


  // Local variables
  Double_t xprof[6], yprof[6];
  Double_t height, xloc, yloc, zloc, alpharot, alphafold;


  // The assembly holding the metallic structure
  TGeoVolumeAssembly *trayStructure = CreateSDDSSDTraysSideC("ITSsupportSDDTrayC");

  // Now the volumes inside it
  // The cooling manifold: four boxes
  // (X and Z are inverted on tray reference system)
  TGeoBBox *coolManifPOM = new TGeoBBox(kSideCCoolManifHalfZ,
		 kSideCCoolManifPOMFrac*kSideCCoolManifHalfY,
					kSideCCoolManifHalfX);

  TGeoBBox *coolManifSteel = new TGeoBBox(kSideCCoolManifHalfZ,
		 kSideCCoolManifSteelFrac*kSideCCoolManifHalfY,
					  kSideCCoolManifHalfX);

  TGeoBBox *coolManifWater = new TGeoBBox(kSideCCoolManifHalfZ,
		 kSideCCoolManifWaterFrac*kSideCCoolManifHalfY,
					  kSideCCoolManifHalfX);

  TGeoBBox *coolManifAl = new TGeoBBox(kSideCCoolManifHalfZ,
		 kSideCCoolManifAlFrac*kSideCCoolManifHalfY,
				       kSideCCoolManifHalfX);

  // The cooling tubes: three Xtru's
  alpharot = kSideCFoldAngle*TMath::DegToRad();

  TGeoXtru *coolTubesPUR = new TGeoXtru(2);

  height = kSideCCoolTubesHigh*kSideCCoolTubesPURFrac;

  xprof[0] = 2*kSideCCoolManifHalfZ;
  yprof[0] = 2*kSideCHalfThick + kSideCCoolTubesTrans;
  xprof[1] = kSideCLength1;
  yprof[1] = yprof[0];
  xprof[2] = xprof[1] + kSideCLength2*TMath::Cos(alpharot);
  yprof[2] = yprof[1] + kSideCLength2*TMath::Sin(alpharot);
  xprof[3] = xprof[2] - height*TMath::Sin(alpharot);
  yprof[3] = yprof[2] + height*TMath::Cos(alpharot);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      height, xprof[4], yprof[4]);
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + height;

  coolTubesPUR->DefinePolygon(6, xprof, yprof);
  coolTubesPUR->DefineSection(0,-kSideCCoolManifHalfX);
  coolTubesPUR->DefineSection(1, kSideCCoolManifHalfX);

  TGeoXtru *coolTubesWater = new TGeoXtru(2);

  height = kSideCCoolTubesHigh*kSideCCoolTubesWaterFrac;

  xprof[0] = coolTubesPUR->GetX(5);
  yprof[0] = coolTubesPUR->GetY(5);
  xprof[1] = coolTubesPUR->GetX(4);
  yprof[1] = coolTubesPUR->GetY(4);
  xprof[2] = coolTubesPUR->GetX(3);
  yprof[2] = coolTubesPUR->GetY(3);
  xprof[3] = xprof[2] - height*TMath::Sin(alpharot);
  yprof[3] = yprof[2] + height*TMath::Cos(alpharot);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      height, xprof[4], yprof[4]);
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + height;

  coolTubesWater->DefinePolygon(6, xprof, yprof);
  coolTubesWater->DefineSection(0,-kSideCCoolManifHalfX);
  coolTubesWater->DefineSection(1, kSideCCoolManifHalfX);

  TGeoXtru *coolTubesAir = new TGeoXtru(2);

  height = kSideCCoolTubesHigh*kSideCCoolTubesAirFrac;

  xprof[0] = coolTubesWater->GetX(5);
  yprof[0] = coolTubesWater->GetY(5);
  xprof[1] = coolTubesWater->GetX(4);
  yprof[1] = coolTubesWater->GetY(4);
  xprof[2] = coolTubesWater->GetX(3);
  yprof[2] = coolTubesWater->GetY(3);
  xprof[3] = xprof[2] - height*TMath::Sin(alpharot);
  yprof[3] = yprof[2] + height*TMath::Cos(alpharot);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      height, xprof[4], yprof[4]);
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + height;

  coolTubesAir->DefinePolygon(6, xprof, yprof);
  coolTubesAir->DefineSection(0,-kSideCCoolManifHalfX);
  coolTubesAir->DefineSection(1, kSideCCoolManifHalfX);

  // The optical fiber connectors: three boxes
  // (X and Z are inverted on tray reference system)
  TGeoBBox *optConnPBT = new TGeoBBox(kSideCOptConnHalfZ,
		 kSideCOptConnPBTFrac*kSideCCoolManifHalfY,
				      kSideCOptConnHalfX);

  TGeoBBox *optConnSteel = new TGeoBBox(kSideCOptConnHalfZ,
		 kSideCOptConnSteelFrac*kSideCCoolManifHalfY,
					kSideCOptConnHalfX);

  TGeoBBox *optConnAl = new TGeoBBox(kSideCOptConnHalfZ,
		 kSideCOptConnAlFrac*kSideCCoolManifHalfY,
				     kSideCOptConnHalfX);

  // The optical fibers: a Xtru
  TGeoXtru *opticalFibs = new TGeoXtru(2);

  xprof[0] = 2*kSideCOptConnHalfZ;
  yprof[0] = 2*kSideCHalfThick + kSideCOptFibsTrans;
  xprof[1] = kSideCLength1;
  yprof[1] = yprof[0];
  xprof[2] = xprof[1] + kSideCLength2*TMath::Cos(alpharot);
  yprof[2] = yprof[1] + kSideCLength2*TMath::Sin(alpharot);
  xprof[3] = xprof[2] - kSideCOptFibsHigh*TMath::Sin(alpharot);
  yprof[3] = yprof[2] + kSideCOptFibsHigh*TMath::Cos(alpharot);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      kSideCOptFibsHigh, xprof[4], yprof[4]);
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + kSideCOptFibsHigh;

  opticalFibs->DefinePolygon(6, xprof, yprof);
  opticalFibs->DefineSection(0,-kSideCOptFibsWide/2);
  opticalFibs->DefineSection(1, kSideCOptFibsWide/2);

  // The input cables: five boxes
  // (X and Z are inverted on tray reference system)
  TGeoBBox *inputCabsCu = new TGeoBBox(kSideCInputCablesLen/2,
		   kSideCInputCablesCu*kSideCInputCablesHigh/2,
				       kSideCInputCablesWide/2);

  TGeoBBox *inputCabsPlast = new TGeoBBox(kSideCInputCablesLen/2,
		   kSideCInputCablesPlast*kSideCInputCablesHigh/2,
					  kSideCInputCablesWide/2);

  TGeoBBox *inputCabsAl = new TGeoBBox(kSideCInputCablesLen/2,
		   kSideCInputCablesAl*kSideCInputCablesHigh/2,
				       kSideCInputCablesWide/2);

  TGeoBBox *inputCabsKapton = new TGeoBBox(kSideCInputCablesLen/2,
		   kSideCInputCablesKapton*kSideCInputCablesHigh/2,
					   kSideCInputCablesWide/2);

  TGeoBBox *inputCabsPOLYAX = new TGeoBBox(kSideCInputCablesLen/2,
		   kSideCInputCablesPOLYAX*kSideCInputCablesHigh/2,
					   kSideCInputCablesWide/2);

  // The output cables: five Xtru
  TGeoXtru *outputCabsCu = new TGeoXtru(2);

  height = kSideCOutputCablesCu*kSideCOutputCablesHigh;

  xprof[0] = coolTubesAir->GetX(5) + kSideCOutputCablesX0;
  yprof[0] = coolTubesAir->GetY(5);
  xprof[1] = coolTubesAir->GetX(4);
  yprof[1] = coolTubesAir->GetY(4);
  xprof[2] = coolTubesAir->GetX(3);
  yprof[2] = coolTubesAir->GetY(3);
  xprof[3] = xprof[2] - height*TMath::Sin(alpharot);
  yprof[3] = yprof[2] + height*TMath::Cos(alpharot);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      height, xprof[4], yprof[4]);
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + height;

  outputCabsCu->DefinePolygon(6, xprof, yprof);
  outputCabsCu->DefineSection(0,-kSideCOutputCablesWide/2);
  outputCabsCu->DefineSection(1, kSideCOutputCablesWide/2);

  TGeoXtru *outputCabsPlast = new TGeoXtru(2);

  height = kSideCOutputCablesPlast*kSideCOutputCablesHigh;

  xprof[0] = outputCabsCu->GetX(5);
  yprof[0] = outputCabsCu->GetY(5);
  xprof[1] = outputCabsCu->GetX(4);
  yprof[1] = outputCabsCu->GetY(4);
  xprof[2] = outputCabsCu->GetX(3);
  yprof[2] = outputCabsCu->GetY(3);
  xprof[3] = xprof[2] - height*TMath::Sin(alpharot);
  yprof[3] = yprof[2] + height*TMath::Cos(alpharot);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      height, xprof[4], yprof[4]);
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + height;

  outputCabsPlast->DefinePolygon(6, xprof, yprof);
  outputCabsPlast->DefineSection(0,-kSideCOutputCablesWide/2);
  outputCabsPlast->DefineSection(1, kSideCOutputCablesWide/2);

  TGeoXtru *outputCabsAl = new TGeoXtru(2);

  height = kSideCOutputCablesAl*kSideCOutputCablesHigh;

  xprof[0] = outputCabsPlast->GetX(5);
  yprof[0] = outputCabsPlast->GetY(5);
  xprof[1] = outputCabsPlast->GetX(4);
  yprof[1] = outputCabsPlast->GetY(4);
  xprof[2] = outputCabsPlast->GetX(3);
  yprof[2] = outputCabsPlast->GetY(3);
  xprof[3] = xprof[2] - height*TMath::Sin(alpharot);
  yprof[3] = yprof[2] + height*TMath::Cos(alpharot);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      height, xprof[4], yprof[4]);
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + height;

  outputCabsAl->DefinePolygon(6, xprof, yprof);
  outputCabsAl->DefineSection(0,-kSideCOutputCablesWide/2);
  outputCabsAl->DefineSection(1, kSideCOutputCablesWide/2);

  TGeoXtru *outputCabsKapton = new TGeoXtru(2);

  height = kSideCOutputCablesKapton*kSideCOutputCablesHigh;

  xprof[0] = outputCabsAl->GetX(5);
  yprof[0] = outputCabsAl->GetY(5);
  xprof[1] = outputCabsAl->GetX(4);
  yprof[1] = outputCabsAl->GetY(4);
  xprof[2] = outputCabsAl->GetX(3);
  yprof[2] = outputCabsAl->GetY(3);
  xprof[3] = xprof[2] - height*TMath::Sin(alpharot);
  yprof[3] = yprof[2] + height*TMath::Cos(alpharot);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      height, xprof[4], yprof[4]);
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + height;

  outputCabsKapton->DefinePolygon(6, xprof, yprof);
  outputCabsKapton->DefineSection(0,-kSideCOutputCablesWide/2);
  outputCabsKapton->DefineSection(1, kSideCOutputCablesWide/2);

  TGeoXtru *outputCabsPOLYAX = new TGeoXtru(2);

  height = kSideCOutputCablesPOLYAX*kSideCOutputCablesHigh;

  xprof[0] = outputCabsKapton->GetX(5);
  yprof[0] = outputCabsKapton->GetY(5);
  xprof[1] = outputCabsKapton->GetX(4);
  yprof[1] = outputCabsKapton->GetY(4);
  xprof[2] = outputCabsKapton->GetX(3);
  yprof[2] = outputCabsKapton->GetY(3);
  xprof[3] = xprof[2] - height*TMath::Sin(alpharot);
  yprof[3] = yprof[2] + height*TMath::Cos(alpharot);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      height, xprof[4], yprof[4]);
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + height;

  outputCabsPOLYAX->DefinePolygon(6, xprof, yprof);
  outputCabsPOLYAX->DefineSection(0,-kSideCOutputCablesWide/2);
  outputCabsPOLYAX->DefineSection(1, kSideCOutputCablesWide/2);

  // The PCB boards: five boxes
  // (X and Z are inverted on tray reference system)
  TGeoBBox *pcbBoardsCu = new TGeoBBox(kSideCPCBBoardsHalfZ,
		     kSideCPCBBoardsCu*kSideCPCBBoardsHalfY,
				       kSideCPCBBoardsHalfX);

  TGeoBBox *pcbBoardsEpoxy = new TGeoBBox(kSideCPCBBoardsHalfZ,
		     kSideCPCBBoardsEpoxy*kSideCPCBBoardsHalfY,
					  kSideCPCBBoardsHalfX);

  TGeoBBox *pcbBoardsPlast = new TGeoBBox(kSideCPCBBoardsHalfZ,
		     kSideCPCBBoardsPlast*kSideCPCBBoardsHalfY,
					  kSideCPCBBoardsHalfX);

  TGeoBBox *pcbBoardsSteel = new TGeoBBox(kSideCPCBBoardsHalfZ,
		     kSideCPCBBoardsSteel*kSideCPCBBoardsHalfY,
					  kSideCPCBBoardsHalfX);

  TGeoBBox *pcbBoardsPPS = new TGeoBBox(kSideCPCBBoardsHalfZ,
		     kSideCPCBBoardsPPS*kSideCPCBBoardsHalfY,
					kSideCPCBBoardsHalfX);


  // We have all shapes: now create the real volumes
  TGeoMedium *medPOM    = mgr->GetMedium("ITS_POLYOXYMETHYLENE$");
  TGeoMedium *medSteel  = mgr->GetMedium("ITS_INOX$");
  TGeoMedium *medWater  = mgr->GetMedium("ITS_WATER$");
  TGeoMedium *medAl     = mgr->GetMedium("ITS_ALUMINUM$");
  TGeoMedium *medCu     = mgr->GetMedium("ITS_COPPER$");
  TGeoMedium *medPUR    = mgr->GetMedium("ITS_POLYURETHANE$");
  TGeoMedium *medPOLYAX = mgr->GetMedium("ITS_POLYAX$");
  TGeoMedium *medKapton = mgr->GetMedium("ITS_SDDKAPTON (POLYCH2)$");
  TGeoMedium *medAir    = mgr->GetMedium("ITS_AIR$");
  TGeoMedium *medPBT    = mgr->GetMedium("ITS_PBT$");
  TGeoMedium *medOptFib = mgr->GetMedium("ITS_SDD OPTICFIB$");
  TGeoMedium *medPPS    = mgr->GetMedium("ITS_PPS$");
  TGeoMedium *medEpoxy  = mgr->GetMedium("ITS_EPOXY$");

  TGeoVolume *pomCoolManif = new TGeoVolume("ITSsuppSDDSideCCoolManifPOM",
					    coolManifPOM, medPOM);

  pomCoolManif->SetVisibility(kTRUE);
  pomCoolManif->SetLineColor(kRed); // Red
  pomCoolManif->SetLineWidth(1);
  pomCoolManif->SetFillColor(pomCoolManif->GetLineColor());
  pomCoolManif->SetFillStyle(4000); // 0% transparent

  TGeoVolume *steelCoolManif = new TGeoVolume("ITSsuppSDDSideCCoolManifSteel",
					      coolManifSteel, medSteel);

  steelCoolManif->SetVisibility(kTRUE);
  steelCoolManif->SetLineColor(kBlue); // Blue
  steelCoolManif->SetLineWidth(1);
  steelCoolManif->SetFillColor(steelCoolManif->GetLineColor());
  steelCoolManif->SetFillStyle(4000); // 0% transparent

  TGeoVolume *waterCoolManif = new TGeoVolume("ITSsuppSDDSideCCoolManifWater",
					      coolManifWater, medWater);

  waterCoolManif->SetVisibility(kTRUE);
  waterCoolManif->SetLineColor(33); // Light Blue
  waterCoolManif->SetLineWidth(1);
  waterCoolManif->SetFillColor(waterCoolManif->GetLineColor());
  waterCoolManif->SetFillStyle(4000); // 0% transparent

  TGeoVolume *alCoolManif = new TGeoVolume("ITSsuppSDDSideCCoolManifAl",
					   coolManifAl, medAl);

  alCoolManif->SetVisibility(kTRUE);
  alCoolManif->SetLineColor(6); // Purple
  alCoolManif->SetLineWidth(1);
  alCoolManif->SetFillColor(alCoolManif->GetLineColor());
  alCoolManif->SetFillStyle(4000); // 0% transparent

  TGeoVolume *purCoolTubes = new TGeoVolume("ITSsuppSDDSideCCoolTubesPUR",
					    coolTubesPUR, medPUR);

  purCoolTubes->SetVisibility(kTRUE);
  purCoolTubes->SetLineColor(kRed); // Red
  purCoolTubes->SetLineWidth(1);
  purCoolTubes->SetFillColor(purCoolTubes->GetLineColor());
  purCoolTubes->SetFillStyle(4000); // 0% transparent

  TGeoVolume *waterCoolTubes = new TGeoVolume("ITSsuppSDDSideCCoolTubesWater",
					      coolTubesWater, medWater);

  waterCoolTubes->SetVisibility(kTRUE);
  waterCoolTubes->SetLineColor(33); // Light Blue
  waterCoolTubes->SetLineWidth(1);
  waterCoolTubes->SetFillColor(waterCoolTubes->GetLineColor());
  waterCoolTubes->SetFillStyle(4000); // 0% transparent

  TGeoVolume *airCoolTubes = new TGeoVolume("ITSsuppSDDSideCCoolTubesAir",
					    coolTubesAir, medAir);

  airCoolTubes->SetVisibility(kTRUE);
  airCoolTubes->SetLineColor(41);
  airCoolTubes->SetLineWidth(1);
  airCoolTubes->SetFillColor(airCoolTubes->GetLineColor());
  airCoolTubes->SetFillStyle(4000); // 0% transparent

  TGeoVolume *pbtOptConn = new TGeoVolume("ITSsuppSDDSideCOptConnPBT",
					  optConnPBT, medPBT);

  pbtOptConn->SetVisibility(kTRUE);
  pbtOptConn->SetLineColor(kRed); // Red
  pbtOptConn->SetLineWidth(1);
  pbtOptConn->SetFillColor(pbtOptConn->GetLineColor());
  pbtOptConn->SetFillStyle(4000); // 0% transparent

  TGeoVolume *steelOptConn = new TGeoVolume("ITSsuppSDDSideCOptConnSteel",
					    optConnSteel, medSteel);

  steelOptConn->SetVisibility(kTRUE);
  steelOptConn->SetLineColor(kBlue); // Blue
  steelOptConn->SetLineWidth(1);
  steelOptConn->SetFillColor(steelOptConn->GetLineColor());
  steelOptConn->SetFillStyle(4000); // 0% transparent

  TGeoVolume *alOptConn = new TGeoVolume("ITSsuppSDDSideCOptConnAl",
					 optConnAl, medAl);

  alOptConn->SetVisibility(kTRUE);
  alOptConn->SetLineColor(6); // Purple
  alOptConn->SetLineWidth(1);
  alOptConn->SetFillColor(alOptConn->GetLineColor());
  alOptConn->SetFillStyle(4000); // 0% transparent

  TGeoVolume *optFibs = new TGeoVolume("ITSsuppSDDSideCOptFibs",
				       opticalFibs, medOptFib);

  optFibs->SetVisibility(kTRUE);
  optFibs->SetLineColor(kOrange+2); // Orange
  optFibs->SetLineWidth(1);
  optFibs->SetFillColor(optFibs->GetLineColor());
  optFibs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *cuInputCabs = new TGeoVolume("ITSsuppSDDSideCInputCabsCu",
					   inputCabsCu, medCu);

  cuInputCabs->SetVisibility(kTRUE);
  cuInputCabs->SetLineColor(kBlack); // Black
  cuInputCabs->SetLineWidth(1);
  cuInputCabs->SetFillColor(cuInputCabs->GetLineColor());
  cuInputCabs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *plastInputCabs = new TGeoVolume("ITSsuppSDDSideCInputCabsPlast",
					      inputCabsPlast, medPUR);

  plastInputCabs->SetVisibility(kTRUE);
  plastInputCabs->SetLineColor(kRed); // Red
  plastInputCabs->SetLineWidth(1);
  plastInputCabs->SetFillColor(plastInputCabs->GetLineColor());
  plastInputCabs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *alInputCabs = new TGeoVolume("ITSsuppSDDSideCInputCabsAl",
					   inputCabsAl, medAl);

  alInputCabs->SetVisibility(kTRUE);
  alInputCabs->SetLineColor(6); // Purple
  alInputCabs->SetLineWidth(1);
  alInputCabs->SetFillColor(alInputCabs->GetLineColor());
  alInputCabs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *kaptonInputCabs = new TGeoVolume("ITSsuppSDDSideCInputCabsKapton",
					       inputCabsKapton, medKapton);

  kaptonInputCabs->SetVisibility(kTRUE);
  kaptonInputCabs->SetLineColor(14); // 
  kaptonInputCabs->SetLineWidth(1);
  kaptonInputCabs->SetFillColor(kaptonInputCabs->GetLineColor());
  kaptonInputCabs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *polyaxInputCabs = new TGeoVolume("ITSsuppSDDSideCInputCabsPOLYAX",
					       inputCabsPOLYAX, medPOLYAX);

  polyaxInputCabs->SetVisibility(kTRUE);
  polyaxInputCabs->SetLineColor(34); // 
  polyaxInputCabs->SetLineWidth(1);
  polyaxInputCabs->SetFillColor(polyaxInputCabs->GetLineColor());
  polyaxInputCabs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *cuOutputCabs = new TGeoVolume("ITSsuppSDDSideCOutputCabsCu",
					    outputCabsCu, medCu);

  cuOutputCabs->SetVisibility(kTRUE);
  cuOutputCabs->SetLineColor(kBlack); // Black
  cuOutputCabs->SetLineWidth(1);
  cuOutputCabs->SetFillColor(cuOutputCabs->GetLineColor());
  cuOutputCabs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *plastOutputCabs = new TGeoVolume("ITSsuppSDDSideCOutputCabsPlast",
					       outputCabsPlast, medPUR);

  plastOutputCabs->SetVisibility(kTRUE);
  plastOutputCabs->SetLineColor(kRed); // Red
  plastOutputCabs->SetLineWidth(1);
  plastOutputCabs->SetFillColor(plastOutputCabs->GetLineColor());
  plastOutputCabs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *alOutputCabs = new TGeoVolume("ITSsuppSDDSideCOutputCabsAl",
					    outputCabsAl, medAl);

  alOutputCabs->SetVisibility(kTRUE);
  alOutputCabs->SetLineColor(6); // Purple
  alOutputCabs->SetLineWidth(1);
  alOutputCabs->SetFillColor(alOutputCabs->GetLineColor());
  alOutputCabs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *kaptonOutputCabs = new TGeoVolume("ITSsuppSDDSideCOutputCabsKapton",
						outputCabsKapton, medKapton);

  kaptonOutputCabs->SetVisibility(kTRUE);
  kaptonOutputCabs->SetLineColor(14); // 
  kaptonOutputCabs->SetLineWidth(1);
  kaptonOutputCabs->SetFillColor(kaptonOutputCabs->GetLineColor());
  kaptonOutputCabs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *polyaxOutputCabs = new TGeoVolume("ITSsuppSDDSideCOutputCabsPOLYAX",
						outputCabsPOLYAX, medPOLYAX);

  polyaxOutputCabs->SetVisibility(kTRUE);
  polyaxOutputCabs->SetLineColor(34); // 
  polyaxOutputCabs->SetLineWidth(1);
  polyaxOutputCabs->SetFillColor(polyaxOutputCabs->GetLineColor());
  polyaxOutputCabs->SetFillStyle(4000); // 0% transparent

  TGeoVolume *cuPCBBoards = new TGeoVolume("ITSsuppSDDSideCPCBBoardsCu",
					   pcbBoardsCu, medCu);

  cuPCBBoards->SetVisibility(kTRUE);
  cuPCBBoards->SetLineColor(kBlack); // Black
  cuPCBBoards->SetLineWidth(1);
  cuPCBBoards->SetFillColor(cuPCBBoards->GetLineColor());
  cuPCBBoards->SetFillStyle(4000); // 0% transparent

  TGeoVolume *epoxyPCBBoards = new TGeoVolume("ITSsuppSDDSideCPCBBoardsEpoxy",
					      pcbBoardsEpoxy, medEpoxy);

  epoxyPCBBoards->SetVisibility(kTRUE);
  epoxyPCBBoards->SetLineColor(22); //
  epoxyPCBBoards->SetLineWidth(1);
  epoxyPCBBoards->SetFillColor(epoxyPCBBoards->GetLineColor());
  epoxyPCBBoards->SetFillStyle(4000); // 0% transparent

  TGeoVolume *plastPCBBoards = new TGeoVolume("ITSsuppSDDSideCPCBBoardsPlast",
					      pcbBoardsPlast, medPUR);

  plastPCBBoards->SetVisibility(kTRUE);
  plastPCBBoards->SetLineColor(kRed); // Red
  plastPCBBoards->SetLineWidth(1);
  plastPCBBoards->SetFillColor(plastPCBBoards->GetLineColor());
  plastPCBBoards->SetFillStyle(4000); // 0% transparent

  TGeoVolume *steelPCBBoards = new TGeoVolume("ITSsuppSDDSideCPCBBoardsSteel",
					      pcbBoardsSteel, medSteel);

  steelPCBBoards->SetVisibility(kTRUE);
  steelPCBBoards->SetLineColor(kBlue); // Blue
  steelPCBBoards->SetLineWidth(1);
  steelPCBBoards->SetFillColor(steelPCBBoards->GetLineColor());
  steelPCBBoards->SetFillStyle(4000); // 0% transparent

  TGeoVolume *ppsPCBBoards = new TGeoVolume("ITSsuppSDDSideCPCBBoardsPPS",
					    pcbBoardsPPS, medPPS);

  ppsPCBBoards->SetVisibility(kTRUE);
  ppsPCBBoards->SetLineColor(kGreen); // Green
  ppsPCBBoards->SetLineWidth(1);
  ppsPCBBoards->SetFillColor(ppsPCBBoards->GetLineColor());
  ppsPCBBoards->SetFillStyle(4000); // 0% transparent


  // Now fill the tray
  xloc = coolManifPOM->GetDX();
  yloc = 2*kSideCHalfThick + coolManifPOM->GetDY();
  trayStructure->AddNode(pomCoolManif, 1,
			 new TGeoTranslation( xloc, yloc, 0) );

  yloc += coolManifPOM->GetDY() + coolManifSteel->GetDY();
  trayStructure->AddNode(steelCoolManif, 1,
			 new TGeoTranslation( xloc, yloc, 0) );

  yloc += coolManifSteel->GetDY() + coolManifWater->GetDY();
  trayStructure->AddNode(waterCoolManif, 1,
			 new TGeoTranslation( xloc, yloc, 0) );

  yloc += coolManifWater->GetDY() + coolManifAl->GetDY();
  trayStructure->AddNode(alCoolManif, 1,
			 new TGeoTranslation( xloc, yloc, 0) );

  xloc = inputCabsCu->GetDX();
  yloc += coolManifWater->GetDY() + inputCabsCu->GetDY()
        + kSideCInputCablesTrans;
  trayStructure->AddNode(cuInputCabs, 1,
			 new TGeoTranslation( xloc, yloc, 0) );

  yloc += inputCabsCu->GetDY() + inputCabsPlast->GetDY();
  trayStructure->AddNode(plastInputCabs, 1,
			 new TGeoTranslation( xloc, yloc, 0) );

  yloc += inputCabsPlast->GetDY() + inputCabsAl->GetDY();
  trayStructure->AddNode(alInputCabs, 1,
			 new TGeoTranslation( xloc, yloc, 0) );

  yloc += inputCabsAl->GetDY() + inputCabsKapton->GetDY();
  trayStructure->AddNode(kaptonInputCabs, 1,
			 new TGeoTranslation( xloc, yloc, 0) );

  yloc += inputCabsKapton->GetDY() + inputCabsPOLYAX->GetDY();
  trayStructure->AddNode(polyaxInputCabs, 1,
			 new TGeoTranslation( xloc, yloc, 0) );

  trayStructure->AddNode(purCoolTubes  , 1, 0);
  trayStructure->AddNode(waterCoolTubes, 1, 0);
  trayStructure->AddNode(airCoolTubes  , 1, 0);

  xloc = optConnPBT->GetDX();
  yloc = 2*kSideCHalfThick + optConnPBT->GetDY();
  zloc = coolManifPOM->GetDZ() + optConnPBT->GetDZ();
  trayStructure->AddNode(pbtOptConn, 1,
			 new TGeoTranslation( xloc, yloc, zloc) );
  trayStructure->AddNode(pbtOptConn, 2,
			 new TGeoTranslation( xloc, yloc,-zloc) );

  yloc += optConnPBT->GetDY() + optConnSteel->GetDY();
  trayStructure->AddNode(steelOptConn, 1,
			 new TGeoTranslation( xloc, yloc, zloc) );
  trayStructure->AddNode(steelOptConn, 2,
			 new TGeoTranslation( xloc, yloc,-zloc) );

  yloc += optConnSteel->GetDY() + optConnAl->GetDY();
  trayStructure->AddNode(alOptConn, 1,
			 new TGeoTranslation( xloc, yloc, zloc) );
  trayStructure->AddNode(alOptConn, 2,
			 new TGeoTranslation( xloc, yloc,-zloc) );

  trayStructure->AddNode(optFibs, 1,
			 new TGeoTranslation( 0, 0, zloc) );
  trayStructure->AddNode(optFibs, 2,
			 new TGeoTranslation( 0, 0,-zloc) );

  trayStructure->AddNode(cuOutputCabs    , 1, 0);
  trayStructure->AddNode(plastOutputCabs , 1, 0);
  trayStructure->AddNode(alOutputCabs    , 1, 0);
  trayStructure->AddNode(kaptonOutputCabs, 1, 0);
  trayStructure->AddNode(polyaxOutputCabs, 1, 0);

  xloc = kXShiftBarCool + kBarCoolRmax + pcbBoardsCu->GetDX();
  yloc = outputCabsPOLYAX->GetY(5) + pcbBoardsCu->GetDY();
  trayStructure->AddNode(cuPCBBoards, 1,
			 new TGeoTranslation( xloc, yloc , 0) );

  yloc += pcbBoardsCu->GetDY() + pcbBoardsEpoxy->GetDY();
  trayStructure->AddNode(epoxyPCBBoards, 1,
			 new TGeoTranslation( xloc, yloc , 0) );

  yloc += pcbBoardsEpoxy->GetDY() + pcbBoardsPlast->GetDY();
  trayStructure->AddNode(plastPCBBoards, 1,
			 new TGeoTranslation( xloc, yloc , 0) );

  yloc += pcbBoardsPlast->GetDY() + pcbBoardsSteel->GetDY();
  trayStructure->AddNode(steelPCBBoards, 1,
			 new TGeoTranslation( xloc, yloc , 0) );

  yloc += pcbBoardsSteel->GetDY() + pcbBoardsPPS->GetDY();
  trayStructure->AddNode(ppsPCBBoards, 1,
			 new TGeoTranslation( xloc, yloc , 0) );


  // Finally put everything in the mother volume
  alphafold = kSideCFoldAngle;

  for (Int_t jt = 0; jt < kNumTraySideC; jt++) {
    alpharot = kTraySideCAlphaRot[jt];
    xloc = kTraySideCRPos*SinD(alpharot);
    yloc = kTraySideCRPos*CosD(alpharot);
    moth->AddNode(trayStructure,jt+1,
		       new TGeoCombiTrans(-xloc, yloc, kTraySideCZPos,
		       new TGeoRotation("",-90.+alpharot,-90.,90.+alphafold)));
  }


  return;
}


//______________________________________________________________________
void AliITSv11GeometrySupport::SSDCableTraysSideA(TGeoVolume *moth,
					    const TGeoManager *mgr){
//
// Creates the SSD cable trays which are outside the ITS support cones
// but still inside the TPC on Side A
// (part of this code is taken or anyway inspired to ServicesCableSupport
// method of AliITSv11GeometrySupport.cxx,v 1.9 2007/06/06)
//
// Input:
//         moth : the TGeoVolume owing the volume structure
//         mgr  : the GeoManager (default gGeoManager)
// Output:
//
// Created:         ???       Bjorn S. Nilsen
// Updated:      30 Dec 2009  Mario Sitta
//
// Technical data are taken from AutoCAD drawings, L.Simonetti technical
// drawings and other (oral) information given by F.Tosello and
// Ton van den Brink
// Cables and cooling tubes are approximated with proper materials and
// rectangular cross sections, always preserving the total material budget.
//

  // Dimensions and positions of the A-Side Cable Trays
  // (parts of 0872/G/D)
  const Double_t kTrayARTrans            =  408.35 *fgkmm;
  const Double_t kTrayAZTrans            = 1011.00 *fgkmm;
  const Double_t kForwardSideYTrans      =   12.00 *fgkmm;//!!!TO BE CHECKED!!!
  const Double_t kCoversYTrans           =    2.00 *fgkmm;
  const Double_t kTrayAZRot              = (180-169.5);// Degrees
  const Double_t kTrayAFirstRotAng       =   22.00;    // Degrees
  const Double_t kTrayASecondRotAng      =   15.00;    // Degrees

  const Double_t kTrayTotalHeight        =   52.00 *fgkmm;
  const Double_t kTrayHeighToBend        =   32.00 *fgkmm;
  const Double_t kTrayWidth              =  130.00 *fgkmm;
  const Double_t kTrayThick              =    2.00 *fgkmm;

  const Double_t kTrayBendAngle          =   22.00 *TMath::DegToRad();

  const Double_t kForwardTrayTotalLen    =  853.00 *fgkmm;
  const Double_t kForwardTrayFirstLen    =  350.00 *fgkmm;
  const Double_t kForwardTrayFirstHeight =   47.00 *fgkmm;
  const Double_t kForwardCoverLen        =  420.00 *fgkmm;

  const Double_t kForwardSideLength      = kForwardTrayFirstLen;//!!!TO BE CHECKED!!!
  const Double_t kForwardSideHeight      =   90.00 *fgkmm;//!!!TO BE CHECKED!!!
  const Double_t kForwardSideThick       =    1.00 *fgkmm;//!!!TO BE CHECKED!!!
  const Double_t kForwardCoverHeight     =   10.00 *fgkmm;//!!!TO BE CHECKED!!!

  const Double_t kExternalTrayTotalLen   = 1200.00 *fgkmm;
  const Double_t kExternalCoverLen       = kExternalTrayTotalLen;
  const Double_t kExternalCoverThick     =    5.00 *fgkmm;

  const Int_t    kForwardTrayNpoints     =   16;

  const Double_t kServicesWidth          =  100.00 *fgkmm;
  const Double_t kCopperHeight           =   11.20 *fgkmm;// 1120 mm^2
  const Double_t kCablePlasticHeight     =   11.50 *fgkmm;// 1150 mm^2
  const Double_t kCoolingWaterHeight     =    2.65 *fgkmm;//  265 mm^2
  const Double_t kPoliUrethaneHeight     =    4.62 *fgkmm;//  462 mm^2


  // Local variables
  Double_t xprof[kForwardTrayNpoints], yprof[kForwardTrayNpoints];
  Double_t xloc, yloc, zloc, alpharot, totalhi;


  // The two tray components as assemblies
  TGeoVolumeAssembly *cableTrayAForw =
    new TGeoVolumeAssembly("ITSsupportSSDTrayAForw");
  TGeoVolumeAssembly *cableTrayAExt =
    new TGeoVolumeAssembly("ITSsupportSSDTrayAExt");
  

  // First create all needed shapes

  // The first part of the forward tray (part of 0872/G/D/07): a Xtru
  TGeoXtru *forwTrayPart1 = new TGeoXtru(2);

  xprof[3] = kTrayWidth/2;
  yprof[3] = kForwardTrayFirstHeight;
  xprof[2] = xprof[3] - kTrayThick;
  yprof[2] = yprof[3];
  xprof[4] = xprof[3];
  yprof[4] = kTrayTotalHeight - kTrayHeighToBend;
  xprof[5] = xprof[4] - yprof[4]*TMath::Tan(kTrayBendAngle);
  yprof[5] = 0;

  InsidePoint( xprof[3], yprof[3], xprof[4], yprof[4], xprof[5], yprof[5],
	      -kTrayThick, xprof[1], yprof[1]);

  xprof[6] = -xprof[5];
  yprof[6] =  yprof[5];

  InsidePoint( xprof[4], yprof[4], xprof[5], yprof[5], xprof[6], yprof[6],
	      -kTrayThick, xprof[0], yprof[0]);

  // We did the right side, now reflex on the left side
  for (Int_t jp = 0; jp < 6; jp++) {
    xprof[6+jp] = -xprof[5-jp];
    yprof[6+jp] =  yprof[5-jp];
  }

  // And now the actual Xtru
  forwTrayPart1->DefinePolygon(12, xprof, yprof);
  forwTrayPart1->DefineSection(0, 0);
  forwTrayPart1->DefineSection(1, kForwardTrayFirstLen);

  // The second part of the forward tray (part of 0872/G/D/07): a Xtru
  TGeoXtru *forwTrayPart2 =
    CreateSDDSSDTraysSideA(kForwardTrayTotalLen - kForwardTrayFirstLen,
			   kTrayTotalHeight);

  // The external tray (as 0872/G/D/03): a Xtru with same profile
  TGeoXtru *externalTray = CreateSDDSSDTraysSideA(kExternalTrayTotalLen,
						  kTrayTotalHeight);

  // The side wall of the forward tray: a BBox
  TGeoBBox *forwSide = new TGeoBBox(kForwardSideThick/2,
				    kForwardSideHeight/2,
				    kForwardSideLength/2);

  // The side cover over the walls: a Xtru
  TGeoXtru *forwSideCover = new TGeoXtru(2);
  forwSideCover->SetName("ITSsuppSSDForwCover");

  xprof[0] = kTrayWidth/2 + 2*kForwardSideThick;
  yprof[0] = kForwardCoverHeight;
  xprof[1] = xprof[0];
  yprof[1] = 0;
  xprof[2] = xprof[1] - kForwardSideThick;
  yprof[2] = yprof[1];
  xprof[3] = xprof[2];
  yprof[3] = yprof[0] - kForwardSideThick;

  // We did the right side, now reflex on the left side
  for (Int_t jp = 0; jp < 4; jp++) {
    xprof[4+jp] = -xprof[3-jp];
    yprof[4+jp] =  yprof[3-jp];
  }

  forwSideCover->DefinePolygon(8, xprof, yprof);
  forwSideCover->DefineSection(0, 0);
  forwSideCover->DefineSection(1, kForwardSideLength);

  // The forward and external covers: two Composite Shape's
  TGeoCompositeShape *forwardCover = CreateTrayAForwardCover(kForwardCoverLen);

  TGeoCompositeShape *externCover = CreateTrayAExternalCover(kExternalCoverLen);

  // The cable copper inside the forward tray: a BBox
  TGeoBBox *forwCopper = new TGeoBBox(kServicesWidth/2,
				      kCopperHeight/2,
				      kForwardTrayTotalLen/2);

  // The cable copper inside the forward tray: a Xtru
  TGeoXtru *extCopper = new TGeoXtru(2);
  extCopper->SetName("ITSsuppSSDExtTrayCopper");

  totalhi = kTrayTotalHeight + kExternalCoverThick - kCoversYTrans
	  - kTrayThick;

  xprof[0] = -totalhi*TanD(kTrayAZRot);
  yprof[0] = kTrayThick;
  xprof[1] = kExternalTrayTotalLen;
  yprof[1] = yprof[0];
  xprof[2] = xprof[1];
  yprof[2] = yprof[1] + kCopperHeight;
  totalhi -= kCopperHeight;
  xprof[3] = -totalhi*TanD(kTrayAZRot);
  yprof[3] = yprof[2];

  extCopper->DefinePolygon(4, xprof, yprof);
  extCopper->DefineSection(0, 0);
  extCopper->DefineSection(1, kServicesWidth);

  // The cable plastic inside the forward tray: a BBox
  TGeoBBox *forwPlastic = new TGeoBBox(kServicesWidth/2,
				       kCablePlasticHeight/2,
				       kForwardTrayTotalLen/2);

  // The cable plastic inside the forward tray: a Xtru
  TGeoXtru *extPlastic = new TGeoXtru(2);
  extPlastic->SetName("ITSsuppSSDExtTrayPlastic");

  totalhi = kTrayTotalHeight + kExternalCoverThick - kCoversYTrans
	  - kTrayThick - kCopperHeight;

  xprof[0] = -totalhi*TanD(kTrayAZRot);
  yprof[0] = kTrayThick;
  xprof[1] = kExternalTrayTotalLen;
  yprof[1] = yprof[0];
  xprof[2] = xprof[1];
  yprof[2] = yprof[1] + kCablePlasticHeight;
  totalhi -= kCablePlasticHeight;
  xprof[3] = -totalhi*TanD(kTrayAZRot);
  yprof[3] = yprof[2];

  extPlastic->DefinePolygon(4, xprof, yprof);
  extPlastic->DefineSection(0, 0);
  extPlastic->DefineSection(1, kServicesWidth);

  // The cooling water inside the forward tray: a BBox
  TGeoBBox *forwWater = new TGeoBBox(kServicesWidth/2,
				     kCoolingWaterHeight/2,
				     kForwardTrayTotalLen/2);

  // The cooling water inside the forward tray: a Xtru
  TGeoXtru *extWater = new TGeoXtru(2);
  extWater->SetName("ITSsuppSSDExtTrayWater");

  totalhi = kTrayTotalHeight + kExternalCoverThick - kCoversYTrans
	  - kTrayThick - kCopperHeight - kCablePlasticHeight;

  xprof[0] = -totalhi*TanD(kTrayAZRot);
  yprof[0] = kTrayThick;
  xprof[1] = kExternalTrayTotalLen;
  yprof[1] = yprof[0];
  xprof[2] = xprof[1];
  yprof[2] = yprof[1] + kCoolingWaterHeight;
  totalhi -= kCoolingWaterHeight;
  xprof[3] = -totalhi*TanD(kTrayAZRot);
  yprof[3] = yprof[2];

  extWater->DefinePolygon(4, xprof, yprof);
  extWater->DefineSection(0, 0);
  extWater->DefineSection(1, kServicesWidth);

  // The polyurethane inside the forward tray: a BBox
  TGeoBBox *forwPUR = new TGeoBBox(kServicesWidth/2,
				   kPoliUrethaneHeight/2,
				   kForwardTrayTotalLen/2);

  // The poliurethane inside the forward tray: a Xtru
  TGeoXtru *extPUR = new TGeoXtru(2);
  extPUR->SetName("ITSsuppSSDExtTrayPUR");

  totalhi = kTrayTotalHeight + kExternalCoverThick - kCoversYTrans
	  - kTrayThick - kCopperHeight - kCablePlasticHeight
	  - kCoolingWaterHeight;

  xprof[0] = -totalhi*TanD(kTrayAZRot);
  yprof[0] = kTrayThick;
  xprof[1] = kExternalTrayTotalLen;
  yprof[1] = yprof[0];
  xprof[2] = xprof[1];
  yprof[2] = yprof[1] + kPoliUrethaneHeight;
  totalhi -= kPoliUrethaneHeight;
  xprof[3] = -totalhi*TanD(kTrayAZRot);
  yprof[3] = yprof[2];

  extPUR->DefinePolygon(4, xprof, yprof);
  extPUR->DefineSection(0, 0);
  extPUR->DefineSection(1, kServicesWidth);


  // We have all shapes: now create the real volumes
  TGeoMedium *medAl    = mgr->GetMedium("ITS_ALUMINUM$");
  TGeoMedium *medAntic = mgr->GetMedium("ITS_ANTICORODAL$");
  TGeoMedium *medCu    = mgr->GetMedium("ITS_COPPER$");
  TGeoMedium *medFEP   = mgr->GetMedium("ITS_SSD FEP$");
  TGeoMedium *medH2O   = mgr->GetMedium("ITS_WATER$");
  TGeoMedium *medPUR   = mgr->GetMedium("ITS_POLYURETHANE$");

  TGeoVolume *forwTrayFirst = new TGeoVolume("ITSsuppSSDSideAForwTrayFirst",
					     forwTrayPart1, medAl);

  forwTrayFirst->SetVisibility(kTRUE);
  forwTrayFirst->SetLineColor(6); // Purple
  forwTrayFirst->SetLineWidth(1);
  forwTrayFirst->SetFillColor(forwTrayFirst->GetLineColor());
  forwTrayFirst->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwTraySecond = new TGeoVolume("ITSsuppSSDSideAForwTraySecond",
					      forwTrayPart2, medAl);

  forwTraySecond->SetVisibility(kTRUE);
  forwTraySecond->SetLineColor(6); // Purple
  forwTraySecond->SetLineWidth(1);
  forwTraySecond->SetFillColor(forwTraySecond->GetLineColor());
  forwTraySecond->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwTraySide = new TGeoVolume("ITSsuppSSDSideAForwTraySide",
					    forwSide, medAl);

  forwTraySide->SetVisibility(kTRUE);
  forwTraySide->SetLineColor(6); // Purple
  forwTraySide->SetLineWidth(1);
  forwTraySide->SetFillColor(forwTraySide->GetLineColor());
  forwTraySide->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwTraySideCover = new TGeoVolume("ITSsuppSSDSideAForwTraySideCover",
					    forwSideCover, medAl);

  forwTraySideCover->SetVisibility(kTRUE);
  forwTraySideCover->SetLineColor(6); // Purple
  forwTraySideCover->SetLineWidth(1);
  forwTraySideCover->SetFillColor(forwTraySideCover->GetLineColor());
  forwTraySideCover->SetFillStyle(4000); // 0% transparent

  TGeoVolume *externalTraySSD = new TGeoVolume("ITSsuppSSDSideAExternalTray",
					       externalTray, medAl);

  externalTraySSD->SetVisibility(kTRUE);
  externalTraySSD->SetLineColor(6); // Purple
  externalTraySSD->SetLineWidth(1);
  externalTraySSD->SetFillColor(externalTraySSD->GetLineColor());
  externalTraySSD->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwardTrayCover = new TGeoVolume("ITSsuppSSDSideAForwTrayCover",
						forwardCover, medAntic);

  forwardTrayCover->SetVisibility(kTRUE);
  forwardTrayCover->SetLineColor(kMagenta+1); // Purple
  forwardTrayCover->SetLineWidth(1);
  forwardTrayCover->SetFillColor(forwardTrayCover->GetLineColor());
  forwardTrayCover->SetFillStyle(4000); // 0% transparent

  TGeoVolume *externTrayCover = new TGeoVolume("ITSsuppSSDSideAExtTrayCover",
					       externCover, medAntic);

  externTrayCover->SetVisibility(kTRUE);
  externTrayCover->SetLineColor(kMagenta+1); // Purple
  externTrayCover->SetLineWidth(1);
  externTrayCover->SetFillColor(externTrayCover->GetLineColor());
  externTrayCover->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwCableCu = new TGeoVolume("ITSsuppSSDSideAForwCableCu",
					   forwCopper, medCu);

  forwCableCu->SetVisibility(kTRUE);
  forwCableCu->SetLineColor(kRed); // Red
  forwCableCu->SetLineWidth(1);
  forwCableCu->SetFillColor(forwCableCu->GetLineColor());
  forwCableCu->SetFillStyle(4000); // 0% transparent

  TGeoVolume *extCableCu = new TGeoVolume("ITSsuppSSDSideAExtCableCu",
					  extCopper, medCu);

  extCableCu->SetVisibility(kTRUE);
  extCableCu->SetLineColor(kRed); // Red
  extCableCu->SetLineWidth(1);
  extCableCu->SetFillColor(extCableCu->GetLineColor());
  extCableCu->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwCableFEP = new TGeoVolume("ITSsuppSSDSideAForwCableFEP",
					    forwPlastic, medFEP);

  forwCableFEP->SetVisibility(kTRUE);
  forwCableFEP->SetLineColor(kYellow); // Yellow
  forwCableFEP->SetLineWidth(1);
  forwCableFEP->SetFillColor(forwCableFEP->GetLineColor());
  forwCableFEP->SetFillStyle(4000); // 0% transparent

  TGeoVolume *extCableFEP = new TGeoVolume("ITSsuppSSDSideAExtCableFEP",
					   extPlastic, medFEP);

  extCableFEP->SetVisibility(kTRUE);
  extCableFEP->SetLineColor(kYellow); // Yellow
  extCableFEP->SetLineWidth(1);
  extCableFEP->SetFillColor(extCableFEP->GetLineColor());
  extCableFEP->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwTrayWater = new TGeoVolume("ITSsuppSSDSideAForwTrayWater",
					     forwWater, medH2O);

  forwTrayWater->SetVisibility(kTRUE);
  forwTrayWater->SetLineColor(kBlue); // Blue
  forwTrayWater->SetLineWidth(1);
  forwTrayWater->SetFillColor(forwTrayWater->GetLineColor());
  forwTrayWater->SetFillStyle(4000); // 0% transparent

  TGeoVolume *extTrayWater = new TGeoVolume("ITSsuppSSDSideAExtTrayWater",
					    extWater, medH2O);

  extTrayWater->SetVisibility(kTRUE);
  extTrayWater->SetLineColor(kBlue); // Blue
  extTrayWater->SetLineWidth(1);
  extTrayWater->SetFillColor(extTrayWater->GetLineColor());
  extTrayWater->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwPolyUr = new TGeoVolume("ITSsuppSSDSideAForwPolyUr",
					  forwPUR, medPUR);

  forwPolyUr->SetVisibility(kTRUE);
  forwPolyUr->SetLineColor(kGray); // Gray
  forwPolyUr->SetLineWidth(1);
  forwPolyUr->SetFillColor(forwPolyUr->GetLineColor());
  forwPolyUr->SetFillStyle(4000); // 0% transparent

  TGeoVolume *extPolyUr = new TGeoVolume("ITSsuppSSDSideAExtPolyUr",
					 extPUR, medPUR);

  extPolyUr->SetVisibility(kTRUE);
  extPolyUr->SetLineColor(kGray); // Gray
  extPolyUr->SetLineWidth(1);
  extPolyUr->SetFillColor(extPolyUr->GetLineColor());
  extPolyUr->SetFillStyle(4000); // 0% transparent


  // Now build up the tray
  cableTrayAForw->AddNode(forwTrayFirst, 1, 0);

  cableTrayAForw->AddNode(forwTraySecond, 1,
			new TGeoTranslation(0, 0, kForwardTrayFirstLen) );

  xloc = kTrayWidth/2 + kForwardSideThick/2;
  yloc = kForwardTrayFirstHeight + kForwardSideHeight/2 - kForwardSideYTrans;
  zloc = kForwardSideLength/2;
  cableTrayAForw->AddNode(forwTraySide,1,
			new TGeoTranslation( xloc, yloc, zloc) );
  cableTrayAForw->AddNode(forwTraySide,2,
			new TGeoTranslation(-xloc, yloc, zloc) );

  yloc = kForwardTrayFirstHeight + kForwardSideHeight - kForwardSideYTrans
       - kForwardCoverHeight;
  cableTrayAForw->AddNode(forwTraySideCover,1,
			new TGeoTranslation(0, yloc, 0) );

  yloc = kTrayTotalHeight - kCoversYTrans;
  zloc = kForwardTrayTotalLen - kForwardCoverLen;
  cableTrayAForw->AddNode(forwardTrayCover,1,
			new TGeoTranslation(0, yloc, zloc) );

  yloc = kTrayThick + forwCopper->GetDY();
  zloc = forwCopper->GetDZ();
  cableTrayAForw->AddNode(forwCableCu, 1,
			new TGeoTranslation(0, yloc, zloc) );

  yloc = kTrayThick + kCopperHeight + forwPlastic->GetDY();
  zloc = forwPlastic->GetDZ();
  cableTrayAForw->AddNode(forwCableFEP, 1,
			new TGeoTranslation(0, yloc, zloc) );

  yloc = kTrayThick + kCopperHeight + kCablePlasticHeight + forwWater->GetDY();
  zloc = forwWater->GetDZ();
  cableTrayAForw->AddNode(forwTrayWater, 1,
			new TGeoTranslation(0, yloc, zloc) );

  yloc = kTrayThick + kCopperHeight + kCablePlasticHeight
       + kCoolingWaterHeight + forwPUR->GetDY();
  zloc = forwPUR->GetDZ();
  cableTrayAForw->AddNode(forwPolyUr, 1,
			new TGeoTranslation(0, yloc, zloc) );

  // To simplify following placement in MARS, origin is on top
  totalhi = kTrayTotalHeight + kExternalCoverThick - kCoversYTrans;

  yloc = -totalhi;
  cableTrayAExt->AddNode(externalTraySSD, 1,
			new TGeoTranslation(0, yloc, 0) );

  yloc = -totalhi + kTrayTotalHeight - kCoversYTrans;
  cableTrayAExt->AddNode(externTrayCover,1,
			new TGeoTranslation(0, yloc, 0) );

  xloc = extCopper->GetDZ();
  yloc = -totalhi;
  cableTrayAExt->AddNode(extCableCu,1,
		        new TGeoCombiTrans( xloc, yloc, 0,
		        new TGeoRotation("",-90, 90, 90)        ) );

  xloc = extPlastic->GetDZ();
  yloc = -totalhi + kCopperHeight;
  cableTrayAExt->AddNode(extCableFEP,1,
		        new TGeoCombiTrans( xloc, yloc, 0,
		        new TGeoRotation("",-90, 90, 90)        ) );

  xloc = extWater->GetDZ();
  yloc = -totalhi + kCopperHeight + kCablePlasticHeight;
  cableTrayAExt->AddNode(extTrayWater,1,
		        new TGeoCombiTrans( xloc, yloc, 0,
		        new TGeoRotation("",-90, 90, 90)        ) );

  xloc = extPUR->GetDZ();
  yloc = -totalhi + kCopperHeight + kCablePlasticHeight + kCoolingWaterHeight;
  cableTrayAExt->AddNode(extPolyUr,1,
		        new TGeoCombiTrans( xloc, yloc, 0,
		        new TGeoRotation("",-90, 90, 90)        ) );


  // Finally put everything in the mother volume
  zloc = kTrayAZTrans;
  Double_t zlocext = zloc + kForwardTrayTotalLen;
  Double_t rExtTray = kTrayARTrans + kTrayTotalHeight;

  alpharot = kTrayAFirstRotAng;
  xloc = kTrayARTrans*SinD(alpharot);
  yloc = kTrayARTrans*CosD(alpharot);
  moth->AddNode(cableTrayAForw,1,
			    new TGeoCombiTrans( xloc, yloc, zloc,
			    new TGeoRotation("",-alpharot,0,0)   )   );
  xloc = rExtTray*SinD(alpharot);
  yloc = rExtTray*CosD(alpharot);
  moth->AddNode(cableTrayAExt,1,
			    new TGeoCombiTrans( xloc, yloc, zlocext,
			    new TGeoRotation("",-alpharot,-kTrayAZRot,0)  )  );

  alpharot += 180;
  xloc = kTrayARTrans*SinD(alpharot);
  yloc = kTrayARTrans*CosD(alpharot);
  moth->AddNode(cableTrayAForw,2,
			    new TGeoCombiTrans( xloc, yloc, zloc,
			    new TGeoRotation("",-alpharot,0,0)   )   );
  xloc = rExtTray*SinD(alpharot);
  yloc = rExtTray*CosD(alpharot);
  moth->AddNode(cableTrayAExt,2,
			    new TGeoCombiTrans( xloc, yloc, zlocext,
			    new TGeoRotation("",-alpharot,-kTrayAZRot,0)  )  );

  alpharot = -kTrayAFirstRotAng - 2*kTrayASecondRotAng;
  xloc = kTrayARTrans*SinD(alpharot);
  yloc = kTrayARTrans*CosD(alpharot);
  moth->AddNode(cableTrayAForw,3,
			    new TGeoCombiTrans( xloc, yloc, zloc,
			    new TGeoRotation("",-alpharot,0,0)   )   );
  xloc = rExtTray*SinD(alpharot);
  yloc = rExtTray*CosD(alpharot);
  moth->AddNode(cableTrayAExt,3,
			    new TGeoCombiTrans( xloc, yloc, zlocext,
			    new TGeoRotation("",-alpharot,-kTrayAZRot,0)  )  );

  alpharot += 180;
  xloc = kTrayARTrans*SinD(alpharot);
  yloc = kTrayARTrans*CosD(alpharot);
  moth->AddNode(cableTrayAForw,4,
			    new TGeoCombiTrans( xloc, yloc, zloc,
			    new TGeoRotation("",-alpharot,0,0)   )   );
  xloc = rExtTray*SinD(alpharot);
  yloc = rExtTray*CosD(alpharot);
  moth->AddNode(cableTrayAExt,4,
			    new TGeoCombiTrans( xloc, yloc, zlocext,
			    new TGeoRotation("",-alpharot,-kTrayAZRot,0)  )  );


  return;
}

//______________________________________________________________________
void AliITSv11GeometrySupport::SSDCableTraysSideC(TGeoVolume *moth,
					    const TGeoManager *mgr){
//
// Creates the SSD cable trays which are outside the ITS support cones
// but still inside the TPC on Side C
// (part of this code is taken or anyway inspired to ServicesCableSupport
// method of AliITSv11GeometrySupport.cxx,v 1.9 2007/06/06)
//
// Input:
//         moth : the TGeoVolume owing the volume structure
//         mgr  : the GeoManager (default gGeoManager)
// Output:
//
// Created:         ???       Bjorn S. Nilsen
// Updated:      15 Apr 2010  Mario Sitta
//
// Technical data are taken from AutoCAD drawings and other (oral)
// information given by F.Tosello
//

  // Dimensions and positions of the C-Side Cable Tray elements
  const Int_t    kNumTraySideC           =    4;

  const Double_t kSideCFoldAngle         =    5.00 *fgkDegree;

  const Double_t kServicesWidth          =  100.00 *fgkmm;
  const Double_t kCopperHeight           =   11.20 *fgkmm;// 1120 mm^2
  const Double_t kCablePlasticHeight     =   11.50 *fgkmm;// 1150 mm^2
  const Double_t kCoolingWaterHeight     =    2.65 *fgkmm;//  265 mm^2
  const Double_t kPoliUrethaneHeight     =    4.62 *fgkmm;//  462 mm^2

  // Overall position and rotation of the C-Side Cable Trays
  const Double_t kTraySideCRPos          =   45.30    *fgkcm;
  const Double_t kTraySideCZPos          = -102.40    *fgkcm;
  const Double_t kTraySideCAlphaRot[kNumTraySideC]  = {     23.0,     -59.0,
    /* from Patch panel position */		       180.+23.0, 180.-59.0};


  // Local variables
  Double_t xprof[6], yprof[6];
  Double_t xloc, yloc, alpharot, alphafold;


  // The assembly holding the metallic structure
  TGeoVolumeAssembly *trayStructure =
				CreateSDDSSDTraysSideC("ITSsupportSSDTrayC");

  // The cable copper inside the tray: a Xtru
  TGeoXtru *copper = new TGeoXtru(2);
  copper->SetName("ITSsuppSSDTrayCCopper");

  // Copper lies on the lower plate: get position of its points
  TGeoXtru *lowerplate = (TGeoXtru*)(mgr->GetVolume("ITSsuppTraySideCLower")->GetShape());
  xprof[0] = lowerplate->GetX(5);
  yprof[0] = lowerplate->GetY(5);
  xprof[1] = lowerplate->GetX(4);
  yprof[1] = lowerplate->GetY(4);
  xprof[2] = lowerplate->GetX(3);
  yprof[2] = lowerplate->GetY(3);
  xprof[3] = xprof[2] - kCopperHeight*SinD(kSideCFoldAngle);
  yprof[3] = yprof[2] + kCopperHeight*CosD(kSideCFoldAngle);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      kCopperHeight , xprof[4], yprof[4]);
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + kCopperHeight;

  copper->DefinePolygon(6, xprof, yprof);
  copper->DefineSection(0, -kServicesWidth/2);
  copper->DefineSection(1,  kServicesWidth/2);

  // The cable plastic inside the tray: a Xtru
  TGeoXtru *plastic = new TGeoXtru(2);
  plastic->SetName("ITSsuppSSDTrayCPlastic");

  xprof[0] = copper->GetX(5);
  yprof[0] = copper->GetY(5);
  xprof[1] = copper->GetX(4);
  yprof[1] = copper->GetY(4);
  xprof[2] = copper->GetX(3);
  yprof[2] = copper->GetY(3);
  xprof[3] = xprof[2] - kCablePlasticHeight*SinD(kSideCFoldAngle);
  yprof[3] = yprof[2] + kCablePlasticHeight*CosD(kSideCFoldAngle);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      kCablePlasticHeight , xprof[4], yprof[4]);
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + kCablePlasticHeight;

  plastic->DefinePolygon(6, xprof, yprof);
  plastic->DefineSection(0, -kServicesWidth/2);
  plastic->DefineSection(1,  kServicesWidth/2);

  // The cooling water inside the tray: a Xtru
  TGeoXtru *water = new TGeoXtru(2);
  water->SetName("ITSsuppSSDTrayCWater");

  xprof[0] = plastic->GetX(5);
  yprof[0] = plastic->GetY(5);
  xprof[1] = plastic->GetX(4);
  yprof[1] = plastic->GetY(4);
  xprof[2] = plastic->GetX(3);
  yprof[2] = plastic->GetY(3);
  xprof[3] = xprof[2] - kCoolingWaterHeight*SinD(kSideCFoldAngle);
  yprof[3] = yprof[2] + kCoolingWaterHeight*CosD(kSideCFoldAngle);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      kCoolingWaterHeight , xprof[4], yprof[4]);
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + kCoolingWaterHeight;

  water->DefinePolygon(6, xprof, yprof);
  water->DefineSection(0, -kServicesWidth/2);
  water->DefineSection(1,  kServicesWidth/2);

  // The poliurethane inside the tray: a Xtru
  TGeoXtru *pur = new TGeoXtru(2);
  pur->SetName("ITSsuppSSDTrayCPUR");
  xprof[0] = water->GetX(5);
  yprof[0] = water->GetY(5);
  xprof[1] = water->GetX(4);
  yprof[1] = water->GetY(4);
  xprof[2] = water->GetX(3);
  yprof[2] = water->GetY(3);
  xprof[3] = xprof[2] - kPoliUrethaneHeight*SinD(kSideCFoldAngle);
  yprof[3] = yprof[2] + kPoliUrethaneHeight*CosD(kSideCFoldAngle);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      kPoliUrethaneHeight , xprof[4], yprof[4]);
  xprof[5] = xprof[0];
  yprof[5] = yprof[0] + kPoliUrethaneHeight;

  pur->DefinePolygon(6, xprof, yprof);
  pur->DefineSection(0, -kServicesWidth/2);
  pur->DefineSection(1,  kServicesWidth/2);


  // We have all shapes: now create the real volumes
  TGeoMedium *medCu    = mgr->GetMedium("ITS_COPPER$");
  TGeoMedium *medFEP   = mgr->GetMedium("ITS_SSD FEP$");
  TGeoMedium *medH2O   = mgr->GetMedium("ITS_WATER$");
  TGeoMedium *medPUR   = mgr->GetMedium("ITS_POLYURETHANE$");

  TGeoVolume *copperCable = new TGeoVolume("ITSsuppSSDSideCCableCu",
					   copper, medCu);

  copperCable->SetVisibility(kTRUE);
  copperCable->SetLineColor(kRed); // Red
  copperCable->SetLineWidth(1);
  copperCable->SetFillColor(copperCable->GetLineColor());
  copperCable->SetFillStyle(4000); // 0% transparent

  TGeoVolume *cableFEP = new TGeoVolume("ITSsuppSSDSideCCableFEP",
					plastic, medFEP);

  cableFEP->SetVisibility(kTRUE);
  cableFEP->SetLineColor(kYellow); // Yellow
  cableFEP->SetLineWidth(1);
  cableFEP->SetFillColor(cableFEP->GetLineColor());
  cableFEP->SetFillStyle(4000); // 0% transparent

  TGeoVolume *trayWater = new TGeoVolume("ITSsuppSSDSideCTrayWater",
					 water, medH2O);

  trayWater->SetVisibility(kTRUE);
  trayWater->SetLineColor(kBlue); // Blue
  trayWater->SetLineWidth(1);
  trayWater->SetFillColor(trayWater->GetLineColor());
  trayWater->SetFillStyle(4000); // 0% transparent

  TGeoVolume *trayPolyUr = new TGeoVolume("ITSsuppSSDSideCPolyUr",
					  pur, medPUR);

  trayPolyUr->SetVisibility(kTRUE);
  trayPolyUr->SetLineColor(kGray); // Gray
  trayPolyUr->SetLineWidth(1);
  trayPolyUr->SetFillColor(trayPolyUr->GetLineColor());
  trayPolyUr->SetFillStyle(4000); // 0% transparent


  // Now fill in the tray
  trayStructure->AddNode(copperCable,1,0);
  trayStructure->AddNode(cableFEP,1,0);
  trayStructure->AddNode(trayWater,1,0);
  trayStructure->AddNode(trayPolyUr,1,0);


  // Finally put everything in the mother volume
  alphafold = kSideCFoldAngle;

  for (Int_t jt = 0; jt < kNumTraySideC; jt++) {
    alpharot = kTraySideCAlphaRot[jt];
    xloc = kTraySideCRPos*SinD(alpharot);
    yloc = kTraySideCRPos*CosD(alpharot);
    moth->AddNode(trayStructure,jt+1,
		       new TGeoCombiTrans(-xloc, yloc, kTraySideCZPos,
		       new TGeoRotation("",-90.+alpharot,-90.,90.+alphafold)));
  }


  return;
}

//______________________________________________________________________
void AliITSv11GeometrySupport::CreateSDDForwardTraySideA(TGeoVolumeAssembly *tray,
						   const TGeoManager *mgr){
//
// Creates the forward SDD tray on Side A (0872/G/D/01)
//
// Input:
//         tray : the TGeoVolumeAssembly to put the elements in
//         mgr  : the GeoManager (used only to get the proper material)
//
// Output:
//
// Return:
//
// Created:      08 Jan 2010  Mario Sitta
// Updated:      07 Sep 2010  Mario Sitta
//
// Technical data are taken from AutoCAD drawings, L.Simonetti technical
// drawings and other (oral) information given by F.Tosello
//

  // Dimensions of the A-Side Forward Cable Tray (0872/G/D/01)
  const Double_t kForwardTrayThick        =    2.00 *fgkmm;
  const Double_t kForwardTraySideLength   =  823.00 *fgkmm;
  const Double_t kForwardTrayTailLength   =  212.00 *fgkmm;
  const Double_t kForwardTrayBaseHalfWide =   55.00 *fgkmm;
  const Double_t kForwardTrayNotchLength  =   47.20 *fgkmm;
  const Double_t kForwardTrayNotchHeight  =   25.00 *fgkmm;
  const Double_t kForwardTrayNotchDown    =   10.00 *fgkmm;
  const Double_t kForwardTraySide1Height  =   39.00 *fgkmm;
  const Double_t kForwardTraySide2Height  =   26.00 *fgkmm;
  const Double_t kForwardTraySide2Expand  =   10.50 *fgkmm;
  const Double_t kForwardTraySide3TailLen =  418.00 *fgkmm;
  const Double_t kForwardTraySide3TailHi  =   31.00 *fgkmm;
  const Double_t kForwardTraySide3HeadLen =  425.00 *fgkmm;
  const Double_t kForwardTraySide3HeadHi  =   72.00 *fgkmm;
  const Double_t kForwardTrayHorWingWide  =   10.50 *fgkmm;
  const Double_t kForwardTrayVertWingWide =   15.00 *fgkmm;

  const Int_t    kForwardTraySideNpoints  =    9;


  // Local variables
  Double_t xprof[kForwardTraySideNpoints], yprof[kForwardTraySideNpoints];
  Double_t ylen, zlen;
  Double_t xloc, yloc, zloc;


  // The tray has a very complex shape, so it is made by assembling
  // different elements (with some small simplifications)

  // The tray base: a BBox
  zlen = (kForwardTraySideLength-kForwardTrayTailLength)/2;
  TGeoBBox *trayBase = new TGeoBBox(kForwardTrayBaseHalfWide,
				    kForwardTrayThick/2, zlen);

  // The first part of the side wall: a Xtru
  TGeoXtru *traySide1 = new TGeoXtru(2);

  xprof[0] = 0;
  yprof[0] = kForwardTrayThick;
  xprof[1] = kForwardTraySideLength-kForwardTrayTailLength;
  yprof[1] = yprof[0];
  xprof[2] = kForwardTraySideLength;
  yprof[2] = kForwardTraySide1Height + kForwardTrayThick;
  xprof[3] = 0;
  yprof[3] = yprof[2];

  traySide1->DefinePolygon(4, xprof, yprof);
  traySide1->DefineSection(0, 0);
  traySide1->DefineSection(1, kForwardTrayThick);

  // The second part of the side wall: a Xtru
  TGeoXtru *traySide2 = new TGeoXtru(2);

  xprof[0] = kForwardTrayBaseHalfWide - kForwardTrayThick;
  yprof[0] = traySide1->GetY(2);
  xprof[1] = kForwardTrayBaseHalfWide;
  yprof[1] = yprof[0];
  xprof[2] = xprof[1] + kForwardTraySide2Expand;
  yprof[2] = yprof[1] + kForwardTraySide2Height;
  xprof[3] = xprof[2] - kForwardTrayThick;
  yprof[3] = yprof[2];

  traySide2->DefinePolygon(4, xprof, yprof);
  traySide2->DefineSection(0, 0);
  traySide2->DefineSection(1, kForwardTraySideLength);

  // The third part of the side wall: a Xtru
  TGeoXtru *traySide3 = new TGeoXtru(2);

  xprof[0] = 0;
  yprof[0] = traySide2->GetY(2);
  xprof[1] = kForwardTraySideLength;
  yprof[1] = yprof[0];
  xprof[2] = xprof[1];
  yprof[2] = yprof[1] + kForwardTraySide3TailHi - kForwardTrayThick;
  xprof[3] = xprof[2] - kForwardTraySide3TailLen - kForwardTrayThick;
  yprof[3] = yprof[2];
  xprof[4] = xprof[3];
  yprof[4] = yprof[3] + kForwardTraySide3HeadHi + kForwardTrayThick;
  xprof[5] = xprof[4] - kForwardTraySide3HeadLen;
  yprof[5] = yprof[4];
  xprof[6] = xprof[5];
  yprof[6] = yprof[5] - kForwardTrayNotchHeight;
  xprof[7] = xprof[6] + kForwardTrayNotchLength;
  yprof[7] = yprof[6];
  xprof[8] = xprof[7];
  yprof[8] = yprof[7] - kForwardTrayNotchDown;

  traySide3->DefinePolygon(9, xprof, yprof);
  traySide3->DefineSection(0, 0);
  traySide3->DefineSection(1, kForwardTrayThick);

  // The horizontal wing: a BBox
  TGeoBBox *trayHorWing = new TGeoBBox(kForwardTrayHorWingWide/2,
				       kForwardTrayThick/2,
				       kForwardTraySide3TailLen/2);

  // The vertical wing: a BBox
  ylen = (traySide3->GetY(4) - traySide3->GetY(3))/2;
  TGeoBBox *trayVertWing = new TGeoBBox(kForwardTrayVertWingWide/2,
					ylen, kForwardTrayThick/2);


  // We have all shapes: now create the real volumes
  TGeoMedium *medAl    = mgr->GetMedium("ITS_ALUMINUM$");

  TGeoVolume *forwTrayBase = new TGeoVolume("ITSsuppSDDSideAForwTrayBase",
					    trayBase, medAl);

  forwTrayBase->SetVisibility(kTRUE);
  forwTrayBase->SetLineColor(6); // Purple
  forwTrayBase->SetLineWidth(1);
  forwTrayBase->SetFillColor(forwTrayBase->GetLineColor());
  forwTrayBase->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwTraySide1 = new TGeoVolume("ITSsuppSDDSideAForwTraySide1",
					    traySide1, medAl);

  forwTraySide1->SetVisibility(kTRUE);
  forwTraySide1->SetLineColor(6); // Purple
  forwTraySide1->SetLineWidth(1);
  forwTraySide1->SetFillColor(forwTraySide1->GetLineColor());
  forwTraySide1->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwTraySide2 = new TGeoVolume("ITSsuppSDDSideAForwTraySide2",
					    traySide2, medAl);

  forwTraySide2->SetVisibility(kTRUE);
  forwTraySide2->SetLineColor(6); // Purple
  forwTraySide2->SetLineWidth(1);
  forwTraySide2->SetFillColor(forwTraySide2->GetLineColor());
  forwTraySide2->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwTraySide3 = new TGeoVolume("ITSsuppSDDSideAForwTraySide3",
					    traySide3, medAl);

  forwTraySide3->SetVisibility(kTRUE);
  forwTraySide3->SetLineColor(6); // Purple
  forwTraySide3->SetLineWidth(1);
  forwTraySide3->SetFillColor(forwTraySide3->GetLineColor());
  forwTraySide3->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwTrayHWing = new TGeoVolume("ITSsuppSDDSideAForwTrayHorWing",
					    trayHorWing, medAl);

  forwTrayHWing->SetVisibility(kTRUE);
  forwTrayHWing->SetLineColor(6); // Purple
  forwTrayHWing->SetLineWidth(1);
  forwTrayHWing->SetFillColor(forwTrayHWing->GetLineColor());
  forwTrayHWing->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwTrayVWing = new TGeoVolume("ITSsuppSDDSideAForwTrayVertWing",
					    trayVertWing, medAl);

  forwTrayVWing->SetVisibility(kTRUE);
  forwTrayVWing->SetLineColor(6); // Purple
  forwTrayVWing->SetLineWidth(1);
  forwTrayVWing->SetFillColor(forwTrayVWing->GetLineColor());
  forwTrayVWing->SetFillStyle(4000); // 0% transparent


  // Now build up the tray
  yloc = kForwardTrayThick/2;
  zloc = zlen;
  tray->AddNode(forwTrayBase, 1,
		new TGeoTranslation(0, yloc, zloc) );

  xloc = kForwardTrayBaseHalfWide;
  tray->AddNode(forwTraySide1, 1,
		new TGeoCombiTrans(xloc, 0, 0,
				   new TGeoRotation("",90,-90,-90)));
  xloc = -xloc + kForwardTrayThick;
  tray->AddNode(forwTraySide1, 2,
		new TGeoCombiTrans(xloc, 0, 0,
				   new TGeoRotation("",90,-90,-90)));

  tray->AddNode(forwTraySide2, 1, 0);
  zloc = kForwardTraySideLength;
  tray->AddNode(forwTraySide2, 2,
		new TGeoCombiTrans(0, 0, zloc,
				   new TGeoRotation("",90,-180,-90)));

  xloc = kForwardTrayBaseHalfWide + kForwardTraySide2Expand;
  tray->AddNode(forwTraySide3, 1,
		new TGeoCombiTrans(xloc, 0, 0,
				   new TGeoRotation("",90,-90,-90)));
  xloc = -xloc + kForwardTrayThick;
  tray->AddNode(forwTraySide3, 2,
		new TGeoCombiTrans(xloc, 0, 0,
				   new TGeoRotation("",90,-90,-90)));

  xloc = kForwardTrayBaseHalfWide + kForwardTraySide2Expand
       - kForwardTrayHorWingWide/2;
  yloc = traySide3->GetY(2) + kForwardTrayThick/2;
  zloc = kForwardTraySideLength - trayHorWing->GetDZ();
  tray->AddNode(forwTrayHWing, 1,
		new TGeoTranslation( xloc, yloc, zloc) );
  tray->AddNode(forwTrayHWing, 2,
		new TGeoTranslation(-xloc, yloc, zloc) );

  xloc = kForwardTrayBaseHalfWide + kForwardTraySide2Expand
       - kForwardTrayVertWingWide/2;
  yloc = traySide3->GetY(2) + trayVertWing->GetDY();
  zloc = traySide3->GetX(3) + kForwardTrayThick/2;
  tray->AddNode(forwTrayVWing, 1,
		new TGeoTranslation( xloc, yloc, zloc) );
  tray->AddNode(forwTrayVWing, 2,
		new TGeoTranslation(-xloc, yloc, zloc) );


  return;
}

//______________________________________________________________________
TGeoCompositeShape* AliITSv11GeometrySupport::CreateTrayAForwardCover(const Double_t coverLen){
//
// Creates the forward cover of the SDD and SSD cable trays on Side A
// (0872/G/D/02)
//
// Input:
//             coverLen: the total length of the cover
//
// Output:
//
// Return:     a TGeoCompositeShape for the cover
//
// Created:      03 Jan 2010  Mario Sitta
//
// Technical data are taken from AutoCAD drawings, L.Simonetti technical
// drawings and other (oral) information given by F.Tosello
//

  // Dimensions and positions of the A-Side Cable Tray Forward Cover
  // (0872/G/D/02)
  const Double_t kForwardCoverWide        =  130.00 *fgkmm;
  const Double_t kForwardCoverSideWide    =   10.00 *fgkmm;
  const Double_t kForwardCoverHoleLen     =  160.00 *fgkmm;
  const Double_t kForwardCoverHoleWide    =   90.00 *fgkmm;
  const Double_t kForwardCoverHoleR10     =   10.00 *fgkmm;
  const Double_t kForwardCoverTotalThick  =    5.00 *fgkmm;
  const Double_t kForwardCoverSideThick   =    3.00 *fgkmm;
  const Double_t kForwardCoverInternThick =    2.00 *fgkmm;

  const Double_t kForwardCoverHoleZTrans  =   40.00 *fgkmm;


  // Local variables
  Double_t xprof[16], yprof[16];
  Double_t yloc, zloc;


  // The main shape: a Xtru
  TGeoXtru *forwCoverMain = new TGeoXtru(2);
  forwCoverMain->SetName("ITSsuppForwCoverMain");

  xprof[0] = kForwardCoverWide/2;
  yprof[0] = kForwardCoverTotalThick;
  xprof[1] = xprof[0];
  yprof[1] = yprof[0] - kForwardCoverSideThick;
  xprof[2] = xprof[1] - kForwardCoverSideWide;
  yprof[2] = yprof[1];
  xprof[3] = xprof[2];
  yprof[3] = 0;

  // We did the right side, now reflex on the left side
  for (Int_t jp = 0; jp < 4; jp++) {
    xprof[4+jp] = -xprof[3-jp];
    yprof[4+jp] =  yprof[3-jp];
  }

  // And now the actual Xtru
  forwCoverMain->DefinePolygon(8, xprof, yprof);
  forwCoverMain->DefineSection(0, 0);
  forwCoverMain->DefineSection(1, coverLen);

  // The hole: another Xtru (rounded corners approximated with segments)
  TGeoXtru *forwCoverHole = new TGeoXtru(2);
  forwCoverHole->SetName("ITSsuppForwCoverHole");

  CreateTrayACoverHolesShape(kForwardCoverHoleWide, kForwardCoverHoleLen,
			     kForwardCoverHoleR10 , xprof, yprof);

  // And now the actual Xtru
  forwCoverHole->DefinePolygon(16, xprof, yprof);
  forwCoverHole->DefineSection(0, 0);
  forwCoverHole->DefineSection(1, kForwardCoverTotalThick-kForwardCoverInternThick);

  // Now the proper rototranslation matrices for the two holes
  yloc = kForwardCoverTotalThick-kForwardCoverInternThick-0.01;//Precision fix
  zloc = kForwardCoverHoleZTrans;
  TGeoCombiTrans *mf1 = new TGeoCombiTrans(0, yloc, zloc,
					   new TGeoRotation("", 0, 90, 0) );
  mf1->SetName("mf1");
  mf1->RegisterYourself();

  zloc = coverLen - kForwardCoverHoleZTrans - kForwardCoverHoleLen;
  TGeoCombiTrans *mf2 = new TGeoCombiTrans(0, yloc, zloc,
					   new TGeoRotation("", 0, 90, 0) );
  mf2->SetName("mf2");
  mf2->RegisterYourself();

  // Finally the actual cover shape
  TGeoCompositeShape *cover = new TGeoCompositeShape("ITSsuppForwardCoverMain",
    "ITSsuppForwCoverMain-ITSsuppForwCoverHole:mf1-ITSsuppForwCoverHole:mf2");

  return cover;
}

//______________________________________________________________________
TGeoCompositeShape* AliITSv11GeometrySupport::CreateTrayAExternalCover(const Double_t coverLen){
//
// Creates the external cover of the SDD and SSD cable trays on Side A
// (0872/G/D/04)
//
// Input:
//             coverLen: the total length of the cover
//
// Output:
//
// Return:     a TGeoCompositeShape for the cover
//
// Created:      03 Jan 2010  Mario Sitta
//
// Technical data are taken from AutoCAD drawings, L.Simonetti technical
// drawings and other (oral) information given by F.Tosello
//

  // Dimensions and positions of the A-Side Cable Tray External Cover
  // (0872/G/D/04)
  const Double_t kExternalCoverWide        =  130.00 *fgkmm;
  const Double_t kExternalCoverSideWide    =   10.00 *fgkmm;
  const Double_t kExternalCoverHoleLen1    =  262.00 *fgkmm;
  const Double_t kExternalCoverHoleLen2    =  280.00 *fgkmm;
  const Double_t kExternalCoverHoleLen3    =  205.00 *fgkmm;
  const Double_t kExternalCoverHoleLen4    =   55.00 *fgkmm;
  const Double_t kExternalCoverHoleWide    =   90.00 *fgkmm;
  const Double_t kExternalCoverHoleR10     =   10.00 *fgkmm;
  const Double_t kExternalCoverTotalThick  =    5.00 *fgkmm;
  const Double_t kExternalCoverSideThick   =    3.00 *fgkmm;
  const Double_t kExternalCoverInternThick =    2.00 *fgkmm;

  const Double_t kExternalCoverHole1ZTrans =   28.00 *fgkmm;
  const Double_t kExternalCoverHolesZTrans =   20.00 *fgkmm;


  // Local variables
  Double_t xprof[16], yprof[16];
  Double_t yloc, zloc;


  // The main shape: a Xtru
  TGeoXtru *externCoverMain = new TGeoXtru(2);
  externCoverMain->SetName("ITSsuppExternCoverMain");

  xprof[0] = kExternalCoverWide/2;
  yprof[0] = kExternalCoverTotalThick;
  xprof[1] = xprof[0];
  yprof[1] = yprof[0] - kExternalCoverSideThick;
  xprof[2] = xprof[1] - kExternalCoverSideWide;
  yprof[2] = yprof[1];
  xprof[3] = xprof[2];
  yprof[3] = 0;

  // We did the right side, now reflex on the left side
  for (Int_t jp = 0; jp < 4; jp++) {
    xprof[4+jp] = -xprof[3-jp];
    yprof[4+jp] =  yprof[3-jp];
  }

  // And now the actual Xtru
  externCoverMain->DefinePolygon(8, xprof, yprof);
  externCoverMain->DefineSection(0, 0);
  externCoverMain->DefineSection(1, coverLen);

  // The first hole: a Xtru (rounded corners approximated with segments)
  Double_t holethick = kExternalCoverTotalThick-kExternalCoverInternThick;

  TGeoXtru *extCoverHole1 = new TGeoXtru(2);
  extCoverHole1->SetName("ITSsuppExtCoverHole1");

  CreateTrayACoverHolesShape(kExternalCoverHoleWide, kExternalCoverHoleLen1,
			     kExternalCoverHoleR10 , xprof, yprof);

  extCoverHole1->DefinePolygon(16, xprof, yprof);
  extCoverHole1->DefineSection(0, 0);
  extCoverHole1->DefineSection(1, holethick);

  // The second (and third) hole: another Xtru
  TGeoXtru *extCoverHole2 = new TGeoXtru(2);
  extCoverHole2->SetName("ITSsuppExtCoverHole2");

  CreateTrayACoverHolesShape(kExternalCoverHoleWide, kExternalCoverHoleLen2,
			     kExternalCoverHoleR10 , xprof, yprof);

  extCoverHole2->DefinePolygon(16, xprof, yprof);
  extCoverHole2->DefineSection(0, 0);
  extCoverHole2->DefineSection(1, holethick);

  // The fourth hole: another Xtru
  TGeoXtru *extCoverHole3 = new TGeoXtru(2);
  extCoverHole3->SetName("ITSsuppExtCoverHole3");

  CreateTrayACoverHolesShape(kExternalCoverHoleWide, kExternalCoverHoleLen3,
			     kExternalCoverHoleR10 , xprof, yprof);

  extCoverHole3->DefinePolygon(16, xprof, yprof);
  extCoverHole3->DefineSection(0, 0);
  extCoverHole3->DefineSection(1, holethick);

  // The fifth and last hole: another Xtru
  TGeoXtru *extCoverHole4 = new TGeoXtru(2);
  extCoverHole4->SetName("ITSsuppExtCoverHole4");

  CreateTrayACoverHolesShape(kExternalCoverHoleWide, kExternalCoverHoleLen4,
			     kExternalCoverHoleR10 , xprof, yprof);

  extCoverHole4->DefinePolygon(16, xprof, yprof);
  extCoverHole4->DefineSection(0, 0);
  extCoverHole4->DefineSection(1, holethick);

  // Now the proper rototranslation matrices for the holes
  yloc = kExternalCoverTotalThick - kExternalCoverInternThick-0.01;
  zloc = kExternalCoverHole1ZTrans;
  TGeoCombiTrans *me1 = new TGeoCombiTrans(0, yloc, zloc,
					   new TGeoRotation("", 0, 90, 0) );
  me1->SetName("me1");
  me1->RegisterYourself();

  zloc += (kExternalCoverHoleLen1 + kExternalCoverHolesZTrans);
  TGeoCombiTrans *me2 = new TGeoCombiTrans(0, yloc, zloc,
					   new TGeoRotation("", 0, 90, 0) );
  me2->SetName("me2");
  me2->RegisterYourself();

  zloc += (kExternalCoverHoleLen2 + kExternalCoverHolesZTrans);
  TGeoCombiTrans *me3 = new TGeoCombiTrans(0, yloc, zloc,
					   new TGeoRotation("", 0, 90, 0) );
  me3->SetName("me3");
  me3->RegisterYourself();

  zloc += (kExternalCoverHoleLen2 + kExternalCoverHolesZTrans);
  TGeoCombiTrans *me4 = new TGeoCombiTrans(0, yloc, zloc,
					   new TGeoRotation("", 0, 90, 0) );
  me4->SetName("me4");
  me4->RegisterYourself();

  zloc += (kExternalCoverHoleLen3 + kExternalCoverHolesZTrans);
  TGeoCombiTrans *me5 = new TGeoCombiTrans(0, yloc, zloc,
					   new TGeoRotation("", 0, 90, 0) );
  me5->SetName("me5");
  me5->RegisterYourself();

  // Finally the actual cover shape
  TGeoCompositeShape *cover = new TGeoCompositeShape("ITSsuppExternCoverMain",
    "ITSsuppExternCoverMain-ITSsuppExtCoverHole1:me1-ITSsuppExtCoverHole2:me2-ITSsuppExtCoverHole2:me3-ITSsuppExtCoverHole3:me4-ITSsuppExtCoverHole4:me5");

  return cover;
}

//______________________________________________________________________
void AliITSv11GeometrySupport::CreateTrayACoverHolesShape(const Double_t wide,
			       const Double_t length, const Double_t r10,
			       Double_t *x, Double_t *y){
//
// Creates the proper sequence of X and Y coordinates to determine
// the base XTru polygon for the holes in the SDD and SSD tray covers
// (here the rounded corners are approximated with segments)
//
// Input:
//        wide   : the hole wide
//        length : the hole length
//        r10    : the radius of the rounded corners
//
// Output:
//        x, y : coordinate vectors [16]
//
// Created:      03 Jan 2010  Mario Sitta
//
// Caller must guarantee that x and y have the correct dimensions
// (but being this a private method it's easy to tell)
//

  x[0] = wide/2 - r10;
  y[0] = length;
  x[1] = x[0] + r10*SinD(30);
  y[1] = y[0] - r10*(1 - CosD(30));
  x[2] = x[0] + r10*SinD(60);
  y[2] = y[0] - r10*(1 - CosD(60));
  x[3] = x[0] + r10;
  y[3] = y[0] - r10;
  x[4] = x[3];
  y[4] = r10;
  x[5] = x[4] - r10*(1 - CosD(30));
  y[5] = y[4] - r10*SinD(30);
  x[6] = x[4] - r10*(1 - CosD(60));
  y[6] = y[4] - r10*SinD(60);
  x[7] = x[4] - r10;
  y[7] = 0;

  // We did the right side, now reflex on the left side
  for (Int_t jp = 0; jp < 8; jp++) {
    x[8+jp] = -x[7-jp];
    y[8+jp] =  y[7-jp];
  }

  return;
}

//______________________________________________________________________
TGeoXtru* AliITSv11GeometrySupport::CreateSDDSSDTraysSideA(
					      const Double_t trayLen,
					      const Double_t trayHi){
//
// Creates parts of the SDD and SSD Trays on Side A which are identical
// (0872/G/D/03, part of 0872/G/D/07, 0872/G/C/11)
//
// Input:
//         trayLen : the length of the tray part
//         trayHi  : the height of the tray part
//
// Output:
//
// Return:     a TGeoXtru
//
// Created:      26 Feb 2010  Mario Sitta
//
// Technical data are taken from AutoCAD drawings, L.Simonetti technical
// drawings and other (oral) information given by F.Tosello
//

  // Dimensions and positions of the A-Side Cable Trays
  // (parts of 0872/G/C)
  const Double_t kTrayWidth              =  130.00 *fgkmm;
  const Double_t kTrayWingWidth          =   10.00 *fgkmm;
  const Double_t kTrayHeightToBend       =   20.00 *fgkmm;
  const Double_t kTrayThick              =    2.00 *fgkmm;

  const Double_t kTrayBendAngle          =   22.00 *TMath::DegToRad();

  const Int_t    kTrayNpoints            =   16;

  // Local variables
  Double_t xprof[kTrayNpoints], yprof[kTrayNpoints];


  // The tray shape: a Xtru
  TGeoXtru *trayPart = new TGeoXtru(2);

  xprof[2] = kTrayWidth/2 - kTrayThick;
  yprof[2] = trayHi - kTrayThick;
  xprof[3] = kTrayWidth/2 - kTrayWingWidth;
  yprof[3] = yprof[2];
  xprof[4] = xprof[3];
  yprof[4] = trayHi;
  xprof[5] = kTrayWidth/2;
  yprof[5] = yprof[4];
  xprof[6] = xprof[5];
  yprof[6] = kTrayHeightToBend;
  xprof[7] = xprof[6] - yprof[6]*TMath::Tan(kTrayBendAngle);
  yprof[7] = 0;

  InsidePoint( xprof[5], yprof[5], xprof[6], yprof[6], xprof[7], yprof[7],
	      -kTrayThick, xprof[1], yprof[1]);

  xprof[8] = -xprof[7];
  yprof[8] =  yprof[7];

  InsidePoint( xprof[6], yprof[6], xprof[7], yprof[7], xprof[8], yprof[8],
	      -kTrayThick, xprof[0], yprof[0]);

  // We did the right side, now reflex on the left side
  for (Int_t jp = 0; jp < 8; jp++) {
    xprof[8+jp] = -xprof[7-jp];
    yprof[8+jp] =  yprof[7-jp];
  }

  // And now the actual Xtru
  trayPart->DefinePolygon(kTrayNpoints, xprof, yprof);
  trayPart->DefineSection(0, 0);
  trayPart->DefineSection(1, trayLen);


  return trayPart;
}

//______________________________________________________________________
TGeoVolumeAssembly* AliITSv11GeometrySupport::CreateSDDSSDTraysSideC(
						       const char *trayName,
						       const TGeoManager *mgr){

//
// Creates the SDD and SSD Trays on Side C which are supposedly identical
//
// Input:
//         trayName : the assembly name
//
// Output:
//
// Return:     a TGeoVolumeAssembly
//
// Created:      16 Apr 2010  Mario Sitta
//
// Technical data are taken from AutoCAD drawings and other (oral)
// information given by F.Tosello
//

  const Double_t kSideCHalfThick      =    0.100   *fgkcm;
  const Double_t kSideCFoldAngle      =    5.000   *TMath::DegToRad();

  const Double_t kSideCLength1        =  172.800   *fgkcm;
  const Double_t kSideCLength2        =  189.300   *fgkcm;
  const Double_t kSideCHalfWide       =    6.350   *fgkcm;
  const Double_t kSideCHeight1        =   11.800   *fgkcm;
  const Double_t kSideCHeight2        =    4.300   *fgkcm;
  const Double_t kSideCSideLength1    =   10.800   *fgkcm;
  const Double_t kSideCSideLength2    =   63.800   *fgkcm;
  const Double_t kSideCSideHeight     =    8.800   *fgkcm;
  const Int_t    kNPointsLowerFace    =    6;
  const Int_t    kNPointsLateralFace  =    9;

  const Double_t kSideCWingAHalfLen   =    5.000   *fgkcm;
  const Double_t kSideCWingBHalfLen   =   30.500   *fgkcm;
  const Double_t kSideCWingCHalfLen   =    2.000   *fgkcm;
  const Double_t kSideCWingDHalfLen   =   48.500   *fgkcm;
  const Double_t kSideCWingEHalfLen   =   83.000   *fgkcm;
  const Double_t kSideCWingsHalfWide  =    0.450   *fgkcm;

  const Int_t    kNPointsCoverFace    =   12;

  const Double_t kPlateHalfLen        =    6.000   *fgkcm;
  const Double_t kPlateThick          =    0.600   *fgkcm;
  const Double_t kPlateHeight         =    4.200   *fgkcm;
  const Int_t    kNPointsPlate        =    6;

  const Double_t kBarCoolRmax         =    0.4     *fgkcm;
  const Int_t    kNumBarCool          =    2;
  const Double_t kXShiftBarCool[kNumBarCool] = { 8.7, 13.0 };
  const Double_t kYShiftBarCool[kNumBarCool] = { 8.5,  5.0 };


  // Local variables
  Double_t xprof[12], yprof[12];
  Double_t xloc, yloc, zloc, delta, alpharot;

  // The single C-Side Cable tray as an assembly
  TGeoVolumeAssembly *cableTrayC = new TGeoVolumeAssembly(trayName);

  // First create all needed shapes

  // The Cable Tray lower face: a Xtru
  TGeoXtru *sideCLowerFace = new TGeoXtru(2);

  xprof[0] = 0.;
  yprof[0] = 0.;
  xprof[1] = kSideCLength1;
  yprof[1] = 0.;
  xprof[2] = xprof[1] + kSideCLength2*TMath::Cos(kSideCFoldAngle);
  yprof[2] = yprof[1] + kSideCLength2*TMath::Sin(kSideCFoldAngle);
  xprof[3] = xprof[2] - 2*kSideCHalfThick*TMath::Sin(kSideCFoldAngle);
  yprof[3] = yprof[2] + 2*kSideCHalfThick*TMath::Cos(kSideCFoldAngle);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      2*kSideCHalfThick , xprof[4], yprof[4]);
  xprof[5] = 0.;
  yprof[5] = 2*kSideCHalfThick;

  sideCLowerFace->DefinePolygon(kNPointsLowerFace, xprof, yprof);
  sideCLowerFace->DefineSection(0,-kSideCHalfWide);
  sideCLowerFace->DefineSection(1, kSideCHalfWide);

  // The Cable Tray lateral face: a Xtru
  TGeoXtru *sideCLateralFace = new TGeoXtru(2);

  xprof[0] = 0.;
  yprof[0] = 0.;
  xprof[1] = kSideCLength1;
  yprof[1] = 0.;
  xprof[2] = xprof[1] + kSideCLength2*TMath::Cos(kSideCFoldAngle);
  yprof[2] = yprof[1] + kSideCLength2*TMath::Sin(kSideCFoldAngle);
  xprof[3] = xprof[2] - kSideCHeight2*TMath::Sin(kSideCFoldAngle);
  yprof[3] = yprof[2] + kSideCHeight2*TMath::Cos(kSideCFoldAngle);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      kSideCHeight2, xprof[4], yprof[4]);
  xprof[5] = kSideCSideLength1 + kSideCSideLength2;
  yprof[5] = kSideCHeight2;
  xprof[6] = xprof[5];
  yprof[6] = kSideCSideHeight;
  xprof[7] = kSideCSideLength1;
  yprof[7] = kSideCHeight1;
  xprof[8] = 0;
  yprof[8] = yprof[7];

  sideCLateralFace->DefinePolygon(kNPointsLateralFace, xprof, yprof);
  sideCLateralFace->DefineSection(0,-kSideCHalfThick);
  sideCLateralFace->DefineSection(1, kSideCHalfThick);

  // The lateral wings: four BBox's
  TGeoBBox *sideCLateralWingA = new TGeoBBox(kSideCWingAHalfLen,
					     kSideCHalfThick,
					     kSideCWingsHalfWide);

  TGeoBBox *sideCLateralWingB = new TGeoBBox(kSideCWingBHalfLen,
					     kSideCHalfThick,
					     kSideCWingsHalfWide);

  TGeoBBox *sideCLateralWingC = new TGeoBBox(kSideCHalfThick,    // With these
					     kSideCWingCHalfLen, // X,Y avoid
					     kSideCWingsHalfWide);//rotations

  TGeoBBox *sideCLateralWingD = new TGeoBBox(kSideCWingDHalfLen,
					     kSideCHalfThick,
					     kSideCWingsHalfWide);

  TGeoBBox *sideCLateralWingE = new TGeoBBox(kSideCWingEHalfLen,
					     kSideCHalfThick,
					     kSideCWingsHalfWide);

  // The connecting lower plate: a Xtru
  TGeoXtru *sideCLowerPlate =  new TGeoXtru(2);

  xprof[0] = 0.;
  yprof[0] = 0.;
  xprof[1] = kPlateHalfLen;
  yprof[1] = 0.;
  xprof[2] = xprof[1] + kPlateHalfLen*TMath::Cos(kSideCFoldAngle);
  yprof[2] = kPlateHalfLen*TMath::Sin(kSideCFoldAngle);
  xprof[3] = xprof[2] - kPlateThick*TMath::Sin(kSideCFoldAngle);
  yprof[3] = yprof[2] + kPlateThick*TMath::Cos(kSideCFoldAngle);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      kPlateThick, xprof[4], yprof[4]);
  xprof[5] = 0.;
  yprof[5] = kPlateThick;

  sideCLowerPlate->DefinePolygon(kNPointsPlate, xprof, yprof);
  Double_t zwide = kSideCHalfWide + 2*kSideCHalfThick;
  sideCLowerPlate->DefineSection(0,-zwide);
  sideCLowerPlate->DefineSection(1, zwide);

  // The connecting side plate: a Xtru
  TGeoXtru *sideCLateralPlate = new TGeoXtru(2);

  xprof[0] = 0.;
  yprof[0] = 0.;
  xprof[1] = kPlateHalfLen;
  yprof[1] = 0.;
  xprof[2] = xprof[1] + kPlateHalfLen*TMath::Cos(kSideCFoldAngle);
  yprof[2] = kPlateHalfLen*TMath::Sin(kSideCFoldAngle);
  xprof[3] = xprof[2] - kPlateHeight*TMath::Sin(kSideCFoldAngle);
  yprof[3] = yprof[2] + kPlateHeight*TMath::Cos(kSideCFoldAngle);
  InsidePoint(xprof[0], yprof[0], xprof[1], yprof[1], xprof[2], yprof[2],
	      kPlateHeight, xprof[4], yprof[4]); // Avoid small overlap
  xprof[5] = 0.;
  yprof[5] = kPlateHeight;

  sideCLateralPlate->DefinePolygon(kNPointsPlate, xprof, yprof);
  sideCLateralPlate->DefineSection(0,-kPlateThick/2);
  sideCLateralPlate->DefineSection(1, kPlateThick/2);

  // The bar fixing the cooling tubes: a Tube
  TGeoTube *coolBar = new TGeoTube(0., kBarCoolRmax, kSideCHalfWide);

  // The Cable Tray cover: a (complex) Xtru
  TGeoXtru *sideCCoverFace = new TGeoXtru(2);

  xprof[ 0] = sideCLateralFace->GetX(8);
  yprof[ 0] = sideCLateralFace->GetY(8);
  xprof[ 1] = sideCLateralFace->GetX(7);
  yprof[ 1] = sideCLateralFace->GetY(7);
  xprof[ 2] = sideCLateralFace->GetX(6);
  yprof[ 2] = sideCLateralFace->GetY(6);
  xprof[ 3] = sideCLateralFace->GetX(5);
  yprof[ 3] = sideCLateralFace->GetY(5);
  xprof[ 4] = sideCLateralFace->GetX(4);
  yprof[ 4] = sideCLateralFace->GetY(4);

  xloc = (kSideCLength1 + (kSideCSideLength1+kSideCSideLength2))/2;
  delta  = kSideCLength1 - (xloc + kSideCWingDHalfLen);
  xprof[ 5] = xprof[4]
	    + (delta + 2*kSideCWingEHalfLen)*TMath::Cos(kSideCFoldAngle);
  yprof[ 5] = yprof[4]
	    + (delta + 2*kSideCWingEHalfLen)*TMath::Sin(kSideCFoldAngle);

  xprof[ 6] = xprof[5] - 2*kSideCHalfThick*TMath::Sin(kSideCFoldAngle);
  yprof[ 6] = yprof[5] + 2*kSideCHalfThick*TMath::Cos(kSideCFoldAngle);
  InsidePoint(xprof[3], yprof[3], xprof[4], yprof[4], xprof[5], yprof[5],
	      2*kSideCHalfThick, xprof[7], yprof[7]);
  InsidePoint(xprof[2], yprof[2], xprof[3], yprof[3], xprof[4], yprof[4],
	      2*kSideCHalfThick, xprof[8], yprof[8]);
  xprof[ 9] = xprof[2] + 2*kSideCHalfThick;
  yprof[ 9] = yprof[2] + 2*kSideCHalfThick;
  xprof[10] = xprof[1];
  yprof[10] = yprof[1] + 2*kSideCHalfThick;
  xprof[11] = xprof[0];
  yprof[11] = yprof[0] + 2*kSideCHalfThick;

  sideCCoverFace->DefinePolygon(kNPointsCoverFace, xprof, yprof);
  zloc = kSideCHalfWide + 2*kSideCHalfThick + 2*kSideCWingsHalfWide;
  sideCCoverFace->DefineSection(0,-zloc);
  sideCCoverFace->DefineSection(1, zloc);


  // We have all shapes: now create the real volumes
  TGeoMedium *medAl      = mgr->GetMedium("ITS_ALUMINUM$");

  TGeoVolume *traySideCLowerFace  = new TGeoVolume("ITSsuppTraySideCLower",
						   sideCLowerFace, medAl);

  traySideCLowerFace->SetVisibility(kTRUE);
  traySideCLowerFace->SetLineColor(6); // Purple
  traySideCLowerFace->SetLineWidth(1);
  traySideCLowerFace->SetFillColor(traySideCLowerFace->GetLineColor());
  traySideCLowerFace->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCLateralFace  = new TGeoVolume("ITSsuppTraySideCLateral",
						     sideCLateralFace, medAl);

  traySideCLateralFace->SetVisibility(kTRUE);
  traySideCLateralFace->SetLineColor(6); // Purple
  traySideCLateralFace->SetLineWidth(1);
  traySideCLateralFace->SetFillColor(traySideCLateralFace->GetLineColor());
  traySideCLateralFace->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCLateralWingA =
    new TGeoVolume("ITSsuppTraySideCLateralWingA", sideCLateralWingA,  medAl);

  traySideCLateralWingA->SetVisibility(kTRUE);
  traySideCLateralWingA->SetLineColor(6); // Purple
  traySideCLateralWingA->SetLineWidth(1);
  traySideCLateralWingA->SetFillColor(traySideCLateralWingA->GetLineColor());
  traySideCLateralWingA->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCLateralWingB =
    new TGeoVolume("ITSsuppTraySideCLateralWingB", sideCLateralWingB,  medAl);

  traySideCLateralWingB->SetVisibility(kTRUE);
  traySideCLateralWingB->SetLineColor(6); // Purple
  traySideCLateralWingB->SetLineWidth(1);
  traySideCLateralWingB->SetFillColor(traySideCLateralWingB->GetLineColor());
  traySideCLateralWingB->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCLateralWingC =
    new TGeoVolume("ITSsuppTraySideCLateralWingC", sideCLateralWingC,  medAl);

  traySideCLateralWingC->SetVisibility(kTRUE);
  traySideCLateralWingC->SetLineColor(6); // Purple
  traySideCLateralWingC->SetLineWidth(1);
  traySideCLateralWingC->SetFillColor(traySideCLateralWingC->GetLineColor());
  traySideCLateralWingC->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCLateralWingD =
    new TGeoVolume("ITSsuppTraySideCLateralWingD", sideCLateralWingD,  medAl);

  traySideCLateralWingD->SetVisibility(kTRUE);
  traySideCLateralWingD->SetLineColor(6); // Purple
  traySideCLateralWingD->SetLineWidth(1);
  traySideCLateralWingD->SetFillColor(traySideCLateralWingD->GetLineColor());
  traySideCLateralWingD->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCLateralWingE =
    new TGeoVolume("ITSsuppTraySideCLateralWingE", sideCLateralWingE,  medAl);

  traySideCLateralWingE->SetVisibility(kTRUE);
  traySideCLateralWingE->SetLineColor(6); // Purple
  traySideCLateralWingE->SetLineWidth(1);
  traySideCLateralWingE->SetFillColor(traySideCLateralWingE->GetLineColor());
  traySideCLateralWingE->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCLowerPlate =
    new TGeoVolume("ITSsuppTraySideCLowerPlate", sideCLowerPlate,  medAl);

  traySideCLowerPlate->SetVisibility(kTRUE);
  traySideCLowerPlate->SetLineColor(6); // Purple
  traySideCLowerPlate->SetLineWidth(1);
  traySideCLowerPlate->SetFillColor(traySideCLowerPlate->GetLineColor());
  traySideCLowerPlate->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCLateralPlate =
    new TGeoVolume("ITSsuppTraySideCLateralPlate", sideCLateralPlate,  medAl);

  traySideCLateralPlate->SetVisibility(kTRUE);
  traySideCLateralPlate->SetLineColor(6); // Purple
  traySideCLateralPlate->SetLineWidth(1);
  traySideCLateralPlate->SetFillColor(traySideCLateralPlate->GetLineColor());
  traySideCLateralPlate->SetFillStyle(4000); // 0% transparent

  TGeoVolume *traySideCCoverFace =
    new TGeoVolume("ITSsuppTraySideCCoverFace", sideCCoverFace,  medAl);

  traySideCCoverFace->SetVisibility(kTRUE);
  traySideCCoverFace->SetLineColor(6); // Purple
  traySideCCoverFace->SetLineWidth(1);
  traySideCCoverFace->SetFillColor(traySideCCoverFace->GetLineColor());
  traySideCCoverFace->SetFillStyle(4000); // 0% transparent

  TGeoVolume *coolingTubeBar = new TGeoVolume("ITSsuppTraySideCCoolBar",
					      coolBar, medAl);

  coolingTubeBar->SetVisibility(kTRUE);
  coolingTubeBar->SetLineColor(6); // Purple
  coolingTubeBar->SetLineWidth(1);
  coolingTubeBar->SetFillColor(coolingTubeBar->GetLineColor());
  coolingTubeBar->SetFillStyle(4000); // 0% transparent


  // Now build up the tray
  cableTrayC->AddNode(traySideCLowerFace,1,0);

  zloc = kSideCHalfWide + kSideCHalfThick;
  cableTrayC->AddNode(traySideCLateralFace,1,
			    new TGeoTranslation(0., 0., zloc) );
  cableTrayC->AddNode(traySideCLateralFace,2,
			    new TGeoTranslation(0., 0.,-zloc) );

  xloc = kSideCWingAHalfLen;
  yloc = kSideCHeight1 - kSideCHalfThick;
  zloc = kSideCHalfWide + 2*kSideCHalfThick + kSideCWingsHalfWide;
  cableTrayC->AddNode(traySideCLateralWingA,1,
			    new TGeoTranslation(xloc, yloc, zloc) );
  cableTrayC->AddNode(traySideCLateralWingA,2,
			    new TGeoTranslation(xloc, yloc,-zloc) );

  xloc = kSideCSideLength1 + kSideCSideLength2/2;
  yloc = Yfrom2Points(kSideCSideLength1,kSideCHeight1,
		      kSideCSideLength1+kSideCSideLength2,kSideCSideHeight,
		      xloc) - kSideCHalfThick -0.0012; // Avoid small overlap
  zloc = kSideCHalfWide + 2*kSideCHalfThick + kSideCWingsHalfWide;
  alpharot = (-(kSideCHeight1 - kSideCSideHeight)/kSideCSideLength2 )*
		TMath::RadToDeg();
  cableTrayC->AddNode(traySideCLateralWingB,1,
			    new TGeoCombiTrans(xloc, yloc, zloc,
					new TGeoRotation("",alpharot,0,0) ) );
  cableTrayC->AddNode(traySideCLateralWingB,2,
			    new TGeoCombiTrans(xloc, yloc,-zloc,
					new TGeoRotation("",alpharot,0,0) ) );

  xloc = kSideCSideLength1 + kSideCSideLength2 - kSideCHalfThick;
  yloc = kSideCSideHeight - kSideCWingCHalfLen;
  zloc = kSideCHalfWide + 2*kSideCHalfThick + kSideCWingsHalfWide;
  cableTrayC->AddNode(traySideCLateralWingC,1,
			    new TGeoTranslation(xloc, yloc, zloc) );
  cableTrayC->AddNode(traySideCLateralWingC,2,
			    new TGeoTranslation(xloc, yloc,-zloc) );

  xloc = (kSideCLength1 + (kSideCSideLength1+kSideCSideLength2))/2;
  yloc = kSideCHeight2 - kSideCHalfThick;
  zloc = kSideCHalfWide + 2*kSideCHalfThick + kSideCWingsHalfWide;
  cableTrayC->AddNode(traySideCLateralWingD,1,
			    new TGeoTranslation(xloc, yloc, zloc) );
  cableTrayC->AddNode(traySideCLateralWingD,2,
			    new TGeoTranslation(xloc, yloc,-zloc) );

  delta = kSideCLength1 - (xloc + kSideCWingDHalfLen);
  xloc = kSideCLength1 + delta + kSideCWingEHalfLen;
  yloc = (xloc - kSideCLength1)*TMath::Tan(kSideCFoldAngle) +
	  kSideCHeight2*TMath::Cos(kSideCFoldAngle) - kSideCHalfThick;
  zloc = kSideCHalfWide + 2*kSideCHalfThick + kSideCWingsHalfWide;
  alpharot = kSideCFoldAngle*TMath::RadToDeg();
  cableTrayC->AddNode(traySideCLateralWingE,1,
			    new TGeoCombiTrans(xloc, yloc, zloc,
					new TGeoRotation("",alpharot,0,0) ) );
  cableTrayC->AddNode(traySideCLateralWingE,2,
			    new TGeoCombiTrans(xloc, yloc,-zloc,
					new TGeoRotation("",alpharot,0,0) ) );

  xloc = kSideCLength1 - kPlateHalfLen;
  yloc = -kPlateThick -0.0025; // Avoid small overlap
  cableTrayC->AddNode(traySideCLowerPlate,1,
			    new TGeoTranslation(xloc, yloc, 0.) );

  xloc = kSideCLength1 - kPlateHalfLen;
  yloc = -kPlateThick;
  zloc = kSideCHalfWide + 2*kSideCHalfThick + kPlateThick/2;
  cableTrayC->AddNode(traySideCLateralPlate,1,
			    new TGeoTranslation(xloc, yloc, zloc) );
  cableTrayC->AddNode(traySideCLateralPlate,2,
			    new TGeoTranslation(xloc, yloc,-zloc) );

  for (Int_t jc = 0; jc <kNumBarCool; jc++) {
    xloc = kXShiftBarCool[jc];
    yloc = kYShiftBarCool[jc];
    cableTrayC->AddNode(coolingTubeBar,jc+1,
			      new TGeoTranslation(xloc, yloc, 0.) );
  }

  cableTrayC->AddNode(traySideCCoverFace,1,0);


  // Finally return what we made up

  return cableTrayC;
}

//______________________________________________________________________
void AliITSv11GeometrySupport::ITSTPCSupports(TGeoVolume *moth,
					const TGeoManager *mgr){
//
// Creates the elements suspending the ITS to the TPC and other fixed
// elements used to hook the rails (0872/C and its daughters)
//
//         moth : the TGeoVolume owing the volume structure
//         mgr  : the GeoManager (default gGeoManager)
// Output:
//
// Return:
//
// Created:      28 Oct 2010  Mario Sitta
//
// Technical data are taken from AutoCAD drawings, L.Simonetti technical
// drawings and other (oral) information given by F.Tosello
//

  // Dimensions and positions of the half ring C2/C3 (0872/C/04)
  const Double_t kRingCZPos           =   733.000*fgkmm;

  const Double_t kRingCThick          =    12.000*fgkmm;
  const Double_t kRingCRmin           =   565.000*fgkmm;
  const Double_t kRingCRmax           =   592.000*fgkmm;
  const Double_t kRingCHeight         =   560.000*fgkmm;
  const Double_t kRingCXToInsert      =   515.000*fgkmm;
  const Double_t kRingCYToInsert      =   113.000*fgkmm;

  const Int_t kNumberOfRingPoints     =    23; // N.points to approximate arc

  // Dimensions of the forward upper hook (0872/C/09)
  const Double_t kForwUpHookThick     =    20.000*fgkmm;
  const Double_t kForwUpHookRext      =   590.000*fgkmm;
  const Double_t kForwUpHookRint      =    20.000*fgkmm;
  const Double_t kForwUpHookHiTot     =    89.000*fgkmm;
  const Double_t kForwUpHookHiInt     =    59.000*fgkmm;
  const Double_t kForwUpHookWide      =    96.000*fgkmm;
  const Double_t kForwUpHookHalfBase  =    25.000*fgkmm;
  const Double_t kForwUpHookBaseCut   =    10.000*fgkmm;
  const Double_t kForwUpHookHoleWide  =    25.000*fgkmm;
  const Double_t kForwUpHookHoleHi    =    22.500*fgkmm;
  const Double_t kForwUpHookHoleBase  =     5.000*fgkmm;
  const Double_t kForwUpHookHoleR5    =     5.000*fgkmm;
  const Double_t kForwUpHookHoleY     =     8.000*fgkmm;
  const Double_t kForwUpHookHollowHi  =    35.000*fgkmm;
  const Double_t kForwUpHookHollowWide=     5.000*fgkmm;

  const Int_t kNumberOfForwUpHookPts  =    11;
  const Int_t kNumbOfForwUpHookHolePts=     5;

  // Dimensions of the forward lower hook (0872/C/08)
  const Double_t kForwLwHookThick     =    20.000*fgkmm;
  const Double_t kForwLwHookRext      =   590.000*fgkmm;
  const Double_t kForwLwHookRint      =    20.000*fgkmm;
  const Double_t kForwLwHookHiTot     =    88.500*fgkmm;
  const Double_t kForwLwHookWide      =    96.000*fgkmm;
  const Double_t kForwLwHookHalfBase  =    25.000*fgkmm;
  const Double_t kForwLwHookBaseCut   =    10.000*fgkmm;
  const Double_t kForwLwHookYToHollow =     3.500*fgkmm;
  const Double_t kForwLwHookHoleR     =     7.500*fgkmm;
  const Double_t kForwLwHookHoleIntHi =    35.000*fgkmm;
  const Double_t kForwLwHookHoleYPos  =    13.500*fgkmm;
  const Double_t kForwLwHookHollowHi  =    62.000*fgkmm;
  const Double_t kForwLwHookHollowWide=     5.000*fgkmm;

  const Int_t kNumberOfForwLwHookPts  =    11;
  const Int_t kNumbOfForwLwHookHolePts=     7;

  // Dimensions of the rear upper hook (0872/C/10)
  const Double_t kRearUpHookThick     =    15.000*fgkmm;
  const Double_t kRearUpHookRext      =   590.000*fgkmm;
  const Double_t kRearUpHookRint      =    20.000*fgkmm;
  const Double_t kRearUpHookHiTot     =    53.500*fgkmm;
  const Double_t kRearUpHookHiInt     =    23.500*fgkmm;
  const Double_t kRearUpHookWide      =    96.000*fgkmm;
  const Double_t kRearUpHookHalfBase  =    25.000*fgkmm;
  const Double_t kRearUpHookHoleWide  =    25.000*fgkmm;
  const Double_t kRearUpHookHoleHi    =    22.500*fgkmm;
  const Double_t kRearUpHookHoleBase  =     5.000*fgkmm;
  const Double_t kRearUpHookHoleR5    =     5.000*fgkmm;
  const Double_t kRearUpHookHoleY     =     8.000*fgkmm;

  const Int_t kNumberOfRearUpHookPts  =    10;
  const Int_t kNumbOfRearUpHookHolePts=     5;

  // Dimensions of the forward lower hook (0872/C/11)
  const Double_t kRearLwHookThick     =    20.000*fgkmm;
  const Double_t kRearLwHookRext      =   590.000*fgkmm;
  const Double_t kRearLwHookHiTot     =    30.000*fgkmm;
  const Double_t kRearLwHookWide      =    96.000*fgkmm;

  const Int_t kNumberOfRearLwHookPts  =     3;

  // Dimensions of the rear lower brackets (0872/C/16)
  const Double_t kRearLwBracketThick  =    15.000*fgkmm;
  const Double_t kRearLwBracketHi1    =    42.000*fgkmm;
  const Double_t kRearLwBracketHi2    =    12.000*fgkmm;
  const Double_t kRearLwBracketWide1  =    34.000*fgkmm;
  const Double_t kRearLwBracketWide2  =    10.000*fgkmm;
//  const Double_t kRearLwBracketR5     =     5.000*fgkmm

  // Dimensions of the forward webcam supports (0872/C/V/01-03-04)
  const Double_t kForwWebSStirrDep    =    20.000*fgkmm;
  const Double_t kForwWebSStirrLen1   =    15.000*fgkmm;
  const Double_t kForwWebSStirrLen2   =    55.000*fgkmm;
  const Double_t kForwWebSStirrLen3   =    10.000*fgkmm;
  const Double_t kForwWebSStirrWide1  =    45.000*fgkmm;
  const Double_t kForwWebSStirrWide2  =    38.000*fgkmm;
  const Double_t kForwWebSStirrWide3  =    23.000*fgkmm;
  const Double_t kForwWebTStirrThick  =     5.000*fgkmm;
  const Double_t kForwWebTStirrWide1  =    30.000*fgkmm;
  const Double_t kForwWebTStirrWide2  =    10.000*fgkmm;
  const Double_t kForwWebTStirrTotLen3=    58.500*fgkmm;
  const Double_t kForwWebTStirrTotLen4=    36.000*fgkmm;
  const Double_t kForwWebTStirrLen1   =    10.000*fgkmm;

  // Dimensions of the forward and rear webcam clamps (0872/C/V/02)
  const Double_t kFRWebClampThick     =    10.000*fgkmm;
  const Double_t kFRWebClampExtWide   =    30.000*fgkmm;
  const Double_t kFRWebClampIntWide   =    18.000*fgkmm;
  const Double_t kFRWebClampExtHi     =    22.000*fgkmm;
  const Double_t kFRWebClampIntHi     =    17.000*fgkmm;

  // Dimensions of the webcam itself
  const Double_t kWebcamLength        =    35.000*fgkmm;//ESTIMATED!!!

  // Dimensions of the rear upper webcam supports (0872/C/V/05-06)
  const Double_t kRearUpWebStirrWide  =    71.000*fgkmm;
  const Double_t kRearUpWebStirrDep   =    15.000*fgkmm;
  const Double_t kRearUpWebStirrThick =     5.000*fgkmm;
  const Double_t kRearUpWebStirrH1    =    27.000*fgkmm;
  const Double_t kRearUpWebStirrH2    =    32.000*fgkmm;
  const Double_t kRearUpWebBarLen     =   130.000*fgkmm;
  const Double_t kRearUpWebBarHi      =    20.000*fgkmm;
  const Double_t kRearUpWebBarThick   =     5.000*fgkmm;

  // Dimensions of the upper wheel slides (0872/C/Z/00-01-02)
  const Double_t kUpperSlideTotHeight =    93.500*fgkmm;
  const Double_t kUpperSlideBlockHi   =    62.500*fgkmm;
  const Double_t kUpperSlideWidth     =    36.000*fgkmm;
  const Double_t kUpperSlideTotDepth  =    51.000*fgkmm;
  const Double_t kUpperSlideIntDepth  =    36.000*fgkmm;
  const Double_t kUpperSlideStubHi    =    15.000*fgkmm;
  const Double_t kUpperSlideStubDep   =     8.000*fgkmm;
  const Double_t kUpperSlideWheelHi   =    18.500*fgkmm;
  const Double_t kUpperSlideHoleRout  =    11.000*fgkmm;
  const Double_t kUpperSlideHoleRint1 =     9.000*fgkmm;
  const Double_t kUpperSlideHoleRint2 =    11.500*fgkmm;
  const Double_t kUpperSlideHoleH1    =     7.000*fgkmm;
  const Double_t kUpperSlideHoleH2    =    46.000*fgkmm;
  const Double_t kUpperSlideHoleH3    =     1.100*fgkmm;
  const Double_t kUpperSlideHoleXPos  =    20.000*fgkmm;
  const Double_t kUpperSlidePinRmin   =     4.000*fgkmm;
  const Double_t kUpperSlidePinRmax   =     6.000*fgkmm;
  const Double_t kUpperSlidePinH1     =     7.000*fgkmm;
  const Double_t kUpperSlidePinH2     =    46.000*fgkmm;
  const Double_t kUpperSlidePinH3     =    25.500*fgkmm;

  // Dimensions of the lower wheel slides (0872/C/W/00-01-02-03)
  const Double_t kLowerSlideTotHeight =    80.000*fgkmm;
  const Double_t kLowerSlideBlockHi   =    28.000*fgkmm;
  const Double_t kLowerSlideWidth     =    36.000*fgkmm;
  const Double_t kLowerSlideTotDepth  =    60.000*fgkmm;
  const Double_t kLowerSlideHoleRout  =     9.500*fgkmm;
  const Double_t kLowerSlideHoleRint  =     4.700*fgkmm;
  const Double_t kLowerSlideHoleH1    =    12.000*fgkmm;
  const Double_t kLowerSlideNoseBase  =    40.000*fgkmm;
  const Double_t kLowerSlideNoseBasHi =     6.000*fgkmm;//Computed
  const Double_t kLowerSlideNoseUpWid =    25.000*fgkmm;
  const Double_t kLowerSlideNoseDepth =    10.000*fgkmm;
  const Double_t kLowerSlidePinRmin   =     3.000*fgkmm;
  const Double_t kLowerSlidePinRmax   =     4.000*fgkmm;
  const Double_t kLowerSlidePinH1     =    12.000*fgkmm;
  const Double_t kLowerSlidePinH2     =    10.000*fgkmm;


  // Local variables
  Double_t xprof[2*kNumberOfRingPoints],yprof[2*kNumberOfRingPoints];
  Double_t xpos, ypos, zpos, alpha;
  

  // First create all needed shapes

  // The Supporting Ring (0872/C/04): a really complex Xtru
  // to approximate the arc with a polyline
  TGeoXtru *ringC2C3 = new TGeoXtru(2);

  for (Int_t j=0; j<11; j++) { // The external arc
    xprof[j] = kRingCRmax*SinD(90*j/10);
    yprof[j] = kRingCRmax*CosD(90*j/10);
  }

  xprof[11] = kRingCRmin;
  yprof[11] = yprof[10];

  alpha = TMath::ASin(kRingCYToInsert/kRingCRmin); // Now the insert
  xprof[12] = kRingCRmin*TMath::Cos(alpha/2);
  yprof[12] = kRingCRmin*TMath::Sin(alpha/2);
  xprof[13] = kRingCRmin*TMath::Cos(alpha);
  yprof[13] = kRingCRmin*TMath::Sin(alpha);

  xprof[14] = kRingCXToInsert;
  yprof[14] = yprof[13];

  alpha = TMath::ACos(kRingCXToInsert/kRingCRmin); // The insert ending angle
  xprof[15] = kRingCRmin*TMath::Cos(alpha);
  yprof[15] = kRingCRmin*TMath::Sin(alpha);

  for (Int_t j=7; j>1; j--) { // The internal arc
    xprof[23-j] = kRingCRmin*SinD(90*j/10);
    yprof[23-j] = kRingCRmin*CosD(90*j/10);
  }

  alpha = TMath::ASin(kRingCHeight/kRingCRmin);    // The angle till the notch
  xprof[22] = kRingCRmin*TMath::Cos(alpha);
  yprof[22] = kRingCRmin*TMath::Sin(alpha);

  xprof[23] = xprof[0];
  yprof[23] = yprof[22];

  // We did the right side, now reflex on the left side
  for (Int_t jp = 0; jp < 22; jp++) {
    xprof[24+jp] = -xprof[23-1-jp];
    yprof[24+jp] =  yprof[23-1-jp];
  }

  // wow! now the actual Xtru
  ringC2C3->DefinePolygon(2*kNumberOfRingPoints, xprof, yprof);
  ringC2C3->DefineSection(0, 0);
  ringC2C3->DefineSection(1, kRingCThick);

  // The Forward Upper Hook (0872/C/09): a Composite Shape made of
  // a really complex Xtru to approximate the arc with a polyline,
  // another Xtru for the hole, and a BBox for the hollow
  // The main body
  TGeoXtru *forwUpHookMainBody = new TGeoXtru(2);
  forwUpHookMainBody->SetName("ITSforwUpHookMainBody");

  xprof[ 0] = kForwUpHookHalfBase - kForwUpHookBaseCut;
  yprof[ 0] = kForwUpHookRext - kForwUpHookHiTot;
  xprof[ 1] = kForwUpHookHalfBase;
  yprof[ 1] = yprof[0] + kForwUpHookBaseCut;
  xprof[ 2] = xprof[1];
  yprof[ 2] = yprof[0] + (kForwUpHookHiInt - kForwUpHookRint);
  for (Int_t j=1; j<6; j++) {
    xprof[2+j] = xprof[2] + kForwUpHookRint*(1 - CosD(90*j/5));
    yprof[2+j] = yprof[2] + kForwUpHookRint*SinD(90*j/5);
  }
  xprof[ 8] = kForwUpHookWide/2;
  yprof[ 8] = yprof[7];
  xprof[ 9] = xprof[8];
  alpha = TMath::ASin(0.5*kForwUpHookWide/kForwUpHookRext);
  yprof[ 9] = kForwUpHookRext*TMath::Cos(alpha);
  xprof[10] = kForwUpHookRext*TMath::Sin(alpha/2);
  yprof[10] = kForwUpHookRext*TMath::Cos(alpha/2);
  xprof[11] = 0;
  yprof[11] = kForwUpHookRext;

  // We did the right side, now reflex on the left side
  for (Int_t jp = 0; jp < kNumberOfForwUpHookPts; jp++) {
    xprof[12+jp] = -xprof[10-jp];
    yprof[12+jp] =  yprof[10-jp];
  }

  // Now the actual Xtru
  forwUpHookMainBody->DefinePolygon(2*kNumberOfForwUpHookPts+1, xprof, yprof);
  forwUpHookMainBody->DefineSection(0, 0);
  forwUpHookMainBody->DefineSection(1, kForwUpHookThick);

  // The hole
  TGeoXtru *forwUpHookHole = new TGeoXtru(2);
  forwUpHookHole->SetName("ITSforwUpHookHole");

  xprof[0] = kForwUpHookHoleBase/2;
  yprof[0] = forwUpHookMainBody->GetY(0) + kForwUpHookHoleY;
  xprof[1] = kForwUpHookHoleWide/2;
  yprof[1] = yprof[0] + (xprof[1] - xprof[0]); // Go at 45deg
  xprof[2] = xprof[1];
  yprof[2] = yprof[0] + kForwUpHookHoleHi - kForwUpHookHoleR5;
  xprof[3] = xprof[2] - kForwUpHookHoleR5*(1 - CosD(45));
  yprof[3] = yprof[2] + kForwUpHookHoleR5*SinD(45);
  xprof[4] = xprof[2] - kForwUpHookHoleR5;
  yprof[4] = yprof[0] + kForwUpHookHoleHi;

  // We did the right side, now reflex on the left side
  for (Int_t jp = 0; jp < kNumbOfForwUpHookHolePts; jp++) {
    xprof[5+jp] = -xprof[4-jp];
    yprof[5+jp] =  yprof[4-jp];
  }

  // Now the actual Xtru
  forwUpHookHole->DefinePolygon(2*kNumbOfForwUpHookHolePts, xprof, yprof);
  forwUpHookHole->DefineSection(0, -0.1);
  forwUpHookHole->DefineSection(1, kForwUpHookThick+0.1);

  // The hollow
  TGeoBBox *forwUpHookHollow = new TGeoBBox(2.1 *kForwUpHookHalfBase,
					    0.55*kForwUpHookHollowHi,
					    0.55*kForwUpHookHollowWide);
  forwUpHookHollow->SetName("ITSforwUpHookHollow");

  TGeoTranslation *forwUpHookHollPos = new TGeoTranslation(0.,
		      forwUpHookMainBody->GetY(0) + 0.5*kForwUpHookHollowHi,
		      forwUpHookMainBody->GetZ(1) - 0.5*kForwUpHookHollowWide);
  forwUpHookHollPos->SetName("ITSforwUpHookHollPos");
  forwUpHookHollPos->RegisterYourself();

  // Finally the actual shape: a CompositeShape
  TGeoCompositeShape *forwUpHookShape = new TGeoCompositeShape("ITSforwUpHookMainBody-ITSforwUpHookHole-ITSforwUpHookHollow:ITSforwUpHookHollPos");

  // The Forward Lower Hook (0872/C/08): a Composite Shape made of
  // a really complex Xtru to approximate the arc with a polyline,
  // another Xtru for the hole, and a BBox for the hollow
  // The main body
  TGeoXtru *forwLwHookMainBody = new TGeoXtru(2);
  forwLwHookMainBody->SetName("ITSforwLwHookMainBody");

  xprof[ 0] = kForwLwHookHalfBase - kForwLwHookBaseCut;
  yprof[ 0] = kForwLwHookRext - kForwLwHookHiTot;
  xprof[ 1] = kForwLwHookHalfBase;
  yprof[ 1] = yprof[0] + kForwLwHookBaseCut;
  xprof[ 2] = xprof[1];
  yprof[ 2] = yprof[0] + (kForwLwHookHollowHi - kForwLwHookYToHollow
			  - kForwLwHookRint);
  for (Int_t j=1; j<6; j++) {
    xprof[2+j] = xprof[2] + kForwLwHookRint*(1 - CosD(90*j/5));
    yprof[2+j] = yprof[2] + kForwLwHookRint*SinD(90*j/5);
  }
  xprof[ 8] = kForwLwHookWide/2;
  yprof[ 8] = yprof[7];
  xprof[ 9] = xprof[8];
  alpha = TMath::ASin(0.5*kForwLwHookWide/kForwLwHookRext);
  yprof[ 9] = kForwLwHookRext*TMath::Cos(alpha);
  xprof[10] = kForwLwHookRext*TMath::Sin(alpha/2);
  yprof[10] = kForwLwHookRext*TMath::Cos(alpha/2);
  xprof[11] = 0;
  yprof[11] = kForwLwHookRext;

  // We did the right side, now reflex on the left side
  for (Int_t jp = 0; jp < kNumberOfForwLwHookPts; jp++) {
    xprof[12+jp] = -xprof[10-jp];
    yprof[12+jp] =  yprof[10-jp];
  }

  // Now the actual Xtru
  forwLwHookMainBody->DefinePolygon(2*kNumberOfForwLwHookPts+1, xprof, yprof);
  forwLwHookMainBody->DefineSection(0, 0);
  forwLwHookMainBody->DefineSection(1, kForwLwHookThick);

  // The hole
  TGeoXtru *forwLwHookHole = new TGeoXtru(2);
  forwLwHookHole->SetName("ITSforwLwHookHole");

  xprof[0] = 0;
  yprof[0] = forwLwHookMainBody->GetY(0) + kForwLwHookHoleYPos
	   - kForwLwHookHoleR;
  for (Int_t j=1; j<3; j++) {
    xprof[0+j] = xprof[0] + kForwLwHookHoleR*SinD(90*j/3);
    yprof[0+j] = yprof[0] + kForwLwHookHoleR*(1 - CosD(90*j/3));
  }
  xprof[3] = xprof[0] + kForwLwHookHoleR;
  yprof[3] = yprof[0] + kForwLwHookHoleR;
  xprof[4] = xprof[3];
  yprof[4] = yprof[3] + kForwLwHookHoleIntHi;
  for (Int_t j=1; j<3; j++) {
    xprof[4+j] = xprof[4] - kForwLwHookHoleR*(1 - CosD(90*j/3));
    yprof[4+j] = yprof[4] + kForwLwHookHoleR*SinD(90*j/3);
  }
  xprof[7] = xprof[0];
  yprof[7] = yprof[4] + kForwLwHookHoleR;

  // We did the right side, now reflex on the left side
  for (Int_t jp = 0; jp < kNumbOfForwLwHookHolePts-1; jp++) {
    xprof[8+jp] = -xprof[6-jp];
    yprof[8+jp] =  yprof[6-jp];
  }

  // Now the actual Xtru
  forwLwHookHole->DefinePolygon(2*kNumbOfForwLwHookHolePts, xprof, yprof);
  forwLwHookHole->DefineSection(0, -0.1);
  forwLwHookHole->DefineSection(1, kForwLwHookThick+0.1);

  // The hollow
  TGeoBBox *forwLwHookHollow = new TGeoBBox(2.1 *kForwLwHookHalfBase,
					    0.55*kForwLwHookHollowHi,
					    0.55*kForwLwHookHollowWide);
  forwLwHookHollow->SetName("ITSforwLwHookHollow");

  TGeoTranslation *forwLwHookHollPos = new TGeoTranslation(0.,
		      forwLwHookMainBody->GetY(0) + 0.5*kForwLwHookHollowHi,
		      forwLwHookMainBody->GetZ(1) - 0.5*kForwLwHookHollowWide);
  forwLwHookHollPos->SetName("ITSforwLwHookHollPos");
  forwLwHookHollPos->RegisterYourself();

  // Finally the actual shape: a CompositeShape
  TGeoCompositeShape *forwLwHookShape = new TGeoCompositeShape("ITSforwLwHookMainBody-ITSforwLwHookHole-ITSforwLwHookHollow:ITSforwLwHookHollPos");

  // The Rear Upper Hook (0872/C/10): a Composite Shape made of
  // a really complex Xtru to approximate the arc with a polyline,
  // and another Xtru for the hole
  // The main body
  TGeoXtru *rearUpHookMainBody = new TGeoXtru(2);
  rearUpHookMainBody->SetName("ITSrearUpHookMainBody");

  xprof[0] = kRearUpHookHalfBase;
  yprof[0] = kRearUpHookRext - kRearUpHookHiTot;
  xprof[1] = xprof[0];
  yprof[1] = yprof[0] + (kRearUpHookHiInt - kRearUpHookRint); 
  for (Int_t j=1; j<6; j++) {
    xprof[1+j] = xprof[1] + kRearUpHookRint*(1 - CosD(90*j/5));
    yprof[1+j] = yprof[1] + kRearUpHookRint*SinD(90*j/5);
  }
  xprof[ 7] = kRearUpHookWide/2;
  yprof[ 7] = yprof[5];
  xprof[ 8] = xprof[7];
  alpha = TMath::ASin(0.5*kRearUpHookWide/kRearUpHookRext);
  yprof[ 8] = kRearUpHookRext*TMath::Cos(alpha);
  xprof[ 9] = kRearUpHookRext*TMath::Sin(alpha/2);
  yprof[ 9] = kRearUpHookRext*TMath::Cos(alpha/2);
  xprof[10] = 0;
  yprof[10] = kRearUpHookRext;

  // We did the right side, now reflex on the left side
  for (Int_t jp = 0; jp < kNumberOfRearUpHookPts; jp++) {
    xprof[11+jp] = -xprof[9-jp];
    yprof[11+jp] =  yprof[9-jp];
  }

  // Now the actual Xtru
  rearUpHookMainBody->DefinePolygon(2*kNumberOfRearUpHookPts+1, xprof, yprof);
  rearUpHookMainBody->DefineSection(0, 0);
  rearUpHookMainBody->DefineSection(1, kRearUpHookThick);

  // The hole
  TGeoXtru *rearUpHookHole = new TGeoXtru(2);
  rearUpHookHole->SetName("ITSrearUpHookHole");

  xprof[0] = kRearUpHookHoleBase/2;
  yprof[0] = rearUpHookMainBody->GetY(0) + kRearUpHookHoleY;
  xprof[1] = kRearUpHookHoleWide/2;
  yprof[1] = yprof[0] + (xprof[1] - xprof[0]); // Go at 45deg
  xprof[2] = xprof[1];
  yprof[2] = yprof[0] + kRearUpHookHoleHi - kRearUpHookHoleR5;
  xprof[3] = xprof[2] - kRearUpHookHoleR5*(1 - CosD(45));
  yprof[3] = yprof[2] + kRearUpHookHoleR5*SinD(45);
  xprof[4] = xprof[2] - kRearUpHookHoleR5;
  yprof[4] = yprof[0] + kRearUpHookHoleHi;

  // We did the right side, now reflex on the left side
  for (Int_t jp = 0; jp < kNumbOfRearUpHookHolePts; jp++) {
    xprof[5+jp] = -xprof[4-jp];
    yprof[5+jp] =  yprof[4-jp];
  }

  // Now the actual Xtru
  rearUpHookHole->DefinePolygon(2*kNumbOfRearUpHookHolePts, xprof, yprof);
  rearUpHookHole->DefineSection(0, -0.1);
  rearUpHookHole->DefineSection(1, kRearUpHookThick+0.1);

  // Finally the actual shape: a CompositeShape
  TGeoCompositeShape *rearUpHookShape = new TGeoCompositeShape("ITSrearUpHookMainBody-ITSrearUpHookHole");

  // The Rear Lower Hook (0872/C/11): a Xtru
  TGeoXtru *rearLwHookShape = new TGeoXtru(2);
  rearLwHookShape->SetName("ITSrearLwHookShape");

  xprof[0] = kRearLwHookWide/2;
  yprof[0] = kRearLwHookRext - kRearLwHookHiTot;
  xprof[1] = xprof[0];
  alpha = TMath::ASin(0.5*kRearLwHookWide/kRearLwHookRext);
  yprof[1] = kRearLwHookRext*TMath::Cos(alpha);
  xprof[2] = kRearLwHookRext*TMath::Sin(alpha/2);
  yprof[2] = kRearLwHookRext*TMath::Cos(alpha/2);
  xprof[3] = 0;
  yprof[3] = kRearLwHookRext;

  // We did the right side, now reflex on the left side
  for (Int_t jp = 0; jp < kNumberOfRearLwHookPts; jp++) {
    xprof[4+jp] = -xprof[2-jp];
    yprof[4+jp] =  yprof[2-jp];
  }

  // Now the actual Xtru
  rearLwHookShape->DefinePolygon(2*kNumberOfRearLwHookPts+1, xprof, yprof);
  rearLwHookShape->DefineSection(0, 0);
  rearLwHookShape->DefineSection(1, kRearLwHookThick);

  // The Rear Lower Bracket (0872/C/16): a Xtru
  TGeoXtru *rearLwBrackShape = new TGeoXtru(2);
  rearLwBrackShape->SetName("ITSrearLwBrackShape");

  xprof[0] = 0;
  yprof[0] = 0;
  xprof[1] = xprof[0] + kRearLwBracketWide1 - kRearLwBracketWide2;
  yprof[1] = yprof[0];
  xprof[2] = xprof[1];
  yprof[2] = yprof[0] + kRearLwBracketHi2;
  xprof[3] = xprof[2] - kRearLwBracketWide1;
  yprof[3] = yprof[2];
  xprof[4] = xprof[3];
  yprof[4] = yprof[3] - kRearLwBracketHi1;
  xprof[5] = xprof[0];
  yprof[5] = yprof[4];

  rearLwBrackShape->DefinePolygon(6, xprof, yprof);
  rearLwBrackShape->DefineSection(0,-kRearLwBracketThick/2);
  rearLwBrackShape->DefineSection(1, kRearLwBracketThick/2);

  // The Forward S-shaped Stirrup for the webcam (0872/C/V/01): a Xtru
  TGeoXtru *forwWebSStirrSh = new TGeoXtru(2);

  xprof[0] = 0;
  yprof[0] = 0;
  xprof[1] = xprof[0] + kForwWebSStirrLen1;
  yprof[1] = yprof[0];
  xprof[2] = xprof[1];
  yprof[2] = yprof[1] + kForwWebSStirrWide1;
  xprof[3] = xprof[0] - kForwWebSStirrLen2 + kForwWebSStirrLen3;
  yprof[3] = yprof[2];
  xprof[4] = xprof[3];
  yprof[4] = yprof[3] + kForwWebSStirrWide3;
  xprof[5] = xprof[4] - kForwWebSStirrLen3;
  yprof[5] = yprof[4];
  xprof[6] = xprof[5];
  yprof[6] = yprof[0] + kForwWebSStirrWide2;
  xprof[7] = xprof[0];
  yprof[7] = yprof[6];

  forwWebSStirrSh->DefinePolygon(8, xprof, yprof);
  forwWebSStirrSh->DefineSection(0,-kForwWebSStirrDep/2);
  forwWebSStirrSh->DefineSection(1, kForwWebSStirrDep/2);

  // The Forward T-shaped Stirrups for the webcam (0872/C/V/03-04): two Xtru
  TGeoXtru *forwWebTStirr3Sh = new TGeoXtru(2);

  xprof[0] = -kForwWebTStirrWide2/2;
  yprof[0] = 0;
  xprof[1] = -kForwWebTStirrWide1/2;
  yprof[1] = yprof[0];
  xprof[2] = xprof[1];
  yprof[2] = yprof[1] - kForwWebTStirrLen1;
  xprof[3] =-xprof[2];
  yprof[3] = yprof[2];
  xprof[4] = xprof[3];
  yprof[4] = yprof[1];
  xprof[5] =-xprof[0];
  yprof[5] = yprof[4];
  xprof[6] = xprof[5];
  yprof[6] = kForwWebTStirrTotLen3 - kForwWebTStirrLen1;
  xprof[7] = xprof[0];
  yprof[7] = yprof[6];

  forwWebTStirr3Sh->DefinePolygon(8, xprof, yprof);
  forwWebTStirr3Sh->DefineSection(0, 0);
  forwWebTStirr3Sh->DefineSection(1, kForwWebTStirrThick);

  TGeoXtru *forwWebTStirr4Sh = new TGeoXtru(2);

  yprof[6] = kForwWebTStirrTotLen4 - kForwWebTStirrLen1;
  yprof[7] = yprof[6];

  forwWebTStirr4Sh->DefinePolygon(8, xprof, yprof);
  forwWebTStirr4Sh->DefineSection(0, 0);
  forwWebTStirr4Sh->DefineSection(1, kForwWebTStirrThick);

  // The Forward and Rear clamp for the webcam (0872/C/V/02): a Xtru
  TGeoXtru *frWebClampSh = new TGeoXtru(2);

  xprof[0] = kFRWebClampIntWide/2;
  yprof[0] = kFRWebClampIntHi;
  xprof[1] = xprof[0];
  yprof[1] = 0;
  xprof[2] = kFRWebClampExtWide/2;
  yprof[2] = yprof[1];
  xprof[3] = xprof[2];
  yprof[3] = kFRWebClampExtHi;
  for (Int_t jp = 0; jp < 4; jp++) {
    xprof[4+jp] = -xprof[3-jp];
    yprof[4+jp] =  yprof[3-jp];
  }

  frWebClampSh->DefinePolygon(8, xprof, yprof);
  frWebClampSh->DefineSection(0,-kFRWebClampThick/2);
  frWebClampSh->DefineSection(1, kFRWebClampThick/2);

  // The Rear Upper Stirrup for the webcam (0872/C/V/05): a Xtru
  TGeoXtru *upWebStirrSh = new TGeoXtru(2);

  xprof[0] = 0;
  yprof[0] = 0;
  xprof[1] = xprof[0] - (kRearUpWebStirrWide - 2*kRearUpWebStirrThick);
  yprof[1] = yprof[0];
  xprof[2] = xprof[1];
  yprof[2] = yprof[1] + (kRearUpWebStirrH1 - kRearUpWebStirrThick);
  xprof[3] = xprof[2] - kRearUpWebStirrThick;
  yprof[3] = yprof[2];
  xprof[4] = xprof[3];
  yprof[4] = yprof[3] - kRearUpWebStirrH1;
  xprof[5] = xprof[4] + kRearUpWebStirrWide;
  yprof[5] = yprof[4];
  xprof[6] = xprof[5];
  yprof[6] = yprof[5] + kRearUpWebStirrH2;
  xprof[7] = xprof[0];
  yprof[7] = yprof[6];

  upWebStirrSh->DefinePolygon(8, xprof, yprof);
  upWebStirrSh->DefineSection(0,-kRearUpWebStirrDep/2);
  upWebStirrSh->DefineSection(1, kRearUpWebStirrDep/2);

  // The Rear Upper Bar for the webcam (0872/C/V/06): a BBox
  TGeoBBox *upRearWebBarSh = new TGeoBBox(kRearUpWebBarLen/2,
					  kRearUpWebBarHi/2,
					  kRearUpWebBarThick/2);

  // The Webcam: a BBox
  TGeoBBox *webcamShape = new TGeoBBox(kFRWebClampIntWide/2,
				       kWebcamLength/2,
				       kFRWebClampIntHi/2);

  // The Upper Wheel Slide (0872/C/Z/00-01-02)
  // A mother volume of air (to avoid assembly) contains the Alluminum block
  // (a Composite Shape: a Xtru and a Pcon for the hole) and the Steel pin
  // (a Pcon) (The wheels are approximated as part of the block itself)
  // The Air mother volume
  TGeoXtru *upSlideAirSh = new TGeoXtru(2);
  upSlideAirSh->SetName("ITSupperSlideAirShape");

  xprof[0] = 0;
  yprof[0] = 0;
  xprof[1] = xprof[0];
  yprof[1] = kUpperSlideBlockHi + kUpperSlideStubHi - kUpperSlideWheelHi;
  xprof[2] = xprof[1] - kUpperSlideIntDepth;
  yprof[2] = yprof[1];
  xprof[3] = xprof[2];
  yprof[3] = yprof[2] - kUpperSlideTotHeight;
  xprof[4] = xprof[3] + kUpperSlideTotDepth;
  yprof[4] = yprof[3];
  xprof[5] = xprof[4];
  yprof[5] = yprof[0];

  upSlideAirSh->DefinePolygon(6, xprof, yprof);
  upSlideAirSh->DefineSection(0,-kUpperSlideWidth/2);
  upSlideAirSh->DefineSection(1, kUpperSlideWidth/2);

  // The (filled) Aluminum block: a Xtru
  TGeoXtru *upSlideAluSh = new TGeoXtru(2);
  upSlideAluSh->SetName("ITSupperSlideAluShape");

  xprof[0] = upSlideAirSh->GetX(0);
  yprof[0] = upSlideAirSh->GetY(0);
  xprof[1] = upSlideAirSh->GetX(1);
  yprof[1] = upSlideAirSh->GetY(1);
  xprof[2] = xprof[1] - kUpperSlideStubDep;
  yprof[2] = yprof[1];
  xprof[3] = xprof[2];
  yprof[3] = yprof[2] - kUpperSlideStubHi;
  xprof[4] = upSlideAirSh->GetX(2);
  yprof[4] = yprof[3];
  xprof[5] = xprof[4];
  yprof[5] = yprof[4] - kUpperSlideBlockHi;
  xprof[6] = upSlideAirSh->GetX(5);
  yprof[6] = yprof[5];
  xprof[7] = xprof[6];
  yprof[7] = yprof[0];

  upSlideAluSh->DefinePolygon(8, xprof, yprof);
  upSlideAluSh->DefineSection(0, upSlideAirSh->GetZ(0));
  upSlideAluSh->DefineSection(1, upSlideAirSh->GetZ(1));

  // The cylindrical hole in the block; a Pcon
  TGeoPcon *upSlideHoleSh = new TGeoPcon(0, 360, 10);
  upSlideHoleSh->SetName("ITSupperSlideHoleShape");

  zpos = upSlideAluSh->GetY(5);
  upSlideHoleSh->DefineSection(0, zpos-0.1, 0, kUpperSlideHoleRout);
  zpos += (kUpperSlideBlockHi - kUpperSlideHoleH3 - kUpperSlideHoleH2
	- 2*kUpperSlideHoleH1);
  upSlideHoleSh->DefineSection(1, zpos, 0, kUpperSlideHoleRout);
  upSlideHoleSh->DefineSection(2, zpos, 0, kUpperSlideHoleRint2);
  zpos += kUpperSlideHoleH3;
  upSlideHoleSh->DefineSection(3, zpos, 0, kUpperSlideHoleRint2);
  upSlideHoleSh->DefineSection(4, zpos, 0, kUpperSlideHoleRout);
  zpos += kUpperSlideHoleH1;
  upSlideHoleSh->DefineSection(5, zpos, 0, kUpperSlideHoleRout);
  upSlideHoleSh->DefineSection(6, zpos, 0, kUpperSlideHoleRint1);
  zpos += kUpperSlideHoleH2;
  upSlideHoleSh->DefineSection(7, zpos, 0, kUpperSlideHoleRint1);
  upSlideHoleSh->DefineSection(8, zpos, 0, kUpperSlideHoleRout);
  zpos += kUpperSlideHoleH1;
  upSlideHoleSh->DefineSection(9, zpos+0.1, 0, kUpperSlideHoleRout);

  TGeoCombiTrans *upSlideHolePos = new TGeoCombiTrans(-kUpperSlideHoleXPos,0,0,
				   new TGeoRotation("",0,-90,0) );
  upSlideHolePos->SetName("ITSupperSlideHolePos");
  upSlideHolePos->RegisterYourself();

  // The actual block: a CompositeShape
  TGeoCompositeShape *upSlideBlockSh = new TGeoCompositeShape("ITSupperSlideAluShape-ITSupperSlideHoleShape:ITSupperSlideHolePos");

  // The Steel pin in the block; a Pcon
  TGeoPcon *upSlidePinSh = new TGeoPcon(0, 360, 6);
  upSlidePinSh->SetName("ITSupperSlidePinShape");

  zpos = upSlideAluSh->GetY(5) - (kUpperSlidePinH1 + kUpperSlidePinH2
       + kUpperSlidePinH3 - kUpperSlideBlockHi);
  upSlidePinSh->DefineSection(0, zpos, 0, kUpperSlidePinRmin);
  zpos += kUpperSlidePinH3;
  upSlidePinSh->DefineSection(1, zpos, 0, kUpperSlidePinRmin);
  upSlidePinSh->DefineSection(2, zpos, 0, kUpperSlidePinRmax);
  zpos += kUpperSlidePinH2;
  upSlidePinSh->DefineSection(3, zpos, 0, kUpperSlidePinRmax);
  upSlidePinSh->DefineSection(4, zpos, 0, kUpperSlidePinRmin);
  zpos += kUpperSlidePinH1;
  upSlidePinSh->DefineSection(5, zpos, 0, kUpperSlidePinRmin);

  // The Lower Wheel Slide (0872/C/W/00-01-02-03)
  // A mother volume of air (to avoid assembly) contains the Alluminum block
  // (a Composite Shape: a Xtru and a Pcon for the hole), the Alluminum nose
  // (a Xtru) and the Steel pin (a Pcon)
  // (The wheels are approximated as part of the block itself)
  // The Air mother volume
  TGeoXtru *lwSlideAirSh = new TGeoXtru(2);
  lwSlideAirSh->SetName("ITSlowerSlideAirShape");

  xprof[0] = 0;
  yprof[0] = 0;
  xprof[1] = xprof[0] + kLowerSlideTotDepth/2 - kLowerSlideNoseBase/2;
  yprof[1] = yprof[0];
  xprof[2] = xprof[1];
  yprof[2] = yprof[1] - (kLowerSlideBlockHi + kLowerSlidePinH2);
  xprof[3] = xprof[2] - kLowerSlideTotDepth;
  yprof[3] = yprof[2];
  xprof[4] = xprof[3];
  yprof[4] = yprof[3] + kLowerSlidePinH2 + kLowerSlideTotHeight;
  xprof[5] = xprof[0];
  yprof[5] = yprof[4];

  lwSlideAirSh->DefinePolygon(6, xprof, yprof);
  lwSlideAirSh->DefineSection(0,-kLowerSlideWidth/2);
  lwSlideAirSh->DefineSection(1, kLowerSlideWidth/2);

  // The (filled) Aluminum block: a Xtru
  TGeoXtru *lwSlideAluSh = new TGeoXtru(2);
  lwSlideAluSh->SetName("ITSlowerSlideAluShape");

  xprof[0] = lwSlideAirSh->GetX(0);
  yprof[0] = lwSlideAirSh->GetY(0);
  xprof[1] = lwSlideAirSh->GetX(1);
  yprof[1] = lwSlideAirSh->GetY(1);
  xprof[2] = xprof[1];
  yprof[2] = yprof[1] - kLowerSlideBlockHi;
  xprof[3] = lwSlideAirSh->GetX(3);
  yprof[3] = yprof[2];
  xprof[4] = xprof[3];
  yprof[4] = yprof[3] + kLowerSlideBlockHi;
  xprof[5] = xprof[4] + kLowerSlideTotDepth/2;
  yprof[5] = yprof[4];
  xprof[6] = xprof[5];
  yprof[6] = lwSlideAirSh->GetY(4);
  xprof[7] = xprof[0];
  yprof[7] = yprof[6];

  lwSlideAluSh->DefinePolygon(8, xprof, yprof);
  lwSlideAluSh->DefineSection(0, lwSlideAirSh->GetZ(0));
  lwSlideAluSh->DefineSection(1, lwSlideAirSh->GetZ(1));

  // The cylindrical hole in the block; a Pcon
  TGeoPcon *lwSlideHoleSh = new TGeoPcon(0, 360, 4);
  lwSlideHoleSh->SetName("ITSlowerSlideHoleShape");

  zpos = lwSlideAluSh->GetY(2);
  lwSlideHoleSh->DefineSection(0, zpos-0.1, 0, kLowerSlideHoleRout);
  zpos += kLowerSlideHoleH1;
  lwSlideHoleSh->DefineSection(1, zpos, 0, kLowerSlideHoleRout);
  lwSlideHoleSh->DefineSection(2, zpos, 0, kLowerSlideHoleRint);
  zpos = lwSlideAluSh->GetY(4);
  lwSlideHoleSh->DefineSection(3, zpos, 0, kLowerSlideHoleRint);

  TGeoCombiTrans *lwSlideHolePos = new TGeoCombiTrans(lwSlideAluSh->GetX(5),
						      0, 0,
				   new TGeoRotation("",0,-90,0) );
  lwSlideHolePos->SetName("ITSlowerSlideHolePos");
  lwSlideHolePos->RegisterYourself();

  // The actual block: a CompositeShape
  TGeoCompositeShape *lwSlideBlockSh = new TGeoCompositeShape("ITSlowerSlideAluShape-ITSlowerSlideHoleShape:ITSlowerSlideHolePos");

  // The Aluminum nose: a Xtru
  TGeoXtru *lwSlideNoseSh = new TGeoXtru(2);
  lwSlideNoseSh->SetName("ITSlowerSlideNoseShape");

  xprof[0] = lwSlideAluSh->GetX(5);
  yprof[0] = lwSlideAluSh->GetY(5);
  xprof[1] = xprof[0] - kLowerSlideNoseBase/2;
  yprof[1] = yprof[0];
  xprof[2] = xprof[1];
  yprof[2] = yprof[1] + kLowerSlideNoseBasHi;
  xprof[3] = lwSlideAluSh->GetX(0) - kLowerSlideNoseUpWid;
  yprof[3] = lwSlideAluSh->GetY(6);
  xprof[4] = xprof[0];
  yprof[4] = yprof[3];

  lwSlideNoseSh->DefinePolygon(5, xprof, yprof);
  lwSlideNoseSh->DefineSection(0,-kLowerSlideNoseDepth/2);
  lwSlideNoseSh->DefineSection(1, kLowerSlideNoseDepth/2);

  // The Steel pin in the block; a Pcon
  TGeoPcon *lwSlidePinSh = new TGeoPcon(0, 360, 4);
  lwSlidePinSh->SetName("ITSlowerSlidePinShape");

  zpos = lwSlideAirSh->GetY(2);
  lwSlidePinSh->DefineSection(0, zpos, 0, kLowerSlidePinRmax);
  zpos += kLowerSlidePinH2;
  lwSlidePinSh->DefineSection(1, zpos, 0, kLowerSlidePinRmax);
  lwSlidePinSh->DefineSection(2, zpos, 0, kLowerSlidePinRmin);
  zpos += kLowerSlidePinH1;
  lwSlidePinSh->DefineSection(3, zpos, 0, kLowerSlidePinRmin);


  // We have all shapes: now create the real volumes
  TGeoMedium *medAlcoa   = mgr->GetMedium("ITS_ALUMINUM$"); // To code!!!!!!
  TGeoMedium *medHokotol = mgr->GetMedium("ITS_HOKOTOL$");
  TGeoMedium *medAnticor = mgr->GetMedium("ITS_ANTICORODAL$");
  TGeoMedium *medAisi    = mgr->GetMedium("ITS_AISI304L$");
  TGeoMedium *medAir     = mgr->GetMedium("ITS_AIR$");
  TGeoMedium *medPlexy   = mgr->GetMedium("ITS_PLEXYGLAS$");
  TGeoMedium *medPVC     = mgr->GetMedium("ITS_PVC$");

  TGeoVolume *suppRingC2C3  = new TGeoVolume("ITSTPCsupportRingC2C3",
					     ringC2C3, medAlcoa);

  suppRingC2C3->SetVisibility(kTRUE);
  suppRingC2C3->SetLineColor(6); // Purple
  suppRingC2C3->SetLineWidth(1);
  suppRingC2C3->SetFillColor(suppRingC2C3->GetLineColor());
  suppRingC2C3->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwUpHook  = new TGeoVolume("ITSTPCsupportForwUpHook",
					   forwUpHookShape, medHokotol);

  forwUpHook->SetVisibility(kTRUE);
  forwUpHook->SetLineColor(6); // Purple
  forwUpHook->SetLineWidth(1);
  forwUpHook->SetFillColor(forwUpHook->GetLineColor());
  forwUpHook->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwLwHook  = new TGeoVolume("ITSTPCsupportForwLwHook",
					   forwLwHookShape, medHokotol);

  forwLwHook->SetVisibility(kTRUE);
  forwLwHook->SetLineColor(6); // Purple
  forwLwHook->SetLineWidth(1);
  forwLwHook->SetFillColor(forwLwHook->GetLineColor());
  forwLwHook->SetFillStyle(4000); // 0% transparent

  TGeoVolume *rearUpHook  = new TGeoVolume("ITSTPCsupportRearUpHook",
					   rearUpHookShape, medHokotol);

  rearUpHook->SetVisibility(kTRUE);
  rearUpHook->SetLineColor(6); // Purple
  rearUpHook->SetLineWidth(1);
  rearUpHook->SetFillColor(rearUpHook->GetLineColor());
  rearUpHook->SetFillStyle(4000); // 0% transparent

  TGeoVolume *rearLwHook  = new TGeoVolume("ITSTPCsupportRearLwHook",
					   rearLwHookShape, medAnticor);

  rearLwHook->SetVisibility(kTRUE);
  rearLwHook->SetLineColor(6); // Purple
  rearLwHook->SetLineWidth(1);
  rearLwHook->SetFillColor(rearLwHook->GetLineColor());
  rearLwHook->SetFillStyle(4000); // 0% transparent

  TGeoVolume *rearLwBrack  = new TGeoVolume("ITSTPCsupportRearLwBracket",
					    rearLwBrackShape, medAnticor);

  rearLwBrack->SetVisibility(kTRUE);
  rearLwBrack->SetLineColor(6); // Purple
  rearLwBrack->SetLineWidth(1);
  rearLwBrack->SetFillColor(rearLwBrack->GetLineColor());
  rearLwBrack->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwWebSStirrup  = new TGeoVolume("ITSTPCsupportForwWebSStirrup",
						forwWebSStirrSh, medAnticor);

  forwWebSStirrup->SetVisibility(kTRUE);
  forwWebSStirrup->SetLineColor(6); // Purple
  forwWebSStirrup->SetLineWidth(1);
  forwWebSStirrup->SetFillColor(forwWebSStirrup->GetLineColor());
  forwWebSStirrup->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwWebTStirr3  = new TGeoVolume("ITSTPCsupportForwWebTStirrup3",
					       forwWebTStirr3Sh, medAnticor);

  forwWebTStirr3->SetVisibility(kTRUE);
  forwWebTStirr3->SetLineColor(6); // Purple
  forwWebTStirr3->SetLineWidth(1);
  forwWebTStirr3->SetFillColor(forwWebTStirr3->GetLineColor());
  forwWebTStirr3->SetFillStyle(4000); // 0% transparent

  TGeoVolume *forwWebTStirr4  = new TGeoVolume("ITSTPCsupportForwWebTStirrup4",
					       forwWebTStirr4Sh, medAnticor);

  forwWebTStirr4->SetVisibility(kTRUE);
  forwWebTStirr4->SetLineColor(6); // Purple
  forwWebTStirr4->SetLineWidth(1);
  forwWebTStirr4->SetFillColor(forwWebTStirr4->GetLineColor());
  forwWebTStirr4->SetFillStyle(4000); // 0% transparent

  TGeoVolume *frWebClamp  = new TGeoVolume("ITSTPCsupportForwRearWebClamp",
					   frWebClampSh, medPlexy);

  frWebClamp->SetVisibility(kTRUE);
  frWebClamp->SetLineColor(kAzure);
  frWebClamp->SetLineWidth(1);
  frWebClamp->SetFillColor(frWebClamp->GetLineColor());
  frWebClamp->SetFillStyle(4000); // 0% transparent

  TGeoVolume *upWebStirrup  = new TGeoVolume("ITSTPCsupportUpperWebStirrup",
					     upWebStirrSh, medAnticor);

  upWebStirrup->SetVisibility(kTRUE);
  upWebStirrup->SetLineColor(6); // Purple
  upWebStirrup->SetLineWidth(1);
  upWebStirrup->SetFillColor(upWebStirrup->GetLineColor());
  upWebStirrup->SetFillStyle(4000); // 0% transparent

  TGeoVolume *upRearWebBar  = new TGeoVolume("ITSTPCsupportUpperRearWebBar",
					     upRearWebBarSh, medPlexy);

  upRearWebBar->SetVisibility(kTRUE);
  upRearWebBar->SetLineColor(kAzure);
  upRearWebBar->SetLineWidth(1);
  upRearWebBar->SetFillColor(upRearWebBar->GetLineColor());
  upRearWebBar->SetFillStyle(4000); // 0% transparent

  TGeoVolume *webCam  = new TGeoVolume("ITSTPCsupportWebcam",
				       webcamShape, medPVC);

  webCam->SetVisibility(kTRUE);
  webCam->SetLineColor(kBlack);
  webCam->SetLineWidth(1);
  webCam->SetFillColor(webCam->GetLineColor());
  webCam->SetFillStyle(4000); // 0% transparent

  TGeoVolume *upSlideVol  = new TGeoVolume("ITSTPCsupportUpperSlide",
					   upSlideAirSh, medAir);

  upSlideVol->SetVisibility(kFALSE);

  TGeoVolume *upSlideBlock  = new TGeoVolume("ITSTPCsupportUpperSlideBlock",
					     upSlideBlockSh, medAnticor);

  upSlideBlock->SetVisibility(kTRUE);
  upSlideBlock->SetLineColor(6); // Purple
  upSlideBlock->SetLineWidth(1);
  upSlideBlock->SetFillColor(upSlideBlock->GetLineColor());
  upSlideBlock->SetFillStyle(4000); // 0% transparent

  TGeoVolume *upSlidePin  = new TGeoVolume("ITSTPCsupportUpperSlidePin",
					   upSlidePinSh, medAisi);

  upSlidePin->SetVisibility(kTRUE);
  upSlidePin->SetLineColor(kGray);
  upSlidePin->SetLineWidth(1);
  upSlidePin->SetFillColor(upSlidePin->GetLineColor());
  upSlidePin->SetFillStyle(4000); // 0% transparent

  TGeoVolume *lwSlideVol  = new TGeoVolume("ITSTPCsupportLowerSlide",
					   lwSlideAirSh, medAir);

  lwSlideVol->SetVisibility(kFALSE);

  TGeoVolume *lwSlideBlock  = new TGeoVolume("ITSTPCsupportLowerSlideBlock",
					     lwSlideBlockSh, medAnticor);

  lwSlideBlock->SetVisibility(kTRUE);
  lwSlideBlock->SetLineColor(6); // Purple
  lwSlideBlock->SetLineWidth(1);
  lwSlideBlock->SetFillColor(lwSlideBlock->GetLineColor());
  lwSlideBlock->SetFillStyle(4000); // 0% transparent

  TGeoVolume *lwSlideNose  = new TGeoVolume("ITSTPCsupportLowerSlideNose",
					    lwSlideNoseSh, medAnticor);

  lwSlideNose->SetVisibility(kTRUE);
  lwSlideNose->SetLineColor(6); // Purple
  lwSlideNose->SetLineWidth(1);
  lwSlideNose->SetFillColor(lwSlideNose->GetLineColor());
  lwSlideNose->SetFillStyle(4000); // 0% transparent

  TGeoVolume *lwSlidePin  = new TGeoVolume("ITSTPCsupportLowerSlidePin",
					   lwSlidePinSh, medAisi);

  lwSlidePin->SetVisibility(kTRUE);
  lwSlidePin->SetLineColor(kGray);
  lwSlidePin->SetLineWidth(1);
  lwSlidePin->SetFillColor(lwSlidePin->GetLineColor());
  lwSlidePin->SetFillStyle(4000); // 0% transparent


  // Build up the wheel slides
  upSlideVol->AddNode(upSlideBlock,1,0);
  upSlideVol->AddNode(upSlidePin,  1,
		      new TGeoCombiTrans(-kUpperSlideHoleXPos, 0, 0,
					 new TGeoRotation("",0,-90,0) ) );

  lwSlideVol->AddNode(lwSlideBlock,1,0);
  lwSlideVol->AddNode(lwSlideNose ,1,0);
  lwSlideVol->AddNode(lwSlidePin,  1,
		      new TGeoCombiTrans(lwSlideAluSh->GetX(5), 0, 0,
					 new TGeoRotation("",0,-90,0) ) );


  // Finally put everything in the mother volume
  moth->AddNode(suppRingC2C3,1,
		new TGeoTranslation(0, 0, kRingCZPos) );
  moth->AddNode(suppRingC2C3,2,
		new TGeoCombiTrans( 0, 0,-kRingCZPos,
				   new TGeoRotation("",0.,180.,0.) ) );
  moth->AddNode(suppRingC2C3,3,
		new TGeoCombiTrans( 0, 0, kRingCZPos,
				   new TGeoRotation("",0.,0.,180.) ) );
  moth->AddNode(suppRingC2C3,4,
		new TGeoCombiTrans( 0, 0,-kRingCZPos,
				   new TGeoRotation("",0.,180.,180.) ) );

  zpos = kRingCZPos + kRingCThick;
  moth->AddNode(forwUpHook,1,
		new TGeoTranslation( 0, 0, zpos) );

  zpos = kRingCZPos + kRingCThick;
  moth->AddNode(forwLwHook,1,
		new TGeoCombiTrans( 0, 0, zpos,
				   new TGeoRotation("",0.,0.,180.) ) );

  zpos = kRingCZPos + kRingCThick + kRearUpHookThick;
  moth->AddNode(rearUpHook,1,
		new TGeoTranslation( 0, 0,-zpos) );

  zpos = kRingCZPos + kRingCThick + kRearLwHookThick;
  moth->AddNode(rearLwHook,1,
		new TGeoCombiTrans( 0, 0,-zpos,
				   new TGeoRotation("",0.,0.,180.) ) );

  xpos =  kRearLwHookWide/2 + kRearLwBracketThick/2;
  ypos = -kRingCHeight;
  moth->AddNode(rearLwBrack,1,
		new TGeoCombiTrans( xpos, ypos,-zpos,
				   new TGeoRotation("", 90.,-90.,-90.) ) );
  moth->AddNode(rearLwBrack,2,
		new TGeoCombiTrans(-xpos, ypos,-zpos,
				   new TGeoRotation("", 90.,-90.,-90.) ) );

  xpos = kForwUpHookWide/2;
  ypos = (forwUpHookMainBody->GetY(8) + forwUpHookMainBody->GetY(9))/2;
  zpos = kRingCZPos + kRingCThick;
  moth->AddNode(forwWebSStirrup,1,
		new TGeoCombiTrans( xpos, ypos, zpos,
				   new TGeoRotation("", 0., 90., 0.) ) );
  xpos = kForwLwHookWide/2;
  ypos = (forwLwHookMainBody->GetY(8) + forwLwHookMainBody->GetY(9))/2;
  moth->AddNode(forwWebSStirrup,2,
		new TGeoCombiTrans( xpos,-ypos, zpos,
				   new TGeoRotation("", 0., 90., 0.) ) );

  xpos = kForwUpHookWide/2
	+ (forwWebSStirrSh->GetX(4) + forwWebSStirrSh->GetX(5))/2;
  ypos = (forwUpHookMainBody->GetY(8) + forwUpHookMainBody->GetY(9))/2
	+  forwWebSStirrSh->GetZ(1) - forwWebTStirr3Sh->GetY(7);
  zpos += (forwWebSStirrSh->GetY(4) - forwWebSStirrSh->GetY(0));
  moth->AddNode(forwWebTStirr3,1,
		new TGeoTranslation( xpos, ypos, zpos) );

  ypos -= frWebClampSh->GetZ(1);
  moth->AddNode(frWebClamp,1,
		new TGeoCombiTrans( xpos, ypos, zpos+forwWebTStirr3Sh->GetZ(1),
				   new TGeoRotation("", 0., 90., 0.) ) );

  ypos -= webcamShape->GetDY()/2;
  moth->AddNode(webCam,1,
		new TGeoTranslation( xpos, ypos,
		     zpos+forwWebTStirr3Sh->GetZ(1)+webcamShape->GetDZ()) );

  xpos = kForwLwHookWide/2
	+ (forwWebSStirrSh->GetX(4) + forwWebSStirrSh->GetX(5))/2;
  ypos = (forwLwHookMainBody->GetY(8) + forwLwHookMainBody->GetY(9))/2
	+  forwWebSStirrSh->GetZ(1) - forwWebTStirr4Sh->GetY(7);
  moth->AddNode(forwWebTStirr4,1,
		new TGeoCombiTrans( xpos,-ypos, zpos,
				   new TGeoRotation("", 180., 0., 0.) ) );

  ypos -= frWebClampSh->GetZ(1);
  moth->AddNode(frWebClamp,2,
		new TGeoCombiTrans( xpos,-ypos, zpos+forwWebTStirr4Sh->GetZ(1),
				   new TGeoRotation("", 0., 90., 0.) ) );

  ypos -= webcamShape->GetDY()/2;
  moth->AddNode(webCam,2,
		new TGeoTranslation( xpos,-ypos,
		     zpos+forwWebTStirr4Sh->GetZ(1)+webcamShape->GetDZ()) );

  xpos = kRearUpHookWide/2 + kRearUpWebStirrDep/2;
  ypos = kRingCHeight;
  zpos = kRingCZPos + kRingCThick;
  moth->AddNode(upWebStirrup,1,
		new TGeoCombiTrans( xpos, ypos,-zpos,
				   new TGeoRotation("",-90.,-90., 90.) ) );
  moth->AddNode(upWebStirrup,2,
		new TGeoCombiTrans(-xpos, ypos,-zpos,
				   new TGeoRotation("",-90.,-90., 90.) ) );

  ypos = kRingCHeight + upWebStirrSh->GetY(2) - upRearWebBarSh->GetDY();
  zpos = kRingCZPos + kRingCThick + upWebStirrSh->GetX(3)
       - upRearWebBarSh->GetDZ();
  moth->AddNode(upRearWebBar,1,
		new TGeoTranslation( 0, ypos,-zpos) );

  zpos -= upRearWebBarSh->GetDZ();
  moth->AddNode(frWebClamp,3,
		new TGeoCombiTrans( 0, ypos,-zpos,
				   new TGeoRotation("", 0., 90., 0.) ) );

  ypos -= webcamShape->GetDY()/2;
  zpos -= webcamShape->GetDZ();
  moth->AddNode(webCam,3,
		new TGeoTranslation( 0, ypos,-zpos) );

  xpos = ringC2C3->GetX(14) + kUpperSlideWidth/2;
  ypos = ringC2C3->GetY(14);
  zpos = kRingCZPos + kRingCThick;
  moth->AddNode(upSlideVol,1,
		new TGeoCombiTrans( xpos, ypos, zpos,
				   new TGeoRotation("",-90.,-90., 90.) ) );
  moth->AddNode(upSlideVol,2,
		new TGeoCombiTrans(-xpos, ypos, zpos,
				   new TGeoRotation("",-90.,-90., 90.) ) );
  moth->AddNode(upSlideVol,3,
		new TGeoCombiTrans( xpos, ypos, -zpos,
				   new TGeoRotation("", 90.,-90.,-90.) ) );
  moth->AddNode(upSlideVol,4,
		new TGeoCombiTrans(-xpos, ypos, -zpos,
				   new TGeoRotation("", 90.,-90.,-90.) ) );

  moth->AddNode(lwSlideVol,1,
		new TGeoCombiTrans( xpos,-ypos,-zpos,
				   new TGeoRotation("",-90.,-90.,-90.) ) );
  moth->AddNode(lwSlideVol,2,
		new TGeoCombiTrans(-xpos,-ypos,-zpos,
				   new TGeoRotation("",-90.,-90.,-90.) ) );
  moth->AddNode(lwSlideVol,3,
		new TGeoCombiTrans( xpos,-ypos, zpos,
				   new TGeoRotation("", 90.,-90., 90.) ) );
  moth->AddNode(lwSlideVol,4,
		new TGeoCombiTrans(-xpos,-ypos, zpos,
				   new TGeoRotation("", 90.,-90., 90.) ) );


  return;
}

