
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
//  Flex class for ALICE MFT upgrade
//  This version uses TGeo
//  Authors:
//  F. Manso
// march 2016: main implementaion of the class
//-------------------------------------------------------------------------

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoTrd2.h"
#include "TGeoMatrix.h"
#include "TGeoBBox.h"
#include "TGeoTube.h"
#include "AliLog.h"
#include "AliMFTLadderSegmentation.h"
#include "AliMFTChipSegmentation.h"
#include "AliMFTFlex.h"
#include "AliMFTChip.h"
#include "AliMFTLadder.h"
#include "AliMFTGeometry.h"
#include "AliMFTPlane.h"
#include "TGeoCompositeShape.h"
#include "TGeoBoolNode.h"

ClassImp(AliMFTFlex)

AliMFTFlex::AliMFTFlex():
TNamed(),
fLadderSeg(NULL)
{
  // Constructor
}

AliMFTFlex::~AliMFTFlex() {
}

AliMFTFlex::AliMFTFlex(AliMFTLadderSegmentation *ladder):
TNamed(),
fLadderSeg(ladder)
{
  // Constructor
}


TGeoVolumeAssembly* AliMFTFlex::MakeFlex(Int_t nbsensors, Double_t length)
{
  // Informations from the technical report mft_flex_proto_5chip_v08_laz50p.docx on MFT twiki and private communications

  // For the naming
  AliMFTGeometry * mftGeom = AliMFTGeometry::Instance();
  Int_t idHalfMFT = mftGeom->GetHalfMFTID(fLadderSeg->GetUniqueID());
  Int_t idHalfDisk = mftGeom->GetHalfDiskID(fLadderSeg->GetUniqueID());
  Int_t idLadder = mftGeom->GetLadderID(fLadderSeg->GetUniqueID());

  // First a global pointer for the flex
  TGeoMedium *kMedAir = gGeoManager->GetMedium("MFT_Air$");
  TGeoVolumeAssembly*  flex  = new TGeoVolumeAssembly(Form("flex_%d_%d_%d",idHalfMFT,idHalfDisk,idLadder));

  // Defining one single layer for the strips and the AVDD and DVDD
  TGeoVolume* lines = Make_Lines(nbsensors,length-AliMFTGeometry::kClearance, AliMFTGeometry::kFlexHeight - AliMFTGeometry::kClearance, AliMFTGeometry::kAluThickness);

  // AGND and DGND layers
  TGeoVolume* agnd_dgnd = Make_AGND_DGND(length-AliMFTGeometry::kClearance, AliMFTGeometry::kFlexHeight-AliMFTGeometry::kClearance, AliMFTGeometry::kAluThickness);

  // The others layers
  TGeoVolume* kaptonlayer     = Make_Kapton(length, AliMFTGeometry::kFlexHeight, AliMFTGeometry::kKaptonThickness);
  TGeoVolume* varnishlayerIn  = Make_Varnish(length, AliMFTGeometry::kFlexHeight, AliMFTGeometry::kVarnishThickness,0);
  TGeoVolume* varnishlayerOut = Make_Varnish(length, AliMFTGeometry::kFlexHeight, AliMFTGeometry::kVarnishThickness,1);
    
  // Final flex building
  Double_t zvarnishIn = AliMFTGeometry::kKaptonThickness/2+AliMFTGeometry::kAluThickness+AliMFTGeometry::kVarnishThickness/2-AliMFTGeometry::kGlueThickness;
  Double_t zgnd = AliMFTGeometry::kKaptonThickness/2+AliMFTGeometry::kAluThickness/2-AliMFTGeometry::kGlueThickness;
  Double_t zkaptonlayer = -AliMFTGeometry::kGlueThickness;
  Double_t zlines = -AliMFTGeometry::kKaptonThickness/2-AliMFTGeometry::kAluThickness/2-AliMFTGeometry::kGlueThickness;
  Double_t zvarnishOut = -AliMFTGeometry::kKaptonThickness/2-AliMFTGeometry::kAluThickness-AliMFTGeometry::kVarnishThickness/2-AliMFTGeometry::kGlueThickness;

  //-----------------------------------------------------------------------------------------
  //-------------------------- Adding all layers of the FPC ----------------------------------
  //-----------------------------------------------------------------------------------------
  
  flex->AddNode(varnishlayerIn,  1,  new TGeoTranslation(0., 0., zvarnishIn));    // inside, in front of the cold plate
  flex->AddNode(agnd_dgnd,       1,  new TGeoTranslation(0., 0., zgnd));
  flex->AddNode(kaptonlayer,     1,  new TGeoTranslation(0., 0., zkaptonlayer));
  flex->AddNode(lines,           1,  new TGeoTranslation(0., 0., zlines));
  flex->AddNode(varnishlayerOut, 1,  new TGeoTranslation(0., 0., zvarnishOut));   // outside

  Make_ElectricComponents(flex, nbsensors, length, zvarnishOut);
  //-----------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------

  return flex;
}


void AliMFTFlex::Make_ElectricComponents(TGeoVolumeAssembly*  flex, Int_t nbsensors, Double_t length, Double_t zvarnish) 
{
  // Making and adding all the electric components
  TGeoVolumeAssembly *electric[200];

  // 2 components on the connector side
  Int_t total;

  TGeoRotation *rotation = new TGeoRotation ("rotation", 90., 0., 0.);
  TGeoRotation *rotationpi = new TGeoRotation ("rotationpi", 180., 0., 0.);
  TGeoCombiTrans *transformation0 = new TGeoCombiTrans(length/2 - 0.1, AliMFTGeometry::kFlexHeight/2-0.2, 
						      zvarnish-AliMFTGeometry::kVarnishThickness/2-AliMFTGeometry::kCapacitorDz/2, rotation);
  TGeoCombiTrans *transformation1 = new TGeoCombiTrans(length/2 - 0.1, AliMFTGeometry::kFlexHeight/2-0.6, 
						      zvarnish-AliMFTGeometry::kVarnishThickness/2-AliMFTGeometry::kCapacitorDz/2, rotation);

  for(Int_t id=0; id < 2; id++)electric[id] = Make_ElectricComponent(AliMFTGeometry::kCapacitorDy, AliMFTGeometry::kCapacitorDx, AliMFTGeometry::kCapacitorDz, id);
  flex->AddNode(electric[0], 1, transformation0);
  flex->AddNode(electric[1], 2, transformation1);
  total=2;

   // 2 lines of electric components along the FPC in the middle (4 per sensor)
  for(Int_t id=0; id < 4*nbsensors; id++)electric[id+total] = Make_ElectricComponent(AliMFTGeometry::kCapacitorDy, AliMFTGeometry::kCapacitorDx, 
										     AliMFTGeometry::kCapacitorDz, id+total);
  for(Int_t id=0; id < 2*nbsensors; id++){
    flex->AddNode(electric[id+total], id+1000, new TGeoTranslation(-length/2 + (id+0.5)*AliMFTGeometry::kSensorLength/2, AliMFTGeometry::kFlexHeight/2 - 0.35, 
								   zvarnish - AliMFTGeometry::kVarnishThickness/2 - AliMFTGeometry::kCapacitorDz/2));
    flex->AddNode(electric[id+total+2*nbsensors], id+2000, new TGeoTranslation(-length/2 + (id+0.5)*AliMFTGeometry::kSensorLength/2, 0., 
									       zvarnish - AliMFTGeometry::kVarnishThickness/2 - AliMFTGeometry::kCapacitorDz/2));
  }
  total=total+4*nbsensors;

  // ------- 3 components on the FPC side -------- 
  for(Int_t id=0; id < 3; id++)electric[id+total] = Make_ElectricComponent(AliMFTGeometry::kCapacitorDy, AliMFTGeometry::kCapacitorDx, 
									   AliMFTGeometry::kCapacitorDz, id+total);
  for(Int_t id=0 ; id < 3; id++){
    flex->AddNode(electric[id+total], id+3000, new TGeoTranslation(-length/2+AliMFTGeometry::kSensorLength+(id+1)*0.3-0.6, -AliMFTGeometry::kFlexHeight/2 + 0.2, 
								   zvarnish-AliMFTGeometry::kVarnishThickness/2 - AliMFTGeometry::kCapacitorDz/2));
  }
  total=total+3;

  
  /*
  // The connector of the FPC
  for(Int_t id=0; id < 74; id++)electric[id+total] = Make_ElectricComponent(AliMFTGeometry::kConnectorLength, AliMFTGeometry::kConnectorWidth, 
									    AliMFTGeometry::kConnectorThickness, id+total);
  for(Int_t id=0; id < 37; id++){
    flex->AddNode(electric[id+total], id+100, new TGeoTranslation(length/2+0.15-AliMFTGeometry::kConnectorOffset, id*0.04-AliMFTGeometry::kFlexHeight/2 + 0.1, 
								  zvarnish-AliMFTGeometry::kVarnishThickness/2-AliMFTGeometry::kCapacitorDz/2));
    flex->AddNode(electric[id+total+37], id+200, new TGeoTranslation(length/2-0.15-AliMFTGeometry::kConnectorOffset, id*0.04-AliMFTGeometry::kFlexHeight/2 + 0.1, 
								     zvarnish - AliMFTGeometry::kVarnishThickness/2 - AliMFTGeometry::kCapacitorDz/2));
  }
  total=total+74;
  */
  

  
  //-------------------------- New Connector ----------------------
  TGeoMedium *kMedAlu = gGeoManager->GetMedium("MFT_Alu$");
  TGeoMedium *kMedPeek = gGeoManager->GetMedium("MFT_PEEK$");

  TGeoBBox *connect = new TGeoBBox("connect", AliMFTGeometry::kConnectorLength/2, AliMFTGeometry::kConnectorWidth/2, AliMFTGeometry::kConnectorHeight/2);
  TGeoBBox *remov = new TGeoBBox("remov", AliMFTGeometry::kConnectorLength/2, AliMFTGeometry::kConnectorWidth/2+AliMFTGeometry::kEpsilon, 
				 AliMFTGeometry::kConnectorHeight/2+AliMFTGeometry::kEpsilon);

  TGeoTranslation    *t1= new TGeoTranslation ("t1", AliMFTGeometry::kConnectorThickness, 0., -0.01);
  TGeoSubtraction    *connecto = new TGeoSubtraction(connect, remov, NULL, t1);
  TGeoCompositeShape *connector = new TGeoCompositeShape("connector", connecto);
  TGeoVolume *connectord = new TGeoVolume("connectord", connector, kMedAlu);
  connectord->SetVisibility(kTRUE);
  connectord->SetLineColor(kRed);
  connectord->SetLineWidth(1);
  connectord->SetFillColor(connectord->GetLineColor());
  connectord->SetFillStyle(4000); // 0% transparent

  Double_t interspace = 0.1; // interspace inside the 2 ranges of connector pads
  Double_t step = 0.04;      // interspace between each pad inside the connector
  for(Int_t id=0; id < 37; id++){
    flex->AddNode(connectord,id+total,new TGeoTranslation(length/2+interspace/2+AliMFTGeometry::kConnectorLength/2-AliMFTGeometry::kConnectorOffset, 
							  id*step-AliMFTGeometry::kFlexHeight/2 + 0.1, zvarnish - AliMFTGeometry::kVarnishThickness/2 
							  - AliMFTGeometry::kConnectorHeight/2));
    TGeoCombiTrans *transformationpi = new TGeoCombiTrans(length/2-interspace/2-AliMFTGeometry::kConnectorLength/2-AliMFTGeometry::kConnectorOffset, 
							  id*step-AliMFTGeometry::kFlexHeight/2 + 0.1, zvarnish - AliMFTGeometry::kVarnishThickness/2 - 
							  AliMFTGeometry::kConnectorHeight/2, rotationpi);
  flex->AddNode(connectord,id+total+37, transformationpi);
  }
  

  Double_t boxthickness = 0.05;
  TGeoBBox *boxconnect = new TGeoBBox("boxconnect", (2*AliMFTGeometry::kConnectorThickness+interspace+boxthickness)/2, AliMFTGeometry::kFlexHeight/2-0.04, 
				      AliMFTGeometry::kConnectorHeight/2);
  TGeoBBox *boxremov = new TGeoBBox("boxremov", (2*AliMFTGeometry::kConnectorThickness+interspace)/2, (AliMFTGeometry::kFlexHeight-0.1-step)/2, 
				    AliMFTGeometry::kConnectorHeight/2+0.001);
  TGeoSubtraction *boxconnecto = new TGeoSubtraction(boxconnect, boxremov, NULL, NULL);
  TGeoCompositeShape *boxconnector = new TGeoCompositeShape("boxconnector", boxconnecto);
  TGeoVolume *boxconnectord = new TGeoVolume("boxconnectord", boxconnector, kMedPeek);
  flex->AddNode(boxconnectord,1,new TGeoTranslation(length/2-AliMFTGeometry::kConnectorOffset, -step/2, zvarnish-AliMFTGeometry::kVarnishThickness/2
						    -AliMFTGeometry::kConnectorHeight/2-AliMFTGeometry::kConnectorThickness));
  

  //---------------------------------------------------------------
  
}


TGeoVolumeAssembly* AliMFTFlex::Make_ElectricComponent(Double_t dx, Double_t dy,  Double_t dz, Int_t id)
{
  
  AliMFTGeometry * mftGeom = AliMFTGeometry::Instance();
  Int_t idHalfMFT = mftGeom->GetHalfMFTID(fLadderSeg->GetUniqueID());
  Int_t idHalfDisk = mftGeom->GetHalfDiskID(fLadderSeg->GetUniqueID());
  Int_t idLadder = mftGeom->GetLadderID(fLadderSeg->GetUniqueID());
  //------------------------------------------------------
  TGeoMedium *kmedX7R  = gGeoManager->GetMedium("MFT_X7Rcapacitors$");
  TGeoMedium *kmedX7Rw = gGeoManager->GetMedium("MFT_X7Rweld$");

  TGeoVolumeAssembly* X7R0402 = new TGeoVolumeAssembly(Form("X7R_%d_%d_%d_%d",idHalfMFT,idHalfDisk,idLadder,id));

  TGeoBBox *capacit = new TGeoBBox("capacitor", dx/2, dy/2, dz/2);
  TGeoBBox *weld = new TGeoBBox("weld", (dx/4)/2, dy/2, (dz/2)/2);
  TGeoVolume*  capacitor = new TGeoVolume(Form("capacitor_%d_%d_%d_%d",idHalfMFT,idHalfDisk,idLadder,id), capacit, kmedX7R);
  TGeoVolume*  welding0 = new TGeoVolume(Form("welding0_%d_%d_%d_%d",idHalfMFT,idHalfDisk,idLadder,id), weld, kmedX7Rw);
  TGeoVolume*  welding1 = new TGeoVolume(Form("welding1_%d_%d_%d_%d",idHalfMFT,idHalfDisk,idLadder,id), weld, kmedX7Rw);
  capacitor->SetVisibility(kTRUE);
  capacitor->SetLineColor(kRed);
  capacitor->SetLineWidth(1);
  capacitor->SetFillColor(capacitor->GetLineColor());
  capacitor->SetFillStyle(4000); // 0% transparent

  welding0->SetVisibility(kTRUE);
  welding0->SetLineColor(kGray);
  welding0->SetLineWidth(1);
  welding0->SetFillColor(welding0->GetLineColor());
  welding0->SetFillStyle(4000); // 0% transparent

  welding1->SetVisibility(kTRUE);
  welding1->SetLineColor(kGray);
  welding1->SetLineWidth(1);
  welding1->SetFillColor(welding1->GetLineColor());
  welding1->SetFillStyle(4000); // 0% transparent


  X7R0402->AddNode(capacitor,  1,  new TGeoTranslation(0., 0., 0.)); 
  X7R0402->AddNode(welding0,   1,  new TGeoTranslation( dx/2+(dx/4)/2, 0., (dz/2)/2)); 
  X7R0402->AddNode(welding1,   1,  new TGeoTranslation(-dx/2-(dx/4)/2, 0., (dz/2)/2));
  
  X7R0402->SetVisibility(kTRUE);

  return X7R0402;
  

  //------------------------------------------------------

  /*
  // the medium has to be changed, see ITS capacitors...
  TGeoMedium *kMedCopper = gGeoManager->GetMedium("MFT_Cu$");

  AliMFTGeometry * mftGeom = AliMFTGeometry::Instance();
  Int_t idHalfMFT = mftGeom->GetHalfMFTID(fLadderSeg->GetUniqueID());
  Int_t idHalfDisk = mftGeom->GetHalfDiskID(fLadderSeg->GetUniqueID());
  Int_t idLadder = mftGeom->GetLadderID(fLadderSeg->GetUniqueID());
  
  TGeoVolume* electriccomponent = new TGeoVolume(Form("electric_%d_%d_%d_%d",idHalfMFT,idHalfDisk,idLadder,id), new TGeoBBox("BOX", dy/2, dx/2, dz/2), kMedCopper);
  electriccomponent->SetVisibility(1);
  electriccomponent->SetLineColor(kRed);
  return electriccomponent;
  */
}


TGeoVolume* AliMFTFlex::Make_Lines(Int_t nbsensors, Double_t length, Double_t widthflex,  Double_t thickness)
{
  // One line is built by removing 3 lines of aluminium in the TGeoBBox *layer_def layer. Then one line is made by the 2 remaining aluminium strips. 

  // the initial layer of aluminium
  TGeoBBox *layer_def = new TGeoBBox("layer_def", length/2, widthflex/2, thickness/2);

  // Two holes for fixing and positionning of the FPC on the cold plate
  TGeoTube *hole1 = new TGeoTube("hole1", 0., AliMFTGeometry::kRadiusHole1, thickness/2 + AliMFTGeometry::kEpsilon);
  TGeoTube *hole2 = new TGeoTube("hole2", 0., AliMFTGeometry::kRadiusHole2, thickness/2 + AliMFTGeometry::kEpsilon);

  TGeoTranslation    *t1= new TGeoTranslation ("t1", length/2 - AliMFTGeometry::kHoleShift1, 0., 0.);
  TGeoSubtraction    *layerholesub1 = new TGeoSubtraction(layer_def, hole1, NULL, t1);
  TGeoCompositeShape *layerhole1 = new TGeoCompositeShape("layerhole1", layerholesub1);

  TGeoTranslation    *t2= new TGeoTranslation ("t2", length/2 - AliMFTGeometry::kHoleShift2, 0., 0.);
  TGeoSubtraction    *layerholesub2 = new TGeoSubtraction(layerhole1, hole2, NULL, t2);
  TGeoCompositeShape *layer = new TGeoCompositeShape("layerhole2", layerholesub2);

  TGeoBBox *line[25];
  TGeoTranslation *t[6],*ts[15],*tvdd, *tl[2];
  TGeoSubtraction *layerl[25];
  TGeoCompositeShape *layern[25];
  Int_t istart, istop;
  Int_t kTotalLinesNb=0;
  Int_t kTotalLinesNb1, kTotalLinesNb2;
  Double_t length_line;


  // ----------- two lines along the FPC digital side --------------
  t[0] = new TGeoTranslation ("t0", AliMFTGeometry::kSensorLength/2-AliMFTGeometry::kConnectorOffset/2, -widthflex/2 + 2*AliMFTGeometry::kLineWidth, 0.);    
  line[0]  = new TGeoBBox("line0",  length/2 - AliMFTGeometry::kConnectorOffset/2 - AliMFTGeometry::kSensorLength/2, AliMFTGeometry::kLineWidth/2, 
			  thickness/2 + AliMFTGeometry::kEpsilon);
  layerl[0] = new TGeoSubtraction(layer, line[0], NULL, t[0]);
  layern[0] = new TGeoCompositeShape(Form("layer%d",0), layerl[0]);

  istart = 1; istop = 6;
  for (int iline = istart; iline < istop; iline++){
    t[iline] = new TGeoTranslation (Form("t%d",iline), AliMFTGeometry::kSensorLength/2 - AliMFTGeometry::kConnectorOffset/2, 
				    -widthflex/2 + 2*(iline+1)*AliMFTGeometry::kLineWidth, 0.);
    line[iline]  = new TGeoBBox(Form("line%d",iline),  length/2 - AliMFTGeometry::kConnectorOffset/2 - AliMFTGeometry::kSensorLength/2, AliMFTGeometry::kLineWidth/2, 
				  thickness/2 + AliMFTGeometry::kEpsilon);
    layerl[iline] = new TGeoSubtraction(layern[iline-1], line[iline], NULL, t[iline]);
    layern[iline] = new TGeoCompositeShape(Form("layer%d",iline), layerl[iline]);
    kTotalLinesNb++;
  }

  // ---------  lines for the sensors, one line/sensor -------------
  istart = kTotalLinesNb+1; istop = 6+3*nbsensors;
  for (int iline = istart; iline < istop; iline++){
    length_line=length - AliMFTGeometry::kConnectorOffset - TMath::Nint((iline-6)/3)*AliMFTGeometry::kSensorLength - AliMFTGeometry::kSensorLength/2;
    ts[iline] = new TGeoTranslation (Form("t%d",iline), length/2-length_line/2-AliMFTGeometry::kConnectorOffset, -2*(iline-6)*AliMFTGeometry::kLineWidth+0.5-widthflex/2, 0.);
    line[iline]  = new TGeoBBox(Form("line%d",iline), length_line/2, AliMFTGeometry::kLineWidth/2, thickness/2 + AliMFTGeometry::kEpsilon);
    layerl[iline] = new TGeoSubtraction(layern[iline-1], line[iline], NULL, ts[iline]);
    layern[iline] = new TGeoCompositeShape(Form("layer%d",iline), layerl[iline]);
    kTotalLinesNb++;
  }

  // ---------  an interspace to separate AVDD and DVDD -------------
  kTotalLinesNb++;
  tvdd = new TGeoTranslation ("tvdd", 0., widthflex/2-AliMFTGeometry::kShiftDDGNDline, 0.);    
  line[kTotalLinesNb]  = new TGeoBBox(Form("line%d",kTotalLinesNb),  length/2, 2*AliMFTGeometry::kLineWidth/2, thickness/2 + AliMFTGeometry::kEpsilon);
  layerl[kTotalLinesNb] = new TGeoSubtraction(layern[kTotalLinesNb-1], line[kTotalLinesNb], NULL, tvdd);
  layern[kTotalLinesNb] = new TGeoCompositeShape(Form("layer%d",kTotalLinesNb), layerl[kTotalLinesNb]);
  kTotalLinesNb++;

  // ---------  one line along the FPC analog side -------------  
  istart = kTotalLinesNb; istop = kTotalLinesNb + 2;
  for (int iline = istart; iline < istop; iline++){
    length_line=length - AliMFTGeometry::kConnectorOffset;
    tl[iline-istart] = new TGeoTranslation (Form("tl%d",iline), length/2-length_line/2-AliMFTGeometry::kConnectorOffset, 
					    widthflex/2-AliMFTGeometry::kShiftline-2.*(iline-istart)*AliMFTGeometry::kLineWidth, 0.);
    line[iline]  = new TGeoBBox(Form("line%d",iline), length_line/2, AliMFTGeometry::kLineWidth/2, thickness/2 + AliMFTGeometry::kEpsilon);
    layerl[iline] = new TGeoSubtraction(layern[iline-1], line[iline], NULL, tl[iline-istart]);
    layern[iline] = new TGeoCompositeShape(Form("layer%d",iline), layerl[iline]);
    kTotalLinesNb++;
  }

  AliMFTGeometry * mftGeom = AliMFTGeometry::Instance();
  Int_t idHalfMFT = mftGeom->GetHalfMFTID(fLadderSeg->GetUniqueID());
  Int_t idHalfDisk = mftGeom->GetHalfDiskID(fLadderSeg->GetUniqueID());
  Int_t idLadder = mftGeom->GetLadderID(fLadderSeg->GetUniqueID());

  TGeoMedium *kMedAlu = gGeoManager->GetMedium("MFT_Alu$");

  TGeoVolume *lineslayer = new TGeoVolume(Form("lineslayer_%d_%d_%d",idHalfMFT,idHalfDisk,idLadder), layern[kTotalLinesNb-1], kMedAlu);
  lineslayer->SetVisibility(1);
  lineslayer->SetLineColor(kBlue);
  return lineslayer;
}


TGeoVolume* AliMFTFlex::Make_AGND_DGND(Double_t length, Double_t widthflex,  Double_t thickness)
{  
  // AGND and DGND layers
  TGeoBBox *layer = new TGeoBBox("layer", length/2, widthflex/2, thickness/2);
  TGeoTube *hole1 = new TGeoTube("hole1", 0., AliMFTGeometry::kRadiusHole1, thickness/2 + AliMFTGeometry::kEpsilon);
  TGeoTube *hole2 = new TGeoTube("hole2", 0., AliMFTGeometry::kRadiusHole2, thickness/2 + AliMFTGeometry::kEpsilon);
  
  TGeoTranslation    *t1= new TGeoTranslation ("t1", length/2-AliMFTGeometry::kHoleShift1, 0., 0.);
  TGeoSubtraction    *layerholesub1 = new TGeoSubtraction(layer, hole1, NULL, t1);
  TGeoCompositeShape *layerhole1 = new TGeoCompositeShape("layerhole1", layerholesub1);

  TGeoTranslation    *t2= new TGeoTranslation ("t2", length/2-AliMFTGeometry::kHoleShift2, 0., 0.);
  TGeoSubtraction    *layerholesub2 = new TGeoSubtraction(layerhole1, hole2, NULL, t2);
  TGeoCompositeShape *layerhole2 = new TGeoCompositeShape("layerhole2", layerholesub2);

  //--------------
  TGeoBBox *line[3];
  TGeoTranslation *t[3];
  TGeoCompositeShape *layern[3];
  TGeoSubtraction *layerl[3];
  Double_t length_line;
  length_line=length - AliMFTGeometry::kConnectorOffset;

  // First, the two lines along the FPC side
  t[0] = new TGeoTranslation("t0", length/2-length_line/2-AliMFTGeometry::kConnectorOffset, widthflex/2 - AliMFTGeometry::kShiftline, 0.);
  line[0]  = new TGeoBBox("line0",  length/2 - AliMFTGeometry::kConnectorOffset/2, AliMFTGeometry::kLineWidth/2, thickness/2 + AliMFTGeometry::kEpsilon);
  layerl[0] = new TGeoSubtraction(layerhole2, line[0], NULL, t[0]);
  layern[0] = new TGeoCompositeShape(Form("layer%d",0), layerl[0]);

  t[1] = new TGeoTranslation("t1", length/2-length_line/2-AliMFTGeometry::kConnectorOffset, 
			     widthflex/2 - AliMFTGeometry::kShiftline - 2*AliMFTGeometry::kLineWidth, 0.);
  line[1]  = new TGeoBBox("line1",  length/2 - AliMFTGeometry::kConnectorOffset/2, AliMFTGeometry::kLineWidth/2, 
			  thickness/2 + AliMFTGeometry::kEpsilon);
  layerl[1] = new TGeoSubtraction(layern[0], line[1], NULL, t[1]);
  layern[1] = new TGeoCompositeShape(Form("layer%d",1), layerl[1]);

  // Now the interspace to separate the AGND et DGND --> same interspace compare the AVDD et DVDD
  t[2] = new TGeoTranslation("t2", length/2-length_line/2, widthflex/2 - AliMFTGeometry::kShiftDDGNDline, 0.);
  line[2]  = new TGeoBBox("line2",  length/2 - AliMFTGeometry::kConnectorOffset/2, AliMFTGeometry::kLineWidth, 
			  thickness/2 + AliMFTGeometry::kEpsilon);
  layerl[2] = new TGeoSubtraction(layern[1], line[2], NULL, t[2]);
  layern[2] = new TGeoCompositeShape(Form("layer%d",2), layerl[2]);

  //--------------

  AliMFTGeometry * mftGeom = AliMFTGeometry::Instance();
  Int_t idHalfMFT = mftGeom->GetHalfMFTID(fLadderSeg->GetUniqueID());
  Int_t idHalfDisk = mftGeom->GetHalfDiskID(fLadderSeg->GetUniqueID());
  Int_t idLadder = mftGeom->GetLadderID(fLadderSeg->GetUniqueID());

  TGeoMedium *kMedAlu = gGeoManager->GetMedium("MFT_Alu$");
  TGeoVolume *alulayer = new TGeoVolume(Form("alulayer_%d_%d_%d",idHalfMFT,idHalfDisk,idLadder), layern[2], kMedAlu);
  alulayer->SetVisibility(1);
  alulayer->SetLineColor(kBlue);
  return alulayer;
}


TGeoVolume* AliMFTFlex::Make_Kapton(Double_t length, Double_t widthflex, Double_t thickness)
{
  TGeoBBox *layer = new TGeoBBox("layer", length/2, widthflex/2, thickness/2);
  // Two holes for fixing and positionning of the FPC on the cold plate
  TGeoTube *hole1 = new TGeoTube("hole1", 0., AliMFTGeometry::kRadiusHole1, thickness/2+AliMFTGeometry::kEpsilon);
  TGeoTube *hole2 = new TGeoTube("hole2", 0., AliMFTGeometry::kRadiusHole2, thickness/2+AliMFTGeometry::kEpsilon);
  
  TGeoTranslation    *t1= new TGeoTranslation ("t1", length/2-AliMFTGeometry::kHoleShift1, 0., 0.);
  TGeoSubtraction    *layerholesub1 = new TGeoSubtraction(layer, hole1, NULL, t1);
  TGeoCompositeShape *layerhole1 = new TGeoCompositeShape("layerhole1", layerholesub1);

  TGeoTranslation    *t2= new TGeoTranslation ("t2", length/2-AliMFTGeometry::kHoleShift2, 0., 0.);
  TGeoSubtraction    *layerholesub2 = new TGeoSubtraction(layerhole1, hole2, NULL, t2);
  TGeoCompositeShape *layerhole2 = new TGeoCompositeShape("layerhole2", layerholesub2);

  AliMFTGeometry * mftGeom = AliMFTGeometry::Instance();
  Int_t idHalfMFT = mftGeom->GetHalfMFTID(fLadderSeg->GetUniqueID());
  Int_t idHalfDisk = mftGeom->GetHalfDiskID(fLadderSeg->GetUniqueID());
  Int_t idLadder = mftGeom->GetLadderID(fLadderSeg->GetUniqueID());

  TGeoMedium *kMedKapton = gGeoManager->GetMedium("MFT_Kapton$");
  TGeoVolume *kaptonlayer = new TGeoVolume(Form("kaptonlayer_%d_%d_%d",idHalfMFT,idHalfDisk,idLadder), layerhole2, kMedKapton);
  kaptonlayer->SetVisibility(1);
  kaptonlayer->SetLineColor(kYellow);
  return kaptonlayer;
}


TGeoVolume* AliMFTFlex::Make_Varnish(Double_t length, Double_t widthflex,  Double_t thickness, Int_t iflag)
{
  TGeoBBox *layer = new TGeoBBox("layer", length/2, widthflex/2, thickness/2);
  // Two holes for fixing and positionning of the FPC on the cold plate
  TGeoTube *hole1 = new TGeoTube("hole1", 0., AliMFTGeometry::kRadiusHole1, thickness/2+AliMFTGeometry::kEpsilon);
  TGeoTube *hole2 = new TGeoTube("hole2", 0., AliMFTGeometry::kRadiusHole2, thickness/2+AliMFTGeometry::kEpsilon);
  
  TGeoTranslation    *t1= new TGeoTranslation ("t1", length/2-AliMFTGeometry::kHoleShift1, 0., 0.);
  TGeoSubtraction    *layerholesub1 = new TGeoSubtraction(layer, hole1, NULL, t1);
  TGeoCompositeShape *layerhole1 = new TGeoCompositeShape("layerhole1", layerholesub1);

  TGeoTranslation    *t2= new TGeoTranslation ("t2", length/2-AliMFTGeometry::kHoleShift2, 0., 0.);
  TGeoSubtraction    *layerholesub2 = new TGeoSubtraction(layerhole1, hole2, NULL, t2);
  TGeoCompositeShape *layerhole2 = new TGeoCompositeShape("layerhole2", layerholesub2);

  AliMFTGeometry * mftGeom = AliMFTGeometry::Instance();
  Int_t idHalfMFT = mftGeom->GetHalfMFTID(fLadderSeg->GetUniqueID());
  Int_t idHalfDisk = mftGeom->GetHalfDiskID(fLadderSeg->GetUniqueID());
  Int_t idLadder = mftGeom->GetLadderID(fLadderSeg->GetUniqueID());

  TGeoMedium *kMedVarnish = gGeoManager->GetMedium("MFT_Epoxy$");  // we assume that varnish = epoxy ...
  TGeoVolume *varnishlayer = new TGeoVolume(Form("varnishlayer_%d_%d_%d_%d",idHalfMFT,idHalfDisk,idLadder,iflag), layerhole2, kMedVarnish);
  varnishlayer->SetVisibility(1);
  varnishlayer->SetLineColor(kGreen-1);
  return varnishlayer;
}

