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
  Revision 1.12  2000/06/13 13:06:28  jbarbosa
  Fixed compiling error for HP (multiple declaration)

  Revision 1.11  2000/06/12 15:35:44  jbarbosa
  Cleaned up version.

  Revision 1.10  2000/06/09 14:59:25  jbarbosa
  New default version. No setters needed, no hits.

  Revision 1.9  2000/05/31 08:19:38  jbarbosa
  Fixed bug in StepManager

  Revision 1.8  2000/05/26 17:30:08  jbarbosa
  Cerenkov angle now stored within cerenkov data structure.

  Revision 1.7  2000/05/18 10:31:36  jbarbosa
  Fixed positioning of spacers inside freon.
  Fixed positioning of proximity gap
  inside methane.
  Fixed cut on neutral particles in the StepManager.

  Revision 1.6  2000/04/28 11:51:58  morsch
   Dimensions of arrays hits and Ckov_data corrected.

  Revision 1.5  2000/04/19 13:28:46  morsch
  Major changes in geometry (parametrised), materials (updated) and
  step manager (diagnostics) (JB, AM)

*/



/////////////////////////////////////////////////////////////
//  Manager and hits classes for set: RICH default version //
/////////////////////////////////////////////////////////////

#include <TTUBE.h>
#include <TNode.h> 
#include <TRandom.h> 

#include "AliRICHv0.h"
#include "AliRICHSegmentation.h"
#include "AliRICHResponse.h"
#include "AliRICHSegmentationV0.h"
#include "AliRICHResponseV0.h"
#include "AliRICHGeometry.h"
#include "AliRun.h"
#include "AliMC.h"
#include "iostream.h"
#include "AliCallf77.h"
#include "AliConst.h" 
#include "AliPDG.h" 
#include "TGeant3.h"

ClassImp(AliRICHv0)
    
//___________________________________________
AliRICHv0::AliRICHv0() : AliRICH()
{

// Default constructor

    //fChambers = 0;
}

//___________________________________________
AliRICHv0::AliRICHv0(const char *name, const char *title)
    : AliRICH(name,title)
{
    //
// Version 0
// Default Segmentation, no hits
    AliRICHSegmentationV0* segmentationV0 = new AliRICHSegmentationV0;
//
//  Segmentation parameters
    segmentationV0->SetPadSize(0.84,0.80);
    segmentationV0->SetDAnod(0.84/2);
//
//  Geometry parameters
    AliRICHGeometry* geometry = new AliRICHGeometry;
    geometry->SetGapThickness(8);
    geometry->SetProximityGapThickness(.4);
    geometry->SetQuartzLength(131);
    geometry->SetQuartzWidth(126.2);
    geometry->SetQuartzThickness(.5);
    geometry->SetOuterFreonLength(131);
    geometry->SetOuterFreonWidth(40.3);
    geometry->SetInnerFreonLength(131);
    geometry->SetInnerFreonWidth(40.3);
    geometry->SetFreonThickness(1);
//
//  Response parameters
    AliRICHResponseV0*  responseV0   = new AliRICHResponseV0;
    responseV0->SetSigmaIntegration(5.);
    responseV0->SetChargeSlope(40.);
    responseV0->SetChargeSpread(0.18, 0.18);
    responseV0->SetMaxAdc(1024);
    responseV0->SetAlphaFeedback(0.05);
    responseV0->SetEIonisation(26.e-9);
    responseV0->SetSqrtKx3(0.77459667);
    responseV0->SetKx2(0.962);
    responseV0->SetKx4(0.379);
    responseV0->SetSqrtKy3(0.77459667);
    responseV0->SetKy2(0.962);
    responseV0->SetKy4(0.379);
    responseV0->SetPitch(0.25);
//
//
//    AliRICH *RICH = (AliRICH *) gAlice->GetDetector("RICH"); 
    
    fCkovNumber=0;
    fFreonProd=0;
    Int_t i=0;
    
    fChambers = new TObjArray(kNCH);
    for (i=0; i<kNCH; i++) {
      
      (*fChambers)[i] = new AliRICHChamber();  
      
    }
  
    for (i=0; i<kNCH; i++) {
      SetGeometryModel(i,geometry);
      SetSegmentationModel(i, segmentationV0);
      SetResponseModel(i, responseV0);
      SetNsec(i,1);
      SetDebugLevel(0);
    }
}


//___________________________________________
void AliRICHv0::CreateGeometry()
{
    //
    // Create the geometry for RICH version 1
    //
    // Modified by:  N. Colonna (INFN - BARI, Nicola.Colonna@ba.infn.it) 
    //               R.A. Fini  (INFN - BARI, Rosanna.Fini@ba.infn.it) 
    //               R.A. Loconsole (Bari University, loco@riscom.ba.infn.it) 
    //
    //Begin_Html
    /*
      <img src="picts/AliRICHv1.gif">
    */
    //End_Html
    //Begin_Html
    /*
      <img src="picts/AliRICHv1Tree.gif">
    */
    //End_Html

  AliRICH *pRICH = (AliRICH *) gAlice->GetDetector("RICH"); 
  AliRICHSegmentation*  segmentation;
  AliRICHGeometry*  geometry;
  AliRICHChamber*       iChamber;

  iChamber = &(pRICH->Chamber(0));
  segmentation=iChamber->GetSegmentationModel(0);
  geometry=iChamber->GetGeometryModel();

  Float_t distance;
  distance = geometry->GetFreonThickness()/2 + geometry->GetQuartzThickness() + geometry->GetGapThickness();
  geometry->SetRadiatorToPads(distance);
    
    
    Int_t *idtmed = fIdtmed->GetArray()-999;
    
    Int_t i;
    Float_t zs;
    Int_t idrotm[1099];
    Float_t par[3];
    
    // --- Define the RICH detector 
    //     External aluminium box 
    par[0] = 71.1;
    par[1] = 11.5;                 //Original Settings
    par[2] = 73.15;
    /*par[0] = 73.15;
    par[1] = 11.5;
    par[2] = 71.1;*/
    gMC->Gsvolu("RICH", "BOX ", idtmed[1009], par, 3);
    
    //     Sensitive part of the whole RICH 
    par[0] = 64.8;
    par[1] = 11.5;                 //Original Settings
    par[2] = 66.55;
    /*par[0] = 66.55;
    par[1] = 11.5;
    par[2] = 64.8;*/
    gMC->Gsvolu("SRIC", "BOX ", idtmed[1000], par, 3);
    
    //     Honeycomb 
    par[0] = 63.1;
    par[1] = .188;                 //Original Settings
    par[2] = 66.55;
    /*par[0] = 66.55;
    par[1] = .188;
    par[2] = 63.1;*/
    gMC->Gsvolu("HONE", "BOX ", idtmed[1001], par, 3);
    
    //     Aluminium sheet 
    par[0] = 63.1;
    par[1] = .025;                 //Original Settings
    par[2] = 66.55;
    /*par[0] = 66.5;
    par[1] = .025;
    par[2] = 63.1;*/
    gMC->Gsvolu("ALUM", "BOX ", idtmed[1009], par, 3);
    
    //     Quartz 
    par[0] = geometry->GetQuartzWidth()/2;
    par[1] = geometry->GetQuartzThickness()/2;
    par[2] = geometry->GetQuartzLength()/2;
    /*par[0] = 63.1;
    par[1] = .25;                  //Original Settings
    par[2] = 65.5;*/
    /*par[0] = geometry->GetQuartzWidth()/2;
    par[1] = geometry->GetQuartzThickness()/2;
    par[2] = geometry->GetQuartzLength()/2;*/
    //printf("\n\n\n\n\n\n\n\\n\n\n\n Gap Thickness: %f %f %f\n\n\n\n\n\n\n\n\n\n\n\n\n\n",par[0],par[1],par[2]);
    gMC->Gsvolu("QUAR", "BOX ", idtmed[1002], par, 3);
    
    //     Spacers (cylinders) 
    par[0] = 0.;
    par[1] = .5;
    par[2] = geometry->GetFreonThickness()/2;
    gMC->Gsvolu("SPAC", "TUBE", idtmed[1002], par, 3);
    
    //     Opaque quartz 
    par[0] = 61.95;
    par[1] = .2;                   //Original Settings
    par[2] = 66.5;
    /*par[0] = 66.5;
    par[1] = .2;
    par[2] = 61.95;*/
    gMC->Gsvolu("OQUA", "BOX ", idtmed[1007], par, 3);
  
    //     Frame of opaque quartz
    par[0] = geometry->GetOuterFreonWidth()/2;
    par[1] = geometry->GetFreonThickness()/2;
    par[2] = geometry->GetOuterFreonLength()/2 + 1; 
    /*par[0] = 20.65;
    par[1] = .5;                   //Original Settings
    par[2] = 66.5;*/
    /*par[0] = 66.5;
    par[1] = .5;
    par[2] = 20.65;*/
    gMC->Gsvolu("OQF1", "BOX ", idtmed[1007], par, 3);

    par[0] = geometry->GetInnerFreonWidth()/2;
    par[1] = geometry->GetFreonThickness()/2;
    par[2] = geometry->GetInnerFreonLength()/2 + 1; 
    gMC->Gsvolu("OQF2", "BOX ", idtmed[1007], par, 3);
    
    //     Little bar of opaque quartz 
    par[0] = .275;
    par[1] = geometry->GetQuartzThickness()/2;
    par[2] = geometry->GetInnerFreonLength()/2 - 2.4; 
    /*par[0] = .275;
    par[1] = .25;                   //Original Settings
    par[2] = 63.1;*/
    /*par[0] = 63.1;
    par[1] = .25;
    par[2] = .275;*/
    gMC->Gsvolu("BARR", "BOX ", idtmed[1007], par, 3);
    
    //     Freon 
    par[0] = geometry->GetOuterFreonWidth()/2;
    par[1] = geometry->GetFreonThickness()/2;
    par[2] = geometry->GetOuterFreonLength()/2; 
    /*par[0] = 20.15;
    par[1] = .5;                   //Original Settings
    par[2] = 65.5;*/
    /*par[0] = 65.5;
    par[1] = .5;
    par[2] = 20.15;*/
    gMC->Gsvolu("FRE1", "BOX ", idtmed[1003], par, 3);

    par[0] = geometry->GetInnerFreonWidth()/2;
    par[1] = geometry->GetFreonThickness()/2;
    par[2] = geometry->GetInnerFreonLength()/2; 
    gMC->Gsvolu("FRE2", "BOX ", idtmed[1003], par, 3);
    
    //     Methane 
    par[0] = 64.8;
    par[1] = geometry->GetGapThickness()/2;
    //printf("\n\n\n\n\n\n\n\\n\n\n\n Gap Thickness: %f\n\n\n\n\n\n\n\n\n\n\n\n\n\n",par[1]);
    par[2] = 64.8;
    gMC->Gsvolu("META", "BOX ", idtmed[1004], par, 3);
    
    //     Methane gap 
    par[0] = 64.8;
    par[1] = geometry->GetProximityGapThickness()/2;
    //printf("\n\n\n\n\n\n\n\\n\n\n\n Gap Thickness: %f\n\n\n\n\n\n\n\n\n\n\n\n\n\n",par[1]);
    par[2] = 64.8;
    gMC->Gsvolu("GAP ", "BOX ", idtmed[1008], par, 3);
    
    //     CsI photocathode 
    par[0] = 64.8;
    par[1] = .25;
    par[2] = 64.8;
    gMC->Gsvolu("CSI ", "BOX ", idtmed[1005], par, 3);
    
    //     Anode grid 
    par[0] = 0.;
    par[1] = .001;
    par[2] = 20.;
    gMC->Gsvolu("GRID", "TUBE", idtmed[1006], par, 3);
    
    // --- Places the detectors defined with GSVOLU 
    //     Place material inside RICH 
    gMC->Gspos("SRIC", 1, "RICH", 0., 0., 0., 0, "ONLY");
    
    gMC->Gspos("ALUM", 1, "SRIC", 0., 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 -.05 - .376 -.025, 0., 0, "ONLY");
    gMC->Gspos("HONE", 1, "SRIC", 0., 1.276- geometry->GetGapThickness()/2  - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 -.05 - .188, 0., 0, "ONLY");
    gMC->Gspos("ALUM", 2, "SRIC", 0., 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .025, 0., 0, "ONLY");
    gMC->Gspos("OQUA", 1, "SRIC", 0., 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .2, 0., 0, "ONLY");
    
    AliMatrix(idrotm[1019], 0., 0., 90., 0., 90., 90.);
    
    Int_t nspacers = (Int_t)(TMath::Abs(geometry->GetInnerFreonLength()/14.4));
    //printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n Spacers:%d\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",nspacers); 

    //printf("Nspacers: %d", nspacers);
    
    //for (i = 1; i <= 9; ++i) {
      //zs = (5 - i) * 14.4;                       //Original settings 
    for (i = 0; i < nspacers; i++) {
	zs = (TMath::Abs(nspacers/2) - i) * 14.4;
	gMC->Gspos("SPAC", i, "FRE1", 6.7, 0., zs, idrotm[1019], "ONLY");  //Original settings 
	//gMC->Gspos("SPAC", i, "FRE1", zs, 0., 6.7, idrotm[1019], "ONLY"); 
    }
    //for (i = 10; i <= 18; ++i) {
      //zs = (14 - i) * 14.4;                       //Original settings 
    for (i = nspacers; i < nspacers*2; ++i) {
	zs = (nspacers + TMath::Abs(nspacers/2) - i) * 14.4;
	gMC->Gspos("SPAC", i, "FRE1", -6.7, 0., zs, idrotm[1019], "ONLY"); //Original settings  
	//gMC->Gspos("SPAC", i, "FRE1", zs, 0., -6.7, idrotm[1019], "ONLY");  
    }

    //for (i = 1; i <= 9; ++i) {
      //zs = (5 - i) * 14.4;                       //Original settings 
      for (i = 0; i < nspacers; i++) {
	zs = (TMath::Abs(nspacers/2) - i) * 14.4;
	gMC->Gspos("SPAC", i, "FRE2", 6.7, 0., zs, idrotm[1019], "ONLY");  //Original settings 
	//gMC->Gspos("SPAC", i, "FRE2", zs, 0., 6.7, idrotm[1019], "ONLY");
    }
    //for (i = 10; i <= 18; ++i) {
      //zs = (5 - i) * 14.4;                       //Original settings 
      for (i = nspacers; i < nspacers*2; ++i) {
	zs = (nspacers + TMath::Abs(nspacers/2) - i) * 14.4;
	gMC->Gspos("SPAC", i, "FRE2", -6.7, 0., zs, idrotm[1019], "ONLY");  //Original settings 
	//gMC->Gspos("SPAC", i, "FRE2", zs, 0., -6.7, idrotm[1019], "ONLY");
    }
    
    /*gMC->Gspos("FRE1", 1, "OQF1", 0., 0., 0., 0, "ONLY");
    gMC->Gspos("FRE2", 1, "OQF2", 0., 0., 0., 0, "ONLY");
    gMC->Gspos("OQF1", 1, "SRIC", 31.3, -4.724, 41.3, 0, "ONLY");
    gMC->Gspos("OQF2", 2, "SRIC", 0., -4.724, 0., 0, "ONLY");
    gMC->Gspos("OQF1", 3, "SRIC", -31.3, -4.724, -41.3, 0, "ONLY");
    gMC->Gspos("BARR", 1, "QUAR", -21.65, 0., 0., 0, "ONLY");           //Original settings 
    gMC->Gspos("BARR", 2, "QUAR", 21.65, 0., 0., 0, "ONLY");            //Original settings 
    gMC->Gspos("QUAR", 1, "SRIC", 0., -3.974, 0., 0, "ONLY");
    gMC->Gspos("GAP ", 1, "META", 0., 4.8, 0., 0, "ONLY");
    gMC->Gspos("META", 1, "SRIC", 0., 1.276, 0., 0, "ONLY");
    gMC->Gspos("CSI ", 1, "SRIC", 0., 6.526, 0., 0, "ONLY");*/


    gMC->Gspos("FRE1", 1, "OQF1", 0., 0., 0., 0, "ONLY");
    gMC->Gspos("FRE2", 1, "OQF2", 0., 0., 0., 0, "ONLY");
    gMC->Gspos("OQF1", 1, "SRIC", geometry->GetOuterFreonWidth()/2 + geometry->GetInnerFreonWidth()/2, 1.276 - geometry->GetGapThickness()/2- geometry->GetQuartzThickness() -geometry->GetFreonThickness()/2, 0., 0, "ONLY"); //Original settings (31.3)
    gMC->Gspos("OQF2", 2, "SRIC", 0., 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()/2, 0., 0, "ONLY");          //Original settings 
    gMC->Gspos("OQF1", 3, "SRIC", - (geometry->GetOuterFreonWidth()/2 + geometry->GetInnerFreonWidth()/2), 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()/2, 0., 0, "ONLY");       //Original settings (-31.3)
    gMC->Gspos("BARR", 1, "QUAR", -21.65, 0., 0., 0, "ONLY");           //Original settings 
    gMC->Gspos("BARR", 2, "QUAR", 21.65, 0., 0., 0, "ONLY");            //Original settings 
    gMC->Gspos("QUAR", 1, "SRIC", 0., 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness()/2, 0., 0, "ONLY");
    gMC->Gspos("GAP ", 1, "META", 0., geometry->GetGapThickness()/2 - geometry->GetProximityGapThickness()/2 - 0.0001, 0., 0, "ONLY");
    gMC->Gspos("META", 1, "SRIC", 0., 1.276, 0., 0, "ONLY");
    gMC->Gspos("CSI ", 1, "SRIC", 0., 1.276 + geometry->GetGapThickness()/2 + .25, 0., 0, "ONLY");

    //printf("Position of the gap: %f to %f\n", 1.276 + geometry->GetGapThickness()/2 - geometry->GetProximityGapThickness()/2 - .2, 1.276 + geometry->GetGapThickness()/2 - geometry->GetProximityGapThickness()/2 + .2);
    
    //     Place RICH inside ALICE apparatus 
  
    AliMatrix(idrotm[1000], 90., 0., 70.69, 90., 19.31, -90.);
    AliMatrix(idrotm[1001], 90., -20., 90., 70., 0., 0.);
    AliMatrix(idrotm[1002], 90., 0., 90., 90., 0., 0.);
    AliMatrix(idrotm[1003], 90., 20., 90., 110., 0., 0.);
    AliMatrix(idrotm[1004], 90., 340., 108.2, 70., 18.2, 70.);
    AliMatrix(idrotm[1005], 90., 0., 109.31, 90., 19.31, 90.);
    AliMatrix(idrotm[1006], 90., 20., 108.2, 110., 18.2, 110.);
    
    gMC->Gspos("RICH", 1, "ALIC", 0., 471.9, 165.26,     idrotm[1000], "ONLY");
    gMC->Gspos("RICH", 2, "ALIC", 171., 470., 0.,        idrotm[1001], "ONLY");
    gMC->Gspos("RICH", 3, "ALIC", 0., 500., 0.,          idrotm[1002], "ONLY");
    gMC->Gspos("RICH", 4, "ALIC", -171., 470., 0.,       idrotm[1003], "ONLY");
    gMC->Gspos("RICH", 5, "ALIC", 161.4, 443.4, -165.3,  idrotm[1004], "ONLY");
    gMC->Gspos("RICH", 6, "ALIC", 0., 471.9, -165.3,     idrotm[1005], "ONLY");
    gMC->Gspos("RICH", 7, "ALIC", -161.4, 443.4, -165.3, idrotm[1006], "ONLY");
    
}


//___________________________________________
void AliRICHv0::CreateMaterials()
{
    //
    // *** DEFINITION OF AVAILABLE RICH MATERIALS *** 
    // ORIGIN    : NICK VAN EIJNDHOVEN 
    // Modified by:  N. Colonna (INFN - BARI, Nicola.Colonna@ba.infn.it) 
    //               R.A. Fini  (INFN - BARI, Rosanna.Fini@ba.infn.it) 
    //               R.A. Loconsole (Bari University, loco@riscom.ba.infn.it) 
    //
    Int_t   isxfld = gAlice->Field()->Integ();
    Float_t sxmgmx = gAlice->Field()->Max();
    Int_t i;

    /************************************Antonnelo's Values (14-vectors)*****************************************/
    /*
    Float_t ppckov[14] = { 5.63e-9,5.77e-9,5.9e-9,6.05e-9,6.2e-9,6.36e-9,6.52e-9,
			   6.7e-9,6.88e-9,7.08e-9,7.3e-9,7.51e-9,7.74e-9,8e-9 };
    Float_t rIndexQuarz[14] = { 1.528309,1.533333,
				 1.538243,1.544223,1.550568,1.55777,
				 1.565463,1.574765,1.584831,1.597027,
			       1.611858,1.6277,1.6472,1.6724 };
    Float_t rIndexOpaqueQuarz[14] = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };
    Float_t rIndexMethane[14] = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };
    Float_t rIndexGrid[14] = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };
    Float_t abscoFreon[14] = { 179.0987,179.0987,
				179.0987,179.0987,179.0987,142.92,56.65,13.95,10.43,7.07,2.03,.5773,.33496,0. };
    //Float_t abscoFreon[14] = { 1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,
	//			 1e-5,1e-5,1e-5,1e-5,1e-5 };
    Float_t abscoQuarz[14] = { 64.035,39.98,35.665,31.262,27.527,22.815,21.04,17.52,
				14.177,9.282,4.0925,1.149,.3627,.10857 };
    Float_t abscoOpaqueQuarz[14] = { 1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,
				 1e-5,1e-5,1e-5,1e-5,1e-5 };
    Float_t abscoCsI[14] = { 1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,
			      1e-4,1e-4,1e-4,1e-4 };
    Float_t abscoMethane[14] = { 1e6,1e6,1e6,1e6,1e6,1e6,1e6,1e6,1e6,1e6,1e6,
				  1e6,1e6,1e6 };
    Float_t abscoGrid[14] = { 1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,
			      1e-4,1e-4,1e-4,1e-4 };
    Float_t efficAll[14] = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };
    Float_t efficCsI[14] = { 6e-4,.005,.0075,.01125,.045,.117,.135,.16575,
			      .17425,.1785,.1836,.1904,.1938,.221 };
    Float_t efficGrid[14] = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };
    */
   
    
    /**********************************End of Antonnelo's Values**********************************/
    
    /**********************************Values from rich_media.f (31-vectors)**********************************/
    

    //Photons energy intervals
    Float_t ppckov[26];
    for (i=0;i<26;i++) 
    {
	ppckov[i] = (Float_t(i)*0.1+5.5)*1e-9;
	//printf ("Energy intervals: %e\n",ppckov[i]);
    }
    
    
    //Refraction index for quarz
    Float_t rIndexQuarz[26];
    Float_t  e1= 10.666;
    Float_t  e2= 18.125;
    Float_t  f1= 46.411;
    Float_t  f2= 228.71;
    for (i=0;i<26;i++)
    {
	Float_t ene=ppckov[i]*1e9;
	Float_t a=f1/(e1*e1 - ene*ene);
	Float_t b=f2/(e2*e2 - ene*ene);
	rIndexQuarz[i] = TMath::Sqrt(1. + a + b );
	//printf ("rIndexQuarz: %e\n",rIndexQuarz[i]);
    } 
    
    //Refraction index for opaque quarz, methane and grid
    Float_t rIndexOpaqueQuarz[26];
    Float_t rIndexMethane[26];
    Float_t rIndexGrid[26];
    for (i=0;i<26;i++)
    {
	rIndexOpaqueQuarz[i]=1;
	rIndexMethane[i]=1.000444;
	rIndexGrid[i]=1;
	//printf ("rIndexOpaqueQuarz , etc: %e, %e, %e\n",rIndexOpaqueQuarz[i], rIndexMethane[i], rIndexGrid[i]=1);
    } 
    
    //Absorption index for freon
    Float_t abscoFreon[26] = {179.0987, 179.0987, 179.0987, 179.0987, 179.0987,  179.0987, 179.0987, 179.0987, 
	 		       179.0987, 142.9206, 56.64957, 25.58622, 13.95293, 12.03905, 10.42953, 8.804196, 
			       7.069031, 4.461292, 2.028366, 1.293013, .577267,   .40746,  .334964, 0., 0., 0.};
    
    //Absorption index for quarz
    /*Float_t Qzt [21] = {.0,.0,.005,.04,.35,.647,.769,.808,.829,.844,.853,.858,.869,.887,.903,.902,.902,
	 		.906,.907,.907,.907};
    Float_t Wavl2[] = {150.,155.,160.0,165.0,170.0,175.0,180.0,185.0,190.0,195.0,200.0,205.0,210.0,
	 	       215.0,220.0,225.0,230.0,235.0,240.0,245.0,250.0};		 		 
    Float_t abscoQuarz[31];	     
    for (Int_t i=0;i<31;i++)
    {
	Float_t Xlam = 1237.79 / (ppckov[i]*1e9);
	if (Xlam <= 160) abscoQuarz[i] = 0;
	if (Xlam > 250) abscoQuarz[i] = 1;
	else 
	{
	    for (Int_t j=0;j<21;j++)
	    {
		//printf ("Passed\n");
		if (Xlam > Wavl2[j] && Xlam < Wavl2[j+1])
		{
		    Float_t Dabs = (Qzt[j+1] - Qzt[j])/(Wavl2[j+1] - Wavl2[j]);
		    Float_t Abso = Qzt[j] + Dabs*(Xlam - Wavl2[j]);
		    abscoQuarz[i] = -5.0/(TMath::Log(Abso));
		} 
	    }
	}
	printf ("abscoQuarz: %e abscoFreon: %e for energy: %e\n",abscoQuarz[i],abscoFreon[i],ppckov[i]);
    }*/

    /*Float_t abscoQuarz[31] = {49.64211, 48.41296, 47.46989, 46.50492, 45.13682, 44.47883, 43.1929 , 41.30922, 40.5943 ,
			       39.82956, 38.98623, 38.6247 , 38.43448, 37.41084, 36.22575, 33.74852, 30.73901, 24.25086, 
			       17.94531, 11.88753, 5.99128,  3.83503,  2.36661,  1.53155, 1.30582, 1.08574, .8779708, 
			       .675275, 0., 0., 0.};
    
    for (Int_t i=0;i<31;i++)
    {
	abscoQuarz[i] = abscoQuarz[i]/10;
    }*/

    Float_t abscoQuarz [26] = {105.8, 65.52, 48.58, 42.85, 35.79, 31.262, 28.598, 27.527, 25.007, 22.815, 21.004,
				19.266, 17.525, 15.878, 14.177, 11.719, 9.282, 6.62, 4.0925, 2.601, 1.149, .667, .3627,
				.192, .1497, .10857};
    
    //Absorption index for methane
    Float_t abscoMethane[26];
    for (i=0;i<26;i++) 
    {
	abscoMethane[i]=AbsoCH4(ppckov[i]*1e9); 
	//printf("abscoMethane: %e for energy: %e\n", abscoMethane[i],ppckov[i]*1e9);
    }
    
    //Absorption index for opaque quarz, csi and grid, efficiency for all and grid
    Float_t abscoOpaqueQuarz[26];
    Float_t abscoCsI[26];
    Float_t abscoGrid[26];
    Float_t efficAll[26];
    Float_t efficGrid[26];
    for (i=0;i<26;i++)
    { 
	abscoOpaqueQuarz[i]=1e-5; 
	abscoCsI[i]=1e-4; 
	abscoGrid[i]=1e-4; 
	efficAll[i]=1; 
	efficGrid[i]=1;
	//printf ("All must be 1: %e,  %e,  %e,  %e,  %e\n",abscoOpaqueQuarz[i],abscoCsI[i],abscoGrid[i],efficAll[i],efficGrid[i]);
    } 
    
    //Efficiency for csi 
    
    Float_t efficCsI[26] = {0.000199999995, 0.000600000028, 0.000699999975, 0.00499999989, 0.00749999983, 0.010125,
			     0.0242999997, 0.0405000001, 0.0688500032, 0.105299994, 0.121500008, 0.141749993, 0.157949999,
			     0.162, 0.166050002, 0.167669997, 0.174299985, 0.176789999, 0.179279998, 0.182599992, 0.18592,
			     0.187579989, 0.189239994, 0.190899998, 0.207499996, 0.215799987};
	
    

    //FRESNEL LOSS CORRECTION FOR PERPENDICULAR INCIDENCE AND
    //UNPOLARIZED PHOTONS

    for (i=0;i<26;i++)
    {
	efficCsI[i] = efficCsI[i]/(1.-Fresnel(ppckov[i]*1e9,1.,0)); 
	//printf ("Fresnel result: %e for energy: %e\n",Fresnel(ppckov[i]*1e9,1.,0),ppckov[i]*1e9);
    }
	
    /*******************************************End of rich_media.f***************************************/

  

    
    
    
    Float_t afre[2], agri, amet[2], aqua[2], ahon, zfre[2], zgri, zhon, 
    zmet[2], zqua[2];
    Int_t nlmatfre;
    Float_t densquao;
    Int_t nlmatmet, nlmatqua;
    Float_t wmatquao[2], rIndexFreon[26];
    Float_t aquao[2], epsil, stmin, zquao[2];
    Int_t nlmatquao;
    Float_t radlal, densal, tmaxfd, deemax, stemax;
    Float_t aal, zal, radlgri, densfre, radlhon, densgri, denshon,densqua, densmet, wmatfre[2], wmatmet[2], wmatqua[2];
    
    Int_t *idtmed = fIdtmed->GetArray()-999;
    
    TGeant3 *geant3 = (TGeant3*) gMC;
    
    // --- Photon energy (GeV) 
    // --- Refraction indexes 
    for (i = 0; i < 26; ++i) {
      rIndexFreon[i] = ppckov[i] * .0172 * 1e9 + 1.177;
      //rIndexFreon[i] = 1;
	//printf ("rIndexFreon: %e \n efficCsI: %e for energy: %e\n",rIndexFreon[i], efficCsI[i], ppckov[i]);
    }
            
    // --- Detection efficiencies (quantum efficiency for CsI) 
    // --- Define parameters for honeycomb. 
    //     Used carbon of equivalent rad. lenght 
    
    ahon    = 12.01;
    zhon    = 6.;
    denshon = 2.265;
    radlhon = 18.8;
    
    // --- Parameters to include in GSMIXT, relative to Quarz (SiO2) 
    
    aqua[0]    = 28.09;
    aqua[1]    = 16.;
    zqua[0]    = 14.;
    zqua[1]    = 8.;
    densqua    = 2.64;
    nlmatqua   = -2;
    wmatqua[0] = 1.;
    wmatqua[1] = 2.;
    
    // --- Parameters to include in GSMIXT, relative to opaque Quarz (SiO2) 
    
    aquao[0]    = 28.09;
    aquao[1]    = 16.;
    zquao[0]    = 14.;
    zquao[1]    = 8.;
    densquao    = 2.64;
    nlmatquao   = -2;
    wmatquao[0] = 1.;
    wmatquao[1] = 2.;
    
    // --- Parameters to include in GSMIXT, relative to Freon (C6F14) 
    
    afre[0]    = 12.;
    afre[1]    = 19.;
    zfre[0]    = 6.;
    zfre[1]    = 9.;
    densfre    = 1.7;
    nlmatfre   = -2;
    wmatfre[0] = 6.;
    wmatfre[1] = 14.;
    
    // --- Parameters to include in GSMIXT, relative to methane (CH4) 
    
    amet[0]    = 12.01;
    amet[1]    = 1.;
    zmet[0]    = 6.;
    zmet[1]    = 1.;
    densmet    = 7.17e-4;
    nlmatmet   = -2;
    wmatmet[0] = 1.;
    wmatmet[1] = 4.;
    
    // --- Parameters to include in GSMIXT, relative to anode grid (Cu) 
  
    agri    = 63.54;
    zgri    = 29.;
    densgri = 8.96;
    radlgri = 1.43;
    
    // --- Parameters to include in GSMATE related to aluminium sheet 
    
    aal    = 26.98;
    zal    = 13.;
    densal = 2.7;
    radlal = 8.9;
    
    AliMaterial(1, "Air     $", 14.61, 7.3, .001205, 30420., 67500);
    AliMaterial(6, "HON", ahon, zhon, denshon, radlhon, 0);
    AliMaterial(16, "CSI", ahon, zhon, denshon, radlhon, 0);
    AliMixture(20, "QUA", aqua, zqua, densqua, nlmatqua, wmatqua);
    AliMixture(21, "QUAO", aquao, zquao, densquao, nlmatquao, wmatquao);
    AliMixture(30, "FRE", afre, zfre, densfre, nlmatfre, wmatfre);
    AliMixture(40, "MET", amet, zmet, densmet, nlmatmet, wmatmet);
    AliMixture(41, "METG", amet, zmet, densmet, nlmatmet, wmatmet);
    AliMaterial(11, "GRI", agri, zgri, densgri, radlgri, 0);
    AliMaterial(50, "ALUM", aal, zal, densal, radlal, 0);
    
    tmaxfd = -10.;
    stemax = -.1;
    deemax = -.2;
    epsil  = .001;
    stmin  = -.001;
    
    AliMedium(1, "DEFAULT MEDIUM AIR$", 1, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(2, "HONEYCOMB$", 6, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(3, "QUARZO$", 20, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(4, "FREON$", 30, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(5, "METANO$", 40, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(6, "CSI$", 16, 1, isxfld, sxmgmx,tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(7, "GRIGLIA$", 11, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(8, "QUARZOO$", 21, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(9, "GAP$", 41, 1, isxfld, sxmgmx,tmaxfd, .1, -deemax, epsil, -stmin);
    AliMedium(10, "ALUMINUM$", 50, 1, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
    

    geant3->Gsckov(idtmed[1000], 26, ppckov, abscoMethane, efficAll, rIndexMethane);
    geant3->Gsckov(idtmed[1001], 26, ppckov, abscoMethane, efficAll, rIndexMethane);
    geant3->Gsckov(idtmed[1002], 26, ppckov, abscoQuarz, efficAll,rIndexQuarz);
    geant3->Gsckov(idtmed[1003], 26, ppckov, abscoFreon, efficAll,rIndexFreon);
    geant3->Gsckov(idtmed[1004], 26, ppckov, abscoMethane, efficAll, rIndexMethane);
    geant3->Gsckov(idtmed[1005], 26, ppckov, abscoCsI, efficCsI, rIndexMethane);
    geant3->Gsckov(idtmed[1006], 26, ppckov, abscoGrid, efficGrid, rIndexGrid);
    geant3->Gsckov(idtmed[1007], 26, ppckov, abscoOpaqueQuarz, efficAll, rIndexOpaqueQuarz);
    geant3->Gsckov(idtmed[1008], 26, ppckov, abscoMethane, efficAll, rIndexMethane);
    geant3->Gsckov(idtmed[1009], 26, ppckov, abscoGrid, efficGrid, rIndexGrid);
}

//___________________________________________

Float_t AliRICHv0::Fresnel(Float_t ene,Float_t pdoti, Bool_t pola)
{

    //ENE(EV), PDOTI=COS(INC.ANG.), PDOTR=COS(POL.PLANE ROT.ANG.)
    
    Float_t en[36] = {5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.0,6.1,6.2,
		      6.3,6.4,6.5,6.6,6.7,6.8,6.9,7.0,7.1,7.2,7.3,7.4,7.5,7.6,7.7,
		      7.8,7.9,8.0,8.1,8.2,8.3,8.4,8.5};
     

    Float_t csin[36] = {2.14,2.21,2.33,2.48,2.76,2.97,2.99,2.59,2.81,3.05,
			2.86,2.53,2.55,2.66,2.79,2.96,3.18,3.05,2.84,2.81,2.38,2.11,
			2.01,2.13,2.39,2.73,3.08,3.15,2.95,2.73,2.56,2.41,2.12,1.95,
			1.72,1.53};
      
    Float_t csik[36] = {0.,0.,0.,0.,0.,0.196,0.408,0.208,0.118,0.49,0.784,0.543,
	 		0.424,0.404,0.371,0.514,0.922,1.102,1.139,1.376,1.461,1.253,0.878,
			0.69,0.612,0.649,0.824,1.347,1.571,1.678,1.763,1.857,1.824,1.824,
			1.714,1.498};
    Float_t xe=ene;
    Int_t  j=Int_t(xe*10)-49;
    Float_t cn=csin[j]+((csin[j+1]-csin[j])/0.1)*(xe-en[j]);
    Float_t ck=csik[j]+((csik[j+1]-csik[j])/0.1)*(xe-en[j]);

    //FORMULAE FROM HANDBOOK OF OPTICS, 33.23 OR
    //W.R. HUNTER, J.O.S.A. 54 (1964),15 , J.O.S.A. 55(1965),1197

    Float_t sinin=TMath::Sqrt(1-pdoti*pdoti);
    Float_t tanin=sinin/pdoti;

    Float_t c1=cn*cn-ck*ck-sinin*sinin;
    Float_t c2=4*cn*cn*ck*ck;
    Float_t aO=TMath::Sqrt(0.5*(TMath::Sqrt(c1*c1+c2)+c1));
    Float_t b2=0.5*(TMath::Sqrt(c1*c1+c2)-c1);
    
    Float_t rs=((aO-pdoti)*(aO-pdoti)+b2)/((aO+pdoti)*(aO+pdoti)+b2);
    Float_t rp=rs*((aO-sinin*tanin)*(aO-sinin*tanin)+b2)/((aO+sinin*tanin)*(aO+sinin*tanin)+b2);
    

    //CORRECTION FACTOR FOR SURFACE ROUGHNESS
    //B.J. STAGG  APPLIED OPTICS, 30(1991),4113

    Float_t sigraf=18.;
    Float_t lamb=1240/ene;
    Float_t fresn;
 
    Float_t  rO=TMath::Exp(-(4*TMath::Pi()*pdoti*sigraf/lamb)*(4*TMath::Pi()*pdoti*sigraf/lamb));

    if(pola)
    {
	Float_t pdotr=0.8;                                 //DEGREE OF POLARIZATION : 1->P , -1->S
	fresn=0.5*(rp*(1+pdotr)+rs*(1-pdotr));
    }
    else
	fresn=0.5*(rp+rs);
      
    fresn = fresn*rO;
    return(fresn);
}

//__________________________________________

Float_t AliRICHv0::AbsoCH4(Float_t x)
{

    //KLOSCH,SCH4(9),WL(9),EM(9),ALENGTH(31)
    Float_t sch4[9] = {.12,.16,.23,.38,.86,2.8,7.9,28.,80.};              //MB X 10^22
    //Float_t wl[9] = {153.,152.,151.,150.,149.,148.,147.,146.,145};
    Float_t em[9] = {8.1,8.158,8.212,8.267,8.322,8.378,8.435,8.493,8.55};
    const Float_t kLosch=2.686763E19;                                      // LOSCHMIDT NUMBER IN CM-3
    const Float_t kIgas1=100, kIgas2=0, kOxy=10., kWater=5., kPressure=750.,kTemperature=283.;                                      
    Float_t pn=kPressure/760.;
    Float_t tn=kTemperature/273.16;
    
	
// ------- METHANE CROSS SECTION -----------------
// ASTROPH. J. 214, L47 (1978)
	
    Float_t sm=0;
    if (x<7.75) 
	sm=.06e-22;
    
    if(x>=7.75 && x<=8.1)
    {
	Float_t c0=-1.655279e-1;
	Float_t c1=6.307392e-2;
	Float_t c2=-8.011441e-3;
	Float_t c3=3.392126e-4;
	sm=(c0+c1*x+c2*x*x+c3*x*x*x)*1.e-18;
    }
    
    if (x> 8.1)
    {
	Int_t j=0;
	while (x<=em[j] && x>=em[j+1])
	{
	    j++;
	    Float_t a=(sch4[j+1]-sch4[j])/(em[j+1]-em[j]);
	    sm=(sch4[j]+a*(x-em[j]))*1e-22;
	}
    }
    
    Float_t dm=(kIgas1/100.)*(1.-((kOxy+kWater)/1.e6))*kLosch*pn/tn;
    Float_t abslm=1./sm/dm;
    
//    ------- ISOBUTHANE CROSS SECTION --------------
//     i-C4H10 (ai) abs. length from curves in
//     Lu-McDonald paper for BARI RICH workshop .
//     -----------------------------------------------------------
    
    Float_t ai;
    Float_t absli;
    if (kIgas2 != 0) 
    {
	if (x<7.25)
	    ai=100000000.;
	
	if(x>=7.25 && x<7.375)
	    ai=24.3;
	
	if(x>=7.375)
	    ai=.0000000001;
	
	Float_t si = 1./(ai*kLosch*273.16/293.);                    // ISOB. CRO.SEC.IN CM2
	Float_t di=(kIgas2/100.)*(1.-((kOxy+kWater)/1.e6))*kLosch*pn/tn;
	absli =1./si/di;
    }
    else
	absli=1.e18;
//    ---------------------------------------------------------
//
//       transmission of O2
//
//       y= path in cm, x=energy in eV
//       so= cross section for UV absorption in cm2
//       do= O2 molecular density in cm-3
//    ---------------------------------------------------------
    
    Float_t abslo;
    Float_t so=0;
    if(x>=6.0)
    {
	if(x>=6.0 && x<6.5)
	{
	    so=3.392709e-13 * TMath::Exp(2.864104 *x);
	    so=so*1e-18;
	}
	
	if(x>=6.5 && x<7.0) 
	{
	    so=2.910039e-34 * TMath::Exp(10.3337*x);
	    so=so*1e-18;
	}
	    

	if (x>=7.0) 
	{
	    Float_t a0=-73770.76;
	    Float_t a1=46190.69;
	    Float_t a2=-11475.44;
	    Float_t a3=1412.611;
	    Float_t a4=-86.07027;
	    Float_t a5=2.074234;
	    so= a0+(a1*x)+(a2*x*x)+(a3*x*x*x)+(a4*x*x*x*x)+(a5*x*x*x*x*x);
	    so=so*1e-18;
	}
	
	Float_t dox=(kOxy/1e6)*kLosch*pn/tn;
	abslo=1./so/dox;
    }
    else
	abslo=1.e18;
//     ---------------------------------------------------------
//
//       transmission of H2O
//
//       y= path in cm, x=energy in eV
//       sw= cross section for UV absorption in cm2
//       dw= H2O molecular density in cm-3
//     ---------------------------------------------------------
    
    Float_t abslw;
    
    Float_t b0=29231.65;
    Float_t b1=-15807.74;
    Float_t b2=3192.926;
    Float_t b3=-285.4809;
    Float_t b4=9.533944;
    
    if(x>6.75)
    {    
	Float_t sw= b0+(b1*x)+(b2*x*x)+(b3*x*x*x)+(b4*x*x*x*x);
	sw=sw*1e-18;
	Float_t dw=(kWater/1e6)*kLosch*pn/tn;
	abslw=1./sw/dw;
    }
    else
    	abslw=1.e18;
	    
//    ---------------------------------------------------------
    
    Float_t alength=1./(1./abslm+1./absli+1./abslo+1./abslw);
    return (alength);
}




//___________________________________________

void AliRICHv0::Init()
{

  printf("*********************************** RICH_INIT ***********************************\n");
  printf("*                                                                               *\n");
  printf("*                       AliRICHv0 Default version started                       *\n");
  printf("*                                                                               *\n");

  
  AliRICHSegmentation*  segmentation;
  AliRICHGeometry*  geometry;
  AliRICHResponse*  response;


    // 
    // Initialize Tracking Chambers
    //
    for (Int_t i=1; i<kNCH; i++) {
	//printf ("i:%d",i);
	( (AliRICHChamber*) (*fChambers)[i])->Init();  
    }  
    
    //
    // Set the chamber (sensitive region) GEANT identifier
    
    ((AliRICHChamber*)(*fChambers)[0])->SetGid(1);  
    ((AliRICHChamber*)(*fChambers)[1])->SetGid(2);  
    ((AliRICHChamber*)(*fChambers)[2])->SetGid(3);  
    ((AliRICHChamber*)(*fChambers)[3])->SetGid(4);  
    ((AliRICHChamber*)(*fChambers)[4])->SetGid(5);  
    ((AliRICHChamber*)(*fChambers)[5])->SetGid(6);  
    ((AliRICHChamber*)(*fChambers)[6])->SetGid(7); 

    Float_t pos1[3]={0,471.8999,165.2599};
    Chamber(0).SetChamberTransform(pos1[0],pos1[1],pos1[2],new TRotMatrix("rot993","rot993",90,0,70.69,90,19.30999,-90));

    Float_t pos2[3]={171,470,0};
    Chamber(1).SetChamberTransform(pos2[0],pos2[1],pos2[2],new TRotMatrix("rot994","rot994",90,-20,90,70,0,0));

    Float_t pos3[3]={0,500,0};
    Chamber(2).SetChamberTransform(pos3[0],pos3[1],pos3[2],new TRotMatrix("rot995","rot995",90,0,90,90,0,0));
    
    Float_t pos4[3]={-171,470,0};
    Chamber(3).SetChamberTransform(pos4[0],pos4[1],pos4[2], new TRotMatrix("rot996","rot996",90,20,90,110,0,0));  

    Float_t pos5[3]={161.3999,443.3999,-165.3};
    Chamber(4).SetChamberTransform(pos5[0],pos5[1],pos5[2],new TRotMatrix("rot997","rot997",90,340,108.1999,70,18.2,70));

    Float_t pos6[3]={0., 471.9, -165.3,};
    Chamber(5).SetChamberTransform(pos6[0],pos6[1],pos6[2],new TRotMatrix("rot998","rot998",90,0,109.3099,90,19.30999,90));

    Float_t pos7[3]={-161.399,443.3999,-165.3};
    Chamber(6).SetChamberTransform(pos7[0],pos7[1],pos7[2],new TRotMatrix("rot999","rot999",90,20,108.1999,110,18.2,110));
    
    segmentation=Chamber(0).GetSegmentationModel(0);
    geometry=Chamber(0).GetGeometryModel();
    response=Chamber(0).GetResponseModel();
    
     
    printf("*                            Pads            : %3dx%3d                          *\n",segmentation->Npx(),segmentation->Npy());
    printf("*                            Pad size        : %5.2f x%5.2f mm2                 *\n",segmentation->Dpx(),segmentation->Dpy()); 
    printf("*                            Gap Thickness   : %5.1f cm                         *\n",geometry->GetGapThickness());
    printf("*                            Radiator Width  : %5.1f cm                         *\n",geometry->GetQuartzWidth());
    printf("*                            Radiator Length : %5.1f cm                         *\n",geometry->GetQuartzLength());
    printf("*                            Freon Thickness : %5.1f cm                         *\n",geometry->GetFreonThickness());
    printf("*                            Charge Slope    : %5.1f ADC                        *\n",response->ChargeSlope());
    printf("*                            Feedback Prob.  : %5.2f %%                         *\n",response->AlphaFeedback()*100);
    printf("*                                                                               *\n");
    printf("*                                   Success!                                    *\n");
    printf("*                                                                               *\n");
    printf("*********************************************************************************\n");

}

//___________________________________________
void AliRICHv0::StepManager()
{
    //Dummy step manager
}

  
//___________________________________________
