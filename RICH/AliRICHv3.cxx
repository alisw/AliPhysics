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

#include <Riostream.h>

#include <TBRIK.h>
#include <TGeometry.h>
#include <TLorentzVector.h>
#include <TNode.h>
#include <TParticle.h>
#include <TVector3.h>
#include <TVirtualMC.h>

#include "AliConst.h"
#include "AliMagF.h"
#include "AliPDG.h"
#include "AliRICHGeometry.h"
#include "AliRICHHit.h"
#include "AliRICHResponseV0.h"
#include "AliRICHSegmentationV1.h"
#include "AliRICHv3.h"
#include "AliRun.h"

ClassImp(AliRICHv3)

//______________________________________________________________
//    Implementation of the RICH version 3 with azimuthal rotation


AliRICHv3::AliRICHv3(const char *sName, const char *sTitle)
      	  :AliRICH(sName,sTitle)
{
// The named ctor currently creates a single copy of 
// AliRICHGeometry AliRICHSegmentationV1 AliRICHResponseV0
// and initialises the corresponding models of all 7 chambers with these stuctures.
// Note: all chambers share the single copy of models. MUST be changed later (???).
   cout<<ClassName()<<"::named ctor(sName,sTitle)>\n"; // no way to control it as ctor is called before call to SetDebugXXXX()

   fCkovNumber=fFreonProd=fDebugLevel=0;
   
   AliRICHGeometry       *pRICHGeometry    =new AliRICHGeometry;           // ??? to be moved to AlRICHChamber::named ctor
   AliRICHSegmentationV1 *pRICHSegmentation=new AliRICHSegmentationV1;     // ??? to be moved to AlRICHChamber::named ctor
   AliRICHResponseV0     *pRICHResponse    =new AliRICHResponseV0;         // ??? to be moved to AlRICHChamber::named ctor
     
   fChambers = new TObjArray(kNCH);
   for (Int_t i=0; i<kNCH; i++){    
      fChambers->AddAt(new AliRICHChamber,i); // ??? to be changed to named ctor of AliRICHChamber
      SetGeometryModel(i,pRICHGeometry);
      SetSegmentationModel(i,pRICHSegmentation);
      SetResponseModel(i,pRICHResponse);
      ((AliRICHChamber*)fChambers->At(i))->Init(i); // ??? to be removed     
   }
}//AliRICHv3::ctor(const char *pcName, const char *pcTitle)

AliRICHv3::~AliRICHv3()
{
// Dtor deletes RICH models. In future (???) AliRICHChamber will be responsible for that.
   if(IsDebugStart()) cout<<ClassName()<<"::dtor()>\n";
      
   delete GetChamber(0)->GetGeometryModel();
   delete GetChamber(0)->GetResponseModel();
   delete GetChamber(0)->GetSegmentationModel();
}//AliRICHv3::dtor()


void AliRICHv3::CreateGeometry()
{
// Provides geometry structure for simulation (currently GEANT volumes tree)         
   if(IsDebugStart()) cout<<ClassName()<<"::CreateGeometry()>\n";

  AliRICH *pRICH = (AliRICH *) gAlice->GetDetector("RICH"); 
  AliRICHSegmentationV0*  segmentation;
  AliRICHGeometry*  geometry;
  AliRICHChamber*       iChamber;

  iChamber = &(pRICH->Chamber(0));
  segmentation=(AliRICHSegmentationV0*) iChamber->GetSegmentationModel(0);
  geometry=iChamber->GetGeometryModel();

  Float_t distance;
  distance = geometry->GetFreonThickness()/2 + geometry->GetQuartzThickness() + geometry->GetGapThickness();
  geometry->SetRadiatorToPads(distance);
    
  //Opaque quartz thickness
  Float_t oqua_thickness = .5;
  //CsI dimensions


  Float_t csi_width = segmentation->Npx()*segmentation->Dpx() + segmentation->DeadZone();
  Float_t csi_length = segmentation->Npy()*segmentation->Dpy() + 2*segmentation->DeadZone();
  
  
  Int_t *idtmed = fIdtmed->GetArray()-999;
    
    Int_t i;
    Float_t zs;
    Int_t idrotm[1099];
    Float_t par[3];
    
    // --- Define the RICH detector 
    //     External aluminium box 
    par[0] = 68.8;
    par[1] = 13;                 //Original Settings
    par[2] = 70.86;
    gMC->Gsvolu("RICH", "BOX ", idtmed[1009], par, 3);
    
    //     Air 
    par[0] = 66.3;
    par[1] = 13;                 //Original Settings
    par[2] = 68.35;
    gMC->Gsvolu("SRIC", "BOX ", idtmed[1000], par, 3);
    
    //    Air 2 (cutting the lower part of the box)
    
    par[0] = 1.25;
    par[1] = 3;                 //Original Settings
    par[2] = 70.86;
    gMC->Gsvolu("AIR2", "BOX ", idtmed[1000], par, 3);

    //    Air 3 (cutting the lower part of the box)
    
    par[0] = 66.3;
    par[1] = 3;                 //Original Settings
    par[2] = 1.2505;
    gMC->Gsvolu("AIR3", "BOX ", idtmed[1000], par, 3);
    
    //     Honeycomb 
    par[0] = 66.3;
    par[1] = .188;                 //Original Settings
    par[2] = 68.35;
    gMC->Gsvolu("HONE", "BOX ", idtmed[1001], par, 3);
    
    //     Aluminium sheet 
    par[0] = 66.3;
    par[1] = .025;                 //Original Settings
    par[2] = 68.35;
    /*par[0] = 66.5;
    par[1] = .025;
    par[2] = 63.1;*/
    gMC->Gsvolu("ALUM", "BOX ", idtmed[1009], par, 3);
    
    //     Quartz 
    par[0] = geometry->GetQuartzWidth()/2;
    par[1] = geometry->GetQuartzThickness()/2;
    par[2] = geometry->GetQuartzLength()/2;
    gMC->Gsvolu("QUAR", "BOX ", idtmed[1002], par, 3);
    
    //     Spacers (cylinders) 
    par[0] = 0.;
    par[1] = .5;
    par[2] = geometry->GetFreonThickness()/2;
    gMC->Gsvolu("SPAC", "TUBE", idtmed[1002], par, 3);
    
    //     Feet (freon slabs supports)

    par[0] = .7;
    par[1] = .3;
    par[2] = 1.9;
    gMC->Gsvolu("FOOT", "BOX", idtmed[1009], par, 3);

    //     Opaque quartz 
    par[0] = geometry->GetQuartzWidth()/2;
    par[1] = .2;
    par[2] = geometry->GetQuartzLength()/2;
    gMC->Gsvolu("OQUA", "BOX ", idtmed[1007], par, 3);
  
    //     Frame of opaque quartz
    par[0] = geometry->GetOuterFreonWidth()/2;
    par[1] = geometry->GetFreonThickness()/2;
    par[2] = geometry->GetOuterFreonLength()/2; 
    gMC->Gsvolu("OQF1", "BOX ", idtmed[1007], par, 3);

    par[0] = geometry->GetInnerFreonWidth()/2;
    par[1] = geometry->GetFreonThickness()/2;
    par[2] = geometry->GetInnerFreonLength()/2; 
    gMC->Gsvolu("OQF2", "BOX ", idtmed[1007], par, 3);
    
    
    //     Freon 
    par[0] = geometry->GetOuterFreonWidth()/2 - oqua_thickness;
    par[1] = geometry->GetFreonThickness()/2;
    par[2] = geometry->GetOuterFreonLength()/2 - 2*oqua_thickness; 
    gMC->Gsvolu("FRE1", "BOX ", idtmed[1003], par, 3);

    par[0] = geometry->GetInnerFreonWidth()/2 - oqua_thickness;
    par[1] = geometry->GetFreonThickness()/2;
    par[2] = geometry->GetInnerFreonLength()/2 - 2*oqua_thickness; 
    gMC->Gsvolu("FRE2", "BOX ", idtmed[1003], par, 3);
    
    //     Methane 
    par[0] = csi_width/2;
    par[1] = geometry->GetGapThickness()/2;
    par[2] = csi_length/2;
    gMC->Gsvolu("META", "BOX ", idtmed[1004], par, 3);
    
    //     Methane gap 
    par[0] = csi_width/2;
    par[1] = geometry->GetProximityGapThickness()/2;
    par[2] = csi_length/2;
    gMC->Gsvolu("GAP ", "BOX ", idtmed[1008], par, 3);
    
    //     CsI photocathode 
    par[0] = csi_width/2;
    par[1] = .25;
    par[2] = csi_length/2;
    gMC->Gsvolu("CSI ", "BOX ", idtmed[1005], par, 3);
    
    //     Anode grid 
    par[0] = 0.;
    par[1] = .001;
    par[2] = 20.;
    gMC->Gsvolu("GRID", "TUBE", idtmed[1006], par, 3);

    // Wire supports
    // Bar of metal
    
    par[0] = csi_width/2;
    par[1] = 1.05;
    par[2] = 1.05;
    gMC->Gsvolu("WSMe", "BOX ", idtmed[1009], par, 3);

    // Ceramic pick up (base)
    
    par[0] =  csi_width/2;
    par[1] = .25;
    par[2] = 1.05;
    gMC->Gsvolu("WSG1", "BOX ", idtmed[1010], par, 3);

    // Ceramic pick up (head)

    par[0] = csi_width/2;
    par[1] = .1;
    par[2] = .1;
    gMC->Gsvolu("WSG2", "BOX ", idtmed[1010], par, 3);

    // Aluminium supports for methane and CsI
    // Short bar

    par[0] = csi_width/2;
    par[1] = geometry->GetGapThickness()/2 + .25;
    par[2] = (68.35 - csi_length/2)/2;
    gMC->Gsvolu("SMSH", "BOX", idtmed[1009], par, 3);
    
    // Long bar

    par[0] = (66.3 - csi_width/2)/2;
    par[1] = geometry->GetGapThickness()/2 + .25;
    par[2] = csi_length/2 + 68.35 - csi_length/2;
    gMC->Gsvolu("SMLG", "BOX", idtmed[1009], par, 3);
    
    // Aluminium supports for freon
    // Short bar

    par[0] = geometry->GetQuartzWidth()/2;
    par[1] = .3;
    par[2] = (68.35 - geometry->GetQuartzLength()/2)/2;
    gMC->Gsvolu("SFSH", "BOX", idtmed[1009], par, 3);
    
    // Long bar

    par[0] = (66.3 - geometry->GetQuartzWidth()/2)/2;
    par[1] = .3;
    par[2] = geometry->GetQuartzLength()/2 + 68.35 - geometry->GetQuartzLength()/2;
    gMC->Gsvolu("SFLG", "BOX", idtmed[1009], par, 3);
    
    // PCB backplane
    
    par[0] = csi_width/2;
    par[1] = .25;
    par[2] = csi_length/4 -.5025;
    gMC->Gsvolu("PCB ", "BOX", idtmed[1011], par, 3);

    
    // Backplane supports

    // Aluminium slab
    
    par[0] = 33.15;
    par[1] = 2;
    par[2] = 21.65;
    gMC->Gsvolu("BACK", "BOX", idtmed[1009], par, 3);
    
    // Big hole
    
    par[0] = 9.05;
    par[1] = 2;
    par[2] = 4.4625;
    gMC->Gsvolu("BKHL", "BOX", idtmed[1000], par, 3);

    // Small hole
    
    par[0] = 5.7;
    par[1] = 2;
    par[2] = 4.4625;
    gMC->Gsvolu("BKHS", "BOX", idtmed[1000], par, 3);

    // Place holes inside backplane support

    gMC->Gspos("BKHS", 1, "BACK", .8 + 5.7,0., .6 + 4.4625, 0, "ONLY");
    gMC->Gspos("BKHS", 2, "BACK", -.8 - 5.7,0., .6 + 4.4625, 0, "ONLY");
    gMC->Gspos("BKHS", 3, "BACK", .8 + 5.7,0., -.6 - 4.4625, 0, "ONLY");
    gMC->Gspos("BKHS", 4, "BACK", -.8 - 5.7,0., -.6 - 4.4625, 0, "ONLY");
    gMC->Gspos("BKHS", 5, "BACK", .8 + 5.7,0., .6 + 8.925 + 1.2 + 4.4625, 0, "ONLY");
    gMC->Gspos("BKHS", 6, "BACK", -.8 - 5.7,0., .6 + 8.925 + 1.2 + 4.4625, 0, "ONLY");
    gMC->Gspos("BKHS", 7, "BACK", .8 + 5.7,0., -.6 - 8.925 - 1.2 - 4.4625, 0, "ONLY");
    gMC->Gspos("BKHS", 8, "BACK", -.8 - 5.7,0., -.6 - 8.925 - 1.2 - 4.4625, 0, "ONLY");
    gMC->Gspos("BKHL", 1, "BACK", .8 + 11.4 + 1.6 + 9.05, 0., .6 + 4.4625, 0, "ONLY");
    gMC->Gspos("BKHL", 2, "BACK", -.8 - 11.4 - 1.6 - 9.05, 0., .6 + 4.4625, 0, "ONLY");
    gMC->Gspos("BKHL", 3, "BACK", .8 + 11.4 + 1.6 + 9.05, 0., -.6 - 4.4625, 0, "ONLY");
    gMC->Gspos("BKHL", 4, "BACK", -.8 - 11.4 - 1.6 - 9.05, 0., -.6 - 4.4625, 0, "ONLY");
    gMC->Gspos("BKHL", 5, "BACK", .8 + 11.4+ 1.6 + 9.05, 0., .6 + 8.925 + 1.2 + 4.4625, 0, "ONLY");
    gMC->Gspos("BKHL", 6, "BACK", -.8 - 11.4 - 1.6 - 9.05, 0., .6 + 8.925 + 1.2 + 4.4625, 0, "ONLY");
    gMC->Gspos("BKHL", 7, "BACK", .8 + 11.4 + 1.6 + 9.05, 0., -.6 - 8.925 - 1.2 - 4.4625, 0, "ONLY");
    gMC->Gspos("BKHL", 8, "BACK", -.8 - 11.4 - 1.6 - 9.05, 0., -.6 - 8.925 - 1.2 - 4.4625, 0, "ONLY");

    
  
    // --- Places the detectors defined with GSVOLU 
    //     Place material inside RICH 
    gMC->Gspos("SRIC", 1, "RICH", 0.,0., 0., 0, "ONLY");
    gMC->Gspos("AIR2", 1, "RICH", 66.3 + 1.2505, 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .6 - .05 - .376 -.5 - 3.35, 0., 0, "ONLY");
    gMC->Gspos("AIR2", 2, "RICH", -66.3 - 1.2505, 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .6 - .05 - .376 -.5 - 3.35, 0., 0, "ONLY");
    gMC->Gspos("AIR3", 1, "RICH", 0.,  1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .6 - .05 - .376 -.5 - 3.35, -68.35 - 1.25, 0, "ONLY");
    gMC->Gspos("AIR3", 2, "RICH", 0., 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .6 - .05 - .376 -.5 - 3.35,  68.35 + 1.25, 0, "ONLY");
    
      
    gMC->Gspos("ALUM", 1, "SRIC", 0., 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .6 - .05 - .376 -.025, 0., 0, "ONLY");
    gMC->Gspos("HONE", 1, "SRIC", 0., 1.276- geometry->GetGapThickness()/2  - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .6 - .05 - .188, 0., 0, "ONLY");
    gMC->Gspos("ALUM", 2, "SRIC", 0., 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .6 - .025, 0., 0, "ONLY");
    gMC->Gspos("FOOT", 1, "SRIC", 64.95, 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .3, 36.9, 0, "ONLY");
    gMC->Gspos("FOOT", 2, "SRIC", 21.65, 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .3 , 36.9, 0, "ONLY");
    gMC->Gspos("FOOT", 3, "SRIC", -21.65, 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .3, 36.9, 0, "ONLY");
    gMC->Gspos("FOOT", 4, "SRIC", -64.95, 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .3, 36.9, 0, "ONLY");
    gMC->Gspos("FOOT", 5, "SRIC", 64.95, 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .3, -36.9, 0, "ONLY");
    gMC->Gspos("FOOT", 6, "SRIC", 21.65, 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .3, -36.9, 0, "ONLY");
    gMC->Gspos("FOOT", 7, "SRIC", -21.65, 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .3, -36.9, 0, "ONLY");
    gMC->Gspos("FOOT", 8, "SRIC", -64.95, 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .4 - .3, -36.9, 0, "ONLY");
    gMC->Gspos("OQUA", 1, "SRIC", 0., 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()- .2, 0., 0, "ONLY");
    
    // Supports placing

    // Methane supports
    gMC->Gspos("SMLG", 1, "SRIC", csi_width/2 + (66.3 - csi_width/2)/2, 1.276 + .25, 0., 0, "ONLY");
    gMC->Gspos("SMLG", 2, "SRIC", - csi_width/2 - (66.3 - csi_width/2)/2, 1.276 + .25, 0., 0, "ONLY");
    gMC->Gspos("SMSH", 1, "SRIC", 0., 1.276 + .25, csi_length/2 + (68.35 - csi_length/2)/2, 0, "ONLY");
    gMC->Gspos("SMSH", 2, "SRIC", 0., 1.276 + .25, - csi_length/2 - (68.35 - csi_length/2)/2, 0, "ONLY");

    //Freon supports

    Float_t supp_y = 1.276 - geometry->GetGapThickness()/2- geometry->GetQuartzThickness() -geometry->GetFreonThickness() - .2 + .3; //y position of freon supports

    gMC->Gspos("SFLG", 1, "SRIC", geometry->GetQuartzWidth()/2 + (66.3 - geometry->GetQuartzWidth()/2)/2, supp_y, 0., 0, "ONLY");
    gMC->Gspos("SFLG", 2, "SRIC", - geometry->GetQuartzWidth()/2 - (66.3 - geometry->GetQuartzWidth()/2)/2, supp_y, 0., 0, "ONLY");
    gMC->Gspos("SFSH", 1, "SRIC", 0., supp_y, geometry->GetQuartzLength()/2 + (68.35 - geometry->GetQuartzLength()/2)/2, 0, "ONLY");
    gMC->Gspos("SFSH", 2, "SRIC", 0., supp_y, - geometry->GetQuartzLength()/2 - (68.35 - geometry->GetQuartzLength()/2)/2, 0, "ONLY");
    
    AliMatrix(idrotm[1019], 0., 0., 90., 0., 90., 90.);
    

    Int_t nspacers = 30;
    
    for (i = 0; i < nspacers/3; i++) {
	zs = -11.6/2 + (TMath::Abs(nspacers/6) - i) * 12.2;
	gMC->Gspos("SPAC", i, "FRE1", 10.5, 0., zs, idrotm[1019], "ONLY");  //Original settings 
    }
    
    for (i = nspacers/3; i < (nspacers*2)/3; i++) {
	zs = -11.6/2 + (nspacers/3 + TMath::Abs(nspacers/6) - i) * 12.2;
	gMC->Gspos("SPAC", i, "FRE1", 0, 0., zs, idrotm[1019], "ONLY");  //Original settings 
    }
    
    for (i = (nspacers*2)/3; i < nspacers; ++i) {
	zs = -11.6/2 + ((nspacers*2)/3 + TMath::Abs(nspacers/6) - i) * 12.2;
	gMC->Gspos("SPAC", i, "FRE1", -10.5, 0., zs, idrotm[1019], "ONLY"); //Original settings  
    }

    for (i = 0; i < nspacers/3; i++) {
	zs = -11.6/2 + (TMath::Abs(nspacers/6) - i) * 12.2;
	gMC->Gspos("SPAC", i, "FRE2", 10.5, 0., zs, idrotm[1019], "ONLY");  //Original settings 
    }
    
    for (i = nspacers/3; i < (nspacers*2)/3; i++) {
	zs = -11.6/2 + (nspacers/3 + TMath::Abs(nspacers/6) - i) * 12.2;
	gMC->Gspos("SPAC", i, "FRE2", 0, 0., zs, idrotm[1019], "ONLY");  //Original settings 
    }
    
    for (i = (nspacers*2)/3; i < nspacers; ++i) {
	zs = -11.6/2 + ((nspacers*2)/3 + TMath::Abs(nspacers/6) - i) * 12.2;
	gMC->Gspos("SPAC", i, "FRE2", -10.5, 0., zs, idrotm[1019], "ONLY"); //Original settings  
    }

    
    gMC->Gspos("FRE1", 1, "OQF1", 0., 0., 0., 0, "ONLY");
    gMC->Gspos("FRE2", 1, "OQF2", 0., 0., 0., 0, "ONLY");
    gMC->Gspos("OQF1", 1, "SRIC", geometry->GetOuterFreonWidth()/2 + geometry->GetInnerFreonWidth()/2 + 2, 1.276 - geometry->GetGapThickness()/2- geometry->GetQuartzThickness() -geometry->GetFreonThickness()/2, 0., 0, "ONLY"); //Original settings (31.3)
    gMC->Gspos("OQF2", 2, "SRIC", 0., 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()/2, 0., 0, "ONLY");          //Original settings 
    gMC->Gspos("OQF1", 3, "SRIC", - (geometry->GetOuterFreonWidth()/2 + geometry->GetInnerFreonWidth()/2) - 2, 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness() - geometry->GetFreonThickness()/2, 0., 0, "ONLY");       //Original settings (-31.3)
    gMC->Gspos("QUAR", 1, "SRIC", 0., 1.276 - geometry->GetGapThickness()/2 - geometry->GetQuartzThickness()/2, 0., 0, "ONLY");
    gMC->Gspos("GAP ", 1, "META", 0., geometry->GetGapThickness()/2 - geometry->GetProximityGapThickness()/2 - 0.0001, 0., 0, "ONLY");
    gMC->Gspos("META", 1, "SRIC", 0., 1.276, 0., 0, "ONLY");
    gMC->Gspos("CSI ", 1, "SRIC", 0., 1.276 + geometry->GetGapThickness()/2 + .25, 0., 0, "ONLY");
    printf("CSI pos: %f\n",1.276 + geometry->GetGapThickness()/2 + .25);
   
    // Wire support placing

    gMC->Gspos("WSG2", 1, "GAP ", 0., geometry->GetProximityGapThickness()/2 - .1, 0., 0, "ONLY");
    gMC->Gspos("WSG1", 1, "CSI ", 0., 0., 0., 0, "ONLY");
    gMC->Gspos("WSMe", 1, "SRIC ", 0., 1.276 + geometry->GetGapThickness()/2 + .5 + 1.05, 0., 0, "ONLY");

    // Backplane placing
    
    gMC->Gspos("BACK", 1, "SRIC ", -33.15, 1.276 + geometry->GetGapThickness()/2 + .5 + 2.1 + 2, 43.3, 0, "ONLY");
    gMC->Gspos("BACK", 2, "SRIC ", 33.15, 1.276 + geometry->GetGapThickness()/2 + .5 + 2.1 + 2 , 43.3, 0, "ONLY");
    gMC->Gspos("BACK", 3, "SRIC ", -33.15, 1.276 + geometry->GetGapThickness()/2 + .5 + 2.1 + 2, 0., 0, "ONLY");
    gMC->Gspos("BACK", 4, "SRIC ", 33.15, 1.276 + geometry->GetGapThickness()/2 + .5 + 2.1 + 2, 0., 0, "ONLY");
    gMC->Gspos("BACK", 5, "SRIC ", 33.15, 1.276 + geometry->GetGapThickness()/2 + .5 + 2.1 + 2, -43.3, 0, "ONLY");
    gMC->Gspos("BACK", 6, "SRIC ", -33.15, 1.276 + geometry->GetGapThickness()/2 + .5 + 2.1 + 2, -43.3, 0, "ONLY");

    // PCB placing
    
    gMC->Gspos("PCB ", 1, "SRIC ", 0.,  1.276 + geometry->GetGapThickness()/2 + .5 + 1.05, csi_width/4 + .5025 + 2.5, 0, "ONLY");
    gMC->Gspos("PCB ", 2, "SRIC ", 0.,  1.276 + geometry->GetGapThickness()/2 + .5 + 1.05, -csi_width/4 - .5025 - 2.5, 0, "ONLY");

// Place chambers into mother volume ALIC
           
   Double_t dOffset        = geometry->GetOffset() - geometry->GetGapThickness()/2;  // distance from center of mother volume ALIC to methane
   
   Double_t dAlpha         = geometry->GetAlphaAngle(); // angle between centers of chambers - y-z plane
   Double_t dAlphaRad      = dAlpha*kDegrad;
   
   Double_t dBeta          = geometry->GetBetaAngle();   // angle between center of chambers - y-x plane
   Double_t dBetaRad       = dBeta*kDegrad;
   
   Double_t dRotAngle      = geometry->GetRotationAngle();     // the whole RICH is to be rotated in x-y plane + means clockwise rotation 
   Double_t dRotAngleRad   = dRotAngle*kDegrad;
    
   
   TRotMatrix *pRotMatrix; // tmp pointer
   
   TVector3 vector(0,dOffset,0); // Position of chamber 2 without rotation
    
// Chamber 0  standalone (no other chambers in this row) 
   pRotMatrix = new TRotMatrix("rot993","rot993", 0., 0., 0.,0.,0.,0.);
   const Double_t* r   = pRotMatrix->SetAngles(90., 0., 90.-dAlpha , 90.,  dAlpha, -90.);
   Double_t* rr  = RotateXY(r, -dRotAngleRad);
   AliMatrix(idrotm[1000], rr[0], rr[1], rr[2], rr[3], rr[4], rr[5]);
   pRotMatrix->SetAngles(rr[0], rr[1], rr[2], rr[3], rr[4], rr[5]);

   vector.SetXYZ(0,dOffset,0);  vector.RotateX(dAlphaRad); 
   vector.RotateZ(-dRotAngleRad);
   
   gMC->Gspos("RICH",1,"ALIC",vector.X(),vector.Y(),vector.Z(),idrotm[1000], "ONLY");           
   Chamber(0).SetChamberTransform(vector.X(),vector.Y(),vector.Z(),pRotMatrix);
  if(GetDebug()) Info("CreateGeometry 0","%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f",rr[0],rr[1],rr[2],rr[3],rr[4],rr[5]);
  if(GetDebug()) Info("CreateGeometry 0","x=%8.3f y=%8.3f z=%8.3f",vector.X(),vector.Y(),vector.Z());   
// Chamber 1   
   pRotMatrix = new TRotMatrix("rot994","rot994", 0., 0., 0.,0.,0.,0.);
   r   = pRotMatrix->SetAngles(90., -dBeta, 90., 90.-dBeta,  0., 0.);
   rr  = RotateXY(r, -dRotAngleRad);
   AliMatrix(idrotm[1001], rr[0], rr[1], rr[2], rr[3], rr[4], rr[5]);
   pRotMatrix->SetAngles(rr[0], rr[1], rr[2], rr[3], rr[4], rr[5]);
   vector.SetXYZ(0,dOffset,0);  vector.RotateZ(-dBetaRad); 
   vector.RotateZ(-dRotAngleRad);
   
   gMC->Gspos("RICH",2,"ALIC",vector.X(),vector.Y(),vector.Z(),idrotm[1001], "ONLY");           
   Chamber(1).SetChamberTransform(vector.X(),vector.Y(),vector.Z(),pRotMatrix);
  if(GetDebug()) Info("CreateGeometry 1","%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f",rr[0],rr[1],rr[2],rr[3],rr[4],rr[5]);
  if(GetDebug()) Info("CreateGeometry 1","x=%8.3f y=%8.3f z=%8.3f",vector.X(),vector.Y(),vector.Z());
// Chamber 2   the top one with no Alpha-Beta rotation
   pRotMatrix = new TRotMatrix("rot995","rot995", 0., 0., 0.,0.,0.,0.);
   r   = pRotMatrix->SetAngles(90., 0., 90., 90.,  0., 0.);
   rr  = RotateXY(r, -dRotAngleRad);
   AliMatrix(idrotm[1002], rr[0], rr[1], rr[2], rr[3], rr[4], rr[5]);
   pRotMatrix->SetAngles(rr[0], rr[1], rr[2], rr[3], rr[4], rr[5]);
   vector.SetXYZ(0,dOffset,0);
   vector.RotateZ(-dRotAngleRad);
   gMC->Gspos("RICH",3,"ALIC",vector.X(),vector.Y(),vector.Z(),idrotm[1002], "ONLY");           
   Chamber(2).SetChamberTransform(vector.X(),vector.Y(),vector.Z(),pRotMatrix);
  if(GetDebug()) Info("CreateGeometry 2","%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f",rr[0],rr[1],rr[2],rr[3],rr[4],rr[5]);
  if(GetDebug()) Info("CreateGeometry 2","x=%8.3f y=%8.3f z=%8.3f",vector.X(),vector.Y(),vector.Z());   
// Chamber 3
   pRotMatrix = new TRotMatrix("rot996","rot996", 0., 0., 0.,0.,0.,0.);
   r   = pRotMatrix->SetAngles(90., dBeta, 90., 90.+dBeta,  0., 0.);
   rr  = RotateXY(r, -dRotAngleRad);
   AliMatrix(idrotm[1003], rr[0], rr[1], rr[2], rr[3], rr[4], rr[5]);
   pRotMatrix->SetAngles(rr[0], rr[1], rr[2], rr[3], rr[4], rr[5]);
   vector.SetXYZ(0,dOffset,0);  vector.RotateZ(dBetaRad); 
   vector.RotateZ(-dRotAngleRad);
   
   gMC->Gspos("RICH",4,"ALIC",vector.X(),vector.Y(),vector.Z(),idrotm[1003], "ONLY");           
   Chamber(3).SetChamberTransform(vector.X(),vector.Y(),vector.Z(),pRotMatrix);
  if(GetDebug()) Info("CreateGeometry 3","%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f",rr[0],rr[1],rr[2],rr[3],rr[4],rr[5]);
  if(GetDebug()) Info("CreateGeometry 3","x=%8.3f y=%8.3f z=%8.3f",vector.X(),vector.Y(),vector.Z());
// Chamber 4   
   pRotMatrix = new TRotMatrix("rot997","rot997", 0., 0., 0.,0.,0.,0.);
   r   = pRotMatrix->SetAngles(90., 360.-dBeta, 108.2, 90.-dBeta,  18.2, 90.-dBeta);
   rr  = RotateXY(r, -dRotAngleRad);
   AliMatrix(idrotm[1004], rr[0], rr[1], rr[2], rr[3], rr[4], rr[5]);
   pRotMatrix->SetAngles(rr[0], rr[1], rr[2], rr[3], rr[4], rr[5]);
   vector.SetXYZ(0,dOffset,0);  vector.RotateZ(-dBetaRad); vector.RotateX(-dAlphaRad); 
   vector.RotateZ(-dRotAngleRad);
   
   gMC->Gspos("RICH",5,"ALIC",vector.X(),vector.Y(),vector.Z(),idrotm[1004], "ONLY");
   Chamber(4).SetChamberTransform(vector.X(),vector.Y(),vector.Z(),pRotMatrix);
  if(GetDebug()) Info("CreateGeometry 4","%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f",rr[0],rr[1],rr[2],rr[3],rr[4],rr[5]);
  if(GetDebug()) Info("CreateGeometry 4","x=%8.3f y=%8.3f z=%8.3f",vector.X(),vector.Y(),vector.Z());
// Chamber 5   
   pRotMatrix = new TRotMatrix("rot998","rot998", 0., 0., 0.,0.,0.,0.);
   r   = pRotMatrix->SetAngles(90., 0., 90.+dAlpha, 90.,  dAlpha, 90.);
   rr  = RotateXY(r, -dRotAngleRad);
   AliMatrix(idrotm[1005], rr[0], rr[1], rr[2], rr[3], rr[4], rr[5]);
   pRotMatrix->SetAngles(rr[0], rr[1], rr[2], rr[3], rr[4], rr[5]);   
   vector.SetXYZ(0,dOffset,0); vector.RotateX(-dAlphaRad); 
   vector.RotateZ(-dRotAngleRad);
      
   gMC->Gspos("RICH",6,"ALIC",vector.X(),vector.Y(),vector.Z(),idrotm[1005], "ONLY");           
   Chamber(5).SetChamberTransform(vector.X(),vector.Y(),vector.Z(),pRotMatrix);
  if(GetDebug()) Info("CreateGeometry 5","%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f",rr[0],rr[1],rr[2],rr[3],rr[4],rr[5]);
  if(GetDebug()) Info("CreateGeometry 5","x=%8.3f y=%8.3f z=%8.3f",vector.X(),vector.Y(),vector.Z());
// Chamber 6          
   pRotMatrix = new TRotMatrix("rot999","rot999", 0., 0., 0.,0.,0.,0.);
   r   = pRotMatrix->SetAngles(90., dBeta, 108.2, 90.+dBeta,  18.2, 90.+dBeta);
   rr  = RotateXY(r, -dRotAngleRad);
   AliMatrix(idrotm[1006], rr[0], rr[1], rr[2], rr[3], rr[4], rr[5]);
   pRotMatrix->SetAngles(rr[0], rr[1], rr[2], rr[3], rr[4], rr[5]);
   vector.SetXYZ(0,dOffset,0);  vector.RotateZ(dBetaRad); vector.RotateX(-dAlphaRad); 
   vector.RotateZ(-dRotAngleRad);
      
   gMC->Gspos("RICH",7,"ALIC",vector.X(),vector.Y(),vector.Z(),idrotm[1006], "ONLY");
   Chamber(6).SetChamberTransform(vector.X(),vector.Y(),vector.Z(),pRotMatrix);
  if(GetDebug()) Info("CreateGeometry 6","%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f",rr[0],rr[1],rr[2],rr[3],rr[4],rr[5]);
  if(GetDebug()) Info("CreateGeometry 6","x=%8.3f y=%8.3f z=%8.3f",vector.X(),vector.Y(),vector.Z());
      
}//void AliRICHv3::CreateGeometry()
//______________________________________________________________________________
void AliRICHv3::Init()
{
// Makes nothing for a while   
   if(IsDebugStart()) cout<<ClassName()<<"::Init()>\n";
    
}
//______________________________________________________________________________
void AliRICHv3::BuildGeometry()    
{                                      
// Provides geometry structure for event display (ROOT TNode tree)
   
   if(IsDebugStart()) cout<<ClassName()<<"::BuildGeometry()>\n";
  
    TNode *node, *subnode, *top;
    
    const int kColorRICH = kRed;
    //
    top=gAlice->GetGeometry()->GetNode("alice");

    AliRICH *pRICH = (AliRICH *) gAlice->GetDetector("RICH"); 
    AliRICHChamber*       iChamber;
    AliRICHGeometry*  geometry;
 
    iChamber = &(pRICH->Chamber(0));
    AliRICHSegmentationV1* segmentation=(AliRICHSegmentationV1*) iChamber->GetSegmentationModel(0);
    geometry=iChamber->GetGeometryModel();
    
    new TBRIK("S_RICH","S_RICH","void",71.09999,11.5,73.15);

    Float_t padplane_width = segmentation->GetPadPlaneWidth();
    Float_t padplane_length = segmentation->GetPadPlaneLength();


    new TBRIK("PHOTO","PHOTO","void", padplane_width/2,.1,padplane_length/2);

// Chamber 0             
    top->cd();
    node = new TNode("RICH1","RICH1","S_RICH",Chamber(0).GetX(),Chamber(0).GetY(),Chamber(0).GetZ(),"rot993");
    node->SetLineColor(kColorRICH);
    node->cd();
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    fNodes->Add(node);

// Chamber 1
    top->cd(); 
    node = new TNode("RICH2","RICH2","S_RICH",Chamber(1).GetX(),Chamber(1).GetY(),Chamber(1).GetZ(),"rot994");
    node->SetLineColor(kColorRICH);
    node->cd();
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    fNodes->Add(node);

// Chamber 2
    top->cd();
    node = new TNode("RICH3","RICH3","S_RICH",Chamber(2).GetX(),Chamber(2).GetY(),Chamber(2).GetZ(),"rot995");
    node->SetLineColor(kColorRICH);
    node->cd();
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    fNodes->Add(node);
    
// Chamber 3
    top->cd();
    node = new TNode("RICH4","RICH4","S_RICH",Chamber(3).GetX(),Chamber(3).GetY(),Chamber(3).GetZ(),"rot996");
    node->SetLineColor(kColorRICH);
    node->cd();
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    fNodes->Add(node);

// Chamber 4
    top->cd();
    node = new TNode("RICH5","RICH5","S_RICH",Chamber(4).GetX(),Chamber(4).GetY(),Chamber(4).GetZ(),"rot997");
    node->SetLineColor(kColorRICH);
    node->cd();
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    fNodes->Add(node);

// Chamber 5
    top->cd();
    node = new TNode("RICH6","RICH6","S_RICH",Chamber(5).GetX(),Chamber(5).GetY(),Chamber(5).GetZ(),"rot998");
    node->SetLineColor(kColorRICH);
    fNodes->Add(node);node->cd();
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);

// Chamber 6
    top->cd();
    node = new TNode("RICH7","RICH7","S_RICH",Chamber(6).GetX(),Chamber(6).GetY(),Chamber(6).GetZ(),"rot999");
    node->SetLineColor(kColorRICH);
    node->cd();
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,padplane_length/2 + segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",padplane_width + segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",0,5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    subnode = new TNode("PHOTO1","PHOTO1","PHOTO",-padplane_width - segmentation->DeadZone(),5,-padplane_length/2 - segmentation->DeadZone()/2,"");
    subnode->SetLineColor(kGreen);
    fNodes->Add(subnode);
    fNodes->Add(node); 
    
}//AliRICHv3::BuildGeometry()
//______________________________________________________________________________
Double_t* AliRICHv3::RotateXY(const Double_t* r, Double_t a)
{
    // Rotatation in xy-plane
    // by angle a
    // The resulting rotation matrix is given back in the G3 notation. 
    Double_t* rr = new Double_t[6];
    Double_t m[9];
    Int_t i,j,k;
    
    for (i = 0; i < 3; i++) {
	j = 3*i;
	m[j]   = r[j] * TMath::Cos(a) - r[j+1] * TMath::Sin(a);
	m[j+1] = r[j] * TMath::Sin(a) + r[j+1] * TMath::Cos(a);
	m[j+2] = r[j+2];
    }
    
    for (i = 0; i < 3; i++) {
	    j = 3*i;
	    k = 2*i;
	    rr[k]    = TMath::ACos(m[j+2])        * kRaddeg;
	    rr[k+1]  = TMath::ATan2(m[j+1], m[j]) * kRaddeg;
    }
    return rr;
}//Double_t* AliRICHv3::RotateXY(const Double_t* r, Double_t a)
