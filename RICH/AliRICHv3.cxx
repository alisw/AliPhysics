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


#include <Riostream.h>

#include <TBRIK.h>
#include <TGeometry.h>
#include <TLorentzVector.h>
#include <TNode.h>
#include <TParticle.h>
#include <TVector3.h>
#include <TVirtualMC.h>
#include <TPDGCode.h> //for kNuetron
#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>

#include "AliConst.h"
#include "AliMagF.h"
#include "AliPDG.h"
#include "AliRICHGeometry.h"
#include "AliRICHResponseV0.h"
#include "AliRICHSegmentationV1.h"
#include "AliRICHv3.h"
#include "AliRun.h"
#include "AliRICHRecHit3D.h"
#include "AliRICHRawCluster.h"
#include "AliRICHDigit.h"
#include "AliRICHRecHit1D.h"


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
  if(GetDebug())Info("named ctor","Start.");

   fCkovNumber=fFreonProd=0;
   
   AliRICHGeometry       *pRICHGeometry    =new AliRICHGeometry;           // ??? to be moved to AlRICHChamber::named ctor
   AliRICHSegmentationV1 *pRICHSegmentation=new AliRICHSegmentationV1;     // ??? to be moved to AlRICHChamber::named ctor
   AliRICHResponseV0     *pRICHResponse    =new AliRICHResponseV0;         // ??? to be moved to AlRICHChamber::named ctor
     
   for (Int_t i=1; i<=kNCH; i++){    
      SetGeometryModel(i,pRICHGeometry);
      SetSegmentationModel(i,pRICHSegmentation);
      SetResponseModel(i,pRICHResponse);
      C(i)->Init(i); // ??? to be removed     
   }
  if(GetDebug())Info("named ctor","Stop.");
}//AliRICHv3::ctor(const char *pcName, const char *pcTitle)

AliRICHv3::~AliRICHv3()
{
// Dtor deletes RICH models. In future (???) AliRICHChamber will be responsible for that.
   if(GetDebug()) cout<<ClassName()<<"::dtor()>\n";
      
   if(fChambers) {
     AliRICHChamber *ch =C(1); 
     if(ch) {
       delete ch->GetGeometryModel();
       delete ch->GetResponseModel();
       delete ch->GetSegmentationModel();
     }
   }
}//AliRICHv3::dtor()


void AliRICHv3::CreateGeometry()
{
// Provides geometry structure for simulation (currently GEANT volumes tree)         
   if(GetDebug()) cout<<ClassName()<<"::CreateGeometry()>\n";

  AliRICH *pRICH = (AliRICH *) gAlice->GetDetector("RICH"); 
  AliRICHSegmentationV0*  segmentation;
  AliRICHGeometry*  geometry;
  AliRICHChamber*       iChamber;

  iChamber = &(pRICH->Chamber(0));
  segmentation=(AliRICHSegmentationV0*) iChamber->GetSegmentationModel();
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
{//Makes nothing for a while   
  if(GetDebug())Info("Init","Start.");
  if(GetDebug())Info("Init","Stop.");    
}
//______________________________________________________________________________
void AliRICHv3::BuildGeometry()    
{//Provides geometry structure for event display (ROOT TNode tree)
  if(GetDebug())Info("BuildGeometry","Start.");
  
    TNode *node, *subnode, *top;
    
    const int kColorRICH = kRed;
    //
    top=gAlice->GetGeometry()->GetNode("alice");

    AliRICH *pRICH = (AliRICH *) gAlice->GetDetector("RICH"); 
    AliRICHChamber*       iChamber;
    AliRICHGeometry*  geometry;
 
    iChamber = &(pRICH->Chamber(0));
    AliRICHSegmentationV1* segmentation=(AliRICHSegmentationV1*) iChamber->GetSegmentationModel();
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
  if(GetDebug())Info("BuildGeometry","Stop.");    
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
//______________________________________________________________________________
void AliRICHv3::StepManager()
{//Full Step Manager

    Int_t          copy, id;
    static Int_t   idvol;
    static Int_t   vol[2];
    Int_t          ipart;
    static Float_t hits[22];
    static Float_t ckovData[19];
    TLorentzVector position;
    TLorentzVector momentum;
    Float_t        pos[3];
    Float_t        mom[4];
    Float_t        localPos[3];
    Float_t        localMom[4];
    Float_t        localTheta,localPhi;
    Float_t        theta,phi;
    Float_t        destep, step;
    Double_t        ranf[2];
    Float_t        coscerenkov;
    static Float_t eloss, xhit, yhit, tlength;
    const  Float_t kBig=1.e10;
       
    TClonesArray &lhits = *fHits;
    TParticle *current = (TParticle*)(*gAlice->Particles())[gAlice->GetCurrentTrackNumber()];

 //if (current->Energy()>1)
   //{
        
    // Only gas gap inside chamber
    // Tag chambers and record hits when track enters 
    
 
    id=gMC->CurrentVolID(copy);
    idvol = copy-1;
    Float_t cherenkovLoss=0;
    //gAlice->KeepTrack(gAlice->GetCurrentTrackNumber());
    
    gMC->TrackPosition(position);
    pos[0]=position(0);
    pos[1]=position(1);
    pos[2]=position(2);
    //bzero((char *)ckovData,sizeof(ckovData)*19);
    ckovData[1] = pos[0];                 // X-position for hit
    ckovData[2] = pos[1];                 // Y-position for hit
    ckovData[3] = pos[2];                 // Z-position for hit
    ckovData[6] = 0;                      // dummy track length
    //ckovData[11] = gAlice->GetCurrentTrackNumber();
    
    //printf("\n+++++++++++\nTrack: %d\n++++++++++++\n",gAlice->GetCurrentTrackNumber());

    //AliRICH *RICH = (AliRICH *) gAlice->GetDetector("RICH"); 
    
    /********************Store production parameters for Cerenkov photons************************/ 
//is it a Cerenkov photon? 
    if (gMC->TrackPid() == 50000050) { 

      //if (gMC->VolId("GAP ")==gMC->CurrentVolID(copy))
        //{                    
	  Float_t ckovEnergy = current->Energy();
	  //energy interval for tracking
	  if  (ckovEnergy > 5.6e-09 && ckovEnergy < 7.8e-09 )       
	    //if (ckovEnergy > 0)
	    {
	      if (gMC->IsTrackEntering()){        //is track entering?
		//printf("Track entered (1)\n");
		if (gMC->VolId("FRE1")==gMC->CurrentVolID(copy) || gMC->VolId("FRE2")==gMC->CurrentVolID(copy))
		  {                                                          //is it in freo?
		    if (gMC->IsNewTrack()){                          //is it the first step?
		      //printf("I'm in!\n");
		      Int_t mother = current->GetFirstMother(); 
		      
		      //printf("Second Mother:%d\n",current->GetSecondMother());
		      
		      ckovData[10] = mother;
		      ckovData[11] = gAlice->GetCurrentTrackNumber();
		      ckovData[12] = 1;             //Media where photon was produced 1->Freon, 2->Quarz
		      //printf("Produced in FREO\n");
		      fCkovNumber++;
		      fFreonProd=1;
		      //printf("Index: %d\n",fCkovNumber);
		    }    //first step question
		  }        //freo question
		
		if (gMC->IsNewTrack()){                                  //is it first step?
		  if (gMC->VolId("QUAR")==gMC->CurrentVolID(copy))             //is it in quarz?
		    {
		      ckovData[12] = 2;
		      //printf("Produced in QUAR\n");
		    }    //quarz question
		}        //first step question
		
		//printf("Before %d\n",fFreonProd);
	      }   //track entering question
	      
	      if (ckovData[12] == 1)                                        //was it produced in Freon?
		//if (fFreonProd == 1)
		{
		  if (gMC->IsTrackEntering()){                                     //is track entering?
		    //printf("Track entered (2)\n");
		    //printf("Current volume (should be META): %s\n",gMC->CurrentVolName());
		    //printf("VolId: %d, CurrentVolID: %d\n",gMC->VolId("META"),gMC->CurrentVolID(copy));
		    if (gMC->VolId("META")==gMC->CurrentVolID(copy))                //is it in gap?      
		      {
			//printf("Got in META\n");
			gMC->TrackMomentum(momentum);
			mom[0]=momentum(0);
			mom[1]=momentum(1);
			mom[2]=momentum(2);
			mom[3]=momentum(3);
			
			gMC->Gmtod(mom,localMom,2);
			Float_t cophi = TMath::Cos(TMath::ATan2(localMom[0], localMom[1]));
			Float_t t = (1. - .025 / cophi) * (1. - .05 /  cophi);
			/**************** Photons lost in second grid have to be calculated by hand************/ 
			gMC->GetRandom()->RndmArray(1,ranf);
			if (ranf[0] > t) {
			  gMC->StopTrack();
			  ckovData[13] = 5;
			  AddCerenkov(gAlice->GetCurrentTrackNumber(),vol,ckovData);
			  //printf("Added One (1)!\n");
			  //printf("Lost one in grid\n");
			}
			/**********************************************************************************/
		      }    //gap
		    
		    //printf("Current volume (should be CSI) (1): %s\n",gMC->CurrentVolName());
		    //printf("VolId: %d, CurrentVolID: %d\n",gMC->VolId("CSI "),gMC->CurrentVolID(copy));
		    if (gMC->VolId("CSI ")==gMC->CurrentVolID(copy))             //is it in csi?      
		      {
			//printf("Got in CSI\n");
			gMC->TrackMomentum(momentum);
			mom[0]=momentum(0);
			mom[1]=momentum(1);
			mom[2]=momentum(2);
			mom[3]=momentum(3);

			gMC->Gmtod(mom,localMom,2);
			/********* Photons lost by Fresnel reflection have to be calculated by hand********/ 
			/***********************Cerenkov phtons (always polarised)*************************/
			Double_t localTc = localMom[0]*localMom[0]+localMom[2]*localMom[2];
			Double_t localRt = TMath::Sqrt(localTc);
			localTheta   = Float_t(TMath::ATan2(localRt,Double_t(localMom[1])));
			Double_t cotheta = TMath::Abs(cos(localTheta));
			Float_t t = Fresnel(ckovEnergy*1e9,cotheta,1);
			    gMC->GetRandom()->RndmArray(1,ranf);
			    if (ranf[0] < t) {
			      gMC->StopTrack();
			      ckovData[13] = 6;
			      AddCerenkov(gAlice->GetCurrentTrackNumber(),vol,ckovData);
				
			      //printf("Added One (2)!\n");
			      //printf("Lost by Fresnel\n");
			    }
			    /**********************************************************************************/
		      }
		  } //track entering?
		  
		  
		  /********************Evaluation of losses************************/
		  /******************still in the old fashion**********************/
		  
		  TArrayI procs;
		  Int_t i1 = gMC->StepProcesses(procs);            //number of physics mechanisms acting on the particle
		  for (Int_t i = 0; i < i1; ++i) {
		    //        Reflection loss 
		    if (procs[i] == kPLightReflection) {        //was it reflected
		      ckovData[13]=10;
		      if (gMC->VolId("FRE1")==gMC->CurrentVolID(copy) || gMC->VolId("FRE2")==gMC->CurrentVolID(copy)) 
			ckovData[13]=1;
		      if (gMC->CurrentVolID(copy) == gMC->VolId("QUAR")) 
			ckovData[13]=2;
		      //gMC->StopTrack();
		      //AddCerenkov(gAlice->GetCurrentTrackNumber(),vol,ckovData);
		    } //reflection question
		     
		    //        Absorption loss 
		    else if (procs[i] == kPLightAbsorption) {              //was it absorbed?
		      //printf("Got in absorption\n");
		      ckovData[13]=20;
		      if (gMC->VolId("FRE1")==gMC->CurrentVolID(copy) || gMC->VolId("FRE2")==gMC->CurrentVolID(copy)) 
			ckovData[13]=11;
		      if (gMC->CurrentVolID(copy) == gMC->VolId("QUAR")) 
			ckovData[13]=12;
		      if (gMC->CurrentVolID(copy) == gMC->VolId("META")) 
			ckovData[13]=13;
		      if (gMC->CurrentVolID(copy) == gMC->VolId("GAP ")) 
			ckovData[13]=13;
		      
		      if (gMC->CurrentVolID(copy) == gMC->VolId("SRIC")) 
			ckovData[13]=15;
		      
		      //        CsI inefficiency 
		      if (gMC->CurrentVolID(copy) == gMC->VolId("CSI ")) {
			ckovData[13]=16;
		      }
		      gMC->StopTrack();
		      AddCerenkov(gAlice->GetCurrentTrackNumber(),vol,ckovData);
		      //printf("Added One (3)!\n");
		      //printf("Added cerenkov %d\n",fCkovNumber);
		    } //absorption question 
		    
		    
		    //        Photon goes out of tracking scope 
		    else if (procs[i] == kPStop) {                 //is it below energy treshold?
		      ckovData[13]=21;
		      gMC->StopTrack();
		      AddCerenkov(gAlice->GetCurrentTrackNumber(),vol,ckovData);
		      //printf("Added One (4)!\n");
		    }	// energy treshold question	    
		  }  //number of mechanisms cycle
		  /**********************End of evaluation************************/
		} //freon production question
	    } //energy interval question
	//}//inside the proximity gap question
    } //cerenkov photon question
      
    /**************************************End of Production Parameters Storing*********************/ 
    
    
    /*******************************Treat photons that hit the CsI (Ckovs and Feedbacks)************/ 
    
    if (gMC->TrackPid() == 50000050 || gMC->TrackPid() == 50000051) {
      //printf("Cerenkov\n");
      
      //if (gMC->TrackPid() == 50000051)
	//printf("Tracking a feedback\n");
      
      if (gMC->VolId("CSI ")==gMC->CurrentVolID(copy))
	{
	  //printf("Current volume (should be CSI) (2): %s\n",gMC->CurrentVolName());
	  //printf("VolId: %d, CurrentVolID: %d\n",gMC->VolId("CSI "),gMC->CurrentVolID(copy));
	  //printf("Got in CSI\n");
	  //printf("Tracking a %d\n",gMC->TrackPid());
	  if (gMC->Edep() > 0.){
		gMC->TrackPosition(position);
		gMC->TrackMomentum(momentum);
		pos[0]=position(0);
		pos[1]=position(1);
		pos[2]=position(2);
		mom[0]=momentum(0);
		mom[1]=momentum(1);
		mom[2]=momentum(2);
		mom[3]=momentum(3);
		Double_t tc = mom[0]*mom[0]+mom[1]*mom[1];
		Double_t rt = TMath::Sqrt(tc);
		theta   = Float_t(TMath::ATan2(rt,Double_t(mom[2])))*kRaddeg;
		phi     = Float_t(TMath::ATan2(Double_t(mom[1]),Double_t(mom[0])))*kRaddeg;
		
		gMC->CurrentVolOffID(2,copy);
		vol[0]=copy;
		idvol=vol[0]-1;
		

		gMC->Gmtod(pos,localPos,1);

		//Chamber(idvol).GlobaltoLocal(pos,localPos);
                                                                    
		gMC->Gmtod(mom,localMom,2);

		//Chamber(idvol).GlobaltoLocal(mom,localMom);
		
		gMC->CurrentVolOffID(2,copy);
		vol[0]=copy;
		idvol=vol[0]-1;

		//Int_t sector=((AliRICHChamber*) (*fChambers)[idvol])
			//->Sector(localPos[0], localPos[2]);
		//printf("Sector:%d\n",sector);

		/*if (gMC->TrackPid() == 50000051){
		  fFeedbacks++;
		  printf("Feedbacks:%d\n",fFeedbacks);
		}*/	
		
        //PH		((AliRICHChamber*) (*fChambers)[idvol])
		((AliRICHChamber*)fChambers->At(idvol))
		    ->SigGenInit(localPos[0], localPos[2], localPos[1]);
		if(idvol<kNCH) {	
		    ckovData[0] = gMC->TrackPid();        // particle type
		    ckovData[1] = pos[0];                 // X-position for hit
		    ckovData[2] = pos[1];                 // Y-position for hit
		    ckovData[3] = pos[2];                 // Z-position for hit
		    ckovData[4] = theta;                      // theta angle of incidence
		    ckovData[5] = phi;                      // phi angle of incidence 
		    ckovData[8] = (Float_t) fNsdigits;      // first sdigit
		    ckovData[9] = -1;                       // last pad hit
		    ckovData[13] = 4;                       // photon was detected
		    ckovData[14] = mom[0];
		    ckovData[15] = mom[1];
		    ckovData[16] = mom[2];
		    
		    destep = gMC->Edep();
		    gMC->SetMaxStep(kBig);
		    cherenkovLoss  += destep;
		    ckovData[7]=cherenkovLoss;
		    
		    ckovData[17] = Hits2SDigits(localPos[0],localPos[2],cherenkovLoss,idvol,kPhoton);//for photons in CsI 
		    		    
		    if (fNsdigits > (Int_t)ckovData[8]) {
			ckovData[8]= ckovData[8]+1;
			ckovData[9]= (Float_t) fNsdigits;
		    }

		    
		    //TClonesArray *Hits = RICH->Hits();
		    AliRICHhit *mipHit =  (AliRICHhit*) (fHits->UncheckedAt(0));
		    if (mipHit)
		      {
			mom[0] = current->Px();
			mom[1] = current->Py();
			mom[2] = current->Pz();
			Float_t mipPx = mipHit->MomX();
			Float_t mipPy = mipHit->MomY();
			Float_t mipPz = mipHit->MomZ();
			
			Float_t r = mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2];
			Float_t rt = TMath::Sqrt(r);
			Float_t mipR = mipPx*mipPx + mipPy*mipPy + mipPz*mipPz;	
			Float_t mipRt = TMath::Sqrt(mipR);
			if ((rt*mipRt) > 0)
			  {
			    coscerenkov = (mom[0]*mipPx + mom[1]*mipPy + mom[2]*mipPz)/(rt*mipRt);
			  }
			else
			  {
			    coscerenkov = 0;
			  }
			Float_t cherenkov = TMath::ACos(coscerenkov);
			ckovData[18]=cherenkov;
		      }
		    //if (sector != -1)
		    //{
		    AddHit(gAlice->GetCurrentTrackNumber(),vol,ckovData);
		    AddCerenkov(gAlice->GetCurrentTrackNumber(),vol,ckovData);
		    //printf("Added One (5)!\n");
		    //}
		}
	    }
	}
    }
    
    /***********************************************End of photon hits*********************************************/
    

    /**********************************************Charged particles treatment*************************************/

    else if (gMC->TrackCharge()){
//If MIP
	/*if (gMC->IsTrackEntering())
	  {                
	    hits[13]=20;//is track entering?
	  }*/
	if (gMC->VolId("FRE1")==gMC->CurrentVolID(copy) || gMC->VolId("FRE2")==gMC->CurrentVolID(copy))
	  {
	    gMC->TrackMomentum(momentum);
	    mom[0]=momentum(0);
	    mom[1]=momentum(1);
	    mom[2]=momentum(2);
	    mom[3]=momentum(3);
	    hits [19] = mom[0];
	    hits [20] = mom[1];
	    hits [21] = mom[2];
	    fFreonProd=1;
	  }

	if (gMC->VolId("GAP ")== gMC->CurrentVolID(copy)) {//is in GAP?
// Get current particle id (ipart), track position (pos)  and momentum (mom)
	    
	    gMC->CurrentVolOffID(3,copy);
	    vol[0]=copy;
	    idvol=vol[0]-1;

	    //Int_t sector=((AliRICHChamber*) (*fChambers)[idvol])
			//->Sector(localPos[0], localPos[2]);
	    //printf("Sector:%d\n",sector);
	    
	    gMC->TrackPosition(position);
	    gMC->TrackMomentum(momentum);
	    pos[0]=position(0);
	    pos[1]=position(1);
	    pos[2]=position(2);
	    mom[0]=momentum(0);
	    mom[1]=momentum(1);
	    mom[2]=momentum(2);
	    mom[3]=momentum(3);

	    gMC->Gmtod(pos,localPos,1);
	    
	    //Chamber(idvol).GlobaltoLocal(pos,localPos);
                                                                    
	    gMC->Gmtod(mom,localMom,2);

	    //Chamber(idvol).GlobaltoLocal(mom,localMom);
	    
	    ipart  = gMC->TrackPid();
	    //
	    // momentum loss and steplength in last step
	    destep = gMC->Edep();
	    step   = gMC->TrackStep();
  
	    //
	    // record hits when track enters ...
	    if( gMC->IsTrackEntering()) {
//		gMC->SetMaxStep(fMaxStepGas);
		Double_t tc = mom[0]*mom[0]+mom[1]*mom[1];
		Double_t rt = TMath::Sqrt(tc);
		theta   = Float_t(TMath::ATan2(rt,Double_t(mom[2])))*kRaddeg;
		phi     = Float_t(TMath::ATan2(Double_t(mom[1]),Double_t(mom[0])))*kRaddeg;
		

		Double_t localTc = localMom[0]*localMom[0]+localMom[2]*localMom[2];
		Double_t localRt = TMath::Sqrt(localTc);
		localTheta   = Float_t(TMath::ATan2(localRt,Double_t(localMom[1])))*kRaddeg;                       
		localPhi     = Float_t(TMath::ATan2(Double_t(localMom[2]),Double_t(localMom[0])))*kRaddeg;    
		
		hits[0] = Float_t(ipart);         // particle type
		hits[1] = localPos[0];                 // X-position for hit
		hits[2] = localPos[1];                 // Y-position for hit
		hits[3] = localPos[2];                 // Z-position for hit
		hits[4] = localTheta;                  // theta angle of incidence
		hits[5] = localPhi;                    // phi angle of incidence 
		hits[8] = (Float_t) fNsdigits;    // first sdigit
		hits[9] = -1;                     // last pad hit
		hits[13] = fFreonProd;           // did id hit the freon?
		hits[14] = mom[0];
		hits[15] = mom[1];
		hits[16] = mom[2];
		hits[18] = 0;               // dummy cerenkov angle

		tlength = 0;
		eloss   = 0;
		fFreonProd = 0;
	
		Chamber(idvol).LocaltoGlobal(localPos,hits+1);
	   
		
		//To make chamber coordinates x-y had to pass localPos[0], localPos[2]
		xhit    = localPos[0];
		yhit    = localPos[2];
		// Only if not trigger chamber
		if(idvol<kNCH) {
		    //
		    //  Initialize hit position (cursor) in the segmentation model 
          //PH		    ((AliRICHChamber*) (*fChambers)[idvol])
		    ((AliRICHChamber*)fChambers->At(idvol))
			->SigGenInit(localPos[0], localPos[2], localPos[1]);
		}
	    }
	    
	    // 
	    // Calculate the charge induced on a pad (disintegration) in case 
	    //
	    // Mip left chamber ...
	    if( gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared()){
		gMC->SetMaxStep(kBig);
		eloss   += destep;
		tlength += step;
		
				
		// Only if not trigger chamber
		if(idvol<kNCH) {
		  if (eloss > 0) 
		    {
		      if(gMC->TrackPid() == kNeutron)
			printf("\n\n\n\n\n Neutron Making Pad Hit!!! \n\n\n\n");
		      hits[17] = Hits2SDigits(xhit,yhit,eloss,idvol,kMip); //for MIP 
		    }
		}
		
		hits[6]=tlength;
		hits[7]=eloss;
		if (fNsdigits > (Int_t)hits[8]) {
		    hits[8]= hits[8]+1;
		    hits[9]= (Float_t) fNsdigits;
		}
		
		//if(sector !=-1)
		new(lhits[fNhits++]) AliRICHhit(fIshunt,gAlice->GetCurrentTrackNumber(),vol,hits);
		eloss = 0; 
		//
		// Check additional signal generation conditions 
		// defined by the segmentation
		// model (boundary crossing conditions) 
	    }else if(((AliRICHChamber*)fChambers->At(idvol))->SigGenCond(localPos[0], localPos[2], localPos[1])){
		((AliRICHChamber*)fChambers->At(idvol))->SigGenInit(localPos[0], localPos[2], localPos[1]);
		if (eloss > 0) 
		  {
		    if(gMC->TrackPid() == kNeutron)
		      printf("\n\n\n\n\n Neutron Making Pad Hit!!! \n\n\n\n");
		    hits[17] = Hits2SDigits(xhit,yhit,eloss,idvol,kMip);//for n
		  }
		xhit     = localPos[0];
		yhit     = localPos[2]; 
		eloss    = destep;
		tlength += step ;
		//
		// nothing special  happened, add up energy loss
	    } else {        
		eloss   += destep;
		tlength += step ;
	    }
	}//is in GAP?
      }//is MIP?
    /*************************************************End of MIP treatment**************************************/
}//void AliRICHv3::StepManager()
//__________________________________________________________________________________________________
Int_t AliRICHv3::Hits2SDigits(Float_t xhit,Float_t yhit,Float_t eloss, Int_t idvol, ResponseType res)
{//calls the charge disintegration method of the current chamber and adds all generated sdigits to the list of digits
   
   Int_t iChamber=kBad,iPadX=kBad,iPadY=kBad,iAdc=kBad,iTrack=kBad;
   Float_t list[4][500];
   Int_t iNdigits;
  ((AliRICHChamber*)fChambers->At(idvol))->DisIntegration(eloss, xhit, yhit, iNdigits, list, res);
    Int_t ic=0;
    
  for(Int_t i=0; i<iNdigits; i++) {
    if(Int_t(list[0][i]) > 0) {
	    ic++;
	    iAdc = Int_t(list[0][i]);
	    iPadX = Int_t(list[1][i]);
	    iPadY = Int_t(list[2][i]);
	    iChamber = Int_t(list[3][i]);
	    AddSDigit(iChamber,iPadX,iPadY,iAdc,iTrack);
	}
  }
    
  if(fLoader->TreeS()){
    fLoader->TreeS()->Fill();
    fLoader->WriteSDigits("OVERWRITE");
  }
   return iNdigits;
}//Int_t AliRICHv3::Hits2SDigits(Float_t xhit,Float_t yhit,Float_t eloss, Int_t idvol, ResponseType res)
//__________________________________________________________________________________________________
void AliRICHv3::DiagnosticsFE(Int_t evNumber1,Int_t evNumber2)
{
  
  Int_t NpadX = 162;                 // number of pads on X
  Int_t NpadY = 162;                 // number of pads on Y
  
  Int_t Pad[162][162];
  for (Int_t i=0;i<NpadX;i++) {
    for (Int_t j=0;j<NpadY;j++) {
      Pad[i][j]=0;
    }
  }
  
  //  Create some histograms

  TH1F *pionspectra1 = new TH1F("pionspectra1","Pion Spectra",200,-4,2);
  TH1F *pionspectra2 = new TH1F("pionspectra2","Pion Spectra",200,-4,2);
  TH1F *pionspectra3 = new TH1F("pionspectra3","Pion Spectra",200,-4,2);
  TH1F *protonspectra1 = new TH1F("protonspectra1","Proton Spectra",200,-4,2);
  TH1F *protonspectra2 = new TH1F("protonspectra2","Proton Spectra",200,-4,2);
  TH1F *protonspectra3 = new TH1F("protonspectra3","Proton Spectra",200,-4,2);
  TH1F *kaonspectra1 = new TH1F("kaonspectra1","Kaon Spectra",100,-4,2);
  TH1F *kaonspectra2 = new TH1F("kaonspectra2","Kaon Spectra",100,-4,2);
  TH1F *kaonspectra3 = new TH1F("kaonspectra3","Kaon Spectra",100,-4,2);
  TH1F *electronspectra1 = new TH1F("electronspectra1","Electron Spectra",100,-4,2);
  TH1F *electronspectra2 = new TH1F("electronspectra2","Electron Spectra",100,-4,2);
  TH1F *electronspectra3 = new TH1F("electronspectra3","Electron Spectra",100,-4,2);
  TH1F *muonspectra1 = new TH1F("muonspectra1","Muon Spectra",100,-4,2);
  TH1F *muonspectra2 = new TH1F("muonspectra2","Muon Spectra",100,-4,2);
  TH1F *muonspectra3 = new TH1F("muonspectra3","Muon Spectra",100,-4,2);
  TH1F *neutronspectra1 = new TH1F("neutronspectra1","Neutron Spectra",100,-4,2);
  TH1F *neutronspectra2 = new TH1F("neutronspectra2","Neutron Spectra",100,-4,2);
  TH1F *neutronspectra3 = new TH1F("neutronspectra2","Neutron Spectra",100,-4,2);
  TH1F *chargedspectra1 = new TH1F("chargedspectra1","Charged particles above 1 GeV Spectra",100,-1,3);
  TH1F *chargedspectra2 = new TH1F("chargedspectra2","Charged particles above 1 GeV Spectra",100,-1,3);
  TH1F *chargedspectra3 = new TH1F("chargedspectra2","Charged particles above 1 GeV Spectra",100,-1,3);
  TH1F *pionptspectrafinal = new TH1F("pionptspectrafinal","Primary Pions Transverse Momenta at HMPID",20,0,5);
  TH1F *pionptspectravertex = new TH1F("pionptspectravertex","Primary Pions Transverse Momenta at vertex",20,0,5);
  TH1F *kaonptspectrafinal = new TH1F("kaonptspectrafinal","Primary Kaons Transverse Momenta at HMPID",20,0,5);
  TH1F *kaonptspectravertex = new TH1F("kaonptspectravertex","Primary Kaons Transverse Momenta at vertex",20,0,5);
  //TH1F *hitsPhi = new TH1F("hitsPhi","Distribution of phi angle of incidence",100,-180,180);
  TH1F *hitsTheta = new TH1F("hitsTheta","Distribution of Theta angle of incidence, all tracks",100,0,50);
  TH1F *hitsTheta500MeV = new TH1F("hitsTheta500MeV","Distribution of Theta angle of incidence, 0.5-1 GeV primary tracks",100,0,50);
  TH1F *hitsTheta1GeV = new TH1F("hitsTheta1GeV","Distribution of Theta angle of incidence, 1-2 GeV primary tracks",100,0,50);
  TH1F *hitsTheta2GeV = new TH1F("hitsTheta2GeV","Distribution of Theta angle of incidence, 2-3 GeV primary tracks",100,0,50);
  TH1F *hitsTheta3GeV = new TH1F("hitsTheta3GeV","Distribution of Theta angle of incidence, >3 GeV primary tracks",100,0,50);
  TH2F *production = new TH2F("production","Mother production vertices",100,-300,300,100,0,600);
   
   
   

//   Start loop over events 

  Int_t pion=0, kaon=0, proton=0, electron=0, positron=0, neutron=0, highneutrons=0, muon=0;
  Int_t chargedpions=0,primarypions=0,highprimarypions=0,chargedkaons=0,primarykaons=0,highprimarykaons=0;
  Int_t photons=0, primaryphotons=0, highprimaryphotons=0;
  TRandom* random=0;

   for (int nev=0; nev<= evNumber2; nev++) {
       Int_t nparticles = gAlice->GetEvent(nev);
       

       if (nev < evNumber1) continue;
       if (nparticles <= 0) return;
       
// Get pointers to RICH detector and Hits containers
       
       AliRICH *pRICH = (AliRICH *) gAlice->GetDetector("RICH");
     
       TTree *treeH = TreeH();
       Int_t ntracks =(Int_t) treeH->GetEntries();
            
// Start loop on tracks in the hits containers
       
       for (Int_t track=0; track<ntracks;track++) {
	   printf ("Processing Track: %d\n",track);
	   gAlice->ResetHits();
	   treeH->GetEvent(track);
	   	   	   
	   for(AliRICHhit* mHit=(AliRICHhit*)pRICH->FirstHit(-1); 
	       mHit;
	       mHit=(AliRICHhit*)pRICH->NextHit()) 
	     {
	       //Int_t nch  = mHit->fChamber;              // chamber number
	       //Float_t x  = mHit->X();                    // x-pos of hit
	       //Float_t y  = mHit->Z();                    // y-pos
	       //Float_t z  = mHit->Y();
	       //Float_t phi = mHit->Phi();                 //Phi angle of incidence
	       Float_t theta = mHit->Theta();             //Theta angle of incidence
	       Float_t px = mHit->MomX();
	       Float_t py = mHit->MomY();
	       Int_t index = mHit->Track();
	       Int_t particle = (Int_t)(mHit->Particle());    
	       Float_t R;
	       Float_t PTfinal;
	       Float_t PTvertex;

	      TParticle *current = gAlice->Particle(index);
	      
	      //Float_t energy=current->Energy(); 

	      R=TMath::Sqrt(current->Vx()*current->Vx() + current->Vy()*current->Vy());
	      PTfinal=TMath::Sqrt(px*px + py*py);
	      PTvertex=TMath::Sqrt(current->Px()*current->Px() + current->Py()*current->Py());
	      
	      

	      if (TMath::Abs(particle) < 10000000)
		{
		  hitsTheta->Fill(theta,(float) 1);
		  if (R<5)
		    {
		      if (PTvertex>.5 && PTvertex<=1)
			{
			  hitsTheta500MeV->Fill(theta,(float) 1);
			}
		      if (PTvertex>1 && PTvertex<=2)
			{
			  hitsTheta1GeV->Fill(theta,(float) 1);
			}
		      if (PTvertex>2 && PTvertex<=3)
			{
			  hitsTheta2GeV->Fill(theta,(float) 1);
			}
		      if (PTvertex>3)
			{
			  hitsTheta3GeV->Fill(theta,(float) 1);
			}
		    }
		  
		}

	      //if (nch == 3)
		//{
	      
	      if (TMath::Abs(particle) < 50000051)
		{
		  //if (TMath::Abs(particle) == 50000050 || TMath::Abs(particle) == 2112)
		  if (TMath::Abs(particle) == 2112 || TMath::Abs(particle) == 50000050)
		    {
		      //gMC->Rndm(&random, 1);
		      if (random->Rndm() < .1)
			production->Fill(current->Vz(),R,(float) 1);
		      if (TMath::Abs(particle) == 50000050)
			//if (TMath::Abs(particle) > 50000000)
			{
			  photons +=1;
			  if (R<5)
			    {
			      primaryphotons +=1;
			      if (current->Energy()>0.001)
				highprimaryphotons +=1;
			    }
			}	
		      if (TMath::Abs(particle) == 2112)
			{
			  neutron +=1;
			  if (current->Energy()>0.0001)
			    highneutrons +=1;
			}
		    }
		  if (TMath::Abs(particle) < 50000000)
		    {
		      production->Fill(current->Vz(),R,(float) 1);
		    }
		  //mip->Fill(x,y,(float) 1);
		}
	      
	      if (TMath::Abs(particle)==211 || TMath::Abs(particle)==111)
		{
		  if (R<5)
		    {
		      pionptspectravertex->Fill(PTvertex,(float) 1);
		      pionptspectrafinal->Fill(PTfinal,(float) 1);
		    }
		}
	      
	      if (TMath::Abs(particle)==321 || TMath::Abs(particle)==130 || TMath::Abs(particle)==310 
		  || TMath::Abs(particle)==311)
		{
		  if (R<5)
		    {
		      kaonptspectravertex->Fill(PTvertex,(float) 1);
		      kaonptspectrafinal->Fill(PTfinal,(float) 1);
		    }
		}
	      
	      
	      if (TMath::Abs(particle)==211 || TMath::Abs(particle)==111)
		{
		  pionspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (current->Vx()>5 && current->Vy()>5 && current->Vz()>5)
		    pionspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>250 && R<450)
		    {
		      pionspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		    }
		  pion +=1;
		  if (TMath::Abs(particle)==211)
		    {
		      chargedpions +=1;
		      if (R<5)
			{
			  primarypions +=1;
			  if (current->Energy()>1)
			    highprimarypions +=1;
			}
		    }	
		}
	      if (TMath::Abs(particle)==2212)
		{
		  protonspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  //ptspectra->Fill(Pt,(float) 1);
		  if (current->Vx()>5 && current->Vy()>5 && current->Vz()>5)
		    protonspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>250 && R<450)
		    protonspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  proton +=1;
		}
	      if (TMath::Abs(particle)==321 || TMath::Abs(particle)==130 || TMath::Abs(particle)==310 
		  || TMath::Abs(particle)==311)
		{
		  kaonspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  //ptspectra->Fill(Pt,(float) 1);
		  if (current->Vx()>5 && current->Vy()>5 && current->Vz()>5)
		    kaonspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>250 && R<450)
		    kaonspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  kaon +=1;
		  if (TMath::Abs(particle)==321)
		    {
		      chargedkaons +=1;
		      if (R<5)
			{
			  primarykaons +=1;
			  if (current->Energy()>1)
			    highprimarykaons +=1;
			}
		    }
		}
	      if (TMath::Abs(particle)==11)
		{
		  electronspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  //ptspectra->Fill(Pt,(float) 1);
		  if (current->Vx()>5 && current->Vy()>5 && current->Vz()>5)
		    electronspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>250 && R<450)
		    electronspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (particle == 11)
		    electron +=1;
		  if (particle == -11)
		    positron +=1;
		}
	      if (TMath::Abs(particle)==13)
		{
		  muonspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  //ptspectra->Fill(Pt,(float) 1);
		  if (current->Vx()>5 && current->Vy()>5 && current->Vz()>5)
		    muonspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>250 && R<450)
		    muonspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  muon +=1;
		}
	      if (TMath::Abs(particle)==2112)
		{
		  neutronspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  //ptspectra->Fill(Pt,(float) 1);
		  if (current->Vx()>5 && current->Vy()>5 && current->Vz()>5)
		    neutronspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		  if (R>250 && R<450)
		    {
		      neutronspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		    }
		  neutron +=1;
		}
	      if(TMath::Abs(particle)==211 || TMath::Abs(particle)==2212 || TMath::Abs(particle)==321)
		{
		  if (current->Energy()-current->GetCalcMass()>1)
		    {
		      chargedspectra1->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		      if (current->Vx()>5 && current->Vy()>5 && current->Vz()>5)
			chargedspectra2->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		      if (R>250 && R<450)
			chargedspectra3->Fill(TMath::Log10(current->Energy() - current->GetCalcMass()),(float) 1);
		    }
		}
	      // Fill the histograms
	      //Nh1+=nhits;
	      //h->Fill(x,y,(float) 1);
	      //}
	      //}
	   }          
	   
       }
       
   }
   //   }

   TStyle *mystyle=new TStyle("Plain","mystyle");
   mystyle->SetPalette(1,0);
   mystyle->cd();
   
   //Create canvases, set the view range, show histograms

    TCanvas *c2 = new TCanvas("c2","Angles of incidence",150,150,100,150);
    c2->Divide(2,2);
    //c2->SetFillColor(42);
    
    c2->cd(1);
    hitsTheta500MeV->SetFillColor(5);
    hitsTheta500MeV->Draw();
    c2->cd(2);
    hitsTheta1GeV->SetFillColor(5);
    hitsTheta1GeV->Draw();
    c2->cd(3);
    hitsTheta2GeV->SetFillColor(5);
    hitsTheta2GeV->Draw();
    c2->cd(4);
    hitsTheta3GeV->SetFillColor(5);
    hitsTheta3GeV->Draw();
    
            
   
    TCanvas *c15 = new TCanvas("c15","Mothers Production Vertices",50,50,600,600);
    c15->cd();
    production->SetFillColor(42);
    production->SetXTitle("z (m)");
    production->SetYTitle("R (m)");
    production->Draw();

    TCanvas *c10 = new TCanvas("c10","Pt Spectra",50,50,600,700);
    c10->Divide(2,2);
    c10->cd(1);
    pionptspectravertex->SetFillColor(5);
    pionptspectravertex->SetXTitle("Pt (GeV)");
    pionptspectravertex->Draw();
    c10->cd(2);
    pionptspectrafinal->SetFillColor(5);
    pionptspectrafinal->SetXTitle("Pt (GeV)");
    pionptspectrafinal->Draw();
    c10->cd(3);
    kaonptspectravertex->SetFillColor(5);
    kaonptspectravertex->SetXTitle("Pt (GeV)");
    kaonptspectravertex->Draw();
    c10->cd(4);
    kaonptspectrafinal->SetFillColor(5);
    kaonptspectrafinal->SetXTitle("Pt (GeV)");
    kaonptspectrafinal->Draw();
   
  
   TCanvas *c16 = new TCanvas("c16","Particles Spectra II",150,150,600,350);
   c16->Divide(2,1);
   
   c16->cd(1);
   //TCanvas *c13 = new TCanvas("c13","Electron Spectra",400,10,600,700);
   electronspectra1->SetFillColor(5);
   electronspectra1->SetXTitle("log(GeV)");
   electronspectra2->SetFillColor(46);
   electronspectra2->SetXTitle("log(GeV)");
   electronspectra3->SetFillColor(10);
   electronspectra3->SetXTitle("log(GeV)");
   //c13->SetLogx();
   electronspectra1->Draw();
   electronspectra2->Draw("same");
   electronspectra3->Draw("same");
   
   c16->cd(2);
   //TCanvas *c14 = new TCanvas("c14","Muon Spectra",400,10,600,700);
   muonspectra1->SetFillColor(5);
   muonspectra1->SetXTitle("log(GeV)");
   muonspectra2->SetFillColor(46);
   muonspectra2->SetXTitle("log(GeV)");
   muonspectra3->SetFillColor(10);
   muonspectra3->SetXTitle("log(GeV)");
   //c14->SetLogx();
   muonspectra1->Draw();
   muonspectra2->Draw("same");
   muonspectra3->Draw("same");
   
   //c16->cd(3);
   //TCanvas *c16 = new TCanvas("c16","Neutron Spectra",400,10,600,700);
   //neutronspectra1->SetFillColor(42);
   //neutronspectra1->SetXTitle("log(GeV)");
   //neutronspectra2->SetFillColor(46);
   //neutronspectra2->SetXTitle("log(GeV)");
   //neutronspectra3->SetFillColor(10);
   //neutronspectra3->SetXTitle("log(GeV)");
   //c16->SetLogx();
   //neutronspectra1->Draw();
   //neutronspectra2->Draw("same");
   //neutronspectra3->Draw("same");

   TCanvas *c9 = new TCanvas("c9","Particles Spectra",150,150,600,700);
   //TCanvas *c9 = new TCanvas("c9","Pion Spectra",400,10,600,700);
   c9->Divide(2,2);
   
   c9->cd(1);
   pionspectra1->SetFillColor(5);
   pionspectra1->SetXTitle("log(GeV)");
   pionspectra2->SetFillColor(46);
   pionspectra2->SetXTitle("log(GeV)");
   pionspectra3->SetFillColor(10);
   pionspectra3->SetXTitle("log(GeV)");
   //c9->SetLogx();
   pionspectra1->Draw();
   pionspectra2->Draw("same");
   pionspectra3->Draw("same");
   
   c9->cd(2);
   //TCanvas *c10 = new TCanvas("c10","Proton Spectra",400,10,600,700);
   protonspectra1->SetFillColor(5);
   protonspectra1->SetXTitle("log(GeV)");
   protonspectra2->SetFillColor(46);
   protonspectra2->SetXTitle("log(GeV)");
   protonspectra3->SetFillColor(10);
   protonspectra3->SetXTitle("log(GeV)");
   //c10->SetLogx();
   protonspectra1->Draw();
   protonspectra2->Draw("same");
   protonspectra3->Draw("same");
   
   c9->cd(3);
   //TCanvas *c11 = new TCanvas("c11","Kaon Spectra",400,10,600,700); 
   kaonspectra1->SetFillColor(5);
   kaonspectra1->SetXTitle("log(GeV)");
   kaonspectra2->SetFillColor(46);
   kaonspectra2->SetXTitle("log(GeV)");
   kaonspectra3->SetFillColor(10);
   kaonspectra3->SetXTitle("log(GeV)");
   //c11->SetLogx();
   kaonspectra1->Draw();
   kaonspectra2->Draw("same");
   kaonspectra3->Draw("same");
   
   c9->cd(4);
   //TCanvas *c12 = new TCanvas("c12","Charged Particles Spectra",400,10,600,700);
   chargedspectra1->SetFillColor(5);
   chargedspectra1->SetXTitle("log(GeV)");
   chargedspectra2->SetFillColor(46);
   chargedspectra2->SetXTitle("log(GeV)");
   chargedspectra3->SetFillColor(10);
   chargedspectra3->SetXTitle("log(GeV)");
   //c12->SetLogx();
   chargedspectra1->Draw();
   chargedspectra2->Draw("same");
   chargedspectra3->Draw("same");
   


   printf("*****************************************\n");
   printf("* Particle                   *  Counts  *\n");
   printf("*****************************************\n");

   printf("* Pions:                     *   %4d   *\n",pion);
   printf("* Charged Pions:             *   %4d   *\n",chargedpions);
   printf("* Primary Pions:             *   %4d   *\n",primarypions);
   printf("* Primary Pions (p>1GeV/c):  *   %4d   *\n",highprimarypions);
   printf("* Kaons:                     *   %4d   *\n",kaon);
   printf("* Charged Kaons:             *   %4d   *\n",chargedkaons);
   printf("* Primary Kaons:             *   %4d   *\n",primarykaons);
   printf("* Primary Kaons (p>1GeV/c):  *   %4d   *\n",highprimarykaons);
   printf("* Muons:                     *   %4d   *\n",muon);
   printf("* Electrons:                 *   %4d   *\n",electron);
   printf("* Positrons:                 *   %4d   *\n",positron);
   printf("* Protons:                   *   %4d   *\n",proton);
   printf("* All Charged:               *   %4d   *\n",(chargedpions+chargedkaons+muon+electron+positron+proton));
   printf("*****************************************\n");
   //printf("* Photons:                   *   %3.1f   *\n",photons); 
   //printf("* Primary Photons:           *   %3.1f   *\n",primaryphotons);
   //printf("* Primary Photons (p>1MeV/c):*   %3.1f   *\n",highprimaryphotons);
   //printf("*****************************************\n");
   //printf("* Neutrons:                  *   %3.1f   *\n",neutron);
   //printf("* Neutrons (p>100keV/c):     *   %3.1f   *\n",highneutrons);
   //printf("*****************************************\n");

   if (gAlice->TreeD())
     {
       gAlice->TreeD()->GetEvent(0);
   
       Float_t occ[7]; 
       Float_t sum=0;
       Float_t mean=0; 
       printf("\n*****************************************\n");
       printf("* Chamber   * Digits      * Occupancy   *\n");
       printf("*****************************************\n");
       
       for (Int_t ich=0;ich<7;ich++)
	 {
	   TClonesArray *Digits = DigitsAddress(ich);    //  Raw clusters branch
	   Int_t ndigits = Digits->GetEntriesFast();
	   occ[ich] = Float_t(ndigits)/(160*144);
	   sum += Float_t(ndigits)/(160*144);
	   printf("*   %d      *    %d      *   %3.1f%%     *\n",ich,ndigits,occ[ich]*100);
	 }
       mean = sum/7;
       printf("*****************************************\n");
       printf("* Mean occupancy          *   %3.1f%%     *\n",mean*100);
       printf("*****************************************\n");
     }
 
  printf("\nEnd of analysis\n");
   
}//void AliRICHv3::DiagnosticsFE(Int_t evNumber1,Int_t evNumber2)
//__________________________________________________________________________________________________
void AliRICHv3::DiagnosticsSE(Int_t diaglevel,Int_t evNumber1,Int_t evNumber2)
{

AliRICH *pRICH  = (AliRICH*)gAlice->GetDetector("RICH");
   AliRICHSegmentationV0*  segmentation;
   AliRICHChamber*       chamber;
   
   chamber = &(pRICH->Chamber(0));
   segmentation=(AliRICHSegmentationV0*) chamber->GetSegmentationModel();

   Int_t NpadX = segmentation->Npx();                 // number of pads on X
   Int_t NpadY = segmentation->Npy();                 // number of pads on Y
    
   Int_t xmin= -NpadX/2;  
   Int_t xmax=  NpadX/2;
   Int_t ymin= -NpadY/2;
   Int_t ymax=  NpadY/2;

   Float_t PTfinal = 0;
   Int_t pionCount = 0;
   Int_t kaonCount = 0;
   Int_t protonCount = 0;
   
   TH2F *feedback = 0;
   TH2F *mip = 0;
   TH2F *cerenkov = 0;
   TH2F *h = 0;
   TH1F *hitsX = 0;
   TH1F *hitsY = 0;

   TH2F *hc0 = new TH2F("hc0","Zoom on center of central chamber",150,-25,25,150,-45,5);

   if (diaglevel == 1)
     {
       printf("Single Ring Hits\n");
       feedback = new TH2F("feedback","Feedback hit distribution",150,-20,20,150,-35,5);
       mip = new TH2F("mip","Mip hit distribution",150,-20,20,150,-35,5);
       cerenkov = new TH2F("cerenkov","Cerenkov hit distribution",150,-20,20,150,-35,5);
       h = new TH2F("h","Detector hit distribution",150,-20,20,150,-35,5);
       hitsX = new TH1F("hitsX","Distribution of hits along x-axis",150,-50,50);
       hitsY = new TH1F("hitsY","Distribution of hits along z-axis",150,-50,50);
     }       
   else
     {
       printf("Full Event Hits\n");
       
       feedback = new TH2F("feedback","Feedback hit distribution",150,-300,300,150,-300,300);
       mip = new TH2F("mip","Mip hit distribution",150,-300,300,150,-300,300);
       cerenkov = new TH2F("cerenkov","Cerenkov hit distribution",150,-300,300,150,-300,300);
       h = new TH2F("h","Detector hit distribution",150,-300,300,150,-300,300); 
       hitsX = new TH1F("digitsX","Distribution of hits along x-axis",200,-300,300);
       hitsY = new TH1F("digitsY","Distribution of hits along z-axis",200,-300,300);
     }
   


   TH2F *hc1 = new TH2F("hc1","Chamber 1 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *hc2 = new TH2F("hc2","Chamber 2 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *hc3 = new TH2F("hc3","Chamber 3 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *hc4 = new TH2F("hc4","Chamber 4 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *hc5 = new TH2F("hc5","Chamber 5 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *hc6 = new TH2F("hc6","Chamber 6 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
   TH2F *hc7 = new TH2F("hc7","Chamber 7 signal distribution",NpadX,xmin,xmax,NpadY,ymin,ymax);
      
   TH1F *Clcharge = new TH1F("Clcharge","Cluster Charge Distribution",500,0.,500.);
   TH1F *ckovangle = new TH1F("ckovangle","Cerenkov angle per photon",100,.35,.8);
   TH1F *hckphi = new TH1F("hckphi","Cerenkov phi angle per photon",620,-3.1,3.1);
   TH1F *mother = new TH1F("mother","Cerenkovs per Mip",75,0.,75.);
   TH1F *radius = new TH1F("radius","Mean distance to Mip",100,0.,20.);
   TH1F *phspectra1 = new TH1F("phspectra1","Detected Photon Spectra",200,5.,10.);
   TH1F *phspectra2 = new TH1F("phspectra2","Produced Photon Spectra",200,5.,10.);
   TH1F *totalphotonstrack = new TH1F("totalphotonstrack","Produced Photons per Mip",100,200,700.);
   TH1F *totalphotonsevent = new TH1F("totalphotonsevent","Produced Photons per Mip",100,200,700.);
   //TH1F *feedbacks = new TH1F("feedbacks","Produced Feedbacks per Mip",50,0.5,50.);
   TH1F *padnumber = new TH1F("padnumber","Number of pads per cluster",50,-0.5,50.);
   TH1F *padsev = new TH1F("padsev","Number of pads hit per MIP",50,0.5,100.);
   TH1F *clusev = new TH1F("clusev","Number of clusters per MIP",50,0.5,50.);
   TH1F *photev = new TH1F("photev","Number of detected photons per MIP",50,0.5,50.);
   TH1F *feedev = new TH1F("feedev","Number of feedbacks per MIP",50,0.5,50.);
   TH1F *padsmip = new TH1F("padsmip","Number of pads per event inside MIP region",50,0.5,50.);
   TH1F *padscl = new TH1F("padscl","Number of pads per event from cluster count",50,0.5,100.);
   TH1F *pionspectra = new TH1F("pionspectra","Pion Spectra",200,.5,10.);
   TH1F *protonspectra = new TH1F("protonspectra","Proton Spectra",200,.5,10.);
   TH1F *kaonspectra = new TH1F("kaonspectra","Kaon Spectra",100,.5,10.);
   TH1F *chargedspectra = new TH1F("chargedspectra","Charged particles above 1 GeV Spectra",100,.5,10.);
   TH1F *hitsPhi = new TH1F("hitsPhi","Distribution of phi angle of incidence",50,0,360);
   TH1F *hitsTheta = new TH1F("hitsTheta","Distribution of theta angle of incidence",50,0,15);
   TH1F *Omega1D = new TH1F("omega","Reconstructed Cerenkov angle per track",50,.5,1);
   TH1F *Theta = new TH1F("theta","Reconstructed theta incidence angle per track",100,0,15);
   TH1F *Phi = new TH1F("phi","Reconstructed phi incidence per track",100,0,360);
   TH1F *Omega3D = new TH1F("omega","Reconstructed Cerenkov angle per track",100,.35,.8);
   TH1F *PhotonCer = new TH1F("photoncer","Reconstructed Cerenkov angle per photon",100,.35,.8);
   TH2F *PadsUsed = new TH2F("padsused","Pads Used for Reconstruction",100,-30,30,100,-30,30);
   TH1F *MeanRadius = new TH1F("radius","Mean Radius for reconstructed track",100,0.,20.);
   TH2F *identification = new TH2F("identification","Particle Identification",100,1,5,100,0,.8);
   TH1F *OriginalOmega = new TH1F("Original Omega","Cerenkov angle per track",100,.35,.8);
   TH1F *OriginalPhi = new TH1F("Original Phi","Distribution of phi angle of incidence per track",100,0,360);
   TH1F *OriginalTheta = new TH1F("Original Theta","Distribution of theta angle per track",100,0,15);
   TH1F *OmegaError = new TH1F("Omega Error","Difference between original an reconstructed cerenkov angle",100,0,.2);
   TH1F *PhiError = new TH1F("Phi Error","Difference between original an reconstructed phi angle",100,0,360);
   TH1F *ThetaError = new TH1F("Theta Error","Difference between original an reconstructed phi angle",100,0,15);


//   Start loop over events 

   Int_t Nh=0;
   Int_t pads=0;
   Int_t Nh1=0;
   Int_t mothers[80000];
   Int_t mothers2[80000];
   Float_t mom[3];
   Int_t nraw=0;
   Int_t phot=0;
   Int_t feed=0;
   Int_t padmip=0;
   Float_t x=0,y=0;

   Float_t chiSquareOmega = 0;
   Float_t chiSquareTheta = 0;
   Float_t chiSquarePhi = 0;

   Float_t recEffEvent = 0;
   Float_t recEffTotal = 0;

   Float_t trackglob[3];
   Float_t trackloc[3];

   
   for (Int_t i=0;i<100;i++) mothers[i]=0;

   for (int nev=0; nev<= evNumber2; nev++) {
       Int_t nparticles = gAlice->GetEvent(nev);
       

       //cout<<"nev  "<<nev<<endl;
       printf ("\n**********************************\nProcessing Event: %d\n",nev);
       //cout<<"nparticles  "<<nparticles<<endl;
       printf ("Particles       : %d\n\n",nparticles);
       if (nev < evNumber1) continue;
       if (nparticles <= 0) return;
       
// Get pointers to RICH detector and Hits containers
       

       TTree *TH = TreeH(); 
       Stat_t ntracks = TH->GetEntries();

       // Start loop on tracks in the hits containers
       //Int_t Nc=0;
       for (Int_t track=0; track<ntracks;track++) {
	   
	 printf ("\nProcessing Track: %d\n",track);
	 gAlice->ResetHits();
	 TH->GetEvent(track);
	 Int_t nhits = pRICH->Hits()->GetEntriesFast();
	 if (nhits) Nh+=nhits;
	 printf("Hits            : %d\n",nhits);
	 for(AliRICHhit* mHit=(AliRICHhit*)pRICH->FirstHit(-1); 
	     mHit;
	     mHit=(AliRICHhit*)pRICH->NextHit()) 
	   {
	     Int_t nch  = mHit->Chamber();              // chamber number
	     trackglob[0] = mHit->X();                 // x-pos of hit
	     trackglob[1] = mHit->Y();
	     trackglob[2] = mHit->Z();                 // y-pos of hit
	     //x  = mHit->X();                           // x-pos of hit
	     //y  = mHit->Z();                           // y-pos
	     Float_t phi = mHit->Phi();                 //Phi angle of incidence
	     Float_t theta = mHit->Theta();             //Theta angle of incidence
	     Int_t index = mHit->Track();
	     Int_t particle = (Int_t)(mHit->Particle());        
	     //Int_t freon = (Int_t)(mHit->fLoss);    
	     Float_t px = mHit->MomX();
	     Float_t py = mHit->MomY();
	     
	     if (TMath::Abs(particle) < 10000000)
	       {
		 PTfinal=TMath::Sqrt(px*px + py*py);
	       }
	
	     chamber = &(pRICH->Chamber(nch-1));
	     
	     
	     chamber->GlobaltoLocal(trackglob,trackloc);
	     
	     chamber->LocaltoGlobal(trackloc,trackglob);
	     
       
	     x=trackloc[0];
	     y=trackloc[2];
	     
	     hitsX->Fill(x,(float) 1);
	     hitsY->Fill(y,(float) 1);
	       
	      
	      TParticle *current = (TParticle*)gAlice->Particle(index);

	      hitsTheta->Fill(theta,(float) 1);
	     
	      if (current->GetPdgCode() < 10000000)
		{
		  mip->Fill(x,y,(float) 1);
		  hitsPhi->Fill(TMath::Abs(phi),(float) 1);
		}
	      
	      if (TMath::Abs(particle)==211 || TMath::Abs(particle)==111)
		{
		  pionspectra->Fill(current->Energy() - current->GetCalcMass(),(float) 1);
		}
	      if (TMath::Abs(particle)==2212)
		{
		  protonspectra->Fill(current->Energy() - current->GetCalcMass(),(float) 1);
		}
	      if (TMath::Abs(particle)==321 || TMath::Abs(particle)==130 || TMath::Abs(particle)==310 
		  || TMath::Abs(particle)==311)
		{
		  kaonspectra->Fill(current->Energy() - current->GetCalcMass(),(float) 1);
		}
	      if(TMath::Abs(particle)==211 || TMath::Abs(particle)==2212 || TMath::Abs(particle)==321)
		{
		  if (current->Energy() - current->GetCalcMass()>1)
		    chargedspectra->Fill(current->Energy() - current->GetCalcMass(),(float) 1);
		}
	      //printf("Hits:%d\n",hit);
	      //printf ("Chamber number:%d x:%f y:%f\n",nch,x,y);
              // Fill the histograms
	      Nh1+=nhits;
	      h->Fill(x,y,(float) 1);
		  //}
              //}
	   }
	   
	   Int_t ncerenkovs = pRICH->Cerenkovs()->GetEntriesFast();
	   //if (current->GetPdgCode() < 50000051 && current->GetPdgCode() > 50000040)
	   //totalphotonsevent->Fill(ncerenkovs,(float) 1);

 	   if (ncerenkovs) {
	     printf("Cerenkovs       : %d\n",ncerenkovs);
	     totalphotonsevent->Fill(ncerenkovs,(float) 1);
	     for (Int_t hit=0;hit<ncerenkovs;hit++) {
	       AliRICHCerenkov* cHit = (AliRICHCerenkov*) pRICH->Cerenkovs()->UncheckedAt(hit);
	       Int_t nchamber = cHit->fChamber;     // chamber number
	       Int_t index =    cHit->Track();
	       //Int_t pindex =   (Int_t)(cHit->fIndex);
	       trackglob[0] = cHit->X();                 // x-pos of hit
	       trackglob[1] = cHit->Y();
	       trackglob[2] = cHit->Z();                 // y-pos of hit
	       //Float_t cx  =      cHit->X();                // x-position
	       //Float_t cy  =      cHit->Z();                // y-position
	       Int_t cmother =  cHit->fCMother;      // Index of mother particle
	       Int_t closs =    (Int_t)(cHit->fLoss);           // How did the particle get lost? 
	       Float_t cherenkov = cHit->fCerenkovAngle;   //production cerenkov angle
	       
	       chamber = &(pRICH->Chamber(nchamber-1));
	     
	       //printf("Nch:%d\n",nch);
	       
	       chamber->GlobaltoLocal(trackglob,trackloc);
	     
	       chamber->LocaltoGlobal(trackloc,trackglob);
	     
       
	       Float_t cx=trackloc[0];
	       Float_t cy=trackloc[2];
	       
	       //printf ("Cerenkov hit number %d/%d, X:%f, Y:%f\n",hit,ncerenkovs,cx,cy); 


	       //printf("Particle:%9d\n",index);
		 		 
	       TParticle *current = (TParticle*)gAlice->Particle(index);
	       Float_t energyckov = current->Energy();
	       
	       if (current->GetPdgCode() == 50000051)
		 {
		   if (closs==4)
		     {
		       feedback->Fill(cx,cy,(float) 1);
		       feed++;
		     }
		 }
	       if (current->GetPdgCode() == 50000050)
		 {
		   
		   if (closs !=4)
		     {
		       phspectra2->Fill(energyckov*1e9,(float) 1);
		     }
		       
		   if (closs==4)
		     {
		       cerenkov->Fill(cx,cy,(float) 1); 
		       
		       //printf ("Cerenkov hit number %d/%d, X:%d, Y:%d\n",hit,ncerenkovs,cx,cy); 
		       
		       //TParticle *MIP = (TParticle*)gAlice->Particle(cmother);
		       AliRICHhit* mipHit = (AliRICHhit*) pRICH->Hits()->UncheckedAt(0);
		       mom[0] = current->Px();
		       mom[1] = current->Py();
		       mom[2] = current->Pz();
		       //mom[0] = cHit->fMomX;
		       // mom[1] = cHit->fMomZ;
		       //mom[2] = cHit->fMomY;
		       //Float_t energymip = MIP->Energy();
		       //Float_t Mip_px = mipHit->fMomFreoX;
		       //Float_t Mip_py = mipHit->fMomFreoY;
		       //Float_t Mip_pz = mipHit->fMomFreoZ;
		       //Float_t Mip_px = MIP->Px();
		       //Float_t Mip_py = MIP->Py();
		       //Float_t Mip_pz = MIP->Pz();
		       
		       
		       
		       //Float_t r = mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2];
		       //Float_t rt = TMath::Sqrt(r);
		       //Float_t Mip_r = Mip_px*Mip_px + Mip_py*Mip_py + Mip_pz*Mip_pz;	
		       //Float_t Mip_rt = TMath::Sqrt(Mip_r);
		       //Float_t coscerenkov = (mom[0]*Mip_px + mom[1]*Mip_py + mom[2]*Mip_pz)/(rt*Mip_rt+0.0000001);
		       //Float_t cherenkov = TMath::ACos(coscerenkov);
		       ckovangle->Fill(cherenkov,(float) 1);                           //Cerenkov angle calculus
		       //printf("Cherenkov: %f\n",cherenkov);
		       Float_t ckphi=TMath::ATan2(mom[0], mom[2]);
		       hckphi->Fill(ckphi,(float) 1);
		       
		       
		       //Float_t mix = MIP->Vx();
		       //Float_t miy = MIP->Vy();
		       Float_t mx = mipHit->X();
		       Float_t my = mipHit->Z();
		       //printf("FX %e, FY %e, VX %e, VY %e\n",cx,cy,mx,my);
		       Float_t dx = trackglob[0] - mx;
		       Float_t dy = trackglob[2] - my;
		       //printf("Dx:%f, Dy:%f\n",dx,dy);
		       Float_t final_radius = TMath::Sqrt(dx*dx+dy*dy);
		       //printf("Final radius:%f\n",final_radius);
		       radius->Fill(final_radius,(float) 1);
		       
		       phspectra1->Fill(energyckov*1e9,(float) 1);
		       phot++;
		     }
		   for (Int_t nmothers=0;nmothers<=ntracks;nmothers++){
		     if (cmother == nmothers){
		       if (closs == 4)
			 mothers2[cmother]++;
		       mothers[cmother]++;
		     }
		   } 
		 }
	     }
	   }
	   

	   if(gAlice->TreeR())
	     {
	       Int_t nent=(Int_t)gAlice->TreeR()->GetEntries();
	       gAlice->TreeR()->GetEvent(nent-1);
	       TClonesArray *Rawclusters = pRICH->RawClustAddress(2);    //  Raw clusters branch
	       //printf ("Rawclusters:%p",Rawclusters);
	       Int_t nrawclusters = Rawclusters->GetEntriesFast();
	       	       
	       if (nrawclusters) {
		 printf("Raw Clusters    : %d\n",nrawclusters);
		 for (Int_t hit=0;hit<nrawclusters;hit++) {
		   AliRICHRawCluster* rcHit = (AliRICHRawCluster*) pRICH->RawClustAddress(2)->UncheckedAt(hit);
		   //Int_t nchamber = rcHit->fChamber;     // chamber number
		   //Int_t nhit = cHit->fHitNumber;        // hit number
		   Int_t qtot = rcHit->fQ;                 // charge
		   Float_t fx  =  rcHit->fX;                 // x-position
		   Float_t fy  =  rcHit->fY;                 // y-position
		   //Int_t type = rcHit->fCtype;             // cluster type ?   
		   Int_t mult = rcHit->fMultiplicity;      // How many pads form the cluster
		   pads += mult;
		   if (qtot > 0) {
		     //printf ("fx: %d, fy: %d\n",fx,fy);
		     if (fx>(x-4) && fx<(x+4)  && fy>(y-4) && fy<(y+4)) {
		       //printf("There %d \n",mult);
		       padmip+=mult;
		     } else {
		       padnumber->Fill(mult,(float) 1);
		       nraw++;
		       if (mult<4) Clcharge->Fill(qtot,(float) 1);
		     }
		     
		   }
		 }
	       }
	       
	       
	       TClonesArray *RecHits1D = pRICH->RecHitsAddress1D(2);
	       Int_t nrechits1D = RecHits1D->GetEntriesFast();
	       //printf (" nrechits:%d\n",nrechits);
	       
	       if(nrechits1D)
		 {
		   for (Int_t hit=0;hit<nrechits1D;hit++) {
		     AliRICHRecHit1D* recHit1D = (AliRICHRecHit1D*) pRICH->RecHitsAddress1D(2)->UncheckedAt(hit);
		     Float_t r_omega = recHit1D->fOmega;                  // Cerenkov angle
		     Float_t *cer_pho = recHit1D->fCerPerPhoton;        // Cerenkov angle per photon
		     Int_t *padsx = recHit1D->fPadsUsedX;           // Pads Used fo reconstruction (x)
		     Int_t *padsy = recHit1D->fPadsUsedY;           // Pads Used fo reconstruction (y)
		     Int_t goodPhotons = recHit1D->fGoodPhotons;    // Number of pads used for reconstruction
		     
		     Omega1D->Fill(r_omega,(float) 1);
		     
		     for (Int_t i=0; i<goodPhotons; i++)
		       {
			 PhotonCer->Fill(cer_pho[i],(float) 1);
			 PadsUsed->Fill(padsx[i],padsy[i],1);
			 //printf("Angle:%f, pad: %d %d\n",cer_pho[i],padsx[i],padsy[i]);
		       }
		     
		     //printf("Omega: %f, Theta: %f, Phi: %f\n",r_omega,r_theta,r_phi);
		   }
		 }

	       
	       TClonesArray *RecHits3D = pRICH->RecHitsAddress3D(2);
	       Int_t nrechits3D = RecHits3D->GetEntriesFast();
	       //printf (" nrechits:%d\n",nrechits);
	       
	       if(nrechits3D)
		 {
		   recEffEvent = 0;
		   
		   //for (Int_t hit=0;hit<nrechits3D;hit++) {
		   AliRICHRecHit3D* recHit3D = (AliRICHRecHit3D*) pRICH->RecHitsAddress3D(2)->UncheckedAt(track);
		   Float_t r_omega    = recHit3D->fOmega;                  // Cerenkov angle
		   Float_t r_theta    = recHit3D->fTheta;                  // Theta angle of incidence
		   Float_t r_phi      = recHit3D->fPhi;                    // Phi angle if incidence
		   Float_t meanradius = recHit3D->fMeanRadius;              // Mean radius for reconstructed point
		   Float_t originalOmega = recHit3D->fOriginalOmega;       // Real Cerenkov angle
		   Float_t originalTheta = recHit3D->fOriginalTheta;       // Real incidence angle
		   Float_t originalPhi = recHit3D->fOriginalPhi;           // Real azimuthal angle
		   
		   
		   //correction to track cerenkov angle
		   originalOmega = (Float_t) ckovangle->GetMean();
		   
		   if(diaglevel == 4)
		     {
		       printf("\nMean cerenkov angle: %f\n", originalOmega);
		       printf("Reconstructed cerenkov angle: %f\n",r_omega);
		     }
		   
		   Float_t omegaError = TMath::Abs(originalOmega - r_omega);
		   Float_t thetaError = TMath::Abs(originalTheta - r_theta);
		   Float_t phiError   = TMath::Abs(originalPhi - r_phi);
		   
		   
		   if(TMath::Abs(omegaError) < 0.015)
		     recEffEvent += 1;
		   
		   Omega3D->Fill(r_omega,(float) 1);
		   Theta->Fill(r_theta*180/TMath::Pi(),(float) 1);
		   Phi->Fill(r_phi*180/TMath::Pi()-180,(float) 1);
		   MeanRadius->Fill(meanradius,(float) 1);
		   identification->Fill(PTfinal, r_omega,1);
		   OriginalOmega->Fill(originalOmega, (float) 1);
		   OriginalTheta->Fill(originalTheta, (float) 1);
		   OriginalPhi->Fill(TMath::Abs(originalPhi), (float) 1);
		   OmegaError->Fill(omegaError, (float) 1);
		   ThetaError->Fill(thetaError, (float) 1);
		   PhiError->Fill(phiError, (float) 1);
		   
		   recEffEvent = recEffEvent;
		   recEffTotal += recEffEvent;
		   
		   Float_t pioncer = acos(sqrt((.139*.139+PTfinal*PTfinal)/(PTfinal*PTfinal*1.285*1.285)));
		   Float_t kaoncer = acos(sqrt((.439*.439+PTfinal*PTfinal)/(PTfinal*PTfinal*1.285*1.285)));
		   Float_t protoncer = acos(sqrt((.938*.938+PTfinal*PTfinal)/(PTfinal*PTfinal*1.285*1.285)));

		   Float_t piondist = TMath::Abs(r_omega - pioncer);
		   Float_t kaondist = TMath::Abs(r_omega - kaoncer);
		   Float_t protondist = TMath::Abs(r_omega - protoncer);

		   if(diaglevel == 4)
		     {
		       if(pioncer<r_omega)
			 {
			   printf("Identified as a PION!\n");
			   pionCount += 1;
			 }
		       if(kaoncer<r_omega && pioncer>r_omega)
			 {
			   if(kaondist>piondist)
			     {
			       printf("Identified as a PION!\n");
			       pionCount += 1;
			     }
			   else
			     {
			       printf("Identified as a KAON!\n");
			       kaonCount += 1;
			     }
			 }			 }
		       if(protoncer<r_omega && kaoncer>r_omega)
			 {
			   if(kaondist>protondist)
			     {
			       printf("Identified as a PROTON!\n");
			       protonCount += 1;
			     }
			   else
			     {
			       printf("Identified as a KAON!\n");
			       pionCount += 1;
			     }
			 }
		       if(protoncer>r_omega)
			 {
			   printf("Identified as a PROTON!\n");
			   protonCount += 1;
			 }

		       printf("\nReconstruction efficiency: %5.2f%%\n", recEffEvent*100);
		 }
	     }
       }
   
       
       for (Int_t nmothers=0;nmothers<ntracks;nmothers++){
	 totalphotonstrack->Fill(mothers[nmothers],(float) 1);
	 mother->Fill(mothers2[nmothers],(float) 1);
       }
       
       clusev->Fill(nraw,(float) 1);
       photev->Fill(phot,(float) 1);
       feedev->Fill(feed,(float) 1);
       padsmip->Fill(padmip,(float) 1);
       padscl->Fill(pads,(float) 1);
       phot = 0;
       feed = 0;
       pads = 0;
       nraw=0;
       padmip=0;
       
       
       
       gAlice->ResetDigits();
       gAlice->TreeD()->GetEvent(0);
       
       if (diaglevel < 4)
	 {
	   
	   
	   TClonesArray *Digits  = pRICH->DigitsAddress(2);
	   Int_t ndigits = Digits->GetEntriesFast();
	   printf("Digits          : %d\n",ndigits);
	   padsev->Fill(ndigits,(float) 1);
	   for (Int_t hit=0;hit<ndigits;hit++) {
	     AliRICHDigit* dHit = (AliRICHDigit*) Digits->UncheckedAt(hit);
	     Int_t qtot = dHit->Signal();                // charge
	     Int_t ipx  = dHit->PadX();               // pad number on X
	     Int_t ipy  = dHit->PadY();               // pad number on Y
	     //printf("%d, %d\n",ipx,ipy);
	     if( ipx<=100 && ipy <=100) hc0->Fill(ipx,ipy,(float) qtot);
	   }
	 }
       
       if (diaglevel == 5)
	 {
	   for (Int_t ich=0;ich<7;ich++)
	     {
	       TClonesArray *Digits = pRICH->DigitsAddress(ich);    //  Raw clusters branch
	       Int_t ndigits = Digits->GetEntriesFast();
	       //printf("Digits:%d\n",ndigits);
	       padsev->Fill(ndigits,(float) 1); 
	       if (ndigits) {
		 for (Int_t hit=0;hit<ndigits;hit++) {
		   AliRICHDigit* dHit = (AliRICHDigit*) Digits->UncheckedAt(hit);
		   Int_t qtot = dHit->Signal();                // charge
		   Int_t ipx  = dHit->PadX();               // pad number on X
		   Int_t ipy  = dHit->PadY();               // pad number on Y
		   if( ipx<=100 && ipy <=100 && ich==2) hc0->Fill(ipx,ipy,(float) qtot);
		   if( ipx<=162 && ipy <=162 && ich==0) hc1->Fill(ipx,ipy,(float) qtot);
		   if( ipx<=162 && ipy <=162 && ich==1) hc2->Fill(ipx,ipy,(float) qtot);
		   if( ipx<=162 && ipy <=162 && ich==2) hc3->Fill(ipx,ipy,(float) qtot);
		   if( ipx<=162 && ipy <=162 && ich==3) hc4->Fill(ipx,ipy,(float) qtot);
		   if( ipx<=162 && ipy <=162 && ich==4) hc5->Fill(ipx,ipy,(float) qtot);
		   if( ipx<=162 && ipy <=162 && ich==5) hc6->Fill(ipx,ipy,(float) qtot);
		   if( ipx<=162 && ipy <=162 && ich==6) hc7->Fill(ipx,ipy,(float) qtot);
		 }
	       }
	     }
	 }
   }
   
   if(diaglevel == 4)
     {

       Stat_t omegaE;
       Stat_t thetaE;
       Stat_t phiE;
       
       Stat_t omegaO;
       Stat_t thetaO;
       Stat_t phiO;
       
       for(Int_t i=0;i<99;i++)
	 {
	   omegaE = OriginalOmega->GetBinContent(i);
	   if(omegaE != 0)
	     {
	       omegaO = Omega3D->GetBinContent(i);
	       chiSquareOmega += (TMath::Power(omegaE,2) - TMath::Power(omegaO,2))/omegaO;
	     }

  	   thetaE = OriginalTheta->GetBinContent(i);
	   if(thetaE != 0)
	     {
	       thetaO = Theta->GetBinContent(i);
	       chiSquareTheta += (TMath::Power(thetaE,2) - TMath::Power(thetaO,2))/thetaO;
	     }

	   phiE = OriginalPhi->GetBinContent(i);
	   if(phiE != 0)
	     {
	       phiO = Phi->GetBinContent(i);
	       chiSquarePhi += (TMath::Power(phiE,2) - TMath::Power(phiO,2))/phiO;
	     }
	 }

       

       printf("\nChi square test values:   Omega - %f\n", chiSquareOmega);
       printf("                          Theta - %f\n", chiSquareTheta);
       printf("                          Phi   - %f\n", chiSquarePhi);
       
       printf("\nKolmogorov test values:   Omega - %5.4f\n", Omega3D->KolmogorovTest(OriginalOmega));
       printf("                          Theta - %5.4f\n", Theta->KolmogorovTest(OriginalTheta));
       printf("                          Phi   - %5.4f\n", Phi->KolmogorovTest(OriginalPhi));

       recEffTotal = recEffTotal/evNumber2;
       printf("\nTotal reconstruction efficiency: %5.2f%%\n", recEffTotal*100);
       printf("\n Pions: %d\n Kaons: %d\n Protons:%d\n",pionCount, kaonCount, protonCount);

     }
   
   
   //Create canvases, set the view range, show histograms

   TCanvas *c1 = 0;
   TCanvas *c2 = 0;
   TCanvas *c3 = 0;
   TCanvas *c4 = 0;
   TCanvas *c5 = 0;
   TCanvas *c6 = 0;
   TCanvas *c7 = 0;
   TCanvas *c8 = 0;
   TCanvas *c9 = 0;
   TCanvas *c10 = 0;
   TCanvas *c11 = 0;
   TCanvas *c12 = 0;
   TCanvas *c13 = 0;

   
   TStyle *mystyle=new TStyle("Plain","mystyle");
   mystyle->SetPalette(1,0);
   mystyle->SetFuncColor(2);
   mystyle->SetDrawBorder(0);
   mystyle->SetTitleBorderSize(0);
   mystyle->SetOptFit(1111);
   mystyle->cd();

   
   TClonesArray *RecHits3D = pRICH->RecHitsAddress3D(2);
   Int_t nrechits3D = RecHits3D->GetEntriesFast();
   TClonesArray *RecHits1D = pRICH->RecHitsAddress1D(2);
   Int_t nrechits1D = RecHits1D->GetEntriesFast();

  switch(diaglevel)
     {
     case 1:
       
       c1 = new TCanvas("c1","Alice RICH digits",50,50,300,350);
       hc0->SetXTitle("ix (npads)");
       hc0->Draw("colz");
	
       c2 = new TCanvas("c2","Hits per type",100,100,600,700);
       c2->Divide(2,2);
       //c4->SetFillColor(42);

       c2->cd(1);
       feedback->SetXTitle("x (cm)");
       feedback->SetYTitle("y (cm)");
       feedback->Draw("colz");
       
       c2->cd(2);
       //mip->SetFillColor(5);
       mip->SetXTitle("x (cm)");
       mip->SetYTitle("y (cm)");
       mip->Draw("colz");
       
       c2->cd(3);
       //cerenkov->SetFillColor(5);
       cerenkov->SetXTitle("x (cm)");
       cerenkov->SetYTitle("y (cm)"); 
       cerenkov->Draw("colz");
       
       c2->cd(4);
       //h->SetFillColor(5);
       h->SetXTitle("x (cm)");
       h->SetYTitle("y (cm)");
       h->Draw("colz");

       c3 = new TCanvas("c3","Hits distribution",150,150,600,350);
       c3->Divide(2,1);
       //c10->SetFillColor(42);
       
       c3->cd(1);
       hitsX->SetFillColor(5);
       hitsX->SetXTitle("(cm)");
       hitsX->Draw();
       
       c3->cd(2);
       hitsY->SetFillColor(5);
       hitsY->SetXTitle("(cm)");
       hitsY->Draw();
       
      
       break;
     case 2:
       
       c4 = new TCanvas("c4","Photon Spectra",50,50,600,350);
       c4->Divide(2,1);
       
       c4->cd(1);
       phspectra2->SetFillColor(5);
       phspectra2->SetXTitle("energy (eV)");
       phspectra2->Draw();
       c4->cd(2);
       phspectra1->SetFillColor(5);
       phspectra1->SetXTitle("energy (eV)");
       phspectra1->Draw();
       
       c5 = new TCanvas("c5","Particles Spectra",100,100,600,700);
       c5->Divide(2,2);
       
       c5->cd(1);
       pionspectra->SetFillColor(5);
       pionspectra->SetXTitle("(GeV)");
       pionspectra->Draw();
       
       c5->cd(2);
       protonspectra->SetFillColor(5);
       protonspectra->SetXTitle("(GeV)");
       protonspectra->Draw();
       
       c5->cd(3);
       kaonspectra->SetFillColor(5);
       kaonspectra->SetXTitle("(GeV)");
       kaonspectra->Draw();
       
       c5->cd(4);
       chargedspectra->SetFillColor(5);
       chargedspectra->SetXTitle("(GeV)");
       chargedspectra->Draw();

       break;
       
     case 3:

       
       if(gAlice->TreeR())
	 {
	   c6=new TCanvas("c6","Clusters Statistics",50,50,600,700);
	   c6->Divide(2,2);
	   
	   c6->cd(1);
	   Clcharge->SetFillColor(5);
	   Clcharge->SetXTitle("ADC counts");
	   if (evNumber2>10)
	     {
	       Clcharge->Fit("expo");
	     }
	   Clcharge->Draw();
	   
	   c6->cd(2);
	   padnumber->SetFillColor(5);
	   padnumber->SetXTitle("(counts)");
	   padnumber->Draw();
	   
	   c6->cd(3);
	   clusev->SetFillColor(5);
	   clusev->SetXTitle("(counts)");
	   if (evNumber2>10)
	     {
	       clusev->Fit("gaus");
	       //gaus->SetLineColor(2);
	       //gaus->SetLineWidth(3);
	     }
	   clusev->Draw();
	   
	   c6->cd(4);
	   padsmip->SetFillColor(5);
	   padsmip->SetXTitle("(counts)");
	   padsmip->Draw(); 
	 }
       
       if(evNumber2<1)
	 {
	   c11 = new TCanvas("c11","Cherenkov per Mip",400,10,600,700);
	   mother->SetFillColor(5);
	   mother->SetXTitle("counts");
	   mother->Draw();
	 }

       c7 = new TCanvas("c7","Production Statistics",100,100,600,700);
       c7->Divide(2,2);
       //c7->SetFillColor(42);
       
       c7->cd(1);
       totalphotonsevent->SetFillColor(5);
       totalphotonsevent->SetXTitle("Photons (counts)");
       if (evNumber2>10)
	   {
	     totalphotonsevent->Fit("gaus");
	     //gaus->SetLineColor(2);
	     //gaus->SetLineWidth(3);
	   }
       totalphotonsevent->Draw();
       
       c7->cd(2);
       photev->SetFillColor(5);
       photev->SetXTitle("(counts)");
       if (evNumber2>10)
	 {
	   photev->Fit("gaus");
	   //gaus->SetLineColor(2);
	   //gaus->SetLineWidth(3);
	 }
       photev->Draw();
       
       c7->cd(3);
       feedev->SetFillColor(5);
       feedev->SetXTitle("(counts)");
       if (evNumber2>10)
	 {
	   feedev->Fit("gaus");
	 }
       feedev->Draw();

       c7->cd(4);
       padsev->SetFillColor(5);
       padsev->SetXTitle("(counts)");
       if (evNumber2>10)
	 {
	   padsev->Fit("gaus");
	 }
       padsev->Draw();

       break;

     case 4:
       

       if(nrechits3D)
	 {
	   c8 = new TCanvas("c8","3D reconstruction of Phi angle",50,50,300,1050);
	   c8->Divide(1,3);
	   //c2->SetFillColor(42);
	   
	   
	   // data per hit
	   c8->cd(1);
	   hitsPhi->SetFillColor(5);
	   if (evNumber2>10)
	     hitsPhi->Fit("gaus");
	   hitsPhi->Draw();
	   
	    //data per track
	   c8->cd(2);
	   OriginalPhi->SetFillColor(5);
	   if (evNumber2>10)
	     OriginalPhi->Fit("gaus");
	   OriginalPhi->Draw();

	   //recontructed data
	   c8->cd(3);
	   Phi->SetFillColor(5);
	   if (evNumber2>10)
	     Phi->Fit("gaus");
	   Phi->Draw();

	   c9 = new TCanvas("c9","3D reconstruction of theta angle",75,75,300,1050);
	   c9->Divide(1,3);

	   // data per hit
	   c9->cd(1);
	   hitsTheta->SetFillColor(5);
	   if (evNumber2>10)
	     hitsTheta->Fit("gaus");
	   hitsTheta->Draw();
	   
	   //data per track
	   c9->cd(2);
	   OriginalTheta->SetFillColor(5);
	   if (evNumber2>10)
	     OriginalTheta->Fit("gaus");
	   OriginalTheta->Draw();

	   //recontructed data
	   c9->cd(3);
	   Theta->SetFillColor(5);
	   if (evNumber2>10)
	     Theta->Fit("gaus");
	   Theta->Draw();

	   c10 = new TCanvas("c10","3D reconstruction of cherenkov angle",100,100,300,1050);
	   c10->Divide(1,3);

	   // data per hit
	   c10->cd(1);
	   ckovangle->SetFillColor(5);
	   ckovangle->SetXTitle("angle (radians)");
	   if (evNumber2>10)
	     ckovangle->Fit("gaus");
	   ckovangle->Draw();
	   
	   //data per track
	   c10->cd(2);
	   OriginalOmega->SetFillColor(5);
	   OriginalOmega->SetXTitle("angle (radians)");
	   if (evNumber2>10)
	     OriginalOmega->Fit("gaus");
	   OriginalOmega->Draw();

	   //recontructed data
	   c10->cd(3);
	   Omega3D->SetFillColor(5);
	   Omega3D->SetXTitle("angle (radians)");
	   if (evNumber2>10)
	     Omega3D->Fit("gaus");
	   Omega3D->Draw(); 


	   c11 = new TCanvas("c11","3D reconstruction of mean radius",125,125,300,700);
	   c11->Divide(1,2);

	   // data per hit
	   c11->cd(1);
	   radius->SetFillColor(5);
	   radius->SetXTitle("radius (cm)");
	   radius->Draw();

	   //recontructed data
	   c11->cd(2);
	   MeanRadius->SetFillColor(5);
	   MeanRadius->SetXTitle("radius (cm)");
	   MeanRadius->Draw();

	   
	   c12 = new TCanvas("c12","Cerenkov angle vs. Momentum",150,150,550,350);

	   c12->cd(1);
	   identification->SetFillColor(5);
	   identification->SetXTitle("Momentum (GeV/c)");
	   identification->SetYTitle("Cherenkov angle (radians)");
	   
	   TF1 *pionplot = new TF1("pion","acos(sqrt((.139*.139+x*x)/(x*x*1.285*1.285)))",1,5);
	   TF1 *kaonplot = new TF1("kaon","acos(sqrt((.439*.439+x*x)/(x*x*1.285*1.285)))",1,5);
	   TF1 *protonplot = new TF1("proton","acos(sqrt((.938*.938+x*x)/(x*x*1.285*1.285)))",1,5);
	   
	   identification->Draw();

	   pionplot->SetLineColor(5);
	   pionplot->Draw("same");

	   kaonplot->SetLineColor(4);
	   kaonplot->Draw("same");

	   protonplot->SetLineColor(3);
	   protonplot->Draw("same");

	   c13 = new TCanvas("c13","Reconstruction Errors",200,200,900,350);
	   c13->Divide(3,1);

	   c13->cd(1);
	   PhiError->SetFillColor(5);
	   if (evNumber2>10)
	     PhiError->Fit("gaus");
	   PhiError->Draw();
	   c13->cd(2);
	   ThetaError->SetFillColor(5);
	   if (evNumber2>10)
	     ThetaError->Fit("gaus");
	   ThetaError->Draw();
	   c13->cd(3);
	   OmegaError->SetFillColor(5);
	   OmegaError->SetXTitle("angle (radians)");
	   if (evNumber2>10)
	     OmegaError->Fit("gaus");
	   OmegaError->Draw();
	   
	 }
       
       if(nrechits1D)
	 {
	   c9 = new TCanvas("c9","1D Reconstruction",100,100,1100,700);
	   c9->Divide(3,2);
	   //c5->SetFillColor(42);
	   
	   c9->cd(1);
	   ckovangle->SetFillColor(5);
	   ckovangle->SetXTitle("angle (radians)");
	   ckovangle->Draw();
	   
	   c9->cd(2);
	   radius->SetFillColor(5);
	   radius->SetXTitle("radius (cm)");
	   radius->Draw();
	   
	   c9->cd(3);
	   hc0->SetXTitle("pads");
	   hc0->Draw("box"); 
	   
	   c9->cd(5);
	   Omega1D->SetFillColor(5);
	   Omega1D->SetXTitle("angle (radians)");
	   Omega1D->Draw();
	   
	   c9->cd(4);
	   PhotonCer->SetFillColor(5);
	   PhotonCer->SetXTitle("angle (radians)");
	   PhotonCer->Draw();
	   
	   c9->cd(6);
	   PadsUsed->SetXTitle("pads");
	   PadsUsed->Draw("box"); 
	 }
       
       break;
       
     case 5:
       
       printf("Drawing histograms.../n");

       c10 = new TCanvas("c10","Alice RICH digits",50,50,1200,700);
       c1->Divide(4,2);
       
       c10->cd(1);
       hc1->SetXTitle("ix (npads)");
       hc1->Draw("box");
       c10->cd(2);
       hc2->SetXTitle("ix (npads)");
       hc2->Draw("box");
       c10->cd(3);
       hc3->SetXTitle("ix (npads)");
       hc3->Draw("box");
       c10->cd(4);
       hc4->SetXTitle("ix (npads)");
       hc4->Draw("box");
       c10->cd(5);
       hc5->SetXTitle("ix (npads)");
       hc5->Draw("box");
       c10->cd(6);
       hc6->SetXTitle("ix (npads)");
       hc6->Draw("box");
       c10->cd(7);
       hc7->SetXTitle("ix (npads)");
       hc7->Draw("box");
       c10->cd(8);
       hc0->SetXTitle("ix (npads)");
       hc0->Draw("box");
       c11 = new TCanvas("c11","Hits per type",100,100,600,700);
       c11->Divide(2,2);
       
       c11->cd(1);
       feedback->SetXTitle("x (cm)");
       feedback->SetYTitle("y (cm)");
       feedback->Draw();
       
       c11->cd(2);
       mip->SetXTitle("x (cm)");
       mip->SetYTitle("y (cm)");
       mip->Draw();
       
       c11->cd(3);
       cerenkov->SetXTitle("x (cm)");
       cerenkov->SetYTitle("y (cm)"); 
       cerenkov->Draw();
       
       c11->cd(4);
       h->SetXTitle("x (cm)");
       h->SetYTitle("y (cm)");
       h->Draw();

       c12 = new TCanvas("c12","Hits distribution",150,150,600,350);
       c12->Divide(2,1);
       
       c12->cd(1);
       hitsX->SetFillColor(5);
       hitsX->SetXTitle("(cm)");
       hitsX->Draw();
       
       c12->cd(2);
       hitsY->SetFillColor(5);
       hitsY->SetXTitle("(cm)");
       hitsY->Draw();
       
       break;
       
     }
       

   printf("\nEnd of analysis\n");
   printf("**********************************\n");
}//void AliRICHv3::DiagnosticsSE(Int_t diaglevel,Int_t evNumber1,Int_t evNumber2)

//__________________________________________________________________________________________________
void AliRICHv3::MakeBranch(Option_t* option)
{//Create Tree branches for the RICH.
  if(GetDebug())Info("MakeBranch","Start with option= %s.",option);
    
  const Int_t kBufferSize = 4000;
  char branchname[20];
      
   
  const char *cH = strstr(option,"H");
  const char *cD = strstr(option,"D");
  const char *cR = strstr(option,"R");
  const char *cS = strstr(option,"S");


  if(cH&&TreeH()){
    if(!fHits) fHits=new TClonesArray("AliRICHhit",1000  );
    if(!fCerenkovs) fCerenkovs  = new TClonesArray("AliRICHCerenkov",1000);
    MakeBranchInTree(TreeH(),"RICHCerenkov", &fCerenkovs, kBufferSize, 0) ;

    if(!fSDigits) fSDigits    = new TClonesArray("AliRICHdigit",100000);
    MakeBranchInTree(TreeH(),"RICHSDigits", &fSDigits, kBufferSize, 0) ;
  }     
  AliDetector::MakeBranch(option);//this is after cH because we need to guarantee that fHits array is created
      
  if(cS&&fLoader->TreeS()){  
    if(!fSDigits) fSDigits=new TClonesArray("AliRICHdigit",100000);
    MakeBranchInTree(fLoader->TreeS(),"RICH",&fSDigits,kBufferSize,0) ;
  }
   
  int i;
  if (cD&&fLoader->TreeD()){
    if(!fDchambers){
      fDchambers=new TObjArray(kNCH);    // one branch for digits per chamber
      for(i=0;i<kNCH;i++){ 
        fDchambers->AddAt(new TClonesArray("AliRICHDigit",10000), i); 
      }       
    }
    for (i=0; i<kNCH ;i++) 
      {
        sprintf(branchname,"%sDigits%d",GetName(),i+1);	
        MakeBranchInTree(fLoader->TreeD(),branchname, &((*fDchambers)[i]), kBufferSize, 0);
      }
   }

  if (cR&&gAlice->TreeR()){//one branch for raw clusters per chamber
    Int_t i;
    if (fRawClusters == 0x0 ) 
     {
       fRawClusters = new TObjArray(kNCH);
       for (i=0; i<kNCH ;i++) 
         {
           fRawClusters->AddAt(new TClonesArray("AliRICHRawCluster",10000), i); 
         }
     }
     
    if (fRecHits1D == 0x0) 
     {
        fRecHits1D = new TObjArray(kNCH);
        for (i=0; i<kNCH ;i++) 
         {
          fRecHits1D->AddAt(new TClonesArray("AliRICHRecHit1D",1000), i);
         }
     }

    if (fRecHits3D == 0x0) 
     {
        fRecHits3D = new TObjArray(kNCH);
        for (i=0; i<kNCH ;i++) 
         {
          fRecHits3D->AddAt(new TClonesArray("AliRICHRecHit3D",1000), i);
         }
     }
       
    for (i=0; i<kNCH ;i++){
       sprintf(branchname,"%sRawClusters%d",GetName(),i+1);      
       MakeBranchInTree(gAlice->TreeR(),branchname, &((*fRawClusters)[i]), kBufferSize, 0);
       sprintf(branchname,"%sRecHits1D%d",GetName(),i+1);
       MakeBranchInTree(fLoader->TreeR(),branchname, &((*fRecHits1D)[i]), kBufferSize, 0);
       sprintf(branchname,"%sRecHits3D%d",GetName(),i+1);  
       MakeBranchInTree(fLoader->TreeR(),branchname, &((*fRecHits3D)[i]), kBufferSize, 0);
     }
   }//if (cR && gAlice->TreeR())
  if(GetDebug())Info("MakeBranch","Stop.");   
}
//______________________________________________________________________________
void AliRICHv3::SetTreeAddress()
{//Set branch address for the Hits and Digits Tree.
  if(GetDebug())Info("SetTreeAddress","Start.");
  
  char branchname[20];
  Int_t i;

    
  TBranch *branch;
  TTree *treeH = fLoader->TreeH();
  TTree *treeD = fLoader->TreeD();
  TTree *treeR = fLoader->TreeR();
  TTree *treeS = fLoader->TreeS();
    
  if(treeH){
    if(GetDebug())Info("SetTreeAddress","tree H is requested.");
    if(fHits==0x0) fHits=new TClonesArray("AliRICHhit",1000); 
    
    branch = treeH->GetBranch("RICHCerenkov");
    if(branch){
      if (fCerenkovs == 0x0) fCerenkovs  = new TClonesArray("AliRICHCerenkov",1000); 
        branch->SetAddress(&fCerenkovs);
    }
       
      branch = treeH->GetBranch("RICHSDigits");
      if (branch) 
       {
         if (fSDigits == 0x0) fSDigits    = new TClonesArray("AliRICHdigit",100000);
         branch->SetAddress(&fSDigits);
       }
  }//if(treeH)
 
   //this is after TreeH because we need to guarantee that fHits array is created
  AliDetector::SetTreeAddress();
    
  if(treeS){
    if(GetDebug())Info("SetTreeAddress","tree S is requested.");
    branch = treeS->GetBranch("RICH");
    if(branch){
      if(!fSDigits) fSDigits=new TClonesArray("AliRICHdigit",100000);
      branch->SetAddress(&fSDigits);
    }
  }
    
    
  if(treeD){
    if(GetDebug())Info("SetTreeAddress","tree D is requested.");

      if (fDchambers == 0x0) 
        {
           fDchambers = new TObjArray(kNCH);
           for (i=0; i<kNCH ;i++) 
             {
               fDchambers->AddAt(new TClonesArray("AliRICHDigit",10000), i); 
             }
        }
      
      for (i=0; i<kNCH; i++) {
        sprintf(branchname,"%sDigits%d",GetName(),i+1);
        if (fDchambers) {
           branch = treeD->GetBranch(branchname);
           if (branch) branch->SetAddress(&((*fDchambers)[i]));
        }
      }
    }
    
  if(treeR){
    if(GetDebug())Info("SetTreeAddress","tree R is requested.");

    if (fRawClusters == 0x0 ) 
     {
       fRawClusters = new TObjArray(kNCH);
       for (i=0; i<kNCH ;i++) 
         {
           fRawClusters->AddAt(new TClonesArray("AliRICHRawCluster",10000), i); 
         }
     }
     
    if (fRecHits1D == 0x0) 
     {
        fRecHits1D = new TObjArray(kNCH);
        for (i=0; i<kNCH ;i++) 
         {
          fRecHits1D->AddAt(new TClonesArray("AliRICHRecHit1D",1000), i);
         }
     }

    if (fRecHits3D == 0x0) 
     {
        fRecHits3D = new TObjArray(kNCH);
        for (i=0; i<kNCH ;i++) 
         {
          fRecHits3D->AddAt(new TClonesArray("AliRICHRecHit3D",1000), i);
         }
     }
    
    for (i=0; i<kNCH; i++) {
	  sprintf(branchname,"%sRawClusters%d",GetName(),i+1);
	  if (fRawClusters) {
	      branch = treeR->GetBranch(branchname);
	      if (branch) branch->SetAddress(&((*fRawClusters)[i]));
	  }
    }
      
    for (i=0; i<kNCH; i++) {
	sprintf(branchname,"%sRecHits1D%d",GetName(),i+1);
	if (fRecHits1D) {
	  branch = treeR->GetBranch(branchname);
	  if (branch) branch->SetAddress(&((*fRecHits1D)[i]));
	  }
     }
      
     for (i=0; i<kNCH; i++) {
	sprintf(branchname,"%sRecHits3D%d",GetName(),i+1);
	if (fRecHits3D) {
	  branch = treeR->GetBranch(branchname);
	  if (branch) branch->SetAddress(&((*fRecHits3D)[i]));
	  }
      } 
      
  }//if(treeR)
  if(GetDebug())Info("SetTreeAddress","Stop.");
}//void AliRICHv3::SetTreeAddress()
