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



////////////////////////////////////////////////////////
//  Manager and hits classes for set:RICH version 0    //
/////////////////////////////////////////////////////////

#include <TTUBE.h>
#include <TNode.h> 
#include <TRandom.h> 

#include "AliRICHv0.h"
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
    //fChambers = 0;
}

//___________________________________________
AliRICHv0::AliRICHv0(const char *name, const char *title)
    : AliRICH(name,title)
{
    fCkov_number=0;
    fFreon_prod=0;

    fChambers = new TObjArray(7);
    for (Int_t i=0; i<7; i++) {
    
	(*fChambers)[i] = new AliRICHChamber();  
	
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

  AliRICH *RICH = (AliRICH *) gAlice->GetDetector("RICH"); 
  AliRICHSegmentation*  segmentation;
  AliRICHGeometry*  geometry;
  AliRICHChamber*       iChamber;

  iChamber = &(RICH->Chamber(0));
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

    printf("Nspacers: %d", nspacers);
    
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

    printf("Position of the gap: %f to %f\n", 1.276 + geometry->GetGapThickness()/2 - geometry->GetProximityGapThickness()/2 - .2, 1.276 + geometry->GetGapThickness()/2 - geometry->GetProximityGapThickness()/2 + .2);
    
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
    Int_t   ISXFLD = gAlice->Field()->Integ();
    Float_t SXMGMX = gAlice->Field()->Max();
    Int_t i;

    /************************************Antonnelo's Values (14-vectors)*****************************************/
    /*
    Float_t ppckov[14] = { 5.63e-9,5.77e-9,5.9e-9,6.05e-9,6.2e-9,6.36e-9,6.52e-9,
			   6.7e-9,6.88e-9,7.08e-9,7.3e-9,7.51e-9,7.74e-9,8e-9 };
    Float_t rindex_quarz[14] = { 1.528309,1.533333,
				 1.538243,1.544223,1.550568,1.55777,
				 1.565463,1.574765,1.584831,1.597027,
			       1.611858,1.6277,1.6472,1.6724 };
    Float_t rindex_quarzo[14] = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };
    Float_t rindex_methane[14] = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };
    Float_t rindex_gri[14] = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };
    Float_t absco_freon[14] = { 179.0987,179.0987,
				179.0987,179.0987,179.0987,142.92,56.65,13.95,10.43,7.07,2.03,.5773,.33496,0. };
    //Float_t absco_freon[14] = { 1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,
	//			 1e-5,1e-5,1e-5,1e-5,1e-5 };
    Float_t absco_quarz[14] = { 64.035,39.98,35.665,31.262,27.527,22.815,21.04,17.52,
				14.177,9.282,4.0925,1.149,.3627,.10857 };
    Float_t absco_quarzo[14] = { 1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,
				 1e-5,1e-5,1e-5,1e-5,1e-5 };
    Float_t absco_csi[14] = { 1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,
			      1e-4,1e-4,1e-4,1e-4 };
    Float_t absco_methane[14] = { 1e6,1e6,1e6,1e6,1e6,1e6,1e6,1e6,1e6,1e6,1e6,
				  1e6,1e6,1e6 };
    Float_t absco_gri[14] = { 1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,
			      1e-4,1e-4,1e-4,1e-4 };
    Float_t effic_all[14] = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };
    Float_t effic_csi[14] = { 6e-4,.005,.0075,.01125,.045,.117,.135,.16575,
			      .17425,.1785,.1836,.1904,.1938,.221 };
    Float_t effic_gri[14] = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };
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
    Float_t rindex_quarz[26];
    Float_t  e1= 10.666;
    Float_t  e2= 18.125;
    Float_t  f1= 46.411;
    Float_t  f2= 228.71;
    for (i=0;i<26;i++)
    {
	Float_t ene=ppckov[i]*1e9;
	Float_t a=f1/(e1*e1 - ene*ene);
	Float_t b=f2/(e2*e2 - ene*ene);
	rindex_quarz[i] = TMath::Sqrt(1. + a + b );
	//printf ("Rindex_quarz: %e\n",rindex_quarz[i]);
    } 
    
    //Refraction index for opaque quarz, methane and grid
    Float_t rindex_quarzo[26];
    Float_t rindex_methane[26];
    Float_t rindex_gri[26];
    for (i=0;i<26;i++)
    {
	rindex_quarzo[i]=1;
	rindex_methane[i]=1.000444;
	rindex_gri[i]=1;
	//printf ("Rindex_quarzo , etc: %e, %e, %e\n",rindex_quarzo[i], rindex_methane[i], rindex_gri[i]=1);
    } 
    
    //Absorption index for freon
    Float_t absco_freon[26] = {179.0987, 179.0987, 179.0987, 179.0987, 179.0987,  179.0987, 179.0987, 179.0987, 
	 		       179.0987, 142.9206, 56.64957, 25.58622, 13.95293, 12.03905, 10.42953, 8.804196, 
			       7.069031, 4.461292, 2.028366, 1.293013, .577267,   .40746,  .334964, 0., 0., 0.};
    
    //Absorption index for quarz
    /*Float_t Qzt [21] = {.0,.0,.005,.04,.35,.647,.769,.808,.829,.844,.853,.858,.869,.887,.903,.902,.902,
	 		.906,.907,.907,.907};
    Float_t Wavl2[] = {150.,155.,160.0,165.0,170.0,175.0,180.0,185.0,190.0,195.0,200.0,205.0,210.0,
	 	       215.0,220.0,225.0,230.0,235.0,240.0,245.0,250.0};		 		 
    Float_t absco_quarz[31];	     
    for (Int_t i=0;i<31;i++)
    {
	Float_t Xlam = 1237.79 / (ppckov[i]*1e9);
	if (Xlam <= 160) absco_quarz[i] = 0;
	if (Xlam > 250) absco_quarz[i] = 1;
	else 
	{
	    for (Int_t j=0;j<21;j++)
	    {
		//printf ("Passed\n");
		if (Xlam > Wavl2[j] && Xlam < Wavl2[j+1])
		{
		    Float_t Dabs = (Qzt[j+1] - Qzt[j])/(Wavl2[j+1] - Wavl2[j]);
		    Float_t Abso = Qzt[j] + Dabs*(Xlam - Wavl2[j]);
		    absco_quarz[i] = -5.0/(TMath::Log(Abso));
		} 
	    }
	}
	printf ("Absco_quarz: %e Absco_freon: %e for energy: %e\n",absco_quarz[i],absco_freon[i],ppckov[i]);
    }*/

    /*Float_t absco_quarz[31] = {49.64211, 48.41296, 47.46989, 46.50492, 45.13682, 44.47883, 43.1929 , 41.30922, 40.5943 ,
			       39.82956, 38.98623, 38.6247 , 38.43448, 37.41084, 36.22575, 33.74852, 30.73901, 24.25086, 
			       17.94531, 11.88753, 5.99128,  3.83503,  2.36661,  1.53155, 1.30582, 1.08574, .8779708, 
			       .675275, 0., 0., 0.};
    
    for (Int_t i=0;i<31;i++)
    {
	absco_quarz[i] = absco_quarz[i]/10;
    }*/

    Float_t absco_quarz [26] = {105.8, 65.52, 48.58, 42.85, 35.79, 31.262, 28.598, 27.527, 25.007, 22.815, 21.004,
				19.266, 17.525, 15.878, 14.177, 11.719, 9.282, 6.62, 4.0925, 2.601, 1.149, .667, .3627,
				.192, .1497, .10857};
    
    //Absorption index for methane
    Float_t absco_methane[26];
    for (i=0;i<26;i++) 
    {
	absco_methane[i]=AbsoCH4(ppckov[i]*1e9); 
	//printf("Absco_methane: %e for energy: %e\n", absco_methane[i],ppckov[i]*1e9);
    }
    
    //Absorption index for opaque quarz, csi and grid, efficiency for all and grid
    Float_t absco_quarzo[26];
    Float_t absco_csi[26];
    Float_t absco_gri[26];
    Float_t effic_all[26];
    Float_t effic_gri[26];
    for (i=0;i<26;i++)
    { 
	absco_quarzo[i]=1e-5; 
	absco_csi[i]=1e-4; 
	absco_gri[i]=1e-4; 
	effic_all[i]=1; 
	effic_gri[i]=1;
	//printf ("All must be 1: %e,  %e,  %e,  %e,  %e\n",absco_quarzo[i],absco_csi[i],absco_gri[i],effic_all[i],effic_gri[i]);
    } 
    
    //Efficiency for csi 
    
    Float_t effic_csi[26] = {0.000199999995, 0.000600000028, 0.000699999975, 0.00499999989, 0.00749999983, 0.010125,
			     0.0242999997, 0.0405000001, 0.0688500032, 0.105299994, 0.121500008, 0.141749993, 0.157949999,
			     0.162, 0.166050002, 0.167669997, 0.174299985, 0.176789999, 0.179279998, 0.182599992, 0.18592,
			     0.187579989, 0.189239994, 0.190899998, 0.207499996, 0.215799987};
	
    

    //FRESNEL LOSS CORRECTION FOR PERPENDICULAR INCIDENCE AND
    //UNPOLARIZED PHOTONS

    for (i=0;i<26;i++)
    {
	effic_csi[i] = effic_csi[i]/(1.-Fresnel(ppckov[i]*1e9,1.,0)); 
	//printf ("Fresnel result: %e for energy: %e\n",Fresnel(ppckov[i]*1e9,1.,0),ppckov[i]*1e9);
    }
	
    /*******************************************End of rich_media.f***************************************/

  

    
    
    
    Float_t afre[2], agri, amet[2], aqua[2], ahon, zfre[2], zgri, zhon, 
    zmet[2], zqua[2];
    Int_t nlmatfre;
    Float_t densquao;
    Int_t nlmatmet, nlmatqua;
    Float_t wmatquao[2], rindex_freon[26];
    Float_t aquao[2], epsil, stmin, zquao[2];
    Int_t nlmatquao;
    Float_t radlal, densal, tmaxfd, deemax, stemax;
    Float_t aal, zal, radlgri, densfre, radlhon, densgri, denshon,densqua, densmet, wmatfre[2], wmatmet[2], wmatqua[2];
    
    Int_t *idtmed = fIdtmed->GetArray()-999;
    
    TGeant3 *geant3 = (TGeant3*) gMC;
    
    // --- Photon energy (GeV) 
    // --- Refraction indexes 
    for (i = 0; i < 26; ++i) {
	rindex_freon[i] = ppckov[i] * .0172 * 1e9 + 1.177;
	//printf ("Rindex_freon: %e \n Effic_csi: %e for energy: %e\n",rindex_freon[i], effic_csi[i], ppckov[i]);
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
    
    AliMedium(1, "DEFAULT MEDIUM AIR$", 1, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(2, "HONEYCOMB$", 6, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(3, "QUARZO$", 20, 1, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(4, "FREON$", 30, 1, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(5, "METANO$", 40, 1, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(6, "CSI$", 16, 1, ISXFLD, SXMGMX,tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(7, "GRIGLIA$", 11, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(8, "QUARZOO$", 21, 1, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
    AliMedium(9, "GAP$", 41, 1, ISXFLD, SXMGMX,tmaxfd, .1, -deemax, epsil, -stmin);
    AliMedium(10, "ALUMINUM$", 50, 1, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
    

    geant3->Gsckov(idtmed[1000], 26, ppckov, absco_methane, effic_all, rindex_methane);
    geant3->Gsckov(idtmed[1001], 26, ppckov, absco_methane, effic_all, rindex_methane);
    geant3->Gsckov(idtmed[1002], 26, ppckov, absco_quarz, effic_all,rindex_quarz);
    geant3->Gsckov(idtmed[1003], 26, ppckov, absco_freon, effic_all,rindex_freon);
    geant3->Gsckov(idtmed[1004], 26, ppckov, absco_methane, effic_all, rindex_methane);
    geant3->Gsckov(idtmed[1005], 26, ppckov, absco_csi, effic_csi, rindex_methane);
    geant3->Gsckov(idtmed[1006], 26, ppckov, absco_gri, effic_gri, rindex_gri);
    geant3->Gsckov(idtmed[1007], 26, ppckov, absco_quarzo, effic_all, rindex_quarzo);
    geant3->Gsckov(idtmed[1008], 26, ppckov, absco_methane, effic_all, rindex_methane);
    geant3->Gsckov(idtmed[1009], 26, ppckov, absco_gri, effic_gri, rindex_gri);
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

    //LOSCH,SCH4(9),WL(9),EM(9),ALENGTH(31)
    Float_t sch4[9] = {.12,.16,.23,.38,.86,2.8,7.9,28.,80.};              //MB X 10^22
    //Float_t wl[9] = {153.,152.,151.,150.,149.,148.,147.,146.,145};
    Float_t em[9] = {8.1,8.158,8.212,8.267,8.322,8.378,8.435,8.493,8.55};
    const Float_t losch=2.686763E19;                                      // LOSCHMIDT NUMBER IN CM-3
    const Float_t igas1=100, igas2=0, oxy=10., wat=5., pre=750.,tem=283.;                                      
    Float_t pn=pre/760.;
    Float_t tn=tem/273.16;
    
	
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
    
    Float_t dm=(igas1/100.)*(1.-((oxy+wat)/1.e6))*losch*pn/tn;
    Float_t abslm=1./sm/dm;
    
//    ------- ISOBUTHANE CROSS SECTION --------------
//     i-C4H10 (ai) abs. length from curves in
//     Lu-McDonald paper for BARI RICH workshop .
//     -----------------------------------------------------------
    
    Float_t ai;
    Float_t absli;
    if (igas2 != 0) 
    {
	if (x<7.25)
	    ai=100000000.;
	
	if(x>=7.25 && x<7.375)
	    ai=24.3;
	
	if(x>=7.375)
	    ai=.0000000001;
	
	Float_t si = 1./(ai*losch*273.16/293.);                    // ISOB. CRO.SEC.IN CM2
	Float_t di=(igas2/100.)*(1.-((oxy+wat)/1.e6))*losch*pn/tn;
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
	
	Float_t dox=(oxy/1e6)*losch*pn/tn;
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
	Float_t dw=(wat/1e6)*losch*pn/tn;
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
    printf("\n\n\n Start Init for version 0 - CPC chamber type \n\n\n");
    
    // 
    // Initialize Tracking Chambers
    //
    for (Int_t i=1; i<7; i++) {
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
    
    printf("\n\n\n Finished Init for version 0 - CPC chamber type\n\n\n");
}

//___________________________________________
void AliRICHv0::StepManager()
{
    Int_t          copy, id;
    static Int_t   idvol;
    static Int_t   vol[2];
    Int_t          ipart;
    static Float_t hits[18];
    static Float_t Ckov_data[19];
    TLorentzVector Position;
    TLorentzVector Momentum;
    Float_t        pos[3];
    Float_t        mom[4];
    Float_t        Localpos[3];
    Float_t        Localmom[4];
    Float_t        Localtheta,Localphi;
    Float_t        theta,phi;
    Float_t        destep, step;
    Float_t        ranf[2];
    Int_t          NPads;
    Float_t        coscerenkov;
    static Float_t eloss, xhit, yhit, tlength;
    const  Float_t big=1.e10;
       
    TClonesArray &lhits = *fHits;
    TGeant3 *geant3 = (TGeant3*) gMC;
    TParticle *current = (TParticle*)(*gAlice->Particles())[gAlice->CurrentTrack()];

 //if (current->Energy()>1)
   //{
        
    // Only gas gap inside chamber
    // Tag chambers and record hits when track enters 
    
    idvol=-1;
    id=gMC->CurrentVolID(copy);
    Float_t cherenkov_loss=0;
    //gAlice->KeepTrack(gAlice->CurrentTrack());
    
    gMC->TrackPosition(Position);
    pos[0]=Position(0);
    pos[1]=Position(1);
    pos[2]=Position(2);
    Ckov_data[1] = pos[0];                 // X-position for hit
    Ckov_data[2] = pos[1];                 // Y-position for hit
    Ckov_data[3] = pos[2];                 // Z-position for hit
    //Ckov_data[11] = gAlice->CurrentTrack();

    AliRICH *RICH = (AliRICH *) gAlice->GetDetector("RICH"); 
    
    /********************Store production parameters for Cerenkov photons************************/ 
//is it a Cerenkov photon? 
    if (gMC->TrackPid() == 50000050) {          

      //if (gMC->VolId("GAP ")==gMC->CurrentVolID(copy))
        //{                    
	  Float_t Ckov_energy = current->Energy();
	  //energy interval for tracking
	  if  (Ckov_energy > 5.6e-09 && Ckov_energy < 7.8e-09 )       
	    //if (Ckov_energy > 0)
	    {
	      if (gMC->IsTrackEntering()){                                     //is track entering?
		if (gMC->VolId("FRE1")==gMC->CurrentVolID(copy) || gMC->VolId("FRE2")==gMC->CurrentVolID(copy))
		  {                                                          //is it in freo?
		    if (geant3->Gctrak()->nstep<1){                          //is it the first step?
		      Int_t mother = current->GetFirstMother(); 
		      
		      //printf("Second Mother:%d\n",current->GetSecondMother());
		      
		      Ckov_data[10] = mother;
		      Ckov_data[11] = gAlice->CurrentTrack();
		      Ckov_data[12] = 1;             //Media where photon was produced 1->Freon, 2->Quarz
		      fCkov_number++;
		      fFreon_prod=1;
		      //printf("Index: %d\n",fCkov_number);
		    }    //first step question
		  }        //freo question
		
		if (geant3->Gctrak()->nstep<1){                                  //is it first step?
		  if (gMC->VolId("QUAR")==gMC->CurrentVolID(copy))             //is it in quarz?
		    {
		      Ckov_data[12] = 2;
		    }    //quarz question
		}        //first step question
		
		//printf("Before %d\n",fFreon_prod);
	      }   //track entering question
	      
	      if (Ckov_data[12] == 1)                                        //was it produced in Freon?
		//if (fFreon_prod == 1)
		{
		  if (gMC->IsTrackEntering()){                                     //is track entering?
		    //printf("Got in");
		    if (gMC->VolId("META")==gMC->CurrentVolID(copy))                //is it in gap?      
		      {
			//printf("Got in\n");
			gMC->TrackMomentum(Momentum);
			mom[0]=Momentum(0);
			mom[1]=Momentum(1);
			mom[2]=Momentum(2);
			mom[3]=Momentum(3);
			// Z-position for hit
			
			
			/**************** Photons lost in second grid have to be calculated by hand************/ 
			
			Float_t cophi = TMath::Cos(TMath::ATan2(mom[0], mom[1]));
			Float_t t = (1. - .025 / cophi) * (1. - .05 /  cophi);
			gMC->Rndm(ranf, 1);
			//printf("grid calculation:%f\n",t);
			if (ranf[0] > t) {
			  geant3->StopTrack();
			  Ckov_data[13] = 5;
			  AddCerenkov(gAlice->CurrentTrack(),vol,Ckov_data);
			  //printf("Lost one in grid\n");
			}
			/**********************************************************************************/
		      }    //gap
		    
		    if (gMC->VolId("CSI ")==gMC->CurrentVolID(copy))             //is it in csi?      
		      {
			gMC->TrackMomentum(Momentum);
			mom[0]=Momentum(0);
			mom[1]=Momentum(1);
			mom[2]=Momentum(2);
			mom[3]=Momentum(3);
			
			/********* Photons lost by Fresnel reflection have to be calculated by hand********/ 
			/***********************Cerenkov phtons (always polarised)*************************/
			
			Float_t cophi = TMath::Cos(TMath::ATan2(mom[0], mom[1]));
			Float_t t = Fresnel(Ckov_energy*1e9,cophi,1);
			gMC->Rndm(ranf, 1);
			if (ranf[0] < t) {
			  geant3->StopTrack();
			  Ckov_data[13] = 6;
			  AddCerenkov(gAlice->CurrentTrack(),vol,Ckov_data);
			  //printf("Lost by Fresnel\n");
			}
			/**********************************************************************************/
		      }
		  } //track entering?
		  
		  
		  /********************Evaluation of losses************************/
		  /******************still in the old fashion**********************/
		  
		  Int_t i1 = geant3->Gctrak()->nmec;            //number of physics mechanisms acting on the particle
		  for (Int_t i = 0; i < i1; ++i) {
		    //        Reflection loss 
		    if (geant3->Gctrak()->lmec[i] == 106) {        //was it reflected
		      Ckov_data[13]=10;
		      if (gMC->VolId("FRE1")==gMC->CurrentVolID(copy) || gMC->VolId("FRE2")==gMC->CurrentVolID(copy)) 
			Ckov_data[13]=1;
		      if (gMC->CurrentVolID(copy) == gMC->VolId("QUAR")) 
			Ckov_data[13]=2;
		      geant3->StopTrack();
		      AddCerenkov(gAlice->CurrentTrack(),vol,Ckov_data);
		    } //reflection question
		    
		    
		    //        Absorption loss 
		    else if (geant3->Gctrak()->lmec[i] == 101) {              //was it absorbed?
		      Ckov_data[13]=20;
		      if (gMC->VolId("FRE1")==gMC->CurrentVolID(copy) || gMC->VolId("FRE2")==gMC->CurrentVolID(copy)) 
			Ckov_data[13]=11;
		      if (gMC->CurrentVolID(copy) == gMC->VolId("QUAR")) 
			Ckov_data[13]=12;
		      if (gMC->CurrentVolID(copy) == gMC->VolId("META")) 
			Ckov_data[13]=13;
		      if (gMC->CurrentVolID(copy) == gMC->VolId("GAP ")) 
			Ckov_data[13]=13;
		      
		      if (gMC->CurrentVolID(copy) == gMC->VolId("SRIC")) 
			Ckov_data[13]=15;
		      
		      //        CsI inefficiency 
		      if (gMC->CurrentVolID(copy) == gMC->VolId("CSI ")) {
			Ckov_data[13]=16;
		      }
		      geant3->StopTrack();
		      AddCerenkov(gAlice->CurrentTrack(),vol,Ckov_data);
		      //printf("Added cerenkov %d\n",fCkov_number);
		    } //absorption question 
		    
		    
		    //        Photon goes out of tracking scope 
		    else if (geant3->Gctrak()->lmec[i] == 30) {                 //is it below energy treshold?
		      Ckov_data[13]=21;
		      geant3->StopTrack();
		      AddCerenkov(gAlice->CurrentTrack(),vol,Ckov_data);
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
	if (gMC->VolId("CSI ")==gMC->CurrentVolID(copy))
	{
	    
	  if (gMC->Edep() > 0.){
		gMC->TrackPosition(Position);
		gMC->TrackMomentum(Momentum);
		pos[0]=Position(0);
		pos[1]=Position(1);
		pos[2]=Position(2);
		mom[0]=Momentum(0);
		mom[1]=Momentum(1);
		mom[2]=Momentum(2);
		mom[3]=Momentum(3);
		Double_t tc = mom[0]*mom[0]+mom[1]*mom[1];
		Double_t rt = TMath::Sqrt(tc);
		theta   = Float_t(TMath::ATan2(rt,Double_t(mom[2])))*kRaddeg;
		phi     = Float_t(TMath::ATan2(Double_t(mom[1]),Double_t(mom[0])))*kRaddeg;
		gMC->Gmtod(pos,Localpos,1);                                                                    
		gMC->Gmtod(mom,Localmom,2);
		
		gMC->CurrentVolOffID(2,copy);
		vol[0]=copy;
		idvol=vol[0]-1;

		//Int_t sector=((AliRICHChamber*) (*fChambers)[idvol])
			//->Sector(Localpos[0], Localpos[2]);
		//printf("Sector:%d\n",sector);

		/*if (gMC->TrackPid() == 50000051){
		  fFeedbacks++;
		  printf("Feedbacks:%d\n",fFeedbacks);
		}*/	
		
		((AliRICHChamber*) (*fChambers)[idvol])
		    ->SigGenInit(Localpos[0], Localpos[2], Localpos[1]);
		if(idvol<7) {	
		    Ckov_data[0] = gMC->TrackPid();        // particle type
		    Ckov_data[1] = pos[0];                 // X-position for hit
		    Ckov_data[2] = pos[1];                 // Y-position for hit
		    Ckov_data[3] = pos[2];                 // Z-position for hit
		    Ckov_data[4] = theta;                      // theta angle of incidence
		    Ckov_data[5] = phi;                      // phi angle of incidence 
		    Ckov_data[8] = (Float_t) fNPadHits;      // first padhit
		    Ckov_data[9] = -1;                       // last pad hit
		    Ckov_data[13] = 4;                       // photon was detected
		    Ckov_data[14] = mom[0];
		    Ckov_data[15] = mom[1];
		    Ckov_data[16] = mom[2];
		    
		    destep = gMC->Edep();
		    gMC->SetMaxStep(big);
		    cherenkov_loss  += destep;
		    Ckov_data[7]=cherenkov_loss;
		    
		    NPads = MakePadHits(Localpos[0],Localpos[2],cherenkov_loss,idvol,cerenkov);
		    if (fNPadHits > (Int_t)Ckov_data[8]) {
			Ckov_data[8]= Ckov_data[8]+1;
			Ckov_data[9]= (Float_t) fNPadHits;
		    }

		    Ckov_data[17] = NPads;
		    //printf("Npads:%d",NPads);
		    
		    //TClonesArray *Hits = RICH->Hits();
		    AliRICHHit *mipHit =  (AliRICHHit*) (fHits->UncheckedAt(0));
		    if (mipHit)
		      {
			mom[0] = current->Px();
			mom[1] = current->Py();
			mom[2] = current->Pz();
			Float_t Mip_px = mipHit->fMomX;
			Float_t Mip_py = mipHit->fMomY;
			Float_t Mip_pz = mipHit->fMomZ;
			
			Float_t r = mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2];
			Float_t rt = TMath::Sqrt(r);
			Float_t Mip_r = Mip_px*Mip_px + Mip_py*Mip_py + Mip_pz*Mip_pz;	
			Float_t Mip_rt = TMath::Sqrt(Mip_r);
			if ((rt*Mip_rt) > 0)
			  {
			    coscerenkov = (mom[0]*Mip_px + mom[1]*Mip_py + mom[2]*Mip_pz)/(rt*Mip_rt);
			  }
			else
			  {
			    coscerenkov = 0;
			  }
			Float_t cherenkov = TMath::ACos(coscerenkov);
			Ckov_data[18]=cherenkov;
		      }
		    //if (sector != -1)
		    //{
		    AddHit(gAlice->CurrentTrack(),vol,Ckov_data);
		    AddCerenkov(gAlice->CurrentTrack(),vol,Ckov_data);
		    //}
		}
	    }
	}
    }
    
    /***********************************************End of photon hits*********************************************/
    

    /**********************************************Charged particles treatment*************************************/

    else if (gMC->TrackCharge())
    //else if (1 == 1)
      {
//If MIP
	/*if (gMC->IsTrackEntering())
	  {                
	    hits[13]=20;//is track entering?
	  }*/
	if (gMC->VolId("FRE1")==gMC->CurrentVolID(copy) || gMC->VolId("FRE2")==gMC->CurrentVolID(copy))
	  {
	    fFreon_prod=1;
	  }

	if (gMC->VolId("GAP ")== gMC->CurrentVolID(copy)) {
// Get current particle id (ipart), track position (pos)  and momentum (mom)
	    
	    gMC->CurrentVolOffID(3,copy);
	    vol[0]=copy;
	    idvol=vol[0]-1;

	    //Int_t sector=((AliRICHChamber*) (*fChambers)[idvol])
			//->Sector(Localpos[0], Localpos[2]);
	    //printf("Sector:%d\n",sector);
	    
	    gMC->TrackPosition(Position);
	    gMC->TrackMomentum(Momentum);
	    pos[0]=Position(0);
	    pos[1]=Position(1);
	    pos[2]=Position(2);
	    mom[0]=Momentum(0);
	    mom[1]=Momentum(1);
	    mom[2]=Momentum(2);
	    mom[3]=Momentum(3);
	    gMC->Gmtod(pos,Localpos,1);                                                                    
	    gMC->Gmtod(mom,Localmom,2);
	    
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
		

		Double_t Localtc = Localmom[0]*Localmom[0]+Localmom[2]*Localmom[2];
		Double_t Localrt = TMath::Sqrt(Localtc);
		Localtheta   = Float_t(TMath::ATan2(Localrt,Double_t(Localmom[1])))*kRaddeg;                       
		Localphi     = Float_t(TMath::ATan2(Double_t(Localmom[2]),Double_t(Localmom[0])))*kRaddeg;    
		
		hits[0] = Float_t(ipart);         // particle type
		hits[1] = Localpos[0];                 // X-position for hit
		hits[2] = Localpos[1];                 // Y-position for hit
		hits[3] = Localpos[2];                 // Z-position for hit
		hits[4] = Localtheta;                  // theta angle of incidence
		hits[5] = Localphi;                    // phi angle of incidence 
		hits[8] = (Float_t) fNPadHits;    // first padhit
		hits[9] = -1;                     // last pad hit
		hits[13] = fFreon_prod;           // did id hit the freon?
		hits[14] = mom[0];
		hits[15] = mom[1];
		hits[16] = mom[2];

		tlength = 0;
		eloss   = 0;
		fFreon_prod = 0;
	
		Chamber(idvol).LocaltoGlobal(Localpos,hits+1);
	   
		
		//To make chamber coordinates x-y had to pass LocalPos[0], LocalPos[2]
		xhit    = Localpos[0];
		yhit    = Localpos[2];
		// Only if not trigger chamber
		if(idvol<7) {
		    //
		    //  Initialize hit position (cursor) in the segmentation model 
		    ((AliRICHChamber*) (*fChambers)[idvol])
			->SigGenInit(Localpos[0], Localpos[2], Localpos[1]);
		}
	    }
	    
	    // 
	    // Calculate the charge induced on a pad (disintegration) in case 
	    //
	    // Mip left chamber ...
	    if( gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared()){
		gMC->SetMaxStep(big);
		eloss   += destep;
		tlength += step;
		
				
		// Only if not trigger chamber
		if(idvol<7) {
		  if (eloss > 0) 
		    {
		      if(gMC->TrackPid() == kNeutron)
			printf("\n\n\n\n\n Neutron Making Pad Hit!!! \n\n\n\n");
		      NPads = MakePadHits(xhit,yhit,eloss,idvol,mip);
		      hits[17] = NPads;
		      //printf("Npads:%d",NPads);
		    }
		}
		
		hits[6]=tlength;
		hits[7]=eloss;
		if (fNPadHits > (Int_t)hits[8]) {
		    hits[8]= hits[8]+1;
		    hits[9]= (Float_t) fNPadHits;
		}
		
		//if(sector !=-1)
		new(lhits[fNhits++]) AliRICHHit(fIshunt,gAlice->CurrentTrack(),vol,hits);
		eloss = 0; 
		//
		// Check additional signal generation conditions 
		// defined by the segmentation
		// model (boundary crossing conditions) 
	    } else if 
		(((AliRICHChamber*) (*fChambers)[idvol])
		 ->SigGenCond(Localpos[0], Localpos[2], Localpos[1]))
	    {
		((AliRICHChamber*) (*fChambers)[idvol])
		    ->SigGenInit(Localpos[0], Localpos[2], Localpos[1]);
		if (eloss > 0) 
		  {
		    if(gMC->TrackPid() == kNeutron)
		      printf("\n\n\n\n\n Neutron Making Pad Hit!!! \n\n\n\n");
		    NPads = MakePadHits(xhit,yhit,eloss,idvol,mip);
		    hits[17] = NPads;
		    //printf("Npads:%d",NPads);
		  }
		xhit     = Localpos[0];
		yhit     = Localpos[2]; 
		eloss    = destep;
		tlength += step ;
		//
		// nothing special  happened, add up energy loss
	    } else {        
		eloss   += destep;
		tlength += step ;
	    }
	}
      }
    /*************************************************End of MIP treatment**************************************/
   //}
}

  
//___________________________________________
Int_t AliRICH::MakePadHits(Float_t xhit,Float_t yhit,Float_t eloss, Int_t idvol, Response_t res)
{
//
//  Calls the charge disintegration method of the current chamber and adds
//  the simulated cluster to the root treee 
//
    Int_t clhits[7];
    Float_t newclust[6][500];
    Int_t nnew;
    
//
//  Integrated pulse height on chamber
    
    clhits[0]=fNhits+1;
    
    ((AliRICHChamber*) (*fChambers)[idvol])->DisIntegration(eloss, xhit, yhit, nnew, newclust, res);
    Int_t ic=0;
    
//
//  Add new clusters
    for (Int_t i=0; i<nnew; i++) {
	if (Int_t(newclust[3][i]) > 0) {
	    ic++;
// Cathode plane
	    clhits[1] = Int_t(newclust[5][i]);
//  Cluster Charge
	    clhits[2] = Int_t(newclust[0][i]);
//  Pad: ix
	    clhits[3] = Int_t(newclust[1][i]);
//  Pad: iy 
	    clhits[4] = Int_t(newclust[2][i]);
//  Pad: charge
	    clhits[5] = Int_t(newclust[3][i]);
//  Pad: chamber sector
	    clhits[6] = Int_t(newclust[4][i]);
	    
	    AddPadHit(clhits);
	}
    }
return nnew;
}
