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
  Revision 1.21  2000/07/21 10:21:07  morsch
  fNrawch   = 0; and  fNrechits = 0; in the default constructor.

  Revision 1.20  2000/07/10 15:28:39  fca
  Correction of the inheritance scheme

  Revision 1.19  2000/06/30 16:29:51  dibari
  Added kDebugLevel variable to control output size on demand

  Revision 1.18  2000/06/12 15:15:46  jbarbosa
  Cleaned up version.

  Revision 1.17  2000/06/09 14:58:37  jbarbosa
  New digitisation per particle type

  Revision 1.16  2000/04/19 12:55:43  morsch
  Newly structured and updated version (JB, AM)

*/


////////////////////////////////////////////////
//  Manager and hits classes for set:RICH     //
////////////////////////////////////////////////

#include <TBRIK.h>
#include <TTUBE.h>
#include <TNode.h> 
#include <TRandom.h> 
#include <TObject.h>
#include <TVector.h>
#include <TObjArray.h>
#include <TArrayF.h>
#include <TFile.h>
#include <TParticle.h>
#include <iostream.h>

#include "AliRICH.h"
#include "AliRICHSegmentation.h"
#include "AliRICHHit.h"
#include "AliRICHCerenkov.h"
#include "AliRICHPadHit.h"
#include "AliRICHDigit.h"
#include "AliRICHTransientDigit.h"
#include "AliRICHRawCluster.h"
#include "AliRICHRecHit.h"
#include "AliRICHHitMapA1.h"
#include "AliRICHClusterFinder.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliConst.h"
#include "AliPDG.h"
#include "AliPoints.h"
#include "AliCallf77.h" 
#include "TGeant3.h"


// Static variables for the pad-hit iterator routines
static Int_t sMaxIterPad=0;
static Int_t sCurIterPad=0;
static TClonesArray *fClusters2;
static TClonesArray *fHits2;
static TTree *TrH1;
 
ClassImp(AliRICH)
    
//___________________________________________
AliRICH::AliRICH()
{
// Default constructor for RICH manager class

    fIshunt     = 0;
    fHits       = 0;
    fPadHits    = 0;
    fNPadHits   = 0;
    fNcerenkovs = 0;
    fDchambers  = 0;
    fCerenkovs  = 0;
    fNdch       = 0;
    fNrawch   = 0;
    fNrechits = 0;
}

//___________________________________________
AliRICH::AliRICH(const char *name, const char *title)
    : AliDetector(name,title)
{
//Begin_Html
/*
  <img src="gif/alirich.gif">
*/
//End_Html
    
    fHits       = new TClonesArray("AliRICHHit",1000  );
    gAlice->AddHitList(fHits);
    fPadHits    = new TClonesArray("AliRICHPadHit",100000);
    fCerenkovs  = new TClonesArray("AliRICHCerenkov",1000);
    gAlice->AddHitList(fCerenkovs);
    //gAlice->AddHitList(fHits);
    fNPadHits   = 0;
    fNcerenkovs = 0;
    fIshunt     = 0;
    
    fNdch      = new Int_t[kNCH];
    
    fDchambers = new TObjArray(kNCH);

    fRecHits = new TObjArray(kNCH);
    
    Int_t i;
   
    for (i=0; i<kNCH ;i++) {
	(*fDchambers)[i] = new TClonesArray("AliRICHDigit",10000); 
	fNdch[i]=0;
    }

    fNrawch      = new Int_t[kNCH];
    
    fRawClusters = new TObjArray(kNCH);
    //printf("Created fRwClusters with adress:%p",fRawClusters);

    for (i=0; i<kNCH ;i++) {
      (*fRawClusters)[i] = new TClonesArray("AliRICHRawCluster",10000); 
      fNrawch[i]=0;
    }

    fNrechits      = new Int_t[kNCH];
    
    for (i=0; i<kNCH ;i++) {
	(*fRecHits)[i] = new TClonesArray("AliRICHRecHit",1000); 
    }
    //printf("Created fRecHits with adress:%p",fRecHits);

        
    SetMarkerColor(kRed);
}

AliRICH::AliRICH(const AliRICH& RICH)
{
// Copy Constructor
}


//___________________________________________
AliRICH::~AliRICH()
{

// Destructor of RICH manager class

    fIshunt  = 0;
    delete fHits;
    delete fPadHits;
    delete fCerenkovs;
}

//___________________________________________
void AliRICH::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{

//  
// Adds a hit to the Hits list
//

    TClonesArray &lhits = *fHits;
    new(lhits[fNhits++]) AliRICHHit(fIshunt,track,vol,hits);
}
//_____________________________________________________________________________
void AliRICH::AddCerenkov(Int_t track, Int_t *vol, Float_t *cerenkovs)
{

//
// Adds a RICH cerenkov hit to the Cerenkov Hits list
//

    TClonesArray &lcerenkovs = *fCerenkovs;
    new(lcerenkovs[fNcerenkovs++]) AliRICHCerenkov(fIshunt,track,vol,cerenkovs);
    //printf ("Done for Cerenkov %d\n\n\n\n",fNcerenkovs);
}
//___________________________________________
void AliRICH::AddPadHit(Int_t *clhits)
{

//
// Add a RICH pad hit to the list
//

    TClonesArray &lPadHits = *fPadHits;
    new(lPadHits[fNPadHits++]) AliRICHPadHit(clhits);
} 
//_____________________________________________________________________________
void AliRICH::AddDigits(Int_t id, Int_t *tracks, Int_t *charges, Int_t *digits)
{

  //
  // Add a RICH digit to the list
  //

    TClonesArray &ldigits = *((TClonesArray*)(*fDchambers)[id]);
    new(ldigits[fNdch[id]++]) AliRICHDigit(tracks,charges,digits);
}

//_____________________________________________________________________________
void AliRICH::AddRawCluster(Int_t id, const AliRICHRawCluster& c)
{
    //
    // Add a RICH digit to the list
    //

    TClonesArray &lrawcl = *((TClonesArray*)(*fRawClusters)[id]);
    new(lrawcl[fNrawch[id]++]) AliRICHRawCluster(c);
}

//_____________________________________________________________________________
void AliRICH::AddRecHit(Int_t id, Float_t *rechit, Float_t *photons, Int_t *padsx, Int_t* padsy)
{
  
  //
  // Add a RICH reconstructed hit to the list
  //

    TClonesArray &lrec = *((TClonesArray*)(*fRecHits)[id]);
    new(lrec[fNrechits[id]++]) AliRICHRecHit(id,rechit,photons,padsx,padsy);
}

//___________________________________________
void AliRICH::BuildGeometry()
    
{
  
  //
  // Builds a TNode geometry for event display
  //
    TNode *node, *top;
    
    const int kColorRICH = kGreen;
    //
    top=gAlice->GetGeometry()->GetNode("alice");
    
    
    new TBRIK("S_RICH","S_RICH","void",71.09999,11.5,73.15);
    
    top->cd();
    Float_t pos1[3]={0,471.8999,165.2599};
    //Chamber(0).SetChamberTransform(pos1[0],pos1[1],pos1[2],
    new TRotMatrix("rot993","rot993",90,0,70.69,90,19.30999,-90);
    node = new TNode("RICH1","RICH1","S_RICH",pos1[0],pos1[1],pos1[2],"rot993");
    

    node->SetLineColor(kColorRICH);
    fNodes->Add(node);
    top->cd();
    
    Float_t pos2[3]={171,470,0};
    //Chamber(1).SetChamberTransform(pos2[0],pos2[1],pos2[2],
    new TRotMatrix("rot994","rot994",90,-20,90,70,0,0);
    node = new TNode("RICH2","RICH2","S_RICH",pos2[0],pos2[1],pos2[2],"rot994");
    
    
    node->SetLineColor(kColorRICH);
    fNodes->Add(node);
    top->cd();
    Float_t pos3[3]={0,500,0};
    //Chamber(2).SetChamberTransform(pos3[0],pos3[1],pos3[2],
    new TRotMatrix("rot995","rot995",90,0,90,90,0,0);
    node = new TNode("RICH3","RICH3","S_RICH",pos3[0],pos3[1],pos3[2],"rot995");
    

    node->SetLineColor(kColorRICH);
    fNodes->Add(node);
    top->cd();
    Float_t pos4[3]={-171,470,0};
    //Chamber(3).SetChamberTransform(pos4[0],pos4[1],pos4[2], 
    new TRotMatrix("rot996","rot996",90,20,90,110,0,0);  
    node = new TNode("RICH4","RICH4","S_RICH",pos4[0],pos4[1],pos4[2],"rot996");
    

    node->SetLineColor(kColorRICH);
    fNodes->Add(node);
    top->cd();
    Float_t pos5[3]={161.3999,443.3999,-165.3};
    //Chamber(4).SetChamberTransform(pos5[0],pos5[1],pos5[2],
    new TRotMatrix("rot997","rot997",90,340,108.1999,70,18.2,70);
    node = new TNode("RICH5","RICH5","S_RICH",pos5[0],pos5[1],pos5[2],"rot997");
    
    node->SetLineColor(kColorRICH);
    fNodes->Add(node);
    top->cd();
    Float_t pos6[3]={0., 471.9, -165.3,};
    //Chamber(5).SetChamberTransform(pos6[0],pos6[1],pos6[2],
    new TRotMatrix("rot998","rot998",90,0,109.3099,90,19.30999,90);
    node = new TNode("RICH6","RICH6","S_RICH",pos6[0],pos6[1],pos6[2],"rot998");
    
    
    node->SetLineColor(kColorRICH);
    fNodes->Add(node);
    top->cd();
    Float_t pos7[3]={-161.399,443.3999,-165.3};
    //Chamber(6).SetChamberTransform(pos7[0],pos7[1],pos7[2],
    new TRotMatrix("rot999","rot999",90,20,108.1999,110,18.2,110);
    node = new TNode("RICH7","RICH7","S_RICH",pos7[0],pos7[1],pos7[2],"rot999");
    node->SetLineColor(kColorRICH);
    fNodes->Add(node); 
    
}

//___________________________________________
void AliRICH::CreateGeometry()
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
void AliRICH::CreateMaterials()
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

Float_t AliRICH::Fresnel(Float_t ene,Float_t pdoti, Bool_t pola)
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
Float_t AliRICH::AbsoCH4(Float_t x)
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
Int_t AliRICH::DistancetoPrimitive(Int_t , Int_t )
{

// Default value

    return 9999;
}

//___________________________________________
void AliRICH::MakeBranch(Option_t* option)
{
  // Create Tree branches for the RICH.
    
    const Int_t kBufferSize = 4000;
    char branchname[20];
    
    
    AliDetector::MakeBranch(option);
    sprintf(branchname,"%sCerenkov",GetName());
    if (fCerenkovs   && gAlice->TreeH()) {
	gAlice->TreeH()->Branch(branchname,&fCerenkovs, kBufferSize);
	printf("Making Branch %s for Cerenkov Hits\n",branchname);
    }
    
    sprintf(branchname,"%sPadHits",GetName());
    if (fPadHits   && gAlice->TreeH()) {
	gAlice->TreeH()->Branch(branchname,&fPadHits, kBufferSize);
	printf("Making Branch %s for PadHits\n",branchname);
    }
    
// one branch for digits per chamber
    Int_t i;
    
    for (i=0; i<kNCH ;i++) {
	sprintf(branchname,"%sDigits%d",GetName(),i+1);
	
	if (fDchambers   && gAlice->TreeD()) {
	    gAlice->TreeD()->Branch(branchname,&((*fDchambers)[i]), kBufferSize);
	    printf("Making Branch %s for digits in chamber %d\n",branchname,i+1);
	}	
    }

// one branch for raw clusters per chamber
  for (i=0; i<kNCH ;i++) {
      sprintf(branchname,"%sRawClusters%d",GetName(),i+1);
      
      if (fRawClusters   && gAlice->TreeR()) {
	 gAlice->TreeR()->Branch(branchname,&((*fRawClusters)[i]), kBufferSize);
	 printf("Making Branch %s for raw clusters in chamber %d\n",branchname,i+1);
      }	
  }

  // one branch for rec hits per chamber
  for (i=0; i<kNCH ;i++) {
    sprintf(branchname,"%sRecHits%d",GetName(),i+1);
    
    if (fRecHits   && gAlice->TreeR()) {
      gAlice->TreeR()->Branch(branchname,&((*fRecHits)[i]), kBufferSize);
      printf("Making Branch %s for rec. hits in chamber %d\n",branchname,i+1);
    }	
  }
}

//___________________________________________
void AliRICH::SetTreeAddress()
{
  // Set branch address for the Hits and Digits Tree.
  char branchname[20];
  Int_t i;

    AliDetector::SetTreeAddress();
    
    TBranch *branch;
    TTree *treeH = gAlice->TreeH();
    TTree *treeD = gAlice->TreeD();
    TTree *treeR = gAlice->TreeR();
    
    if (treeH) {
	if (fPadHits) {
	    branch = treeH->GetBranch("RICHPadHits");
	    if (branch) branch->SetAddress(&fPadHits);
	}
	if (fCerenkovs) {
	    branch = treeH->GetBranch("RICHCerenkov");
	    if (branch) branch->SetAddress(&fCerenkovs);
	}
    }
    
    if (treeD) {
	for (int i=0; i<kNCH; i++) {
	    sprintf(branchname,"%sDigits%d",GetName(),i+1);
	    if (fDchambers) {
		branch = treeD->GetBranch(branchname);
		if (branch) branch->SetAddress(&((*fDchambers)[i]));
	    }
	}
    }
  if (treeR) {
      for (i=0; i<kNCH; i++) {
	  sprintf(branchname,"%sRawClusters%d",GetName(),i+1);
	  if (fRawClusters) {
	      branch = treeR->GetBranch(branchname);
	      if (branch) branch->SetAddress(&((*fRawClusters)[i]));
	  }
      }
      
      for (i=0; i<kNCH; i++) {
	sprintf(branchname,"%sRecHits%d",GetName(),i+1);
	if (fRecHits) {
	  branch = treeR->GetBranch(branchname);
	  if (branch) branch->SetAddress(&((*fRecHits)[i]));
	  }
      }
      
  }
}
//___________________________________________
void AliRICH::ResetHits()
{
  // Reset number of clusters and the cluster array for this detector
    AliDetector::ResetHits();
    fNPadHits   = 0;
    fNcerenkovs = 0;
    if (fPadHits)  fPadHits->Clear();
    if (fCerenkovs) fCerenkovs->Clear();
}


//____________________________________________
void AliRICH::ResetDigits()
{
  //
  // Reset number of digits and the digits array for this detector
  //
    for ( int i=0;i<kNCH;i++ ) {
	if ((*fDchambers)[i])   (*fDchambers)[i]->Clear();
	if (fNdch)  fNdch[i]=0;
    }
}

//____________________________________________
void AliRICH::ResetRawClusters()
{
  //
  // Reset number of raw clusters and the raw clust array for this detector
  //
    for ( int i=0;i<kNCH;i++ ) {
	if ((*fRawClusters)[i])    ((TClonesArray*)(*fRawClusters)[i])->Clear();
	if (fNrawch)  fNrawch[i]=0;
    }
}

//____________________________________________
void AliRICH::ResetRecHits()
{
  //
  // Reset number of raw clusters and the raw clust array for this detector
  //
  
  for ( int i=0;i<kNCH;i++ ) {
	if ((*fRecHits)[i])    ((TClonesArray*)(*fRecHits)[i])->Clear();
	if (fNrechits)  fNrechits[i]=0;
    }
}

//___________________________________________
void   AliRICH::SetGeometryModel(Int_t id, AliRICHGeometry *geometry)
{

//
// Setter for the RICH geometry model
//


    ((AliRICHChamber*) (*fChambers)[id])->GeometryModel(geometry);
}

//___________________________________________
void   AliRICH::SetSegmentationModel(Int_t id, AliRICHSegmentation *segmentation)
{

//
// Setter for the RICH segmentation model
//

    ((AliRICHChamber*) (*fChambers)[id])->SegmentationModel(segmentation);
}

//___________________________________________
void   AliRICH::SetResponseModel(Int_t id, AliRICHResponse *response)
{

//
// Setter for the RICH response model
//

    ((AliRICHChamber*) (*fChambers)[id])->ResponseModel(response);
}

void   AliRICH::SetReconstructionModel(Int_t id, AliRICHClusterFinder *reconst)
{

//
// Setter for the RICH reconstruction model (clusters)
//

    ((AliRICHChamber*) (*fChambers)[id])->ReconstructionModel(reconst);
}

void   AliRICH::SetNsec(Int_t id, Int_t nsec)
{

//
// Sets the number of padplanes
//

    ((AliRICHChamber*) (*fChambers)[id])->SetNsec(nsec);
}


//___________________________________________
void AliRICH::StepManager()
{

// Full Step Manager

    Int_t          copy, id;
    static Int_t   idvol;
    static Int_t   vol[2];
    Int_t          ipart;
    static Float_t hits[18];
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
    Float_t        ranf[2];
    Int_t          nPads;
    Float_t        coscerenkov;
    static Float_t eloss, xhit, yhit, tlength;
    const  Float_t kBig=1.e10;
       
    TClonesArray &lhits = *fHits;
    TGeant3 *geant3 = (TGeant3*) gMC;
    TParticle *current = (TParticle*)(*gAlice->Particles())[gAlice->CurrentTrack()];

 //if (current->Energy()>1)
   //{
        
    // Only gas gap inside chamber
    // Tag chambers and record hits when track enters 
    
    idvol=-1;
    id=gMC->CurrentVolID(copy);
    Float_t cherenkovLoss=0;
    //gAlice->KeepTrack(gAlice->CurrentTrack());
    
    gMC->TrackPosition(position);
    pos[0]=position(0);
    pos[1]=position(1);
    pos[2]=position(2);
    bzero(ckovData,sizeof(ckovData)*19);
    ckovData[1] = pos[0];                 // X-position for hit
    ckovData[2] = pos[1];                 // Y-position for hit
    ckovData[3] = pos[2];                 // Z-position for hit
    //ckovData[11] = gAlice->CurrentTrack();

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
	      if (gMC->IsTrackEntering()){                                     //is track entering?
		if (gMC->VolId("FRE1")==gMC->CurrentVolID(copy) || gMC->VolId("FRE2")==gMC->CurrentVolID(copy))
		  {                                                          //is it in freo?
		    if (geant3->Gctrak()->nstep<1){                          //is it the first step?
		      Int_t mother = current->GetFirstMother(); 
		      
		      //printf("Second Mother:%d\n",current->GetSecondMother());
		      
		      ckovData[10] = mother;
		      ckovData[11] = gAlice->CurrentTrack();
		      ckovData[12] = 1;             //Media where photon was produced 1->Freon, 2->Quarz
		      fCkovNumber++;
		      fFreonProd=1;
		      //printf("Index: %d\n",fCkovNumber);
		    }    //first step question
		  }        //freo question
		
		if (geant3->Gctrak()->nstep<1){                                  //is it first step?
		  if (gMC->VolId("QUAR")==gMC->CurrentVolID(copy))             //is it in quarz?
		    {
		      ckovData[12] = 2;
		    }    //quarz question
		}        //first step question
		
		//printf("Before %d\n",fFreonProd);
	      }   //track entering question
	      
	      if (ckovData[12] == 1)                                        //was it produced in Freon?
		//if (fFreonProd == 1)
		{
		  if (gMC->IsTrackEntering()){                                     //is track entering?
		    //printf("Got in");
		    if (gMC->VolId("META")==gMC->CurrentVolID(copy))                //is it in gap?      
		      {
			//printf("Got in\n");
			gMC->TrackMomentum(momentum);
			mom[0]=momentum(0);
			mom[1]=momentum(1);
			mom[2]=momentum(2);
			mom[3]=momentum(3);
			// Z-position for hit
			
			
			/**************** Photons lost in second grid have to be calculated by hand************/ 
			
			Float_t cophi = TMath::Cos(TMath::ATan2(mom[0], mom[1]));
			Float_t t = (1. - .025 / cophi) * (1. - .05 /  cophi);
			gMC->Rndm(ranf, 1);
			//printf("grid calculation:%f\n",t);
			if (ranf[0] > t) {
			  //geant3->StopTrack();
			  ckovData[13] = 5;
			  AddCerenkov(gAlice->CurrentTrack(),vol,ckovData);
			  //printf("Lost one in grid\n");
			}
			/**********************************************************************************/
		      }    //gap
		    
		    if (gMC->VolId("CSI ")==gMC->CurrentVolID(copy))             //is it in csi?      
		      {
			gMC->TrackMomentum(momentum);
			mom[0]=momentum(0);
			mom[1]=momentum(1);
			mom[2]=momentum(2);
			mom[3]=momentum(3);
			
			/********* Photons lost by Fresnel reflection have to be calculated by hand********/ 
			/***********************Cerenkov phtons (always polarised)*************************/
			
			Float_t cophi = TMath::Cos(TMath::ATan2(mom[0], mom[1]));
			Float_t t = Fresnel(ckovEnergy*1e9,cophi,1);
			gMC->Rndm(ranf, 1);
			if (ranf[0] < t) {
			  //geant3->StopTrack();
			  ckovData[13] = 6;
			  AddCerenkov(gAlice->CurrentTrack(),vol,ckovData);
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
		      ckovData[13]=10;
		      if (gMC->VolId("FRE1")==gMC->CurrentVolID(copy) || gMC->VolId("FRE2")==gMC->CurrentVolID(copy)) 
			ckovData[13]=1;
		      if (gMC->CurrentVolID(copy) == gMC->VolId("QUAR")) 
			ckovData[13]=2;
		      //geant3->StopTrack();
		      AddCerenkov(gAlice->CurrentTrack(),vol,ckovData);
		    } //reflection question
		    
		    
		    //        Absorption loss 
		    else if (geant3->Gctrak()->lmec[i] == 101) {              //was it absorbed?
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
		      //geant3->StopTrack();
		      AddCerenkov(gAlice->CurrentTrack(),vol,ckovData);
		      //printf("Added cerenkov %d\n",fCkovNumber);
		    } //absorption question 
		    
		    
		    //        Photon goes out of tracking scope 
		    else if (geant3->Gctrak()->lmec[i] == 30) {                 //is it below energy treshold?
		      ckovData[13]=21;
		      //geant3->StopTrack();
		      AddCerenkov(gAlice->CurrentTrack(),vol,ckovData);
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
		gMC->Gmtod(pos,localPos,1);                                                                    
		gMC->Gmtod(mom,localMom,2);
		
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
		
		((AliRICHChamber*) (*fChambers)[idvol])
		    ->SigGenInit(localPos[0], localPos[2], localPos[1]);
		if(idvol<kNCH) {	
		    ckovData[0] = gMC->TrackPid();        // particle type
		    ckovData[1] = pos[0];                 // X-position for hit
		    ckovData[2] = pos[1];                 // Y-position for hit
		    ckovData[3] = pos[2];                 // Z-position for hit
		    ckovData[4] = theta;                      // theta angle of incidence
		    ckovData[5] = phi;                      // phi angle of incidence 
		    ckovData[8] = (Float_t) fNPadHits;      // first padhit
		    ckovData[9] = -1;                       // last pad hit
		    ckovData[13] = 4;                       // photon was detected
		    ckovData[14] = mom[0];
		    ckovData[15] = mom[1];
		    ckovData[16] = mom[2];
		    
		    destep = gMC->Edep();
		    gMC->SetMaxStep(kBig);
		    cherenkovLoss  += destep;
		    ckovData[7]=cherenkovLoss;
		    
		    nPads = MakePadHits(localPos[0],localPos[2],cherenkovLoss,idvol,kCerenkov);
		    if (fNPadHits > (Int_t)ckovData[8]) {
			ckovData[8]= ckovData[8]+1;
			ckovData[9]= (Float_t) fNPadHits;
		    }

		    ckovData[17] = nPads;
		    //printf("nPads:%d",nPads);
		    
		    //TClonesArray *Hits = RICH->Hits();
		    AliRICHHit *mipHit =  (AliRICHHit*) (fHits->UncheckedAt(0));
		    if (mipHit)
		      {
			mom[0] = current->Px();
			mom[1] = current->Py();
			mom[2] = current->Pz();
			Float_t mipPx = mipHit->fMomX;
			Float_t mipPy = mipHit->fMomY;
			Float_t mipPz = mipHit->fMomZ;
			
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
		    AddHit(gAlice->CurrentTrack(),vol,ckovData);
		    AddCerenkov(gAlice->CurrentTrack(),vol,ckovData);
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
	    fFreonProd=1;
	  }

	if (gMC->VolId("GAP ")== gMC->CurrentVolID(copy)) {
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
	    gMC->Gmtod(mom,localMom,2);
	    
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
		hits[8] = (Float_t) fNPadHits;    // first padhit
		hits[9] = -1;                     // last pad hit
		hits[13] = fFreonProd;           // did id hit the freon?
		hits[14] = mom[0];
		hits[15] = mom[1];
		hits[16] = mom[2];

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
		    ((AliRICHChamber*) (*fChambers)[idvol])
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
		      nPads = MakePadHits(xhit,yhit,eloss,idvol,kMip);
		      hits[17] = nPads;
		      //printf("nPads:%d",nPads);
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
		 ->SigGenCond(localPos[0], localPos[2], localPos[1]))
	    {
		((AliRICHChamber*) (*fChambers)[idvol])
		    ->SigGenInit(localPos[0], localPos[2], localPos[1]);
		if (eloss > 0) 
		  {
		    if(gMC->TrackPid() == kNeutron)
		      printf("\n\n\n\n\n Neutron Making Pad Hit!!! \n\n\n\n");
		    nPads = MakePadHits(xhit,yhit,eloss,idvol,kMip);
		    hits[17] = nPads;
		    //printf("Npads:%d",NPads);
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
	}
      }
    /*************************************************End of MIP treatment**************************************/
   //}
}

void AliRICH::FindClusters(Int_t nev,Int_t lastEntry)
{

//
// Loop on chambers and on cathode planes
//
    for (Int_t icat=1;icat<2;icat++) {
	gAlice->ResetDigits();
	gAlice->TreeD()->GetEvent(1); // spurious +1 ...
	for (Int_t ich=0;ich<kNCH;ich++) {
	  AliRICHChamber* iChamber=(AliRICHChamber*) (*fChambers)[ich];
	  TClonesArray *pRICHdigits  = this->DigitsAddress(ich);
	  if (pRICHdigits == 0)	      
	      continue;
	  //
	  // Get ready the current chamber stuff
	  //
	  AliRICHResponse* response = iChamber->GetResponseModel();
	  AliRICHSegmentation*  seg = iChamber->GetSegmentationModel();
	  AliRICHClusterFinder* rec = iChamber->GetReconstructionModel();
	  if (seg) {	  
	      rec->SetSegmentation(seg);
	      rec->SetResponse(response);
	      rec->SetDigits(pRICHdigits);
	      rec->SetChamber(ich);
	      if (nev==0) rec->CalibrateCOG(); 
	      rec->FindRawClusters();
	  }  
	  TClonesArray *fRch;
	  fRch=RawClustAddress(ich);
	  fRch->Sort();
	} // for ich

	gAlice->TreeR()->Fill();
	TClonesArray *fRch;
	for (int i=0;i<kNCH;i++) {
	    fRch=RawClustAddress(i);
	    int nraw=fRch->GetEntriesFast();
	    printf ("Chamber %d, raw clusters %d\n",i,nraw);
	}
	
	ResetRawClusters();
	
    } // for icat
    
    char hname[30];
    sprintf(hname,"TreeR%d",nev);
    gAlice->TreeR()->Write(hname,kOverwrite,0);
    gAlice->TreeR()->Reset();
    
    //gObjectTable->Print();
}


//______________________________________________________________________________
void AliRICH::Streamer(TBuffer &R__b)
{
    // Stream an object of class AliRICH.
    AliRICHChamber       *iChamber;
    AliRICHSegmentation  *segmentation;
    AliRICHResponse      *response;
    TClonesArray         *digitsaddress;
    TClonesArray         *rawcladdress;
    TClonesArray         *rechitaddress;
      
    if (R__b.IsReading()) {
	Version_t R__v = R__b.ReadVersion(); if (R__v) { }
	AliDetector::Streamer(R__b);
	R__b >> fNPadHits;
	R__b >> fPadHits;   // diff
	R__b >> fNcerenkovs;
	R__b >> fCerenkovs; // diff
	R__b >> fDchambers;
	R__b >> fRawClusters;
	R__b >> fRecHits;  //diff
	R__b >> fDebugLevel;  //diff
	R__b.ReadArray(fNdch);
	R__b.ReadArray(fNrawch);
	R__b.ReadArray(fNrechits);
//
	R__b >> fChambers;
// Stream chamber related information
	for (Int_t i =0; i<kNCH; i++) {
	    iChamber=(AliRICHChamber*) (*fChambers)[i];
	    iChamber->Streamer(R__b);
	    segmentation=iChamber->GetSegmentationModel();
	    segmentation->Streamer(R__b);
	    response=iChamber->GetResponseModel();
	    response->Streamer(R__b);	  
	    rawcladdress=(TClonesArray*) (*fRawClusters)[i];
	    rawcladdress->Streamer(R__b);
	    rechitaddress=(TClonesArray*) (*fRecHits)[i];
	    rechitaddress->Streamer(R__b);
	    digitsaddress=(TClonesArray*) (*fDchambers)[i];
	    digitsaddress->Streamer(R__b);
	}
      R__b >> fDebugLevel;
      R__b >> fCkovNumber;
      R__b >> fCkovQuarz;
      R__b >> fCkovGap;
      R__b >> fCkovCsi;
      R__b >> fLostRfreo;
      R__b >> fLostRquar;
      R__b >> fLostAfreo;
      R__b >> fLostAquarz;
      R__b >> fLostAmeta;
      R__b >> fLostCsi;
      R__b >> fLostWires;
      R__b >> fFreonProd;
      R__b >> fMipx;
      R__b >> fMipy;
      R__b >> fFeedbacks;
      R__b >> fLostFresnel;
      
    } else {
	R__b.WriteVersion(AliRICH::IsA());
	AliDetector::Streamer(R__b);
	R__b << fNPadHits;
	R__b << fPadHits; // diff
	R__b << fNcerenkovs;
	R__b << fCerenkovs; // diff
	R__b << fDchambers;
	R__b << fRawClusters;
	R__b << fRecHits; //diff
	R__b << fDebugLevel; //diff
	R__b.WriteArray(fNdch, kNCH);
	R__b.WriteArray(fNrawch, kNCH);
	R__b.WriteArray(fNrechits, kNCH);
	//
	R__b << fChambers;
//  Stream chamber related information
	for (Int_t i =0; i<kNCH; i++) {
	    iChamber=(AliRICHChamber*) (*fChambers)[i];
	    iChamber->Streamer(R__b);
	    segmentation=iChamber->GetSegmentationModel();
	    segmentation->Streamer(R__b);
	    response=iChamber->GetResponseModel();
	    response->Streamer(R__b);
	    rawcladdress=(TClonesArray*) (*fRawClusters)[i];
	    rawcladdress->Streamer(R__b);
	    rechitaddress=(TClonesArray*) (*fRecHits)[i];
	    rechitaddress->Streamer(R__b);
	    digitsaddress=(TClonesArray*) (*fDchambers)[i];
	    digitsaddress->Streamer(R__b);
	}
      R__b << fDebugLevel;
      R__b << fCkovNumber;
      R__b << fCkovQuarz;
      R__b << fCkovGap;
      R__b << fCkovCsi;
      R__b << fLostRfreo;
      R__b << fLostRquar;
      R__b << fLostAfreo;
      R__b << fLostAquarz;
      R__b << fLostAmeta;
      R__b << fLostCsi;
      R__b << fLostWires;
      R__b << fFreonProd;
      R__b << fMipx;
      R__b << fMipy;
      R__b << fFeedbacks;
      R__b << fLostFresnel;
    }
}
AliRICHPadHit* AliRICH::FirstPad(AliRICHHit*  hit,TClonesArray *clusters ) 
{
//
    // Initialise the pad iterator
    // Return the address of the first padhit for hit
    TClonesArray *theClusters = clusters;
    Int_t nclust = theClusters->GetEntriesFast();
    if (nclust && hit->fPHlast > 0) {
	sMaxIterPad=Int_t(hit->fPHlast);
	sCurIterPad=Int_t(hit->fPHfirst);
	return (AliRICHPadHit*) clusters->UncheckedAt(sCurIterPad-1);
    } else {
	return 0;
    }
    
}

AliRICHPadHit* AliRICH::NextPad(TClonesArray *clusters) 
{

  // Iterates over pads
  
    sCurIterPad++;
    if (sCurIterPad <= sMaxIterPad) {
	return (AliRICHPadHit*) clusters->UncheckedAt(sCurIterPad-1);
    } else {
	return 0;
    }
}


void AliRICH::Digitise(Int_t nev, Int_t flag, Option_t *option,Text_t *filename)
{
    // keep galice.root for signal and name differently the file for 
    // background when add! otherwise the track info for signal will be lost !

    static Bool_t first=kTRUE;
    static TFile *pFile;
    char *addBackground = strstr(option,"Add");

    FILE* points; //these will be the digits...

    points=fopen("points.dat","w");

    AliRICHChamber*       iChamber;
    AliRICHSegmentation*  segmentation;

    Int_t digitse=0;
    Int_t trk[50];
    Int_t chtrk[50];  
    TObjArray *list=new TObjArray;
    static TClonesArray *pAddress=0;
    if(!pAddress) pAddress=new TClonesArray("TVector",1000);
    Int_t digits[5]; 
    
    AliRICH *pRICH = (AliRICH *) gAlice->GetDetector("RICH");
    AliRICHHitMap* pHitMap[10];
    Int_t i;
    for (i=0; i<10; i++) {pHitMap[i]=0;}
    if (addBackground ) {
	if(first) {
	    fFileName=filename;
	    cout<<"filename"<<fFileName<<endl;
	    pFile=new TFile(fFileName);
	    cout<<"I have opened "<<fFileName<<" file "<<endl;
	    fHits2     = new TClonesArray("AliRICHHit",1000  );
	    fClusters2 = new TClonesArray("AliRICHPadHit",10000);
	    first=kFALSE;
	}
	pFile->cd();
	// Get Hits Tree header from file
	if(fHits2) fHits2->Clear();
	if(fClusters2) fClusters2->Clear();
	if(TrH1) delete TrH1;
	TrH1=0;

	char treeName[20];
	sprintf(treeName,"TreeH%d",nev);
	TrH1 = (TTree*)gDirectory->Get(treeName);
	if (!TrH1) {
	    printf("ERROR: cannot find Hits Tree for event:%d\n",nev);
	}
	// Set branch addresses
	TBranch *branch;
	char branchname[20];
	sprintf(branchname,"%s",GetName());
	if (TrH1 && fHits2) {
	    branch = TrH1->GetBranch(branchname);
	    if (branch) branch->SetAddress(&fHits2);
	}
	if (TrH1 && fClusters2) {
	    branch = TrH1->GetBranch("RICHCluster");
	    if (branch) branch->SetAddress(&fClusters2);
	}
    }
    
    AliRICHHitMap* hm;
    Int_t countadr=0;
    Int_t counter=0;
    for (i =0; i<kNCH; i++) {
      iChamber=(AliRICHChamber*) (*fChambers)[i];
      segmentation=iChamber->GetSegmentationModel(1);
      pHitMap[i] = new AliRICHHitMapA1(segmentation, list);
    }
    //
    //   Loop over tracks
    //
    
    TTree *treeH = gAlice->TreeH();
    Int_t ntracks =(Int_t) treeH->GetEntries();
    for (Int_t track=0; track<ntracks; track++) {
      gAlice->ResetHits();
      treeH->GetEvent(track);
      //
      //   Loop over hits
      for(AliRICHHit* mHit=(AliRICHHit*)pRICH->FirstHit(-1); 
	  mHit;
	  mHit=(AliRICHHit*)pRICH->NextHit()) 
	{
	  
	  digitse=0;
	  
	  Int_t   nch   = mHit->fChamber-1;  // chamber number
	  if (nch >kNCH) continue;
	  iChamber = &(pRICH->Chamber(nch));
	  
	  TParticle *current = (TParticle*)(*gAlice->Particles())[track];
	  
	  Int_t particle = current->GetPdgCode();
	  
	  //printf("Flag:%d\n",flag);
	  //printf("Track:%d\n",track);
	  //printf("Particle:%d\n",particle);
	  
	  if (flag == 0)
	    digitse=1;
	  
	  if (flag == 1) 
	    if(TMath::Abs(particle) == 211 || TMath::Abs(particle) == 111)
	      digitse=1;
	  
	  if (flag == 2)
	    if(TMath::Abs(particle)==321 || TMath::Abs(particle)==130 || TMath::Abs(particle)==310 
	       || TMath::Abs(particle)==311)
	      digitse=1;
	  
	  if (flag == 3 && TMath::Abs(particle)==2212)
	    digitse=1;
	  
	  if (flag == 4 && TMath::Abs(particle)==13)
	    digitse=1;
	  
	  if (flag == 5 && TMath::Abs(particle)==11)
	    digitse=1;
	  
	  if (flag == 6 && TMath::Abs(particle)==2112)
	    digitse=1;
	  
	  
	  //printf ("Particle: %d, Flag: %d, Digitse: %d\n",particle,flag,digitse); 
	  
	  
	  if (digitse)
	    {
	      
	      //
	      // Loop over pad hits
	      for (AliRICHPadHit* mPad=
		     (AliRICHPadHit*)pRICH->FirstPad(mHit,fPadHits);
		   mPad;
		   mPad=(AliRICHPadHit*)pRICH->NextPad(fPadHits))
		{
		  Int_t cathode  = mPad->fCathode;    // cathode number
		  Int_t ipx      = mPad->fPadX;       // pad number on X
		  Int_t ipy      = mPad->fPadY;       // pad number on Y
		  Int_t iqpad    = mPad->fQpad;       // charge per pad
		  //
		  //
		  //printf("X:%d, Y:%d, Q:%d\n",ipx,ipy,iqpad);
		  
		  Float_t thex, they;
		  segmentation=iChamber->GetSegmentationModel(cathode);
		  segmentation->GetPadCxy(ipx,ipy,thex,they);
		  new((*pAddress)[countadr++]) TVector(2);
		  TVector &trinfo=*((TVector*) (*pAddress)[countadr-1]);
		  trinfo(0)=(Float_t)track;
		  trinfo(1)=(Float_t)iqpad;
		  
		  digits[0]=ipx;
		  digits[1]=ipy;
		  digits[2]=iqpad;
		  
		  AliRICHTransientDigit* pdigit;
		  // build the list of fired pads and update the info
		  if (!pHitMap[nch]->TestHit(ipx, ipy)) {
		    list->AddAtAndExpand(new AliRICHTransientDigit(nch,digits),counter);
		    pHitMap[nch]->SetHit(ipx, ipy, counter);
		    counter++;
		    pdigit=(AliRICHTransientDigit*)list->At(list->GetLast());
		    // list of tracks
		    TObjArray *trlist=(TObjArray*)pdigit->TrackList();
		    trlist->Add(&trinfo);
		  } else {
		    pdigit=(AliRICHTransientDigit*) pHitMap[nch]->GetHit(ipx, ipy);
		    // update charge
		    (*pdigit).fSignal+=iqpad;
		    // update list of tracks
		    TObjArray* trlist=(TObjArray*)pdigit->TrackList();
		    Int_t lastEntry=trlist->GetLast();
		    TVector *ptrkP=(TVector*)trlist->At(lastEntry);
		    TVector &ptrk=*ptrkP;
		    Int_t lastTrack=Int_t(ptrk(0));
		    Int_t lastCharge=Int_t(ptrk(1));
		    if (lastTrack==track) {
		      lastCharge+=iqpad;
		      trlist->RemoveAt(lastEntry);
		      trinfo(0)=lastTrack;
		      trinfo(1)=lastCharge;
		      trlist->AddAt(&trinfo,lastEntry);
		    } else {
		      trlist->Add(&trinfo);
		    }
		    // check the track list
		    Int_t nptracks=trlist->GetEntriesFast();
		    if (nptracks > 2) {
		      printf("Attention - tracks:  %d (>2)\n",nptracks);
		      //printf("cat,nch,ix,iy %d %d %d %d  \n",icat+1,nch,ipx,ipy);
		      for (Int_t tr=0;tr<nptracks;tr++) {
			TVector *pptrkP=(TVector*)trlist->At(tr);
			TVector &pptrk=*pptrkP;
			trk[tr]=Int_t(pptrk(0));
			chtrk[tr]=Int_t(pptrk(1));
		      }
		    } // end if nptracks
		  } //  end if pdigit
		} //end loop over clusters
	    }// track type condition
	} // hit loop
    } // track loop
    
    // open the file with background
    
    if (addBackground ) {
      ntracks =(Int_t)TrH1->GetEntries();
      //printf("background - icat,ntracks1  %d %d\n",icat,ntracks);
      //printf("background - Start loop over tracks \n");     
      //
      //   Loop over tracks
      //
      for (Int_t trak=0; trak<ntracks; trak++) {
	if (fHits2)       fHits2->Clear();
	if (fClusters2)   fClusters2->Clear();
	TrH1->GetEvent(trak);
	//
	//   Loop over hits
	AliRICHHit* mHit;
	for(int j=0;j<fHits2->GetEntriesFast();++j) 
	  {
	    mHit=(AliRICHHit*) (*fHits2)[j];
	    Int_t   nch   = mHit->fChamber-1;  // chamber number
	    if (nch >6) continue;
	    iChamber = &(pRICH->Chamber(nch));
	    Int_t rmin = (Int_t)iChamber->RInner();
	    Int_t rmax = (Int_t)iChamber->ROuter();
	    //
	    // Loop over pad hits
	    for (AliRICHPadHit* mPad=
		   (AliRICHPadHit*)pRICH->FirstPad(mHit,fClusters2);
		 mPad;
		 mPad=(AliRICHPadHit*)pRICH->NextPad(fClusters2))
	      {
		Int_t cathode  = mPad->fCathode;    // cathode number
		Int_t ipx      = mPad->fPadX;       // pad number on X
		Int_t ipy      = mPad->fPadY;       // pad number on Y
		Int_t iqpad    = mPad->fQpad;       // charge per pad
		
		Float_t thex, they;
		segmentation=iChamber->GetSegmentationModel(cathode);
		segmentation->GetPadCxy(ipx,ipy,thex,they);
		Float_t rpad=TMath::Sqrt(thex*thex+they*they);
		if (rpad < rmin || iqpad ==0 || rpad > rmax) continue;
		new((*pAddress)[countadr++]) TVector(2);
		TVector &trinfo=*((TVector*) (*pAddress)[countadr-1]);
		trinfo(0)=-1;  // tag background
		trinfo(1)=-1;
		digits[0]=ipx;
		digits[1]=ipy;
		digits[2]=iqpad;
		if (trak <4 && nch==0)
		  printf("bgr - pHitMap[nch]->TestHit(ipx, ipy),trak %d %d\n",
			 pHitMap[nch]->TestHit(ipx, ipy),trak);
		AliRICHTransientDigit* pdigit;
		// build the list of fired pads and update the info
		if (!pHitMap[nch]->TestHit(ipx, ipy)) {
		  list->AddAtAndExpand(new AliRICHTransientDigit(nch,digits),counter);
		  
		  pHitMap[nch]->SetHit(ipx, ipy, counter);
		  counter++;
		  printf("bgr new elem in list - counter %d\n",counter);
		  
		  pdigit=(AliRICHTransientDigit*)list->At(list->GetLast());
		  // list of tracks
		  TObjArray *trlist=(TObjArray*)pdigit->TrackList();
		  trlist->Add(&trinfo);
		} else {
		  pdigit=(AliRICHTransientDigit*) pHitMap[nch]->GetHit(ipx, ipy);
		  // update charge
		  (*pdigit).fSignal+=iqpad;
		  // update list of tracks
		  TObjArray* trlist=(TObjArray*)pdigit->TrackList();
		  Int_t lastEntry=trlist->GetLast();
		  TVector *ptrkP=(TVector*)trlist->At(lastEntry);
		  TVector &ptrk=*ptrkP;
		  Int_t lastTrack=Int_t(ptrk(0));
		  if (lastTrack==-1) {
		    continue;
		  } else {
		    trlist->Add(&trinfo);
		  }
		  // check the track list
		  Int_t nptracks=trlist->GetEntriesFast();
		  if (nptracks > 0) {
		    for (Int_t tr=0;tr<nptracks;tr++) {
		      TVector *pptrkP=(TVector*)trlist->At(tr);
		      TVector &pptrk=*pptrkP;
		      trk[tr]=Int_t(pptrk(0));
		      chtrk[tr]=Int_t(pptrk(1));
		    }
		  } // end if nptracks
		} //  end if pdigit
	      } //end loop over clusters
	  } // hit loop
      } // track loop
	    TTree *fAli=gAlice->TreeK();
	    if (fAli) pFile =fAli->GetCurrentFile();
	    pFile->cd();
    } // if Add	
    
    Int_t tracks[10];
    Int_t charges[10];
    //cout<<"Start filling digits \n "<<endl;
    Int_t nentries=list->GetEntriesFast();
    //printf(" \n \n nentries %d \n",nentries);
    
    // start filling the digits
    
    for (Int_t nent=0;nent<nentries;nent++) {
      AliRICHTransientDigit *address=(AliRICHTransientDigit*)list->At(nent);
      if (address==0) continue; 
      
      Int_t ich=address->fChamber;
      Int_t q=address->fSignal; 
      iChamber=(AliRICHChamber*) (*fChambers)[ich];
      AliRICHResponse * response=iChamber->GetResponseModel();
      Int_t adcmax= (Int_t) response->MaxAdc();
      
      
      // add white noise and do zero-suppression and signal truncation (new electronics,old electronics gaus 1.2,0.2)
      //printf("Treshold: %d\n",iChamber->fTresh->GetHitIndex(address->fPadX,address->fPadY));
      Int_t pedestal = iChamber->fTresh->GetHitIndex(address->fPadX,address->fPadY);

      //printf("Pedestal:%d\n",pedestal);
      //Int_t pedestal=0;
      Float_t treshold = (pedestal + 4*1.7);
      
      Float_t meanNoise = gRandom->Gaus(1.7, 0.25);
      Float_t noise     = gRandom->Gaus(0, meanNoise);
      q+=(Int_t)(noise + pedestal);
      //q+=(Int_t)(noise);
      //          magic number to be parametrised !!! 
      if ( q <= treshold) 
	{
	  q = q - pedestal;
	  continue;
	}
      q = q - pedestal;
      if ( q >= adcmax) q=adcmax;
      digits[0]=address->fPadX;
      digits[1]=address->fPadY;
      digits[2]=q;
      
      TObjArray* trlist=(TObjArray*)address->TrackList();
      Int_t nptracks=trlist->GetEntriesFast();
      
      // this was changed to accomodate the real number of tracks
      if (nptracks > 10) {
	cout<<"Attention - tracks > 10 "<<nptracks<<endl;
	nptracks=10;
      }
      if (nptracks > 2) {
	printf("Attention - tracks > 2  %d \n",nptracks);
	//printf("cat,ich,ix,iy,q %d %d %d %d %d \n",
	//icat,ich,digits[0],digits[1],q);
      }
      for (Int_t tr=0;tr<nptracks;tr++) {
	TVector *ppP=(TVector*)trlist->At(tr);
	TVector &pp  =*ppP;
	tracks[tr]=Int_t(pp(0));
	charges[tr]=Int_t(pp(1));
      }      //end loop over list of tracks for one pad
      if (nptracks < 10 ) {
	for (Int_t t=nptracks; t<10; t++) {
	  tracks[t]=0;
	  charges[t]=0;
	}
      }
      //write file
      if (ich==2)
	fprintf(points,"%4d,      %4d,      %4d\n",digits[0],digits[1],digits[2]);
      
      // fill digits
      pRICH->AddDigits(ich,tracks,charges,digits);
    }	
    gAlice->TreeD()->Fill();
    
    list->Delete();
    for(Int_t ii=0;ii<kNCH;++ii) {
      if (pHitMap[ii]) {
	hm=pHitMap[ii];
	delete hm;
	pHitMap[ii]=0;
      }
    }
    
    //TTree *TD=gAlice->TreeD();
    //Stat_t ndig=TD->GetEntries();
    //cout<<"number of digits  "<<ndig<<endl;
    TClonesArray *fDch;
    for (int k=0;k<kNCH;k++) {
      fDch= pRICH->DigitsAddress(k);
      int ndigit=fDch->GetEntriesFast();
      printf ("Chamber %d digits %d \n",k,ndigit);
    }
    pRICH->ResetDigits();
    char hname[30];
    sprintf(hname,"TreeD%d",nev);
    gAlice->TreeD()->Write(hname,kOverwrite,0);
    
    // reset tree
    //    gAlice->TreeD()->Reset();
    delete list;
    pAddress->Clear();
    // gObjectTable->Print();
}

AliRICH& AliRICH::operator=(const AliRICH& rhs)
{
// Assignment operator
    return *this;
    
}

Int_t AliRICH::MakePadHits(Float_t xhit,Float_t yhit,Float_t eloss, Int_t idvol, ResponseType res)
{
//
//  Calls the charge disintegration method of the current chamber and adds
//  the simulated cluster to the root treee 
//
    Int_t clhits[kNCH];
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

