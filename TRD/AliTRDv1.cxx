///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Transition Radiation Detector version 1 -- detailed simulation           //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliTRDv1Class.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TRandom.h>
#include <TVector.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TPGON.h> 

#include "GParticle.h"
#include "AliTRDv1.h"
#include "AliRun.h"
#include "AliConst.h"
#include "AliMC.h"
 
ClassImp(AliTRDv1)

//_____________________________________________________________________________
AliTRDv1::AliTRDv1(const char *name, const char *title) 
         :AliTRD(name, title) 
{
  //
  // Standard constructor for the Transition Radiation Detector version 1
  //
  fIdSens1 = fIdSens2 = fIdSens3 = 0;
}
 
//_____________________________________________________________________________
void AliTRDv1::CreateGeometry()
{
  //
  // Create the geometry for the Transition Radiation Detector version 1
  // --- The coarse geometry of the TRD, that can be used for background 
  //     studies. This version leaves the space in front of the PHOS and 
  //     HMPID empty. 
  // -- Author :    Nick van Eijndhoven (CERN)   24/09/90 
  //
  //Begin_Html
  /*
    <img src="gif/AliTRDv1.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="gif/AliTRDv1Tree.gif">
  */
  //End_Html

  Float_t xpos, ypos, zpos, f;
  Int_t idmat[5];
  Float_t widma, theoc, widmi, tanzr;
  Float_t par_ic[4], par_oc[11], phisec, par_mo[10], par_fr[4], par_su[10];
  
  Int_t *idtmed = gAlice->Idtmed();
  
  // --- Name Conventions : 
  //        TRD       --> Mother TRD volume                       (Air) 
  //        UTRL(S)   --> Long (short) subdetector-type           (Al) 
  //        UTSL(S)   --> Sectors of a subdetector                (Al) 
  //        UTFI(O/S) --> Inner part of the detector frame        (Air) 
  //        UTCI(O/S) --> Frames of the inner and outer chambers  (C) 
  //        UTII(O/S) --> Inner part of the chambers              (Air) 
  //        UTMI(O/S) --> Modules in the chambers                 (Air) 
  //        UT1I(O/S) --> Radiator layer                          (CO2) 
  //        UT2I(O/S) --> Polyethylene layer                      (PE) 
  //        UT3I(O/S) --> Mylar layer                             (Mylar) 
  //        UT4I(O/S) --> Xe/C02 layer                            (Xe/C02) 
  //        UT5I(O/S) --> Cu layer (pads/sensitive)               (Cu) 
  //        UT6I(O/S) --> Kapton layer                            (Kapton) 
  //        UT7I(O/S) --> NOMEX layer                             (C) 
  //        UT8I(O/S) --> Readout layer                           (Al) 
  
  // --- Contains geometry information 
  
  // --- Number of sectors in the full detector 
  // --- Number of modules in each sector 
  // --- z-Coordinates of the TRD-frame 
  // --- r-Coordinates of the TRD-frame 
  // --- Thickness of the aluminium of the support frame 
  // --- Thickness of the interior of the support frame 
  // --- Thickness of the carbon chamber frame 
  // --- Thickness and z-position of the PE-layer in the radiator 
  // --- Thickness and z-position of the radiator 
  // --- Thickness and z-position of the mylar-layer 
  // --- Thickness and z-position of the Xe/C02-layer 
  // --- Thickness and z-position of the Cu-layer (Pads) 
  // --- Thickness and z-position of the kapton-layer 
  // --- Thickness and z-position of the NOMEX-layer 
  //     Simple C-layer for the time being 
  // --- Thickness and z-position of the readout-layer 
  // --- Parameter for the arrays 
  // --- Number of subdetector-types 
  //--- Number of sectors in the first subdetector-type (full theta coverage)
  //--- Number of sectors in the second subdetector-type (with hole for PHOS)
  //************************************************************************
  
  //     Definition of Volumes 
  
  //************************************************************************
  
  const Int_t nsec1 = 5; //Number of sectors in the first subdetector-type
  const Int_t nsec2 = 5; //Number of sectors in the second subdetector-type
  
  AliMC* pMC = AliMC::GetMC();
  
  phisec = 360./nsect;  //The phi-angle of the sectors 
  widmi  = rmin*TMath::Sin(kPI/nsect);
  widma  = rmax*TMath::Sin(kPI/nsect);
  // --- Definition of the Mother volume for the TRD (Al) 
  par_mo[0] = 0.;
  par_mo[1] = 360.;
  par_mo[2] = nsect;
  par_mo[3] = 2.;
  par_mo[4] = -zmax1;
  par_mo[5] = rmin;
  par_mo[6] = rmax;
  par_mo[7] = zmax1;
  par_mo[8] = rmin;
  par_mo[9] = rmax;
  pMC->Gsvolu("TRD ", "PGON", idtmed[1301], par_mo, 10);
  // --- Definition of the 1st subdetector-type (full theta-coverage) (Al) 
  par_su[0] = 120.;
  par_su[1] = nsec1*phisec;
  par_su[2] = nsec1;
  par_su[3] = 2.;
  par_su[4] = -zmax1;
  par_su[5] = rmin;
  par_su[6] = rmax;
  par_su[7] = 0.;
  par_su[8] = rmin;
  par_su[9] = rmax;
  pMC->Gsvolu("UTRL", "PGON", idtmed[1300], par_su, 10);
  pMC->Gsdvn("UTSL", "UTRL", nsec1, 2);
  // --- Definition of the 2nd subdetector-type (hole for PHOS) 
  par_su[0] = 220.;
  par_su[1] = nsec2*phisec;
  par_su[2] = nsec2;
  par_su[3] = 2.;
  par_su[4] = -zmax1;
  par_su[5] = rmin;
  par_su[6] = rmax;
  par_su[7] = -zmax1/2;
  par_su[8] = rmin;
  par_su[9] = rmax;
  pMC->Gsvolu("UTRS", "PGON", idtmed[1300], par_su, 10);
  pMC->Gsdvn("UTSS", "UTRS", nsec2, 2);
  // --- Definition of the inner part of the detector frame (Air) 
  par_fr[0] = widmi;
  par_fr[1] = widma;
  par_fr[2] = zmax1/4 - alfram2/2;
  par_fr[3] = (rmax-rmin)/2;
  pMC->Gsvolu("UTFI", "TRD1", idtmed[1301], par_fr, 4);
  pMC->Gsvolu("UTFO", "TRD1", idtmed[1301], par_fr, 4);
  pMC->Gsvolu("UTFS", "TRD1", idtmed[1301], par_fr, 4);
  // --- Calculate the shape-parameter for the outer chambers 
  tanzr = (zmax1-zmax2)/(rmax-rmin);
  theoc = -kRaddeg*TMath::ATan(tanzr / 2.);
  // --- The carbon frame of the outer chambers 
  par_oc[0] = (rmax-rmin)/2;
  par_oc[1] = theoc;
  par_oc[2] = 90.;
  par_oc[3] = zmax2/2 - zmax1/4 -alfram2/2;
  par_oc[4] = widmi - (inframe+alfram1)/2;
  par_oc[5] = widmi - (inframe+alfram1)/2;
  par_oc[6] = 0.;
  par_oc[7] = zmax1/4 - alfram2/2;
  par_oc[8] = widma - (inframe+alfram1)/2;
  par_oc[9] = widma - (inframe+alfram1)/2;
  par_oc[10] = 0.;
  pMC->Gsvolu("UTCO", "TRAP", idtmed[1306], par_oc, 11);
  // --- The inner part of the outer chambers (Air) 
  par_oc[3] -= ccframe;
  par_oc[4] -= ccframe;
  par_oc[5] -= ccframe;
  par_oc[7] -= ccframe;
  par_oc[8] -= ccframe;
  par_oc[9] -= ccframe;
  pMC->Gsvolu("UTIO", "TRAP", idtmed[1301], par_oc, 11);
  // --- Definition of the six modules within each outer chamber 
  pMC->Gsdvn("UTMO", "UTIO", nmodul, 3);
  // --- Definition of the layers of each outer chamber 
  par_oc[1] = theoc;
  par_oc[2] = 90.;
  par_oc[3] = -1.;
  par_oc[4] = -1.;
  par_oc[5] = -1.;
  par_oc[6] =  0.;
  par_oc[7] = -1.;
  par_oc[8] = -1.;
  par_oc[9] = -1.;
  par_oc[10] = 0.;
  // --- Radiator layer 
  par_oc[0] = rathick/2;
  pMC->Gsvolu("UT1O", "TRAP", idtmed[1311], par_oc, 11);
  // --- Polyethylene layer 
  par_oc[0] = pethick/2;
  pMC->Gsvolu("UT2O", "TRAP", idtmed[1302], par_oc, 11);
  // --- Mylar layer 
  par_oc[0] = mythick/2;
  pMC->Gsvolu("UT3O", "TRAP", idtmed[1307], par_oc, 11);
  // --- Xe/CO2 layer 
  par_oc[0] = xethick/2;
  pMC->Gsvolu("UT4O", "TRAP", idtmed[1308], par_oc, 11);
  // --- Cu layer 
  par_oc[0] = cuthick/2;
  pMC->Gsvolu("UT5O", "TRAP", idtmed[1304], par_oc, 11);
  // --- Kapton layer 
  par_oc[0] = kathick/2;
  pMC->Gsvolu("UT6O", "TRAP", idtmed[1310], par_oc, 11);
  // --- NOMEX layer 
  par_oc[0] = nothick/2;
  pMC->Gsvolu("UT7O", "TRAP", idtmed[1309], par_oc, 11);
  // --- Read out layer 
  par_oc[0] = rothick/2;
  pMC->Gsvolu("UT8O", "TRAP", idtmed[1305], par_oc, 11);
  // --- The carbon frame of the chambers in the short sectors 
  par_oc[0] = (rmax-rmin)/2;
  par_oc[1] = theoc;
  par_oc[2] = 90.;
  par_oc[3] = zmax2/2 - zmax1/4 -alfram2/2;
  par_oc[4] = widmi - (inframe+alfram1)/2;
  par_oc[5] = widmi - (inframe+alfram1)/2;
  par_oc[6] = 0.;
  par_oc[7] = zmax1/4 - alfram2/2;
  par_oc[8] = widma - (inframe+alfram1)/2;
  par_oc[9] = widma - (inframe+alfram1)/2;
  par_oc[10] = 0.;
  pMC->Gsvolu("UTCS", "TRAP", idtmed[1306], par_oc, 11);
  // --- The inner part of the chambers in the short sectors (Air) 
  par_oc[3] -= ccframe;
  par_oc[4] -= ccframe;
  par_oc[5] -= ccframe;
  par_oc[7] -= ccframe;
  par_oc[8] -= ccframe;
  par_oc[9] -= ccframe;
  pMC->Gsvolu("UTIS", "TRAP", idtmed[1301], par_oc, 11);
  //--- Definition of the six modules within each chamber of the short sectors
  pMC->Gsdvn("UTMS", "UTIS", 6, 3);
  // --- Definition of the layers of each chamber in the short sectors 
  par_oc[1] = theoc;
  par_oc[2] = 90.;
  par_oc[3] = -1.;
  par_oc[4] = -1.;
  par_oc[5] = -1.;
  par_oc[6] =  0.;
  par_oc[7] = -1.;
  par_oc[8] = -1.;
  par_oc[9] = -1.;
  par_oc[10] = 0.;
  // --- Radiator layer 
  par_oc[0] = rathick/2;
  pMC->Gsvolu("UT1S", "TRAP", idtmed[1311], par_oc, 11);
  // --- Polyethylene layer 
  par_oc[0] = pethick/2;
  pMC->Gsvolu("UT2S", "TRAP", idtmed[1302], par_oc, 11);
  // --- Mylar layer 
  par_oc[0] = mythick/2;
  pMC->Gsvolu("UT3S", "TRAP", idtmed[1307], par_oc, 11);
  // --- Xe/CO2 layer 
  par_oc[0] = xethick/2;
  pMC->Gsvolu("UT4S", "TRAP", idtmed[1308], par_oc, 11);
  // --- Cu layer 
  par_oc[0] = cuthick/2;
  pMC->Gsvolu("UT5S", "TRAP", idtmed[1304], par_oc, 11);
  // --- Kapton layer 
  par_oc[0] = kathick/2;
  pMC->Gsvolu("UT6S", "TRAP", idtmed[1310], par_oc, 11);
  // --- NOMEX layer 
  par_oc[0] = nothick/2;
  pMC->Gsvolu("UT7S", "TRAP", idtmed[1309], par_oc, 11);
  // --- Read out layer 
  par_oc[0] = rothick/2;
  pMC->Gsvolu("UT8S", "TRAP", idtmed[1305], par_oc, 11);
  // --- The carbon frame of the inner chambers 
  par_ic[0] = widmi - (inframe+alfram1)/2;
  par_ic[1] = widma - (inframe+alfram1)/2;
  par_ic[2] = zmax1/4 - alfram2/2;
  par_ic[3] = (rmax-rmin)/2;
  pMC->Gsvolu("UTCI", "TRD1", idtmed[1306], par_ic, 4);
  // --- The inner part of the inner chambers (Air) 
  par_ic[0] -= ccframe;
  par_ic[1] -= ccframe;
  par_ic[2] -= ccframe;
  pMC->Gsvolu("UTII", "TRD1", idtmed[1301], par_ic, 4);
  // --- Definition of the six modules within each outer chamber 
  pMC->Gsdvn("UTMI", "UTII", nmodul, 3);
  // --- Definition of the layers of each inner chamber 
  par_ic[0] = -1.;
  par_ic[1] = -1.;
  par_ic[2] = -1.;
  // --- Radiator layer 
  par_ic[3] = rathick/2;
  pMC->Gsvolu("UT1I", "TRD1", idtmed[1311], par_ic, 4);
  // --- Polyethylene layer 
  par_ic[3] = pethick/2;
  pMC->Gsvolu("UT2I", "TRD1", idtmed[1302], par_ic, 4);
  // --- Mylar layer 
  par_ic[3] = mythick/2;
  pMC->Gsvolu("UT3I", "TRD1", idtmed[1307], par_ic, 4);
  // --- Xe/CO2 layer 
  par_ic[3] = xethick/2;
  pMC->Gsvolu("UT4I", "TRD1", idtmed[1308], par_ic, 4);
  // --- Cu layer 
  par_ic[3] = cuthick/2;
  pMC->Gsvolu("UT5I", "TRD1", idtmed[1304], par_ic, 4);
  // --- Kapton layer 
  par_ic[3] = kathick/2;
  pMC->Gsvolu("UT6I", "TRD1", idtmed[1310], par_ic, 4);
  // --- NOMEX layer 
  par_ic[3] = nothick/2;
  pMC->Gsvolu("UT7I", "TRD1", idtmed[1309], par_ic, 4);
  // --- Read out layer 
  par_ic[3] = rothick/2;
  pMC->Gsvolu("UT8I", "TRD1", idtmed[1305], par_ic, 4);
  //************************************************************************
  
  //     Positioning of Volumes 
  
  //************************************************************************
  // --- The rotation matrices 
  AliMatrix(idmat[0], 90., 180.,  90., 90.,   0., 0.);
  AliMatrix(idmat[1], 90.,   0.,  90., 90., 180., 0.);
  AliMatrix(idmat[2], 90., 180.,  90., 90., 180., 0.);
  AliMatrix(idmat[3], 90.,  90., 180.,  0.,  90., 0.);
  AliMatrix(idmat[4], 90.,  90.,   0.,  0.,  90., 0.);
  // --- Position of the layers in a TRD module 
  f = TMath::Tan(theoc * kDegrad);
  pMC->Gspos("UT8O", 1, "UTMO", 0., f*rozpos, rozpos, 0, "ONLY");
  pMC->Gspos("UT7O", 1, "UTMO", 0., f*nozpos, nozpos, 0, "ONLY");
  pMC->Gspos("UT6O", 1, "UTMO", 0., f*kazpos, kazpos, 0, "ONLY");
  pMC->Gspos("UT5O", 1, "UTMO", 0., f*cuzpos, cuzpos, 0, "ONLY");
  pMC->Gspos("UT4O", 1, "UTMO", 0., f*xezpos, xezpos, 0, "ONLY");
  pMC->Gspos("UT3O", 1, "UTMO", 0., f*myzpos, myzpos, 0, "ONLY");
  pMC->Gspos("UT1O", 1, "UTMO", 0., f*razpos, razpos, 0, "ONLY");
  pMC->Gspos("UT2O", 1, "UT1O", 0., f*pezpos, pezpos, 0, "ONLY");
  
  pMC->Gspos("UT8S", 1, "UTMS", 0., f*rozpos, rozpos, 0, "ONLY");
  pMC->Gspos("UT7S", 1, "UTMS", 0., f*nozpos, nozpos, 0, "ONLY");
  pMC->Gspos("UT6S", 1, "UTMS", 0., f*kazpos, kazpos, 0, "ONLY");
  pMC->Gspos("UT5S", 1, "UTMS", 0., f*cuzpos, cuzpos, 0, "ONLY");
  pMC->Gspos("UT4S", 1, "UTMS", 0., f*xezpos, xezpos, 0, "ONLY");
  pMC->Gspos("UT3S", 1, "UTMS", 0., f*myzpos, myzpos, 0, "ONLY");
  pMC->Gspos("UT1S", 1, "UTMS", 0., f*razpos, razpos, 0, "ONLY");
  pMC->Gspos("UT2S", 1, "UT1S", 0., f*pezpos, pezpos, 0, "ONLY");
  
  pMC->Gspos("UT8I", 1, "UTMI", 0., 0., rozpos, 0, "ONLY");
  pMC->Gspos("UT7I", 1, "UTMI", 0., 0., nozpos, 0, "ONLY");
  pMC->Gspos("UT6I", 1, "UTMI", 0., 0., kazpos, 0, "ONLY");
  pMC->Gspos("UT5I", 1, "UTMI", 0., 0., cuzpos, 0, "ONLY");
  pMC->Gspos("UT4I", 1, "UTMI", 0., 0., xezpos, 0, "ONLY");
  pMC->Gspos("UT3I", 1, "UTMI", 0., 0., myzpos, 0, "ONLY");
  pMC->Gspos("UT1I", 1, "UTMI", 0., 0., razpos, 0, "ONLY");
  pMC->Gspos("UT2I", 1, "UT1I", 0., 0., pezpos, 0, "ONLY");
  // --- Position of the inner part of the chambers 
  pMC->Gspos("UTII", 1, "UTCI", 0., 0., 0., 0, "ONLY");
  pMC->Gspos("UTIO", 1, "UTCO", 0., 0., 0., 0, "ONLY");
  pMC->Gspos("UTIS", 1, "UTCS", 0., 0., 0., 0, "ONLY");
  // --- Position of the chambers in the support frame 
  xpos = 0.;
  ypos = (zmax1-zmax2)/4;
  zpos = 0.;
  pMC->Gspos("UTCO", 1, "UTFO", xpos, ypos, zpos, 0, "ONLY");
  xpos = 0.;
  ypos = (zmax1-zmax2)/4;
  zpos = 0.;
  pMC->Gspos("UTCS", 1, "UTFS", xpos, ypos, zpos, 0, "ONLY");
  xpos = 0.;
  ypos = 0.;
  zpos = 0.;
  pMC->Gspos("UTCI", 1, "UTFI", xpos, ypos, zpos, 0, "ONLY");
  // --- Position of the inner part of the frame in the sectors 
  xpos = (rmax+rmin)/2;
  ypos = 0;
  zpos = -zmax1*3/4;
  pMC->Gspos("UTFO", 1, "UTSL", xpos, ypos, zpos, idmat[4], "ONLY");
  xpos = (rmax+rmin)/2;
  ypos = 0;
  zpos = -zmax1*3/4;
  pMC->Gspos("UTFS", 1, "UTSS", xpos, ypos, zpos, idmat[4], "ONLY");
  xpos = (rmax+rmin)/2;
  ypos = 0.;
  zpos = -zmax1/4;
  pMC->Gspos("UTFI", 1, "UTSL", xpos, ypos, zpos, idmat[4], "ONLY");
  // --- Position of the subdetectors in the mother volume 
  pMC->Gspos("UTRL", 1, "TRD ", 0., 0., 0., 0,        "ONLY");
  pMC->Gspos("UTRL", 2, "TRD ", 0., 0., 0., idmat[0], "ONLY");
  pMC->Gspos("UTRL", 3, "TRD ", 0., 0., 0., idmat[1], "ONLY");
  pMC->Gspos("UTRL", 4, "TRD ", 0., 0., 0., idmat[2], "ONLY");
  pMC->Gspos("UTRS", 1, "TRD ", 0., 0., 0., 0,        "ONLY");
  pMC->Gspos("UTRS", 2, "TRD ", 0., 0., 0., idmat[1], "ONLY");
  // --- Position of TRD mother volume in ALICE experiment 
  pMC->Gspos("TRD ", 1, "ALIC", 0., 0., 0., 0, "ONLY");
}
 
//_____________________________________________________________________________
void AliTRDv1::DrawModule()
{
  //
  // Draw a shaded view of the Transition Radiation Detector version 1
  //

  AliMC* pMC = AliMC::GetMC();
  
  // Set everything unseen
  pMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  pMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  pMC->Gsatt("TRD","SEEN",0);
  pMC->Gsatt("UTRL","SEEN",0);
  pMC->Gsatt("UTSL","SEEN",0);
  pMC->Gsatt("UTRS","SEEN",0);
  pMC->Gsatt("UTSS","SEEN",0);
  pMC->Gsatt("UTFI","SEEN",0);
  pMC->Gsatt("UTFO","SEEN",0);
  pMC->Gsatt("UTFS","SEEN",0);
  pMC->Gsatt("UTCO","SEEN",0);
  pMC->Gsatt("UTIO","SEEN",0);
  pMC->Gsatt("UTMO","SEEN",0);
  pMC->Gsatt("UT1O","SEEN",1);
  pMC->Gsatt("UT4O","SEEN",1);
  pMC->Gsatt("UTCS","SEEN",0);
  pMC->Gsatt("UTIS","SEEN",0);
  pMC->Gsatt("UTMS","SEEN",0);
  pMC->Gsatt("UT1S","SEEN",1);
  pMC->Gsatt("UT4S","SEEN",1);
  pMC->Gsatt("UTCI","SEEN",0);
  pMC->Gsatt("UTII","SEEN",0);
  pMC->Gsatt("UTMI","SEEN",0);
  pMC->Gsatt("UT1I","SEEN",1);
  pMC->Gsatt("UT4I","SEEN",1);
  //
  pMC->Gdopt("hide", "on");
  pMC->Gdopt("shad", "on");
  pMC->Gsatt("*", "fill", 7);
  pMC->SetClipBox(".");
  pMC->SetClipBox("*", 0, 2000, -2000, 2000, -2000, 2000);
  pMC->DefaultRange();
  pMC->Gdraw("alic", 40, 30, 0, 12, 9.4, .021, .021);
  pMC->Gdhead(1111, "Transition Radiation Detector Version 1");
  pMC->Gdman(18, 4, "MAN");
}

//_____________________________________________________________________________
void AliTRDv1::CreateMaterials()
{
  //
  // Create materials for the Transition Radiation Detector version 1
  //
   printf("TRD: Fast simulation with coarse geometry\n");
   AliTRD::CreateMaterials();
}

//_____________________________________________________________________________
void AliTRDv1::Init() 
{
  //
  // Initialise the Transition Radiation Detector after the geometry is built
  //
  AliTRD::Init();
  AliMC* pMC = AliMC::GetMC();
  fIdSens1 = pMC->VolId("UT5I");
  fIdSens2 = pMC->VolId("UT5O");
  fIdSens3 = pMC->VolId("UT5S");
}

//_____________________________________________________________________________
void AliTRDv1::StepManager() 
{
  //
  // Called at every step in the Transition Radiation Detector
  //
  Int_t         vol[3]; 
  Int_t         icopy1, icopy5, icopy6, idSens, icSens; 
  
  Float_t       hits[4];
  
  TClonesArray &lhits = *fHits;

  AliMC* pMC = AliMC::GetMC();
  
  // Use only charged tracks and count them only once per volume
  if(pMC->TrackCharge() && pMC->TrackEntering()) {
    
    // Check on sensitive volume
    idSens = pMC->CurrentVol(0,icSens);

    // Long sectors
    if      ((idSens == fIdSens1) || (idSens == fIdSens2)) {
      
      pMC->CurrentVolOff(1,0,icopy1);
      pMC->CurrentVolOff(5,0,icopy5);
      pMC->CurrentVolOff(6,0,icopy6);
      
      // The sector number
      if ((icopy6 == 1) || (icopy6 == 3)) 
	vol[0] = icopy5;
      else 
	vol[0] = 16 - icopy5;
      
      // The chamber number
      if (idSens == fIdSens2) { 
        if (icopy6 < 3)
          vol[1] = 4;
        else
          vol[1] = 1;
      }         
      else {
        if (icopy6 < 3)
          vol[1] = 3;
        else
          vol[1] = 2;
      }
      
      // The plane number
      vol[2] = icopy1;    
      
      pMC->TrackPosition(hits);
      hits[3] = 0;
      
      new(lhits[fNhits++]) AliTRDhit(fIshunt,gAlice->CurrentTrack(),vol,hits);
      
    }
    // Short sectors
    else if (idSens == fIdSens3) {
      
      pMC->CurrentVolOff(1,0,icopy1);
      pMC->CurrentVolOff(5,0,icopy5);
      pMC->CurrentVolOff(6,0,icopy6);
      
      // The sector number
      vol[0] = icopy5 + 5;
      
      // The chamber number
      if (icopy6 == 1)
        vol[1] = 4;
      else
        vol[1] = 1;
      
      // The plane number
      vol[2] = icopy1;
      
      pMC->TrackPosition(hits);        
      hits[3] = 0;
      
      new(lhits[fNhits++]) AliTRDhit(fIshunt,gAlice->CurrentTrack(),vol,hits);
      
    }   
  } 
}
