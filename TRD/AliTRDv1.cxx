///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Transition Radiation Detector version 1 -- coarse simulation             //
//  This version has two detector arms, leaving the space in front of the    //
//  HMPID and PHOS empty                                                     //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliTRDv1Class.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TRandom.h>
#include <TVector.h>

#include "AliTRDv1.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliConst.h"
 
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
  // --- Author :  Christoph Blume (GSI) 18/5/99 
  //
  // --- Volume names : 
  //     TRD       --> Mother TRD volume                                     (Air)
  //     UTRD      --> The detector arms                                     (Al)
  //     UTRS      --> Sectors of the sub-detector                           (Al)
  //     UTRI      --> Inner part of the detector frame                      (Air) 
  //     UTCI(N,O) --> Frames of the inner, neighbouring and outer chambers  (C) 
  //     UTII(N,O) --> Inner part of the chambers                            (Air) 
  //     UTMI(N,O) --> Modules in the chambers                               (Air) 
  //     UT0I(N,O) --> Radiator seal                                         (G10)
  //     UT1I(N,O) --> Radiator                                              (CO2)
  //     UT2I(N,O) --> Polyethylene of radiator                              (PE)
  //     UT3I(N,O) --> Entrance window                                       (Mylar)
  //     UT4I(N,O) --> Gas volume (sensitive)                                (Xe/Isobutane)
  //     UT5I(N,O) --> Pad plane                                             (Cu)
  //     UT6I(N,O) --> Support structure                                     (G10)
  //     UT7I(N,O) --> FEE + signal lines                                    (Cu)
  //     UT8I(N,O) --> Polyethylene of cooling device                        (PE)
  //     UT9I(N,O) --> Cooling water                                         (Water)
  //
  //Begin_Html
  /*
    <img src="picts/AliTRDv1.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliTRDv1Tree.gif">
  */
  //End_Html

  Float_t xpos, ypos, zpos, f;
  Int_t   idmat[5];

  const Int_t nparmo = 10;
  const Int_t nparar = 10;
  const Int_t nparfr =  4;
  const Int_t nparic =  4;
  const Int_t nparnc =  4;
  const Int_t nparoc = 11;

  Float_t par_mo[nparmo];
  Float_t par_ar[nparar];
  Float_t par_fr[nparfr];
  Float_t par_ic[nparic];
  Float_t par_nc[nparnc];
  Float_t par_oc[nparoc];
  
  Int_t *idtmed = fIdtmed->GetArray()-1299;
  
  //////////////////////////////////////////////////////////////////////////
  //     Definition of Volumes   
  //////////////////////////////////////////////////////////////////////////

  // Definition of the mother volume for the TRD (Air) 
  par_mo[0] =   0.;
  par_mo[1] = 360.;
  par_mo[2] = nsect;
  par_mo[3] = 2.;
  par_mo[4] = -zmax1;
  par_mo[5] = rmin;
  par_mo[6] = rmax;
  par_mo[7] =  zmax1;
  par_mo[8] = rmin;
  par_mo[9] = rmax;
  gMC->Gsvolu("TRD ", "PGON", idtmed[1302-1], par_mo, nparmo);
  
  Float_t phisec = 360. / nsect;   
  // Definition of the two detector arms (Al) 
  par_ar[0] = 120.;
  par_ar[1] = narmsec * phisec;
  par_ar[2] = narmsec;
  par_ar[3] = 2.;
  par_ar[4] = -zmax1;
  par_ar[5] = rmin;
  par_ar[6] = rmax;
  par_ar[7] =  zmax1;
  par_ar[8] = rmin;
  par_ar[9] = rmax;
  gMC->Gsvolu("UTRD", "PGON", idtmed[1301-1], par_ar, nparar);
  gMC->Gsdvn("UTRS", "UTRD", narmsec, 2);

  // The minimal width of a sector in rphi-direction
  Float_t widmi = rmin * TMath::Tan(kPI/nsect);
  // The maximal width of a sector in rphi-direction
  Float_t widma = rmax * TMath::Tan(kPI/nsect);
  // The total thickness of the spaceframe (Al + Air)
  Float_t frame = widmi - (widpl1 / 2);

  // Definition of the inner part of the detector frame (Air) 
  par_fr[0] = widmi - alframe / 2;
  par_fr[1] = widma - alframe / 2;
  par_fr[2] = zmax1;
  par_fr[3] = (rmax - rmin) / 2;
  gMC->Gsvolu("UTRI", "TRD1", idtmed[1302-1], par_fr, nparfr); 

  // 
  // The outer chambers
  //

  // Calculate some shape-parameter
  Float_t tanzr = (zmax1 - zmax2) / (rmax - rmin);
  Float_t theoc = -kRaddeg * TMath::ATan(tanzr / 2);

  // The carbon frame (C)
  par_oc[0]  = (rmax - rmin) / 2;
  par_oc[1]  = theoc;
  par_oc[2]  = 90.;
  par_oc[3]  = (zmax2 - zlenn - zleni/2)   / 2;
  par_oc[4]  = widmi - frame;
  par_oc[5]  = widmi - frame;
  par_oc[6]  = 0.;
  par_oc[7]  = (zmax1 - zlenn - zleni/2)   / 2;
  par_oc[8]  = widma - frame;
  par_oc[9]  = widma - frame;
  par_oc[10] = 0.;
  gMC->Gsvolu("UTCO", "TRAP", idtmed[1307-1], par_oc, nparoc);

  // The inner part (Air) 
  par_oc[3] -= ccframe;
  par_oc[4] -= ccframe;
  par_oc[5] -= ccframe; 
  par_oc[7] -= ccframe;
  par_oc[8] -= ccframe;
  par_oc[9] -= ccframe;
  gMC->Gsvolu("UTIO", "TRAP", idtmed[1302-1], par_oc, nparoc);

  // Definition of the six modules within each chamber 
  gMC->Gsdvn("UTMO", "UTIO", nmodul, 3);

  // Definition of the layers of each chamber 
  par_oc[1]  =  theoc;
  par_oc[2]  =  90.;
  par_oc[3]  = -1.;
  par_oc[4]  = -1.;
  par_oc[5]  = -1.;
  par_oc[6]  =  0.;
  par_oc[7]  = -1.;
  par_oc[8]  = -1.;
  par_oc[9]  = -1.;
  par_oc[10] =  0.;
  // G10 layer (radiator layer)
  par_oc[0] = sethick / 2;
  gMC->Gsvolu("UT0O", "TRAP", idtmed[1313-1], par_oc, nparoc);
  // CO2 layer (radiator)
  par_oc[0] = rathick / 2;
  gMC->Gsvolu("UT1O", "TRAP", idtmed[1312-1], par_oc, nparoc);
  // PE layer (radiator)
  par_oc[0] = pethick / 2;
  gMC->Gsvolu("UT2O", "TRAP", idtmed[1303-1], par_oc, nparoc);
  // Mylar layer (entrance window + HV cathode) 
  par_oc[0] = mythick / 2;
  gMC->Gsvolu("UT3O", "TRAP", idtmed[1308-1], par_oc, nparoc);
  // Xe/Isobutane layer (gasvolume)
  par_oc[0] = xethick / 2;
  gMC->Gsvolu("UT4O", "TRAP", idtmed[1309-1], par_oc, nparoc);
  // Cu layer (pad plane)
  par_oc[0] = cuthick / 2;
  gMC->Gsvolu("UT5O", "TRAP", idtmed[1305-1], par_oc, nparoc);
  // G10 layer (support structure)
  par_oc[0] = suthick / 2;
  gMC->Gsvolu("UT6O", "TRAP", idtmed[1313-1], par_oc, nparoc);
  // Cu layer (FEE + signal lines)
  par_oc[0] = fethick / 2;
  gMC->Gsvolu("UT7O", "TRAP", idtmed[1305-1], par_oc, nparoc);
  // PE layer (cooling devices)
  par_oc[0] = cothick / 2;
  gMC->Gsvolu("UT8O", "TRAP", idtmed[1303-1], par_oc, nparoc);
  // Water layer (cooling)
  par_oc[0] = wathick / 2;
  gMC->Gsvolu("UT9O", "TRAP", idtmed[1314-1], par_oc, nparoc);

  //
  // The neighbouring chambers
  //

  // The carbon frame (C) 
  par_nc[0] = widmi - frame;
  par_nc[1] = widma - frame;
  par_nc[2] = zlenn / 2;
  par_nc[3] = (rmax - rmin) / 2;
  gMC->Gsvolu("UTCN", "TRD1", idtmed[1307-1], par_nc, nparnc);

  // The inner part (Air) 
  par_nc[0] -= ccframe;
  par_nc[1] -= ccframe;
  par_nc[2] -= ccframe;
  gMC->Gsvolu("UTIN", "TRD1", idtmed[1302-1], par_nc, nparnc);

  // Definition of the six modules within each outer chamber 
  gMC->Gsdvn("UTMN", "UTIN", nmodul, 3);

  // Definition of the layers of each chamber 
  par_nc[0] = -1.;
  par_nc[1] = -1.;
  par_nc[2] = -1.;
  // G10 layer (radiator layer)
  par_nc[3] = sethick / 2;
  gMC->Gsvolu("UT0N", "TRD1", idtmed[1313-1], par_nc, nparnc);
  // CO2 layer (radiator)
  par_nc[3] = rathick / 2;
  gMC->Gsvolu("UT1N", "TRD1", idtmed[1312-1], par_nc, nparnc);
  // PE layer (radiator)
  par_nc[3] = pethick / 2;
  gMC->Gsvolu("UT2N", "TRD1", idtmed[1303-1], par_nc, nparnc);
  // Mylar layer (entrance window + HV cathode) 
  par_nc[3] = mythick / 2;
  gMC->Gsvolu("UT3N", "TRD1", idtmed[1308-1], par_nc, nparnc);
  // Xe/Isobutane layer (gasvolume)
  par_nc[3] = xethick / 2;
  gMC->Gsvolu("UT4N", "TRD1", idtmed[1309-1], par_nc, nparnc);
  // Cu layer (pad plane)
  par_nc[3] = cuthick / 2;
  gMC->Gsvolu("UT5N", "TRD1", idtmed[1305-1], par_nc, nparnc);
  // G10 layer (support structure)
  par_nc[3] = suthick / 2;
  gMC->Gsvolu("UT6N", "TRD1", idtmed[1313-1], par_nc, nparnc);
  // Cu layer (FEE + signal lines)
  par_nc[3] = fethick / 2;
  gMC->Gsvolu("UT7N", "TRD1", idtmed[1305-1], par_nc, nparnc);
  // PE layer (cooling devices)
  par_nc[3] = cothick / 2;
  gMC->Gsvolu("UT8N", "TRD1", idtmed[1303-1], par_nc, nparnc);
  // Water layer (cooling)
  par_nc[3] = wathick / 2;
  gMC->Gsvolu("UT9N", "TRD1", idtmed[1314-1], par_nc, nparnc);

  //
  // The inner chamber
  //

  // The carbon frame (C) 
  par_ic[0] = widmi - frame;
  par_ic[1] = widma - frame;
  par_ic[2] = zleni / 2;
  par_ic[3] = (rmax - rmin) / 2;
  gMC->Gsvolu("UTCI", "TRD1", idtmed[1307-1], par_ic, nparic);

  // The inner part (Air) 
  par_ic[0] -= ccframe;
  par_ic[1] -= ccframe;
  par_ic[2] -= ccframe;
  gMC->Gsvolu("UTII", "TRD1", idtmed[1302-1], par_ic, nparic);

  // Definition of the six modules within each outer chamber 
  gMC->Gsdvn("UTMI", "UTII", nmodul, 3);

  // Definition of the layers of each inner chamber 
  par_ic[0] = -1.;
  par_ic[1] = -1.;
  par_ic[2] = -1.;
  // G10 layer (radiator layer)
  par_ic[3] = sethick / 2;
  gMC->Gsvolu("UT0I", "TRD1", idtmed[1313-1], par_ic, nparic);
  // CO2 layer (radiator)
  par_ic[3] = rathick / 2;
  gMC->Gsvolu("UT1I", "TRD1", idtmed[1312-1], par_ic, nparic);
  // PE layer (radiator)
  par_ic[3] = pethick / 2;
  gMC->Gsvolu("UT2I", "TRD1", idtmed[1303-1], par_ic, nparic);
  // Mylar layer (entrance window + HV cathode) 
  par_ic[3] = mythick / 2;
  gMC->Gsvolu("UT3I", "TRD1", idtmed[1308-1], par_ic, nparic);
  // Xe/Isobutane layer (gasvolume)
  par_ic[3] = xethick / 2;
  gMC->Gsvolu("UT4I", "TRD1", idtmed[1309-1], par_ic, nparic);
  // Cu layer (pad plane)
  par_ic[3] = cuthick / 2;
  gMC->Gsvolu("UT5I", "TRD1", idtmed[1305-1], par_ic, nparic);
  // G10 layer (support structure)
  par_ic[3] = suthick / 2;
  gMC->Gsvolu("UT6I", "TRD1", idtmed[1313-1], par_ic, nparic);
  // Cu layer (FEE + signal lines)
  par_ic[3] = fethick / 2;
  gMC->Gsvolu("UT7I", "TRD1", idtmed[1305-1], par_ic, nparic);
  // PE layer (cooling devices)
  par_ic[3] = cothick / 2;
  gMC->Gsvolu("UT8I", "TRD1", idtmed[1303-1], par_ic, nparic);
  // Water layer (cooling)
  par_ic[3] = wathick / 2;
  gMC->Gsvolu("UT9I", "TRD1", idtmed[1314-1], par_ic, nparic);

  //////////////////////////////////////////////////////////////////////////
  //     Positioning of Volumes 
  //////////////////////////////////////////////////////////////////////////

  // The rotation matrices 
  AliMatrix(idmat[0],  90., 180.,  90.,  90.,   0.,   0.);
  AliMatrix(idmat[1],  90.,  90., 180.,   0.,  90.,   0.);
  AliMatrix(idmat[2],  90., 180.,  90., 270.,   0.,   0.);

  // Position of the layers in a TRD module 
  f = TMath::Tan(theoc * kDegrad);
  gMC->Gspos("UT9O", 1, "UTMO", 0., f*wazpos, wazpos, 0, "ONLY");
  gMC->Gspos("UT8O", 1, "UTMO", 0., f*cozpos, cozpos, 0, "ONLY");
  gMC->Gspos("UT7O", 1, "UTMO", 0., f*fezpos, fezpos, 0, "ONLY");
  gMC->Gspos("UT6O", 1, "UTMO", 0., f*suzpos, suzpos, 0, "ONLY");
  gMC->Gspos("UT5O", 1, "UTMO", 0., f*cuzpos, cuzpos, 0, "ONLY");
  gMC->Gspos("UT4O", 1, "UTMO", 0., f*xezpos, xezpos, 0, "ONLY");
  gMC->Gspos("UT3O", 1, "UTMO", 0., f*myzpos, myzpos, 0, "ONLY");
  gMC->Gspos("UT1O", 1, "UTMO", 0., f*razpos, razpos, 0, "ONLY");
  gMC->Gspos("UT0O", 1, "UTMO", 0., f*sezpos, sezpos, 0, "ONLY");
  gMC->Gspos("UT2O", 1, "UT1O", 0., f*pezpos, pezpos, 0, "ONLY");

  gMC->Gspos("UT9N", 1, "UTMN", 0.,       0., wazpos, 0, "ONLY");
  gMC->Gspos("UT8N", 1, "UTMN", 0.,       0., cozpos, 0, "ONLY");
  gMC->Gspos("UT7N", 1, "UTMN", 0.,       0., fezpos, 0, "ONLY");
  gMC->Gspos("UT6N", 1, "UTMN", 0.,       0., suzpos, 0, "ONLY");
  gMC->Gspos("UT5N", 1, "UTMN", 0.,       0., cuzpos, 0, "ONLY");
  gMC->Gspos("UT4N", 1, "UTMN", 0.,       0., xezpos, 0, "ONLY");
  gMC->Gspos("UT3N", 1, "UTMN", 0.,       0., myzpos, 0, "ONLY");
  gMC->Gspos("UT1N", 1, "UTMN", 0.,       0., razpos, 0, "ONLY");
  gMC->Gspos("UT0N", 1, "UTMN", 0.,       0., sezpos, 0, "ONLY");
  gMC->Gspos("UT2N", 1, "UT1N", 0.,       0., pezpos, 0, "ONLY");

  gMC->Gspos("UT9I", 1, "UTMI", 0.,       0., wazpos, 0, "ONLY");
  gMC->Gspos("UT8I", 1, "UTMI", 0.,       0., cozpos, 0, "ONLY");
  gMC->Gspos("UT7I", 1, "UTMI", 0.,       0., fezpos, 0, "ONLY");
  gMC->Gspos("UT6I", 1, "UTMI", 0.,       0., suzpos, 0, "ONLY");
  gMC->Gspos("UT5I", 1, "UTMI", 0.,       0., cuzpos, 0, "ONLY");
  gMC->Gspos("UT4I", 1, "UTMI", 0.,       0., xezpos, 0, "ONLY");
  gMC->Gspos("UT3I", 1, "UTMI", 0.,       0., myzpos, 0, "ONLY");
  gMC->Gspos("UT1I", 1, "UTMI", 0.,       0., razpos, 0, "ONLY");
  gMC->Gspos("UT0I", 1, "UTMI", 0.,       0., sezpos, 0, "ONLY");
  gMC->Gspos("UT2I", 1, "UT1I", 0.,       0., pezpos, 0, "ONLY");

  // Position of the inner part of the chambers 
  xpos = 0.;
  ypos = 0.;
  zpos = 0.;
  gMC->Gspos("UTII", 1, "UTCI", xpos, ypos, zpos, 0, "ONLY");
  gMC->Gspos("UTIN", 1, "UTCN", xpos, ypos, zpos, 0, "ONLY");
  gMC->Gspos("UTIO", 1, "UTCO", xpos, ypos, zpos, 0, "ONLY");

  // Position of the chambers in the support frame 
  xpos = 0.;
  ypos = ((zmax1 + zmax2) / 2 + zlenn + zleni / 2) / 2;
  zpos = 0.;
  gMC->Gspos("UTCO", 1, "UTRI", xpos, ypos, zpos, idmat[2], "ONLY");
  gMC->Gspos("UTCO", 2, "UTRI", xpos,-ypos, zpos,       0 , "ONLY");
  xpos = 0.;
  ypos = (zlenn + zleni) / 2;
  zpos = 0.;
  gMC->Gspos("UTCN", 1, "UTRI", xpos, ypos, zpos,       0 , "ONLY");
  gMC->Gspos("UTCN", 2, "UTRI", xpos,-ypos, zpos,       0 , "ONLY");
  xpos = 0.;
  ypos = 0.;
  zpos = 0.;
  gMC->Gspos("UTCI", 1, "UTRI", xpos, ypos, zpos,       0 , "ONLY");

  // Position of the inner part of the detector frame
  xpos = (rmax + rmin) / 2;
  ypos = 0.;
  zpos = 0.;
  gMC->Gspos("UTRI", 1, "UTRS", xpos, ypos, zpos, idmat[1], "ONLY");

  // Position of the two arms of the detector
  xpos = 0.;
  ypos = 0.;
  zpos = 0.;
  gMC->Gspos("UTRD", 1, "TRD ", xpos, ypos, zpos,        0, "ONLY");
  gMC->Gspos("UTRD", 2, "TRD ", xpos, ypos, zpos, idmat[0], "ONLY");

  // Position of TRD mother volume in ALICE experiment 
  xpos = 0.;
  ypos = 0.;
  zpos = 0.;
  gMC->Gspos("TRD ", 1, "ALIC", xpos, ypos, zpos,        0, "ONLY");

}
 
//_____________________________________________________________________________
void AliTRDv1::DrawModule()
{
  //
  // Draw a shaded view of the Transition Radiation Detector version 1
  //

  // Set everything unseen
  gMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  gMC->Gsatt("TRD" ,"SEEN",0);
  gMC->Gsatt("UTRD","SEEN",0);
  gMC->Gsatt("UTRS","SEEN",0);
  gMC->Gsatt("UTRI","SEEN",0);
  gMC->Gsatt("UTCO","SEEN",0);
  gMC->Gsatt("UTIO","SEEN",0);
  gMC->Gsatt("UTMO","SEEN",0);
  gMC->Gsatt("UTCN","SEEN",0);
  gMC->Gsatt("UTIN","SEEN",0);
  gMC->Gsatt("UTMN","SEEN",0);
  gMC->Gsatt("UTCI","SEEN",0);
  gMC->Gsatt("UTII","SEEN",0);
  gMC->Gsatt("UTMI","SEEN",0);
  gMC->Gsatt("UT1O","SEEN",1);
  gMC->Gsatt("UT4O","SEEN",1);
  gMC->Gsatt("UT1N","SEEN",1);
  gMC->Gsatt("UT4N","SEEN",1);
  gMC->Gsatt("UT1I","SEEN",1);
  gMC->Gsatt("UT4I","SEEN",1);
  //
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 2000, -2000, 2000, -2000, 2000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 40, 30, 0, 12, 9.4, .021, .021);
  gMC->Gdhead(1111, "Transition Radiation Detector Version 1");
  gMC->Gdman(18, 4, "MAN");
}

//_____________________________________________________________________________
void AliTRDv1::CreateMaterials()
{
  //
  // Create materials for the Transition Radiation Detector version 1
  //
  AliTRD::CreateMaterials();
}

//_____________________________________________________________________________
void AliTRDv1::Init() 
{
  //
  // Initialise the Transition Radiation Detector after the geometry is built
  //
  AliTRD::Init();
  // Retrieve the numeric identifier of the sensitive volumes (gas volume)
  fIdSens1 = gMC->VolId("UT4I");
  fIdSens2 = gMC->VolId("UT4N");
  fIdSens3 = gMC->VolId("UT4O");
}

//_____________________________________________________________________________
void AliTRDv1::StepManager() 
{
  //
  // Procedure called at every step in the TRD
  //

  Int_t         vol[3]; 
  Int_t         icopy1, icopy2;
  Int_t         idSens, icSens; 
  
  Float_t       hits[4];
  
  TClonesArray &lhits = *fHits;

  // Use only charged tracks and count them only once per volume
  if (gMC->TrackCharge() && gMC->TrackExiting()) {
    
    // Check on sensitive volume
    idSens = gMC->CurrentVol(0,icSens);

    // Check on sensitive volume
    idSens = gMC->CurrentVol(0,icSens);
    if ((idSens == fIdSens1) || 
        (idSens == fIdSens2) ||
        (idSens == fIdSens3)) { 
      
      // The sector number
      gMC->CurrentVolOff(5,0,icopy1);
      gMC->CurrentVolOff(6,0,icopy2);
      if (icopy2 == 1)
        vol[0] =     icopy1;
      else
        vol[0] = 6 - icopy1 + 5;
      
      // The chamber number 
      //   1: outer left
      //   2: neighbouring left
      //   3: inner
      //   4: neighbouring right
      //   5: outer right
      gMC->CurrentVolOff(3,0,icopy1);
      if      (idSens == fIdSens3)
        vol[1] = 4 * icopy1 - 3; 
      else if (idSens == fIdSens2)
        vol[1] = 2 * icopy1;
      else 
        vol[1] = 3;
      
      // The plane number
      gMC->CurrentVolOff(1,0,icopy1);
      vol[2] = icopy1;

      if (fSensSelect) {
        Int_t addthishit = 1;
        if ((fSensPlane)   && (vol[2] != fSensPlane  )) addthishit = 0;
        if ((fSensChamber) && (vol[1] != fSensChamber)) addthishit = 0;
        if ((fSensSector)  && (vol[0] != fSensSector )) addthishit = 0;
        if (addthishit) {
          gMC->TrackPosition(hits);
          hits[3] = 0;
          new(lhits[fNhits++]) AliTRDhit(fIshunt,gAlice->CurrentTrack(),vol,hits);
	}
      }
      else {      
        gMC->TrackPosition(hits);
        hits[3] = 0;
        new(lhits[fNhits++]) AliTRDhit(fIshunt,gAlice->CurrentTrack(),vol,hits);
      }

    }

  } 

}
