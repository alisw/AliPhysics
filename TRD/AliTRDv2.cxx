///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Transition Radiation Detector version 2 -- detailed simulation           //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliTRDv2Class.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TVector.h>

#include "AliTRDv2.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliConst.h"

ClassImp(AliTRDv2)

//_____________________________________________________________________________
AliTRDv2::AliTRDv2(const char *name, const char *title) 
         :AliTRD(name, title) 
{
  //
  // Standard constructor for Transition Radiation Detector version 2
  //
  for (Int_t icham = 0; icham < ncham; ++icham) {
    fIdSensI[icham] = 0;
    fIdSensN[icham] = 0;
    fIdSensO[icham] = 0;
  }
  fDeltaE = NULL;

  SetBufferSize(128000);
}

AliTRDv2::~AliTRDv2()
{
   if (fDeltaE)  delete fDeltaE;
}
 
//_____________________________________________________________________________
void AliTRDv2::CreateGeometry()
{
  //
  // Create geometry for the Transition Radiation Detector version 2
  // This version covers the full azimuth. 
  // --- Author :  Christoph Blume (GSI) 20/5/99 
  //
  // --- Volume names : 
  //       TRD         --> Mother TRD volume                              (Al) 
  //       UTRS        --> Sectors of the sub-detector                    (Al)
  //       UTRI        --> Inner part of the detector frame               (Air) 
  //     The chambers 
  //       UCI1-6      --> The frame of the inner chambers                (C) 
  //       UCN1-6      --> The frame of the neighbouring chambers         (C) 
  //       UCO1-6      --> The frame of the outer chambers                (C) 
  //       UII1-6      --> The inner part of the inner chambers           (Air) 
  //       UIN1-6      --> The inner part of the neighbouring chambers    (Air) 
  //       UIO1-6      --> The inner part of the outer chambers           (Air) 
  //     The layers inside a chamber 
  //       UT0I(N,O)   --> Radiator seal                                  (G10)
  //       UT1I(N,O)   --> Radiator                                       (CO2)
  //       UT2I(N,O)   --> Polyethylene of radiator                       (PE)
  //       UT3I(N,O)   --> Entrance window                                (Mylar)
  //       UXI(N,O)1-6 --> Gas volume (sensitive)                         (Xe/Isobutane)
  //       UT5I(N,O)   --> Pad plane                                      (Cu)
  //       UT6I(N,O)   --> Support structure                              (G10)
  //       UT7I(N,O)   --> FEE + signal lines                             (Cu)
  //       UT8I(N,O)   --> Polyethylene of cooling device                 (PE)
  //       UT9I(N,O)   --> Cooling water                                  (Water)
  //
  //Begin_Html
  /*
    <img src="picts/AliTRDv2.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliTRDv2Tree.gif">
  */
  //End_Html
  
  Float_t xpos, ypos, zpos;
  Int_t   idmat[2];

  const Int_t nparmo = 10;
  const Int_t nparfr =  4;
  const Int_t nparch =  3;
  const Int_t nparic =  4;
  const Int_t nparnc =  4;
  const Int_t nparoc = 11;

  Float_t par_mo[nparmo];
  Float_t par_fr[nparfr];
  Float_t par_ch[nparch];
  Float_t par_ic[nparic];
  Float_t par_nc[nparnc];
  Float_t par_oc[nparoc];

  Int_t icham;

  Int_t *idtmed = gAlice->Idtmed();
  
  AliMC* pMC = AliMC::GetMC();

  //////////////////////////////////////////////////////////////////////////
  //     Definition of Volumes 
  //////////////////////////////////////////////////////////////////////////
  
  // Definition of the mother volume for the TRD (Al) 
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
  pMC->Gsvolu("TRD ", "PGON", idtmed[1301-1], par_mo, nparmo);
  pMC->Gsdvn("UTRS", "TRD ", nsect, 2);

  // The minimal width of a sector in rphi-direction
  Float_t widmi = rmin * TMath::Sin(kPI/nsect);
  // The maximal width of a sector in rphi-direction
  Float_t widma = rmax * TMath::Sin(kPI/nsect);
  // The total thickness of the spaceframe (Al + Air)
  Float_t frame = widmi - (widpl1 / 2);

  // Definition of the inner part of the detector frame (Air) 
  par_fr[0] = widmi - alframe / 2.;
  par_fr[1] = widma - alframe / 2.;
  par_fr[2] = zmax1;
  par_fr[3] = (rmax - rmin) / 2;
  pMC->Gsvolu("UTRI", "TRD1", idtmed[1302-1], par_fr, nparfr); 

  // Some parameter for the chambers 
  Float_t lendifc = (zmax1 - zmax2) / nmodul;
  Float_t heightc = (rmax  - rmin ) / nmodul;
  Float_t widdifc = (widma - widmi) / nmodul;

  // Definition of the chambers 
  Char_t ctagc[5], ctagi[5];
  for (icham = 1; icham <= ncham; ++icham) {

    // Carbon frame of the inner chambers (C) 
    par_ch[0] = widmi + (icham-1) * widdifc - frame;
    par_ch[1] = zleni   / 2.;
    par_ch[2] = heightc / 2.;
    sprintf(ctagc,"UCI%1d",icham);
    pMC->Gsvolu(ctagc, "BOX ", idtmed[1307-1], par_ch, nparch);
    // Inner part of the inner chambers (Air) 
    par_ch[0] -= ccframe;
    par_ch[1] -= ccframe;
    sprintf(ctagc,"UII%1d",icham);
    pMC->Gsvolu(ctagc, "BOX ", idtmed[1302-1], par_ch, nparch);

    // Carbon frame of the neighbouring chambers (C) 
    par_ch[0] = widmi + (icham-1) * widdifc - frame;
    par_ch[1] = zlenn   / 2.;
    par_ch[2] = heightc / 2.;
    sprintf(ctagc,"UCN%1d",icham);
    pMC->Gsvolu(ctagc, "BOX ", idtmed[1307-1], par_ch, nparch);
    // Inner part of the neighbouring chambers (Air) 
    par_ch[0] -= ccframe;
    par_ch[1] -= ccframe;
    sprintf(ctagc,"UIN%1d",icham);
    pMC->Gsvolu(ctagc, "BOX ", idtmed[1302-1], par_ch, nparch);

    // Carbon frame of the outer chambers (C) 
    par_ch[0] = widmi + (icham-1) * widdifc - frame;
    par_ch[1] = (icham - 6) * lendifc / 2. + zleno   / 2.;
    par_ch[2] = heightc / 2.;
    sprintf(ctagc,"UCO%1d",icham);
    pMC->Gsvolu(ctagc, "BOX ", idtmed[1307-1], par_ch, nparch);
    // Inner part of the outer chambers (Air) 
    par_ch[0] -= ccframe;
    par_ch[1] -= ccframe;
    sprintf(ctagc,"UIO%1d",icham);
    pMC->Gsvolu(ctagc, "BOX ", idtmed[1302-1], par_ch, nparch);

  }

  // Definition of the layers in each inner chamber 
  par_ic[0] = -1.;
  par_ic[1] = -1.;
  // G10 layer (radiator layer)
  par_ic[2] = sethick / 2;
  pMC->Gsvolu("UT0I", "BOX ", idtmed[1313-1], par_ic, nparic);
  // CO2 layer (radiator)
  par_ic[2] = rathick / 2;
  pMC->Gsvolu("UT1I", "BOX ", idtmed[1312-1], par_ic, nparic);
  // PE layer (radiator)
  par_ic[2] = pethick / 2;
  pMC->Gsvolu("UT2I", "BOX ", idtmed[1303-1], par_ic, nparic);
  // Mylar layer (entrance window + HV cathode) 
  par_ic[2] = mythick / 2;
  pMC->Gsvolu("UT3I", "BOX ", idtmed[1308-1], par_ic, nparic);
  // Xe/Isobutane layer (gasvolume) 
  par_ic[2] = xethick / 2.;
  for (icham = 1; icham <= 6; ++icham) {
    sprintf(ctagc,"UXI%1d",icham);
    pMC->Gsvolu(ctagc, "BOX ", idtmed[1309-1], par_ic, nparic);
  }
  // Cu layer (pad plane)
  par_ic[2] = cuthick / 2;
  pMC->Gsvolu("UT5I", "BOX ", idtmed[1305-1], par_ic, nparic);
  // G10 layer (support structure)
  par_ic[2] = suthick / 2;
  pMC->Gsvolu("UT6I", "BOX ", idtmed[1313-1], par_ic, nparic);
  // Cu layer (FEE + signal lines)
  par_ic[2] = fethick / 2;
  pMC->Gsvolu("UT7I", "BOX ", idtmed[1305-1], par_ic, nparic);
  // PE layer (cooling devices)
  par_ic[2] = cothick / 2;
  pMC->Gsvolu("UT8I", "BOX ", idtmed[1303-1], par_ic, nparic);
  // Water layer (cooling)
  par_ic[2] = wathick / 2;
  pMC->Gsvolu("UT9I", "BOX ", idtmed[1314-1], par_ic, nparic);

  // Definition of the layers in each neighbouring chamber 
  par_nc[0] = -1.;
  par_nc[1] = -1.;
  // G10 layer (radiator layer)
  par_nc[2] = sethick / 2;
  pMC->Gsvolu("UT0N", "BOX ", idtmed[1313-1], par_nc, nparnc);
  // CO2 layer (radiator)
  par_nc[2] = rathick / 2;
  pMC->Gsvolu("UT1N", "BOX ", idtmed[1312-1], par_nc, nparnc);
  // PE layer (radiator)
  par_nc[2] = pethick / 2;
  pMC->Gsvolu("UT2N", "BOX ", idtmed[1303-1], par_nc, nparnc);
  // Mylar layer (entrance window + HV cathode) 
  par_nc[2] = mythick / 2;
  pMC->Gsvolu("UT3N", "BOX ", idtmed[1308-1], par_nc, nparnc);
  // Xe/Isobutane layer (gasvolume) 
  par_nc[2] = xethick / 2.;
  for (icham = 1; icham <= 6; ++icham) {
    sprintf(ctagc,"UXN%1d",icham);
    pMC->Gsvolu(ctagc, "BOX ", idtmed[1309-1], par_nc, nparnc);
  }
  // Cu layer (pad plane)
  par_nc[2] = cuthick / 2;
  pMC->Gsvolu("UT5N", "BOX ", idtmed[1305-1], par_nc, nparnc);
  // G10 layer (support structure)
  par_nc[2] = suthick / 2;
  pMC->Gsvolu("UT6N", "BOX ", idtmed[1313-1], par_nc, nparnc);
  // Cu layer (FEE + signal lines)
  par_nc[2] = fethick / 2;
  pMC->Gsvolu("UT7N", "BOX ", idtmed[1305-1], par_nc, nparnc);
  // PE layer (cooling devices)
  par_nc[2] = cothick / 2;
  pMC->Gsvolu("UT8N", "BOX ", idtmed[1303-1], par_nc, nparnc);
  // Water layer (cooling)
  par_nc[2] = wathick / 2;
  pMC->Gsvolu("UT9N", "BOX ", idtmed[1314-1], par_nc, nparnc);

  // Definition of the layers in each outer chamber 
  par_oc[0] = -1.;
  par_oc[1] = -1.;
  // G10 layer (radiator layer)
  par_oc[2] = sethick / 2;
  pMC->Gsvolu("UT0O", "BOX ", idtmed[1313-1], par_oc, nparoc);
  // CO2 layer (radiator)
  par_oc[2] = rathick / 2;
  pMC->Gsvolu("UT1O", "BOX ", idtmed[1312-1], par_oc, nparoc);
  // PE layer (radiator)
  par_oc[2] = pethick / 2;
  pMC->Gsvolu("UT2O", "BOX ", idtmed[1303-1], par_oc, nparoc);
  // Mylar layer (entrance window + HV cathode) 
  par_oc[2] = mythick / 2;
  pMC->Gsvolu("UT3O", "BOX ", idtmed[1308-1], par_oc, nparoc);
  // Xe/Isobutane layer (gasvolume) 
  par_oc[2] = xethick / 2.;
  for (icham = 1; icham <= 6; ++icham) {
    sprintf(ctagc,"UXO%1d",icham);
    pMC->Gsvolu(ctagc, "BOX ", idtmed[1309-1], par_oc, nparoc);
  }
  // Cu layer (pad plane)
  par_oc[2] = cuthick / 2;
  pMC->Gsvolu("UT5O", "BOX ", idtmed[1305-1], par_oc, nparoc);
  // G10 layer (support structure)
  par_oc[2] = suthick / 2;
  pMC->Gsvolu("UT6O", "BOX ", idtmed[1313-1], par_oc, nparoc);
  // Cu layer (FEE + signal lines)
  par_oc[2] = fethick / 2;
  pMC->Gsvolu("UT7O", "BOX ", idtmed[1305-1], par_oc, nparoc);
  // PE layer (cooling devices)
  par_oc[2] = cothick / 2;
  pMC->Gsvolu("UT8O", "BOX ", idtmed[1303-1], par_oc, nparoc);
  // Water layer (cooling)
  par_oc[2] = wathick / 2;
  pMC->Gsvolu("UT9O", "BOX ", idtmed[1314-1], par_oc, nparoc);

  //////////////////////////////////////////////////////////////////////////
  //     Positioning of Volumes   
  //////////////////////////////////////////////////////////////////////////

  // The rotation matrices 
  AliMatrix(idmat[0],  90.,  90., 180.,   0.,  90.,   0.);
  AliMatrix(idmat[1],  90.,  90.,   0.,   0.,  90.,   0.);

  // Position of the layers in a chamber 
  pMC->Gspos("UT2I", 1, "UT1I", 0., 0., pezpos, 0, "ONLY");
  pMC->Gspos("UT2N", 1, "UT1N", 0., 0., pezpos, 0, "ONLY");
  pMC->Gspos("UT2O", 1, "UT1O", 0., 0., pezpos, 0, "ONLY");
  for (icham = 1; icham <= ncham; ++icham) {
    // The inner chambers 
    sprintf(ctagi,"UII%1d",icham);
    sprintf(ctagc,"UXI%1d",icham);
    pMC->Gspos("UT9I", icham, ctagi, 0., 0., wazpos, 0, "ONLY");
    pMC->Gspos("UT8I", icham, ctagi, 0., 0., cozpos, 0, "ONLY");
    pMC->Gspos("UT7I", icham, ctagi, 0., 0., fezpos, 0, "ONLY");
    pMC->Gspos("UT6I", icham, ctagi, 0., 0., suzpos, 0, "ONLY");
    pMC->Gspos("UT5I", icham, ctagi, 0., 0., cuzpos, 0, "ONLY");
    pMC->Gspos(ctagc ,     1, ctagi, 0., 0., xezpos, 0, "ONLY");
    pMC->Gspos("UT3I", icham, ctagi, 0., 0., myzpos, 0, "ONLY");
    pMC->Gspos("UT1I", icham, ctagi, 0., 0., razpos, 0, "ONLY");
    pMC->Gspos("UT0I", icham, ctagi, 0., 0., sezpos, 0, "ONLY");
    // The neighbouring chambers 
    sprintf(ctagi,"UIN%1d",icham);
    sprintf(ctagc,"UXN%1d",icham);
    pMC->Gspos("UT9N", icham, ctagi, 0., 0., wazpos, 0, "ONLY");
    pMC->Gspos("UT8N", icham, ctagi, 0., 0., cozpos, 0, "ONLY");
    pMC->Gspos("UT7N", icham, ctagi, 0., 0., fezpos, 0, "ONLY");
    pMC->Gspos("UT6N", icham, ctagi, 0., 0., suzpos, 0, "ONLY");
    pMC->Gspos("UT5N", icham, ctagi, 0., 0., cuzpos, 0, "ONLY");
    pMC->Gspos(ctagc ,     1, ctagi, 0., 0., xezpos, 0, "ONLY");
    pMC->Gspos("UT3N", icham, ctagi, 0., 0., myzpos, 0, "ONLY");
    pMC->Gspos("UT1N", icham, ctagi, 0., 0., razpos, 0, "ONLY");
    pMC->Gspos("UT0N", icham, ctagi, 0., 0., sezpos, 0, "ONLY");
    // The outer chambers 
    sprintf(ctagi,"UIO%1d",icham);
    sprintf(ctagc,"UXO%1d",icham);
    pMC->Gspos("UT9O", icham, ctagi, 0., 0., wazpos, 0, "ONLY");
    pMC->Gspos("UT8O", icham, ctagi, 0., 0., cozpos, 0, "ONLY");
    pMC->Gspos("UT7O", icham, ctagi, 0., 0., fezpos, 0, "ONLY");
    pMC->Gspos("UT6O", icham, ctagi, 0., 0., suzpos, 0, "ONLY");
    pMC->Gspos("UT5O", icham, ctagi, 0., 0., cuzpos, 0, "ONLY");
    pMC->Gspos(ctagc ,     1, ctagi, 0., 0., xezpos, 0, "ONLY");
    pMC->Gspos("UT3O", icham, ctagi, 0., 0., myzpos, 0, "ONLY");
    pMC->Gspos("UT1O", icham, ctagi, 0., 0., razpos, 0, "ONLY");
    pMC->Gspos("UT0O", icham, ctagi, 0., 0., sezpos, 0, "ONLY");
  }

  // Position of the inner part of the chambers in the carbon-frames 
  for (icham = 1; icham <= ncham; ++icham) {
    xpos = 0.;
    ypos = 0.;
    zpos = 0.;
    // The inner chambers 
    sprintf(ctagi,"UII%1d",icham);
    sprintf(ctagc,"UCI%1d",icham);
    pMC->Gspos(ctagi, 1, ctagc, xpos, ypos, zpos, 0, "ONLY");
    // The neighbouring chambers 
    sprintf(ctagi,"UIN%1d",icham);
    sprintf(ctagc,"UCN%1d",icham);
    pMC->Gspos(ctagi, 1, ctagc, xpos, ypos, zpos, 0, "ONLY");
    // The outer chambers 
    sprintf(ctagi,"UIO%1d",icham);
    sprintf(ctagc,"UCO%1d",icham);
    pMC->Gspos(ctagi, 1, ctagc, xpos, ypos, zpos, 0, "ONLY");
  }

  // Position of the chambers in the full TRD-setup 
  for (icham = 1; icham <= ncham; ++icham) {
    // The inner chambers
    xpos = 0.;
    ypos = 0.;
    zpos = (icham-0.5) * heightc - (rmax - rmin) / 2;
    sprintf(ctagc,"UCI%1d",icham);
    pMC->Gspos(ctagc, 1, "UTRI", xpos, ypos, zpos, 0, "ONLY");
    // The neighbouring chambers
    xpos = 0.;
    ypos = (zleni + zlenn) / 2.;
    zpos = (icham-0.5) * heightc - (rmax - rmin) / 2;
    sprintf(ctagc,"UCN%1d",icham);
    pMC->Gspos(ctagc, 1, "UTRI", xpos, ypos, zpos, 0, "ONLY");
    ypos = -ypos;
    sprintf(ctagc,"UCN%1d",icham);
    pMC->Gspos(ctagc, 2, "UTRI", xpos, ypos, zpos, 0, "ONLY");
    // The outer chambers 
    xpos = 0.;
    ypos = (zleni / 2. + zlenn + zmax2 + (icham-1) * lendifc) / 2.;
    zpos = (icham-0.5) * heightc - (rmax-rmin)/2;
    sprintf(ctagc,"UCO%1d",icham);
    pMC->Gspos(ctagc, 1, "UTRI", xpos, ypos, zpos, 0, "ONLY");
    ypos = -ypos;
    sprintf(ctagc,"UCO%1d",icham);
    pMC->Gspos(ctagc, 2, "UTRI", xpos, ypos, zpos, 0, "ONLY");
  }

  // Position of the inner part of the detector frame
  xpos = (rmax + rmin) / 2;
  ypos = 0.;
  zpos = 0.;
  pMC->Gspos("UTRI", 1, "UTRS", xpos, ypos, zpos, idmat[0], "ONLY");

  // Position of the TRD mother volume in the ALICE experiment 
  xpos = 0.;
  ypos = 0.;
  zpos = 0.;
  pMC->Gspos("TRD ", 1, "ALIC", xpos, ypos, zpos,        0, "ONLY");

}

//_____________________________________________________________________________
void AliTRDv2::DrawModule()
{
  //
  // Draw a shaded view of the Transition Radiation Detector version 1
  //

  AliMC* pMC = AliMC::GetMC();
  
  // Set everything unseen
  pMC->Gsatt("*", "seen", -1);
   
  // Set ALIC mother transparent
  pMC->Gsatt("ALIC","SEEN",0);
  
  // Set the volumes visible
  pMC->Gsatt("TRD ","SEEN",0);
  pMC->Gsatt("UTRS","SEEN",0);
  pMC->Gsatt("UTRI","SEEN",0);
  Char_t ctag[5];
  for (Int_t icham = 0; icham < ncham; ++icham) {
    sprintf(ctag,"UCI%1d",icham+1);
    pMC->Gsatt(ctag,"SEEN",0);
    sprintf(ctag,"UCN%1d",icham+1);
    pMC->Gsatt(ctag,"SEEN",0);
    sprintf(ctag,"UCO%1d",icham+1);
    pMC->Gsatt(ctag,"SEEN",0);
    sprintf(ctag,"UII%1d",icham+1);
    pMC->Gsatt(ctag,"SEEN",0);
    sprintf(ctag,"UIN%1d",icham+1);
    pMC->Gsatt(ctag,"SEEN",0);
    sprintf(ctag,"UIO%1d",icham+1);
    pMC->Gsatt(ctag,"SEEN",0);
    sprintf(ctag,"UXI%1d",icham+1);
    pMC->Gsatt(ctag,"SEEN",1);
    sprintf(ctag,"UXN%1d",icham+1);
    pMC->Gsatt(ctag,"SEEN",1);
    sprintf(ctag,"UXO%1d",icham+1);
    pMC->Gsatt(ctag,"SEEN",1);
  }
  pMC->Gsatt("UT1I","SEEN",1);
  pMC->Gsatt("UT1N","SEEN",1);
  pMC->Gsatt("UT1O","SEEN",1);

  pMC->Gdopt("hide", "on");
  pMC->Gdopt("shad", "on");
  pMC->Gsatt("*", "fill", 7);
  pMC->SetClipBox(".");
  pMC->SetClipBox("*", 0, 2000, -2000, 2000, -2000, 2000);
  pMC->DefaultRange();
  pMC->Gdraw("alic", 40, 30, 0, 12, 9.4, .021, .021);
  pMC->Gdhead(1111, "Transition Radiation Detector Version 2");
  pMC->Gdman(18, 4, "MAN");
  pMC->Gdopt("hide", "off");
}

//_____________________________________________________________________________
void AliTRDv2::CreateMaterials()
{
  //
  // Create materials for the Transition Radiation Detector version 2
  //
  AliTRD::CreateMaterials();
}

//_____________________________________________________________________________
void AliTRDv2::Init() 
{
  //
  // Initialise Transition Radiation Detector after geometry has been built
  //

  // First ionization potential (eV) for the gas mixture (90% Xe + 10% CO2)
  const Float_t kPoti = 12.1;
  // Maximum energy (50 keV);
  const Float_t kEend = 50000.0;

  AliTRD::Init();

  AliMC* pMC = AliMC::GetMC();

  // Get the sensitive volumes
  Char_t ctag[5];
  for (Int_t icham = 0; icham < ncham; ++icham) {
    sprintf(ctag,"UXI%1d",icham+1);
    fIdSensI[icham] = pMC->VolId(ctag);
    sprintf(ctag,"UXN%1d",icham+1);
    fIdSensN[icham] = pMC->VolId(ctag);
    sprintf(ctag,"UXO%1d",icham+1);
    fIdSensO[icham] = pMC->VolId(ctag);
  }

  Float_t Poti = TMath::Log(kPoti);
  Float_t Eend = TMath::Log(kEend);

  // Ermilova distribution for the delta-ray spectrum
  fDeltaE  = new TF1("deltae",Ermilova,Poti,Eend,0);

}

//_____________________________________________________________________________
void AliTRDv2::StepManager()
{
  //
  // Called at every step in the Transition Radiation Detector version 2
  //

  Int_t          idSens, icSens, id;
  Int_t          iPla, iCha, iSec;
  Int_t          iOut;
  Int_t          vol[3]; 
  Int_t          iPid;

  const Double_t kBig = 1.0E+12;

  Float_t        hits[4];
  Float_t        mom[4];
  Float_t        random[1];
  Float_t        charge;
  Float_t        aMass;

  Double_t       pTot;
  Double_t       qTot;
  Double_t       eDelta;
  Double_t       betaGamma, pp;

  TClonesArray  &lhits = *fHits;

  AliMC* pMC = AliMC::GetMC();

  // Ionization energy
  const Float_t kWion    = 22.04;
  // Maximum energy for e+ e- g for the step-size calculation
  const Float_t kPTotMax = 0.002;
  // Plateau value of the energy-loss for electron in xenon
  // taken from: Allison + Comb, Ann. Rev. Nucl. Sci. (1980), 30, 253
  //const Double_t kPlateau = 1.70;
  // the averaged value (26/3/99)
  const Float_t kPlateau = 1.55;
  // dN1/dx|min for the gas mixture (90% Xe + 10% CO2)
  const Float_t kPrim    = 48.0;
  // First ionization potential (eV) for the gas mixture (90% Xe + 10% CO2)
  const Float_t kPoti    = 12.1;

  // Set the maximum step size to a very large number for all 
  // neutral particles and those outside the driftvolume
  pMC->SetMaxStep(kBig); 

  // Use only charged tracks 
  if (( pMC->TrackCharge()   ) &&
      (!pMC->TrackStop()     ) && 
      (!pMC->TrackDisappear())) {

    // Find the sensitive volume
    idSens = pMC->CurrentVol(0,icSens);
    iPla   = 0;
    iOut   = 0;
    for (Int_t icham = 0; icham < ncham; ++icham) {
      if (idSens == fIdSensI[icham]) {
        iOut = 0;
        iPla = icham + 1;
      }
      if (idSens == fIdSensN[icham]) {
        iOut = 1;
        iPla = icham + 1;
      }
      if (idSens == fIdSensO[icham]) {
        iOut = 2;
        iPla = icham + 1;
      }
    }

    // Inside a sensitive volume?
    if (iPla) {

      // Calculate the energy of the delta-electrons
      eDelta = TMath::Exp(fDeltaE->GetRandom()) - kPoti;
      eDelta = TMath::Max(eDelta,0.0);

      // The number of secondary electrons created
      qTot = (Double_t) ((Int_t) (eDelta / kWion) + 1);

      // The sector number
      id = pMC->CurrentVolOff(4,0,iSec);

      // The chamber number
      //   1: outer left
      //   2: neighbouring left
      //   3: inner
      //   4: neighbouring right
      //   5: outer right
      id = pMC->CurrentVolOff(2,0,iCha);
      if (iCha == 1) 
        iCha = 3 + iOut;
      else
        iCha = 3 - iOut;

      vol[0]  = iSec;
      vol[1]  = iCha;
      vol[2]  = iPla;

      // Check on selected volumes
      Int_t addthishit = 1;
      if (fSensSelect) {
        if ((fSensPlane)   && (vol[2] != fSensPlane  )) addthishit = 0;
        if ((fSensChamber) && (vol[1] != fSensChamber)) addthishit = 0;
        if ((fSensSector)  && (vol[0] != fSensSector )) addthishit = 0;
      }

      if (addthishit) {

        // Add this hit
        pMC->TrackPosition(hits);
        hits[3] = qTot;
        new(lhits[fNhits++]) AliTRDhit(fIshunt,gAlice->CurrentTrack(),vol,hits);

        // The energy loss according to Bethe Bloch
        pMC->TrackMomentum(mom);
        pTot = mom[3];
        iPid = pMC->TrackPid();
        if ( (iPid >  3) ||
	    ((iPid <= 3) && (pTot < kPTotMax))) {
          aMass     = pMC->TrackMass();
          betaGamma = pTot / aMass;
          pp        = kPrim * BetheBloch(betaGamma);
	  // Take charge > 1 into account
          charge = pMC->TrackCharge();
          if (TMath::Abs(charge) > 1) pp = pp * charge*charge;
        }
        // Electrons above 20 Mev/c are at the plateau
        else {
          pp = kPrim * kPlateau;
        }
      
        // Calculate the maximum step size for the next tracking step
        if (pp > 0) {
          do 
            pMC->Rndm(random,1);
          while ((random[0] == 1.) || (random[0] == 0.));
          pMC->SetMaxStep( - TMath::Log(random[0]) / pp);
	}

      }
      else {
        // set step size to maximal value
        pMC->SetMaxStep(kBig); 
      }

    }

  }

}

//_____________________________________________________________________________
Double_t AliTRDv2::BetheBloch(Double_t bg) 
{
  //
  // Parametrization of the Bethe-Bloch-curve
  // The parametrization is the same as for the TPC and is taken from Lehrhaus.
  //

  // The parameters have been adjusted to Xe-data found in:
  // Allison & Cobb, Ann. Rev. Nucl. Sci. (1980), 30, 253
  //const Double_t kP1 = 0.76176E-1;
  //const Double_t kP2 = 10.632;
  //const Double_t kP3 = 3.17983E-6;
  //const Double_t kP4 = 1.8631;
  //const Double_t kP5 = 1.9479;

  // This parameters have been adjusted to averaged values from GEANT
  const Double_t kP1 = 7.17960e-02;
  const Double_t kP2 = 8.54196;
  const Double_t kP3 = 1.38065e-06;
  const Double_t kP4 = 5.30972;
  const Double_t kP5 = 2.83798;

  if (bg > 0) {
    Double_t yy = bg / TMath::Sqrt(1. + bg*bg);
    Double_t aa = TMath::Power(yy,kP4);
    Double_t bb = TMath::Power((1./bg),kP5);
             bb = TMath::Log(kP3 + bb);
    return ((kP2 - aa - bb)*kP1 / aa);
  }
  else
    return 0;

}

//_____________________________________________________________________________
Double_t Ermilova(Double_t *x, Double_t *)
{
  //
  // Calculates the delta-ray energy distribution according to Ermilova
  // Logarithmic scale !
  //

  Double_t energy;
  Double_t dpos;
  Double_t dnde;

  Int_t    pos1, pos2;

  const Int_t nV = 31;

  Float_t vxe[nV] = { 2.3026, 2.9957, 3.4012, 3.6889, 3.9120  
                    , 4.0943, 4.2485, 4.3820, 4.4998, 4.6052
                    , 4.7005, 5.0752, 5.2983, 5.7038, 5.9915
                    , 6.2146, 6.5221, 6.9078, 7.3132, 7.6009
                    , 8.0064, 8.5172, 8.6995, 8.9872, 9.2103
                    , 9.4727, 9.9035,10.3735,10.5966,10.8198
                    ,11.5129 };

  Float_t vye[nV] = { 80.0  , 31.0  , 23.3  , 21.1  , 21.0
                    , 20.9  , 20.8  , 20.0  , 16.0  , 11.0
                    ,  8.0  ,  6.0  ,  5.2  ,  4.6  ,  4.0
                    ,  3.5  ,  3.0  ,  1.4  ,  0.67 ,  0.44
                    ,  0.3  ,  0.18 ,  0.12 ,  0.08 ,  0.056
                    ,  0.04 ,  0.023,  0.015,  0.011,  0.01
		    ,  0.004 };

  energy = x[0];

  // Find the position 
  pos1 = pos2 = 0;
  dpos = 0;
  do {
    dpos = energy - vxe[pos2++];
  } 
  while (dpos > 0);
  pos2--; 
  if (pos2 > nV) pos2 = nV;
  pos1 = pos2 - 1;

  // Differentiate between the sampling points
  dnde = (vye[pos1] - vye[pos2]) / (vxe[pos2] - vxe[pos1]);

  return dnde;

}
