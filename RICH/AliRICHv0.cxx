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
Revision 1.3  1999/09/29 09:24:29  fca
Introduction of the Copyright and cvs Log

*/

/////////////////////////////////////////////////////////
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
#include "TGeant3.h"

ClassImp(AliRICHv0)
    
//___________________________________________
AliRICHv0::AliRICHv0() : AliRICH()
{
    fChambers = 0;
}

//___________________________________________
AliRICHv0::AliRICHv0(const char *name, const char *title)
    : AliRICH(name,title)
{
    
    fChambers = new TObjArray(7);
    for (Int_t i=0; i<7; i++) {
    
	(*fChambers)[i] = new AliRICHchamber();  
	
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
    
    
    Int_t *idtmed = fIdtmed->GetArray()-999;
    
    Int_t i;
    Float_t zs;
    Int_t idrotm[1099];
    Float_t par[3];
    
    // --- Define the RICH detector 
    //     External aluminium box 
    par[0] = 71.1;
    par[1] = 11.5;
    par[2] = 73.15;
    gMC->Gsvolu("RICH", "BOX ", idtmed[1009], par, 3);
    
    //     Sensitive part of the whole RICH 
    par[0] = 64.8;
    par[1] = 11.5;
    par[2] = 66.55;
    gMC->Gsvolu("SRIC", "BOX ", idtmed[1000], par, 3);
    
    //     Honeycomb 
    par[0] = 63.1;
    par[1] = .188;
    par[2] = 66.55;
    gMC->Gsvolu("HONE", "BOX ", idtmed[1001], par, 3);
    
    //     Aluminium sheet 
    par[0] = 63.1;
    par[1] = .025;
    par[2] = 66.55;
    gMC->Gsvolu("ALUM", "BOX ", idtmed[1009], par, 3);
    
    //     Quartz 
    par[0] = 63.1;
    par[1] = .25;
    par[2] = 65.5;
    gMC->Gsvolu("QUAR", "BOX ", idtmed[1002], par, 3);
    
    //     Spacers (cylinders) 
    par[0] = 0.;
    par[1] = .5;
    par[2] = .5;
    gMC->Gsvolu("SPAC", "TUBE", idtmed[1002], par, 3);
    
    //     Opaque quartz 
    par[0] = 61.95;
    par[1] = .2;
    par[2] = 66.5;
    gMC->Gsvolu("OQUA", "BOX ", idtmed[1007], par, 3);
  
    //     Frame of opaque quartz 
    par[0] = 20.65;
    par[1] = .5;
    par[2] = 66.5;
    gMC->Gsvolu("OQUF", "BOX ", idtmed[1007], par, 3);
    
    //     Little bar of opaque quartz 
    par[0] = 63.1;
    par[1] = .25;
    par[2] = .275;
    gMC->Gsvolu("BARR", "BOX ", idtmed[1007], par, 3);
    
    //     Freon 
    par[0] = 20.15;
    par[1] = .5;
    par[2] = 65.5;
    gMC->Gsvolu("FREO", "BOX ", idtmed[1003], par, 3);
    
    //     Methane 
    par[0] = 64.8;
    par[1] = 5.;
    par[2] = 64.8;
    gMC->Gsvolu("META", "BOX ", idtmed[1004], par, 3);
    
    //     Methane gap 
    par[0] = 64.8;
    par[1] = .2;
    par[2] = 64.8;
    gMC->Gsvolu("GAP ", "BOX ", idtmed[1008], par, 3);
    
    //     CsI photocathode 
    par[0] = 64.8;
    par[1] = .25;
    par[2] = 64.8;
    gMC->Gsvolu("CSI ", "BOX ", idtmed[1005], par, 3);
    
    //     Anode grid 
    par[0] = 0.;
    par[1] = .0025;
    par[2] = 20.;
    gMC->Gsvolu("GRID", "TUBE", idtmed[1006], par, 3);
    
    // --- Places the detectors defined with GSVOLU 
    //     Place material inside RICH 
    gMC->Gspos("SRIC", 1, "RICH", 0., 0., 0., 0, "ONLY");
    
    gMC->Gspos("ALUM", 1, "SRIC", 0., -6.075, 0., 0, "ONLY");
    gMC->Gspos("HONE", 1, "SRIC", 0., -5.862, 0., 0, "ONLY");
    gMC->Gspos("ALUM", 2, "SRIC", 0., -5.649, 0., 0, "ONLY");
    gMC->Gspos("OQUA", 1, "SRIC", 0., -5.424, 0., 0, "ONLY");
    
    AliMatrix(idrotm[1019], 0., 0., 90., 0., 90., 90.);
    
    for (i = 1; i <= 9; ++i) {
	zs = (5 - i) * 14.4;
	gMC->Gspos("SPAC", i, "FREO", 6.7, 0., zs, idrotm[1019], "ONLY");
    }
    for (i = 10; i <= 18; ++i) {
	zs = (14 - i) * 14.4;
	gMC->Gspos("SPAC", i, "FREO", -6.7, 0., zs, idrotm[1019], "ONLY");
    }
    
    gMC->Gspos("FREO", 1, "OQUF", 0., 0., 0., 0, "ONLY");
    gMC->Gspos("OQUF", 1, "SRIC", 41.3, -4.724, 0., 0, "ONLY");
    gMC->Gspos("OQUF", 2, "SRIC", 0., -4.724, 0., 0, "ONLY");
    gMC->Gspos("OQUF", 3, "SRIC", -41.3, -4.724, 0., 0, "ONLY");
    gMC->Gspos("BARR", 1, "QUAR", 0., 0., -21.65, 0, "ONLY");
    gMC->Gspos("BARR", 2, "QUAR", 0., 0., 21.65, 0, "ONLY");
    gMC->Gspos("QUAR", 1, "SRIC", 0., -3.974, 0., 0, "ONLY");
    gMC->Gspos("GAP ", 1, "META", 0., 4.8, 0., 0, "ONLY");
    gMC->Gspos("META", 1, "SRIC", 0., 1.276, 0., 0, "ONLY");
    gMC->Gspos("CSI ", 1, "SRIC", 0., 6.526, 0., 0, "ONLY");
    
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
				179.0987,179.0987,179.0987,35.7,12.54,5.92,4.92,3.86,1.42,.336,.134,0. };
    Float_t absco_quarz[14] = { 20.126,16.27,13.49,11.728,9.224,8.38,7.44,7.17,
				6.324,4.483,1.6,.323,.073,0. };
    Float_t absco_quarzo[14] = { 1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,
				 1e-5,1e-5,1e-5,1e-5,1e-5 };
    Float_t absco_csi[14] = { 1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,
			      1e-4,1e-4,1e-4,1e-4 };
    Float_t absco_methane[14] = { 1e6,1e6,1e6,1e6,1e6,1e6,1e6,1e6,1e6,1e6,1e6,
				  1e6,1e6,1e6 };
    Float_t absco_gri[14] = { 1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,
			      1e-4,1e-4,1e-4,1e-4 };
    Float_t effic_all[14] = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };
    Float_t effic_csi[14] = { 4.74e-4,.00438,.009,.0182,.0282,.0653,.1141,.163,
			      .2101,.2554,.293,.376,.3861,.418 };
    Float_t effic_gri[14] = { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. };
    
    Float_t afre[2], agri, amet[2], aqua[2], ahon, zfre[2], zgri, zhon, 
    zmet[2], zqua[2];
    Int_t nlmatfre;
    Float_t densquao;
    Int_t nlmatmet, nlmatqua;
    Float_t wmatquao[2], rindex_freon[14];
    Int_t i;
    Float_t aquao[2], epsil, stmin, zquao[2];
    Int_t nlmatquao;
    Float_t radlal, densal, tmaxfd, deemax, stemax;
    Float_t aal, zal, radlgri, densfre, radlhon, densgri, denshon,densqua, densmet, wmatfre[2], wmatmet[2], wmatqua[2];
    
    Int_t *idtmed = fIdtmed->GetArray()-999;
    
    TGeant3 *geant3 = (TGeant3*) gMC;
    
    // --- Photon energy (GeV) 
    // --- Refraction indexes 
    for (i = 0; i < 14; ++i) {
	rindex_freon[i] = ppckov[i] * .01095 * 1e9 + 1.2177;
    }
    // need to be changed 
    
    // --- Absorbtion lenghts (in cm) 
    //      DATA ABSCO_QUARZ / 
    //     &    5 * 1000000., 
    //     &    29.85,    7.34,     4.134,    1.273,    0.722, 
    //     &    0.365, 0.365, 0.365, 0.  / 
    // need to be changed 
    
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
    
    
    //     Switch on delta-ray production in the methane and freon gaps 
    
    gMC->Gstpar(idtmed[1002], "LOSS", 1.);
    gMC->Gstpar(idtmed[1003], "LOSS", 1.);
    gMC->Gstpar(idtmed[1004], "LOSS", 1.);
    gMC->Gstpar(idtmed[1008], "LOSS", 1.);
    gMC->Gstpar(idtmed[1005], "LOSS", 1.);
    gMC->Gstpar(idtmed[1002], "HADR", 1.);
    gMC->Gstpar(idtmed[1003], "HADR", 1.);
    gMC->Gstpar(idtmed[1004], "HADR", 1.);
    gMC->Gstpar(idtmed[1008], "HADR", 1.);
    gMC->Gstpar(idtmed[1005], "HADR", 1.);
    gMC->Gstpar(idtmed[1002], "DCAY", 1.);
    gMC->Gstpar(idtmed[1003], "DCAY", 1.);
    gMC->Gstpar(idtmed[1004], "DCAY", 1.);
    gMC->Gstpar(idtmed[1008], "DCAY", 1.);
    gMC->Gstpar(idtmed[1005], "DCAY", 1.);
    geant3->Gsckov(idtmed[1000], 14, ppckov, absco_methane, effic_all, rindex_methane);
    geant3->Gsckov(idtmed[1001], 14, ppckov, absco_methane, effic_all, rindex_methane);
    geant3->Gsckov(idtmed[1002], 14, ppckov, absco_quarz, effic_all,rindex_quarz);
    geant3->Gsckov(idtmed[1003], 14, ppckov, absco_freon, effic_all,rindex_freon);
    geant3->Gsckov(idtmed[1004], 14, ppckov, absco_methane, effic_all, rindex_methane);
    geant3->Gsckov(idtmed[1005], 14, ppckov, absco_csi, effic_csi, rindex_methane);
    geant3->Gsckov(idtmed[1006], 14, ppckov, absco_gri, effic_gri, rindex_gri);
    geant3->Gsckov(idtmed[1007], 14, ppckov, absco_quarzo, effic_all, rindex_quarzo);
    geant3->Gsckov(idtmed[1008], 14, ppckov, absco_methane, effic_all, rindex_methane);
    geant3->Gsckov(idtmed[1009], 14, ppckov, absco_gri, effic_gri, rindex_gri);
}

//___________________________________________

void AliRICHv0::Init()
{
    printf("\n\n\n Start Init for version 0 - CPC chamber type \n\n\n");
    
    // 
    // Initialize Tracking Chambers
    //
    for (Int_t i=0; i<7; i++) {
	( (AliRICHchamber*) (*fChambers)[i])->Init();
    }
    
    //
    // Set the chamber (sensitive region) GEANT identifier
    
    ((AliRICHchamber*)(*fChambers)[0])->SetGid(1);
    ((AliRICHchamber*)(*fChambers)[1])->SetGid(2);
    ((AliRICHchamber*)(*fChambers)[2])->SetGid(3);
    ((AliRICHchamber*)(*fChambers)[3])->SetGid(4);
    ((AliRICHchamber*)(*fChambers)[4])->SetGid(5);
    ((AliRICHchamber*)(*fChambers)[5])->SetGid(6);
    ((AliRICHchamber*)(*fChambers)[6])->SetGid(7);
    
    printf("\n\n\n Finished Init for version 0 - CPC chamber type\n\n\n");
}

//___________________________________________
void AliRICHv0::StepManager()
{
    Int_t          copy, id;
    static Int_t   idvol;
    static Int_t   vol[2];
    Int_t          ipart;
    static Float_t hits[10];
    TLorentzVector Position;
    TLorentzVector Momentum;
    Float_t        pos[3];
    Float_t        mom[4];
    Float_t        Localpos[3];
    Float_t        Localmom[4];
    Float_t        Localtheta,Localphi;
    Float_t        theta,phi;
    Float_t        destep, step;
    static Float_t eloss, xhit, yhit, tlength;
    const  Float_t big=1.e10;
    
    TClonesArray &lhits = *fHits;
    TClonesArray &lcerenkovs = *fCerenkovs;
    
    // Only gas gap inside chamber
    // Tag chambers and record hits when track enters 
    
    idvol=-1;
    id=gMC->CurrentVolID(copy);
    Float_t cherenkov_loss=0.00001;
    
    //        Treat photons produced in Freon and Quartz 
    if (gMC->TrackPid() == 50000050 ) {
	if (gMC->IsTrackEntering()){
	    if (gMC->VolId("FREO")==gMC->CurrentVolID(copy) || gMC->VolId("QUAR")==gMC->CurrentVolID(copy)){
		//printf("GOT ONE! Type:%d    \n",gMC->TrackPid());
	    }
	}
    }
    
  
    if (gMC->TrackPid() == 50000050 ) {
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

	      ((AliRICHchamber*) (*fChambers)[idvol])
		->SigGenInit(Localpos[0], Localpos[2], Localpos[1]);
	      if(idvol<7) {
		hits[0] = 50000050;               // particle type
		hits[1] = pos[0];                 // X-position for hit
		hits[2] = pos[1];                 // Y-position for hit
		hits[3] = pos[2];                 // Z-position for hit
		hits[4] = theta;                      // theta angle of incidence
		hits[5] = phi;                      // phi angle of incidence 
		hits[8] = (Float_t) fNclusters;   // first padhit
		hits[9] = -1;                     // last pad hit
		
		MakePadHits(Localpos[0],Localpos[2],cherenkov_loss,idvol,cerenkov);
		if (fNclusters > (Int_t)hits[8]) {
		  hits[8]= hits[8]+1;
		  hits[9]= (Float_t) fNclusters;
		}
		  
		AddHit(gAlice->CurrentTrack(),vol,hits);
		new(lcerenkovs[fNcerenkovs++]) AliRICHCerenkov(fIshunt,gAlice->CurrentTrack(),vol,hits);
		
	      }
	    }
	  }
	
	
// 
// treat charged particles
    } else if (gMC->TrackCharge())
    {
	
//If MIP
	if (gMC->VolId("GAP ")== gMC->CurrentVolID(copy)) {
// Get current particle id (ipart), track position (pos)  and momentum (mom) 

	  gMC->CurrentVolOffID(3,copy);
	  vol[0]=copy;
	  idvol=vol[0]-1;
		
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
		gMC->SetMaxStep(fMaxStepGas);
		Double_t tc = mom[0]*mom[0]+mom[1]*mom[1];
		Double_t rt = TMath::Sqrt(tc);
		theta   = Float_t(TMath::ATan2(rt,Double_t(mom[2])))*kRaddeg;
		phi     = Float_t(TMath::ATan2(Double_t(mom[1]),Double_t(mom[0])))*kRaddeg;
		
		Localtheta   = Float_t(TMath::ATan2(rt,Double_t(Localmom[2])))*kRaddeg;                       
		Localphi     = Float_t(TMath::ATan2(Double_t(Localmom[1]),Double_t(Localmom[0])))*kRaddeg;    
		
		hits[0] = Float_t(ipart);         // particle type
		hits[1] = pos[0];                 // X-position for hit
		hits[2] = pos[1];                 // Y-position for hit
		hits[3] = pos[2];                 // Z-position for hit
		hits[4] = theta;                  // theta angle of incidence
		hits[5] = phi;                    // phi angle of incidence 
		hits[8] = (Float_t) fNclusters;   // first padhit
		hits[9] = -1;                     // last pad hit
		// phi angle of incidence
		tlength = 0;
		eloss   = 0;
		
		Chamber(idvol).LocaltoGlobal(Localpos,hits+1);
		
		//To make chamber coordinates x-y had to pass LocalPos[0], LocalPos[2]
		xhit    = Localpos[0];
		yhit    = Localpos[2];
		// Only if not trigger chamber
		if(idvol<7) {
		    //
		    //  Initialize hit position (cursor) in the segmentation model 
		    ((AliRICHchamber*) (*fChambers)[idvol])
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
		    if (eloss > 0) MakePadHits(xhit,yhit,eloss,idvol,mip);
		}
		
		hits[6]=tlength;
		hits[7]=eloss;
		if (fNclusters > (Int_t)hits[8]) {
		    hits[8]= hits[8]+1;
		    hits[9]= (Float_t) fNclusters;
		}

		new(lhits[fNhits++]) 
		    AliRICHhit(fIshunt,gAlice->CurrentTrack(),vol,hits);
		eloss = 0; 
		//
		// Check additional signal generation conditions 
		// defined by the segmentation
		// model (boundary crossing conditions) 
	    } else if 
		(((AliRICHchamber*) (*fChambers)[idvol])
		 ->SigGenCond(Localpos[0], Localpos[2], Localpos[1]))
	    {
		((AliRICHchamber*) (*fChambers)[idvol])
		    ->SigGenInit(Localpos[0], Localpos[2], Localpos[1]);
		if (eloss > 0) MakePadHits(xhit,yhit,eloss,idvol,mip);
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
}

  
//___________________________________________
void AliRICH::MakePadHits(Float_t xhit,Float_t yhit,Float_t eloss, Int_t idvol, Response_t res)
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
    
    ((AliRICHchamber*) (*fChambers)[idvol])->DisIntegration(eloss, xhit, yhit, nnew, newclust, res);
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
	    
	    AddCluster(clhits);
	}
    }
}

ClassImp(AliRICHchamber)	
    
    AliRICHchamber::AliRICHchamber() 
{
    fSegmentation = new TObjArray(2);
    fResponse= new TObjArray(2);
    fnsec=1;
    frMin=0.1;
    frMax=140;
}

//  
//  Get reference to response model
AliRICHresponse* AliRICHchamber::GetResponseModel(Response_t res)
{
    if (res==mip) {
	return (AliRICHresponse*) (*fResponse)[0];
    } else if (res==cerenkov) {
	return (AliRICHresponse*) (*fResponse)[1];
    }
    return (AliRICHresponse*) (*fResponse)[0];
}

// Configure response model
void   AliRICHchamber::ResponseModel(Response_t res, AliRICHresponse* thisResponse)
{
    
    if (res==mip) {
	(*fResponse)[0]=thisResponse;
    } else if (res==cerenkov) {
	(*fResponse)[1]=thisResponse;
    }
}

void AliRICHchamber::Init()
{
    
    ((AliRICHsegmentation *) (*fSegmentation)[0])->Init(this);
    if (fnsec==2) {
	((AliRICHsegmentation *) (*fSegmentation)[1])->Init(this);
    }
}

void AliRICHchamber::LocaltoGlobal(Float_t pos[3],Float_t Localpos[3])
{

    Double_t *fMatrix;
    fMatrix =  fChamberMatrix->GetMatrix();
    Localpos[0]=pos[0]*fMatrix[0]+pos[1]*fMatrix[3]+pos[2]*fMatrix[6];
    Localpos[1]=pos[0]*fMatrix[1]+pos[1]*fMatrix[4]+pos[2]*fMatrix[7];
    Localpos[2]=pos[0]*fMatrix[2]+pos[1]*fMatrix[5]+pos[2]*fMatrix[8];
    Localpos[0]+=fChamberTrans[0];
    Localpos[1]+=fChamberTrans[1];
    Localpos[2]+=fChamberTrans[2];
}


void AliRICHchamber::DisIntegration(Float_t eloss, Float_t xhit, Float_t yhit,
				    Int_t& nnew,Float_t newclust[6][500],Response_t res) 
{
//    
//  Generates pad hits (simulated cluster) 
//  using the segmentation and the response model
    
    Float_t dx, dy;
    Float_t local[3];
    Float_t source[3];
    Int_t Nfp=0;
    //    Float_t dummy=0;
    
    
    //
    // Width of the integration area
    //
    dx=((AliRICHresponse*) (*fResponse)[0])->Nsigma()*((AliRICHresponse*) (*fResponse)[0])->ChwX();
    dy=((AliRICHresponse*) (*fResponse)[0])->Nsigma()*((AliRICHresponse*) (*fResponse)[0])->ChwY();
    //
    // Get pulse height from energy loss
    Float_t qtot;
    
    if (res==mip) {
	qtot = ((AliRICHresponse*) (*fResponse)[0])->IntPH(eloss);
	
	local[0]=xhit;
	//Z position of the wires relative to the RICH mother volume
	local[1]=6.026;
	local[2]=yhit;
	//Generate feedback photons
	Nfp  = ((AliRICHresponse*) (*fResponse)[0])->FeedBackPhotons(source,qtot);
	//printf("\nFeedbacks (Mip)     :%d",Nfp);
    } else if (res==cerenkov) {
	qtot = ((AliRICHresponse*) (*fResponse)[1])->IntPH();
	local[0]=xhit;
	local[1]=6.026;
	local[2]=yhit;
	Nfp  = ((AliRICHresponse*) (*fResponse)[1])->FeedBackPhotons(source,qtot);
	//printf("\nFeedbacks (Cerenkov):%d",Nfp);
    }

    //
    // Loop Over Pads
    
    Float_t qcheck=0, qp;
    nnew=0;
    for (Int_t i=1; i<=fnsec; i++) {
	qcheck=0;
	AliRICHsegmentation * segmentation=(AliRICHsegmentation *) (*fSegmentation)[i-1];
	for (segmentation->FirstPad(xhit, yhit, dx, dy); 
	     segmentation->MorePads(); 
	     segmentation->NextPad()) 
	{
	    if (res==mip) {
		qp= ((AliRICHresponse*) (*fResponse)[0])->IntXY(segmentation);
	    }
	    if (res==cerenkov) {
		qp= ((AliRICHresponse*) (*fResponse)[0])->IntXY(segmentation);
	    }
	    
	    qp= TMath::Abs(qp);

	    if (qp > 1.e-4) {
		qcheck+=qp;
		//
		// --- store signal information
		newclust[0][nnew]=qtot;
		newclust[1][nnew]=segmentation->Ix();
		newclust[2][nnew]=segmentation->Iy();
		newclust[3][nnew]=qp * qtot;
		newclust[4][nnew]=segmentation->ISector();
		newclust[5][nnew]=(Float_t) i;
		nnew++;
		
		
	    }
	} // Pad loop
    } // Cathode plane loop
}



















