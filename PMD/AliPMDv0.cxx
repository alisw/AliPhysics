///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Photon Multiplicity Detector Version 1                                   //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliPMDv0Class.gif">
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include "AliPMDv0.h"
#include "AliRun.h"
#include "AliMC.h" 
#include "AliConst.h" 
 
static Float_t smod2[3], smod3[3], smod4[3];
static Int_t maxbox, kdet;
static Float_t thgas,thmin,thmax,zdist,zdist1,thlow,
  thhigh,edge;
static Int_t numqu;
static Float_t xbox[40][40], ybox[40][40];
static Int_t pindex[40][40];
 
ClassImp(AliPMDv0)
//_____________________________________________________________________________
AliPMDv0::AliPMDv0() : AliPMD()
{
  //
  // Default constructor 
  //
  fMedSens=0;
}
 
//_____________________________________________________________________________
AliPMDv0::AliPMDv0(const char *name, const char *title)
  : AliPMD(name,title)
{
  //
  // Standard constructor
  //
  fMedSens=0;
}

//_____________________________________________________________________________
void AliPMDv0::Coordnew()
{
  //
  // Find coordinates for pad geometry
  //
  // Author Y.P. Viyogi, VECC Calcutta
  //

  Float_t th1, th2, dbox, dist;
  //Float_t xoff[40][40], yoff[40][40];
  Int_t i, j, nbox;
  Float_t rlow;
  Int_t xoff1[3], yoff1[3], l;
  Float_t rhigh, dmax, hole;
  Int_t kk, nhol;
  Float_t rr, xx, yy;
  
  th1 = thmin * kPI / 180;
  th2 = thmax * kPI / 180;
  /* ESTIMATES FOR OCTAGON */
  dist = zdist * TMath::Tan(th2);
  /* ***  04.06.97 Fixed Module size of 6 cm, 0 mm boundary. */
  /* ***  variable pad sizes of 0.3 mm, 0.5  mm, 1.0 mm and 1.2 mm */
  dbox = edge * 2 + 24;
  maxbox = Int_t(dist / dbox + .5);
  dmax= maxbox * dbox;
  /* NOW GET THE HOLE SIZE ETC. */
  hole = zdist * TMath::Tan(th1);
  nhol = Int_t(hole / dbox + .5);
  hole = nhol * dbox;
  
  rlow = zdist * TMath::Tan(thlow * kPI / 180);
  rhigh = zdist * TMath::Tan(thhigh * kPI / 180);
  for (i = 1; i <= 40; ++i) {
    for (j = 1; j <= 40; ++j) {
      //index[j][i] = 0;
      //xoff[j][i] = 0;
      //yoff[j][i] = 0;
      xbox[j][i] = 0;
      /* L5: */
      ybox[j][i] = 0;
    }
  }
  
  // NOW START PLACING THE BOXES IN VARIOUS LAYERS, START FROM THE CENTRE 
  
  yy = dbox / 2;
  for(i=0;i<3;i++) yoff1[i]=0;
  nbox = 0;
  //        PRINT*,'MAXBOX=',MAXBOX 
  for (i = 1; i <= maxbox; ++i) {
    xx = dbox / 2;
    for(j=0;j<3;j++) xoff1[j]=0;
	for (j = 1; j <= maxbox; ++j) {
	  rr = sqrt(xx*xx+yy*yy);
	  if (rr >= hole && rr <= dmax) {
	    //  BOX CAN BE FITTED 
	    //index[j][i] = 2;
	    //if (rr < rlow) index[j][i] = 1;
	    //else if (rr > rhigh) index[j][i] = 3;
	    xbox[j][i] = xx;
	    ybox[j][i] = yy;
	      ++nbox;
	      //xoff[j][i] = xoff1[index[j][i] - 1];
	      //yoff[j][i] = yoff1[index[j][i] - 1];
	  }
	  if (kdet == 1) kk = 1; else kk = 0;
	  for (l = 1; l <= 3; ++l)
	    xoff1[l - 1] += fNumPads[l + kk - 1];
	  xx += dbox;
	}
	
	if (kdet == 1) kk = 1; else kk=0;
	
	for (l = 1; l <= 3; ++l)
	  yoff1[l - 1] += fNumPads[l + kk - 1];
	yy += dbox;
  }
}

//_____________________________________________________________________________
void AliPMDv0::Coordinates()
{
  //
  //  SUBROUTINE TO COMPUTE THE X- AND Y- COORDINATES OF THE BOXES 
  //  WHICH CAN FIT INTO THE CIRCULAR REGION BETWEEN THE GIVEN ANGLES. 
  //  INPUT : ZDIST, THMIN, THMAX, PADSIZE (FOR INSIDE and OUTSIDE PMD). 
  //  ALL DIMENSIONS IN CM. 
  // -- Author :	Y.P. VIYOGI, 10/05/1996. 

  Float_t hole, dmax, dbox;
  Int_t nhol;
  Float_t dist;
  Int_t nbox;
  Float_t rlow;
  Int_t i, j;
  Float_t rhigh, rr, xx, yy, th1, th2;
  
  th1 = thmin*kPI/180;
  th2 = thmax*kPI/180;
  // ESTIMATES FOR OCTAGON 
  dist = zdist * TMath::Tan(th2);
  // ***  04.06.97 Fixed Module size of 24 cm, 3 mm boundary. 
  // ***  variable pad sizes of 8 mm, 10 mm, 12mm and 15 mm 
  dbox   = edge*2 + 24.;
  maxbox = Int_t(dist / dbox + .5);
  dmax   = maxbox*dbox;
  // NOW GET THE HOLE SIZE ETC. 
  hole = zdist * TMath::Tan(th1);
  nhol = Int_t(hole / dbox + .5);
  hole = nhol * dbox;
  
  rlow  = zdist * TMath::Tan(thlow*kPI/180);
  rhigh = zdist * TMath::Tan(thhigh*kPI/180);
  for (i = 0; i < 40; ++i) {
    for (j = 0; j < 40; ++j) {
      pindex[j][i] = 0;
      xbox[j][i]   = 0;
      ybox[j][i]   = 0;
    }
  }
  
  //  NOW START PLACING THE BOXES IN VARIOUS LAYERS, START FROM THE CENTRE 
  yy   = dbox / 2;
    nbox = 0;
    for (i = 0; i < maxbox; ++i) {
      xx = dbox / 2;
      for (j = 0; j < maxbox; ++j) {
	rr = TMath::Sqrt(xx*xx + yy*yy);
	if (rr >= hole && rr <= dmax) {  //  BOX CAN BE FITTED 
	  pindex[j][i] = 2;
	  if (rr < rlow)  pindex[j][i] = 1;
	  if (rr > rhigh) pindex[j][i] = 3;
	  xbox[j][i] = xx;
	  ybox[j][i] = yy;
	  ++nbox;
	}
	xx += dbox;
      }
      yy += dbox;
    }
}
 
//_____________________________________________________________________________
void AliPMDv0::CreateGeometry()
{
  //
  // Create geometry for Photon Multiplicity Detector Version 1
  //
  //Begin_Html
  /*
    <img src="gif/AliPMDv0.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="gif/AliPMDv0Tree.gif">
  */
  //End_Html
  CreatePads();
  CreateInside();
}
 
//_____________________________________________________________________________
void AliPMDv0::CreateInside()
{
  //
  // Create inside of Pads
  //
  // -- Author :     Y.P. VIYOGI, 07/05/1996. 
  // -- Modified:    P.V.K.S.Baba(JU), 15-12-97. 
  
  Float_t sipmd[3] = { 300.,300.,5. };
  
  Int_t i2;
  
  Float_t xiqa[4], yiqa[4];
  Int_t inum2, inum3, inum4, i, j, k;
  Float_t siqad[4];
  Float_t zd, xd, yd, xp, yp, zp;
  Int_t idrotm[100];
  
  Int_t *idtmed = gAlice->Idtmed();    
  
  //  VOLUMES Names : begining with D for all PMD volumes, 
  // The names of SIZE variables begin with S and have more meaningful
  // characters as shown below. 
  
  // 		VOLUME 	SIZE	MEDIUM	: 	REMARKS 
  // 		------	-----	------	: --------------------------- 
  
  // 		DPMD	SIPMD	AIR	: INSIDE PMD  and its SIZE 
  
  
  
  // *** Define the  DPMD   Volume and fill with air *** 

  AliMC* pMC = AliMC::GetMC();
  
  pMC->Gsvolu("DPMD", "BOX ", idtmed[698], sipmd, 3);
  
  // *** Define DIQU Volume and fill with air 
  siqad[0] = sipmd[0] / 2. - 1.;
  siqad[1] = sipmd[1] / 2. - 1.;
  siqad[2] = sipmd[2];
  pMC->Gsvolu("DIQU","BOX ", idtmed[698], siqad, 3);
  pMC->Gsatt("DIQU", "SEEN", 1);
  
  
  // --- Place the modules in INSIDE PMD (DPMD) 
  // --- FIRST CALCULATE THE COORDINATES OF THE MODULES WHICH CAN BE 
  // --- ACCOMODATED. 
  
  kdet = 1;
  Coordinates();
  
  //inum = 0;
  zd   = 0.;
  AliMatrix(idrotm[1], 90., 0.,   90.,  90., 0., 0.);
  AliMatrix(idrotm[2], 90., 180., 90.,  90., 0., 0.);
  AliMatrix(idrotm[3], 90., 180., 90., 270., 0., 0.);
  AliMatrix(idrotm[4], 90., 0.,   90., 270., 0., 0.);
  // ****  Filling the DIQU Vol. (One Quadrant) 
  inum2 = 0;
  inum3 = 0;
    inum4 = 0;
    for (i = 0; i < maxbox; ++i) {
      i2 = maxbox;
      for (j = 0; j < i2; ++j) {
	if (xbox[j][i] <= 0 && ybox[j][i] <= 0) continue;
	xd = xbox[j][i] - siqad[0];
	yd = ybox[j][i] - siqad[1];
	if (pindex[j][i] == 1) {
	  ++inum2;
	  pMC->Gsposp("DM11", inum2, "DIQU", xd, yd, zd, 0, "ONLY", smod2, 3);
	}
	if (pindex[j][i] == 2) {
	  ++inum3;
	  pMC->Gsposp("DM12", inum3, "DIQU", xd, yd, zd, 0, "ONLY", smod3, 3);
	}
	if (pindex[j][i] == 3) {
	  ++inum4;
	  pMC->Gsposp("DM13", inum4, "DIQU", xd, yd, zd, 0, "ONLY", smod4, 3);
	}
      }
    }
    xiqa[0] = siqad[0];
    xiqa[1] = -siqad[0];
    xiqa[2] = xiqa[1];
    xiqa[3] = xiqa[0];
    yiqa[0] = siqad[0];
    yiqa[1] = yiqa[0];
    yiqa[2] = -siqad[0];
    yiqa[3] = yiqa[2];
    i2      = numqu;
    for (k = 1; k <= i2; ++k) {
      pMC->Gsposp("DIQU", k, "DPMD", xiqa[k-1], yiqa[k-1], zd, idrotm[k], "ONLY", siqad, 3);
    }
    
    // --- Place the DPMD in ALICE with front edge 6.0m from vertex  --- 
    xp = 0.;
    yp = 0.;
    zp = zdist1;
    pMC->Gspos("DPMD", 1, "ALIC", xp, yp, zp, 0, "ONLY");
    
}

//_____________________________________________________________________________
void AliPMDv0::CreatePads()
{
  //
  // Create the geometry of the pads
  // *** DEFINITION OF THE GEOMETRY OF THE PMD  *** 
  // *** DIFFERENT PADS WITH SIZES 8 MM, 10 MM, 12 MM AND 15 MM SQUARE 
  // -- Author :     Y.P. VIYOGI, 04/06/1997. 
  // -- Modified:    P.V.K.S.Baba(JU), 13-12-97. 
  
  AliMC* pMC = AliMC::GetMC();
  
  Int_t npad2;
  Float_t /* scpv1[3], */ scpv2[3] /*, scpv3[3], scpv4[3] */;
  Float_t  spsw1[3], spsw2[3];//, spsw3[3], spsw4[3];
  Float_t  sw[3], xc, yc, zc;
  Float_t sfe[3];
  Float_t spb[3], pad1, pad2, pad3, pad4;
  //  VOLUMES Names : begining with D for all PMD volumes, 
  
  //     DM11 : MODULE TYPE 
  
  // The names of SIZE variables begin with S and have more meaningful
  // characters as shown below. 
  
  // 		VOLUME 	SIZE	MEDIUM	: 	REMARKS 
  // 		------	-----	------	: --------------------------- 
  
  // 		DPPB	SPB	PB	: PB Converter and its SIZE 
  // 		DPFE	SFE	FE	: FE Support Plate and its SIZE 
  
  //               DW11    SPSW3   G10     : PRESHOWER 
  //               DV11    SCPV3   G10     : CPV 
  //     ****************** VOLUME TREE ****************** 
  
  // 			DM11 (Module) 
  // 			       | 
  // 			       | 
  // 	------------------------------------------------- 
  //       |             |               |                 | 
  //       |             |               |                 | 
  //    DV11( CPV)      DPFE            DPPB              DW11(Preshower) 
  //    ************************************************************ 
  
  
  
  Int_t *idtmed = gAlice->Idtmed();
  
  thgas  = fPar[2];
  thmin  = fIn[0];
  thmax  = fIn[1];
  zdist1  = fIn[2];
  zdist  = TMath::Abs(zdist1);
  thlow  = fIn[3];
  thhigh = fIn[4];
  edge   = fGeo[1];
  numqu  = Int_t(fGeo[2]);
  
  pad1  = fPadSize[0];
  pad2  = fPadSize[1];
  pad3  = fPadSize[2];
  pad4  = fPadSize[3];
  npad2 = Int_t(24/fPadSize[1]);
  
  spsw2[0] = (npad2 * pad2)/2 + edge;
  spsw2[1] = spsw2[0];
  spsw2[2] = (thgas + .4) / 2;
  scpv2[0] = spsw2[0];
  scpv2[1] = spsw2[1];
  scpv2[2] = spsw2[2];
// The modules (DW11 and DV11 are filed with gas, G10 plate is ignored)
  pMC->Gsvolu("DW11","BOX ", idtmed[604], spsw2, 3);
  pMC->Gsatt("DW11", "SEEN", 1);
  pMC->Gsvolu("DV11","BOX ", idtmed[604], spsw2, 3);
  pMC->Gsatt("DV11", "SEEN", 1);
  
  // --- DEFINE MODULES, IRON, TUNGSTEN AND LEAD VOLUMES 
  
  
  spb[0] = spsw1[0];
  spb[1] = spsw1[1];
  spb[2] = .75;
  pMC->Gsvolu("DPPB","BOX ", idtmed[600], spb, 3);
  pMC->Gsatt("DPPB", "SEEN", 1);
  
  sw[0] = spsw1[0];
  sw[1] = spsw1[1];
  sw[2] = 0.9/2.;
  pMC->Gsvolu("DPW ","BOX ", idtmed[600], sw, 3);
  pMC->Gsatt("DPW ", "SEEN", 1);
  
  sfe[0] = spsw1[0];
  sfe[1] = spsw1[1];
  sfe[2] = 0.6/2.;
  pMC->Gsvolu("DPFE","BOX ", idtmed[605], sfe, 3);
  pMC->Gsatt("DPFE", "SEEN", 1);
  
  smod2[0] = spsw2[0];
  smod2[1] = smod2[0];
  smod2[2] = spsw2[2] + sfe[2] + spb[2] + scpv2[2];
  pMC->Gsvolu("DM11", "BOX ", idtmed[698], smod2, 3);
  
  // ---  place gas box (as CPV), iron support, lead converter and gas box 
  // ---  (preshower) in the module 
  xc = 0.;
  yc = 0.;
  // --- First the CPV box 
  zc = -(spsw2[2] + sfe[2] + spb[2] + spsw2[2]) + spsw2[2];
  pMC->Gspos("DV11", 1, "DM11", xc, yc, zc, 0, "ONLY");
  // --- Then iron support plate 
  zc = zc + sfe[2] + spsw2[2];
  pMC->Gspos("DPFE", 1, "DM11", xc, yc, zc, 0, "ONLY");
  // --- Then lead converter plate 
  zc = zc + sfe[2] + spb[2];
  pMC->Gspos("DPPB", 1, "DM11", xc, yc, zc, 0, "ONLY");
  // --- Lastly the preshower box 
  zc = zc + spb[2] + spsw2[2];
  pMC->Gspos("DW11", 1, "DM11", xc, yc, zc, 0, "ONLY");
  
}
 
//_____________________________________________________________________________
void AliPMDv0::DrawModule()
{
  //
  // Draw a shaded view of the Photon Multiplicity Detector
  //

  AliMC* pMC = AliMC::GetMC();
  
  pMC->Gsatt("*", "seen", -1);
  pMC->Gsatt("alic", "seen", 0);
  //
  // Set the visibility of the components
  // 
  pMC->Gsatt("DW11","seen",0);
  pMC->Gsatt("DV11","seen",0);
  pMC->Gsatt("DPPB","seen",1);
  pMC->Gsatt("DPW ","seen",1); 
  pMC->Gsatt("DPFE","seen",1);
  pMC->Gsatt("DM11","seen",1);
  pMC->Gsatt("DPMD","seen",0);
  pMC->Gsatt("DIQU","seen",0);
  //
  pMC->Gdopt("hide", "on");
  pMC->Gdopt("shad", "on");
  pMC->Gsatt("*", "fill", 7);
  pMC->SetClipBox(".");
  pMC->SetClipBox("*", 0, 3000, -3000, 3000, -6000, 6000);
  pMC->DefaultRange();
  pMC->Gdraw("alic", 40, 30, 0, 22, 15.5, .04, .04);
  pMC->Gdhead(1111, "Photon Multiplicity Detector Version 1");
  pMC->Gdman(17, 5, "MAN");
  pMC->Gdopt("hide", "off");
}

//_____________________________________________________________________________
void AliPMDv0::CreateMaterials()
{
  //
  // Create materials for the PMD version 1
  //
  // ORIGIN    : Y. P. VIYOGI 
  //
  
  AliMC* pMC = AliMC::GetMC();
  
  // --- The Argon- CO2 mixture --- 
  Float_t ag[2] = { 39.95 };
  Float_t zg[2] = { 18. };
  Float_t wg[2] = { .8,.2 };
  Float_t dar   = .001782;   // --- Ar density in g/cm3 --- 
  // --- CO2 --- 
  Float_t ac[2] = { 12.,16. };
  Float_t zc[2] = { 6.,8. };
  Float_t wc[2] = { 1.,2. };
  Float_t dc    = .001977;
  Float_t dco   = .002;  // --- CO2 density in g/cm3 ---
  
  Float_t absl, radl, a, d, z;
  Float_t dg;
  Float_t x0ar;
  Float_t buf[1];
  Int_t nbuf;
  
  Int_t *idtmed = gAlice->Idtmed();
  Int_t isxfld = gAlice->Field()->Integ();
  Float_t sxmgmx = gAlice->Field()->Max();
  
  // --- Define the various materials for GEANT --- 
  AliMaterial(1, "Pb    $", 207.19, 82., 11.35, .56, 18.5);
  x0ar = 19.55 / dar;
  AliMaterial(2, "Argon$", 39.95, 18., dar, x0ar, 6.5e4);
  AliMixture(3, "CO2  $", ac, zc, dc, -2, wc);
  AliMaterial(4, "Al   $", 26.98, 13., 2.7, 8.9, 18.5);
  AliMaterial(6, "Fe   $", 55.85, 26., 7.87, 1.76, 18.5);
  AliMaterial(7, "W    $", 183.85, 74., 19.3, .35, 10.3);
  AliMaterial(8, "G10  $", 20., 10., 1.7, 19.4, 999);
  AliMaterial(9, "SILIC$", 28.09, 14., 2.33, 9.36, 45.);
  AliMaterial(10, "Be   $", 9.01, 4., 1.848, 35.3, 36.7);
  AliMaterial(15, "Cu   $", 63.54, 29., 8.96, 1.43, 15.);
  AliMaterial(16, "C    $", 12.01, 6., 2.265, 18.8, 49.9);
  
  AliMaterial(96, "MYLAR$", 8.73, 4.55, 1.39, 28.7, 62.);
  AliMaterial(97, "CONCR$", 20., 10., 2.5, 10.7, 40.);
  AliMaterial(98, "Vacum$", 1e-9, 1e-9, 1e-9, 1e16, 1e16);
  AliMaterial(99, "Air  $", 14.61, 7.3, .0012, 30420., 67500.);
  
  // 	define gas-mixtures 
  
  char namate[21];
  pMC->Gfmate((*fIdmate)[3], namate, a, z, d, radl, absl, buf, nbuf);
  ag[1] = a;
  zg[1] = z;
  dg = (dar * 4 + dco) / 5;
  AliMixture(5, "ArCO2$", ag, zg, dg, 2, wg);
  
  // Define tracking media 
  AliMedium(601, "Pb conv.$", 1,  0, 0, isxfld, sxmgmx, 1., .1, .01, .1);
  AliMedium(607, "W  conv.$", 7,  0, 0, isxfld, sxmgmx, 1., .1, .01, .1);
  AliMedium(608, "G10plate$", 8,  0, 0, isxfld, sxmgmx, 1., .1, .01, .1);
  AliMedium(604, "Al      $", 4,  0, 0, isxfld, sxmgmx, .1,  .1, .01, .1);
  AliMedium(606, "Fe      $", 6,  0, 0, isxfld, sxmgmx, .1,  .1, .01, .1);
  AliMedium(605, "ArCO2   $", 5,  1, 0, isxfld, sxmgmx, .1,  .1, .1,  .1);
  AliMedium(609, "SILICON $", 9,  1, 0, isxfld, sxmgmx, .1,  .1, .1,  .1);
  AliMedium(610, "Be      $", 10, 0, 0, isxfld, sxmgmx, .1,  .1, .01, .1);
  AliMedium(698, "Vacuum  $", 98, 0, 0, isxfld, sxmgmx, 1., .1, .1,  10);
  AliMedium(699, "Air gaps$", 99, 0, 0, isxfld, sxmgmx, 1., .1, .1,  .1);
  AliMedium(615, "Cu      $", 15, 0, 0, isxfld, sxmgmx, .1,  .1, .01, .1);
  AliMedium(616, "C       $", 16, 0, 0, isxfld, sxmgmx, .1,  .1, .01, .1);
  
  // --- Generate explicitly delta rays in the iron, aluminium and lead --- 
  pMC->Gstpar(idtmed[600], "LOSS", 3.);
  pMC->Gstpar(idtmed[600], "DRAY", 1.);
  
  pMC->Gstpar(idtmed[603], "LOSS", 3.);
  pMC->Gstpar(idtmed[603], "DRAY", 1.);
  
  pMC->Gstpar(idtmed[604], "LOSS", 3.);
  pMC->Gstpar(idtmed[604], "DRAY", 1.);
  
  pMC->Gstpar(idtmed[605], "LOSS", 3.);
  pMC->Gstpar(idtmed[605], "DRAY", 1.);
  
  pMC->Gstpar(idtmed[606], "LOSS", 3.);
  pMC->Gstpar(idtmed[606], "DRAY", 1.);
  
  pMC->Gstpar(idtmed[607], "LOSS", 3.);
  pMC->Gstpar(idtmed[607], "DRAY", 1.);
  
  // --- Energy cut-offs in the Pb and Al to gain time in tracking --- 
  // --- without affecting the hit patterns --- 
  pMC->Gstpar(idtmed[600], "CUTGAM", 1e-4);
  pMC->Gstpar(idtmed[600], "CUTELE", 1e-4);
  pMC->Gstpar(idtmed[600], "CUTNEU", 1e-4);
  pMC->Gstpar(idtmed[600], "CUTHAD", 1e-4);
  pMC->Gstpar(idtmed[605], "CUTGAM", 1e-4);
  pMC->Gstpar(idtmed[605], "CUTELE", 1e-4);
  pMC->Gstpar(idtmed[605], "CUTNEU", 1e-4);
  pMC->Gstpar(idtmed[605], "CUTHAD", 1e-4);
  pMC->Gstpar(idtmed[606], "CUTGAM", 1e-4);
  pMC->Gstpar(idtmed[606], "CUTELE", 1e-4);
  pMC->Gstpar(idtmed[606], "CUTNEU", 1e-4);
  pMC->Gstpar(idtmed[606], "CUTHAD", 1e-4);
  pMC->Gstpar(idtmed[603], "CUTGAM", 1e-4);
  pMC->Gstpar(idtmed[603], "CUTELE", 1e-4);
  pMC->Gstpar(idtmed[603], "CUTNEU", 1e-4);
  pMC->Gstpar(idtmed[603], "CUTHAD", 1e-4);
  pMC->Gstpar(idtmed[609], "CUTGAM", 1e-4);
  pMC->Gstpar(idtmed[609], "CUTELE", 1e-4);
  pMC->Gstpar(idtmed[609], "CUTNEU", 1e-4);
  pMC->Gstpar(idtmed[609], "CUTHAD", 1e-4);
  
  // --- Prevent particles stopping in the gas due to energy cut-off --- 
  pMC->Gstpar(idtmed[604], "CUTGAM", 1e-5);
  pMC->Gstpar(idtmed[604], "CUTELE", 1e-5);
  pMC->Gstpar(idtmed[604], "CUTNEU", 1e-5);
  pMC->Gstpar(idtmed[604], "CUTHAD", 1e-5);
  pMC->Gstpar(idtmed[604], "CUTMUO", 1e-5);
}

//_____________________________________________________________________________
void AliPMDv0::Init()
{
  //
  // Initialises PMD detector after it has been built
  //
  Int_t i;
  kdet=1;
  //
  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" PMD_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  printf("                 PMD simulation package initialised\n");
  printf(" parameters of pmd\n");
  printf("%6d %10.2f %10.2f %10.2f %10.2f %10.2f\n",kdet,thmin,thmax,zdist,thlow,thhigh);
  //
  for(i=0;i<80;i++) printf("*");
  printf("\n");
  //
  Int_t *idtmed = gAlice->Idtmed();
  fMedSens=idtmed[605-1];
}

//_____________________________________________________________________________
void AliPMDv0::StepManager()
{
  //
  // Called at each step in the PMD
  //
  Int_t   copy;
  Float_t hits[4], destep;
  Float_t center[3] = {0,0,0};
  Int_t   vol[5];
  Text_t namep[5];
  
  AliMC* pMC=AliMC::GetMC();
  if(pMC->GetMedium() == fMedSens && (destep = pMC->Edep())) {
    
    pMC->CurrentVol(namep, copy);
    vol[0]=copy;
    pMC->CurrentVolOff(1,namep,copy);
    vol[1]=copy;
    pMC->CurrentVolOff(2,namep,copy);
    vol[2]=copy;
    if(strncmp(namep,"DW11",4))vol[2]=1;
    if(strncmp(namep,"DV11",4))vol[2]=2;
    pMC->CurrentVolOff(3,namep,copy);
    vol[3]=copy;
    pMC->CurrentVolOff(4,namep,copy);
    vol[4]=copy;
    pMC->Gdtom(center,hits,1);
    hits[3] = destep*1e9; //Number in eV
    AddHit(gAlice->CurrentTrack(), vol, hits);
  }
}
