///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Photon Multiplicity Detector                                             //
//  This class contains the basic functions for the Photon Multiplicity      //
//  Detector. Functions specific to one particular geometry are              //
//  contained in the derived classes                                         //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliPMDClass.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>The responsible person for this module is
<a href="mailto:sub@vecdec.veccal.ernet.in">Subhasis Chattopadhyay</a>.
</font>
<pre>
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TBRIK.h>
#include <TNode.h>
#include "AliPMD.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliConst.h" 
 
static Float_t smod1[3], smod2[3], smod3[3], smod4[3];
static Int_t maxbox, kdet;
//static Float_t pmdin,pmdout,wafer;
static Float_t thgas,thcell,thmin,thmax,zdist,thlow,
  thhigh,edge;
static Int_t numqu;
static Float_t xbox[40][40], ybox[40][40];
static Int_t pindex[40][40];

ClassImp(AliPMD)
 
//_____________________________________________________________________________
AliPMD::AliPMD()
{
  //
  // Default constructor
  //
  fIshunt = 0;
}
 
//_____________________________________________________________________________
AliPMD::AliPMD(const char *name, const char *title)
  : AliDetector(name,title)
{
  //
  // Default constructor
  //

  // 
  // Allocate the array of hits
  fHits   = new TClonesArray("AliPMDhit",  405);
  
  fIshunt =  1;
  
  fPar[0] = 1;
  fPar[1] = 1;
  fPar[2] = 0.8;
  fPar[3] = 0.02;
  fIn[0]  = 6;
  fIn[1]  = 20;
  fIn[2]  = 600;
  fIn[3]  = 27;
  fIn[4]  = 27;
  fGeo[0] = 0;
  fGeo[1] = 0.2;
  fGeo[2] = 4;
  fPadSize[0] = 0.8;
  fPadSize[1] = 1.0;
  fPadSize[2] = 1.2;
  fPadSize[3] = 1.5;
}

//_____________________________________________________________________________
void AliPMD::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a PMD hit
  //
  TClonesArray &lhits = *fHits;
  AliPMDhit *newcell, *curcell;
  //  printf("PMD++ Adding energy %f, prim %d, vol %d %d %d %d\n",
  //	 hits[3],gAlice->GetPrimary(track-1),vol[0],vol[1],vol[2],vol[3]);
  newcell = new AliPMDhit(fIshunt, track, vol, hits);
  Int_t i;
  for (i=0; i<fNhits; i++) {
    //
    // See if this cell has already been hit
    curcell=(AliPMDhit*) lhits[i];
    if (*curcell==*newcell) {
      //      printf("Cell with same numbers found\n") ; curcell->Print();
      *curcell = *curcell+*newcell;
      //      printf("Cell after addition\n") ; curcell->Print();
      delete newcell;
      return;
    }
  }
  new(lhits[fNhits++]) AliPMDhit(newcell);
  delete newcell;
}
 
//_____________________________________________________________________________
void AliPMD::BuildGeometry()
{
  //
  // Build simple ROOT TNode geometry for event display
  //

  TNode *Node, *Top;
  const int kColorPMD  = kRed;

  //
  Top=gAlice->GetGeometry()->GetNode("alice");

  // PMD
  new TBRIK("S_PMD","PMD box","void",300,300,5);
  Top->cd();
  Node = new TNode("PMD","PMD","S_PMD",0,0,600,"");
  Node->SetLineColor(kColorPMD);
  fNodes->Add(Node);
}

//_____________________________________________________________________________
Int_t AliPMD::DistancetoPrimitive(Int_t , Int_t )
{
  //
  // Distance from mouse to detector on the screen
  // dummy routine
  //
   return 9999;
}
 
//_____________________________________________________________________________
void AliPMD::SetPAR(Float_t p1, Float_t p2, Float_t p3,Float_t p4)
{
  //
  // Set PMD parameters
  //
  fPar[0] = p1;
  fPar[1] = p2;
  fPar[2] = p3;
  fPar[3] = p4;
}
 
//_____________________________________________________________________________
void AliPMD::SetIN(Float_t p1, Float_t p2, Float_t p3,Float_t p4,Float_t p5)
{
  //
  // Set PMD parameters
  //
  fIn[0] = p1;
  fIn[1] = p2;
  fIn[2] = p3;
  fIn[3] = p4;
  fIn[4] = p5;
}
 
//_____________________________________________________________________________
void AliPMD::SetGEO(Float_t p1, Float_t p2, Float_t p3)
{
  //
  // Set geometry parameters
  //
  fGeo[0] = p1;
  fGeo[1] = p2;
  fGeo[2] = p3;
}
 
//_____________________________________________________________________________
void AliPMD::SetPadSize(Float_t p1, Float_t p2, Float_t p3,Float_t p4)
{
  //
  // Set pad size
  //
  fPadSize[0] = p1;
  fPadSize[1] = p2;
  fPadSize[2] = p3;
  fPadSize[3] = p4;
}
 
//_____________________________________________________________________________
void AliPMD::StepManager()
{
  //
  // Called at every step in PMD
  //
}

 
//_____________________________________________________________________________
void AliPMD::Undulation(char *undul, Float_t pitch, Float_t thick,
			Float_t zundul,	Float_t rundul, char (*cone)[5])
{
  //
  // RUNDUL   : Internal radius of the undulated chamber
  // THICK    : material thickness
  // PITCH    : one-QUARTER wave of undulation (cm)
  // ZUNDUL   : half length (cm)
  //
  // The undulated structure is desgned as a superposition of eight CONES
  // of suitable sizes, where the inner/outer radius of the cone increases,
  // then decreases, each half of the wave is assumed to be a semicircle,
  // which allows to calculate the thickness and the radii of the cone, by
  // dividing the semicircle into 4 parts of equal arc length.
  // Thus apear the constants 0.293 and 0.707.
  //

  const Float_t const1 = .293;
  const Float_t const2 = .707;
  
  AliMC* pMC = AliMC::GetMC();
  
  // Local variables 
  Int_t j, nwave;
  Float_t dcone1[5], dcone2[5], dcone3[5], dcone4[5], dcone5[5],
    dcone6[5], dcone7[5], dcone8[5];
  Float_t xc, yc, zc, dundul[3];
  Int_t *idtmed = gAlice->Idtmed();
  
  // Function Body 
  
  dcone1[0] = const1 * pitch / 2;
  dcone1[1] = rundul;
  dcone1[2] = dcone1[1] + thick;
  dcone1[3] = dcone1[1] + const2 * pitch;
  dcone1[4] = dcone1[3] + thick;
  
  dcone2[0] = const2 * pitch / 2;
  dcone2[1] = dcone1[3];
  dcone2[2] = dcone1[4];
  dcone2[3] = dcone2[1] + const1 * pitch;
  dcone2[4] = dcone2[3] + thick;

  dcone3[0] = dcone2[0];
  dcone3[1] = dcone2[3];
  dcone3[2] = dcone2[4];
  dcone3[3] = dcone2[1];
  dcone3[4] = dcone2[2];
  
  dcone4[0] = dcone1[0];
  dcone4[1] = dcone1[3];
  dcone4[2] = dcone1[4];
  dcone4[3] = dcone1[1];
  dcone4[4] = dcone1[2];
  
  dcone5[0] = dcone1[0];
  dcone5[1] = dcone1[1] - thick;
  dcone5[2] = dcone1[1];
  dcone5[3] = dcone5[1] - const2 * pitch;
  dcone5[4] = dcone5[3] + thick;
  
  dcone6[0] = dcone2[0];
  dcone6[1] = dcone5[3];
  dcone6[2] = dcone5[4];
  dcone6[3] = dcone6[1] - const1 * pitch;
  dcone6[4] = dcone6[3] + thick;
  
  dcone7[0] = dcone6[0];
  dcone7[1] = dcone6[3];
  dcone7[2] = dcone6[4];
  dcone7[3] = dcone5[3];
  dcone7[4] = dcone5[4];
  
  dcone8[0] = dcone5[0];
  dcone8[1] = dcone7[3];
  dcone8[2] = dcone7[4];
  dcone8[3] = dcone5[1];
  dcone8[4] = dcone5[2];
  
  pMC->Gsvolu(cone[0], "CONE", idtmed[606-1], dcone1, 5);
  pMC->Gsvolu(cone[1], "CONE", idtmed[606-1], dcone2, 5);
  pMC->Gsvolu(cone[2], "CONE", idtmed[606-1], dcone3, 5);
  pMC->Gsvolu(cone[3], "CONE", idtmed[606-1], dcone4, 5);
  pMC->Gsvolu(cone[4], "CONE", idtmed[606-1], dcone5, 5);
  pMC->Gsvolu(cone[5], "CONE", idtmed[606-1], dcone6, 5);
  pMC->Gsvolu(cone[6], "CONE", idtmed[606-1], dcone7, 5);
  pMC->Gsvolu(cone[7], "CONE", idtmed[606-1], dcone8, 5);
  pMC->Gsatt(cone[0], "SEEN", 0);
  pMC->Gsatt(cone[1], "SEEN", 0);
  pMC->Gsatt(cone[2], "SEEN", 0);
  pMC->Gsatt(cone[3], "SEEN", 0);
  pMC->Gsatt(cone[4], "SEEN", 0);
  pMC->Gsatt(cone[5], "SEEN", 0);
  pMC->Gsatt(cone[6], "SEEN", 0);
  pMC->Gsatt(cone[7], "SEEN", 0);
  
  // DEFINE AN IMAGINARY TUBE VOLUME FOR UNDULATED CHAMBER, FILL WITH VACUUM
    
  nwave = Int_t (zundul / (pitch * 2) + .1);
  dundul[2] = pitch * 2 * nwave;
  dundul[1] = rundul + pitch + thick * 2;
  //
  dundul[0] = 1e-4;
  pMC->Gsvolu(undul, "TUBE", idtmed[698-1], dundul, 3);
  
  xc = 0;
  yc = 0;
  zc = -dundul[2] + dcone1[0];
  for (j = 1; j <= nwave; ++j) {
    pMC->Gspos(cone[0], j, undul, xc, yc, zc, 0, "ONLY");
    zc = zc + dcone1[0] + dcone2[0];
    pMC->Gspos(cone[1], j, undul, xc, yc, zc, 0, "ONLY");
    zc = zc + dcone2[0] + dcone3[0];
    pMC->Gspos(cone[2], j, undul, xc, yc, zc, 0, "ONLY");
    zc = zc + dcone3[0] + dcone4[0];
    pMC->Gspos(cone[3], j, undul, xc, yc, zc, 0, "ONLY");
    zc = zc + dcone4[0] + dcone5[0];
    pMC->Gspos(cone[4], j, undul, xc, yc, zc, 0, "ONLY");
    zc = zc + dcone5[0] + dcone6[0];
    pMC->Gspos(cone[5], j, undul, xc, yc, zc, 0, "ONLY");
    zc = zc + dcone6[0] + dcone7[0];
    pMC->Gspos(cone[6], j, undul, xc, yc, zc, 0, "ONLY");
    zc = zc + dcone7[0] + dcone8[0];
    pMC->Gspos(cone[7], j, undul, xc, yc, zc, 0, "ONLY");
    zc = zc + dcone8[0] + dcone1[0];
  }
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Photon Multiplicity Detector Version 1                                   //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliPMDv1Class.gif">
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

ClassImp(AliPMDv1)
 
//_____________________________________________________________________________
AliPMDv1::AliPMDv1() : AliPMD()
{
  //
  // Default constructor 
  //
  fMedSens=0;
}
 
//_____________________________________________________________________________
AliPMDv1::AliPMDv1(const char *name, const char *title)
  : AliPMD(name,title)
{
  //
  // Standard constructor
  //
  fMedSens=0;
}

//_____________________________________________________________________________
void AliPMDv1::Coordnew()
{
  //
  // Find coordinates for pad geometry
  //
  // Author Y.P. Viyogi, VECC Calcutta
  //

  Float_t th1, th2, dbox, dist;
  //Float_t xoff[40][40], yoff[40][40];
  Int_t i, j, nbox;
  Int_t xoff1[3], yoff1[3], l;
  Float_t dmax, hole;
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
  
  //rlow = zdist * TMath::Tan(thlow * kPI / 180);
  //rhigh = zdist * TMath::Tan(thhigh * kPI / 180);
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
void AliPMDv1::Coordinates()
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
void AliPMDv1::CreateGeometry()
{
  //
  // Create geometry for Photon Multiplicity Detector Version 1
  //
  //Begin_Html
  /*
    <img src="gif/AliPMDv1.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="gif/AliPMDv1Tree.gif">
  */
  //End_Html
  CreatePads();
  CreateInside();
}
 
//_____________________________________________________________________________
void AliPMDv1::CreateInside()
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
    zp = zdist;
    pMC->Gspos("DPMD", 1, "ALIC", xp, yp, zp, 0, "ONLY");
    
}

//_____________________________________________________________________________
void AliPMDv1::CreatePads()
{
  //
  // Create the geometry of the pads
  // *** DEFINITION OF THE GEOMETRY OF THE PMD  *** 
  // *** DIFFERENT PADS WITH SIZES 8 MM, 10 MM, 12 MM AND 15 MM SQUARE 
  // -- Author :     Y.P. VIYOGI, 04/06/1997. 
  // -- Modified:    P.V.K.S.Baba(JU), 13-12-97. 
  
  AliMC* pMC = AliMC::GetMC();
  
  Int_t npad1, npad2, npad3, npad4;
  Float_t spad1[3], spad2[3], spad3[3], spad4[3];
  Float_t scpv1[3], scpv2[3], scpv3[3], scpv4[3];
  Int_t i, j;
  Float_t sstr1[3], spsw1[3], sstr2[3], spsw2[3], sstr3[3], spsw3[3],
    sstr4[3], spsw4[3];
  Float_t xa, ya, za, xb, yb, zb, xc, sw[3], yc, zc;
  Float_t sfe[3];
  Float_t spb[3], pad1, pad2, pad3, pad4;
  //  VOLUMES Names : begining with D for all PMD volumes, 
  //     DMO1 : MODULE TYPE 1 ( 8 MM PADS) 
  //     DM11 : MODULE TYPE 2 (10 MM PADS) 
  //     DM12 : MODULE TYPE 3 (12 MM PADS) 
  //     DM13 : MODULE TYPE 4 (15 MM PADS) 
  
  // The names of SIZE variables begin with S and have more meaningful
  // characters as shown below. 
  
  // 		VOLUME 	SIZE	MEDIUM	: 	REMARKS 
  // 		------	-----	------	: --------------------------- 
  
  // 		DPPB	SPB	PB	: PB Converter and its SIZE 
  // 		DPFE	SFE	FE	: FE Support Plate and its SIZE 
  
  //               DP11    SPAD2   GAS     : PAD TYPE 2 (10 MM) 
  //               DP12    SPAD2   GAS     : PAD TYPE 2 FOR CPV(10 MM) 
  //               DS11    SSTR2   FE      : STRIP OF IRON 
  //               DW11    SPSW2   G10     : PRESHOWER 
  //               DV11    SCPV2   G10     : CPV 
  
  //               DP13    SPAD3   GAS     : PAD TYPE 3 (12 MM) 
  //               DP14    SPAD3   GAS     : PAD TYPE 3 FOR CPV(12 MM) 
  //               DS12    SSTR3   FE      : STRIP OF IRON 
  //               DW12    SPSW3   G10     : PRESHOWER 
  //               DV12    SCPV3   G10     : CPV 
  
  //               DP15    SPAD4   GAS     : PAD TYPE 4 (15 MM) 
  //               DP16    SPAD4   GAS     : PAD TYPE 4 FOR CPV(15 MM) 
  //               DS13    SSTR4   FE      : STRIP OF IRON 
  //               DW13    SPSW4   G10     : PRESHOWER 
  //               DV13    SCPV4   G10     : CPV 
  
  //     ****************** VOLUME TREE ****************** 
  
  // 			DM11 (Module) 
  // 			       | 
  // 			       | 
  // 	------------------------------------------------- 
  //       |             |               |                 | 
  //       |             |               |                 | 
  //    DV11( CPV)      DPFE            DPPB              DW11(Preshower) 
  //       |					        | 
  // 	|                                               | 
  //    DS12(Strip)		                       DS11(Strip) 
  //       |                                               | 
  //       |					        | 
  //    DP12(Pads)                                        DP11(Pads) 
  
  //    ************************************************************ 
  
  // --- The above  gives the Volume Tree. PAD is a gas cell of size 
  // --- given by PADSIZE in the input cards.  STRIP is a collection of 
  //--- PADs in a row.   STRIPs are positioned in the PRESHOWER BOX. This is
  //--- then placed in the MODULE. The PSW and the MODULE have the same size
  // --- ; Lead converter, Iron support plate are also placed 
  // --- in the MODULE. 
  
  
  //        DATA PAD1,PAD2,PAD3,PAD4/4*0.8/ 
  //        DATA NPAD1,NPAD2,NPAD3,NPAD4/4*30/ 
  
  Int_t *idtmed = gAlice->Idtmed();
  
  // **** PAD SIZE 8 MM 
  //  pmdin  = fPar[0];
  //  pmdout = fPar[1];
  thgas  = fPar[2];
  thcell = fPar[3];
  thmin  = fIn[0];
  thmax  = fIn[1];
  zdist  = fIn[2];
  thlow  = fIn[3];
  thhigh = fIn[4];
  //  wafer  = fGeo[0];
  edge   = fGeo[1];
  numqu  = Int_t(fGeo[2]);
  
  
  // *BABA 
  pad1  = fPadSize[0];
  pad2  = fPadSize[1];
  pad3  = fPadSize[2];
  pad4  = fPadSize[3];
  npad1 = Int_t(24/fPadSize[0]);
  npad2 = Int_t(24/fPadSize[1]);
  npad3 = Int_t(24/fPadSize[2]);
  npad4 = Int_t(24/fPadSize[3]);
  // *BABA 
  spad1[0] = (pad1 - thcell) / 2.;
  spad1[1] = spad1[0];
  spad1[2] = thgas / 2;
  pMC->Gsvolu("DP21", "BOX ", idtmed[604], spad1, 3);
  pMC->Gsatt("DP21", "SEEN", 1);
  pMC->Gsvolu("DP22", "BOX ", idtmed[604], spad1, 3);
  pMC->Gsatt("DP22", "SEEN", 1);
  
  sstr1[0] = npad1*pad1/2;
  sstr1[1] = pad1/2;
  sstr1[2] = thgas/2;
  pMC->Gsvolu("DS21", "BOX ", idtmed[605], sstr1, 3);
  pMC->Gsatt("DS21", "SEEN", 1);
  pMC->Gsvolu("DS22", "BOX ", idtmed[605], sstr1, 3);
  pMC->Gsatt("DS22", "SEEN", 1);
  
  spsw1[0] = sstr1[0] + edge;
  spsw1[1] = spsw1[0];
  spsw1[2] = (thgas + .4) / 2;
  // 2 mm G10 Plate cover (NMATE = 808) 
  scpv1[0] = spsw1[0];
  scpv1[1] = spsw1[1];
  scpv1[2] = spsw1[2];
  pMC->Gsvolu("DW21", "BOX ", idtmed[607], spsw1, 3);
  pMC->Gsatt("DW21", "SEEN", 1);
  pMC->Gsvolu("DV21", "BOX ", idtmed[607], spsw1, 3);
  pMC->Gsatt("DV21", "SEEN", 1);
  
  // --- place  pads in a strip 
  xa = (-npad1 + 1.) * pad1 / 2.;
  ya = 0.;
  za = 0.;
  for (i = 1; i <= npad1; ++i) {
    pMC->Gsposp("DP21", i, "DS21", xa, ya, za, 0, "ONLY", spad1, 3);
    pMC->Gsposp("DP22", i, "DS22", xa, ya, za, 0, "ONLY", spad1, 3);
    xa += pad1;
  }
  // --- place  strips in the PRESHOWER AND CPV boxes 
  xb = 0.;
  yb = (-npad1 + 1.) * pad1 / 2.;
  zb = 0.;
  for (j = 1; j <= npad1; ++j) {
    pMC->Gsposp("DS21", j, "DW21", xb, yb, zb, 0, "ONLY", sstr1, 3);
    pMC->Gsposp("DS22", j, "DV21", xb, yb, zb, 0, "ONLY", sstr1, 3);
    yb += pad1;
  }
  
  // **** PAD SIZE 10 MM 
  
  spad2[0] = (pad2 - thcell) / 2.;
  spad2[1] = spad2[0];
  spad2[2] = thgas / 2;
  pMC->Gsvolu("DP11", "BOX ", idtmed[604], spad2, 3);
  pMC->Gsatt("DP11", "SEEN", 1);
  pMC->Gsvolu("DP12", "BOX ", idtmed[604], spad2, 3);
  pMC->Gsatt("DP12", "SEEN", 1);
  
  sstr2[0] = npad2 * pad2 / 2;
  sstr2[1] = pad2 / 2;
  sstr2[2] = thgas / 2;
  pMC->Gsvolu("DS11", "BOX ", idtmed[605], sstr2, 3);
  pMC->Gsatt("DS11", "SEEN", 1);
  pMC->Gsvolu("DS12", "BOX ", idtmed[605], sstr2, 3);
  pMC->Gsatt("DS12", "SEEN", 1);
  
  spsw2[0] = sstr2[0] + edge;
  spsw2[1] = spsw2[0];
  spsw2[2] = (thgas + .4) / 2;
  // 2 mm G10 Plate cover (NMATE = 808) 
  scpv2[0] = spsw2[0];
  scpv2[1] = spsw2[1];
  scpv2[2] = spsw2[2];
  pMC->Gsvolu("DW11","BOX ", idtmed[607], spsw2, 3);
  pMC->Gsatt("DW11", "SEEN", 1);
  pMC->Gsvolu("DV11","BOX ", idtmed[607], spsw2, 3);
  pMC->Gsatt("DV11", "SEEN", 1);
  
  // --- place  pads in a strip 
  xa = (-npad2 + 1.) * pad2 / 2.;
  ya = 0.;
  za = 0.;
  for (i = 1; i <= npad2; ++i) {
    pMC->Gsposp("DP11", i, "DS11", xa, ya, za, 0, "ONLY", spad2, 3);
    pMC->Gsposp("DP12", i, "DS12", xa, ya, za, 0, "ONLY", spad2, 3);
    xa += pad2;
  }
  // --- place  strips in the PRESHOWER AND CPV boxes 
  xb = 0.;
  yb = (-npad2 + 1.) * pad2 / 2.;
  zb = 0.;
  for (j = 1; j <= npad2; ++j) {
    pMC->Gsposp("DS11", j, "DW11", xb, yb, zb, 0, "ONLY", sstr2, 3);
    pMC->Gsposp("DS12", j, "DV11", xb, yb, zb, 0, "ONLY", sstr2, 3);
    yb += pad2;
  }
  
  // **** PAD SIZE 12 MM 
  
  spad3[0] = (pad3 - thcell) / 2.;
  spad3[1] = spad3[0];
  spad3[2] = thgas / 2;
  pMC->Gsvolu("DP13", "BOX ", idtmed[604], spad3, 3);
  pMC->Gsatt("DP13", "SEEN", 1);
  pMC->Gsvolu("DP14", "BOX ", idtmed[604], spad3, 3);
  pMC->Gsatt("DP14", "SEEN", 1);
  
  sstr3[0] = npad3 * pad3 / 2;
  sstr3[1] = pad3 / 2;
  sstr3[2] = thgas / 2;
  pMC->Gsvolu("DS13", "BOX ", idtmed[605], sstr3, 3);
  pMC->Gsatt("DS13", "SEEN", 1);
  pMC->Gsvolu("DS14", "BOX ", idtmed[605], sstr3, 3);
  pMC->Gsatt("DS14", "SEEN", 1);
  
  spsw3[0] = sstr3[0] + edge;
  spsw3[1] = spsw3[0];
  spsw3[2] = (thgas + .4) / 2;
  // 2 mm G10 Plate cover (NMATE = 808) 
  scpv3[0] = spsw3[0];
  scpv3[1] = spsw3[1];
  scpv3[2] = spsw3[2];
  pMC->Gsvolu("DW12","BOX ", idtmed[607], spsw3, 3);
  pMC->Gsatt("DW12", "SEEN", 1);
  pMC->Gsvolu("DV12","BOX ", idtmed[607], spsw3, 3);
  pMC->Gsatt("DV12", "SEEN", 1);
  
  // --- place  pads in a strip 
  xa = (-npad3 + 1.) * pad3 / 2.;
  ya = 0.;
  za = 0.;
  for (i = 1; i <= npad3; ++i) {
    pMC->Gsposp("DP13", i, "DS13", xa, ya, za, 0, "ONLY", spad3, 3);
    pMC->Gsposp("DP14", i, "DS14", xa, ya, za, 0, "ONLY", spad3, 3);
    xa += pad3;
  }
  // --- place  strips in the PRESHOWER AND CPV boxes 
  xb = 0.;
  yb = (-npad3 + 1.) * pad3 / 2.;
  zb = 0.;
  for (j = 1; j <= npad3; ++j) {
    pMC->Gsposp("DS13", j, "DW12", xb, yb, zb, 0, "ONLY", sstr3, 3);
    pMC->Gsposp("DS14", j, "DV12", xb, yb, zb, 0, "ONLY", sstr3, 3);
    yb += pad3;
  }
  
  // **** PAD SIZE 15 MM 
  
  spad4[0] = (pad4 - thcell) / 2.;
  spad4[1] = spad4[0];
  spad4[2] = thgas / 2;
  pMC->Gsvolu("DP15","BOX ", idtmed[604], spad4, 3);
  pMC->Gsatt("DP15", "SEEN", 1);
  pMC->Gsvolu("DP16","BOX ", idtmed[604], spad4, 3);
  pMC->Gsatt("DP16", "SEEN", 1);
  
  sstr4[0] = npad4 * pad4 / 2;
  sstr4[1] = pad4 / 2;
  sstr4[2] = thgas / 2;
  pMC->Gsvolu("DS15","BOX ", idtmed[605], sstr4, 3);
  pMC->Gsatt("DS15", "SEEN", 1);
  pMC->Gsvolu("DS16","BOX ", idtmed[605], sstr4, 3);
  pMC->Gsatt("DS16", "SEEN", 1);
  
  spsw4[0] = sstr4[0] + edge;
  spsw4[1] = spsw4[0];
  spsw4[2] = (thgas + .4) / 2;
  // 2 mm G10 Plate cover (NMATE = 808) 
  scpv4[0] = spsw4[0];
  scpv4[1] = spsw4[1];
  scpv4[2] = spsw4[2];
  pMC->Gsvolu("DW13","BOX ", idtmed[607], spsw4, 3);
  pMC->Gsatt("DW13", "SEEN", 1);
  pMC->Gsvolu("DV13","BOX ", idtmed[607], spsw4, 3);
  pMC->Gsatt("DV13", "SEEN", 1);
  
  // --- place  pads in a strip 
  xa = (-npad4 + 1.) * pad4 / 2.;
  ya = 0.;
  za = 0.;
  for (i = 1; i <= npad4; ++i) {
    pMC->Gsposp("DP15", i, "DS15", xa, ya, za, 0, "ONLY", spad4, 3);
    pMC->Gsposp("DP16", i, "DS16", xa, ya, za, 0, "ONLY", spad4, 3);
    xa += pad4;
  }
  // --- place  strips in the PRESHOWER AND CPV boxes 
  xb = 0.;
  yb = (-npad4 + 1.) * pad4 / 2.;
  zb = 0.;
  for (j = 1; j <= npad4; ++j) {
    pMC->Gsposp("DS15", j, "DW13", xb, yb, zb, 0, "ONLY", sstr4, 3);
    pMC->Gsposp("DS16", j, "DV13", xb, yb, zb, 0, "ONLY", sstr4, 3);
    yb += pad4;
  }
  
  
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
  
  smod1[0] = spsw1[0];
  smod1[1] = smod1[0];
  // 	SMOD1(3)=SPSW1(3)+SFE(3)+SW (3)+SCPV1(3) 
  smod1[2] = spsw1[2] + sfe[2] + spb[2] + scpv1[2];
  pMC->Gsvolu("DM21", "BOX ", idtmed[698], smod1, 3);
  
  smod2[0] = spsw2[0];
  smod2[1] = smod2[0];
  smod2[2] = spsw2[2] + sfe[2] + spb[2] + scpv2[2];
  pMC->Gsvolu("DM11", "BOX ", idtmed[698], smod2, 3);
  
  smod3[0] = spsw3[0];
  smod3[1] = smod3[0];
  smod3[2] = spsw3[2] + sfe[2] + spb[2] + scpv3[2];
  pMC->Gsvolu("DM12", "BOX ", idtmed[698], smod3, 3);
  
  smod4[0] = spsw4[0];
  smod4[1] = smod4[0];
  smod4[2] = spsw4[2] + sfe[2] + spb[2] + scpv4[2];
  pMC->Gsvolu("DM13", "BOX ", idtmed[698], smod4, 3);
  
  // **** MODULE TYPE 1 : ALWAYS WITH TUNSGTEN CONVERTER 
  
  // *** try with PB once 8.6.97 
  
  // ---  place gas box (as CPV), iron support, lead converter and gas box 
  // ---  (preshower) in the module 
  xc = 0.;
  yc = 0.;
  // --- First the CPV box 
  zc = -(spsw1[2] + sfe[2] + spb[2] + spsw1[2]) + spsw1[2];
  pMC->Gspos("DV21", 1, "DM21", xc, yc, zc, 0, "ONLY");
  // --- Then iron support plate 
  zc = zc + sfe[2] + spsw1[2];
  pMC->Gspos("DPFE", 1, "DM21", xc, yc, zc, 0, "ONLY");
  // --- Then  converter plate 
  zc = zc + sfe[2] + spb[2];
  pMC->Gspos("DPPB", 1, "DM21", xc, yc, zc, 0, "ONLY");
  // --- Lastly the preshower box 
  zc = zc + spb[2] + spsw1[2];
  pMC->Gspos("DW21", 1, "DM21", xc, yc, zc, 0, "ONLY");
  
  // **** MODULE TYPE 2 
  
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
  
  
  // **** MODULE TYPE 3 
  
  // ---  place gas box (as CPV), iron support, lead converter and gas box 
  // ---  (preshower) in the module 
  xc = 0.;
  yc = 0.;
  // --- First the CPV box 
  zc = -(spsw3[2] + sfe[2] + spb[2] + spsw3[2]) + spsw3[2];
  pMC->Gspos("DV12", 1, "DM12", xc, yc, zc, 0, "ONLY");
  // --- Then iron support plate 
  zc = zc + sfe[2] + spsw3[2];
  pMC->Gspos("DPFE", 1, "DM12", xc, yc, zc, 0, "ONLY");
  // --- Then lead converter plate 
  zc = zc + sfe[2] + spb[2];
  pMC->Gspos("DPPB", 1, "DM12", xc, yc, zc, 0, "ONLY");
  // --- Lastly the preshower box 
  zc = zc + spb[2] + spsw3[2];
  pMC->Gspos("DW12", 1, "DM12", xc, yc, zc, 0, "ONLY");
  

  // **** MODULE TYPE 4 
  
  // ---  place gas box (as CPV), iron support, lead converter and gas box 
  // ---  (preshower) in the module 
  xc = 0.;
  yc = 0.;
  // --- First the CPV box 
  zc = -(spsw4[2] + sfe[2] + spb[2] + spsw4[2]) + spsw4[2];
  pMC->Gspos("DV13", 1, "DM13", xc, yc, zc, 0, "ONLY");
  // --- Then iron support plate 
  zc = zc + sfe[2] + spsw4[2];
  pMC->Gspos("DPFE", 1, "DM13", xc, yc, zc, 0, "ONLY");
  // --- Then lead converter plate 
  zc = zc + sfe[2] + spb[2];
  pMC->Gspos("DPPB", 1, "DM13", xc, yc, zc, 0, "ONLY");
  // --- Lastly the preshower box 
  zc = zc + spb[2] + spsw4[2];
  pMC->Gspos("DW13", 1, "DM13", xc, yc, zc, 0, "ONLY");
  
}
 
//_____________________________________________________________________________
void AliPMDv1::DrawDetector()
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
  pMC->Gsatt("DP21","seen",0);
  pMC->Gsatt("DP22","seen",0);
  pMC->Gsatt("DS21","seen",1);
  pMC->Gsatt("DS22","seen",1);
  pMC->Gsatt("DW21","seen",0);
  pMC->Gsatt("DV21","seen",0);
  pMC->Gsatt("DP11","seen",0);
  pMC->Gsatt("DP12","seen",0);
  pMC->Gsatt("DS11","seen",1);
  pMC->Gsatt("DS12","seen",1);
  pMC->Gsatt("DW11","seen",0);
  pMC->Gsatt("DV11","seen",0);
  pMC->Gsatt("DP13","seen",0);
  pMC->Gsatt("DP14","seen",0);
  pMC->Gsatt("DS13","seen",1);
  pMC->Gsatt("DS14","seen",1);
  pMC->Gsatt("DW12","seen",0); 
  pMC->Gsatt("DV12","seen",0);
  pMC->Gsatt("DP15","seen",0);
  pMC->Gsatt("DP16","seen",0);
  pMC->Gsatt("DS15","seen",1);
  pMC->Gsatt("DS16","seen",1);
  pMC->Gsatt("DW13","seen",0);
  pMC->Gsatt("DV13","seen",0);
  pMC->Gsatt("DPPB","seen",1);
  pMC->Gsatt("DPW ","seen",1); 
  pMC->Gsatt("DPFE","seen",1);
  pMC->Gsatt("DM21","seen",1);
  pMC->Gsatt("DM11","seen",1);
  pMC->Gsatt("DM12","seen",1);
  pMC->Gsatt("DM13","seen",1);
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
void AliPMDv1::CreateMaterials()
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
void AliPMDv1::Init()
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
void AliPMDv1::StepManager()
{
  //
  // Called at each step in the PMD
  //
  Int_t   copy;
  Float_t hits[4], destep;
  Float_t center[3] = {0,0,0};
  Int_t   vol[4];
  
  AliMC* pMC=AliMC::GetMC();
  if(pMC->GetMedium() == fMedSens && (destep = pMC->Edep())) {
    
    // THIS PART MUST BE CHECKED
    pMC->CurrentVol(0, copy);
    vol[0]=copy;
    pMC->CurrentVolOff(1,0,copy);
    vol[1]=copy;
    pMC->CurrentVolOff(2,0,copy);
    vol[2]=copy;
    pMC->CurrentVolOff(1,0,copy);
    vol[3]=copy;
    pMC->Gdtom(center,hits,1);
    hits[3] = destep*1e9; //Number in eV
    AddHit(gAlice->CurrentTrack(), vol, hits);
  }
}

ClassImp(AliPMDhit)
  
//_____________________________________________________________________________
AliPMDhit::AliPMDhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
  AliHit(shunt, track)
{
  //
  // Add a PMD hit
  //
  Int_t i;
  for (i=0;i<4;i++) fVolume[i] = vol[i];
  fX=hits[0];
  fY=hits[1];
  fZ=hits[2];
  fEnergy=hits[3];
}
