///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Of Flight                                                           //
//  This class contains the basic functions for the Time Of Flight           //
//  detector. Functions specific to one particular geometry are              //
//  contained in the derived classes                                         //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliTOFClass.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTOF.h"
#include <TNode.h>
#include <TTUBE.h>
#include <TBRIK.h>
#include "AliRun.h"
#include "AliConst.h"
 
ClassImp(AliTOF)
 
//_____________________________________________________________________________
AliTOF::AliTOF()
{
  //
  // Default constructor
  //
  fIshunt   = 0;
}
 
//_____________________________________________________________________________
AliTOF::AliTOF(const char *name, const char *title)
       : AliDetector(name,title)
{
  //
  // AliTOF standard constructor
  // 
  fHits   = new TClonesArray("AliTOFhit",  405);
  //
  fIshunt     =  0;
  //
  SetMarkerColor(7);
  SetMarkerStyle(2);
  SetMarkerSize(0.4);
}

//_____________________________________________________________________________
void AliTOF::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a TOF hit
  //
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliTOFhit(fIshunt,track,vol,hits);
}
 
//_____________________________________________________________________________
void AliTOF::BuildGeometry()
{
  //
  // Build TOF ROOT geometry for the ALICE event viewver
  //
  TNode *Node, *Top;
  const int kColorTOF  = 27;
  //
  // Find top TNODE
  Top=gAlice->GetGeometry()->GetNode("alice");
  //
  // Define rotation matrixes
  new TRotMatrix("rot501","rot501",90,-18.94737,90,71.05263,0,0);
  new TRotMatrix("rot502","rot502",90,-37.89474,90,52.10526,0,0);
  new TRotMatrix("rot503","rot503",90,-56.84211,90,33.15789,0,0);
  new TRotMatrix("rot504","rot504",90,-75.78947,90,14.21053,0,0);
  new TRotMatrix("rot505","rot505",90,-94.73685,90,-4.736847,0,0);
  new TRotMatrix("rot506","rot506",90,-113.6842,90,-23.68421,0,0);
  new TRotMatrix("rot507","rot507",90,-132.6316,90,-42.63158,0,0);
  new TRotMatrix("rot508","rot508",90,-151.5789,90,-61.57895,0,0);
  new TRotMatrix("rot509","rot509",90,-170.5263,90,-80.52632,0,0);
  new TRotMatrix("rot510","rot510",90,-189.4737,90,-99.47369,0,0);
  new TRotMatrix("rot511","rot511",90,-208.4211,90,-118.4211,0,0);
  new TRotMatrix("rot512","rot512",90,-227.3684,90,-137.3684,0,0);
  new TRotMatrix("rot513","rot513",90,-246.3158,90,-156.3158,0,0);
  new TRotMatrix("rot514","rot514",90,-265.2632,90,-175.2632,0,0);
  new TRotMatrix("rot515","rot515",90,-284.2105,90,-194.2105,0,0);
  new TRotMatrix("rot516","rot516",90,-303.1579,90,-213.1579,0,0);
  new TRotMatrix("rot517","rot517",90,-322.1053,90,-232.1053,0,0);
  new TRotMatrix("rot518","rot518",90,-341.0526,90,-251.0526,0,0);
  new TRotMatrix("rot519","rot519",90,-360,90,-270,0,0);
  //
  // Position the different copies
  //  const Float_t rtof=366;
  // changed by Federico Carminati. TOF people should really look at this
  const Float_t rtof=381;
  const Int_t ntof=19;
  const Float_t angle=2*kPI/ntof;
  Float_t ang;
  const Float_t xtof = rtof*TMath::Sin(kPI/ntof);
  //
  // Define TOF basic volume
  new TBRIK("S_TOF1","TOF box","void",xtof,6.,175.);
  //
  Top->cd();
  ang=2.75*angle;
  Node = new TNode("FTO11","FTO11","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),175,"rot502");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO111","FTO111","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-175,"rot502");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=1.75*angle;
  Node = new TNode("FTO12","FTO12","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),175,"rot503");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO112","FTO112","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-175,"rot503");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();  
  ang=0.75*angle;
  Node = new TNode("FTO13","FTO13","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),175,"rot504");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node);
  //
  Top->cd();
  Node = new TNode("FTO113","FTO113","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-175,"rot504");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=-0.25*angle;
  Node = new TNode("FTO14","FTO14","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),175,"rot505");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO114","FTO114","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-175,"rot505");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=-1.25*angle;
  Node = new TNode("FTO15","FTO15","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),175,"rot506");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO115","FTO115","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-175,"rot506");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();  
  ang=kPI+1.25*angle;
  Node = new TNode("FTO16","FTO16","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),175,"rot513");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO116","FTO116","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-175,"rot513");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=kPI+0.25*angle;
  Node = new TNode("FTO17","FTO17","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),175,"rot514");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO117","FTO117","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-175,"rot514");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=kPI-0.75*angle;
  Node = new TNode("FTO18","FTO18","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),175,"rot515");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO118","FTO118","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-175,"rot515");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=kPI-1.75*angle;
  Node = new TNode("FTO19","FTO19","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),175,"rot516");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO119","FTO119","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-175,"rot516");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=kPI-2.75*angle;
  Node = new TNode("FTO110","FTO110","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),175,"rot517");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO120","FTO120","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-175,"rot517");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node);
  //
  // Define second TOF volume
  new TBRIK("S_TOF2","TOF box","void",xtof,6.,100.);
  //
  // Position the volume
  Top->cd();
  ang=-2.25*angle;
  Node = new TNode("FTO21","FTO21","S_TOF2",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),250,"rot507");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO27","FTO27","S_TOF2",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-250,"rot507");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=-3.25*angle;
  Node = new TNode("FTO22","FTO22","S_TOF2",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),250,"rot508");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO28","FTO28","S_TOF2",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-250,"rot508");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=-4.25*angle;
  Node = new TNode("FTO23","FTO23","S_TOF2",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),250,"rot509");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO29","FTO29","S_TOF2",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-250,"rot509");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=kPI+4.25*angle;
  Node = new TNode("FTO24","FTO24","S_TOF2",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),250,"rot510");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO210","FTO210","S_TOF2",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-250,"rot510");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=kPI+3.25*angle;
  Node = new TNode("FTO25","FTO25","S_TOF2",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),250,"rot511");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO211","FTO211","S_TOF2",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-250,"rot511");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=kPI+2.25*angle;
  Node = new TNode("FTO26","FTO26","S_TOF2",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),250,"rot512");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO212","FTO212","S_TOF2",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-250,"rot512");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  // Define third TOF volume
  new TBRIK("S_TOF3","TOF box","void",xtof,6.,75.);
  //
  // Position it
  Top->cd();
  ang=3.75*angle;
  Node = new TNode("FTO31","FTO31","S_TOF3",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),275,"rot501");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO34","FTO34","S_TOF3",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-275.,"rot501");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=kPI-3.75*angle;
  Node = new TNode("FTO32","FTO32","S_TOF3",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),275,"rot518");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO35","FTO35","S_TOF3",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-275,"rot518");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=kPI/2;
  Node = new TNode("FTO33","FTO33","S_TOF3",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),275.,"rot519");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO36","FTO36","S_TOF3",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-275.,"rot519");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
}

//_____________________________________________________________________________
void AliTOF::CreateGeometry()
{
  //
  // Common geometry code for version 2 and version 3 of the TOF
  //
  //Begin_Html
  /*
    <img src="gif/AliTOFv23.gif">
  */
  //End_Html
  //

  AliMC* pMC = AliMC::GetMC();
  
  const Double_t kPi=TMath::Pi();
  const Double_t kDegrad=kPi/180;
  const Double_t kRaddeg=180/kPi;
  //
  Float_t fil_rich;
  Int_t lmax;
  Float_t xtof, ytof, fil_step, phos_phi, phos_min, phos_max;
  Int_t lmax1;
  Float_t zcor1, zcor2, zcor3;
  Float_t ztof0, ztof1, ztof2, xcor2, ycor2;
  Int_t i;
  Float_t dx, dz, zl, xm, ym, phos_r;
  Int_t idrotm[101];
  Float_t phos_x;
  Float_t rp2, zm0, zm1, zm2;
  Float_t par[10], fil_min, fil_max;
  Float_t fil1;
  //
  Int_t *idtmed = gAlice->Idtmed();
  //
  phos_x = 214.6;
  phos_r = 467.;
  //phos_z = 260.;
  //rich_z = 472.5;
  xtof = 120.;
  ytof = 12.;
  ztof0 = 350.;
  ztof1 = 200.;
  ztof2 = 150.;
  //
  // frame thick along Z axis
  dz = 0.;
  //
  // frame thick along X axis
  dx = 0.;
  //
  // barrel length along Z axis
  zl = 720.;
  //
  // PHOS openings 
  fil_rich = 30.;
  phos_phi = TMath::ATan(phos_x / (phos_r * 2.));
  phos_min = (kPi - phos_phi * 4.) * kRaddeg;
  phos_max = (phos_phi * 4. + kPi) * kRaddeg;
  //
  // barrel radius in module contact point 
  par[0] = 370;
  par[1] = 390;
  par[2] = zl / 2.;
  pMC->Gsvolu("FBAR", "TUBE", idtmed[500], par, 3);
  //
  // --- Set module unseen ---
  pMC->Gspos("FBAR", 1, "ALIC", 0., 0., 0., 0, "ONLY");
  pMC->Gsatt("FBAR", "SEEN", 0);
  //
  // Number of TOF-block 
  lmax = 19;
  //
  // New size YTOF
  par[0] = xtof / 2.;
  par[1] = ytof / 2.;
  par[2] = ztof0 / 2.;
  pMC->Gsvolu("FTO1", "BOX ", idtmed[506], par, 3);
  pMC->Gsatt("FTO1", "SEEN", -2);
  par[2] = ztof1 / 2.;
  pMC->Gsvolu("FTO2", "BOX ", idtmed[506], par, 3);
  pMC->Gsatt("FTO2", "SEEN", -2);
  par[2] = ztof2 / 2.;
  pMC->Gsvolu("FTO3", "BOX ", idtmed[506], par, 3);
  pMC->Gsatt("FTO3", "SEEN", -2);
  //
  // Subtraction of TOF module boundaries 
  xm = xtof - dx * 2.;
  ym = ytof;
  zm0 = ztof0 - dz * 2.;
  zm1 = ztof1 - dz * 2.;
  zm2 = ztof2 - dz * 2.;
  //  
  // TOF module internal definitions
  //  
  TOFpc(xm, ym, zm0, zm1, zm2);
  //
  //rp1 = 382.;
  rp2 = 381;
  fil_step = 360. / lmax;
  fil_min = phos_min - fil_step * .5;
  fil_max = phos_max + fil_step * .5;
  zcor1 = 175.;
  zcor2 = ztof0 - ztof1 / 2.;
  zcor3 = ztof0 - ztof2 / 2.;
  for (i = 1; i <= lmax; ++i) {
    fil1 = fil_step * i;
    //xcor1 = rp1 * TMath::Sin(fil1 * kDegrad) + dx * TMath::Cos(fil1 * kDegrad);
    //ycor1 = rp1 * TMath::Cos(fil1 * kDegrad) - dx * TMath::Sin(fil1 * kDegrad);
    xcor2 = rp2 * TMath::Sin(fil1 * kDegrad);
    ycor2 = rp2 * TMath::Cos(fil1 * kDegrad);
    lmax1 = i + lmax;
    AliMatrix(idrotm[i], 90., -fil1, 90., 90. - fil1, 0., 0.);
    if (fil1 >= fil_min  && fil1 <= fil_max) {
      pMC->Gspos("FTO2", i, "FBAR", xcor2, ycor2, zcor2, idrotm[i], "ONLY");
      pMC->Gspos("FTO2", lmax1, "FBAR", xcor2, ycor2, -zcor2, idrotm[i], "ONLY");
    } else if (fil1 <= fil_rich || fil1 >= 360. - fil_rich) {
      par[2] = ztof2 / 2.;
      pMC->Gspos("FTO3", i, "FBAR", xcor2, ycor2, zcor3, idrotm[i], "ONLY");
      pMC->Gspos("FTO3", lmax1, "FBAR", xcor2, ycor2, -zcor3, idrotm[i], "ONLY");
    } else {
      par[2] = ztof0 / 2.;
      pMC->Gspos("FTO1", i, "FBAR", xcor2, ycor2, zcor1, idrotm[i], "ONLY");
      pMC->Gspos("FTO1", lmax1, "FBAR", xcor2, ycor2, -zcor1, idrotm[i], "ONLY");
    }
  }
}

//_____________________________________________________________________________
void AliTOF::DrawModule()
{
  //
  // Draw a shaded view of the common part of the TOF geometry
  // for versions 2 and 3
  //

  AliMC* pMC = AliMC::GetMC();
  
  // Set everything unseen
  pMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  pMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  pMC->Gsatt("FBAR","SEEN",0);
  pMC->Gsatt("FTO1","SEEN",1);
  pMC->Gsatt("FTO2","SEEN",1);
  pMC->Gsatt("FTO3","SEEN",1);
  //
  pMC->Gdopt("hide", "on");
  pMC->Gdopt("shad", "on");
  pMC->Gsatt("*", "fill", 7);
  pMC->SetClipBox(".");
  pMC->SetClipBox("*", 0, 1000, -1000, 1000, -1000, 1000);
  pMC->DefaultRange();
  pMC->Gdraw("alic", 40, 30, 0, 12, 9.5, .02, .02);
  pMC->Gdhead(1111, "Time Of Flight");
  pMC->Gdman(18, 4, "MAN");
  pMC->Gdopt("hide","off");
}

//_____________________________________________________________________________
void AliTOF::CreateMaterials()
{
  //
  // Defines TOF materials for all versions
  // Authors :   Maxim Martemianov, Boris Zagreev (ITEP)   18/09/98 
  //
  Int_t   ISXFLD = gAlice->Field()->Integ();
  Float_t SXMGMX = gAlice->Field()->Max();
  //
  Float_t ac[2]   = { 12.,16. };
  Float_t zc[2]   = { 6.,8. };
  Float_t wc[2]   = { 1.,2. };
  Float_t ag10[4] = { 12.,1.,16.,28. };
  Float_t zg10[4] = { 6.,1.,8.,14. };
  Float_t wmatg10[4] = { .259,.288,.248,.205 };
  Float_t adme[5] = { 12.,1.,16.,19.,79. };
  Float_t zdme[5] = { 6.,1.,8.,9.,35. };
  Float_t wmatdme[5] = { .4056,.0961,.2562,.1014,.1407 };
  Float_t aal[2] = { 27.,16. };
  Float_t zal[2] = { 13.,8. };
  Float_t wmatal[2] = { 2.,3. };
  //
  Int_t nlmatdme;
  Float_t epsil, stmin, dc, densg10, densal, deemax, stemax;
  Float_t densdme;
  Int_t nlmatg10, nlmatal;
  //
  // --- Vacuum
  AliMaterial(0, "Vacuum$", 1e-16, 1e-16, 1e-16, 1e16, 1e16);
  // --- Air
  AliMaterial(1, "Air$",14.61,7.3,0.001205,30423.24,67500.);
  // --- CO2 
  dc = .001977;
  AliMixture(7, "CO2$", ac, zc, dc, -2, wc);
  // --- G10 
  densg10  = 1.7;
  nlmatg10 = -4;
  AliMixture(5, "G10$", ag10, zg10, densg10, nlmatg10, wmatg10);
  // --- DME 
  densdme  = .00205;
  nlmatdme = 5;
  AliMixture(6, "DME ", adme, zdme, densdme, nlmatdme, wmatdme);
  // ---- ALUMINA (AL203) 
  densal  = 2.3;
  nlmatal = -2;
  AliMixture(8, "ALUMINA$", aal, zal, densal, nlmatal, wmatal);
  //       Previous data 
  //       EPSIL  =  0.1   ! Tracking precision, 
  //       STEMAX = 0.1      ! Maximum displacement for multiple scattering
  //       DEEMAX = 0.1    ! Maximum fractional energy loss, DLS 
  //       STMIN  = 0.1 
  //       New data from 
  epsil  = .001;
  stemax = -1.;
  deemax = -.3;
  stmin  = -.8;
  AliMedium(500, "Vacuum  $", 0, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(501, "Air $", 0, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(505, "G10$", 5, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(506, "DME$", 6, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(507, "CO2$", 7, 0, ISXFLD, SXMGMX, 10., -.01, -.1, .01, -.01);
  AliMedium(508, "ALUMINA$", 8, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(510, "DME$", 6, 1, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
}

//_____________________________________________________________________________
Int_t AliTOF::DistancetoPrimitive(Int_t , Int_t )
{
  //
  // Returns distance from mouse pointer to detector, default version
  //
  return 9999;
}
 
//_____________________________________________________________________________
void AliTOF::Init()
{
  //
  // Initialise TOF detector after it has been built
  //
  Int_t i;
  AliMC *pMC=AliMC::GetMC();
  //
  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" TOF_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  //
  // Set id of TOF sensitive volume
  fIdSens=pMC->VolId("FPG2");
  //
  for(i=0;i<80;i++) printf("*");
  printf("\n");
}

 
ClassImp(AliTOFhit)
 
//___________________________________________
AliTOFhit::AliTOFhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
  AliHit(shunt, track)
{
  //
  // Store a TOF hit
  //
  Int_t i;
  for (i=0;i<3;i++) fVolume[i] = vol[i];
  //
  // Position
  fX=hits[0];
  fY=hits[1];
  fZ=hits[2];
  //
  // Momentum
  fPx=hits[3];
  fPy=hits[4];
  fPz=hits[5];
  fPmom=hits[6];
  //
  // Time Of Flight
  fTof=hits[7];
}
 
 
