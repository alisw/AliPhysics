///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Of Flight                                                           //
//  This class contains the basic functions for the Time Of Flight           //
//  detector. Functions specific to one particular geometry are              //
//  contained in the derived classes                                         //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliTOFClass.gif">
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
#include <stdlib.h>
 
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
  //
  // Check that FRAME is there otherwise we have no place where to
  // put TOF
  AliModule* FRAME=gAlice->GetModule("FRAME");
  if(!FRAME) {
    Error("Ctor","TOF needs FRAME to be present\n");
    exit(1);
  }
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
  new TRotMatrix("rot501","rot501",90,-20,90,90-20,0,0);
  new TRotMatrix("rot502","rot502",90,-40,90,90-40,0,0);
  new TRotMatrix("rot503","rot503",90,-60,90,90-60,0,0);
  new TRotMatrix("rot504","rot504",90,-80,90,90-80,0,0);
  new TRotMatrix("rot505","rot505",90,-100,90,90-100,0,0);
  new TRotMatrix("rot506","rot506",90,-120,90,90-120,0,0);
  new TRotMatrix("rot507","rot507",90,-140,90,90-140,0,0);
  new TRotMatrix("rot508","rot508",90,-160,90,90-160,0,0);
  new TRotMatrix("rot509","rot509",90,-180,90,90-180,0,0);
  new TRotMatrix("rot510","rot510",90,-200,90,90-200,0,0);
  new TRotMatrix("rot511","rot511",90,-220,90,90-220,0,0);
  new TRotMatrix("rot512","rot512",90,-240,90,90-240,0,0);
  new TRotMatrix("rot513","rot513",90,-260,90,90-260,0,0);
  new TRotMatrix("rot514","rot514",90,-280,90,90-280,0,0);
  new TRotMatrix("rot515","rot515",90,-300,90,90-300,0,0);
  new TRotMatrix("rot516","rot516",90,-320,90,90-320,0,0);
  new TRotMatrix("rot517","rot517",90,-340,90,90-340,0,0);
  new TRotMatrix("rot518","rot518",90,-360,90,90-360,0,0);
  //
  // Position the different copies
  const Float_t rtof=(399+370)/2;
  const Int_t ntof=18;
  const Float_t angle=2*kPI/ntof;
  Float_t ang;
  //
  // Define TOF basic volume
  new TBRIK("S_TOF1","TOF box","void",130/2,29/2,190.);
  //
  // Position it
  Top->cd();
  ang=2.5*angle;
  Node = new TNode("FTO002","FTO02","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),190,"rot502");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO102","FTO102","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-190,"rot502");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=1.5*angle;
  Node = new TNode("FTO003","FTO003","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),190,"rot503");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO103","FTO103","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-190,"rot503");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();  
  ang=0.5*angle;
  Node = new TNode("FTO004","FTO004","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),190,"rot504");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node);
  //
  Top->cd();
  Node = new TNode("FTO104","FTO104","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-190,"rot504");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=-0.5*angle;
  Node = new TNode("FTO005","FTO005","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),190,"rot505");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO105","FTO105","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-190,"rot505");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=-1.5*angle;
  Node = new TNode("FTO006","FTO006","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),190,"rot506");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO106","FTO106","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-190,"rot506");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();  
  ang=kPI+1.5*angle;
  Node = new TNode("FTO012","FTO012","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),190,"rot512");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO112","FTO112","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-190,"rot512");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=kPI+0.5*angle;
  Node = new TNode("FTO013","FTO013","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),190,"rot513");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO113","FTO113","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-190,"rot513");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=kPI-0.5*angle;
  Node = new TNode("FTO014","FTO04","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),190,"rot514");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO114","FTO114","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-190,"rot514");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=kPI-1.5*angle;
  Node = new TNode("FTO015","FTO015","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),190,"rot515");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO115","FTO115","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-190,"rot515");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=kPI-2.5*angle;
  Node = new TNode("FTO016","FTO016","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),190,"rot516");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO116","FTO116","S_TOF1",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-190,"rot516");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  // Define second TOF volume
  new TBRIK("S_TOF2","TOF box","void",130/2,29/2,170.);
  //
  // Position the volume
  Top->cd();
  ang=-2.5*angle;
  Node = new TNode("FTO007","FTO007","S_TOF2",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),(2*190-170),"rot507");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO107","FTO107","S_TOF2",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-(2*190-170),"rot507");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=-3.5*angle;
  Node = new TNode("FTO008","FTO008","S_TOF2",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),(2*190-170),"rot508");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO108","FTO108","S_TOF2",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-(2*190-170),"rot508");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=-kPI/2;
  Node = new TNode("FTO009","FTO009","S_TOF2",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),(2*190-170),"rot509");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO109","FTO109","S_TOF2",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-(2*190-170),"rot509");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=kPI+3.5*angle;
  Node = new TNode("FTO010","FTO010","S_TOF2",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),(2*190-170),"rot510");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO110","FTO110","S_TOF2",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-(2*190-170),"rot510");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=kPI+2.5*angle;
  Node = new TNode("FTO011","FTO011","S_TOF2",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),(2*190-170),"rot511");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO111","FTO111","S_TOF2",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-(2*190-170),"rot511");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  // Define third TOF volume
  new TBRIK("S_TOF3","TOF box","void",130/2.,29/2,75.);
  //
  // Position it
  Top->cd();
  ang=3.5*angle;
  Node = new TNode("FTO001","FTO001","S_TOF3",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),(2*190-75),"rot501");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO101","FTO101","S_TOF3",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-(2*190-75),"rot501");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=kPI-3.5*angle;
  Node = new TNode("FTO017","FTO017","S_TOF3",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),(2*190-75),"rot517");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO117","FTO117","S_TOF3",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-(2*190-75),"rot517");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  //
  Top->cd();
  ang=kPI/2;
  Node = new TNode("FTO018","FTO018","S_TOF3",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),(2*190-75),"rot518");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
  //
  Top->cd();
  Node = new TNode("FTO118","FTO118","S_TOF3",rtof*TMath::Cos(ang),rtof*TMath::Sin(ang),-(2*190-75),"rot518");
  Node->SetLineColor(kColorTOF);
  fNodes->Add(Node); 
}

//_____________________________________________________________________________
void AliTOF::CreateGeometry()
{
  //
  // Common geometry code 
  //
  //Begin_Html
  /*
    <img src="picts/AliTOFv23.gif">
  */
  //End_Html
  //

  const Double_t kPi=TMath::Pi();
  const Double_t kDegrad=kPi/180;
  //
  Int_t lmax;
  Float_t xtof, ytof, fil_step;
  Float_t zcor1, zcor2, zcor3;
  Float_t ztof0, ztof1, ztof2;
  Float_t zl, rmin, rmax, xm, ym, dwall;
  Int_t idrotm[18];
  Float_t zm0, zm1, zm2;
  Float_t par[10];
  //
  Int_t *idtmed = fIdtmed->GetArray()-499;
  //
  // barrel iner radius 
  rmin = 370.;
  // barrel outer radius 
  rmax = rmin+29;
  // barrel length along Z axis
  zl = (rmin+2/*distance to sencetive layer*/+7/2)*2;
  //
  // frame inbetween TOF modules
  dwall = 4.;
  // Sizes of TOF module with its support etc..
  xtof = 2 * (rmin*TMath::Tan(10*kDegrad)-dwall/2-.5);
  ytof = rmax-rmin;
  ztof0 = zl/2;
  // Is it full coverage version (3) or not
  if (IsVersion() != 3) {
  ztof1 = ztof0-rmax*TMath::Tan(7.8*kDegrad); // minus Z size of PHOS
  ztof2 = ztof0-rmax*TMath::Tan(54.34/2*kDegrad); // minus Z size of HMPID;
  } else {
  ztof1 = ztof0;
  ztof2 = ztof0;
  }
   // Number of TOF-modules 
  lmax = 18;
  //
/*
  //Some imitation of TRD
  par[0] = 281;
  par[1] = 350.282;
  par[2] = zl/2;
  gMC->Gsvolu("FTRD", "TUBE", idtmed[510], par, 3);
  gMC->Gspos("FTRD", 1, "ALIC", 0., 0., 0., 0, "ONLY");

  par[0] = 0.;
  par[1] = 360.;
  par[2] = lmax;
  par[3] = 2.;
  par[4] = -zl/2;
  par[5] = rmin;
  par[6] = rmax;
  par[7] = zl/2;
  par[8] = rmin;
  par[9] = rmax;
  gMC->Gsvolu("FBAR", "PGON", idtmed[500], par, 10);
  gMC->Gspos("FBAR", 1, "ALIC", 0., 0., 0., 0, "ONLY");
*/
  //
  // TOF size  (CO2)
  par[0] = xtof / 2.;
  par[1] = ytof / 2.;
  par[2] = ztof0 / 2.;
  gMC->Gsvolu("FTO1", "BOX ", idtmed[506], par, 3);
  par[2] = ztof1 / 2.;
  gMC->Gsvolu("FTO2", "BOX ", idtmed[506], par, 3);
  par[2] = ztof2 / 2.;
  gMC->Gsvolu("FTO3", "BOX ", idtmed[506], par, 3);
/*
  // Frame wall
  par[0]=dwall/2.;
  par[1]=(rmax-rmin)/2.;
  par[2]=ztof0/2.;
  gMC->Gsvolu("FFR1", "BOX ", idtmed[508], par, 3);
  gMC->Gsatt("FFR1", "SEEN", -2);
  par[2]=ztof1/2.;
  gMC->Gsvolu("FFR2", "BOX ", idtmed[508], par, 3);
  gMC->Gsatt("FFR2", "SEEN", -2);
  par[2]=ztof2/2.;
  gMC->Gsvolu("FFR2", "BOX ", idtmed[508], par, 3);
  gMC->Gsatt("FFR2", "SEEN", -2);
*/  
  //
  // Subtraction the distanse to TOF module boundaries 
  xm = xtof -(.5 +.5)*2;
  ym = ytof;
  zm0 = ztof0;
  zm1 = ztof1;
  zm2 = ztof2;
  //  
/////////////// TOF module internal definitions //////////////
  TOFpc(xm, ym, zm0, zm1, zm2);
/////////////////////////////////////////////////////////////
  //
  // Position of modules
  fil_step = 360./lmax;
  zcor1 = ztof0/2;
  zcor2 = ztof0 - ztof1 / 2.;
  zcor3 = ztof0 - ztof2 / 2.;
/*
  for (i = 1; i <= lmax; ++i) {
    fil1 = fil_step * i;
    xcor2 = (rmin+rmax)/2 * TMath::Sin(fil1 * kDegrad);
    ycor2 = (rmin+rmax)/2 * TMath::Cos(fil1 * kDegrad);
    lmax1 = i + lmax;
    AliMatrix(idrotm[i], 90., -fil1, 90., 90. -fil1, 0., 0.);
    if (i>=7 && i<=11) { // free space for PHOS
      //    if (fil1 >= 180-50  && fil1 <= 180+50) {
      gMC->Gspos("FTO2", i, "FBAR", xcor2, ycor2, zcor2, idrotm[i], "ONLY");
      gMC->Gspos("FTO2", lmax1, "FBAR", xcor2, ycor2, -zcor2, idrotm[i], "ONLY");
    } else if (i>=17 || i==1) { // free space for RICH
      //    } else if (fil1 <= 30 || fil1 >= 360. - 30) {
      gMC->Gspos("FTO3", i, "FBAR", xcor2, ycor2, zcor3, idrotm[i], "ONLY");
      gMC->Gspos("FTO3", lmax1, "FBAR", xcor2, ycor2, -zcor3, idrotm[i], "ONLY");
    } else {
      gMC->Gspos("FTO1", i, "FBAR", xcor2, ycor2, zcor1, idrotm[i], "ONLY");
      gMC->Gspos("FTO1", lmax1, "FBAR", xcor2, ycor2, -zcor1, idrotm[i], "ONLY");
    }
  }
*/
      AliMatrix(idrotm[0], 90., 0., 0., 0., 90, -90.);
      gMC->Gspos("FTO2", 1, "BTO2", 0, zcor2, 0, idrotm[0], "ONLY");
      gMC->Gspos("FTO2", 2, "BTO2", 0, -zcor2, 0, idrotm[0], "ONLY");

      gMC->Gspos("FTO3", 1, "BTO3", 0, zcor3, 0, idrotm[0], "ONLY");
      gMC->Gspos("FTO3", 2, "BTO3", 0, -zcor3, 0, idrotm[0], "ONLY");

      gMC->Gspos("FTO1", 1, "BTO1", 0, zcor1, 0, idrotm[0], "ONLY");
      gMC->Gspos("FTO1", 2, "BTO1", 0, -zcor1, 0, idrotm[0], "ONLY");
}

//_____________________________________________________________________________
void AliTOF::DrawModule()
{
  //
  // Draw a shaded view of the common part of the TOF geometry
  //

   cout << " Drawing of AliTOF"<< endl; 
  // Set everything unseen
  gMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  gMC->Gsatt("FBAR","SEEN",0);
  gMC->Gsatt("FTO1","SEEN",1);
  gMC->Gsatt("FTO2","SEEN",1);
  gMC->Gsatt("FTO3","SEEN",1);
  gMC->Gsatt("FBT1","SEEN",1);
  gMC->Gsatt("FBT2","SEEN",1);
  gMC->Gsatt("FBT3","SEEN",1);
  gMC->Gsatt("FLT1","SEEN",1);
  gMC->Gsatt("FLT2","SEEN",1);
  gMC->Gsatt("FLT3","SEEN",1);
  //
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 1000, -1000, 1000, -1000, 1000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 40, 30, 0, 12, 9.5, .02, .02);
  gMC->Gdhead(1111, "Time Of Flight");
  gMC->Gdman(18, 4, "MAN");
  gMC->Gdopt("hide","off");
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
  //--- Quartz (SiO2) 
  Float_t   aq[2] = { 28.0855,15.9994 };
  Float_t   zq[2] = { 14.,8. };
  Float_t   wq[2] = { 1.,2. };
  Float_t   dq = 2.20;
  Int_t nq = -2;
  // --- Freon
  Float_t afre[2]  = { 12.011,18.9984032 };
  Float_t zfre[2]  = { 6.,9. };
  Float_t wfre[2]  = { 5.,12. };
  Float_t densfre  = 1.5;
  Int_t nfre = -2;
  // --- CO2 
  Float_t ac[2]   = { 12.,16. };
  Float_t zc[2]   = { 6.,8. };
  Float_t wc[2]   = { 1.,2. };
  Float_t dc = .001977;
  Int_t nc = -2;
   // For mylar (C5H4O2) 
  Float_t amy[3] = { 12., 1., 16. };
  Float_t zmy[3] = {  6., 1.,  8. };
  Float_t wmy[3] = {  5., 4.,  2. };
  Float_t dmy    = 1.39;
  Int_t nmy = -3;
 // For polyethilene (CH2) for honeycomb!!!!
  Float_t ape[2] = { 12., 1. };
  Float_t zpe[2] = {  6., 1. };
  Float_t wpe[2] = {  1., 2. };
  Float_t dpe    = 0.935*0.479; //To have 1%X0 for 1cm as for honeycomb
  Int_t npe = -2;
  // --- G10 
  Float_t ag10[4] = { 12.,1.,16.,28. };
  Float_t zg10[4] = { 6.,1.,8.,14. };
  Float_t wmatg10[4] = { .259,.288,.248,.205 };
  Float_t densg10  = 1.7;
  Int_t nlmatg10 = -4;
  // --- DME 
  Float_t adme[5] = { 12.,1.,16.,19.,79. };
  Float_t zdme[5] = { 6.,1.,8.,9.,35. };
  Float_t wmatdme[5] = { .4056,.0961,.2562,.1014,.1407 };
  Float_t densdme  = .00205;
  Int_t nlmatdme = 5;
  // ---- ALUMINA (AL203) 
  Float_t aal[2] = { 27.,16. };
  Float_t zal[2] = { 13.,8. };
  Float_t wmatal[2] = { 2.,3. };
  Float_t densal  = 2.3;
  Int_t nlmatal = -2;
  // -- Water
  Float_t awa[2] = {  1., 16. };
  Float_t zwa[2] = {  1.,  8. };
  Float_t wwa[2] = {  2.,  1. };
  Float_t dwa    = 1.0;
  Int_t nwa = -2;
  //
  //
  //AliMaterial(0, "Vacuum$", 1e-16, 1e-16, 1e-16, 1e16, 1e16);
  AliMaterial(1, "Air$",14.61,7.3,0.001205,30423.24,67500.);
  AliMaterial(2, "Cu $",  63.54, 29.0, 8.96, 1.43, 14.8);
  AliMaterial(3, "C  $",  12.01,  6.0, 2.265,18.8, 74.4);
  AliMixture(4, "Polyethilene$", ape, zpe, dpe, npe, wpe);
  AliMixture(5, "G10$", ag10, zg10, densg10, nlmatg10, wmatg10);
  AliMixture(6, "DME ", adme, zdme, densdme, nlmatdme, wmatdme);
  AliMixture(7, "CO2$", ac, zc, dc, nc, wc);
  AliMixture(8, "ALUMINA$", aal, zal, densal, nlmatal, wmatal);
  AliMaterial(9, "Al $", 26.98, 13., 2.7, 8.9, 37.2);
  // (TRD simulation) thickness = 69.282cm/18.8cm = 3.685 X/X0
  //  AliMaterial(10, "C-TRD$", 12.01, 6., 2.265*18.8/69.282*10.2/100, 18.8, 74.4); // for 10.2% 
  AliMaterial(10, "C-TRD$", 12.01, 6., 2.265*18.8/69.282*15./100, 18.8, 74.4); // for 15%
  //  AliMaterial(10, "C-TRD$", 12.01, 6., 2.265*18.8/69.282*20./100, 18.8, 74.4); // for 20%
  AliMixture(11, "Mylar$",  amy, zmy, dmy, nmy, wmy);
  AliMixture(12, "Freon$",  afre, zfre, densfre, nfre, wfre);
  AliMixture(13, "Quartz$", aq, zq, dq, nq, wq);
  AliMixture(14, "Water$",  awa, zwa, dwa, nwa, wwa);

  Float_t epsil, stmin, deemax, stemax;
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
  //  AliMedium(0, "Vacuum  $", 0, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(1, "Air$", 1, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(2, "Cu $", 2, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(3, "C  $", 3, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(4, "Pol$", 4, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(5, "G10$", 5, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(6, "DME$", 6, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(7, "CO2$", 7, 0, ISXFLD, SXMGMX, 10., -.01, -.1, .01, -.01);
  AliMedium(8, "ALUMINA$", 8, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(9, "Al Frame$", 9, 0, ISXFLD, SXMGMX, 10, stemax, deemax, epsil, stmin);
  AliMedium(10, "DME-S$", 6, 1, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(11, "C-TRD$", 10, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(12, "Myl$", 11, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(13, "Fre$", 12, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(14, "Fre-S$", 12, 1, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(15, "Glass$", 13, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
  AliMedium(16, "Water$", 14, 0, ISXFLD, SXMGMX, 10., stemax, deemax, epsil, stmin);
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
  //
  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" TOF_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  cout << "TOF version " << IsVersion() <<" initialized" << endl;
  //
  // Set id of TOF sensitive volume
  if (IsVersion() !=0) fIdSens=gMC->VolId("FPG0");
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
 
 
