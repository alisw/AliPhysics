///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Inner Traking System                                                     //
//  This class contains the base procedures for the Inner Tracking System    //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliITSClass.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>The responsible person for this module is
<a href="mailto:roberto.barbera@ct.infn.it">Roberto Barbera</a>.
</font>
<pre>
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
 
#include <TMath.h>
#include <TRandom.h>
#include <TVector.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TTUBE.h>

#include "AliITS.h"
#include "AliRun.h"

ClassImp(AliITS)
 
//_____________________________________________________________________________
AliITS::AliITS() : AliDetector()
{
  //
  // Default initialiser for ITS
  //
  fIshunt   = 0;
  fEuclidOut  =  0;
}
 
//_____________________________________________________________________________
AliITS::AliITS(const char *name, const char *title)
  : AliDetector(name,title)
{
  //
  // Default initialiser for ITS
  //
  
  fHits   = new TClonesArray("AliITShit", 1560);
  fDigits   = new TClonesArray("AliITSdigit",1000);
  
  fIshunt     =  0;
  fEuclidOut  =  0;
  
  fIdSens1 = fIdSens2 = fIdSens3 = fIdSens4 = fIdSens5 = fIdSens6 = 0;
  
  SetMarkerColor(kRed);
}

//_____________________________________________________________________________
AliITS::~AliITS()
{
  //
  // Default distructor for ITS
  //
  delete fHits;
  delete fDigits;
}

//_____________________________________________________________________________
void AliITS::AddDigit(Int_t *tracks, Int_t *digits)
{
  //
  // Add an ITS Digit
  //
  TClonesArray &ldigits = *fDigits;
  new(ldigits[fNdigits++]) AliITSdigit(tracks,digits);
}

//_____________________________________________________________________________
void AliITS::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add an ITS hit
  //
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliITShit(fIshunt,track,vol,hits);
}
 
//_____________________________________________________________________________
void AliITS::BuildGeometry()
{
  //
  // Build ITS TNODE geometry for event display
  //
  TNode *Node, *Top;
  const int kColorITS=kYellow;
  //
  Top=gAlice->GetGeometry()->GetNode("alice");

  new TTUBE("S_layer1","Layer1 of ITS","void",3.9,3.9+0.05475,12.25);
  Top->cd();
  Node = new TNode("Layer1","Layer1","S_layer1",0,0,0,"");
  Node->SetLineColor(kColorITS);
  fNodes->Add(Node);

  new TTUBE("S_layer2","Layer2 of ITS","void",7.6,7.6+0.05475,16.3);
  Top->cd();
  Node = new TNode("Layer2","Layer2","S_layer2",0,0,0,"");
  Node->SetLineColor(kColorITS);
  fNodes->Add(Node);

  new TTUBE("S_layer3","Layer3 of ITS","void",14,14+0.05288,21.1);
  Top->cd();
  Node = new TNode("Layer3","Layer3","S_layer3",0,0,0,"");
  Node->SetLineColor(kColorITS);
  fNodes->Add(Node);

  new TTUBE("S_layer4","Layer4 of ITS","void",24,24+0.05288,29.6);
  Top->cd();
  Node = new TNode("Layer4","Layer4","S_layer4",0,0,0,"");
  Node->SetLineColor(kColorITS);
  fNodes->Add(Node);

  new TTUBE("S_layer5","Layer5 of ITS","void",40,40+0.05382,45.1);
  Top->cd();
  Node = new TNode("Layer5","Layer5","S_layer5",0,0,0,"");
  Node->SetLineColor(kColorITS);
  fNodes->Add(Node);

  new TTUBE("S_layer6","Layer6 of ITS","void",45,45+0.05382,50.4);
  Top->cd();
  Node = new TNode("Layer6","Layer6","S_layer6",0,0,0,"");
  Node->SetLineColor(kColorITS);
  fNodes->Add(Node);
}
 
//_____________________________________________________________________________
void AliITS::CreateMaterials()
{
  //
  // Create ITS materials
  //

  Float_t awat[2]  = { 1.00794,15.9994 };
  Float_t zwat[2]  = { 1.,8. };
  Float_t wwat[2]  = { 2.,1. };
  Float_t denswat  = 1.;
  Float_t afre[2]  = { 12.011,18.9984032 };
  Float_t zfre[2]  = { 6.,9. };
  Float_t wfre[2]  = { 5.,12. };
  Float_t densfre  = 1.5;
  //     94.4% Al2O3 , 2.8% SiO2 , 2.3% MnO , 0.5% Cr2O3 
  Float_t acer[5]  = { 26.981539,15.9994,28.0855,54.93805,51.9961 };
  Float_t zcer[5]  = { 13.,8.,14.,25.,	    24. };
  Float_t wcer[5]  = { .49976,1.01233,.01307,	    .01782,.00342 };
  Float_t denscer  = 3.6;
  //     60% SiO2 , 40% G10FR4 
  Float_t apcb[3]  = { 28.0855,15.9994,17.749 };
  Float_t zpcb[3]  = { 14.,8.,8.875 };
  Float_t wpcb[3]  = { .28,.32,.4 };
  Float_t denspcb  = 1.8;
  Float_t apoly[2] = { 12.01,1. };
  Float_t zpoly[2] = { 6.,1. };
  Float_t wpoly[2] = { .33,.67 };
  Float_t zserv[4] = { 1.,6.,26.,29. };
  Float_t aserv[4] = { 1.,12.,55.8,63.5 };
  Float_t wserv[4] = { .014,.086,.42,.48 };
  
  Int_t ISXFLD   = gAlice->Field()->Integ();
  Float_t SXMGMX = gAlice->Field()->Max();
  
  
  // --- Define the various materials for GEANT --- 
  
  //    200-224 --> Silicon Pixel Detectors (detectors, chips, buses, cooling,..)
  
  AliMaterial(0, "SPD Si$",      28.0855, 14., 2.33, 9.36, 999);
  AliMaterial(1, "SPD Si chip$", 28.0855, 14., 2.33, 9.36, 999);
  AliMaterial(2, "SPD Si bus$",  28.0855, 14., 2.33, 9.36, 999);
  AliMaterial(3, "SPD C$",       12.011,   6., 2.265,18.8, 999);
  // v. dens 
  AliMaterial(4, "SPD Air$",    14.61, 7.3, .001205, 30423., 999);
  AliMaterial(5, "SPD Vacuum$", 1e-16, 1e-16, 1e-16, 1e16, 1e16);
  AliMaterial(6, "SPD Al$",     26.981539, 13., 2.6989, 8.9, 999);
  AliMixture( 7, "SPD Water $", awat, zwat, denswat, -2, wwat);
  AliMixture( 8, "SPD Freon$",  afre, zfre, densfre, -2, wfre);
  // ** 
  AliMedium(200, "SPD Si$",      0, 1, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(201, "SPD Si chip$", 1, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(202, "SPD Si bus$",  2, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(203, "SPD C$",       3, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(204, "SPD Air$",     4, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(205, "SPD Vacuum$",  5, 0, ISXFLD, SXMGMX, 10.,  1., .1, .1,    10.);
  AliMedium(206, "SPD Al$",      6, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(207, "SPD Water $",  7, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(208, "SPD Freon$",   8, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  
  //    225-249 --> Silicon Drift Detectors (detectors, chips, buses, cooling,..)
  
  AliMaterial(25, "SDD Si$",      28.0855, 14., 2.33,  9.36, 999);
  AliMaterial(26, "SDD Si chip$", 28.0855, 14., 2.33,  9.36, 999);
  AliMaterial(27, "SDD Si bus$",  28.0855, 14., 2.33,  9.36, 999);
  AliMaterial(28, "SDD C$",       12.011,   6., 2.265,18.8,  999);
  // v. dens 
  AliMaterial(29, "SDD Air$",     14.61, 7.3, .001205, 30423., 999);
  AliMaterial(30, "SDD Vacuum$",  1e-16, 1e-16, 1e-16, 1e16,  1e16);
  AliMaterial(31, "SDD Al$",      26.981539, 13., 2.6989, 8.9, 999);
  // After a call with ratios by number (negative number of elements), 
  // the ratio array is changed to the ratio by weight, so all successive 
  // calls with the same array must specify the number of elements as 
  // positive 
  AliMixture(32, "SDD Water $", awat, zwat, denswat, 2, wwat);
  // After a call with ratios by number (negative number of elements), 
  // the ratio array is changed to the ratio by weight, so all successive 
  // calls with the same array must specify the number of elements as 
  // positive 
  AliMixture( 33, "SDD Freon$", afre, zfre, densfre, 2, wfre);
  AliMixture( 34, "SDD PCB$",   apcb, zpcb, denspcb, 3, wpcb);
  AliMaterial(35, "SDD Copper$", 63.546, 29., 8.96, 1.43, 999);
  AliMixture( 36, "SDD Ceramics$", acer, zcer, denscer, -5, wcer);
  AliMaterial(37, "SDD Kapton$", 12.011, 6., 1.3, 31.27, 999);
  // ** 
  // check A and Z 
  AliMedium(225, "SDD Si$",      25, 1, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(226, "SDD Si chip$", 26, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(227, "SDD Si bus$",  27, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(228, "SDD C$",       28, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(229, "SDD Air$",     29, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(230, "SDD Vacuum$",  30, 0, ISXFLD, SXMGMX, 10.,  1., .1, .1,    10.);
  AliMedium(231, "SDD Al$",      31, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(232, "SDD Water $",  32, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(233, "SDD Freon$",   33, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(234, "SDD PCB$",     34, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(235, "SDD Copper$",  35, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(236, "SDD Ceramics$",36, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(237, "SDD Kapton$",  37, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  
  //    250-274 --> Silicon Strip Detectors (detectors, chips, buses, cooling,..)
  
  AliMaterial(50, "SSD Si$",      28.0855, 14., 2.33, 9.36, 999.);
  AliMaterial(51, "SSD Si chip$", 28.0855, 14., 2.33, 9.36, 999.);
  AliMaterial(52, "SSD Si bus$",  28.0855, 14., 2.33, 9.36, 999.);
  AliMaterial(53, "SSD C$",       12.011,   6., 2.265,18.8, 999.);
  // v. dens 
  AliMaterial(54, "SSD Air$",     14.61, 7.3, .001205, 30423., 999);
  AliMaterial(55, "SSD Vacuum$",  1e-16, 1e-16, 1e-16, 1e16, 1e16);
  AliMaterial(56, "SSD Al$",      26.981539, 13., 2.6989, 8.9, 999);
  // After a call with ratios by number (negative number of elements), 
  // the ratio array is changed to the ratio by weight, so all successive 
  // calls with the same array must specify the number of elements as 
  // positive 
  AliMixture(57, "SSD Water $", awat, zwat, denswat, 2, wwat);
  // After a call with ratios by number (negative number of elements), 
  // the ratio array is changed to the ratio by weight, so all successive 
  // calls with the same array must specify the number of elements as 
  // positive 
  AliMixture(58, "SSD Freon$", afre, zfre, densfre, 2, wfre);
  AliMixture(59, "SSD PCB$",   apcb, zpcb, denspcb, 3, wpcb);
  AliMaterial(60, "SSD Copper$", 63.546, 29., 8.96, 1.43, 999.);
  // After a call with ratios by number (negative number of elements), 
  // the ratio array is changed to the ratio by weight, so all successive 
  // calls with the same array must specify the number of elements as 
  // positive 
  AliMixture( 61, "SSD Ceramics$", acer, zcer, denscer, 5, wcer);
  AliMaterial(62, "SSD Kapton$", 12.011, 6., 1.3, 31.27, 999.);
  // check A and Z 
  AliMaterial(63, "SDD G10FR4$", 17.749, 8.875, 1.8, 21.822, 999.);
  // ** 
  AliMedium(250, "SSD Si$",      50, 1, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(251, "SSD Si chip$", 51, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(252, "SSD Si bus$",  52, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(253, "SSD C$",       53, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(254, "SSD Air$",     54, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(255, "SSD Vacuum$",  55, 0, ISXFLD, SXMGMX, 10.,  1., .1, .1,    10.);
  AliMedium(256, "SSD Al$",      56, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(257, "SSD Water $",  57, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(258, "SSD Freon$",   58, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(259, "SSD PCB$",     59, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(260, "SSD Copper$",  60, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(261, "SSD Ceramics$",61, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(262, "SSD Kapton$",  62, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(263, "SSD G10FR4$",  63, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  
  //     275-299 --> General (end-caps, frames, cooling, cables, etc.) 
  
  AliMaterial(75, "GEN C$", 12.011, 6., 2.265, 18.8, 999.);
  // verify density 
  AliMaterial(76, "GEN Air$", 14.61, 7.3, .001205, 30423., 999);
  AliMaterial(77, "GEN Vacuum$", 1e-16, 1e-16, 1e-16, 1e16, 1e16);
  AliMixture( 78, "GEN POLYETHYL$", apoly, zpoly, .95, -2, wpoly);
  AliMixture( 79, "GEN SERVICES$",  aserv, zserv, 4.68, 4, wserv);
  AliMaterial(80, "GEN Copper$", 63.546, 29., 8.96, 1.43, 999.);
  // After a call with ratios by number (negative number of elements), 
  // the ratio array is changed to the ratio by weight, so all successive 
  // calls with the same array must specify the number of elements as 
  // positive 
  AliMixture(81, "GEN Water $", awat, zwat, denswat, 2, wwat);
  // ** 
  AliMedium(275, "GEN C$",         75, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(276, "GEN Air$",       76, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(277, "GEN Vacuum$",    77, 0, ISXFLD, SXMGMX, 10., .1,  .1, .1,    10.);
  AliMedium(278, "GEN POLYETHYL$", 78, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(279, "GEN SERVICES$",  79, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(280, "GEN Copper$",    80, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(281, "GEN Water $",    81, 0, ISXFLD, SXMGMX, 10., .01, .1, .003, .003);
}

//_____________________________________________________________________________
Int_t AliITS::DistancetoPrimitive(Int_t , Int_t )
{
  //
  // Distance from mouse to ITS on the screen. Dummy routine
  //
  return 9999;
}

//_____________________________________________________________________________
void AliITS::Init()
{
  //
  // Initialise ITS after it has been built
  //
  Int_t i;
  AliMC* pMC = AliMC::GetMC();
  //
  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" ITS_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  //
  //
  fIdSens1=pMC->VolId("ITS1");
  fIdSens2=pMC->VolId("ITS2");
  fIdSens3=pMC->VolId("ITS3");
  fIdSens4=pMC->VolId("ITS4");
  fIdSens5=pMC->VolId("ITS5");
  fIdSens6=pMC->VolId("ITS6");
  //
  for(i=0;i<80;i++) printf("*");
  printf("\n");
}

//_____________________________________________________________________________
void AliITS::MakeBranch(Option_t* option)
{
  //
  // Create Tree branches for the ITS.
  //
  Int_t buffersize = 4000;
  char branchname[10];
  sprintf(branchname,"%s",GetName());

  AliDetector::MakeBranch(option);

  char *D = strstr(option,"D");

  if (fDigits   && gAlice->TreeD() && D) {
    gAlice->TreeD()->Branch(branchname,&fDigits, buffersize);
    printf("Making Branch %s for digits\n",branchname);
  }	
}
 
//_____________________________________________________________________________
void AliITS::SetEUCLID(Bool_t euclid)
{
  //
  // set euclid output flag
  //
  fEuclidOut=euclid;
}
 
ClassImp(AliITShit)
 
//_____________________________________________________________________________
AliITShit::AliITShit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
  AliHit(shunt, track)
{
  //
  // Create ITS hit
  //
  fLayer      = vol[0];   // Layer number
  fLadder     = vol[2];   // Ladder number
  fDet        = vol[1];   // Detector number   
  fX          = hits[0];
  fY          = hits[1];
  fZ          = hits[2];
  fPx         = hits[3];
  fPy         = hits[4];
  fPz         = hits[5];
  fDestep     = hits[6];
}

ClassImp(AliITSdigit)
 
//_____________________________________________________________________________
AliITSdigit::AliITSdigit(Int_t *tracks, Int_t *digits):
  AliDigit(tracks)
{
  //
  // Create ITS digit
  //
  fEvent      = digits[0];
  fLayer      = digits[1];
  fDet        = digits[2];
  fNoverl     = digits[3];
}
