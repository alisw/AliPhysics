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
Revision 1.8  1999/09/29 09:24:19  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//
//      An overview of the basic philosophy of the ITS code development
// and analysis is show in the figure below.
//Begin_Html
/*
<img src="picts/ITS/ITS_Analysis_schema.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>Roberto Barbera is in charge of the ITS Offline code (1999).
<a href="mailto:roberto.barbera@ct.infn.it">Roberto Barbera</a>.
</font>
<pre>
*/
//End_Html
//
//  AliITS. Inner Traking System base class.
//  This class contains the base procedures for the Inner Tracking System
//
//Begin_Html
/*
<img src="picts/ITS/AliITS_Class_Diagram.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This show the class diagram of the different elements that are part of
the AliITS class.
</font>
<pre>
*/
//End_Html
//
// Version: 0
// Written by Rene Brun, Federico Carminati, and Roberto Barbera
//
// Version: 1
// Modified and documented by Bjorn S. Nilsen
// July 11 1999
//
// AliITS is the general base class for the ITS. Also see AliDetector for
// futher information.
//
///////////////////////////////////////////////////////////////////////////////
 
#include <TMath.h>
#include <TRandom.h>
#include <TVector.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TTUBE.h>

#include "AliITSmodule.h"
#include "AliDetector.h"
#include "AliITS.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "AliITShit.h"
#include "AliITSdigit.h"
#include "AliRun.h"

ClassImp(AliITS)
 
//_____________________________________________________________________________
AliITS::AliITS() {
  //
  // Default initialiser for ITS
  //     The default constructor of the AliITS class. In addition to
  // creating the AliITS class it zeros the variables fIshunt (a member
  // of AliDetector class), fEuclidOut, and fIdN, and zeros the pointers
  // fITSpoints, fIdSens, and fIdName.
  //
  fITSpoints  = 0;
  fIshunt     = 0;
  fEuclidOut  = 0;
  fIdN        = 0;
  fIdName     = 0;
  fIdSens     = 0;
  fITSmodules = 0;

}
//_____________________________________________________________________________
AliITS::AliITS(const char *name, const char *title):AliDetector(name,title){
  //
  // Default initialiser for ITS
  //     The constructor of the AliITS class. In addition to creating the
  // AliITS class, it allocates memory for the TClonesArrays fHits and
  // fDigits, and for the TObjArray fITSpoints. It also zeros the variables
  // fIshunt (a member of AliDetector class), fEuclidOut, and fIdN, and zeros
  // the pointers fIdSens and fIdName. To help in displaying hits via the ROOT
  // macro display.C AliITS also sets the marker color to red. The variables
  // passes with this constructor, const char *name and *title, are used by
  // the constructor of AliDetector class. See AliDetector class for a
  // description of these parameters and its constructor functions.
  //

  fHits       = new TClonesArray("AliITShit", 1560);
  fDigits     = new TClonesArray("AliITSdigit",1000);
  fITSpoints  = new TObjArray();
  fITSmodules = 0; //new AliITSmodules();

  fIshunt     = 0;
  fEuclidOut  = 0;
  fIdN        = 0;
  fIdName     = 0;
  fIdSens     = 0;

  SetMarkerColor(kRed);

}

//_____________________________________________________________________________
AliITS::~AliITS(){
  //
  // Default distructor for ITS
  //     The default destructor of the AliITS class. In addition to deleting
  // the AliITS class it deletes the memory pointed to by the fHits, fDigits,
  // fIdSens, fIdName, and fITSpoints.
  //
  delete fHits;
  delete fDigits;
  if(fIdName!=0) delete[] fIdName;
  if(fIdSens!=0) delete[] fIdSens;
  delete fITSmodules;
  if(fITSpoints!=0) delete fITSpoints;
}

//_____________________________________________________________________________
void AliITS::AddDigit(Int_t *tracks, Int_t *digits){
  //
  // Add an ITS Digit
  //     The function to add information to the AliITSdigits class. See the
  // AliITSdigits class for a full description. This function allocates the
  // necessary new space for the digits information and passes the pointers
  // *track and *digits to the AliITSdigits constructor function.
  //
  TClonesArray &ldigits = *fDigits;
  new(ldigits[fNdigits++]) AliITSdigit(tracks,digits);
}

Int_t AliITS::AddDigit(AliITSdigit* d) {

   fDigits->Add(d);
   fNdigits = fDigits->GetEntriesFast();
   return fNdigits;
}

//_____________________________________________________________________________
void AliITS::AddHit(Int_t track, Int_t *vol, Float_t *hits){
  //
  // Add an ITS hit
  //     The function to add information to the AliITShit class. See the
  // AliITShit class for a full description. This function allocates the
  // necessary new space for the hit information and passes the variable
  // track, and the pointers *vol and *hits to the AliITShit constructor
  // function.
  //
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliITShit(fIshunt,track,vol,hits);
}
//_____________________________________________________________________________
void AliITS::BuildGeometry(){
  //
  // Build ITS TNODE geometry for event display
  //     This function builds a simple ITS geometry used by the ROOT macro
  // display.C. In general the geometry as coded is wrong.
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
void AliITS::CreateMaterials(){
  //
  // Create ITS materials
  //     This function defines the default materials used in the Geant
  // Monte Carlo simulations. In general it is automatically replaced by
  // the CreatMaterials routine defined in AliITSv?. Should the function
  // CreateMaterials not exist for the geometry version you are using this
  // one is used. See the definition found in AliITSv5 or the other routine
  // for a complete definition.
  //
  // Water H2O
  Float_t awat[2]  = { 1.00794,15.9994 };
  Float_t zwat[2]  = { 1.,8. };
  Float_t wwat[2]  = { 2.,1. };
  Float_t denswat  = 1.;
  // Freon
  Float_t afre[2]  = { 12.011,18.9984032 };
  Float_t zfre[2]  = { 6.,9. };
  Float_t wfre[2]  = { 5.,12. };
  Float_t densfre  = 1.5;
  // Ceramics
  //     94.4% Al2O3 , 2.8% SiO2 , 2.3% MnO , 0.5% Cr2O3 
  Float_t acer[5]  = { 26.981539,15.9994,28.0855,54.93805,51.9961 };
  Float_t zcer[5]  = { 13.,8.,14.,25.,	    24. };
  Float_t wcer[5]  = { .49976,1.01233,.01307,	    .01782,.00342 };
  Float_t denscer  = 3.6;
  //
  //     60% SiO2 , 40% G10FR4 
  // PC board
  Float_t apcb[3]  = { 28.0855,15.9994,17.749 };
  Float_t zpcb[3]  = { 14.,8.,8.875 };
  Float_t wpcb[3]  = { .28,.32,.4 };
  Float_t denspcb  = 1.8;
  // POLYETHYL
  Float_t apoly[2] = { 12.01,1. };
  Float_t zpoly[2] = { 6.,1. };
  Float_t wpoly[2] = { .33,.67 };
  // SERVICES
  Float_t zserv[4] = { 1.,6.,26.,29. };
  Float_t aserv[4] = { 1.,12.,55.8,63.5 };
  Float_t wserv[4] = { .014,.086,.42,.48 };
  
  Int_t  ISXFLD  = gAlice->Field()->Integ();
  Float_t SXMGMX = gAlice->Field()->Max();
  
  
  // --- Define the various materials for GEANT --- 
  
  //  200-224 --> Silicon Pixel Detectors (detectors, chips, buses, cooling,..)
  
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
  AliMedium(0, "SPD Si$",      0, 1,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(1, "SPD Si chip$", 1, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(2, "SPD Si bus$",  2, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(3, "SPD C$",       3, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(4, "SPD Air$",     4, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(5, "SPD Vacuum$",  5, 0,ISXFLD,SXMGMX, 10.,1.00, .1, .100,10.00);
  AliMedium(6, "SPD Al$",      6, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(7, "SPD Water $",  7, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(8, "SPD Freon$",   8, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  
  //  225-249 --> Silicon Drift Detectors (detectors, chips, buses, cooling,..)
  
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
  AliMedium(25, "SDD Si$",      25, 1,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(26, "SDD Si chip$", 26, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(27, "SDD Si bus$",  27, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(28, "SDD C$",       28, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(29, "SDD Air$",     29, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(30, "SDD Vacuum$",  30, 0,ISXFLD,SXMGMX, 10.,1.00, .1, .100,10.00);
  AliMedium(31, "SDD Al$",      31, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(32, "SDD Water $",  32, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(33, "SDD Freon$",   33, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(34, "SDD PCB$",     34, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(35, "SDD Copper$",  35, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(36, "SDD Ceramics$",36, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(37, "SDD Kapton$",  37, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  
  //  250-274 --> Silicon Strip Detectors (detectors, chips, buses, cooling,..)
  
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
  AliMedium(50, "SSD Si$",      50, 1,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(51, "SSD Si chip$", 51, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(52, "SSD Si bus$",  52, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(53, "SSD C$",       53, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(54, "SSD Air$",     54, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(55, "SSD Vacuum$",  55, 0,ISXFLD,SXMGMX, 10.,1.00, .1, .100,10.00);
  AliMedium(56, "SSD Al$",      56, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(57, "SSD Water $",  57, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(58, "SSD Freon$",   58, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(59, "SSD PCB$",     59, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(60, "SSD Copper$",  60, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(61, "SSD Ceramics$",61, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(62, "SSD Kapton$",  62, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(63, "SSD G10FR4$",  63, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  
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
  AliMedium(75,"GEN C$",        75, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(76,"GEN Air$",      76, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(77,"GEN Vacuum$",   77, 0,ISXFLD,SXMGMX, 10., .10, .1, .100,10.00);
  AliMedium(78,"GEN POLYETHYL$",78, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(79,"GEN SERVICES$", 79, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(80,"GEN Copper$",   80, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
  AliMedium(81,"GEN Water $",   81, 0,ISXFLD,SXMGMX, 10., .01, .1, .003, .003);
}

//_____________________________________________________________________________
Int_t AliITS::DistancetoPrimitive(Int_t , Int_t ){
  //
  // Distance from mouse to ITS on the screen. Dummy routine
  //     A dummy routine used by the ROOT macro display.C to allow for the
  // use of the mouse (pointing device) in the macro. In general this should
  // never be called. If it is it returns the number 9999 for any value of
  // x and y.
  //
  return 9999;
}

//_____________________________________________________________________________
void AliITS::Init(){
  //
  // Initialise ITS after it has been built
  //     This routine initializes the AliITS class. It is intended to be called
  // from the Init function in AliITSv?. Besides displaying a banner
  // indicating that it has been called it initializes the array fIdSens.
  // Therefore it should be called after a call to CreateGeometry.
  //
  Int_t i;
  //
  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" ITS_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  //
  //
  for(i=0;i<fIdN;i++) fIdSens[i] = gMC->VolId(fIdName[i]);
  //
  for(i=0;i<80;i++) printf("*");
  printf("\n");
}

//_____________________________________________________________________________
void AliITS::MakeBranch(Option_t* option){
  //
  // Create Tree branches for the ITS.
  // Creates the TTree branch where the class AliITS is kept.
  //
  Int_t buffersize = 4000;
  char branchname[10];
  sprintf(branchname,"%s",GetName());

  AliDetector::MakeBranch(option);

  char *D = strstr(option,"D");

  if (fDigits   && gAlice->TreeD() && D) {
    gAlice->TreeD()->Branch(branchname,&fDigits, buffersize);
    printf("Making Branch %s for digits\n",branchname);
  } // end if
}

//____________________________________________________________________________
void AliITS::Streamer(TBuffer &R__b){
   // Stream an object of class AliITS.
    Int_t i,j,l;

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); 
      if (R__v == 1) {
	  AliDetector::Streamer(R__b);
	  R__b >> fITSgeom;
//        R__b >> fITSmodules; //We do not write out modules so don't read them
	  R__b >> fITSpoints;
	  R__b >> fEuclidOut;
	  R__b >> fIdN;
	  if(fIdSens!=0) delete[] fIdSens;
	  if(fIdName!=0) delete[] fIdName;
	  fIdSens = new Int_t[fIdN];
	  fIdName = new char*[fIdN];
	  for(i=0;i<fIdN;i++) R__b >> fIdSens[i];
	  for(i=0;i<fIdN;i++){
	      R__b >> l;
	      fIdName[i] = new char[l+1]; // add room for null character.
	      for(j=0;j<l;j++) R__b >> fIdName[i][j];
	      fIdName[i][l] = '\0'; // Null terminate this string.
	  } // end for i
	  R__b >> fMajorVersion;
	  R__b >> fMinorVersion;
      } // end if (R__v)
   } else {
      R__b.WriteVersion(AliITS::IsA());
      AliDetector::Streamer(R__b);
      R__b << fITSgeom;
//    R__b << fITSmodules; //We don't want to write out the modules class.
      R__b << fITSpoints;
      R__b << fEuclidOut;
      R__b << fIdN;
      for(i=0;i<fIdN;i++) R__b <<fIdSens[i];
      for(i=0;i<fIdN;i++){
	  l = strlen(fIdName[i]);
	  R__b << l;
	  for(j=0;j<l;j++) R__b << fIdName[i][j];
      } // end for i
      R__b << fMajorVersion;
      R__b << fMinorVersion;
   }
}
