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
Revision 1.20  2000/06/08 18:32:57  cblume
Make code compliant to coding conventions

Revision 1.19  2000/06/07 16:25:37  cblume
Try to remove compiler warnings on Sun and HP

Revision 1.18  2000/05/08 16:17:27  cblume
Merge TRD-develop

Revision 1.17.2.1  2000/05/08 14:28:59  cblume
Introduced SetPHOShole() and SetRICHhole(). AliTRDrecPoint container is now a TObjArray

Revision 1.17  2000/02/28 19:10:26  cblume
Include the new TRD classes

Revision 1.16.2.2  2000/02/28 17:53:24  cblume
Introduce TRD geometry classes

Revision 1.16.2.1  2000/02/28 17:04:19  cblume
Include functions and data members for AliTRDrecPoint

Revision 1.16  2000/01/19 17:17:35  fca
Introducing a list of lists of hits -- more hits allowed for detector now

Revision 1.15  1999/11/02 17:04:25  fca
Small syntax change for HP compiler

Revision 1.14  1999/11/02 16:57:02  fca
Avoid non ansi warnings on HP compilers

Revision 1.13  1999/11/02 16:35:56  fca
New version of TRD introduced

Revision 1.12  1999/11/01 20:41:51  fca
Added protections against using the wrong version of FRAME

Revision 1.11  1999/09/29 09:24:34  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Transition Radiation Detector                                            //
//  This class contains the basic functions for the Transition Radiation     //
//  Detector.                                                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>

#include <TMath.h>
#include <TNode.h>
#include <TPGON.h> 

#include "AliTRD.h"
#include "AliRun.h"
#include "AliConst.h"
#include "AliTRDdigitizer.h"
#include "AliTRDclusterizer.h"
#include "AliTRDgeometryHole.h"
#include "AliTRDgeometryFull.h"
#include "AliTRDrecPoint.h"
 
ClassImp(AliTRD)
 
//_____________________________________________________________________________
AliTRD::AliTRD()
{
  //
  // Default constructor
  //

  fIshunt      = 0;
  fGasMix      = 0;
  fHits        = 0;
  fDigits      = 0;

  fRecPoints   = 0;
  fNRecPoints  = 0;

  fGeometry    = 0;

}
 
//_____________________________________________________________________________
AliTRD::AliTRD(const char *name, const char *title)
       : AliDetector(name,title)
{
  //
  // Standard constructor for the TRD
  //

  // Check that FRAME is there otherwise we have no place where to
  // put TRD
  AliModule* frame = gAlice->GetModule("FRAME");
  if (!frame) {
    Error("Ctor","TRD needs FRAME to be present\n");
    exit(1);
  } 

  // Define the TRD geometry according to the FRAME geometry
  if      (frame->IsVersion() == 0) {
    // Geometry with hole
    fGeometry = new AliTRDgeometryHole();
  }
  else if (frame->IsVersion() == 1) {
    // Geometry without hole
    fGeometry = new AliTRDgeometryFull();
  }
  else {
    Error("Ctor","Could not find valid FRAME version\n");
    exit(1);
  }

  // Allocate the hit array
  fHits       = new TClonesArray("AliTRDhit"     ,405);
  gAlice->AddHitList(fHits);

  // Allocate the digits array
  fDigits     = 0;

  // Allocate the rec point array
  fRecPoints  = new TObjArray(400);
  fNRecPoints = 0;
   
  fIshunt = 0;
  fGasMix = 0;

  SetMarkerColor(kWhite);   

}

//_____________________________________________________________________________
AliTRD::AliTRD(const AliTRD &trd)
{
  //
  // Copy constructor
  //

  ((AliTRD &) trd).Copy(*this);

}

//_____________________________________________________________________________
AliTRD::~AliTRD()
{
  //
  // TRD destructor
  //

  fIshunt = 0;

  delete fGeometry;
  delete fHits;
  delete fRecPoints;

}

//_____________________________________________________________________________
void AliTRD::AddRecPoint(Float_t *pos, Int_t *digits, Int_t det, Float_t amp)
{
  //
  // Add a reconstructed point for the TRD
  //
  
  AliTRDrecPoint *recPoint = new AliTRDrecPoint();
  TVector3        posVec(pos[0],pos[1],pos[2]);
  recPoint->SetLocalPosition(posVec);
  recPoint->SetDetector(det);
  recPoint->SetEnergy(amp);
  for (Int_t iDigit = 0; iDigit < 3; iDigit++) {
    recPoint->AddDigit(digits[iDigit]);
  }

  fRecPoints->Add(recPoint);

}

//_____________________________________________________________________________
void AliTRD::AddDigit(Int_t *digits, Int_t *amp)
{
  //
  // Add a digit for the TRD
  //

  TClonesArray &ldigits = *fDigits;
  new(ldigits[fNdigits++]) AliTRDdigit(kFALSE,digits,amp);

}

//_____________________________________________________________________________
void AliTRD::AddHit(Int_t track, Int_t *det, Float_t *hits)
{
  //
  // Add a hit for the TRD
  //

  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliTRDhit(fIshunt,track,det,hits);

}

//_____________________________________________________________________________
void AliTRD::BuildGeometry()
{
  //
  // Create the ROOT TNode geometry for the TRD
  //

  TNode *node, *top;
  TPGON *pgon;
  const Int_t kColorTRD = 46;
  
  // Find the top node alice
  top = gAlice->GetGeometry()->GetNode("alice");
  
  pgon = new TPGON("S_TRD","TRD","void",0,360,kNsect,4);
  Float_t ff    = TMath::Cos(kDegrad * 180 / kNsect);
  Float_t rrmin = kRmin / ff;
  Float_t rrmax = kRmax / ff;
  pgon->DefineSection(0,-kZmax1,rrmax,rrmax);
  pgon->DefineSection(1,-kZmax2,rrmin,rrmax);
  pgon->DefineSection(2, kZmax2,rrmin,rrmax);
  pgon->DefineSection(3, kZmax1,rrmax,rrmax);
  top->cd();
  node = new TNode("TRD","TRD","S_TRD",0,0,0,"");
  node->SetLineColor(kColorTRD);
  fNodes->Add(node);

}
 
//_____________________________________________________________________________
void AliTRD::Copy(TObject &trd)
{
  //
  // Copy function
  //

  ((AliTRD &) trd).fGasMix     = fGasMix;
  ((AliTRD &) trd).fGeometry   = fGeometry;
  ((AliTRD &) trd).fRecPoints  = fRecPoints;
  ((AliTRD &) trd).fNRecPoints = fNRecPoints;

  //AliDetector::Copy(trd);

}

//_____________________________________________________________________________
void AliTRD::CreateGeometry()
{
  //
  // Creates the volumes for the TRD chambers
  //

  // Check that FRAME is there otherwise we have no place where to put the TRD
  AliModule* frame = gAlice->GetModule("FRAME");
  if (!frame) {
    printf(" The TRD needs the FRAME to be defined first\n");
    return;
  }

  fGeometry->CreateGeometry(fIdtmed->GetArray() - 1299);

}
 
//_____________________________________________________________________________
void AliTRD::CreateMaterials()
{
  //
  // Create the materials for the TRD
  // Origin Y.Foka
  //

  Int_t   isxfld = gAlice->Field()->Integ();
  Float_t sxmgmx = gAlice->Field()->Max();
  
  // For polyethilene (CH2) 
  Float_t ape[2] = { 12., 1. };
  Float_t zpe[2] = {  6., 1. };
  Float_t wpe[2] = {  1., 2. };
  Float_t dpe    = 0.95;

  // For mylar (C5H4O2) 
  Float_t amy[3] = { 12., 1., 16. };
  Float_t zmy[3] = {  6., 1.,  8. };
  Float_t wmy[3] = {  5., 4.,  2. };
  Float_t dmy    = 1.39;

  // For CO2 
  Float_t aco[2] = { 12., 16. };
  Float_t zco[2] = {  6.,  8. };
  Float_t wco[2] = {  1.,  2. };
  Float_t dco    = 0.001977;

  // For water
  Float_t awa[2] = {  1., 16. };
  Float_t zwa[2] = {  1.,  8. };
  Float_t wwa[2] = {  2.,  1. };
  Float_t dwa    = 1.0;

  // For isobutane (C4H10)
  Float_t ais[2] = { 12.,  1. };
  Float_t zis[2] = {  6.,  1. };
  Float_t wis[2] = {  4., 10. };
  Float_t dis    = 0.00267;

  // For Xe/CO2-gas-mixture 
  // Xe-content of the Xe/CO2-mixture (90% / 10%) 
  Float_t fxc    = .90;
  // Xe-content of the Xe/Isobutane-mixture (97% / 3%) 
  Float_t fxi    = .97;
  Float_t dxe    = .005858;
  
  // General tracking parameter
  Float_t tmaxfd = -10.;
  Float_t stemax = -1e10;
  Float_t deemax = -0.1;
  Float_t epsil  =  1e-4;
  Float_t stmin  = -0.001;
  
  Float_t absl, radl, d, buf[1];
  Float_t agm[2], dgm, zgm[2], wgm[2];
  Int_t   nbuf;
  
  //////////////////////////////////////////////////////////////////////////
  //     Define Materials 
  //////////////////////////////////////////////////////////////////////////

  AliMaterial( 1, "Al $",  26.98, 13.0, 2.7     ,     8.9 ,    37.2);
  AliMaterial( 2, "Air$",  14.61,  7.3, 0.001205, 30420.0 , 67500.0);
  AliMaterial( 4, "Xe $", 131.29, 54.0, dxe     ,  1447.59,     0.0);
  AliMaterial( 5, "Cu $",  63.54, 29.0, 8.96    ,     1.43,    14.8);
  AliMaterial( 6, "C  $",  12.01,  6.0, 2.265   ,    18.8 ,    74.4);
  AliMaterial(12, "G10$",  20.00, 10.0, 1.7     ,    19.4 ,   999.0);

  // Mixtures 
  AliMixture(3, "Polyethilene$",   ape, zpe, dpe, -2, wpe);
  AliMixture(7, "Mylar$",          amy, zmy, dmy, -3, wmy);
  AliMixture(8, "CO2$",            aco, zco, dco, -2, wco);
  AliMixture(9, "Isobutane$",      ais, zis, dis, -2, wis);
  AliMixture(13,"Water$",          awa, zwa, dwa, -2, wwa);

  // Gas mixtures
  Char_t namate[21];
  // Xe/CO2-mixture
  // Get properties of Xe 
  gMC->Gfmate((*fIdmate)[4], namate, agm[0], zgm[0], d, radl, absl, buf, nbuf);
  // Get properties of CO2 
  gMC->Gfmate((*fIdmate)[8], namate, agm[1], zgm[1], d, radl, absl, buf, nbuf);
  // Create gas mixture 
  wgm[0] = fxc;
  wgm[1] = 1. - fxc;
  dgm    = wgm[0] * dxe + wgm[1] * dco;
  AliMixture(10, "Gas mixture 1$", agm, zgm, dgm,  2, wgm);
  // Xe/Isobutane-mixture
  // Get properties of Xe 
  gMC->Gfmate((*fIdmate)[4], namate, agm[0], zgm[0], d, radl, absl, buf, nbuf);
  // Get properties of Isobutane
  gMC->Gfmate((*fIdmate)[9], namate, agm[1], zgm[1], d, radl, absl, buf, nbuf);
  // Create gas mixture 
  wgm[0] = fxi;
  wgm[1] = 1. - fxi;
  dgm    = wgm[0] * dxe + wgm[1] * dis;
  AliMixture(11, "Gas mixture 2$", agm, zgm, dgm,  2, wgm);
 
  //////////////////////////////////////////////////////////////////////////
  //     Tracking Media Parameters 
  //////////////////////////////////////////////////////////////////////////

  // Al Frame 
  AliMedium(1, "Al Frame$",   1, 0, isxfld, sxmgmx
                , tmaxfd, stemax, deemax, epsil, stmin);
  // Air 
  AliMedium(2, "Air$",        2, 0, isxfld, sxmgmx
                , tmaxfd, stemax, deemax, epsil, stmin);
  // Polyethilene 
  AliMedium(3, "Radiator$",   3, 0, isxfld, sxmgmx
                , tmaxfd, stemax, deemax, epsil, stmin);
  // Xe 
  AliMedium(4, "Xe$",         4, 1, isxfld, sxmgmx
                , tmaxfd, stemax, deemax, epsil, stmin);
  // Cu pads 
  AliMedium(5, "Padplane$",   5, 1, isxfld, sxmgmx
                , tmaxfd, stemax, deemax, epsil, stmin);
  // Fee + cables 
  AliMedium(6, "Readout$",    1, 0, isxfld, sxmgmx
                , tmaxfd, stemax, deemax, epsil, stmin);
  // C frame 
  AliMedium(7, "C Frame$",    6, 0, isxfld, sxmgmx
                , tmaxfd, stemax, deemax, epsil, stmin);
  // Mylar foils 
  AliMedium(8, "Mylar$",      7, 0, isxfld, sxmgmx
                , tmaxfd, stemax, deemax, epsil, stmin);
  if (fGasMix == 1) {
    // Gas-mixture (Xe/CO2) 
    AliMedium(9, "Gas-mix$",   10, 1, isxfld, sxmgmx
                  , tmaxfd, stemax, deemax, epsil, stmin);
  }
  else {
    // Gas-mixture (Xe/Isobutane) 
    AliMedium(9, "Gas-mix$",   11, 1, isxfld, sxmgmx
                  , tmaxfd, stemax, deemax, epsil, stmin);
  }
  // Nomex-honeycomb (use carbon for the time being) 
  AliMedium(10, "Nomex$",      6, 0, isxfld, sxmgmx
                , tmaxfd, stemax, deemax, epsil, stmin);
  // Kapton foils (use Mylar for the time being) 
  AliMedium(11, "Kapton$",     7, 0, isxfld, sxmgmx
                , tmaxfd, stemax, deemax, epsil, stmin);
  // Gas-filling of the radiator 
  AliMedium(12, "CO2$",        8, 0, isxfld, sxmgmx
                , tmaxfd, stemax, deemax, epsil, stmin);
  // G10-plates
  AliMedium(13, "G10-plates$",12, 0, isxfld, sxmgmx
                , tmaxfd, stemax, deemax, epsil, stmin);
  // Cooling water
  AliMedium(14, "Water$",     13, 0, isxfld, sxmgmx
                , tmaxfd, stemax, deemax, epsil, stmin);

}

//_____________________________________________________________________________
void AliTRD::DrawModule()
{
  //
  // Draw a shaded view of the Transition Radiation Detector version 0
  //

  // Set everything unseen
  gMC->Gsatt("*"   ,"SEEN",-1);
  
  // Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN", 0);
  
  // Set the volumes visible
  if (fGeometry->IsVersion() == 0) {
    gMC->Gsatt("B071","SEEN", 0);
    gMC->Gsatt("B074","SEEN", 0);
    gMC->Gsatt("B075","SEEN", 0);
    gMC->Gsatt("B077","SEEN", 0);
    gMC->Gsatt("BTR1","SEEN", 0);
    gMC->Gsatt("BTR2","SEEN", 0);
    gMC->Gsatt("BTR3","SEEN", 0);
    gMC->Gsatt("TRD1","SEEN", 0);
    gMC->Gsatt("TRD2","SEEN", 0);
    gMC->Gsatt("TRD3","SEEN", 0);
  }
  else {
    gMC->Gsatt("B071","SEEN", 0);
    gMC->Gsatt("B074","SEEN", 0);
    gMC->Gsatt("B075","SEEN", 0);
    gMC->Gsatt("B077","SEEN", 0);
    gMC->Gsatt("BTR1","SEEN", 0);
    gMC->Gsatt("BTR2","SEEN", 0);
    gMC->Gsatt("BTR3","SEEN", 0);
    gMC->Gsatt("TRD1","SEEN", 0);
    if (fGeometry->GetPHOShole())
      gMC->Gsatt("TRD2","SEEN", 0);
    if (fGeometry->GetRICHhole())
      gMC->Gsatt("TRD3","SEEN", 0);
  }
  gMC->Gsatt("UCII","SEEN", 0);
  gMC->Gsatt("UCIM","SEEN", 0);
  gMC->Gsatt("UCIO","SEEN", 0);
  gMC->Gsatt("UL02","SEEN", 1);
  gMC->Gsatt("UL05","SEEN", 1);
  gMC->Gsatt("UL06","SEEN", 1);
  
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 2000, -2000, 2000, -2000, 2000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 40, 30, 0, 12, 9.4, .021, .021);
  gMC->Gdhead(1111, "Transition Radiation Detector");
  gMC->Gdman(18, 4, "MAN");

}

//_____________________________________________________________________________
Int_t AliTRD::DistancetoPrimitive(Int_t , Int_t )
{
  //
  // Distance between the mouse and the TRD detector on the screen
  // Dummy routine
  
  return 9999;

}
 
//_____________________________________________________________________________
void AliTRD::Init()
{
  //
  // Initialize the TRD detector after the geometry has been created
  //

  Int_t i;

  printf("\n");
  for (i = 0; i < 35; i++) printf("*");
  printf(" TRD_INIT ");
  for (i = 0; i < 35; i++) printf("*");
  printf("\n");
  printf("\n");

  if      (fGeometry->IsVersion() == 0) {
    printf("          Geometry for spaceframe with holes initialized.\n\n");
  }
  else if (fGeometry->IsVersion() == 1) {
    printf("          Geometry for spaceframe without holes initialized.\n");
    if (fGeometry->GetPHOShole())
      printf("          Leave space in front of PHOS free.\n");
    if (fGeometry->GetRICHhole())
      printf("          Leave space in front of RICH free.\n");
    printf("\n");
  }

  if (fGasMix == 1)
    printf("          Gas Mixture: 90%% Xe + 10%% CO2\n\n");
  else
    printf("          Gas Mixture: 97%% Xe + 3%% Isobutane\n\n");

}

//_____________________________________________________________________________
void AliTRD::MakeBranch(Option_t* option)
{
  //
  // Create Tree branches for the TRD digits and cluster.
  //

  Int_t  buffersize = 4000;
  Char_t branchname[15];

  AliDetector::MakeBranch(option);

  Char_t *r = strstr(option,"R");
  sprintf(branchname,"%srecPoints",GetName());
  if (fRecPoints && gAlice->TreeR() && r) {
    gAlice->TreeR()->Branch(branchname,fRecPoints->IsA()->GetName()
                           ,&fRecPoints,buffersize,0);
    printf("* AliTRD::MakeBranch * Making Branch %s for points in TreeR\n",branchname);
  }

}

//_____________________________________________________________________________
void AliTRD::ResetRecPoints()
{
  //
  // Reset number of reconstructed points and the point array
  //

  fNRecPoints = 0;
  if (fRecPoints) fRecPoints->Delete();

}

//_____________________________________________________________________________
void AliTRD::SetTreeAddress()
{
  //
  // Set the branch addresses for the trees.
  //

  Char_t branchname[15];

  AliDetector::SetTreeAddress();

  TBranch *branch;
  TTree   *treeR = gAlice->TreeR();

  if (treeR) {
    sprintf(branchname,"%srecPoints",GetName());
    if (fRecPoints) {
      branch = treeR->GetBranch(branchname);
      if (branch) {
        branch->SetAddress(&fRecPoints);
      }
    }
  }

}

//_____________________________________________________________________________
void AliTRD::SetGasMix(Int_t imix)
{
  //
  // Defines the gas mixture (imix=0:  Xe/Isobutane imix=1: Xe/CO2)
  //
  
  if ((imix < 0) || (imix > 1)) {
    printf("Wrong input value: %d\n",imix);
    printf("Use standard setting\n");
    fGasMix = 0;
    return;
  }

  fGasMix = imix;

}

//_____________________________________________________________________________
AliTRD &AliTRD::operator=(const AliTRD &trd)
{
  //
  // Assignment operator
  //

  if (this != &trd) ((AliTRD &) trd).Copy(*this);
  return *this;

}
