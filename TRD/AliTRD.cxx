///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Transition Radiation Detector                                            //
//  This class contains the basic functions for the Transition Radiation     //
//  detector. Functions specific to one particular geometry are              //
//  contained in the derived classes                                         //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliTRDClass.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TRandom.h>
#include <TVector.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TObjectTable.h>
#include <TBRIK.h>
#include <TPGON.h> 

#include "GParticle.h"
#include "AliTRD.h"
#include "AliRun.h"
#include <iostream.h>
#include <fstream.h>

#include "TGeant3.h"
 
ClassImp(AliTRD)
 
//_____________________________________________________________________________
AliTRD::AliTRD()
{
  //
  // Default constructor
  //
  fIshunt   = 0;
}
 
//_____________________________________________________________________________
AliTRD::AliTRD(const char *name, const char *title)
       : AliDetector(name,title)
{
  //
  // Standard constructor for the TRD
  //

  //
  // Allocate the hit array
 
  fHits   = new TClonesArray("AliTRDhit",  405);
  
  fIshunt = 0;
  
  SetMarkerColor(kWhite);   
}
 
//_____________________________________________________________________________
void AliTRD::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a hit for the TRD
  //
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliTRDhit(fIshunt,track,vol,hits);
}

//_____________________________________________________________________________
void AliTRD::BuildGeometry()
{
  //
  // Create the ROOT TNode geometry for the TRD
  //
  TNode *Node, *Top;
  TPGON *pgon;
  const Int_t kColorTRD = 46;
  //
  // Find the top node alice
  Top=gAlice->GetGeometry()->GetNode("alice");
  //
  pgon = new TPGON("S_TRD","TRD","void",0,360,18,2);
  pgon->DefineSection(0,-351.2,285.335,355.686);
  pgon->DefineSection(1, 351.2,285.335,355.686);
  Top->cd();
  Node = new TNode("TRD","TRD","S_TRD",0,0,0,"");
  Node->SetLineColor(kColorTRD);
  fNodes->Add(Node);
}
 
//_____________________________________________________________________________
void AliTRD::CreateMaterials()
{
  //
  // Create the materials for the TRD
  // Origin Y.Foka
  //

  AliMC* pMC = AliMC::GetMC();
  
  Int_t   ISXFLD = gAlice->Field()->Integ();
  Float_t SXMGMX = gAlice->Field()->Max();
  
  // --- For polyethilene (CH2) 
  Float_t ape[2] = { 12.,1. };
  Float_t zpe[2] = { 6.,1. };
  Float_t wpe[2] = { .33,.67 };
  Float_t dpe    = .95;
  // --- For mylar (C5H4O2) 
  Float_t amy[3] = { 12.,1.,16. };
  Float_t zmy[3] = { 6.,1.,8. };
  Float_t wmy[3] = { 5.,4.,2. };
  Float_t dmy    = 1.39;
  // --- For CO2 
  Float_t aco[2] = { 12.,16. };
  Float_t zco[2] = { 6.,8. };
  Float_t wco[2] = { 1.,2. };
  Float_t dco    = .001977;
  // --- For Xe/CO2-gas-mixture 
  // --- Xe-content of the gas-mixture (90% / 10%) 
  Float_t fxe    = .9;
  Float_t dxe    = .005858;
  
  // --- General tracking parameter
  Float_t tmaxfd = -10.;
  Float_t stemax = -1e10;
  Float_t deemax = -.1;
  Float_t epsil  = 1e-4;
  Float_t stmin  = -.001;
  
  Float_t absl, radl, d, buf[1];
  Float_t agm[2], dgm, zgm[2], wgm[2];
  Int_t nbuf;
  
  //     Define Materials 
  
  // ***********************************************************************
  AliMaterial(1, "Al $", 26.98, 13., 2.7, 8.9, 37.2);
  AliMaterial(2, "Air$", 14.61, 7.3, .001205, 30420., 67500);
  AliMaterial(4, "Xe $", 131.29, 54., dxe, 1447.59, 0);
  AliMaterial(5, "Cu $", 63.54, 29., 8.96, 1.43, 14.8);
  AliMaterial(6, "C  $", 12.01, 6., 2.265, 18.8, 74.4);
  // --- Mixtures 
  AliMixture(3, "Polyethilene$", ape, zpe, dpe, -2, wpe);
  AliMixture(7, "Mylar$",        amy, zmy, dmy, -3, wmy);
  AliMixture(8, "CO2$",          aco, zco, dco, -2, wco);
  char namate[21];
  // --- Get properties of Xe 
  pMC->Gfmate((*fIdmate)[4], namate, agm[0], zgm[0], d, radl, absl, buf, nbuf);
  // --- Get properties of CO2 
  pMC->Gfmate((*fIdmate)[8], namate, agm[1], zgm[1], d, radl, absl, buf, nbuf);
  // --- Create gas mixture 
  wgm[0] = fxe;
  wgm[1] = 1. - fxe;
  dgm    = wgm[0] * dxe + wgm[1] * dco;
  AliMixture(9, "Gas mixture$", agm, zgm, dgm, 2, wgm);
  // ***********************************************************************
  
  //     Tracking Media Parameters 
  
  // ***********************************************************************
  // --- Al Frame 
  AliMedium(1301, "Al Frame$", 1, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  // --- Air 
  AliMedium(1302, "Air$",      2, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  // --- Polyethilene 
  AliMedium(1303, "Radiator$", 3, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  // --- Xe 
  AliMedium(1304, "Xe$",       4, 1, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  // --- Cu pads 
  AliMedium(1305, "Padplane$", 5, 1, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  // --- Fee + cables 
  AliMedium(1306, "Readout$",  1, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  // --- C frame 
  AliMedium(1307, "C Frame$",  6, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  // --- Mylar foils 
  AliMedium(1308, "Mylar$",    7, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  // --- Gas-mixture 
  AliMedium(1309, "Gas-mix$",  9, 1, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  // --- Nomex-honeycomb (use carbon for the time being) 
  AliMedium(1310, "Nomex$",    6, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  // --- Kapton foils (use Mylar for the time being) 
  AliMedium(1311, "Kapton$",   7, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  // --- Gas-filling of the radiator 
  AliMedium(1312, "CO2$",      8, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
}

//_____________________________________________________________________________
Int_t AliTRD::DistancetoPrimitive(Int_t , Int_t )
{
  //
  // Distance between the mouse and the TRD detector on the screen
  // Dummy routine
  //
   return 9999;
}
 
//_____________________________________________________________________________
void AliTRD::Init()
{
  //
  // Initialise the TRD detector after the geometry has been created
  //
  Int_t i;
  //
  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" TRD_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  //
  // Here the TRD initialisation code (if any!)
  for(i=0;i<80;i++) printf("*");
  printf("\n");
}

ClassImp(AliTRDhit)
 
//_____________________________________________________________________________
AliTRDhit::AliTRDhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
  AliHit(shunt, track)
{
  //
  // Create a TRD hit
  //

  //
  // Store volume hierarchy
  fSector  = vol[0]; 
  fChamber = vol[1];
  fPlane   = vol[2];
  //
  // Store position and charge
  fX       = hits[0];
  fY       = hits[1];
  fZ       = hits[2];
  fQ       = hits[3];
}
