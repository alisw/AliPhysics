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
*/

/////////////////////////////////////////////////////////
//  Manager and hits classes for set:PHOS version 0    //
/////////////////////////////////////////////////////////
 
// --- ROOT system ---
#include "TH1.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TBRIK.h"
#include "TNode.h"

// --- galice header files ---
#include "AliPHOSv0.h"
#include "AliRun.h"
#include "AliMC.h" 

ClassImp(AliPHOSv0)

//______________________________________________________________________________


AliPHOSv0::AliPHOSv0()
{
  fIdSens=0;
}
 
//______________________________________________________________________________

AliPHOSv0::AliPHOSv0(const char *name, const char *title)
          : AliPHOS(name, title)
{
  fIdSens=0;
}
 
//___________________________________________
void AliPHOSv0::Init()
{
  fIdSens=gMC->VolId("PXTL");
}

//___________________________________________
void AliPHOSv0::CreateGeometry()
{
// *** DEFINITION OF THE -0.25<y<0.25 TILTED GEOMETRY OF THE PHOS *** 
// ORIGIN    : NICK VAN EIJNDHOVEN 

    Float_t pphi;
    Float_t r, dptcb[3], dpair[3], dphos[3], dpucp[3], dpasp[3];
    Float_t dpxtl[3];
    Float_t yo;
    Int_t idrotm[99];
    Float_t xp1, yp1, xp2, yp2;
    
    Int_t *idtmed = fIdtmed->GetArray()-699;

// --- Dimensions of PbWO4 crystal --- 
      const Float_t XTL_X=2.2;
      const Float_t XTL_Y=18.;
      const Float_t XTL_Z=2.2;
// --- Tyvek wrapper thickness 
      const Float_t PAP_THICK=0.01;
// --- Polystyrene Foam Outer Cover dimensions --- 
      const Float_t FOC_X=214.6;
      const Float_t FOC_Y=80.;
      const Float_t FOC_Z=260.;
// --- Inner AIR volume dimensions --- 
      const Float_t AIR_X=206.;
      const Float_t AIR_Y=66.;
      const Float_t AIR_Z=244.;
// --- Tyvek Crystal Block dimensions --- 
      const Float_t TCB_X=198.;
      const Float_t TCB_Y=25.;
      const Float_t TCB_Z=234.;
// --- Upper Cooling Plate thickness --- 
      const Float_t UCP_Y=0.06;
// --- Al Support Plate thickness --- 
      const Float_t ASP_Y=10.;
//--- Distance from IP to Foam Outer Cover top plate (needs to be 447.) ---
      const Float_t FOC_R=467.;
//--- Distance from IP to Crystal Block top Surface (needs to be 460.) ---
      const Float_t CBS_R=480.;

// --- Dimensions of volumes --- 


// --- Define PHOS box volume, fill with Polystyrene foam --- 
    dphos[0] = FOC_X/2.;
    dphos[1] = FOC_Y/2.;
    dphos[2] = FOC_Z/2.;
    gMC->Gsvolu("PHOS", "BOX ", idtmed[703], dphos, 3);

// --- Define air-filled box, place inside PHOS --- 
    dpair[0] = AIR_X/2.;
    dpair[1] = AIR_Y/2.;
    dpair[2] = AIR_Z/2.;
    gMC->Gsvolu("PAIR", "BOX ", idtmed[798], dpair, 3);
    gMC->Gspos("PAIR", 1, "PHOS", 0., 0., 0., 0, "ONLY");

// --- Define Upper Cooling Panel --- 
// --- place it right behind upper foam plate --- 
    dpucp[0] = TCB_X/2.;
    dpucp[1] = UCP_Y/2.;
    dpucp[2] = TCB_Z/2.;
    gMC->Gsvolu("PUCP", "BOX ", idtmed[701], dpucp, 3);
    yo = (AIR_Y-UCP_Y)/2.;
    gMC->Gspos("PUCP", 1, "PAIR", 0., yo, 0., 0, "ONLY");

// --- Define Crystal Block, fill with Tyvek, position inside PAIR --- 
    dptcb[0] = TCB_X/2.;
    dptcb[1] = TCB_Y/2.;
    dptcb[2] = TCB_Z/2.;
    gMC->Gsvolu("PTCB", "BOX ", idtmed[702], dptcb, 3);
// --- Divide PTCB in X and Z directions -- 
    gMC->Gsdvn("PSEC", "PTCB", 11, 1);
    gMC->Gsdvn("PMOD", "PSEC", 13, 3);
    gMC->Gsdvn("PSTR", "PMOD", 8, 1);
    gMC->Gsdvn("PCEL", "PSTR", 8, 3);
    yo = (FOC_Y-TCB_Y)/2. -(CBS_R-FOC_R);
    gMC->Gspos("PTCB", 1, "PAIR", 0., yo, 0., 0, "ONLY");

// --- Define PbWO4 crystal volume, place inside PCEL --- 
    dpxtl[0] = XTL_X/2.;
    dpxtl[1] = XTL_Y/2.;
    dpxtl[2] = XTL_Z/2.;
    gMC->Gsvolu("PXTL", "BOX ", idtmed[699], dpxtl, 3);
    yo = (TCB_Y-XTL_Y)/2. - PAP_THICK;
    gMC->Gspos("PXTL", 1, "PCEL", 0., yo, 0., 0, "ONLY");

// --- Define Al Support Plate, position it inside PAIR --- 
// --- right beneath PTCB --- 
    dpasp[0] = AIR_X/2.;
    dpasp[1] = ASP_Y/2.;
    dpasp[2] = AIR_Z/2.;
    gMC->Gsvolu("PASP", "BOX ", idtmed[701], dpasp, 3);
    yo = (FOC_Y-ASP_Y)/2. - (CBS_R-FOC_R+TCB_Y);
    gMC->Gspos("PASP", 1, "PAIR", 0., yo, 0., 0, "ONLY");

// --- Divide in X and Z direction (same way as PTCB) --- 
    gMC->Gsdvn("PCMO", "PCSE", 13, 3);
    gMC->Gsdvn("PCST", "PCMO", 8, 1);
    gMC->Gsdvn("PCCE", "PCST", 8, 3);

// --- Position various PHOS units in ALICE setup --- 
// --- PHOS itself first --- 
    r     = FOC_R+FOC_Y/2.;
    pphi  = TMath::ATan(FOC_X/(2.*FOC_R));
    xp1   = -r * TMath::Sin(pphi * 3.);
    yp1   = -r * TMath::Cos(pphi * 3.);
    xp2   = -r * TMath::Sin(pphi);
    yp2   = -r * TMath::Cos(pphi);
    pphi *= 180/kPI;
    AliMatrix(idrotm[0], 90.,-3*pphi, 90., 90-3*pphi, 0., 0.);
    AliMatrix(idrotm[1], 90.,  -pphi, 90., 90-pphi,   0., 0.);
    AliMatrix(idrotm[2], 90.,   pphi, 90., 90+pphi,   0., 0.);
    AliMatrix(idrotm[3], 90., 3*pphi, 90., 90+3*pphi, 0., 0.);
    gMC->Gspos("PHOS", 1, "ALIC", xp1, yp1, 0., idrotm[0], "ONLY");
    gMC->Gspos("PHOS", 2, "ALIC", xp2, yp2, 0., idrotm[1], "ONLY");
    gMC->Gspos("PHOS", 3, "ALIC",-xp2, yp2, 0., idrotm[2], "ONLY");
    gMC->Gspos("PHOS", 4, "ALIC",-xp1, yp1, 0., idrotm[3], "ONLY");

// --- Set modules seen without tree for drawings --- 
    gMC->Gsatt("PMOD", "SEEN", -2);
    gMC->Gsatt("PCMO", "SEEN", -2);
}
 
//___________________________________________
void AliPHOSv0::CreateMaterials()
{
// *** DEFINITION OF AVAILABLE PHOS MATERIALS *** 
// ORIGIN    : NICK VAN EIJNDHOVEN 

    Int_t   ISXFLD = gAlice->Field()->Integ();
    Float_t SXMGMX = gAlice->Field()->Max();
    
// --- The PbWO4 crystals --- 
    Float_t ax[3] = { 207.19,183.85,16. };
    Float_t zx[3] = { 82.,74.,8. };
    Float_t wx[3] = { 1.,1.,4. };
    Float_t dx    = 8.28;
// --- The polysterene scintillator (CH) --- 
    Float_t ap[2] = { 12.011,1.00794 };
    Float_t zp[2] = { 6.,1. };
    Float_t wp[2] = { 1.,1. };
    Float_t dp    = 1.032;
// --- Tyvek (CnH2n) 
    Float_t at[2] = { 12.011,1.00794 };
    Float_t zt[2] = { 6.,1. };
    Float_t wt[2] = { 1.,2. };
    Float_t dt    = .331;
// --- Polystyrene foam --- 
    Float_t af[2] = { 12.011,1.00794 };
    Float_t zf[2] = { 6.,1. };
    Float_t wf[2] = { 1.,1. };
    Float_t df    = .3;

    Int_t *idtmed = fIdtmed->GetArray()-699;
    
    AliMixture( 0, "PbWO4$",       ax, zx, dx, -3, wx);
    AliMixture( 1, "Polystyrene$", ap, zp, dp, -2, wp);
    AliMaterial(2, "Al$",          26.98, 13., 2.7, 8.9, 999.);
// ---                                Absorption length^ is ignored --- 
    AliMixture( 3, "Tyvek$", at, zt, dt, -2, wt);
    AliMixture( 4, "Foam$",  af, zf, df, -2, wf);
    AliMaterial(9, "Air$", 14.61, 7.3, .001205, 30420., 67500);

    AliMedium(0, "PHOS Xtal    $", 0, 1, ISXFLD, SXMGMX, 10., .1, .1, .1, .1);
    AliMedium(2, "Al parts     $", 2, 0, ISXFLD, SXMGMX, 10., .1, .1, .001, .001);
    AliMedium(3, "Tyvek wrapper$", 3, 0, ISXFLD, SXMGMX, 10., .1, .1, .001, .001);
    AliMedium(4, "Polyst. foam $", 4, 0, ISXFLD, SXMGMX, 10., .1, .1, .1, .1);
    AliMedium(99, "Air          $", 9, 0, ISXFLD, SXMGMX, 10., 1., .1, .1, 10.);

// --- Generate explicitly delta rays in aluminium parts --- 
    gMC->Gstpar(idtmed[701], "LOSS", 3.);
    gMC->Gstpar(idtmed[701], "DRAY", 1.);
}

void AliPHOSv0::StepManager()
{

  TClonesArray &lhits = *fHits;
  TLorentzVector p;
  Int_t copy, i;
  Int_t vol[5];
  Float_t hits[4];
  if(gMC->CurrentVolID(copy) == fIdSens) {
    //
    //We are in the sensitive volume
    for(i=0;i<4;i++) {
      gMC->CurrentVolOffID(i+1,copy);
      vol[4-i]=copy;
    }
    gMC->CurrentVolOffID(7,copy);
    vol[0]=copy;
    gMC->TrackPosition(p);
    for(i=0;i<3;++i) hits[i]=p[i];
    hits[3]=gMC->Edep();
    new(lhits[fNhits++]) AliPHOShit(fIshunt,gAlice->CurrentTrack(),vol,hits);
  }
}
