////////////////////////////////////////////////
//  Manager and hits classes for set:CPV      //
//                                            //
//  Author: Yuri Kharlov, IHEP, Protvino      //
//  e-mail: Yuri.Kharlov@cern.ch              //
//  Last modified: 18 September 1999          //
////////////////////////////////////////////////
 
// --- ROOT system ---
#include "TH1.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TBRIK.h"
#include "TNode.h"
#include "TMath.h"

// --- Standard library ---
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// --- galice header files ---
#include "AliCPV.h"
#include "AliRun.h"

extern "C" void cpvrec_ (Float_t *x1, Float_t *x2, Float_t *x3, Int_t *x4, 
                         Float_t *x5, Float_t *x6, Int_t *x7, Float_t *x8);

//==============================================================================
//                              AliCPVExactHit
//==============================================================================

ClassImp(AliCPVExactHit)

//______________________________________________________________________________

AliCPVExactHit::AliCPVExactHit(TLorentzVector p, Float_t *xy, Int_t ipart)
{
//
// Creates an exact CPV hit object
//

  Int_t i;
  for (i=0; i<2; i++) fXYhit[i]  = xy[i];
  fMomentum  = p;
  fIpart     = ipart;
}

//______________________________________________________________________________

void AliCPVExactHit::Print()
{
// Print exact hit

  printf("CPV hit: p  = (%f, %f, %f, %f) GeV,\n",
	 fMomentum.Px(),fMomentum.Py(),fMomentum.Pz(),fMomentum.E());
  printf("         xy = (%f, %f) cm,   ipart = %d\n",
	 fXYhit[0],fXYhit[1],fIpart);
}
 
//==============================================================================
//                              AliCPVHit
//==============================================================================

ClassImp(AliCPVHit)

//______________________________________________________________________________

AliCPVHit::AliCPVHit(Float_t *xy)
{
//
// Creates a CPV hit object
//

  Int_t i;
  for (i=0; i<2; i++) fXYhit[i]  = xy[i];
}

//______________________________________________________________________________

void AliCPVHit::Print()
{
// Print reconstructed hit

  printf("CPV hit: xy = (%f, %f) cm\n",
	 fXYhit[0],fXYhit[1]);
}
 
//==============================================================================
//                              AliCPVCradle
//==============================================================================

ClassImp(AliCPVCradle)

//______________________________________________________________________________

AliCPVCradle::AliCPVCradle(void) {
// Default constructor

  fNhitsExact         = 0;
  fNhitsReconstructed = 0;

  // Allocate the array of exact and reconstructed hits
  fHitExact          = new TClonesArray("AliCPVExactHit", 100);
  fHitReconstructed  = new TClonesArray("AliCPVHit",      100);
}

//______________________________________________________________________________

AliCPVCradle::AliCPVCradle(Int_t   Geometry  ,
                           Float_t PadZSize  ,
                           Float_t PadPhiSize,
                           Float_t Radius    ,
                           Float_t Thickness ,
                           Float_t Angle     ,
                           Int_t   Nz        ,
                           Int_t   Nphi      )
{
// Set geometry parameters and allocate space for the hit storage

  fPadZSize   = PadZSize;
  fPadPhiSize = PadPhiSize;
  fRadius     = Radius;
  fThickness  = Thickness;
  fAngle      = Angle;
  fNz         = Nz;
  fNphi       = Nphi;
  fNhitsExact         = 0;
  fNhitsReconstructed = 0;

  // Allocate the array of exact and reconstructed hits
  fHitExact          = new TClonesArray("AliCPVExactHit", 100);
  fHitReconstructed  = new TClonesArray("AliCPVHit",      100);
}

//______________________________________________________________________________

AliCPVCradle::~AliCPVCradle(void)
{
// Default destructor

  Clear();
}

//______________________________________________________________________________

void AliCPVCradle::Clear(Option_t *opt="")
{
// Clear hit information

  fHitExact         -> Clear(opt);
  fHitReconstructed -> Clear(opt);
  fNhitsExact         = 0;
  fNhitsReconstructed = 0;
}

//______________________________________________________________________________

void AliCPVCradle::AddHit(TLorentzVector p, Float_t *xy, Int_t ipart)
{
// Add this hit to the hits list in CPV detector.

  TClonesArray &lhits = *fHitExact;
  new(lhits[fNhitsExact++]) AliCPVExactHit(p,xy,ipart);
}

//______________________________________________________________________________

void AliCPVCradle::Print(Option_t *opt)
{
// Print AliCPVCradle information.
// 
// options:  'g' - print cradle geometry
//           'p' - print generated hits in the cradle
//           'r' - print reconstructed hits in the cradle

  if (strcmp(opt,"g")==0) {
    printf ("CPVCradle: size = %f x %f cm\n",fPadZSize*fNz,
	                                     fPadPhiSize*fNphi);
  }
  else if (strcmp(opt,"p")==0) {
    printf ("CPV cradle has %d generated hits\n",fNhitsExact);
    TIter next(fHitExact);
    AliCPVExactHit *hit0;
    while ( (hit0 = (AliCPVExactHit*)next()) ) hit0 -> Print();
  }
  else if (strcmp(opt,"r")==0) {
    printf ("CPV cradle has %d reconstructed hits\n",fNhitsReconstructed);
    TIter next(fHitReconstructed);
    AliCPVHit *hit0;
    while ( (hit0 = (AliCPVHit*)next()) ) hit0 -> Print();
  }
}

//______________________________________________________________________________

void AliCPVCradle::Reconstruction(Float_t min_distance, Float_t min_signal)
{
// This function:
// 1) Takes on input the array of generated (exact) hits in the CPV cradle,
// 2) Founds the pad response according the known pad-response function,
// 3) Solves the inverse problem to find the array of reconstructed hits

  Float_t DzCPV=fPadZSize*fNz/2,
          DyCPV=fPadPhiSize*fNphi/2,
          DxCPV=fThickness/2;
  Float_t Xin[5000][2],
          Pin[5000][4],
          Xout[5000][2];
  fHitReconstructed->Clear();   
  fNhitsReconstructed=0;
  if (fHitExact) {
    TIter next(fHitExact);
    AliCPVExactHit *hit;
    TClonesArray &lhits = *fHitReconstructed;
    for(int i=0;i<fNhitsExact;i++) {
      hit = (AliCPVExactHit*)next();
      Xin[i][0]=hit->fXYhit[0];
      Xin[i][1]=hit->fXYhit[1];
      Pin[i][0]=hit->fMomentum.Pz();
      Pin[i][1]=hit->fMomentum.Py();
      Pin[i][2]=hit->fMomentum.Px();
      Pin[i][3]=hit->fMomentum.E();
    }
    cpvrec_(&DzCPV,&DyCPV,&DxCPV,&fNhitsExact,&Xin[0][0],&Pin[0][0],
	    &fNhitsReconstructed,&Xout[0][0]);
    for(int i=0;i<fNhitsReconstructed;i++) new(lhits[i]) AliCPVHit(Xout[i]);
  }
}

//==============================================================================
//                              AliCPV
//==============================================================================

ClassImp(AliCPV)

//______________________________________________________________________________

AliCPV::~AliCPV(void)
{
  fCradles->Delete();
  delete fCradles;
}

//______________________________________________________________________________

AliCPV::AliCPV() :
  fDebugLevel            (0)
{
  fNCradles   = 4;           // Number of cradles

  if( NULL==(fCradles=new TClonesArray("AliCPVCradle",fNCradles)) )
  {
    Error("AliCPV","Can not create fCradles");
    exit(1);
  }
}
 
//______________________________________________________________________________

AliCPV::AliCPV(const char *name, const char *title)
       : AliDetector (name,title),
         fDebugLevel (0)
{
//Begin_Html
/*
<img src="picts/aliCPV.gif">
*/
//End_Html
 
  SetMarkerColor(kGreen);
  SetMarkerStyle(2);
  SetMarkerSize(0.4);
  
  fNCradles   = 4;           // Number of cradles
  fPadZSize   = 2.26;        // Pad size along beam       [cm]
  fPadPhiSize = 1.13;        // Pad size across beam      [cm]
  fRadius     = 455.;        // Distance of CPV from IP   [cm]
  fThickness  = 1.;          // CPV thickness             [cm]
  fAngle      = 27.0;        // Angle between CPV cradles [deg]
  fNz         = 52;          // Number of pads along beam
  fNphi       = 176;         // Number of pads across beam

  if( NULL==(fCradles=new TClonesArray("AliCPVCradle",fNCradles)) )
  {
    Error("AliCPV","Can not create fCradles");
    exit(1);
  }
}

//______________________________________________________________________________

void AliCPV::SetGeometry (Int_t ncradles, Int_t nz, Int_t nphi, Float_t angle)
{
// Set CPV geometry which differs from the default one

  fNCradles   = ncradles;     // Number of cradles
  fNz         = nz;           // Number of pads along beam
  fNphi       = nphi;         // Number of pads across beam
  fAngle      = angle;        // Angle between CPV cradles [deg]
}

//______________________________________________________________________________

void AliCPV::Init()
{
  Int_t i;
  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" CPV_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  for(i=0;i<80;i++) printf("*");
  printf("\n");
}

//______________________________________________________________________________

void AliCPV::BuildGeometry()
{
//
// Build simple ROOT TNode geometry for event display
//

  TNode *Node, *Top;

  const int kColorCPV = kGreen;
  //
  Top=gAlice->GetGeometry()->GetNode("alice");


  // CPV
  Float_t pphi=fAngle;
  new TRotMatrix("rotCPV1","rotCPV1",90,-3*pphi,90,90-3*pphi,0,0);
  new TRotMatrix("rotCPV2","rotCPV2",90,-  pphi,90,90-  pphi,0,0);
  new TRotMatrix("rotCPV3","rotCPV3",90,   pphi,90,90+  pphi,0,0);
  new TRotMatrix("rotCPV4","rotCPV4",90, 3*pphi,90,90+3*pphi,0,0);
  new TBRIK("S_CPV","CPV box","void",99.44,0.50,58.76);
  Top->cd();
  Node = new TNode("CPV1","CPV1","S_CPV",-317.824921,-395.014343,0,"rot988");
  Node->SetLineColor(kColorCPV);
  fNodes->Add(Node);
  Top->cd();
  Node = new TNode("CPV2","CPV2","S_CPV",-113.532333,-494.124908,0,"rot989");
  fNodes->Add(Node);
  Node->SetLineColor(kColorCPV);
  Top->cd();
  Node = new TNode("CPV3","CPV3","S_CPV", 113.532333,-494.124908,0,"rot990");
  Node->SetLineColor(kColorCPV);
  fNodes->Add(Node);
  Top->cd();
  Node = new TNode("CPV4","CPV4","S_CPV", 317.824921,-395.014343,0,"rot991");
  Node->SetLineColor(kColorCPV);
  fNodes->Add(Node);
}
 
//______________________________________________________________________________
void AliCPV::CreateMaterials()
{
// *** DEFINITION OF AVAILABLE CPV MATERIALS *** 

// CALLED BY : CPV_MEDIA 
// ORIGIN    : NICK VAN EIJNDHOVEN 

  Int_t *idtmed = fIdtmed->GetArray()-1999;

  Int_t   ISXFLD = gAlice->Field()->Integ();
  Float_t SXMGMX = gAlice->Field()->Max();
    
// --- The polysterene scintillator (CH) --- 
  Float_t ap[2] = { 12.011,1.00794 };
  Float_t zp[2] = { 6.,1. };
  Float_t wp[2] = { 1.,1. };
  Float_t dp    = 1.032;
  
  AliMixture ( 1, "Polystyrene$",    ap, zp, dp, -2, wp);
  AliMedium  ( 1, "CPV scint. $", 1, 1, ISXFLD, SXMGMX, 10., .1, .1, .1, .1);

  // --- Generate explicitly delta rays in the CPV media --- 
  gMC->Gstpar(idtmed[2000], "LOSS", 3.);
  gMC->Gstpar(idtmed[2000], "DRAY", 1.);
}

//______________________________________________________________________________

void AliCPV::AddCPVCradles()
{
// Create array of CPV cradles

  TClonesArray &lcradle = *fCradles;
  for(Int_t i=0; i<fNCradles; i++) {
    new(lcradle[i]) AliCPVCradle( IsVersion(),
                                  fPadZSize  ,
                                  fPadPhiSize,
                                  fRadius    ,
                                  fThickness ,
                                  fAngle     ,
                                  fNz        ,
                                  fNphi     );
  }
}

//______________________________________________________________________________

void AliCPV::Reconstruction(Float_t min_distance, Float_t min_signal)
{
// Performs reconstruction for all CPV cradles

  for( Int_t i=0; i<fNCradles; i++ )
    GetCradle(i).Reconstruction(min_distance, min_signal);
}

//______________________________________________________________________________

void AliCPV::ResetDigits(void)
{
  AliDetector::ResetDigits();

  for( int i=0; i<fNCradles; i++ )
    ((AliCPVCradle*)(*fCradles)[i]) -> Clear();
}

//______________________________________________________________________________

void AliCPV::FinishEvent(void)
{
// Called at the end of each 'galice' event.

}

//______________________________________________________________________________

void AliCPV::FinishRun(void)
{
}

//______________________________________________________________________________

void AliCPV::Print(Option_t *opt)
{
// Print CPV information.
// For each AliCPVCradle the function AliCPVCradle::Print(opt) is called.

  AliCPV &CPV = *(AliCPV *)this;     // Removing 'const'...

  if (strcmp(opt,"g")==0) {
    for( Int_t i=0; i<fNCradles; i++ ) {
      printf("CPV cradle %d of %d\n",i+1, fNCradles);
      CPV.GetCradle(i).Print(opt);
      printf( "-------------------------------------------------------------\n");
    }
  }
  else if (strcmp(opt,"p")==0) {
    for( Int_t i=0; i<fNCradles; i++ ) {
      if ( (CPV.GetCradle(i).GetHitExact())->GetEntries()!=0 ) {
	printf("CPV cradle %d of %d\n",i+1, fNCradles);
	CPV.GetCradle(i).Print(opt);
	printf( "-------------------------------------------------------------\n");
      }
    }
  }
  else if (strcmp(opt,"r")==0) {
    for( Int_t i=0; i<fNCradles; i++ ) {
      if ( (CPV.GetCradle(i).GetHitReconstructed())->GetEntries()!=0 ) {
	printf("CPV cradle %d of %d\n",i+1, fNCradles);
	CPV.GetCradle(i).Print(opt);
	printf( "-------------------------------------------------------------\n");
      }
    }
  }
}
