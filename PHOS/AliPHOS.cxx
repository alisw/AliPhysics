////////////////////////////////////////////////
//  Manager and hits classes for set:PHOS     //
////////////////////////////////////////////////
 
// --- ROOT system ---
#include "TH1.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TBRIK.h"
#include "TNode.h"

// --- Standard library ---
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// --- galice header files ---
#include "AliPHOS.h"
#include "AliRun.h"
#include "AliMC.h" 
#include "TGeant3.h"

//______________________________________________________________________________


ClassImp(AliPHOS)

//______________________________________________________________________________

AliPHOS::~AliPHOS(void)
{
  fCradles->Delete();
  delete fCradles;
}

//______________________________________________________________________________

AliPHOS::AliPHOS() :
         fDebugLevel            (0),
         fTreePHOS              (NULL),
         fBranchNameOfCradles   ("AliPHOSCradles"),
         fTreeName              ("PHOS")
{
   fIshunt   = 0;

  if( NULL==(fCradles=new TObjArray) )
  {
    Error("AliPHOS","Can not create fCradles");
    exit(1);
  }
  DefPars();
}
 
//______________________________________________________________________________

AliPHOS::AliPHOS(const char *name, const char *title)
       : AliDetector            (name,title),
         fDebugLevel            (0),
         fTreePHOS              (NULL),
         fBranchNameOfCradles   ("AliPHOSCradles"),
         fTreeName              ("PHOS")
{
//Begin_Html
/*
<img src="gif/aliphos.gif">
*/
//End_Html
 
   fHits   = new TClonesArray("AliPHOShit",  405);
 
   fIshunt     =  0;

   SetMarkerColor(kGreen);
   SetMarkerStyle(2);
   SetMarkerSize(0.4);

  if( NULL==(fCradles=new TObjArray) ) {
     Error("AliPHOS","Can not create fCradles");
     exit(1);
  }
  DefPars();
}

//______________________________________________________________________________

void AliPHOS::DefPars()
{ 
      PHOSflags[0]=0;
      PHOSflags[1]=1;
      PHOSflags[2]=0;
      PHOSflags[3]=0;
      PHOSflags[4]=0;
      PHOSflags[5]=0;
      PHOSflags[6]=0;
      PHOSflags[7]=0;
      PHOSflags[8]=0;
      PHOScell[0]=2.2;
      PHOScell[1]=18.;
      PHOScell[2]=0.01;
      PHOScell[3]=0.01;
      PHOScell[4]=1.0;
      PHOScell[5]=0.1;
      PHOScell[6]=0.;
      PHOScell[7]=0.;
      PHOScell[8]=0.;
      PHOSradius=460.;
      PHOSsize[0]=104;
      PHOSsize[1]=88;
      PHOSsize[2]=4;
      PHOScradlesA=0.;
      PHOSCPV[0]=1.;
      PHOSCPV[1]=2.;
      PHOSCPV[2]=0.;
      PHOSCPV[3]=0.;
      PHOSCPV[4]=0.;
      PHOSCPV[5]=0.;
      PHOSCPV[6]=0.;
      PHOSCPV[7]=0.;
      PHOSCPV[8]=0.;
      PHOSextra[0]=0.001;
      PHOSextra[1]=6.95;
      PHOSextra[2]=4.;
      PHOSextra[3]=5.;
      PHOSextra[4]=2.;
      PHOSextra[5]=0.06;
      PHOSextra[6]=10.;
      PHOSextra[7]=3.;
      PHOSextra[8]=1.;
      PHOSTXW[0]=209.;
      PHOSTXW[1]=71.;
      PHOSTXW[2]=250.;
      PHOSAIR[0]=206.;
      PHOSAIR[1]=66.;
      PHOSAIR[2]=244.;
      PHOSFTI[0]=214.6;
      PHOSFTI[1]=80.;
      PHOSFTI[2]=260.;
      PHOSFTI[3]=467.;
}
//______________________________________________________________________________

void AliPHOS::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliPHOShit(fIshunt,track,vol,hits);
}
 
//___________________________________________
void AliPHOS::BuildGeometry()
{

  TNode *Node, *Top;

  const int kColorPHOS = kRed;
  //
  Top=gAlice->GetGeometry()->GetNode("alice");


  // PHOS
  Float_t pphi=12.9399462;
  new TRotMatrix("rot988","rot988",90,-3*pphi,90,90-3*pphi,0,0);
  new TRotMatrix("rot989","rot989",90,-  pphi,90,90-  pphi,0,0);
  new TRotMatrix("rot990","rot990",90,   pphi,90,90+  pphi,0,0);
  new TRotMatrix("rot991","rot991",90, 3*pphi,90,90+3*pphi,0,0);
  new TBRIK("S_PHOS","PHOS box","void",107.3,40,130);
  Top->cd();
  Node = new TNode("PHOS1","PHOS1","S_PHOS",-317.824921,-395.014343,0,"rot988");
  Node->SetLineColor(kColorPHOS);
  fNodes->Add(Node);
  Top->cd();
  Node = new TNode("PHOS2","PHOS2","S_PHOS",-113.532333,-494.124908,0,"rot989");
  fNodes->Add(Node);
  Node->SetLineColor(kColorPHOS);
  Top->cd();
  Node = new TNode("PHOS3","PHOS3","S_PHOS", 113.532333,-494.124908,0,"rot990");
  Node->SetLineColor(kColorPHOS);
  fNodes->Add(Node);
  Top->cd();
  Node = new TNode("PHOS4","PHOS4","S_PHOS", 317.824921,-395.014343,0,"rot991");
  Node->SetLineColor(kColorPHOS);
  fNodes->Add(Node);
}
 
//___________________________________________
void AliPHOS::CreateMaterials()
{
// *** DEFINITION OF AVAILABLE PHOS MATERIALS *** 

// CALLED BY : PHOS_MEDIA 
// ORIGIN    : NICK VAN EIJNDHOVEN 


  AliMC* pMC = AliMC::GetMC();

    Int_t   ISXFLD = gAlice->Field()->Integ();
    Float_t SXMGMX = gAlice->Field()->Max();
    
// --- The PbWO4 crystals --- 
    Float_t ax[3] = { 207.19,183.85,16. };
    Float_t zx[3] = { 82.,74.,8. };
    Float_t wx[3] = { 1.,1.,4. };
    Float_t dx    = 8.28;
// --- Stainless Steel --- 
    Float_t as[5] = { 55.847,12.011,51.9961,58.69,28.0855 };
    Float_t zs[5] = { 26.,6.,24.,28.,14. };
    Float_t ws[5] = { .6392,8e-4,.2,.14,.02 };
    Float_t ds    = 8.;
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
    Float_t df    = .12;
//--- Foam thermo insulation (actual chemical composition unknown yet!) ---
    Float_t ati[2] = { 12.011,1.00794 };
    Float_t zti[2] = { 6.,1. };
    Float_t wti[2] = { 1.,1. };
    Float_t dti    = .1;
// --- Textolit (actual chemical composition unknown yet!) --- 
    Float_t atx[2] = { 12.011,1.00794 };
    Float_t ztx[2] = { 6.,1. };
    Float_t wtx[2] = { 1.,1. };
    Float_t dtx    = 1.83;

    Int_t *idtmed = gAlice->Idtmed();
    

    AliMixture(  0, "PbWO4$",          ax, zx, dx, -3, wx);
    AliMixture(  1, "Polystyrene$",    ap, zp, dp, -2, wp);
    AliMaterial( 2, "Al$",             26.98, 13., 2.7, 8.9, 999);
// ---                                Absorption length^ is ignored --- 
    AliMixture(  3, "Tyvek$",           at, zt, dt, -2, wt);
    AliMixture(  4, "Foam$",            af, zf, df, -2, wf);
    AliMixture(  5, "Stainless Steel$", as, zs, ds, 5, ws);
    AliMaterial( 6, "Si$",              28.09, 14., 2.33, 9.36, 42.3);
    AliMixture(  7, "Thermo Insul.$",   ati, zti, dti, -2, wti);
    AliMixture(  8, "Textolit$",        atx, ztx, dtx, -2, wtx);
    AliMaterial(99, "Air$",             14.61, 7.3, .001205, 30420., 67500);

    AliMedium(700, "PHOS Xtal    $", 0, 1, ISXFLD, SXMGMX, 10., .1, .1, .1, .1);
    AliMedium(701, "CPV scint.   $", 1, 1, ISXFLD, SXMGMX, 10., .1, .1, .1, .1);
    AliMedium(702, "Al parts     $", 2, 0, ISXFLD, SXMGMX, 10., .1, .1, .001, .001);
    AliMedium(703, "Tyvek wrapper$", 3, 0, ISXFLD, SXMGMX, 10., .1, .1, .001, .001);
    AliMedium(704, "Polyst. foam $", 4, 0, ISXFLD, SXMGMX, 10., .1, .1, .1, .1);
    AliMedium(705, "Steel cover  $", 5, 0, ISXFLD, SXMGMX, 10., .1, .1, 1e-4, 1e-4);
    AliMedium(706, "Si PIN       $", 6, 0, ISXFLD, SXMGMX, 10., .1, .1, .01, .01);
    AliMedium(707, "Thermo Insul.$", 7, 0, ISXFLD, SXMGMX, 10., .1, .1, .1, .1);
    AliMedium(708, "Textolit     $", 8, 0, ISXFLD, SXMGMX, 10., .1, .1, .1, .1);
    AliMedium(799, "Air          $",99, 0, ISXFLD, SXMGMX, 10., 1., .1, .1, 10);

// --- Generate explicitly delta rays in the steel cover --- 
    pMC->Gstpar(idtmed[704], "LOSS", 3.);
    pMC->Gstpar(idtmed[704], "DRAY", 1.);
// --- and in aluminium parts --- 
    pMC->Gstpar(idtmed[701], "LOSS", 3.);
    pMC->Gstpar(idtmed[701], "DRAY", 1.);
}
 
//______________________________________________________________________________

void AliPHOS::AddPHOSCradles()
{
  Int_t i;
  for(i=0;i<GetCradlesAmount();i++) {
    
    int n = fCradles->GetEntries();
    fCradles->Add(new AliPHOSCradle( IsVersion(),            // geometry.
				     GetCrystalSideSize    (),
				     GetCrystalLength      (),
				     GetWrapThickness      (),
				     GetAirThickness       (),
				     GetPIN_SideSize       (),
				     GetPIN_Length         (),
				     GetRadius             (),
				     GetCPV_Thickness      (),
				     GetCPV_PHOS_Distance  (),
				     GetNz                 (),
				     GetNphi               (),
				     GetCradleAngle        (i)));
    
    if( n+1 != fCradles->GetEntries() || NULL == fCradles->At(n) )
      {
	cout << "  Can not create or add AliPHOSCradle.\n";
	exit(1);
      }
  }
}

//______________________________________________________________________________

Int_t AliPHOS::DistancetoPrimitive(Int_t , Int_t )
{
   return 9999;
}
 
//___________________________________________
void AliPHOS::Init()
{
  Int_t i;
  //
  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" PHOS_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  //
  // Here the ABSO initialisation code (if any!)
  for(i=0;i<80;i++) printf("*");
  printf("\n");
}

//______________________________________________________________________________

void AliPHOS::MakeBranch(Option_t *)
{
// ROOT output initialization to ROOT file.
// 
// AliDetector::MakeBranch()  is always called.
//
// There will be also special tree "PHOS" with one branch "AliPHOSCradles"
// if it was set next flag in the galice card file:
//  * PHOSflags:    YES: X<>0   NO: X=0
//  * PHOSflags(1) : -----X.  Create branch for TObjArray of AliPHOSCradle
//     Examples:
//     PHOSflags      1.
//     PHOSflags 636301.
// In that case special bit CradlesBranch_Bit will be set for AliPHOS

  AliDetector::MakeBranch();
  
  int i;
  float t = GetPHOS_flag(0)/10;
  i = (int) t;
  i = (int) ((t-i)*10);
  if( !i )
    return;

  SetBit(CradlesBranch_Bit);

  if( NULL==(fTreePHOS=new TTree(fTreeName.Data(),"PHOS events tree")) )
  {
    Error("MakeBranch","Can not create TTree");
    exit(1);
  }

  if( NULL==fTreePHOS->GetCurrentFile() )
  {
    Error("MakeBranch","There is no opened ROOT file");
    exit(1);
  }

  // Create a new branch in the current Root Tree.

  if( NULL==fTreePHOS->Branch(fBranchNameOfCradles.Data(),"TObjArray",&fCradles,4000,0) )
  {
    Error("MakeBranch","Can not create branch");
    exit(1);
  }

  printf("The branch %s has been created\n",fBranchNameOfCradles.Data());
}

//______________________________________________________________________________

void AliPHOS::SetTreeAddress(void)
{
// ROOT input initialization.
//
// AliDetector::SetTreeAddress()  is always called.
//
// If CradlesBranch_Bit is set (see AliPHOS::MakeBranch) than fTreePHOS is
// initilized.

  AliDetector::SetTreeAddress();

  if( !TestBit(CradlesBranch_Bit) )
    return;

  if( NULL==(fTreePHOS=(TTree*)gDirectory->Get((char*)(fTreeName.Data()))  ) )
  {
    Error("Can not find Tree \"%s\"\n",fTreeName.Data());
    exit(1);
  }

  TBranch *branch = fTreePHOS->GetBranch(fBranchNameOfCradles.Data());
  if( NULL==branch )
  {
    Error("SetTreeAddress","Can not find branch %s in TTree:%s",fBranchNameOfCradles.Data(),fTreeName.Data());
    exit(1);
  }

  branch->SetAddress(&fCradles);
}

//______________________________________________________________________________

AliPHOSCradle *AliPHOS::GetCradleOfTheParticle(const Hep3Vector &p,const Hep3Vector &v) const
{
// For a given direction 'p' and source point 'v' returns pointer to AliPHOSCradle
// in that direction or NULL if AliPHOSCradle was not found.

  for( int m=0; m<fCradles->GetEntries(); m++ )
  {
    AliPHOS *PHOS = (AliPHOS *)this;     // Removing 'const'...
    AliPHOSCradle *cradle = (AliPHOSCradle *)PHOS->fCradles->operator[](m);

    float x,y,l;
    const float d = cradle->GetRadius()-cradle->GetCPV_PHOS_Distance()-cradle->GetCPV_Thikness();
    cradle->GetXY(p,v,d,x,y,l);

    if( l>0 && fabs(x)<cradle->GetNz  ()*cradle->GetCellSideSize()/2 
            && fabs(y)<cradle->GetNphi()*cradle->GetCellSideSize()/2 )
      return cradle;
  }

  return NULL;
}

//______________________________________________________________________________

void AliPHOS::Reconstruction(Float_t signal_step, UInt_t min_signal_reject)
{
// Call AliPHOSCradle::Reconstruction(Float_t signal_step, UInt_t min_signal_reject)
// for all AliPHOSCradles.

  for( int i=0; i<fCradles->GetEntries(); i++ )
    GetCradle(i).Reconstruction(signal_step,min_signal_reject);
}

//______________________________________________________________________________

void AliPHOS::ResetDigits(void)
{
  AliDetector::ResetDigits();

  for( int i=0; i<fCradles->GetEntries(); i++ )
    ((AliPHOSCradle*)(*fCradles)[i]) -> Clear();
}

//______________________________________________________________________________

void AliPHOS::FinishEvent(void)
{
// Called at the end of each 'galice' event.

  if( NULL!=fTreePHOS )
    fTreePHOS->Fill();
}

//______________________________________________________________________________

void AliPHOS::FinishRun(void)
{
}

//______________________________________________________________________________

void AliPHOS::Print(Option_t *opt)
{
// Print PHOS information.
// For each AliPHOSCradle the function AliPHOSCradle::Print(opt) is called.

  AliPHOS &PHOS = *(AliPHOS *)this;     // Removing 'const'...

  for( int i=0; i<fCradles->GetEntries(); i++ )
  {
    printf("PHOS cradle %d from %d\n",i+1, fCradles->GetEntries());
    PHOS.GetCradle(i).Print(opt);
    printf( "---------------------------------------------------\n");
  }
}

//______________________________________________________________________________
void AliPHOS::SetFlags(Float_t p1,Float_t p2,Float_t p3,Float_t p4,
                       Float_t p5,Float_t p6,Float_t p7,Float_t p8,Float_t p9)
{
  PHOSflags[0]=p1;
  PHOSflags[1]=p2;
  PHOSflags[2]=p3;
  PHOSflags[3]=p4;
  PHOSflags[4]=p5;
  PHOSflags[5]=p6;
  PHOSflags[6]=p7;
  PHOSflags[7]=p8;
  PHOSflags[8]=p9;
}

//______________________________________________________________________________
void AliPHOS::SetCell(Float_t p1,Float_t p2,Float_t p3,Float_t p4,
                       Float_t p5,Float_t p6,Float_t p7,Float_t p8,Float_t p9)
{
  PHOScell[0]=p1;
  PHOScell[1]=p2;
  PHOScell[2]=p3;
  PHOScell[3]=p4;
  PHOScell[4]=p5;
  PHOScell[5]=p6;
  PHOScell[6]=p7;
  PHOScell[7]=p8;
  PHOScell[8]=p9;
}

//______________________________________________________________________________
void AliPHOS::SetRadius(Float_t radius)
{
   PHOSradius=radius;
}

//______________________________________________________________________________
void AliPHOS::SetCradleSize(Int_t nz, Int_t nphi, Int_t ncradles)
{
   PHOSsize[0]=nz;
   PHOSsize[1]=nphi;
   PHOSsize[2]=ncradles;
}

//______________________________________________________________________________
void AliPHOS::SetCradleA(Float_t angle)
{
   PHOScradlesA=angle;
}

//______________________________________________________________________________
void AliPHOS::SetCPV(Float_t p1,Float_t p2,Float_t p3,Float_t p4,
                     Float_t p5,Float_t p6,Float_t p7,Float_t p8,Float_t p9)
{
   PHOSCPV[0] = p1;
   PHOSCPV[1] = p2;
   PHOSCPV[2] = p3;
   PHOSCPV[3] = p4;
   PHOSCPV[4] = p5;
   PHOSCPV[5] = p6;
   PHOSCPV[6] = p7;
   PHOSCPV[7] = p8;
   PHOSCPV[8] = p9;
}

//______________________________________________________________________________
void AliPHOS::SetExtra(Float_t p1,Float_t p2,Float_t p3,Float_t p4,
                       Float_t p5,Float_t p6,Float_t p7,Float_t p8,Float_t p9)
{
   PHOSextra[0] = p1;
   PHOSextra[1] = p2;
   PHOSextra[2] = p3;
   PHOSextra[3] = p4;
   PHOSextra[4] = p5;
   PHOSextra[5] = p6;
   PHOSextra[6] = p7;
   PHOSextra[7] = p8;
   PHOSextra[8] = p9;
}

//______________________________________________________________________________
void AliPHOS::SetTextolitWall(Float_t dx, Float_t dy, Float_t dz)
{
   PHOSTXW[0] = dx;
   PHOSTXW[1] = dy;
   PHOSTXW[2] = dz;
}

//______________________________________________________________________________
void AliPHOS::SetInnerAir(Float_t dx, Float_t dy, Float_t dz)
{
   PHOSAIR[0] = dx;
   PHOSAIR[1] = dy;
   PHOSAIR[2] = dz;
}

//______________________________________________________________________________
void AliPHOS::SetFoam(Float_t dx, Float_t dy, Float_t dz, Float_t dr)
{
   PHOSFTI[0] = dx;
   PHOSFTI[1] = dy;
   PHOSFTI[2] = dz;
   PHOSFTI[3] = dr;
}

ClassImp(AliPHOSCradle)

//______________________________________________________________________________

AliPHOSCradle::AliPHOSCradle(void) {}

//______________________________________________________________________________

AliPHOSCradle::AliPHOSCradle( int   Geometry           ,
                              float CrystalSideSize    ,
                              float CrystalLength      ,
                              float WrapThickness      ,
                              float AirThickness       ,
                              float PIN_SideSize       ,
                              float PIN_Length         ,
                              float Radius             ,
                              float CPV_Thickness      ,
                              float CPV_PHOS_Distance  ,
                              int   Nz                 ,
                              int   Nphi               ,
                              float Angle              ) :
    fGeometry                   (Geometry),
//  fCellEnergy                 (),
//  fChargedTracksInPIN         (),
//  fCPV_hitsX                  (),
//  fCPV_hitsY                  (),
    fCrystalSideSize            (CrystalSideSize),
    fCrystalLength              (CrystalLength),
    fWrapThickness              (WrapThickness),
    fAirThickness               (AirThickness),
    fPIN_SideSize               (PIN_SideSize),
    fPIN_Length                 (PIN_Length),
    fRadius                     (Radius),
    fCPV_PHOS_Distance          (CPV_PHOS_Distance),
    fCPV_Thickness              (CPV_Thickness),
    fNz                         (Nz),
    fNphi                       (Nphi),
    fPhi                        (Angle)
{
        fCellEnergy         = TH2F("CellE","Energy deposition in a cells",fNz,0,fNz,fNphi,0,fNphi);
        fCellEnergy           .SetDirectory(0);
        fChargedTracksInPIN = TH2S("PINCtracks","Amount of charged tracks in PIN",fNz,0,fNz,fNphi,0,fNphi);
        fChargedTracksInPIN   .SetDirectory(0);
}

//______________________________________________________________________________

void AliPHOSCradle::Clear(Option_t *)
{
// Clear digit. information.

  fCellEnergy              .Reset();
  fChargedTracksInPIN      .Reset();
  GetParticles()           .Delete();
  GetParticles()           .Compress();
  GetGammasReconstructed() .Delete();
  GetGammasReconstructed() .Compress();

  fCPV_hitsX.Set(0);
  fCPV_hitsY.Set(0);
}

//______________________________________________________________________________

void AliPHOSCradle::AddCPVHit(float x,float y)
{
// Add this hit to the hits list in CPV detector.

  TArrayF a(fCPV_hitsX.GetSize()+1);
  
  memcpy(a.GetArray(),fCPV_hitsX.GetArray(),sizeof(Float_t)*fCPV_hitsX.GetSize());
  a[fCPV_hitsX.GetSize()] = x;
  fCPV_hitsX = a;

  // It must be:   fCPV_hitsX.GetSize() == fCPV_hitsY.GetSize()

  memcpy(a.GetArray(),fCPV_hitsY.GetArray(),sizeof(Float_t)*fCPV_hitsY.GetSize());
  a[fCPV_hitsY.GetSize()] = y;
  fCPV_hitsY = a;
}

//______________________________________________________________________________

void AliPHOSCradle::GetXY(const Hep3Vector &p,const Hep3Vector &v,float R,float &x,float &y,float &l) const
{
// This function calculates hit position (x,y) in the CRADLE cells plain from particle in
// the direction given by 'p' (not required to be normalized) and start point
// given by 3-vector 'v'. So the particle trajectory is   t(l) = v + p*l
// were 'l' is a number (distance from 'v' to CRADLE cells plain) and 't' is resulting
// three-vector of trajectory point.
// 
// After the call to this function user should test that l>=0 (the particle HITED the
// plain) and (x,y) are in the region of CRADLE:
// 
// Example:
//   AliPHOSCradle cradle(......);
//   Hep3Vector p(....), v(....);
//   Float_t x,y,l;
//   cradle.GetXY(p,v,x,y,l);
//   if( l<0 || fabs(x)>cradle.GetNz()  *cradle.GetCellSideSize()/2
//           || fabs(y)>cradle.GetNphi()*cradle.GetCellSideSize()/2 )
//     cout << "Outside the CRADLE.\n";

  // We have to create three vectors:
  //    s  - central point on the PHOS surface
  //    n1 - first vector in CRADLE plain
  //    n2 - second vector in CRADLE plain
  // This three vectors are orthonormalized.

  double phi = fPhi/180*M_PI;
  Hep3Vector      n1(   0        ,   0        , 1 ),   // Z direction (X)
                  n2(  -sin(phi) ,   cos(phi) , 0 ),   // around beam (Y)
                  s ( R*cos(phi) , R*sin(phi) , 0 );   // central point

  const double l1_min = 1e-2;
  double l1,
         p_n1 = p.dot(n1),        // dot() - scalar product.
         p_n2 = p.dot(n2),
         v_n1 = v.dot(n1),
         v_n2 = v.dot(n2),
         s_n1 = s.dot(n1), // 0
         s_n2 = s.dot(n2); // 0
  
  if      ( fabs(l1=p.x()-n1.x()*p_n1-n2.x()*p_n2)>l1_min )
    { l = (-v.x()+s.x()+n1.x()*(v_n1-s_n1)+n2.x()*(v_n2-s_n2))/l1; }
  else if ( fabs(l1=p.y()-n1.y()*p_n1-n2.y()*p_n2)>l1_min )
    { l = (-v.y()+s.y()+n1.y()*(v_n1-s_n1)+n2.y()*(v_n2-s_n2))/l1; }
  else if ( fabs(l1=p.z()-n1.z()*p_n1-n2.z()*p_n2)>l1_min )
    { l = (-v.z()+s.z()+n1.z()*(v_n1-s_n1)+n2.z()*(v_n2-s_n2))/l1; }

//         double lx = (-v.x()+s.x()+n1.x()*(v.dot(n1)-s.dot(n1))+n2.x()*(v.dot(n2)-s.dot(n2)))/
//                     (p.x()-n1.x()*p.dot(n1)-n2.x()*p.dot(n2)),
//                ly = (-v.y()+s.y()+n1.y()*(v.dot(n1)-s.dot(n1))+n2.y()*(v.dot(n2)-s.dot(n2)))/
//                     (p.y()-n1.y()*p.dot(n1)-n2.y()*p.dot(n2)),
//                lz = (-v.z()+s.z()+n1.z()*(v.dot(n1)-s.dot(n1))+n2.z()*(v.dot(n2)-s.dot(n2)))/
//                     (p.z()-n1.z()*p.dot(n1)-n2.z()*p.dot(n2));
//         cout.form("x: %g %g %g %g\n",lx,-v.x()+s.x()+n1.x()*(v.dot(n1)-s.dot(n1))+n2.x()*(v.dot(n2)-s.dot(n2)),p.x()-n1.x()*p.dot(n1)-n2.x()*p.dot(n2));
//         cout.form("y: %g %g %g %g\n",lx,-v.y()+s.y()+n1.y()*(v.dot(n1)-s.dot(n1))+n2.y()*(v.dot(n2)-s.dot(n2)),p.y()-n1.y()*p.dot(n1)-n2.y()*p.dot(n2));
//         cout.form("z: %g %g %g %g\n",lx,-v.z()+s.z()+n1.z()*(v.dot(n1)-s.dot(n1))+n2.z()*(v.dot(n2)-s.dot(n2)),p.z()-n1.z()*p.dot(n1)-n2.z()*p.dot(n2));
//         cout.form("lx,ly,lz =   %g,%g,%g\n",lx,ly,lz);

  x = p_n1*l + v_n1 - s_n1;
  y = p_n2*l + v_n2 - s_n2;
}

//______________________________________________________________________________

void AliPHOSCradle::Print(Option_t *opt)
{
// Print AliPHOSCradle information.
// 
// options:  'd' - print energy deposition for EVERY cell
//           'p' - print particles list that hit the cradle
//           'r' - print list of reconstructed particles

  AliPHOSCradle *cr = (AliPHOSCradle *)this;     // Removing 'const'...

  printf("AliPHOSCradle:  Nz=%d  Nphi=%d, fPhi=%f, E=%g, CPV hits amount = %d\n",fNz,fNphi,fPhi,
       cr->fCellEnergy.GetSumOfWeights(),fCPV_hitsX.GetSize());

  if( NULL!=strchr(opt,'d') )
  {
    printf("\n\nCells Energy (in MeV):\n\n   |");
    for( int x=0; x<fNz; x++ )
      printf(" %4d|",x+1);
    printf("\n");

    for( int y=fNphi-1; y>=0; y-- )
    {
      printf("%3d|",y+1);
      for( int x=0; x<fNz; x++ )
        printf("%6d",(int)(cr->fCellEnergy.GetBinContent(cr->fCellEnergy.GetBin(x,y))*1000));
      printf("\n");
    }
    printf("\n");
  }

  if( NULL!=strchr(opt,'p') )
  {
    printf("This cradle was hit by %d particles\n",
         ((AliPHOSCradle*)this)->GetParticles().GetEntries());
    TObjArray &p=((AliPHOSCradle*)this)->GetParticles();
    for( int i=0; i<p.GetEntries(); i++ )
      ((AliPHOSgamma*)(p[i]))->Print();
  }

  if( NULL!=strchr(opt,'p') )
  {
    printf("Amount of reconstructed gammas is %d\n",
         ((AliPHOSCradle*)this)->GetGammasReconstructed().GetEntries());

    TObjArray &p=((AliPHOSCradle*)this)->GetGammasReconstructed();
    for( int i=0; i<p.GetEntries(); i++ )
      ((AliPHOSgamma*)(p[i]))->Print();
  }
}

//______________________________________________________________________________

void AliPHOSCradle::Distortion(const TH2F *Noise, const TH2F *Stochastic, const TH2F *Calibration)
{
// This function changes histogram of cell energies fCellEnergy on the base of input
// histograms Noise, Stochastic, Calibration. The histograms must have
// size Nz x Nphi. 

  //////////////////////////////////
  // Testing the histograms size. //
  //////////////////////////////////
  
  if( fNz!=fCellEnergy.GetNbinsX() || fNphi!=fCellEnergy.GetNbinsY() )
  {
    printf      ("Bad size of CellEnergy!   Must be:   Nz x Nphi = %d x %d\n"
                 "but size of CellEnergy is:  %d x %d\n",
                 fNz,fNphi,fCellEnergy.GetNbinsX(),fCellEnergy.GetNbinsY());
    exit(1);
  }

  if( fNz!=fChargedTracksInPIN.GetNbinsX() || fNphi!=fChargedTracksInPIN.GetNbinsY() )
  {
    printf      ("Bad size of ChargedTracksInPIN!   Must be:   Nz x Nphi = %d x %d\n"
                 "but size of ChargedTracksInPIN is:  %d x %d\n",
                 fNz,fNphi,fChargedTracksInPIN.GetNbinsX(),fChargedTracksInPIN.GetNbinsY());
    exit(1);
  }

  if( NULL!=Noise && (fNz!=Noise->GetNbinsX() || fNphi!=Noise->GetNbinsX()) )
  {
    printf      ("Bad size of Noise!   Must be:   Nz x Nphi = %d x %d\n"
                 "but size of Noise is:  %d x %d\n",
                 fNz,fNphi,fChargedTracksInPIN.GetNbinsX(),fChargedTracksInPIN.GetNbinsY());
    exit(1);
  }

  if( NULL!=Stochastic && (fNz!=Stochastic->GetNbinsX() || fNphi!=Stochastic->GetNbinsX()) )
  {
    printf      ("Bad size of Stochastic!   Must be:   Nz x Nphi = %d x %d\n"
                 "but size of Stochastic is:  %d x %d\n",
                 fNz,fNphi,fChargedTracksInPIN.GetNbinsX(),fChargedTracksInPIN.GetNbinsY());
    exit(1);
  }

  if( NULL!=Calibration && (fNz!=Calibration->GetNbinsX() || fNphi!=Calibration->GetNbinsX()) )
  {
    printf      ("Bad size of Calibration!   Must be:   Nz x Nphi = %d x %d\n"
                 "but size of Calibration is:  %d x %d\n",
                 fNz,fNphi,fChargedTracksInPIN.GetNbinsX(),fChargedTracksInPIN.GetNbinsY());
    exit(1);
  }

  ////////////////////
  // Do distortion! //
  ////////////////////

  for( int y=0; y<fNphi; y++ )
    for( int x=0; x<fNz; x++ )
    {
      const int n = fCellEnergy.GetBin(x,y);   // Bin number
      static TRandom r;
    
      Float_t   E_old=fCellEnergy.GetBinContent(n),   E_new=E_old;

      if( NULL!=Stochastic )
        E_new   = r.Gaus(E_old,sqrt(E_old)*GetDistortedValue(Stochastic,n));

      if( NULL!=Calibration )
        E_new  *=  GetDistortedValue(Calibration,n);

      if( NULL!=Noise )
        E_new  +=  GetDistortedValue(Noise,n);

      fCellEnergy.SetBinContent(n,E_new);
    }
}

////////////////////////////////////////////////////////////////////////////////

TH2F* AliPHOSCradle::CreateHistForDistortion(const char *name, const char *title,
                                             Int_t Nx, Int_t Ny,
                                             Float_t MU_mu,    Float_t MU_sigma,
                                             Float_t SIGMA_mu, Float_t SIGMA_sigma)
{
// Create (new TH2F(...)) histogram with information (for every bin) that will
// be used for VALUE creation.
// Two values will be created for each bin:
// MU    = TRandom::Gaus(MU_mu,MU_sigma)
// and
// SIGMA = TRandom::Gaus(SIGMA_mu,SIGMA_sigma)
// The VALUE in a particluar bin will be equal
// VALUE = TRandom::Gaus(MU,SIGMA)
// 
// Do not forget to delete the histogram at the end of the work.

  TH2F *h = new TH2F( name,title, Nx,1,Nx, Ny,1,Ny );
  if( h==NULL )
  {
    Error("CreateHistForDistortion","Can not create the histogram");
    exit(1);
  }
  h->SetDirectory(0);

  for( int y=0; y<Ny; y++ )
    for( int x=0; x<Nx; x++ )
    {
      const int n = h->GetBin(x,y);
      h->SetBinContent(n,r.Gaus(   MU_mu,   MU_sigma));
      h->SetBinError  (n,r.Gaus(SIGMA_mu,SIGMA_sigma));
    }

  return h;
}

////////////////////////////////////////////////////////////////////////////////

Float_t AliPHOSCradle::GetDistortedValue(const TH2F *h, UInt_t n)
{
  return r.Gaus(((TH2F*)h)->GetBinContent(n),n);
}

////////////////////////////////////////////////////////////////////////////////
//______________________________________________________________________________

#ifdef WIN32
  #define common_for_event_storing COMMON_FOR_EVENT_STORING
#else
  #define common_for_event_storing common_for_event_storing_
#endif

extern "C" struct
{
  enum { crystals_matrix_amount_max=4, crystals_in_matrix_amount_max=40000 };

  // Event-independent information
  UShort_t      crystals_matrix_amount_PHOS,
                crystal_matrix_type,
                amount_of_crystals_on_Z,
                amount_of_crystals_on_PHI;
  Float_t       radius,
                crystal_size,
                crystal_length,
                matrix_coordinate_Z             [crystals_matrix_amount_max],
                matrix_coordinate_PHI           [crystals_matrix_amount_max];
  UInt_t        event_number;
  UShort_t      crystals_amount_with_amplitudes [crystals_matrix_amount_max],
                crystals_amplitudes_Iad         [crystals_matrix_amount_max]
                                                [crystals_in_matrix_amount_max][2];
} common_for_event_storing;

//       integer*4 crystals_amount_max,crystals_in_matrix_amount_max,
//      +          crystals_matrix_amount_max
//       parameter (crystals_matrix_amount_max=4)
//       parameter (crystals_in_matrix_amount_max=40000)
//       parameter (crystals_amount_max =crystals_matrix_amount_max*
//      +                                crystals_in_matrix_amount_max)
// 
// * All units are in GeV, cm, radian
//       real       crystal_amplitudes_unit, radius_unit,
//      +           crystal_size_unit, crystal_length_unit,
//      +           matrix_coordinate_Z_unit, matrix_coordinate_PHI_unit
//       integer    crystal_amplitudes_in_units_min
//       parameter (crystal_amplitudes_in_units_min        = 1)
//       parameter (crystal_amplitudes_unit                = 0.001 ) ! 1.0  MeV
//       parameter (radius_unit                            = 0.1   ) ! 0.1  cm
//       parameter (crystal_size_unit                      = 0.01  ) ! 0.01 cm
//       parameter (crystal_length_unit                    = 0.01  ) ! 0.01 cm
//       parameter (matrix_coordinate_Z_unit               = 0.1   ) ! 0.1  cm
//       parameter (matrix_coordinate_PHI_unit             = 1e-4  ) ! 1e-4 radian
// 
//       integer*2 crystals_matrix_amount_PHOS, crystal_matrix_type,
//      +          amount_of_crystals_on_Z, amount_of_crystals_on_PHI,
//      +          crystals_amount_with_amplitudes, crystals_amplitudes_Iad
//       integer*4 event_number
// 
//       real      radius, crystal_size, crystal_length,
//      +          matrix_coordinate_Z, matrix_coordinate_PHI
// 
//       real      crystals_amplitudes, crystals_energy_total
//       integer   event_file_unit_number
// 
//       common /common_for_event_storing/
//      + ! Event-independent information
//      +        crystals_matrix_amount_PHOS,
//      +        crystal_matrix_type,
//      +        amount_of_crystals_on_Z,
//      +        amount_of_crystals_on_PHI,
//      +        radius,
//      +        crystal_size,
//      +        crystal_length,
//      +        matrix_coordinate_Z     (crystals_matrix_amount_max),
//      +        matrix_coordinate_PHI   (crystals_matrix_amount_max),
//      +
//      + ! Event-dependent information
//      +        event_number,
//      +        crystals_amount_with_amplitudes
//      +                                (crystals_matrix_amount_max),
//      +        crystals_amplitudes_Iad (2,crystals_in_matrix_amount_max,
//      +                                 crystals_matrix_amount_max),
//      +        
//      + ! These information don't store in data file
//      +        crystals_amplitudes     (crystals_amount_max),
//      +        crystals_energy_total,
//      +        event_file_unit_number


// 	parameter (NGp=1000,nsps=10,nvertmax=1000)
//         COMMON /GAMMA/KG,MW(ngp),ID(ngp),JD(ngp),E(ngp),E4(ngp),
//      ,  XW(ngp),YW(ngp),ES(nsps,ngp),ET(nsps,ngp),ISsd(ngp),
//      ,  IGDEV(ngp),ZGDEV(ngp),sigexy(3,ngp),Emimx(2,nsps,ngp),
//      ,  kgfix,igfix(ngp),cgfix(3,ngp),sgfix(3,ngp),hiw(ngp),
//      ,  wsw(nsps,ngp),h1w(ngp),h0w(ngp),raxay(5,ngp),
//      ,  sigmaes0(nsps,ngp),dispeces(nsps,ngp),
//      ,  igamvert(ngp)


#ifdef WIN32
#define rcgamma RCGAMMA
#else
#define rcgamma rcgamma_
#endif

extern "C" struct
{
  enum {NGP=1000, nsps=10, nvertmax=1000};
  int   recons_gammas_amount, mw[NGP],ID[NGP],JD[NGP];
  float E[NGP], E4[NGP], XW[NGP], YW[NGP], ES[NGP][nsps],ET[NGP][nsps],ISsd[NGP],
        igdev[NGP],Zgdev[NGP];
//      sigexy(3,ngp),Emimx(2,nsps,ngp),
//   ,  kgfix,igfix(ngp),cgfix(3,ngp),sgfix(3,ngp),hiw(ngp),
//   ,  wsw(nsps,ngp),h1w(ngp),h0w(ngp),raxay(5,ngp),
//   ,  sigmaes0(nsps,ngp),dispeces(nsps,ngp),
//   ,  igamvert(ngp)
} rcgamma;

#ifdef WIN32
#define reconsfirst RECONSFIRST
#define type_of_call _stdcall
#else
#define reconsfirst reconsfirst_
#define type_of_call
#endif

extern "C" void type_of_call reconsfirst(const float &,const float &);

void AliPHOSCradle::Reconstruction(Float_t signal_step, UInt_t min_signal_reject)
{
// Call of PHOS reconstruction program.
// signal_step=0.001  GeV (1MeV)
// min_signal_reject = 15 or 30 MeV

  common_for_event_storing.event_number                       = 0;  // We do not know event number?
  common_for_event_storing.crystals_matrix_amount_PHOS        = 1;
  common_for_event_storing.crystal_matrix_type                = 1; // 1 - rectangular
  common_for_event_storing.amount_of_crystals_on_Z            = fNz;
  common_for_event_storing.amount_of_crystals_on_PHI          = fNphi;

  common_for_event_storing.radius                             = fRadius;
  common_for_event_storing.crystal_size                       = GetCellSideSize();
  common_for_event_storing.crystal_length                     = fCrystalLength;

  common_for_event_storing.matrix_coordinate_Z            [0] = 0;
  common_for_event_storing.matrix_coordinate_PHI          [0] = fPhi;

  #define  k    common_for_event_storing.crystals_amount_with_amplitudes[0] 
  k=0;

  for( int y=0; y<fNphi; y++ )
    for( int x=0; x<fNz; x++ )
    {
      UInt_t    n       = fCellEnergy.GetBin(x,y);
      UInt_t    signal  = (int) (fCellEnergy.GetBinContent(n)/signal_step);
      if( signal>=min_signal_reject )
      {
        common_for_event_storing.crystals_amplitudes_Iad[0][k][0] = signal;
        common_for_event_storing.crystals_amplitudes_Iad[0][k][1] = x + y*fNz;
        k++;
      }
    }
  #undef  k

  GetGammasReconstructed().Delete();
  GetGammasReconstructed().Compress();

  const float   stochastic_term   = 0.03,        // per cents over sqrt(E);  E in GeV
                electronic_noise  = 0.01;        // GeV
  reconsfirst(stochastic_term,electronic_noise); // Call of reconstruction program.
  

  for( int i=0; i<rcgamma.recons_gammas_amount; i++ )
  {
//     new (GetGammasReconstructed().UncheckedAt(i) ) AliPHOSgamma;
//     AliPHOSgamma &g = *(AliPHOSgamma*)(GetGammasReconstructed().UncheckedAt(i));

    AliPHOSgamma *gggg = new AliPHOSgamma;
    if( NULL==gggg )
    {
      Error("Reconstruction","Can not create AliPHOSgamma");
      exit(1);
    }

    GetGammasReconstructed().Add(gggg);
    AliPHOSgamma &g=*gggg;
    
    Float_t thetta, alpha, betta, R=fRadius+rcgamma.Zgdev[i]/10;

    g.fX      = rcgamma.YW[i]/10;
    g.fXsigma = 0.2;
    g.fY      = rcgamma.XW[i]/10;
    g.fYsigma = 0.2;
    g.fE      = rcgamma.E [i];
    g.fEsigma = 0.01*sqrt(rcgamma.E[i])+0.05;

    thetta      = atan(g.fX/R);

    alpha = atan(g.fY/R);
    betta = fPhi/180*M_PI + alpha;

    g.fPx = g.fE * cos(thetta) * cos(betta);
    g.fPy = g.fE * cos(thetta) * sin(betta);
    g.fPz = g.fE * sin(thetta);
  }
}

//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________

ClassImp(AliPHOSgamma)

//______________________________________________________________________________

void AliPHOSgamma::Print(Option_t *)
{
  float mass = fE*fE - fPx*fPx - fPy*fPy - fPz*fPz;

  if( mass>=0 )
    mass =  sqrt( mass);
  else
    mass = -sqrt(-mass);

  printf("XY=(%+7.2f,%+7.2f)  (%+7.2f,%+7.2f,%+7.2f;%7.2f)  mass=%8.4f\n",
          fX,fY,fPx,fPy,fPz,fE,mass);
}

//______________________________________________________________________________

AliPHOSgamma &AliPHOSgamma::operator=(const AliPHOSgamma &g)
{
  fX           = g.fX;
  fXsigma      = g.fXsigma;
  fY           = g.fY;
  fYsigma      = g.fYsigma;
  fE           = g.fE;
  fEsigma      = g.fEsigma;
  fPx          = g.fPx;
  fPy          = g.fPy;
  fPz          = g.fPz;

  return *this;
}

//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________

ClassImp(AliPHOShit)

//______________________________________________________________________________

AliPHOShit::AliPHOShit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
AliHit(shunt, track)
{
   Int_t i;
   for (i=0;i<5;i++) fVolume[i] = vol[i];
   fX       = hits[0];
   fY       = hits[1];
   fZ       = hits[2];
   fELOS    = hits[3];
}
 
//______________________________________________________________________________
