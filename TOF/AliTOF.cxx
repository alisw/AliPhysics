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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Of Flight                                                           //
//  This class contains the basic functions for the Time Of Flight           //
//  detector. Functions specific to one particular geometry are              //
//  contained in the derived classes                                         //
//
//  VERSIONE WITH 5 SYMMETRIC MODULES ALONG Z AXIS
//  ==============================================
//  
//  VERSION WITH HOLES FOR PHOS AND TRD IN SPACEFRAME WITH HOLES
//
//  Volume sensibile : FPAD
//
//
//
// Begin_Html
/*
<img src="picts/AliTOFClass.gif">
*/
//End_Html
//             
//
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <Riostream.h>
#include <Rstrstream.h>
#include <stdlib.h>

#include <TBRIK.h>
#include <TFile.h>
#include <TFolder.h>
#include <TNode.h>
#include <TObject.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTask.h>
#include <TTree.h>
#include <TVirtualMC.h>

#include "AliConst.h"
#include "AliLoader.h"
#include "AliMagF.h"
#include "AliRun.h"
#include "AliTOF.h"
#include "AliTOFRawDigit.h"
#include "AliTOFRawSector.h"
#include "AliTOFRoc.h"
#include "AliTOFSDigit.h"
#include "AliTOFdigit.h"
#include "AliTOFhit.h"
#include "AliTOFhitT0.h"
#include "AliMC.h"
#include "AliTOFDigitizer.h"
 
ClassImp(AliTOF)
 
//_____________________________________________________________________________
AliTOF::AliTOF()
{
  //
  // Default constructor
  //
  fFGeom = 0x0;
  fDTask = 0x0;
  fReTask = 0x0;
  fIshunt   = 0;
  fSDigits       = 0 ;
  fNSDigits = 0;
  fDigits        = 0 ;
  fReconParticles = 0x0;
  fName="TOF";
  fMerger = 0x0;
  fTZero = kFALSE;
  fTOFGeometry = 0;
}
 
//_____________________________________________________________________________
AliTOF::AliTOF(const char *name, const char *title, Option_t *option)
       : AliDetector(name,title)
{
  //
  // AliTOF standard constructor
  // 
  // Here are fixed some important parameters
  //

  // Initialization of hits, sdigits and digits array
  // added option for time zero analysis
  fFGeom = 0x0; //skowron
  fDTask = 0x0;
  fReTask= 0x0;
  fReconParticles= 0x0;
  fMerger = 0x0;
  fTOFGeometry = 0;

  if (strstr(option,"tzero")){
    fHits   = new TClonesArray("AliTOFhitT0",  1000);
    fTZero = kTRUE;
    cout << "tzero option requires AliTOFv4T0 as TOF version (check Your Config.C)" << endl;
  }else{
    fHits   = new TClonesArray("AliTOFhit",  1000);
    fTZero = kFALSE;
  }
  if (gAlice==0) {
     Fatal("AliTOF","gAlice==0 !");
  }
  if (gAlice->GetMCApp()->GetHitLists())
     gAlice->GetMCApp()->AddHitList(fHits);
  else Error("AliTOF","gAlice->GetHitLists()==0");

  fIshunt  = 0;
  fSDigits       = new TClonesArray("AliTOFSDigit",  1000);
  fDigits        = new TClonesArray("AliTOFdigit",  1000);
  fNSDigits = 0;

  fFGeom = 0x0;
  fDTask = 0x0;
  fReTask = 0x0;
  fReconParticles = 0x0;
  fMerger = 0x0;

  //
  // Digitization parameters
  //
  // (Transfer Functions to be inserted here)
  //
  SetMarkerColor(7);
  SetMarkerStyle(2);
  SetMarkerSize(0.4);

// General Geometrical Parameters
  fNTof    =  18;  // number of sectors
  fRmax    = 399.0;//cm 
  fRmin    = 370.0;//cm
  fZlenC   = 177.5;//cm length of module C
  fZlenB   = 141.0;//cm length of module B
  fZlenA   = 106.0;//cm length of module A
  fZtof    = 371.5;//cm total semi-length of TOF detector

// Strip Parameters
  fStripLn = 122.0;//cm  Strip Length
  fSpace   =   5.5;//cm  Space Beetween the strip and the bottom of the plate 
  fDeadBndZ=   1.5;//cm  Dead Boundaries of a Strip along Z direction (width)
  fDeadBndX=   1.0;//cm  Dead Boundaries of a Strip along X direction (length)
  fXpad    =   2.5;//cm  X size of a pad
  fZpad    =   3.5;//cm  Z size of a pad
  fGapA    =   4.; //cm  Gap beetween tilted strip in A-type plate
  fGapB    =   6.; //cm  Gap beetween tilted strip in B-type plate
  fOverSpc =  15.3;//cm Space available for sensitive layers in radial direction
  fNpadX   =  48;  // Number of pads in a strip along the X direction
  fNpadZ   =   2;  // Number of pads in a strip along the Z direction
  fPadXStr = fNpadX*fNpadZ; //Number of pads per strip
  fNStripA = 15; // number of strips in A type module 
  fNStripB = 19; // number of strips in B type module
  fNStripC = 20; // number of strips in C type module
 
// Physical performances
  fTimeRes = 100.;//ps
  fChrgRes = 100.;//pC

// DAQ characteristics
// cfr. TOF-TDR pag. 105 for Glossary
// TARODA : TOF-ALICE Read Out and Data Acquisition system
  fPadXSector = 8928; // number of pad per sector -with no holes-
                      // ((15+2*19+2*20)*(48*2))
  fNRoc       = 14;   // number of Roc (Read Out Controller) (TARODA)
  fNFec       = 32;   // number of Fec (Front-End electronic Card)
                      // (TARODA)
  fNTdc       = 32;   // number of Tdc (Time to Digital Converter)
  fNPadXRoc   = (Int_t)fPadXSector/fNRoc; // number of pads for each ROC
  
}

//_____________________________________________________________________________
void AliTOF::CreateTOFFolders()
{
  // create the ALICE TFolder
  // create the ALICE TTasks
  // create the ALICE main TFolder
  // to be done by AliRun

  TFolder * alice = new TFolder();
  alice->SetNameTitle("FPAlice", "Alice Folder") ;
  gROOT->GetListOfBrowsables()->Add(alice) ;

  TFolder * aliceF  = alice->AddFolder("folders", "Alice memory Folder") ;
  //  make it the owner of the objects that it contains
  aliceF->SetOwner() ;
  // geometry folder
  TFolder * geomF = aliceF->AddFolder("Geometry", "Geometry objects") ;
  TFolder * aliceT  = alice->AddFolder("tasks", "Alice tasks Folder") ;   
  //  make it the owner of the objects that it contains
  aliceT->SetOwner() ;

  TTask * aliceDi = new TTask("(S)Digitizer", "Alice SDigitizer & Digitizer") ;
  aliceT->Add(aliceDi);

  TTask * aliceRe = new TTask("Reconstructioner", "Alice Reconstructioner") ;
  aliceT->Add(aliceRe);

  char * tempo = new char[80] ;

  // creates the TOF Digitizer and adds it to alice main (S)Digitizer task
  sprintf(tempo, "%sDigitizers container",GetName() ) ;
  fDTask = new TTask(GetName(), tempo);
  aliceDi->Add(fDTask) ;

  // creates the TOF reconstructioner and adds it to alice main Reconstructioner task
  sprintf(tempo, "%sReconstructioner container",GetName() ) ;
  fReTask = new TTask(GetName(), tempo);
  aliceRe->Add(fReTask) ;

  delete [] tempo ;
 
  // creates the TOF geometry  folder
  geomF->AddFolder("TOF", "Geometry for TOF") ;
}

//_____________________________________________________________________________
AliTOF::~AliTOF()
{
  // dtor:
  // it remove also the alice folder 
  // and task that TOF creates instead of AliRun
  /* PH Temporarily commented because of problems
  TFolder * alice = (TFolder*)gROOT->GetListOfBrowsables()->FindObject("FPAlice") ;
  delete alice;
  alice = 0;
  */
  if (fHits)
    {
      fHits->Delete ();
      delete fHits;
      fHits = 0;
    }
  if (fDigits)
    {
      fDigits->Delete ();
      delete fDigits;
      fDigits = 0;
    }
  if (fSDigits)
    {
      fSDigits->Delete();
      delete fSDigits;
      fSDigits = 0;
    }
  if (fReconParticles)
    {
      fReconParticles->Delete ();
      delete fReconParticles;
      fReconParticles = 0;
    }

}

//_____________________________________________________________________________
void AliTOF::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a TOF hit
  // new with placement used
  //
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliTOFhit(fIshunt, track, vol, hits);
}

//_____________________________________________________________________________
void AliTOF::AddT0Hit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a TOF hit
  // new with placement used
  //
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliTOFhitT0(fIshunt, track, vol, hits);
}

//_____________________________________________________________________________
void AliTOF::AddDigit(Int_t *tracks, Int_t *vol, Float_t *digits)
{
  //
  // Add a TOF digit
  // new with placement used
  //
  TClonesArray &ldigits = *fDigits;
  new (ldigits[fNdigits++]) AliTOFdigit(tracks, vol, digits);
}

//___________________________________________
void AliTOF::AddSDigit(Int_t tracknum, Int_t *vol, Float_t *digits)
{
     
//
// Add a TOF sdigit
//
        
  TClonesArray &lSDigits = *fSDigits;   
  new(lSDigits[fNSDigits++]) AliTOFSDigit(tracknum, vol, digits);
}

//_____________________________________________________________________________
void AliTOF::SetTreeAddress ()
{
  // Set branch address for the Hits and Digits Tree.
  char branchname[30];
  
  if (fLoader->TreeH())
   {
     if (fHits == 0x0)
      {
        if (fTZero) fHits   = new TClonesArray("AliTOFhitT0", 1000);
        else fHits   = new TClonesArray("AliTOFhit", 1000);
      }
   }
  AliDetector::SetTreeAddress ();

  TBranch *branch;
  TTree *treeD = fLoader->TreeD();

  if (treeD)
    {
      branch = treeD->GetBranch (branchname);
      if (branch)
       {
         if (fDigits == 0x0) fDigits = new TClonesArray("AliTOFdigit",  1000); 
         branch->SetAddress (&fDigits);
       }
    }

  if (fLoader->TreeS () )
    {
      branch = fLoader->TreeS ()->GetBranch ("TOF");
      if (branch) {
        if (fSDigits == 0x0) fSDigits = new TClonesArray("AliTOFSDigit",  1000);
        branch->SetAddress (&fSDigits);
      }
    }

  if (fLoader->TreeR() && fReconParticles) //I do not know where this array is created - skowron
    {
      branch = fLoader->TreeR()->GetBranch("TOF"); 
      if (branch) 
       {
         branch->SetAddress(&fReconParticles) ;
       }
    }   
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
  const Double_t kDegrad=kPi/180.;
  //
  Float_t xTof, yTof, wall;

  // frame inbetween TOF modules
  wall = 4.;//cm

  // Sizes of TOF module with its support etc..
  xTof = 2.*(fRmin*TMath::Tan(10*kDegrad)-wall/2-.5);
  yTof = fRmax-fRmin;

//  TOF module internal definitions 
  TOFpc(xTof, yTof, fZlenC, fZlenB, fZlenA, fZtof);
}

//_____________________________________________________________________________
void AliTOF::DrawModule() const
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
  gMC->Gsatt("FTOA","SEEN",1);
  gMC->Gsatt("FTOB","SEEN",1);
  gMC->Gsatt("FTOC","SEEN",1);
  gMC->Gsatt("FLTA","SEEN",1);
  gMC->Gsatt("FLTB","SEEN",1);
  gMC->Gsatt("FLTC","SEEN",1);
  gMC->Gsatt("FSTR","SEEN",1);
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
  // Revision: F. Pierella 18-VI-2002
  //

  Int_t   isxfld = gAlice->Field()->Integ();
  Float_t sxmgmx = gAlice->Field()->Max();
  //
  //--- Quartz (SiO2) to simulate float glass
  //    density tuned to have correct float glass 
  //    radiation length
  Float_t   aq[2] = { 28.0855,15.9994 };
  Float_t   zq[2] = { 14.,8. };
  Float_t   wq[2] = { 1.,2. };
  Float_t   dq = 2.55; // std value: 2.2
  Int_t nq = -2;

  // --- Freon C2F4H2 (TOF-TDR pagg.)
  // Geant Manual CONS110-1, pag. 43 (Geant, Detector Description and Simulation Tool)
  Float_t afre[3]  = {12.011,18.998,1.007};
  Float_t zfre[3]  = { 6., 9., 1.}; 
  Float_t wfre[3]  = { 2., 4., 2.};
  Float_t densfre  = 0.00375;   
// http://www.fi.infn.it/sezione/prevprot/gas/freon.html
  Int_t nfre = -3; 
/*
  //-- Isobutane quencher C4H10 (5% in the sensitive mixture)
  Float_t aiso[2]  = {12.011,1.007};
  Float_t ziso[2]  = { 6.,  1.};
  Float_t wiso[2]  = { 4., 10.};
  Float_t densiso  = .......;  // (g/cm3) density
  Int_t nfre = -2; // < 0 i.e. proportion by number of atoms of each kind
  //-- SF6 (5% in the sensitive mixture)
  Float_t asf[3]  = {32.066,18.998};
  Float_t zsf[3]  = { 16., 9.};
  Float_t wsf[3]  = {  1., 6.}; 
  Float_t denssf  = .....;   // (g/cm3) density
  Int_t nfre = -2; // < 0 i.e. proportion by number of atoms of each kind
*/
  // --- CO2 
  Float_t ac[2]   = {12.,16.};
  Float_t zc[2]   = { 6., 8.};
  Float_t wc[2]   = { 1., 2.};
  Float_t dc = .001977;
  Int_t nc = -2;
   // For mylar (C5H4O2) 
  Float_t amy[3] = { 12., 1., 16. };
  Float_t zmy[3] = {  6., 1.,  8. };
  Float_t wmy[3] = {  5., 4.,  2. };
  Float_t dmy    = 1.39;
  Int_t nmy = -3;
 // For polyethilene (CH2) - honeycomb -
  Float_t ape[2] = { 12., 1. };
  Float_t zpe[2] = {  6., 1. };
  Float_t wpe[2] = {  1., 2. };
  Float_t dpe    = 0.935*0.479; //To have 1%X0 for 1cm as for honeycomb
  Int_t npe = -2;
  // --- G10 
  Float_t ag10[4] = { 12.,1.,16.,28. };
  Float_t zg10[4] = {  6.,1., 8.,14. };
  Float_t wmatg10[4] = { .259,.288,.248,.205 };
  Float_t densg10  = 1.7;
  Int_t nlmatg10 = -4;

  // plexiglass CH2=C(CH3)CO2CH3
  Float_t aplex[3] = { 12.,1.,16.};
  Float_t zplex[3] = {  6.,1., 8.};
  Float_t wmatplex[3] = {5.,8.,2.};
  Float_t densplex  =1.16;
  Int_t nplex = -3;

  // ---- ALUMINA (AL203) 
  Float_t aal[2] = { 27.,16.};
  Float_t zal[2] = { 13., 8.};
  Float_t wmatal[2] = { 2.,3. };
  Float_t densal  = 2.3;
  Int_t nlmatal = -2;
  // -- Water
  Float_t awa[2] = {  1., 16. };
  Float_t zwa[2] = {  1.,  8. };
  Float_t wwa[2] = {  2.,  1. };
  Float_t dwa    = 1.0;
  Int_t nwa = -2;

// stainless steel
  Float_t asteel[4] = { 55.847,51.9961,58.6934,28.0855 };
  Float_t zsteel[4] = { 26.,24.,28.,14. };
  Float_t wsteel[4] = { .715,.18,.1,.005 };

  //
  //AliMaterial(0, "Vacuum$", 1e-16, 1e-16, 1e-16, 1e16, 1e16);
  AliMaterial( 1, "Air$",14.61,7.3,0.001205,30423.24,67500.);
  AliMaterial( 2, "Cu $",  63.54, 29.0, 8.96, 1.43, 14.8);
  AliMaterial( 3, "C  $",  12.01,  6.0, 2.265,18.8, 74.4);
  AliMixture ( 4, "Polyethilene$", ape, zpe, dpe, npe, wpe);
  AliMixture ( 5, "G10$", ag10, zg10, densg10, nlmatg10, wmatg10);
  AliMixture ( 6, "PLE$", aplex, zplex, densplex, nplex, wmatplex);
  AliMixture ( 7, "CO2$", ac, zc, dc, nc, wc);
  AliMixture ( 8, "ALUMINA$", aal, zal, densal, nlmatal, wmatal);
  AliMaterial( 9, "Al $", 26.98, 13., 2.7, 8.9, 37.2);
  AliMaterial(10, "C-TRD$", 12.01, 6., 2.265*18.8/69.282*15./100, 18.8, 74.4); // for 15%
  AliMixture (11, "Mylar$",  amy, zmy, dmy, nmy, wmy);
  AliMixture (12, "Freon$",  afre, zfre, densfre, nfre, wfre);
  AliMixture (13, "Glass$", aq, zq, dq, nq, wq);
  AliMixture (14, "Water$",  awa, zwa, dwa, nwa, wwa);
  AliMixture (15, "STAINLESS STEEL$", asteel, zsteel, 7.88, 4, wsteel);

  Float_t epsil, stmin, deemax, stemax;
 
  //   Previous data
  //       EPSIL  =  0.1   ! Tracking precision, 
  //       STEMAX = 0.1      ! Maximum displacement for multiple scattering
  //       DEEMAX = 0.1    ! Maximum fractional energy loss, DLS 
  //       STMIN  = 0.1 
  //
  //   New data  
  epsil  = .001;  // Tracking precision,
  stemax = -1.;   // Maximum displacement for multiple scattering
  deemax = -.3;   // Maximum fractional energy loss, DLS
  stmin  = -.8;

  AliMedium( 1, "Air$"  ,  1, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 2, "Cu $"  ,  2, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 3, "C  $"  ,  3, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 4, "Pol$"  ,  4, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 5, "G10$"  ,  5, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 6, "PLE$"  ,  6, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 7, "CO2$"  ,  7, 0, isxfld, sxmgmx, 10., -.01, -.1, .01, -.01);
  AliMedium( 8,"ALUMINA$", 8, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium( 9,"Al Frame$",9, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(10, "DME-S$",  6, 1, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(11, "C-TRD$", 10, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(12, "Myl$"  , 11, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(13, "Fre$"  , 12, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(14, "Fre-S$", 12, 1, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(15, "Glass$", 13, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(16, "Water$", 14, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
  AliMedium(17, "STEEL$", 15, 0, isxfld, sxmgmx, 10., stemax, deemax, epsil, stmin);
}

//_____________________________________________________________________________
Int_t AliTOF::DistancetoPrimitive(Int_t , Int_t ) const
{
  //
  // Returns distance from mouse pointer to detector, default version
  //
  return 9999;
}

//___________________________________________
void AliTOF::ResetHits ()
{
  // Reset number of clusters and the cluster array for this detector
  AliDetector::ResetHits ();
}

//____________________________________________
void AliTOF::ResetDigits ()
{
  //
  // Reset number of digits and the digits array for this detector
  AliDetector::ResetDigits ();
  //
} 
//____________________________________________
void AliTOF::ResetSDigits ()
{
  //
  // Reset number of sdigits and the sdigits array for this detector
  fNSDigits = 0;
  //fSDigits = 0x0;
  //
} 
//_____________________________________________________________________________
void AliTOF::Init()
{
  //
  // Initialise TOF detector after it has been built
  //
  // Set id of TOF sensitive volume
  if (IsVersion() !=0) fIdSens=gMC->VolId("FPAD");
  //
}

//____________________________________________________________________________
void AliTOF::MakeBranch(Option_t* option)
{
 //
 // Initializes the Branches of the TOF inside the 
 // trees written for each event. 
 // AliDetector::MakeBranch initializes just the 
 // Branch inside TreeH. Here we add the branches in 
 // TreeD, TreeS and TreeR.
 //
  const char *oH = strstr(option,"H");
  if (fLoader->TreeH() && oH)
   {
     if (fHits == 0x0)
      {
        if (fTZero) fHits   = new TClonesArray("AliTOFhitT0", 1000);
        else fHits   = new TClonesArray("AliTOFhit", 1000);
      }
   }
  
  AliDetector::MakeBranch(option);

  Int_t buffersize = 4000;
  Char_t branchname[10];
  sprintf(branchname,"%s",GetName());
  
  const char *oD = strstr(option,"D");
  const char *oS = strstr(option,"S");
  const char *oR = strstr(option,"R");

  if (fLoader->TreeD() && oD){
    if (fDigits == 0x0) fDigits = new TClonesArray("AliTOFdigit",  1000); 
    MakeBranchInTree(fLoader->TreeD(), branchname, &fDigits,buffersize, 0) ;
  }

  if (fLoader->TreeS() && oS){
    if (fSDigits == 0x0) fSDigits = new TClonesArray("AliTOFSDigit",  1000);
    MakeBranchInTree(fLoader->TreeS(), branchname, &fSDigits,buffersize, 0) ;
  }

  if (fReconParticles && fLoader->TreeR() && oR){
    MakeBranchInTree(fLoader->TreeR(), branchname, &fReconParticles,buffersize, 0) ;
  }
}

//____________________________________________________________________________
void AliTOF::Makehits(Bool_t hits) 
{
// default argument used, see AliTOF.h
// Enable/Disable the writing of the TOF-hits branch 
// on TreeH
// by default : enabled for TOFv1, v2, v3, v4
//              disabled for TOFv0
// 
   if (hits &&  (IsVersion()!=0))
      fIdSens = gMC->VolId("FPAD");
   else
      cout << "Option for writing the TOF-hits branch on TreeH: disabled" << endl;
}

//____________________________________________________________________________
void AliTOF::FinishEvent()
{
// do nothing
}

//___________________________________________
void AliTOF::SDigits2Digits()
{
//
// Generate digits performing merging
//
  /*
    int nparticles = gAlice->GetNtrack();
    cout << "Particles       :" <<nparticles<<endl;
    if (nparticles > 0 ) {
      
      AliTOF::Hits2Digits();
      
    }
  */
  cout<<"AliTOF::SDigits2Digits"<<endl; 
    if (fMerger) {
      fMerger->Init();
      cout<<"AliTOF::SDigits2Digits Init"<<endl; 
      fMerger->Digitise();
      cout<<"AliTOF::SDigits2Digits Digitise() "<<endl; 
     }
}

//---------------------------------------------------------------------
void   AliTOF::SetMerger(AliTOFMerger* merger)
{
// Set pointer to merger
    fMerger = merger;
}

//---------------------------------------------------------------------
AliTOFMerger*  AliTOF::Merger()
{
// Return pointer to merger
    return fMerger;
}

//---------------------------------------------------------------------
void AliTOF::Hits2SDigits()
{
//
// Use the TOF SDigitizer to make TOF SDigits
//

//  cout<<"AliTOF::Hits2SDigits> start...\n";
  
  AliRunLoader * rl = fLoader->GetRunLoader();
  AliTOFSDigitizer sd((rl->GetFileName()).Data());
  if (GetDebug()) sd.Print("");

  sd.Exec("") ;

}

//---------------------------------------------------------------------
void AliTOF::Hits2SDigits(Int_t evNumber1, Int_t evNumber2)
{
//
// Use the TOF SDigitizer to make TOF SDigits
//

  if ((evNumber2-evNumber1)==1) 
      cout << "<AliTOF::Hits2SDigits>: I am making sdigits for the " << evNumber1 << "th event \n";
  else if ((evNumber2-evNumber1)>1)
      cout << "<AliTOF::Hits2SDigits>: I am making sdigits for the events from the " 
           << evNumber1 << "th to the " << evNumber2-1 << "th \n";
 
  AliRunLoader * rl = fLoader->GetRunLoader();
  AliTOFSDigitizer sd((rl->GetFileName()).Data(),evNumber1,evNumber2) ;
  if (GetDebug()) sd.Print("");

  sd.Exec("") ;

}

//___________________________________________________________________________
AliDigitizer* AliTOF::CreateDigitizer(AliRunDigitizer* manager) const
{
  return new AliTOFDigitizer(manager);
}

//___________________________________________________________________________
Bool_t AliTOF::CheckOverlap(Int_t* vol, Float_t* digit,Int_t Track)
{
//
// Checks if 2 or more hits belong to the same pad.
// In this case the data assigned to the digit object
// are the ones of the first hit in order of Time.
// 2 hits from the same track on the same pad are collected.
// Called only by Hits2SDigits.
// This procedure has to be optimized in the next TOF release.
//

        Bool_t overlap = kFALSE;
        Int_t  vol2[5];

        for (Int_t ndig=0; ndig<fSDigits->GetEntries(); ndig++){
	   AliTOFdigit* currentDigit = (AliTOFdigit*)(fSDigits->UncheckedAt(ndig));
           currentDigit->GetLocation(vol2);
           Bool_t idem= kTRUE;
	   // check on digit volume
           for (Int_t i=0;i<=4;i++){
	       if (!idem) break;
               if (vol[i]!=vol2[i]) idem=kFALSE;}

           if (idem){  // same pad fired
 	      Float_t tdc2 = digit[0];
              Float_t tdc1 = currentDigit->GetTdc();

	      // we separate two digits on the same pad if
	      // they are separated in time by at least 25 ns
	      // remember that tdc time is given in ps

              if (TMath::Abs(tdc1-tdc2)<25000){
		  // in case of overlap we take the earliest
		  if (tdc1>tdc2){
                   currentDigit->SetTdc(tdc2); 
                   currentDigit->SetAdc(digit[1]);
		  }
		  else {
		   currentDigit->SetTdc(tdc1);
		   currentDigit->SetAdc(digit[1]);
		  }
                  currentDigit->AddTrack(Track); // add track number in the track array
                  overlap = kTRUE;
		  return overlap;
              } else 
		overlap= kFALSE;

           } // close if (idem) -> two digits on the same TOF pad

        } // end loop on existing sdigits

        return overlap;
}


//____________________________________________________________________________
void AliTOF::Digits2Raw(Int_t evNumber)
{
//
// Starting from digits, writes the 
// Raw Data objects, i.e. a 
// TClonesArray of 18 AliTOFRawSector objects
//

  TTree* tD;

  // do nothing if no particles
  Int_t nparticles = gAlice->GetEvent(evNumber); 
  if (nparticles <= 0) return;

  tD = gAlice->TreeD();
  
  TClonesArray* tofdigits = this->Digits();
  Int_t ndigits = tofdigits->GetEntriesFast();

  TClonesArray* rawsectors = new TClonesArray("AliTOFRawSector",fNTof+2); 

  for (Int_t isect=1;isect<=fNTof;isect++){
     AliTOFRawSector* currentSector = (AliTOFRawSector*)rawsectors->UncheckedAt(isect);
     TClonesArray* rocData = (TClonesArray*)currentSector->GetRocData();

     for (Int_t digit=0; digit<ndigits; digit++){
        AliTOFdigit* currentDigit = (AliTOFdigit*)tofdigits->UncheckedAt(digit);
        Int_t sector = currentDigit->GetSector();
        if (sector==isect){
	    Int_t   pad    = currentDigit -> GetTotPad();
	    Int_t   roc    = (Int_t)(pad/fNPadXRoc)-1;
	    if (roc>=fNRoc) printf("Wrong n. of ROC ! Roc = %i",roc);
            Int_t   padRoc = (Int_t) pad%fNPadXRoc;
	    Int_t   fec    = (Int_t)(padRoc/fNFec)-1;
            Int_t   tdc    = (Int_t)(padRoc%fNFec)-1;
            Float_t time   = currentDigit->GetTdc();
            Float_t charge = currentDigit->GetAdc();
	    AliTOFRoc* currentROC = (AliTOFRoc*)rocData->UncheckedAt(roc);
	    Int_t error    = 0;
            currentROC->AddItem(fec, tdc, error, charge, time);
	} // close if (sector==isect) i.e. end loop on digits for the current sector
     } // end loop on TOF digits
     
     UInt_t totSize=16,rocSize=0;
     UInt_t rocHead[14],rocChek[14];
     UInt_t globalCheckSum=0;

     for (UInt_t iRoc = 1; iRoc<(UInt_t)fNRoc; iRoc++){
        AliTOFRoc* currentRoc = (AliTOFRoc*)rocData->UncheckedAt(iRoc); 
	rocSize  = currentRoc->GetItems()*2+1;
	totSize += rocSize*4;
	if (rocSize>=TMath::Power(2,16)) rocSize=0;
	rocHead[iRoc]   = iRoc<<28;
	rocHead[iRoc]  += rocSize;
	rocChek[iRoc]   = currentRoc->GetCheckSum();
        Int_t headCheck = currentRoc->BitCount(rocHead[iRoc]);
	globalCheckSum += headCheck;
	globalCheckSum += rocChek[iRoc];
     }
     
     AliTOFRoc* dummyRoc = new AliTOFRoc();
     totSize *= 4;
     if (totSize>=TMath::Power(2,24)) totSize=0;
     UInt_t header = totSize;
     UInt_t sectId = ((UInt_t)isect)<<24;
     header += sectId;
     globalCheckSum += dummyRoc->BitCount(header);
     currentSector->SetGlobalCS(globalCheckSum);
     currentSector->SetHeader(header);
  }  
}
 
//____________________________________________________________________________
void AliTOF::Raw2Digits(Int_t evNumber)
{
//
//  Converts Raw Data objects into digits objects.
//  We schematize the raw data with a 
//  TClonesArray of 18 AliTOFRawSector objects
//

  TTree    *tD;
  Int_t    vol[5];
  Int_t    tracks[3];
  Float_t  digit[2];
 
  tracks[0]=0;
  tracks[1]=0;
  tracks[2]=0;
 
  Int_t nparticles = gAlice->GetEvent(evNumber); 
  if (nparticles <= 0) return;

  tD = gAlice->TreeD();
  
  TClonesArray* rawsectors = new TClonesArray("AliTOFRawSector",fNTof+2);
  
  for(Int_t nSec=1; nSec<=fNTof; nSec++){
     AliTOFRawSector* currentSector = (AliTOFRawSector*)rawsectors->UncheckedAt(nSec);
     TClonesArray* rocData = (TClonesArray*)currentSector->GetRocData();
     for(Int_t nRoc=1; nRoc<=14; nRoc++){
        AliTOFRoc* currentRoc = (AliTOFRoc*)rocData->UncheckedAt(nRoc);
        Int_t currentItems = currentRoc->GetItems();
        for(Int_t item=1; item<currentItems; item++){ 
           Int_t nPad = currentRoc->GetTotPad(item);        
	   vol[0] = nSec;
	   Int_t nStrip = (Int_t)(nPad/fPadXStr)+1;
	   Int_t nPlate = 5;
	   if (nStrip<=fNStripC+2*fNStripB+fNStripA) nPlate = 4;
	   if (nStrip<=fNStripC+fNStripB+fNStripA)   nPlate = 3;
	   if (nStrip<=fNStripC+fNStripB)            nPlate = 2;
	   if (nStrip<=fNStripC)                     nPlate = 1;
	   vol[1] = nPlate;
	   switch (nPlate){
	   case 1: break;
	   case 2: nStrip -= (fNStripC);
	           break;
	   case 3: nStrip -= (fNStripC+fNStripB);
	           break;
	   case 4: nStrip -= (fNStripC+fNStripB+fNStripA);
	           break;
	   case 5: nStrip -= (fNStripC+2*fNStripB+fNStripA);
	           break;
	   }
           vol[2] = nStrip;
           Int_t pad = nPad%fPadXStr;
	   if (pad==0) pad=fPadXStr;
	   Int_t nPadX=0, nPadZ=0;
	   (pad>fNpadX)? nPadX -= fNpadX : nPadX = pad ;
	   vol[3] = nPadX;
	   (pad>fNpadX)? nPadZ = 2 : nPadZ = 1 ;
	   vol[4] = nPadZ;
	   UInt_t error=0;
	   Float_t tdc = currentRoc->GetTime(item,error);
	   if (!error) digit[0]=tdc;
	   digit[1] = currentRoc->GetCharge(item);
	   AddDigit(tracks,vol,digit);
        }
     }
  }
  tD->Fill();
  tD->Write(0,TObject::kOverwrite);
} 

////////////////////////////////////////////////////////////////////////
void AliTOF::RecreateSDigitsArray() {
//
// delete TClonesArray fSDigits and create it again
//  needed for backward compatability with PPR test production
//
  delete fSDigits;
  fSDigits       = new TClonesArray("AliTOFSDigit",  1000);
}
////////////////////////////////////////////////////////////////////////
void AliTOF::CreateSDigitsArray() {
//
// create TClonesArray fSDigits
//  needed for backward compatability with PPR test production
//
  fSDigits       = new TClonesArray("AliTOFSDigit",  1000);
}

