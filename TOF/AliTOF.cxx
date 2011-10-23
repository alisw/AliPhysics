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
//                                                                           //
//  VERSIONE WITH 5 SYMMETRIC MODULES ALONG Z AXIS                           //
//  ============================================================             //
//                                                                           //
//  VERSION WITH HOLES FOR PHOS AND TRD IN SPACEFRAME WITH HOLES             //
//                                                                           //
//  Volume sensibile : FPAD                                                  //
//                                                                           //
//                                                                           //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
// Begin_Html
/*
<img src="picts/AliTOFClass.gif">
*/
//End_Html

#include <TClonesArray.h>
#include <TFile.h>
#include <TFolder.h>
#include <TROOT.h>
#include <TTree.h>
#include <TVirtualMC.h>
#include <TStopwatch.h>

#include "AliConst.h"
#include "AliLoader.h"
#include "AliLog.h"
#include "AliMC.h"
#include "AliRun.h"
#include "AliDAQ.h"
#include "AliRawReader.h"

//#include "AliTOFDDLRawData.h"
#include "AliTOFDigitizer.h"
#include "AliTOFdigit.h"
#include "AliTOFhitT0.h"
#include "AliTOFhit.h"
#include "AliTOFGeometry.h"
#include "AliTOFSDigitizer.h"
#include "AliTOFSDigit.h"
#include "AliTOF.h"
#include "AliTOFrawData.h"
#include "AliTOFRawStream.h"

class AliTOFcluster;

extern TROOT *gROOT;
extern TVirtualMC *gMC;

extern AliRun *gAlice;

ClassImp(AliTOF)

//_____________________________________________________________________________
AliTOF::AliTOF():
  fFGeom(0x0),
  fSDigits(0x0),
  fNSDigits(0),
  fReconParticles(0x0),
  fIdSens(-1),
  fTZero(kFALSE),
  fTOFHoles(kTRUE),
  fTOFGeometry(0x0),
  fTOFRawWriter(AliTOFDDLRawData())
{
  //
  // Default constructor
  //

  //by default all sectors switched on
  for (Int_t ii=0; ii<18; ii++) fTOFSectors[ii]=0;

  fDigits = 0;
  fIshunt   = 0;
  fName = "TOF";

}
 
//_____________________________________________________________________________
AliTOF::AliTOF(const char *name, const char *title, Option_t *option)
       : 
  AliDetector(name,title),
  fFGeom(0x0),
  fSDigits(0x0),
  fNSDigits(0),
  fReconParticles(0x0),
  fIdSens(-1),
  fTZero(kFALSE),
  fTOFHoles(kTRUE),
  fTOFGeometry(0x0),
  fTOFRawWriter(AliTOFDDLRawData())
{
  //
  // AliTOF standard constructor
  // 
  // Here are fixed some important parameters
  //

  // Initialization of hits, sdigits and digits array
  // added option for time zero analysis
  //skowron
  fTOFGeometry = new AliTOFGeometry();

  //by default all sectors switched on
  for (Int_t ii=0; ii<18; ii++) fTOFSectors[ii]=0;

  if (strstr(option,"tzero")){
    fHits   = new TClonesArray("AliTOFhitT0",  1000);
    fTZero = kTRUE;
    //    AliWarning("tzero option requires AliTOFv4T0/AliTOFv5T0 as TOF version (check Your Config.C)");
  }else{
    fHits   = new TClonesArray("AliTOFhit",  1000);
    fTZero = kFALSE;
  }
  if (gAlice==0) {
     AliFatal("gAlice==0 !");
  }

  AliMC *mcApplication = (AliMC*)gAlice->GetMCApp();

  if (mcApplication->GetHitLists())
     mcApplication->AddHitList(fHits);
  else AliError("gAlice->GetHitLists()==0");

  fIshunt  = 0;
  fSDigits = new TClonesArray("AliTOFSDigit", 1000);
  fDigits  = new TClonesArray("AliTOFdigit",  1000);

  //
  // Digitization parameters
  //
  // (Transfer Functions to be inserted here)
  //
  //PH  SetMarkerColor(7);
  //PH  SetMarkerStyle(2);
  //PH  SetMarkerSize(0.4);

// Strip Parameters
  //fGapA    =   4.; //cm  Gap beetween tilted strip in A-type plate
  //fGapB    =   6.; //cm  Gap beetween tilted strip in B-type plate

  // Physical performances
  //fTimeRes = 100.;//ps
  //fChrgRes = 100.;//pC

}

//____________________________________________________________________________
void AliTOF::SetTOFSectors(Int_t * const sectors)
{
  // Setter for partial/full TOF configuration

  for(Int_t isec=0;isec<18;isec++){
    fTOFSectors[isec]=sectors[isec];
  }
}
//____________________________________________________________________________
void AliTOF::GetTOFSectors(Int_t *sectors) const
{
  // Getter for partial/full TOF configuration

  for(Int_t isec=0;isec<18;isec++){
    sectors[isec]=fTOFSectors[isec];
  }
}

//_____________________________________________________________________________
void AliTOF::CreateTOFFolders()
{
  // create the ALICE TFolder
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
 
  // creates the TOF geometry  folder
  geomF->AddFolder("TOF", "Geometry for TOF") ;
}

//_____________________________________________________________________________
AliTOF::~AliTOF()
{
  // dtor:
  // it remove also the alice folder 
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
void AliTOF::AddDigit(Int_t *tracks, Int_t *vol, Int_t *digits)
{
  //
  // Add a TOF digit
  // new with placement used
  //
  TClonesArray &ldigits = *fDigits;
  new (ldigits[fNdigits++]) AliTOFdigit(tracks, vol, digits);
}

//_____________________________________________________________________________
void AliTOF::AddSDigit(Int_t tracknum, Int_t *vol, Int_t *digits)
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

  if (fLoader->TreeS () )
    {
      branch = fLoader->TreeS ()->GetBranch ("TOF");
      if (branch) {
        if (fSDigits == 0x0) fSDigits = new TClonesArray("AliTOFSDigit",  1000);
        branch->SetAddress (&fSDigits);
      }
    }

  if (fLoader->TreeR() ) 
    {
      branch = fLoader->TreeR()->GetBranch("TOF"); 
      if (branch) 
       {
	 if (fReconParticles == 0x0) fReconParticles = new TClonesArray("AliTOFcluster",  1000);
         branch->SetAddress(&fReconParticles);
       }
    }

  /*
  if (fLoader->TreeR() && fReconParticles) //I do not know where this array is created - skowron
    {
      branch = fLoader->TreeR()->GetBranch("TOF"); 
      if (branch) 
       {
         branch->SetAddress(&fReconParticles) ;
       }
    }
  */
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

  Float_t xTof, yTof;

  if (IsVersion()==8) {

    xTof = 124.5;//fTOFGeometry->StripLength()+2.*(0.3+0.03); // cm,  x-dimension of FTOA volume
    yTof = fTOFGeometry->Rmax()-fTOFGeometry->Rmin(); // cm,  y-dimension of FTOA volume
    Float_t zTof = fTOFGeometry->ZlenA();             // cm,  z-dimension of FTOA volume
    
    //  TOF module internal definitions
    TOFpc(xTof, yTof, zTof);

  } else if (IsVersion()==7) {

    xTof = 124.5;//fTOFGeometry->StripLength()+2.*(0.3+0.03); // cm,  x-dimension of FTOA volume
    yTof = fTOFGeometry->Rmax()-fTOFGeometry->Rmin(); // cm,  y-dimension of FTOA volume
    Float_t zTof = fTOFGeometry->ZlenA();             // cm,  z-dimension of FTOA volume
    
    //  TOF module internal definitions
    TOFpc(xTof, yTof, zTof, fTOFGeometry->ZlenB());

  } else {

    Float_t wall = 4.;//cm // frame inbetween TOF modules

    // Sizes of TOF module with its support etc..
    xTof = 2.*(fTOFGeometry->Rmin()*TMath::Tan(10.*kDegrad)-wall/2.-0.5);
    yTof = fTOFGeometry->Rmax()-fTOFGeometry->Rmin();

    //  TOF module internal definitions 
    TOFpc(xTof, yTof, fTOFGeometry->ZlenC(), fTOFGeometry->ZlenB(), fTOFGeometry->ZlenA(), fTOFGeometry->MaxhZtof());
  }

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

  /*
  // Save the geometry
  TDirectory* saveDir = gDirectory;
  AliRunLoader::Instance()->CdGAFile();
  fTOFGeometry->Write("TOFGeometry");
  saveDir->cd();
  */
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
  const Int_t kSize=10;
  Char_t branchname[kSize];
  snprintf(branchname,kSize,"%s",GetName());
  
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

  if (fLoader->TreeR() && oR){
    if (fReconParticles == 0x0) fReconParticles = new TClonesArray("AliTOFcluster",  1000);
    MakeBranchInTree(fLoader->TreeR(), branchname, &fReconParticles,buffersize, 0) ;
  }

  /*
  if (fReconParticles && fLoader->TreeR() && oR){
    MakeBranchInTree(fLoader->TreeR(), branchname, &fReconParticles,buffersize, 0) ;
  }
  */
}

//____________________________________________________________________________
void AliTOF::Makehits(Bool_t hits) 
{
// default argument used, see AliTOF.h
// Enable/Disable the writing of the TOF-hits branch 
// on TreeH
// by default :  enabled for TOFv1, v2, v3, v4, v5
//              disabled for TOFv0
// 
   if (hits &&  (IsVersion()!=0))
      fIdSens = gMC->VolId("FPAD");
   else
      AliInfo("Option for writing the TOF-hits branch on TreeH: disabled");
}

//____________________________________________________________________________
void AliTOF::FinishEvent()
{
// do nothing
}

//____________________________________________________________________________
void AliTOF::Hits2SDigits()
{
//
// Use the TOF SDigitizer to make TOF SDigits
//

//  AliInfo("Start...");
  
  AliRunLoader * rl = fLoader->GetRunLoader();
  AliDebug(2,"Initialized runLoader");
  AliTOFSDigitizer sd((rl->GetFileName()).Data());
  AliDebug(2,"Initialized TOF sdigitizer");
  //ToAliDebug(1, sd.Print(""));
  //AliInfo("ToAliDebug");

  //sd.Digitize("all") ;
  sd.Digitize("partial") ;

  AliDebug(2,"I am sorting from AliTOF class");

}

//____________________________________________________________________________
void AliTOF::Hits2SDigits(Int_t evNumber1, Int_t evNumber2)
{
//
// Use the TOF SDigitizer to make TOF SDigits
//

  if ((evNumber2-evNumber1)==1) 
    AliDebug(1, Form("I am making sdigits for the %dth event", evNumber1));
  if ((evNumber2-evNumber1)>1)
    AliDebug(1, Form("I am making sdigits for the events from the %dth to the %dth", evNumber1, evNumber2-1));
 
  AliRunLoader * rl = fLoader->GetRunLoader();
  AliTOFSDigitizer sd((rl->GetFileName()).Data(),evNumber1,evNumber2) ;
  ToAliDebug(1, sd.Print(""));

  sd.Digitize("") ;

}

//___________________________________________________________________________
AliDigitizer* AliTOF::CreateDigitizer(AliDigitizationInput* digInput) const
{
  AliDebug(2,"I am creating the TOF digitizer");
  return new AliTOFDigitizer(digInput);
}

//___________________________________________________________________________
Bool_t AliTOF::CheckOverlap(const Int_t * const vol,
			    Int_t* digit,Int_t Track)
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
      Int_t tdc2 = digit[0];
      Int_t tdc1 = currentDigit->GetTdc();

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
void AliTOF::Digits2Raw()
{
//
// Starting from the TOF digits, writes the Raw Data objects
//

  TStopwatch stopwatch;
  stopwatch.Start();

  fLoader->LoadDigits();

  TTree* digits = fLoader->TreeD();
  if (!digits) {
    AliError("no digits tree");
    return;
  }
  
  //AliTOFDDLRawData rawWriter;
  fTOFRawWriter.Clear();
  fTOFRawWriter.SetVerbose(0);
  if (fTOFRawWriter.GetPackedAcquisitionMode()) {
    if(fTOFRawWriter.GetMatchingWindow()>8192)
      AliWarning(Form("You are running in packing mode and the matching window is %.2f ns, i.e. greater than 199.8848 ns",
		      fTOFRawWriter.GetMatchingWindow()*AliTOFGeometry::TdcBinWidth()*1.e-03));
  }
  
  AliDebug(1,"Formatting raw data for TOF");
  digits->GetEvent(0);
  fTOFRawWriter.RawDataTOF(digits->GetBranch("TOF"));  

  fLoader->UnloadDigits();
  
  AliDebug(1, Form("Execution time to write TOF raw data : R:%.2fs C:%.2fs",
		   stopwatch.RealTime(),stopwatch.CpuTime()));

}

//____________________________________________________________________________
void AliTOF::RecreateSDigitsArray() {
//
// delete TClonesArray fSDigits and create it again
//  needed for backward compatability with PPR test production
//
  fSDigits->Clear();
}
//____________________________________________________________________________
void AliTOF::CreateSDigitsArray() {
//
// create TClonesArray fSDigits
//  needed for backward compatability with PPR test production
//
  fSDigits       = new TClonesArray("AliTOFSDigit",  1000);
}
//____________________________________________________________________________
Bool_t AliTOF::Raw2SDigits(AliRawReader* rawReader)
{
  //
  // Converts raw data to sdigits for TOF
  //

  TStopwatch stopwatch;
  stopwatch.Start();

  if(!GetLoader()->TreeS()) {MakeTree("S");  MakeBranch("S");}
  //TClonesArray &aSDigits = *fSDigits;

  AliTOFRawStream tofRawStream = AliTOFRawStream();
  tofRawStream.Raw2SDigits(rawReader, fSDigits);

  GetLoader()->TreeS()->Fill(); GetLoader()->WriteSDigits("OVERWRITE");//write out sdigits
  Int_t nSDigits = fSDigits->GetEntries();

  ResetSDigits();

  AliDebug(1, Form("Got %d TOF sdigits", nSDigits));
  AliDebug(1, Form("Execution time to read TOF raw data and fill TOF sdigit tree : R:%.2fs C:%.2fs",
		   stopwatch.RealTime(),stopwatch.CpuTime()));

  return kTRUE;

}

//____________________________________________________________________________
void AliTOF::Raw2Digits(AliRawReader* rawReader)
{
  //
  // Converts raw data to digits for TOF
  //

  TStopwatch stopwatch;
  stopwatch.Start();

  if(!GetLoader()->TreeD()) {MakeTree("D");  MakeBranch("D");}
  //TClonesArray &aDigits = *fDigits;

  AliTOFRawStream tofRawStream = AliTOFRawStream();
  tofRawStream.Raw2Digits(rawReader, fDigits);

  GetLoader()->TreeD()->Fill(); GetLoader()->WriteDigits("OVERWRITE");//write out digits
  Int_t nDigits = fDigits->GetEntries();

  ResetDigits();

  AliDebug(1, Form("Got %d TOF digits", nDigits));
  AliDebug(1, Form("Execution time to read TOF raw data and fill TOF digit tree : R:%.2fs C:%.2fs",
		   stopwatch.RealTime(),stopwatch.CpuTime()));

}
