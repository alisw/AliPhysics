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
//  Creates and handles digits from TRD hits                                 //
//  Author: C. Blume (C.Blume@gsi.de)                                        //
//                                                                           //
//  The following effects are included:                                      //
//      - Diffusion                                                          //
//      - ExB effects                                                        //
//      - Gas gain including fluctuations                                    //
//      - Pad-response (simple Gaussian approximation)                       //
//      - Time-response                                                      //
//      - Electronics noise                                                  //
//      - Electronics gain                                                   //
//      - Digitization                                                       //
//      - ADC threshold                                                      //
//  The corresponding parameter can be adjusted via the various              //
//  Set-functions. If these parameters are not explicitly set, default       //
//  values are used (see Init-function).                                     //
//  As an example on how to use this class to produce digits from hits       //
//  have a look at the macro hits2digits.C                                   //
//  The production of summable digits is demonstrated in hits2sdigits.C      //
//  and the subsequent conversion of the s-digits into normal digits is      //
//  explained in sdigits2digits.C.                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>

#include <TMath.h>
#include <TVector.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TList.h>
#include <TTask.h>

#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliConfig.h"
#include "AliMagF.h"
#include "AliRunDigitizer.h"
#include "AliRunLoader.h"
#include "AliLoader.h"

#include "AliTRD.h"
#include "AliTRDhit.h"
#include "AliTRDdigitizer.h"
#include "AliTRDdataArrayI.h"
#include "AliTRDdataArrayF.h"
#include "AliTRDsegmentArray.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDgeometry.h"
#include "AliTRDparameter.h"

ClassImp(AliTRDdigitizer)

//_____________________________________________________________________________
AliTRDdigitizer::AliTRDdigitizer()
{
  //
  // AliTRDdigitizer default constructor
  //

  fRunLoader          = 0;
  fDigitsManager      = 0;
  fSDigitsManager     = 0;
  fSDigitsManagerList = 0;
  fTRD                = 0;
  fGeo                = 0;
  fPar                = 0;
  fEvent              = 0;
  fMasks              = 0;
  fCompress           = kTRUE;
  fDebug              = 0;
  fSDigits            = kFALSE;
  fSDigitsScale       = 0.0;
  fMergeSignalOnly    = kFALSE;
  fSimpleSim          = kFALSE;
  fSimpleDet          = 0;
 
}

//_____________________________________________________________________________
AliTRDdigitizer::AliTRDdigitizer(const Text_t *name, const Text_t *title)
                :AliDigitizer(name,title)
{
  //
  // AliTRDdigitizer constructor
  //

  fRunLoader          = 0;

  //NewIO: These data members probably are not needed anymore
  fDigitsManager      = 0;
  fSDigitsManager     = 0;
  fSDigitsManagerList = 0;
  fTRD                = 0;
  fGeo                = 0;
  fPar                = 0;
  //End NewIO comment
  fEvent              = 0;
  fMasks              = 0;
  fCompress           = kTRUE;
  fDebug              = 0;
  fSDigits            = kFALSE;
  fSDigitsScale       = 100.; // For the summable digits
  fMergeSignalOnly    = kFALSE;
  fSimpleSim          = kFALSE;
  fSimpleDet          = 0;
 

}

//_____________________________________________________________________________
AliTRDdigitizer::AliTRDdigitizer(AliRunDigitizer *manager
                                , const Text_t *name, const Text_t *title)
                :AliDigitizer(manager,name,title)
{
  //
  // AliTRDdigitizer constructor
  //

  fRunLoader          = 0;
  fDigitsManager      = 0;
  fSDigitsManager     = 0;
  fSDigitsManagerList = 0;
  fTRD                = 0;
  fGeo                = 0;
  fPar                = 0;
  fEvent              = 0;
  fMasks              = 0;
  fCompress           = kTRUE;
  fDebug              = 0;
  fSDigits            = kFALSE;
  fSDigitsScale       = 100.; // For the summable digits
  fMergeSignalOnly    = kFALSE;
  fSimpleSim          = kFALSE;
  fSimpleDet          = 0;
 

}

//_____________________________________________________________________________
AliTRDdigitizer::AliTRDdigitizer(AliRunDigitizer *manager)
                :AliDigitizer(manager,"AliTRDdigitizer","TRD digitizer")
{
  //
  // AliTRDdigitizer constructor
  //


  fRunLoader          = 0;
  fDigitsManager      = 0;
  fSDigitsManager     = 0;
  fSDigitsManagerList = 0;
  fTRD                = 0;
  fGeo                = 0;
  fPar                = 0;
  fEvent              = 0;
  fMasks              = 0;
  fCompress           = kTRUE;
  fDebug              = 0;
  fSDigits            = kFALSE;
  fSDigitsScale       = 100.;  // For the summable digits
  fMergeSignalOnly    = kFALSE;
  fSimpleSim          = kFALSE;
  fSimpleDet          = 0;


}

//_____________________________________________________________________________
AliTRDdigitizer::AliTRDdigitizer(const AliTRDdigitizer &d):AliDigitizer(d)
{
  //
  // AliTRDdigitizer copy constructor
  //

  ((AliTRDdigitizer &) d).Copy(*this);

}

//_____________________________________________________________________________
AliTRDdigitizer::~AliTRDdigitizer()
{
  //
  // AliTRDdigitizer destructor
  //


  if (fDigitsManager) {
    delete fDigitsManager;
    fDigitsManager = 0;
  }

  fSDigitsManager = 0;

  if (fSDigitsManagerList) {
    fSDigitsManagerList->Delete();
    delete fSDigitsManagerList;
    fSDigitsManagerList = 0;
  }

  if (fMasks) {
    delete [] fMasks;
    fMasks = 0;
  }

}

//_____________________________________________________________________________
AliTRDdigitizer &AliTRDdigitizer::operator=(const AliTRDdigitizer &d)
{
  //
  // Assignment operator
  //

  if (this != &d) ((AliTRDdigitizer &) d).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDdigitizer::Copy(TObject &d)
{
  //
  // Copy function
  //

  ((AliTRDdigitizer &) d).fRunLoader          = 0;
  ((AliTRDdigitizer &) d).fDigitsManager      = 0;
  ((AliTRDdigitizer &) d).fSDigitsManager     = 0;
  ((AliTRDdigitizer &) d).fSDigitsManagerList = 0;
  ((AliTRDdigitizer &) d).fTRD                = 0;
  ((AliTRDdigitizer &) d).fGeo                = 0;
  ((AliTRDdigitizer &) d).fPar                = 0;
  ((AliTRDdigitizer &) d).fEvent              = 0;
  ((AliTRDdigitizer &) d).fMasks              = 0;
  ((AliTRDdigitizer &) d).fCompress           = fCompress;
  ((AliTRDdigitizer &) d).fDebug              = fDebug  ;
  ((AliTRDdigitizer &) d).fSDigits            = fSDigits;
  ((AliTRDdigitizer &) d).fSDigitsScale       = fSDigitsScale;
  ((AliTRDdigitizer &) d).fMergeSignalOnly    = fMergeSignalOnly;
  ((AliTRDdigitizer &) d).fSimpleSim          = fSimpleSim;
  ((AliTRDdigitizer &) d).fSimpleDet          = fSimpleDet;
                                       
}

//_____________________________________________________________________________
void AliTRDdigitizer::Exec(Option_t* option)
{
  //
  // Executes the merging
  //

  Int_t iInput;

  AliTRDdigitsManager *sdigitsManager;

  TString optionString = option;
  if (optionString.Contains("deb")) {
    fDebug = 1;
    if (optionString.Contains("2")) {
      fDebug = 2;
    }
    printf("<AliTRDdigitizer::Exec> ");
    printf("Called with debug option %d\n",fDebug);
  }

  // The AliRoot file is already connected by the manager
  AliRunLoader* inrl;
  
  if (gAlice) 
   {
    if (fDebug > 0) {
      printf("<AliTRDdigitizer::Exec> ");
      printf("AliRun object found on file.\n");
    }
   }
  else {
    inrl = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(0));
    inrl->LoadgAlice();
    gAlice = inrl->GetAliRun();
    if (!gAlice)
     {
       printf("<AliTRDdigitizer::Exec> ");
       printf("Could not find AliRun object.\n");
       return;
     }
  }
                                                                           
  Int_t nInput = fManager->GetNinputs();
  fMasks = new Int_t[nInput];
  for (iInput = 0; iInput < nInput; iInput++) {
    fMasks[iInput] = fManager->GetMask(iInput);
  }

  // Initialization

  AliRunLoader* orl = AliRunLoader::GetRunLoader(fManager->GetOutputFolderName());
  if (InitDetector()) {
    AliLoader* ogime = orl->GetLoader("TRDLoader");

    TTree* tree = 0;
    if (fSDigits)
      { 
	//if we produce SDigits
	tree = ogime->TreeS();
	if (!tree)
	  {
	    ogime->MakeTree("S");
	    tree = ogime->TreeS();
	  }
      }
    else
      {//if we produce Digits
	tree = ogime->TreeD();
	if (!tree)
	  {
	    ogime->MakeTree("D");
	    tree = ogime->TreeD();
	  }
      }
    MakeBranch(tree);
  }
 
  for (iInput = 0; iInput < nInput; iInput++) {

    if (fDebug > 0) {
      printf("<AliTRDdigitizer::Exec> ");
      printf("Add input stream %d\n",iInput);
    }

    // check if the input tree exists
    inrl = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(iInput));
    AliLoader* gime = inrl->GetLoader("TRDLoader");

    TTree * treees =  gime->TreeS();
    if (treees == 0x0) 
     {
      if (gime->LoadSDigits())
       {
         Error("Exec","Error Occured while loading S. Digits for input %d.",iInput);
         return;
       }
      treees =  gime->TreeS();
     }
    
    if (treees == 0x0) {
      printf("<AliTRDdigitizer::Exec> ");
      printf("Input stream %d does not exist\n",iInput);
      return;
    } 

    // Read the s-digits via digits manager
    sdigitsManager = new AliTRDdigitsManager();
    sdigitsManager->SetDebug(fDebug);
    sdigitsManager->SetSDigits(kTRUE);
    
    AliRunLoader* rl = AliRunLoader::GetRunLoader(fManager->GetInputFolderName(iInput));
    AliLoader* gimme = rl->GetLoader("TRDLoader");
    if (!gimme->TreeS()) gimme->LoadSDigits();
    sdigitsManager->ReadDigits(gimme->TreeS());

    // Add the s-digits to the input list 
    AddSDigitsManager(sdigitsManager);

  }

  // Convert the s-digits to normal digits
  if (fDebug > 0) {
    printf("<AliTRDdigitizer::Exec> ");
    printf("Do the conversion\n");
  }
  SDigits2Digits();

  // Store the digits
  if (fDebug > 0) {
    printf("<AliTRDdigitizer::Exec> ");
    printf("Write the digits\n");
  }
  
  fDigitsManager->WriteDigits();

  //Write parameters
  orl->CdGAFile();
  if (!gFile->Get("TRDparameter")) GetParameter()->Write();

  if (fDebug > 0) {
    printf("<AliTRDdigitizer::Exec> ");
    printf("Done\n");
  }

  DeleteSDigitsManager();

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::Open(const Char_t *file, Int_t nEvent)
{
  //
  // Opens a ROOT-file with TRD-hits and reads in the hit-tree
  //

  // Connect the AliRoot file containing Geometry, Kine, and Hits
  

  TString evfoldname = AliConfig::GetDefaultEventFolderName();
  fRunLoader = AliRunLoader::GetRunLoader(evfoldname);
  if (!fRunLoader)
    fRunLoader = AliRunLoader::Open(file,AliConfig::GetDefaultEventFolderName(),
				    "UPDATE");
  
  if (!fRunLoader)
   {
     Error("Open","Can not open session for file %s.",file);
     return kFALSE;
   }
   
  if (!fRunLoader->GetAliRun()) fRunLoader->LoadgAlice();
  gAlice = fRunLoader->GetAliRun();
  
  if (gAlice) {
    if (fDebug > 0) {
      printf("<AliTRDdigitizer::Open> ");
      printf("AliRun object found on file.\n");
    }
  }
  else {
    printf("<AliTRDdigitizer::Open> ");
    printf("Could not find AliRun object.\n");
    return kFALSE;
  }

  fEvent = nEvent;

  // Import the Trees for the event nEvent in the file
  fRunLoader->GetEvent(fEvent);
  
  AliLoader* loader = fRunLoader->GetLoader("TRDLoader");
  if (!loader)
   {
     Error("Open","Can not get TRD loader from Run Loader");
     return kFALSE;
   }
  
  if (InitDetector()) {
    TTree* tree = 0;
    if (fSDigits)
     { 
     //if we produce SDigits
       tree = loader->TreeS();
       if (!tree)
        {
         loader->MakeTree("S");
         tree = loader->TreeS();
        }
     }
    else
     {//if we produce Digits
       if (!tree)
        {
         loader->MakeTree("D");
         tree = loader->TreeD();
        }
     }
    return MakeBranch(tree);
  }
  else {
    return kFALSE;
  }

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::InitDetector()
{
  //
  // Sets the pointer to the TRD detector and the geometry
  //

  // Get the pointer to the detector class and check for version 1
  fTRD = (AliTRD *) gAlice->GetDetector("TRD");
  if (!fTRD) {
    printf("<AliTRDdigitizer::InitDetector> ");
    printf("No TRD module found\n");
    exit(1);
  }
  if (fTRD->IsVersion() != 1) {
    printf("<AliTRDdigitizer::InitDetector> ");
    printf("TRD must be version 1 (slow simulator).\n");
    exit(1);
  }

  // Get the geometry
  fGeo = fTRD->GetGeometry();
  if (fDebug > 0) {
    printf("<AliTRDdigitizer::InitDetector> ");
    printf("Geometry version %d\n",fGeo->IsVersion());
  }

  // Create a digits manager
  delete fDigitsManager;
  fDigitsManager = new AliTRDdigitsManager();
  fDigitsManager->SetSDigits(fSDigits);
  fDigitsManager->CreateArrays();
  fDigitsManager->SetEvent(fEvent);
  fDigitsManager->SetDebug(fDebug);

  // The list for the input s-digits manager to be merged
  if (fSDigitsManagerList) {
    fSDigitsManagerList->Delete();
  } else {
    fSDigitsManagerList = new TList();
  }

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::MakeBranch(TTree* tree) const
{
  // 
  // Create the branches for the digits array
  //

  return fDigitsManager->MakeBranch(tree);

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::MakeDigits()
{
  //
  // Creates digits.
  //

  ///////////////////////////////////////////////////////////////
  // Parameter 
  ///////////////////////////////////////////////////////////////

  // Converts number of electrons to fC
  const Double_t kEl2fC  = 1.602E-19 * 1.0E15; 

  ///////////////////////////////////////////////////////////////

  // Number of pads included in the pad response
  const Int_t kNpad  = 3;

  // Number of track dictionary arrays
  const Int_t kNDict = AliTRDdigitsManager::kNDict;

  // Half the width of the amplification region
  const Float_t kAmWidth = AliTRDgeometry::AmThick() / 2.;

  Int_t   iRow, iCol, iTime, iPad;
  Int_t   iDict  = 0;
  Int_t   nBytes = 0;

  Int_t   totalSizeDigits = 0;
  Int_t   totalSizeDict0  = 0;
  Int_t   totalSizeDict1  = 0;
  Int_t   totalSizeDict2  = 0;

  Int_t   timeTRDbeg = 0;
  Int_t   timeTRDend = 1;

  Float_t pos[3];
  Float_t rot[3];
  Float_t xyz[3];
  Float_t padSignal[kNpad];
  Float_t signalOld[kNpad];

  AliTRDdataArrayF *signals = 0;
  AliTRDdataArrayI *digits  = 0;
  AliTRDdataArrayI *dictionary[kNDict];

  // Create a default parameter class if none is defined
  if (!fPar) {
    fPar = new AliTRDparameter("TRDparameter","Standard TRD parameter");
    if (fDebug > 0) {
      printf("<AliTRDdigitizer::MakeDigits> ");
      printf("Create the default parameter object\n");
    }
  }

  // Create a container for the amplitudes
  AliTRDsegmentArray *signalsArray 
                     = new AliTRDsegmentArray("AliTRDdataArrayF"
                                             ,AliTRDgeometry::Ndet());

  if (fPar->TRFOn()) {
    timeTRDbeg = ((Int_t) (-fPar->GetTRFlo() / fPar->GetTimeBinSize())) - 1;
    timeTRDend = ((Int_t) ( fPar->GetTRFhi() / fPar->GetTimeBinSize())) - 1;
    if (fDebug > 0) {
      printf("<AliTRDdigitizer::MakeDigits> ");
      printf("Sample the TRF between -%d and %d\n",timeTRDbeg,timeTRDend);
    }
  }

  Float_t elAttachProp = fPar->GetElAttachProp() / 100.; 

  if (!fGeo) {
    printf("<AliTRDdigitizer::MakeDigits> ");
    printf("No geometry defined\n");
    return kFALSE;
  }

  if (fDebug > 0) {
    printf("<AliTRDdigitizer::MakeDigits> ");
    printf("Start creating digits.\n");
  }

  AliLoader* gimme = fRunLoader->GetLoader("TRDLoader");
  if (!gimme->TreeH()) gimme->LoadHits();
  TTree* hitTree = gimme->TreeH();
  if (hitTree == 0x0)
    {
      Error("MakeDigits","Can not get TreeH");
      return kFALSE;
    }
  fTRD->SetTreeAddress();
  
  // Get the number of entries in the hit tree
  // (Number of primary particles creating a hit somewhere)
  Int_t nTrack = 1;
  if (!fSimpleSim) {
    nTrack = (Int_t) hitTree->GetEntries();
    if (fDebug > 0) {
      printf("<AliTRDdigitizer::MakeDigits> ");
      printf("Found %d primary particles\n",nTrack);
    } 
  }

  Int_t detectorOld = -1;
  Int_t countHits   =  0; 

  if (fDebug > 0) {
    printf("<AliTRDdigitizer::MakeDigits> ");
    printf("Driftvelocity = %.2f, Sampling = %.2fns\n",fPar->GetDriftVelocity(),
	   fPar->GetTimeBinSize() / fPar->GetDriftVelocity() * 1000.);
    printf("<AliTRDdigitizer::MakeDigits> ");
    printf("Gain = %d, Noise = %d\n",(Int_t)fPar->GetGasGain(),(Int_t)fPar->GetNoise());
    if (fPar->TimeStructOn()) {
      printf("<AliTRDdigitizer::MakeDigits> ");
      printf("Time Structure of drift cells implemented.\n");
    } else {
      printf("<AliTRDdigitizer::MakeDigits> ");
      printf("Constant drift velocity in drift cells.\n");      
    }
  }

  // Loop through all entries in the tree
  for (Int_t iTrack = 0; iTrack < nTrack; iTrack++) {

    if (!fSimpleSim) {   
      gAlice->ResetHits();
      nBytes += hitTree->GetEvent(iTrack);
    }

    // Loop through the TRD hits
    Int_t iHit = 0;
    AliTRDhit *hit = (AliTRDhit *) fTRD->FirstHit(-1);
    while (hit) {
 
      countHits++;
      iHit++;

              pos[0]      = hit->X();
              pos[1]      = hit->Y();
              pos[2]      = hit->Z();
      Float_t q           = hit->GetCharge();
      Int_t   track       = hit->Track();
      Int_t   detector    = hit->GetDetector();
      Int_t   plane       = fGeo->GetPlane(detector);
      Int_t   sector      = fGeo->GetSector(detector);
      Int_t   chamber     = fGeo->GetChamber(detector);
      Int_t   nRowMax     = fPar->GetRowMax(plane,chamber,sector);
      Int_t   nColMax     = fPar->GetColMax(plane);
      Int_t   nTimeMax    = fPar->GetTimeMax();
      Int_t   nTimeBefore = fPar->GetTimeBefore();
      Int_t   nTimeAfter  = fPar->GetTimeAfter();
      Int_t   nTimeTotal  = fPar->GetTimeTotal();
      Float_t row0        = fPar->GetRow0(plane,chamber,sector);
      Float_t col0        = fPar->GetCol0(plane);
      Float_t time0       = fPar->GetTime0(plane);
      Float_t rowPadSize  = fPar->GetRowPadSize(plane,chamber,sector);
      Float_t colPadSize  = fPar->GetColPadSize(plane);
      Float_t timeBinSize = fPar->GetTimeBinSize();
      Float_t divideRow   = 1.0 / rowPadSize;
      Float_t divideCol   = 1.0 / colPadSize;
      Float_t divideTime  = 1.0 / timeBinSize;

      if (fDebug > 1) {
        printf("Analyze hit no. %d ",iHit);
        printf("-----------------------------------------------------------\n");
        hit->Dump();
        printf("plane = %d, sector = %d, chamber = %d\n"
              ,plane,sector,chamber);
        printf("nRowMax = %d, nColMax = %d, nTimeMax = %d\n" 
              ,nRowMax,nColMax,nTimeMax);
        printf("nTimeBefore = %d, nTimeAfter = %d, nTimeTotal = %d\n"
	      ,nTimeBefore,nTimeAfter,nTimeTotal);
        printf("row0 = %f, col0 = %f, time0 = %f\n"
              ,row0,col0,time0);
        printf("rowPadSize = %f, colPadSize = %f, timeBinSize = %f\n"
	       ,rowPadSize,colPadSize,timeBinSize); 
      }
       
      // Don't analyze test hits and switched off detectors
      if ((CheckDetector(plane,chamber,sector)) &&
          (((Int_t) q) != 0)) {

        if (detector != detectorOld) {

          if (fDebug > 1) {
            printf("<AliTRDdigitizer::MakeDigits> ");
            printf("Get new container. New det = %d, Old det = %d\n"
                  ,detector,detectorOld);
	  }
          // Compress the old one if enabled
          if ((fCompress) && (detectorOld > -1)) {
            if (fDebug > 1) {
              printf("<AliTRDdigitizer::MakeDigits> ");
              printf("Compress the old container ...");
	    }
            signals->Compress(1,0);
            for (iDict = 0; iDict < kNDict; iDict++) {
              dictionary[iDict]->Compress(1,0);
	    }
            if (fDebug > 1) printf("done\n");
	  }
	  // Get the new container
          signals = (AliTRDdataArrayF *) signalsArray->At(detector);
          if (signals->GetNtime() == 0) {
            // Allocate a new one if not yet existing
            if (fDebug > 1) {
              printf("<AliTRDdigitizer::MakeDigits> ");
              printf("Allocate a new container ... ");
	    }
            signals->Allocate(nRowMax,nColMax,nTimeTotal);
	  }
          else if (fSimpleSim) {
            // Clear an old one for the simple simulation
            if (fDebug > 1) {
              printf("<AliTRDdigitizer::MakeDigits> ");
              printf("Clear a old container ... ");
            }
            signals->Clear();
          }
          else {
	    // Expand an existing one
            if (fCompress) {
              if (fDebug > 1) {
                printf("<AliTRDdigitizer::MakeDigits> ");
                printf("Expand an existing container ... ");
	      }
              signals->Expand();
	    }
	  }
	  // The same for the dictionary
          if (!fSimpleSim) {       
            for (iDict = 0; iDict < kNDict; iDict++) {       
              dictionary[iDict] = fDigitsManager->GetDictionary(detector,iDict);
              if (dictionary[iDict]->GetNtime() == 0) {
                dictionary[iDict]->Allocate(nRowMax,nColMax,nTimeTotal);
	      }
              else {
                if (fCompress) dictionary[iDict]->Expand();
	      }
	    }
          }      
          if (fDebug > 1) printf("done\n");
          detectorOld = detector;
        }

        // Rotate the sectors on top of each other
        if (fSimpleSim) {
          rot[0] = pos[0];
          rot[1] = pos[1];
          rot[2] = pos[2];
        }
        else {
          fGeo->Rotate(detector,pos,rot);
	}

        // The driftlength. It is negative if the hit is in the 
        // amplification region.
        Float_t driftlength = time0 - rot[0];

        // Take also the drift in the amplification region into account
        Float_t driftlengthL = TMath::Abs(driftlength + kAmWidth);
        if (fPar->ExBOn()) driftlengthL /= TMath::Sqrt(fPar->GetLorentzFactor());

        // Loop over all electrons of this hit
        // TR photons produce hits with negative charge
        Int_t nEl = ((Int_t) TMath::Abs(q));
        for (Int_t iEl = 0; iEl < nEl; iEl++) {

          xyz[0] = rot[0];
          xyz[1] = rot[1];
          xyz[2] = rot[2];

	  // Stupid patch to take care of TR photons that are absorbed
	  // outside the chamber volume. A real fix would actually need
	  // a more clever implementation of the TR hit generation
          Float_t zz = xyz[2] - row0;
          if ((zz < 0.0) || (zz > rowPadSize*nRowMax)) {
            if (iEl == 0) {
              printf("<AliTRDdigitizer::MakeDigits> ");
              printf("Hit outside of sensitive volume, row (Q = %d)\n",((Int_t) q));
	    }
            continue;
	  }
          Int_t tt = ((Int_t) (10*(driftlength + 2.0*kAmWidth)));
          if (tt < 0) {
            if (iEl == 0) {
              printf("<AliTRDdigitizer::MakeDigits> ");
              printf("Hit outside of sensitive volume, time (Q = %d)\n",((Int_t) q));
	    }
            continue;
	  }

          // Electron attachment
          if (fPar->ElAttachOn()) {
            if (gRandom->Rndm() < (driftlengthL * elAttachProp)) 
              continue;
          }

          // Apply the diffusion smearing
          if (fPar->DiffusionOn()) {
            if (!(fPar->Diffusion(driftlengthL,xyz))) continue;
	  }

          // Apply E x B effects (depends on drift direction)
          if (fPar->ExBOn()) { 
            if (!(fPar->ExB(driftlength+kAmWidth,xyz))) continue;   
	  }

	  // Apply the drift time correction
	  if (fPar->TimeStructOn()) {
	    // Get z-position with respect to anode wire:
	    Float_t Z  =  xyz[2] - row0 + fPar->GetAnodeWireOffset();
	    Z -= ((Int_t)(2*Z))/2.;
	    if (Z>0.25)   Z  = 0.5-Z;
	    //if (Z>0.248) cout << Z << ": " << time0 - xyz[0] + kAmWidth << " -> ";
            if (!(fPar->TimeStruct(driftlength+2*kAmWidth,Z,xyz))) continue;
	    //if (Z>0.248) cout << time0 - xyz[0] + kAmWidth << endl;
	  }

          // The electron position after diffusion and ExB in pad coordinates 
          // The pad row (z-direction)
          Float_t rowDist   = xyz[2] - row0;
          Int_t   rowE      = ((Int_t) (rowDist * divideRow));
          if ((rowE < 0) || (rowE >= nRowMax)) continue;   
          Float_t rowOffset = ((((Float_t) rowE) + 0.5) * rowPadSize) - rowDist;

          // The pad column (rphi-direction)
          Float_t col0tilt  = fPar->Col0Tilted(col0,rowOffset,plane);
          Float_t colDist   = xyz[1] - col0tilt;
          Int_t   colE      = ((Int_t) (colDist * divideCol));
          if ((colE < 0) || (colE >= nColMax)) continue;   
          Float_t colOffset = ((((Float_t) colE) + 0.5) * colPadSize) - colDist;    

          // The time bin (negative for hits in the amplification region)
	  // In the amplification region the electrons drift from both sides
	  // to the middle (anode wire plane)
          Float_t timeDist   = time0 - xyz[0];
          Float_t timeOffset = 0;
          Int_t   timeE      = 0;
          if (timeDist > 0) {
	    // The time bin
            timeE      = ((Int_t) (timeDist * divideTime));
            // The distance of the position to the middle of the timebin
            timeOffset = ((((Float_t) timeE) + 0.5) * timeBinSize) - timeDist;
	  }
          else {
	    // Difference between half of the amplification gap width and
	    // the distance to the anode wire
            Float_t anodeDist = kAmWidth - TMath::Abs(timeDist + kAmWidth);
            // The time bin
            timeE      = -1 * (((Int_t ) (anodeDist * divideTime)) + 1);
            // The distance of the position to the middle of the timebin
            timeOffset = ((((Float_t) timeE) + 0.5) * timeBinSize) + anodeDist;
	  }
 
          // Apply the gas gain including fluctuations
          Float_t ggRndm = 0.0;
          do {
            ggRndm = gRandom->Rndm();
	  } while (ggRndm <= 0);
          Int_t signal = (Int_t) (-fPar->GetGasGain() * TMath::Log(ggRndm));

          // Apply the pad response 
          if (fPar->PRFOn()) {
  	    // The distance of the electron to the center of the pad 
	    // in units of pad width
            Float_t dist = - colOffset * divideCol;
            if (!(fPar->PadResponse(signal,dist,plane,padSignal))) continue;
	  }
	  else {
            padSignal[0] = 0.0;
            padSignal[1] = signal;
            padSignal[2] = 0.0;
	  }

	  // Sample the time response inside the drift region
	  // + additional time bins before and after.
          // The sampling is done always in the middle of the time bin
          for (Int_t iTimeBin = TMath::Max(timeE-timeTRDbeg,        -nTimeBefore) 
  	            ;iTimeBin < TMath::Min(timeE+timeTRDend,nTimeMax+nTimeAfter ) 
        	    ;iTimeBin++) {

     	    // Apply the time response
            Float_t timeResponse = 1.0;
            Float_t crossTalk    = 0.0;
            Float_t time         = (iTimeBin - timeE) * timeBinSize + timeOffset;
            if (fPar->TRFOn()) {
              timeResponse = fPar->TimeResponse(time);
	    }
            if (fPar->CTOn()) {
              crossTalk    = fPar->CrossTalk(time);
            }

            signalOld[0] = 0.0;
            signalOld[1] = 0.0;
            signalOld[2] = 0.0;

            for (iPad = 0; iPad < kNpad; iPad++) {

              Int_t colPos = colE + iPad - 1;
              if (colPos <        0) continue;
              if (colPos >= nColMax) break;

              // Add the signals
              // Note: The time bin number is shifted by nTimeBefore to avoid negative
              // time bins. This has to be subtracted later.
              Int_t iCurrentTimeBin = iTimeBin + nTimeBefore;
              signalOld[iPad]  = signals->GetDataUnchecked(rowE,colPos,iCurrentTimeBin);
              if( colPos != colE ) {
                signalOld[iPad] += padSignal[iPad] * (timeResponse + crossTalk);
              } 
              else {
                signalOld[iPad] += padSignal[iPad] * timeResponse;
              }
              signals->SetDataUnchecked(rowE,colPos,iCurrentTimeBin,signalOld[iPad]);

              // Store the track index in the dictionary
              // Note: We store index+1 in order to allow the array to be compressed
              if ((signalOld[iPad] > 0) && (!fSimpleSim)) { 
                for (iDict = 0; iDict < kNDict; iDict++) {
                  Int_t oldTrack = dictionary[iDict]->GetDataUnchecked(rowE
                                                                      ,colPos
                                                                      ,iCurrentTimeBin);
                  if (oldTrack == track+1) break;
                  if (oldTrack ==       0) {
                    dictionary[iDict]->SetDataUnchecked(rowE,colPos,iCurrentTimeBin,track+1);
                    break;
                  }
                }
              }

	    } // Loop: pads

	  } // Loop: time bins

        } // Loop: electrons of a single hit

      } // If: detector and test hit

      hit = (AliTRDhit *) fTRD->NextHit();   

    } // Loop: hits of one primary track

  } // Loop: primary tracks

  if (fDebug > 0) {
    printf("<AliTRDdigitizer::MakeDigits> ");
    printf("Finished analyzing %d hits\n",countHits);
  }

  // The coupling factor
  Float_t coupling = fPar->GetPadCoupling() 
                   * fPar->GetTimeCoupling();

  // The conversion factor
  Float_t convert  = kEl2fC
                   * fPar->GetChipGain();

  // Loop through all chambers to finalize the digits
  Int_t iDetBeg = 0;
  Int_t iDetEnd = AliTRDgeometry::Ndet();
  if (fSimpleSim) {
    iDetBeg = fSimpleDet;
    iDetEnd = iDetBeg + 1;
  }
  for (Int_t iDet = iDetBeg; iDet < iDetEnd; iDet++) {

    Int_t plane       = fGeo->GetPlane(iDet);
    Int_t sector      = fGeo->GetSector(iDet);
    Int_t chamber     = fGeo->GetChamber(iDet);
    Int_t nRowMax     = fPar->GetRowMax(plane,chamber,sector);
    Int_t nColMax     = fPar->GetColMax(plane);
    Int_t nTimeMax    = fPar->GetTimeMax();
    Int_t nTimeTotal  = fPar->GetTimeTotal();

    Double_t *inADC  = new Double_t[nTimeTotal];
    Double_t *outADC = new Double_t[nTimeTotal];

    if (fDebug > 0) {
      printf("<AliTRDdigitizer::MakeDigits> ");
      printf("Digitization for chamber %d\n",iDet);
    }

    // Add a container for the digits of this detector
    digits = fDigitsManager->GetDigits(iDet);        
    // Allocate memory space for the digits buffer
    if (digits->GetNtime() == 0) {
      digits->Allocate(nRowMax,nColMax,nTimeTotal);
    }
    else if (fSimpleSim) {
      digits->Clear();
    }
 
    // Get the signal container
    signals = (AliTRDdataArrayF *) signalsArray->At(iDet);
    if (signals->GetNtime() == 0) {
      // Create missing containers
      signals->Allocate(nRowMax,nColMax,nTimeTotal);      
    }
    else {
      // Expand the container if neccessary
      if (fCompress) signals->Expand();
    }
    // Create the missing dictionary containers
    if (!fSimpleSim) {    
      for (iDict = 0; iDict < kNDict; iDict++) {       
        dictionary[iDict] = fDigitsManager->GetDictionary(iDet,iDict);
        if (dictionary[iDict]->GetNtime() == 0) {
          dictionary[iDict]->Allocate(nRowMax,nColMax,nTimeTotal);
        }
      } 
    }

    Int_t nDigits = 0;

    // Don't create noise in detectors that are switched off
    if (CheckDetector(plane,chamber,sector)) {

      // Create the digits for this chamber
      for (iRow  = 0; iRow  <  nRowMax;   iRow++ ) {
        for (iCol  = 0; iCol  <  nColMax;   iCol++ ) {

	  // Create summable digits
          if (fSDigits) {

            for (iTime = 0; iTime < nTimeTotal; iTime++) {         
              Float_t signalAmp = signals->GetDataUnchecked(iRow,iCol,iTime);
              signalAmp *= fSDigitsScale;
              signalAmp  = TMath::Min(signalAmp,(Float_t) 1.0e9);
              Int_t adc  = (Int_t) signalAmp;
              if (adc > 0) nDigits++;
              digits->SetDataUnchecked(iRow,iCol,iTime,adc);
	    }

	  }
	  // Create normal digits
          else {

            for (iTime = 0; iTime < nTimeTotal; iTime++) {         
              Float_t signalAmp = signals->GetDataUnchecked(iRow,iCol,iTime);
              // Pad and time coupling
              signalAmp *= coupling;
              // Add the noise, starting from minus ADC baseline in electrons
              Double_t baselineEl = fPar->GetADCbaseline() * (fPar->GetADCinRange()
                                                           / fPar->GetADCoutRange()) 
                                                           / convert;
              signalAmp  = TMath::Max((Double_t) gRandom->Gaus(signalAmp,fPar->GetNoise())
                                     ,-baselineEl);
              // Convert to mV
              signalAmp *= convert;
              // Add ADC baseline in mV
              signalAmp += fPar->GetADCbaseline() * (fPar->GetADCinRange()
                                                   / fPar->GetADCoutRange());
 	      // Convert to ADC counts. Set the overflow-bit fADCoutRange if the 
	      // signal is larger than fADCinRange
              Int_t adc  = 0;
              if (signalAmp >= fPar->GetADCinRange()) {
                adc = ((Int_t) fPar->GetADCoutRange());
	      }
              else {
                adc = ((Int_t) (signalAmp * (fPar->GetADCoutRange() 
                                           / fPar->GetADCinRange())));
	      }
              inADC[iTime]  = adc;
              outADC[iTime] = adc;
	    }

	    // Apply the tail cancelation via the digital filter
            if (fPar->TCOn()) {
              DeConvExp(inADC,outADC,nTimeTotal,fPar->GetTCnexp());
	    }

            for (iTime = 0; iTime < nTimeTotal; iTime++) {   
              // Store the amplitude of the digit if above threshold
              if (outADC[iTime] > fPar->GetADCthreshold()) {
                if (fDebug > 2) {
                  printf("  iRow = %d, iCol = %d, iTime = %d, adc = %f\n"
                        ,iRow,iCol,iTime,outADC[iTime]);
	        }
                nDigits++;
                digits->SetDataUnchecked(iRow,iCol,iTime,((Int_t) outADC[iTime]));
  	      }
	    }

	  }

        }
      }

    }

    // Compress the arrays
    if (!fSimpleSim) {  
      digits->Compress(1,0);
      for (iDict = 0; iDict < kNDict; iDict++) {
        dictionary[iDict]->Compress(1,0);
      }

      totalSizeDigits += digits->GetSize();
      totalSizeDict0  += dictionary[0]->GetSize();
      totalSizeDict1  += dictionary[1]->GetSize();
      totalSizeDict2  += dictionary[2]->GetSize();

      Float_t nPixel = nRowMax * nColMax * nTimeMax;
      if (fDebug > 0) {
        printf("<AliTRDdigitizer::MakeDigits> ");
        printf("Found %d digits in detector %d (%3.0f).\n"
              ,nDigits,iDet
              ,100.0 * ((Float_t) nDigits) / nPixel);
      } 

      if (fCompress) signals->Compress(1,0);

    }

    delete [] inADC;
    delete [] outADC;

  }

  if (signalsArray) {
    delete signalsArray;
    signalsArray = 0;
  }

  if (fDebug > 0) {
    printf("<AliTRDdigitizer::MakeDigits> ");
    printf("Total number of analyzed hits = %d\n",countHits);
    if (!fSimpleSim) {    
      printf("<AliTRDdigitizer::MakeDigits> ");
      printf("Total digits data size = %d, %d, %d, %d\n",totalSizeDigits
                                                        ,totalSizeDict0
                                                        ,totalSizeDict1
                                                        ,totalSizeDict2);        
    }
  }

  return kTRUE;

}

//_____________________________________________________________________________
void AliTRDdigitizer::AddSDigitsManager(AliTRDdigitsManager *man)
{
  //
  // Add a digits manager for s-digits to the input list.
  //

  fSDigitsManagerList->Add(man);

}

//_____________________________________________________________________________
void AliTRDdigitizer::DeleteSDigitsManager()
{
  //
  // Removes digits manager from the input list.
  //

  fSDigitsManagerList->Delete();

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::ConvertSDigits()
{
  //
  // Converts s-digits to normal digits
  //

  // Number of track dictionary arrays
  const Int_t    kNDict = AliTRDdigitsManager::kNDict;

  // Converts number of electrons to fC
  const Double_t kEl2fC = 1.602E-19 * 1.0E15; 

  Int_t iDict = 0;
  Int_t iRow;
  Int_t iCol;
  Int_t iTime;

  if (!fPar) {    
    fPar = new AliTRDparameter("TRDparameter","Standard parameter");
    if (fDebug > 0) {
      printf("<AliTRDdigitizer::ConvertSDigits> ");
      printf("Create the default parameter object\n");
    }
  }

  Double_t sDigitsScale = 1.0 / GetSDigitsScale();
  Double_t noise        = fPar->GetNoise();
  Double_t padCoupling  = fPar->GetPadCoupling();
  Double_t timeCoupling = fPar->GetTimeCoupling();
  Double_t chipGain     = fPar->GetChipGain();
  Double_t coupling     = padCoupling * timeCoupling;
  Double_t convert      = kEl2fC * chipGain;
  Double_t adcInRange   = fPar->GetADCinRange();
  Double_t adcOutRange  = fPar->GetADCoutRange();
  Int_t    adcThreshold = fPar->GetADCthreshold();
  Int_t    adcBaseline  = fPar->GetADCbaseline();   

  AliTRDdataArrayI *digitsIn;
  AliTRDdataArrayI *digitsOut;
  AliTRDdataArrayI *dictionaryIn[kNDict];
  AliTRDdataArrayI *dictionaryOut[kNDict];

  // Loop through the detectors
  for (Int_t iDet = 0; iDet < AliTRDgeometry::Ndet(); iDet++) {

    if (fDebug > 0) {
      printf("<AliTRDdigitizer::ConvertSDigits> ");
      printf("Convert detector %d to digits.\n",iDet);
    }

    Int_t plane      = fGeo->GetPlane(iDet);
    Int_t sector     = fGeo->GetSector(iDet);
    Int_t chamber    = fGeo->GetChamber(iDet);
    Int_t nRowMax    = fPar->GetRowMax(plane,chamber,sector);
    Int_t nColMax    = fPar->GetColMax(plane);
    Int_t nTimeTotal = fPar->GetTimeTotal();

    Double_t *inADC  = new Double_t[nTimeTotal];
    Double_t *outADC = new Double_t[nTimeTotal];

    digitsIn  = fSDigitsManager->GetDigits(iDet);
    digitsIn->Expand();
    digitsOut = fDigitsManager->GetDigits(iDet);
    digitsOut->Allocate(nRowMax,nColMax,nTimeTotal);
    for (iDict = 0; iDict < kNDict; iDict++) {
      dictionaryIn[iDict]  = fSDigitsManager->GetDictionary(iDet,iDict);
      dictionaryIn[iDict]->Expand();
      dictionaryOut[iDict] = fDigitsManager->GetDictionary(iDet,iDict);
      dictionaryOut[iDict]->Allocate(nRowMax,nColMax,nTimeTotal);
    }

    for (iRow  = 0; iRow  <  nRowMax;   iRow++ ) {
      for (iCol  = 0; iCol  <  nColMax;   iCol++ ) {

        for (iTime = 0; iTime < nTimeTotal; iTime++) {         
          Double_t signal = (Double_t) digitsIn->GetDataUnchecked(iRow,iCol,iTime);
          signal *= sDigitsScale;
          // Pad and time coupling
          signal *= coupling;
          // Add the noise, starting from minus ADC baseline in electrons
          Double_t baselineEl = adcBaseline * (adcInRange / adcOutRange) / convert;
          signal  = TMath::Max((Double_t) gRandom->Gaus(signal,noise),-baselineEl);
          // Convert to mV
          signal *= convert;
          // add ADC baseline in mV
          signal += adcBaseline * (adcInRange / adcOutRange);
	  // Convert to ADC counts. Set the overflow-bit adcOutRange if the 
	  // signal is larger than adcInRange
          Int_t adc  = 0;
          if (signal >= adcInRange) {
            adc = ((Int_t) adcOutRange);
	  }
          else {
            adc = ((Int_t) (signal * (adcOutRange / adcInRange)));
	  }
          inADC[iTime]  = adc;
          outADC[iTime] = adc;
	}

	// Apply the tail cancelation via the digital filter
        if (fPar->TCOn()) {
          DeConvExp(inADC,outADC,nTimeTotal,fPar->GetTCnexp());
	}

        for (iTime = 0; iTime < nTimeTotal; iTime++) {   
          // Store the amplitude of the digit if above threshold
          if (outADC[iTime] > adcThreshold) {
            digitsOut->SetDataUnchecked(iRow,iCol,iTime,((Int_t) outADC[iTime]));
  	    // Copy the dictionary
            for (iDict = 0; iDict < kNDict; iDict++) { 
              Int_t track = dictionaryIn[iDict]->GetDataUnchecked(iRow,iCol,iTime);
              dictionaryOut[iDict]->SetDataUnchecked(iRow,iCol,iTime,track);
	    }
	  }
	}

      }
    }

    if (fCompress) {
      digitsIn->Compress(1,0);
      digitsOut->Compress(1,0);
      for (iDict = 0; iDict < kNDict; iDict++) {
        dictionaryIn[iDict]->Compress(1,0);
        dictionaryOut[iDict]->Compress(1,0);
      }
    }

    delete [] inADC;
    delete [] outADC;

  }    

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::MergeSDigits()
{
  //
  // Merges the input s-digits:
  //   - The amplitude of the different inputs are summed up.
  //   - Of the track IDs from the input dictionaries only one is
  //     kept for each input. This works for maximal 3 different merged inputs.
  //

  // Number of track dictionary arrays
  const Int_t kNDict = AliTRDdigitsManager::kNDict;

  if (!fPar) {
    fPar = new AliTRDparameter("TRDparameter","Standard parameter");
    if (fDebug > 0) {
      printf("<AliTRDdigitizer::MergeSDigits> ");
      printf("Create the default parameter object\n");
    }
  }

  Int_t iDict = 0;
  Int_t jDict = 0;

  AliTRDdataArrayI *digitsA;
  AliTRDdataArrayI *digitsB;
  AliTRDdataArrayI *dictionaryA[kNDict];
  AliTRDdataArrayI *dictionaryB[kNDict];

  // Get the first s-digits
  fSDigitsManager = (AliTRDdigitsManager *) fSDigitsManagerList->First();
  if (!fSDigitsManager) return kFALSE;

  // Loop through the other sets of s-digits
  AliTRDdigitsManager *mergeSDigitsManager;
  mergeSDigitsManager = (AliTRDdigitsManager *) 
                        fSDigitsManagerList->After(fSDigitsManager);

  if (fDebug > 0) {
    if (mergeSDigitsManager) {
      printf("<AliTRDdigitizer::MergeSDigits> ");
      printf("Merge %d input files.\n",fSDigitsManagerList->GetSize());
    }
    else {
      printf("<AliTRDdigitizer::MergeSDigits> ");
      printf("Only one input file.\n");
    }
  }

  Int_t iMerge = 0;
  while (mergeSDigitsManager) {

    iMerge++;

    // Loop through the detectors
    for (Int_t iDet = 0; iDet < AliTRDgeometry::Ndet(); iDet++) {

      Int_t plane      = fGeo->GetPlane(iDet);
      Int_t sector     = fGeo->GetSector(iDet);
      Int_t chamber    = fGeo->GetChamber(iDet);
      Int_t nRowMax    = fPar->GetRowMax(plane,chamber,sector);
      Int_t nColMax    = fPar->GetColMax(plane);
      Int_t nTimeTotal = fPar->GetTimeTotal();

      // Loop through the pixels of one detector and add the signals
      digitsA = fSDigitsManager->GetDigits(iDet);
      digitsB = mergeSDigitsManager->GetDigits(iDet);
      digitsA->Expand();
      digitsB->Expand();
      for (iDict = 0; iDict < kNDict; iDict++) {
        dictionaryA[iDict] = fSDigitsManager->GetDictionary(iDet,iDict);
        dictionaryB[iDict] = mergeSDigitsManager->GetDictionary(iDet,iDict);
        dictionaryA[iDict]->Expand();
        dictionaryB[iDict]->Expand();
      }

      // Merge only detectors that contain a signal
      Bool_t doMerge = kTRUE;
      if (fMergeSignalOnly) {
        if (digitsA->GetOverThreshold(0) == 0) {
          doMerge = kFALSE;
	}
      }

      if (doMerge) {

        if (fDebug > 0) {
          printf("<AliTRDdigitizer::MergeSDigits> ");
          printf("Merge detector %d of input no.%d\n",iDet,iMerge+1);
        }

        for (Int_t iRow  = 0; iRow  <  nRowMax;   iRow++ ) {
          for (Int_t iCol  = 0; iCol  <  nColMax;   iCol++ ) {
            for (Int_t iTime = 0; iTime < nTimeTotal; iTime++) {         

	      // Add the amplitudes of the summable digits 
              Int_t ampA = digitsA->GetDataUnchecked(iRow,iCol,iTime);
              Int_t ampB = digitsB->GetDataUnchecked(iRow,iCol,iTime);
              ampA += ampB;
              digitsA->SetDataUnchecked(iRow,iCol,iTime,ampA);

	     // Add the mask to the track id if defined.
              for (iDict = 0; iDict < kNDict; iDict++) {
                Int_t trackB = dictionaryB[iDict]->GetDataUnchecked(iRow,iCol,iTime);
                if ((fMasks) && (trackB > 0)) {
                  for (jDict = 0; jDict < kNDict; jDict++) { 
                    Int_t trackA = dictionaryA[iDict]->GetDataUnchecked(iRow,iCol,iTime);
                    if (trackA == 0) {
                      trackA = trackB + fMasks[iMerge];
                      dictionaryA[iDict]->SetDataUnchecked(iRow,iCol,iTime,trackA);
		    }
	  	  }
	        }
	      }

	    }
	  }
        }

      }

      if (fCompress) {
        digitsA->Compress(1,0);
        digitsB->Compress(1,0);
        for (iDict = 0; iDict < kNDict; iDict++) {
          dictionaryA[iDict]->Compress(1,0);
          dictionaryB[iDict]->Compress(1,0);
        }
      }

    }    

    // The next set of s-digits
    mergeSDigitsManager = (AliTRDdigitsManager *) 
                          fSDigitsManagerList->After(mergeSDigitsManager);

  }

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::SDigits2Digits()
{
  //
  // Merges the input s-digits and converts them to normal digits
  //

  if (!MergeSDigits()) return kFALSE;

  return ConvertSDigits();

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::CheckDetector(Int_t plane, Int_t chamber, Int_t sector)
{
  //
  // Checks whether a detector is enabled
  //

  if (fSimpleSim) return kTRUE; 

  if ((fTRD->GetSensChamber() >=       0) &&
      (fTRD->GetSensChamber() != chamber)) return kFALSE;
  if ((fTRD->GetSensPlane()   >=       0) &&
      (fTRD->GetSensPlane()   !=   plane)) return kFALSE;
  if ( fTRD->GetSensSector()  >=       0) {
    Int_t sens1 = fTRD->GetSensSector();
    Int_t sens2 = sens1 + fTRD->GetSensSectorRange();
    sens2 -= ((Int_t) (sens2 / AliTRDgeometry::Nsect())) 
           * AliTRDgeometry::Nsect();
    if (sens1 < sens2) {
      if ((sector < sens1) || (sector >= sens2)) return kFALSE;
    }
    else {
      if ((sector < sens1) && (sector >= sens2)) return kFALSE;
    }
  }

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::WriteDigits() const
{
  //
  // Writes out the TRD-digits and the dictionaries
  //

  //Write parameters
  fRunLoader->CdGAFile();
  if (!gFile->Get("TRDparameter")) GetParameter()->Write();

  // Store the digits and the dictionary in the tree
  return fDigitsManager->WriteDigits();

}

//_____________________________________________________________________________
void AliTRDdigitizer::DeConvExp(Double_t *source, Double_t *target
                              , Int_t n, Int_t nexp) 
{
  //
  // Does the deconvolution by the digital filter.
  //
  // Author:        Marcus Gutfleisch, KIP Heidelberg
  // Optimized for: New TRF from Venelin Angelov, simulated with CADENCE
  //                Pad-ground capacitance = 25 pF
  //                Pad-pad cross talk capacitance = 6 pF
  //                For 10 MHz digitization, corresponding to 20 time bins
  //                in the drift region
  //

  Double_t rates[2];
  Double_t coefficients[2];

  /* initialize (coefficient = alpha, rates = lambda) */
  
  if( nexp == 1 ) {
    rates[0] = 0.466998;
    /* no rescaling */
    coefficients[0] = 1.0;
  }
  if( nexp == 2 ) {
    rates[0] = 0.8988162;
    coefficients[0] = 0.11392069;
    rates[1] = 0.3745688;
    coefficients[1] = 0.8860793;
    /* no rescaling */
    Float_t sumc = coefficients[0]+coefficients[1];
    coefficients[0] /= sumc;
    coefficients[1] /= sumc;
  }
      
  Int_t i, k;
  Double_t reminder[2];
  Double_t correction, result;

  /* attention: computation order is important */
  correction=0.0;
  for ( k = 0; k < nexp; k++ ) reminder[k]=0.0;
    
  for ( i = 0; i < n; i++ ) {
    result = ( source[i] - correction );    /* no rescaling */
    target[i] = result;
    
    for ( k = 0; k < nexp; k++ ) reminder[k] = rates[k] 
                             * ( reminder[k] + coefficients[k] * result);
      
    correction=0.0;
    for ( k = 0; k < nexp; k++ ) correction += reminder[k];
  }
  
}

//_____________________________________________________________________________
void AliTRDdigitizer::InitOutput(Int_t iEvent)
{
  //
  // Initializes the output branches
  //

  fEvent = iEvent;
   
  if (!fRunLoader)
   {
     Error("InitOutput","Run Loader is NULL");
     return;  
   }
  AliLoader* loader = fRunLoader->GetLoader("TRDLoader");
  if (!loader)
   {
     Error("Open","Can not get TRD loader from Run Loader");
     return;
   }

  TTree* tree = 0;
  
  if (fSDigits)
   { 
   //if we produce SDigits
    tree = loader->TreeS();
    if (!tree)
     {
      loader->MakeTree("S");
      tree = loader->TreeS();
     }
   }
  else
   {//if we produce Digits
     tree = loader->TreeD();
     if (!tree)
      {
       loader->MakeTree("D");
       tree = loader->TreeD();
      }
   }
  fDigitsManager->SetEvent(iEvent);
  fDigitsManager->MakeBranch(tree);

}
