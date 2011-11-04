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

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Creates and handles digits from TRD hits                              //
//                                                                        //
//  Authors: C. Blume (blume@ikf.uni-frankfurt.de)                        //
//           C. Lippmann                                                  //
//           B. Vulpescu                                                  //
//                                                                        //
//  The following effects are included:                                   //
//      - Diffusion                                                       //
//      - ExB effects                                                     //
//      - Gas gain including fluctuations                                 //
//      - Pad-response (simple Gaussian approximation)                    //
//      - Time-response                                                   //
//      - Electronics noise                                               //
//      - Electronics gain                                                //
//      - Digitization                                                    //
//      - Zero suppression                                                //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TGeoManager.h>
#include <TList.h>
#include <TMath.h>
#include <TRandom.h>
#include <TTree.h>

#include "AliRun.h"
#include "AliMC.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliConfig.h"
#include "AliDigitizationInput.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliLog.h"

#include "AliTRD.h"
#include "AliTRDhit.h"
#include "AliTRDdigitizer.h"
#include "AliTRDarrayDictionary.h"
#include "AliTRDarrayADC.h"
#include "AliTRDarraySignal.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"
#include "AliTRDcalibDB.h"
#include "AliTRDSimParam.h"
#include "AliTRDCommonParam.h"
#include "AliTRDfeeParam.h"
#include "AliTRDmcmSim.h"
#include "AliTRDdigitsParam.h"

#include "Cal/AliTRDCalROC.h"
#include "Cal/AliTRDCalDet.h"
#include "Cal/AliTRDCalOnlineGainTableROC.h"

ClassImp(AliTRDdigitizer)

//_____________________________________________________________________________
AliTRDdigitizer::AliTRDdigitizer()
  :AliDigitizer()
  ,fRunLoader(0)
  ,fDigitsManager(0)
  ,fSDigitsManager(0)
  ,fSDigitsManagerList(0)
  ,fTRD(0)
  ,fGeo(0)
  ,fMcmSim(new AliTRDmcmSim)
  ,fEvent(0)
  ,fMasks(0)
  ,fCompress(kTRUE)
  ,fSDigits(kFALSE)
  ,fMergeSignalOnly(kFALSE)
{
  //
  // AliTRDdigitizer default constructor
  //
  
}

//_____________________________________________________________________________
AliTRDdigitizer::AliTRDdigitizer(const Text_t *name, const Text_t *title)              
  :AliDigitizer(name,title)
  ,fRunLoader(0)
  ,fDigitsManager(0)
  ,fSDigitsManager(0)
  ,fSDigitsManagerList(0)
  ,fTRD(0)
  ,fGeo(0)
  ,fMcmSim(new AliTRDmcmSim)
  ,fEvent(0)
  ,fMasks(0)
  ,fCompress(kTRUE)
  ,fSDigits(kFALSE)
  ,fMergeSignalOnly(kFALSE)
{
  //
  // AliTRDdigitizer constructor
  //

}

//_____________________________________________________________________________
AliTRDdigitizer::AliTRDdigitizer(AliDigitizationInput* digInput
                               , const Text_t *name, const Text_t *title)
  :AliDigitizer(digInput,name,title)
  ,fRunLoader(0)
  ,fDigitsManager(0)
  ,fSDigitsManager(0)
  ,fSDigitsManagerList(0)
  ,fTRD(0)
  ,fGeo(0)
  ,fMcmSim(new AliTRDmcmSim)
  ,fEvent(0)
  ,fMasks(0)
  ,fCompress(kTRUE)
  ,fSDigits(kFALSE)
  ,fMergeSignalOnly(kFALSE)
{
  //
  // AliTRDdigitizer constructor
  //

}

//_____________________________________________________________________________
AliTRDdigitizer::AliTRDdigitizer(AliDigitizationInput* digInput)
  :AliDigitizer(digInput,"AliTRDdigitizer","TRD digitizer")
  ,fRunLoader(0)
  ,fDigitsManager(0)
  ,fSDigitsManager(0)
  ,fSDigitsManagerList(0)
  ,fTRD(0)
  ,fGeo(0)
  ,fMcmSim(new AliTRDmcmSim)
  ,fEvent(0)
  ,fMasks(0)
  ,fCompress(kTRUE)
  ,fSDigits(kFALSE)
  ,fMergeSignalOnly(kFALSE)
{
  //
  // AliTRDdigitizer constructor
  //

}

//_____________________________________________________________________________
AliTRDdigitizer::AliTRDdigitizer(const AliTRDdigitizer &d)
  :AliDigitizer(d)
  ,fRunLoader(0)
  ,fDigitsManager(0)
  ,fSDigitsManager(0)
  ,fSDigitsManagerList(0)
  ,fTRD(0)
  ,fGeo(0)
  ,fMcmSim(new AliTRDmcmSim)
  ,fEvent(0)
  ,fMasks(0)
  ,fCompress(d.fCompress)
  ,fSDigits(d.fSDigits)
  ,fMergeSignalOnly(d.fMergeSignalOnly)
{
  //
  // AliTRDdigitizer copy constructor
  //

}

//_____________________________________________________________________________
AliTRDdigitizer::~AliTRDdigitizer()
{
  //
  // AliTRDdigitizer destructor
  //

  if (fDigitsManager) {
    delete fDigitsManager;
    fDigitsManager      = 0;
  }

  if (fSDigitsManager) {
    // s-digitsmanager will be deleted via list
    fSDigitsManager     = 0;
  }

  if (fSDigitsManagerList) {
    fSDigitsManagerList->Delete();
    delete fSDigitsManagerList;
    fSDigitsManagerList = 0;
  }

  if (fMasks) {
    delete [] fMasks;
    fMasks       = 0;
  }

  if (fMcmSim) {
    delete fMcmSim;
    fMcmSim = 0;
  }

  if (fGeo) {
    delete fGeo;
    fGeo = 0;
  }

}

//_____________________________________________________________________________
AliTRDdigitizer &AliTRDdigitizer::operator=(const AliTRDdigitizer &d)
{
  //
  // Assignment operator
  //

  if (this != &d) {
    ((AliTRDdigitizer &) d).Copy(*this);
  }

  return *this;

}

//_____________________________________________________________________________
void AliTRDdigitizer::Copy(TObject &d) const
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
  ((AliTRDdigitizer &) d).fEvent              = 0;
  ((AliTRDdigitizer &) d).fMasks              = 0;
  ((AliTRDdigitizer &) d).fCompress           = fCompress;
  ((AliTRDdigitizer &) d).fSDigits            = fSDigits;
  ((AliTRDdigitizer &) d).fMergeSignalOnly    = fMergeSignalOnly;

}

//_____________________________________________________________________________
void AliTRDdigitizer::Digitize(const Option_t* option)
{
  //
  // Executes the merging
  //

  Int_t iInput;

  AliTRDdigitsManager *sdigitsManager;

  TString optionString = option;
  if (optionString.Contains("deb")) {
    AliLog::SetClassDebugLevel("AliTRDdigitizer",1);
    AliInfo("Called with debug option");
  }

  // The AliRoot file is already connected by the manager
  AliRunLoader *inrl = 0x0;
  
  if (gAlice) {
    AliDebug(1,"AliRun object found on file.");
  }
  else {
    inrl = AliRunLoader::GetRunLoader(fDigInput->GetInputFolderName(0));
    inrl->LoadgAlice();
    gAlice = inrl->GetAliRun();
    if (!gAlice) {
      AliError("Could not find AliRun object.");
      return;
    }
  }
                                                                           
  Int_t nInput = fDigInput->GetNinputs();
  fMasks       = new Int_t[nInput];
  for (iInput = 0; iInput < nInput; iInput++) {
    fMasks[iInput] = fDigInput->GetMask(iInput);
  }

  //
  // Initialization
  //

  AliRunLoader *orl = AliRunLoader::GetRunLoader(fDigInput->GetOutputFolderName());

  if (InitDetector()) {

    AliLoader *ogime = orl->GetLoader("TRDLoader");

    TTree *tree = 0;
    if (fSDigits) { 
      // If we produce SDigits
      tree = ogime->TreeS();
      if (!tree) {
	ogime->MakeTree("S");
	tree = ogime->TreeS();
      }
    }
    else {
      // If we produce Digits
      tree = ogime->TreeD();
      if (!tree) {
	ogime->MakeTree("D");
        tree = ogime->TreeD();
      }
    }

    MakeBranch(tree);

  }
 
  for (iInput = 0; iInput < nInput; iInput++) {

    AliDebug(1,Form("Add input stream %d",iInput));

    // Check if the input tree exists
    inrl = AliRunLoader::GetRunLoader(fDigInput->GetInputFolderName(iInput));
    AliLoader *gime = inrl->GetLoader("TRDLoader");

    TTree *treees = gime->TreeS();
    if (treees == 0x0) {
      if (gime->LoadSDigits()) {
        AliError(Form("Error Occured while loading S. Digits for input %d.",iInput));
        return;
      }
      treees = gime->TreeS();
    }
    
    if (treees == 0x0) {
      AliError(Form("Input stream %d does not exist",iInput));
      return;
    } 

    // Read the s-digits via digits manager
    sdigitsManager = new AliTRDdigitsManager();
    sdigitsManager->SetSDigits(kTRUE);
    
    AliRunLoader *rl = AliRunLoader::GetRunLoader(fDigInput->GetInputFolderName(iInput));
    AliLoader *gimme = rl->GetLoader("TRDLoader");
    if (!gimme->TreeS()) 
      {
	gimme->LoadSDigits();
      }

    sdigitsManager->ReadDigits(gimme->TreeS());
   
    // Add the s-digits to the input list 
    AddSDigitsManager(sdigitsManager);

  }

  // Convert the s-digits to normal digits
  AliDebug(1,"Do the conversion");
  SDigits2Digits();

  // Store the digits
  AliDebug(1,"Write the digits");
  fDigitsManager->WriteDigits();

  // Write parameters
  orl->CdGAFile();

  // Clean up
  DeleteSDigitsManager();

  AliDebug(1,"Done");

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::Open(const Char_t *file, Int_t nEvent)
{
  //
  // Opens a ROOT-file with TRD-hits and reads in the hit-tree
  //
  // Connect the AliRoot file containing Geometry, Kine, and Hits
  //  

  TString evfoldname = AliConfig::GetDefaultEventFolderName();

  fRunLoader = AliRunLoader::GetRunLoader(evfoldname);
  if (!fRunLoader) {
    fRunLoader = AliRunLoader::Open(file,evfoldname,"UPDATE");
  }  
  if (!fRunLoader) {
    AliError(Form("Can not open session for file %s.",file));
    return kFALSE;
  }
   
  if (!fRunLoader->GetAliRun()) {
    fRunLoader->LoadgAlice();
  }
  gAlice = fRunLoader->GetAliRun();
  
  if (gAlice) {
    AliDebug(1,"AliRun object found on file.");
  }
  else {
    AliError("Could not find AliRun object.");
    return kFALSE;
  }

  fEvent = nEvent;

  AliLoader *loader = fRunLoader->GetLoader("TRDLoader");
  if (!loader) {
    AliError("Can not get TRD loader from Run Loader");
    return kFALSE;
  }
  
  if (InitDetector()) {
    TTree *tree = 0;
    if (fSDigits) { 
      // If we produce SDigits
      tree = loader->TreeS();
      if (!tree) {
        loader->MakeTree("S");
        tree = loader->TreeS();
      }
    }
    else {
      // If we produce Digits
      tree = loader->TreeD();
      if (!tree) {
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
Bool_t AliTRDdigitizer::Open(AliRunLoader * const runLoader, Int_t nEvent)
{
  //
  // Opens a ROOT-file with TRD-hits and reads in the hit-tree
  //
  // Connect the AliRoot file containing Geometry, Kine, and Hits
  //  

  fRunLoader = runLoader;
  if (!fRunLoader) {
    AliError("RunLoader does not exist");
    return kFALSE;
  }
   
  if (!fRunLoader->GetAliRun()) {
    fRunLoader->LoadgAlice();
  }
  gAlice = fRunLoader->GetAliRun();
  
  if (gAlice) {
    AliDebug(1,"AliRun object found on file.");
  }
  else {
    AliError("Could not find AliRun object.");
    return kFALSE;
  }

  fEvent = nEvent;

  AliLoader *loader = fRunLoader->GetLoader("TRDLoader");
  if (!loader) {
    AliError("Can not get TRD loader from Run Loader");
    return kFALSE;
  }
  
  if (InitDetector()) {
    TTree *tree = 0;
    if (fSDigits) { 
      // If we produce SDigits
      tree = loader->TreeS();
      if (!tree) {
        loader->MakeTree("S");
        tree = loader->TreeS();
      }
    }
    else {
      // If we produce Digits
      tree = loader->TreeD();
      if (!tree) {
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
    AliFatal("No TRD module found");
    exit(1);
  }
  if (fTRD->IsVersion() != 1) {
    AliFatal("TRD must be version 1 (slow simulator)");
    exit(1);
  }

  // Get the geometry
  fGeo = new AliTRDgeometry();

  // Create a digits manager
  if (fDigitsManager) {
    delete fDigitsManager;
  }
  fDigitsManager = new AliTRDdigitsManager();
  fDigitsManager->SetSDigits(fSDigits);
  fDigitsManager->CreateArrays();
  fDigitsManager->SetEvent(fEvent);

  // The list for the input s-digits manager to be merged
  if (fSDigitsManagerList) {
    fSDigitsManagerList->Delete();
  } 
  else {
    fSDigitsManagerList = new TList();
  }

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::MakeBranch(TTree *tree) const
{
  // 
  // Create the branches for the digits array
  //

  return fDigitsManager->MakeBranch(tree);

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
Bool_t AliTRDdigitizer::MakeDigits()
{
  //
  // Creates digits.
  //

  AliDebug(1,"Start creating digits");

  if (!fGeo) {
    AliError("No geometry defined");
    return kFALSE;
  }

  AliTRDcalibDB *calibration = AliTRDcalibDB::Instance();
  if (!calibration) {
    AliFatal("Could not get calibration object");
    return kFALSE;
  }

  const Int_t kNdet  = AliTRDgeometry::Ndet();

  Float_t **hits = new Float_t*[kNdet];
  Int_t    *nhit = new Int_t[kNdet];

  AliTRDarraySignal *signals = 0x0;

  // Check the number of time bins from simParam against OCDB,
  // if OCDB value is not supposed to be used.
  // As default, the value from OCDB is taken
  if (AliTRDSimParam::Instance()->GetNTBoverwriteOCDB()) {
    if (calibration->GetNumberOfTimeBinsDCS() != AliTRDSimParam::Instance()->GetNTimeBins()) {
      AliWarning(Form("Number of time bins is different to OCDB value [SIM=%d, OCDB=%d]"
                     ,AliTRDSimParam::Instance()->GetNTimeBins()
                     ,calibration->GetNumberOfTimeBinsDCS()));
    }
    // Save the values for the raw data headers
    fDigitsManager->GetDigitsParam()->SetNTimeBinsAll(AliTRDSimParam::Instance()->GetNTimeBins());
  }
  else {
    // Save the values for the raw data headers
    fDigitsManager->GetDigitsParam()->SetNTimeBinsAll(calibration->GetNumberOfTimeBinsDCS());
  }

  // Save the values for the raw data headers
  fDigitsManager->GetDigitsParam()->SetADCbaselineAll(AliTRDSimParam::Instance()->GetADCbaseline());
 
  // Sort all hits according to detector number
  if (!SortHits(hits,nhit)) {
    AliError("Sorting hits failed");
    delete [] hits;
    delete [] nhit;
    return kFALSE;
  }

  // Loop through all detectors
  for (Int_t det = 0; det < kNdet; det++) {

    // Detectors that are switched off, not installed, etc.
    if ((!calibration->IsChamberNoData(det))    &&
        ( fGeo->ChamberInGeometry(det))         &&
        (nhit[det] > 0)) {

      signals = new AliTRDarraySignal();
	  
      // Convert the hits of the current detector to detector signals
      if (!ConvertHits(det,hits[det],nhit[det],signals)) {
	AliError(Form("Conversion of hits failed for detector=%d",det));
        delete [] hits;
        delete [] nhit;
        delete signals;
        signals = 0x0;
        return kFALSE;
      }

      // Convert the detector signals to digits or s-digits
      if (!ConvertSignals(det,signals)) {
	AliError(Form("Conversion of signals failed for detector=%d",det));
        delete [] hits;
        delete [] nhit;
        delete signals;
        signals = 0x0;
	return kFALSE;
      }

      // Delete the signals array
      delete signals;
      signals = 0x0;

    } // if: detector status

    delete [] hits[det];

  } // for: detector

  if (!fSDigits) {
    if (AliDataLoader *trklLoader 
          = AliRunLoader::Instance()->GetLoader("TRDLoader")->GetDataLoader("tracklets")) {
      if (trklLoader->Tree())
        trklLoader->WriteData("OVERWRITE");
    }
  }
    
  delete [] hits;
  delete [] nhit;

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::SortHits(Float_t **hits, Int_t *nhit)
{
  //
  // Read all the hits and sorts them according to detector number
  // in the output array <hits>.
  //

  AliDebug(1,"Start sorting hits");

  const Int_t kNdet = AliTRDgeometry::Ndet();
  // Size of the hit vector
  const Int_t kNhit = 6;

  Float_t *xyz      = 0;
  Int_t    nhitTrk  = 0;

  Int_t   *lhit     = new Int_t[kNdet];

  for (Int_t det = 0; det < kNdet; det++) {
    lhit[det] = 0;
    nhit[det] = 0;
    hits[det] = 0;
  }

  AliLoader *gimme = fRunLoader->GetLoader("TRDLoader");
  if (!gimme->TreeH()) {
    gimme->LoadHits();
  }
  TTree *hitTree = gimme->TreeH();
  if (hitTree == 0x0) {
    AliError("Can not get TreeH");
    delete [] lhit;
    return kFALSE;
  }
  fTRD->SetTreeAddress();

  // Get the number of entries in the hit tree
  // (Number of primary particles creating a hit somewhere)
  Int_t nTrk = (Int_t) hitTree->GetEntries();
  AliDebug(1,Form("Found %d tracks",nTrk));

  // Loop through all the tracks in the tree
  for (Int_t iTrk = 0; iTrk < nTrk; iTrk++) {

    gAlice->GetMCApp()->ResetHits();
    hitTree->GetEvent(iTrk);

    if (!fTRD->Hits()) {
      AliError(Form("No hits array for track = %d",iTrk));
      continue;
    }

    // Number of hits for this track
    nhitTrk = fTRD->Hits()->GetEntriesFast();

    Int_t hitCnt = 0;
    // Loop through the TRD hits
    AliTRDhit *hit = (AliTRDhit *) fTRD->FirstHit(-1);
    while (hit) {

      hitCnt++;

      // Don't analyze test hits
      if (((Int_t) hit->GetCharge()) != 0) {

        Int_t   trk  = hit->Track();
        Int_t   det  = hit->GetDetector();
        Int_t   q    = hit->GetCharge();
        Float_t x    = hit->X();
        Float_t y    = hit->Y();
        Float_t z    = hit->Z();
        Float_t time = hit->GetTime();

        if (nhit[det] == lhit[det]) {
          // Inititialization of new detector
          xyz        = new Float_t[kNhit*(nhitTrk+lhit[det])];
          if (hits[det]) {
            memcpy(xyz,hits[det],sizeof(Float_t)*kNhit*lhit[det]);
            delete [] hits[det];
	  }
          lhit[det] += nhitTrk;
          hits[det]  = xyz;
	}
        else {
          xyz        = hits[det];
	}
        xyz[nhit[det]*kNhit+0] = x;
        xyz[nhit[det]*kNhit+1] = y;
        xyz[nhit[det]*kNhit+2] = z;
        xyz[nhit[det]*kNhit+3] = q;
        xyz[nhit[det]*kNhit+4] = trk;
        xyz[nhit[det]*kNhit+5] = time;
        nhit[det]++;

      } // if: charge != 0

      hit = (AliTRDhit *) fTRD->NextHit();   

    } // for: hits of one track

  } // for: tracks

  delete [] lhit;

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::ConvertHits(Int_t det
                                  , const Float_t * const hits
                                  , Int_t nhit
                                  , AliTRDarraySignal *signals)
{
  //
  // Converts the detectorwise sorted hits to detector signals
  //

  AliDebug(1,Form("Start converting hits for detector=%d (nhits=%d)",det,nhit));

  // Number of pads included in the pad response
  const Int_t kNpad      = 3;
  // Number of track dictionary arrays
  const Int_t kNdict     = AliTRDdigitsManager::kNDict;
  // Size of the hit vector
  const Int_t kNhit      = 6;

  // Width of the amplification region
  const Float_t kAmWidth = AliTRDgeometry::AmThick();
  // Width of the drift region
  const Float_t kDrWidth = AliTRDgeometry::DrThick();
  // Drift + amplification region 
  const Float_t kDrMin   =          - 0.5 * kAmWidth;
  const Float_t kDrMax   = kDrWidth + 0.5 * kAmWidth;
  
  Int_t    iPad          = 0;
  Int_t    dict          = 0;
  Int_t    timeBinTRFend = 1;

  Double_t pos[3];
  Double_t loc[3];
  Double_t padSignal[kNpad];
  Double_t signalOld[kNpad];

  AliTRDarrayDictionary *dictionary[kNdict];

  AliTRDSimParam    *simParam    = AliTRDSimParam::Instance();
  AliTRDCommonParam *commonParam = AliTRDCommonParam::Instance();
  AliTRDcalibDB     *calibration = AliTRDcalibDB::Instance();

  if (!commonParam) {
    AliFatal("Could not get common parameterss");
    return kFALSE;
  }
  if (!simParam) {
    AliFatal("Could not get simulation parameters");
    return kFALSE;
  }  
  if (!calibration) {
    AliFatal("Could not get calibration object");  
    return kFALSE;
  }

  // Get the detector wise calibration objects
  AliTRDCalROC       *calVdriftROC      = 0;
  Float_t             calVdriftDetValue = 0.0;
  const AliTRDCalDet *calVdriftDet      = calibration->GetVdriftDet();  
  AliTRDCalROC       *calT0ROC          = 0;
  Float_t             calT0DetValue     = 0.0;
  const AliTRDCalDet *calT0Det          = calibration->GetT0Det();  
  Double_t            calExBDetValue    = 0.0;
  const AliTRDCalDet *calExBDet         = calibration->GetExBDet();

  if (simParam->TRFOn()) {
    timeBinTRFend = ((Int_t) (simParam->GetTRFhi() 
                  * commonParam->GetSamplingFrequency())) - 1;
  }

  Int_t   nTimeTotal   = fDigitsManager->GetDigitsParam()->GetNTimeBins(det);
  Float_t samplingRate = commonParam->GetSamplingFrequency();
  Float_t elAttachProp = simParam->GetElAttachProp() / 100.0; 

  AliTRDpadPlane *padPlane = fGeo->GetPadPlane(det);
  Int_t   layer   = fGeo->GetLayer(det);         //update
  Float_t row0    = padPlane->GetRow0ROC();
  Int_t   nRowMax = padPlane->GetNrows();
  Int_t   nColMax = padPlane->GetNcols();

  // Create a new array for the signals
  signals->Allocate(nRowMax,nColMax,nTimeTotal);

  // Create a new array for the dictionary
  for (dict = 0; dict < kNdict; dict++) {       
    dictionary[dict] = (AliTRDarrayDictionary *) fDigitsManager->GetDictionary(det,dict);
    dictionary[dict]->Allocate(nRowMax,nColMax,nTimeTotal);
  }      

  // Loop through the hits in this detector
  for (Int_t hit = 0; hit < nhit; hit++) {

    pos[0]          = hits[hit*kNhit+0];
    pos[1]          = hits[hit*kNhit+1];
    pos[2]          = hits[hit*kNhit+2];
    Float_t q       = hits[hit*kNhit+3];
    Float_t hittime = hits[hit*kNhit+5];
    Int_t   track   = ((Int_t) hits[hit*kNhit+4]);

    Int_t   inDrift = 1;

    // Find the current volume with the geo manager
    gGeoManager->SetCurrentPoint(pos);
    gGeoManager->FindNode();      
    if (strstr(gGeoManager->GetPath(),"/UK")) {
      inDrift = 0;
    }

    // Get the calibration objects
    calVdriftROC      = calibration->GetVdriftROC(det);
    calVdriftDetValue = calVdriftDet->GetValue(det);
    calT0ROC          = calibration->GetT0ROC(det);
    calT0DetValue     = calT0Det->GetValue(det);
    calExBDetValue    = calExBDet->GetValue(det);

    // Go to the local coordinate system:
    // loc[0] - col  direction in amplification or driftvolume
    // loc[1] - row  direction in amplification or driftvolume
    // loc[2] - time direction in amplification or driftvolume
    gGeoManager->MasterToLocal(pos,loc);
    if (inDrift) {
      // Relative to middle of amplification region
      loc[2] = loc[2] - kDrWidth/2.0 - kAmWidth/2.0;
    } 

    // The driftlength [cm] (w/o diffusion yet !).
    // It is negative if the hit is between pad plane and anode wires.
    Double_t driftlength = -1.0 * loc[2];

    // Stupid patch to take care of TR photons that are absorbed
    // outside the chamber volume. A real fix would actually need
    // a more clever implementation of the TR hit generation
    if (q < 0.0) {
      if ((loc[1] < padPlane->GetRowEndROC()) ||
          (loc[1] > padPlane->GetRow0ROC())) {
        continue;
      }
      if ((driftlength < kDrMin) ||
          (driftlength > kDrMax)) {
        continue;
      }
    }

    // Get row and col of unsmeared electron to retrieve drift velocity
    // The pad row (z-direction)
    Int_t    rowE         = padPlane->GetPadRowNumberROC(loc[1]);
    if (rowE < 0) {
      continue;
    }
    Double_t rowOffset    = padPlane->GetPadRowOffsetROC(rowE,loc[1]);
    // The pad column (rphi-direction)
    Double_t offsetTilt   = padPlane->GetTiltOffset(rowOffset);
    Int_t    colE         = padPlane->GetPadColNumber(loc[0]+offsetTilt);
    if (colE < 0) {
      continue;	  
    }
    Double_t colOffset    = 0.0;

    // Normalized drift length
    Float_t  driftvelocity  = calVdriftDetValue * calVdriftROC->GetValue(colE,rowE);
    Double_t absdriftlength = TMath::Abs(driftlength);
    if (commonParam->ExBOn()) {
      absdriftlength /= TMath::Sqrt(1.0 / (1.0 + calExBDetValue*calExBDetValue));
    }

    // Loop over all electrons of this hit
    // TR photons produce hits with negative charge
    Int_t nEl = ((Int_t) TMath::Abs(q));
    for (Int_t iEl = 0; iEl < nEl; iEl++) {

      // Now the real local coordinate system of the ROC
      // column direction: locC
      // row direction:    locR 
      // time direction:   locT
      // locR and locC are identical to the coordinates of the corresponding
      // volumina of the drift or amplification region.
      // locT is defined relative to the wire plane (i.e. middle of amplification
      // region), meaning locT = 0, and is negative for hits coming from the
      // drift region. 
      Double_t locC = loc[0];
      Double_t locR = loc[1];
      Double_t locT = loc[2];

      // Electron attachment
      if (simParam->ElAttachOn()) {
        if (gRandom->Rndm() < (absdriftlength * elAttachProp)) {
          continue;
        }
      }
          
      // Apply the diffusion smearing
      if (simParam->DiffusionOn()) {
        if (!(Diffusion(driftvelocity,absdriftlength,calExBDetValue,locR,locC,locT))) {
          continue;
	}
      }

      // Apply E x B effects (depends on drift direction)
      if (commonParam->ExBOn()) {
        locC = locC + calExBDetValue * driftlength;
      }

      // The electron position after diffusion and ExB in pad coordinates.
      // The pad row (z-direction)
      rowE       = padPlane->GetPadRowNumberROC(locR);
      if (rowE < 0) continue;
      rowOffset  = padPlane->GetPadRowOffsetROC(rowE,locR);

      // The pad column (rphi-direction)
      offsetTilt = padPlane->GetTiltOffset(rowOffset);
      colE       = padPlane->GetPadColNumber(locC+offsetTilt);
      if (colE < 0) continue;         
      colOffset  = padPlane->GetPadColOffset(colE,locC+offsetTilt);
	  
      // Also re-retrieve drift velocity because col and row may have changed
      driftvelocity = calVdriftDetValue * calVdriftROC->GetValue(colE,rowE);
      Float_t t0    = calT0DetValue     + calT0ROC->GetValue(colE,rowE);

      // Convert the position to drift time [mus], using either constant drift velocity or
      // time structure of drift cells (non-isochronity, GARFIELD calculation).
      // Also add absolute time of hits to take pile-up events into account properly
      Double_t drifttime;
      if (simParam->TimeStructOn()) {
	// Get z-position with respect to anode wire
        Double_t zz  =  row0 - locR + padPlane->GetAnodeWireOffset();
	zz -= ((Int_t)(2 * zz)) / 2.0;
        if (zz > 0.25) {
          zz  = 0.5 - zz;
        }
        // Use drift time map (GARFIELD)
        drifttime = commonParam->TimeStruct(driftvelocity,0.5*kAmWidth-1.0*locT,zz)
                  + hittime;
      } 
      else {
	// Use constant drift velocity
        drifttime = TMath::Abs(locT) / driftvelocity
                  + hittime;
      }

      // Apply the gas gain including fluctuations
      Double_t ggRndm = 0.0;
      do {
        ggRndm = gRandom->Rndm();
      } while (ggRndm <= 0);
      Double_t signal = -(simParam->GetGasGain()) * TMath::Log(ggRndm);

      // Apply the pad response 
      if (simParam->PRFOn()) {
        // The distance of the electron to the center of the pad 
        // in units of pad width
        Double_t dist = (colOffset - 0.5*padPlane->GetColSize(colE))
                      / padPlane->GetColSize(colE);
        // This is a fixed parametrization, i.e. not dependent on
        // calibration values !
        if (!(calibration->PadResponse(signal,dist,layer,padSignal))) continue;
      }
      else {
        padSignal[0] = 0.0;
        padSignal[1] = signal;
        padSignal[2] = 0.0;
      }

      // The time bin (always positive), with t0 distortion
      Double_t timeBinIdeal = drifttime * samplingRate + t0;
      // Protection 
      if (TMath::Abs(timeBinIdeal) > 2*nTimeTotal) {
        timeBinIdeal = 2 * nTimeTotal;
      }
      Int_t    timeBinTruncated = ((Int_t) timeBinIdeal);
      // The distance of the position to the middle of the timebin
      Double_t timeOffset       = ((Float_t) timeBinTruncated 
                                + 0.5 - timeBinIdeal) / samplingRate;
          
      // Sample the time response inside the drift region
      // + additional time bins before and after.
      // The sampling is done always in the middle of the time bin
      for (Int_t iTimeBin = TMath::Max(timeBinTruncated,0)
          ;iTimeBin < TMath::Min(timeBinTruncated+timeBinTRFend,nTimeTotal)
	  ;iTimeBin++) {

        // Apply the time response
        Double_t timeResponse = 1.0;
        Double_t crossTalk    = 0.0;
        Double_t time         = (iTimeBin - timeBinTruncated) / samplingRate + timeOffset;

        if (simParam->TRFOn()) {
          timeResponse = simParam->TimeResponse(time);
        }
        if (simParam->CTOn()) {
          crossTalk    = simParam->CrossTalk(time);
        }

        signalOld[0] = 0.0;
        signalOld[1] = 0.0;
        signalOld[2] = 0.0;

        for (iPad = 0; iPad < kNpad; iPad++) {

          Int_t colPos = colE + iPad - 1;
          if (colPos <        0) continue;
          if (colPos >= nColMax) break;

          // Add the signals
          signalOld[iPad] = signals->GetData(rowE,colPos,iTimeBin);

          if (colPos != colE) {
	    // Cross talk added to non-central pads
            signalOld[iPad] += padSignal[iPad] 
	                     * (timeResponse + crossTalk);
          } 
          else {
	    // W/o cross talk at central pad
            signalOld[iPad] += padSignal[iPad] 
	                     * timeResponse;
          }

          signals->SetData(rowE,colPos,iTimeBin,signalOld[iPad]);

          // Store the track index in the dictionary
          // Note: We store index+1 in order to allow the array to be compressed
	  // Note2: Taking out the +1 in track
          if (signalOld[iPad] > 0.0) { 
            for (dict = 0; dict < kNdict; dict++) {
	      Int_t oldTrack = dictionary[dict]->GetData(rowE,colPos,iTimeBin); 
	      if (oldTrack == track) break;
	      if (oldTrack ==  -1 ) {
		dictionary[dict]->SetData(rowE,colPos,iTimeBin,track);
	        break;
	      }
            }
          }

	} // Loop: pads

      } // Loop: time bins

    } // Loop: electrons of a single hit

  } // Loop: hits

  AliDebug(2,Form("Finished analyzing %d hits",nhit));

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::ConvertSignals(Int_t det, AliTRDarraySignal *signals)
{
  //
  // Convert signals to digits
  //

  AliDebug(1,Form("Start converting the signals for detector %d",det));

  if (fSDigits) {
    // Convert the signal array to s-digits
    if (!Signal2SDigits(det,signals)) {
      return kFALSE;
    }
  }
  else {
    // Convert the signal array to digits
    if (!Signal2ADC(det,signals)) {
      return kFALSE;
    }
    // Run digital processing for digits
    RunDigitalProcessing(det);
  }

  // Compress the arrays
  CompressOutputArrays(det);   

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::Signal2ADC(Int_t det, AliTRDarraySignal *signals)
{
  //
  // Converts the sampled electron signals to ADC values for a given chamber
  //

  AliDebug(1,Form("Start converting signals to ADC values for detector=%d",det));

  AliTRDcalibDB     *calibration = AliTRDcalibDB::Instance();
  if (!calibration) {
    AliFatal("Could not get calibration object");
    return kFALSE;
  }

  AliTRDSimParam    *simParam    = AliTRDSimParam::Instance();
  if (!simParam) {
    AliFatal("Could not get simulation parameters");
    return kFALSE;
  }  

  // Converts number of electrons to fC
  const Double_t kEl2fC = 1.602e-19 * 1.0e15; 

  // Coupling factor
  Double_t coupling     = simParam->GetPadCoupling() 
                        * simParam->GetTimeCoupling();
  // Electronics conversion factor
  Double_t convert      = kEl2fC 
                        * simParam->GetChipGain();
  // ADC conversion factor
  Double_t adcConvert   = simParam->GetADCoutRange()
                        / simParam->GetADCinRange();
  // The electronics baseline in mV
  Double_t baseline     = simParam->GetADCbaseline() 
                        / adcConvert;
  // The electronics baseline in electrons
  Double_t baselineEl   = baseline
                        / convert;

  Int_t row  = 0;
  Int_t col  = 0;
  Int_t time = 0;

  Int_t nRowMax    = fGeo->GetPadPlane(det)->GetNrows();
  Int_t nColMax    = fGeo->GetPadPlane(det)->GetNcols();
  Int_t nTimeTotal = fDigitsManager->GetDigitsParam()->GetNTimeBins(det);
  if (fSDigitsManager->GetDigitsParam()->GetNTimeBins(det)) {
    nTimeTotal = fSDigitsManager->GetDigitsParam()->GetNTimeBins(det);
  }
  else {
    AliFatal("Could not get number of time bins");
    return kFALSE;
  }

  // The gain factor calibration objects
  const AliTRDCalDet *calGainFactorDet      = calibration->GetGainFactorDet();  
  AliTRDCalROC       *calGainFactorROC      = 0x0;
  Float_t             calGainFactorDetValue = 0.0;

  AliTRDarrayADC     *digits = 0x0;

  if (!signals) {
    AliError(Form("Signals array for detector %d does not exist\n",det));
    return kFALSE;
  }
  if (signals->HasData()) {
    // Expand the container if neccessary
    signals->Expand();   
  }
  else {
    // Create missing containers
    signals->Allocate(nRowMax,nColMax,nTimeTotal);      
  }

  // Get the container for the digits of this detector
  if (fDigitsManager->HasSDigits()) {
    AliError("Digits manager has s-digits");
    return kFALSE;
  }

  digits = (AliTRDarrayADC *) fDigitsManager->GetDigits(det);
  // Allocate memory space for the digits buffer
  if (!digits->HasData()) {
    digits->Allocate(nRowMax,nColMax,nTimeTotal);
  }

  // Get the calibration objects
  calGainFactorROC      = calibration->GetGainFactorROC(det);
  calGainFactorDetValue = calGainFactorDet->GetValue(det);

  // Get the online gain factors
  //AliTRDCalOnlineGainTableROC *onlineGainFactorROC 
  //  = calibration->GetOnlineGainTableROC(det);

  // Create the digits for this chamber
  for (row  = 0; row  <  nRowMax; row++ ) {
    for (col  = 0; col  <  nColMax; col++ ) {

      // Check whether pad is masked
      // Bridged pads are not considered yet!!!
      if (calibration->IsPadMasked(det,col,row) || 
          calibration->IsPadNotConnected(det,col,row)) {
        continue;
      }

      // The gain factors
      Float_t padgain = calGainFactorDetValue 
                      * calGainFactorROC->GetValue(col,row);
      if (padgain <= 0) {
        AliError(Form("Not a valid gain %f, %d %d %d",padgain,det,col,row));
      }

      for (time = 0; time < nTimeTotal; time++) {

	// Get the signal amplitude
        Float_t signalAmp = signals->GetData(row,col,time);
        // Pad and time coupling
        signalAmp *= coupling;
	// Gain factors
	signalAmp *= padgain;

        // Add the noise, starting from minus ADC baseline in electrons
        signalAmp  = TMath::Max((Double_t) gRandom->Gaus(signalAmp,simParam->GetNoise())
                               ,-baselineEl);

        // Convert to mV
        signalAmp *= convert;
        // Add ADC baseline in mV
        signalAmp += baseline;

 	// Convert to ADC counts. Set the overflow-bit fADCoutRange if the
	// signal is larger than fADCinRange
        Short_t adc  = 0;
        if (signalAmp >= simParam->GetADCinRange()) {
          adc = ((Short_t) simParam->GetADCoutRange());
	}
        else {
          adc = TMath::Nint(signalAmp * adcConvert);
	}

        // Saving all digits
	digits->SetData(row,col,time,adc);

      } // for: time

    } // for: col
  } // for: row

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::Signal2SDigits(Int_t det, AliTRDarraySignal *signals)
{
  //
  // Converts the sampled electron signals to s-digits
  //

  AliDebug(1,Form("Start converting signals to s-digits for detector=%d",det));

  AliTRDcalibDB     *calibration = AliTRDcalibDB::Instance();
  if (!calibration) {
    AliFatal("Could not get calibration object");
    return kFALSE;
  }

  Int_t row  = 0;
  Int_t col  = 0;
  Int_t time = 0;

  Int_t nRowMax    = fGeo->GetPadPlane(det)->GetNrows();
  Int_t nColMax    = fGeo->GetPadPlane(det)->GetNcols();
  Int_t nTimeTotal = fDigitsManager->GetDigitsParam()->GetNTimeBins(det);

  // Get the container for the digits of this detector
  if (!fDigitsManager->HasSDigits()) {
    AliError("Digits manager has no s-digits");
    return kFALSE;
  }

  AliTRDarraySignal *digits = (AliTRDarraySignal *) fDigitsManager->GetSDigits(det);
  // Allocate memory space for the digits buffer
  if (!digits->HasData()) {
    digits->Allocate(nRowMax,nColMax,nTimeTotal);
  }

  // Create the sdigits for this chamber
  for (row  = 0; row  <  nRowMax; row++ ) {
    for (col  = 0; col  <  nColMax; col++ ) {
      for (time = 0; time < nTimeTotal; time++) {         
        digits->SetData(row,col,time,signals->GetData(row,col,time));
      } // for: time
    } // for: col
  } // for: row
  
  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::Digits2SDigits(AliTRDdigitsManager * const manDig
                                     , AliTRDdigitsManager * const manSDig)
{
  //
  // Converts digits into s-digits. Needed for embedding into real data.
  //

  AliDebug(1,"Start converting digits to s-digits");

  if (!fGeo) {
    fGeo = new AliTRDgeometry();
  }

  AliTRDcalibDB     *calibration = AliTRDcalibDB::Instance();
  if (!calibration) {
    AliFatal("Could not get calibration object");
    return kFALSE;
  }

  AliTRDSimParam    *simParam    = AliTRDSimParam::Instance();
  if (!simParam) {
    AliFatal("Could not get simulation parameters");
    return kFALSE;
  }  

  // Converts number of electrons to fC
  const Double_t kEl2fC = 1.602e-19 * 1.0e15; 

  // Coupling factor
  Double_t coupling     = simParam->GetPadCoupling() 
                        * simParam->GetTimeCoupling();
  // Electronics conversion factor
  Double_t convert      = kEl2fC 
                        * simParam->GetChipGain();
  // ADC conversion factor
  Double_t adcConvert   = simParam->GetADCoutRange()
                        / simParam->GetADCinRange();
  // The electronics baseline in mV
  Double_t baseline     = simParam->GetADCbaseline() 
                        / adcConvert;
  // The electronics baseline in electrons
  //Double_t baselineEl   = baseline
  //                      / convert;

  // The gainfactor calibration objects
  //const AliTRDCalDet *calGainFactorDet      = calibration->GetGainFactorDet();  
  //AliTRDCalROC       *calGainFactorROC      = 0;
  //Float_t             calGainFactorDetValue = 0.0;

  Int_t row  = 0;
  Int_t col  = 0;
  Int_t time = 0;

  for (Int_t det = 0; det < AliTRDgeometry::Ndet(); det++) {

    Int_t nRowMax    = fGeo->GetPadPlane(det)->GetNrows();
    Int_t nColMax    = fGeo->GetPadPlane(det)->GetNcols();
    Int_t nTimeTotal = manDig->GetDigitsParam()->GetNTimeBins(det);

    // Get the calibration objects
    //calGainFactorROC      = calibration->GetGainFactorROC(det);
    //calGainFactorDetValue = calGainFactorDet->GetValue(det);

    // Get the digits
    AliTRDarrayADC *digits = (AliTRDarrayADC *) manDig->GetDigits(det);

    if (!manSDig->HasSDigits()) {
      AliError("SDigits manager has no s-digits");
      return kFALSE;
    }
    // Get the s-digits
    AliTRDarraySignal     *sdigits = (AliTRDarraySignal *)     manSDig->GetSDigits(det);
    AliTRDarrayDictionary *tracks0 = (AliTRDarrayDictionary *) manSDig->GetDictionary(det,0);
    AliTRDarrayDictionary *tracks1 = (AliTRDarrayDictionary *) manSDig->GetDictionary(det,1);
    AliTRDarrayDictionary *tracks2 = (AliTRDarrayDictionary *) manSDig->GetDictionary(det,2);
    // Allocate memory space for the digits buffer
    sdigits->Allocate(nRowMax,nColMax,nTimeTotal);
    tracks0->Allocate(nRowMax,nColMax,nTimeTotal);
    tracks1->Allocate(nRowMax,nColMax,nTimeTotal);
    tracks2->Allocate(nRowMax,nColMax,nTimeTotal);

    // Keep the digits param
    manSDig->GetDigitsParam()->SetNTimeBinsAll(manDig->GetDigitsParam()->GetNTimeBins(0));
    manSDig->GetDigitsParam()->SetADCbaselineAll(manDig->GetDigitsParam()->GetADCbaseline(0));

    if (digits->HasData()) {

      digits->Expand();

      // Create the sdigits for this chamber
      for (row  = 0; row  <  nRowMax; row++ ) {
        for (col  = 0; col  <  nColMax; col++ ) {

          // The gain factors
          //Float_t padgain = calGainFactorDetValue 
          //                * calGainFactorROC->GetValue(col,row);

          for (time = 0; time < nTimeTotal; time++) {

            Short_t  adcVal = digits->GetData(row,col,time);
            Double_t signal = (Double_t) adcVal;
            // ADC -> signal in mV
            signal /= adcConvert;
            // Subtract baseline in mV
            signal -= baseline;
            // Signal in mV -> signal in #electrons
            signal /= convert;
            // Gain factor
            //signal /= padgain; // Not needed for real data
            // Pad and time coupling
            signal /= coupling;

            sdigits->SetData(row,col,time,signal);
            tracks0->SetData(row,col,time,0);
            tracks1->SetData(row,col,time,0);
            tracks2->SetData(row,col,time,0);

          } // for: time

        } // for: col
      } // for: row
  
    } // if: has data

    sdigits->Compress(0);
    tracks0->Compress();
    tracks1->Compress();
    tracks2->Compress();

    // No compress just remove
    manDig->RemoveDigits(det);
    manDig->RemoveDictionaries(det);      

  } // for: det

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::SDigits2Digits()
{
  //
  // Merges the input s-digits and converts them to normal digits
  //

  if (!MergeSDigits()) {
    return kFALSE;
  }

  return ConvertSDigits();

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

  AliTRDSimParam    *simParam    = AliTRDSimParam::Instance();
  if (!simParam) {
    AliFatal("Could not get simulation parameters");
    return kFALSE;
  }
  
  AliTRDcalibDB     *calibration = AliTRDcalibDB::Instance();
  if (!calibration) {
    AliFatal("Could not get calibration object");
    return kFALSE;
  }
  
  Int_t iDict = 0;
  Int_t jDict = 0;

  AliTRDarraySignal     *digitsA;
  AliTRDarraySignal     *digitsB;
  AliTRDarrayDictionary *dictionaryA[kNDict];
  AliTRDarrayDictionary *dictionaryB[kNDict];

  AliTRDdigitsManager   *mergeSDigitsManager = 0x0;
  // Get the first s-digits
  fSDigitsManager = (AliTRDdigitsManager *) fSDigitsManagerList->First();
  if (!fSDigitsManager) { 
    AliError("No SDigits manager");
    return kFALSE;
  }

  // Loop through the other sets of s-digits
  mergeSDigitsManager = (AliTRDdigitsManager *) fSDigitsManagerList->After(fSDigitsManager);

  if (mergeSDigitsManager) {
    AliDebug(1,Form("Merge %d input files.",fSDigitsManagerList->GetSize()));
  }
  else {
    AliDebug(1,"Only one input file.");
  }
  
  Int_t iMerge = 0;

  while (mergeSDigitsManager) {

    iMerge++;

    // Loop through the detectors
    for (Int_t iDet = 0; iDet < AliTRDgeometry::Ndet(); iDet++) {

      Int_t nTimeTotal = fSDigitsManager->GetDigitsParam()->GetNTimeBins(iDet);
      if (mergeSDigitsManager->GetDigitsParam()->GetNTimeBins(iDet) != nTimeTotal) {
        AliError(Form("Mismatch in the number of time bins [%d,%d] in detector %d"
                     ,nTimeTotal
  		     ,mergeSDigitsManager->GetDigitsParam()->GetNTimeBins(iDet)
		     ,iDet));
        return kFALSE;
      }

      Int_t nRowMax = fGeo->GetPadPlane(iDet)->GetNrows();
      Int_t nColMax = fGeo->GetPadPlane(iDet)->GetNcols();
	  
      // Loop through the pixels of one detector and add the signals
      digitsA = (AliTRDarraySignal *) fSDigitsManager->GetSDigits(iDet);    
      digitsB = (AliTRDarraySignal *) mergeSDigitsManager->GetSDigits(iDet); 
      digitsA->Expand();  
      if (!digitsA->HasData()) continue;
      digitsB->Expand();    
      if (!digitsB->HasData()) continue;
	  
      for (iDict = 0; iDict < kNDict; iDict++) {
	dictionaryA[iDict] = (AliTRDarrayDictionary *) fSDigitsManager->GetDictionary(iDet,iDict);
	dictionaryB[iDict] = (AliTRDarrayDictionary *) mergeSDigitsManager->GetDictionary(iDet,iDict);
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
	      
	AliDebug(1,Form("Merge detector %d of input no.%d",iDet,iMerge+1));
	      
	for (Int_t iRow  = 0; iRow  <  nRowMax;   iRow++ ) {
	  for (Int_t iCol  = 0; iCol  <  nColMax;   iCol++ ) {
	    for (Int_t iTime = 0; iTime < nTimeTotal; iTime++) {
         	
	      // Add the amplitudes of the summable digits 
	      Float_t ampA = digitsA->GetData(iRow,iCol,iTime);
	      Float_t ampB = digitsB->GetData(iRow,iCol,iTime);
	      ampA += ampB;
	      digitsA->SetData(iRow,iCol,iTime,ampA);

	      // Add the mask to the track id if defined.
	      for (iDict = 0; iDict < kNDict; iDict++) {
		Int_t trackB = dictionaryB[iDict]->GetData(iRow,iCol,iTime);
	        if ((fMasks) && (trackB > 0))  {
		  for (jDict = 0; jDict < kNDict; jDict++) { 
		    Int_t trackA = dictionaryA[iDict]->GetData(iRow,iCol,iTime); 
		    if (trackA == 0) {
		      trackA = trackB + fMasks[iMerge];
		      dictionaryA[iDict]->SetData(iRow,iCol,iTime,trackA);  
		    } // if:  track A == 0
		  } // for: jDict
		} // if:  fMasks and trackB > 0
	      } // for: iDict

	    } // for: iTime
	  } // for: iCol
        } // for: iRow

      } // if:  doMerge

      mergeSDigitsManager->RemoveDigits(iDet);
      mergeSDigitsManager->RemoveDictionaries(iDet);
  
      if (fCompress) {
        digitsA->Compress(0); 
        for (iDict = 0; iDict < kNDict; iDict++) {                                     
	  dictionaryA[iDict]->Compress();
        }
      }
      
    } // for: detectors    
      
    // The next set of s-digits
    mergeSDigitsManager = (AliTRDdigitsManager *) fSDigitsManagerList->After(mergeSDigitsManager);
    
  } // while: mergeDigitsManagers
  
  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::ConvertSDigits()
{
  //
  // Converts s-digits to normal digits
  //

  AliTRDarraySignal *digitsIn = 0x0;

  if (!fSDigitsManager->HasSDigits()) {
    AliError("No s-digits in digits manager");
    return kFALSE;
  }

  // Loop through the detectors
  for (Int_t det = 0; det < AliTRDgeometry::Ndet(); det++) {

    // Get the merged s-digits (signals)
    digitsIn = (AliTRDarraySignal *) fSDigitsManager->GetSDigits(det);
    if (!digitsIn->HasData()) {
      AliDebug(2,Form("No digits for det=%d",det));
      continue;
    }
    
    // Convert the merged sdigits to digits
    if (!Signal2ADC(det,digitsIn)) {
      continue;
    }

    // Copy the dictionary information to the output array
    if (!CopyDictionary(det)) {
      continue;
    }

    // Delete 
    fSDigitsManager->RemoveDigits(det);
    fSDigitsManager->RemoveDictionaries(det);

    // Run digital processing
    RunDigitalProcessing(det);

    // Compress the arrays
    CompressOutputArrays(det);

  } // for: detector numbers

  if (AliDataLoader *trklLoader = AliRunLoader::Instance()->GetLoader("TRDLoader")->GetDataLoader("tracklets")) {
    if (trklLoader->Tree())
      trklLoader->WriteData("OVERWRITE");
  }

  // Save the values for the raw data headers
  if (AliTRDSimParam::Instance()->GetNTBoverwriteOCDB()) {
    fDigitsManager->GetDigitsParam()->SetNTimeBinsAll(AliTRDSimParam::Instance()->GetNTimeBins());
  }
  else {
    fDigitsManager->GetDigitsParam()->SetNTimeBinsAll(AliTRDcalibDB::Instance()->GetNumberOfTimeBinsDCS());
  }
  fDigitsManager->GetDigitsParam()->SetADCbaselineAll(AliTRDSimParam::Instance()->GetADCbaseline());

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::CopyDictionary(Int_t det)
{
  //
  // Copies the dictionary information from the s-digits arrays
  // to the output arrays
  //

  AliTRDcalibDB     *calibration = AliTRDcalibDB::Instance();
  if (!calibration) {
    AliFatal("Could not get calibration object");
    return kFALSE;
  }

  AliDebug(1,Form("Start copying dictionaries for detector=%d",det));

  const Int_t kNDict = AliTRDdigitsManager::kNDict;
  AliTRDarrayDictionary *dictionaryIn[kNDict];
  AliTRDarrayDictionary *dictionaryOut[kNDict];

  Int_t nRowMax    = fGeo->GetPadPlane(det)->GetNrows();
  Int_t nColMax    = fGeo->GetPadPlane(det)->GetNcols();
  Int_t nTimeTotal = fSDigitsManager->GetDigitsParam()->GetNTimeBins(det);

  Int_t row  = 0;
  Int_t col  = 0;
  Int_t time = 0;
  Int_t dict = 0;

  for (dict = 0; dict < kNDict; dict++) {

    dictionaryIn[dict]  = (AliTRDarrayDictionary *) fSDigitsManager->GetDictionary(det,dict);
    dictionaryIn[dict]->Expand();
    dictionaryOut[dict] = (AliTRDarrayDictionary *) fDigitsManager->GetDictionary(det,dict);
    dictionaryOut[dict]->Allocate(nRowMax,nColMax,nTimeTotal);

    for (row = 0; row < nRowMax; row++) {
      for (col = 0; col < nColMax; col++) {
        for (time = 0; time < nTimeTotal; time++) {
	  Int_t track = dictionaryIn[dict]->GetData(row,col,time);
	  dictionaryOut[dict]->SetData(row,col,time,track);
	} // for: time
      } // for: col
    } // for: row
    
  } // for: dictionaries
  
  return kTRUE;

}

//_____________________________________________________________________________
void AliTRDdigitizer::CompressOutputArrays(Int_t det)
{
  //
  // Compress the output arrays
  //

  const Int_t kNDict = AliTRDdigitsManager::kNDict;
  AliTRDarrayDictionary *dictionary = 0x0;

  if (fCompress) {

    if (!fSDigits) {
      AliTRDarrayADC *digits = 0x0;  
      digits = (AliTRDarrayADC *) fDigitsManager->GetDigits(det);
      digits->Compress();
    }

    if (fSDigits) {
      AliTRDarraySignal *digits = 0x0; 
      digits = (AliTRDarraySignal *) fDigitsManager->GetSDigits(det);
      digits->Compress(0);
    }

    for (Int_t dict = 0; dict < kNDict; dict++) {
      dictionary = (AliTRDarrayDictionary *) fDigitsManager->GetDictionary(det,dict);
      dictionary->Compress();
    }

  }

}

//_____________________________________________________________________________
Bool_t AliTRDdigitizer::WriteDigits() const
{
  //
  // Writes out the TRD-digits and the dictionaries
  //

  // Write parameters
  fRunLoader->CdGAFile();

  // Store the digits and the dictionary in the tree
  return fDigitsManager->WriteDigits();

}

//_____________________________________________________________________________
void AliTRDdigitizer::InitOutput(Int_t iEvent)
{
  //
  // Initializes the output branches
  //

  fEvent = iEvent;
   
  if (!fRunLoader) {
    AliError("Run Loader is NULL");
    return;  
  }

  AliLoader *loader = fRunLoader->GetLoader("TRDLoader");
  if (!loader) {
    AliError("Can not get TRD loader from Run Loader");
    return;
  }

  TTree *tree = 0;
  
  if (fSDigits) { 
    // If we produce SDigits
    tree = loader->TreeS();
    if (!tree) {
      loader->MakeTree("S");
      tree = loader->TreeS();
    }
  }
  else {
    // If we produce Digits
    tree = loader->TreeD();
    if (!tree) {
      loader->MakeTree("D");
      tree = loader->TreeD();
    }
  }
  fDigitsManager->SetEvent(iEvent);
  fDigitsManager->MakeBranch(tree);

}
  
//_____________________________________________________________________________
Int_t AliTRDdigitizer::Diffusion(Float_t vdrift, Double_t absdriftlength
                               , Double_t exbvalue
                               , Double_t &lRow, Double_t &lCol, Double_t &lTime)
{
  //
  // Applies the diffusion smearing to the position of a single electron.
  // Depends on absolute drift length.
  //
  
  Float_t diffL = 0.0;
  Float_t diffT = 0.0;

  if (AliTRDCommonParam::Instance()->GetDiffCoeff(diffL,diffT,vdrift)) {

    Float_t driftSqrt = TMath::Sqrt(absdriftlength);
    Float_t sigmaT    = driftSqrt * diffT;
    Float_t sigmaL    = driftSqrt * diffL;
    lRow  = gRandom->Gaus(lRow ,sigmaT);
    if (AliTRDCommonParam::Instance()->ExBOn()) {
      lCol  = gRandom->Gaus(lCol ,sigmaT * 1.0 / (1.0 + exbvalue*exbvalue));
      lTime = gRandom->Gaus(lTime,sigmaL * 1.0 / (1.0 + exbvalue*exbvalue));
    }
    else {
      lCol  = gRandom->Gaus(lCol ,sigmaT);
      lTime = gRandom->Gaus(lTime,sigmaL);
    }

    return 1;

  }
  else {

    return 0;

  }

}
  
//_____________________________________________________________________________
void AliTRDdigitizer::RunDigitalProcessing(Int_t det)
{
  //
  // Run the digital processing in the TRAP
  //

  AliTRDfeeParam *feeParam = AliTRDfeeParam::Instance();

  AliTRDarrayADC *digits = fDigitsManager->GetDigits(det);
  if (!digits)
    return;

  //Call the methods in the mcm class using the temporary array as input  
  for(Int_t rob = 0; rob < digits->GetNrow() / 2; rob++)
  {
    for(Int_t mcm = 0; mcm < 16; mcm++)
    {
      fMcmSim->Init(det, rob, mcm); 
      fMcmSim->SetDataByPad(digits, fDigitsManager);
      fMcmSim->Filter();
      if (feeParam->GetTracklet()) {
        fMcmSim->Tracklet();
        fMcmSim->StoreTracklets();
      }
      fMcmSim->ZSMapping();
      fMcmSim->WriteData(digits);
    }
  }

}

