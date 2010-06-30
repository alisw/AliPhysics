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

// $Id$

#include "AliMUONPedestalEventGenerator.h"

#include "AliCodeTimer.h"
#include "AliDAQ.h"
#include "AliHeader.h"
#include "AliLog.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONRawWriter.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONVStore.h"
#include "AliMpCathodType.h"
#include "AliMpConstants.h"
#include "AliMpDEStore.h"
#include "AliMpDetElement.h"
#include "AliMpPlaneType.h"
#include "AliRawDataHeaderSim.h"
#include "AliRunLoader.h"
#include <TClonesArray.h>
#include <TMath.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TStopwatch.h>
#include <TSystem.h>

//-----------------------------------------------------------------------------
/// \class AliMUONPedestalEventGenerator
///
/// Generate simulated pedestal events for MUON TRK, to be able to e.g. test
/// online calibration routines.
///
/// The pedestals themselves are taken from the CDB. What we get from the CDB
/// is, per channel, the mean and the sigma of the pedestal. We then use
/// those informations to randomly get the pedestals for each channel, for
/// each event (picking in a gaus(mean,sigma)).
///
/// Output can be just digits, or digits + raw (ddl), or digits + raw (ddl)
/// + raw (date files, one per LDC), depending of ctor and MakeDDL() method.
///
/// \author L. Aphecetche
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONPedestalEventGenerator)
/// \endcond

Int_t AliMUONPedestalEventGenerator::fgCounter(0);

//std::streambuf* RedirectTo(std::ostream& what, std::ostream& to)
//{
//  std::streambuf* old = what.rdbuf();
//  
//  std::streambuf* psbuf = to.rdbuf();
//  what.rdbuf(psbuf);
//  
//  return old;
//}

//_____________________________________________________________________________
AliMUONPedestalEventGenerator::AliMUONPedestalEventGenerator(Int_t runNumber,
                                                             Int_t nevents,
                                                             const char* filename)
: TTask("AliMUONPedestalEventGenerator","Generate fake pedestal events"), 
fCalibrationData(new AliMUONCalibrationData(runNumber)),
fDateFileName(filename),
fGAliceFileName("galice.root"),
fMakeDDL(kTRUE),
fLoader(0x0),
fPedestals(fCalibrationData->Pedestals()),
fDigitStore(0x0),
fRawWriter(0x0)
{
  /// Will generate pedestals according to (mean,sigma)s found in CDB
  /// for run runNumber.
  /// Will generate nevents events
  /// If filename is != "", it will be the basename of the output LDC files
  ///
  if (!gSystem->IsAbsoluteFileName(fGAliceFileName)) 
  {
    char* absFileName = gSystem->ConcatFileName(gSystem->WorkingDirectory(),
                                                fGAliceFileName);
    fGAliceFileName = absFileName;
    delete[] absFileName;
  }
  
  AliRunLoader* runLoader = LoadRun("recreate");
  
  runLoader->SetNumberOfEventsPerFile(nevents);
  
  if (!runLoader)
  {
    AliError("Could not create RunLoader");
    return;
  }
  
  // Initialize event headers.
  runLoader->MakeTree("E");

  for ( Int_t iEvent = 0; iEvent < nevents; ++iEvent )
  {
    runLoader->SetEventNumber(iEvent);
    runLoader->GetHeader()->Reset(runNumber,iEvent);
    runLoader->TreeE()->Fill();
  }
  runLoader->WriteHeader("OVERWRITE");
  runLoader->CdGAFile();
  runLoader->Write(0, TObject::kOverwrite);  

  delete runLoader;
  fLoader = 0x0;
}


//_____________________________________________________________________________
AliMUONPedestalEventGenerator::~AliMUONPedestalEventGenerator()
{
  /// dtor
  delete fCalibrationData;
  AliInfo(Form("make a digit counter %d",fgCounter));
  delete fDigitStore;
  delete fRawWriter;
}

//_____________________________________________________________________________
Bool_t 
AliMUONPedestalEventGenerator::ConvertRawFilesToDate()
{
  /// convert raw data DDL files to DATE files with the program "dateStream".
  /// we make one file per LDC
  
  AliCodeTimerAuto("",0)
  
  AliInfo("Converting raw to date");
  
  const Int_t kIDet = AliDAQ::DetectorID("MUONTRK");
  
  const Int_t kNLDCs = 5;//TMath::CeilNint(AliDAQ::NumberOfLdcs(kIDet));
  
  char* path = gSystem->Which(gSystem->Getenv("PATH"), "dateStream");
  if (!path) 
  {
    AliError("the program dateStream was not found");
    return kFALSE;
  }
  
  delete[] path;
  
  AliRunLoader* runLoader = LoadRun("read");
  if (!runLoader) return kFALSE;
  
  AliInfo(Form("converting raw data DDL files to DATE files %s", fDateFileName.Data()));
  FILE** pipe = new FILE*[kNLDCs];
  
  for ( Int_t iFile = 0; iFile < kNLDCs; ++iFile)
  {
    char command[256];
    // Note the option -s. It is used in order to avoid
    // the generation of SOR/EOR events.
    sprintf(command, "dateStream -c -D -o %s.LDC%d -# %d -C", 
            fDateFileName.Data(), iFile, runLoader->GetNumberOfEvents());
    pipe[iFile] = gSystem->OpenPipe(command, "w");
  }
  
  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); ++iEvent) 
  {
    Float_t ldc = 0;
    Int_t prevLDC = -1;

    for (Int_t iDDL = 0; iDDL < AliDAQ::NumberOfDdls(kIDet); ++iDDL) 
    {        
      Int_t ddlID = AliDAQ::DdlID(kIDet,iDDL);
      Int_t ldcID = Int_t(ldc + 0.0001);
      ldc += AliDAQ::NumberOfLdcs(kIDet) / AliDAQ::NumberOfDdls(kIDet);
      
      char rawFileName[256];
      sprintf(rawFileName, "raw%d/%s", 
              iEvent, AliDAQ::DdlFileName(kIDet,iDDL));
      
      // check existence and size of raw data file
      FILE* file = fopen(rawFileName, "rb");
      if (!file) continue;
      fseek(file, 0, SEEK_END);
      unsigned long size = ftell(file);
      fclose(file);
      if (!size) continue;
      
      if (ldcID != prevLDC) {
        fprintf(pipe[ldcID], " LDC Id %d\n", ldcID);
        prevLDC = ldcID;
      }
      fprintf(pipe[ldcID], "  Equipment Id %d Payload %s\n", ddlID, rawFileName);
    }
  }
  
  Int_t result(0);
  
  for ( Int_t iFile = 0; iFile < kNLDCs; ++iFile)
  {
    result += gSystem->ClosePipe(pipe[iFile]);
  }
  
  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); ++iEvent) 
  {
    char command[256];
    sprintf(command, "rm -r raw%d", iEvent);
    gSystem->Exec(command);
  }
  
  delete [] pipe;
  delete runLoader;
  fLoader=0x0;
  return (result == 0);
}

//_____________________________________________________________________________
AliMUONVDigitStore*
AliMUONPedestalEventGenerator::DigitStore()
{
/// Return digt container; create it if it does not exist

  if (!fDigitStore) fDigitStore = AliMUONVDigitStore::Create("AliMUONDigitStoreV2R");
  return fDigitStore;
}

//_____________________________________________________________________________
void
AliMUONPedestalEventGenerator::Exec(Option_t*)
{  
  /// Main steering method
  
  AliCodeTimerAuto("",0)
  
  if (!fPedestals)
  {
    AliError("No pedestal store. Cannot proceed.");
    return;
  }
  
  AliRunLoader* runLoader = LoadRun("update");
  
  Int_t nevents = runLoader->GetNumberOfEvents();
    
  for ( Int_t i = 0; i < nevents ; ++i )
  {
    runLoader->GetEvent(i);
    
    fLoader->MakeDigitsContainer();
    TTree* treeD = fLoader->TreeD();  
    if (!treeD)
    {
      AliError(Form("Could not get TreeD for event %d",i));
      continue;
    }
    
    DigitStore()->Connect(*treeD);
        
    GenerateDigits(*(DigitStore()));

    // Fill the output treeD
    treeD->Fill();
    
    // Write to the output tree(D).
    // Please note that as GlobalTrigger, LocalTrigger and Digits are in the same
    // tree (=TreeD) in different branches, this WriteDigits in fact writes all of 
    // the 3 branches.

    AliCodeTimerStart("WriteDigits")
    fLoader->WriteDigits("OVERWRITE");
    AliCodeTimerStop("WriteDigits")
    
    fLoader->UnloadDigits();
    
    if ( fMakeDDL )
    {
      AliCodeTimerAuto("Digits2Raw",1);
      Digits2Raw(i);
    }
  }
    
  runLoader->WriteRunLoader("OVERWRITE");
  delete runLoader;
  fLoader = 0x0;
    
  // Finally, if instructed to do so, convert DDL files to DATE file(s)
  if ( fMakeDDL && fDateFileName.Length() > 0 ) 
  {
    AliCodeTimerAuto("ConvertRawFilesToDate",1)
    Bool_t dateOutput = ConvertRawFilesToDate();
    if (!dateOutput) 
    {
      AliError("DATE output failed. Exiting.");
      return;
    }    
  }
}

//_____________________________________________________________________________
void
AliMUONPedestalEventGenerator::Digits2Raw(Int_t event)
{
  /// Converts digits (from MUON.Digits.root file) to Raw DDL ascii files.
  
  AliCodeTimerAuto("",0)
  
  if (!fRawWriter) 
  {
      AliRawDataHeaderSim header;
      fRawWriter = new AliMUONRawWriter;
      fRawWriter->SetHeader(header);
  }
  
  // Generate RAW data from the digits
  // Be carefull to create&change to the correct directory first...
  
  TString baseDir = gSystem->WorkingDirectory();
  
  char dirName[256];
  sprintf(dirName, "raw%d", event);
  gSystem->MakeDirectory(dirName);
  if (!gSystem->ChangeDirectory(dirName)) 
  {
    AliError(Form("couldn't change to directory %s", dirName));
    return;
  }
  
  fRawWriter->Digits2Raw(DigitStore(),0);
  
  gSystem->ChangeDirectory(baseDir);
}

//_____________________________________________________________________________
void
AliMUONPedestalEventGenerator::GenerateDigits(AliMUONVDigitStore& digitStore)
{  
  /// Generate digits (where ADC is set to pedestal value) for all MUON TRK
  /// and for 1 event.
  
  AliCodeTimerAuto("",0)

  digitStore.Clear();
  
  Int_t ngenerated(0);
  Int_t nmanus(0);
  TIter next(fPedestals->CreateIterator());
  AliMUONVCalibParam* pedestals;
  
  while ( ( pedestals = static_cast<AliMUONVCalibParam*>(next())) )
  {
    Int_t detElemId = pedestals->ID0();
    Int_t manuId = pedestals->ID1();
    
    AliMpDetElement* de = AliMpDEStore::Instance()->GetDetElement(detElemId);
    AliMp::PlaneType planeType = AliMp::kBendingPlane;
    if ( manuId & AliMpConstants::ManuMask(AliMp::kNonBendingPlane) ) 
    {
      planeType = AliMp::kNonBendingPlane;
    }
    AliMp::CathodType cathode = de->GetCathodType(planeType);
    
    ++nmanus;
        
    for ( Int_t manuChannel = 0; manuChannel < pedestals->Size(); ++manuChannel )
    {
      Float_t mean = pedestals->ValueAsFloat(manuChannel,0);
      if (mean == AliMUONVCalibParam::InvalidFloatValue())
      {
        // This is a poor's man way of knowing if that channel really exists.
        // Better and safer way (but much slower too) would be to check pad existence
        // using AliMpVSegmentation::PadByLocation(manuId,manuChannel)
        continue;
      }
      else if ( mean < 1 || mean >  4095 ) 
      {
        AliFatal(Form("Got an invalid mean pedestal value for DE %d Manu %d"
                      " channel %d : mean = %e",detElemId,manuId,manuChannel,
                      mean));
      }
      else
      {
        Float_t sigma = pedestals->ValueAsFloat(manuChannel,1);
        
        if ( sigma < 0 ) 
        {
          AliWarning(Form("Got a negative sigma pedestal value for DE %d Manu %d"
                          " channel %d : sigma = %e, will use Abs()=%e",
                          detElemId,manuId,manuChannel,
                          sigma,-sigma));
          sigma = -sigma;
        }
        
        AliMUONVDigit* d = digitStore.Add(detElemId,manuId,manuChannel,
                                          cathode,
                                          AliMUONVDigitStore::kIgnore);
        
        Float_t ped = -1;
        while ( ped <= 0 ) 
        {
          ped = gRandom->Gaus(mean,sigma);
        }
        Int_t pedADC = TMath::FloorNint(ped);

        d->SetADC(pedADC);
        d->SetCharge(ped);
        // we do not set the remaining parts of the digit, as in principle
        // this is all we need : manuId, manuChannel and ADC, as far as
        // real data is concerned.
        ++fgCounter;
        ++ngenerated;
      }
    }
  }
  AliDebug(1,Form("ngenerated=%d nmanus=%d",ngenerated,nmanus));
}

//_____________________________________________________________________________
AliRunLoader*
AliMUONPedestalEventGenerator::LoadRun(const char* mode)
{
  /// Get access to AliRunLoader object
  while (AliRunLoader::Instance()) 
  {
    AliDebug(1,Form("Deleting AliRunLoader %p",AliRunLoader::Instance()));
    delete AliRunLoader::Instance();
  }
  
  AliRunLoader* runLoader = 
    AliRunLoader::Open(fGAliceFileName,AliConfig::GetDefaultEventFolderName(), 
                       mode);

  AliDebug(1,Form("AliRunLoader(%s)=%p",mode,runLoader));
  
  if (!runLoader) 
  {
    AliError("No run loader found in file galice.root");
  }
    
  TString smode(mode);
  smode.ToUpper();
  
  if (smode.Contains("RECREATE"))
  {
    AliInfo("Creating folder structure");
    AliConfig::Instance()
    ->CreateDetectorFolders(runLoader->GetEventFolder(), 
                            "MUON", "MUON");
    fLoader = new AliLoader("MUON",runLoader->GetEventFolder());
    runLoader->AddLoader(fLoader);
  }
  
  fLoader = static_cast<AliLoader*>(runLoader->GetDetectorLoader("MUON"));
    
  return runLoader;
}
