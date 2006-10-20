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

#include "AliDAQ.h"
#include "AliHeader.h"
#include "AliLog.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONData.h"
#include "AliMUONDigit.h"
#include "AliMUONRawWriter.h"
#include "AliMUONVCalibParam.h"
#include "AliMpIntPair.h"
#include "AliMpManuList.h"
#include "AliRunLoader.h"

#include "TClonesArray.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TMath.h"

///
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
///

/// \cond CLASSIMP
ClassImp(AliMUONPedestalEventGenerator)
/// \endcond

Int_t AliMUONPedestalEventGenerator::fgCounter(0);

//_____________________________________________________________________________
AliMUONPedestalEventGenerator::AliMUONPedestalEventGenerator(Int_t runNumber,
                                                             Int_t nevents,
                                                             const char* filename)
: TTask("AliMUONPedestalEventGenerator","Generate fake pedestal events"), 
fManuList(AliMpManuList::ManuList()),
fCalibrationData(new AliMUONCalibrationData(runNumber)),
fDateFileName(filename),
fGAliceFileName("galice.root"),
fMakeDDL(kTRUE)
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
  
  AliConfig::Instance()
    ->CreateDetectorFolders(runLoader->GetEventFolder(), 
                            "MUON", "MUON");
  
  AliLoader* loader = new AliLoader("MUON",runLoader->GetEventFolder());
  
	runLoader->AddLoader(loader);
  
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

  assert(runLoader->GetNumberOfEvents()==nevents);
  assert(runLoader->GetHeader()->GetRun()==runNumber);
  
  delete runLoader;
}


//_____________________________________________________________________________
AliMUONPedestalEventGenerator::~AliMUONPedestalEventGenerator()
{
  /// dtor
  delete fManuList;
  delete fCalibrationData;
  AliInfo(Form("make a digit counter %d",fgCounter));
}

//_____________________________________________________________________________
Bool_t 
AliMUONPedestalEventGenerator::ConvertRawFilesToDate()
{
  /// convert raw data DDL files to DATE files with the program "dateStream".
  /// we make one file per LDC
  
  const Int_t kIDet = AliDAQ::DetectorID("MUONTRK");
  
  const Int_t kNLDCs = TMath::CeilNint(AliDAQ::NumberOfLdcs(kIDet));
  
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
    sprintf(command, "dateStream -s -D -o %s.LDC%d -# %d -C", 
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
  
  delete [] pipe;
  delete runLoader;
  return (result == 0);
}

//_____________________________________________________________________________
void
AliMUONPedestalEventGenerator::Exec(Option_t*)
{  
  /// Main steering method
  
  TStopwatch timer;
  timer.Start(kTRUE);
  
  TStopwatch writeTimer;
  writeTimer.Start(kTRUE); writeTimer.Stop();
  TStopwatch fillTimer;
  fillTimer.Start(kTRUE); fillTimer.Stop();
  TStopwatch resetTimer;
  resetTimer.Start(kTRUE); resetTimer.Stop();
  TStopwatch getEventTimer;
  getEventTimer.Start(kTRUE); getEventTimer.Stop();
  TStopwatch branchTimer;
  branchTimer.Start(kTRUE); branchTimer.Stop();
  TStopwatch generateDigitsTimer;
  generateDigitsTimer.Start(kTRUE); generateDigitsTimer.Stop();
  TStopwatch digits2RawTimer;
  digits2RawTimer.Start(kTRUE); digits2RawTimer.Stop();
  TStopwatch convertRawFilesToDateTimer;
  convertRawFilesToDateTimer.Start(kTRUE); convertRawFilesToDateTimer.Stop();
  TStopwatch getTimer;

  getTimer.Start(kTRUE); 
  AliMUONData* data = GetDataAccess("update");
  AliRunLoader* runLoader = data->GetLoader()->GetRunLoader();  
  getTimer.Stop();
  
  for ( Int_t i = 0; i < runLoader->GetNumberOfEvents(); ++i )
  {
    getEventTimer.Start(kFALSE);
    runLoader->GetEvent(i);
    getEventTimer.Stop();
    
    branchTimer.Start(kFALSE);
    if ( data->TreeD() == 0x0 )
    {
      AliDebug(1,"Calling MakeDigitsContainer");
      data->GetLoader()->MakeDigitsContainer();
    }
    data->MakeBranch("D,GLT");
    data->SetTreeAddress("D,GLT");
    branchTimer.Stop();
    
    generateDigitsTimer.Start(kFALSE);
    GenerateDigits(data);
    generateDigitsTimer.Stop();
    
    fillTimer.Start(kTRUE);
    // Fill the output treeD
    data->Fill("D,GLT");
    fillTimer.Stop();
    
    // Write to the output tree(D).
    // Please note that as GlobalTrigger, LocalTrigger and Digits are in the same
    // tree (=TreeD) in different branches, this WriteDigits in fact writes all of 
    // the 3 branches.
    
    writeTimer.Start(kFALSE);
    data->GetLoader()->WriteDigits("OVERWRITE");
    writeTimer.Stop();
    
    // Finally, we clean up after ourselves.
    resetTimer.Start(kTRUE);
    data->ResetDigits();
    data->ResetTrigger();
    data->GetLoader()->UnloadDigits();
    resetTimer.Stop();
  }
  
  TStopwatch endTimer;
  endTimer.Start(kTRUE);
  runLoader->WriteRunLoader("OVERWRITE");
  delete data;
  delete runLoader;
  endTimer.Stop();
  
  if ( fMakeDDL ) 
  {
    // Now convert the digits.root file(s) to DDL files
    digits2RawTimer.Start(kFALSE);
    Digits2Raw();
    digits2RawTimer.Stop();
  }
  
  // Finally, if instructed to do so, convert DDL files to DATE file(s)
  if ( fMakeDDL && fDateFileName.Length() > 0 ) 
  {
    convertRawFilesToDateTimer.Start(kFALSE);
    AliDebug(1,"Converting to DATE file(s)");
    
    Bool_t dateOutput = ConvertRawFilesToDate();
    convertRawFilesToDateTimer.Stop();
    if (!dateOutput) 
    {
      AliError("DATE output failed. Aborting.");
      return;
    }    
  }
  
  AliInfo(Form("Execution time Exec : R:%.2fs C:%.2fs",
               timer.RealTime(),timer.CpuTime()));

  AliInfo(Form("  Execution time branch : R:%.2fs C:%.2fs",
               branchTimer.RealTime(),branchTimer.CpuTime()));
  AliInfo(Form("  Execution time getEvent : R:%.2fs C:%.2fs",
               getEventTimer.RealTime(),getEventTimer.CpuTime()));
  AliInfo(Form("  Execution time fill digits : R:%.2fs C:%.2fs",
               fillTimer.RealTime(),fillTimer.CpuTime()));
  AliInfo(Form("  Execution time write digits : R:%.2fs C:%.2fs",
               writeTimer.RealTime(),writeTimer.CpuTime()));
  AliInfo(Form("  Execution time reset digits : R:%.2fs C:%.2fs",
               resetTimer.RealTime(),resetTimer.CpuTime()));
  AliInfo(Form("  Execution time for GenerateDigits : R:%.2fs C:%.2fs",
               generateDigitsTimer.RealTime(),generateDigitsTimer.CpuTime()));
  AliInfo(Form("  Execution time for Digits2Raw : R:%.2fs C:%.2fs",
               digits2RawTimer.RealTime(),digits2RawTimer.CpuTime()));
  AliInfo(Form("  Execution time for ConvertRawFilesToDate : R:%.2fs C:%.2fs",
               convertRawFilesToDateTimer.RealTime(),convertRawFilesToDateTimer.CpuTime()));
  AliInfo(Form("  Execution time for get : R:%.2fs C:%.2fs",
               getTimer.RealTime(),getTimer.CpuTime()));
  AliInfo(Form("  Execution time for end : R:%.2fs C:%.2fs",
               endTimer.RealTime(),endTimer.CpuTime()));
}

//_____________________________________________________________________________
void
AliMUONPedestalEventGenerator::Digits2Raw()
{
  /// Converts digits (from MUON.Digits.root file) to Raw DDL ascii files.
  
  AliMUONData* data = GetDataAccess("read");
  AliRunLoader* runLoader = data->GetLoader()->GetRunLoader();
  AliDebug(1,Form("runLoader=%p",runLoader));
  
  AliMUONRawWriter rawWriter(data);
  
  // Generate RAW data from the digits
  // Be carefull to create&change to the correct directory first...

  TString baseDir = gSystem->WorkingDirectory();

  for ( Int_t i = 0; i < runLoader->GetNumberOfEvents(); ++i )
  {
    runLoader->GetEvent(i);
    
    AliDebug(1,Form("processing event %d", i));
    
    char dirName[256];
    sprintf(dirName, "raw%d", i);
    gSystem->MakeDirectory(dirName);
    if (!gSystem->ChangeDirectory(dirName)) 
    {
      AliError(Form("couldn't change to directory %s", dirName));
      return;
    }

    rawWriter.Digits2Raw();
    
    gSystem->ChangeDirectory(baseDir);
  }
  
  delete data;
  delete runLoader;
  
  AliDebug(1,Form("DDL files written successfully"));    
}

//_____________________________________________________________________________
void
AliMUONPedestalEventGenerator::GenerateDigits(AliMUONData* data)
{  
  /// Generate digits (where ADC is set to pedestal value) for all MUON TRK
  /// and for 1 event.
  
  TIter next(fManuList);
  AliMpIntPair* p;
  Int_t ngenerated(0);
  
  while ( ( p = static_cast<AliMpIntPair*>(next())) )
  {
    Int_t detElemId = p->GetFirst();

    Int_t chamber = detElemId/100-1;
        
    TClonesArray* pdigits = data->Digits(chamber);
    if (!pdigits)
    {
      AliError(Form("No digits for chamber %d",1));
      continue;
    }
    TClonesArray& digits = *pdigits;
    
    Int_t manuId = p->GetSecond();
    
    AliMUONVCalibParam* pedestals = fCalibrationData->Pedestals(detElemId,manuId);
    if (!pedestals)
    {
      AliError(Form("Could not find pedestals for (DE,MANU)=(%d,%d)",detElemId,
                    manuId));
      return;
    }
    for ( Int_t manuChannel = 0; manuChannel < pedestals->Size(); ++manuChannel )
    {
      Float_t mean = pedestals->ValueAsFloat(manuChannel,0);
      if (mean == AliMUONVCalibParam::InvalidFloatValue())
      {
        // This is a poor's man way of knowing if that channel really exists.
        // Better and safer way (but much slower too) would be to check pad existence
        // using AliMpVSegmentation::PadByLocation(AliMpIntPair(manuId,manuChannel))
        continue;
      }
      else
      {
        Float_t sigma = pedestals->ValueAsFloat(manuChannel,1);
        AliMUONDigit d;
        Float_t ped = gRandom->Gaus(mean,sigma);
	Int_t pedADC = TMath::FloorNint(ped);
        d.SetElectronics(manuId, manuChannel);
        d.SetADC(pedADC);
        d.SetSignal(ped);
        d.SetDetElemId(detElemId); 
        // we do not set the remaining parts of the digit, as in principle
        // this is all we need : manuId, manuChannel and ADC, as far as
        // real data is concerned.
        new (digits[digits.GetLast()+1]) AliMUONDigit(d);
        ++fgCounter;
        ++ngenerated;
      }
    }
  }
  AliDebug(1,Form("ngenerated=%d",ngenerated));
}

//_____________________________________________________________________________
AliMUONData*
AliMUONPedestalEventGenerator::GetDataAccess(const char* mode)
{
  /// Get the pointer to AliMUONData object
  AliRunLoader* runLoader = LoadRun(mode);
  if (!runLoader)
  {
    AliError("Could not get RunLoader");
    return 0x0;
  }
  AliLoader* loader = static_cast<AliLoader*>(runLoader->GetLoader("MUONLoader"));
  if (!loader)
  {
    AliError("Could not get MuonLoader");
    return 0x0;
  }
  
  return new AliMUONData(loader,"MUON","MUONData");
}

//_____________________________________________________________________________
AliRunLoader*
AliMUONPedestalEventGenerator::LoadRun(const char* mode)
{
  /// Get access to AliRunLoader object
  while (AliRunLoader::GetRunLoader()) 
  {
    AliDebug(1,Form("Deleting AliRunLoader %p",AliRunLoader::GetRunLoader()));
    delete AliRunLoader::GetRunLoader();
  }
  
  AliRunLoader* runLoader = 
    AliRunLoader::Open(fGAliceFileName,AliConfig::GetDefaultEventFolderName(), 
                       mode);
  
  AliDebug(1,Form("AliRunLoader(%s)=%p",mode,runLoader));
           
  if (!runLoader) 
  {
    AliError("No run loader found in file galice.root");
  }
  return runLoader;
}
