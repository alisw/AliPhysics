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

/* $Id: */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for running the Quality Assurance Checker
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliLog.h"
#include "AliModule.h" 
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliQualAss.h"
#include "AliQualAssChecker.h"

#include <TObjArray.h>
#include <TStopwatch.h> 
#include <TString.h> 


ClassImp(AliQualAssChecker)

  TFile * AliQualAssChecker::fgOutFile = 0x0 ; 
  TString AliQualAssChecker::fgOutDir("/RUN/") ; 
  TString AliQualAssChecker::fgOutName("QA.root") ; 
  TString AliQualAssChecker::fgRefDir("/QA/Ref/") ; 
  TString AliQualAssChecker::fgRefName("QA.root") ; 

//_____________________________________________________________________________
AliQualAssChecker::AliQualAssChecker(const char* name, const char* title) :
  TNamed(name, title),
  fGAliceFileName("galice.root"),
  fStopOnError(kFALSE)
{
  // create simulation object with default parameters

  SetGAliceFile("galice.root");
}

//_____________________________________________________________________________
AliQualAssChecker::AliQualAssChecker(const AliQualAssChecker& qac) :
  TNamed(qac),
  fGAliceFileName(qac.fGAliceFileName),
  fStopOnError(qac.fStopOnError)
{
// copy constructor
}

//_____________________________________________________________________________
AliQualAssChecker& AliQualAssChecker::operator = (const AliQualAssChecker& qac)
{
// assignment operator

  this->~AliQualAssChecker();
  new(this) AliQualAssChecker(qac);
  return *this;
}

//_____________________________________________________________________________
AliQualAssChecker::~AliQualAssChecker()
{
// clean up
  fgOutFile->Close() ; 
  delete fgOutFile ;
}

//_____________________________________________________________________________
TFile * AliQualAssChecker:: GetOutFile() 
{
  // Check if file to store QA exists, if not create it

  if (fgOutFile) { 
    if (fgOutFile->IsOpen()){
      fgOutFile->Close() ; 
      fgOutFile = 0x0 ; 
    }
  }   
  fgOutName.Prepend(fgOutDir) ;
  TString opt("") ; 
  if ( !gSystem->AccessPathName(fgOutName) )
    opt = "UPDATE" ; 
  else 
    opt = "NEW" ; 
  fgOutFile = TFile::Open(fgOutName, opt) ;   
  
  return fgOutFile ; 
}

//_____________________________________________________________________________
AliRunLoader* AliQualAssChecker::LoadRun(const char* mode) const
{
// delete existing run loaders, open a new one and load gAlice
  while (AliRunLoader::GetRunLoader()) 
    delete AliRunLoader::GetRunLoader();
  AliRunLoader* runLoader = 
    AliRunLoader::Open(fGAliceFileName.Data(), 
		       AliConfig::GetDefaultEventFolderName(), mode);
  if (!runLoader) {
    AliError(Form("no run loader found in file %s", fGAliceFileName.Data()));
    return NULL;
  }
  runLoader->LoadgAlice();
  gAlice = runLoader->GetAliRun();
  if (!gAlice) {
    AliError(Form("no gAlice object found in file %s", 
                  fGAliceFileName.Data()));
    return NULL;
  }
  return runLoader;
}

//_____________________________________________________________________________
Bool_t AliQualAssChecker::Run()
{
 // run the Quality Assurance Checker for all tasks Hits, SDigits, Digits, recpoints, tracksegments, recparticles and ESDs

  Bool_t rv = kFALSE ; 

  TStopwatch stopwatch;
  stopwatch.Start();

  AliRunLoader* runLoader = LoadRun("READ");
  if (!runLoader) 
    return rv ;

  TObjArray* detArray = runLoader->GetAliRun()->Detectors();
  for (Int_t iDet = 0; iDet < detArray->GetEntriesFast(); iDet++) {
    AliModule* det = (AliModule*) detArray->At(iDet);
    if (!det || !det->IsActive()) 
      continue;
      AliInfo(Form("QA checking %s", det->GetName()));
      det->CheckQA();
    }

  delete runLoader;

  AliInfo(Form("Execution time for QA: R:%.2fs C:%.2fs", stopwatch.RealTime(),stopwatch.CpuTime()));

  return rv;
}

//_____________________________________________________________________________
void AliQualAssChecker::SetGAliceFile(const char* fileName)
{
// set the name of the galice file
// the path is converted to an absolute one if it is relative

  fGAliceFileName = fileName;
  if (!gSystem->IsAbsoluteFileName(fGAliceFileName)) {
    char* absFileName = gSystem->ConcatFileName(gSystem->WorkingDirectory(),
						fGAliceFileName);
    fGAliceFileName = absFileName;
    delete[] absFileName;
  }

  AliDebug(2, Form("galice file name set to %s", fileName));
}

//_____________________________________________________________________________
void AliQualAssChecker::SetOutDir(const char * outDir)
{
  // Set the root directory where to store the QA status object

  fgOutDir.Prepend(outDir) ; 
  AliInfo(Form("QA results are in  %s", fgOutDir.Data())) ;
  if ( fgOutDir.Contains("local://")) 
    fgOutDir.ReplaceAll("local:/", "") ;  
}

//_____________________________________________________________________________
void AliQualAssChecker::SetRefDir(const char * refDir)
{
  // Set the root directory of reference data

  fgRefDir.Prepend(refDir) ; 
  fgRefDir.Append(fgRefName) ; 
  AliInfo(Form("Reference data are taken from %s", fgRefDir.Data())) ;
  if ( fgRefDir.Contains("local://")) 
    fgRefDir.ReplaceAll("local:/", "") ; 
}





