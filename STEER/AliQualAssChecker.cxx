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
#include "AliQualAss.h"
#include "AliQualAssChecker.h"
#include "AliQualAssCheckerBase.h"

#include <TKey.h>
#include <TObjArray.h>
#include <TPluginManager.h> 
#include <TROOT.h>
#include <TStopwatch.h> 
#include <TString.h> 
#include <TSystem.h> 

ClassImp(AliQualAssChecker)
  TFile   * AliQualAssChecker::fgQAResultFile = 0x0 ;  
  TString   AliQualAssChecker::fgQAResultDirName = "local://RUN/";  
  TString   AliQualAssChecker::fgQAResultFileName = "QA.root" ; 

//_____________________________________________________________________________
AliQualAssChecker::AliQualAssChecker(const char* name, const char* title) :
  TNamed(name, title),
  fDataFile(0x0), 
  fRefDirName("/QA/Ref/"), 
  fRefName("QA.root"), 
  fFoundDetectors(".")
{
  // ctor: initialise checkers and open the data file   
  for (Int_t det = 0 ; det < AliQualAss::kNDET ; det++) 
    fCheckers[det] = NULL ; 
  
  GetDataFile() ; 
}

//_____________________________________________________________________________
AliQualAssChecker::AliQualAssChecker(const AliQualAssChecker& qac) :
  TNamed(qac),
  fDataFile(qac.fDataFile), 
  fRefDirName(qac.fRefDirName), 
  fRefName(qac.fRefName), 
  fFoundDetectors(qac.fFoundDetectors)
{
  // copy constructor
  
  for (Int_t det = 0 ; det < AliQualAss::kNDET ; det++) 
    fCheckers[det] = NULL ; 
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
}

//_____________________________________________________________________________
TFile * AliQualAssChecker:: GetDataFile()
{
  // Open if necessary the Data file and return its pointer

  if (!fDataFile) 
    fDataFile =  TFile::Open(AliQualAss::GetDataName()) ;
  if (!fDataFile) 
    AliFatal(Form("QA Data File %s does not exist", AliQualAss::GetDataName() )) ; 
  return fDataFile ; 
}

//_____________________________________________________________________________
TFile * AliQualAssChecker:: GetQAResultFile() 
{
  // Check if file to store QA exists, if not create it

  if (fgQAResultFile) { 
    if (fgQAResultFile->IsOpen()){
      fgQAResultFile->Close() ; 
      fgQAResultFile = 0x0 ; 
    }
  }   
  if ( fgQAResultFileName.Contains("local://")) 
    fgQAResultFileName.ReplaceAll("local:/", "") ;
  
  TString opt("") ; 
  if ( !gSystem->AccessPathName(fgQAResultFileName) )
    opt = "UPDATE" ; 
  else 
    opt = "NEW" ; 
  fgQAResultFile = TFile::Open(fgQAResultFileName, opt) ;   
      
  return fgQAResultFile ; 
}

//_____________________________________________________________________________
  AliQualAssCheckerBase * AliQualAssChecker::GetDetQualAssChecker(Int_t det)
{
  // Gets the Quality Assurance checker for the detector specified by its name
  
  if (fCheckers[det]) 
    return fCheckers[det];

 TString detName(AliQualAss::GetDetName(det)) ; 

  AliInfo(Form("Retrieving QA checker for %s", detName.Data())) ; 
  TPluginManager* pluginManager = gROOT->GetPluginManager() ;
  TString qacName = "Ali" + detName + "QualAssChecker" ;

  AliQualAssCheckerBase * qac = NULL ;
  // first check if a plugin is defined for the quality assurance checker
  TPluginHandler* pluginHandler = pluginManager->FindHandler("AliQualAssChecker", detName.Data());
  // if not, add a plugin for it
  if (!pluginHandler) {
    //AliInfo(Form("defining plugin for %s", qacName.Data()));
    TString libs = gSystem->GetLibraries();
 
   if (libs.Contains("lib" + detName + "base.so") || (gSystem->Load("lib" + detName + "base.so") >= 0))
      pluginManager->AddHandler("AliQualAssChecker", detName, qacName, detName + "qac", qacName + "()");
    else 
      pluginManager->AddHandler("AliQualAssChecker", detName, qacName, detName, qacName + "()");

   pluginHandler = pluginManager->FindHandler("AliQualAssChecker", detName);

  if (pluginHandler && (pluginHandler->LoadPlugin() == 0)) 
    qac = (AliQualAssCheckerBase *) pluginHandler->ExecPlugin(0);
  
  if (qac) 
    fCheckers[det] = qac ; 
  }

 return qac ; 
}


//_____________________________________________________________________________
TDirectory * AliQualAssChecker::GetRefSubDir(const char * det, const char * task)     
{ 
  // Opens and returns the file with the reference data 
  TFile * f = TFile::Open(fRefDirName, "READ") ;
  if (!f) 
    AliFatal(Form("Cannot find reference file %s", fRefDirName.Data())) ; 
  TDirectory * rv = NULL ; 
  rv = f->GetDirectory(det) ; 
  if (!rv) {
    AliWarning(Form("Directory %s not found in %d", det, fRefDirName.Data())) ; 
  } else {
    rv = rv->GetDirectory(task) ; 
    if (!rv) 
      AliWarning(Form("Directory %s/%s not found in %s", det, task, fRefDirName.Data())) ; 
  }  
  return rv ; 
}

//_____________________________________________________________________________
Bool_t AliQualAssChecker::Run()
{
  // run the Quality Assurance Checker for all tasks Hits, SDigits, Digits, recpoints, tracksegments, recparticles and ESDs

  Bool_t rv = kFALSE ; 
  
  TStopwatch stopwatch;
  stopwatch.Start();

  //search for all detectors QA directories
  TList * detKeyList = GetDataFile()->GetListOfKeys() ; 
  TIter nextd(detKeyList) ; 
  TKey * detKey ; 
  while ( (detKey = dynamic_cast<TKey *>(nextd()) ) ) {
    AliInfo(Form("Found %s", detKey->GetName())) ;
    //Check which detector
    TString detName ; 
    TString detNameQA(detKey->GetName()) ; 
    Int_t det ; 
    for ( det = 0; det < AliQualAss::kNDET ; det++) {
      detName = AliQualAss::GetDetName(det) ; 
      if (detNameQA.Contains(detName)) {
	fFoundDetectors+=detName ; 
	fFoundDetectors+="." ;		
	break ; 
      }
    } 
    TDirectory * detDir = GetDataFile()->GetDirectory(detKey->GetName()) ; 
    TList * taskKeyList = detDir->GetListOfKeys() ;
    TIter nextt(taskKeyList) ; 
    TKey * taskKey ; 
    // now search for the tasks dir
    while ( (taskKey = static_cast<TKey *>(nextt()) ) ) {
      TString taskName( taskKey->GetName() ) ; 
      AliInfo(Form("Found %s", taskName.Data())) ;
      TDirectory * taskDir = detDir->GetDirectory(taskName.Data()) ; 
      taskDir->cd() ; 
      AliQualAssCheckerBase * qac = GetDetQualAssChecker(det) ; 
      if (qac)
 	AliInfo(Form("QA checker found for %s", detName.Data())) ; 
      if (!qac)
 	AliFatal(Form("QA checker not found for %s", detName.Data())) ; 
      AliQualAss::ALITASK index = AliQualAss::kNULLTASK ; 
      if ( taskName == AliQualAss::GetTaskName(AliQualAss::kHITS) ) 
	index = AliQualAss::kSIM ; 
      if ( taskName == AliQualAss::GetTaskName(AliQualAss::kSDIGITS) ) 
	index = AliQualAss::kSIM ; 
      if ( taskName == AliQualAss::GetTaskName(AliQualAss::kDIGITS) ) 
	index = AliQualAss::kSIM ; 
      if ( taskName == AliQualAss::GetTaskName(AliQualAss::kRECPOINTS) ) 
	index = AliQualAss::kREC ; 
      if ( taskName == AliQualAss::GetTaskName(AliQualAss::kTRACKSEGMENTS) ) 
	index = AliQualAss::kREC ; 
      if ( taskName == AliQualAss::GetTaskName(AliQualAss::kRECPARTICLES) ) 
	index = AliQualAss::kREC ; 
      if ( taskName == AliQualAss::GetTaskName(AliQualAss::kESDS) ) 
	index = AliQualAss::kESD ; 
      qac->Init(AliQualAss::DETECTORINDEX(det)) ; 
      qac->SetRefandData(GetRefSubDir(detNameQA.Data(), taskName.Data()), taskDir) ; 
      qac->Run(index) ; 
    }
 }
  AliInfo("QA performed for following detectors:") ; 
  for ( Int_t det = 0; det < AliQualAss::kNDET; det++) {
    if (fFoundDetectors.Contains(AliQualAss::GetDetName(det))) {
      printf("%s, ",AliQualAss::GetDetName(det)) ; 
      fFoundDetectors.ReplaceAll(AliQualAss::GetDetName(det), "") ; 
    }	
  }
  printf("\n") ; 
  rv = kTRUE ; 

  return rv ; 
  
}

//_____________________________________________________________________________
void AliQualAssChecker::SetQAResultDirName(const char * name)
{
  // Set the root directory where to store the QA status object

  fgQAResultDirName.Prepend(name) ; 
  AliInfo(Form("QA results are in  %s", fgQAResultDirName.Data())) ;
  if ( fgQAResultDirName.Contains("local://")) 
    fgQAResultDirName.ReplaceAll("local:/", "") ;
  fgQAResultFileName.Prepend(fgQAResultDirName) ;
}

//_____________________________________________________________________________
void AliQualAssChecker::SetRefDirName(const char * name)
{
  // Set the root directory of reference data

  fRefDirName.Prepend(name) ; 
  fRefDirName.Append(fRefName) ; 
  AliInfo(Form("Reference data are taken from %s", fRefDirName.Data())) ;
  if ( fRefDirName.Contains("local://")) 
    fRefDirName.ReplaceAll("local:/", "") ; 
}





