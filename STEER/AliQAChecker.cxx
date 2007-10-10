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
#include "AliQA.h"
#include "AliQAChecker.h"
#include "AliQACheckerBase.h"

#include <TKey.h>
#include <TObjArray.h>
#include <TPluginManager.h> 
#include <TROOT.h>
#include <TStopwatch.h> 
#include <TString.h> 
#include <TSystem.h> 
#include <TList.h>

ClassImp(AliQAChecker)
  AliQAChecker * AliQAChecker::fgQAChecker = 0x0 ;
  TFile   * AliQAChecker::fgQAResultFile        = 0x0 ;  
  TString   AliQAChecker::fgQAResultDirName     = "local://RUN/";  
  TString   AliQAChecker::fgQAResultFileName    = "QA.root" ; 

//_____________________________________________________________________________
AliQAChecker::AliQAChecker(const char* name, const char* title) :
  TNamed(name, title),
  fDataFile(0x0), 
  fRefDirName("./Ref/"), 
  fRefName("QA.root"), 
  fFoundDetectors(".")
{
  // ctor: initialise checkers and open the data file   
  for (Int_t det = 0 ; det < AliQA::kNDET ; det++) 
    fCheckers[det] = NULL ; 
  
  GetDataFile() ; 
  fRefDirName.Append(fRefName) ; 
}

//_____________________________________________________________________________
AliQAChecker::AliQAChecker(const AliQAChecker& qac) :
  TNamed(qac),
  fDataFile(qac.fDataFile), 
  fRefDirName(qac.fRefDirName), 
  fRefName(qac.fRefName), 
  fFoundDetectors(qac.fFoundDetectors)
{
  // copy constructor
  
  for (Int_t det = 0 ; det < AliQA::kNDET ; det++) 
    fCheckers[det] = NULL ; 
}

//_____________________________________________________________________________
AliQAChecker& AliQAChecker::operator = (const AliQAChecker& qac)
{
// assignment operator

  this->~AliQAChecker();
  new(this) AliQAChecker(qac);
  return *this;
}

//_____________________________________________________________________________
AliQAChecker::~AliQAChecker()
{
// clean up
  delete [] fCheckers ; 
  fgQAResultFile->Close() ; 
}

//_____________________________________________________________________________
TFile * AliQAChecker:: GetDataFile()
{
  // Open if necessary the Data file and return its pointer

  if (!fDataFile) 
	if  (gSystem->AccessPathName(AliQA::GetDataName()))
     fDataFile =  TFile::Open(AliQA::GetDataName()) ;
  return fDataFile ; 
}

//_____________________________________________________________________________
TFile * AliQAChecker:: GetQAResultFile() 
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
  AliQACheckerBase * AliQAChecker::GetDetQAChecker(Int_t det)
{
  // Gets the Quality Assurance checker for the detector specified by its name
  
  if (fCheckers[det]) 
    return fCheckers[det];

 TString detName(AliQA::GetDetName(det)) ; 

  AliInfo(Form("Retrieving QA checker for %s", detName.Data())) ; 
  TPluginManager* pluginManager = gROOT->GetPluginManager() ;
  TString qacName = "Ali" + detName + "QAChecker" ;

  AliQACheckerBase * qac = NULL ;
  // first check if a plugin is defined for the quality assurance checker
  TPluginHandler* pluginHandler = pluginManager->FindHandler("AliQAChecker", detName.Data());
  // if not, add a plugin for it
  if (!pluginHandler) {
    //AliInfo(Form("defining plugin for %s", qacName.Data()));
    TString libs = gSystem->GetLibraries();
 
   if (libs.Contains("lib" + detName + "base.so") || (gSystem->Load("lib" + detName + "base.so") >= 0))
      pluginManager->AddHandler("AliQAChecker", detName, qacName, detName + "qac", qacName + "()");
    else 
      pluginManager->AddHandler("AliQAChecker", detName, qacName, detName, qacName + "()");

   pluginHandler = pluginManager->FindHandler("AliQAChecker", detName);

  if (pluginHandler && (pluginHandler->LoadPlugin() == 0)) 
    qac = (AliQACheckerBase *) pluginHandler->ExecPlugin(0);
  
  if (qac) 
    fCheckers[det] = qac ; 
  }

 return qac ; 
}


//_____________________________________________________________________________
TDirectory * AliQAChecker::GetRefSubDir(const char * det, const char * task)     
{ 
  // Opens and returns the file with the reference data 
  TFile * f = TFile::Open(fRefDirName, "READ") ;
  if (!f) 
    AliError(Form("Cannot find reference file %s", fRefDirName.Data())) ; 
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
AliQAChecker * AliQAChecker::Instance()
{
	// returns unique instance of the checker
  if ( ! fgQAChecker ) 
   fgQAChecker = new AliQAChecker() ; 
  return fgQAChecker ;  
}

//_____________________________________________________________________________
Bool_t AliQAChecker::Run()
{
  // run the Quality Assurance Checker for all tasks Hits, SDigits, Digits, recpoints, tracksegments, recparticles and ESDs
  // starting from data in file

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
    for ( det = 0; det < AliQA::kNDET ; det++) {
      detName = AliQA::GetDetName(det) ; 
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
      AliQACheckerBase * qac = GetDetQAChecker(det) ; 
      if (qac)
		AliInfo(Form("QA checker found for %s", detName.Data())) ; 
      if (!qac)
		AliFatal(Form("QA checker not found for %s", detName.Data())) ; 
      AliQA::ALITASK index = AliQA::kNULLTASK ; 
      if ( taskName == AliQA::GetTaskName(AliQA::kHITS) ) 
		index = AliQA::kSIM ; 
      if ( taskName == AliQA::GetTaskName(AliQA::kSDIGITS) ) 
		index = AliQA::kSIM ; 
      if ( taskName == AliQA::GetTaskName(AliQA::kDIGITS) ) 
		index = AliQA::kSIM ; 
      if ( taskName == AliQA::GetTaskName(AliQA::kRECPOINTS) ) 
		index = AliQA::kREC ; 
      if ( taskName == AliQA::GetTaskName(AliQA::kTRACKSEGMENTS) ) 
		index = AliQA::kREC ; 
      if ( taskName == AliQA::GetTaskName(AliQA::kRECPARTICLES) ) 
		index = AliQA::kREC ; 
      if ( taskName == AliQA::GetTaskName(AliQA::kESDS) ) 
		index = AliQA::kESD ; 
      qac->Init(AliQA::DETECTORINDEX(det)) ; 
	  TDirectory * refDir = GetRefSubDir(detNameQA.Data(), taskName.Data()) ;
	  if ( refDir ) { 
		qac->SetRefandData(refDir, taskDir) ; 
		qac->Run(index) ; 
	  }
    }
 }
  AliInfo("QA performed for following detectors:") ; 
  for ( Int_t det = 0; det < AliQA::kNDET; det++) {
    if (fFoundDetectors.Contains(AliQA::GetDetName(det))) {
      printf("%s, ",AliQA::GetDetName(det)) ; 
      fFoundDetectors.ReplaceAll(AliQA::GetDetName(det), "") ; 
    }	
  }
  printf("\n") ; 
  rv = kTRUE ; 

  return rv ; 
  
}

//_____________________________________________________________________________
Bool_t AliQAChecker::Run(AliQA::DETECTORINDEX det, AliQA::TASKINDEX task, TList * list)
{
  // run the Quality Assurance Checker for detector det, for task task starting from data in list

  AliQACheckerBase * qac = GetDetQAChecker(det) ; 
  if (qac)
    AliInfo(Form("QA checker found for %s", AliQA::GetDetName(det).Data())) ;
  if (!qac)
	AliFatal(Form("QA checker not found for %s", AliQA::GetDetName(det).Data())) ; 
  
  AliQA::ALITASK index = AliQA::kNULLTASK ; 
  if ( task == AliQA::kRAWS ) 
		index = AliQA::kRAW ; 
  else if ( task == AliQA::kHITS ) 
		index = AliQA::kSIM ; 
  else if ( task == AliQA::kSDIGITS ) 
		index = AliQA::kSIM ; 
  else if ( task == AliQA::kDIGITS ) 
		index = AliQA::kSIM ; 
  else if ( task == AliQA::kRECPOINTS ) 
		index = AliQA::kREC ; 
  else if ( task == AliQA::kTRACKSEGMENTS ) 
		index = AliQA::kREC ; 
  else if ( task == AliQA::kRECPARTICLES ) 
		index = AliQA::kREC ; 
  else if ( task == AliQA::kESDS ) 
		index = AliQA::kESD ; 

  qac->Init(det) ; 
  qac->SetRefandData(GetRefSubDir(AliQA::GetDetName(det).Data(), AliQA::GetTaskName(task).Data())) ; 
  qac->Run(index, list) ; 
  
  return kTRUE ; 
  
}

//_____________________________________________________________________________
void AliQAChecker::SetQAResultDirName(const char * name)
{
  // Set the root directory where to store the QA status object

  fgQAResultDirName.Prepend(name) ; 
  AliInfo(Form("QA results are in  %s", fgQAResultDirName.Data())) ;
  if ( fgQAResultDirName.Contains("local://")) 
    fgQAResultDirName.ReplaceAll("local:/", "") ;
  fgQAResultFileName.Prepend(fgQAResultDirName) ;
}

//_____________________________________________________________________________
void AliQAChecker::SetRefDirName(const char * name)
{
  // Set the root directory of reference data

  fRefDirName.Prepend(name) ; 
  fRefDirName.Append(fRefName) ; 
  AliInfo(Form("Reference data are taken from %s", fRefDirName.Data())) ;
  if ( fRefDirName.Contains("local://")) 
    fRefDirName.ReplaceAll("local:/", "") ; 
}





