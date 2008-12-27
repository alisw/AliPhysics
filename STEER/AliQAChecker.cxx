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

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliRunInfo.h" 
#include "AliLog.h"
#include "AliModule.h" 
#include "AliQA.h"
#include "AliQAChecker.h"
#include "AliQACheckerBase.h"
#include "AliCorrQAChecker.h"
#include "AliGlobalQAChecker.h"
#include "AliGRPObject.h"

#include <TKey.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TPluginManager.h> 
#include <TROOT.h>
#include <TStopwatch.h> 
#include <TString.h> 
#include <TSystem.h> 
#include <TList.h>
#include <TNtupleD.h>

ClassImp(AliQAChecker)
  AliQAChecker * AliQAChecker::fgQAChecker = 0x0 ;

//_____________________________________________________________________________
AliQAChecker::AliQAChecker(const char* name, const char* title) :
  TNamed(name, title),
  fDataFile(0x0), 
  fRunInfo(0x0), 
  fRunInfoOwner(kFALSE), 
  fRefFile(0x0), 
  fFoundDetectors("."), 
  fEventSpecie(AliRecoParam::kDefault) 
{
  // ctor: initialise checkers and open the data file   
  for (Int_t det = 0 ; det < AliQA::kNDET ; det++) 
    fCheckers[det] = NULL ; 
}

//_____________________________________________________________________________
AliQAChecker::AliQAChecker(const AliQAChecker& qac) :
  TNamed(qac),
  fDataFile(qac.fDataFile), 
  fRunInfo(qac.fRunInfo), 
  fRunInfoOwner(kFALSE),   
  fRefFile(qac.fRefFile), 
  fFoundDetectors(qac.fFoundDetectors),
  fEventSpecie(qac.fEventSpecie)
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
  if (fRunInfo)
    delete fRunInfo ; 
  delete [] fCheckers ; 
  AliQA::Close() ; 
}

//_____________________________________________________________________________
  AliQACheckerBase * AliQAChecker::GetDetQAChecker(Int_t det)
{
	// Gets the Quality Assurance checker for the detector specified by its name
	
	if (fCheckers[det]) 
    return fCheckers[det];

	AliQACheckerBase * qac = NULL ;

	TString detName(AliQA::GetDetName(det)) ; 
	
	if (det == AliQA::kGLOBAL) {
		qac = new AliGlobalQAChecker() ; 
	} else if (det == AliQA::kCORR) {
		qac = new AliCorrQAChecker() ; 
	} else {
		AliDebug(1, Form("Retrieving QA checker for %s", detName.Data())) ; 
		TPluginManager* pluginManager = gROOT->GetPluginManager() ;
		TString qacName = "Ali" + detName + "QAChecker" ;

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
  
		}
	}
	if (qac) 
		fCheckers[det] = qac ;
	
	return qac ; 
}
 
//_____________________________________________________________________________
void AliQAChecker::GetRefSubDir(const char * det, const char * task, TDirectory *& dirFile, TObjArray **& dirOCDB)     
{ 
  // Opens and returns the file with the reference data 
	
  dirFile = NULL ; 
  TString refStorage(AliQA::GetQARefStorage()) ; 
  if (refStorage.Contains(AliQA::GetLabLocalFile())) {	
    refStorage.ReplaceAll(AliQA::GetLabLocalFile(), "") ; 
    refStorage += AliQA::GetQARefFileName() ;
    if ( fRefFile ) 
      if ( fRefFile->IsOpen() ) 
					fRefFile->Close() ; 
    fRefFile = TFile::Open(refStorage.Data()) ; 
    if (!fRefFile) { 
      AliError(Form("Cannot find reference file %s", refStorage.Data())) ; 
      dirFile = NULL ; 
    }
    dirFile = fRefFile->GetDirectory(det) ; 
    if (!dirFile) {
      AliWarning(Form("Directory %s not found in %d", det, refStorage.Data())) ; 
    } else {
			dirFile = dirFile->GetDirectory(task) ; 
      if (!dirFile) 
				AliWarning(Form("Directory %s/%s not found in %s", det, task, refStorage.Data())) ; 
    }  
  } else if (refStorage.Contains(AliQA::GetLabLocalOCDB()) || refStorage.Contains(AliQA::GetLabAliEnOCDB())) {	
    AliCDBManager* man = AliCDBManager::Instance() ;
    for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
      if ( !AliQA::Instance()->IsEventSpecieSet(specie) ) 
        continue ; 
      //if ( strcmp(AliQA::GetRefDataDirName(), "") == 0 ) { // the name of the last level of the directory is not set (EventSpecie)
        // Get it from RunInfo
        //if (!fRunInfo)  // not yet set, get the info from GRP
        //  LoadRunInfoFromGRP() ; 
      AliQA::SetQARefDataDirName(specie) ;
      //}
      if ( ! man->GetLock() ) { 
        man->SetDefaultStorage(AliQA::GetQARefStorage()) ; 
        man->SetSpecificStorage("*", AliQA::GetQARefStorage()) ;
      }
      char * detOCDBDir = Form("%s/%s/%s", det, AliQA::GetRefOCDBDirName(), AliQA::GetRefDataDirName()) ; 
      AliCDBEntry * entry = man->Get(detOCDBDir, man->GetRun()) ;
      if (entry) {
        dirOCDB = new TObjArray*[AliRecoParam::kNSpecies] ;	
        TList * listDetQAD = dynamic_cast<TList *>(entry->GetObject()) ;
        TIter next(listDetQAD) ;
        TObjArray * ar ; 
        while ( ar = (TObjArray*)next() )
          if ( listDetQAD ) 
          dirOCDB[specie] = dynamic_cast<TObjArray *>(listDetQAD->FindObject(Form("%s/%s", task, AliRecoParam::GetEventSpecieName(specie)))) ; 
      }
    }
  }
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
void AliQAChecker::LoadRunInfoFromGRP() 
{
  AliCDBManager* man = AliCDBManager::Instance() ;
  AliCDBEntry* entry = man->Get(AliQA::GetGRPPath().Data());
  AliGRPObject* grpObject = 0x0;
  if (entry) {

	  TMap* m = dynamic_cast<TMap*>(entry->GetObject());  // old GRP entry

	  if (m) {
	    AliInfo("It is a map");
	    //m->Print();
	    grpObject = new AliGRPObject();
	         grpObject->ReadValuesFromMap(m);
    }

    else {
	    AliInfo("It is a new GRP object");
        grpObject = dynamic_cast<AliGRPObject*>(entry->GetObject());  // new GRP entry
    }

    entry->SetOwner(0);
    AliCDBManager::Instance()->UnloadFromCache("GRP/GRP/Data");
  }

  if (!grpObject) {
     AliFatal("No GRP entry found in OCDB!");
  }

  TString lhcState = grpObject->GetLHCState();
  if (lhcState==AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the LHC state ! Using UNKNOWN");
    lhcState = "UNKNOWN";
  }

  TString beamType = grpObject->GetBeamType();
  if (beamType==AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the beam type ! Using UNKNOWN");
    beamType = "UNKNOWN";
  }

  Float_t beamEnergy = grpObject->GetBeamEnergy();
  if (beamEnergy==AliGRPObject::GetInvalidFloat()) {
    AliError("GRP/GRP/Data entry:  missing value for the beam energy ! Using 0");
    beamEnergy = 0;
  }

  TString runType = grpObject->GetRunType();
  if (runType==AliGRPObject::GetInvalidString()) {
    AliError("GRP/GRP/Data entry:  missing value for the run type ! Using UNKNOWN");
    runType = "UNKNOWN";
  }

  Int_t activeDetectors = grpObject->GetDetectorMask();
  if (activeDetectors==AliGRPObject::GetInvalidInt()) {
    AliError("GRP/GRP/Data entry:  missing value for the detector mask ! Using 1074790399");
    activeDetectors = 1074790399;
  }

  fRunInfo = new AliRunInfo(lhcState, beamType, beamEnergy, runType, activeDetectors);

  fRunInfoOwner = kTRUE ; 

  // set the event specie
  fEventSpecie = AliRecoParam::kDefault ;
  if (strcmp(runType,"PHYSICS")) {
    // Not a physics run, the event specie is set to kCalib
    fEventSpecie = AliRecoParam::kCalib ;
    return;
  }
  if (strcmp(lhcState,"STABLE_BEAMS") == 0) {
    // Heavy ion run (any beam tha is not pp, the event specie is set to kHighMult
    fEventSpecie = AliRecoParam::kHighMult ;
    if ((strcmp(beamType,"p-p") == 0) ||
        (strcmp(beamType,"p-")  == 0) ||
        (strcmp(beamType,"-p")  == 0) ||
        (strcmp(beamType,"P-P") == 0) ||
        (strcmp(beamType,"P-")  == 0) ||
        (strcmp(beamType,"-P")  == 0)) {
      // Proton run, the event specie is set to kLowMult
      fEventSpecie = AliRecoParam::kLowMult ;
    }
    else if (strcmp(beamType,"-") == 0) {
      // No beams, we assume cosmic data
      fEventSpecie = AliRecoParam::kCosmic ;
    }
    else if (strcmp(beamType,"UNKNOWN") == 0) {
      // No LHC beam information is available, we use the default event specie
      fEventSpecie = AliRecoParam::kDefault ;
    }
  }
}

//_____________________________________________________________________________
Bool_t AliQAChecker::Run(const char * fileName)
{
  // run the Quality Assurance Checker for all tasks Hits, SDigits, Digits, recpoints, tracksegments, recparticles and ESDs
  // starting from data in file  
  TStopwatch stopwatch;
  stopwatch.Start();

  //search for all detectors QA directories
  TList * detKeyList = AliQA::GetQADataFile(fileName)->GetListOfKeys() ; 
  TIter nextd(detKeyList) ; 
  TKey * detKey ; 
  while ( (detKey = dynamic_cast<TKey *>(nextd()) ) ) {
    AliDebug(1, Form("Found %s", detKey->GetName())) ;
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
    TDirectory * detDir = AliQA::GetQADataFile(fileName)->GetDirectory(detKey->GetName()) ; 
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
      AliQA::ALITASK_t index = AliQA::kNULLTASK ; 
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
      qac->Init(AliQA::DETECTORINDEX_t(det)) ; 
      
      TDirectory * refDir     = NULL ; 
      TObjArray ** refOCDBDir = NULL ;	
      GetRefSubDir(detNameQA.Data(), taskName.Data(), refDir, refOCDBDir) ;
		  qac->SetRefandData(refDir, refOCDBDir, taskDir) ;
		  qac->Run(index) ; 
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
  return kTRUE ; 
}

//_____________________________________________________________________________
Bool_t AliQAChecker::Run(AliQA::DETECTORINDEX_t det, AliQA::TASKINDEX_t task, TObjArray ** list)
{
	// run the Quality Assurance Checker for detector det, for task task starting from data in list

	AliQACheckerBase * qac = GetDetQAChecker(det) ; 
	if (qac)
		AliDebug(1, Form("QA checker found for %s", AliQA::GetDetName(det).Data())) ;
	if (!qac)
		AliError(Form("QA checker not found for %s", AliQA::GetDetName(det).Data())) ; 
  
	AliQA::ALITASK_t index = AliQA::kNULLTASK ; 
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

	TDirectory * refDir     = NULL ; 
	TObjArray ** refOCDBDir = NULL  ;	
  qac->Init(det) ; 
  GetRefSubDir(AliQA::GetDetName(det), AliQA::GetTaskName(task), refDir, refOCDBDir) ;
  qac->SetRefandData(refDir, refOCDBDir) ; 
  qac->Run(index, list) ; 
	return kTRUE ; 
}

//_____________________________________________________________________________
Bool_t AliQAChecker::Run(AliQA::DETECTORINDEX_t det, AliQA::TASKINDEX_t task, TNtupleD ** list)
{
	// run the Quality Assurance Checker for detector det, for task task starting from data in list
  
	AliQACheckerBase * qac = GetDetQAChecker(det) ; 
	if (qac)
		AliDebug(1, Form("QA checker found for %s", AliQA::GetDetName(det).Data())) ;
	if (!qac)
		AliError(Form("QA checker not found for %s", AliQA::GetDetName(det).Data())) ; 
  
	AliQA::ALITASK_t index = AliQA::kNULLTASK ; 
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
  
	TDirectory * refDir     = NULL ; 
	TObjArray ** refOCDBDir = NULL ;	
  qac->Init(det) ; 
  GetRefSubDir(AliQA::GetDetName(det), AliQA::GetTaskName(task), refDir, refOCDBDir) ;
  qac->SetRefandData(refDir, refOCDBDir) ; 
  qac->Run(index, list) ; 
  return kTRUE ; 
}
