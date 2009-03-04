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

/* $Id: AliQAManager.cxx 30894 2009-02-05 13:46:48Z schutz $ */
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for running the QA makers                                           //
//                                                                           //
//   AliQAManager qas;                                                //
//   qas.Run(AliQA::kRAWS, rawROOTFileName);                                 //
//   qas.Run(AliQA::kHITS);                                                  //
//   qas.Run(AliQA::kSDIGITS);                                               //
//   qas.Run(AliQA::kDIGITS);                                                //
//   qas.Run(AliQA::kRECPOINTS);                                             //
//   qas.Run(AliQA::kESDS);                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TKey.h>
#include <TFile.h>
#include <TFileMerger.h>
#include <TPluginManager.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliCDBId.h"
#include "AliCDBMetaData.h"
#include "AliCodeTimer.h"
#include "AliCorrQADataMakerRec.h"
#include "AliDetectorRecoParam.h"
#include "AliESDEvent.h"
#include "AliGeomManager.h"
#include "AliGlobalQADataMaker.h"
#include "AliHeader.h"
#include "AliLog.h"
#include "AliModule.h"
#include "AliQA.h"
#include "AliQADataMakerRec.h"
#include "AliQADataMakerSim.h"
#include "AliQAManager.h" 
#include "AliRawReaderDate.h"
#include "AliRawReaderFile.h"
#include "AliRawReaderRoot.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliRunTag.h"

ClassImp(AliQAManager) 
AliQAManager* AliQAManager::fgQAInstance = 0x0;

//_____________________________________________________________________________
AliQAManager::AliQAManager() :
  AliCDBManager(), 
  fCurrentEvent(0),   
  fCycleSame(kFALSE),
  fDetectors("ALL"), 
  fDetectorsW("ALL"), 
  fESD(NULL), 
  fESDTree(NULL),
  fGAliceFileName(""), 
  fFirstEvent(0),        
  fMaxEvents(0),   
  fMode(""), 
  fNumberOfEvents(999999), 
  fRecoParam(),
  fRunNumber(0), 
  fRawReader(NULL), 
  fRawReaderDelete(kTRUE), 
  fRunLoader(NULL), 
  fTasks(""), 
  fEventSpecie(AliRecoParam::kDefault)
{
	// default ctor
	fMaxEvents = fNumberOfEvents ; 
	for (UInt_t iDet = 0; iDet < fgkNDetectors; iDet++) {
		if (IsSelected(AliQA::GetDetName(iDet))) {
			fLoader[iDet]      = NULL ;
			fQADataMaker[iDet] = NULL ;
			fQACycles[iDet]    = 999999 ;
      fQAWriteExpert[iDet] = kTRUE ;
		}
	}	
}

//_____________________________________________________________________________
AliQAManager::AliQAManager(const Char_t * mode, const Char_t* gAliceFilename) :
  AliCDBManager(), 
  fCurrentEvent(0),  
	fCycleSame(kFALSE),
	fDetectors("ALL"), 
	fDetectorsW("ALL"), 
	fESD(NULL), 
	fESDTree(NULL),
	fGAliceFileName(gAliceFilename), 
	fFirstEvent(0),        
	fMaxEvents(0),   
  fMode(mode), 
	fNumberOfEvents(999999), 
  fRecoParam(),
	fRunNumber(0), 
	fRawReader(NULL), 
	fRawReaderDelete(kTRUE), 
	fRunLoader(NULL), 
  fTasks(""), 
  fEventSpecie(AliRecoParam::kDefault)
{
	// default ctor
	fMaxEvents = fNumberOfEvents ; 
	for (UInt_t iDet = 0; iDet < fgkNDetectors; iDet++) {
		if (IsSelected(AliQA::GetDetName(iDet))) {
			fLoader[iDet]      = NULL ;
			fQADataMaker[iDet] = NULL ;
			fQACycles[iDet]    = 999999 ;
      fQAWriteExpert[iDet] = kTRUE ;
		}
	}	  
  // set the default storage
  TString sto(AliQA::GetQARefStorage()) ;
    if (sto.IsNull()) 
      AliFatal("QA reference storage not set do: AliQA::SetQARefStorage") ; 
  SetDefaultStorage(sto.Data());
}

//_____________________________________________________________________________
AliQAManager::AliQAManager(const AliQAManager & qas) : 
  AliCDBManager(), 
	fCurrentEvent(qas.fCurrentEvent),  
	fCycleSame(kFALSE),
	fDetectors(qas.fDetectors), 
	fDetectorsW(qas.fDetectorsW), 
	fESD(NULL), 
	fESDTree(NULL), 
	fGAliceFileName(qas.fGAliceFileName), 
	fFirstEvent(qas.fFirstEvent),        
	fMaxEvents(qas.fMaxEvents),    
	fMode(qas.fMode), 
	fNumberOfEvents(qas.fNumberOfEvents), 
	fRecoParam(),		
	fRunNumber(qas.fRunNumber), 
	fRawReader(NULL), 
	fRawReaderDelete(kTRUE), 
	fRunLoader(NULL), 
  fTasks(qas.fTasks),
  fEventSpecie(qas.fEventSpecie)
{
	// cpy ctor
	for (UInt_t iDet = 0; iDet < fgkNDetectors; iDet++) {
		fLoader[iDet]         = qas.fLoader[iDet] ;
		fQADataMaker[iDet]    = qas.fQADataMaker[iDet] ;
		fQACycles[iDet]       = qas.fQACycles[iDet] ;	
    fQAWriteExpert[iDet] = qas.fQAWriteExpert[iDet] ;
  }
  // set the default storage
    TString sto(AliQA::GetQARefStorage()) ;
    if (sto.IsNull()) 
      AliFatal("QA reference storage not set do: AliQA::SetQARefStorage") ; 
    SetDefaultStorage(sto.Data());
}

//_____________________________________________________________________________
AliQAManager & AliQAManager::operator = (const AliQAManager & qas) 
{
	// assignment operator
  this->~AliQAManager() ;
  new(this) AliQAManager(qas) ;
  return *this ;
}

//_____________________________________________________________________________
AliQAManager::~AliQAManager() 
{
	// dtor
	for (UInt_t iDet = 0; iDet < fgkNDetectors; iDet++) {
		if (IsSelected(AliQA::GetDetName(iDet))) {
		  fLoader[iDet] = NULL;
		  if (fQADataMaker[iDet]) {
			  (fQADataMaker[iDet])->Finish() ; 
				delete fQADataMaker[iDet] ;
		  }
		}
	}
	if (fRawReaderDelete) { 
		fRunLoader = NULL ;
		delete fRawReader ;
		fRawReader = NULL ;
	}
}

//_____________________________________________________________________________
Bool_t AliQAManager::DoIt(const AliQA::TASKINDEX_t taskIndex)
{
	// Runs all the QA data Maker for every detector
		
	Bool_t rv = kFALSE ;
    // Fill QA data in event loop 
	for (UInt_t iEvent = fFirstEvent ; iEvent < (UInt_t)fMaxEvents ; iEvent++) {
		fCurrentEvent++ ; 
		// Get the event
		if ( iEvent%10 == 0  ) 
			AliInfo(Form("processing event %d", iEvent));
		if ( taskIndex == AliQA::kRAWS ) {
			if ( !fRawReader->NextEvent() )
				break ;
		} else if ( taskIndex == AliQA::kESDS ) {
			if ( fESDTree->GetEntry(iEvent) == 0 )
				break ;
		} else {
			if ( fRunLoader->GetEvent(iEvent) != 0 )
				break ;
		}
		// loop  over active loaders
		for (UInt_t iDet = 0; iDet < fgkNDetectors ; iDet++) {
			if (IsSelected(AliQA::GetDetName(iDet))) {
				AliQADataMaker * qadm = GetQADataMaker(iDet) ;
				if (!qadm) continue; // This detector doesn't have any QA (for example, HLT)
				if ( qadm->IsCycleDone() ) {
          qadm->EndOfCycle(taskIndex) ;
				}
				TTree * data = NULL ; 
				AliLoader* loader = GetLoader(qadm->GetUniqueID());
				switch (taskIndex) {
					case AliQA::kNULLTASKINDEX : 
						break ; 
					case AliQA::kRAWS :
						qadm->Exec(taskIndex, fRawReader) ; 
						break ; 
					case AliQA::kHITS :
            if( loader ) {
              loader->LoadHits() ; 
              data = loader->TreeH() ; 
              if ( ! data ) {
                AliWarning(Form(" Hit Tree not found for  %s", AliQA::GetDetName(iDet))) ; 
                break ; 
              } 
            } 
            qadm->Exec(taskIndex, data) ;
						break ;
						case AliQA::kSDIGITS :
            if( loader ) {      
              loader->LoadSDigits() ; 
              data = loader->TreeS() ; 
              if ( ! data ) {
                AliWarning(Form(" SDigit Tree not found for  %s", AliQA::GetDetName(iDet))) ; 
                break ; 
              } 
            }
            qadm->Exec(taskIndex, data) ; 
						break; 
						case AliQA::kDIGITS :
            if( loader ) {      
              loader->LoadDigits() ; 
              data = loader->TreeD() ; 
              if ( ! data ) {
                AliWarning(Form(" Digit Tree not found for  %s", AliQA::GetDetName(iDet))) ; 
                break ; 
              } 
            }
            qadm->Exec(taskIndex, data) ;
						break; 
						case AliQA::kRECPOINTS :
            if( loader ) {      
              loader->LoadRecPoints() ; 
              data = loader->TreeR() ; 
              if (!data) {
                AliWarning(Form("RecPoints not found for %s", AliQA::GetDetName(iDet))) ; 
                break ; 
              } 
            }
            qadm->Exec(taskIndex, data) ; 
            break; 
						case AliQA::kTRACKSEGMENTS :
						break; 
						case AliQA::kRECPARTICLES :
						break; 
						case AliQA::kESDS :
						qadm->Exec(taskIndex, fESD) ;
						break; 
						case AliQA::kNTASKINDEX :
						break; 
				} //task switch
			}
		} // detector loop
    Increment() ; 
	} // event loop	
	// Save QA data for all detectors
	rv = Finish(taskIndex) ;
	
	if ( taskIndex == AliQA::kRAWS ) 
		fRawReader->RewindEvents() ;

	return rv ; 
}

//_____________________________________________________________________________
Bool_t AliQAManager::Finish(const AliQA::TASKINDEX_t taskIndex) 
{
	// write output to file for all detectors
	for (UInt_t iDet = 0; iDet < fgkNDetectors ; iDet++) {
		if (IsSelected(AliQA::GetDetName(iDet))) {
			AliQADataMaker * qadm = GetQADataMaker(iDet) ;
			if (qadm) 
        qadm->EndOfCycle(taskIndex) ;
		}
	}
	return kTRUE ; 
}

//_____________________________________________________________________________
TObjArray * AliQAManager::GetFromOCDB(AliQA::DETECTORINDEX_t det, AliQA::TASKINDEX_t task, const char * year) const 
{
	// Retrieve the list of QA data for a given detector and a given task 
	TObjArray * rv = NULL ;
	if ( !strlen(AliQA::GetQARefStorage()) ) { 
		AliError("No storage defined, use AliQA::SetQARefStorage") ; 
		return NULL ; 
	}	
	if ( ! IsDefaultStorageSet() ) {
		TString tmp(AliQA::GetQARefDefaultStorage()) ; 
		tmp.Append(year) ; 
		tmp.Append("/") ; 
		Instance()->SetDefaultStorage(tmp.Data()) ; 		
		Instance()->SetSpecificStorage(Form("%s/*", AliQA::GetQAName()), AliQA::GetQARefStorage()) ;
	}
	TString detOCDBDir(Form("%s/%s/%s", AliQA::GetQAName(), AliQA::GetDetName((Int_t)det), AliQA::GetRefOCDBDirName())) ; 
	AliInfo(Form("Retrieving reference data from %s/%s for %s", AliQA::GetQARefStorage(), detOCDBDir.Data(), AliQA::GetTaskName(task).Data())) ; 
	AliCDBEntry* entry = QAManager()->Get(detOCDBDir.Data(), 0) ; //FIXME 0 --> Run Number
	TList * listDetQAD = dynamic_cast<TList *>(entry->GetObject()) ;
	if ( listDetQAD ) 
		rv = dynamic_cast<TObjArray *>(listDetQAD->FindObject(AliQA::GetTaskName(task))) ; 
	return rv ; 
}

//_____________________________________________________________________________
AliLoader * AliQAManager::GetLoader(Int_t iDet)
{
	// get the loader for a detector

	if ( !fRunLoader || iDet == AliQA::kCORR) 
		return NULL ; 
	
	TString detName = AliQA::GetDetName(iDet) ;
    fLoader[iDet] = fRunLoader->GetLoader(detName + "Loader");
	if (fLoader[iDet]) 
		return fLoader[iDet] ;
	
	// load the QA data maker object
	TPluginManager* pluginManager = gROOT->GetPluginManager() ;
	TString loaderName = "Ali" + detName + "Loader" ;

	AliLoader * loader = NULL ;
	// first check if a plugin is defined for the quality assurance data maker
	TPluginHandler* pluginHandler = pluginManager->FindHandler("AliLoader", detName) ;
	// if not, add a plugin for it
	if (!pluginHandler) {
		AliDebug(1, Form("defining plugin for %s", loaderName.Data())) ;
		TString libs = gSystem->GetLibraries() ;
		if (libs.Contains("lib" + detName + "base.so") || (gSystem->Load("lib" + detName + "base.so") >= 0)) {
			pluginManager->AddHandler("AliQADataMaker", detName, loaderName, detName + "loader", loaderName + "()") ;
		} else {
			pluginManager->AddHandler("AliLoader", detName, loaderName, detName, loaderName + "()") ;
		}
		pluginHandler = pluginManager->FindHandler("AliLoader", detName) ;
	}
	if (pluginHandler && (pluginHandler->LoadPlugin() == 0)) {
		loader = (AliLoader *) pluginHandler->ExecPlugin(0) ;
	}
	if (loader) 
		fLoader[iDet] = loader ;
	return loader ;
}

//_____________________________________________________________________________
AliQA * AliQAManager::GetQA(UInt_t run, UInt_t evt) 
{
// retrieves the QA object stored in a file named "Run{run}.Event{evt}_1.ESD.tag.root"  
  char * fileName = Form("Run%d.Event%d_1.ESD.tag.root", run, evt) ; 
  TFile * tagFile = TFile::Open(fileName) ;
  if ( !tagFile ) {
    AliError(Form("File %s not found", fileName)) ;
    return NULL ; 
  }
  TTree * tagTree = dynamic_cast<TTree *>(tagFile->Get("T")) ; 
  if ( !tagTree ) {
    AliError(Form("Tree T not found in %s", fileName)) ; 
    tagFile->Close() ; 
    return NULL ; 
  }
  AliRunTag * tag = new AliRunTag ; 
  tagTree->SetBranchAddress("AliTAG", &tag) ; 
  tagTree->GetEntry(evt) ; 
  AliQA * qa = AliQA::Instance(tag->GetQALength(), tag->GetQA(), tag->GetESLength(), tag->GetEventSpecies()) ; 
  tagFile->Close() ; 
  return qa ; 
}

//_____________________________________________________________________________
AliQADataMaker * AliQAManager::GetQADataMaker(const Int_t iDet) 
{
	// get the quality assurance data maker for a detector
	
	if (fQADataMaker[iDet]) {
    fQADataMaker[iDet]->SetEventSpecie(fEventSpecie) ; 
		return fQADataMaker[iDet] ;
  }
	
	AliQADataMaker * qadm = NULL ;
	
	if (iDet == AliQA::kGLOBAL) { //Global QA
		qadm = new AliGlobalQADataMaker();
		qadm->SetName(AliQA::GetDetName(iDet));
		qadm->SetUniqueID(iDet);
		fQADataMaker[iDet] = qadm;
    qadm->SetEventSpecie(fEventSpecie) ; 
		return qadm;
	}

	if (iDet == AliQA::kCORR) { //the data maker for correlations among detectors
    qadm = new AliCorrQADataMakerRec(fQADataMaker) ; 
		qadm->SetName(AliQA::GetDetName(iDet));
		qadm->SetUniqueID(iDet);
		fQADataMaker[iDet] = qadm;
    qadm->SetEventSpecie(fEventSpecie) ; 
		return qadm;
  }

	// load the QA data maker object
	TPluginManager* pluginManager = gROOT->GetPluginManager() ;
	TString detName = AliQA::GetDetName(iDet) ;
	TString tmp(fMode) ; 
	if (tmp.Contains("sim")) 
		tmp.ReplaceAll("s", "S") ; 
	else if (tmp.Contains("rec")) 
		tmp.ReplaceAll("r", "R") ; 

	TString qadmName = "Ali" + detName + "QADataMaker" + tmp ;

	// first check if a plugin is defined for the quality assurance data maker
	TPluginHandler* pluginHandler = pluginManager->FindHandler("AliQADataMaker", detName) ;
	// if not, add a plugin for it
	if (!pluginHandler) {
		AliDebug(1, Form("defining plugin for %s", qadmName.Data())) ;
		TString libs = gSystem->GetLibraries() ;
		if (libs.Contains("lib" + detName + fMode + ".so") || (gSystem->Load("lib" + detName + fMode + ".so") >= 0)) {
			pluginManager->AddHandler("AliQADataMaker", detName, qadmName, detName + "qadm", qadmName + "()") ;
		} else {
			pluginManager->AddHandler("AliQADataMaker", detName, qadmName, detName, qadmName + "()") ;
		}
		pluginHandler = pluginManager->FindHandler("AliQADataMaker", detName) ;
	}
	if (pluginHandler && (pluginHandler->LoadPlugin() == 0)) {
		qadm = (AliQADataMaker *) pluginHandler->ExecPlugin(0) ;
	}
	if (qadm) {
		qadm->SetName(AliQA::GetDetName(iDet));
		qadm->SetUniqueID(iDet);
		fQADataMaker[iDet] = qadm ;
    qadm->SetEventSpecie(fEventSpecie) ; 
	}

  return qadm ;
}

//_____________________________________________________________________________
void  AliQAManager::EndOfCycle(TObjArray * detArray) 
{
	// End of cycle QADataMakers 
	
	for (UInt_t iDet = 0; iDet < fgkNDetectors ; iDet++) {
		if (IsSelected(AliQA::GetDetName(iDet))) {
			AliQADataMaker * qadm = GetQADataMaker(iDet) ;
			if (!qadm) 
				continue ;	
			// skip non active detectors
			if (detArray) {
				AliModule* det = static_cast<AliModule*>(detArray->FindObject(AliQA::GetDetName(iDet))) ;
				if (!det || !det->IsActive())  
					continue ;
			}
			for (UInt_t taskIndex = 0; taskIndex < AliQA::kNTASKINDEX; taskIndex++) {
				if ( fTasks.Contains(Form("%d", taskIndex)) ) 
					qadm->EndOfCycle(AliQA::GetTaskIndex(AliQA::GetTaskName(taskIndex))) ;
			}
			qadm->Finish();
		}
	}
}

//_____________________________________________________________________________
void  AliQAManager::EndOfCycle(TString detectors) 
{
	// End of cycle QADataMakers 
	
	for (UInt_t iDet = 0; iDet < fgkNDetectors ; iDet++) {
		if (IsSelected(AliQA::GetDetName(iDet))) {
			AliQADataMaker * qadm = GetQADataMaker(iDet) ;
			if (!qadm) 
				continue ;	
			// skip non active detectors
      if (!detectors.Contains(AliQA::GetDetName(iDet))) 
        continue ;
   		for (UInt_t taskIndex = 0; taskIndex < AliQA::kNTASKINDEX; taskIndex++) {
				if ( fTasks.Contains(Form("%d", taskIndex)) ) 
					qadm->EndOfCycle(AliQA::GetTaskIndex(AliQA::GetTaskName(taskIndex))) ;
			}
			qadm->Finish();
		}
	}
}

//_____________________________________________________________________________
void AliQAManager::Increment()
{
  // Increments the cycle counter for all QA Data Makers
 	for (UInt_t iDet = 0; iDet < fgkNDetectors ; iDet++) {
		if (IsSelected(AliQA::GetDetName(iDet))) {
			AliQADataMaker * qadm = GetQADataMaker(iDet) ;
			if (qadm) 
        qadm->Increment() ;
    }
  }
}
  
//_____________________________________________________________________________
Bool_t AliQAManager::InitQA(const AliQA::TASKINDEX_t taskIndex, const  char * input )
{
	// Initialize the event source and QA data makers
	
	fTasks += Form("%d", taskIndex) ; 

	if (taskIndex == AliQA::kRAWS) { 
		if (!fRawReader) {
		        fRawReader = AliRawReader::Create(input);
		}
		if ( ! fRawReader ) 
			return kFALSE ; 
		fRawReaderDelete = kTRUE ; 
		fRawReader->NextEvent() ; 
		fRunNumber = fRawReader->GetRunNumber() ; 
		SetRun(fRunNumber) ; 
		fRawReader->RewindEvents();
		fNumberOfEvents = 999999 ;
		if ( fMaxEvents < 0 ) 
			fMaxEvents = fNumberOfEvents ; 
		} else if (taskIndex == AliQA::kESDS) {
			fTasks = AliQA::GetTaskName(AliQA::kESDS) ; 
      if (!gSystem->AccessPathName("AliESDs.root")) { // AliESDs.root exists
        TFile * esdFile = TFile::Open("AliESDs.root") ;
        fESDTree = dynamic_cast<TTree *> (esdFile->Get("esdTree")) ; 
        if ( !fESDTree ) {
          AliError("esdTree not found") ; 
          return kFALSE ; 
        } else {
          fESD     = new AliESDEvent() ;
          fESD->ReadFromTree(fESDTree) ;
          fESDTree->GetEntry(0) ; 
          fRunNumber = fESD->GetRunNumber() ; 
          fNumberOfEvents = fESDTree->GetEntries() ;
          if ( fMaxEvents < 0 ) 
            fMaxEvents = fNumberOfEvents ; 
        }
      } else {
        AliError("AliESDs.root not found") ; 
        return kFALSE ; 
      }			
    } else {
      if ( !InitRunLoader() ) { 
        AliWarning("No Run Loader not found") ; 
      } else {
        fNumberOfEvents = fRunLoader->GetNumberOfEvents() ;
        if ( fMaxEvents < 0 ) 
          fMaxEvents = fNumberOfEvents ; 
      }
    }

  // Get Detectors 
  TObjArray* detArray = NULL ; 
	if (fRunLoader) // check if RunLoader exists 
		if ( fRunLoader->GetAliRun() ) { // check if AliRun exists in gAlice.root
			detArray = fRunLoader->GetAliRun()->Detectors() ;
			fRunNumber = fRunLoader->GetHeader()->GetRun() ; 
		}

	// Initialize all QA data makers for all detectors
	fRunNumber = AliCDBManager::Instance()->GetRun() ; 
	if ( !  AliGeomManager::GetGeometry() ) 
		AliGeomManager::LoadGeometry() ; 
	
	InitQADataMaker(fRunNumber, detArray) ; //, fCycleSame, kTRUE, detArray) ; 
	return kTRUE ; 
}

//_____________________________________________________________________________
void  AliQAManager::InitQADataMaker(UInt_t run, TObjArray * detArray) 
{
	// Initializes The QADataMaker for all active detectors and for all active tasks 
	for (UInt_t iDet = 0; iDet < fgkNDetectors ; iDet++) {
		if (IsSelected(AliQA::GetDetName(iDet))) {
			AliQADataMaker * qadm = GetQADataMaker(iDet) ;
			if (!qadm) {
				AliError(Form("AliQADataMaker not found for %s", AliQA::GetDetName(iDet))) ; 
				fDetectorsW.ReplaceAll(AliQA::GetDetName(iDet), "") ; 
			} else {
        if (fQAWriteExpert[iDet])
          qadm->SetWriteExpert() ; 
				AliDebug(1, Form("Data Maker found for %s", qadm->GetName())) ; 
				// skip non active detectors
				if (detArray) {
					AliModule* det = static_cast<AliModule*>(detArray->FindObject(AliQA::GetDetName(iDet))) ;
					if (!det || !det->IsActive())  
						continue ;
				}
				if (fQAWriteExpert[iDet]) qadm->SetWriteExpert() ; 
	      // Set default reco params
        Bool_t sameCycle = kFALSE ; 
				for (UInt_t taskIndex = 0; taskIndex < AliQA::kNTASKINDEX; taskIndex++) {
					if ( fTasks.Contains(Form("%d", taskIndex)) ) {
						qadm->Init(AliQA::GetTaskIndex(AliQA::GetTaskName(taskIndex)), GetQACycles(qadm->GetUniqueID())) ;
            qadm->StartOfCycle(AliQA::GetTaskIndex(AliQA::GetTaskName(taskIndex)), run,  sameCycle) ;
            sameCycle = kTRUE ;
					}
				}
			}
		}
	}
}


//_____________________________________________________________________________
Bool_t AliQAManager::InitRunLoader()
{
	// get or create the run loader
	if (fRunLoader) {
		fCycleSame = kTRUE ; 
	} else {
		if (!gSystem->AccessPathName(fGAliceFileName.Data())) { // galice.root exists
			// load all base libraries to get the loader classes
			TString libs = gSystem->GetLibraries() ;
			for (UInt_t iDet = 0; iDet < fgkNDetectors; iDet++) {
				if (!IsSelected(AliQA::GetDetName(iDet))) 
					continue ; 
				TString detName = AliQA::GetDetName(iDet) ;
				if (detName == "HLT") 
					continue;
				if (libs.Contains("lib" + detName + "base.so")) 
					continue;
				gSystem->Load("lib" + detName + "base.so");
			}
			fRunLoader = AliRunLoader::Open(fGAliceFileName.Data());
			if (!fRunLoader) {
				AliError(Form("no run loader found in file %s", fGAliceFileName.Data()));
				return kFALSE;
			}
			fRunLoader->CdGAFile();
			if (fRunLoader->LoadgAlice() == 0) {
				gAlice = fRunLoader->GetAliRun();
			}

			if (!gAlice) {
				AliError(Form("no gAlice object found in file %s", fGAliceFileName.Data()));
				return kFALSE;
			}

		} else {               // galice.root does not exist
			AliError(Form("the file %s does not exist", fGAliceFileName.Data()));
			return kFALSE;
		}
	}

	if (!fRunNumber) { 
		fRunLoader->LoadHeader();
		fRunNumber = fRunLoader->GetHeader()->GetRun() ; 
	}
	return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliQAManager::IsSelected(const char * det) 
{
  // check whether detName is contained in detectors
	// if yes, it is removed from detectors
	
	Bool_t rv = kFALSE;
	const TString detName(det) ;
  // always activates Correlation
  if ( detName.Contains(AliQA::GetDetName(AliQA::kCORR))) {
    rv = kTRUE ; 
  } else {
    // check if all detectors are selected
    if (fDetectors.Contains("ALL")) {
      fDetectors = "ALL";
      rv = kTRUE;
    } else if ((fDetectors.CompareTo(detName) == 0) ||
               fDetectors.BeginsWith(detName+" ") ||
               fDetectors.EndsWith(" "+detName) ||
               fDetectors.Contains(" "+detName+" ")) {
      rv = kTRUE;
    }
  }
	return rv ;
}

//_____________________________________________________________________________
Bool_t AliQAManager::Merge(const Int_t runNumber) const
{
	// Merge data from all the cycles from all detectors in one single file per run
	// Merge the QA results from all the data chunks in one run 
 Bool_t rv = MergeData(runNumber) ; 
 rv *= MergeResults(runNumber) ;
 return rv ; 
}
	
	
//_____________________________________________________________________________
Bool_t AliQAManager::MergeData(const Int_t runNumber) const
{
	// Merge all the cycles from all detectors in one single file per run
	TString cmd ;
	if (runNumber == -1) 
		cmd = Form(".! ls *%s*.%d.root > tempo.txt", AliQA::GetQADataFileName(), runNumber) ; 
	else 
		cmd = Form(".! ls *%s*.*.root > tempo.txt", AliQA::GetQADataFileName()) ; 
	gROOT->ProcessLine(cmd.Data()) ;
	ifstream in("tempo.txt") ; 
	const Int_t runMax = 10 ;  
	TString file[AliQA::kNDET*runMax] ;
	
	Int_t index = 0 ; 
	while ( 1 ) {
		in >> file[index] ; 
		if ( !in.good() ) 
			break ; 
		AliInfo(Form("index = %d file = %s", index, (file[index]).Data())) ; 
		index++ ;
	}
	
	if ( index == 0 ) { 
		AliError(Form("run number %d not found", runNumber)) ; 
		return kFALSE ; 
	}

  TFileMerger merger ; 
  TString outFileName(Form("Merged.%s.Data.%d.root",AliQA::GetQADataFileName(),runNumber)); 
  merger.OutputFile(outFileName.Data()) ; 
  for (Int_t ifile = 0 ; ifile < index-1 ; ifile++) {
    TString pattern(Form("%s.%d.", AliQA::GetQADataFileName(), runNumber)); 
    TString tmp(file[ifile]) ; 
    if (tmp.Contains(pattern)) {
      merger.AddFile(tmp) ; 
    }
	}
  merger.Merge() ; 
	return kTRUE ; 
}

//_____________________________________________________________________________
Bool_t AliQAManager::MergeResults(const Int_t runNumber) const
{
	// Merge the QA result from all the data chunks in a run 
	TString cmd ;
	cmd = Form(".! ls %s*.root > tempo.txt", AliQA::GetQADataFileName()) ; 
	gROOT->ProcessLine(cmd.Data()) ;
	ifstream in("tempo.txt") ; 
	const Int_t chunkMax = 100 ;  
	TString fileList[chunkMax] ;
	
	Int_t index = 0 ; 
	while ( 1 ) {
		TString file ; 
		in >> fileList[index] ; 
		if ( !in.good() ) 
			break ; 
		AliInfo(Form("index = %d file = %s", index, (fileList[index].Data()))) ; 
		index++ ;
	}
	
	if ( index == 0 ) { 
		AliError("No QA Result File found") ; 
		return kFALSE ; 
	}
	
	TFileMerger merger ; 
	TString outFileName(Form("Merged.%s.Result.%d.root", AliQA::GetQADataFileName(), runNumber));		
	merger.OutputFile(outFileName.Data()) ; 
	for (Int_t ifile = 0 ; ifile < index ; ifile++) {
		TString file = fileList[ifile] ; 
		merger.AddFile(file) ; 
	}
	merger.Merge() ; 
	
	return kTRUE ; 
}

//_____________________________________________________________________________
void AliQAManager::Reset(const Bool_t sameCycle)
{
	// Reset the default data members

	for (UInt_t iDet = 0; iDet < fgkNDetectors ; iDet++) {
		if (IsSelected(AliQA::GetDetName(iDet))) {
			AliQADataMaker * qadm = GetQADataMaker(iDet);
			qadm->Reset();
		}
	} 
	if (fRawReaderDelete) { 
		delete fRawReader ;
		fRawReader      = NULL ;
	}

	fCycleSame      = sameCycle ; 
	fESD            = NULL ; 
	fESDTree        = NULL ; 
	//fFirst          = kTRUE ;   
	fNumberOfEvents = 999999 ;  
}

//_____________________________________________________________________________
AliQAManager * AliQAManager::QAManager(const Char_t * mode, TMap *entryCache, Int_t run) 
{
  // returns AliQAManager instance (singleton)
  
	if (!fgQAInstance) {
    fgQAInstance = new AliQAManager(mode) ;  
    if (!entryCache)
		  fgQAInstance->Init();
		else
		  fgQAInstance->InitFromCache(entryCache,run);
  }
	return fgQAInstance;
}

//_____________________________________________________________________________
TString AliQAManager::Run(const char * detectors, AliRawReader * rawReader, const Bool_t sameCycle) 
{
	//Runs all the QA data Maker for Raws only
	
	fCycleSame       = sameCycle ;
	fRawReader       = rawReader ;
	fDetectors       = detectors ; 
	fDetectorsW      = detectors ; 	
	
	AliCDBManager* man = AliCDBManager::Instance() ; 

	if ( man->GetRun() == -1 ) {// check if run number not set previously and set it from raw data
		rawReader->NextEvent() ; 
		man->SetRun(fRawReader->GetRunNumber()) ;
		rawReader->RewindEvents() ;
	}	
	
	if (!fCycleSame) 
    if ( !InitQA(AliQA::kRAWS) ) 
      return kFALSE ; 
  fRawReaderDelete = kFALSE ; 

	DoIt(AliQA::kRAWS) ; 
	return 	fDetectorsW ;
}

//_____________________________________________________________________________
TString AliQAManager::Run(const char * detectors, const char * fileName, const Bool_t sameCycle) 
{
	//Runs all the QA data Maker for Raws only

	fCycleSame       = sameCycle ;
	fDetectors       = detectors ; 
	fDetectorsW      = detectors ; 	
	
	AliCDBManager* man = AliCDBManager::Instance() ; 
	if ( man->GetRun() == -1 ) { // check if run number not set previously and set it from AliRun
		AliRunLoader * rl = AliRunLoader::Open("galice.root") ;
		if ( ! rl ) {
			AliFatal("galice.root file not found in current directory") ; 
		} else {
			rl->CdGAFile() ; 
			rl->LoadgAlice() ;
			if ( ! rl->GetAliRun() ) {
				AliFatal("AliRun not found in galice.root") ;
			} else {
				rl->LoadHeader() ;
				man->SetRun(rl->GetHeader()->GetRun());
			}
		}
	}
	
	if (!fCycleSame) 
    if ( !InitQA(AliQA::kRAWS, fileName) ) 
      return kFALSE ; 
	
	DoIt(AliQA::kRAWS) ; 
	return 	fDetectorsW ;
}

//_____________________________________________________________________________
TString AliQAManager::Run(const char * detectors, const AliQA::TASKINDEX_t taskIndex, Bool_t const sameCycle, const  char * fileName ) 
{
	// Runs all the QA data Maker for every detector
	
	fCycleSame       = sameCycle ;
	fDetectors       = detectors ; 
	fDetectorsW      = detectors ; 		
	
	AliCDBManager* man = AliCDBManager::Instance() ; 	
	if ( man->GetRun() == -1 ) { // check if run number not set previously and set it from AliRun
		AliRunLoader * rl = AliRunLoader::Open("galice.root") ;
		if ( ! rl ) {
			AliFatal("galice.root file not found in current directory") ; 
		} else {
			rl->CdGAFile() ; 
			rl->LoadgAlice() ;
			if ( ! rl->GetAliRun() ) {
				AliInfo("AliRun not found in galice.root") ;
			} else {
				rl->LoadHeader() ;
				man->SetRun(rl->GetHeader()->GetRun()) ;
			}
		}
	}
	

	if ( taskIndex == AliQA::kNULLTASKINDEX) { 
		for (UInt_t task = 0; task < AliQA::kNTASKINDEX; task++) {
			if ( fTasks.Contains(Form("%d", task)) ) {
        if (!fCycleSame)
          if ( !InitQA(AliQA::GetTaskIndex(AliQA::GetTaskName(task)), fileName) ) 
            return kFALSE ;
        DoIt(AliQA::GetTaskIndex(AliQA::GetTaskName(task))) ;
			}
		}
	} else {
    if (! fCycleSame )
      if ( !InitQA(taskIndex, fileName) ) 
        return kFALSE ; 
      DoIt(taskIndex) ; 
  } 		
	
	return fDetectorsW ;

}

//_____________________________________________________________________________
void AliQAManager::RunOneEvent(AliRawReader * rawReader) 
{
	//Runs all the QA data Maker for Raws only and on one event only (event loop done by calling method)
  if ( ! rawReader ) 
    return ; 
	AliCodeTimerAuto("") ;
  if (fTasks.Contains(Form("%d", AliQA::kRAWS))){
    for (UInt_t iDet = 0; iDet < fgkNDetectors; iDet++) {
      if (!IsSelected(AliQA::GetDetName(iDet))) 
        continue;
      AliQADataMaker *qadm = GetQADataMaker(iDet);  
      if (!qadm) 
        continue;
      if ( qadm->IsCycleDone() ) {
        qadm->EndOfCycle() ;
      }
      AliCodeTimerStart(Form("running RAW quality assurance data maker for %s", AliQA::GetDetName(iDet))); 
      qadm->SetEventSpecie(fEventSpecie) ; 
			qadm->Exec(AliQA::kRAWS, rawReader) ;
      AliCodeTimerStop(Form("running RAW quality assurance data maker for %s", AliQA::GetDetName(iDet)));
		}
  }
}

//_____________________________________________________________________________
void AliQAManager::RunOneEvent(AliESDEvent *& esd) 
{
	//Runs all the QA data Maker for ESDs only and on one event only (event loop done by calling method)
	
  AliCodeTimerAuto("") ;
  if (fTasks.Contains(Form("%d", AliQA::kESDS))) {
    for (UInt_t iDet = 0; iDet < fgkNDetectors; iDet++) {
      if (!IsSelected(AliQA::GetDetName(iDet))) 
        continue;
      AliQADataMaker *qadm = GetQADataMaker(iDet);  
      if (!qadm) 
        continue;
      if ( qadm->IsCycleDone() ) {
        qadm->EndOfCycle() ;
      }
      AliCodeTimerStart(Form("running ESD quality assurance data maker for %s", AliQA::GetDetName(iDet)));
			qadm->Exec(AliQA::kESDS, esd) ;
      AliCodeTimerStop(Form("running ESD quality assurance data maker for %s", AliQA::GetDetName(iDet)));
		}
	}
}

//_____________________________________________________________________________
void AliQAManager::RunOneEventInOneDetector(Int_t det, TTree * tree) 
{
	// Runs all the QA data Maker for ESDs only and on one event only (event loop done by calling method)
	AliCodeTimerAuto("") ;
  if (fTasks.Contains(Form("%d", AliQA::kRECPOINTS))) {
    if (IsSelected(AliQA::GetDetName(det))) {
      AliQADataMaker *qadm = GetQADataMaker(det);  
      if (qadm) { 
        if ( qadm->IsCycleDone() ) {
          qadm->EndOfCycle() ;
        }
        AliCodeTimerStart(Form("running RecPoints quality assurance data maker for %s", AliQA::GetDetName(det)));
        qadm->Exec(AliQA::kRECPOINTS, tree) ;
        AliCodeTimerStop(Form("running RecPoints quality assurance data maker for %s", AliQA::GetDetName(det)));
      }
    }
  }
}

//_____________________________________________________________________________
Bool_t AliQAManager::Save2OCDB(const Int_t runNumber, AliRecoParam::EventSpecie_t es, const char * year, const char * detectors) const
{
	// take the locasl QA data merge into a single file and save in OCDB 
	Bool_t rv = kTRUE ; 
	TString tmp(AliQA::GetQARefStorage()) ; 
	if ( tmp.IsNull() ) { 
		AliError("No storage defined, use AliQA::SetQARefStorage") ; 
		return kFALSE ; 
	}
	if ( !(tmp.Contains(AliQA::GetLabLocalOCDB()) || tmp.Contains(AliQA::GetLabAliEnOCDB())) ) {
		AliError(Form("%s is a wrong storage, use %s or %s", AliQA::GetQARefStorage(), AliQA::GetLabLocalOCDB().Data(), AliQA::GetLabAliEnOCDB().Data())) ; 
		return kFALSE ; 
	}
	TString sdet(detectors) ; 
	sdet.ToUpper() ;
	TFile * inputFile ; 
	if ( sdet.Contains("ALL") ) {
		rv = Merge(runNumber) ; 
		if ( ! rv )
			return kFALSE ; 
		TString inputFileName(Form("Merged.%s.Data.%d.root", AliQA::GetQADataFileName(), runNumber)) ; 
		inputFile = TFile::Open(inputFileName.Data()) ; 
		rv = SaveIt2OCDB(runNumber, inputFile, year, es) ; 
	} else {
		for (Int_t index = 0; index < AliQA::kNDET; index++) {
			if (sdet.Contains(AliQA::GetDetName(index))) {
				TString inputFileName(Form("%s.%s.%d.root", AliQA::GetDetName(index), AliQA::GetQADataFileName(), runNumber)) ; 
				inputFile = TFile::Open(inputFileName.Data()) ; 			
				rv *= SaveIt2OCDB(runNumber, inputFile, year, es) ; 
			}
		}
	}
	return rv ; 
}

//_____________________________________________________________________________
Bool_t AliQAManager::SaveIt2OCDB(const Int_t runNumber, TFile * inputFile, const char * year, AliRecoParam::EventSpecie_t es) const
{
	// reads the TH1 from file and adds it to appropriate list before saving to OCDB
	Bool_t rv = kTRUE ;
	AliInfo(Form("Saving TH1s in %s to %s", inputFile->GetName(), AliQA::GetQARefStorage())) ; 
	if ( ! IsDefaultStorageSet() ) {
		TString tmp( AliQA::GetQARefStorage() ) ; 
		if ( tmp.Contains(AliQA::GetLabLocalOCDB()) ) 
			Instance()->SetDefaultStorage(AliQA::GetQARefStorage()) ;
		else {
			TString tmp1(AliQA::GetQARefDefaultStorage()) ; 
			tmp1.Append(year) ; 
			tmp1.Append("?user=alidaq") ; 
			Instance()->SetDefaultStorage(tmp1.Data()) ; 
		}
	}
	Instance()->SetSpecificStorage("*", AliQA::GetQARefStorage()) ; 
	if(GetRun() < 0) 
		Instance()->SetRun(runNumber);

	AliCDBMetaData mdr ;
	mdr.SetResponsible("yves schutz");

	for ( Int_t detIndex = 0 ; detIndex < AliQA::kNDET ; detIndex++) {
		TDirectory * detDir = inputFile->GetDirectory(AliQA::GetDetName(detIndex)) ; 
		if ( detDir ) {
			AliInfo(Form("Entering %s", detDir->GetName())) ;
      AliQA::SetQARefDataDirName(es) ;
			TString detOCDBDir(Form("%s/%s/%s", AliQA::GetDetName(detIndex), AliQA::GetRefOCDBDirName(), AliQA::GetRefDataDirName())) ; 
			AliCDBId idr(detOCDBDir.Data(), runNumber, AliCDBRunRange::Infinity())  ;
			TList * listDetQAD = new TList() ;
			TString listName(Form("%s QA data Reference", AliQA::GetDetName(detIndex))) ; 
			mdr.SetComment(Form("%s QA stuff", AliQA::GetDetName(detIndex)));
			listDetQAD->SetName(listName) ; 
			TList * taskList = detDir->GetListOfKeys() ; 
			TIter nextTask(taskList) ; 
			TKey * taskKey ; 
			while ( (taskKey = dynamic_cast<TKey*>(nextTask())) ) {
				TDirectory * taskDir = detDir->GetDirectory(taskKey->GetName()) ; 
        TDirectory * esDir   = taskDir->GetDirectory(AliRecoParam::GetEventSpecieName(es)) ; 
				AliInfo(Form("Saving %s", esDir->GetName())) ; 
				TObjArray * listTaskQAD = new TObjArray(100) ; 
				listTaskQAD->SetName(Form("%s/%s", taskKey->GetName(), AliRecoParam::GetEventSpecieName(es))) ;
				listDetQAD->Add(listTaskQAD) ; 
				TList * histList = esDir->GetListOfKeys() ; 
				TIter nextHist(histList) ; 
				TKey * histKey ; 
				while ( (histKey = dynamic_cast<TKey*>(nextHist())) ) {
					TObject * odata = esDir->Get(histKey->GetName()) ; 
					if ( !odata ) {
						AliError(Form("%s in %s/%s returns a NULL pointer !!", histKey->GetName(), detDir->GetName(), taskDir->GetName())) ;
					} else {
            if ( AliQA::GetExpert() == histKey->GetName() ) {
              TDirectory * expertDir   = esDir->GetDirectory(histKey->GetName()) ; 
              TList * expertHistList = expertDir->GetListOfKeys() ; 
              TIter nextExpertHist(expertHistList) ; 
              TKey * expertHistKey ; 
              while ( (expertHistKey = dynamic_cast<TKey*>(nextExpertHist())) ) {
                TObject * expertOdata = expertDir->Get(expertHistKey->GetName()) ; 
                if ( !expertOdata ) {
                  AliError(Form("%s in %s/%s/Expert returns a NULL pointer !!", expertHistKey->GetName(), detDir->GetName(), taskDir->GetName())) ;
                } else {
                  AliInfo(Form("Adding %s", expertHistKey->GetName())) ;
                  if ( expertOdata->IsA()->InheritsFrom("TH1") ) {
                    AliInfo(Form("Adding %s", expertHistKey->GetName())) ;
                    TH1 * hExpertdata = static_cast<TH1*>(expertOdata) ; 
                    listTaskQAD->Add(hExpertdata) ; 
                  }                  
                }                
              }
            }
						AliInfo(Form("Adding %s", histKey->GetName())) ;
						if ( odata->IsA()->InheritsFrom("TH1") ) {
							AliInfo(Form("Adding %s", histKey->GetName())) ;
							TH1 * hdata = static_cast<TH1*>(odata) ; 
							listTaskQAD->Add(hdata) ; 
						}
					}
				}
			}
			Instance()->Put(listDetQAD, idr, &mdr) ;
		}
	}
	return rv ; 
}	

//_____________________________________________________________________________
void AliQAManager::SetEventSpecie(AliRecoParam::EventSpecie_t es) 
{
  // set the current event specie and inform AliQA that this event specie has been encountered
  fEventSpecie = es ;
  AliQA::Instance()->SetEventSpecie(es) ; 
}

//_____________________________________________________________________________
void AliQAManager::SetRecoParam(const Int_t det, const AliDetectorRecoParam *par) 
{
  // Set custom reconstruction parameters for a given detector
  // Single set of parameters for all the events
  GetQADataMaker(det)->SetRecoParam(par) ; 
}


