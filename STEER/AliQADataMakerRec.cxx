/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notifce   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


/* $Id$ */

//
//  Base Class
//  Produces the data needed to calculate the quality assurance for Reconstruction
//  All data must be mergeable objects.
//  Y. Schutz CERN July 2007
//

// --- ROOT system ---
#include <TFile.h>
#include <TTree.h>
#include <TNtupleD.h>
#include <TObjArray.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#include "AliDetectorRecoParam.h"
#include "AliCDBManager.h"

#include "AliLog.h"
#include "AliQADataMakerRec.h"
#include "AliESDEvent.h"
#include "AliRawReader.h"

ClassImp(AliQADataMakerRec)
             
//____________________________________________________________________________ 
AliQADataMakerRec::AliQADataMakerRec(const char * name, const char * title) : 
  AliQADataMaker(name, title), 
  fESDsQAList(NULL), 
  fRawsQAList(NULL), 
  fRecPointsQAList(NULL),
  fCorrNt(NULL), 
  fRecoParam(NULL) 
{
  // ctor
	fDetectorDirName = GetName() ; 
}

//____________________________________________________________________________ 
AliQADataMakerRec::AliQADataMakerRec(const AliQADataMakerRec& qadm) :
  AliQADataMaker(qadm.GetName(), qadm.GetTitle()), 
  fESDsQAList(qadm.fESDsQAList),
  fRawsQAList(qadm.fRawsQAList),
  fRecPointsQAList(qadm.fRecPointsQAList),
  fCorrNt(qadm.fCorrNt),  
  fRecoParam(qadm.fRecoParam) 
{
  //copy ctor
 	SetName(qadm.GetName()) ; 
	SetTitle(qadm.GetTitle()) ; 
	fDetectorDirName = GetName() ; 
}

//____________________________________________________________________________ 
AliQADataMakerRec::~AliQADataMakerRec()
{
	//dtor: delete the TObjArray and thei content
	if ( fESDsQAList ) {
    for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
      if ( fESDsQAList[specie] ) {
        if ( fESDsQAList[specie]->IsOwner() ) 
          fESDsQAList[specie]->Delete() ;     
      }
    }
    delete[] fESDsQAList ;
	}
	if ( fRawsQAList ) {
    for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
      if ( fRawsQAList[specie] ) {
        if ( fRawsQAList[specie]->IsOwner() ) 
          fRawsQAList[specie]->Delete() ;
      }
    }
    delete[] fRawsQAList ;
  }
	if ( fRecPointsQAList ) {
    for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
      if ( fRecPointsQAList[specie] ) {
        if ( fRecPointsQAList[specie]->IsOwner() ) 
          fRecPointsQAList[specie]->Delete() ;
      }
    }
		delete[] fRecPointsQAList ; 
	}
}

//__________________________________________________________________
AliQADataMakerRec& AliQADataMakerRec::operator = (const AliQADataMakerRec& qadm )
{
  // Assignment operator.
  this->~AliQADataMakerRec();
  new(this) AliQADataMakerRec(qadm);
  return *this;
}

//____________________________________________________________________________
void AliQADataMakerRec::EndOfCycle() 
{
  // Finishes a cycle of QA for all the tasks
  EndOfCycle(AliQA::kRAWS) ; 
  EndOfCycle(AliQA::kRECPOINTS) ; 
  EndOfCycle(AliQA::kESDS) ; 
  ResetCycle() ; 
}

//____________________________________________________________________________
void AliQADataMakerRec::EndOfCycle(AliQA::TASKINDEX_t task) 
{
	// Finishes a cycle of QA 
	
	TObjArray ** list = NULL ; 
	
	if ( task == AliQA::kRAWS )     
		list = fRawsQAList ; 
	else if ( task == AliQA::kRECPOINTS ) 
		list = fRecPointsQAList ; 
	else if ( task == AliQA::kESDS )
		list = fESDsQAList ; 

 
	if ( ! list && ! fCorrNt ) 
    return ; 
  //DefaultEndOfDetectorCycle(task) ;
	EndOfDetectorCycle(task, list) ;
	TDirectory * subDir = NULL ;
	if (fDetectorDir) 
		subDir = fDetectorDir->GetDirectory(AliQA::GetTaskName(task)) ; 
	if ( subDir ) {
		subDir->cd() ; 
    for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
      TDirectory * eventSpecieDir = subDir->GetDirectory(AliRecoParam::GetEventSpecieName(specie)) ;
      if (eventSpecieDir) {
        eventSpecieDir->cd() ;    
        if (list[specie]) {
          TIter next(list[specie]) ; 
          TObject * obj ; 
          while( (obj = next()) ) {
            if (!obj->TestBit(AliQA::GetExpertBit()))
              obj->Write() ;
          }
          if (WriteExpert()) {
            TDirectory * expertDir = eventSpecieDir->GetDirectory(AliQA::GetExpert()) ; 
            if ( expertDir ) { // Write only if requested
              expertDir->cd() ;
              next.Reset() ; 
              while( (obj = next()) ) {
                if (!obj->TestBit(AliQA::GetExpertBit()))
                  continue ; 
                obj->Write() ;
              }      
            }
          }
        }
        if ( !fCorrNt )
          continue ; 
        if (fCorrNt[specie] && AliQA::GetDetIndex(GetName()) == AliQA::kCORR) {
          eventSpecieDir->cd() ; 
          fCorrNt[specie]->Write() ; 
        }
      }
    }
    fOutput->Save() ; 
	}
}
 
//____________________________________________________________________________
void AliQADataMakerRec::Exec(AliQA::TASKINDEX_t task, TObject * data) 
{ 
  // creates the quality assurance data for the various tasks (Hits, SDigits, Digits, ESDs)
	
	if ( task == AliQA::kRAWS ) {
		AliDebug(1, "Processing Raws QA") ; 
		AliRawReader * rawReader = dynamic_cast<AliRawReader *>(data) ; 
		if (rawReader) 
			MakeRaws(rawReader) ;
		else
		AliInfo("Raw data are not processed") ;     
	} else if ( task == AliQA::kRECPOINTS ) {
		AliDebug(1, "Processing RecPoints QA") ; 
		TTree * tree = dynamic_cast<TTree *>(data) ; 
		if (tree) {
			MakeRecPoints(tree) ; 
		} else {
			AliWarning("data are not a TTree") ; 
		}
	} else if ( task == AliQA::kESDS ) {
		AliDebug(1, "Processing ESDs QA") ; 
		AliESDEvent * esd = dynamic_cast<AliESDEvent *>(data) ; 
		if (esd) 
			MakeESDs(esd) ;
		else 
			AliError("Wrong type of esd container") ; 
	}  
}

//____________________________________________________________________________ 
TObjArray **  AliQADataMakerRec::Init(AliQA::TASKINDEX_t task, Int_t cycles)
{
  // general intialisation
  InitRecoParams() ;
	TObjArray ** rv = NULL ; 
  
	if (cycles > 0)
		SetCycle(cycles) ;  
	
	if ( task == AliQA::kRAWS ) {
		if (! fRawsQAList ) { 
      fRawsQAList = new TObjArray *[AliRecoParam::kNSpecies] ; 
      for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
        fRawsQAList[specie] = new TObjArray(100) ;	 
        fRawsQAList[specie]->SetName(Form("%s_%s_%s", GetName(), AliQA::GetTaskName(task).Data(), AliRecoParam::GetEventSpecieName(specie))) ;
      }
			InitRaws() ;
		}
		rv = fRawsQAList ;
	} else if ( task == AliQA::kRECPOINTS ) {
		if ( ! fRecPointsQAList ) {
      fRecPointsQAList = new TObjArray *[AliRecoParam::kNSpecies] ; 
      for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
        fRecPointsQAList[specie] = new TObjArray(100) ; 
        fRecPointsQAList[specie]->SetName(Form("%s_%s_%s", GetName(), AliQA::GetTaskName(task).Data(), AliRecoParam::GetEventSpecieName(specie))) ; 
      }
      InitRecPoints() ;
		}
		rv = fRecPointsQAList ;
	} else if ( task == AliQA::kESDS ) {
		if ( ! fESDsQAList ) {
      fESDsQAList = new TObjArray *[AliRecoParam::kNSpecies] ; 
      for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
        fESDsQAList[specie] = new TObjArray(100) ;
        fESDsQAList[specie]->SetName(Form("%s_%s", GetName(), AliQA::GetTaskName(task).Data(), AliRecoParam::GetEventSpecieName(specie))) ; 
      }
			InitESDs() ;
		}
		rv = fESDsQAList ;
	}
	return rv ; 
}

//____________________________________________________________________________ 
void AliQADataMakerRec::Init(AliQA::TASKINDEX_t task, TObjArray ** list, Int_t run, Int_t cycles)
{
  // Intialisation by passing the list of QA data booked elsewhere
  
  InitRecoParams() ;
  fRun = run ;
	if (cycles > 0)
		SetCycle(cycles) ;  
	
	if ( task == AliQA::kRAWS ) {
		fRawsQAList = list ;	 
	} else if ( task == AliQA::kRECPOINTS ) {
		fRecPointsQAList = list ; 
	} else if ( task == AliQA::kESDS ) {
		fESDsQAList = list ; 
	}
}

//____________________________________________________________________________
void AliQADataMakerRec::InitRecoParams() 
{
  if (!fRecoParam) {
    AliInfo(Form("Loading reconstruction parameter objects for detector %s", GetName()));
    AliCDBPath path(GetName(),"Calib","RecoParam");
    AliCDBEntry *entry=AliCDBManager::Instance()->Get(path.GetPath());
    if(!entry) {
      fRecoParam = NULL ; 
      AliWarning(Form("Couldn't find RecoParam entry in OCDB for detector %s",GetName()));
    }
    else {
      TObject * recoParamObj = entry->GetObject() ; 
      if (dynamic_cast<AliDetectorRecoParam*>(recoParamObj)) {
        // The detector has only onse set of reco parameters
        // Registering it in AliRecoParam
        AliInfo(Form("Single set of reconstruction parameters found for detector %s",GetName()));
        dynamic_cast<AliDetectorRecoParam*>(recoParamObj)->SetAsDefault();
        fRecoParam = dynamic_cast<AliDetectorRecoParam*>(recoParamObj) ;
      } else { 
        AliError(Form("No valid RecoParam object found in the OCDB for detector %s",GetName()));
      }
    }
  }
}

//____________________________________________________________________________
void AliQADataMakerRec::StartOfCycle(Int_t run) 
{
  // Finishes a cycle of QA for all the tasks
  Bool_t samecycle = kFALSE ; 
  StartOfCycle(AliQA::kRAWS,      run, samecycle) ;
  samecycle = kTRUE ; 
  StartOfCycle(AliQA::kRECPOINTS, run, samecycle) ; 
  StartOfCycle(AliQA::kESDS,      run, samecycle) ; 
}

//____________________________________________________________________________
void AliQADataMakerRec::StartOfCycle(AliQA::TASKINDEX_t task, Int_t run, const Bool_t sameCycle) 
{ 
  // Finishes a cycle of QA data acquistion
  if ( run > 0 ) 
    fRun = run ; 
	if ( !sameCycle || fCurrentCycle == -1) {
		ResetCycle() ;
		if (fOutput) 
			fOutput->Close() ; 
		fOutput = AliQA::GetQADataFile(GetName(), fRun) ; 	
	}	
	AliInfo(Form(" Run %d Cycle %d task %s file %s", 
				 fRun, fCurrentCycle, AliQA::GetTaskName(task).Data(), fOutput->GetName() )) ;

	fDetectorDir = fOutput->GetDirectory(GetDetectorDirName()) ; 
	if (!fDetectorDir)
		fDetectorDir = fOutput->mkdir(GetDetectorDirName()) ; 

	TDirectory * subDir = fDetectorDir->GetDirectory(AliQA::GetTaskName(task)) ; 
	if (!subDir)
		subDir = fDetectorDir->mkdir(AliQA::GetTaskName(task)) ;  
  
  for ( Int_t specie = AliRecoParam::kDefault ; specie < AliRecoParam::kNSpecies ; specie++ ) {
    TDirectory * eventSpecieDir = subDir->GetDirectory(AliRecoParam::GetEventSpecieName(specie)) ; 
    if (!eventSpecieDir) 
      eventSpecieDir = subDir->mkdir(AliRecoParam::GetEventSpecieName(specie)) ; 
    TDirectory * expertDir = eventSpecieDir->GetDirectory(AliQA::GetExpert()) ; 
    if (!expertDir)
      expertDir = eventSpecieDir->mkdir(AliQA::GetExpert()) ; 
  } 
	StartOfDetectorCycle() ; 
}
