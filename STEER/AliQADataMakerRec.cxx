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

// --- Standard library ---

// --- AliRoot header files ---
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
  fObject(NULL), 
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
  fObject(qadm.fObject),  
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
		if ( fESDsQAList->IsOwner() ) 
			fESDsQAList->Delete() ;     
		delete fESDsQAList ;     
	}
	if ( fRawsQAList ) {
		if ( fRawsQAList->IsOwner() ) 
			fRawsQAList->Delete() ;
		delete fRawsQAList ;
	}
	if ( fRecPointsQAList ) {
		if ( fRecPointsQAList->IsOwner() ) 
			fRecPointsQAList->Delete() ; 
		delete fRecPointsQAList ; 
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
void AliQADataMakerRec::EndOfCycle(AliQA::TASKINDEX_t task) 
{
	// Finishes a cycle of QA data acquistion
	
	TObjArray * list = NULL ; 
	
	if ( task == AliQA::kRAWS )     
		list = fRawsQAList ; 
	else if ( task == AliQA::kRECPOINTS ) 
		list = fRecPointsQAList ; 
	else if ( task == AliQA::kESDS )
		list = fESDsQAList ; 

	//DefaultEndOfDetectorCycle(task) ;
	EndOfDetectorCycle(task, list) ;
	TDirectory * subDir = NULL ;
	if (fDetectorDir) 
		subDir = fDetectorDir->GetDirectory(AliQA::GetTaskName(task)) ; 
	if ( subDir ) {
		subDir->cd() ; 
		if (list) 
			list->Write() ;
    if (fObject)
      fObject->Write() ; 
	}
	//Finish() ; 
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
		AliError("Wrong data type") ;     
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
TObjArray *  AliQADataMakerRec::Init(AliQA::TASKINDEX_t task, Int_t run, Int_t cycles)
{
  // general intialisation
	
	TObjArray * rv = NULL ; 
  
	fRun = run ;
	if (cycles > 0)
		SetCycle(cycles) ;  
	
	if ( task == AliQA::kRAWS ) {
		if (! fRawsQAList ) { 
			fRawsQAList = new TObjArray(100) ;	 
      fRawsQAList->SetName(Form("%s_%s", GetName(), AliQA::GetTaskName(task).Data())) ; 
			InitRaws() ;
		}
		rv = fRawsQAList ;
	} else if ( task == AliQA::kRECPOINTS ) {
		if ( ! fRecPointsQAList ) {
			fRecPointsQAList = new TObjArray(100) ; 
      fRecPointsQAList->SetName(Form("%s_%s", GetName(), AliQA::GetTaskName(task).Data())) ; 
			InitRecPoints() ;
		}
		rv = fRecPointsQAList ;
	} else if ( task == AliQA::kESDS ) {
		if ( ! fESDsQAList ) {
			fESDsQAList = new TObjArray(100) ;
      fESDsQAList->SetName(Form("%s_%s", GetName(), AliQA::GetTaskName(task).Data())) ; 
			InitESDs() ;
		}
		rv = fESDsQAList ;
	}
	
	return rv ; 
}

//____________________________________________________________________________ 
void AliQADataMakerRec::Init(AliQA::TASKINDEX_t task, TObjArray * list, Int_t run, Int_t cycles)
{
  // Intialisation by passing the list of QA data booked elsewhere
  
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
void AliQADataMakerRec::StartOfCycle(AliQA::TASKINDEX_t task, const Bool_t sameCycle) 
{ 
  // Finishes a cycle of QA data acquistion
	if ( !sameCycle || fCurrentCycle == -1) {
		ResetCycle() ;
		if (fOutput) 
			fOutput->Close() ; 
		fOutput = AliQA::GetQADataFile(GetName(), fRun, fCurrentCycle) ; 	
	}	
	AliInfo(Form(" Run %d Cycle %d task %s file %s", 
				 fRun, fCurrentCycle, AliQA::GetTaskName(task).Data(), fOutput->GetName() )) ;

	fDetectorDir = fOutput->GetDirectory(GetDetectorDirName()) ; 
	if (!fDetectorDir)
		fDetectorDir = fOutput->mkdir(GetDetectorDirName()) ; 

	TDirectory * subDir = fDetectorDir->GetDirectory(AliQA::GetTaskName(task)) ; 
	if (!subDir)
		subDir = fDetectorDir->mkdir(AliQA::GetTaskName(task)) ;  
	subDir->cd() ; 

	TObjArray * list = NULL ; 
  
  if ( task == AliQA::kRAWS ) 
	  list = fRawsQAList ; 
  else if ( task == AliQA::kRECPOINTS)  
	  list = fRecPointsQAList ;
  else if ( task == AliQA::kESDS )  
	  list = fESDsQAList ;
	
// Should be the choice of detectors
//	TIter next(list) ;
//	TH1 * h ; 
//	while ( (h = dynamic_cast<TH1 *>(next())) )
//		h->Reset() ;  

	StartOfDetectorCycle() ; 
}
