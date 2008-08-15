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

//
//  Base Class
//  Produces the data needed to calculate the quality assurance. 
//  All data must be mergeable objects.
//  Y. Schutz CERN July 2007
//

// --- ROOT system ---
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQADataMakerSim.h"

ClassImp(AliQADataMakerSim)
             
//____________________________________________________________________________ 
AliQADataMakerSim::AliQADataMakerSim(const char * name, const char * title) : 
  AliQADataMaker(name, title), 
  fDigitsQAList(0x0), 
  fHitsQAList(0x0),
  fSDigitsQAList(0x0)
{
	// ctor
	fDetectorDirName = GetName() ; 
}

//____________________________________________________________________________ 
AliQADataMakerSim::AliQADataMakerSim(const AliQADataMakerSim& qadm) :
  AliQADataMaker(qadm.GetName(), qadm.GetTitle()), 
  fDigitsQAList(qadm.fDigitsQAList),
  fHitsQAList(qadm.fHitsQAList),
  fSDigitsQAList(qadm.fSDigitsQAList) 
{
  //copy ctor
  fDetectorDirName = GetName() ; 
}

//____________________________________________________________________________ 
AliQADataMakerSim::~AliQADataMakerSim()
{
	//dtor: delete the TObjArray and thei content
	if ( fDigitsQAList ) { 
		if ( fDigitsQAList->IsOwner() )
			fDigitsQAList->Delete() ;     
		delete fDigitsQAList ;     
	}
	if ( fHitsQAList ) {
		if ( fHitsQAList->IsOwner() ) 
			fHitsQAList->Delete() ;
		delete fHitsQAList ;
	}
	if ( fSDigitsQAList ) { 
		if ( fSDigitsQAList->IsOwner() ) 
			fSDigitsQAList->Delete() ; 
		delete fSDigitsQAList ; 
	}
}

//__________________________________________________________________
AliQADataMakerSim& AliQADataMakerSim::operator = (const AliQADataMakerSim& qadm )
{
  // Assignment operator.
  this->~AliQADataMakerSim();
  new(this) AliQADataMakerSim(qadm);
  return *this;
}

//____________________________________________________________________________
void AliQADataMakerSim::EndOfCycle(AliQA::TASKINDEX_t task) 
{ 
  // Finishes a cycle of QA data acquistion
	TObjArray * list = 0x0 ; 
	
	if ( task == AliQA::kHITS ) 
		list = fHitsQAList ; 
	else if ( task == AliQA::kSDIGITS )
		list = fSDigitsQAList ; 
	else if ( task == AliQA::kDIGITS ) 
		list = fDigitsQAList ; 
	
	EndOfDetectorCycle(task, list) ; 
	TDirectory * subDir = fDetectorDir->GetDirectory(AliQA::GetTaskName(task)) ; 
	if (subDir) { 
		subDir->cd() ; 
		list->Write() ; 
	}
}
 
//____________________________________________________________________________
void AliQADataMakerSim::Exec(AliQA::TASKINDEX_t task, TObject * data) 
{ 
  // creates the quality assurance data for the various tasks (Hits, SDigits, Digits, ESDs)
    
	if ( task == AliQA::kHITS ) {  
		AliDebug(1, "Processing Hits QA") ; 
		TClonesArray * arr = dynamic_cast<TClonesArray *>(data) ; 
		if (arr) { 
			MakeHits(arr) ;
		} else {
			TTree * tree = dynamic_cast<TTree *>(data) ; 
			if (tree) {
				MakeHits(tree) ; 
			} else {
				AliWarning("data are neither a TClonesArray nor a TTree") ; 
			}
		}
	} else if ( task == AliQA::kSDIGITS ) {
		AliDebug(1, "Processing SDigits QA") ; 
		TClonesArray * arr = dynamic_cast<TClonesArray *>(data) ; 
		if (arr) { 
			MakeSDigits(arr) ;
		} else {
			TTree * tree = dynamic_cast<TTree *>(data) ; 
			if (tree) {
				MakeSDigits(tree) ; 
			} else {
				AliWarning("data are neither a TClonesArray nor a TTree") ; 
			}
		}
	} else if ( task == AliQA::kDIGITS ) {
		AliDebug(1, "Processing Digits QA") ; 
		TClonesArray * arr = dynamic_cast<TClonesArray *>(data) ; 
		if (arr) { 
			MakeDigits(arr) ;
		} else {
			TTree * tree = dynamic_cast<TTree *>(data) ; 
			if (tree) {
				MakeDigits(tree) ; 
			} else {
				AliWarning("data are neither a TClonesArray nor a TTree") ; 
			}
		}
	}
}

//____________________________________________________________________________ 
TObjArray *  AliQADataMakerSim::Init(AliQA::TASKINDEX_t task, Int_t run, Int_t cycles)
{
  // general intialisation
	
	fRun = run ;
	if (cycles > 0)
		SetCycle(cycles) ;  
	TObjArray * rv = NULL ; 
	if ( task == AliQA::kHITS ) {
		if ( ! fHitsQAList ) {
			fHitsQAList = new TObjArray(100) ;	 
   fHitsQAList->SetName(Form("%s/%s", GetName(), AliQA::GetTaskName(task).Data())) ; 
			InitHits() ;
		}
		rv = fHitsQAList ;
	} else if ( task == AliQA::kSDIGITS ) {
		if ( ! fSDigitsQAList ) {
			fSDigitsQAList = new TObjArray(100) ; 
   fSDigitsQAList->SetName(Form("%s/%s", GetName(), AliQA::GetTaskName(task).Data())) ; 
			InitSDigits() ;
		}
		rv = fSDigitsQAList ;
   } else if ( task == AliQA::kDIGITS ) {
	   if ( ! fDigitsQAList ) {
		   fDigitsQAList = new TObjArray(100) ;
     fDigitsQAList->SetName(Form("%s/%s", GetName(), AliQA::GetTaskName(task).Data())) ; 
		   InitDigits() ;
	   }
	   rv =  fDigitsQAList ;
   }
  
	return rv ; 
}

//____________________________________________________________________________ 
void AliQADataMakerSim::Init(AliQA::TASKINDEX_t task, TObjArray * list, Int_t run, Int_t cycles)
{
  // Intialisation by passing the list of QA data booked elsewhere
  
	fRun = run ;
	if (cycles > 0)
		SetCycle(cycles) ;  
	
	if ( task == AliQA::kHITS ) {
		fHitsQAList = list ;	 
	} else if ( task == AliQA::kSDIGITS) {
		fSDigitsQAList = list ; 
	} else if ( task == AliQA::kDIGITS ) {
		fDigitsQAList = list ; 
	} 
}

//____________________________________________________________________________
void AliQADataMakerSim::StartOfCycle(AliQA::TASKINDEX_t task, const Bool_t sameCycle) 
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
	
	TObjArray * list = 0x0 ; 
  
	if ( task == AliQA::kHITS )  
		list = fHitsQAList ;
	else if ( task == AliQA::kSDIGITS )  
		list = fSDigitsQAList ;
	else  if ( task == AliQA::kDIGITS ) 
		list = fDigitsQAList ;
	
// Should be the choice of detectors
//	TIter next(list) ;
//	TH1 * h ; 
//	while ( (h = dynamic_cast<TH1 *>(next())) )
//		h->Reset() ;  
//
	StartOfDetectorCycle() ; 
}
