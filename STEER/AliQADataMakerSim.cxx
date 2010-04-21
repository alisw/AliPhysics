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
  fDigitsQAList(NULL), 
  fHitsQAList(NULL),
  fSDigitsQAList(NULL),  
  fHitsArray(NULL),
  fSDigitsArray(NULL)
{
	// ctor
	fDetectorDirName = GetName() ; 
}

//____________________________________________________________________________ 
AliQADataMakerSim::AliQADataMakerSim(const AliQADataMakerSim& qadm) :
  AliQADataMaker(qadm.GetName(), qadm.GetTitle()), 
  fDigitsQAList(qadm.fDigitsQAList),
  fHitsQAList(qadm.fHitsQAList),
  fSDigitsQAList(qadm.fSDigitsQAList),  
  fHitsArray(NULL),
  fSDigitsArray(NULL)
{
  //copy ctor
  fDetectorDirName = GetName() ; 
}

//____________________________________________________________________________ 
AliQADataMakerSim::~AliQADataMakerSim()
{
	//dtor: delete the TObjArray and thei content
	if ( fDigitsQAList ) { 
    for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
			fDigitsQAList[specie]->Delete() ;
    }
		delete[] fDigitsQAList ;
  }
	if ( fHitsQAList ) {
    for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
			fHitsQAList[specie]->Delete() ;
    }
   	delete[] fHitsQAList ;
  }
	if ( fSDigitsQAList ) { 
    for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
			fSDigitsQAList[specie]->Delete() ; 
    }
 		delete[] fSDigitsQAList ;
  }
  if (fHitsArray) {
    fHitsArray->Clear() ; 
    delete fHitsArray ;
  }
  if (fSDigitsArray) {
    fSDigitsArray->Clear() ; 
    delete fSDigitsArray ;
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
void AliQADataMakerSim::EndOfCycle() 
{ 
  // Finishes a cycle of QA for all tasks
  EndOfCycle(AliQAv1::kHITS) ; 
  EndOfCycle(AliQAv1::kSDIGITS) ; 
  EndOfCycle(AliQAv1::kDIGITS) ;
  ResetCycle() ; 
}

//____________________________________________________________________________
void AliQADataMakerSim::EndOfCycle(AliQAv1::TASKINDEX_t task) 
{ 
  // Finishes a cycle of QA data acquistion
	TObjArray ** list = NULL ; 
	
	if ( task == AliQAv1::kHITS ) 
		list = fHitsQAList ; 
	else if ( task == AliQAv1::kSDIGITS )
		list = fSDigitsQAList ; 
	else if ( task == AliQAv1::kDIGITS ) 
		list = fDigitsQAList ; 
  
  if ( ! list ) 
    return ; 
	EndOfDetectorCycle(task, list) ; 
  fDetectorDir = fOutput->GetDirectory(GetDetectorDirName()) ;
	if (!fDetectorDir) 
    fDetectorDir = fOutput->mkdir(GetDetectorDirName()) ; 
  TDirectory * subDir = fDetectorDir->GetDirectory(AliQAv1::GetTaskName(task)) ; 
  if (!subDir)
    subDir = fDetectorDir->mkdir(AliQAv1::GetTaskName(task)) ;  
  subDir->cd() ; 
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    if (! AliQAv1::Instance(AliQAv1::GetDetIndex(GetName()))->IsEventSpecieSet(AliRecoParam::ConvertIndex(specie)) ) 
      continue ;
    if (list[specie]->GetEntries() != 0 ) {
      TDirectory * eventSpecieDir = subDir->GetDirectory(AliRecoParam::GetEventSpecieName(specie)) ;
      if (!eventSpecieDir) 
        eventSpecieDir = subDir->mkdir(AliRecoParam::GetEventSpecieName(specie)) ; 
      eventSpecieDir->cd() ; 
      TIter next(list[specie]) ; 
      TObject * obj ; 
      while ( (obj = next()) )  {
        if (!obj->TestBit(AliQAv1::GetExpertBit()))
          obj->Write() ;
      }
      if (WriteExpert()) {
        TDirectory * expertDir = eventSpecieDir->GetDirectory(AliQAv1::GetExpert()) ; 
        if (!expertDir) 
          expertDir = eventSpecieDir->mkdir(AliQAv1::GetExpert()) ; 
        expertDir->cd() ;
        next.Reset() ; 
        while ( (obj = next()) ) {
          if (!obj->TestBit(AliQAv1::GetExpertBit()))
            continue ; 
          obj->Write() ;
        }      
      }
    }
    fOutput->Save() ; 
  }
}

//____________________________________________________________________________
void AliQADataMakerSim::Exec(AliQAv1::TASKINDEX_t task, TObject * data) 
{ 
  // creates the quality assurance data for the various tasks (Hits, SDigits, Digits, ESDs)
  
	if ( task == AliQAv1::kHITS ) {  
		AliDebug(AliQAv1::GetQADebugLevel(), "Processing Hits QA") ; 
 		if (strcmp(data->ClassName(), "TClonesArray") == 0) { 
      fHitsArray = static_cast<TClonesArray *>(data) ; 
			MakeHits() ;
		} else if (strcmp(data->ClassName(), "TTree") == 0) {
			TTree * tree = static_cast<TTree *>(data) ; 
      MakeHits(tree) ; 
    } else {
      AliWarning("data are neither a TClonesArray nor a TTree") ; 
    }
	} else if ( task == AliQAv1::kSDIGITS ) {
		AliDebug(AliQAv1::GetQADebugLevel(), "Processing SDigits QA") ; 
		if (strcmp(data->ClassName(), "TClonesArray") == 0) { 
      fSDigitsArray = static_cast<TClonesArray *>(data) ; 
			MakeSDigits() ;
		} else if (strcmp(data->ClassName(), "TTree") == 0) {
			TTree * tree = static_cast<TTree *>(data) ; 
      MakeSDigits(tree) ; 
    } else {
      AliWarning("data are neither a TClonesArray nor a TTree") ; 
    }
 	} else if ( task == AliQAv1::kDIGITS ) {
		AliDebug(AliQAv1::GetQADebugLevel(), "Processing Digits QA") ; 
		if (strcmp(data->ClassName(), "TClonesArray") == 0) { 
      fDigitsArray = static_cast<TClonesArray *>(data) ; 
			MakeDigits() ;
		} else if (strcmp(data->ClassName(), "TTree") == 0)  {
			TTree * tree = static_cast<TTree *>(data) ; 
      MakeDigits(tree) ; 
    } else {
      AliWarning("data are neither a TClonesArray nor a TTree") ; 
    }
  }
}

//____________________________________________________________________________ 
TObjArray **  AliQADataMakerSim::Init(AliQAv1::TASKINDEX_t task, Int_t cycles)
{
  // general intialisation
	
	if (cycles > 0)
		SetCycle(cycles) ;  
	TObjArray ** rv = NULL ; 
	if ( task == AliQAv1::kHITS ) {
		if ( ! fHitsQAList ) {
      fHitsQAList = new TObjArray *[AliRecoParam::kNSpecies] ; 
      for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
        fHitsQAList[specie] = new TObjArray(AliQAv1::GetMaxQAObj()) ;	 
        fHitsQAList[specie]->SetName(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(task).Data(), AliRecoParam::GetEventSpecieName(specie))) ;
      }
		}
		rv = fHitsQAList ;
	} else if ( task == AliQAv1::kSDIGITS ) {
		if ( ! fSDigitsQAList ) {
      fSDigitsQAList = new TObjArray *[AliRecoParam::kNSpecies] ; 
      for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
        fSDigitsQAList[specie] = new TObjArray(AliQAv1::GetMaxQAObj()) ; 
        fSDigitsQAList[specie]->SetName(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(task).Data(), AliRecoParam::GetEventSpecieName(specie))) ; 
      }
		}
		rv = fSDigitsQAList ;
   } else if ( task == AliQAv1::kDIGITS ) {
	   if ( ! fDigitsQAList ) {
       fDigitsQAList = new TObjArray *[AliRecoParam::kNSpecies] ; 
       for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {    
         fDigitsQAList[specie] = new TObjArray(AliQAv1::GetMaxQAObj()) ;
         fDigitsQAList[specie]->SetName(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(task).Data(), AliRecoParam::GetEventSpecieName(specie))) ;
       }
	   }
	   rv =  fDigitsQAList ;
   }
  
	return rv ; 
} 

//____________________________________________________________________________ 
void AliQADataMakerSim::Init(AliQAv1::TASKINDEX_t task, TObjArray ** list, Int_t run, Int_t cycles)
{
  // Intialisation by passing the list of QA data booked elsewhere
  
	fRun = run ;
	if (cycles > 0)
		SetCycle(cycles) ;  
	
	if ( task == AliQAv1::kHITS ) {
		fHitsQAList = list ;	 
	} else if ( task == AliQAv1::kSDIGITS) {
		fSDigitsQAList = list ; 
	} else if ( task == AliQAv1::kDIGITS ) {
		fDigitsQAList = list ; 
	} 
}

//____________________________________________________________________________ 
void AliQADataMakerSim::ResetDetector(AliQAv1::TASKINDEX_t task)
{
    // default reset that resets all the QA objects.
    // to be overloaded by detectors, if necessary
  
  TObjArray ** list = NULL ; 
  if ( task == AliQAv1::kHITS ) {
		list = fHitsQAList ;	 
	} else if ( task == AliQAv1::kSDIGITS ) {
		list = fSDigitsQAList ; 
	} else if ( task == AliQAv1::kDIGITS ) {
		list = fDigitsQAList ; 
	}
    //list was not initialized, skip
  if (!list) 
    return ; 
  
  for (int spec = 0; spec < AliRecoParam::kNSpecies; spec++) {
    if (!AliQAv1::Instance()->IsEventSpecieSet(AliRecoParam::ConvertIndex(spec)))
      continue;
    TIter next(list[spec]) ; 
    TH1 * histo = NULL ; 
    while ( (histo = dynamic_cast<TH1*> (next())) ) {
      histo->Reset() ;
    }
  }
}
  
//____________________________________________________________________________
void AliQADataMakerSim::StartOfCycle(Int_t run) 
{ 
  // Finishes a cycle of QA for all tasks
  Bool_t samecycle = kFALSE ; 
  StartOfCycle(AliQAv1::kHITS,    run, samecycle) ;
  samecycle = kTRUE ; 
  StartOfCycle(AliQAv1::kSDIGITS, run, samecycle) ;
  StartOfCycle(AliQAv1::kDIGITS,  run, samecycle) ;
}

//____________________________________________________________________________
void AliQADataMakerSim::StartOfCycle(AliQAv1::TASKINDEX_t task, Int_t run, const Bool_t sameCycle) 
{ 
  // Finishes a cycle of QA data acquistion
  if ( run > 0 ) 
    fRun = run ; 
	if ( !sameCycle || fCurrentCycle == -1) {
		ResetCycle() ;
	if (fOutput) 
		fOutput->Close() ; 
	fOutput = AliQAv1::GetQADataFile(GetName(), fRun) ; 	
	}	

	AliDebug(AliQAv1::GetQADebugLevel(), Form(" Run %d Cycle %d task %s file %s", 
				 fRun, fCurrentCycle, AliQAv1::GetTaskName(task).Data(), fOutput->GetName() )) ;

	//fDetectorDir = fOutput->GetDirectory(GetDetectorDirName()) ; 
//	if (!fDetectorDir)
//		fDetectorDir = fOutput->mkdir(GetDetectorDirName()) ; 
//
//	TDirectory * subDir = fDetectorDir->GetDirectory(AliQAv1::GetTaskName(task)) ; 
//	if (!subDir)
//		subDir = fDetectorDir->mkdir(AliQAv1::GetTaskName(task)) ;  
//  
//  for ( Int_t index = AliRecoParam::kDefault ; index < AliRecoParam::kNSpecies ; index++ ) {
//    TDirectory * eventSpecieDir = subDir->GetDirectory(AliRecoParam::GetEventSpecieName(index)) ; 
//    if (!eventSpecieDir) 
//      eventSpecieDir = subDir->mkdir(AliRecoParam::GetEventSpecieName(index)) ; 
//    TDirectory * expertDir = eventSpecieDir->GetDirectory(AliQAv1::GetExpert()) ; 
//    if (!expertDir) 
//      expertDir = eventSpecieDir->mkdir(AliQAv1::GetExpert()) ; 
//   }   
	StartOfDetectorCycle() ; 
}
