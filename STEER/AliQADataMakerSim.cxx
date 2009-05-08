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
#include <TCanvas.h>
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
    for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
      if ( fDigitsQAList[specie]->IsOwner() )
			fDigitsQAList[specie]->Delete() ;
    }
		delete[] fDigitsQAList ;     
	}
	if ( fHitsQAList ) {
    for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
      if ( fHitsQAList[specie]->IsOwner() ) 
			fHitsQAList[specie]->Delete() ;
    }
		delete[] fHitsQAList ;
	}
	if ( fSDigitsQAList ) { 
    for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
      if ( fSDigitsQAList[specie]->IsOwner() ) 
			fSDigitsQAList[specie]->Delete() ; 
    }
		delete[] fSDigitsQAList ; 
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
  TDirectory * subDir = NULL ;
	if (fDetectorDir) 
    subDir = fDetectorDir->GetDirectory(AliQAv1::GetTaskName(task)) ; 
	if (subDir) { 
		subDir->cd() ; 
    for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
      TDirectory * eventSpecieDir = subDir->GetDirectory(AliRecoParam::GetEventSpecieName(specie)) ;
      if (eventSpecieDir) {
        eventSpecieDir->cd() ; 
        TIter next(list[specie]) ; 
        TObject * obj ; 
        while ( (obj = next()) )  {
          if (!obj->TestBit(AliQAv1::GetExpertBit()))
            obj->Write() ;
        }
        if (WriteExpert()) {
          TDirectory * expertDir = eventSpecieDir->GetDirectory(AliQAv1::GetExpert()) ; 
          if ( expertDir ) {
            expertDir->cd() ;
            next.Reset() ; 
            while ( (obj = next()) ) {
              if (!obj->TestBit(AliQAv1::GetExpertBit()))
                continue ; 
            obj->Write() ;
            }      
          }
        }
      }
    }
    fOutput->Save() ; 
  }
  MakeImage(task) ; 
}
 
//____________________________________________________________________________
void AliQADataMakerSim::Exec(AliQAv1::TASKINDEX_t task, TObject * data) 
{ 
  // creates the quality assurance data for the various tasks (Hits, SDigits, Digits, ESDs)
    
	if ( task == AliQAv1::kHITS ) {  
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
	} else if ( task == AliQAv1::kSDIGITS ) {
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
	} else if ( task == AliQAv1::kDIGITS ) {
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
void AliQADataMakerSim::MakeImage(AliQAv1::TASKINDEX_t task)
{
  // create a drawing of detetor defined histograms
  TObjArray ** list = NULL ;  
  switch (task) {
    case AliQAv1::kRAWS:
      break;
    case AliQAv1::kHITS:
      list = fHitsQAList ;
      break;
    case AliQAv1::kSDIGITS:
      list = fSDigitsQAList ;
      break;  
    case AliQAv1::kDIGITS:
      list = fDigitsQAList ;
      break;  
    case AliQAv1::kRECPOINTS:
      break;
    case AliQAv1::kTRACKSEGMENTS:
      break;
    case AliQAv1::kRECPARTICLES:
      break;
    case AliQAv1::kESDS:
      break;
    case AliQAv1::kNTASKINDEX:
      break;
    default:
    break;
  }
  if ( !list) {
    AliFatal("data not initialized, call AliQADataMaker::Init"); 
  return ; 
  }
  TIter next(list[0]) ;  
  TH1 * hdata = NULL ; 
  Int_t nImages = 0 ;
  while ( (hdata=dynamic_cast<TH1 *>(next())) ) {
    if ( hdata->TestBit(AliQAv1::GetImageBit()) )
      nImages++; 
  }
  if ( nImages == 0 ) {
    AliInfo(Form("No histogram will be plotted for %s %s\n", GetName(), AliQAv1::GetTaskName(task).Data())) ;  
  } else {
    AliInfo(Form("%d histograms will be plotted for %s %s\n", nImages, GetName(), AliQAv1::GetTaskName(task).Data())) ;  
    Double_t w  = 1000 ;
    Double_t h  = 1000 ;
    for (Int_t esIndex = 0 ; esIndex < AliRecoParam::kNSpecies ; esIndex++) {
      TCanvas * canvasQA = new TCanvas(Form("QA_%s_%s_%s", 
                                            GetName(), 
                                            AliQAv1::GetTaskName(task).Data(), 
                                            AliRecoParam::GetEventSpecieName(esIndex)), 
                                       Form("QA control plots for det=%s task=%s eventspecie=%s", 
                                            GetName(), 
                                            AliQAv1::GetTaskName(task).Data(), 
                                            AliRecoParam::GetEventSpecieName(esIndex)), 
                                       w, h) ;
      canvasQA->SetWindowSize(w + (w - canvasQA->GetWw()), h + (h - canvasQA->GetWh())) ;
      Int_t nx = TMath::Sqrt(nImages) ; 
      Int_t ny = nx  ; 
      if ( nx < TMath::Sqrt(nImages)) 
        ny++ ; 
      canvasQA->Divide(nx, ny) ; 
      TIter nexthist(list[esIndex]) ; 
      TH1* hist = NULL ;
      Int_t npad = 1 ; 
      canvasQA->cd(npad) ; 
      while ( (hist=dynamic_cast<TH1*>(nexthist())) ) {
        if(hist->TestBit(AliQAv1::GetImageBit())) {
          hist->Draw() ; 
          canvasQA->cd(++npad) ; 
        }
      }
      canvasQA->Print() ; 
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
        fHitsQAList[specie] = new TObjArray(100) ;	 
        fHitsQAList[specie]->SetName(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(task).Data(), AliRecoParam::GetEventSpecieName(specie))) ;
      }
			InitHits() ;
		}
		rv = fHitsQAList ;
	} else if ( task == AliQAv1::kSDIGITS ) {
		if ( ! fSDigitsQAList ) {
      fSDigitsQAList = new TObjArray *[AliRecoParam::kNSpecies] ; 
      for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
        fSDigitsQAList[specie] = new TObjArray(100) ; 
        fSDigitsQAList[specie]->SetName(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(task).Data(), AliRecoParam::GetEventSpecieName(specie))) ; 
      }
      InitSDigits() ;
		}
		rv = fSDigitsQAList ;
   } else if ( task == AliQAv1::kDIGITS ) {
	   if ( ! fDigitsQAList ) {
       fDigitsQAList = new TObjArray *[AliRecoParam::kNSpecies] ; 
       for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {    
         fDigitsQAList[specie] = new TObjArray(100) ;
         fDigitsQAList[specie]->SetName(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(task).Data(), AliRecoParam::GetEventSpecieName(specie))) ;
       }
		   InitDigits() ;
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

	AliInfo(Form(" Run %d Cycle %d task %s file %s", 
				 fRun, fCurrentCycle, AliQAv1::GetTaskName(task).Data(), fOutput->GetName() )) ;

	fDetectorDir = fOutput->GetDirectory(GetDetectorDirName()) ; 
	if (!fDetectorDir)
		fDetectorDir = fOutput->mkdir(GetDetectorDirName()) ; 

	TDirectory * subDir = fDetectorDir->GetDirectory(AliQAv1::GetTaskName(task)) ; 
	if (!subDir)
		subDir = fDetectorDir->mkdir(AliQAv1::GetTaskName(task)) ;  
  
  for ( Int_t index = AliRecoParam::kDefault ; index < AliRecoParam::kNSpecies ; index++ ) {
    TDirectory * eventSpecieDir = subDir->GetDirectory(AliRecoParam::GetEventSpecieName(index)) ; 
    if (!eventSpecieDir) 
      eventSpecieDir = subDir->mkdir(AliRecoParam::GetEventSpecieName(index)) ; 
    TDirectory * expertDir = eventSpecieDir->GetDirectory(AliQAv1::GetExpert()) ; 
    if (!expertDir) 
      expertDir = eventSpecieDir->mkdir(AliQAv1::GetExpert()) ; 
   }   
	StartOfDetectorCycle() ; 
}
