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
#include <TCanvas.h>
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
  fDigitsQAList(NULL),
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
  fDigitsQAList(qadm.fDigitsQAList),
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
	if ( fDigitsQAList ) {
    for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
      if ( fDigitsQAList[specie] ) {
        if ( fDigitsQAList[specie]->IsOwner() ) 
          fDigitsQAList[specie]->Delete() ;
      }
    }
		delete[] fDigitsQAList ; 
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
  EndOfCycle(AliQAv1::kRAWS) ; 
  EndOfCycle(AliQAv1::kDIGITSR) ; 
  EndOfCycle(AliQAv1::kRECPOINTS) ; 
  EndOfCycle(AliQAv1::kESDS) ; 
  ResetCycle() ; 
}

//____________________________________________________________________________
void AliQADataMakerRec::EndOfCycle(AliQAv1::TASKINDEX_t task) 
{
	// Finishes a cycle of QA 
	
	TObjArray ** list = NULL ; 
	
	if ( task == AliQAv1::kRAWS )     
		list = fRawsQAList ; 
	else if ( task == AliQAv1::kDIGITSR ) 
		list = fDigitsQAList ; 
	else if ( task == AliQAv1::kRECPOINTS ) 
		list = fRecPointsQAList ; 
	else if ( task == AliQAv1::kESDS )
		list = fESDsQAList ; 

 
	if ( ! list && ! fCorrNt ) 
    return ; 
  //DefaultEndOfDetectorCycle(task) ;
	EndOfDetectorCycle(task, list) ;
	TDirectory * subDir = NULL ;
	if (fDetectorDir) 
		subDir = fDetectorDir->GetDirectory(AliQAv1::GetTaskName(task)) ; 
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
            if (!obj->TestBit(AliQAv1::GetExpertBit()))
              obj->Write() ;
          }
          if (WriteExpert()) {
            TDirectory * expertDir = eventSpecieDir->GetDirectory(AliQAv1::GetExpert()) ; 
            if ( expertDir ) { // Write only if requested
              expertDir->cd() ;
              next.Reset() ; 
              while( (obj = next()) ) {
                if (!obj->TestBit(AliQAv1::GetExpertBit()))
                  continue ; 
                obj->Write() ;
              }      
            }
          }
        }
        if ( !fCorrNt )
          continue ; 
        if (fCorrNt[specie] && AliQAv1::GetDetIndex(GetName()) == AliQAv1::kCORR) {
          eventSpecieDir->cd() ; 
          fCorrNt[specie]->Write() ; 
        }
      }
    }
    fOutput->Save() ; 
	}
  if ( AliDebugLevel()  == AliQAv1::GetQADebugLevel() ) 
    MakeImage(task) ; 
}

//____________________________________________________________________________ 
void AliQADataMakerRec::MakeImage(AliQAv1::TASKINDEX_t task)
{
  // create a drawing of detetor defined histograms
  TObjArray ** list = NULL ;  
  switch (task) {
    case AliQAv1::kRAWS:
      list = fRawsQAList ; 
      break;
    case AliQAv1::kHITS:
      break;
    case AliQAv1::kSDIGITS:
      break;  
    case AliQAv1::kDIGITS:
      break;  
    case AliQAv1::kDIGITSR:
      list = fDigitsQAList ; 
      break;  
    case AliQAv1::kRECPOINTS:
      list = fRecPointsQAList ; 
      break;
    case AliQAv1::kTRACKSEGMENTS:
      break;
    case AliQAv1::kRECPARTICLES:
      break;
    case AliQAv1::kESDS:
      list = fESDsQAList ; 
      break;
    case AliQAv1::kNTASKINDEX:
      break;
    default:
      break;
  }
  if ( !list) {
    AliError("data not initialized, call AliQADataMaker::Init"); 
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
    AliWarning(Form("No histogram will be plotted for %s %s\n", GetName(), AliQAv1::GetTaskName(task).Data())) ;  
  } else {
    AliDebug(AliQAv1::GetQADebugLevel(), Form("%d histograms will be plotted for %s %s\n", nImages, GetName(), AliQAv1::GetTaskName(task).Data())) ;  
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
void AliQADataMakerRec::Exec(AliQAv1::TASKINDEX_t task, TObject * data) 
{ 
  // creates the quality assurance data for the various tasks (Hits, SDigits, Digits, ESDs)
	
	if ( task == AliQAv1::kRAWS ) {
		AliDebug(AliQAv1::GetQADebugLevel(), "Processing Raws QA") ; 
		AliRawReader * rawReader = dynamic_cast<AliRawReader *>(data) ; 
		if (rawReader) 
			MakeRaws(rawReader) ;
		else
      AliDebug(AliQAv1::GetQADebugLevel(), "Raw data are not processed") ;     
	} else if ( task == AliQAv1::kDIGITSR ) {
		AliDebug(AliQAv1::GetQADebugLevel(), "Processing Digits QA") ; 
		TTree * tree = dynamic_cast<TTree *>(data) ; 
		if (tree) {
			MakeDigits(tree) ; 
		} else {
			AliWarning("data are not a TTree") ; 
		}
	} else if ( task == AliQAv1::kRECPOINTS ) {
		AliDebug(AliQAv1::GetQADebugLevel(), "Processing RecPoints QA") ; 
		TTree * tree = dynamic_cast<TTree *>(data) ; 
		if (tree) {
			MakeRecPoints(tree) ; 
		} else {
			AliWarning("data are not a TTree") ; 
		}
	} else if ( task == AliQAv1::kESDS ) {
		AliDebug(AliQAv1::GetQADebugLevel(), "Processing ESDs QA") ; 
		AliESDEvent * esd = dynamic_cast<AliESDEvent *>(data) ; 
		if (esd) 
			MakeESDs(esd) ;
		else 
			AliError("Wrong type of esd container") ; 
	}  
}

//____________________________________________________________________________ 
TObjArray **  AliQADataMakerRec::Init(AliQAv1::TASKINDEX_t task, Int_t cycles)
{
  // general intialisation
  InitRecoParams() ;
	TObjArray ** rv = NULL ; 
  
	if (cycles > 0)
		SetCycle(cycles) ;  
	
	if ( task == AliQAv1::kRAWS ) {
		if (! fRawsQAList ) { 
      fRawsQAList = new TObjArray *[AliRecoParam::kNSpecies] ; 
      for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
        fRawsQAList[specie] = new TObjArray(100) ;	 
        fRawsQAList[specie]->SetName(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(task).Data(), AliRecoParam::GetEventSpecieName(specie))) ;
      }
			InitRaws() ;
		}
		rv = fRawsQAList ;
	} else if ( task == AliQAv1::kDIGITSR ) {
		if ( ! fDigitsQAList ) {
      fDigitsQAList = new TObjArray *[AliRecoParam::kNSpecies] ; 
      for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
        fDigitsQAList[specie] = new TObjArray(100) ; 
        fDigitsQAList[specie]->SetName(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(task).Data(), AliRecoParam::GetEventSpecieName(specie))) ; 
      }
      InitDigits() ;
		}
		rv = fDigitsQAList ;
	} else if ( task == AliQAv1::kRECPOINTS ) {
		if ( ! fRecPointsQAList ) {
      fRecPointsQAList = new TObjArray *[AliRecoParam::kNSpecies] ; 
      for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
        fRecPointsQAList[specie] = new TObjArray(100) ; 
        fRecPointsQAList[specie]->SetName(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(task).Data(), AliRecoParam::GetEventSpecieName(specie))) ; 
      }
      InitRecPoints() ;
		}
		rv = fRecPointsQAList ;
	} else if ( task == AliQAv1::kESDS ) {
		if ( ! fESDsQAList ) {
      fESDsQAList = new TObjArray *[AliRecoParam::kNSpecies] ; 
      for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
        fESDsQAList[specie] = new TObjArray(100) ;
        fESDsQAList[specie]->SetName(Form("%s_%s", GetName(), AliQAv1::GetTaskName(task).Data(), AliRecoParam::GetEventSpecieName(specie))) ; 
      }
			InitESDs() ;
		}
		rv = fESDsQAList ;
	}
	return rv ; 
}

//____________________________________________________________________________ 
void AliQADataMakerRec::Init(AliQAv1::TASKINDEX_t task, TObjArray ** list, Int_t run, Int_t cycles)
{
  // Intialisation by passing the list of QA data booked elsewhere
  
  InitRecoParams() ;
  fRun = run ;
	if (cycles > 0)
		SetCycle(cycles) ;  
	
	if ( task == AliQAv1::kRAWS ) {
		fRawsQAList = list ;	 
	} else if ( task == AliQAv1::kDIGITSR ) {
		fDigitsQAList = list ; 
	} else if ( task == AliQAv1::kRECPOINTS ) {
		fRecPointsQAList = list ; 
	} else if ( task == AliQAv1::kESDS ) {
		fESDsQAList = list ; 
	}
}

//____________________________________________________________________________
void AliQADataMakerRec::InitRecoParams() 
{
  if (!fRecoParam) {
    AliDebug(AliQAv1::GetQADebugLevel(), Form("Loading reconstruction parameter objects for detector %s", GetName()));
    AliCDBPath path(GetName(),"Calib","RecoParam");
    AliCDBEntry *entry=AliCDBManager::Instance()->Get(path.GetPath());
    if(!entry) {
      fRecoParam = NULL ; 
      AliWarning(Form("Couldn't find RecoParam entry in OCDB for detector %s",GetName()));
    }
    else {
      TObject * recoParamObj = entry->GetObject() ; 
      if (dynamic_cast<TObjArray*>(recoParamObj)) {
        // The detector has only one set of reco parameters
        AliDebug(AliQAv1::GetQADebugLevel(), Form("Array of reconstruction parameters found for detector %s",GetName()));
        TObjArray *recoParamArray = dynamic_cast<TObjArray*>(recoParamObj) ;
        for (Int_t iRP=0; iRP<recoParamArray->GetEntriesFast(); iRP++) {
          fRecoParam = dynamic_cast<AliDetectorRecoParam*>(recoParamArray->At(iRP)) ;
          if (fRecoParam->IsDefault()) break;
        }
      }
      else if (dynamic_cast<AliDetectorRecoParam*>(recoParamObj)) {
        // The detector has only onse set of reco parameters
        // Registering it in AliRecoParam
        AliDebug(AliQAv1::GetQADebugLevel(), Form("Single set of reconstruction parameters found for detector %s",GetName()));
        dynamic_cast<AliDetectorRecoParam*>(recoParamObj)->SetAsDefault();
        fRecoParam = dynamic_cast<AliDetectorRecoParam*>(recoParamObj) ;
      } else { 
        AliError(Form("No valid RecoParam object found in the OCDB for detector %s",GetName()));
      }
    }
    AliCDBManager::Instance()->UnloadFromCache(path.GetPath());
  }
}

//____________________________________________________________________________
void AliQADataMakerRec::StartOfCycle(Int_t run) 
{
  // Finishes a cycle of QA for all the tasks
  Bool_t samecycle = kFALSE ; 
  StartOfCycle(AliQAv1::kRAWS,      run, samecycle) ;
  samecycle = kTRUE ; 
  StartOfCycle(AliQAv1::kDIGITSR,   run, samecycle) ; 
  StartOfCycle(AliQAv1::kRECPOINTS, run, samecycle) ; 
  StartOfCycle(AliQAv1::kESDS,      run, samecycle) ; 
}

//____________________________________________________________________________
void AliQADataMakerRec::StartOfCycle(AliQAv1::TASKINDEX_t task, Int_t run, const Bool_t sameCycle) 
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

	fDetectorDir = fOutput->GetDirectory(GetDetectorDirName()) ; 
	if (!fDetectorDir)
		fDetectorDir = fOutput->mkdir(GetDetectorDirName()) ; 

	TDirectory * subDir = fDetectorDir->GetDirectory(AliQAv1::GetTaskName(task)) ; 
	if (!subDir)
		subDir = fDetectorDir->mkdir(AliQAv1::GetTaskName(task)) ;  
  
  for ( Int_t specie = AliRecoParam::kDefault ; specie < AliRecoParam::kNSpecies ; specie++ ) {
    TDirectory * eventSpecieDir = subDir->GetDirectory(AliRecoParam::GetEventSpecieName(specie)) ; 
    if (!eventSpecieDir) 
      eventSpecieDir = subDir->mkdir(AliRecoParam::GetEventSpecieName(specie)) ; 
    TDirectory * expertDir = eventSpecieDir->GetDirectory(AliQAv1::GetExpert()) ; 
    if (!expertDir)
      expertDir = eventSpecieDir->mkdir(AliQAv1::GetExpert()) ; 
  } 
	StartOfDetectorCycle() ; 
}
