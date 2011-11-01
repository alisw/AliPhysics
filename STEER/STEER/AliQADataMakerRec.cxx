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
#include "AliQAManager.h"
#include "AliESDEvent.h"
#include "AliRawReader.h"

ClassImp(AliQADataMakerRec)
             
//____________________________________________________________________________ 
AliQADataMakerRec::AliQADataMakerRec(const char * name, const char * title) : 
AliQADataMaker(name, title), 
  fDigitsQAList(NULL),
  fESDsQAList(NULL), 
  fRawsQAList(NULL), 
  fRecPointsQAList(NULL),
  fCorrNt(NULL), 
  fRecoParam(NULL),
  fRecPointsArray(NULL)
{
  // ctor
  fDetectorDirName = GetName() ; 
}

//____________________________________________________________________________ 
AliQADataMakerRec::AliQADataMakerRec(const AliQADataMakerRec& qadm) :
  AliQADataMaker(qadm.GetName(), qadm.GetTitle()), 
  fDigitsQAList(qadm.fDigitsQAList),
  fESDsQAList(qadm.fESDsQAList),
  fRawsQAList(qadm.fRawsQAList),
  fRecPointsQAList(qadm.fRecPointsQAList),
  fCorrNt(qadm.fCorrNt),  
  fRecoParam(qadm.fRecoParam),
  fRecPointsArray(NULL)
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
	fESDsQAList[specie]->Delete() ;     
      }
    }
    delete[] fESDsQAList ;
  }
  if ( fRawsQAList ) {
    for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
      if ( fRawsQAList[specie] ) {
        fRawsQAList[specie]->Delete() ;
      }
    }
    delete[] fRawsQAList ;
  }
  if ( fDigitsQAList ) {
    for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
      if ( fDigitsQAList[specie] ) {
        fDigitsQAList[specie]->Delete() ;
      }
    }
    delete[] fDigitsQAList ; 
  }
  if ( fRecPointsQAList ) {
    for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
      if ( fRecPointsQAList[specie] ) {
        fRecPointsQAList[specie]->Delete() ;
      }
    }
    delete[] fRecPointsQAList ; 
  }
  if (fRecPointsArray) {
    fRecPointsArray->Clear() ; 
    delete fRecPointsArray ; 
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
  
  if (! AliQAManager::QAManager(AliQAv1::kRECMODE)->IsSaveData()) return; 

  fDetectorDir = fOutput->GetDirectory(GetDetectorDirName()) ; 
  if (!fDetectorDir) fDetectorDir = fOutput->mkdir(GetDetectorDirName()) ; 
  TDirectory * subDir = fDetectorDir->GetDirectory(AliQAv1::GetTaskName(task)) ; 
  if (!subDir) subDir = fDetectorDir->mkdir(AliQAv1::GetTaskName(task)) ;  
  subDir->cd();
  // 
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) { // skip Default
    if (! AliQAv1::Instance(AliQAv1::GetDetIndex(GetName()))->IsEventSpecieSet(AliRecoParam::ConvertIndex(specie)) || AliRecoParam::ConvertIndex(specie) == AliRecoParam::kDefault) continue ; 
    TDirectory * eventSpecieDir = subDir->GetDirectory(AliRecoParam::GetEventSpecieName(specie)) ;
    if (!eventSpecieDir) eventSpecieDir = subDir->mkdir(AliRecoParam::GetEventSpecieName(specie)) ; 
    eventSpecieDir->cd();    
    if (list) {
      if (list[specie]) {
        TIter next(list[specie]) ; 
        TObject * obj ; 
        while( (obj = next()) ) {
          if (!obj->TestBit(AliQAv1::GetExpertBit())) obj->Write(); // RS: Note, this can be also TObjArray of clones
        }
        if (WriteExpert()) {
          TDirectory * expertDir = eventSpecieDir->GetDirectory(AliQAv1::GetExpert()) ; 
          if (!expertDir)
            expertDir = eventSpecieDir->mkdir(AliQAv1::GetExpert()) ; 
          expertDir->cd() ;
          next.Reset() ; 
          while( (obj = next()) ) {
            if (!obj->TestBit(AliQAv1::GetExpertBit())) continue; 
            obj->Write();  // RS: Note, this can be also TObjArray of clones
          }      
        }
      }
    } 
    else if ( fCorrNt ) {
      if (fCorrNt[specie] && AliQAv1::GetDetIndex(GetName()) == AliQAv1::kCORR) {
        if (fCorrNt[specie]->GetNvar() != 0) {
          eventSpecieDir->cd() ; 
          fCorrNt[specie]->Write() ; 
        }
      }
      fOutput->Save() ; 
    }
  }
}

//____________________________________________________________________________
void AliQADataMakerRec::Exec(AliQAv1::TASKINDEX_t task, TObject * data) 
{ 
  // creates the quality assurance data for the various tasks (Hits, SDigits, Digits, ESDs)
	
  if ( task == AliQAv1::kRAWS ) {
    AliDebug(AliQAv1::GetQADebugLevel(), "Processing Raws QA") ; 
    AliRawReader * rawReader = static_cast<AliRawReader *>(data) ; 
    if (rawReader) 
      MakeRaws(rawReader) ;
    else
      AliDebug(AliQAv1::GetQADebugLevel(), "Raw data are not processed") ;     
  } else if ( task == AliQAv1::kDIGITSR ) {
    AliDebug(AliQAv1::GetQADebugLevel(), "Processing Digits QA") ; 
    TTree * tree = static_cast<TTree *>(data) ; 
    if (strcmp(tree->ClassName(), "TTree") == 0) {
      MakeDigits(tree) ; 
    } else {
      AliWarning("data are not a TTree") ; 
    }
  } else if ( task == AliQAv1::kRECPOINTS ) {
    AliDebug(AliQAv1::GetQADebugLevel(), "Processing RecPoints QA") ; 
    TTree * tree = static_cast<TTree *>(data) ; 
    if (strcmp(tree->ClassName(), "TTree") == 0) {
      MakeRecPoints(tree) ; 
    } else {
      AliWarning("data are not a TTree") ; 
    }
  } else if ( task == AliQAv1::kESDS ) {
    AliDebug(AliQAv1::GetQADebugLevel(), "Processing ESDs QA") ; 
    AliESDEvent * esd = static_cast<AliESDEvent *>(data) ; 
    if (strcmp(esd->ClassName(), "AliESDEvent") == 0) 
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
        fRawsQAList[specie] = new TObjArray(AliQAv1::GetMaxQAObj()) ;	 
        fRawsQAList[specie]->SetName(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(task).Data(), AliRecoParam::GetEventSpecieName(specie))) ;
      }
    }
    rv = fRawsQAList ;
  } else if ( task == AliQAv1::kDIGITSR ) {
    if ( ! fDigitsQAList ) {
      fDigitsQAList = new TObjArray *[AliRecoParam::kNSpecies] ; 
      for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
        fDigitsQAList[specie] = new TObjArray(AliQAv1::GetMaxQAObj()) ; 
        fDigitsQAList[specie]->SetName(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(task).Data(), AliRecoParam::GetEventSpecieName(specie))) ; 
      }
    }
    rv = fDigitsQAList ;
  } else if ( task == AliQAv1::kRECPOINTS ) {
    if ( ! fRecPointsQAList ) {
      fRecPointsQAList = new TObjArray *[AliRecoParam::kNSpecies] ; 
      for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
        fRecPointsQAList[specie] = new TObjArray(AliQAv1::GetMaxQAObj()) ; 
        fRecPointsQAList[specie]->SetName(Form("%s_%s_%s", GetName(), AliQAv1::GetTaskName(task).Data(), AliRecoParam::GetEventSpecieName(specie))) ; 
      }
    }
    rv = fRecPointsQAList ;
  } else if ( task == AliQAv1::kESDS ) {
    if ( ! fESDsQAList ) {
      fESDsQAList = new TObjArray *[AliRecoParam::kNSpecies] ; 
      for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
        fESDsQAList[specie] = new TObjArray(AliQAv1::GetMaxQAObj()) ;
        fESDsQAList[specie]->SetName(Form("%s_%s", GetName(), AliQAv1::GetTaskName(task).Data())); //, AliRecoParam::GetEventSpecieName(specie))) ; 
      }
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
  // Get the recoparam form the OCDB 
  if (!fRecoParam) {
    AliDebug(AliQAv1::GetQADebugLevel(), Form("Loading reconstruction parameter objects for detector %s", GetName()));
    AliCDBPath path(GetName(),"Calib","RecoParam");
    AliCDBEntry *entry=AliCDBManager::Instance()->Get(path.GetPath());    
    if(!entry) {
      fRecoParam = NULL ; 
      AliDebug(AliQAv1::GetQADebugLevel(), Form("Couldn't find RecoParam entry in OCDB for detector %s",GetName()));
    }
    else {
      //      entry->SetOwner(kTRUE);
      TObject * recoParamObj = entry->GetObject() ; 
      if ( strcmp(recoParamObj->ClassName(), "TObjArray") == 0 ) {
        // The detector has only one set of reco parameters
        AliDebug(AliQAv1::GetQADebugLevel(), Form("Array of reconstruction parameters found for detector %s",GetName()));
        TObjArray *recoParamArray = static_cast<TObjArray*>(recoParamObj) ;
	//        recoParamArray->SetOwner(kTRUE);
        for (Int_t iRP=0; iRP<recoParamArray->GetEntriesFast(); iRP++) {
          fRecoParam = static_cast<AliDetectorRecoParam*>(recoParamArray->At(iRP)) ;
          if (!fRecoParam) 
            break ; 
          else if (fRecoParam->IsDefault()) 
            break ; 
        }
      }
      else if (recoParamObj->InheritsFrom("AliDetectorRecoParam")) {
        // The detector has only one set of reco parameters
        // Registering it in AliRecoParam
        AliDebug(AliQAv1::GetQADebugLevel(), Form("Single set of reconstruction parameters found for detector %s",GetName()));
        fRecoParam = static_cast<AliDetectorRecoParam*>(recoParamObj) ;
        static_cast<AliDetectorRecoParam*>(recoParamObj)->SetAsDefault();
      } else { 
        AliError(Form("No valid RecoParam object found in the OCDB for detector %s",GetName()));
      }
    }
  }
}

//____________________________________________________________________________ 
void AliQADataMakerRec::ResetDetector(AliQAv1::TASKINDEX_t task)
{
  // default reset that resets all the QA objects.
  // to be overloaded by detectors, if necessary

  TObjArray ** list = NULL ; 
  if ( task == AliQAv1::kRAWS ) {
    list = fRawsQAList ;	 
  } else if ( task == AliQAv1::kDIGITSR ) {
    list = fDigitsQAList ; 
  } else if ( task == AliQAv1::kRECPOINTS ) {
    list = fRecPointsQAList ; 
  } else if ( task == AliQAv1::kESDS ) {
    list = fESDsQAList ; 
  }
  //list was not initialized, skip
  if (!list) return; 
  
  for (int spec = 0; spec < AliRecoParam::kNSpecies; spec++) {
    if (!AliQAv1::Instance()->IsEventSpecieSet(AliRecoParam::ConvertIndex(spec))) continue;
    TIter next(list[spec]) ; 
    TObject *obj = NULL; 
    while ( (obj = next()) ) {
      if (obj->TestBit(AliQAv1::GetClonedBit())) { // this is array of cloned histos
	TObjArray* arr = (TObjArray*)obj;
	for (int ih=arr->GetEntriesFast();ih--;) {
	  TH1* histo = (TH1*)arr->At(ih); 
	  if (!histo) continue;
	  histo->Reset("ICE");
	  histo->ResetStats();
	}
      }
      else {
	((TH1*)obj)->Reset("ICE");
	((TH1*)obj)->ResetStats();
      }
    }
    ResetEvCountCycle(AliRecoParam::ConvertIndex(spec));
    ResetEvCountTotal(AliRecoParam::ConvertIndex(spec));
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
  if ( run > 0 ) fRun = run ; 
  if ( !sameCycle || fCurrentCycle == -1) {
    ResetCycle() ;
    if (fOutput) 
      fOutput->Close() ; 
    if (AliQAManager::QAManager(AliQAv1::kRECMODE)->IsSaveData())
      fOutput = AliQAv1::GetQADataFile(GetName(), fRun) ; 	
  }	
  AliDebug(AliQAv1::GetQADebugLevel(), Form(" Run %d Cycle %d task %s file %s", 
					    fRun, fCurrentCycle, AliQAv1::GetTaskName(task).Data(), fOutput->GetName() )) ;

  //	fDetectorDir = fOutput->GetDirectory(GetDetectorDirName()) ; 
  //	if (!fDetectorDir)
  //		fDetectorDir = fOutput->mkdir(GetDetectorDirName()) ; 
  //  
  //	TDirectory * subDir = fDetectorDir->GetDirectory(AliQAv1::GetTaskName(task)) ; 
  //	if (!subDir)
  //		subDir = fDetectorDir->mkdir(AliQAv1::GetTaskName(task)) ;  
  //  
  //  for ( Int_t specie = AliRecoParam::kDefault ; specie < AliRecoParam::kNSpecies ; specie++ ) {
  //    TDirectory * eventSpecieDir = subDir->GetDirectory(AliRecoParam::GetEventSpecieName(specie)) ; 
  //    if (!eventSpecieDir) 
  //      eventSpecieDir = subDir->mkdir(AliRecoParam::GetEventSpecieName(specie)) ; 
  //    TDirectory * expertDir = eventSpecieDir->GetDirectory(AliQAv1::GetExpert()) ; 
  //    if (!expertDir)
  //      expertDir = eventSpecieDir->mkdir(AliQAv1::GetExpert()) ; 
  //  } 
  StartOfDetectorCycle();
  ResetEvCountCycle();
}

//____________________________________________________________________________
void AliQADataMakerRec::ClonePerTrigClass(AliQAv1::TASKINDEX_t task)
{
  // clone the histos of the array corresponding to task
  switch (task) 
    {
    case AliQAv1::kRAWS      : ClonePerTrigClassL(fRawsQAList, task);      break;
    case AliQAv1::kDIGITS    : ClonePerTrigClassL(fDigitsQAList, task);    break;
    case AliQAv1::kRECPOINTS : ClonePerTrigClassL(fRecPointsQAList, task); break;
    case AliQAv1::kESDS      : ClonePerTrigClassL(fESDsQAList, task);      break;
    default : AliError(Form("Task %s is invalid in this context", AliQAv1::GetTaskName(task).Data() )); break;
    }
  //
}
