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

/*
  Base Class
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.
  Y. Schutz CERN July 2007
*/

// --- ROOT system ---
#include <TSystem.h> 
#include <TFile.h>
#include <TList.h> 
#include <TTree.h>
#include <TClonesArray.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQADataMaker.h"
#include "AliQAChecker.h"
#include "AliESDEvent.h"
#include "AliRawReader.h"

ClassImp(AliQADataMaker)
             
//____________________________________________________________________________ 
AliQADataMaker::AliQADataMaker(const char * name, const char * title) : 
  TNamed(name, title), 
  fOutput(0x0),
  fDetectorDir(0x0),
  fDetectorDirName(""), 
  fDigitsQAList(0x0), 
  fESDsQAList(0x0), 
  fHitsQAList(0x0),
  fRawsQAList(0x0), 
  fRecPointsQAList(0x0), 
  fSDigitsQAList(0x0), 
  fCurrentCycle(-1), 
  fCycle(9999999), 
  fCycleCounter(0), 
  fRun(0)
{
  // ctor
  fDetectorDirName = GetName() ; 
}

//____________________________________________________________________________ 
AliQADataMaker::AliQADataMaker(const AliQADataMaker& qadm) :
  TNamed(qadm.GetName(), qadm.GetTitle()),
  fOutput(qadm.fOutput),
  fDetectorDir(qadm.fDetectorDir),
  fDetectorDirName(qadm.fDetectorDirName),
  fDigitsQAList(qadm.fDigitsQAList),
  fESDsQAList(qadm.fESDsQAList),
  fHitsQAList(qadm.fHitsQAList),
  fRawsQAList(qadm.fRecPointsQAList),
  fRecPointsQAList(qadm.fRecPointsQAList),
  fSDigitsQAList(qadm.fSDigitsQAList), 
  fCurrentCycle(qadm.fCurrentCycle), 
  fCycle(qadm.fCycle), 
  fCycleCounter(qadm.fCycleCounter), 
  fRun(qadm.fRun)
{
  //copy ctor
  fDetectorDirName = GetName() ; 
}

//__________________________________________________________________
AliQADataMaker& AliQADataMaker::operator = (const AliQADataMaker& qadm )
{
  // Assignment operator.
  this->~AliQADataMaker();
  new(this) AliQADataMaker(qadm);
  return *this;
}

//____________________________________________________________________________
void AliQADataMaker::EndOfCycle(AliQA::TASKINDEX task) 
{ 
  // Finishes a cycle of QA data acquistion
  
 TList * list = 0x0 ; 
  
 switch (task) { 
  
  case AliQA::kRAWS:    
	list = fRawsQAList ; 
  break ; 

  case AliQA::kHITS:
	list = fHitsQAList ; 
  break ; 

  case AliQA::kSDIGITS:
 	list = fSDigitsQAList ; 
  break ; 
    
  case AliQA::kDIGITS:
 	list = fDigitsQAList ; 
  break ;  
 
   case AliQA::kRECPOINTS:
	list = fRecPointsQAList ; 
   break ;  

   case AliQA::kTRACKSEGMENTS:
   break ;  
  
   case AliQA::kRECPARTICLES:
   break ;  
    
   case AliQA::kESDS:
	list = fESDsQAList ; 
   break ;  
  }	
  
 EndOfDetectorCycle(task, list) ; 
 TDirectory * subDir = fDetectorDir->GetDirectory(AliQA::GetTaskName(task)) ; 
 subDir->cd() ; 
 list->Write() ; 
}
 
//____________________________________________________________________________
void AliQADataMaker::Exec(AliQA::TASKINDEX task, TObject * data) 
{ 
  // creates the quality assurance data for the various tasks (Hits, SDigits, Digits, ESDs)
    
	switch (task) { 
  
		case AliQA::kRAWS:
		{
			AliDebug(1, "Processing Raws QA") ; 
			AliRawReader * rawReader = dynamic_cast<AliRawReader *>(data) ; 
			if (rawReader) 
				MakeRaws(rawReader) ;
			else
			AliError("Wrong data type") ;     
			break ; 
		}
		case AliQA::kHITS:
		{  
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
			break ; 
		}
		case AliQA::kSDIGITS:
		{
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
			break ; 
		}  
		case AliQA::kDIGITS:
		{
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
			break ;  
		}
		case AliQA::kRECPOINTS:
		{
			AliDebug(1, "Processing RecPoints QA") ; 
			TTree * tree = dynamic_cast<TTree *>(data) ; 
			if (tree) {
				MakeRecPoints(tree) ; 
			} else {
				AliWarning("data are not a TTree") ; 
			}
			break ;  
		}
		case AliQA::kTRACKSEGMENTS:
			AliInfo("Processing Track Segments QA: not existing anymore") ; 
		//       MakeTrackSegments(ts) ;
		break ;  
  
		case AliQA::kRECPARTICLES:
			AliInfo("Processing RecParticles QA: not existing anymore") ; 
			//       MakeRecParticles(recpar) ;
		break ;  
		case AliQA::kESDS:
		{
			AliDebug(1, "Processing ESDs QA") ; 
			AliESDEvent * esd = dynamic_cast<AliESDEvent *>(data) ; 
			if (esd) 
				MakeESDs(esd) ;
			else 
				AliError("Wrong type of esd container") ; 
			break ;
		}  
	}	  
}

//____________________________________________________________________________ 
void AliQADataMaker::Finish() const 
{ 
  // write to the output File
  fOutput->Close() ; 
} 

//____________________________________________________________________________ 
TList *  AliQADataMaker::Init(AliQA::TASKINDEX task, Int_t run, Int_t cycles)
{
  // general intialisation
  
  fRun = run ;
  if (cycles > 0)
    SetCycle(cycles) ;  
	
  switch (task) {
  case AliQA::kRAWS: 
   {
	fRawsQAList = new TList() ;	 
    InitRaws() ;
	return fRawsQAList ;
    break ; 
   }
  case AliQA::kHITS: 
   {
	fHitsQAList = new TList() ;	 
    InitHits() ;
	return fHitsQAList ;
    break ; 
   }
  case AliQA::kSDIGITS: 
   {
	fSDigitsQAList = new TList() ; 
    InitSDigits() ;
	return fSDigitsQAList ;
    break ; 
   }
  case AliQA::kDIGITS: 
   {
	fDigitsQAList = new TList(); 
	InitDigits() ;
	return fDigitsQAList ;
	break ; 
   }	  
  case AliQA::kRECPOINTS: 
   {
	fRecPointsQAList = new TList() ; 
    InitRecPoints() ;
	return fRecPointsQAList ;
    break ; 
  }
  case AliQA::kTRACKSEGMENTS: 
//  InitTrackSegments() ;
    break ; 
    
  case AliQA::kRECPARTICLES: 
//    InitRecParticles() ;
    break ; 
    
  case AliQA::kESDS: 
   {
	fESDsQAList = new TList() ; 
	InitESDs() ;
	return fRecPointsQAList ;
    break ; 
   }
  }  
  return 0x0 ; 
}

//____________________________________________________________________________ 
void AliQADataMaker::Init(AliQA::TASKINDEX task, TList * list, Int_t run, Int_t cycles)
{
  // Intialisation by passing the list of QA data booked elsewhere
  
  fRun = run ;
  if (cycles > 0)
    SetCycle(cycles) ;  
	
  switch (task) {
  case AliQA::kRAWS: 
   {
	fRawsQAList = list ;	 
    break ; 
   }
  case AliQA::kHITS: 
   {
	fHitsQAList = list ;	 
    break ; 
   }
  case AliQA::kSDIGITS: 
   {
	fSDigitsQAList = list ; 
    break ; 
   }
  case AliQA::kDIGITS: 
   {
	fDigitsQAList = list ; 
	break ; 
   }	  
  case AliQA::kRECPOINTS: 
   {
	fRecPointsQAList = list ; 
    break ; 
  }
  case AliQA::kTRACKSEGMENTS: 
    break ; 
    
  case AliQA::kRECPARTICLES: 
    break ; 
    
  case AliQA::kESDS: 
   {
	fESDsQAList = list ; 
    break ; 
   }
  }  
}

//____________________________________________________________________________
void AliQADataMaker::Reset() 
{ 
  // Resets defaut value of data members 
  fCurrentCycle = -1 ;  
  fCycleCounter = 0 ; 
}

//____________________________________________________________________________
void AliQADataMaker::StartOfCycle(AliQA::TASKINDEX task, const Bool_t sameCycle) 
{ 
  // Finishes a cycle of QA data acquistion
 
 if ( !sameCycle ) {
	ResetCycle() ;
	if (fOutput) 
		fOutput->Close() ; 
	fOutput = AliQA::GetQADMOutFile(GetName(), fRun, fCurrentCycle) ; 	
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

  TList * list = 0x0 ; 
  
  switch (task) { 
  case AliQA::kRAWS: 
	list = fRawsQAList ; 
    break ; 

  case AliQA::kHITS: 
	list = fHitsQAList ; 
    break ; 
  
  case AliQA::kSDIGITS: 
	list = fSDigitsQAList ;
    break ; 

  case AliQA::kDIGITS: 
	list = fDigitsQAList ;
	break ; 
	  
  case AliQA::kRECPOINTS: 
	list = fRecPointsQAList ;
	break ; 

  case AliQA::kTRACKSEGMENTS: 
    break ; 
    
  case AliQA::kRECPARTICLES: 
    break ; 
    
  case AliQA::kESDS: 
  	list = fESDsQAList ;
    break ; 
  }  

 TIter next(list) ;
 TH1 * h ; 
 while ( (h = dynamic_cast<TH1 *>(next())) )
    h->Reset() ;  

 StartOfDetectorCycle() ; 

}
