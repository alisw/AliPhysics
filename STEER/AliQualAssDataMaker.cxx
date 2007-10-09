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
#include "AliQualAssDataMaker.h"
#include "AliESDEvent.h"
#include "AliRawReader.h"

ClassImp(AliQualAssDataMaker)
  
TString AliQualAssDataMaker::fDetectorDirName("") ;

           
//____________________________________________________________________________ 
AliQualAssDataMaker::AliQualAssDataMaker(const char * name, const char * title) : 
  TNamed(name, title), 
  fOutput(0x0),
  fDetectorDir(0x0),
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
AliQualAssDataMaker::AliQualAssDataMaker(const AliQualAssDataMaker& qadm) :
  TNamed(qadm.GetName(), qadm.GetTitle()),
  fOutput(qadm.fOutput),
  fDetectorDir(qadm.fDetectorDir),
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
AliQualAssDataMaker& AliQualAssDataMaker::operator = (const AliQualAssDataMaker& qadm )
{
  // Equal operator.
  this->~AliQualAssDataMaker();
  new(this) AliQualAssDataMaker(qadm);
  return *this;
}

//____________________________________________________________________________
void AliQualAssDataMaker::EndOfCycle(AliQualAss::TASKINDEX task) 
{ 
  // Finishes a cycle of QA data acquistion
 
 EndOfDetectorCycle() ; 
 TDirectory * subDir = fDetectorDir->GetDirectory(AliQualAss::GetTaskName(task)) ; 
 
 switch (task) { 
  
  case AliQualAss::kRAWS:
    subDir->cd() ; 
	fRawsQAList->Write() ; 
  break ; 

  case AliQualAss::kHITS:
    subDir->cd() ; 
	fHitsQAList->Write() ; 
  break ; 

  case AliQualAss::kSDIGITS:
    subDir->cd() ; 
	fSDigitsQAList->Write() ; 
  break ; 
    
  case AliQualAss::kDIGITS:
    subDir->cd() ; 
	fDigitsQAList->Write() ; 
  break ;  
 
   case AliQualAss::kRECPOINTS:
    subDir->cd() ; 
	fRecPointsQAList->Write() ; 
   break ;  

   case AliQualAss::kTRACKSEGMENTS:
   break ;  
  
   case AliQualAss::kRECPARTICLES:
   break ;  
    
   case AliQualAss::kESDS:
    subDir->cd() ; 
	fESDsQAList->Write() ; 
   break ;  
  }	
}
 
//____________________________________________________________________________
void AliQualAssDataMaker::Exec(AliQualAss::TASKINDEX task, TObject * data) 
{ 
  // creates the quality assurance data for the various tasks (Hits, SDigits, Digits, ESDs)
    
  switch (task) { 
  
  case AliQualAss::kRAWS:
  {
    AliInfo("Processing Raws QA") ; 
	AliRawReader * rawReader = dynamic_cast<AliRawReader *>(data) ; 
    if (rawReader) 
	  MakeRaws(rawReader) ;
	 else
	  AliError("Wrong data type") ;     
    break ; 
  }
  case AliQualAss::kHITS:
  {  
	AliInfo("Processing Hits QA") ; 
	TClonesArray * hits = dynamic_cast<TClonesArray *>(data) ; 
    if (hits) 
     MakeHits(hits) ;
	else 
     AliError("Wrong type of hits container") ;
    break ; 
  }
  case AliQualAss::kSDIGITS:
  {
    AliInfo("Processing SDigits QA") ; 
    TClonesArray * sdigits = dynamic_cast<TClonesArray *>(data) ; 
	if (sdigits) 
      MakeSDigits(sdigits) ;
	 else    
      AliError("Wrong type of sdigits container") ; 
    break ; 
  }  
  case AliQualAss::kDIGITS:
  {
    TClonesArray * digits = dynamic_cast<TClonesArray *>(data) ; 
    if (digits) 
	  MakeDigits(digits) ;
	 else 
      AliError("Wrong type of digits container") ; 
    break ;  
  }
  case AliQualAss::kRECPOINTS:
  {
     AliInfo("Processing RecPoints QA") ; 
     TTree * recpoints = dynamic_cast<TTree *>(data) ; 
    if (recpoints) 
      MakeRecPoints(recpoints) ;
    else 
      AliError("Wrong type of recpoints container") ; 
    break ;  
  }
   case AliQualAss::kTRACKSEGMENTS:
    AliInfo("Processing Track Segments QA: not existing anymore") ; 
//     TTree * ts = dynamic_cast<TTree *>(data) ; 
//     if (ts) 
//       MakeTrackSegments(ts) ;
//     else 
//       AliError("Wrong type of track segments container") ; 
    break ;  
  
    case AliQualAss::kRECPARTICLES:
    AliInfo("Processing RecParticles QA: not existing anymore") ; 
//     TTree * recpar = dynamic_cast<TTree *>(data) ; 
//     if (recpar) 
//       MakeRecParticles(recpar) ;
//     else 
//       AliError("Wrong type of recparticles container") ; 
    break ;  
    
  case AliQualAss::kESDS:
   {
    AliInfo("Processing ESDs QA") ; 
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
void AliQualAssDataMaker::Finish(AliQualAss::TASKINDEX) const 
{ 
  // write to the output File
  fOutput->Close() ; 
} 

//____________________________________________________________________________ 
TList *  AliQualAssDataMaker::Init(AliQualAss::TASKINDEX task, Int_t run, Int_t cycles)
{
  // general intialisation
  
  fRun = run ;
  if (cycles > 0)
    SetCycle(cycles) ;  
	
  switch (task) {
  case AliQualAss::kRAWS: 
   {
	fRawsQAList = new TList() ;	 
    InitRaws() ;
	return fRawsQAList ;
    break ; 
   }
  case AliQualAss::kHITS: 
   {
	fHitsQAList = new TList() ;	 
    InitHits() ;
	return fHitsQAList ;
    break ; 
   }
  case AliQualAss::kSDIGITS: 
   {
	fSDigitsQAList = new TList() ; 
    InitSDigits() ;
	return fSDigitsQAList ;
    break ; 
   }
  case AliQualAss::kDIGITS: 
   {
	fDigitsQAList = new TList(); 
	InitDigits() ;
	return fDigitsQAList ;
	break ; 
   }	  
  case AliQualAss::kRECPOINTS: 
   {
	fRecPointsQAList = new TList ; 
    InitRecPoints() ;
	return fRecPointsQAList ;
    break ; 
  }
  case AliQualAss::kTRACKSEGMENTS: 
//  InitTrackSegments() ;
    break ; 
    
  case AliQualAss::kRECPARTICLES: 
//    InitRecParticles() ;
    break ; 
    
  case AliQualAss::kESDS: 
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
void AliQualAssDataMaker::StartOfCycle(AliQualAss::TASKINDEX task, Option_t * sameCycle) 
{ 
  // Finishes a cycle of QA data acquistion
 
 if ( (strcmp(sameCycle, "new") == 0) )  {
   ResetCycle() ;
   if (fOutput) 
	fOutput->Close() ; 
   fOutput = AliQualAss::GetQADMOutFile(GetName(), fRun, fCurrentCycle) ; 	
 }
    	
 AliInfo(Form(" Run %d Cycle %d task %s file %s", 
	fRun, fCurrentCycle, AliQualAss::GetTaskName(task).Data(), fOutput->GetName() )) ;

 fDetectorDir = fOutput->GetDirectory(GetDetectorDirName()) ; 
 if (!fDetectorDir)
   fDetectorDir = fOutput->mkdir(GetDetectorDirName()) ; 

 TDirectory * subDir = fDetectorDir->GetDirectory(AliQualAss::GetTaskName(task)) ; 
 if (!subDir)
   subDir = fDetectorDir->mkdir(AliQualAss::GetTaskName(task)) ;  
 subDir->cd() ; 

  TList * list = 0x0 ; 
  
  switch (task) { 
  case AliQualAss::kRAWS: 
	list = fRawsQAList ; 
    break ; 

  case AliQualAss::kHITS: 
	list = fHitsQAList ; 
    break ; 
  
  case AliQualAss::kSDIGITS: 
	list = fSDigitsQAList ;
    break ; 

  case AliQualAss::kDIGITS: 
	list = fDigitsQAList ;
	break ; 
	  
  case AliQualAss::kRECPOINTS: 
	list = fRecPointsQAList ;
	break ; 

  case AliQualAss::kTRACKSEGMENTS: 
    break ; 
    
  case AliQualAss::kRECPARTICLES: 
    break ; 
    
  case AliQualAss::kESDS: 
  	list = fESDsQAList ;
    break ; 
  }  

 TIter next(list) ;
 TH1 * h ; 
 while ( (h = dynamic_cast<TH1 *>(next())) )
    h->Reset() ;  

 StartOfDetectorCycle() ; 

}
