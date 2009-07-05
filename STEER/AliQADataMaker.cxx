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
#include <TROOT.h> 
#include <TSystem.h> 
#include <TFile.h>
#include <TList.h> 
#include <TTree.h>
#include <TClonesArray.h>
#include <TParameter.h>
#include <TH1K.h>
#include <TH2C.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TH2I.h>
#include <TH3C.h>
#include <TH3D.h>
#include <TH3F.h>
#include <TH3I.h>
#include <TH3S.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQADataMaker.h"
#include "AliQAChecker.h"
#include "AliESDEvent.h"
#include "AliRawReader.h"
#include "AliDetectorRecoParam.h"


ClassImp(AliQADataMaker)
             
//____________________________________________________________________________ 
AliQADataMaker::AliQADataMaker(const Char_t * name, const Char_t * title) : 
  TNamed(name, title), 
  fOutput(0x0),
  fDetectorDir(0x0),
  fDetectorDirName(""), 
  fCurrentCycle(0), 
  fCycle(9999999), 
  fCycleCounter(0), 
  fWriteExpert(kFALSE),
  fParameterList(new TList*[AliRecoParam::kNSpecies]), 
  fRun(0), 
  fEventSpecie(AliRecoParam::kDefault), 
  fDigitsArray(NULL) 
{
  // ctor
  fDetectorDirName = GetName() ; 
    for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
      fParameterList[specie] = NULL ; 
    }
}

//____________________________________________________________________________ 
AliQADataMaker::AliQADataMaker(const AliQADataMaker& qadm) :
  TNamed(qadm.GetName(), qadm.GetTitle()),
  fOutput(qadm.fOutput),
  fDetectorDir(qadm.fDetectorDir),
  fDetectorDirName(qadm.fDetectorDirName),
  fCurrentCycle(qadm.fCurrentCycle), 
  fCycle(qadm.fCycle), 
  fCycleCounter(qadm.fCycleCounter), 
  fWriteExpert(qadm.fWriteExpert),
  fParameterList(qadm.fParameterList),  
  fRun(qadm.fRun), 
  fEventSpecie(qadm.fEventSpecie),
  fDigitsArray(NULL) 
{
  //copy ctor
  fDetectorDirName = GetName() ; 
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    fParameterList[specie] = qadm.fParameterList[specie] ; 
  //  fImage[specie] = qadm.fImage[specie] ; 
  }
}

//____________________________________________________________________________ 
AliQADataMaker::~AliQADataMaker()
{
  for (Int_t esIndex = 0 ; esIndex < AliRecoParam::kNSpecies ; esIndex++) {
    if (fParameterList[esIndex] )
      delete fParameterList[esIndex] ; 
  }
  delete[] fParameterList ; 

  if (fDigitsArray) {
    fDigitsArray->Clear() ; 
    delete fDigitsArray ;
  }
}

//____________________________________________________________________________
Int_t AliQADataMaker::Add2List(TH1 * hist, const Int_t index, TObjArray ** list, const Bool_t expert, const Bool_t image, const Bool_t saveForCorr) 
{ 
	// Set histograms memory resident and add to the list
	// Maximm allowed is 10000
  
  Int_t rv = -1 ; 
  TClass * classType = hist->Class() ;
  TString className(classType->GetName()) ; 
  if( ! className.BeginsWith("T") && ! classType->InheritsFrom("TH1") ) {
    AliError(Form("QA data Object must be a generic ROOT object and derive fom TH1 and not %s", className.Data())) ; 
	} else if ( index > AliQAv1::GetMaxQAObj() ) {
		AliError(Form("Max number of authorized QA objects is %d", AliQAv1::GetMaxQAObj())) ; 
  } else {
    hist->SetDirectory(0) ; 
    if (expert) 
      hist->SetBit(AliQAv1::GetExpertBit()) ;
    if (image) 
      hist->SetBit(AliQAv1::GetImageBit()) ;  
    const Char_t * name = Form("%s_%s", AliRecoParam::GetEventSpecieName(fEventSpecie), hist->GetName()) ;
    hist->SetName(name) ; 
    if(saveForCorr) {  
      const Char_t * cname = Form("%s_%s", list[AliRecoParam::AConvert(AliRecoParam::kDefault)]->GetName(), hist->GetName()) ;  
      TParameter<double> * p = new TParameter<double>(cname, 9999.9999) ;
      if ( fParameterList[AliRecoParam::AConvert(fEventSpecie)] == NULL )
        fParameterList[AliRecoParam::AConvert(fEventSpecie)] = new TList() ; 
      fParameterList[AliRecoParam::AConvert(fEventSpecie)]->Add(p) ;
    }
    list[AliRecoParam::AConvert(fEventSpecie)]->AddAtAndExpand(hist, index) ; 
    rv = list[AliRecoParam::AConvert(fEventSpecie)]->GetLast() ;
  }
  return rv ; 
}

//____________________________________________________________________________
TH1 *  AliQADataMaker::CloneMe(TH1 * hist, Int_t specie) const  
{
  // clones a histogram 
  const Char_t * name = Form("%s_%s", AliRecoParam::GetEventSpecieName(specie), hist->GetName()) ;
  TH1 * hClone = static_cast<TH1 *>(hist->Clone(name)) ; 
  if ( hist->TestBit(AliQAv1::GetExpertBit()) )
    hClone->SetBit(AliQAv1::GetExpertBit()) ; 
  if ( hist->TestBit(AliQAv1::GetImageBit()) )
    hClone->SetBit(AliQAv1::GetImageBit()) ; 
  return hClone ; 
}

//____________________________________________________________________________
void AliQADataMaker::DefaultEndOfDetectorCycle(AliQAv1::TASKINDEX_t task) 
{
	// this method must be oveloaded by detectors
	// sets the QA result to Fatal
	AliQAv1::Instance(AliQAv1::GetDetIndex(GetName())) ;
	AliQAv1 * qa = AliQAv1::Instance(task) ;
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
    qa->Set(AliQAv1::kFATAL, specie) ; 
	AliQAv1::GetQAResultFile()->cd() ; 
	qa->Write(AliQAv1::GetQAName(), kWriteDelete) ;   
	AliQAv1::GetQAResultFile()->Close() ; 
}

//____________________________________________________________________________ 
void AliQADataMaker::Finish() const 
{ 
	// write to the output File
	if (fOutput) 
		fOutput->Close() ; 
} 

//____________________________________________________________________________ 
TObject * AliQADataMaker::GetData(TObjArray ** list, const Int_t index)  
{ 
	// Returns the QA object at index. Limit is AliQAv1::GetMaxQAObj() 
  if ( ! list ) {
		AliError("Data list is NULL !!") ; 
		return NULL ; 		
	}

  SetEventSpecie(fEventSpecie) ;  
  Int_t esindex = AliRecoParam::AConvert(fEventSpecie) ; 
  TH1 * histClone = NULL ; 
  TObjArray * arr = list[esindex] ; 
	if (arr) {
    if ( ! arr->GetEntriesFast() ) {
      // Initializes the histograms 
      TString arrName(arr->GetName()) ; 
      if (arrName.Contains(AliQAv1::GetTaskName(AliQAv1::kRAWS)))
        InitRaws() ; 
      else if (arrName.Contains(AliQAv1::GetTaskName(AliQAv1::kHITS)))
        InitHits() ; 
      else if (arrName.Contains(AliQAv1::GetTaskName(AliQAv1::kSDIGITS)))
        InitSDigits() ; 
      else if (arrName.Contains(AliQAv1::GetTaskName(AliQAv1::kDIGITS)))
        InitDigits() ; 
      else if (arrName.Contains(AliQAv1::GetTaskName(AliQAv1::kDIGITSR)))
        InitDigits() ; 
      else if (arrName.Contains(AliQAv1::GetTaskName(AliQAv1::kRECPOINTS)))
        InitRecPoints() ; 
      else if (arrName.Contains(AliQAv1::GetTaskName(AliQAv1::kESDS)))
        InitESDs() ; 
    }
		if ( index > AliQAv1::GetMaxQAObj() ) {
			AliError(Form("Max number of authorized QA objects is %d", AliQAv1::GetMaxQAObj())) ; 
		} else {
      if ( arr->At(index) )  {
        histClone = static_cast<TH1*>(arr->At(index)) ; 
      } 	
    }
  }
  return histClone ; 		
}

//____________________________________________________________________________ 
TObjArray*  AliQADataMaker::Init(AliQAv1::TASKINDEX_t task, AliRecoParam::EventSpecie_t es, Int_t cycles)
{
  // Initialializes and  returns the QAData list for a given event specie
  TObjArray ** ar = Init(task, cycles) ; 
  return ar[AliRecoParam::AConvert(es)] ;  
}

//____________________________________________________________________________ 
Bool_t AliQADataMaker::IsValidEventSpecie(Int_t eventSpecieIndex, TObjArray ** list)
{
  // check if event specie was present in current run or 
  // if histograms of this event specie have been created
  if (! AliQAv1::Instance(AliQAv1::GetDetIndex(GetName()))->IsEventSpecieSet(AliRecoParam::ConvertIndex(eventSpecieIndex)) || ! list[eventSpecieIndex]->GetEntriesFast() )
    return kFALSE ;
  else
    return kTRUE ;
}
