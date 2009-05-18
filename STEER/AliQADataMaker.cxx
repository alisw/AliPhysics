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
AliQADataMaker::AliQADataMaker(const char * name, const char * title) : 
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
  fEventSpecie(AliRecoParam::kDefault)
{
  // ctor
  fDetectorDirName = GetName() ; 
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
    fParameterList[specie] = NULL ; 
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
  fEventSpecie(qadm.fEventSpecie)
{
  //copy ctor
  fDetectorDirName = GetName() ; 
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
	} else if ( index > 10000 ) {
		AliError("Max number of authorized QA objects is 10000") ; 
  } else {    
    if (expert) 
      hist->SetBit(AliQAv1::GetExpertBit()) ;
    if (image) 
      hist->SetBit(AliQAv1::GetImageBit()) ;  
    TH1 * histClone[AliRecoParam::kNSpecies] ; 
    for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
      histClone[specie] = CloneMe(hist, specie) ; 
      histClone[specie]->SetDirectory(0) ; 
      list[specie]->AddAtAndExpand(histClone[specie], index) ; 
      if(saveForCorr) {  
        char * name = Form("%s_%s", list[AliRecoParam::AConvert(AliRecoParam::kDefault)]->GetName(), hist->GetName()) ;  
        TParameter<double> * p = new TParameter<double>(name, 9999.9999) ;
        if ( fParameterList[specie] == NULL )
          fParameterList[specie] = new TList() ; 
        fParameterList[specie]->Add(p) ;
      }
    }
    rv = list[AliRecoParam::kDefault]->GetLast() ;
  }
  delete hist ; 
  return rv ; 
}

//____________________________________________________________________________
TH1 *  AliQADataMaker::CloneMe(TH1 * hist, Int_t specie) const  
{
  // clones a histogram 
  char * name = Form("%s_%s", AliRecoParam::GetEventSpecieName(specie), hist->GetName()) ;
  TH1 * hClone = dynamic_cast<TH1 *>(hist->Clone(name)) ; 
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
	// Returns the QA object at index. Limit is 100. 
  if ( ! list ) {
		AliError("Data list is NULL !!") ; 
		return NULL ; 		
	}

  SetEventSpecie(fEventSpecie) ;  
  if ( GetRecoParam() ) {
    if ( AliRecoParam::Convert(GetRecoParam()->GetEventSpecie()) != AliRecoParam::kDefault) {
      SetEventSpecie(GetRecoParam()->GetEventSpecie()) ; 
    } else { 
      AliError(Form("Event Specie from RecoParam of %s is = %d\n", GetName(), fEventSpecie));
    }
  }
	if (list[AliRecoParam::AConvert(fEventSpecie)]) {
		if ( index > 10000 ) {
			AliError("Max number of authorized QA objects is 10000") ; 
			return NULL ; 
		} else {
      Int_t esindex = AliRecoParam::AConvert(fEventSpecie) ; 
      return list[esindex]->At(index) ; 
		} 	
  } else {
		AliError("Data list is NULL !!") ; 
		return NULL ; 		
	}
}

