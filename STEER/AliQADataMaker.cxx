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
  fCurrentCycle(qadm.fCurrentCycle), 
  fCycle(qadm.fCycle), 
  fCycleCounter(qadm.fCycleCounter), 
  fRun(qadm.fRun)
{
  //copy ctor
  fDetectorDirName = GetName() ; 
}

//____________________________________________________________________________
Int_t AliQADataMaker::Add2List(TH1 * hist, const Int_t index, TObjArray * list) 
{ 
	// Set histograms memory resident and add to the list
	// Maximm allowed is 10000
	if ( index > 10000 ) {
		AliError("Max number of authorized QA objects is 10000") ; 
		return -1 ; 
	} else {
		hist->SetDirectory(0) ; 
		list->AddAtAndExpand(hist, index) ; 
		return list->GetLast() ;
	}
}

//____________________________________________________________________________ 
void AliQADataMaker::Finish() const 
{ 
  // write to the output File
  fOutput->Close() ; 
} 

//____________________________________________________________________________ 
TObject * AliQADataMaker::GetData(TObjArray * list, const Int_t index)  
{ 
	// Returns the QA object at index. Limit is 100. 
	if ( index > 10000 ) {
		AliError("Max number of authorized QA objects is 10000") ; 
		return NULL; 
	} else {
		return list->At(index) ; 
	}
} 

//____________________________________________________________________________
void AliQADataMaker::Reset(const Bool_t sameCycle) 
{ 
  // Resets defaut value of data members 
	if (!sameCycle) {
		fCycleCounter = 0 ; 
	}
}
