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

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Analysis task for providing various flow informations
// author: D.J. Kim(dong.jo.kim@cern.ch)
// ALICE Group University of Jyvaskyla
// Finland
//
// Fill the analysis containers for ESD or AOD
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////
#include <AliAnalysisManager.h>
#include <AliJBaseTrack.h>
#include "AliJHSCTask.h"

//______________________________________________________________________________
AliJHSCTask::AliJHSCTask() :   
	AliAnalysisTaskSE("JHSCTask"),
	fJCatalystTask(NULL),
	fJCatalystTaskName("JCatalystTask"),
    fFirstEvent(kTRUE),
	fIsMC(kFALSE),
	fTwoMultiAna(NULL),
	fOutput(NULL)
{
}

//______________________________________________________________________________
AliJHSCTask::AliJHSCTask(const char *name):
	AliAnalysisTaskSE(name), 
	fJCatalystTask(NULL),
	fJCatalystTaskName("JCatalystTask"),
    fFirstEvent(kTRUE),
	fIsMC(kFALSE),
	fTwoMultiAna(NULL),
	fOutput(NULL)
{
	// Constructor
	AliInfo("---- AliJHSCTask Constructor ----");
	DefineOutput (1, TList::Class());
}

//____________________________________________________________________________
AliJHSCTask::AliJHSCTask(const AliJHSCTask& ap) :
	AliAnalysisTaskSE(ap.GetName()), 
	fJCatalystTask(ap.fJCatalystTask),
	fJCatalystTaskName(ap.fJCatalystTaskName),
    fFirstEvent(ap.fFirstEvent),
	fIsMC(ap.fIsMC),
	fTwoMultiAna(ap.fTwoMultiAna),
	fOutput( ap.fOutput )
{ 

	AliInfo("----DEBUG AliJHSCTask COPY ----");

}

//_____________________________________________________________________________
AliJHSCTask& AliJHSCTask::operator = (const AliJHSCTask& ap)
{
	// assignment operator

	AliInfo("----DEBUG AliJHSCTask operator= ----");
	this->~AliJHSCTask();
	new(this) AliJHSCTask(ap);
	return *this;
}

//______________________________________________________________________________
AliJHSCTask::~AliJHSCTask()
{
	// destructor 
	//delete fOutput;
	delete fTwoMultiAna;

}

//________________________________________________________________________

void AliJHSCTask::UserCreateOutputObjects()
{  
	//=== create the jcorran outputs objects
	if(fDebug > 1) printf("AliJHSCTask::UserCreateOutPutData() \n");
	//=== Get AnalysisManager
	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
	cout<< "AliJCatalystTask Name = " << fJCatalystTaskName << endl;
	fJCatalystTask = (AliJCatalystTask*)(man->GetTask( fJCatalystTaskName ));
	fTwoMultiAna = new AliAnalysisAnaTwoMultiCorrelations("TwoMultiCorrelations",kFALSE);
	OpenFile(1);

	fTwoMultiAna->UserCreateOutputObjects();
	fFirstEvent = kTRUE;
	PostData(1, fTwoMultiAna->GetTList());
}

//______________________________________________________________________________
void AliJHSCTask::UserExec(Option_t* /*option*/) 
{
	// Processing of one event
	if(fDebug > 5) cout << "------- AliJHSCTask Exec-------"<<endl;
	if(!((Entry()-1)%1))  AliInfo(Form(" Processing event # %lld",  Entry())); 

	if( fJCatalystTask->GetJCatalystEntry() != fEntry ) return;
	if( fFirstEvent ) {
	    fFirstEvent = kFALSE;
	}
	if(fJCatalystTask->GetCentrality()>100. || fJCatalystTask->GetInputList()->GetEntriesFast()<1) return;
	if(fDebug > 5) cout << Form("Nch = %d, cent=%.0f",fJCatalystTask->GetInputList()->GetEntriesFast(), fJCatalystTask->GetCentrality()) << endl;
	fTwoMultiAna->SetInputList( fJCatalystTask->GetInputList() );
	fTwoMultiAna->SetEventCentrality( fJCatalystTask->GetCentrality() );
	fTwoMultiAna->UserExec("");
}

//______________________________________________________________________________
void AliJHSCTask::Init()
{
	// Intialisation of parameters
	AliInfo("Doing initialization") ; 

}

//______________________________________________________________________________
void AliJHSCTask::Terminate(Option_t *)
{
	// Processing when the event loop is ended
	cout<<"AliJHSCTask Analysis DONE !!"<<endl; 
}
