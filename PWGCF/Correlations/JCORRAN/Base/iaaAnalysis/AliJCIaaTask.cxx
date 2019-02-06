/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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
//______________________________________________________________________________
// Analysis task for high pt particle correlations 
// author: Jussi Viinikainen, D. J.~Kim
// ALICE Group University of Jyvaskyla 
// Finland 
//////////////////////////////////////////////////////////////////////////////

#include <AliAnalysisTaskSE.h>
#include <AliAODHandler.h>
#include <AliAnalysisManager.h>
#include <AliLog.h>

#include "AliJCIaaTask.h"
#include "../AliJCatalystTask.h"
#include "../AliJCard.h"

//______________________________________________________________________________
AliJCIaaTask::AliJCIaaTask() :   
	AliAnalysisTaskSE("AliJCIaaTask"),
	fJCatalystTask(NULL),
	fJCatalystTaskName("JCatalystTask"),
	fIaaAna(0x0),
	fOutput(NULL)
{
	DefineOutput (1, TDirectory::Class());
}

//______________________________________________________________________________
AliJCIaaTask::AliJCIaaTask(const char *name, TString inputformat):
	AliAnalysisTaskSE(name),
	fJCatalystTask(NULL),
	fJCatalystTaskName("JCatalystTask"),
	fIaaAna(0x0),
	fOutput(NULL)
{
	// Constructor
	AliInfo("---- AliJCIaaTask Constructor ----");

	JUNUSED(inputformat);
	DefineOutput (1, TDirectory::Class());
}

//____________________________________________________________________________
AliJCIaaTask::AliJCIaaTask(const AliJCIaaTask& ap) :
	AliAnalysisTaskSE(ap.GetName()),
	fJCatalystTask(ap.fJCatalystTask),
	fJCatalystTaskName(ap.fJCatalystTaskName),
	fIaaAna( ap.fIaaAna ),
	fOutput( ap.fOutput )
{ 

	AliInfo("----DEBUG AliJCIaaTask COPY ----");

}

//_____________________________________________________________________________
AliJCIaaTask& AliJCIaaTask::operator = (const AliJCIaaTask& ap)
{
	// assignment operator

	AliInfo("----DEBUG AliJCIaaTask operator= ----");
	this->~AliJCIaaTask();
	new(this) AliJCIaaTask(ap);
	return *this;
}

//______________________________________________________________________________
AliJCIaaTask::~AliJCIaaTask()
{
	// destructor 

	delete fIaaAna;

}

//________________________________________________________________________

void AliJCIaaTask::UserCreateOutputObjects()
{  
	//=== create the jcorran outputs objects
	if(fDebug > 1) printf("AliJCIaaTask::UserCreateOutPutData() \n");

	//=== Get AnalysisManager
	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();

	fJCatalystTask = (AliJCatalystTask*)(man->GetTask( fJCatalystTaskName ));

	OpenFile(1);
	fOutput = gDirectory;
	fOutput->cd();

	// Order should be kept for correct functionality
	fIaaAna->UserCreateOutputObjects();
	fIaaAna->SetTrackList(fJCatalystTask->GetInputList());
	fIaaAna->GetCard()->WriteCard(fOutput);

	PostData(1, fOutput);

}

//______________________________________________________________________________
void AliJCIaaTask::UserExec(Option_t* /*option*/) 
{

	// Processing of one event
	if(fDebug > 5) std::cout << "------- AliJCIaaTask Exec-------"<<endl;
	if( fJCatalystTask->GetJCatalystEntry() != fEntry ) return;
	// Event selection is done in fJCatalystTask
	// need to put cent/vtx/run#
	fIaaAna->SetRunNumber(fJCatalystTask->GetRunNumber());
	fIaaAna->SetCentrality(fJCatalystTask->GetCentrality());
	fIaaAna->SetZVertex(fJCatalystTask->GetZVertex());
	fIaaAna->UserExec();
	PostData(1, fOutput );


	if(fDebug > 5) std::cout << "\t------- End UserExec "<<endl;
}

//______________________________________________________________________________
void AliJCIaaTask::Init()
{
	// Intialisation of parameters
	AliInfo("Doing initialization") ; 
	fIaaAna->Init();
}

//______________________________________________________________________________
void AliJCIaaTask::Terminate(Option_t *)
{

	std::cout<<"AliJCIaaTask Analysis DONE !!"<<endl;

}
