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

#include "AliJJtTask.h"
#include "../AliJCatalystTask.h"
#include "../AliJCard.h"

//______________________________________________________________________________
AliJJtTask::AliJJtTask() :   
	AliAnalysisTaskSE("AliJJtTask"),
	fJCatalystTask(NULL),
	fJCatalystTaskName("JCatalystTask"),
	fJtAna(0x0),
	fOutput(NULL)
{
	DefineOutput (1, TDirectory::Class());
}

//______________________________________________________________________________
AliJJtTask::AliJJtTask(const char *name, TString inputformat):
	AliAnalysisTaskSE(name),
	fJCatalystTask(NULL),
	fJCatalystTaskName("JCatalystTask"),
	fJtAna(0x0),
	fOutput(NULL)
{
	// Constructor
	AliInfo("---- AliJJtTask Constructor ----");

	JUNUSED(inputformat);
	DefineOutput (1, TDirectory::Class());
}

//____________________________________________________________________________
AliJJtTask::AliJJtTask(const AliJJtTask& ap) :
	AliAnalysisTaskSE(ap.GetName()),
	fJCatalystTask(ap.fJCatalystTask),
	fJCatalystTaskName(ap.fJCatalystTaskName),
	fJtAna( ap.fJtAna ),
	fOutput( ap.fOutput )
{ 

	AliInfo("----DEBUG AliJJtTask COPY ----");

}

//_____________________________________________________________________________
AliJJtTask& AliJJtTask::operator = (const AliJJtTask& ap)
{
	// assignment operator

	AliInfo("----DEBUG AliJJtTask operator= ----");
	this->~AliJJtTask();
	new(this) AliJJtTask(ap);
	return *this;
}

//______________________________________________________________________________
AliJJtTask::~AliJJtTask()
{
	// destructor 

	delete fJtAna;

}

//________________________________________________________________________

void AliJJtTask::UserCreateOutputObjects()
{  
	//=== create the jcorran outputs objects
	if(fDebug > 1) printf("AliJJtTask::UserCreateOutPutData() \n");

	//=== Get AnalysisManager
	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();

	fJCatalystTask = (AliJCatalystTask*)(man->GetTask( fJCatalystTaskName ));

	OpenFile(1);
	fOutput = gDirectory;
	fOutput->cd();

	// Order should be kept for correct functionality
	fJtAna->UserCreateOutputObjects();
	fJtAna->SetTrackList(fJCatalystTask->GetInputList());
	fJtAna->GetCard()->WriteCard(fOutput);

	PostData(1, fOutput);

}

//______________________________________________________________________________
void AliJJtTask::UserExec(Option_t* /*option*/) 
{

	// Processing of one event
	if(fDebug > 5) std::cout << "------- AliJJtTask Exec-------"<<endl;
	if( fJCatalystTask->GetJCatalystEntry() != fEntry ) return;
	// Event selection is done in fJCatalystTask
	// need to put cent/vtx/run#
	fJtAna->SetRunNumber(fJCatalystTask->GetRunNumber());
	fJtAna->SetCentrality(fJCatalystTask->GetCentrality());
	fJtAna->SetZVertex(fJCatalystTask->GetZVertex());
	fJtAna->UserExec();
	PostData(1, fOutput );


	if(fDebug > 5) std::cout << "\t------- End UserExec "<<endl;
}

//______________________________________________________________________________
void AliJJtTask::Init()
{
	// Intialisation of parameters
	AliInfo("Doing initialization") ; 
	fJtAna->Init();
}

//______________________________________________________________________________
void AliJJtTask::Terminate(Option_t *)
{

	std::cout<<"AliJJtTask Analysis DONE !!"<<endl;

}
