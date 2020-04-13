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
// author: Marton Vargyas, Jussi Viinikainen
// ALICE Group University of Jyvaskyla 
// Finland 
// Fill the analysis containers for ESD or AOD
// Adapted for AliAnalysisTaskSE and AOD objects  
//////////////////////////////////////////////////////////////////////////////

#include <AliAnalysisTaskSE.h>
#include <AliAODHandler.h>
#include <AliAnalysisManager.h>
#include <AliLog.h>

#include "AliJDiHadronIaaTask.h"
#include "../AliJEventHeader.h"
#include "../AliJRunHeader.h"
#include "../AliJCORRANTask.h"
#include "../AliJCard.h"

//______________________________________________________________________________
AliJDiHadronIaaTask::AliJDiHadronIaaTask() :
	AliAnalysisTaskSE("AliJDiHadronIaaTaskTask"),
	fFilterTask(NULL),
	fFilterTaskName(""),
	fIaaAnalysis(0x0),
	fOutput(NULL)
{
	DefineOutput (1, TDirectory::Class());
}

//______________________________________________________________________________
AliJDiHadronIaaTask::AliJDiHadronIaaTask(const char *name, TString inputformat):
	AliAnalysisTaskSE(name),
	fFilterTask(NULL),
	fFilterTaskName(""),
	fIaaAnalysis(0x0),
	fOutput(NULL)
{
	// Constructor
	AliInfo("---- AliJDiHadronIaaTask Constructor ----");

	JUNUSED(inputformat);
	DefineOutput (1, TDirectory::Class());
}

//____________________________________________________________________________
AliJDiHadronIaaTask::AliJDiHadronIaaTask(const AliJDiHadronIaaTask& ap) :
	AliAnalysisTaskSE(ap.GetName()),
	fFilterTask(ap.fFilterTask),
	fFilterTaskName(ap.fFilterTaskName),
	fIaaAnalysis( ap.fIaaAnalysis ),
	fOutput( ap.fOutput )
{ 

	AliInfo("----DEBUG AliJDiHadronIaaTask COPY ----");

}

//_____________________________________________________________________________
AliJDiHadronIaaTask& AliJDiHadronIaaTask::operator = (const AliJDiHadronIaaTask& ap)
{
	// assignment operator

	AliInfo("----DEBUG AliJDiHadronIaaTask operator= ----");
	this->~AliJDiHadronIaaTask();
	new(this) AliJDiHadronIaaTask(ap);
	return *this;
}

//______________________________________________________________________________
AliJDiHadronIaaTask::~AliJDiHadronIaaTask()
{
	// destructor

	delete fIaaAnalysis;

}

//________________________________________________________________________

void AliJDiHadronIaaTask::UserCreateOutputObjects()
{  
	//=== create the jcorran outputs objects
	if(fDebug > 1) printf("AliJDiHadronIaaTask::UserCreateOutPutData() \n");

	//=== Get AnalysisManager
	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();

	fFilterTask = (AliJCORRANTask*)(man->GetTask( fFilterTaskName ));

	OpenFile(1);
	fOutput = gDirectory;
	fOutput->cd();

	// Order should be kept for correct functionality
	fIaaAnalysis->SetRunHeader(fFilterTask->GetFilter()->GetAliJRunHeader()); // caution RunHeader must be taken from the filter since it is set in the task macro
	fIaaAnalysis->UserCreateOutputObjects();
	fIaaAnalysis->SetHeaderList(fFilterTask->GetFilter()->GetHeaderList());
	fIaaAnalysis->SetTrackList(fFilterTask->GetFilter()->GetTrackList());
	fIaaAnalysis->SetMCTrackList(fFilterTask->GetFilter()->GetMCTrackList());
//	fIaaAnalysis->GetCard()->WriteCard(fOutput);

	PostData(1, fOutput);

	cout << "Add(fAliRunHeader) in UserCreateObject() ======= " << endl;

}

//______________________________________________________________________________
void AliJDiHadronIaaTask::UserExec(Option_t* /*option*/)
{

	// Processing of one event
	if(fDebug > 5) cout << "------- AliJDiHadronIaaTask Exec-------"<<endl;

	if( fFilterTask->GetFilterEntry() != fEntry ) return;

	if( fFilterTask->GetFilter()->GetEventSuccess() ){
		fIaaAnalysis->UserExec();
		PostData(1, fOutput );
	}


	if(fDebug > 5) cout << "\t------- End UserExec "<<endl;
}

//______________________________________________________________________________
void AliJDiHadronIaaTask::Init()
{
	// Intialisation of parameters
	AliInfo("Doing initialization") ;
	fIaaAnalysis->Init();
}

//______________________________________________________________________________
void AliJDiHadronIaaTask::Terminate(Option_t *)
{

	cout<<"AliJDiHadronIaaTask Analysis DONE !!"<<endl;

}
