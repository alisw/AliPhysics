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
// author: Jussi Viinikainen
// ALICE Group University of Jyvaskyla 
// Finland 
// Fill the analysis containers for ESD or AOD
// Adapted for AliAnalysisTaskSE and AOD objects  
//////////////////////////////////////////////////////////////////////////////

#include <AliAnalysisTaskSE.h>
#include <AliAODHandler.h>
#include <AliAnalysisManager.h>
#include <AliLog.h>

#include "AliJDiHadronJtTask.h"
#include "../AliJEventHeader.h"
#include "../AliJRunHeader.h"
#include "../AliJCORRANTask.h"
#include "../AliJCard.h"

//______________________________________________________________________________
AliJDiHadronJtTask::AliJDiHadronJtTask() :   
  AliAnalysisTaskSE("AliJDiHadronJtTaskTask"),
  fFilterTask(NULL),
  fFilterTaskName(""),
  fJtAnalysis(0x0),
  fOutput(NULL)
{
  DefineOutput (1, TDirectory::Class());
}

//______________________________________________________________________________
AliJDiHadronJtTask::AliJDiHadronJtTask(const char *name, TString inputformat):
  AliAnalysisTaskSE(name),
  fFilterTask(NULL),
  fFilterTaskName(""),
  fJtAnalysis(0x0),
  fOutput(NULL)
{
  // Constructor
  AliInfo("---- AliJDiHadronJtTask Constructor ----");

  JUNUSED(inputformat);
  DefineOutput (1, TDirectory::Class());
}

//____________________________________________________________________________
AliJDiHadronJtTask::AliJDiHadronJtTask(const AliJDiHadronJtTask& ap) :
  AliAnalysisTaskSE(ap.GetName()),
  fFilterTask(ap.fFilterTask),
  fFilterTaskName(ap.fFilterTaskName),
  fJtAnalysis( ap.fJtAnalysis ),
  fOutput( ap.fOutput )
{ 

  AliInfo("----DEBUG AliJDiHadronJtTask COPY ----");

}

//_____________________________________________________________________________
AliJDiHadronJtTask& AliJDiHadronJtTask::operator = (const AliJDiHadronJtTask& ap)
{
  // assignment operator

  AliInfo("----DEBUG AliJDiHadronJtTask operator= ----");
  this->~AliJDiHadronJtTask();
  new(this) AliJDiHadronJtTask(ap);
  return *this;
}

//______________________________________________________________________________
AliJDiHadronJtTask::~AliJDiHadronJtTask()
{
  // destructor 

   delete fJtAnalysis;

}

//________________________________________________________________________

void AliJDiHadronJtTask::UserCreateOutputObjects()
{  
  //=== create the jcorran outputs objects
  if(fDebug > 1) printf("AliJDiHadronJtTask::UserCreateOutPutData() \n");
  
  //=== Get AnalysisManager
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  // no need!!!!
  //if(!man->GetOutputEventHandler()) {
  //  Fatal("UserCreateOutputObjects", "This task needs an AOD handler");
  //  return;
  //}

   fFilterTask = (AliJCORRANTask*)(man->GetTask( fFilterTaskName ));

   OpenFile(1);
   fOutput = gDirectory;
   fOutput->cd();

	// Order should be kept for correct functionality
   fJtAnalysis->SetRunHeader(fFilterTask->GetFilter()->GetAliJRunHeader()); // caution RunHeader must be taken from the filter since it is set in the task macro
   fJtAnalysis->UserCreateOutputObjects();
   fJtAnalysis->SetHeaderList(fFilterTask->GetFilter()->GetHeaderList());
   fJtAnalysis->SetTrackList(fFilterTask->GetFilter()->GetTrackList());
   fJtAnalysis->SetMCTrackList(fFilterTask->GetFilter()->GetMCTrackList());
   fJtAnalysis->GetCard()->WriteCard(fOutput);

   PostData(1, fOutput);

   cout << "Add(fAliRunHeader) in UserCreateObject() ======= " << endl;

}

//______________________________________________________________________________
void AliJDiHadronJtTask::UserExec(Option_t* /*option*/) 
{

	// Processing of one event
	if(fDebug > 5) cout << "------- AliJDiHadronJtTask Exec-------"<<endl;

	if( fFilterTask->GetFilterEntry() != fEntry ) return;

	if( fFilterTask->GetFilter()->GetEventSuccess() ){
		fJtAnalysis->UserExec();
		PostData(1, fOutput );
	}


	if(fDebug > 5) cout << "\t------- End UserExec "<<endl;
}

//______________________________________________________________________________
void AliJDiHadronJtTask::Init()
{
	// Intialisation of parameters
	AliInfo("Doing initialization") ; 
	fJtAnalysis->Init();
}

//______________________________________________________________________________
void AliJDiHadronJtTask::Terminate(Option_t *)
{

	cout<<"AliJDiHadronJtTask Analysis DONE !!"<<endl;

}
