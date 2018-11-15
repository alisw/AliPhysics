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
// Analysis task for high pt particle correlations w.r.t EP 
// author: D. J.~Kim
// ALICE Group University of Jyvaskyla 
// Finland 
//////////////////////////////////////////////////////////////////////////////

#include <AliAnalysisTaskSE.h>
#include <AliAODHandler.h>
#include <AliAnalysisManager.h>
#include <AliLog.h>

#include <AliJCard.h>
#include "AliJCIaaEPTask.h"
//#include "AliJFlowBaseTask.h"

//______________________________________________________________________________
AliJCIaaEPTask::AliJCIaaEPTask() :   
	AliAnalysisTaskSE("AliJCIaaEPTask"),
	fJFlowBaseTask(NULL),
	fJFlowBaseTaskName("JFlowBaseTask"),
	fEPDetID(0),
	fIaaAna(0x0),
	fOutput(NULL)
{
	DefineOutput (1, TDirectory::Class());
}

//______________________________________________________________________________
AliJCIaaEPTask::AliJCIaaEPTask(const char *name, TString inputformat):
	AliAnalysisTaskSE(name),
	fJFlowBaseTask(NULL),
	fJFlowBaseTaskName("JFlowBaseTask"),
	fEPDetID(0),
	fIaaAna(0x0),
	fOutput(NULL)
{
	// Constructor
	AliInfo("---- AliJCIaaEPTask Constructor ----");

	JUNUSED(inputformat);
	DefineOutput (1, TDirectory::Class());
}

//____________________________________________________________________________
AliJCIaaEPTask::AliJCIaaEPTask(const AliJCIaaEPTask& ap) :
	AliAnalysisTaskSE(ap.GetName()),
	fJFlowBaseTask(ap.fJFlowBaseTask),
	fJFlowBaseTaskName(ap.fJFlowBaseTaskName),
	fIaaAna( ap.fIaaAna ),
	fOutput( ap.fOutput )
{ 

	AliInfo("----DEBUG AliJCIaaEPTask COPY ----");

}

//_____________________________________________________________________________
AliJCIaaEPTask& AliJCIaaEPTask::operator = (const AliJCIaaEPTask& ap)
{
	// assignment operator

	AliInfo("----DEBUG AliJCIaaEPTask operator= ----");
	this->~AliJCIaaEPTask();
	new(this) AliJCIaaEPTask(ap);
	return *this;
}

//______________________________________________________________________________
AliJCIaaEPTask::~AliJCIaaEPTask()
{
	// destructor 

	delete fIaaAna;

}

//________________________________________________________________________

void AliJCIaaEPTask::UserCreateOutputObjects()
{  
	//=== create the jcorran outputs objects
	if(fDebug > 1) printf("AliJCIaaEPTask::UserCreateOutPutData() \n");

	//=== Get AnalysisManager
	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();

	fJFlowBaseTask = (AliJFlowBaseTask*)(man->GetTask( fJFlowBaseTaskName ));

	OpenFile(1);
	fOutput = gDirectory;
	fOutput->cd();

	// Order should be kept for correct functionality
	fIaaAna->UserCreateOutputObjects();
	fIaaAna->GetCard()->WriteCard(fOutput);

	PostData(1, fOutput);

}

//______________________________________________________________________________
void AliJCIaaEPTask::UserExec(Option_t* /*option*/) 
{

	// Processing of one event
	if(fDebug > 5) std::cout << "------- AliJCIaaEPTask Exec-------"<<endl;
	if( fJFlowBaseTask->GetJCatalystTask()->GetJCatalystEntry() != fEntry ) return;
	if(fDebug > 5 ) std::cout << fJFlowBaseTask->GetJCatalystTask()->GetJCatalystTaskName() << endl;
	if(fJFlowBaseTask->IsMC()) {
		fIaaAna->SetTrackList(fJFlowBaseTask->GetJCatalystTask()->GetInputListALICE());
	} else {
		fIaaAna->SetTrackList(fJFlowBaseTask->GetJCatalystTask()->GetInputList());
	}
	// Event selection is done in fJFlowBaseTask
	// need to put cent/vtx/run#
	int iH = 0; // 0:2nd 1:3rd
	double EP2 = -999;
	if(fJFlowBaseTask->IsMC()) {
		TComplex Qvector;
		Qvector = fJFlowBaseTask->GetQvectorsEP(fEPDetID, iH);
		EP2 = Qvector.Theta()/double(iH+2);
	} else {
		EP2 = fJFlowBaseTask->GetEventPlaneALICE(fEPDetID,iH);
	}
	if(fDebug > 5) printf("entry, EP2 = %d, %.2f\n", fEntry,EP2);
	if(fDebug > 5) printf("runN,cent,zvtx = %d, %.1f, %.1f\n",fJFlowBaseTask->GetJCatalystTask()->GetRunNumber(),fJFlowBaseTask->GetJCatalystTask()->GetCentrality(),fJFlowBaseTask->GetJCatalystTask()->GetZVertex());
	fIaaAna->SetRunNumber(fJFlowBaseTask->GetJCatalystTask()->GetRunNumber());
	fIaaAna->SetCentrality(fJFlowBaseTask->GetJCatalystTask()->GetCentrality());
	fIaaAna->SetZVertex(fJFlowBaseTask->GetJCatalystTask()->GetZVertex());
	fIaaAna->SetEventPlane(EP2);
	fIaaAna->SetEPmin(fEPmin);
	fIaaAna->SetEPmax(fEPmax);
	fIaaAna->UserExec();
	PostData(1, fOutput );


	if(fDebug > 5) std::cout << "\t------- End UserExec "<<endl;
}

//______________________________________________________________________________
void AliJCIaaEPTask::Init()
{
	// Intialisation of parameters
	AliInfo("Doing initialization") ; 
	fIaaAna->Init();
}

//______________________________________________________________________________
void AliJCIaaEPTask::Terminate(Option_t *)
{

	std::cout<<"AliJCIaaEPTask Analysis DONE !!"<<endl;

}
