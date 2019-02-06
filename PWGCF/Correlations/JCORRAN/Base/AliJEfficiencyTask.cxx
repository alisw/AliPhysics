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

//______________________________________________________________________________
// Analysis task for high pt particle correlations 
// author: R.Diaz, J. Rak,  D.J. Kim
// ALICE Group University of Jyvaskyla 
// Finland 
// Fill the analysis containers for ESD or AOD
// Adapted for AliAnalysisTaskSE and AOD objects  
//////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskSE.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"

#include "AliJEfficiencyTask.h" 
#include "AliJEventHeader.h"
#include "AliJRunHeader.h"
#include "AliJCORRANTask.h"

//______________________________________________________________________________
AliJEfficiencyTask::AliJEfficiencyTask() :   
    AliAnalysisTaskSE("JEfficiencyTask"),
	fFilterTask(NULL),
	fFilterTaskName(""),
    fEfficiencyScanner(0x0),
    fEffHistDir(0x0)
{
  DefineOutput (1, TDirectory::Class());
}

//______________________________________________________________________________
AliJEfficiencyTask::AliJEfficiencyTask(const char *name, TString inputformat):
    AliAnalysisTaskSE(name), 
	fFilterTask(NULL),
	fFilterTaskName(""),
    fEfficiencyScanner(0x0),
    fEffHistDir(0x0)
{
  // Constructor
  AliInfo("---- AliJEfficiencyTask Constructor ----");

  JUNUSED(inputformat);
  DefineOutput (1, TDirectory::Class());
  fEfficiencyScanner = new AliJEfficiencyScanner( Form("%sEffScanner",name ));
}

//____________________________________________________________________________
AliJEfficiencyTask::AliJEfficiencyTask(const AliJEfficiencyTask& ap) :
    AliAnalysisTaskSE(ap.GetName()), 
	fFilterTask(ap.fFilterTask),
	fFilterTaskName(ap.fFilterTaskName),
    fEfficiencyScanner( ap.fEfficiencyScanner ),
    fEffHistDir( ap.fEffHistDir )
{ 

  AliInfo("----DEBUG AliJEfficiencyTask COPY ----");

}

//_____________________________________________________________________________
AliJEfficiencyTask& AliJEfficiencyTask::operator = (const AliJEfficiencyTask& ap)
{
  // assignment operator

  AliInfo("----DEBUG AliJEfficiencyTask operator= ----");
  this->~AliJEfficiencyTask();
  new(this) AliJEfficiencyTask(ap);
  return *this;
}

//______________________________________________________________________________
AliJEfficiencyTask::~AliJEfficiencyTask()
{
  // destructor 

   delete fEfficiencyScanner;

}

//________________________________________________________________________

void AliJEfficiencyTask::UserCreateOutputObjects()
{  
  //=== create the jcorran outputs objects
  if(fDebug > 1) printf("AliJEfficiencyTask::UserCreateOutPutData() \n");
  
   AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
   fFilterTask = (AliJCORRANTask*)(man->GetTask( fFilterTaskName ));

   OpenFile(1);

   fEfficiencyScanner->SetJRunHeader( fFilterTask->GetFilter()->GetAliJRunHeader()); // caution RunHeader must be taken from the filter since it is set in the task macro
   fEfficiencyScanner->SetJTrackList( fFilterTask->GetFilter()->GetTrackList()     );
   fEfficiencyScanner->SetJMCTrackList( fFilterTask->GetFilter()->GetMCTrackList()    );
   fEffHistDir = gDirectory;//->mkdir("EffHist"); // no need
   fEffHistDir->cd();
   fEfficiencyScanner->UserCreateOutputObjects();
   PostData( 1, fEffHistDir );

   cout << "Add(fAliRunHeader) in UserCreateObject() ======= " << endl;

}

//______________________________________________________________________________
void AliJEfficiencyTask::UserExec(Option_t* /*option*/) 
{

	// Processing of one event
	if(fDebug > 5) cout << "------- AliJEfficiencyTask Exec-------"<<endl;
	if(!((Entry()-1)%100))  AliInfo(Form(" Processing event # %lld",  Entry())); 
	if( fFilterTask->GetFilterEntry() != fEntry ) return;

	if( fFilterTask->GetFilter()->GetEventSuccess() ){
		fEfficiencyScanner->SetJEventHeader( (AliJEventHeader*) fFilterTask->GetFilter()->GetHeaderList()->At(0) );
		fEfficiencyScanner->UserExec("");
		PostData(1, fEffHistDir );
	}


	if(fDebug > 5) cout << "\t------- End UserExec "<<endl;
}

//______________________________________________________________________________
void AliJEfficiencyTask::Init()
{
	// Intialisation of parameters
	AliInfo("Doing initialization") ; 
	fEfficiencyScanner->Init();
}

//______________________________________________________________________________
void AliJEfficiencyTask::Terminate(Option_t *)
{

	//fEfficiencyScanner->Terminate();
	//OpenFile(1);
	//fEffHistDir->Write();
	cout<<"JEfficiency Analysis DONE !!"<<endl; 

}
