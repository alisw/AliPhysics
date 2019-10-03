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

#include "AliJDiHadronCorrTask.h" 
#include "AliJEventHeader.h"
#include "AliJRunHeader.h"
#include "AliJCORRANTask.h"
#include "AliJCard.h"

//______________________________________________________________________________
AliJDiHadronCorrTask::AliJDiHadronCorrTask() :   
    AliAnalysisTaskSE("AliJDiHadronCorrTaskTask"),
	fFilterTask(NULL),
	fFilterTaskName(""),
    fJCORRAN(0x0),
	fOutput(NULL)
{
  DefineOutput (1, TDirectory::Class());
}

//______________________________________________________________________________
AliJDiHadronCorrTask::AliJDiHadronCorrTask(const char *name, TString inputformat):
    AliAnalysisTaskSE(name), 
	fFilterTask(NULL),
	fFilterTaskName(""),
    fJCORRAN(0x0),
	fOutput(NULL)
{
  // Constructor
  AliInfo("---- AliJDiHadronCorrTask Constructor ----");

  JUNUSED(inputformat);
  DefineOutput (1, TDirectory::Class());
}

//____________________________________________________________________________
AliJDiHadronCorrTask::AliJDiHadronCorrTask(const AliJDiHadronCorrTask& ap) :
    AliAnalysisTaskSE(ap.GetName()), 
	fFilterTask(ap.fFilterTask),
	fFilterTaskName(ap.fFilterTaskName),
    fJCORRAN( ap.fJCORRAN ),
	fOutput( ap.fOutput )
{ 

  AliInfo("----DEBUG AliJDiHadronCorrTask COPY ----");

}

//_____________________________________________________________________________
AliJDiHadronCorrTask& AliJDiHadronCorrTask::operator = (const AliJDiHadronCorrTask& ap)
{
  // assignment operator

  AliInfo("----DEBUG AliJDiHadronCorrTask operator= ----");
  this->~AliJDiHadronCorrTask();
  new(this) AliJDiHadronCorrTask(ap);
  return *this;
}

//______________________________________________________________________________
AliJDiHadronCorrTask::~AliJDiHadronCorrTask()
{
  // destructor 

   delete fJCORRAN;

}

//________________________________________________________________________

void AliJDiHadronCorrTask::UserCreateOutputObjects()
{  
  //=== create the jcorran outputs objects
  if(fDebug > 1) printf("AliJDiHadronCorrTask::UserCreateOutPutData() \n");
  
  //=== Get AnalysisManager
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if(!man->GetOutputEventHandler()) {
    Fatal("UserCreateOutputObjects", "This task needs an AOD handler");
    return;
  }

   fFilterTask = (AliJCORRANTask*)(man->GetTask( fFilterTaskName ));

   OpenFile(1);
   fOutput = gDirectory;//->mkdir("JDiHadronCorr");
   fOutput->cd();

	// Order should be kept
	// TODO : better implementation
   //bool orignalTH1AdddirectoryStatus=TH1::AddDirectoryStatus();
   // TODO Why? 
   //TH1::AddDirectory(kTRUE);
   //if( !orignalTH1AdddirectoryStatus ) cout<<"DEBUG :: TH1::AddDirectory is turned on"<<endl;
   fJCORRAN->SetRunHeader( fFilterTask->GetJRunHeader() );//TODO
   fJCORRAN->UserCreateOutputObjects();
   fJCORRAN->SetHeaderList(fFilterTask->GetFilter()->GetHeaderList() );
   fJCORRAN->SetTrackList( fFilterTask->GetFilter()->GetTrackList()     );
   fJCORRAN->SetMCTrackList( fFilterTask->GetFilter()->GetMCTrackList()    );
   fJCORRAN->GetCard()->WriteCard(fOutput);

   PostData( 1, fOutput );
   //TH1::AddDirectory( orignalTH1AdddirectoryStatus );
   //cout<<"DEBYG :: TH1::AddDirectory get orignal Value = "<<( orignalTH1AdddirectoryStatus?"True":"False" )<<endl;

   cout << "Add(fAliRunHeader) in UserCreateObject() ======= " << endl;

}

//______________________________________________________________________________
void AliJDiHadronCorrTask::UserExec(Option_t* /*option*/) 
{

	// Processing of one event
	if(fDebug > 5) cout << "------- AliJDiHadronCorrTask Exec-------"<<endl;
	//if(!((Entry()-1)%100))  cout << Form(" Processing event # %lld",  Entry()) << endl; 
	if( fFilterTask->GetFilterEntry() != fEntry ) return;

	//cout << "AliJDiHadronCorrTask::UserExec fFilterTask->GetFilter()->GetEventSuccess() = " << fFilterTask->GetFilter()->GetEventSuccess() << endl;
	if( fFilterTask->GetFilter()->GetEventSuccess() ){
		fJCORRAN->UserExec();
		PostData(1, fOutput );
	}


	if(fDebug > 5) cout << "\t------- End UserExec "<<endl;
}

//______________________________________________________________________________
void AliJDiHadronCorrTask::Init()
{
	// Intialisation of parameters
	AliInfo("Doing initialization") ; 
	fJCORRAN->Init();
}

//______________________________________________________________________________
void AliJDiHadronCorrTask::Terminate(Option_t *)
{
	/*
	   fJCORRAN->Terminate();
	   OpenFile(1);
	   fOutput->cd();

	   cout<<"# Write Data "<<endl;
	//write Trigg ID and Assoc ID
	fOutput->Write(0,TObject::kOverwrite);
	PostData(1,fOutput);
	*/
	cout<<"AliJDiHadronCorrTask Analysis DONE !!"<<endl; 

}
