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
#include "AliJSoundsTask.h"

//______________________________________________________________________________
AliJSoundsTask::AliJSoundsTask() :   
	AliAnalysisTaskSE("JFlowBaseTask"),
	fJCatalystTask(NULL),
	fJCatalystTaskName("JCatalystTask"),
	fEfficiency(NULL),
        fFirstEvent(kTRUE),
	fIsMC(kTRUE),
	fhistos(NULL),
	fCBin(-1),
	fOutput(NULL)
{
}

//______________________________________________________________________________
AliJSoundsTask::AliJSoundsTask(const char *name, TString inputformat):
	AliAnalysisTaskSE(name), 
	fJCatalystTask(NULL),
	fJCatalystTaskName("JCatalystTask"),
	fEfficiency(NULL),
        fFirstEvent(kTRUE),
	fIsMC(kTRUE),
	fhistos(NULL),
	fCBin(-1),
	fOutput(NULL)
{
	// Constructor
	AliInfo("---- AliJSoundsTask Constructor ----");
	DefineOutput (1, TDirectory::Class());
}

//____________________________________________________________________________
AliJSoundsTask::AliJSoundsTask(const AliJSoundsTask& ap) :
	AliAnalysisTaskSE(ap.GetName()), 
	fJCatalystTask(ap.fJCatalystTask),
	fJCatalystTaskName(ap.fJCatalystTaskName),
	fEfficiency(ap.fEfficiency),
        fFirstEvent(ap.fFirstEvent),
	fIsMC(ap.fIsMC),
	fhistos(ap.fhistos),
	fCBin(ap.fCBin),
	fOutput( ap.fOutput )
{ 

	AliInfo("----DEBUG AliJSoundsTask COPY ----");

}

//_____________________________________________________________________________
AliJSoundsTask& AliJSoundsTask::operator = (const AliJSoundsTask& ap)
{
	// assignment operator

	AliInfo("----DEBUG AliJSoundsTask operator= ----");
	this->~AliJSoundsTask();
	new(this) AliJSoundsTask(ap);
	return *this;
}

//______________________________________________________________________________
AliJSoundsTask::~AliJSoundsTask()
{
	// destructor 
	delete fOutput;
	delete fhistos;

}

//________________________________________________________________________

void AliJSoundsTask::UserCreateOutputObjects()
{  
	//=== create the jcorran outputs objects
	if(fDebug > 1) printf("AliJSoundsTask::UserCreateOutPutData() \n");
	//=== Get AnalysisManager
	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();

	fJCatalystTask = (AliJCatalystTask*)(man->GetTask( fJCatalystTaskName ));

	OpenFile(1);
	fOutput = gDirectory;
	fOutput->cd();

	fhistos = new AliJFlowHistos();
	fhistos->CreateEventTrackHistos();

	fhistos->fHMG->Print();

	fEfficiency = new AliJEfficiency();
	fEfficiency->SetMode(1 ) ; // 1: priod should work for you
	fEfficiency->SetDataPath( "alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data" );

	fFirstEvent = kTRUE;

	PostData(1, fOutput);

}

//______________________________________________________________________________
void AliJSoundsTask::UserExec(Option_t* /*option*/) 
{

	// Processing of one event
	if(fDebug > 5) cout << "------- AliJSoundsTask Exec-------"<<endl;
	if(!((Entry()-1)%1000))  AliInfo(Form(" Processing event # %lld",  Entry())); 
	if( fJCatalystTask->GetJCatalystEntry() != fEntry ) return;
	fCBin = AliJFlowHistos::GetCentralityClass(fJCatalystTask->GetCentrality());
	if(fCBin == -1)
		return;

	if( fFirstEvent ) {
	    fEfficiency->SetRunNumber( fJCatalystTask->GetRunNumber() );
   	    fEfficiency->Load();
	    fFirstEvent = kFALSE;
	}

	TClonesArray *fInputList = (TClonesArray*)fJCatalystTask->GetInputList();
	BookHistos(fInputList);
}

//______________________________________________________________________________
void AliJSoundsTask::Init()
{
	// Intialisation of parameters
	AliInfo("Doing initialization") ; 

}

//______________________________________________________________________________
void AliJSoundsTask::Terminate(Option_t *)
{
	// Processing when the event loop is ended
	cout<<"AliJSoundsTask Analysis DONE !!"<<endl; 
}

//______________________________________________________________________________
void AliJSoundsTask::BookHistos(TClonesArray *inList) {

	int noTracks = inList->GetEntries();
	for(int itrack=0;itrack<noTracks; itrack++){
		AliJBaseTrack *trk = (AliJBaseTrack*)inList->At(itrack);
		double phi = trk->Phi();
		double eta = trk->Eta();
		double pt =  trk->Pt();
		Double_t effCorr = fEfficiency->GetCorrection( pt, fJCatalystTask->GetEffFilterBit(), fJCatalystTask->GetCentrality());
		fhistos->fh_eta[fCBin][0]->Fill(eta);
		fhistos->fh_phi[fCBin][0]->Fill(phi);
		fhistos->fh_pt[fCBin][0]->Fill(pt); // uncorrected
		fhistos->fh_pt[fCBin][1]->Fill(pt,1./effCorr); //corrected
	}
}
