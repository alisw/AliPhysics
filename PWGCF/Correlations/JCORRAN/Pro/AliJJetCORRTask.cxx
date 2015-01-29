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
// Analysis task for jet correlation Analysis
// author: B.S Chang D.J. Kim
// ALICE Group University of Jyvaskyla 
// Finland 
// Fill the analysis containers for ESD or AOD
// Adapted for AliAnalysisTaskSE and AOD objects  
//////////////////////////////////////////////////////////////////////////////

#include "AliAODEvent.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliJJetCORRTask.h" 
#include "AliJEventHeader.h"
#include "AliJCard.h"
#include "AliJRunTable.h"
#include "AliJEfficiency.h" 
#include "AliJHistos.h"

//______________________________________________________________________________
AliJJetCORRTask::AliJJetCORRTask() :   
    AliAnalysisTaskSE("AliJJetCORRTaskTask"),
	fJetTask(NULL),
	fJetTaskName(""),
    fJJetCORRAnalysis(0x0),
	fOutput(NULL),
    fCard(NULL),
    fAnaUtils(NULL),
    fRunTable(NULL),
    fFirstEvent(kTRUE),
    cBin(-1),
    zBin(-1),
	zVert(-999),
	fDebugMode(0),
	fTargetJetIndex(0),
	fevt(-1)
{
	DefineOutput (1, TDirectory::Class());
}

//______________________________________________________________________________
AliJJetCORRTask::AliJJetCORRTask(const char *name, TString inputformat):
	AliAnalysisTaskSE(name), 
	fJetTask(NULL),
	fJetTaskName(""),
	fJJetCORRAnalysis(0x0),
	fOutput(NULL),
	fCard(NULL),
	fAnaUtils(NULL),
	fRunTable(NULL),
	fFirstEvent(kTRUE),
	cBin(-1),
	zBin(-1),
	zVert(-999),
	fDebugMode(0),
	fTargetJetIndex(0),
	fevt(-1)
{
	// Constructor
	AliInfo("---- AliJJetCORRTask Constructor ----");

	JUNUSED(inputformat);
	DefineOutput (1, TDirectory::Class());
}

//____________________________________________________________________________
AliJJetCORRTask::AliJJetCORRTask(const AliJJetCORRTask& ap) :
	AliAnalysisTaskSE(ap.GetName()), 
	fJetTask(ap.fJetTask),
	fJetTaskName(ap.fJetTaskName),
	fJJetCORRAnalysis( ap.fJJetCORRAnalysis ),
	fOutput( ap.fOutput ),
	fCard( ap.fCard ),
	fAnaUtils(ap.fAnaUtils),
	fRunTable(ap.fRunTable),
	fFirstEvent(kTRUE),
	cBin(-1),
	zBin(-1),
	zVert(-999),
	fDebugMode(0),
	fTargetJetIndex(0),
	fevt(-1)
{ 

	AliInfo("----DEBUG AliJJetCORRTask COPY ----");

}

//_____________________________________________________________________________
AliJJetCORRTask& AliJJetCORRTask::operator = (const AliJJetCORRTask& ap)
{
	// assignment operator

	AliInfo("----DEBUG AliJJetCORRTask operator= ----");
	this->~AliJJetCORRTask();
	new(this) AliJJetCORRTask(ap);
	return *this;
}

//______________________________________________________________________________
AliJJetCORRTask::~AliJJetCORRTask()
{
	// destructor 

	delete fJJetCORRAnalysis;
	delete fAnaUtils;

}

//________________________________________________________________________

void AliJJetCORRTask::UserCreateOutputObjects()
{  
	//=== create the jcorran outputs objects
	if(fDebugMode > 1) printf("AliJJetCORRTask::UserCreateOutPutData() \n");

	fAnaUtils = new AliAnalysisUtils();
	fAnaUtils->SetUseOutOfBunchPileUp( kTRUE );

	//=== Get AnalysisManager
	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
	fJetTask = (AliJJetTask*)(man->GetTask( fJetTaskName));

	OpenFile(1);
	fOutput = gDirectory;//->mkdir("JDiHadronCorr");
	fOutput->cd();

	fJJetCORRAnalysis = new AliJJetCORRAnalysis(fCard);
	fJJetCORRAnalysis->UserCreateOutputObjects();
	TNamed jetFinderName("JetFinder", fJetTask->GetJetFinderString()[fTargetJetIndex].Data() );
	jetFinderName.Write();
	fCard->WriteCard( gDirectory );
	

	PostData( 1, fOutput );


	cout << fJetTaskName <<"\t NumberOfJetCollections = "<< fJetTask->GetNumberOfJetCollections() << endl;

	// Add jets into the fJJetCORRAnalysis-
	for( int ij=0;ij< fJetTask->GetNumberOfJetCollections();ij++ ){
		fJJetCORRAnalysis->AddJets( fJetTask->GetAliJJetList( ij ) );
	}

	fFirstEvent = kTRUE;
	fevt = -1;

	cout << "Add(fAliRunHeader) in UserCreateObject() ======= " << endl;
}

//______________________________________________________________________________
void AliJJetCORRTask::UserExec(Option_t* /*option*/) 
{

	// Processing of one event
	if(fDebugMode > 5) cout << "------- AliJJetCORRTask Exec-------"<<endl;
	if( fJetTask->GetTaskEntry() != fEntry ) return; // Make sure if we loop over the same event

	fevt++;
	if(fevt % 1000 == 0) cout << "AliJJetCORRTask:: Numer of event scanned = "<< fevt << endl;
	
	fJJetCORRAnalysis->ClearBeforeEvent();

	// Main loop
	// Called for each event
	AliVEvent *event = InputEvent();
	if(!event) return;

	//---------------------------------------------------------------
	// check if the event was triggered or not and vertex

	AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(event);
	if(!aodEvent) return;

	if( fFirstEvent ) {
		fRunTable = & AliJRunTable::GetSpecialInstance();
		fRunTable->SetRunNumber( aodEvent->GetRunNumber() );
		fJJetCORRAnalysis->GetAliJEfficiency()->SetRunNumber( aodEvent->GetRunNumber() );
		fJJetCORRAnalysis->GetAliJEfficiency()->Load();
		fFirstEvent = kFALSE;
	}

	if(!IsGoodEvent( event )) return; // zBin is set there
	if(fDebugMode>5) cout << "zvtx = " << zVert << endl;

	// centrality
	float fcent = -999;
	if(fRunTable->IsHeavyIon() || fRunTable->IsPA()){
		AliCentrality *cent = event->GetCentrality();
		if( ! cent ) return;
		fcent = cent->GetCentralityPercentile("V0M");
	} else {
		fcent = -1;
	}
	cBin = fCard->GetBin(kCentrType, fcent);;

	if(cBin<0) return;

	fJJetCORRAnalysis->SetCentralityBin( cBin );
	fJJetCORRAnalysis->SetCentrality( fcent );
	fJJetCORRAnalysis->SetZVertexBin( zBin );
	fJJetCORRAnalysis->SetZVertex( zVert );
	fJJetCORRAnalysis->SetTargetJetIndex( fTargetJetIndex );

	fJJetCORRAnalysis->UserExec();

	PostData(1, fOutput );

	if(fDebugMode > 5) cout << "\t------- End UserExec "<<endl;
}

//______________________________________________________________________________
void AliJJetCORRTask::Init()
{
	// Intialisation of parameters
	AliInfo("Doing initialization") ; 
	//fJJetCORRAnalysis->Init();
}

//______________________________________________________________________________
void AliJJetCORRTask::Terminate(Option_t *)
{
	cout<<"AliJJetCORRTask Analysis DONE !!"<<endl; 

}

//________________________________________________________________________
bool AliJJetCORRTask::IsGoodEvent(AliVEvent *event) {

	if(fRunTable->IsPP() && fAnaUtils->IsPileUpEvent(event)) {
		return kFALSE;
	} else {
		Bool_t triggeredEventMB = kFALSE; //init

		fJJetCORRAnalysis->GetAliJHistos()->fhEvents->Fill( 0 );

		Bool_t triggerkMB = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & ( AliVEvent::kMB );

		if( triggerkMB ){
			triggeredEventMB = kTRUE;  //event triggered as minimum bias
			fJJetCORRAnalysis->GetAliJHistos()->fhEvents->Fill( 1 );
		}
		//--------------------------------------------------------------
		// check reconstructed vertex
		int ncontributors = 0;
		Bool_t goodRecVertex = kFALSE;
		const AliVVertex *vtx = event->GetPrimaryVertex();
		if(vtx){
			ncontributors = vtx->GetNContributors();
			if(ncontributors > 0){
				zVert = vtx->GetZ();
				fJJetCORRAnalysis->GetAliJHistos()->fhEvents->Fill( 2 );
				if(fCard->VertInZRange(zVert)) {
					goodRecVertex = kTRUE;
					fJJetCORRAnalysis->GetAliJHistos()->fhEvents->Fill( 3 );
					zBin  = fCard->GetBin(kZVertType, zVert);
				}
			}
		}
		return goodRecVertex;
	}
	//---------------------------------
}
