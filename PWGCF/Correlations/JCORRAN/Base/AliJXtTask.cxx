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

#include "TChain.h"
#include "TFile.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODHandler.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliAODEvent.h" 
#include "AliJTrack.h"
#include "AliJMCTrack.h"
#include "AliJEventHeader.h"
#include "AliJHistManager.h"
#include "AliInputEventHandler.h"

#include "AliJXtTask.h"
#include "AliJXtAnalysis.h"
#include "AliJRunTable.h"
#include "AliJEfficiency.h"

//______________________________________________________________________________
AliJXtTask::AliJXtTask():
	fetaRange(0.8),
	fzvertexRange(10.0),
	fDebugLevel(0),
	fEvtNum(0),
	fTriggerBit(0),
	fTrackFilterBit(0),
	fEffFilterBit(0),
	frunNumber(-1),
	fEffMode(0),
	fSQRTS(7000.0),
	fFirstEvent(kTRUE),
	fIsMC(kFALSE),
	fIsPP(kTRUE),
	fRunTable(0x0),
	fEfficiency(0x0),
	fEfficiencyIsolated(0x0),
	fInputList(0x0),
	fXtAna(0x0),
	fOutput()
{
  // default constructor
}

//______________________________________________________________________________
AliJXtTask::AliJXtTask(const char *name, int CollisionCandidates, Bool_t IsMC):
	fetaRange(0.8),
	fzvertexRange(10.0),
	fDebugLevel(0),
	fEvtNum(0),
	fTriggerBit(0),
	fEffFilterBit(0),
	frunNumber(-1),
	fEffMode(0),
	fSQRTS(7000.0),
	fFirstEvent(kTRUE),
	fIsPP(kTRUE),
	fRunTable(0x0),
	fEfficiency(0x0),
	fEfficiencyIsolated(0x0),
	fInputList(0),
	AliAnalysisTaskSE(name), 
	fXtAna(0x0),
	fOutput()
{
  DefineOutput(1, TDirectory::Class());
  fTrackFilterBit = CollisionCandidates;
  fTaskName = name;
  fIsMC = IsMC;
}

//____________________________________________________________________________
AliJXtTask::AliJXtTask(const AliJXtTask& ap) :
	fetaRange(ap.fetaRange),
	fzvertexRange(ap.fzvertexRange),
	fDebugLevel(ap.fDebugLevel),
	fEvtNum(ap.fEvtNum),
	fTriggerBit(ap.fTriggerBit),
	fTrackFilterBit(ap.fTrackFilterBit),
	fEffFilterBit(ap.fEffFilterBit),
	frunNumber(ap.frunNumber),
	fEffMode(ap.fEffMode),
	fSQRTS(ap.fSQRTS),
	fFirstEvent(ap.fFirstEvent),
	fIsMC(ap.fIsMC),
	fIsPP(ap.fIsPP),
	fRunTable(ap.fRunTable),
	fEfficiency(ap.fEfficiency),
        fEfficiencyIsolated(ap.fEfficiencyIsolated),
	fInputList(ap.fInputList),
	AliAnalysisTaskSE(ap.GetName()),
	fXtAna(ap.fXtAna),
	fOutput(ap.fOutput)
{ 
  AliInfo("----DEBUG AliJXtTask COPY ----");
}

//_____________________________________________________________________________
AliJXtTask& AliJXtTask::operator = (const AliJXtTask& ap)
{
  // assignment operator
  AliInfo("----DEBUG AliJXtTask operator= ----");
  this->~AliJXtTask();
  new(this) AliJXtTask(ap);
  return *this;
}

//______________________________________________________________________________
AliJXtTask::~AliJXtTask()
{
	// destructor 
	delete fXtAna;
	delete fInputList;
	delete fOutput;
	delete fRunTable;
	delete fEfficiency;
	delete fEfficiencyIsolated;
}

//________________________________________________________________________
void AliJXtTask::UserCreateOutputObjects()
{ 
  fXtAna =  new AliJXtAnalysis( fTaskName );
  fXtAna->SetDebugLevel(fDebugLevel);
  fXtAna->SetEtaRange(fetaRange);
  fXtAna->SetZVertexRange(fzvertexRange);
  
  cout << "INFO: set etaRange = " << fetaRange << " and z-vertex range = " << fzvertexRange << endl;
  
  fInputList = new TClonesArray("AliJBaseTrack" , 2500);
  fInputList->SetOwner(kTRUE);

  OpenFile(1);
  fOutput = gDirectory;
  fOutput->cd();
  fXtAna->UserCreateOutputObjects(); 
  
  // currently, accept only MB events in the xT analysis
  fTriggerBit = AliVEvent::kMB;
  
  //cout << "======= DEBUG: UserCreateOutputObjects finished in the Task" << endl;
  
  PostData(1, fOutput);
}

//______________________________________________________________________________
void AliJXtTask::UserExec(Option_t* /*option*/) 
{
	//==== RunTable
	if( fFirstEvent ) {
	  // run table
	  frunNumber = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetEvent()->GetRunNumber();
	  if( frunNumber < 0 ){ cout << "ERROR: Task did not find run number!" << endl; exit(1); }
	  fRunTable = & AliJRunTable::GetSpecialInstance();
	  fRunTable->SetRunNumber( frunNumber );
	  fSQRTS = fRunTable->GetBeamEnergy(fRunTable->GetPeriod());
	  fIsPP = fRunTable->IsPP();
	  // Efficiency - no isolation
	  fEffMode = 1;
	  fEfficiency = new AliJEfficiency();
	  fEfficiency->SetMode( fEffMode ) ; // 0:NoEff 1:Period 2:RunNum 3:Auto
	  fEfficiency->SetRunNumber( frunNumber );
	  fEfficiency->SetDataPath( "alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data" );
	  fEfficiency->Load();
	  // Efficiency - isolated
	  fEfficiencyIsolated = new AliJEfficiency();
	  fEfficiencyIsolated->SetMode( fEffMode );
	  fEfficiencyIsolated->SetRunNumber( frunNumber );
	  fEfficiencyIsolated->SetDataPath( "alien:///alice/cern.ch/user/s/srasanen/efficiency/isol-rel10-cone04/" );
	  fEfficiencyIsolated->Load();
	  // update analysis
	  fXtAna->SetEfficiency( fEfficiency );
	  fXtAna->SetEfficiencyIsolated( fEfficiencyIsolated );
	  fXtAna->SetTrackFilterBit( fTrackFilterBit );
	  GetEfficiencyFilterBit( fTrackFilterBit );
	  fXtAna->SetEfficiencyFilterBit( fEffFilterBit );
	  fXtAna->SetSqrts( fSQRTS );
	  fFirstEvent = kFALSE;
	  //cout << "======= DEBUG: The first event configuration finished in the Task" << endl;
        }
	
	// Processing of one event
	if(!((Entry()-1)%100))  AliInfo(Form(" Processing event # %lld",  Entry()));
	
         // load current event and save track, event info
	AliAODEvent *currentEvent = dynamic_cast<AliAODEvent*>(InputEvent());
	fEvtNum++;
	
	// clear previous event 
	fInputList->Clear();
	
	float fCent = -999;
	double fvertex[3] = {-999,-999,-999};	

	// routine "IsGoodEvent()" checks is the event triggered and if it has a good vertex 
	if( IsGoodEvent( currentEvent )){
			//cout << "======= DEBUG: enter IsGoodEvent " << endl;
			ReadAODTracks( currentEvent, fInputList ) ; // read tracklist
			ReadVertexInfo( currentEvent, fvertex); // read vertex info
			// Set here centrality bin always to 0 for pp data
			fCent = fIsPP ? 0 : ReadAODCentrality( currentEvent, "V0M"  ); 
			Int_t Ntracks = fInputList->GetEntriesFast();

			// Analysis Part 
			fXtAna->Init();
			fXtAna->SetInputList( fInputList ); 
			fXtAna->SetEventCentrality( fCent );
			fXtAna->SetEventVertex( fvertex );
			fXtAna->UserExec(""); // analysis code 
	}
}

//______________________________________________________________________________
void AliJXtTask::ReadVertexInfo ( AliAODEvent *aod, double*  fvtx )
{
	fvtx[0] = aod->GetPrimaryVertex()->GetX();
	fvtx[1] = aod->GetPrimaryVertex()->GetY();
	fvtx[2] = aod->GetPrimaryVertex()->GetZ();
}
//______________________________________________________________________________
float AliJXtTask::ReadAODCentrality( AliAODEvent *aod, TString Trig )
{
	AliCentrality *cent = aod->GetCentrality();
	if(!cent){ Printf("No Cent event"); return -999; };
	return cent->GetCentralityPercentile(Trig.Data()); 
}
//______________________________________________________________________________
void AliJXtTask::Init()
{
	// initiationg all parameters and variables

	// Intialisation of parameters
	AliInfo("Doing initialization") ; 
//
}
//______________________________________________________________________________
void AliJXtTask::ReadAODTracks( AliAODEvent *aod , TClonesArray *TrackList)
{
	// 
	// NOTE: in task level no eta cut and very loose vertex cut!
	//
	//aod->Print();
	if( fIsMC == kTRUE ){  // how to get a flag to check  MC or not ! 

		TClonesArray *mcArray = (TClonesArray*) aod->FindListObject(AliAODMCParticle::StdBranchName());
		if(!mcArray){ Printf("Error not a proper MC event"); };  // check mc array
		
		Int_t nt = mcArray->GetEntriesFast();
		Int_t ntrack =0;
		for( int it=0; it< nt ; it++){
			AliAODMCParticle *track = (AliAODMCParticle*)mcArray->At(it);
			if(!track) { Error("ReadEventAODMC", "Could not receive particle %d",(int) it); continue; };
			if( track->IsPhysicalPrimary() ){
				Int_t pdg = track->GetPdgCode();
				Char_t ch = (Char_t) track->Charge();
				Int_t label = track->GetLabel();
				AliJBaseTrack *itrack = new ((*TrackList)[ntrack++])AliJBaseTrack;
				itrack->SetLabel( label );
				itrack->SetParticleType( pdg);
				itrack->SetPxPyPzE( track->Px(), track->Py(), track->Pz(), track->E() );
				itrack->SetCharge(ch) ;
			}
		}
	} // read mc track done.
	else if( fIsMC == kFALSE){
		Int_t nt = aod->GetNumberOfTracks();
		Int_t ntrack =0;
		for( int it=0; it<nt ; it++){
			AliAODTrack *track = dynamic_cast<AliAODTrack*>(aod->GetTrack(it));
			if(!track) { Error("ReadEventAOD", "Could not receive partice %d", (int) it); continue; };
			if(track->TestFilterBit( fTrackFilterBit )){ // hybrid cut
				AliJBaseTrack *itrack = new( (*TrackList)[ntrack++])AliJBaseTrack; 
				itrack->SetID( TrackList->GetEntriesFast() ) ;
				itrack->SetPxPyPzE( track->Px(), track->Py(), track->Pz(), track->E() );
				itrack->SetParticleType(kJHadron);
				itrack->SetCharge(track->Charge() );
				itrack->SetStatus(track->GetStatus() );
			}
		}
	} //read aod reco track done.
}


//______________________________________________________________________________
Bool_t AliJXtTask::IsGoodEvent( AliAODEvent *event ){
	
	Bool_t Event_status = kFALSE;
	
	// First, check the trigger
	Bool_t isTriggered = kFALSE;
	isTriggered = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & fTriggerBit);
	if( !isTriggered ) return Event_status;
	
	// in task-level require only that there is a vertex with at least 
	// one contributor. Require only a loose z-vertex cut
	AliVVertex *vtx = event->GetPrimaryVertex();
	if(vtx){
	  if( vtx->GetNContributors() > 0 ){
	    double zVert = vtx->GetZ();
	    if( zVert > -30.0 && zVert < 30.0) Event_status = kTRUE;
	  }
	}
	return Event_status;
}
//______________________________________________________________________________
void AliJXtTask::GetEfficiencyFilterBit(int inputTrackCut ){
	fEffFilterBit = 0; // default - corresponds to TPConly
	if( inputTrackCut == 128  ) fEffFilterBit = 0; // TPConly
	if( inputTrackCut == 1024 ) fEffFilterBit = 1; // RAA
	if( inputTrackCut == 32   ) fEffFilterBit = 2; // Global tight DCA
	if( inputTrackCut == 16   ) fEffFilterBit = 3; // Global DCA
	if( inputTrackCut == 96   ) fEffFilterBit = 4; // Global SDD
	if( inputTrackCut == 768  ) fEffFilterBit = 5; // Hybrid
	if( inputTrackCut == 272  ) fEffFilterBit = 6; // Hybrid at AOD86
	return;
}
//______________________________________________________________________________
void AliJXtTask::Terminate(Option_t *)
{
	// Processing when the event loop is ended
	cout<<"AliJXtTask Analysis DONE !!"<<endl; 
}

