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
// Analysis task for di-jet Analysis
// author: B.K Kim,  D.J. Kim
// ALICE Group University of Jyvaskyla 
// Finland 
// Fill the analysis containers for ESD or AOD
// Adapted for AliAnalysisTaskSE and AOD objects  
//////////////////////////////////////////////////////////////////////////////

#include "AliAODEvent.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODHandler.h"
#include "AliParticleContainer.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"
#include "AliJMCTrack.h"
#include "AliJJetJtTask.h" 
#include "AliJCard.h"
#include "AliJetContainer.h"
#include "AliJEfficiency.h"

//______________________________________________________________________________
AliJJetJtTask::AliJJetJtTask() :   
    AliAnalysisTaskSE("AliJJetJtTaskTask"),
	fJetTask(NULL),
	fMCJetTask(NULL),
	fJetTaskName(""),
	fMCJetTaskName(""),
    fJJetJtAnalysis(0x0),
	fOutput(NULL),
    fCard(NULL),
    fFirstEvent(kTRUE),
    cBin(-1),
    zBin(-1),
    fDoMC(0),
    NRandom(1),
    moveJet(0),
    zVert(-999),
    fAnaUtils(NULL),
    fRunTable(NULL),
  fJMCTracks(NULL),
    fEventHist(NULL)
    
{
  DefineOutput (1, TDirectory::Class());
}

//______________________________________________________________________________
AliJJetJtTask::AliJJetJtTask(const char *name, TString inputformat):
    AliAnalysisTaskSE(name), 
	fJetTask(NULL),
	fMCJetTask(NULL),
	fJetTaskName(""),
	fMCJetTaskName(""),
    fJJetJtAnalysis(0x0),
    fOutput(NULL),
    fCard(NULL),
    fFirstEvent(kTRUE),
    cBin(-1),
    zBin(-1),
    fDoMC(0),
    NRandom(1),
    moveJet(0),
    zVert(-999),
    fAnaUtils(NULL),
    fRunTable(NULL),
  fJMCTracks(NULL),
    fEventHist(NULL)
{
  // Constructor
  AliInfo("---- AliJJetJtTask Constructor ----");

  JUNUSED(inputformat);
  DefineOutput (1, TDirectory::Class());
}

//____________________________________________________________________________
AliJJetJtTask::AliJJetJtTask(const AliJJetJtTask& ap) :
    AliAnalysisTaskSE(ap.GetName()), 
	fJetTask(ap.fJetTask),
	fJetTaskName(ap.fJetTaskName),
	fMCJetTask(ap.fMCJetTask),
	fMCJetTaskName(ap.fMCJetTaskName),
    fJJetJtAnalysis( ap.fJJetJtAnalysis ),
	fOutput( ap.fOutput ),
    fCard(ap.fCard),
    fFirstEvent(ap.fFirstEvent),
    cBin(-1),
    zBin(-1),
    NRandom(1),
    moveJet(0),
    zVert(-999),
  fJMCTracks(ap.fJMCTracks),
    fAnaUtils(ap.fAnaUtils),
    fRunTable(ap.fRunTable),
    fEventHist(ap.fEventHist)
{ 

  AliInfo("----DEBUG AliJJetJtTask COPY ----");

}

//_____________________________________________________________________________
AliJJetJtTask& AliJJetJtTask::operator = (const AliJJetJtTask& ap)
{
  // assignment operator

  AliInfo("----DEBUG AliJJetJtTask operator= ----");
  this->~AliJJetJtTask();
  new(this) AliJJetJtTask(ap);
  return *this;
}

//______________________________________________________________________________
AliJJetJtTask::~AliJJetJtTask()
{
  // destructor 

    delete fJJetJtAnalysis;

}

//________________________________________________________________________

void AliJJetJtTask::UserCreateOutputObjects()
{  
  //=== create the jcorran outputs objects
  if(fDebug > 1) printf("AliJJetJtTask::UserCreateOutPutData() \n");
  
  //=== Get AnalysisManager
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  //OpenFile(1);
  //if(!man->GetOutputEventHandler()) {
  //  Fatal("UserCreateOutputObjects", "This task needs an AOD handler");
  //  return;
  //}


   OpenFile(1);
   fOutput = gDirectory;//->mkdir("JDiHadronCorr");
   fOutput->cd();    

   



   fJetTask = (AliJJetTask*)(man->GetTask( fJetTaskName));


   fJJetJtAnalysis = new AliJJetJtAnalysis(fCard);
   fJJetJtAnalysis->SetJetFinderName(fJetTask->GetJetFinderString());
   fJJetJtAnalysis->SetNumberOfJetFinders(fJetTask->GetNumberOfJetCollections());
   fJJetJtAnalysis->SetNrandom(NRandom);
   fJJetJtAnalysis->SetMoveJet(moveJet);
   fJJetJtAnalysis->SetMC(fDoMC);
   fJJetJtAnalysis->UserCreateOutputObjects();

   fCard->WriteCard(gDirectory);
   fEventHist = new TH1D("EventHist","event numbers",10,0,10);




   PostData( 1, fOutput );


   for( int ij=0;ij< fJetTask->GetNumberOfJetCollections();ij++ ){
	   fJJetJtAnalysis->AddJets( fJetTask->GetAliJJetList( ij ),fJetTask->GetTrackOrMCParticle(ij) );
   }
   fJJetJtAnalysis->SetJTracks(fJetTask->GetJTracks());
   if(fDoMC){
     fJMCTracks = fJetTask->GetMCJTracks();
	   fJJetJtAnalysis->SetMCJTracks(fJetTask->GetMCJTracks());
   }
   fFirstEvent = kTRUE;
}

//______________________________________________________________________________
void AliJJetJtTask::UserExec(Option_t* /*option*/) 
{
	// Processing of one event
	if( fJetTask->GetTaskEntry() != fEntry ) return;
	fJJetJtAnalysis->ClearBeforeEvent();

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
		fJJetJtAnalysis->GetAliJEfficiency()->SetRunNumber( aodEvent->GetRunNumber() );
		fJJetJtAnalysis->GetAliJEfficiency()->Load();
		fFirstEvent = kFALSE;
	}

	if(!IsGoodEvent( event )) return; // zBin is set there

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

	fJJetJtAnalysis->SetCentralityBin( cBin );
	fJJetJtAnalysis->SetCentrality( fcent );
	fJJetJtAnalysis->SetZVertexBin( zBin );
	fJJetJtAnalysis->SetZVertex( zVert );


  if(fDoMC){
    AliMCParticleContainer * mcTracksCont = fJetTask->GetMCParticleContainer("mcparticles");
    if(mcTracksCont){
      TClonesArray * jets = new TClonesArray("AliJJet",1000); // just alias for AliJJet array
      for(int itrack = 0; itrack<mcTracksCont->GetNParticles(); itrack++){
        AliAODMCParticle *track = static_cast<AliAODMCParticle*>(mcTracksCont->GetParticle(itrack));
        if(TMath::Abs(track->Eta()) > 1.0) continue;
        if(track->Pt() > 5){
          if(track->GetNDaughters() > 1){
            /*cout << "Mother: " << endl;
            cout << "Track " << itrack << " Px: " << track->Px() << " Py: " << track->Py() << " Pz: " << track->Pz() << " charge: " << track->Charge() << endl;
            cout << "pT " << track->Pt() << " eta: " << track->Eta() << endl;
            track->Print();
            cout << "Daughters: " << endl;*/
            AliJJet *jet = new AliJJet(track->Px(),track->Py(),track->Pz(), track->E(),track->GetLabel(),0,0);
            FindDaughters(jet,track,mcTracksCont);
            //cout << "Jet pT: " << jet->Pt() << " Eta: " << jet->Eta() << " Phi: " << jet->Phi() << " Constituents: " << jet->GetNConstituents() << endl;
          }
          // If .... new (jets[iJet]) AliJJet(jet->Px(),jet->Py(), jet->Pz(), jet->E(), jet->GetLabel(),0,0);
        }
      }
    }
  }

  fJJetJtAnalysis->UserExec();
  PostData(1, fOutput );
  //PostData(1, fJJetJtAnalysis->GetHistogramsDirectory());

  if(fDebug > 5) cout << "\t------- End UserExec "<<endl;

}

void AliJJetJtTask::FindDaughters(AliJJet *jet, AliAODMCParticle *track, AliMCParticleContainer *mcTracksCont){
  
  for(int id = track->GetFirstDaughter(); id <= track->GetLastDaughter() ; id++){
    AliAODMCParticle *daughter = static_cast<AliAODMCParticle*>(mcTracksCont->GetParticle(id));
    if(daughter->GetNDaughters() > 0){
      FindDaughters(jet,daughter,mcTracksCont);
    }
    if(daughter->IsPrimary() && daughter->Charge() != 0){
      /*cout << "Add constituent" << endl;
      cout << "Track " << id << " Px: " << daughter->Px() << " Py: " << daughter->Py() << " Pz: " << daughter->Pz() << " charge: " << daughter->Charge() << endl;
      cout << "pT " << daughter->Pt() << " eta: " << daughter->Eta() << endl;
      daughter->Print();*/
      AliJMCTrack * particle = (AliJMCTrack*)fJMCTracks->At(id);
      //cout << "Add Constituent: " << id << " pT: " << particle->Pt() << " Eta: " << particle->Eta() << " Phi: " << particle->Phi() << endl;
      //jet->AddConstituent(&fJetTask->GetMCJTracks()[id]);
      jet->AddConstituent(particle);
    }
  }

}

//______________________________________________________________________________
void AliJJetJtTask::Init()
{
  // Intialisation of parameters
  AliInfo("Doing initialization") ; 
  //fJJetJtAnalysis->Init();
}


void AliJJetJtTask::FinishTaskOutput(){

  //TFile::Open("a.root","recreate");
  OpenFile(1);
  fJJetJtAnalysis->WriteHistograms();

}


//________________________________________________________________________
bool AliJJetJtTask::IsGoodEvent(AliVEvent *event) {

  if(fRunTable->IsPP() && fAnaUtils->IsPileUpEvent(event)) {
    return kFALSE;
  } else {
    Bool_t triggeredEventMB = kFALSE; //init

    fEventHist->Fill( 0 );

    Bool_t triggerkMB = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & ( AliVEvent::kMB );

    if( triggerkMB ){
      triggeredEventMB = kTRUE;  //event triggered as minimum bias
      fEventHist->Fill( 1 );
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
        //cout<<zVert<<endl;
        fEventHist->Fill( 2 );
        if(fCard->VertInZRange(zVert)) {
          goodRecVertex = kTRUE;
          fEventHist->Fill( 3 );
          zBin  = fCard->GetBin(kZVertType, zVert);
        }
      }
    }
    return goodRecVertex;
  }
  //---------------------------------
}


//______________________________________________________________________________
void AliJJetJtTask::Terminate(Option_t *)
{
  //fJJetJtAnalysis->WriteHistograms();

}
