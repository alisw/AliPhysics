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
// Analysis task for jT Analysis
// author: B.K Kim,  D.J. Kim, T. Snellman
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
#include "AliMultSelection.h"
#include "AliAnalysisManager.h"
#include "AliJHistManager.h"
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
  fSelector(""),
  fJJetJtAnalysis(0x0),
  fOutput(NULL),
  fCard(NULL),
  fFirstEvent(kTRUE),
  cBin(-1),
  zBin(-1),
  fDoMC(0),
  fSide(0),
  fDoLog(0),
  NRandom(1),
  moveJet(0),
  fCentCut(100),
  zVert(-999),
  fAnaUtils(NULL),
  fRunTable(NULL),
  fJMCTracks(NULL),
  fEventHist()

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
  fSelector(""),
  fJJetJtAnalysis(0x0),
  fOutput(NULL),
  fCard(NULL),
  fFirstEvent(kTRUE),
  cBin(-1),
  zBin(-1),
  fDoMC(0),
  fSide(0),
  fDoLog(0),
  NRandom(1),
  moveJet(0),
  fCentCut(100),
  zVert(-999),
  fAnaUtils(NULL),
  fRunTable(NULL),
  fJMCTracks(NULL),
  fEventHist()
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
  fSelector(ap.fSelector),
  fJJetJtAnalysis( ap.fJJetJtAnalysis ),
  fOutput( ap.fOutput ),
  fCard(ap.fCard),
  fFirstEvent(ap.fFirstEvent),
  cBin(ap.cBin),
  zBin(ap.zBin),
  fSide(ap.fSide),
  NRandom(ap.NRandom),
  moveJet(ap.moveJet),
  fCentCut(ap.fCentCut),
  zVert(ap.zVert),
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
  delete fAnaUtils;

}

//________________________________________________________________________
void AliJJetJtTask::UserCreateOutputObjects()
{  
  //=== create the jcorran outputs objects
  if(fDebug > 1) printf("AliJJetJtTask::UserCreateOutPutData() \n");

  fAnaUtils = new AliAnalysisUtils();
  fAnaUtils->SetUseOutOfBunchPileUp(kTRUE);

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
  if(fDoLog) fJJetJtAnalysis->SetLog(fDoLog);
  fJJetJtAnalysis->SetLeadingJets(fLeadingJets);
  fJJetJtAnalysis->SetMaxDeltaRCorr(fmaxDeltaRCorr);
  fJJetJtAnalysis->SetSide(fSide);
  fJJetJtAnalysis->SetnR(fJetTask->GetnR());
  fJJetJtAnalysis->Setnkt(fJetTask->Getnkt());
  fJJetJtAnalysis->UserCreateOutputObjects();

  fCard->WriteCard(gDirectory);
  fEventHist << TH1D("EventHist","event numbers",10,0,10) << "END";




  PostData( 1, fOutput );


  for( int ij=0;ij< fJetTask->GetNumberOfJetCollections();ij++ ){
    fJJetJtAnalysis->AddJets( fJetTask->GetAliJJetList( ij ),fJetTask->GetTrackOrMCParticle(ij), fJetTask->GetConeSize(ij));
  }
  fJJetJtAnalysis->SetJTracks(fJetTask->GetJTracks());
  fJJetJtAnalysis->SetIncludeFullJets(fJetTask->GetIncludeFullJets());
  if(fDoMC){
    fJMCTracks = fJetTask->GetMCJTracks();
    fJJetJtAnalysis->SetMCJTracks(fJetTask->GetMCJTracks());
  }
  fFirstEvent = kTRUE;
}


//______________________________________________________________________________
/// Primary loop
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
	if(fDebug > 6) cout << fRunTable->GetPeriodName() << endl;
	sel = (AliMultSelection*) InputEvent() -> FindListObject("MultSelection");
	if (sel) {
		if(!fSelector.IsNull()){ //If centrality selector is set in wagon configuration, otherwise default to V0A
			fcent = sel->GetMultiplicityPercentile(fSelector.Data());
		}else{
			fcent = sel->GetMultiplicityPercentile("V0A");
		}
	}
	else {
		if(fDebug > 2) cout << "Sel not found" << endl;
		fcent = -1;
	}

  if(fcent > fCentCut){
    if(fDebug > 2) cout << "Skip event, Centrality was " << fcent << " and cut is " << fCentCut << endl;
    return;
  }

  cBin = fCard->GetBin(kCentrType, fcent);;

  //if(cBin<0) return; //TODO CHECK

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

/// Obsolete?
/// Find daughters of given track
/// If daughters are primary tracks, add them as constituents to the original jet
/// If they are not primary then keep searching for daughters
/// \param jet Pointer to jet, i.e. the original track
/// \param track Find daughters of this track
/// \param mcTracksCont Monte carlo track container
/// 
void AliJJetJtTask::FindDaughters(AliJJet *jet, AliAODMCParticle *track, AliMCParticleContainer *mcTracksCont){

  for(int id = track->GetDaughterFirst(); id <= track->GetDaughterLast() ; id++){
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
/// Function to initialize the parameters
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
    // Make sure if V0OR or whayt

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
