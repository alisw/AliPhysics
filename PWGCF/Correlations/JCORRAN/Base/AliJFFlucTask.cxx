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
//#include "TList.h"
//#include "TTree.h"
#include "TFile.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODHandler.h"
#include "AliAODMCParticle.h"
#include "AliJFFlucTask.h" 
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliAODEvent.h" 
#include "AliJTrack.h"
#include "AliJMCTrack.h"
//#include "AliJPhoton.h"
#include "AliJEventHeader.h"
#include "AliJHistManager.h"
#include "AliInputEventHandler.h"

//______________________________________________________________________________
AliJFFlucTask::AliJFFlucTask():
	fInputList(0),
    AliAnalysisTaskSE(), 
    fFFlucAna(0x0),
	fOutput(),
	h_ratio(0)
{
 fEvtNum=0;
 fFilterBit = 0;
 fEta_min = 0;
 fEta_max = 0;
//  DefineOutput(1, TDirectory::Class());
}


//______________________________________________________________________________
AliJFFlucTask::AliJFFlucTask(const char *name, int CollisionCandidates, Bool_t IsMC, Bool_t IsExcludeWeakDecay):
	fInputList(0),
    AliAnalysisTaskSE(name), 
    fFFlucAna(0x0),
	fOutput(),
	h_ratio(0)
{
 DefineOutput(1, TDirectory::Class());
  fEvtNum=0;
 fFilterBit = 0;
 fEta_min = 0;
 fEta_max = 0;
 fTaskName = name;
}

//____________________________________________________________________________
AliJFFlucTask::AliJFFlucTask(const AliJFFlucTask& ap) :
	fInputList(ap.fInputList),
	AliAnalysisTaskSE(ap.GetName()), 
    fFFlucAna(ap.fFFlucAna),
	fOutput(ap.fOutput)
{ 
  AliInfo("----DEBUG AliJFFlucTask COPY ----");
}

//_____________________________________________________________________________
AliJFFlucTask& AliJFFlucTask::operator = (const AliJFFlucTask& ap)
{
  // assignment operator
  AliInfo("----DEBUG AliJFFlucTask operator= ----");
  this->~AliJFFlucTask();
  new(this) AliJFFlucTask(ap);
  return *this;
}

//______________________________________________________________________________
AliJFFlucTask::~AliJFFlucTask()
{
		// destructor 
		delete fFFlucAna;
		delete fInputList;
		delete fOutput;
		delete h_ratio;
}

//________________________________________________________________________
void AliJFFlucTask::UserCreateOutputObjects()
{ 
  fFFlucAna =  new AliJFFlucAnalysis( fTaskName );
  fFFlucAna->SetDebugLevel(fDebugLevel); 
  //=== create the jcorran outputs objects

  //=== Get AnalysisManager
  /*
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if(!man->GetOutputEventHandler()) {
    Fatal("UserCreateOutputObjects", "This task needs an AOD handler");
    return;
  }
  */
  fInputList = new TClonesArray("AliJBaseTrack" , 2500);
  fInputList->SetOwner(kTRUE);


  OpenFile(1);
  fOutput = gDirectory;
  fOutput->cd();
  fFFlucAna->UserCreateOutputObjects(); 

  PostData(1, fOutput);
}

//______________________________________________________________________________
void AliJFFlucTask::UserExec(Option_t* /*option*/) 
{
	// Processing of one event
	if(!((Entry()-1)%100))  AliInfo(Form(" Processing event # %lld",  Entry()));
	// initiaizing variables from last event 
	fInputList->Clear();
	
	float fCent = -999;
	double fvertex[3] = {-999,-999,-999};	

	fEvtNum++;
//	if(fEvtNum % 100 == 0){ cout << "evt : " << fEvtNum <<endl ;}

	// load current event and save track, event info 
	AliAODEvent *currentEvent = dynamic_cast<AliAODEvent*>(InputEvent());
	fCent = ReadAODCentrality( currentEvent, "V0M"  ) ; 

	if( IsGoodEvent( currentEvent )){

			ReadAODTracks( currentEvent, fInputList ) ; // read tracklist
			ReadVertexInfo( currentEvent, fvertex); // read vertex info
			Int_t Ntracks = fInputList->GetEntriesFast();

			// Analysis Part 
			fFFlucAna->Init();
			fFFlucAna->SetInputList( fInputList ); 
			fFFlucAna->SetEventCentrality( fCent );
			fFFlucAna->SetEventVertex( fvertex );
			fFFlucAna->SetEtaRange( fEta_min, fEta_max ) ;
			fFFlucAna->UserExec(""); // doing some analysis here. 
			// 
	}
}

//______________________________________________________________________________
void AliJFFlucTask::ReadVertexInfo ( AliAODEvent *aod, double*  fvtx )
{
	fvtx[0] = aod->GetPrimaryVertex()->GetX();
	fvtx[1] = aod->GetPrimaryVertex()->GetY();
	fvtx[2] = aod->GetPrimaryVertex()->GetZ();
}
//______________________________________________________________________________
float AliJFFlucTask::ReadAODCentrality( AliAODEvent *aod, TString Trig )
{
	AliCentrality *cent = aod->GetCentrality();
	if(!cent){ Printf("No Cent event"); return -999; };
	return cent->GetCentralityPercentile(Trig.Data()); 
}
//______________________________________________________________________________
void AliJFFlucTask::Init()
{
	// initiationg all parameters and variables

	// Intialisation of parameters
	AliInfo("Doing initialization") ; 
//
}
//______________________________________________________________________________
void AliJFFlucTask::ReadAODTracks( AliAODEvent *aod , TClonesArray *TrackList)
{
		//aod->Print();
	if( IsMC == kTRUE ){  // how to get a flag to check  MC or not ! 

		TClonesArray *mcArray = (TClonesArray*) aod->FindListObject(AliAODMCParticle::StdBranchName());
		if(!mcArray){ Printf("Error not a proper MC event"); };  // check mc array
		
		Int_t nt = mcArray->GetEntriesFast();
		Int_t ntrack =0;
		for( int it=0; it< nt ; it++){
				AliAODMCParticle *track = (AliAODMCParticle*)mcArray->At(it);
				if(!track) { Error("ReadEventAODMC", "Could not receive particle %d",(int) it); continue; };
				if( track->IsPhysicalPrimary() ){
						// insert AMTP weak decay switch here
						if(IsExcludeWeakDecay == kTRUE){
								//cout << "finding weak decaying particle ... " << endl;
								Bool_t kExcludeParticle = kFALSE;
								Int_t gMotherIndex = track->GetMother(); // check and ask about this to DJ changed to mother from firstmother
								if(gMotherIndex != -1) {
//								cout << "this mother is " << gMotherIndex << endl;
										AliAODMCParticle *motherParticle = (AliAODMCParticle*)mcArray->At(gMotherIndex);
//										cout << "mother pdg code is " << motherParticle->GetPdgCode() << endl;
										if(motherParticle) {
												if(IsThisAWeakDecayingParticle(motherParticle)){
//													cout << "!!!!!!! FOUND!!!!!!!" <<endl;
													kExcludeParticle = kTRUE;
												}
										}
								}
								//Exclude from the analysis decay products of weakly decaying particles
								if(kExcludeParticle) continue;
						} // weak decay particles are excluded
						double track_abs_eta = TMath::Abs( track->Eta() );
						double Pt = track->Pt(); 
						if( fPt_min > 0){
							if( track->Pt() < fPt_min || track->Pt() > fPt_max ) continue ; // pt cut
						}
						Int_t pdg = track->GetPdgCode();
						Char_t ch = (Char_t) track->Charge();
						Int_t label = track->GetLabel();
						AliJBaseTrack *itrack = new ((*TrackList)[ntrack++])AliJBaseTrack;
						itrack->SetLabel( label );
						itrack->SetParticleType( pdg);
						itrack->SetPxPyPzE( track->Px(), track->Py(), track->Pz(), track->E() );
						itrack->SetCharge(ch) ;
						//}; no eta cut in task file
				}
		}
	} // read mc track done.
	else if( IsMC == kFALSE){
		Int_t nt = aod->GetNumberOfTracks();
		Int_t ntrack =0;
		for( int it=0; it<nt ; it++){
				AliAODTrack *track = dynamic_cast<AliAODTrack*>(aod->GetTrack(it));
				if(!track) { Error("ReadEventAOD", "Could not receive partice %d", (int) it); continue; };
				if(track->TestFilterBit( fFilterBit )){ // hybrid cut
						double Pt = track->Pt(); 
						if( fPt_min > 0){
							if( track->Pt() < fPt_min || track->Pt() > fPt_max ) continue ; // pt cut
						}

//						double track_abs_eta = TMath::Abs( track->Eta() );
//						if( track_abs_eta > fEta_min && track_abs_eta < fEta_max){ // eta cut
								AliJBaseTrack *itrack = new( (*TrackList)[ntrack++])AliJBaseTrack;
								//itrack->SetID( track->GetID() );
								itrack->SetID( TrackList->GetEntriesFast() ) ;
								itrack->SetPxPyPzE( track->Px(), track->Py(), track->Pz(), track->E() );
								itrack->SetParticleType(kJHadron);
								itrack->SetCharge(track->Charge() );
								itrack->SetStatus(track->GetStatus() );
								//						} no eta cut in task file
				}
		}
	} //read aod reco track done.
}
//______________________________________________________________________________
Bool_t AliJFFlucTask::IsGoodEvent( AliAODEvent *event){

		//event selection here!
		Bool_t Event_status = kFALSE;

		// check vertex 
		AliVVertex *vtx = event->GetPrimaryVertex();
		if(vtx){
				if( vtx->GetNContributors() > 0 ){
						double zVert = vtx->GetZ();
						if( zVert > -10 && zVert < 10) Event_status = kTRUE;
				}
		}

		// event cent flatting  --- do it only when IsCentFlat is true
		if( Event_status== kTRUE && IsCentFlat == kTRUE ){
				float centrality = ReadAODCentrality( event, "V0M");
				double cent_flat_ratio = h_ratio->GetBinContent( (h_ratio->GetXaxis()->FindBin(centrality))) ;
				if (gRandom->Uniform(0, cent_flat_ratio) > 1 ){
				//	cout << "this event is excluded with random pro" << endl;
					Event_status = kFALSE;
				}
				
		}

		return Event_status;
		/*
		//	Bool_t Trigger_kCentral = kFALSE; // init
		//	Trigger_kCentral =(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))
		->IsEventSelected() & AliVEvent::kCentral);
		//	return Trigger_kCentral;
		//	*/
}
//
//______________________________________________________________________________
void AliJFFlucTask::Terminate(Option_t *)
{
		//    fFFlucAna->Terminate("");
		// Processing when the event loop is ended
		cout<<"AliJFFlucTask Analysis DONE !!"<<endl; 
}


//______________________________________________________________________________
Bool_t AliJFFlucTask::IsThisAWeakDecayingParticle(AliAODMCParticle *thisGuy)
{
		// In order to prevent analyzing daughters from weak decays 
		// - AMPT does not only strong decays, so IsPhysicalPrimary does not catch it
		Int_t pdgcode = TMath::Abs( thisGuy->GetPdgCode() );
		Int_t myWeakParticles[7] = { 3322, 3312, 3222, // Xi0 Xi+- Sigma-+
				3122, 3112, // Lambda0 Sigma+-
				130, 310 // K_L0 K_S0
		};
		Bool_t found = kFALSE;
		for(Int_t i=0; i!=7; ++i)
				if( myWeakParticles[i] == pdgcode ) {
						found = kTRUE;
						break;
				}
		return found;
}
//______________________________________________________________________________
void AliJFFlucTask::SetIsCentFlat( Bool_t isCentFlat )
{
		cout << "Setting to flatting Centrality with LHC11h data.. " << endl;
		if( IsMC == kTRUE ) Printf("this is not LHC11h data");
		IsCentFlat = isCentFlat;

		if(IsCentFlat == kTRUE){
				// ratio from ExtractCentRatio.C // do not change this manually	
				h_ratio = new TH1D("h_ratio","",60,0,60);
				h_ratio->SetBinContent(1,1);
				h_ratio->SetBinContent(2,1.04705);
				h_ratio->SetBinContent(3,1.06116);
				h_ratio->SetBinContent(4,1.06961);
				h_ratio->SetBinContent(5,1.06418);
				h_ratio->SetBinContent(6,1.32022);
				h_ratio->SetBinContent(7,1.33321);
				h_ratio->SetBinContent(8,1.29442);
				h_ratio->SetBinContent(9,1.22629);
				h_ratio->SetBinContent(10,1);
				h_ratio->SetBinContent(11,2.53223);
				h_ratio->SetBinContent(12,1.43706);
				h_ratio->SetBinContent(13,1.04602);
				h_ratio->SetBinContent(14,1.02736);
				h_ratio->SetBinContent(15,1.02187);
				h_ratio->SetBinContent(16,1);
				h_ratio->SetBinContent(17,1.02364);
				h_ratio->SetBinContent(18,1.00895);
				h_ratio->SetBinContent(19,1.00213);
				h_ratio->SetBinContent(20,1.03277);
				h_ratio->SetBinContent(21,1.02793);
				h_ratio->SetBinContent(22,1.05164);
				h_ratio->SetBinContent(23,1.03922);
				h_ratio->SetBinContent(24,1.04664);
				h_ratio->SetBinContent(25,1.01487);
				h_ratio->SetBinContent(26,1.01047);
				h_ratio->SetBinContent(27,1.03473);
				h_ratio->SetBinContent(28,1.04211);
				h_ratio->SetBinContent(29,1.00174);
				h_ratio->SetBinContent(30,1);
				h_ratio->SetBinContent(31,1.04261);
				h_ratio->SetBinContent(32,1.01121);
				h_ratio->SetBinContent(33,1);
				h_ratio->SetBinContent(34,1.00863);
				h_ratio->SetBinContent(35,1.02587);
				h_ratio->SetBinContent(36,1.00558);
				h_ratio->SetBinContent(37,1.02542);
				h_ratio->SetBinContent(38,1.00648);
				h_ratio->SetBinContent(39,1.03365);
				h_ratio->SetBinContent(40,1.0166);
				h_ratio->SetBinContent(41,1.01691);
				h_ratio->SetBinContent(42,1.00623);
				h_ratio->SetBinContent(43,1);
				h_ratio->SetBinContent(44,1.02026);
				h_ratio->SetBinContent(45,1.0043);
				h_ratio->SetBinContent(46,1.03178);
				h_ratio->SetBinContent(47,1.00555);
				h_ratio->SetBinContent(48,1.01814);
				h_ratio->SetBinContent(49,1.00046);
				h_ratio->SetBinContent(50,1.00857);
				h_ratio->SetBinContent(51,13.4844);
				h_ratio->SetBinContent(52,12.6347);
				h_ratio->SetBinContent(53,11.2659);
				h_ratio->SetBinContent(54,8.80011);
				h_ratio->SetBinContent(55,5.4266);
				h_ratio->SetBinContent(56,2.38822);
				h_ratio->SetBinContent(57,1.12582);
				h_ratio->SetBinContent(58,1.01973);
				h_ratio->SetBinContent(59,1);
				h_ratio->SetBinContent(60,1.03022);
		}
}
//
//

