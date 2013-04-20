// -------------------------------------------------------------------------
//  INFO
// -------------------------------------------------------------------------

#include <iostream>

// Basic Includes
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THn.h"
#include "TFile.h"
#include "TChain.h"
#include "TObject.h"
#include "TObjArray.h"

// Manager/ Handler
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

// Event pool includes.
#include "AliEventPoolManager.h"

// PID includes.
#include "AliPIDResponse.h"

// AOD includes.
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliAODVertex.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"

// Additional includes.
#include "AliTrackDiHadronPID.h"
#include "AliAODTrackCutsDiHadronPID.h"
#include "AliAODEventCutsDiHadronPID.h"
#include "AliHistToolsDiHadronPID.h"

// AnalysisTask headers.
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskDiHadronPID.h"

using namespace std;

ClassImp(AliAnalysisTaskDiHadronPID);

// -------------------------------------------------------------------------
AliAnalysisTaskDiHadronPID::AliAnalysisTaskDiHadronPID():
	AliAnalysisTaskSE(),
	fPIDResponse(0x0),
	fEventCuts(0x0),
	fTrackCutsTrigger(0x0),
	fTrackCutsAssociated(0x0),
	fPoolMgr(0x0),
	fTriggerTracks(0x0),
	fAssociatedTracks(0x0),
	fCurrentAODEvent(0x0),
	fOutputList(0x0),
	fPtSpectrum(0x0),
	fCorrelations(0x0),
	fMixedEvents(0x0),
	fTOFhistos(0x0),
	fNDEtaBins(32),
	fNDPhiBins(32),	
	fMinNEventsForMixing(5),
	fPoolTrackDepth(2000),
	fPoolSize(1000),
	fMixEvents(kTRUE),
	fMixTriggers(kFALSE),
	fCalculateTOFmismatch(kTRUE),
	fT0Fill(0x0),
	fLvsEta(0x0),
	fLvsEtaProjections(0x0),		
	fDebug(0)

{

	//
	// Default Constructor.
	//

	if (fDebug > 0) {AliInfo("AliAnalysisTaskDiHadronPID Default Constructor.");}		

	for (Int_t iPtClass = 0; iPtClass < 5; iPtClass++) {
			for (Int_t iSpecies = 0; iSpecies < 3; iSpecies++) {
			fCorrelationsTOF[iPtClass][iSpecies] = 0x0;
			fCorrelationsTOFTPC[iPtClass][iSpecies] = 0x0;
		}
	}


}

// -------------------------------------------------------------------------
AliAnalysisTaskDiHadronPID::AliAnalysisTaskDiHadronPID(const char* name):
	AliAnalysisTaskSE(name),
	fPIDResponse(0x0),
	fEventCuts(0x0),
	fTrackCutsTrigger(0x0),
	fTrackCutsAssociated(0x0),
	fPoolMgr(0x0),
	fTriggerTracks(0x0),
	fAssociatedTracks(0x0),	
	fCurrentAODEvent(0x0),
	fOutputList(0x0),
	fPtSpectrum(0x0),
	fCorrelations(0x0),
	fMixedEvents(0x0),
	fTOFhistos(0x0),	
	fNDEtaBins(32),
	fNDPhiBins(32),
	fMinNEventsForMixing(5),
	fPoolTrackDepth(2000),
	fPoolSize(1000),
	fMixEvents(kTRUE),	
	fMixTriggers(kFALSE),
	fCalculateTOFmismatch(kTRUE),
	fT0Fill(0x0),
	fLvsEta(0x0),
	fLvsEtaProjections(0x0),						
	fDebug(0) 

{

	//
	// Named Constructor.
	//

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	for (Int_t iPtClass = 0; iPtClass < 5; iPtClass++) {
			for (Int_t iSpecies = 0; iSpecies < 3; iSpecies++) {
			fCorrelationsTOF[iPtClass][iSpecies] = 0x0;
			fCorrelationsTOFTPC[iPtClass][iSpecies] = 0x0;
		}
	}

	DefineInput(0,TChain::Class());
	DefineOutput(1,TList::Class());

}

// -------------------------------------------------------------------------
AliAnalysisTaskDiHadronPID::~AliAnalysisTaskDiHadronPID() {;

	//
	// Destructor.
	//

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

}

// -------------------------------------------------------------------------
void AliAnalysisTaskDiHadronPID::UserCreateOutputObjects() {

	//
	// Create Output objects.
	//

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// --- BEGIN: Initialization on the worker nodes ---
	AliAnalysisManager* manager = AliAnalysisManager::GetAnalysisManager();
	AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*> (manager->GetInputEventHandler());

	// Getting the pointer to the PID response object.
	fPIDResponse = inputHandler->GetPIDResponse();	

	// Not very neat - only set up for 0-5% analysis.
	Int_t nCentralityBins  = 5;
	Double_t centralityBins[] = {0.,1.,2.,3.,4.,5.};

	Int_t nZvtxBins  = 7;
	Double_t vertexBins[] = {-7.,-5.,-3.,-1.,1.,3.,5.,7.};

	fPoolMgr = new AliEventPoolManager(fPoolSize, fPoolTrackDepth, nCentralityBins, (Double_t*) centralityBins, nZvtxBins, (Double_t*) vertexBins);
    // --- END ---

	// Create the output list.
	fOutputList = new TList();
	fOutputList->SetOwner(kTRUE);

	// Creating all requested histograms locally.
	fEventCuts->CreateHistos();
	fOutputList->Add(fEventCuts);

	fTrackCutsTrigger->CreateHistos();
	fOutputList->Add(fTrackCutsTrigger);

	fTrackCutsAssociated->CreateHistos();
	fOutputList->Add(fTrackCutsAssociated);

	// Get the pT axis for the PID histograms.
	Double_t* ptaxis = fTrackCutsAssociated->GetPtAxisPID();
	Int_t nptbins = fTrackCutsAssociated->GetNPtBinsPID();

	// Create Pt spectrum histogram.
	fPtSpectrum = new TH1F("fPtSpectrum","p_{T} Spectrum;p_{T} (GeV/c);Count",nptbins,ptaxis);
	fOutputList->Add(fPtSpectrum);

	// Create unidentified correlations histogram.
	fCorrelations = AliHistToolsDiHadronPID::MakeHist3D("fCorrelations","Correlations;#Delta#phi;#Delta#eta;p_{T} (GeV/c)",
		fNDPhiBins,-TMath::Pi()/2.,3.*TMath::Pi()/2.,
		fNDEtaBins,-1.6,1.6,
		nptbins, ptaxis);
	fOutputList->Add(fCorrelations);

	// Create unidentified mixed events histogram.
	fMixedEvents = AliHistToolsDiHadronPID::MakeHist3D("fMixedEvents","Mixed Events;#Delta#phi;#Delta#eta;p_{T} (GeV/c)",
		fNDPhiBins,-TMath::Pi()/2.,3.*TMath::Pi()/2.,
		fNDEtaBins,-1.6,1.6,
		nptbins, ptaxis);
	fOutputList->Add(fMixedEvents);

	// Create TOF correlations histograms (DPhi,DEta,pt,TOF).
	fTOFhistos = new TObjArray(15);
	fTOFhistos->SetName("CorrelationsTOF");
	fTOFhistos->SetOwner(kTRUE);

	Int_t nbins[4] = {fNDPhiBins,fNDEtaBins,0,0};
	Double_t min[4] = {-TMath::Pi()/2.,-1.6,0.,0.};
	Double_t max[4] = {3.*TMath::Pi()/2.,1.6,0.,0.};

	for (Int_t iPtClass = 0; iPtClass < 5; iPtClass++) {
		
		nbins[2] = fTrackCutsAssociated->GetNPtBinsPID(iPtClass);
		min[2] = fTrackCutsAssociated->GetPtClassMin(iPtClass);
		max[2] = fTrackCutsAssociated->GetPtClassMax(iPtClass);

		for (Int_t iSpecies = 0; iSpecies < 3; iSpecies++) {

			nbins[3] = fTrackCutsAssociated->GetNTOFbins(iPtClass,iSpecies);
			min[3] = fTrackCutsAssociated->GetTOFmin(iPtClass,iSpecies);
			max[3] = fTrackCutsAssociated->GetTOFmax(iPtClass,iSpecies);

			fCorrelationsTOF[iPtClass][iSpecies] = new THnF(
				Form("fCorrelationsTOF_%i_%i",iPtClass,iSpecies),
				Form("CorrelationsTOF_%i_%i",iPtClass,iSpecies),
				4,nbins,min,max);

			fTOFhistos->Add(fCorrelationsTOF[iPtClass][iSpecies]);

		}
	}

	fOutputList->Add(fTOFhistos);

	// TODO: Create TOF/TPC correlations histogram.

	// Load external TOF histograms if flag is set.
	if (fCalculateTOFmismatch) {LoadExtMismatchHistos();}

	PostData(1,fOutputList);

}

// -------------------------------------------------------------------------
void AliAnalysisTaskDiHadronPID::LocalInit() {

	//
	// Initialize on the client. (or on my computer?? - I think so...)
	//

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}
 	
}

// -------------------------------------------------------------------------
void AliAnalysisTaskDiHadronPID::UserExec(Option_t*) {

	//
	// Main Loop.
	//

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Input Current Event.
	fCurrentAODEvent = dynamic_cast<AliAODEvent*>(InputEvent());
	if (!fCurrentAODEvent) AliFatal("No Event Found!");

	if (!fEventCuts->IsSelected(fCurrentAODEvent)) {return;}

	// Fill the global tracks array. - NOT NEEDED I THINK, since we're not using
	// bit 1<<7 for the associated tracks!

	// Let the track cut objects know that a new event will start.
	fTrackCutsTrigger->StartNewEvent();
	fTrackCutsAssociated->StartNewEvent();

	// Create arrays for trigger/associated particles.
	fTriggerTracks = new TObjArray(350);
	fTriggerTracks->SetOwner(kTRUE);

	fAssociatedTracks = new TObjArray(3500);
	fAssociatedTracks->SetOwner(kTRUE);

	for (Int_t iTrack = 0; iTrack < fCurrentAODEvent->GetNumberOfTracks(); iTrack++) {

		AliAODTrack* track = (AliAODTrack*)fCurrentAODEvent->GetTrack(iTrack);
		AliTrackDiHadronPID* pidtrack = new AliTrackDiHadronPID(track,0x0,0x0,fPIDResponse);
		pidtrack->ForgetAboutPointers();
		pidtrack->SetDebugLevel(fDebug);

		Double_t rndhittime = -1.e21;
		if (fCalculateTOFmismatch) rndhittime = GenerateRandomHit(pidtrack->Eta());

		// Fill p_T spectrum.
		fPtSpectrum->Fill(pidtrack->Pt());

		// Fill the trigger/associated tracks array.
		if (fTrackCutsTrigger->IsSelectedData(pidtrack,rndhittime)) {fTriggerTracks->AddLast(pidtrack);}
		else if (fTrackCutsAssociated->IsSelectedData(pidtrack,rndhittime)) {fAssociatedTracks->AddLast(pidtrack);} 
		else {delete pidtrack;}

	}

	// Fill Correlation histograms.
	for (Int_t iTrigger = 0; iTrigger < fTriggerTracks->GetEntriesFast(); iTrigger++) {
		AliTrackDiHadronPID* triggertrack = (AliTrackDiHadronPID*)fTriggerTracks->At(iTrigger);

		for (Int_t iAssociated = 0; iAssociated < fAssociatedTracks->GetEntriesFast(); iAssociated++) {
			AliTrackDiHadronPID* associatedtrack = (AliTrackDiHadronPID*)fAssociatedTracks->At(iAssociated);

			Double_t DPhi = triggertrack->Phi() - associatedtrack->Phi();
			if (DPhi < -TMath::Pi()/2.) {DPhi += 2.*TMath::Pi();}
			if (DPhi > 3.*TMath::Pi()/2.) {DPhi -= 2.*TMath::Pi();}

			Double_t DEta = triggertrack->Eta() - associatedtrack->Eta();
			fCorrelations->Fill(DPhi,DEta,associatedtrack->Pt());

			Double_t tofhistfill[4] = {DPhi,DEta,associatedtrack->Pt(),-999.};

			for (Int_t iSpecies = 0; iSpecies < 3; iSpecies++) {
				tofhistfill[3] = associatedtrack->GetTOFsignalMinusExpected(iSpecies);

				for (Int_t iPtClass = 0; iPtClass < 5; iPtClass++) {

					// prevent under/over-flow bins to be filled.
					Double_t ptmin = fTrackCutsAssociated->GetPtClassMin(iPtClass);
					Double_t ptmax = fTrackCutsAssociated->GetPtClassMax(iPtClass);
					Double_t apt = associatedtrack->Pt();

					if ( (apt > ptmin) && (apt < ptmax) ) {
						fCorrelationsTOF[iPtClass][iSpecies]->Fill(tofhistfill);
					}

				}
			}
		}
	}

	cout<<"Triggers: "<<fTriggerTracks->GetEntriesFast()<<" Associateds: "<<fAssociatedTracks->GetEntriesFast()<<endl;	

	// Determine centrality of current event.
	TString centralityestimator = fEventCuts->GetCentralityEstimator();
	AliCentrality* currentcentrality = fCurrentAODEvent->GetCentrality();
	Float_t percentile = currentcentrality->GetCentralityPercentile(centralityestimator.Data());

	// Determine vtxz of current event.
	AliAODVertex* currentprimaryvertex = fCurrentAODEvent->GetPrimaryVertex();
	Double_t vtxz = currentprimaryvertex->GetZ();

	AliEventPool* poolin = fPoolMgr->GetEventPool(percentile, vtxz); 
	if (!poolin) {AliFatal(Form("No pool found for centrality = %f, vtxz = %f", percentile, vtxz));}
	// TObjArray* fGlobalTracksArray; 

	// Give a print out of the pool manager's contents.
	if (fDebug > 0) PrintPoolManagerContents();

	// Mix events if there are enough events in the pool.
	if (poolin->GetCurrentNEvents() >= fMinNEventsForMixing) {
		{cout << "Mixing Events." << endl;}

		// Loop over all events in the event pool.
		for (Int_t iMixEvent = 0; iMixEvent < poolin->GetCurrentNEvents(); iMixEvent++) {
	    	TObjArray* mixtracks = poolin->GetEvent(iMixEvent);

	    	// Mix either the triggers or the associateds.
	    	if (fMixTriggers) {

				// Loop over all associateds in this event.
				for (Int_t iAssociated = 0; iAssociated < fAssociatedTracks->GetEntriesFast(); iAssociated++) {
					AliTrackDiHadronPID* associatedtrack = (AliTrackDiHadronPID*)fAssociatedTracks->At(iAssociated);

		    		// Loop over all mixed tracks.
		    		for (Int_t iMixTrack = 0; iMixTrack < mixtracks->GetEntriesFast(); iMixTrack++) {
		    			AliTrackDiHadronPID* mixtrack = (AliTrackDiHadronPID*)mixtracks->At(iMixTrack);
							
						Double_t DPhi = mixtrack->Phi() - associatedtrack->Phi();
						if (DPhi < -TMath::Pi()/2.) {DPhi += 2.*TMath::Pi();}
						if (DPhi > 3.*TMath::Pi()/2.) {DPhi -= 2.*TMath::Pi();}

						Double_t DEta = mixtrack->Eta() - associatedtrack->Eta();
						fMixedEvents->Fill(DPhi,DEta,associatedtrack->Pt());

		    		}
		   		}

		   	} else {

				// Loop over all triggers in this event.
				for (Int_t iTrigger = 0; iTrigger < fTriggerTracks->GetEntriesFast(); iTrigger++) {
					AliTrackDiHadronPID* triggertrack = (AliTrackDiHadronPID*)fTriggerTracks->At(iTrigger);

		    		// Loop over all mixed tracks.
		    		for (Int_t iMixTrack = 0; iMixTrack < mixtracks->GetEntriesFast(); iMixTrack++) {
		    			AliTrackDiHadronPID* mixtrack = (AliTrackDiHadronPID*)mixtracks->At(iMixTrack);
							
						Double_t DPhi = triggertrack->Phi() - mixtrack->Phi();
						if (DPhi < -TMath::Pi()/2.) {DPhi += 2.*TMath::Pi();}
						if (DPhi > 3.*TMath::Pi()/2.) {DPhi -= 2.*TMath::Pi();}

						Double_t DEta = triggertrack->Eta() - mixtrack->Eta();
						fMixedEvents->Fill(DPhi,DEta,mixtrack->Pt());

		    		}
		   		}

		   	} // End if  	
	   	}
	}	

	// Update the event pool.
	AliEventPool* poolout = fPoolMgr->GetEventPool(percentile, vtxz); // Get the buffer associated with the current centrality and z-vtx
	if (!poolout) AliFatal(Form("No pool found for centrality = %f, vtx_z = %f", percentile, vtxz));

	// Q: is it a problem that the fAssociatedTracks array can be bigger than the number of tracks inside?
	if (fMixTriggers) {
		poolout->UpdatePool(fTriggerTracks);
		fAssociatedTracks->Delete();
		delete fAssociatedTracks;
	}
	else {
		poolout->UpdatePool(fAssociatedTracks);
		fTriggerTracks->Delete();
		delete fTriggerTracks;
	}

	fTriggerTracks = 0x0;
	fAssociatedTracks = 0x0;

	// Tell the track cut object that the event is done.
	fTrackCutsTrigger->EventIsDone(0);
	fTrackCutsAssociated->EventIsDone(0);

	PostData(1,fOutputList);

}

// -------------------------------------------------------------------------
Bool_t AliAnalysisTaskDiHadronPID::LoadExtMismatchHistos() {

	//
	// Attempting to load a root file containing information needed
	// to generate random TOF hits.
 	//

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Opening external TOF file.
	if (fDebug > 0) cout<<"Trying to open TOFmismatchHistos.root ..."<<endl;
	TFile* fin = 0x0;
	fin = TFile::Open("alien:///alice/cern.ch/user/m/mveldhoe/rootfiles/TOFmismatchHistos.root");
	if (!fin) {
		AliWarning("Couln't open TOFmismatchHistos, will not calculate mismatches...");
		fCalculateTOFmismatch = kFALSE;
		return kFALSE;
	}

	// Check if the required histograms are present.
	TH1F* tmp1 = (TH1F*)fin->Get("hNewT0Fill");
	if (!tmp1) {
		AliWarning("Couln't find hNewT0Fill, will not calculate mismatches...");
		fCalculateTOFmismatch = kFALSE;
		return kFALSE;	
	}
	TH2F* tmp2 = (TH2F*)fin->Get("hLvsEta");
	if (!tmp2) {
		AliWarning("Couln't find hLvsEta, will not calculate mismatches...");
		fCalculateTOFmismatch = kFALSE;
		return kFALSE;	
	}	

	// Make a deep copy of the files in the histogram.
	fT0Fill = (TH1F*)tmp1->Clone("fT0Fill");
	fLvsEta = (TH2F*)tmp2->Clone("fLvsEta");

	// Close the external file.
	AliInfo("Closing external file.");
	fin->Close();

	// Creating a TObjArray for LvsEta projections.
	const Int_t nbinseta = fLvsEta->GetNbinsX();
	fLvsEtaProjections = new TObjArray(nbinseta);
	fLvsEtaProjections->SetOwner(kTRUE);

	// Making the projections needed (excluding underflow/ overflow).
	for (Int_t iEtaBin = 1; iEtaBin < (nbinseta + 1); iEtaBin++) {
		TH1F* tmp = (TH1F*)fLvsEta->ProjectionY(Form("LvsEtaProjection_%i",iEtaBin),iEtaBin,iEtaBin);
		tmp->SetDirectory(0);
		fLvsEtaProjections->AddAt(tmp,iEtaBin);
	}

	return kTRUE;

}

// -------------------------------------------------------------------------
Double_t AliAnalysisTaskDiHadronPID::GenerateRandomHit(Double_t eta) {

	//
	// Returns a random TOF time.
	//

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	// Default (error) value:
	Double_t rndhittime = -1.e21;

	// TOF mismatch flag is not turned on.
	if (!fCalculateTOFmismatch) {
		AliFatal("Called GenerateRandomHit() method, but flag fCalculateTOFmismatch not set.");
		return rndhittime;
	}

	// TOF doesn't extend much further than 0.8.
	if (TMath::Abs(eta) > 0.8) {
		if (fDebug) {AliInfo("Tried to get a random hit for a track with eta > 0.8.");}
		return rndhittime;
	}

	// Finding the bin of the eta.
	TAxis* etaAxis = fLvsEta->GetXaxis();
	Int_t etaBin = etaAxis->FindBin(eta);
	//cout<<"Eta: "<<eta<<" bin: "<<etaBin<<endl;
	const TH1F* lengthDistribution = (const TH1F*)fLvsEtaProjections->At(etaBin);

	if (!lengthDistribution) {
		AliFatal("length Distribution not found.");
		return rndhittime;
	}

	Double_t currentRndLength = lengthDistribution->GetRandom(); // in cm.

	// Similar to Roberto's code.
	Double_t currentRndTime = currentRndLength / (TMath::C() * 1.e2 / 1.e12);
	Double_t t0fill = -1.26416e+04;
	rndhittime = fT0Fill->GetRandom() - t0fill + currentRndTime;

	return rndhittime;

}

// -------------------------------------------------------------------------
void AliAnalysisTaskDiHadronPID::PrintPoolManagerContents() {

	//
	// Prints out the current contents of the event pool manager.
	//

	// Determine the number of pools in the pool manager.
	AliEventPool* poolin = fPoolMgr->GetEventPool(0,0);
	Int_t NPoolsCentrality = 0;
	while (poolin) {
		NPoolsCentrality++;
		poolin = fPoolMgr->GetEventPool(NPoolsCentrality,0);
	} 

	poolin = fPoolMgr->GetEventPool(0,0);
	Int_t NPoolsVtxZ = 0;	
	while (poolin) {
		NPoolsVtxZ++;
		poolin = fPoolMgr->GetEventPool(0,NPoolsVtxZ);
	} 

	// Loop over all Pools in the matrix of the pool manager.
	cout<<" Pool manager contents: (Nevt,NTrack)"<<endl;
	for (Int_t iCentrality = 0; iCentrality < NPoolsCentrality; iCentrality++) {
		cout<<Form("Centrality Bin: %2i --> ", iCentrality);

		for (Int_t iVtxZ = 0; iVtxZ < NPoolsVtxZ; iVtxZ++) {

			poolin = fPoolMgr->GetEventPool(iCentrality, iVtxZ);

			cout<<Form("(%2i,%4i) ",poolin->GetCurrentNEvents(), poolin->NTracksInPool());

		}

		cout<<endl;
	}

}

// -------------------------------------------------------------------------
void AliAnalysisTaskDiHadronPID::Terminate(Option_t*) {;

	//
	// Called when task is done.
	//

	if (fDebug > 0) {cout << Form("File: %s, Line: %i, Function: %s",__FILE__,__LINE__,__func__) << endl;}

	delete fT0Fill;
	fT0Fill = 0x0;
	delete fLvsEta;
	fLvsEta = 0x0;
	delete fLvsEtaProjections;
	fLvsEtaProjections = 0x0;

}
