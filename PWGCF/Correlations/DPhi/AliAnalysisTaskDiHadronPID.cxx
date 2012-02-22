// ----------------------------------------------------------------------------
// This class makes di-hadron correlations, with TOF and TPC signals for
// the associated particles. It runs on AOD049 and AOD73 productions.
// ----------------------------------------------------------------------------
// Author: Misha Veldhoen (misha.veldhoen@cern.ch)
// Start: July 21st, 2011.
// Last edit: Feb. 17th 2012. (v08)
// ----------------------------------------------------------------------------
//

#include <iostream>
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TFile.h"

#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODVertex.h"
//#include "AliAODPid.h"

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"
//#include "AliTOFPIDResponse.h"

#include "AliAnalysisTaskDiHadronPID.h"

using namespace std;

ClassImp(AliAnalysisTaskDiHadronPID);

//_____________________________________________________________________________
AliAnalysisTaskDiHadronPID::AliAnalysisTaskDiHadronPID():
	AliAnalysisTaskSE(),
	fPIDResponse(0x0),
	fAODEvent(0x0),
	fAODHeader(0x0),
	fAODVertex(0x0),
	fAODTrack(0x0),
	fPIDPartners(0x0),
	fCentrality(0x0),
	fVertexZ(0x0),
    fTrackCuts(0x0),
	fPtSpectrum(0x0),
	fAssociatedDistribution(0x0),
	fTPCnSigmaProton(0x0),
	fTOFnSigmaProton(0x0),
	fTPCnSigmaPion(0x0),
	fTOFnSigmaPion(0x0),
	fTPCnSigmaKaon(0x0),
	fTOFnSigmaKaon(0x0),
	fTPCSignal(0x0),
	fTOFSignal(0x0),
	fDiHadron(0x0),
	fMixedEvents(0x0),
	fHistoList(0x0),
	fVerbose(kFALSE),
	fMask(0),	
	fCalculateMixedEvents(kFALSE),
	fTrigBufferIndex(0),	
	fTrigBufferSize(0)

{
	//
	// Default Constructor.
	//
	
	// Trigger buffer.
	for(Int_t i=0; i<25000; i++) {
		for(Int_t j=0; j<4; j++) {
			fTrigBuffer[i][j]=0;
		}				
	}	

	
	fMask = 1<<7;
	
    // The identified di-hadron correlations.
	for (Int_t i = 0; i < 3; i++) {
		for (Int_t j = 0; j < 10; j++) {
            fDiHadronTPC[i][j]=0x0;
            fDiHadronTOF[i][j]=0x0;
			fDiHadronTPCTOF[i][j]=0x0;
		}
	}
    
}

//_____________________________________________________________________________
AliAnalysisTaskDiHadronPID::AliAnalysisTaskDiHadronPID(const char *name):
	AliAnalysisTaskSE(name),
	fPIDResponse(0x0),
	fAODEvent(0x0),
	fAODHeader(0x0),
	fAODVertex(0x0),
	fAODTrack(0x0),
	fPIDPartners(0x0),
	fCentrality(0x0),
	fVertexZ(0x0),
    fTrackCuts(0x0),
	fPtSpectrum(0x0),
    fAssociatedDistribution(0x0),
	fTPCnSigmaProton(0x0),
	fTOFnSigmaProton(0x0),
	fTPCnSigmaPion(0x0),
	fTOFnSigmaPion(0x0),
	fTPCnSigmaKaon(0x0),
	fTOFnSigmaKaon(0x0),
	fTPCSignal(0x0),
	fTOFSignal(0x0),
	fDiHadron(0x0),
	fMixedEvents(0x0),
	fHistoList(0x0),
	fVerbose(kFALSE),
	fMask(0),
	fCalculateMixedEvents(kFALSE),
	fTrigBufferIndex(0),	
	fTrigBufferSize(0)

{
	//
	// Named Constructor.
	//
	
	DefineInput(0, TChain::Class());
	DefineOutput(1, TList::Class());
	
	// Trigger buffer.
	for(Int_t i=0; i<25000; i++) {
		for(Int_t j=0; j<4; j++) {
			fTrigBuffer[i][j]=0;
		}				
	}	
	
	fMask = 1<<7;
    
    // The identified di-hadron correlations.
	for (Int_t i = 0; i < 3; i++) {
		for (Int_t j = 0; j < 10; j++) {
            fDiHadronTPC[i][j]=0x0;
            fDiHadronTOF[i][j]=0x0;
			fDiHadronTPCTOF[i][j]=0x0;
		}
	}

}

//_____________________________________________________________________________
AliAnalysisTaskDiHadronPID::~AliAnalysisTaskDiHadronPID() {

	//
	// Destructor.
	//
	
    if(fPIDPartners) {
        delete fPIDPartners;
        fPIDPartners=0;
    }

}

//_____________________________________________________________________________
void AliAnalysisTaskDiHadronPID::UserCreateOutputObjects() 

{
	//
	// Output objects.
	//
    
	// The Analysis Manager.
	AliAnalysisManager* manager = AliAnalysisManager::GetAnalysisManager();
	AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*> (manager->GetInputEventHandler());

	if(!inputHandler) {
	  AliFatal("Error getting AliInputEventHandler");
	}

	// Pointers to PID Response objects.	
	fPIDResponse = inputHandler->GetPIDResponse(); 
	cout << "PID Response object: " << fPIDResponse << endl;

	// Create the output of the task.	
	fHistoList = new TList();
	fHistoList->SetOwner(kTRUE); 
	
	// Ranges in dPhi, dEta, and the number of bins for di-hadron correlations.
	Int_t binsHisto[4] = {36,25};
	Double_t minHisto[4] = {-TMath::Pi()/2.,-1.8};
	Double_t maxHisto[4] = {3.*TMath::Pi()/2,1.8};
	
    // --- EVENT SAMPLE PLOTS (PER EVENT) ---
    
	// Centrality Histogram.
	fCentrality = new TH1F("fCentrality","Centrality;Centrality;N_{evt}",10,0,100);
	fHistoList->Add(fCentrality);
	
	// Vertex Z Histogram.
	fVertexZ = new TH1F("fVertexZ","Vertex Z position;z (cm);N_{evt}",30,-15,15);
	fHistoList->Add(fVertexZ);
	
	// --- UNIDENTIFIED SPECTRA ---
	
	fPtSpectrum = new TH1F("fPtSpectrum","p_{T} Spectrum Unidentified;p_{T}",100,0,50);
	fHistoList->Add(fPtSpectrum);
	
	fAssociatedDistribution = new TH3F("fAssociatedDistribution","Associated Distribution;#phi;#eta;p_{T_assoc}",binsHisto[0],0.,2.*TMath::Pi(),binsHisto[1],-0.9,0.9,10,0.,5.);
    fAssociatedDistribution->Sumw2();
    fHistoList->Add(fAssociatedDistribution);
    
    // --- TRACK CUTS ---
    
    fTrackCuts = new TH1F("fTrackCuts","Track Cuts;Cut;Count",4,0,4);
    (fTrackCuts->GetXaxis())->SetBinLabel(1,"Standard Cuts");
    (fTrackCuts->GetXaxis())->SetBinLabel(2,"+ TPCpid");
    (fTrackCuts->GetXaxis())->SetBinLabel(3,"+ TOFpid");
    (fTrackCuts->GetXaxis())->SetBinLabel(4,"+ no TOFmismatch");
    fHistoList->Add(fTrackCuts);
	
    // --- QA PLOTS PID ---
    
	fTPCnSigmaProton = new TH2F("fTPCnSigmaProton","n#sigma plot for the TPC (Protons);p (GeV/c);n#sigma",100,0.,5.,100,-10.,10.);
	fHistoList->Add(fTPCnSigmaProton);
	
	fTOFnSigmaProton = new TH2F("fTOFnSigmaProton","n#sigma plot for the TOF (Protons);p (GeV/c);n#sigma",100,0.,5.,100,-10.,10.);
	fHistoList->Add(fTOFnSigmaProton);

	fTPCnSigmaPion = new TH2F("fTPCnSigmaPion","n#sigma plot for the TPC (Pions);p (GeV/c);n#sigma",100,0.,5.,100,-10.,10.);
	fHistoList->Add(fTPCnSigmaPion);
	
	fTOFnSigmaPion = new TH2F("fTOFnSigmaPion","n#sigma plot for the TOF (Pions);p (GeV/c);n#sigma",100,0.,5.,100,-10.,10.);
	fHistoList->Add(fTOFnSigmaPion);
	
	fTPCnSigmaKaon = new TH2F("fTPCnSigmaKaon","n#sigma plot for the TPC (Kaons);p (GeV/c);n#sigma",100,0.,5.,100,-10.,10.);
	fHistoList->Add(fTPCnSigmaKaon);
	
	fTOFnSigmaKaon = new TH2F("fTOFnSigmaKaon","n#sigma plot for the TOF (Kaons);p (GeV/c);n#sigma",100,0.,5.,100,-10.,10.);
	fHistoList->Add(fTOFnSigmaKaon);	
	
	// --- PID SIGNALS ---
	
	fTPCSignal = new TH3F("fTPCSignal","TPC Signal;#eta;p_{T} GeV/c;dE/dx",25,-.9,.9,100.,0.,5.,150,0.,300.);
	fHistoList->Add(fTPCSignal);
	
	fTOFSignal = new TH3F("fTOFSignal","TOF Signal;#eta;p_{T} GeV/c;t",25,-.9,.9,100.,0.,5.,150,10000.,25000.);
	fHistoList->Add(fTOFSignal);
	
	// --- UNIDENTIFIED DI-HADRON CORRELATIONS & MIXED EVENTS ---
		
	// Di Hadron Correlations, unidentified. (Dphi,Deta,p_T_assoc)
	fDiHadron = new TH3F("fDiHadron","Di-Hadron Correlations;#Delta#phi;#Delta#eta;p_{T,assoc}",binsHisto[0],minHisto[0],maxHisto[0],binsHisto[1],minHisto[1],maxHisto[1],10,0.,5.);
	fHistoList->Add(fDiHadron);
		
    if (fCalculateMixedEvents) {
        fMixedEvents = new TH3F("fMixedEvents","Mixed Events;#Delta#phi;#Delta#eta;p_{T_assoc}",binsHisto[0],minHisto[0],maxHisto[0],binsHisto[1],minHisto[1],maxHisto[1],10,0.,5.);
        fMixedEvents->Sumw2();
        fHistoList->Add(fMixedEvents);
    }
		
    // --- DI-HADRON CORRELATIONS WITH TPC AND TOF SIGNALS ---
    
	// Di Hadron Correlations with two PID signals, for each p_T bin and particle species separately. (i.e. 30 histo's) 
	// Axes: {Dphi, Deta, TPC signal, TOF signal}
    TString basenameTPC("fDiHadronTPC");
    TString basenameTOF("fDiHadronTOF");
	TString basenameTPCTOF("fDiHadronTPCTOF");
	TString basetitle("Di-Hadron correlation");
	TString finalname, finaltitle;
	TString species[3] = {"Pion","Kaon","Proton"};
	TString ptbins[11] = {"0.0","0.5","1.0","1.5","2.0","2.5","3.0","3.5","4.0","4.5","5.0"};
	
	Int_t binsTPC[3][10] = {{100,100,100,100,100,100,100,100,100,100},
                            {100,100,100,100,100,100,100,100,100,100},
                            {100,100,100,100,100,100,100,100,100,100}};
		
	Double_t minTPC[3][10] =    {{   -50., -50.,-30.,-30.,-30.,-30.,-30.,-30.,-30.,-30.},
                                 {	-100., -50.,-20.,-20.,-20.,-20.,-20.,-20.,-20.,-20.},
                                 {  -200.,-150.,-50.,-30.,-20.,-20.,-20.,-20.,-20.,-20.}};
	Double_t maxTPC[3][10] =    {{ 150., 100., 50., 30., 25., 25., 25., 25., 25., 25.},
                                 {  50., 100., 40., 30., 30., 30., 30., 30., 30., 30.},
                                 { 100.,  50., 50., 30., 30., 30., 30., 30., 30., 30.}};
	
	Int_t binsTOF[3][10] = {{100,100,100,100,100,100,100,100,100,100},
                            {100,100,100,100,100,100,100,100,100,100},
                            {100,100,100,100,100,100,100,100,100,100}};
	
	Double_t minTOF[3][10] =    {{-2000.,-2000.,-1000., -500., -500., -500., -500.,	-500., -500., -500.},
                                 {-2000.,-2000.,-2000.,-1000.,-1000.,-1000.,-1000.,-1000.,-1000.,-1000.},
                                 {-2000.,-2000.,-6000.,-3000.,-2000.,-1500.,-1500.,-1500.,-1500.,-1500.}};
	Double_t maxTOF[3][10] =    {{ 2000., 2000., 5000., 3000., 2000., 1500., 1500., 1500., 1500., 1500.},
                                 { 2000., 2000., 5000., 2000., 1500., 1500., 1500., 1500., 1500., 1500.},
                                 { 2000., 2000., 2000., 1000., 1000., 1000., 1000., 1000., 1000., 1000.}};
	
	// Recall that AliPID::kPion = 2, AliPID::kKaon = 3, AliPID::kProton = 4.
	for (Int_t iSpecies = 0; iSpecies < 3; iSpecies++) {
				
		for (Int_t iPtBin = 0; iPtBin < 10; iPtBin++) {
			
			// setting the variables for each histogram.
			finaltitle = basetitle;
			(((((finaltitle += " (") += species[iSpecies]) += ") ") += ptbins[iPtBin]) += " < P_t < ") += ptbins[iPtBin+1];
            finaltitle+=";#Delta#phi;#Delta#eta;TPC signal;TOF signal;";
            
			binsHisto[2] = binsTPC[iSpecies][iPtBin];
			binsHisto[3] = binsTOF[iSpecies][iPtBin];
			minHisto[2] = minTPC[iSpecies][iPtBin];
			minHisto[3] = minTOF[iSpecies][iPtBin];
			maxHisto[2] = maxTPC[iSpecies][iPtBin];
			maxHisto[3] = maxTOF[iSpecies][iPtBin];
                
            // Make the di-hadron correlations with different pid signals.
            finalname = basenameTPC;
			(((finalname += "_") += species[iSpecies]) += "_") += iPtBin;
            fDiHadronTPC[iSpecies][iPtBin] = new TH3F(finalname,finaltitle,binsHisto[0],minHisto[0],maxHisto[0],binsHisto[1],minHisto[1],maxHisto[1],binsHisto[2],minHisto[2],maxHisto[2]);
            fHistoList->Add(fDiHadronTPC[iSpecies][iPtBin]);

            finalname = basenameTOF;
			(((finalname += "_") += species[iSpecies]) += "_") += iPtBin;
            fDiHadronTOF[iSpecies][iPtBin] = new TH3F(finalname,finaltitle,binsHisto[0],minHisto[0],maxHisto[0],binsHisto[1],minHisto[1],maxHisto[1],binsHisto[3],minHisto[3],maxHisto[3]);
            fHistoList->Add(fDiHadronTOF[iSpecies][iPtBin]);

            finalname = basenameTPCTOF;
			(((finalname += "_") += species[iSpecies]) += "_") += iPtBin;
            fDiHadronTPCTOF[iSpecies][iPtBin] = new THnSparseF(finalname,finaltitle,4,binsHisto,minHisto,maxHisto);
            fHistoList->Add(fDiHadronTPCTOF[iSpecies][iPtBin]);
            
		}
	}
    
	PostData(1, fHistoList);
	
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskDiHadronPID::SelectTrack(AliAODTrack *track, Int_t cuts)

{
	
    // This selects tracks with three different kinds of track cuts:
    //
    //  Case 0: 
    //   - Tracks must pass the filterbit,
    //   - Tracks must have eta < 0.9.
    //  Case 1:
    //   - TPCpid hit is demanded.
    //  Case 2:
    //   - TOFpid hit is demanded.
    //  Case 3:
    //   - no TOFmismatch is demanded.
    //
    
    ULong_t status;
    
    switch (cuts) {
        case 0:
            //cout<<"test0"<<endl;
            if (!(track->TestFilterMask(fMask))) {return kFALSE;}
            if (!(TMath::Abs(track->Eta())<0.9)) {return kFALSE;}
            break;
        case 1:
            //cout<<"test1"<<endl;
            status=GetTrackPartner(track)->GetStatus();  
            if (!((status&AliAODTrack::kTPCpid)==AliAODTrack::kTPCpid)) {return kFALSE;}
            break;
        case 2:
            //cout<<"test2"<<endl;
            status=GetTrackPartner(track)->GetStatus();  
            if (!((status&AliAODTrack::kTOFpid)==AliAODTrack::kTOFpid)) {return kFALSE;}
            break;
        case 3:
            //cout<<"test3"<<endl;
            status=GetTrackPartner(track)->GetStatus();  
            if ((status&AliAODTrack::kTOFmismatch)==AliAODTrack::kTOFmismatch) {return kFALSE;}
            break;
        default:
            return kFALSE;
            break;
    }

    return kTRUE;
    
}

//____________________________________________________________________________
Bool_t AliAnalysisTaskDiHadronPID::SelectEvent(AliAODVertex* vertex)

{
	
	//
	// Event Selection.
	//
	
	Double_t primVtx[3];
	vertex->GetXYZ(primVtx);
	if (TMath::Sqrt(primVtx[0]*primVtx[0] + primVtx[1]*primVtx[1])>1. || TMath::Abs(primVtx[2])>10.) {
		cout << "AliAnalysisTaskDiHadronPID::SelectEvent: Vertex Out of Range." << endl;
        cout << "AliAnalysisTaskDiHadronPID::SelectEvent: Event not selected." << endl;
		return kFALSE;
	}
	cout << "AliAnalysisTaskDiHadronPID::SelectEvent: Vertex is OK." << endl;
    
    // We also wish to make a 0-10% centrality cut.
    if (fAODHeader->GetCentrality()>10.) {
        cout << "AliAnalysisTaskDiHadronPID::SelectEvent: Non-central event." << endl;
        cout << "AliAnalysisTaskDiHadronPID::SelectEvent: Event not selected." << endl;
		return kFALSE;
    }
	cout << "AliAnalysisTaskDiHadronPID::SelectEvent: Central Event." << endl;
    
    cout << "AliAnalysisTaskDiHadronPID::SelectEvent: Event selected." << endl;
	return kTRUE;
	
}

//_____________________________________________________________________________
void AliAnalysisTaskDiHadronPID::FillPIDPartnersArray() {

	// Initialize the mapping for corresponding PID tracks. (see
	// GetTrackPartner(AliAODTrack* track)). 
	//
	
	if (!fAODEvent) {
		cout << "ERROR in CreatePIDPartersArray(): fAODEvent not set." << endl;
		return;
	}
	
	if (!fMask) {
		cout << "ERROR in CreatePIDPartersArray(): fMask not set." << endl;
		return;
	}
	
	fPIDPartners = new TObjArray();
	AliAODTrack* track = 0x0;
		
	for (Int_t iTrack = 0; iTrack < fAODEvent->GetNumberOfTracks(); iTrack++) {
		
		track = fAODEvent->GetTrack(iTrack);
		
		// cout << "Track: " << iTrack << " FilterMaskPass: "<< track->TestFilterMask(fMask) << " Track ID: " << track->GetID() << endl;
		
		// I.e., if it does NOT pass the filtermask.
		if (!(track->TestFilterMask(fMask))) {
            fPIDPartners->AddAtAndExpand(track,track->GetID());
            if (track->GetID()<1) cout<<"Track ID: "<<track->GetID()<<" Partner ID: "<<(-track->GetID()-1)<<endl;
		}
        
	}
	
}

//_____________________________________________________________________________
AliAODTrack* AliAnalysisTaskDiHadronPID::GetTrackPartner(AliAODTrack* track) {
	
	//
	// Returns the "parner track" of track, which contains the pid information.
	//
	
	AliAODTrack* partner = 0x0;
	    
    partner = (AliAODTrack*)(fPIDPartners->At(-track->GetID()-1));
	    
	if (!partner&&fVerbose) cout<<"GetTrackPartner: No Partner found!"<<endl;
	
	return partner;
	
}

//_____________________________________________________________________________
Double_t AliAnalysisTaskDiHadronPID::PhiRange(Double_t DPhi)

{
	//
	// Puts the argument in the range [-pi/2,3 pi/2].
	//
	
	if (DPhi < -TMath::Pi()/2) DPhi += 2*TMath::Pi();
	if (DPhi > 3*TMath::Pi()/2) DPhi -= 2*TMath::Pi();	

	return DPhi;
	
}

//_____________________________________________________________________________
void AliAnalysisTaskDiHadronPID::UserExec(Option_t *)

{
	//
	// UserExec.
	//
	
	// Input the event.
	fAODEvent = dynamic_cast<AliAODEvent*>(InputEvent());
	
	if (!fAODEvent) {
		cout << "ERROR: No AliAODEvent pointer could be created." << endl;
		return;
	}
	
	// Get the event header.
	fAODHeader = fAODEvent->GetHeader();
	
	if (!fAODHeader) {
		cout << "ERROR: No AliAODHeader pointer could be created."<<endl;
		return;
	}
	
	// Get the event vertex.
	fAODVertex = fAODEvent->GetPrimaryVertex();
	
	if (!fAODVertex) {
		cout << "ERROR: No AliAODVertex pointer could be created." << endl;
		return;
	}
		
	// Display basic event information.
	cout << endl;
	cout << "Event centrality: " << fAODHeader->GetCentrality() << endl;
	cout << "Event Vertex Z: " << fAODVertex->GetZ() << endl;
	cout << "Event tracks in AOD: " << fAODEvent->GetNumberOfTracks() << endl;
	cout << endl;

    if (!SelectEvent(fAODVertex)) return;
	// Fill the TObjArray which holds PID partners.
	FillPIDPartnersArray();
	
	// Filling Centrality/ VertexZ Hisogram.
	fCentrality->Fill(fAODHeader->GetCentrality());
	fVertexZ->Fill(fAODVertex->GetZ());
	
	TObjArray *triggers		= new TObjArray();
	TObjArray *associateds	= new TObjArray();
	
	// 1. Identifying sets of triggers and associateds. OK!
	// 2. Filling Pt distribution (UnID). OK!
	// 3. Making the nSigma plots for TPC and TOF.
	    
	for (Int_t iTrack = 0; iTrack < fAODEvent->GetNumberOfTracks(); iTrack++) {
        
		fAODTrack = fAODEvent->GetTrack(iTrack);

        // Skip track if there is no pointer, or if it doesn't pass the strandard track cuts.
		if ((!fAODTrack)||(!SelectTrack(fAODTrack,0))) continue;
        
        if (fAODTrack->GetID()>-1) cout<<"Track ID: "<<fAODTrack->GetID()<<endl;
        
        // Make the pT spectrum with only standard track cuts.
        fPtSpectrum->Fill(fAODTrack->Pt());  
        
		if (fAODTrack->Pt() > 5.0) {
			if (fVerbose) cout << "AliAnalysisTaskDiHadronPID::UserExec: Trigger found!" << endl;
			triggers->AddLast(fAODTrack);
		}
		
		if (fAODTrack->Pt() < 5.0) {
        
            fTrackCuts->AddBinContent(1);
        
            // Now we demand a TPC hit.            
            if (!SelectTrack(fAODTrack,1)) continue;
            //cout<<"Passed Track cuts 1"<<endl;
            fTrackCuts->AddBinContent(2);
                        
			// the associateds array contains tracks pT < 5.0 GeV/c && TPC hit.
			associateds->AddLast(fAODTrack);
            if (fAODTrack->GetID()>-1) cout<<"Assoc. Track ID: "<<fAODTrack->GetID()<<endl;

			// Make the nSigma plots.
			Double_t mom, nSigma;

			mom = GetTrackPartner(fAODTrack)->GetTPCmomentum();
			nSigma = fPIDResponse->NumberOfSigmasTPC(GetTrackPartner(fAODTrack),AliPID::kProton);
			fTPCnSigmaProton->Fill(mom,nSigma);
			nSigma = fPIDResponse->NumberOfSigmasTPC(GetTrackPartner(fAODTrack),AliPID::kPion);
			fTPCnSigmaPion->Fill(mom,nSigma);
            nSigma = fPIDResponse->NumberOfSigmasTPC(GetTrackPartner(fAODTrack),AliPID::kKaon);
			fTPCnSigmaKaon->Fill(mom,nSigma);

			fTPCSignal->Fill(fAODTrack->Eta(),fAODTrack->Pt(),GetTrackPartner(fAODTrack)->GetTPCsignal());
			
            if (!SelectTrack(fAODTrack,2)) continue;
            fTrackCuts->AddBinContent(3);
            if (!SelectTrack(fAODTrack,3)) continue;
            fTrackCuts->AddBinContent(4);
            
			mom = GetTrackPartner(fAODTrack)->P();
			nSigma = fPIDResponse->NumberOfSigmasTOF(GetTrackPartner(fAODTrack),AliPID::kProton);
			fTOFnSigmaProton->Fill(mom,nSigma);			
            nSigma = fPIDResponse->NumberOfSigmasTOF(GetTrackPartner(fAODTrack),AliPID::kPion);
			fTOFnSigmaPion->Fill(mom,nSigma);	
            nSigma = fPIDResponse->NumberOfSigmasTOF(GetTrackPartner(fAODTrack),AliPID::kKaon);
			fTOFnSigmaKaon->Fill(mom,nSigma);
            
			fTOFSignal->Fill(fAODTrack->Eta(),fAODTrack->Pt(),GetTrackPartner(fAODTrack)->GetTOFsignal());

		}
				
	}
	
    //cout<<"Done with the first loop."<<endl;
    
	// 1. Making the di-hadron correlations. (to do: set this up a bit nicer!)
    // {Dphi, Deta, TPC signal, TOF signal}
	Double_t histoFill[4];
	AliAODTrack* currentTrigger = 0x0;
	AliAODTrack* currentAssociated = 0x0;
	AliAODTrack* currentPIDPartner = 0x0;
	//AliAODPid* currentPIDObject = 0x0;
	
	AliTPCPIDResponse& TPCPIDResponse = fPIDResponse->GetTPCResponse();
	//AliTOFPIDResponse& TOFPIDResponse = fPIDResponse->GetTOFResponse();
        
	for (Int_t iTrig = 0; iTrig < triggers->GetEntriesFast(); iTrig++){
	
		currentTrigger = (AliAODTrack*)(triggers->At(iTrig));
		
		for (Int_t iAssoc = 0; iAssoc < associateds->GetEntriesFast(); iAssoc++) {
		
			currentAssociated = (AliAODTrack*)(associateds->At(iAssoc));
						
			histoFill[0] = PhiRange(currentTrigger->Phi() - currentAssociated->Phi());
			histoFill[1] = currentTrigger->Eta() - currentAssociated->Eta();
			Double_t pt = currentAssociated->Pt();

			// Be aware that there may be a caveat here when Pt = 5.00000000
			const Int_t ptbin = (Int_t)(2*currentAssociated->Pt());
			
            fDiHadron->Fill(histoFill[0],histoFill[1],pt);
			
			currentPIDPartner = GetTrackPartner(currentAssociated);
            //currentPIDObject = currentPIDPartner->GetDetPid();
			
			if (currentPIDPartner/*&&currentPIDObject*/) {
				
				Double_t TPCmom = currentPIDPartner->GetTPCmomentum();
                Double_t TPCsignal = currentPIDPartner->GetTPCsignal();
                Double_t TOFsignal = -999.;
                Double_t expectedTOFsignalKaon=0,expectedTOFsignalPion=0,expectedTOFsignalProton=0;
                Double_t times[AliPID::kSPECIES]; // For the expected time per particle species. 
            
				if (SelectTrack(currentAssociated,2)) {
                    TOFsignal = currentPIDPartner->GetTOFsignal();
                    
                    //currentPIDObject->GetIntegratedTimes(times);
                    currentPIDPartner->GetIntegratedTimes(times);
                    
                    expectedTOFsignalPion = times[AliPID::kPion];
                    expectedTOFsignalKaon = times[AliPID::kKaon];
                    expectedTOFsignalProton = times[AliPID::kProton]; 
                }
                
				Double_t expectedTPCsignalPion = TPCPIDResponse.GetExpectedSignal(TPCmom,AliPID::kPion);
                Double_t expectedTPCsignalKaon = TPCPIDResponse.GetExpectedSignal(TPCmom,AliPID::kKaon);
				Double_t expectedTPCsignalProton = TPCPIDResponse.GetExpectedSignal(TPCmom,AliPID::kProton);
				
                histoFill[2] = TPCsignal - expectedTPCsignalPion;
                fDiHadronTPC[0][ptbin]->Fill(histoFill[0],histoFill[1],histoFill[2]);                
                if (SelectTrack(currentAssociated,2)&&SelectTrack(currentAssociated,3)) {
                    histoFill[3] = TOFsignal - expectedTOFsignalPion;
                    fDiHadronTOF[0][ptbin]->Fill(histoFill[0],histoFill[1],histoFill[3]);
                    fDiHadronTPCTOF[0][ptbin]->Fill(histoFill);
                }

                histoFill[2] = TPCsignal - expectedTPCsignalKaon;
                fDiHadronTPC[1][ptbin]->Fill(histoFill[0],histoFill[1],histoFill[2]);                
                if (SelectTrack(currentAssociated,2)&&SelectTrack(currentAssociated,3)) {
                    histoFill[3] = TOFsignal - expectedTOFsignalKaon;
                    fDiHadronTOF[1][ptbin]->Fill(histoFill[0],histoFill[1],histoFill[3]);
                    fDiHadronTPCTOF[1][ptbin]->Fill(histoFill);
                }
                
                histoFill[2] = TPCsignal - expectedTPCsignalProton;
                fDiHadronTPC[2][ptbin]->Fill(histoFill[0],histoFill[1],histoFill[2]);                
                if (SelectTrack(currentAssociated,2)&&SelectTrack(currentAssociated,3)) {
                    histoFill[3] = TOFsignal - expectedTOFsignalProton;
                    fDiHadronTOF[2][ptbin]->Fill(histoFill[0],histoFill[1],histoFill[3]);
                    fDiHadronTPCTOF[2][ptbin]->Fill(histoFill);
                }
                                
                fAssociatedDistribution->Fill(currentAssociated->Phi(),currentAssociated->Eta(),currentAssociated->Pt());
			}
		}
	}
	
    if (fCalculateMixedEvents) {
        
        
        // Loop over the trigger buffer.
        if (fVerbose) cout << "AliAnalysisTaskDiHadronPID::UserExec: Mixing the events with "<<fTrigBufferSize<<" triggers from the buffer." <<endl;
        if (fVerbose) cout << "AliAnalysisTaskDiHadronPID::UserExec: Buffer size: "<<fTrigBufferIndex<<endl;
        
        for (Int_t iTrig=0;iTrig<fTrigBufferSize;iTrig++) {
            
            // Check if the trigger and the associated have a reconstructed
            // vertext no further than 2cm apart.
            
            // fTrigBuffer[i][0] = z
            // fTrigBuffer[i][1] = phi
            // fTrigBuffer[i][2] = eta
            // fTrigBuffer[i][3] = p_t
            
            if (TMath::Abs(fTrigBuffer[iTrig][0]-fAODVertex->GetZ())<2.) {
                
                cout<<"AliAnalysisTaskDiHadronPID::UserExec: Mixing with trigger Z: "<<fTrigBuffer[iTrig][0]<<", Pt: "<<fTrigBuffer[iTrig][3]<<endl;
                
                for (Int_t iAssoc = 0; iAssoc < associateds->GetEntriesFast(); iAssoc++) {
                    
                    currentAssociated = (AliAODTrack*)(associateds->At(iAssoc));
                    currentPIDPartner = GetTrackPartner(currentAssociated);
                    //currentPIDObject = currentPIDPartner->GetDetPid();
                    
                    if (currentPIDPartner/*&&currentPIDObject*/) {
                        
                        Double_t DPhi = PhiRange(fTrigBuffer[iTrig][1] - currentAssociated->Phi());
                        Double_t DEta = fTrigBuffer[iTrig][2] - currentAssociated->Eta();
                        Double_t Ptassoc = currentAssociated->Pt();
                        
                        fMixedEvents->Fill(DPhi,DEta,Ptassoc);
                    }
                }
            }
        }
    
        // Copy the triggers from the current event into the buffer.
        if (fAODVertex->GetZ()<10.) {
        
            if (fVerbose) cout<<"AliAnalysisTaskDiHadronPID::UserExec: Copying "<<triggers->GetEntriesFast()<<" triggers to the buffer with vertex z = "<<fAODVertex->GetZ()<<endl;
        
            for (Int_t iTrig = 0; iTrig<triggers->GetEntriesFast(); iTrig++) {
            
                currentTrigger = (AliAODTrack*)(triggers->At(iTrig));
                if (fVerbose) cout<<"AliAnalysisTaskDiHadronPID::UserExec: Trigger pt = "<<currentTrigger->Pt()<<endl;
            
                fTrigBuffer[fTrigBufferIndex][0] = fAODVertex->GetZ();
                fTrigBuffer[fTrigBufferIndex][1] = currentTrigger->Phi();
                fTrigBuffer[fTrigBufferIndex][2] = currentTrigger->Eta();
                fTrigBuffer[fTrigBufferIndex][3] = currentTrigger->Pt();
                fTrigBufferIndex++;
                if (fTrigBufferSize<25000) {fTrigBufferSize++;}
                if (fTrigBufferIndex==25000) {fTrigBufferIndex=0;}
            }
        }
    }        
        
    if (fVerbose) cout<<"AliAnalysisTaskDiHadronPID::UserExec: Trigger buffer index: "<<fTrigBufferIndex<<", and size: "<<fTrigBufferSize<<endl;
    
	delete triggers;
	delete associateds;
	 
	PostData(1,fHistoList);
	
}

//_____________________________________________________________________________
void AliAnalysisTaskDiHadronPID::Terminate(Option_t *)

{
	//
	// Terminate.
    //
    
}

