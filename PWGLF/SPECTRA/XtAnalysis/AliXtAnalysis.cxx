/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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

/**************************************************************************
 *     AliXtAnalysis:
 * This class constructs inclusive and isolated (based on charged tracks)
 * invariant spectra. Isolated case uses a different efficiency where the
 * contamination from cases where non-isolated particle appears as an
 * isolated one due to finite detector efficiency is taken into account.
 *
 * contact: Sami Räsänen
 *          University of Jyväskylä, Finland 
 *          sami.s.rasanen@jyu.fi
**************************************************************************/

// general + root classes
#include <TList.h>
#include <TChain.h>
#include <TObjArray.h>

// AliRoot classes
#include <AliAnalysisManager.h>
#include <AliInputEventHandler.h>
#include <AliAnalysisUtils.h>
#include <AliAODTrack.h>
#include <AliAODEvent.h>

// Jyväskylä classes
#include <AliJCard.h>  // this class at PWGCF / JCORRAN
#include "AliJXtHistos.h"
//#include "AliIsolatedEfficiency.h"

// This analysis
#include "AliXtAnalysis.h"
 
ClassImp(AliXtAnalysis)

//________________________________________________________________________
AliXtAnalysis::AliXtAnalysis() 
      : AliAnalysisTaskSE(),
        fOutputList(0x0),
        fAnaUtils(0x0),
        fCard(0x0),
        fHistDir(0x0),
        fHistos(0x0),
        fhEvents(0x0),
	//fEfficiency(0x0),
	fChargedList(0x0),
	fIsolatedChargedList(0x0),
	fCentBin(-1),
	fZBin(-1),
	fZVert(999999.),
	fevt(0),
	fDebugMode(0)
{
// Constructor
}
//________________________________________________________________________
AliXtAnalysis::AliXtAnalysis(const char *name, const char * cardname) 
      :	AliAnalysisTaskSE(name),
	fOutputList(0x0), 
	fAnaUtils(0x0),
   	fCard(0x0),
	fHistDir(0x0),
	fHistos(0x0),
   	fhEvents(0x0),
	//fEfficiency(0x0),
	fChargedList(0x0),
	fIsolatedChargedList(0x0),
	fCentBin(-1),
	fZBin(-1),
	fZVert(999999.),
	fevt(0),
	fDebugMode(0)
{
	// All parameters of the analysis
	fCard = new AliJCard(cardname);
	
	cout << "debug: card created to address " << fCard << endl;
	
	fCard->PrintOut();
	
	fAnaUtils = new AliAnalysisUtils();
	fAnaUtils->SetUseOutOfBunchPileUp( kTRUE );
	fAnaUtils->SetUseSPDCutInMultBins( kTRUE);
	
	// Define input and output slots here
	// Input slot #0 works with a TChain
	DefineInput(0, TChain::Class());
	// Output slot #0 writes into a TH1 container
	DefineOutput(1, TList::Class());
	// JHistos into TDirectory
	DefineOutput(2, TDirectory::Class());
}
//________________________________________________________________________
AliXtAnalysis::AliXtAnalysis(const AliXtAnalysis& a)
      :	AliAnalysisTaskSE(a.GetName()),
	fOutputList(a.fOutputList),
	fAnaUtils(a.fAnaUtils),
	fCard(a.fCard),
	fHistDir(a.fHistDir),
	fHistos(a.fHistos),
	fhEvents(a.fhEvents),
	//fEfficiency(a.fEfficiency),
	fChargedList(a.fChargedList),
	fIsolatedChargedList(a.fIsolatedChargedList),
	fCentBin(-1),
	fZBin(-1),
	fZVert(999999.),
	fevt(0),
	fDebugMode(0)
{
    //copy constructor
    fhEvents = (TH1D*) a.fhEvents->Clone(a.fhEvents->GetName());
}
//________________________________________________________________________
AliXtAnalysis& AliXtAnalysis::operator = (const AliXtAnalysis& ap){
	// assignment operator
	this->~AliXtAnalysis();
	new(this) AliXtAnalysis(ap);
	return *this;
}
//________________________________________________________________________
AliXtAnalysis::~AliXtAnalysis() {
    // destructor
    delete fOutputList;
    delete [] fAnaUtils;
    delete [] fhEvents;
    delete fHistos;
    //delete fEfficiency;
    delete fCard;
    delete fChargedList;
    delete fIsolatedChargedList;
}
//________________________________________________________________________
void AliXtAnalysis::UserCreateOutputObjects(){
    // Create histograms
    // Called once
    
    cout<<"\n=============== CARD =============="<<endl;
    fCard->PrintOut();
    cout<<"===================================\n"<<endl;
    
    fhEvents = new TH1D("hEvents","events passing cuts", 11, -0.5, 10.5 );
    fhEvents->SetXTitle( "0 - all, 1 - pileup, 2 - kMB selected, 3 - has vertex, 4 - good vertex, 5 - MB + good vertex, 6 - MB + good vertex + centrality" );
    fhEvents->Sumw2();
    
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    //for any MC track filled with MC pt
    fOutputList->Add(fhEvents);

    fChargedList = new TObjArray;
    fChargedList->SetOwner(kTRUE);
    fIsolatedChargedList = new TObjArray;
    
    PostData(1, fOutputList);
    
    bool orignalTH1AdddirectoryStatus=TH1::AddDirectoryStatus();
    TH1::AddDirectory(kTRUE);
    if( !orignalTH1AdddirectoryStatus ) cout<<"DEBUG : TH1::AddDirectory is turned on"<<endl;
    //TFile * file2 = OpenFile(2);
    
    fHistDir=gDirectory;
    fHistos = new AliJXtHistos(fCard);
    fHistos->CreateXtHistos();
    
    PostData(2, fHistDir);
    TH1::AddDirectory( orignalTH1AdddirectoryStatus );
    cout<<"DEBUG : TH1::AddDirectory get orignal Value = "<<( orignalTH1AdddirectoryStatus?"True":"False" )<<endl;
    
    //fEfficiency = new AliJEfficiency();
    
    fevt = 0;
}
//________________________________________________________________________
void AliXtAnalysis::UserExec(Option_t *) {
    
    // Main loop - called for each event
    
    fevt++;
    if(fevt % 10000 == 0) cout << "Number of event scanned = "<< fevt << endl;
    
    Int_t filterBit = AliAODTrack::kTrkGlobal; //fCard->Get("TrackFilterBit");
    
    AliVEvent *event = InputEvent();
    if(!event) return;
    
    // check if the event was triggered or not and vertex
    if( !IsGoodEvent( event ) ) return; // fZBin is set there
    fhEvents->Fill( 5 );
    if(fDebugMode > 0) cout << "zvtx = " << fZVert << endl;
    
    AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(event);
    if(!aodEvent) return;
    
    // centrality - in case of pPb and PbPb, dictated by DataType in the card
    double fcent;
    if( fCard->Get("DataType") > 0.5 ){
        AliCentrality *cent = event->GetCentrality();
        if( ! cent ) return;
        fcent = cent->GetCentralityPercentile("V0M");
        fCentBin = fCard->GetBin(kCentrType, fcent);;
        //cout <<"Centrality = "<< fcent <<"\t"<< fCentBin << endl;
    }else{
        fcent = 0.0;
        fCentBin = 0;
    }
    
    if(fCentBin<0) return;
    fhEvents->Fill( 6 );
    Int_t nt = aodEvent->GetNumberOfTracks();
    
    fHistos->FillCentralityHistos(fcent, fCentBin);
    
    // clear them up for every event
    fChargedList->Clear();
    fIsolatedChargedList->Clear();
    
    // Is there a better way than hard code this into card???
    double sqrts = fCard->Get("Sqrts");
    
    // Isolation method: 0 = relative, 1 = absolute
    double isolMethod = fCard->Get("IsolationMethod");
    
    // Minimum isolation pT that is considered
    double minIsolationPt = fCard->Get("MinimumIsolatedPt");
    
    // Dummy numbers; either of these is updated based on selected method (other is void)
    double isolationFraction  = 99999.;
    double isolationThreshold = 99999.;
    
    if( isolMethod > 0.5 ){
        isolationThreshold = fCard->Get("IsolationThreshold");  // this is fixed threshold
    }else{
        isolationFraction = fCard->Get("IsolationFraction");  // fraction in pT, threshold updated later in a loop
    }
    
    // Isolation radius
    double isolR = fCard->Get("IsolationRadius");
    
    // Get good tracks to the list
    for(Int_t it = 0; it < nt; it++) {
        AliAODTrack *track = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(it));
        if( !track ) continue;
        if( !track->TestFilterBit(filterBit) ) continue;
        double eta = track->Eta();
        if( !fCard->IsInEtaRange( eta ) ) continue; // Do not use fiducial cut here
        AliAODTrack *acceptedTrack = new AliAODTrack( *track );
        double pT = acceptedTrack->Pt();
        double xT = 2.*pT/sqrts;
        double phi = acceptedTrack->Phi();
        
        // TODO:
        double effCorr = 1.0;
        
        // Fill histos
        fHistos->FillInclusiveHistograms(pT, xT, eta, phi, effCorr, fCentBin);
        
        // Add an accepted track into list of inclusive charged particles
        fChargedList->Add( acceptedTrack );
    }
    
    // Check isolation of the found good tracks
    Int_t nAcc = fChargedList->GetEntriesFast();
    
    // Get eta range for fiducial cut
    double etaRange = fCard->Get("EtaRange");
    
    for( Int_t it = 0; it < nAcc; it++ ){
        AliAODTrack *track = (AliAODTrack*)fChargedList->At(it);
        // Fiducial eta cut for isolated particles
        double eta = track->Eta();
        if( abs(eta) > etaRange - isolR ) continue;
        
        // To check the isolation with cone, use TLorentzVector -routines
        TLorentzVector lvTrack = TLorentzVector( track->Px(), track->Py(), track->Pz(), track->P() );
        
        // Set isolation threshold, if relative isolation
        double pT = lvTrack.Pt();
        if( pT < minIsolationPt ) continue;
        if( isolMethod < 0.5 ) isolationThreshold = pT * isolationFraction; // here: relative isolation
        
        // Isolation check
        double sum = 0.0;
        for( Int_t ia = 0; ia < nAcc; ia++ ){
            if( ia == it ) continue; // Do not count the particle itself
            AliAODTrack *assocTrack = (AliAODTrack*)fChargedList->At(ia);
            TLorentzVector lvAssoc = TLorentzVector( assocTrack->Px(), assocTrack->Py(), assocTrack->Pz(), assocTrack->P() );
            if( lvTrack.DeltaR( lvAssoc ) < isolR ) sum += lvAssoc.Pt();
        }
        
        // Cone activity
        fHistos->FillInclusiveConeActivities(pT,sum);  // efficiency correction?
        
        // If the particle was isolated, fill histograms and the list
        if( sum < isolationThreshold ){
            
            double xT = 2.*pT/sqrts;
            
            if( fDebugMode > 0 ){
                cout << "DEBUG: isolated particle found " << endl;
                cout << "fCentBin = " << fCentBin << ", pT = " << pT << ", xT = " << xT << " and eta = " << eta << endl;
            }
            
            fHistos->FillIsolatedConeActivities(pT, sum);  // efficiency correction?
            
            // TODO:
            double effCorr = 1.0;

            // Fill histos
            fHistos->FillIsolatedHistograms(pT, xT, eta, track->Phi(), effCorr, fCentBin);
            
            // Add tracks to the isolated list
            AliAODTrack *acceptedIsolatedTrack = new AliAODTrack( *track );
            fIsolatedChargedList->Add( acceptedIsolatedTrack );
        }
    }
    
    // For further analysis, there is now lists of all charged and isolated charged available, will be usefull.
    
    return;  // END
}
//________________________________________________________________________
void AliXtAnalysis::Terminate(Option_t *)
{
    cout<<"Successfully finished"<<endl;
}
//________________________________________________________________________
bool AliXtAnalysis::IsGoodEvent(AliVEvent *event) {
    
    // This function checks that the event has the requested trigger
    // and that the z-vertex is in given range
    
    
    fhEvents->Fill( 0 );
    
    if(fAnaUtils->IsPileUpEvent(event)) {
        return kFALSE;
    } else {
        Bool_t triggeredEventMB = kFALSE; //init
        
        fhEvents->Fill( 1 );
        
        Bool_t triggerkMB = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & ( AliVEvent::kMB );
        
        if( triggerkMB ){
            triggeredEventMB = kTRUE;  //event triggered as minimum bias
            fhEvents->Fill( 2 );
        }
        //--------------------------------------------------------------
        // check reconstructed vertex
        int ncontributors = 0;
        Bool_t goodRecVertex = kFALSE;
        const AliVVertex *vtx = event->GetPrimaryVertex();
        if(vtx){
            fhEvents->Fill( 3 );
            fZVert = vtx->GetZ();
            fHistos->FillRawVertexHisto(fZVert);
            ncontributors = vtx->GetNContributors();
            if(ncontributors > 0){
                if(fCard->VertInZRange(fZVert)) {
                    goodRecVertex = kTRUE;
                    fhEvents->Fill( 4 );
                    fHistos->FillAcceptedVertexHisto(fZVert, fCentBin);
                    fZBin  = fCard->GetBin(kZVertType, fZVert);
                }
            }
        }
        return triggerkMB && goodRecVertex;
    }
    //---------------------------------
}
