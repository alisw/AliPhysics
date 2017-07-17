#include <TChain.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TVectorT.h>


#include "AliJBaseTrack.h"
#include "AliJHSInterplayTask.h"
#include "AliJCorrelations.h"
#include "AliJEventPool.h"
#include "AliJEbePercentile.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliESDtrackCuts.h"
#include "AliVVertex.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliGenHijingEventHeader.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h" 
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliESDVertex.h"
#include "AliVParticle.h"
#include "AliCentrality.h" 
#include "AliEventplane.h"

#include "AliJHistos.h"
#include "AliJEbeHistos.h"
#include "AliJJet.h"

#include "AliMultSelection.h"
//#include "AliJRunTable.h"


ClassImp(AliJHSInterplayTask)

AliJHSInterplayTask::AliJHSInterplayTask() 
	: AliAnalysisTaskSE(), fOutput(0x0),
	fAnaUtils(NULL),
	fFFlucAna(NULL),
	fJetTask(NULL),
	fJetTaskName(""),
	fJetSel(0),
	fCard(NULL),
	fHistos(NULL),
	fEfficiency(NULL),
	fHadronSelectionCut(0),
	fInputList(NULL),
	fInputListSpectra(NULL),
	fVnMethod(0),
	fESMethod(0),
	fInputListFlow(NULL),
	fFirstEvent(kTRUE),
	cBin(-1),
	zBin(-1),
	zVert(-999),
	fevt(0),
	fDebugMode(0),
	trkfilterBit(0),
	fRunTable(0),
	fRandom(NULL),
	IsMC(kFALSE),
	IsKinematicOnly(kFALSE)
{
	pfOutlierLowCut = new TF1("fLowCut","[0]+[1]*x - 5.*([2]+[3]*x+[4]*x*x+[5]*x*x*x)",0,100);
	pfOutlierHighCut = new TF1("fHighCut","[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)",0,100);
	
	pfOutlierLowCut->SetParameters(0.0157497, 0.973488, 0.673612, 0.0290718, -0.000546728, 5.82749e-06);
	pfOutlierHighCut->SetParameters(0.0157497, 0.973488, 0.673612, 0.0290718, -0.000546728, 5.82749e-06);
}
//________________________________________________________________________
	AliJHSInterplayTask::AliJHSInterplayTask(const char *name) 
: AliAnalysisTaskSE(name), 
	fOutput(0x0), 
	fAnaUtils(0x0),
	fFFlucAna(0x0),
	fJetTask(0x0),
	fJetTaskName(""),
	fJetSel(0),
	fCard(0x0),
	fHistos(0x0),
	fEfficiency(0x0),
	fHadronSelectionCut(0),
	fInputList(0x0),
	fInputListSpectra(0x0),
	fVnMethod(0),
	fESMethod(0),
	fInputListFlow(0x0),
	fFirstEvent(kTRUE),
	cBin(-1),
	zBin(-1),
	zVert(-999),
	fevt(0),
	fDebugMode(0),
	trkfilterBit(0),
	fRunTable(0),
	fRandom(0x0),
	IsMC(kFALSE),
	IsKinematicOnly(kFALSE),
	fPtHardMin(0),
	fPtHardMax(0),
	TagThisEvent(kFALSE)
{

	pfOutlierLowCut = new TF1("fLowCut","[0]+[1]*x - 5.*([2]+[3]*x+[4]*x*x+[5]*x*x*x)",0,100);
	pfOutlierHighCut = new TF1("fHighCut","[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)",0,100);
	
	pfOutlierLowCut->SetParameters(0.0157497, 0.973488, 0.673612, 0.0290718, -0.000546728, 5.82749e-06);
	pfOutlierHighCut->SetParameters(0.0157497, 0.973488, 0.673612, 0.0290718, -0.000546728, 5.82749e-06);


	// Define input and output slots here
	// Input slot #0 works with a TChain
	DefineInput(0, TChain::Class());
	// Output slot #0 writes into a TH1 container
	DefineOutput(1, TDirectory::Class());
}

//________________________________________________________________________
AliJHSInterplayTask::AliJHSInterplayTask(const AliJHSInterplayTask& a):
	AliAnalysisTaskSE(a.GetName()),
	fOutput(a.fOutput),
	fAnaUtils(a.fAnaUtils),
	fFFlucAna(a.fFFlucAna),
	fJetTask(a.fJetTask),
	fJetTaskName(a.fJetTaskName),
	fCard(a.fCard),
	fHistos(a.fHistos),
	fEfficiency(a.fEfficiency),
	fHadronSelectionCut(a.fHadronSelectionCut),
	fInputList(a.fInputList),
	fInputListSpectra(a.fInputListSpectra),
	fVnMethod(a.fVnMethod),
	fESMethod(a.fVnMethod),
	fInputListFlow(a.fInputListFlow),
	fFirstEvent(kTRUE),
	cBin(-1),
	zBin(-1),
	zVert(-999),
	fevt(-1),
	fDebugMode(0), 
	trkfilterBit(a.trkfilterBit),
	fRunTable(a.fRunTable),
	fRandom(a.fRandom),
	IsMC(a.IsMC),
	IsKinematicOnly(a.IsKinematicOnly),
	fPtHardMin(a.fPtHardMin),
	fPtHardMax(a.fPtHardMax),
	TagThisEvent(a.TagThisEvent)
{
	pfOutlierLowCut = (TF1*)a.pfOutlierLowCut->Clone();
	pfOutlierHighCut = (TF1*)a.pfOutlierHighCut->Clone();
}
//________________________________________________________________________
AliJHSInterplayTask& AliJHSInterplayTask::operator = (const AliJHSInterplayTask& ap){
	// assignment operator

	JUNUSED(ap);
	this->~AliJHSInterplayTask();
	new(this) AliJHSInterplayTask(ap);
	return *this;
}

Bool_t AliJHSInterplayTask::UserNotify() {

	cout <<"UserNotify  ievt =" << fevt << endl;
	return kTRUE;
}


//________________________________________________________________________
void AliJHSInterplayTask::UserCreateOutputObjects(){
	// Create histograms
	// Called once

	fRandom = new TRandom(0);

	fAnaUtils = new AliAnalysisUtils();
	fAnaUtils->SetUseOutOfBunchPileUp( kTRUE );
	fAnaUtils->SetUseSPDCutInMultBins( kTRUE);


	cout<<"\n=============== CARD =============="<<endl;
	fCard->PrintOut();
	cout<<"===================================\n"<<endl;

	fInputList  = new TClonesArray("AliJBaseTrack",1500);
	fInputListSpectra  = new TClonesArray("AliJBaseTrack",1500);
	fInputListFlow  = new TClonesArray("AliJBaseTrack",1500);
	fInputList->SetOwner(kTRUE);
	fInputListSpectra->SetOwner(kTRUE);
	fInputListFlow->SetOwner(kTRUE);

	OpenFile(1);
	fOutput = gDirectory;
	fOutput->cd();


	fHistos = new AliJHistos(fCard);
	fHistos->CreateEventTrackHistos();
	fHistos->CreateJetHistos();
	fHistos->fHMG->Print();

	fEfficiency = new AliJEfficiency();	
	fEfficiency->SetMode( fCard->Get("EfficiencyMode") ); // 0:NoEff, 1:Period 2:RunNum 3:Auto
	if(IsMC || IsKinematicOnly) fEfficiency->SetMode( 0 );
	fEfficiency->SetDataPath("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data"); // Efficiency root file location local or alien

	fHadronSelectionCut = int ( fCard->Get("HadronSelectionCut"));
	trkfilterBit	= Int_t ( fCard->Get("AODTrackFilterBit"));
	//fVnMethod 	= int (fCard->Get("VnMethod"));
	fESMethod 	= int (fCard->Get("ESMethod"));
	fJetSel		= int (fCard->Get("JetSel"));

	//=== Get AnalysisManager
	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
	// Add JJet
	fJetTask = (AliJJetTask*)(man->GetTask( fJetTaskName));
	// Add a AliJFlucAnalysis
	fFFlucAna = new AliJFFlucAnalysis("JFFlucAnalysis");
	fFFlucAna->SetDebugLevel(1);
	fFFlucAna->SetIsPhiModule( kFALSE );
	fFFlucAna->SetIsSCptdep( kTRUE );
	fFFlucAna->SetSCwithQC( kTRUE );
	fFFlucAna->SetEbEWeight( kTRUE ); // << this is important
	fOutput->cd();
	fFFlucAna->UserCreateOutputObjects();

	fCard->WriteCard(fOutput);
	fHistos->fHMG->WriteConfig();
	fFirstEvent = kTRUE;
	fevt = -1;
	PostData(1, fOutput);
}
//________________________________________________________________________
AliJHSInterplayTask::~AliJHSInterplayTask() {
	delete pfOutlierLowCut;
	delete pfOutlierHighCut;
	delete fOutput; 
	delete fAnaUtils;
	delete fFFlucAna;
	delete fHistos;
	delete fEfficiency;
	delete fInputList;
	delete fInputListSpectra;
	delete fCard;
	delete fRandom;
}

//________________________________________________________________________
void AliJHSInterplayTask::UserExec(Option_t *) {
	// Main loop
	// Called for each event
	fevt++;
	if(fevt % 1000 == 0) cout << "Number of events scanned = "<< fevt << endl;

	// clear them up for every event
	fInputList->Clear();
	fInputListSpectra->Clear();
	fInputListFlow->Clear();

	float fcent = -999;
	float fImpactParameter = -999;
	double fvertex[3] = {-999,-999,-999};
	AliAODEvent* aodEvent;
	AliVEvent *event = NULL;
	AliMCEvent *mcEvent = NULL;

	// load current event and save track, event info
	if(IsKinematicOnly) {
		mcEvent = MCEvent();
		if (!mcEvent) {
			AliError("ERROR: mcEvent not available");
			return;
		}
		Double_t fReactionPlane = 0.;
		AliGenHijingEventHeader* headerH = dynamic_cast<AliGenHijingEventHeader*>(mcEvent->GenEventHeader());
		if (headerH) {
			fReactionPlane = headerH->ReactionPlaneAngle();
			fImpactParameter = headerH->ImpactParameter();
			fcent = GetCentralityFromImpactPar(fImpactParameter);
		}
		AliGenEventHeader *header = mcEvent->GenEventHeader();
		if(!header) return;
		TArrayF gVertexArray;
		header->PrimaryVertex(gVertexArray);
		for(int i=0;i<3;i++) fvertex[i]=gVertexArray.At(i);
		zVert = fvertex[2];
		zBin  = fCard->GetBin(kZVertType, zVert); // Real data IsGoodEvent
		ReadKineTracks( mcEvent, fInputList);
		ReadKineTracks( mcEvent, fInputListSpectra);
	} else { // Kine
		event = InputEvent();
		if(!event) return;

		aodEvent = dynamic_cast<AliAODEvent*>(event);
		if(!aodEvent) return;

		if( fFirstEvent ) {
			fRunTable = & AliJRunTable::GetSpecialInstance();
			fRunTable->SetRunNumber( aodEvent->GetRunNumber() );
			fEfficiency->SetRunNumber( aodEvent->GetRunNumber() );
			fEfficiency->Load();
			fFirstEvent = kFALSE;
		}

		if(!IsGoodEvent( event ))
			return; // zBin is set there
		if(fDebug) cout << "zvtx = " << zVert << endl;

		// centrality 
		if(fRunTable->IsHeavyIon() || fRunTable->IsPA()){
			AliMultSelection *pms = (AliMultSelection*)event->FindListObject("MultSelection");
			if(!pms){
				AliError("MultSelection unavailable.");
				return;
			}

			fcent = pms->GetMultiplicityPercentile("V0M");
			//Float_t cl0cent = pms->GetMultiplicityPercentile("CL0");
		} else {
			fcent = -1;
			cout<<"warning: centrality unavailable";
		}
	}

	//cout<<"cent = "<<fcent<<endl;
	cBin = fCard->GetBin(kCentrType, fcent);
	if(cBin<0) return;

	fHistos->fhZVert[cBin]->Fill(zVert);
	if(fDebug ) cout <<"Centrality = "<< fcent <<"\t"<< cBin << endl;

	fHistos->fhEvents->Fill( 4 );

	int counter = 0;

	// Making the inputlist from AOD
	if( IsMC == kTRUE && IsKinematicOnly == kFALSE ){  // how to get a flag to check  MC or not !
		TClonesArray *mcArray = (TClonesArray*) aodEvent->FindListObject(AliAODMCParticle::StdBranchName());
		if(!mcArray){ Printf("Error not a proper MC event"); return;};  // check mc array

		Int_t nt = mcArray->GetEntriesFast();
		Int_t ntrack =0;
		for( int it=0; it< nt ; it++){
			AliAODMCParticle *mctrack = (AliAODMCParticle*)mcArray->At(it);
			if( !mctrack ) continue;
			if( mctrack->IsPhysicalPrimary() ){
				if(!fCard->IsInEtaRange(mctrack->Eta())) continue;
				AliJBaseTrack* track = new( (*fInputList)[fInputList->GetEntriesFast()] ) AliJBaseTrack;
				track->SetPxPyPzE(mctrack->Px(), mctrack->Py(), mctrack->Pz(), mctrack->E());
				Int_t pdg = mctrack->GetPdgCode();
				Char_t ch = (Char_t) mctrack->Charge();
				Int_t label = mctrack->GetLabel();
				track->SetID(fInputList->GetEntriesFast());
				track->SetParticleType(kJHadron);
				track->SetCharge(ch);
				track->SetLabel( label );
				double ptt = mctrack->Pt();
				double effCorr = 1.;  // here you generate warning if ptt>30
				fHistos->fhTrackingEfficiency[cBin]->Fill( ptt, 1./effCorr );
				track->SetTrackEff( 1./effCorr );
			} // PhysicalPrimary
		} // mcArray
		RegisterList(fInputListSpectra, fInputList, 0.2, 200.);
	} // read mc track done.

	if( IsMC == kFALSE && IsKinematicOnly == kFALSE ){  
		Int_t nt = aodEvent->GetNumberOfTracks();
		for(Int_t it = 0; it < nt; it++) {
			AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(it));
			if( !aodtrack ) continue;
			if(!fCard->IsInEtaRange(aodtrack->Eta())) continue;
			// For Vn selection and correlation
			if(aodtrack->TestFilterBit(trkfilterBit)) { 
				AliJBaseTrack* track = new( (*fInputList)[fInputList->GetEntriesFast()] ) AliJBaseTrack;
				track->SetPxPyPzE(aodtrack->Px(), aodtrack->Py(), aodtrack->Pz(), aodtrack->E());
				track->SetID(fInputList->GetEntriesFast());
				track->SetParticleType(kJHadron);
				track->SetCharge(aodtrack->Charge());
				double ptt = track->Pt();
				double effCorr = 1./fEfficiency->GetCorrection(ptt, fHadronSelectionCut, fcent);  // here you generate warning if ptt>30
				track->SetTrackEff( 1./effCorr );
			}
			// Spectra analysis
			if(aodtrack->TestFilterBit(1024)) {  // Raa Cut for spectra bit 1024 Hardon selection for eff =1  
				AliJBaseTrack* track = new( (*fInputListSpectra)[fInputListSpectra->GetEntriesFast()] ) AliJBaseTrack;
				track->SetPxPyPzE(aodtrack->Px(), aodtrack->Py(), aodtrack->Pz(), aodtrack->E());
				track->SetID(fInputListSpectra->GetEntriesFast());
				track->SetParticleType(kJHadron);
				track->SetCharge(aodtrack->Charge());
				double ptt = track->Pt();
				double effCorr = 1./fEfficiency->GetCorrection(ptt, 1, fcent);  // here you generate warning if ptt>30
				fHistos->fhTrackingEfficiency[cBin]->Fill( ptt, 1./effCorr );
				track->SetTrackEff( 1./effCorr );
			}
		} // end of aodtrack
	}

	// Run the flow analysis here
	if(fDebugMode) cout << "Start of RunEbEFlowAnalysis.. FilterBit (AOD,Ntrack)="<< trkfilterBit << ","<< fInputList->GetEntries() <<endl;
	// Spectra Ana
	if(fDebugMode) cout << "Spectra Ana fInputListSpectra->GetEntriesFast() = "<< fInputListSpectra->GetEntriesFast() << endl;
	if(fDebugMode) cout << " fPtHardMin = "<< fPtHardMin <<" , fPtHardMax = "<< fPtHardMax << endl;

	// find LP in a event
	float pT_max = -1.;     //maximal pT in the event
	int indexLP      = -1;      //index of the leading particle
	for(int i=0;i<fInputListSpectra->GetEntriesFast();i++) {
		AliJBaseTrack *track = (AliJBaseTrack*)fInputListSpectra->At(i);
		double ptt = track->Pt();
		if(ptt > pT_max ){  //Find the leading particle in the event
			indexLP  = i;
			pT_max = ptt;
		}
	}

	if(fESMethod==5) { // di-jet Asymm
		// Make a decision for running the analysis or not
		// analyze the events only if we find a pt > PtHardMin
		TObjArray *fjets = (TObjArray*)fJetTask->GetAliJJetList(fJetSel); // only selected jet
		AliJJet *Ljet = dynamic_cast<AliJJet*>( fjets->At(0) );
		AliJJet *subLjet = dynamic_cast<AliJJet*>( fjets->At(1) );
		double minSubLeadingJetPt = 5.0;
		if(  Ljet && subLjet ) { 
			double ptt = Ljet->Pt();
			fHistos->fhJetPt[cBin]->Fill(ptt);
			double asym = (Ljet->Pt() - subLjet->Pt())/(Ljet->Pt() + subLjet->Pt());
			double InvM = ( *Ljet + *subLjet).M();
			fHistos->fhDiJetAsym[cBin]->Fill( asym );
			fHistos->fhRecoDiJetM[cBin]->Fill( InvM );
			if( ptt > fPtHardMin && ptt < fPtHardMax && subLjet->Pt()>minSubLeadingJetPt && asym > fDiJetAsymMin ) TagThisEvent = kTRUE;
		}
	} else if(fESMethod==4) { // di-jet
		// Make a decision for running the analysis or not
		// analyze the events only if we find a pt > PtHardMin
		TObjArray *fjets = (TObjArray*)fJetTask->GetAliJJetList(fJetSel); // only selected jet
		AliJJet *Ljet = dynamic_cast<AliJJet*>( fjets->At(0) );
		AliJJet *subLjet = dynamic_cast<AliJJet*>( fjets->At(1) );
		double minSubLeadingJetPt = 5.0;
		if(  Ljet && subLjet ) { 
			double ptt = Ljet->Pt();
			fHistos->fhJetPt[cBin]->Fill(ptt);
			if( ptt > fPtHardMin && ptt < fPtHardMax && subLjet->Pt()>minSubLeadingJetPt ) TagThisEvent = kTRUE;
			if( TagThisEvent ) {
				double asym = (Ljet->Pt() - subLjet->Pt())/(Ljet->Pt() + subLjet->Pt());
				double InvM = ( *Ljet + *subLjet).M();
				fHistos->fhDiJetAsym[cBin]->Fill( asym );
				fHistos->fhRecoDiJetM[cBin]->Fill( InvM );
			}
		}
	} else if(fESMethod==3) { // Leading jet
		// Make a decision for running the analysis or not
		// analyze the events only if we find a pt > PtHardMin
		TObjArray *fjets = (TObjArray*)fJetTask->GetAliJJetList(fJetSel); // only selected jet
		AliJJet *Ljet = dynamic_cast<AliJJet*>( fjets->At(0) );
		if (Ljet) {
			double ptt = Ljet->Pt();
			fHistos->fhJetPt[cBin]->Fill(ptt);
			if( ptt > fPtHardMin && ptt < fPtHardMax ) TagThisEvent = kTRUE;
		}
	} else if(fESMethod==2) { // jet
		// Make a decision for running the analysis or not
		// analyze the events only if we find a pt > PtHardMin
		TObjArray *fjets = (TObjArray*)fJetTask->GetAliJJetList(fJetSel); // only selected jet
		for (int ijet = 0; ijet<fjets->GetEntriesFast(); ijet++){
			AliJJet *jet = dynamic_cast<AliJJet*>( fjets->At(ijet) );
			double ptt = jet->Pt();
			fHistos->fhJetPt[cBin]->Fill(ptt);
			if( ptt > fPtHardMin && ptt < fPtHardMax ) TagThisEvent = kTRUE;
		}
	} else if(fESMethod==1) { // Leading particle
		// Make a decision for running the analysis or not
		// analyze the events only if we find a pt > PtHardMin
		if( pT_max > fPtHardMin && pT_max < fPtHardMax ) TagThisEvent = kTRUE;
	} else if (fESMethod==0) {
		// Find whether this event has a partilce at least pt > PtHardMin
		for(int i=0;i<fInputListSpectra->GetEntriesFast();i++) {
			AliJBaseTrack *track = (AliJBaseTrack*)fInputListSpectra->At(i);
			double ptt = track->Pt();
			// Make a decision for running the analysis or not
			// analyze the events only if we find a pt > PtHardMin
			if( ptt > fPtHardMin && ptt < fPtHardMax ) TagThisEvent = kTRUE;
		}
	}

	if(TagThisEvent) {
		fHistos->fhiCentr->Fill(cBin);
		fHistos->fhCentr->Fill(fcent);

		for(int i=0;i<fInputListSpectra->GetEntriesFast();i++) {
			AliJBaseTrack *track = (AliJBaseTrack*)fInputListSpectra->At(i);
			double ptt = track->Pt();
			double etat = track->Eta();
			double effCorr = 1.0/track->GetTrackEff();

			fHistos->fhChargedEta->Fill(etat);
			fHistos->fhChargedPt[cBin]->Fill(ptt, effCorr );
			fHistos->fhChargedPtNoCorr[cBin]->Fill( ptt );
			fHistos->fhChargedEta->Fill(etat, effCorr);
			//fHistos->fhChargedPtJacek[cBin]->Fill(ptt, effCorr );
			fHistos->fhChargedPtJacek[cBin]->Fill(ptt, ptt>0 ? 1./ptt*effCorr : 0); //One CANNOT do 1/ptt here!! First unfold.
			if(track->GetCharge()>0.) fHistos->fhChargedPtJacekPos[cBin]->Fill(ptt, ptt>0 ? 1./ptt*effCorr : 0); //One CANNOT do 1/ptt here!! First unfold.
			if(track->GetCharge()<0.) fHistos->fhChargedPtJacekNeg[cBin]->Fill(ptt, ptt>0 ? 1./ptt*effCorr : 0); //One CANNOT do 1/ptt here!! First unfold.
			if( -0.8<etat && etat<-0.2) fHistos->fhChargedPtJacekEta[cBin][0]->Fill(ptt, ptt>0 ? 1./ptt*effCorr : 0);
			if( -0.2<etat && etat< 0.3) fHistos->fhChargedPtJacekEta[cBin][1]->Fill(ptt, ptt>0 ? 1./ptt*effCorr : 0);
			if(  0.3<etat && etat< 0.8) fHistos->fhChargedPtJacekEta[cBin][2]->Fill(ptt, ptt>0 ? 1./ptt*effCorr : 0);
			fHistos->fhChargedPtFiete->Fill(ptt, effCorr );
		}

		// Run the flow analysis only for PtHardMin Cut
		RegisterList(fInputListFlow, fInputList, 0.2, 5.); // this must be done here input for fFFlucAna and RunFlow..
		double Eta_min=0.4, Eta_max=0.8;
		double vertex[3] = { 0, 0, zVert};
		fFFlucAna->Init();
		fFFlucAna->SetInputList( fInputListFlow );
		fFFlucAna->SetEventCentrality( fcent );
		fFFlucAna->SetEventImpactParameter( GetImpactParFromCentrality(fcent) );
		fFFlucAna->SetEventVertex ( vertex );
		fFFlucAna->SetEtaRange( Eta_min, Eta_max);
		fFFlucAna->UserExec("");
	}

	PostData(1, fOutput);
}

//________________________________________________________________________
void AliJHSInterplayTask::Terminate(Option_t *) 
{

	cout<<"Sucessfully Finished"<<endl;

}


//________________________________________________________________________
bool AliJHSInterplayTask::IsGoodEvent(AliVEvent *event) {

	if(fRunTable->IsPP() && fAnaUtils->IsPileUpEvent(event)) {
		return kFALSE;
	} else {
		Bool_t triggeredEventMB = kFALSE; //init

		fHistos->fhEvents->Fill( 0 );

		AliInputEventHandler *pieh = (AliInputEventHandler*)AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
		Bool_t triggerkMB = pieh->IsEventSelected() & ( AliVEvent::kMB );

		if( triggerkMB ){
			triggeredEventMB = kTRUE;  //event triggered as minimum bias
			fHistos->fhEvents->Fill( 1 );
		}

		int frunNumber = pieh->GetEvent()->GetRunNumber();
		if(frunNumber < 0)
			cout << "ERROR: unknown run number" << endl;
		//AliJRunTable *fRunTable = & AliJRunTable::GetSpecialInstance();
		fRunTable->SetRunNumber( frunNumber );

		if(fRunTable->GetRunNumberToPeriod(frunNumber) == AliJRunTable::kLHC15o){
			const AliVVertex* vtTrc = event->GetPrimaryVertex();
			const AliVVertex* vtSPD = event->GetPrimaryVertexSPD();
			double covTrc[6],covSPD[6];
			vtTrc->GetCovarianceMatrix(covTrc);
			vtSPD->GetCovarianceMatrix(covSPD);
			double dz = vtTrc->GetZ()-vtSPD->GetZ();
			double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
			double errTrc = TMath::Sqrt(covTrc[5]);
			double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
			if(TMath::Abs(dz) > 0.2 || nsigTot > 10 || nsigTrc > 20)
				return kFALSE;
		
			AliMultSelection *pms = (AliMultSelection*)event->FindListObject("MultSelection");
			if(!pms){
				AliError("MultSelection unavailable.");
				return kFALSE;
			}

			Float_t v0mcent = pms->GetMultiplicityPercentile("V0M");
			Float_t cl0cent = pms->GetMultiplicityPercentile("CL0");
			if(cl0cent < pfOutlierLowCut->Eval(v0mcent) || cl0cent > pfOutlierHighCut->Eval(v0mcent))
				return kFALSE;
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
				fHistos->fhEvents->Fill( 2 );
				if(fCard->VertInZRange(zVert)) {
					goodRecVertex = kTRUE; 
					fHistos->fhEvents->Fill( 3 );
					zBin  = fCard->GetBin(kZVertType, zVert); 
				}
			}
		}
		return goodRecVertex;
	}
	//---------------------------------
}


void AliJHSInterplayTask::ReadKineTracks( AliMCEvent *mcEvent, TClonesArray *TrackList)
{
	Bool_t IsExcludeWeakDecay = kTRUE;
	Int_t nt = mcEvent->GetNumberOfPrimaries();
	Int_t ntrack =0;
	for (Int_t it = 0; it < nt; it++) {
		AliMCParticle* track = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(it));
		if(mcEvent->IsPhysicalPrimary(it)) {
			if(!fCard->IsInEtaRange(track->Eta())) continue;
			double Pt = track->Pt();
			TParticle *particle = track->Particle();
			if(!particle) continue;

			if( IsExcludeWeakDecay == kTRUE){
				Bool_t kExcludeParticle = kFALSE;
				Int_t gMotherIndex = particle->GetFirstMother(); //
				if(gMotherIndex != -1){ // -1 means don't have mother.
					//DEBUG( 4,  "this particle has a mother " );
					AliMCParticle* motherParticle= dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(gMotherIndex));
					if(motherParticle){
						if(IsThisAWeakDecayingParticle( motherParticle)){
							kExcludeParticle = kTRUE;
							//DEBUG ( 4, Form("this particle will be removed because it comes from : %d", motherParticle->PdgCode() ));
						}
					}
				}
				if(kExcludeParticle) continue;
			}

			Int_t pdg = particle->GetPdgCode();
			Char_t ch = (Char_t) track->Charge();
			Int_t label = track->GetLabel();
			AliJBaseTrack *itrack = new ((*TrackList)[ntrack++])AliJBaseTrack;
			itrack->SetLabel( label );
			itrack->SetID(TrackList->GetEntriesFast());
			itrack->SetParticleType(kJHadron);
			itrack->SetPxPyPzE( track->Px(), track->Py(), track->Pz(), track->E() );
			itrack->SetCharge(ch) ;
			double effCorr = 1.;  // here you generate warning if ptt>30
			itrack->SetTrackEff( 1./effCorr );
		}
	}
}

//______________________________________________________________________________
Bool_t AliJHSInterplayTask::IsThisAWeakDecayingParticle(AliMCParticle *thisGuy)
{
	// In order to prevent analyzing daughters from weak decays
	// - AMPT does not only strong decays, so IsPhysicalPrimary does not catch it
	Int_t pdgcode = TMath::Abs( thisGuy->PdgCode() );
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

double AliJHSInterplayTask::GetCentralityFromImpactPar(double ip) {
	/*
	   \begin{tabular}{ |c|c c c c c c c c| }
	   \hline
	   Centrality(\%)&0-5      &5-10     &10-20    &20-30    &30-40     &40-50      &50-60      &60-70\\
	   \hline
	   b(fm) AMPT    &0.00-3.72&3.72-5.23&5.23-7.31&7.31-8.88&8.88-10.20&10.20-11.38&11.38-12.47&12.47-13.50\\
	   b(fm) HIJING  &0.00-3.60&3.60-5.09&5.09-7.20&7.20-8.83&8.83-10.20&10.20-11.40&11.40-12.49&12.49-13.49\\
	   b(fm) ALICE   &0.00-3.50&3.50-4.94&4.94-6.98&6.98-    &    -9.88 &9.81-      &     -12.09&12.09-\\
	   \hline
	   \end{tabular}
	   \begin{tablenotes}
	   \url{https://twiki.cern.ch/twiki/bin/viewauth/ALICE/CentStudies}
	   \end{tablenotes}
	 */
	double bmin[10]={0.0,3.72,5.23,7.31,8.88,10.20,11.38,12.47,14.51,100};
	double centmean[10]={2.5,7.5,15,25,35,45,55,65,75,90};
	int iC = -1;
	for(int i=0;i<9;i++){
		if(bmin[i]<ip&&ip<=bmin[i+1]) {iC=i;  break;}
	}
	return centmean[iC];
}

double AliJHSInterplayTask::GetImpactParFromCentrality(double cent) {
	/*
	   \begin{tabular}{ |c|c c c c c c c c| }
	   \hline
	   Centrality(\%)&0-5      &5-10     &10-20    &20-30    &30-40     &40-50      &50-60      &60-70\\
	   \hline
	   b(fm) AMPT    &0.00-3.72&3.72-5.23&5.23-7.31&7.31-8.88&8.88-10.20&10.20-11.38&11.38-12.47&12.47-13.50\\
	   b(fm) HIJING  &0.00-3.60&3.60-5.09&5.09-7.20&7.20-8.83&8.83-10.20&10.20-11.40&11.40-12.49&12.49-13.49\\
	   b(fm) ALICE   &0.00-3.50&3.50-4.94&4.94-6.98&6.98-    &    -9.88 &9.81-      &     -12.09&12.09-\\
	   \hline
	   \end{tabular}
	   \begin{tablenotes}
	   \url{https://twiki.cern.ch/twiki/bin/viewauth/ALICE/CentStudies}
	   \end{tablenotes}
	 */
	double bmin[10]={0.0,3.72,5.23,7.31,8.88,10.20,11.38,12.47,14.51,100};
	double centmin[11]={0,5,10,20,30,40,50,60,70,80,100};
	int iC = -1;
	for(int i=0;i<9;i++){
		if(centmin[i]<cent&&cent<=centmin[i+1]) {iC=i;  break;}
	}
	return (bmin[iC]+bmin[iC+1])/2.;
}

void AliJHSInterplayTask::RegisterList(TClonesArray* listToFill, TClonesArray* listFromToFill,double lpt, double hpt){
	int count = 0;
	for(int i=0;i<listFromToFill->GetEntriesFast();i++) {
		AliJBaseTrack *track = (AliJBaseTrack*)listFromToFill->At(i);
		double ptt = track->Pt();
		if(ptt > lpt && ptt<hpt ) new ((*listToFill)[count++]) AliJBaseTrack(*track);
	}
}
