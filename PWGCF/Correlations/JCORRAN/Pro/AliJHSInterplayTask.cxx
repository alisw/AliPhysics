#include "AliJBaseTrack.h"
#include "AliJHSInterplayTask.h"
#include "AliAnalysisManager.h"
#include "AliJJet.h"
#include "AliJCDijetAna.h"


ClassImp(AliJHSInterplayTask)

AliJHSInterplayTask::AliJHSInterplayTask() 
	: AliAnalysisTaskSE(), fOutput(0x0),
	fJCatalystTask(NULL),
	fJCatalystTaskName("JCatalystTask"),
	fFFlucAna(NULL),
	fJetTask(NULL),
	fJetTaskName("JJetTask"),
	fJetSel(0),
	fJFJTask(NULL),
	fJFJTaskName("JFJTask"),
	fCard(NULL),
	fHistos(NULL),
	fInputListSpectra(NULL),
	fInputListFlow(NULL),
	fVnMethod(0),
	jettask_tagging(1),
	fESMethod(0),
	cBin(-1),
	zBin(-1),
	zVert(-999),
	fevt(0),
	fDebugMode(0),
	fPtHardMin(0),
	fPtHardMax(0),
	flags(0),
	TagThisEvent()
{
	for(int iE=0;iE<kNESE;iE++) TagThisEvent[iE]=kFALSE; // init
}
//________________________________________________________________________
	AliJHSInterplayTask::AliJHSInterplayTask(const char *name) 
: AliAnalysisTaskSE(name), 
	fOutput(0x0), 
	fJCatalystTask(NULL),
	fJCatalystTaskName("JCatalystTask"),
	fFFlucAna(0x0),
	fJetTask(0x0),
	fJetTaskName(""),
	fJetSel(0),
	fJFJTask(NULL),
	fJFJTaskName("JFJTask"),
	fCard(0x0),
	fHistos(0x0),
	fInputListSpectra(0x0),
	fInputListFlow(0x0),
	fVnMethod(0),
	jettask_tagging(1),
	fESMethod(0),
	cBin(-1),
	zBin(-1),
	zVert(-999),
	fevt(0),
	fDebugMode(0),
	fPtHardMin(0),
	fPtHardMax(0),
	flags(0),
	TagThisEvent()
{
	// Define input and output slots here
	for(int iE=0;iE<kNESE;iE++) TagThisEvent[iE]=kFALSE; // init
	// Input slot #0 works with a TChain
	DefineInput(0, TDirectory::Class());
	// Output slot #0 writes into a TH1 container
	DefineOutput(1, TDirectory::Class());
}

//________________________________________________________________________
AliJHSInterplayTask::AliJHSInterplayTask(const AliJHSInterplayTask& a):
	AliAnalysisTaskSE(a.GetName()),
	fOutput(a.fOutput),
	fJCatalystTask(a.fJCatalystTask),
	fJCatalystTaskName(a.fJCatalystTaskName),
	fFFlucAna(a.fFFlucAna),
	fJetTask(a.fJetTask),
	fJetTaskName(a.fJetTaskName),
	fJetSel(a.fJetSel),
    fJFJTask(a.fJFJTask),
	fJFJTaskName(a.fJFJTaskName),	
	fCard(a.fCard),
	fHistos(a.fHistos),
	fInputListSpectra(a.fInputListSpectra),
	fInputListFlow(a.fInputListFlow),
	fVnMethod(a.fVnMethod),
	jettask_tagging(a.jettask_tagging),
	fESMethod(a.fVnMethod),
	cBin(-1),
	zBin(-1),
	zVert(-999),
	fevt(-1),
	fDebugMode(0), 
	fPtHardMin(a.fPtHardMin),
	fPtHardMax(a.fPtHardMax),
	flags(a.flags),
	TagThisEvent()
{

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
	cout<<"\n=============== CARD =============="<<endl;
	fCard->PrintOut();
	cout<<"===================================\n"<<endl;
	fInputListSpectra  = new TClonesArray("AliJBaseTrack",1500);
	fInputListFlow  = new TClonesArray("AliJBaseTrack",1500);
	fInputListSpectra->SetOwner(kTRUE);
	fInputListFlow->SetOwner(kTRUE);

	OpenFile(1);
	fOutput = gDirectory;
	fOutput->cd();

	fHistos = new AliJHistos(fCard);
	fHistos->CreateEventTrackHistos();
	fHistos->CreateJetHistos();
	fHistos->fHMG->Print();

	fESMethod 	= int (fCard->Get("ESMethod"));
	fJetSel		= int (fCard->Get("JetSel"));
	fDiJetAsymMin = int (fCard->Get("DiJetAsymMin"));

	for(int iE=0;iE<kNESE;iE++) TagThisEvent[iE]=kFALSE; // init

	//=== Get AnalysisManager
	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
	// Add JCatalystTask
	fJCatalystTask = (AliJCatalystTask*)(man->GetTask( fJCatalystTaskName ));
	// Add JJet
	fJetTask = (AliJJetTask*)(man->GetTask( fJetTaskName));
	// Add JFJ fastjet standalone task.
	fJFJTask = (AliJFJTask*)(man->GetTask( fJFJTaskName));
	// Add a AliJFlucAnalysis
	fFFlucAna = new AliJFFlucAnalysis("JFFlucAnalysis");
	fFFlucAna->SetBinning(AliJFFlucAnalysis::BINNING_CENT_PbPb);
	if(flags & HSINT_SCPT)
		fFFlucAna->AddFlags(AliJFFlucAnalysis::FLUC_SCPT);
	if(flags & HSINT_EBE_WEIGHTING)  // centality
		fFFlucAna->AddFlags(AliJFFlucAnalysis::FLUC_EBE_WEIGHTING); // only for 18q+r no need for 15o
	if(flags & HSINT_PHI_CORRECTION)
		fFFlucAna->AddFlags(AliJFFlucAnalysis::FLUC_PHI_CORRECTION);
	// fFFlucAna->SetQCEtaCut( -0.8, 0.8, 0.4 ); // not used anymore and need to check with JP

	fOutput->cd();
	//fFFlucAna->SetEffConfig( fJCatalystTask->GetEffMode(), fJCatalystTask->GetEffFilterBit() );
	fFFlucAna->UserCreateOutputObjects();

	//fCard->WriteCard(fOutput);
	fHistos->fHMG->WriteConfig();
	fevt = -1;
	PostData(1, fOutput);
}
//________________________________________________________________________
AliJHSInterplayTask::~AliJHSInterplayTask() {
	delete fOutput; 
	delete fFFlucAna;
	delete fInputListSpectra;
	delete fCard;
}

//________________________________________________________________________
void AliJHSInterplayTask::UserExec(Option_t *) {
	// Main loop
	// Called for each event
	fevt++;
	if(fevt % 1000 == 0) cout << "Number of events scanned = "<< fevt << endl;
	if( fJCatalystTask->GetJCatalystEntry() != fEntry ) return;
	// Event selection is done in fJCatalystTask
	// need to put cent/vtx/run# fJCatalystTask->GetRunNumber() fJCatalystTask->GetCentrality()fJCatalystTask->GetZVertex()
	// clear them up for every event
	zVert = fJCatalystTask->GetZVertex();
	float fcent = fJCatalystTask->GetCentrality();

	TObjArray *fInputList = (TObjArray*)fJCatalystTask->GetInputList();
	fInputListSpectra->Clear();
	fInputListFlow->Clear();

    RegisterList(fInputListSpectra, fInputList, 0., 100.);
    	// Run the flow analysis only for PtHardMin Cut
	RegisterList(fInputListFlow, fInputList, 0.2, 5.); // this must be done here input for fFFlucAna and RunFlow..

	cBin = fCard->GetBin(kCentrType, fcent);
	if(cBin<0) return;

	if(fDebugMode) cout << "Spectra Ana fInputListSpectra->GetEntriesFast() = "<< fInputListSpectra->GetEntriesFast() << endl;
	if(fDebugMode) cout << " fPtHardMin = "<< fPtHardMin <<" , fPtHardMax = "<< fPtHardMax << endl;
	if(fInputListSpectra->GetEntriesFast() < 1) return;
	// find LP in a event
	float pT_max = -1.;     //maximal pT in the event
	for(int i=0;i<fInputListSpectra->GetEntriesFast();i++) {
		AliJBaseTrack *track = (AliJBaseTrack*)fInputListSpectra->At(i);
		double ptt = track->Pt();
		if(ptt > pT_max ){  //Find the leading particle in the event
			pT_max = ptt;
		}
	}
	/*
	for(int iE=0;iE<kNESE;iE++) TagThisEvent[iE]=kFALSE; // init for every events
	for(int iE=0;iE<kNESE;iE++) ESETagging(jettask_tagging, iE, pT_max);
	for(int iE=0;iE<kNESE;iE++) cout << TagThisEvent[iE]<<":";
		cout << "LP pt="<< pT_max << endl;
	*/
	// Now only for selected one
    TagThisEvent[fESMethod]=kFALSE; // init for every events
    ESETagging(jettask_tagging, fESMethod, pT_max); // event tagging
    if(fDebugMode) cout << fESMethod <<"\t" << TagThisEvent[fESMethod] << endl;
	if(TagThisEvent[fESMethod]) {
		fHistos->fhiCentr->Fill(cBin);
		fHistos->fhCentr->Fill(fcent);

		for(int i=0;i<fInputListSpectra->GetEntriesFast();i++) {
			AliJBaseTrack *track = (AliJBaseTrack*)fInputListSpectra->At(i);
			double ptt = track->Pt();
			double etat = track->Eta();
			fHistos->fhChargedEta->Fill(etat);
			fHistos->fhChargedPtNoCorr[cBin]->Fill( ptt );
			fHistos->fhChargedPtJacek[cBin]->Fill(ptt, ptt>0 ? 1./ptt : 0); 
			if(track->GetCharge()>0.) fHistos->fhChargedPtJacekPos[cBin]->Fill(ptt, ptt>0 ? 1./ptt : 0); //One CANNOT do 1/ptt here!! First unfold.
			if(track->GetCharge()<0.) fHistos->fhChargedPtJacekNeg[cBin]->Fill(ptt, ptt>0 ? 1./ptt : 0); //One CANNOT do 1/ptt here!! First unfold.
			if( -0.8<etat && etat<-0.2) fHistos->fhChargedPtJacekEta[cBin][0]->Fill(ptt, ptt>0 ? 1./ptt : 0);
			if( -0.2<etat && etat< 0.3) fHistos->fhChargedPtJacekEta[cBin][1]->Fill(ptt, ptt>0 ? 1./ptt : 0);
			if(  0.3<etat && etat< 0.8) fHistos->fhChargedPtJacekEta[cBin][2]->Fill(ptt, ptt>0 ? 1./ptt : 0);
			fHistos->fhChargedPtFiete->Fill(ptt);
		}
	
		double Eta_min=0.4, Eta_max=0.8;
		double vertex[3] = { 0, 0, zVert};
		fFFlucAna->Init();
		fFFlucAna->SetInputList( fInputListFlow );
		fFFlucAna->SetEventCentrality( fcent );
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

void AliJHSInterplayTask::ESETagging(int itask, int iESE, double pT_max) {
	double Ljetpt = -999;
	double subLjetpt = -999;
	double minSubLeadingJetPt = 5.0;
	double asym = -999;
	double InvM = -999;
    if(itask==1) {
#if !defined(__CINT__) && !defined(__MAKECINT__)
        AliJCDijetAna *fJFJAna = (AliJCDijetAna*)fJFJTask->GetJCDijetAna();
	    vector<vector<fastjet::PseudoJet>> jets = fJFJAna->GetJets();
	    int iBGSubtr = 1; // private AliJCDijetAna::iBGSubtr
	    jets.at(iBGSubtr) = fastjet::sorted_by_pt(jets.at(iBGSubtr));
	    if(fDebugMode) {
	 	   for (unsigned i = 0; i < jets.at(iBGSubtr).size(); i++) {
			  cout << jets.at(iBGSubtr).at(i).pt()<<",";
	    	}
	    	cout << endl;
		}

		vector<fastjet::PseudoJet> psjets = jets.at(iBGSubtr); // only selected jet
		if(psjets.size() == 0) return;
		fastjet::PseudoJet psLjet = jets.at(iBGSubtr).at(0);
		Ljetpt = psLjet.pt();
		if(psjets.size()>1) {
			fastjet::PseudoJet pssubLjet = jets.at(iBGSubtr).at(1);
			subLjetpt = pssubLjet.pt();
			asym = (Ljetpt - subLjetpt)/(Ljetpt + subLjetpt);
			InvM = ( psLjet + pssubLjet).m();
		}
		if(fDebugMode) cout << Form("LPpt=%.1f:LPjet=%.1f:subJet=%.1f:DiJetAsym=%.1f",pT_max,Ljetpt,subLjetpt,asym) << endl;
		switch(iESE) {
		  case 0: // Leading particle
				TagThisEvent[iESE] = kTRUE;
				break;
		  case 1: // Leading particle
				if( pT_max > fPtHardMin && pT_max < fPtHardMax ) TagThisEvent[iESE] = kTRUE;
				break;
		  case 2: // jet
				for (unsigned ijet = 0; ijet<psjets.size(); ijet++){
					double ptt = psjets.at(ijet).pt();
					fHistos->fhJetPt[cBin]->Fill(ptt);
					if( ptt > fPtHardMin && ptt < fPtHardMax ) TagThisEvent[iESE] = kTRUE;
				}
				break;
		  case 3:  // Leading jet	
				fHistos->fhJetPt[cBin]->Fill(Ljetpt);
				if( Ljetpt > fPtHardMin && Ljetpt < fPtHardMax ) TagThisEvent[iESE] = kTRUE;
				break;
		  case 4: // di-jet
				fHistos->fhJetPt[cBin]->Fill(Ljetpt);
				if( Ljetpt > fPtHardMin && Ljetpt < fPtHardMax && subLjetpt > minSubLeadingJetPt ) TagThisEvent[iESE] = kTRUE;
				if( TagThisEvent[iESE] ) {
					fHistos->fhDiJetAsym[cBin]->Fill( asym );
					fHistos->fhRecoDiJetM[cBin]->Fill( InvM );
				}
				break;					
		  case 5: // di-jet Asymm
				fHistos->fhJetPt[cBin]->Fill(Ljetpt);
				fHistos->fhDiJetAsym[cBin]->Fill( asym );
				fHistos->fhRecoDiJetM[cBin]->Fill( InvM );
				if( Ljetpt > fPtHardMin && Ljetpt < fPtHardMax && subLjetpt > minSubLeadingJetPt && asym > fDiJetAsymMin ) TagThisEvent[iESE] = kTRUE;
				break;
		} // end of switch
#endif 
    } // jet task

    if(itask==0) {
		TObjArray *fjets = (TObjArray*)fJetTask->GetAliJJetList(fJetSel); // only selected jet
		AliJJet *Ljet = dynamic_cast<AliJJet*>( fjets->At(0) );
		AliJJet *subLjet = dynamic_cast<AliJJet*>( fjets->At(1) );
		Ljetpt = Ljet->Pt();
		subLjetpt = subLjet->Pt();
		asym = (Ljetpt - subLjetpt)/(Ljetpt + subLjetpt);
		InvM = ( *Ljet + *subLjet).M();
		
		switch(iESE) {
		  case 0: // Leading particle
				TagThisEvent[iESE] = kTRUE;
				break;
		  case 1: // Leading particle
				if( pT_max > fPtHardMin && pT_max < fPtHardMax ) TagThisEvent[iESE] = kTRUE;
				break;
		  case 2: // jet
				for (int ijet = 0; ijet<fjets->GetEntriesFast(); ijet++){
					AliJJet *jet = dynamic_cast<AliJJet*>( fjets->At(ijet) );
					double ptt = jet->Pt();
					fHistos->fhJetPt[cBin]->Fill(ptt);
					if( ptt > fPtHardMin && ptt < fPtHardMax ) TagThisEvent[iESE] = kTRUE;
				}
				break;
		  case 3:  // Leading jet	
				fHistos->fhJetPt[cBin]->Fill(Ljetpt);
				if( Ljetpt > fPtHardMin && Ljetpt < fPtHardMax ) TagThisEvent[iESE] = kTRUE;
				break;
		  case 4: // di-jet
				fHistos->fhJetPt[cBin]->Fill(Ljetpt);
				if( Ljetpt > fPtHardMin && Ljetpt < fPtHardMax && subLjetpt > minSubLeadingJetPt ) TagThisEvent[iESE] = kTRUE;
				if( TagThisEvent[iESE] ) {
					fHistos->fhDiJetAsym[cBin]->Fill( asym );
					fHistos->fhRecoDiJetM[cBin]->Fill( InvM );
				}
				break;					
		  case 5: // di-jet Asymm
				fHistos->fhJetPt[cBin]->Fill(Ljetpt);
				fHistos->fhDiJetAsym[cBin]->Fill( asym );
				fHistos->fhRecoDiJetM[cBin]->Fill( InvM );
				if( Ljetpt > fPtHardMin && Ljetpt < fPtHardMax && subLjetpt > minSubLeadingJetPt && asym > fDiJetAsymMin ) TagThisEvent[iESE] = kTRUE;
				break;
		} // end of switch
    } // jet task
}
//--------------------------------
void AliJHSInterplayTask::RegisterList(TClonesArray* listToFill, TObjArray* listFromToFill,double lpt, double hpt){
	int count = 0;
	for(int i=0;i<listFromToFill->GetEntriesFast();i++) {
		AliJBaseTrack *track = (AliJBaseTrack*)listFromToFill->At(i);
		double ptt = track->Pt();
		if(ptt > lpt && ptt<hpt ) new ((*listToFill)[count++]) AliJBaseTrack(*track);
	}
}
