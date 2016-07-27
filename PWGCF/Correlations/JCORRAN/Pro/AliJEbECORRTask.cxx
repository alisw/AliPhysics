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
#include "AliJEbECORRTask.h"
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


ClassImp(AliJEbECORRTask)

	AliJEbECORRTask::AliJEbECORRTask() 
: AliAnalysisTaskSE(), fOutput(0x0),
	fAnaUtils(NULL),
	fCard(NULL),
	fHistos(NULL),
	fEbeHistos(NULL),
	fEfficiency(NULL),
	fHadronSelectionCut(0),
	fInputList(NULL),
	ftriggList(NULL),
	fassocList(NULL),
	fcorrelations(NULL),
	fassocPool(NULL),
	fEbePercentile(NULL),
	fFirstEvent(kTRUE),
	cBin(-1),
	ebeBin(-1),
	zBin(-1),
	zVert(-999),
	fevt(0),
	fDebugMode(0),
	trkfilterBit(0),
	fEbECentBinBorders(NULL),
	ebePercentileInputFileName(""),
	fRunTable(0),
	fRandom(NULL),
	IsMC(kFALSE),
	fenableCORR(kFALSE)
{
	// Constructor
}
//________________________________________________________________________
AliJEbECORRTask::AliJEbECORRTask(const char *name) 
	: AliAnalysisTaskSE(name), 
	fOutput(0x0), 
	fAnaUtils(0x0),
	fCard(0x0),
	fHistos(0x0),
	fEbeHistos(0x0),
	fEfficiency(0x0),
	fHadronSelectionCut(0),
	fInputList(0x0),
	ftriggList(0x0),
	fassocList(0x0),
	fcorrelations(0x0),
	fassocPool(0x0),
	fEbePercentile(0x0),
	fFirstEvent(kTRUE),
	cBin(-1),
	ebeBin(-1),
	zBin(-1),
	zVert(-999),
	fevt(0),
	fDebugMode(0),
	trkfilterBit(0),
	fEbECentBinBorders(0x0),
	ebePercentileInputFileName(""),
	fRunTable(0),
	fRandom(0x0),
	IsMC(kFALSE),
	fenableCORR(kFALSE) 
{

	// Constructor


	// Define input and output slots here
	// Input slot #0 works with a TChain
	DefineInput(0, TChain::Class());
	// Output slot #0 writes into a TH1 container
	DefineOutput(1, TDirectory::Class());
}

//________________________________________________________________________
AliJEbECORRTask::AliJEbECORRTask(const AliJEbECORRTask& a):
	AliAnalysisTaskSE(a.GetName()),
	fOutput(a.fOutput),
	fAnaUtils(a.fAnaUtils),
	fCard(a.fCard),
	fHistos(a.fHistos),
	fEbeHistos(a.fEbeHistos),
	fEfficiency(a.fEfficiency),
	fHadronSelectionCut(a.fHadronSelectionCut),
	fInputList(a.fInputList),
	ftriggList(a.ftriggList),
	fassocList(a.fassocList),
	fcorrelations(a.fcorrelations),
	fassocPool(a.fassocPool),
	fEbePercentile(a.fEbePercentile),
	fFirstEvent(kTRUE),
	cBin(-1),
	ebeBin(-1),
	zBin(-1),
	zVert(-999),
	fevt(-1),
	fDebugMode(0), 
	trkfilterBit(a.trkfilterBit),
	fEbECentBinBorders(a.fEbECentBinBorders),
	ebePercentileInputFileName(""),
	fRunTable(a.fRunTable),
	fRandom(a.fRandom),
	IsMC(a.IsMC),
	fenableCORR(a.fenableCORR)
{
	//copy constructor
}
//________________________________________________________________________
AliJEbECORRTask& AliJEbECORRTask::operator = (const AliJEbECORRTask& ap){
	// assignment operator

	JUNUSED(ap);
	this->~AliJEbECORRTask();
	new(this) AliJEbECORRTask(ap);
	return *this;
}

Bool_t AliJEbECORRTask::UserNotify() {

	cout <<"UserNotify  ievt =" << fevt << endl;
	return kTRUE;
}


//________________________________________________________________________
void AliJEbECORRTask::UserCreateOutputObjects(){
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
	fInputList->SetOwner(kTRUE);
	ftriggList  = new TClonesArray("AliJBaseTrack",1500);
	fassocList  = new TClonesArray("AliJBaseTrack",1500);

	OpenFile(1);
	fOutput = gDirectory;
	fOutput->cd();

	fHistos = new AliJHistos(fCard);
	fHistos->CreateEventTrackHistos();
	if(fenableCORR) {
		fHistos->CreateAzimuthCorrHistos();
		//fHistos->CreateIAAMoons();
		//fHistos->CreateXEHistos();
	}
	fHistos->CreateXtHistos();
	//fHistos->CreatePairPtCosThetaStar();

	fHistos->fHMG->Print();
	// E-b-E
	fEbeHistos = new AliJEbeHistos(fCard);
	fEbeHistos->CreateUnfoldingHistos();

	fcorrelations = new AliJCorrelations(fCard, fHistos);
	fassocPool   = new AliJEventPool( fCard, fHistos, fcorrelations, kJHadron);
	fEbePercentile = new AliJEbePercentile(fCard, ebePercentileInputFileName);
	fEbECentBinBorders = fCard->GetVector("EbECentBinBorders");

	fEfficiency = new AliJEfficiency();	
	fEfficiency->SetMode( fCard->Get("EfficiencyMode") ); // 0:NoEff, 1:Period 2:RunNum 3:Auto
	if(IsMC) fEfficiency->SetMode( 0 );
	fEfficiency->SetDataPath("alien:///alice/cern.ch/user/d/djkim/legotrain/efficieny/data"); // Efficiency root file location local or alien

	fHadronSelectionCut = int ( fCard->Get("HadronSelectionCut"));
	trkfilterBit	= Int_t ( fCard->Get("AODTrackFilterBit"));

	fCard->WriteCard(fOutput);
	fHistos->fHMG->WriteConfig();
	fFirstEvent = kTRUE;
	fevt = -1;
	PostData(1, fOutput);
}
//________________________________________________________________________
AliJEbECORRTask::~AliJEbECORRTask() {
	delete fOutput; 
	delete fAnaUtils;
	delete fHistos;
	delete fEbeHistos;
	delete fEfficiency;
	delete fInputList;
	delete ftriggList;
	delete fassocList;
	delete fCard;
	delete fRandom;
	delete fcorrelations;
	if( fassocPool ) delete fassocPool;
	delete fEbePercentile;
}

//________________________________________________________________________
void AliJEbECORRTask::UserExec(Option_t *) {
	// Main loop
	// Called for each event
	fevt++;
	if(fevt % 1000 == 0) cout << "Numer of event scanned = "<< fevt << endl;

	AliVEvent *event = InputEvent();
	if(!event) return;

	//---------------------------------------------------------------
	// check if the event was triggered or not and vertex   

	AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(event);
	if(!aodEvent) return;

	if( fFirstEvent ) {
		fRunTable = & AliJRunTable::GetSpecialInstance();
		fRunTable->SetRunNumber( aodEvent->GetRunNumber() );
		fEfficiency->SetRunNumber( aodEvent->GetRunNumber() );
		fEfficiency->Load();
		fFirstEvent = kFALSE;
	}

	if(!IsGoodEvent( event )) return; // zBin is set there
	if(fDebug) cout << "zvtx = " << zVert << endl;

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

	fHistos->fhZVert[cBin]->Fill(zVert);
	if(fDebug ) cout <<"Centrality = "<< fcent <<"\t"<< cBin << endl;

	fHistos->fhEvents->Fill( 4 );
	Int_t nt = aodEvent->GetNumberOfTracks();


	// clear them up for every event
	fInputList->Clear();
	ftriggList->Clear();
	fassocList->Clear();

	int counter = 0;

	// Making the inputlist from AOD
	if( IsMC == kTRUE ){  // how to get a flag to check  MC or not !
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
	} // read mc track done.
	if( IsMC == kFALSE ){  
		for(Int_t it = 0; it < nt; it++) {
			AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(it));
			if( !aodtrack ) continue;
			if(aodtrack->TestFilterBit(trkfilterBit)) { 
				if(!fCard->IsInEtaRange(aodtrack->Eta())) continue;
				AliJBaseTrack* track = new( (*fInputList)[fInputList->GetEntriesFast()] ) AliJBaseTrack;
				track->SetPxPyPzE(aodtrack->Px(), aodtrack->Py(), aodtrack->Pz(), aodtrack->E());
				track->SetID(fInputList->GetEntriesFast());
				track->SetParticleType(kJHadron);
				track->SetCharge(aodtrack->Charge());
				double ptt = track->Pt();
				double effCorr = 1./fEfficiency->GetCorrection(ptt, fHadronSelectionCut, fcent);  // here you generate warning if ptt>30
				fHistos->fhTrackingEfficiency[cBin]->Fill( ptt, 1./effCorr );
				track->SetTrackEff( 1./effCorr );
			}
		} // end of aodtrack
	}

	// Run the flow analysis here
	if(fDebugMode) cout << "Start of RunEbEFlowAnalysis.. FilterBit (AOD,Ntrack)="<< trkfilterBit << ","<< nt <<","<< fInputList->GetEntries() <<endl;

	double ebeCent = RunEbEFlowAnalysis(event, fInputList);

	if(fDebugMode) fEbECentBinBorders->Print();
	if(fDebugMode) cout << "EbECentBinBorders="<< (*fEbECentBinBorders)[1]<<","<< (*fEbECentBinBorders)[2]<<endl;

	if( (*fEbECentBinBorders)[1] > ebeCent || ebeCent > (*fEbECentBinBorders)[2] ) return;

	fHistos->fhCentr->Fill(fcent);
	fHistos->fhiCentr->Fill(cBin);
    //only for v3 selection ih ==3

	if(fDebugMode) cout << "End of RunEbEFlowAnalysis.."<<endl;

	// Trigger and Assoc List
	int noTrigg=0,noAssoc=0;
	for(int i=0;i<fInputList->GetEntriesFast();i++) {
		AliJBaseTrack *track = (AliJBaseTrack*)fInputList->At(i);
		track->SetTriggBin( fCard->GetBin(kTriggType, track->Pt()) );
		track->SetAssocBin( fCard->GetBin(kAssocType, track->Pt()) );

		double ptt = track->Pt();
		double etat = track->Eta();
		double effCorr = 1.0/track->GetTrackEff();

		fHistos->fhChargedEta->Fill(etat);
		fHistos->fhChargedPt[cBin]->Fill(ptt, effCorr );
		fHistos->fhChargedPtNoCorr[cBin]->Fill( ptt );
		fHistos->fhChargedEta->Fill(etat, effCorr);
		//fHistos->fhChargedPtJacek[cBin]->Fill(ptt, effCorr );
		fHistos->fhChargedPtJacek[cBin]->Fill(ptt, ptt>0 ? 1./ptt*effCorr : 0); //One CANNOT do 1/ptt here!! First unfold.
		if( -0.8<etat && etat<-0.2) fHistos->fhChargedPtJacekEta[cBin][0]->Fill(ptt, ptt>0 ? 1./ptt*effCorr : 0);
		if( -0.2<etat && etat< 0.3) fHistos->fhChargedPtJacekEta[cBin][1]->Fill(ptt, ptt>0 ? 1./ptt*effCorr : 0);
		if(  0.3<etat && etat< 0.8) fHistos->fhChargedPtJacekEta[cBin][2]->Fill(ptt, ptt>0 ? 1./ptt*effCorr : 0);
		fHistos->fhChargedPtFiete->Fill(ptt, effCorr );
		if(fenableCORR) {
			if( track->GetTriggBin() > -1 ) { new ((*ftriggList)[noTrigg++]) AliJBaseTrack(*track); }//ftriggList->Add( track );
			if( track->GetAssocBin() > -1 ) { new ((*fassocList)[noAssoc++]) AliJBaseTrack(*track); }//
			if(track->GetAssocBin() > -1){
				int ipta  = track->GetAssocBin();
				double effCorrection = 1.0/track->GetTrackEff();
				fHistos->fhIphiAssoc[cBin][ipta]->Fill( track->Phi(), effCorrection);
				fHistos->fhIetaAssoc[cBin][ipta]->Fill( track->Eta(), effCorrection);
			}
		}
	}

	// correlation loop
	if(fDebugMode) cout << "Start of Correlation Loop noTrigg = "<< ftriggList->GetEntriesFast()<<"\t noAssoc="<<fassocList->GetEntriesFast()<<endl;
	if(fenableCORR) {
		if(fassocList->GetEntriesFast()>0 ) fassocPool->AcceptList(fassocList, fcent, zVert, fInputList->GetEntriesFast(), fevt);
		PlayCorrelation(ftriggList, fassocList);
		fassocPool->Mix(ftriggList, kAzimuthFill, fcent, zVert, fInputList->GetEntriesFast(), fevt);
	}
	PostData(1, fOutput);
}

//________________________________________________________________________
void AliJEbECORRTask::Terminate(Option_t *) 
{

	/*
	   for (int hic = 0;hic < fCard->GetNoOfBins(kCentrType);hic++){
	   ScaleNotEquidistantHisto( fHistos->fhChargedPt[hic], 1);
	   ScaleNotEquidistantHisto( fHistos->fhChargedPtNoCorr[hic], 1);
	   ScaleNotEquidistantHisto( fHistos->fhChargedPtJacek[hic], 1);
	   }
	   OpenFile(2);

	   fcorrelations->PrintOut();
	//fassocPool->PrintOut();
	 */
	cout<<"Sucessfully Finished"<<endl;

}


//________________________________________________________________________
bool AliJEbECORRTask::IsGoodEvent(AliVEvent *event) {

	if(fRunTable->IsPP() && fAnaUtils->IsPileUpEvent(event)) {
		return kFALSE;
	} else {
		Bool_t triggeredEventMB = kFALSE; //init

		fHistos->fhEvents->Fill( 0 );

		Bool_t triggerkMB = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & ( AliVEvent::kMB );

		if( triggerkMB ){
			triggeredEventMB = kTRUE;  //event triggered as minimum bias
			fHistos->fhEvents->Fill( 1 );
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

//________________________________________________________________________
double AliJEbECORRTask::RunEbEFlowAnalysis(AliVEvent *event, TClonesArray* inputlist){

	Double_t qxTot = 0.0, qyTot = 0.0;
	Double_t qxTotA = 0.0, qyTotA = 0.0;
	Double_t qxTotC[kNHarmonics] = {0.0}, qyTotC[kNHarmonics] = {0.0};
	double Qxa[kNHarmonics];
	double Qxb[kNHarmonics];
	double Qya[kNHarmonics];
	double Qyb[kNHarmonics];
	double Qxalt[kNHarmonics];
	double Qyalt[kNHarmonics];
	double Psi[kNHarmonics];
	double PsiA[kNHarmonics];
	double PsiC[kNHarmonics];
	double vnxa[kNHarmonics];
	double vnya[kNHarmonics];
	double vnxb[kNHarmonics];
	double vnyb[kNHarmonics];
	double vnx[kNHarmonics];
	double vny[kNHarmonics];
	double vobsalt[kNHarmonics];

	int firstH = 1;

	for( int ih=1;ih<kNHarmonics;ih++ ){
		fEbeHistos->fhEPCosndPhi[ih]->Reset();
		fEbeHistos->fhEPCosndPhi2[ih]->Reset();
	}
	AliEventplane *ep = event->GetEventplane();

	for( int ih=firstH;ih<kNHarmonics;ih++ ) Psi[ih] = 0;
	if(ep){
		for( int ih=firstH;ih<kNHarmonics;ih++ ){
			Psi[ih] = ep->CalculateVZEROEventPlane(event,10,ih,qxTot,qyTot);
			PsiA[ih] = ep->CalculateVZEROEventPlane(event,8,ih,qxTotA,qyTotA);
			PsiC[ih] = ep->CalculateVZEROEventPlane(event,9,ih,qxTotC[ih],qyTotC[ih]);
			fEbeHistos->fhQvectorV0[cBin][ih]->Fill(TMath::Sqrt(qxTot*qxTot+qyTot*qyTot));
			fEbeHistos->fhQvectorV0A[cBin][ih]->Fill(TMath::Sqrt(qxTotA*qxTotA+qyTotA*qyTotA));
			fEbeHistos->fhQvectorV0C[cBin][ih]->Fill(TMath::Sqrt(qxTotC[ih]*qxTotC[ih]+qyTotC[ih]*qyTotC[ih]));
			fEbeHistos->fhEventPlane[cBin][ih]->Fill(Psi[ih]);
			fEbeHistos->fhQvectorCorrelation[cBin][ih]->Fill(TMath::Sqrt(qxTotA*qxTotA+qyTotA*qyTotA),TMath::Sqrt(qxTotC[ih]*qxTotC[ih]+qyTotC[ih]*qyTotC[ih]));
		}
	}
	for(int ih = 0 ; ih < kNHarmonics; ih++){
		Qxa[ih] = 0;
		Qxb[ih] = 0;
		Qya[ih] = 0;
		Qyb[ih] = 0;
		Qxalt[ih]  = 0;
		Qyalt[ih]  = 0;
	}

	int counterA = 0;
	int counterB = 0;
	int counterAB = 0;
	int counterMulti = 0;
	double lowpTcut = 0.5;// ATLAS paper > 0.5

	for(int i=0;i<inputlist->GetEntriesFast();i++) {
		AliJBaseTrack *track = (AliJBaseTrack*)inputlist->At(i);
		//double w = tk->Pt();
		double w = 1;
		if(track->Pt()>lowpTcut){ 
			counterMulti++;
			if(fRandom->Uniform(1.0)<0.5){
				for(int ih = firstH ; ih < kNHarmonics ; ih++){
					Qxa[ih] = Qxa[ih] + w*TMath::Cos(ih*track->Phi());
					Qya[ih] = Qya[ih] + w*TMath::Sin(ih*track->Phi());
				}
				//counterA+=w;
				counterA++;
			}
			else{
				for(int ih = firstH ; ih < kNHarmonics ; ih++){
					Qxb[ih] = Qxb[ih] + w*TMath::Cos(ih*track->Phi());
					Qyb[ih] = Qyb[ih] + w*TMath::Sin(ih*track->Phi());
				}
				//      counterB+=w;
				counterB++;
			}
			for(int ih = firstH; ih < kNHarmonics; ih++){
				Qxalt[ih] = Qxalt[ih] + w*TMath::Cos(ih*track->Phi());
				Qyalt[ih] = Qyalt[ih] + w*TMath::Sin(ih*track->Phi());
			}
			//counterAB+=w;
			counterAB++;
		}
		for(int ih = firstH; ih < kNHarmonics ; ih++){
			fEbeHistos->fhEPCosndPhi[ih]->Fill( cos(ih*(track->Phi()-Psi[ih])) );
			if(track->Pt()>lowpTcut){
				fEbeHistos->fhEPCosndPhi2[ih]->Fill( cos(ih*(track->Phi()-Psi[ih])) );
			}
			fEbeHistos->fhCounter[cBin][ih]->Fill(track->Pt(),1);
		}
	} // track loop

	for(int ih = firstH ; ih < kNHarmonics ; ih++)  vobsalt[ih] = 0;
	// calculating Vn obs
	if(counterA > 1 && counterB > 1){

		for(int ih = firstH ; ih < kNHarmonics ; ih++){
			vnxa[ih] = Qxa[ih]/counterA;
			vnya[ih] = Qya[ih]/counterA;
			vnxb[ih] = Qxb[ih]/counterB;
			vnyb[ih] = Qyb[ih]/counterB;
			vnx[ih] = Qxalt[ih]/counterAB;
			vny[ih] = Qyalt[ih]/counterAB;
			vobsalt[ih] = TMath::Sqrt(vnx[ih]*vnx[ih]+vny[ih]*vny[ih]);
			fEbeHistos->fhVnObsVector[cBin][ih]->Fill(vobsalt[ih]);
			//cout << "ih: " << ih << " vnalt Average : " << vnalt[ih]->GetMean() << " Current vobs: " << vobsalt[ih] << " Counter: " << counterAB << endl;
			fEbeHistos->fhResponseDist[cBin][ih]->Fill(vnxa[ih]-vnxb[ih],vnya[ih]-vnyb[ih]);
		}
		//multiCount->Fill(counterMulti);
		fEbeHistos->fhMultiCount[cBin]->Fill(counterAB);
	}

	// Vn from EP 
	for(int ih = firstH; ih < kNHarmonics; ih++){
		fEbeHistos->fhVnObsEP[cBin][ih]->Fill(fEbeHistos->fhEPCosndPhi2[ih]->GetMean());
		//VnMulti[ih]->Fill(hEPCosndPhi2[ih]->GetEntries(),hEPCosndPhi2[ih]->GetMean());
	}
	for(int i=0;i<inputlist->GetEntriesFast();i++) {
		AliJBaseTrack *track = (AliJBaseTrack*)inputlist->At(i);
		for(int ih = firstH; ih < kNHarmonics; ih++)
			fEbeHistos->fheCosndPhiPt[cBin][ih]->Fill(track->Pt(), pow( (cos (ih* (track->Phi()-Psi[ih]) - fEbeHistos->fhCosndPhiPt[cBin][ih]->GetBinContent( ( fEbeHistos->fhCosndPhiPt[cBin][ih]->FindBin(track->Pt()))))),2));
	}

	// correlation VnObs vs Qvector from V0C
	for(int ih = firstH; ih < kNHarmonics; ih++){
		fEbeHistos->fhVnObsVsQvectorCorrelation[cBin][ih]->Fill(vobsalt[ih],TMath::Sqrt(qxTotC[ih]*qxTotC[ih]+qyTotC[ih]*qyTotC[ih]));
	}
	// E-b-E vn selection here
	int ih = 3;// for 3 rd harmonics
	double ebeCent = fEbePercentile->GetEbeFlowPercentile( cBin, ih, vobsalt[ih] );//From the data scan observed distribution at this moment, we need matrix to unfloded vn
    ebeCent *=100; // make it to a percentile;

    if( (*fEbECentBinBorders)[1] < ebeCent && ebeCent < (*fEbECentBinBorders)[2] ) {
        fEbeHistos->fhVnObsVectorAfterSelection[cBin][ih]->Fill(vobsalt[ih]);
    }
	return ebeCent;
}

//________________________________________________________________________
void AliJEbECORRTask::PlayCorrelation(TClonesArray *triggList, TClonesArray *assocList){
	int noTrigg = triggList->GetEntries();
	int noAssoc = assocList->GetEntries();

	fHistos->fhAssocMult->Fill(noAssoc);
	//------------------------------------------------------------------
	//==== Correlation Loop
	//------------------------------------------------------------------
	for(int ii=0;ii<noTrigg;ii++){ // trigger loop
		AliJBaseTrack * triggTr = (AliJBaseTrack*)triggList->At(ii);
		int iptt   = triggTr->GetTriggBin();
		if( iptt < 0 ) continue;
		double effCorr = 1.0/triggTr->GetTrackEff();
		fHistos->fhTriggPtBin[cBin][zBin][iptt]->Fill(triggTr->Pt(), effCorr);//inclusive
		for(int jj=0;jj<noAssoc;jj++){ // assoc loop
			AliJBaseTrack  *assocTr = (AliJBaseTrack*)assocList->At(jj);
			fcorrelations->FillAzimuthHistos(kReal, cBin, zBin, triggTr, assocTr);
		}
	} // end of trigg
}

void AliJEbECORRTask::ScaleNotEquidistantHisto(TH1D *hid, const double sc=1){
	// histo scaler
	for(int i=1;i<= hid->GetNbinsX();i++){
		hid->SetBinContent(i,hid->GetBinContent(i)*sc/hid->GetBinWidth(i));
		hid->SetBinError(i,hid->GetBinError(i)*sc/hid->GetBinWidth(i));
	}  
}
