/************************************************************************

   This macro produces: 
   - spherocity distributions (only works for reference estimator, |eta|<0.8) using ESD information. 
   - Aditya Nath Mishra, WIGNER RCP Budapest
   
   Please report bugs to: aditya.nath.mishra@cern.ch
   First version: 21/09/2021

 ************************************************************************/


#include "AliAnalysisTaskS0.h"

// ROOT includes
#include <TList.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2D.h>
#include <TParticle.h>
#include <TFile.h>

// AliRoot includes
#include <AliAnalysisManager.h>
#include <AliAnalysisFilter.h>
#include <AliESDInputHandler.h>
#include <AliESDEvent.h>
#include <AliEventCuts.h>
#include <AliESDVertex.h>
#include <AliLog.h>
#include <AliExternalTrackParam.h>
#include <AliESDtrackCuts.h>
#include <AliESDVZERO.h>
#include <AliAODVZERO.h>
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include <AliStack.h>
#include <TTreeStream.h>
#include <AliHeader.h>
#include <AliAnalysisUtils.h>
#include <AliSpherocityUtils.h>
#include <AliMultiplicity.h>
#include <AliMultSelection.h>
#include <AliAODInputHandler.h> 
#include <AliAODHandler.h> 
#include <AliAODVertex.h>
#include <AliAODTrack.h> 
#include <AliAODPid.h> 
#include <AliAODMCHeader.h> 
#include <AliDataFile.h>

#include <iostream>
using namespace std;

ClassImp(AliAnalysisTaskS0)

	//_____________________________________________________________________________
AliAnalysisTaskS0::AliAnalysisTaskS0():
AliAnalysisTaskSE(),
  fESD(0x0),
  fEventCuts(0x0),
  fMC(0x0),
  fMCStack(0x0),
  fSpheroUtils(0x0),
  fTrackFilter(0x0),
  fAnalysisType("ESD"),
  fAnalysisMC(kFALSE),	
  fnRefGlobal(0x0),
  fS0R(-1),
  fS0T(-1),
  fListOfObjects(0x0),	
  fHistEventCounter(0x0),
  hS0(0x0),
  hS0RefMult(0x0),         
  hpt(0x0),
  hS0pt(0x0),
  hphiso(0x0),
  hRefMult(0x0),
  hS0phi(0x0),
  hetaso(0x0),
  hS0eta(0x0),
  hMultTrue(0x0),
  hS0True(0x0),
  hS0Truept(0x0),
  hS0Truephi(0x0),
  hS0Trueeta(0x0),
  hS0TrueMult(0x0)
{
	// Default constructor (should not be used)
}

//______________________________________________________________________________
AliAnalysisTaskS0::AliAnalysisTaskS0(const char *name):
  AliAnalysisTaskSE(name),
  fESD(0x0),
  fEventCuts(0x0),
  fMC(0x0),
  fMCStack(0x0),
  fSpheroUtils(0x0),
  fTrackFilter(0x0),
  fAnalysisType("ESD"),
  fAnalysisMC(kFALSE),	
  fnRefGlobal(0x0),
  fS0R(-1),
  fS0T(-1),
  fListOfObjects(0x0),	
  fHistEventCounter(0x0),
  hS0(0x0),
  hS0RefMult(0x0),         
  hpt(0x0),
  hS0pt(0x0),
  hphiso(0x0),
  hRefMult(0x0),
  hS0phi(0x0),
  hetaso(0x0),
  hS0eta(0x0),
  hMultTrue(0x0),
  hS0True(0x0),
  hS0Truept(0x0),
  hS0Truephi(0x0),
  hS0Trueeta(0x0),
  hS0TrueMult(0x0)
{
	DefineOutput(1, TList::Class());
}

//________________________________________________________________________

void AliAnalysisTaskS0::Exit(const char *msg) {

	Printf("%s", msg);
	return;
}


//_____________________________________________________________________________
AliAnalysisTaskS0::~AliAnalysisTaskS0()
{
	// Destructor
	// histograms are in the output list and deleted when the output
	// list is deleted by the TSelector dtor
	if (fListOfObjects && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()){
		delete fListOfObjects;
		fListOfObjects = 0x0;
	}

}

//______________________________________________________________________________
void AliAnalysisTaskS0::UserCreateOutputObjects()
{
	// This method is called once per worker node
	// Here we define the output: histograms and debug tree if requested 

	if(!fSpheroUtils){
		fSpheroUtils = new AliSpherocityUtils();
		fSpheroUtils->Init();
	}

	// Definition of trackcuts
	if(!fTrackFilter){	
		fTrackFilter = new AliAnalysisFilter("trackFilter2015");
		SetTrackCuts(fTrackFilter);
	}
	// Helper to obtain the spherocity percentile
	//	fnsoB = fSoBining->GetNbinsX();
	//	fnMultbins = 100;
	

	const Int_t nPtBins      = 58;
	Double_t xBins[nPtBins+1] = {
		0.1,0.12,0.15,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
		0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,
		1.5,1.6,1.7,1.8,1.9,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,
		4.5,5,5.5,6,6.5,7,8,9,10,11,12,13,14,15,16,18,20.0
	};

	OpenFile(1);
	fListOfObjects = new TList();
	fListOfObjects->SetOwner();

	//
	// Histograms
	//  
	if(! fHistEventCounter ){
		fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",10,0,10);

		fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
		fHistEventCounter->GetXaxis()->SetBinLabel(2, "Trigger");//NotinVertexcut");
		fHistEventCounter->GetXaxis()->SetBinLabel(3, "Physics Selection"); //CINT7-B-NOPF-CENT");
		fListOfObjects->Add(fHistEventCounter);
	}

	hS0 = new TH1D("hS0",";Spherocity;count",100,0,1);
	fListOfObjects->Add(hS0);

	hRefMult = new TH1D("hRefMult",";Measured mult. (reference, |#eta|<0.8); Counts", 200, 0, 200);
	fListOfObjects->Add(hRefMult);

	hS0RefMult = new TH2D("hS0RefMult",";Spherocity;Measured mult. (reference, |#eta|<0.8)", 100,0,1, 200, 0, 200);
	fListOfObjects->Add(hS0RefMult);

	hphiso = 0;
	hphiso = new TH1D("hphiso","spherocity; #phi; counts",64,0.0,2*TMath::Pi());
	fListOfObjects->Add(hphiso);
	
	hetaso = 0;
	hetaso = new TH1D("hetaso","spherocity; #eta; counts",40,-1.0,1.0);
	fListOfObjects->Add(hetaso);
	
	hpt = new TH1D("hpt",";#it{p}_{T} (GeV/#it{c});count",nPtBins,xBins);
	fListOfObjects->Add(hpt);

	hS0pt = new TH2D("hS0pt",";Spherocity;#it{p}_{T} (GeV/#it{c})",100,0,1,nPtBins,xBins);
	fListOfObjects->Add(hS0pt);
	hS0phi= new TH2D("hS0phi",";Spherocity;#phi (rad.)",100,0,1,64,0,2.0*TMath::Pi());
	fListOfObjects->Add(hS0phi);
	hS0eta= new TH2D("hS0eta",";Spherocity;#eta",100,0,1,20,-1,1);
	fListOfObjects->Add(hS0eta);

	
	if (fAnalysisMC){
	  
	  hMultTrue = new TH1D("hMultTrue",";True mult. (reference, |#eta|<0.8); Counts", 200, 0, 200);
	  fListOfObjects->Add(hMultTrue);
	  
	  hS0True = new TH1D("hS0True",";True Spherocity;count",100,0,1);
	  fListOfObjects->Add(hS0True);
	  
	  hS0TrueMult = new TH2D("hS0TrueMult",";Spherocity;True mult. (reference, |#eta|<0.8)", 100,0,1, 200, 0, 200);
	  fListOfObjects->Add(hS0TrueMult);
	
	  hS0Truept = new TH2D("hS0Truept",";True Spherocity;#it{p}_{T} (GeV/#it{c})",100,0,1,nPtBins,xBins);
	  fListOfObjects->Add(hS0Truept);

	  hS0Truephi= new TH2D("hS0Truephi",";True Spherocity;#phi (rad.)",100,0,1,64,0,2.0*TMath::Pi());
	  fListOfObjects->Add(hS0Truephi);

	  hS0Trueeta= new TH2D("hS0Trueeta",";True Spherocity;#eta",100,0,1,20,-1,1);
	  fListOfObjects->Add(hS0Trueeta);

	}



	fEventCuts.AddQAplotsToList(fListOfObjects);
	PostData(1, fListOfObjects);

}

//______________________________________________________________________________
void AliAnalysisTaskS0::UserExec(Option_t *)
{

	// -----------------------------------------------------
	//			 InputEvent
	// -----------------------------------------------------

	AliVEvent *event = InputEvent();
	if (!event) {
		Error("UserExec", "Could not retrieve event");
		return;
	}

	// -----------------------------------------------------
	//			 E S D
	// -----------------------------------------------------

	if (fAnalysisType == "ESD"){
		fESD = dynamic_cast<AliESDEvent*>(event);

		if(!fESD){
			Printf("%s:%d ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
			this->Dump();
			return;
		}
	}

	// -----------------------------------------------------
	//			 M C
	// -----------------------------------------------------

	if (fAnalysisMC) {

		//	E S D

		if (fAnalysisType == "ESD"){
			fMC = dynamic_cast<AliMCEvent*>(MCEvent());
			if(!fMC){
				Printf("%s:%d MCEvent not found in Input Manager",(char*)__FILE__,__LINE__);
				this->Dump();
				return;
			}    

			fMCStack = fMC->Stack();

			if(!fMCStack){
				Printf("%s:%d MCStack not found in Input Manager",(char*)__FILE__,__LINE__);
				this->Dump();
				return;
			}    
		}
	}


	fHistEventCounter->Fill(0.5); 	//	All events

	AliAnalysisUtils * utils = new AliAnalysisUtils();
	if (!utils)
	{
		cout<<"------- No AnalysisUtils Object Found --------"<<utils<<endl;
		return;
	}

	// Cuts at event level
	UInt_t fSelectMask= fInputHandler->IsEventSelected();
	Bool_t isINT7selected = fSelectMask&AliVEvent::kINT7;
	if(!isINT7selected)
		return;
	fHistEventCounter->Fill(1.5);

	if (!fEventCuts.AcceptEvent(event)) {
		PostData(1, fListOfObjects);
		return;
	}

	fHistEventCounter->Fill(2.5);

	fnRefGlobal = -1;
	fnRefGlobal = AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8 );

	// Extract the event spherocity using the official cuts 
	//Double_t SOm = -1.0;
	fS0R = fSpheroUtils->GetEventShape( event, hphiso, hetaso );
	hS0->Fill(fS0R);

	hRefMult->Fill(fnRefGlobal);
	hS0RefMult->Fill(fS0R,fnRefGlobal);
	cout<<"S0  ===== "<<fS0R<<"fnRefGlobal ==  "<<fnRefGlobal<<endl;

	MakeESDAnalysis( fESD, 0.8 );


	//	Double_t S0T = -1.0;
	Int_t True08 = -1;
	if (fAnalysisMC){
		True08 = GetMultiplicityParticles(0.8);
		fS0T = fSpheroUtils->GetEventShapeTrue( fMCStack );
		hS0TrueMult->Fill(True08,fS0T);
		hMultTrue->Fill(True08);
		hS0True->Fill(fS0T);
	}

	// Post output data.
	PostData(1, fListOfObjects);

}
//_____________________________________________________________________________
void AliAnalysisTaskS0::MakeESDAnalysis( AliESDEvent *ESDevent, Double_t etaCut ){

	const Int_t nESDTracks = ESDevent->GetNumberOfTracks();

	for(Int_t iT = 0; iT < nESDTracks; iT++) {

		AliESDtrack* esdTrack = ESDevent->GetTrack(iT);

		Double_t eta      = esdTrack->Eta();
		Double_t phi      = esdTrack->Phi();
		//Double_t momentum = esdTrack->P();
		Double_t pt       = esdTrack->Pt();

		if(TMath::Abs(eta) > etaCut) continue;
		//quality cuts
		if(!fTrackFilter->IsSelected(esdTrack))  continue;
		if ( TMath::Abs(pt) < 0.1 ) continue;

		Int_t nClustersTPC = -1;
		nClustersTPC = esdTrack->GetTPCclusters(0);
		Int_t nClustersTPCShared = esdTrack->GetTPCnclsS();
		Float_t fracClustersTPCShared = -1.;
		if (nClustersTPC!=0) {
		  fracClustersTPCShared = Float_t(nClustersTPCShared)/Float_t(nClustersTPC);
		}

		hpt->Fill(pt);
		hS0pt->Fill(fS0R,pt);
		hS0phi->Fill(fS0R,phi);
		hS0eta->Fill(fS0R,eta);
	
	}

}
//________________________________________________________________________
void AliAnalysisTaskS0::MakeESDAnalysisMC(AliMCEvent* fMCEvent, Double_t etaCut)
{
        // VZ
	for ( int iT = 0 ; iT < fMCEvent->GetNumberOfTracks(); iT++ ) // loop over TRUE MC
	{
	  TParticle *mcParticle = fMCEvent->GetTrack(iT)->Particle();
	  
	  if (!mcParticle){
	    cout<<"no mcParticle"<<endl;
	    continue;
	  }

	  // VZ
	  if(!(fMCEvent->IsPhysicalPrimary(iT)))  continue;
	  
	  Double_t eta      = mcParticle->Eta();
	  Double_t phi      = mcParticle->Phi();
	  Double_t pt       = mcParticle->Pt();
	  
	  Int_t partPDG = TMath::Abs(mcParticle->GetPdgCode());
	  if(TMath::Abs(eta) > etaCut) continue;
	  //quality cuts
	  // if(!fTrackFilter->IsSelected(esdTrack))  continue;
	  if ( TMath::Abs(pt) < 0.1 ) continue;

	  hpt->Fill(pt);
	  hS0Truept->Fill(fS0T,pt);
	  hS0Truephi->Fill(fS0T,phi);
	  hS0Trueeta->Fill(fS0T,eta);
	}
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskS0::GetMultiplicityParticles(Double_t etaCut) 
{
	// Fill the special MC histogram with the MC truth info



	Int_t trackmult = 0;
	const Int_t nTracksMC = fMCStack->GetNtrack();

	for (Int_t iTracks = 0; iTracks < nTracksMC; iTracks++) {

		TParticle* trackMC = fMCStack->Particle(iTracks);
		if(!trackMC)
			continue;

		if( !(fMCStack->IsPhysicalPrimary(iTracks)) )
			continue;


		TParticlePDG* pdgPart = trackMC ->GetPDG();
		Double_t chargeMC = pdgPart->Charge();

		if( TMath::Abs(chargeMC) < 0.1 )
			continue;

		if ( TMath::Abs( trackMC->Eta() ) > etaCut )
			continue;

		trackmult++;

	}//MC track loop

	return trackmult;

}
//____________________________________________________________
void AliAnalysisTaskS0::SetTrackCuts(AliAnalysisFilter* fTrackFilter){

	AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts();
	esdTrackCuts->SetMaxFractionSharedTPCClusters(0.4);//
	esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);//
	esdTrackCuts->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);//
	esdTrackCuts->SetMaxChi2PerClusterTPC(4);//
	esdTrackCuts->SetAcceptKinkDaughters(kFALSE);//
	esdTrackCuts->SetRequireTPCRefit(kTRUE);//
	esdTrackCuts->SetRequireITSRefit(kTRUE);//
	esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
			AliESDtrackCuts::kAny);//
	esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");//
	esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);//
	esdTrackCuts->SetMaxDCAToVertexZ(2);//
	esdTrackCuts->SetDCAToVertex2D(kFALSE);//
	esdTrackCuts->SetRequireSigmaToVertex(kFALSE);//
	esdTrackCuts->SetMaxChi2PerClusterITS(36);//
	fTrackFilter->AddCuts(esdTrackCuts);

}
