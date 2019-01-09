
/*

   This macro produces: 
   - spherocity distributions as a function of multiplicity (only works for reference estimator, |eta|<0.8) using ESD information. 
   Antonio Ortiz, ICN-UNAM
   Please report bugs to: aortizve@cern.ch / antonio.ortiz@nucleares.unam.mx 
   First version: 11/09/2018

 */

#include "AliAnalysisTaskUeSpherocity.h"

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

TH1D * hSOMPerc[100];
TH2D * hSOMAux[100];
ClassImp(AliAnalysisTaskUeSpherocity)

	//_____________________________________________________________________________
	AliAnalysisTaskUeSpherocity::AliAnalysisTaskUeSpherocity():
		AliAnalysisTaskSE(),
		fESD(0x0),
		fEventCuts(0x0),
		fMC(0x0),
		fMCStack(0x0),
		fSpheroUtils(0x0),
		fTrackFilter(0x0),
		fAnalysisType("ESD"),
		fAnalysisMC(kFALSE),
		fnMultbins(0x0),
		fnRefGlobal(0x0),
		fSoBining(0x0),
		fnsoB(0x0),
		fbinSom(0x0),
		fListOfObjects(0x0),
		fHistEventCounter(0x0),
		fEvents(0x0), 
		hSOGlobal08(0x0),
		hptRefGlobal08(0x0),
		hSOTrue08(0x0),
		hRMGlobal(0x0),
		hphiso(0x0),
		hetaso(0x0),  
		fn1(0x0)
{
	// Default constructor (should not be used)
	for(Int_t i=0;i<200;++i){
		hSOGlobal[i]=0;
		hptRefGlobal[i]=0;
	}
}

//______________________________________________________________________________
AliAnalysisTaskUeSpherocity::AliAnalysisTaskUeSpherocity(const char *name):
	AliAnalysisTaskSE(name),
	fESD(0x0),
	fEventCuts(0x0),
	fMC(0x0),
	fMCStack(0x0),
	fSpheroUtils(0x0),
	fTrackFilter(0x0),
	fAnalysisType("ESD"),
	fAnalysisMC(kFALSE),
	fnMultbins(0x0),
	fnRefGlobal(0x0),
	fSoBining(0x0),
	fnsoB(0x0),
	fbinSom(0x0),
	fListOfObjects(0x0),
	fHistEventCounter(0x0), 
	fEvents(0x0), 
	hSOGlobal08(0x0),
	hptRefGlobal08(0x0),
	hSOTrue08(0x0),
	hRMGlobal(0x0),
	hphiso(0x0),
	hetaso(0x0), 
	fn1(0x0)
{
	// Output slot #1 writes into a TList
	for(Int_t i=0;i<200;++i){
		hSOGlobal[i]=0;
		hptRefGlobal[i]=0;
	}

	DefineOutput(1, TList::Class());
}

//________________________________________________________________________

void AliAnalysisTaskUeSpherocity::Exit(const char *msg) {

	Printf("%s", msg);
	return;
}


//_____________________________________________________________________________
AliAnalysisTaskUeSpherocity::~AliAnalysisTaskUeSpherocity()
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
void AliAnalysisTaskUeSpherocity::UserCreateOutputObjects()
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
	TFile * finPercent = 0;
	TString nameSoHelper = AliDataFile::GetFileName("PWGMM/spheroLHC15fpass2.root");
	finPercent = TFile::Open(nameSoHelper);
	fnsoB = fSoBining->GetNbinsX();
	fnMultbins = 100;
	for( Int_t i_mult = 0; i_mult < fnMultbins; ++i_mult ){
		hSOMPerc[i_mult] = 0;
		hSOMPerc[i_mult] = (TH1D *)finPercent->Get(Form("hSOMPerc%d",i_mult));
		hSOMAux[i_mult] = 0;
		hSOMAux[i_mult] = (TH2D *)finPercent->Get(Form("hSOMauxMult%d",i_mult));
	}


	const Int_t nPtBins      = 59;
	Double_t xBins[nPtBins+1] = {
		0.01,0.1,0.12,0.15,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
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

	fEvents = new TH1I("fEvents","Number of analyzed events; Events; Counts", 1, 0, 1);
	fListOfObjects->Add(fEvents);


	hSOGlobal08 = new TH2D("hSOGlobal08","spherocity; measured spherocity; measured mult. (reference, |#eta|<0.8)",fnMultbins, 0, fnMultbins, 1000,0.0,1.0);
	fListOfObjects->Add(hSOGlobal08);
	hptRefGlobal08 = new TH2D("hptRefGlobal08",";measured mult. (reference, |#eta|<0.8); #it{p}_{T} (GeV/#it{c})",fnMultbins,0,fnMultbins,nPtBins,xBins);
	fListOfObjects->Add(hptRefGlobal08);

	// Histos vs spherocity
	for(Int_t i=0;i<fnsoB;++i){

		hSOGlobal[i] = 0;
		hSOGlobal[i] = new TH2D(Form("hSOGlobalES%d",i),"spherocity; measured spherocity; measured mult. (reference, |#eta|<0.8)",fnMultbins, 0, fnMultbins, 1000,0.0,1.0);
		fListOfObjects->Add(hSOGlobal[i]);

		hptRefGlobal[i] = new TH2D(Form("hptRefGlobalES%d",i),";measured mult. (reference, |#eta|<0.8); #it{p}_{T} (GeV/#it{c})",fnMultbins,0,fnMultbins,nPtBins,xBins);
		fListOfObjects->Add(hptRefGlobal[i]);

	}



	if (fAnalysisMC){

		hSOTrue08 = new TH2D("hSOTrue08","spherocity; generated spherocity; generated mult. (#it{p}_{T}>0 |#eta|<0.8)",fnMultbins, 0, fnMultbins, 1000,0.0,1.0);
		fListOfObjects->Add(hSOTrue08);
		// Multiplicity response matrices
		hRMGlobal = new TH2D("hRMGlobal","; generated mult. (#it{p}_{T}>0 |#eta|<0.8); measured mult. (reference, |#eta|<0.8)",fnMultbins,0,fnMultbins,fnMultbins,0,fnMultbins);
		fListOfObjects->Add(hRMGlobal);

	}

	hphiso = 0;
	hphiso = new TH1D("hphiso","spherocity; #phi; counts",64,0.0,2*TMath::Pi());
	fListOfObjects->Add(hphiso);

	hetaso = 0;
	hetaso = new TH1D("hetaso","spherocity; #eta; counts",40,-1.0,1.0);
	fListOfObjects->Add(hetaso);


	fn1=new TH1F("fn1","fn1",5001,-1,5000);
	fListOfObjects->Add(fn1);



	fEventCuts.AddQAplotsToList(fListOfObjects);
	PostData(1, fListOfObjects);

}

//______________________________________________________________________________
void AliAnalysisTaskUeSpherocity::UserExec(Option_t *)
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

	Float_t V0MPercentile = -1;
	Float_t V0APercentile = -1;
	Float_t ADMPercentile = -1;
	Bool_t fSPDCvsTCutStatus = kFALSE;

	fSPDCvsTCutStatus = utils->IsSPDClusterVsTrackletBG(fESD);
	AliMultSelection *MultSelection = (AliMultSelection*) fESD -> FindListObject("MultSelection");
	if(MultSelection-> IsEventSelected()){
		V0MPercentile = MultSelection->GetMultiplicityPercentile("V0M",false);
		V0APercentile = MultSelection->GetMultiplicityPercentile("V0A",false);
		ADMPercentile = MultSelection->GetMultiplicityPercentile("ADM",false);
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
	Double_t SOm = -1.0;
	SOm = fSpheroUtils->GetEventShape( event, hphiso, hetaso );
	// spherocity percentil
	Double_t SotPerc = -1;
	SotPerc = GetSpheroPercentile(SOm,fnRefGlobal);
	fbinSom = -1;
	fbinSom = GetSpheroPercentileBin(SotPerc);
	// Correlation between sphericity and multiplicity 
	hSOGlobal08->Fill(fnRefGlobal,SOm);
	if(fbinSom>=0&&fbinSom<fnsoB)
		hSOGlobal[fbinSom]->Fill(fnRefGlobal,SOm);


	MakeESDAnalysis( fESD, 0.8 );


	Double_t SOt = -1.0;
	Int_t True08 = -1;
	if (fAnalysisMC){
		True08 = GetMultiplicityParticles(0.8);
		SOt = fSpheroUtils->GetEventShapeTrue( fMCStack );
		hSOTrue08->Fill(True08,SOt);
		hRMGlobal->Fill(True08,fnRefGlobal);
	}




	// Post output data.
	PostData(1, fListOfObjects);

}
//_____________________________________________________________________________
void AliAnalysisTaskUeSpherocity::MakeESDAnalysis( AliESDEvent *ESDevent, Double_t etaCut ){

	const Int_t nESDTracks = ESDevent->GetNumberOfTracks();

	for(Int_t iT = 0; iT < nESDTracks; iT++) {

		AliESDtrack* esdTrack = ESDevent->GetTrack(iT);

		if(TMath::Abs(esdTrack->Eta()) > etaCut)
			continue;

		//quality cuts
		if(!fTrackFilter->IsSelected(esdTrack))
			continue;

		Int_t nClustersTPC = -1;
		nClustersTPC = esdTrack->GetTPCclusters(0);
		Int_t nClustersTPCShared = esdTrack->GetTPCnclsS();
		Float_t fracClustersTPCShared = -1.;
		if (nClustersTPC!=0) {
			fracClustersTPCShared = Float_t(nClustersTPCShared)/Float_t(nClustersTPC);
		}

		hptRefGlobal08->Fill(fnRefGlobal,esdTrack->Pt());
		if(fbinSom>=0&&fbinSom<fnsoB)
			hptRefGlobal[fbinSom]->Fill(fnRefGlobal,esdTrack->Pt());

	}

}
//_____________________________________________________________________________
Int_t AliAnalysisTaskUeSpherocity::GetMultiplicityParticles(Double_t etaCut) 
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
Double_t AliAnalysisTaskUeSpherocity::GetSpheroPercentile( Double_t valES, Int_t valMult ){

	if(valMult<0)
		return -1.0;

	Double_t porcentaje = -1.0;

	for(Int_t bin = 1; bin <= hSOMPerc[valMult]->GetNbinsX(); ++bin){

		Double_t minEs = 0.0;
		Double_t maxEs = 0.0;
		minEs = hSOMPerc[valMult]->GetBinLowEdge(bin);
		maxEs = minEs + hSOMPerc[valMult]->GetBinWidth(bin);
		if( (valES >= minEs) && (valES < maxEs) )
			porcentaje = hSOMPerc[valMult]->GetBinContent(bin);

	}

	return porcentaje;

}
//___________________________________________________________
Int_t AliAnalysisTaskUeSpherocity::GetSpheroPercentileBin( Double_t valES ){

	Int_t BinSot = -1;

	for( Int_t i_so = 0; i_so < fnsoB; ++i_so ){
		if( valES >= fSoBining->GetBinLowEdge(i_so+1) && valES < fSoBining->GetBinLowEdge(i_so+1)+fSoBining->GetBinWidth(i_so+1) )
			BinSot = i_so;

	}

	return BinSot;
}
//____________________________________________________________
void AliAnalysisTaskUeSpherocity::SetTrackCuts(AliAnalysisFilter* fTrackFilter){

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
