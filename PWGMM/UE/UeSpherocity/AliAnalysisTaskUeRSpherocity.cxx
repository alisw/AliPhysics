
/*

   Antonio Ortiz, ICN-UNAM
   Please report bugs to: aortizve@cern.ch / antonio.ortiz@nucleares.unam.mx 
   First version: 29/07/2019
   rho as a function of R and Nch


 */

#include "AliAnalysisTaskUeRSpherocity.h"

// ROOT includes
#include <TList.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2D.h>
#include <TProfile.h>
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
#include <AliDataFile.h>

#include <iostream>
using namespace std;

const Int_t n_samples = 5;
const Char_t * NameSample[n_samples]={"REF","HM_SPD","HM_V0M","HM_SPD_JETTY","HM_V0M_JETTY"};


const Int_t nPtBins = 66;
Double_t xBins[nPtBins+1] = {
	0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
	0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
	1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
	2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
	4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
	11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
	26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0 };

const Int_t nPtLeadingBins = 6;
Double_t xPtLeading[nPtLeadingBins+1]={1.0,3.0,4.0,6.0,10.0,20.0,100};

const Double_t pi = 3.1415926535897932384626433832795028841971693993751058209749445;

ClassImp(AliAnalysisTaskUeRSpherocity)

	//_____________________________________________________________________________
	AliAnalysisTaskUeRSpherocity::AliAnalysisTaskUeRSpherocity():
		AliAnalysisTaskSE(),
		fESD(0x0),
		fdata_set("k"),
		fEventCuts(0x0),
		fTrackFilter(0x0),
		fSpheroUtils(0x0),
		fAnalysisType("ESD"),
		fListOfObjects(0x0),
		fHistEventCounter(0x0),
		hcounter(0x0),
		hphiso(0x0),
		hetaso(0x0)


{
	// Default constructor (should not be used)
	for( Int_t i_sam= 0; i_sam < n_samples; ++i_sam ){
		hR[i_sam] = 0;
		hPtLeading[i_sam] = 0;
		for(Int_t i_ptl = 0; i_ptl < nPtLeadingBins+1; ++i_ptl){

			hRDPhi[i_ptl][i_sam] = 0;
			hPtVsR[i_ptl][i_sam] = 0;

		}
	}

}

//______________________________________________________________________________
AliAnalysisTaskUeRSpherocity::AliAnalysisTaskUeRSpherocity(const char *name):
	AliAnalysisTaskSE(name),
	fESD(0x0),
	fdata_set("k"),
	fEventCuts(0x0),
	fTrackFilter(0x0),
	fSpheroUtils(0x0),
	fAnalysisType("ESD"),
	fListOfObjects(0x0),
	fHistEventCounter(0x0),
	hcounter(0x0),
	hphiso(0x0),
	hetaso(0x0)
{
	// Output slot #1 writes into a TList
        for( Int_t i_sam= 0; i_sam < n_samples; ++i_sam ){
                hR[i_sam] = 0;
                hPtLeading[i_sam] = 0;
                for(Int_t i_ptl = 0; i_ptl < nPtLeadingBins+1; ++i_ptl){
                        
                        hRDPhi[i_ptl][i_sam] = 0;
                        hPtVsR[i_ptl][i_sam] = 0;
                
                }
        }
	DefineOutput(1, TList::Class());
}

//________________________________________________________________________

void AliAnalysisTaskUeRSpherocity::Exit(const char *msg) {

	Printf("%s", msg);
	return;
}


//_____________________________________________________________________________
AliAnalysisTaskUeRSpherocity::~AliAnalysisTaskUeRSpherocity()
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
void AliAnalysisTaskUeRSpherocity::UserCreateOutputObjects()
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

	hcounter = 0;
	hcounter = new TH1D("hcounter","",12,-0.5,11.5);
	fListOfObjects->Add(hcounter);

	hphiso = 0;
	hphiso = new TH1D("hphiso","spherocity; #phi; counts",64,0.0,2*TMath::Pi());
	fListOfObjects->Add(hphiso);

	hetaso = 0;
	hetaso = new TH1D("hetaso","spherocity; #eta; counts",40,-1.0,1.0);
	fListOfObjects->Add(hetaso);

	for( Int_t i_sam= 0; i_sam < n_samples; ++i_sam ){

		hR[i_sam] = 0;
		hR[i_sam] = new TH1D(Form("hR_%s",NameSample[i_sam]),Form("%s leading h^{+}+h^{-}; #it{R}; counts",NameSample[i_sam]),50,0,5.0);
		fListOfObjects->Add(hR[i_sam]);

		hPtLeading[i_sam] = 0;
		hPtLeading[i_sam] = new TH1D(Form("hPtLeading_%s",NameSample[i_sam]),Form("%s leading h^{+}+h^{-}; #it{p}_{T}^{leading}; counts",NameSample[i_sam]),nPtLeadingBins,xPtLeading);
		fListOfObjects->Add(hPtLeading[i_sam]);

		for(Int_t i = 0; i< nPtLeadingBins+1; ++i){

			hRDPhi[i][i_sam] = 0;
			hRDPhi[i][i_sam] = new TH2D(Form("hRDPhi_%d_%s",i,NameSample[i_sam]),Form("%s leading h^{+}+h^{-}",NameSample[i_sam]),50,0,5.0,62,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);
			fListOfObjects->Add(hRDPhi[i][i_sam]);

			hPtVsR[i][i_sam] = 0;
			hPtVsR[i][i_sam] = new TH2D(Form("hPtVsR_%d_%s",i,NameSample[i_sam]),Form("%s leading h^{+}+h^{-}",NameSample[i_sam]),50,0,5.0,nPtBins,xBins);
			fListOfObjects->Add(hPtVsR[i][i_sam]);


		}
	}

	fEventCuts.AddQAplotsToList(fListOfObjects);

	PostData(1, fListOfObjects);

}

//______________________________________________________________________________
void AliAnalysisTaskUeRSpherocity::UserExec(Option_t *)
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



	fHistEventCounter->Fill(0.5); 	//	All events

	AliAnalysisUtils * utils = new AliAnalysisUtils();
	if (!utils)
	{
		cout<<"------- No AnalysisUtils Object Found --------"<<utils<<endl;
		return;
	}

	//Bool_t fSPDCvsTCutStatus = kFALSE;

	//fSPDCvsTCutStatus = utils->IsSPDClusterVsTrackletBG(fESD);

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

	// Get Leading track
	Int_t index_lider = -1;
	index_lider =IndexLeadingTrack(0.8);
	Double_t leading_pt = -1;
	if(index_lider>=0)
		leading_pt = ((AliESDtrack*)fESD->GetTrack(index_lider))->Pt(); 

	// Get multiplicity at mid-rapidity
	Int_t nchmid = -1;
	nchmid = AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8 );

	// Get multiplicity V0M
	AliMultSelection *MultSelection = (AliMultSelection*) fESD -> FindListObject("MultSelection");
	Double_t v0M = -1;
	v0M = MultSelection->GetMultiplicityPercentile("V0M",false);

	if(leading_pt<1)
		return;

	Double_t SOm = -1.0;
	if(nchmid>=51 || (v0M<0.1&&v0M>=0))
		SOm = fSpheroUtils->GetEventShape( event, hphiso, hetaso );

	MakeAnalysis(0,index_lider,0.8);// reference
	if(nchmid>=51){//HM SPD sample
		MakeAnalysis(1,index_lider,0.8);
		if( (SOm <= 0.522) && (SOm >= 0) )// jetty
			MakeAnalysis(3,index_lider,0.8);
	}
	if(v0M<0.1&&v0M>=0){
		MakeAnalysis(2,index_lider,0.8);
		if( (SOm <= 0.522) && (SOm >= 0) )// jetty
			MakeAnalysis(4,index_lider,0.8);
	}

	//cout<<"ptleading=="<<leading_pt<<" nchmid="<<nchmid<<"  V0M="<<v0M<<"   Som="<<SOm<<endl;

	// Post output data.
	PostData(1, fListOfObjects);

}

Int_t AliAnalysisTaskUeRSpherocity::IndexLeadingTrack( Double_t etaCut ){

	// selection on leading particle
	Int_t index_leading = 0;
	Double_t pt_leading    = 0;

	if(fESD->GetNumberOfTracks()==0)
		return -1;

	for(Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {

		AliESDtrack* esdTrack = fESD->GetTrack(i);

		Double_t eta      = esdTrack->Eta();
		Double_t pt       = esdTrack->Pt();

		if(TMath::Abs(eta) > etaCut)
			continue;

		//quality cuts, standard 2015 track cuts
		if(!fTrackFilter->IsSelected(esdTrack))
			continue;

		if(pt<0.15)
			continue;

		if(pt>pt_leading){
			pt_leading  = pt;
			index_leading  = i;
		}



	}// end loop over tracks
	if(pt_leading==0)
		return -1;

	return index_leading;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskUeRSpherocity::MakeAnalysis( Int_t index_sample, Int_t index_leading, Double_t etaCut ){

	Double_t phi_leading = ((AliESDtrack*)fESD->GetTrack(index_leading))->Phi();
	Double_t eta_leading = ((AliESDtrack*)fESD->GetTrack(index_leading))->Eta();
	Double_t pt_leading = ((AliESDtrack*)fESD->GetTrack(index_leading))->Pt();


	// Selecting the pTleading interval
	Int_t binPtLeading = -1;
	for( Int_t i=0; i<nPtLeadingBins; ++i ){

		if( pt_leading>=xPtLeading[i] && pt_leading<xPtLeading[i+1] ){
			binPtLeading=i+1;
		}
	}

	if(binPtLeading<=0)
		return kFALSE;

	hPtLeading[index_sample]->Fill(pt_leading);

	//cout<<"!!!!!!!!!!!!!!!!!!      ptleading="<<pt_leading<<"  bin="<<binPtLeading<<endl;


	// Fill the relevant histograms
	for(Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {

		AliESDtrack* esdTrack = fESD->GetTrack(i);

		Double_t eta      = esdTrack->Eta();
		Double_t phi      = esdTrack->Phi();
		Double_t pt       = esdTrack->Pt();

		if(i==index_leading)
			continue;

		if(TMath::Abs(eta) > etaCut)
			continue;

		//quality cuts, standard 2015 track cuts
		if(!fTrackFilter->IsSelected(esdTrack))
			continue;

		if(pt<0.15)
			continue;

		Double_t DPhi = DeltaPhi(phi,phi_leading);
		Double_t Deta = TMath::Abs(eta-eta_leading);
		Double_t R = TMath::Sqrt(DPhi*DPhi+Deta*Deta);
		hR[index_sample]->Fill(R);
		hRDPhi[0][index_sample]->Fill(R,DPhi);
		hRDPhi[binPtLeading][index_sample]->Fill(R,DPhi);
		hPtVsR[0][index_sample]->Fill(R,pt);
		hPtVsR[binPtLeading][index_sample]->Fill(R,pt);

	}// end loop over tracks

	return kTRUE;
}
//____________________________________________________________
void AliAnalysisTaskUeRSpherocity::SetTrackCuts(AliAnalysisFilter* fTrackFilter){

	AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts;
	// TPC
	esdTrackCuts->SetMinNCrossedRowsTPC(70);
	esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);

	esdTrackCuts->SetMaxChi2PerClusterTPC(4);
	esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
	esdTrackCuts->SetRequireTPCRefit(kTRUE);
	// ITS
	esdTrackCuts->SetRequireITSRefit(kTRUE);
	esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
			AliESDtrackCuts::kAny);
	// 7*(0.0015+0.0050/pt^1.1)
	esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");

	esdTrackCuts->SetMaxDCAToVertexZ(2);
	esdTrackCuts->SetDCAToVertex2D(kFALSE);
	esdTrackCuts->SetRequireSigmaToVertex(kFALSE);

	esdTrackCuts->SetMaxChi2PerClusterITS(36);

	/*
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
	 */
	fTrackFilter->AddCuts(esdTrackCuts);

}
//______________________________________________________________
Double_t AliAnalysisTaskUeRSpherocity::DeltaPhi(Double_t phia, Double_t phib,
		Double_t rangeMin, Double_t rangeMax)
{
	Double_t dphi = -999;
	Double_t pi = TMath::Pi();

	if (phia < 0)         phia += 2*pi;
	else if (phia > 2*pi) phia -= 2*pi;
	if (phib < 0)         phib += 2*pi;
	else if (phib > 2*pi) phib -= 2*pi;
	dphi = phib - phia;
	if (dphi < rangeMin)      dphi += 2*pi;
	else if (dphi > rangeMax) dphi -= 2*pi;

	return dphi;
}


