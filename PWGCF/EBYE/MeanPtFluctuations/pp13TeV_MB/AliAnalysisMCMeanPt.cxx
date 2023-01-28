#include "Riostream.h"
#include "TChain.h"
#include "TROOT.h"
#include "TMath.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliVEvent.h"
#include "AliAODHeader.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliInputEventHandler.h"
#include "AliHeader.h"
#include "AliMultSelection.h"
#include "AliMultSelectionBase.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCVertex.h"
#include "AliStack.h"
#include "AliLog.h"
#include "AliAnalysisUtils.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliVEventHandler.h"
#include "AliAnalysisMCMeanPt.h"
#include "TFile.h"
#include <iostream>
#include <fstream>

//class AliAnalysisMCMeanPt;
using namespace std;

ClassImp(AliAnalysisMCMeanPt)
//================================================
//--------------Constructor-----------------------
//================================================
AliAnalysisMCMeanPt::AliAnalysisMCMeanPt()
:AliAnalysisTaskSE(),
fAOD(0x0),
fOutputList(0),
fEventCount(0),
fCentMC(0),
fHistVx(0),
fHistVy(0),
fHistVz(0),
fHistVzcut(0),
AODMCTrackArray(0),
fHistPhi_gen(0),
fHistPt_gen(0),
fHistEta_gen(0),
fHistPtgen(0),
hScale_gen(0),
hTracks_gen(0),
hTracksq_gen(0),
hTrackavgpt_gen(0),
hTrackavgptsq_gen(0),
hTrackpair_gen(0),
fHistPhi_rec(0),
fHistPt_rec(0),
fHistEta_rec(0),
fHistPtrec(0),
hScale_rec(0),
hTracks_rec(0),
hTracksq_rec(0),
hTrackavgpt_rec(0),
hTrackavgptsq_rec(0),
hTrackpair_rec(0),
fTHnfcentetaptphi(0),
fTHnfcentetaptphigen(0),
fTreeFB768rec(0x0),
fTreeFB768gen(0x0),
nevt(0),
hTracksumpt_rec(0),
hTracksumptiptj_rec(0),
hTracksumpt_gen(0),
hTracksumptiptj_gen(0),
isample(0),
hrec(0),
hprim(0),
hgen(0),
hprimgen(0),
hsec(0),
hsec_wd(0),
hsec_mat(0),
fTHnfnch(0),
fTHnfnchP(0),
fzvtxcutmin(0),
fzvtxcutmax(10),
ftrackBit(128),
fCutTPCMaxCls(100.),
fCutTPCNCls(80.)
{
    // default constructor, don't allocate memory here!
}
//_____________________________________________________________________________
AliAnalysisMCMeanPt::AliAnalysisMCMeanPt(const char* name)
:AliAnalysisTaskSE(name),
fAOD(0x0),
fOutputList(0),
fEventCount(0),
fCentMC(0),
fHistVx(0),
fHistVy(0),
fHistVz(0),
fHistVzcut(0),
AODMCTrackArray(0),
fHistPhi_gen(0),
fHistPt_gen(0),
fHistEta_gen(0),
fHistPtgen(0),
hScale_gen(0),
hTracks_gen(0),
hTracksq_gen(0),
hTrackavgpt_gen(0),
hTrackavgptsq_gen(0),
hTrackpair_gen(0),
fHistPhi_rec(0),
fHistPt_rec(0),
fHistEta_rec(0),
fHistPtrec(0),
hScale_rec(0),
hTracks_rec(0),
hTracksq_rec(0),
hTrackavgpt_rec(0),
hTrackavgptsq_rec(0),
hTrackpair_rec(0),
fTHnfcentetaptphi(0),
fTHnfcentetaptphigen(0),
fTreeFB768rec(0x0),
fTreeFB768gen(0x0),
nevt(0),
hTracksumpt_rec(0),
hTracksumptiptj_rec(0),
hTracksumpt_gen(0),
hTracksumptiptj_gen(0),
isample(0),
hrec(0),
hprim(0),
hgen(0),
hprimgen(0),
hsec(0),
hsec_wd(0),
hsec_mat(0),
fTHnfnch(0),
fTHnfnchP(0),
fzvtxcutmin(0),
fzvtxcutmax(10),
ftrackBit(128),
fCutTPCMaxCls(100.),
fCutTPCNCls(80.)
{
    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());

}
//==========================================================
//--------------Destructor----------------------------------
//==========================================================
AliAnalysisMCMeanPt::~AliAnalysisMCMeanPt()
{
    //destructor
    if(fOutputList) 			        delete fOutputList;
    if (fCentMC)               		    delete fCentMC;
    if (fEventCount)               		delete fEventCount;

}
//============================================================
//--------------UserCreateOutputObjects-----------------------
//============================================================
void AliAnalysisMCMeanPt::UserCreateOutputObjects()
{

    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);

	fEventCount = new TH1I("fEventCounter", "Events", 10, 0, 10);
	fOutputList->Add(fEventCount);

	const int dim = 6;
	int bins[dim] = {142, 142, 142, 142, 142, 142};
	double mins[dim] = {0., 0., 0., 0., 0., 0.};
	double maxs[dim] = {142., 142., 142., 142., 142., 142.};

	fTHnfnch = new THnSparseD("fTHnfnch","fTHnfnch",dim, bins, mins, maxs);
	fTHnfnch ->GetAxis(0)->SetTitle("N_{truth}^{pt, pri}");
	fTHnfnch ->GetAxis(1)->SetTitle("N_{truth}^{pt, sec}");
	fTHnfnch ->GetAxis(2)->SetTitle("N_{truth}^{pri}");
	fTHnfnch ->GetAxis(3)->SetTitle("N_{truth}^{sec}");
	fTHnfnch ->GetAxis(4)->SetTitle("N_{rec}^{pri}");
	fTHnfnch ->GetAxis(5)->SetTitle("N_{rec}^{sec}");
	fOutputList->Add(fTHnfnch);

    //fEventCuts.AddQAplotsToList(fOutputList);

	fCentMC = new TH1F("fCent","Multiplicity",100,0,100);
    fOutputList->Add(fCentMC);

	//Vertex distributions
    fHistVz = new TH1D("fHistVz","Primary vertex - z (Before cut);V_{z} (cm)",400,-20,20);
    fOutputList->Add(fHistVz);

	fHistVzcut = new TH1D("fHistVzcut","Primary vertex - z (After cut);V_{z} (cm)",400,-20,20);
    fOutputList->Add(fHistVzcut);

	hrec = new TH2D("hrec","Physical tracks", 60, 0., 3., 40, -1., 1.);
	hrec->StatOverflows(kTRUE);
	fOutputList->Add(hrec);

	hprim = new TH2D("hprim","Physical tracks", 60, 0., 3., 40, -1., 1.);
	hprim->StatOverflows(kTRUE);
	fOutputList->Add(hprim);

	hgen = new TH2D("hgen","Physical tracks", 60, 0., 3., 40, -1., 1.);
	hgen->StatOverflows(kTRUE);
	fOutputList->Add(hgen);

	hprimgen = new TH2D("hprimgen","Physical tracks", 60, 0., 3., 40, -1., 1.);
	hprimgen->StatOverflows(kTRUE);
	fOutputList->Add(hprimgen);

	//========gen=============
	hScale_gen = new TH2D("hScale_gen","ncal",31, 1, 31, 140, 2, 142);
	hScale_gen->StatOverflows(kTRUE);
	fOutputList->Add(hScale_gen);

	hTracks_gen = new TH2D("hTracks_gen","ncal",31, 1, 31, 140, 2, 142);
	hTracks_gen->StatOverflows(kTRUE);
	fOutputList->Add(hTracks_gen);

	hTracksq_gen = new TH2D("hTracks_sq_gen","",31, 1, 31, 140, 2, 142);
	hTracksq_gen->StatOverflows(kTRUE);
	fOutputList->Add(hTracksq_gen);

	hTrackavgpt_gen = new TH2D("hTrack_evt_pt_gen"," ",31, 1, 31, 140, 2, 142);
	hTrackavgpt_gen->StatOverflows(kTRUE);
	fOutputList->Add(hTrackavgpt_gen);

	hTrackavgptsq_gen = new TH2D("hTrack_evt_pt_sq_gen"," ",31, 1, 31, 140, 2, 142);
	hTrackavgptsq_gen->StatOverflows(kTRUE);
	fOutputList->Add(hTrackavgptsq_gen);

	hTrackpair_gen = new TH2D("hTrack_npairs_gen"," ",31, 1, 31, 140, 2, 142);
	hTrackpair_gen->StatOverflows(kTRUE);
	fOutputList->Add(hTrackpair_gen);

	hTracksumptiptj_gen = new TH2D("hTrack_sumptiptj_gen","term1",31, 1, 31, 140, 2, 142);
	hTracksumptiptj_gen->StatOverflows(kTRUE);
	fOutputList->Add(hTracksumptiptj_gen);

	//=========rec==================

	hScale_rec = new TH2D("hScale_rec","ncal",31, 1, 31, 140, 2, 142);
	hScale_rec->StatOverflows(kTRUE);
	fOutputList->Add(hScale_rec);

	hTracks_rec = new TH2D("hTracks_rec","ncal",31, 1, 31, 140, 2, 142);
	hTracks_rec->StatOverflows(kTRUE);
	fOutputList->Add(hTracks_rec);

	hTracksq_rec = new TH2D("hTracks_sq_rec","",31, 1, 31, 140, 2, 142);
	hTracksq_rec->StatOverflows(kTRUE);
	fOutputList->Add(hTracksq_rec);

	hTrackavgpt_rec = new TH2D("hTrack_evt_pt_rec"," ",31, 1, 31, 140, 2, 142);
	hTrackavgpt_rec->StatOverflows(kTRUE);
	fOutputList->Add(hTrackavgpt_rec);

	hTrackavgptsq_rec = new TH2D("hTrack_evt_pt_sq_rec"," ",31, 1, 31, 140, 2, 142);
	hTrackavgptsq_rec->StatOverflows(kTRUE);
	fOutputList->Add(hTrackavgptsq_rec);

	hTrackpair_rec = new TH2D("hTrack_npairs_rec"," ",31, 1, 31, 140, 2, 142);
	hTrackpair_rec->StatOverflows(kTRUE);
	fOutputList->Add(hTrackpair_rec);

	hTracksumptiptj_rec = new TH2D("hTrack_sumptiptj_rec","term1",31, 1, 31, 140, 2, 142);
	hTracksumptiptj_rec->StatOverflows(kTRUE);
	fOutputList->Add(hTracksumptiptj_rec);

	//------------------------------

    PostData(1, fOutputList);
}
//============================================================
//----------------------UserExec------------------------------
//============================================================
void AliAnalysisMCMeanPt::UserExec(Option_t *){

    fEventCount->Fill("before cuts",1);
	
	AliVEvent *event = InputEvent();
    if (!event){
   		cout<<"NO EVENT FOUND!"<<endl;
   		return;
   	}

	fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  	if (!fAOD){
    	cout<< "Error: AOD not found " <<endl;
    	return;
    }
	fEventCount -> Fill("AOD check",1);

	AliAODInputHandler *eventHandler = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  	if (!eventHandler) {
    	cout<<"Input handler not found."<<endl;
    	return;
  	}

	//Physics selection
    UInt_t fSelectMask= eventHandler->IsEventSelected();
    Bool_t isINT7selected = fSelectMask& AliVEvent::kINT7;
	//from https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGtoolsEventProp
    if(!isINT7selected)return ;
	fEventCount->Fill("Physics Selection",1);

	int fCentrality = 0. ;
	AliMultSelection *MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    if(!MultSelection) {
    	cout << "AliMultSelection object not found!" << endl;
    	return;
    }
    else{
    	fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
		fCentMC->Fill(fCentrality);
    	}
	fEventCount->Fill("centrality check",1);

  	const AliAODVertex *vertex = (AliAODVertex*)fInputEvent->GetPrimaryVertex();//fAOD->GetPrimaryVertex();
    if(!vertex) return;

	fEventCount->Fill("PrimVtx check",1);

   	if(vertex->GetNContributors() < 1) return;
	double vz = vertex->GetZ();
    fHistVz->Fill(vertex->GetZ());

    if (vz < fzvtxcutmin || vz > fzvtxcutmax) return;        // vertex cut
	fEventCount->Fill("Vz-cut check",1);

	fHistVzcut->Fill(vz);

  	AODMCTrackArray = (TClonesArray*) fAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  	if (AODMCTrackArray == NULL) {
   		cout<<"Error: MC particles branch not found!\n";
   		return;
  	}

   	AliMCEvent* mcEvent = MCEvent();
    if (!mcEvent) {
    	cout<< "Could not retrieve MC event"<<endl;
    	return;
    }

	//oooo00000000ooooo000000000ooooooooooo
	// Reconstructed
	//_____________________________________
  	int nTracks = fAOD->GetNumberOfTracks();
	//printf("%8d\t\t",nTracks);
	
    int trkrec =0;
	double sum_trackpt_rec = 0.;
	double sum_ptiptj_rec =0.;
	double fillnch[6] = {0};

    //fRunNumber = fAOD->GetRunNumber();// cout << fRunNumber  <<endl;

	isample = gRandom->Integer(30)+1;
	int nsrec =0, nsrec_sec =0;
	int i;

    for(i = 0; i < nTracks; i++)
	{

    	AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
    	if(!track) continue;
        if(!AcceptTrack(track)) continue;

		AliAODMCParticle* trackrec = (AliAODMCParticle*)AODMCTrackArray->At(TMath::Abs(track->GetLabel()));
		if(!trackrec) continue;

      	double eTa = track->Eta();
      	double pT = track->Pt();

        if(pT< 0.15 || pT > 2.0) continue;		// not included in track selection

		hrec -> Fill(pT, eTa); // All rec.- including sec.
		//if(trackrec->IsSecondaryFromWeakDecay()) hsec_wd ->Fill(pT, eTa);
		//if(trackrec->IsSecondaryFromMaterial()) hsec_mat ->Fill(pT, eTa);
		//if(trackrec->IsSecondaryFromWeakDecay() || track->IsSecondaryFromMaterial()) hsec ->Fill(pT, eTa);

		nsrec_sec++;

		if(!trackrec->IsPhysicalPrimary()) continue;
		hprim -> Fill(pT, eTa);

		sum_trackpt_rec = sum_trackpt_rec + pT ;
		sum_ptiptj_rec = sum_ptiptj_rec + ( pT * pT ) ;

        trkrec++;
	}	//track1

	fillnch[4] = trkrec;
	fillnch[5] = nsrec_sec;

	if(trkrec < 2) return;
	
    double trkrecsq = trkrec * trkrec;
    double meanptrec = sum_trackpt_rec / trkrec;
	double meanptrecsq = meanptrec * meanptrec ;
	double npairrec = trkrec*(trkrec-1);

	if(npairrec < 2.) return;
    //----filling histograms-----
	hScale_rec			->Fill(isample, trkrec);
	hTracks_rec			->Fill(isample, trkrec, trkrec);
	hTracksq_rec		->Fill(isample, trkrec, trkrecsq);
	hTrackavgpt_rec		->Fill(isample, trkrec, meanptrec);
	hTrackavgptsq_rec	->Fill(isample, trkrec, meanptrecsq);
	hTrackpair_rec		->Fill(isample, trkrec, npairrec);
	//++++++++++++++++++++++++++++++
	double term = sum_trackpt_rec * sum_trackpt_rec - sum_ptiptj_rec; 
	hTracksumptiptj_rec ->Fill(isample, trkrec, (term/npairrec));

	//----gen-----------------
  	int mcTracks = AODMCTrackArray->GetEntries();
	//printf("%8d\n",mcTracks);

    int trkgen =0;
	double sum_trackpt_gen = 0.;
	double sum_ptiptj_gen =0.;

	int nsgen_all=0, nsgen_pt =0;
	int nch_all =0, nch_pt =0;
	int iMC;

    for(iMC = 0; iMC < mcTracks; iMC++)
	{

    	AliAODMCParticle* mcparticle = (AliAODMCParticle*)AODMCTrackArray->At(iMC);
      	if(!mcparticle) continue;

        if(mcparticle->Eta() < -0.8 || mcparticle->Eta() > 0.8) continue;
		if(mcparticle->Charge() == 0) continue;
		double mceTa = mcparticle -> Eta();
		double mcpT  = mcparticle -> Pt();
		nch_all++; // no pt cut + no primary condition

        if( mcpT >= 0.15 &&  mcpT <= 2.0)  nch_pt++;	// pt cut + no primary
		hgen -> Fill(mcpT, mceTa);

		if(!mcparticle->IsPhysicalPrimary()) continue;
		nsgen_all++; // no pt cut + primary

		if( mcpT < 0.15 ||  mcpT > 2.0) continue;

		hprimgen -> Fill(mcpT, mceTa);

		sum_trackpt_gen = sum_trackpt_gen + mcpT ;
		sum_ptiptj_gen = sum_ptiptj_gen + ( mcpT * mcpT) ;

        trkgen++;// pt cut + primary condition

	} //particle1

	fillnch[0] = trkgen;
	fillnch[1] = nch_pt;
	fillnch[2] = nsgen_all;
	fillnch[3] = nch_all;
	fTHnfnch->Fill(fillnch);

	if(trkgen < 2) return;

    double trkgensq 	= trkgen * trkgen  ;
    double meanptgen 	= sum_trackpt_gen / trkgen;
	double meanptgensq 	= meanptgen*meanptgen ;
	double npairgen 	= trkgen*(trkgen-1);

	if(npairgen < 2.) return;

    	//----filling histograms-----
	hScale_gen			->Fill(isample, trkgen);
	hTracks_gen			->Fill(isample, trkgen, trkgen);
	hTracksq_gen		->Fill(isample, trkgen, trkgensq);
	hTrackavgpt_gen		->Fill(isample, trkgen, meanptgen);
	hTrackavgptsq_gen	->Fill(isample, trkgen, meanptgensq);
	hTrackpair_gen		->Fill(isample, trkgen, npairgen);
	//++++++++++++++++++++++++++++++
	double term_gen 	= sum_trackpt_gen * sum_trackpt_gen - sum_ptiptj_gen ;
	hTracksumptiptj_gen ->Fill(isample, trkgen, (term_gen/npairgen));

	//--------------------------------------------------------
	fEventCount ->Fill("Events Analyzed",1);
    PostData(1, fOutputList);
}
//---------------------------------------------------------------------------------------
Bool_t AliAnalysisMCMeanPt::AcceptTrack(AliAODTrack* aodtrack) const{

	if(!aodtrack) return kFALSE;
	double pt = aodtrack->Pt();

	if(pt< 0.15 ) return kFALSE;
	if( aodtrack->Charge() == 0 ) return kFALSE;
	if(!aodtrack->TestFilterBit(ftrackBit)) return kFALSE;   
	if(aodtrack->Eta() < -0.8 || aodtrack->Eta() > 0.8) return kFALSE;
	if(aodtrack->GetTPCCrossedRows() < fCutTPCMaxCls) return kFALSE; 

	return kTRUE;
}
//------------------------------------------------------------------------------------------
void AliAnalysisMCMeanPt::Terminate(Option_t *)
{
    // terminate
    Info("AliAnalysisMCMeanPt"," Task Successfully finished");
    //AliInfo(Form("Found  %d MC events",nevt));
}
//____________________________________________________________________________________________
