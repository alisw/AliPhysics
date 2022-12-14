/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
 // Analysis task for Mean p_T fluctuations in HM pp collisions
 
#include "Riostream.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1I.h"
#include "TROOT.h"
#include "TMath.h"
#include "THn.h"
#include "TTree.h"
#include "TList.h"
#include "TRandom.h"
#include "TObjArray.h"
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
#include "AliVEventHandler.h"
#include "THnSparse.h"
#include "AliAnalysisHMMeanPt.h"

class AliAnalysisHMMeanPt;

#include "TFile.h"
#include <iostream>
#include <fstream>

using namespace std;

ClassImp(AliAnalysisHMMeanPt)

//================================================
//--------------Constructor-----------------------
//================================================
AliAnalysisHMMeanPt::AliAnalysisHMMeanPt()
: AliAnalysisTaskSE(),
fAOD(0x0), 
fOutputList(0), 
fEventCount(0), 
fHistPhi(0), 
fHistPt(0), 
fHistEta(0), 
fCent(0), 
fHistVx(0), 
fHistVy(0), 
fHistVz(0), 
hHighMultRuns(0), 
hScale(0), 
hTracks(0), 
hTrackssq(0), 
hTrackavgpt(0), 
hTrackavgptsq(0), 
hTrackpair(0), 
hTracksumpt(0),
hTracksumptiptj(0),  
fHistPtcal(0),
fRunNumber(0), 
fTreept(0x0),
nevt(0),
fHistSignal(0),
fHistClustersTPC(0),
fHistChi2perNDF(0),
fHistTPCFoundFrac(0),
fHistTPCcrossrows(0),
sample(0),
fHistMag(0),
fTHnetaptphi(0),
fMulttrkV0(0),
fMulttrkV0_meanpt(0),
fMulttrk_pair(0),
fMulttrk_meanpt(0),
fMultV0A(0),
fMultV0C(0),
fMultV0(0),
ftrackBit(128),
vzmin(-10.),
vzmax(10.),
ptmin(0.15),
ptmax(2.0),
TPCrows(70),
htrk(0)
{
    
}
//-------------------------------------------------
AliAnalysisHMMeanPt::AliAnalysisHMMeanPt(const char* name)
:AliAnalysisTaskSE(name),
fAOD(0x0), 
fOutputList(0), 
fEventCount(0), 
fHistPhi(0), 
fHistPt(0), 
fHistEta(0), 
fCent(0), 
fHistVx(0), 
fHistVy(0), 
fHistVz(0), 
hHighMultRuns(0), 
hScale(0), 
hTracks(0), 
hTrackssq(0), 
hTrackavgpt(0), 
hTrackavgptsq(0), 
hTrackpair(0), 
hTracksumpt(0),
hTracksumptiptj(0),
fHistPtcal(0),
fRunNumber(0), 
fTreept(0x0),
nevt(0),
fHistSignal(0),
fHistClustersTPC(0),
fHistChi2perNDF(0),
fHistTPCFoundFrac(0),
fHistTPCcrossrows(0),
sample(0),
fHistMag(0),
fTHnetaptphi(0),
fMulttrkV0(0),
fMulttrkV0_meanpt(0),
fMulttrk_pair(0),
fMulttrk_meanpt(0),
fMultV0A(0),
fMultV0C(0),
fMultV0(0),
ftrackBit(128),
vzmin(-10.),
vzmax(10.),
ptmin(0.15),
ptmax(2.0),
TPCrows(70),
htrk(0)
{
  
    //constructor
    DefineOutput(1, TList::Class());
    //DefineOutput(2, TTree::Class());
    
}
//--------------Destructor----------------------------------

AliAnalysisHMMeanPt::~AliAnalysisHMMeanPt()
{
	if(fOutputList) 			    delete fOutputList;
	//if(fTreept) 			        delete fTreept;	
}
//--------------UserCreateOutputObjects-----------------------

void AliAnalysisHMMeanPt::UserCreateOutputObjects()
{
    // create output objects

    fOutputList = new TList();          
    fOutputList->SetOwner(kTRUE);       
	//===========================	
	fEventCount = new TH1I("fEventCounter", "Events", 10, 0, 10);
	fOutputList->Add(fEventCount);

	//fCent = new TH1F("fCent","Multiplicity",100,0,100);
    //fOutputList->Add(fCent);
    
	fMulttrkV0 = new TH2D("fMulttrkV0","Multiplicity; V0mult; N_{acc}",1000,0,1000,200,0,200);
    fOutputList->Add(fMulttrkV0);
    
	//fMultV0A = new TH1F("fMultV0A","V0A Multiplicity",1000,0,1000);
    //fOutputList->Add(fMultV0A);
    
	//fMultV0C = new TH1F("fMultV0C","V0C Multiplicity",1000,0,1000);
    //fOutputList->Add(fMultV0C);
    
	fMultV0 = new TH1F("fMultV0","V0A+C Multiplicity",2000,0,2000);
    fOutputList->Add(fMultV0);          
    
	fHistMag = new TH1D("MAgnetic feild", " ",5,-10.,10.);
	fOutputList->Add(fHistMag);
	
	float run_first = 280000;
	float run_last =  300000;
	int runs = (run_last - run_first) +1;
    
	hHighMultRuns = new TH1D("HighMultV0 triggered runs", ";run#; ",runs,run_first,run_last);
	hHighMultRuns ->SetBuffer(1000);
	hHighMultRuns->StatOverflows(kTRUE);
	fOutputList->Add(hHighMultRuns);
    
    fHistEta = new TH1F("hTrackEta", "fHistEta", 40, -1., 1.);
    fHistEta->GetXaxis()->SetTitle("#eta ");
    fHistEta->GetYaxis()->SetTitle("Counts");      
	fHistEta->StatOverflows(kTRUE);	fHistEta->Sumw2();      
    fOutputList->Add(fHistEta); 
    
    fHistPt = new TH1F("hTrackPt", "fHistPt", 250, 0, 5);   
    fHistPt->GetXaxis()->SetTitle("p_{T} ");
    fHistPt->GetYaxis()->SetTitle("Counts");        
	fHistPt->StatOverflows(kTRUE);	fHistPt->Sumw2();      
	fOutputList->Add(fHistPt);    
    
    fHistPhi = new TH1F("hTrackPhi","fHistPhi",180,0,TMath::TwoPi());
    fHistPhi->GetXaxis()->SetTitle("#phi ");
    fHistPhi->GetYaxis()->SetTitle("Counts");
   	fHistPhi->StatOverflows(kTRUE);	fHistPhi->Sumw2();      
    fOutputList->Add(fHistPhi);

	htrk = new TH1F("htrk","ncal", 200, 2, 202);
	htrk->StatOverflows(kTRUE);	htrk->Sumw2();
	fOutputList->Add(htrk);

	hScale = new TH2D("hScale","ncal",30, 1, 31, 140, 2, 142);
	hScale->StatOverflows(kTRUE);	hScale->Sumw2();
	fOutputList->Add(hScale);

	hTracks = new TH2D("hTracks","ncal",30, 1, 31, 140, 2, 142);
	hTracks->StatOverflows(kTRUE);	hTracks->Sumw2();
	fOutputList->Add(hTracks);

	hTrackssq = new TH2D("hTrackssq","ncal",30, 1, 31, 140, 2, 142);
	hTrackssq->StatOverflows(kTRUE);	hTrackssq->Sumw2();
	fOutputList->Add(hTrackssq);

	hTrackavgpt = new TH2D("hTrack_evt_pt"," ",30, 1, 31, 140, 2, 142);
	hTrackavgpt->StatOverflows(kTRUE);	hTrackavgpt->Sumw2();
	fOutputList->Add(hTrackavgpt);

	hTrackavgptsq = new TH2D("hTrack_evt_pt_sq"," ",30, 1, 31, 140, 2, 142);
	hTrackavgptsq->StatOverflows(kTRUE);	hTrackavgptsq->Sumw2();
	fOutputList->Add(hTrackavgptsq);

	hTrackpair = new TH2D("hTrack_npairs","",30, 1, 31, 140, 2, 142);
	hTrackpair->StatOverflows(kTRUE);	hTrackpair->Sumw2();
	fOutputList->Add(hTrackpair);

	//hTracksumpt = new TH2D("hTrack_sumpt","term2",30, 1, 31, 140, 2, 142);
	//hTracksumpt->StatOverflows(kTRUE);	hTracksumpt->Sumw2();
	//fOutputList->Add(hTracksumpt);

	hTracksumptiptj = new TH2D("hTrack_sumptiptj","term1",30, 1, 31, 140, 2, 142);
	hTracksumptiptj->StatOverflows(kTRUE);	hTracksumptiptj->Sumw2();
	fOutputList->Add(hTracksumptiptj);

	//fMulttrkV0_meanpt = new TH2D("fMulttrkV0_meanpt","; V0mult; M(p_T)",1000,0,1000,200,0,2);
    //fOutputList->Add(fMulttrkV0_meanpt);
    
	//fMulttrk_meanpt = new TH2D("fMulttrk_meanpt"," ;Multiplicity; M(p_T)",200,0,200,200,0,2);
    //fOutputList->Add(fMulttrk_meanpt);
    
	//fMulttrk_pair = new TH2D("fMulttrk_pair",";Multiplicity; N_{pairs}",200,0,200,1000,0,20000);
    //fOutputList->Add(fMulttrk_pair);  
    
	//--------
	//fHistChi2perNDF = new TH1F("fHistChi2perNDF","#chi^2/ndf",50,0,10);
	//fHistChi2perNDF->StatOverflows(kTRUE);
	//fOutputList->Add(fHistChi2perNDF);

	fHistClustersTPC = new TH1F("fHistClustersTPC","Ncls",160,0,160);
	fHistClustersTPC->StatOverflows(kTRUE);	
	fOutputList->Add(fHistClustersTPC);

	fHistTPCFoundFrac = new TH1F("fHistTPCFoundFrac","Ncls/NCrossRows",120,0,1.2);
	fHistTPCFoundFrac->StatOverflows(kTRUE);	
	fOutputList->Add(fHistTPCFoundFrac);

	fHistTPCcrossrows = new TH1F("fHistTPCcrossrows","NCrossRows",160,0,160);
	fHistTPCcrossrows->StatOverflows(kTRUE);	
	fOutputList->Add(fHistTPCcrossrows);

    fEventCuts.AddQAplotsToList(fOutputList);
    
	//Vertex distributions
    fHistVx = new TH1D("fHistVx","Primary vertex distribution - x coordinate;V_{x} (cm)",200,-20, 20);
    fOutputList->Add(fHistVx);
    
    fHistVy = new TH1D("fHistVy","Primary vertex distribution - y coordinate;V_{y} (cm)",200,-20, 20);
    fOutputList->Add(fHistVy);
    
    fHistVz = new TH1D("fHistVz","Primary vertex distribution - z coordinate;V_{z} (cm)",200,-20,20);
    fOutputList->Add(fHistVz);

    /*fTreept = new TTree(Nameoftree,"making"); 
    fTreept->Branch("nch",&nch,"nch/I");
    fTreept->Branch("sample",&sample,"sample/I");     
    fTreept->Branch("charge",&charge,"charge[nch]/I");
    fTreept->Branch("pt", &pt,"pt[nch]/F");*/

    PostData(1, fOutputList);         
    //PostData(2, fTreept);         

}
//----------------------UserExec------------------------------

void AliAnalysisHMMeanPt::UserExec(Option_t *){
    
//cout<<"starting analysis"<<endl;
    
	fEventCount -> Fill("Event captured", 1);

	fAOD = dynamic_cast<AliAODEvent*>(InputEvent());   
	if(!fAOD) return;  

	fEventCount ->Fill("AOD check",1);

	//===magnetic field====
	float magfld = fAOD->GetMagneticField();
	fHistMag->Fill(magfld);

	field = fAOD->GetMagneticField();

   	//Physics selection & trigger
	UInt_t fSelectMask= fInputHandler->IsEventSelected();
   	Bool_t isINT7selected = fSelectMask& AliVEvent::kHighMultV0; 
	if(!isINT7selected) return;
	//if(!fInputHandler->IsEventSelected()) return;

	fEventCount ->Fill("Physsel check",1);

	AliAODVertex *fPrimaryVtx=(AliAODVertex*)fInputEvent->GetPrimaryVertex();
	if (!fPrimaryVtx) return;    
	fEventCount ->Fill("PrimVtx check",1);
    
    if(fPrimaryVtx->GetNContributors() < 1) return;
    double xv=fPrimaryVtx->GetX();
    double yv=fPrimaryVtx->GetY();
    double zv=fPrimaryVtx->GetZ();
    
    fHistVx->Fill(xv);
    fHistVy->Fill(yv);
    fHistVz->Fill(zv);

    	// vertex cut    
    	//if (TMath::Abs(zv) > 10.) return;   	// vrtxz from Addtask     
    	if (zv< vzmin || zv > vzmax) return;   	// vrtxz from Addtask     
		//cout << vzmin <<"\t"<< vzmax <<endl;

	fEventCount ->Fill("Vzcut check",1);
    
	//Centrality    
	float centrality(0);
    	AliMultSelection *multSelection =static_cast<AliMultSelection*>(fAOD->FindListObject("MultSelection"));
    	if(!multSelection) return;
    	if(multSelection) centrality = multSelection->GetMultiplicityPercentile("V0M");
	fEventCount ->Fill("centrality check",1);
   
	//***************VZERO-AMPLITUDE************************	
	AliAODVZERO *aodV0 = fAOD->GetVZEROData();
	float fV0A_mult = aodV0->GetMTotV0A();
	float fV0C_mult = aodV0->GetMTotV0C();      
	fV0_total = fV0A_mult + fV0C_mult;
	
	fEventCount ->Fill("V0 multiplicity check",1);	
   	//--------------------------------------------------------
   
    fRunNumber = fInputEvent->GetRunNumber();
    hHighMultRuns->Fill(fRunNumber);

    fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kHighMultV0);

    if (!fEventCuts.AcceptEvent(fInputEvent)) return;

    
    //fCent	->Fill(centrality);    
	fMultV0 ->Fill(fV0_total);
	//fMultV0A->Fill(fV0A_mult);
	//fMultV0C->Fill(fV0C_mult);	
	
	fEventCount ->Fill("AliVEvent check",1);
	
//----------end of detector cuts---------        
    int trk =0;	
    double sum_trackpt = 0.;
    double sum_ptiptj =0.; 

    sample = gRandom->Integer(30)+1 ;

    Int_t iTracks(fAOD->GetNumberOfTracks());      
         
    for(Int_t i(0); i < iTracks; i++) {                 
    
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         
        if(!AcceptTrack(track)) continue;

		fHistClustersTPC->Fill(track->GetTPCNcls());
		//fHistChi2perNDF->Fill(track->Chi2perNDF());
		fHistTPCcrossrows->Fill(track->GetTPCCrossedRows());
		fHistTPCFoundFrac->Fill(track->GetTPCFoundFraction());

		fHistEta->Fill(track->Eta());        
        fHistPt->Fill(track->Pt()); 
       	fHistPhi->Fill(track->Phi()); 

		//TPCrows[trk] = track->GetTPCCrossedRows();
		//pt[trk] 	= track->Pt();
		//charge[trk] = track->Charge();

		sum_trackpt = sum_trackpt + track->Pt() ;
		sum_ptiptj  = sum_ptiptj + ( track->Pt() * track->Pt() ) ;
			
		/* ALICE method...
			int trkall =0;		
			for(int j =i+1; j < iTracks; j++){
			if(i == j) continue;
		
			AliAODTrack* track1 = static_cast<AliAODTrack*>(fAOD->GetTrack(j));         
        		if(!AcceptTrack(track1)) continue;
        		trkall++;
        	
			sum_ptiptj = sum_ptiptj + ( track->Pt() * track1->Pt() ) ;
		
			} //track_2
		*/

		trk++;               
		
		}//track_1
	//===tree=========
	if(trk < 2) return;
	
	//nch = trk;
	
	//fTreept->Fill();
			
	double trksq = trk * trk;
    double meanpt = sum_trackpt / trk;
	double meanptsq = meanpt * meanpt ;

	double npair = trk*(trk-1);  
	if(npair < 2) return;  
	
	//cout << sum_trackpt << "\t"<< trk << "\t"<< meanpt << "\t"<< npair << "\t" ;
	//cout << "npair\n";
    
    //----filling histograms-----
    htrk			->Fill(trk);
	hScale			->Fill(sample, trk);
	hTracks			->Fill(sample, trk, trk);
	hTrackssq		->Fill(sample, trk, trksq);
	hTrackavgpt		->Fill(sample, trk, meanpt);
	hTrackavgptsq	->Fill(sample, trk, meanptsq);
	hTrackpair		->Fill(sample, trk, npair);
	//++++++++++++++++++++++++++++++

	double term = sum_trackpt * sum_trackpt - sum_ptiptj;
	//hTracksumpt ->Fill(sample, trk, term);
	term = term/npair;

	hTracksumptiptj ->Fill(sample, trk, term);
	//=============extra plots====
	fMulttrkV0 ->Fill(fV0_total, trk);

	//fMulttrk_pair ->Fill(trk, npair);	
	//fMulttrkV0_meanpt ->Fill(fV0_total, meanpt);	
	//fMulttrk_meanpt ->Fill(trk, meanpt);		
	//==============    
	//nevt++;
	fEventCount ->Fill("Events Analyzed",1);   
	//cout << ftrackBit <<"\t"<< vzmax <<"\t"<< Nameoftree<< endl;

	PostData(1, fOutputList);
	//PostData(2, fTreept);                           
                                                    
}
//---------------------------------------------------------------------------------------

bool AliAnalysisHMMeanPt::AcceptTrack(AliAODTrack* aodtrack) const{
if(!aodtrack) return kFALSE;
Double_t pt = aodtrack->Pt();
    
if(pt < ptmin || pt> ptmax) return kFALSE;
if( aodtrack->Charge() == 0 ) return kFALSE;
if(!aodtrack->TestFilterBit(ftrackBit)) return kFALSE;  // Track Bit
if(aodtrack->Eta() < -0.8 || aodtrack->Eta() > 0.8) return kFALSE;
if(aodtrack->GetTPCCrossedRows() < TPCrows) return kFALSE;

return kTRUE;
}
//--------------------------------------------------------------------------------

void AliAnalysisHMMeanPt::Terminate(Option_t *){
// terminate
// called at the END of the analysis (when all events are processed)
Info("AliAnalysisHMMeanPt"," Task Successfully finished");
AliInfo(Form("Found  %d Analysed events",nevt));

}


