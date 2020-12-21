#include "AliAnalysisTaskChargedFlow.h"

// ROOT includes
#include <TList.h>
#include <TFile.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TMatrixDSym.h>
#include <TClonesArray.h>
#include <TRandom3.h>
#include <TComplex.h>

// AliRoot includes
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisFilter.h"
#include "AliAnalysisUtils.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliESDtrackCuts.h"
#include "AliAODITSsaTrackCuts.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliCentrality.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliVVertex.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliAODMCParticle.h"
#include "AliExternalTrackParam.h"
#include "AliMultSelection.h"
#include "AliInputEventHandler.h"

// STL includes
#include <iostream>
using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskChargedFlow)
//___________________________________________________________________________
AliAnalysisTaskChargedFlow::AliAnalysisTaskChargedFlow():
	AliAnalysisTaskSE(),
	fEtaCut(0.8),
	fVtxCut(10.0),  
	fMinPt(0.2),
	fMaxPt(3.0),
	fSample(10),

	fOutputContainer(0),

	hMult(0), 
	hMultTrigger(0), 
	fVtxAfterCuts(0),

	fPhiDis(0), 
	fEtaDis(0),
	fEtaBefore(0),
	fPtDis(0),
	fPtBefore(0),
	hIsPiKp(0),

	fChsc4242(0),
	fChsc4242Gap0(0),
	fChsc4242_3sub(0),
	fChsc4224_3sub(0),
	fChsc4242_3sub_perm2(0),
	fChsc4224_3sub_perm2(0),
	fChsc4242_3sub_perm3(0),
	fChsc4224_3sub_perm3(0),
	fChsc3232(0),
	fChsc3232Gap0(0),
	fChsc3232_3sub(0),
	fChsc3223_3sub(0),
	fChsc3232_3sub_perm2(0),
	fChsc3223_3sub_perm2(0),
	fChsc3232_3sub_perm3(0),
	fChsc3223_3sub_perm3(0)
{
 
	for(int h=0; h<3; h++)
	{
 
		fChcn2Ntrks1bin[h] 									= 0;
		fChcn2Gap0Ntrks1bin[h] 							= 0;
		fChcn2Gap2Ntrks1bin[h] 							= 0;
		fChcn2Gap4Ntrks1bin[h] 							= 0;
		fChcn2Gap8Ntrks1bin[h] 							= 0;
		fChcn2GapNtrks1bin[h] 							= 0;
		fChcn2Gap14Ntrks1bin[h] 						= 0;

		fChcn2_3subLMNtrks1bin[h] 					= 0;
		fChcn2_3subRMNtrks1bin[h] 					= 0;
		fChcn2_3subLM_perm2Ntrks1bin[h] 		= 0;
		fChcn2_3subLR_perm2Ntrks1bin[h] 		= 0;
		fChcn2_3subRM_perm3Ntrks1bin[h] 		= 0;
		fChcn2_3subRL_perm3Ntrks1bin[h] 		= 0;

		fChcn4Ntrks1bin[h] 									= 0;
		fChcn4Gap0Ntrks1bin[h] 							= 0;
		fChcn4_3subNtrks1bin[h] 						= 0;
		fChcn4_3sub_perm2Ntrks1bin[h]				= 0;
		fChcn4_3sub_perm3Ntrks1bin[h]				= 0;

		fChcn6Ntrks1bin[h] 									= 0;
		fChcn6Gap0Ntrks1bin[h] 							= 0;

		fChcn8Ntrks1bin[h] 									= 0;
		fChcn8Gap0Ntrks1bin[h] 							= 0;

 	}

	for(int k=0; k<fSample+1; k++)
	{			

		fsc4242[k] 													= 0;
		fsc4242Gap0[k] 											= 0;
		fsc4242_3sub[k] 										= 0;
		fsc4224_3sub[k] 										= 0;
		fsc4242_3sub_perm2[k] 						= 0;
		fsc4224_3sub_perm2[k] 						= 0;
		fsc4242_3sub_perm3[k]							= 0;
		fsc4224_3sub_perm3[k]							= 0;

		fsc3232[k] 													= 0;
		fsc3232Gap0[k] 											= 0;
		fsc3232_3sub[k] 										= 0;
		fsc3223_3sub[k] 										= 0;
		fsc3232_3sub_perm2[k] 						= 0;
		fsc3223_3sub_perm2[k] 						= 0;
		fsc3232_3sub_perm3[k] 						= 0;
		fsc3223_3sub_perm3[k] 						= 0;

		for(Int_t h = 0; h < 3; h++)	
		{		

			fcn2Ntrks1bin[h][k] 							= 0;
			fcn2GapNtrks1bin[h][k] 						= 0;
			fcn2Gap0Ntrks1bin[h][k] 					= 0;
			fcn2Gap2Ntrks1bin[h][k] 					= 0;
			fcn2Gap4Ntrks1bin[h][k] 					= 0;
			fcn2Gap8Ntrks1bin[h][k] 					= 0;
			fcn2Gap14Ntrks1bin[h][k] 					= 0;

			fcn2_3subLMNtrks1bin[h][k] 				= 0;
			fcn2_3subRMNtrks1bin[h][k] 				= 0;
			fcn2_3subLM_perm2Ntrks1bin[h][k] = 0;
			fcn2_3subLR_perm2Ntrks1bin[h][k] = 0;
			fcn2_3subRM_perm3Ntrks1bin[h][k] = 0;
			fcn2_3subRL_perm3Ntrks1bin[h][k] = 0;

			fcn4Ntrks1bin[h][k] 							= 0;
			fcn4Gap0Ntrks1bin[h][k] 					= 0;
			fcn4_3subNtrks1bin[h][k] 					= 0;	
			fcn4_3sub_perm2Ntrks1bin[h][k]		= 0;
			fcn4_3sub_perm3Ntrks1bin[h][k]		= 0;

			fcn6Ntrks1bin[h][k] 							= 0;
			fcn6Gap0Ntrks1bin[h][k] 					= 0;

			fcn8Ntrks1bin[h][k] 							= 0;
			fcn8Gap0Ntrks1bin[h][k] 					= 0;

	 }
       
  }
    
}
//______________________________________________________________________________
AliAnalysisTaskChargedFlow::AliAnalysisTaskChargedFlow(const char *name):
  AliAnalysisTaskSE(name),
	fEtaCut(0.8),
	fVtxCut(10.0),  
	fMinPt(0.2),
	fMaxPt(3.0),
	fSample(10),

	fOutputContainer(0),

	hMult(0), 
	hMultTrigger(0), 
	fVtxAfterCuts(0),

	fPhiDis(0), 
	fEtaDis(0),
	fEtaBefore(0),
	fPtDis(0),
	fPtBefore(0),
	hIsPiKp(0),

	fChsc4242(0),
	fChsc4242Gap0(0),
	fChsc4242_3sub(0),
	fChsc4224_3sub(0),
	fChsc4242_3sub_perm2(0),
	fChsc4224_3sub_perm2(0),
	fChsc4242_3sub_perm3(0),
	fChsc4224_3sub_perm3(0),
	fChsc3232(0),
	fChsc3232Gap0(0),
	fChsc3232_3sub(0),
	fChsc3223_3sub(0),
	fChsc3232_3sub_perm2(0),
	fChsc3223_3sub_perm2(0),
	fChsc3232_3sub_perm3(0),
	fChsc3223_3sub_perm3(0)
{

	for(int h=0; h<3; h++)
	{
 
		fChcn2Ntrks1bin[h] 									= 0;
		fChcn2Gap0Ntrks1bin[h] 							= 0;
		fChcn2Gap2Ntrks1bin[h] 							= 0;
		fChcn2Gap4Ntrks1bin[h] 							= 0;
		fChcn2Gap8Ntrks1bin[h] 							= 0;
		fChcn2GapNtrks1bin[h] 							= 0;
		fChcn2Gap14Ntrks1bin[h] 						= 0;

		fChcn2_3subLMNtrks1bin[h] 					= 0;
		fChcn2_3subRMNtrks1bin[h] 					= 0;
		fChcn2_3subLM_perm2Ntrks1bin[h] 		= 0;
		fChcn2_3subLR_perm2Ntrks1bin[h] 		= 0;
		fChcn2_3subRM_perm3Ntrks1bin[h] 		= 0;
		fChcn2_3subRL_perm3Ntrks1bin[h] 		= 0;

		fChcn4Ntrks1bin[h] 									= 0;
		fChcn4Gap0Ntrks1bin[h] 							= 0;
		fChcn4_3subNtrks1bin[h] 						= 0;
		fChcn4_3sub_perm2Ntrks1bin[h]				= 0;
		fChcn4_3sub_perm3Ntrks1bin[h]				= 0;

		fChcn6Ntrks1bin[h] 									= 0;
		fChcn6Gap0Ntrks1bin[h] 							= 0;

		fChcn8Ntrks1bin[h] 									= 0;
		fChcn8Gap0Ntrks1bin[h] 							= 0;

 	}

	for(int k=0; k<fSample+1; k++)
	{			

		fsc4242[k] 													= 0;
		fsc4242Gap0[k] 											= 0;
		fsc4242_3sub[k] 										= 0;
		fsc4224_3sub[k] 										= 0;
		fsc4242_3sub_perm2[k] 						= 0;
		fsc4224_3sub_perm2[k] 						= 0;
		fsc4242_3sub_perm3[k]							= 0;
		fsc4224_3sub_perm3[k]							= 0;

		fsc3232[k] 													= 0;
		fsc3232Gap0[k] 											= 0;
		fsc3232_3sub[k] 										= 0;
		fsc3223_3sub[k] 										= 0;
		fsc3232_3sub_perm2[k] 						= 0;
		fsc3223_3sub_perm2[k] 						= 0;
		fsc3232_3sub_perm3[k] 						= 0;
		fsc3223_3sub_perm3[k] 						= 0;

		for(Int_t h = 0; h < 3; h++)	
		{		

			fcn2Ntrks1bin[h][k] 							= 0;
			fcn2GapNtrks1bin[h][k] 						= 0;
			fcn2Gap0Ntrks1bin[h][k] 					= 0;
			fcn2Gap2Ntrks1bin[h][k] 					= 0;
			fcn2Gap4Ntrks1bin[h][k] 					= 0;
			fcn2Gap8Ntrks1bin[h][k] 					= 0;
			fcn2Gap14Ntrks1bin[h][k] 					= 0;

			fcn2_3subLMNtrks1bin[h][k] 				= 0;
			fcn2_3subRMNtrks1bin[h][k] 				= 0;
			fcn2_3subLM_perm2Ntrks1bin[h][k] 				= 0;
			fcn2_3subLR_perm2Ntrks1bin[h][k] 				= 0;
			fcn2_3subRM_perm3Ntrks1bin[h][k] 				= 0;
			fcn2_3subRL_perm3Ntrks1bin[h][k] 				= 0;

			fcn4Ntrks1bin[h][k] 							= 0;
			fcn4Gap0Ntrks1bin[h][k] 					= 0;
			fcn4_3subNtrks1bin[h][k] 					= 0;	
			fcn4_3sub_perm2Ntrks1bin[h][k]		= 0;
			fcn4_3sub_perm3Ntrks1bin[h][k]		= 0;

			fcn6Ntrks1bin[h][k] 							= 0;
			fcn6Gap0Ntrks1bin[h][k] 					= 0;

			fcn8Ntrks1bin[h][k]			 					= 0;
			fcn8Gap0Ntrks1bin[h][k] 					= 0;

	 }
       
  }
    
	// Output slot #1 writes into a TList
  DefineOutput(1, TList::Class());
	
}

//_____________________________________________________________________________
AliAnalysisTaskChargedFlow::~AliAnalysisTaskChargedFlow()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
  //if (fOutputContainer) 
    //delete fOutputContainer;

}

//______________________________________________________________________________
void AliAnalysisTaskChargedFlow::UserCreateOutputObjects()
{ 
	
  //fOutputContainer = new TList();
  //fOutputContainer->SetOwner();

  // Create histograms
  if(fOutputContainer != nullptr){
    delete fOutputContainer;
    fOutputContainer          = nullptr;
  }
  if(fOutputContainer == nullptr){
    fOutputContainer          = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }

	// range on Xaxis: 
	//		for pp, pPb and for lowM PbPb and XeXe the range is 200 in unit bins
	//		for PbPb and XeXe the range is 5000 in bin width 50 (low M part which won't be correct due to mult. fluctuations will be done from unit bins in separate running)
	int range = 200;
	const int nn = 200;
	int rangeV0 = 300;
	const int nnV0 = 60;

  // Event Histograms
	hMult = new TH2F("hMult", ";number of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
	hMult->Sumw2();
	fOutputContainer->Add(hMult);

	hMultTrigger = new TH1F("hMultTrigger", ";number of tracks; entries", nnV0, 0, rangeV0);
	hMultTrigger->Sumw2();
	fOutputContainer->Add(hMultTrigger);

  fVtxAfterCuts = new TH1F("fVtxAfterCuts", "Vtx distribution (after cuts); Vtx z [cm]; Counts", 120, -30, 30);
	fVtxAfterCuts->Sumw2();
  fOutputContainer->Add(fVtxAfterCuts);

	// Track histograms  
  fPhiDis = new TH2F("fPhiDis", "phi distribution; #phi; V0 mult.", 360, 0, 2*TMath::Pi(), nnV0, 0, rangeV0); //..18 sectors of TPC, 20 bins per each sector
	fPhiDis->Sumw2();
  fOutputContainer->Add(fPhiDis);
    
  fEtaDis = new TH2F("fEtaDis", "eta distribution; #eta; V0 mult.", 200, -2., 2., nnV0, 0, rangeV0);
	fEtaDis->Sumw2();
  fOutputContainer->Add(fEtaDis);
   
  fEtaBefore = new TH2F("fEtaBefore", "eta distribution; #eta; V0 mult.", 400, -6., 6., nnV0, 0, rangeV0);
	fEtaBefore->Sumw2();
  fOutputContainer->Add(fEtaBefore);
   
	fPtDis = new TH2F("fPtDis", "pt distribution; p_{T}; V0 mult.", 100, 0, 10, nnV0, 0, rangeV0);
	fPtDis->Sumw2();
	fOutputContainer->Add(fPtDis);
 
	fPtBefore = new TH2F("fPtBefore", "pt distribution; p_{T}; V0 mult.", 100, 0, 10, nnV0, 0, rangeV0);
	fPtBefore->Sumw2();
	fOutputContainer->Add(fPtBefore);
 
	hIsPiKp = new TH1F("hIsPiKp", "", 1, 0, 1);
	fOutputContainer->Add(hIsPiKp);

	// Physics profiles
	// SC(n,m)
	fChsc4242 = new TProfile2D("fChsc4242", ";# of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
	fChsc4242->Sumw2();
	fOutputContainer->Add(fChsc4242);

	fChsc4242Gap0 = new TProfile2D("fChsc4242Gap0", ";# of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
	fChsc4242Gap0->Sumw2();
	fOutputContainer->Add(fChsc4242Gap0);

	fChsc4242_3sub = new TProfile2D("fChsc4242_3sub", ";# of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
	fChsc4242_3sub->Sumw2();
	fOutputContainer->Add(fChsc4242_3sub);

	fChsc4224_3sub = new TProfile2D("fChsc4224_3sub", ";# of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
	fChsc4224_3sub->Sumw2();
	fOutputContainer->Add(fChsc4224_3sub);

	fChsc4242_3sub_perm2 = new TProfile2D("fChsc4242_3sub_perm2", ";# of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
	fChsc4242_3sub_perm2->Sumw2();
	fOutputContainer->Add(fChsc4242_3sub_perm2);

	fChsc4224_3sub_perm2 = new TProfile2D("fChsc4224_3sub_perm2", ";# of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
	fChsc4224_3sub_perm2->Sumw2();
	fOutputContainer->Add(fChsc4224_3sub_perm2);

	fChsc4242_3sub_perm3 = new TProfile2D("fChsc4242_3sub_perm3", ";# of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
	fChsc4242_3sub_perm3->Sumw2();
	fOutputContainer->Add(fChsc4242_3sub_perm3);

	fChsc4224_3sub_perm3 = new TProfile2D("fChsc4224_3sub_perm3", ";# of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
	fChsc4224_3sub_perm3->Sumw2();
	fOutputContainer->Add(fChsc4224_3sub_perm3);


	fChsc3232 = new TProfile2D("fChsc3232", ";# of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
	fChsc3232->Sumw2();
	fOutputContainer->Add(fChsc3232);

	fChsc3232Gap0 = new TProfile2D("fChsc3232Gap0", ";# of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
	fChsc3232Gap0->Sumw2();
	fOutputContainer->Add(fChsc3232Gap0);

	fChsc3232_3sub = new TProfile2D("fChsc3232_3sub", ";# of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
	fChsc3232_3sub->Sumw2();
	fOutputContainer->Add(fChsc3232_3sub);

	fChsc3223_3sub = new TProfile2D("fChsc3223_3sub", ";# of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
	fChsc3223_3sub->Sumw2();
	fOutputContainer->Add(fChsc3223_3sub);

	fChsc3232_3sub_perm2 = new TProfile2D("fChsc3232_3sub_perm2", ";# of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
	fChsc3232_3sub_perm2->Sumw2();
	fOutputContainer->Add(fChsc3232_3sub_perm2);

	fChsc3223_3sub_perm2 = new TProfile2D("fChsc3223_3sub_perm2", ";# of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
	fChsc3223_3sub_perm2->Sumw2();
	fOutputContainer->Add(fChsc3223_3sub_perm2);

	fChsc3232_3sub_perm3 = new TProfile2D("fChsc3232_3sub_perm3", ";# of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
	fChsc3232_3sub_perm3->Sumw2();
	fOutputContainer->Add(fChsc3232_3sub_perm3);

	fChsc3223_3sub_perm3 = new TProfile2D("fChsc3223_3sub_perm3", ";# of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
	fChsc3223_3sub_perm3->Sumw2();
	fOutputContainer->Add(fChsc3223_3sub_perm3);


	for(int h=0; h<3; h++)
	{
	  fChcn2Ntrks1bin[h] = new TProfile2D(Form("fChc%d2Ntrks1bin", h+2), "<<2>> Re; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
    fChcn2Ntrks1bin[h]->Sumw2();
    fOutputContainer->Add(fChcn2Ntrks1bin[h]);

    fChcn2Gap0Ntrks1bin[h] = new TProfile2D(Form("fChc%d2Gap0Ntrks1bin", h+2), "<<2>> Re; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
    fChcn2Gap0Ntrks1bin[h]->Sumw2();
    fOutputContainer->Add(fChcn2Gap0Ntrks1bin[h]);

    fChcn2Gap2Ntrks1bin[h] = new TProfile2D(Form("fChc%d2Gap2Ntrks1bin", h+2), "<<2>> Re; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
    fChcn2Gap2Ntrks1bin[h]->Sumw2();
    fOutputContainer->Add(fChcn2Gap2Ntrks1bin[h]);

    fChcn2Gap4Ntrks1bin[h] = new TProfile2D(Form("fChc%d2Gap4Ntrks1bin", h+2), "<<2>> Re; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
    fChcn2Gap4Ntrks1bin[h]->Sumw2();
    fOutputContainer->Add(fChcn2Gap4Ntrks1bin[h]);

    fChcn2Gap8Ntrks1bin[h] = new TProfile2D(Form("fChc%d2Gap8Ntrks1bin", h+2), "<<2>> Re; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
    fChcn2Gap8Ntrks1bin[h]->Sumw2();
    fOutputContainer->Add(fChcn2Gap8Ntrks1bin[h]);

    fChcn2GapNtrks1bin[h] = new TProfile2D(Form("fChc%d2GapNtrks1bin", h+2), "<<2>> Re; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
    fChcn2GapNtrks1bin[h]->Sumw2();
    fOutputContainer->Add(fChcn2GapNtrks1bin[h]);

    fChcn2Gap14Ntrks1bin[h] = new TProfile2D(Form("fChc%d2Gap14Ntrks1bin", h+2), "<<2>> Re; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
    fChcn2Gap14Ntrks1bin[h]->Sumw2();
    fOutputContainer->Add(fChcn2Gap14Ntrks1bin[h]);


    fChcn2_3subLMNtrks1bin[h] = new TProfile2D(Form("fChc%d2_3subLMNtrks1bin", h+2), "<<2>> Re for 3-subevent method, left+middle; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
    fChcn2_3subLMNtrks1bin[h]->Sumw2();
    fOutputContainer->Add(fChcn2_3subLMNtrks1bin[h]);

    fChcn2_3subRMNtrks1bin[h] = new TProfile2D(Form("fChc%d2_3subRMNtrks1bin", h+2), "<<2>> Re for 3-subevent method, right+middle; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
    fChcn2_3subRMNtrks1bin[h]->Sumw2();
    fOutputContainer->Add(fChcn2_3subRMNtrks1bin[h]);

    fChcn2_3subLM_perm2Ntrks1bin[h] = new TProfile2D(Form("fChc%d2_3subLM_perm2Ntrks1bin", h+2), "<<2>> Re for 3-subevent method, left+middle; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
    fChcn2_3subLM_perm2Ntrks1bin[h]->Sumw2();
    fOutputContainer->Add(fChcn2_3subLM_perm2Ntrks1bin[h]);

    fChcn2_3subLR_perm2Ntrks1bin[h] = new TProfile2D(Form("fChc%d2_3subLR_perm2Ntrks1bin", h+2), "<<2>> Re for 3-subevent method, right+middle; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
    fChcn2_3subLR_perm2Ntrks1bin[h]->Sumw2();
    fOutputContainer->Add(fChcn2_3subLR_perm2Ntrks1bin[h]);

    fChcn2_3subRL_perm3Ntrks1bin[h] = new TProfile2D(Form("fChc%d2_3subRL_perm3Ntrks1bin", h+2), "<<2>> Re for 3-subevent method, left+middle; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
    fChcn2_3subRL_perm3Ntrks1bin[h]->Sumw2();
    fOutputContainer->Add(fChcn2_3subRL_perm3Ntrks1bin[h]);

    fChcn2_3subRM_perm3Ntrks1bin[h] = new TProfile2D(Form("fChc%d2_3subRM_perm3Ntrks1bin", h+2), "<<2>> Re for 3-subevent method, right+middle; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
    fChcn2_3subRM_perm3Ntrks1bin[h]->Sumw2();
    fOutputContainer->Add(fChcn2_3subRM_perm3Ntrks1bin[h]);


	  fChcn4Ntrks1bin[h] = new TProfile2D(Form("fChc%d4Ntrks1bin", h+2), "<<4>> Re; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
    fChcn4Ntrks1bin[h]->Sumw2();
    fOutputContainer->Add(fChcn4Ntrks1bin[h]);

    fChcn4Gap0Ntrks1bin[h] = new TProfile2D(Form("fChc%d4Gap0Ntrks1bin", h+2), "<<4>> Re; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
    fChcn4Gap0Ntrks1bin[h]->Sumw2();
    fOutputContainer->Add(fChcn4Gap0Ntrks1bin[h]);

		fChcn4_3subNtrks1bin[h] = new TProfile2D(Form("fChc%d4_3subNtrks1bin", h+2), "<<4>> 3-subevent method; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
		fChcn4_3subNtrks1bin[h]->Sumw2();
		fOutputContainer->Add(fChcn4_3subNtrks1bin[h]);
	
		fChcn4_3sub_perm2Ntrks1bin[h] = new TProfile2D(Form("fChc%d4_3sub_perm2Ntrks1bin", h+2), "<<4>> for 3-subevent method", nn, 0, range, nnV0, 0, rangeV0);
		fChcn4_3sub_perm2Ntrks1bin[h]->Sumw2();
		fOutputContainer->Add(fChcn4_3sub_perm2Ntrks1bin[h]);

		fChcn4_3sub_perm3Ntrks1bin[h] = new TProfile2D(Form("fChc%d4_3sub_perm3Ntrks1bin", h+2), "<<4>> for 3-subevent method", nn, 0, range, nnV0, 0, rangeV0);
		fChcn4_3sub_perm3Ntrks1bin[h]->Sumw2();
		fOutputContainer->Add(fChcn4_3sub_perm3Ntrks1bin[h]);


	  fChcn6Ntrks1bin[h] = new TProfile2D(Form("fChc%d6Ntrks1bin", h+2), "<<6>> Re; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
    fChcn6Ntrks1bin[h]->Sumw2();
    fOutputContainer->Add(fChcn6Ntrks1bin[h]);

	  fChcn6Gap0Ntrks1bin[h] = new TProfile2D(Form("fChc%d6Gap0Ntrks1bin", h+2), "<<6>> Gap0 Re; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
    fChcn6Gap0Ntrks1bin[h]->Sumw2();
    fOutputContainer->Add(fChcn6Gap0Ntrks1bin[h]);


	  fChcn8Ntrks1bin[h] = new TProfile2D(Form("fChc%d8Ntrks1bin", h+2), "<<8>>  Re; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
    fChcn8Ntrks1bin[h]->Sumw2();
    fOutputContainer->Add(fChcn8Ntrks1bin[h]);

	  fChcn8Gap0Ntrks1bin[h] = new TProfile2D(Form("fChc%d8Gap0Ntrks1bin", h+2), "<<8>> Gap0 Re; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
    fChcn8Gap0Ntrks1bin[h]->Sumw2();
    fOutputContainer->Add(fChcn8Gap0Ntrks1bin[h]);

	}// harmonics

	for(int j=0; j<fSample+1; j++)
	{

		fsc4242[j] = new TProfile2D(Form("fsc4242_number%d", j), "; Ntrks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
		fsc4242[j]->Sumw2();
		fOutputContainer->Add(fsc4242[j]);

		fsc4242Gap0[j] = new TProfile2D(Form("fsc4242Gap0_number%d", j), "; Ntrks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
		fsc4242Gap0[j]->Sumw2();
		fOutputContainer->Add(fsc4242Gap0[j]);

		fsc4242_3sub[j] = new TProfile2D(Form("fsc4242_3sub_number%d", j), "; Ntrks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
		fsc4242_3sub[j]->Sumw2();
		fOutputContainer->Add(fsc4242_3sub[j]);

		fsc4224_3sub[j] = new TProfile2D(Form("fsc4224_3sub_number%d", j), "; Ntrks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
		fsc4224_3sub[j]->Sumw2();
		fOutputContainer->Add(fsc4224_3sub[j]);

		fsc4242_3sub_perm2[j] = new TProfile2D(Form("fsc4242_3sub_perm2_number%d", j), "; Ntrks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
		fsc4242_3sub_perm2[j]->Sumw2();
		fOutputContainer->Add(fsc4242_3sub_perm2[j]);

		fsc4224_3sub_perm2[j] = new TProfile2D(Form("fsc4224_3sub_perm2_number%d", j), "; Ntrks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
		fsc4224_3sub_perm2[j]->Sumw2();
		fOutputContainer->Add(fsc4224_3sub_perm2[j]);

		fsc4242_3sub_perm3[j] = new TProfile2D(Form("fsc4242_3sub_perm3_number%d", j), "; Ntrks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
		fsc4242_3sub_perm3[j]->Sumw2();
		fOutputContainer->Add(fsc4242_3sub_perm3[j]);

		fsc4224_3sub_perm3[j] = new TProfile2D(Form("fsc4224_3sub_perm3_number%d", j), "; Ntrks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
		fsc4224_3sub_perm3[j]->Sumw2();
		fOutputContainer->Add(fsc4224_3sub_perm3[j]);


		fsc3232[j] = new TProfile2D(Form("fsc3232_number%d", j), "; Ntrks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
		fsc3232[j]->Sumw2();
		fOutputContainer->Add(fsc3232[j]);

		fsc3232Gap0[j] = new TProfile2D(Form("fsc3232Gap0_number%d", j), "; Ntrks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
		fsc3232Gap0[j]->Sumw2();
		fOutputContainer->Add(fsc3232Gap0[j]);

		fsc3232_3sub[j] = new TProfile2D(Form("fsc3232_3sub_number%d", j), "; Ntrks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
		fsc3232_3sub[j]->Sumw2();
		fOutputContainer->Add(fsc3232_3sub[j]);

		fsc3223_3sub[j] = new TProfile2D(Form("fsc3223_3sub_number%d", j), "; Ntrks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
		fsc3223_3sub[j]->Sumw2();
		fOutputContainer->Add(fsc3223_3sub[j]);

		fsc3232_3sub_perm2[j] = new TProfile2D(Form("fsc3232_3sub_perm2_number%d", j), "; Ntrks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
		fsc3232_3sub_perm2[j]->Sumw2();
		fOutputContainer->Add(fsc3232_3sub_perm2[j]);

		fsc3223_3sub_perm2[j] = new TProfile2D(Form("fsc3223_3sub_perm2_number%d", j), "; Ntrks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
		fsc3223_3sub_perm2[j]->Sumw2();
		fOutputContainer->Add(fsc3223_3sub_perm2[j]);

		fsc3232_3sub_perm3[j] = new TProfile2D(Form("fsc3232_3sub_perm3_number%d", j), "; Ntrks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
		fsc3232_3sub_perm3[j]->Sumw2();
		fOutputContainer->Add(fsc3232_3sub_perm3[j]);

		fsc3223_3sub_perm3[j] = new TProfile2D(Form("fsc3223_3sub_perm3_number%d", j), "; Ntrks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
		fsc3223_3sub_perm3[j]->Sumw2();
		fOutputContainer->Add(fsc3223_3sub_perm3[j]);

		for(int h=0; h<3; h++)
		{
	    fcn2Ntrks1bin[h][j] = new TProfile2D(Form("fc%d2_number%dNtrks1bin", h+2, j), "<<2>> Re; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
      fcn2Ntrks1bin[h][j]->Sumw2();
      fOutputContainer->Add(fcn2Ntrks1bin[h][j]);

      fcn2Gap0Ntrks1bin[h][j] = new TProfile2D(Form("fc%d2_number%dGap0Ntrks1bin", h+2, j), "<<2>> Re; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
      fcn2Gap0Ntrks1bin[h][j]->Sumw2();
      fOutputContainer->Add(fcn2Gap0Ntrks1bin[h][j]);

      fcn2Gap2Ntrks1bin[h][j] = new TProfile2D(Form("fc%d2_number%dGap2Ntrks1bin", h+2, j), "<<2>> Re; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
      fcn2Gap2Ntrks1bin[h][j]->Sumw2();
      fOutputContainer->Add(fcn2Gap2Ntrks1bin[h][j]);

      fcn2Gap4Ntrks1bin[h][j] = new TProfile2D(Form("fc%d2_number%dGap4Ntrks1bin", h+2, j), "<<2>> Re; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
      fcn2Gap4Ntrks1bin[h][j]->Sumw2();
      fOutputContainer->Add(fcn2Gap4Ntrks1bin[h][j]);

      fcn2Gap8Ntrks1bin[h][j] = new TProfile2D(Form("fc%d2_number%dGap8Ntrks1bin", h+2, j), "<<2>> Re; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
      fcn2Gap8Ntrks1bin[h][j]->Sumw2();
      fOutputContainer->Add(fcn2Gap8Ntrks1bin[h][j]);

      fcn2GapNtrks1bin[h][j] = new TProfile2D(Form("fc%d2_number%dGapNtrks1bin", h+2, j), "<<2>> Re; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
      fcn2GapNtrks1bin[h][j]->Sumw2();
      fOutputContainer->Add(fcn2GapNtrks1bin[h][j]);

      fcn2Gap14Ntrks1bin[h][j] = new TProfile2D(Form("fc%d2_number%dGap14Ntrks1bin", h+2, j), "<<2>> Re; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
      fcn2Gap14Ntrks1bin[h][j]->Sumw2();
      fOutputContainer->Add(fcn2Gap14Ntrks1bin[h][j]);


      fcn2_3subLMNtrks1bin[h][j] = new TProfile2D(Form("fc%d2_number%d_3subLMNtrks1bin", h+2, j), "<<2>> Re for 3-subevent method left+middle; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
      fcn2_3subLMNtrks1bin[h][j]->Sumw2();
      fOutputContainer->Add(fcn2_3subLMNtrks1bin[h][j]);

      fcn2_3subRMNtrks1bin[h][j] = new TProfile2D(Form("fc%d2_number%d_3subRMNtrks1bin", h+2, j), "<<2>> Re for 3-subevent method right+middle; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
      fcn2_3subRMNtrks1bin[h][j]->Sumw2();
      fOutputContainer->Add(fcn2_3subRMNtrks1bin[h][j]);

      fcn2_3subLM_perm2Ntrks1bin[h][j] = new TProfile2D(Form("fc%d2_number%d_3subLM_perm2Ntrks1bin", h+2, j), "<<2>> Re for 3-subevent method left+middle; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
      fcn2_3subLM_perm2Ntrks1bin[h][j]->Sumw2();
      fOutputContainer->Add(fcn2_3subLM_perm2Ntrks1bin[h][j]);

      fcn2_3subLR_perm2Ntrks1bin[h][j] = new TProfile2D(Form("fc%d2_number%d_3subLR_perm2Ntrks1bin", h+2, j), "<<2>> Re for 3-subevent method right+middle; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
      fcn2_3subLR_perm2Ntrks1bin[h][j]->Sumw2();
      fOutputContainer->Add(fcn2_3subLR_perm2Ntrks1bin[h][j]);

      fcn2_3subRM_perm3Ntrks1bin[h][j] = new TProfile2D(Form("fc%d2_number%d_3subRM_perm3Ntrks1bin", h+2, j), "<<2>> Re for 3-subevent method left+middle; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
      fcn2_3subRM_perm3Ntrks1bin[h][j]->Sumw2();
      fOutputContainer->Add(fcn2_3subRM_perm3Ntrks1bin[h][j]);

      fcn2_3subRL_perm3Ntrks1bin[h][j] = new TProfile2D(Form("fc%d2_number%d_3subRL_perm3Ntrks1bin", h+2, j), "<<2>> Re for 3-subevent method right+middle; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
      fcn2_3subRL_perm3Ntrks1bin[h][j]->Sumw2();
      fOutputContainer->Add(fcn2_3subRL_perm3Ntrks1bin[h][j]);


			fcn4Ntrks1bin[h][j] = new TProfile2D(Form("fc%d4_number%dNtrks1bin", h+2, j), "<<4>>; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
			fcn4Ntrks1bin[h][j]->Sumw2();
			fOutputContainer->Add(fcn4Ntrks1bin[h][j]);
	
			fcn4Gap0Ntrks1bin[h][j] = new TProfile2D(Form("fc%d4Gap0_number%dNtrks1bin", h+2, j), "<<4>>; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
			fcn4Gap0Ntrks1bin[h][j]->Sumw2();
			fOutputContainer->Add(fcn4Gap0Ntrks1bin[h][j]);
	
			fcn4_3subNtrks1bin[h][j] = new TProfile2D(Form("fc%d4_3sub_number%dNtrks1bin", h+2, j), "<<4>> 3-subevent method; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
			fcn4_3subNtrks1bin[h][j]->Sumw2();
			fOutputContainer->Add(fcn4_3subNtrks1bin[h][j]);

			fcn4_3sub_perm2Ntrks1bin[h][j] = new TProfile2D(Form("fc%d4_number%d_3sub_perm2Ntrks1bin", h+2, j), "<<4>> for 3-subevent method", nn, 0, range, nnV0, 0, rangeV0);
			fcn4_3sub_perm2Ntrks1bin[h][j]->Sumw2();
			fOutputContainer->Add(fcn4_3sub_perm2Ntrks1bin[h][j]);

			fcn4_3sub_perm3Ntrks1bin[h][j] = new TProfile2D(Form("fc%d4_number%d_3sub_perm3Ntrks1bin", h+2, j), "<<4>> for 3-subevent method", nn, 0, range, nnV0, 0, rangeV0);
			fcn4_3sub_perm3Ntrks1bin[h][j]->Sumw2();
			fOutputContainer->Add(fcn4_3sub_perm3Ntrks1bin[h][j]);


			fcn6Ntrks1bin[h][j] = new TProfile2D(Form("fc%d6_number%dNtrks1bin", h+2, j), "<<6>>; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
			fcn6Ntrks1bin[h][j]->Sumw2();
			fOutputContainer->Add(fcn6Ntrks1bin[h][j]);
	
			fcn6Gap0Ntrks1bin[h][j] = new TProfile2D(Form("fc%d6Gap0_number%dNtrks1bin", h+2, j), "<<6>> Gap0; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
			fcn6Gap0Ntrks1bin[h][j]->Sumw2();
			fOutputContainer->Add(fcn6Gap0Ntrks1bin[h][j]);
	

			fcn8Ntrks1bin[h][j] = new TProfile2D(Form("fc%d8_number%dNtrks1bin", h+2, j), "<<8>> ; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
			fcn8Ntrks1bin[h][j]->Sumw2();
			fOutputContainer->Add(fcn8Ntrks1bin[h][j]);

			fcn8Gap0Ntrks1bin[h][j] = new TProfile2D(Form("fc%d8Gap0_number%dNtrks1bin", h+2, j), "<<8>> Gap0; # of tracks; V0 mult.", nn, 0, range, nnV0, 0, rangeV0);
			fcn8Gap0Ntrks1bin[h][j]->Sumw2();
			fOutputContainer->Add(fcn8Gap0Ntrks1bin[h][j]);

		}// harmonics

	}// samples

  // Post output data.
  PostData(1, fOutputContainer);
}

//______________________________________________________________________________
void AliAnalysisTaskChargedFlow::UserExec(Option_t *) 
{

	fInputEvent = InputEvent();

	fMCEvent = MCEvent();
	if(fMCEvent == nullptr) return;

	//..filling Vz distribution
	const AliVVertex *vtx = fMCEvent->GetPrimaryVertex();
	float fVtxZ = vtx->GetZ();	  
	if(TMath::Abs(fVtxZ) > fVtxCut) return;
	fVtxAfterCuts->Fill(fVtxZ);

	AnalyzeAOD();

  // Post output data.
	PostData(1, fOutputContainer);

}
//________________________________________________________________________
void AliAnalysisTaskChargedFlow::AnalyzeAOD()
{  

	AliStack *stack = fMCEvent->Stack();
	const int nTracks = stack->GetNtrack();

	double NtrksTrigger = 0;
	double NtrksAfter = 0;
	double NtrksAfterGap0M = 0;
	double NtrksAfterGap0P = 0;
	double NtrksAfterGap1M = 0;
	double NtrksAfterGap1P = 0;
	double NtrksAfterGap2M = 0;
	double NtrksAfterGap2P = 0;
	double NtrksAfterGap4M = 0;
	double NtrksAfterGap4P = 0;
	double NtrksAfterGap8M = 0;
	double NtrksAfterGap8P = 0;
	double NtrksAfterGapM = 0;
	double NtrksAfterGapP = 0;
	double NtrksAfterGap14M = 0;
	double NtrksAfterGap14P = 0;
	double NtrksAfter3subL = 0;
	double NtrksAfter3subM = 0;
	double NtrksAfter3subR = 0;

	double Qcos[20][20] = {0};
	double Qsin[20][20] = {0};
	double QcosGapM[20][20] = {0};
	double QsinGapM[20][20] = {0};
	double QcosGapP[20][20] = {0};
	double QsinGapP[20][20] = {0};

	double QcosGap0M[20][20] = {0};
	double QsinGap0M[20][20] = {0};
	double QcosGap0P[20][20] = {0};
	double QsinGap0P[20][20] = {0};
  
	double QcosGap2M[20][20] = {0};
	double QsinGap2M[20][20] = {0};
	double QcosGap2P[20][20] = {0};
	double QsinGap2P[20][20] = {0};
  
	double QcosGap4M[20][20] = {0};
	double QsinGap4M[20][20] = {0};
	double QcosGap4P[20][20] = {0};
	double QsinGap4P[20][20] = {0};
	
	double QcosGap8M[20][20] = {0};
	double QsinGap8M[20][20] = {0};
	double QcosGap8P[20][20] = {0};
	double QsinGap8P[20][20] = {0};
    
	double QcosGap14M[20][20] = {0};
	double QsinGap14M[20][20] = {0};
	double QcosGap14P[20][20] = {0};
	double QsinGap14P[20][20] = {0};
    
	double QcosSubLeft[20][20] = {0};
	double QsinSubLeft[20][20] = {0};
	double QcosSubMiddle[20][20] = {0};
	double QsinSubMiddle[20][20] = {0};
	double QcosSubRight[20][20] = {0};
	double QsinSubRight[20][20] = {0};
    
  Int_t fBin = 0;
  TRandom3 rr(0);
  Double_t ranNum = 0.;

		if(fSample == 1) fBin = 1;
    
    if(fSample == 10)
		{
    	ranNum = rr.Rndm();
        
    	if (ranNum <= 0.1 && ranNum > 0. ) fBin = 1;
    	else if (ranNum <= 0.2 && ranNum > 0.1) fBin = 2;
    	else if (ranNum <= 0.3 && ranNum > 0.2) fBin = 3;
    	else if (ranNum <= 0.4 && ranNum > 0.3) fBin = 4;
    	else if (ranNum <= 0.5 && ranNum > 0.4) fBin = 5;
    	else if (ranNum <= 0.6 && ranNum > 0.5) fBin = 6;
    	else if (ranNum <= 0.7 && ranNum > 0.6) fBin = 7;
    	else if (ranNum <= 0.8 && ranNum > 0.7) fBin = 8;
    	else if (ranNum <= 0.9 && ranNum > 0.8) fBin = 9;
    	else if (ranNum > 0.9 && ranNum <= 1.) fBin = 10;
    	else {
         cout << ranNum << endl;
         
         cout << "This is very strange!!!" << endl;
         fBin = 11;
    	}

    }

	//..LOOP OVER TRACKS........	
	//........................................
  for(Int_t nt = 0; nt < nTracks; nt++)
  {

		TParticle *particle = (TParticle*)stack->Particle(nt);
		if(!particle) continue;

		//if(!(particle->GetPdgCode())) continue;

		//	check if it is primary by checking if it has mother or not
		bool isPrimary = false;
		if(stack->IsPhysicalPrimary(nt)) isPrimary = true;

		//	charge
		bool isCharged = false;
		int pid = TMath::Abs(particle->GetPdgCode());
		if((pid == 211) || (pid == 321) || (pid == 2212)) {
			hIsPiKp->Fill("is charged pi, K, p", 1.);
			isCharged = true;
		}else {
			hIsPiKp->Fill("other than pi, K, p", 1.);
		}

		double eta = particle->Eta(); 
		double pt = particle->Pt();
		double phi = particle->Phi();

		//	get multiplicity to mimic HM trigger
		if((isPrimary == true) && (isCharged == true) && (pt > 0.) && ((eta > -3.7 && eta < -1.7) || (eta > 2.8 && eta < 5.1)))
			NtrksTrigger += 1;

		//	do selection
		if(isPrimary == false) continue;
		if(isCharged == false) continue;

		fEtaBefore->Fill(eta, NtrksTrigger);
		fPtBefore->Fill(pt, NtrksTrigger);

		if(pt < fMinPt) continue;
		if(pt > fMaxPt) continue;

    if(TMath::Abs(eta) > fEtaCut) continue;

		fPhiDis->Fill(phi, NtrksTrigger);
    fEtaDis->Fill(eta, NtrksTrigger);
		fPtDis->Fill(pt, NtrksTrigger);     
 
		NtrksAfter += 1;

		double weight = 1;
		double weightPt = 1;

			//..calculate Q-vectors
			//..no eta gap
			for(int iharm=0; iharm<10; iharm++)
			{
				for(int ipow=0; ipow<10; ipow++)
				{
					Qcos[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*phi);
					Qsin[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*phi);
				}
			}

      //..Gap > 0.0
      if (eta > 0.)
      {
				NtrksAfterGap0P += 1;
         
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap0P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*phi);
						QsinGap0P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*phi);
					}
				}
 
      }
      if (eta < -0.)
      {
				NtrksAfterGap0M += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap0M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*phi);
						QsinGap0M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*phi);
					}
				}
      }
		
      //..Gap > 0.2
      if (eta > 0.1)
      {
        NtrksAfterGap2P += 1; 
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap2P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*phi);
						QsinGap2P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*phi);
					}
				}
 
      }
      if (eta < -0.1)
      {
				NtrksAfterGap2M += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap2M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*phi);
						QsinGap2M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*phi);
					}
				}
      }
		
      //..Gap > 0.4
      if (eta > 0.2)
      {
				NtrksAfterGap4P += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap4P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*phi);
						QsinGap4P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*phi);
					}
				}
      }
      if (eta < -0.2)
      {
				NtrksAfterGap4M += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap4M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*phi);
						QsinGap4M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*phi);
					}
				}
      }
	 	
      //..Gap > 0.8
      if (eta > 0.4)
      {
				NtrksAfterGap8P += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap8P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*phi);
						QsinGap8P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*phi);
					}
				}
      }
      if (eta < -0.4)
      {
				NtrksAfterGap8M += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap8M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*phi);
						QsinGap8M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*phi);
					}
				}

      }

			//..Gap > 1.0
			if(eta < -0.5)
			{
				NtrksAfterGapM += 1;
				for(int iharm=0; iharm<20; iharm++)
				{
					for(int ipow=0; ipow<20; ipow++)
					{
						QcosGapM[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*phi);
						QsinGapM[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*phi);
					}
				}
			} 
			if(eta > 0.5)
			{
				NtrksAfterGapP += 1;
				for(int iharm=0; iharm<20; iharm++)
				{
					for(int ipow=0; ipow<20; ipow++)
					{
						QcosGapP[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*phi);
						QsinGapP[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*phi);	
					}
				}
			}

      //..Gap > 1.4
      if (eta > 0.7)
      {
				NtrksAfterGap14P += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap14P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*phi);
						QsinGap14P[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*phi);
					}
				}
      }
          
      if (eta < -0.7)
      {
				NtrksAfterGap14M += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosGap14M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*phi);
						QsinGap14M[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*phi);
					}
				}

      }

			//..3-subevent method
			if(eta < -0.4)
			{//..left part
				NtrksAfter3subL += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosSubLeft[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*phi);
						QsinSubLeft[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*phi);
					}
				}
			}
			if(eta >= -0.4 && eta <= 0.4)
			{//..middle part
				NtrksAfter3subM += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosSubMiddle[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*phi);
						QsinSubMiddle[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*phi);
					}
				}
			}
			if(eta > 0.4)
			{//..right part
				NtrksAfter3subR += 1;
				for(int iharm=0; iharm<10; iharm++)
				{
					for(int ipow=0; ipow<10; ipow++)
					{
						QcosSubRight[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Cos(iharm*phi);
						QsinSubRight[iharm][ipow] += TMath::Power(weight*weightPt, ipow)*TMath::Sin(iharm*phi);
					}
				}
			}


  } // end loop of all track

	hMult->Fill(NtrksAfter, NtrksTrigger);
	hMultTrigger->Fill(NtrksTrigger);

	//............................
	//..GENERIC FRAMEWORK RP
	//............................
 
	//..calculate Q-vector for each harmonics n and power p
	for(int iharm=0; iharm<20; iharm++)
	{
		for(int ipow=0; ipow<20; ipow++)
		{
			Qvector[iharm][ipow] = TComplex(Qcos[iharm][ipow], Qsin[iharm][ipow]);
			QvectorM[iharm][ipow] = TComplex(QcosGapM[iharm][ipow], QsinGapM[iharm][ipow]);
			QvectorP[iharm][ipow] = TComplex(QcosGapP[iharm][ipow], QsinGapP[iharm][ipow]);
			Qvector0M[iharm][ipow] = TComplex(QcosGap0M[iharm][ipow], QsinGap0M[iharm][ipow]);
			Qvector0P[iharm][ipow] = TComplex(QcosGap0P[iharm][ipow], QsinGap0P[iharm][ipow]);
			Qvector2M[iharm][ipow] = TComplex(QcosGap2M[iharm][ipow], QsinGap2M[iharm][ipow]);
			Qvector2P[iharm][ipow] = TComplex(QcosGap2P[iharm][ipow], QsinGap2P[iharm][ipow]);
			Qvector4M[iharm][ipow] = TComplex(QcosGap4M[iharm][ipow], QsinGap4M[iharm][ipow]);
			Qvector4P[iharm][ipow] = TComplex(QcosGap4P[iharm][ipow], QsinGap4P[iharm][ipow]);
			Qvector8M[iharm][ipow] = TComplex(QcosGap8M[iharm][ipow], QsinGap8M[iharm][ipow]);
			Qvector8P[iharm][ipow] = TComplex(QcosGap8P[iharm][ipow], QsinGap8P[iharm][ipow]);
			Qvector14M[iharm][ipow] = TComplex(QcosGap14M[iharm][ipow], QsinGap14M[iharm][ipow]);
			Qvector14P[iharm][ipow] = TComplex(QcosGap14P[iharm][ipow], QsinGap14P[iharm][ipow]);
			QvectorSubLeft[iharm][ipow] = TComplex(QcosSubLeft[iharm][ipow], QsinSubLeft[iharm][ipow]);
			QvectorSubRight[iharm][ipow] = TComplex(QcosSubRight[iharm][ipow], QsinSubRight[iharm][ipow]);
			QvectorSubMiddle[iharm][ipow] = TComplex(QcosSubMiddle[iharm][ipow], QsinSubMiddle[iharm][ipow]);
		}
	}
	
	//..calculate 2-particle correlations
	//..................................
	double Dn2 = Two(0, 0).Re();
	double Dn2Gap = TwoGap(0, 0).Re();
	double Dn2Gap0 = TwoGap0(0, 0).Re();
	double Dn2Gap2 = TwoGap2(0, 0).Re();
	double Dn2Gap4 = TwoGap4(0, 0).Re();
	double Dn2Gap8 = TwoGap8(0, 0).Re();
	double Dn2Gap14 = TwoGap14(0, 0).Re();
	double Dn2_3subLM = Two_3SubLM(0, 0).Re();
	double Dn2_3subRM = Two_3SubRM(0, 0).Re();
	double Dn2_3subLM_perm2 = Two_3SubLM_perm2(0, 0).Re();
	double Dn2_3subLR_perm2 = Two_3SubLR_perm2(0, 0).Re();
	double Dn2_3subRM_perm3 = Two_3SubRM_perm3(0, 0).Re();
	double Dn2_3subRL_perm3 = Two_3SubRL_perm3(0, 0).Re();

	if(NtrksAfter > 1 && Dn2 != 0)
	{	

		//..v2{2} = <cos2(phi1 - phi2)>
		TComplex v22 = Two(2, -2);
		double v22Re = v22.Re()/Dn2;
		fChcn2Ntrks1bin[0]->Fill(NtrksAfter, NtrksTrigger, v22Re, Dn2);
		fcn2Ntrks1bin[0][fBin]->Fill(NtrksAfter, NtrksTrigger, v22Re, Dn2);

		//..v3{2} = <cos3(phi1 - phi2)>
		TComplex v32 = Two(3, -3);
		double v32Re = v32.Re()/Dn2;
		fChcn2Ntrks1bin[1]->Fill(NtrksAfter, NtrksTrigger, v32Re, Dn2);
		fcn2Ntrks1bin[1][fBin]->Fill(NtrksAfter, NtrksTrigger, v32Re, Dn2);

		//..v4{2} = <cos4(phi1 - phi2)>
		TComplex v42 = Two(4, -4);
		double v42Re = v42.Re()/Dn2;
		fChcn2Ntrks1bin[2]->Fill(NtrksAfter, NtrksTrigger, v42Re, Dn2);
		fcn2Ntrks1bin[2][fBin]->Fill(NtrksAfter, NtrksTrigger, v42Re, Dn2);

	}	

	if(NtrksAfterGap0M > 0 && NtrksAfterGap0P > 0 && Dn2Gap0 != 0)
	{
		//..v2{2} with eta Gap > 0.
		TComplex v22Gap0 = TwoGap0(2, -2);
		double v22ReGap0 = v22Gap0.Re()/Dn2Gap0;
		fChcn2Gap0Ntrks1bin[0]->Fill(NtrksAfter, NtrksTrigger, v22ReGap0, Dn2Gap0);
		fcn2Gap0Ntrks1bin[0][fBin]->Fill(NtrksAfter, NtrksTrigger, v22ReGap0, Dn2Gap0);

		//..v3{2} with eta Gap > 0.
		TComplex v32Gap0 = TwoGap0(3, -3);
		double v32ReGap0 = v32Gap0.Re()/Dn2Gap0;
		fChcn2Gap0Ntrks1bin[1]->Fill(NtrksAfter, NtrksTrigger, v32ReGap0, Dn2Gap0);
		fcn2Gap0Ntrks1bin[1][fBin]->Fill(NtrksAfter, NtrksTrigger, v32ReGap0, Dn2Gap0);

		//..v4{2} with eta Gap > 0.
		TComplex v42Gap0 = TwoGap0(4, -4);
		double v42ReGap0 = v42Gap0.Re()/Dn2Gap0;
		fChcn2Gap0Ntrks1bin[2]->Fill(NtrksAfter, NtrksTrigger, v42ReGap0, Dn2Gap0);
		fcn2Gap0Ntrks1bin[2][fBin]->Fill(NtrksAfter, NtrksTrigger, v42ReGap0, Dn2Gap0);

	}

	if(NtrksAfterGap2M > 0 && NtrksAfterGap2P > 0 && Dn2Gap2 != 0)
	{
		//..v2{2} with eta Gap > 0.
		TComplex v22Gap2 = TwoGap2(2, -2);
		double v22ReGap2 = v22Gap2.Re()/Dn2Gap2;
		fChcn2Gap2Ntrks1bin[0]->Fill(NtrksAfter, NtrksTrigger, v22ReGap2, Dn2Gap2);
		fcn2Gap2Ntrks1bin[0][fBin]->Fill(NtrksAfter, NtrksTrigger, v22ReGap2, Dn2Gap2);

		//..v3{2} with eta Gap > 0.
		TComplex v32Gap2 = TwoGap2(3, -3);
		double v32ReGap2 = v32Gap2.Re()/Dn2Gap2;
		fChcn2Gap2Ntrks1bin[1]->Fill(NtrksAfter, NtrksTrigger, v32ReGap2, Dn2Gap2);
		fcn2Gap2Ntrks1bin[1][fBin]->Fill(NtrksAfter, NtrksTrigger, v32ReGap2, Dn2Gap2);

		//..v4{2} with eta Gap > 0.
		TComplex v42Gap2 = TwoGap2(4, -4);
		double v42ReGap2 = v42Gap2.Re()/Dn2Gap2;
		fChcn2Gap2Ntrks1bin[2]->Fill(NtrksAfter, NtrksTrigger, v42ReGap2, Dn2Gap2);
		fcn2Gap2Ntrks1bin[2][fBin]->Fill(NtrksAfter, NtrksTrigger, v42ReGap2, Dn2Gap2);

	}

	if(NtrksAfterGap4M > 0 && NtrksAfterGap4P > 0 && Dn2Gap4 != 0)
	{
		//..v2{2} with eta Gap > 0.4
		TComplex v22Gap4 = TwoGap4(2, -2);
		double v22ReGap4 = v22Gap4.Re()/Dn2Gap4;
		fChcn2Gap4Ntrks1bin[0]->Fill(NtrksAfter, NtrksTrigger, v22ReGap4, Dn2Gap4);
		fcn2Gap4Ntrks1bin[0][fBin]->Fill(NtrksAfter, NtrksTrigger, v22ReGap4, Dn2Gap4);

		//..v3{2} with eta Gap > 0.4
		TComplex v32Gap4 = TwoGap4(3, -3);
		double v32ReGap4 = v32Gap4.Re()/Dn2Gap4;
		fChcn2Gap4Ntrks1bin[1]->Fill(NtrksAfter, NtrksTrigger, v32ReGap4, Dn2Gap4);
		fcn2Gap4Ntrks1bin[1][fBin]->Fill(NtrksAfter, NtrksTrigger, v32ReGap4, Dn2Gap4);

		//..v4{2} with eta Gap > 0.4
		TComplex v42Gap4 = TwoGap4(4, -4);
		double v42ReGap4 = v42Gap4.Re()/Dn2Gap4;
		fChcn2Gap4Ntrks1bin[2]->Fill(NtrksAfter, NtrksTrigger, v42ReGap4, Dn2Gap4);
		fcn2Gap4Ntrks1bin[2][fBin]->Fill(NtrksAfter, NtrksTrigger, v42ReGap4, Dn2Gap4);

	}

	if(NtrksAfterGap8M > 0 && NtrksAfterGap8P > 0 && Dn2Gap8 != 0)
	{
		//..v2{2} with eta Gap > 0.8
		TComplex v22Gap8 = TwoGap8(2, -2);
		double v22ReGap8 = v22Gap8.Re()/Dn2Gap8;
		fChcn2Gap8Ntrks1bin[0]->Fill(NtrksAfter, NtrksTrigger, v22ReGap8, Dn2Gap8);
		fcn2Gap8Ntrks1bin[0][fBin]->Fill(NtrksAfter, NtrksTrigger, v22ReGap8, Dn2Gap8);

		//..v3{2} with eta Gap > 0.8
		TComplex v32Gap8 = TwoGap8(3, -3);
		double v32ReGap8 = v32Gap8.Re()/Dn2Gap8;
		fChcn2Gap8Ntrks1bin[1]->Fill(NtrksAfter, NtrksTrigger, v32ReGap8, Dn2Gap8);
		fcn2Gap8Ntrks1bin[1][fBin]->Fill(NtrksAfter, NtrksTrigger, v32ReGap8, Dn2Gap8);

		//..v4{2} with eta Gap > 0.8
		TComplex v42Gap8 = TwoGap8(4, -4);
		double v42ReGap8 = v42Gap8.Re()/Dn2Gap8;
		fChcn2Gap8Ntrks1bin[2]->Fill(NtrksAfter, NtrksTrigger, v42ReGap8, Dn2Gap8);
		fcn2Gap8Ntrks1bin[2][fBin]->Fill(NtrksAfter, NtrksTrigger, v42ReGap8, Dn2Gap8);

	}

	if(NtrksAfterGapM > 0 && NtrksAfterGapP > 0 && Dn2Gap != 0)
	{
		//..v2{2} with eta Gap > 1.0
		TComplex v22Gap = TwoGap(2, -2);
		double v22ReGap = v22Gap.Re()/Dn2Gap;
		fChcn2GapNtrks1bin[0]->Fill(NtrksAfter, NtrksTrigger, v22ReGap, Dn2Gap);
		fcn2GapNtrks1bin[0][fBin]->Fill(NtrksAfter, NtrksTrigger, v22ReGap, Dn2Gap);

		//..v3{2} with eta Gap > 1.0
		TComplex v32Gap = TwoGap(3, -3);
		double v32ReGap = v32Gap.Re()/Dn2Gap;
		fChcn2GapNtrks1bin[1]->Fill(NtrksAfter, NtrksTrigger, v32ReGap, Dn2Gap);
		fcn2GapNtrks1bin[1][fBin]->Fill(NtrksAfter, NtrksTrigger, v32ReGap, Dn2Gap);

		//..v4{2} with eta Gap > 1.0
		TComplex v42Gap = TwoGap(4, -4);
		double v42ReGap = v42Gap.Re()/Dn2Gap;
		fChcn2GapNtrks1bin[2]->Fill(NtrksAfter, NtrksTrigger, v42ReGap, Dn2Gap);
		fcn2GapNtrks1bin[2][fBin]->Fill(NtrksAfter, NtrksTrigger, v42ReGap, Dn2Gap);

	}

	if(NtrksAfterGap14M > 0 && NtrksAfterGap14P > 0 && Dn2Gap14 != 0)
	{
		//..v2{2} with eta Gap > 1.4
		TComplex v22Gap14 = TwoGap14(2, -2);
		double v22ReGap14 = v22Gap14.Re()/Dn2Gap14;
		fChcn2Gap14Ntrks1bin[0]->Fill(NtrksAfter, NtrksTrigger, v22ReGap14, Dn2Gap14);
		fcn2Gap14Ntrks1bin[0][fBin]->Fill(NtrksAfter, NtrksTrigger, v22ReGap14, Dn2Gap14);

		//..v3{2} with eta Gap > 1.4
		TComplex v32Gap14 = TwoGap14(3, -3);
		double v32ReGap14 = v32Gap14.Re()/Dn2Gap14;
		fChcn2Gap14Ntrks1bin[1]->Fill(NtrksAfter, NtrksTrigger, v32ReGap14, Dn2Gap14);
		fcn2Gap14Ntrks1bin[1][fBin]->Fill(NtrksAfter, NtrksTrigger, v32ReGap14, Dn2Gap14);

		//..v4{2} with eta Gap > 1.4
		TComplex v42Gap14 = TwoGap14(4, -4);
		double v42ReGap14 = v42Gap14.Re()/Dn2Gap14;
		fChcn2Gap14Ntrks1bin[2]->Fill(NtrksAfter, NtrksTrigger, v42ReGap14, Dn2Gap14);
		fcn2Gap14Ntrks1bin[2][fBin]->Fill(NtrksAfter, NtrksTrigger, v42ReGap14, Dn2Gap14);

	}


	//..for 3-subevent method, Gap0
	if(NtrksAfter3subL > 0 && NtrksAfter3subM > 0 && Dn2_3subLM != 0)
	{//..left+middle
		TComplex v22_3subLM = Two_3SubLM(2, -2);
		double v22Re_3subLM = v22_3subLM.Re()/Dn2_3subLM;
		fChcn2_3subLMNtrks1bin[0]->Fill(NtrksAfter, NtrksTrigger, v22Re_3subLM, Dn2_3subLM);
		fcn2_3subLMNtrks1bin[0][fBin]->Fill(NtrksAfter, NtrksTrigger, v22Re_3subLM, Dn2_3subLM);

		TComplex v32_3subLM = Two_3SubLM(3, -3);
		double v32Re_3subLM = v32_3subLM.Re()/Dn2_3subLM;
		fChcn2_3subLMNtrks1bin[1]->Fill(NtrksAfter, NtrksTrigger, v32Re_3subLM, Dn2_3subLM);
		fcn2_3subLMNtrks1bin[1][fBin]->Fill(NtrksAfter, NtrksTrigger, v32Re_3subLM, Dn2_3subLM);

		TComplex v42_3subLM = Two_3SubLM(4, -4);
		double v42Re_3subLM = v42_3subLM.Re()/Dn2_3subLM;
		fChcn2_3subLMNtrks1bin[2]->Fill(NtrksAfter, NtrksTrigger, v42Re_3subLM, Dn2_3subLM);
		fcn2_3subLMNtrks1bin[2][fBin]->Fill(NtrksAfter, NtrksTrigger, v42Re_3subLM, Dn2_3subLM);

	}

	if(NtrksAfter3subM > 0 && NtrksAfter3subR > 0 && Dn2_3subRM != 0)
	{//..right+middle
		TComplex v22_3subRM = Two_3SubRM(2, -2);
		double v22Re_3subRM = v22_3subRM.Re()/Dn2_3subRM;
		fChcn2_3subRMNtrks1bin[0]->Fill(NtrksAfter, NtrksTrigger, v22Re_3subRM, Dn2_3subRM);
		fcn2_3subRMNtrks1bin[0][fBin]->Fill(NtrksAfter, NtrksTrigger, v22Re_3subRM, Dn2_3subRM);

		TComplex v32_3subRM = Two_3SubRM(3, -3);
		double v32Re_3subRM = v32_3subRM.Re()/Dn2_3subRM;
		fChcn2_3subRMNtrks1bin[1]->Fill(NtrksAfter, NtrksTrigger, v32Re_3subRM, Dn2_3subRM);
		fcn2_3subRMNtrks1bin[1][fBin]->Fill(NtrksAfter, NtrksTrigger, v32Re_3subRM, Dn2_3subRM);

		TComplex v42_3subRM = Two_3SubRM(4, -4);
		double v42Re_3subRM = v42_3subRM.Re()/Dn2_3subRM;
		fChcn2_3subRMNtrks1bin[2]->Fill(NtrksAfter, NtrksTrigger, v42Re_3subRM, Dn2_3subRM);
		fcn2_3subRMNtrks1bin[2][fBin]->Fill(NtrksAfter, NtrksTrigger, v42Re_3subRM, Dn2_3subRM);

	}

	//..for 3-subevent method, perm 2
	if(NtrksAfter3subL > 0 && NtrksAfter3subM > 0 && Dn2_3subLM_perm2 != 0)
	{//..left+middle
		TComplex v22_3subLM_perm2 = Two_3SubLM_perm2(2, -2);
		double v22Re_3subLM_perm2 = v22_3subLM_perm2.Re()/Dn2_3subLM_perm2;
		fChcn2_3subLM_perm2Ntrks1bin[0]->Fill(NtrksAfter, NtrksTrigger, v22Re_3subLM_perm2, Dn2_3subLM_perm2);
		fcn2_3subLM_perm2Ntrks1bin[0][fBin]->Fill(NtrksAfter, NtrksTrigger, v22Re_3subLM_perm2, Dn2_3subLM_perm2);

		TComplex v32_3subLM_perm2 = Two_3SubLM_perm2(3, -3);
		double v32Re_3subLM_perm2 = v32_3subLM_perm2.Re()/Dn2_3subLM_perm2;
		fChcn2_3subLM_perm2Ntrks1bin[1]->Fill(NtrksAfter, NtrksTrigger, v32Re_3subLM_perm2, Dn2_3subLM_perm2);
		fcn2_3subLM_perm2Ntrks1bin[1][fBin]->Fill(NtrksAfter, NtrksTrigger, v32Re_3subLM_perm2, Dn2_3subLM_perm2);

		TComplex v42_3subLM_perm2 = Two_3SubLM_perm2(4, -4);
		double v42Re_3subLM_perm2 = v42_3subLM_perm2.Re()/Dn2_3subLM_perm2;
		fChcn2_3subLM_perm2Ntrks1bin[2]->Fill(NtrksAfter, NtrksTrigger, v42Re_3subLM_perm2, Dn2_3subLM_perm2);
		fcn2_3subLM_perm2Ntrks1bin[2][fBin]->Fill(NtrksAfter, NtrksTrigger, v42Re_3subLM_perm2, Dn2_3subLM_perm2);

	}

	if(NtrksAfter3subL > 0 && NtrksAfter3subR > 0 && Dn2_3subLR_perm2 != 0)
	{//..left+right
		TComplex v22_3subLR_perm2 = Two_3SubLR_perm2(2, -2);
		double v22Re_3subLR_perm2 = v22_3subLR_perm2.Re()/Dn2_3subLR_perm2;
		fChcn2_3subLR_perm2Ntrks1bin[0]->Fill(NtrksAfter, NtrksTrigger, v22Re_3subLR_perm2, Dn2_3subLR_perm2);
		fcn2_3subLR_perm2Ntrks1bin[0][fBin]->Fill(NtrksAfter, NtrksTrigger, v22Re_3subLR_perm2, Dn2_3subLR_perm2);

		TComplex v32_3subLR_perm2 = Two_3SubLR_perm2(3, -3);
		double v32Re_3subLR_perm2 = v32_3subLR_perm2.Re()/Dn2_3subLR_perm2;
		fChcn2_3subLR_perm2Ntrks1bin[1]->Fill(NtrksAfter, NtrksTrigger, v32Re_3subLR_perm2, Dn2_3subLR_perm2);
		fcn2_3subLR_perm2Ntrks1bin[1][fBin]->Fill(NtrksAfter, NtrksTrigger, v32Re_3subLR_perm2, Dn2_3subLR_perm2);

		TComplex v42_3subLR_perm2 = Two_3SubLR_perm2(4, -4);
		double v42Re_3subLR_perm2 = v42_3subLR_perm2.Re()/Dn2_3subLR_perm2;
		fChcn2_3subLR_perm2Ntrks1bin[2]->Fill(NtrksAfter, NtrksTrigger, v42Re_3subLR_perm2, Dn2_3subLR_perm2);
		fcn2_3subLR_perm2Ntrks1bin[2][fBin]->Fill(NtrksAfter, NtrksTrigger, v42Re_3subLR_perm2, Dn2_3subLR_perm2);

	}

	//..for 3-subevent method, perm 3
	if(NtrksAfter3subL > 0 && NtrksAfter3subR > 0 && Dn2_3subRL_perm3 != 0)
	{//..left+right
		TComplex v22_3subRL_perm3 = Two_3SubRL_perm3(2, -2);
		double v22Re_3subRL_perm3 = v22_3subRL_perm3.Re()/Dn2_3subRL_perm3;
		fChcn2_3subRL_perm3Ntrks1bin[0]->Fill(NtrksAfter, NtrksTrigger, v22Re_3subRL_perm3, Dn2_3subRL_perm3);
		fcn2_3subRL_perm3Ntrks1bin[0][fBin]->Fill(NtrksAfter, NtrksTrigger, v22Re_3subRL_perm3, Dn2_3subRL_perm3);

		TComplex v32_3subRL_perm3 = Two_3SubRL_perm3(3, -3);
		double v32Re_3subRL_perm3 = v32_3subRL_perm3.Re()/Dn2_3subRL_perm3;
		fChcn2_3subRL_perm3Ntrks1bin[1]->Fill(NtrksAfter, NtrksTrigger, v32Re_3subRL_perm3, Dn2_3subRL_perm3);
		fcn2_3subRL_perm3Ntrks1bin[1][fBin]->Fill(NtrksAfter, NtrksTrigger, v32Re_3subRL_perm3, Dn2_3subRL_perm3);

		TComplex v42_3subRL_perm3 = Two_3SubRL_perm3(4, -4);
		double v42Re_3subRL_perm3 = v42_3subRL_perm3.Re()/Dn2_3subRL_perm3;
		fChcn2_3subRL_perm3Ntrks1bin[2]->Fill(NtrksAfter, NtrksTrigger, v42Re_3subRL_perm3, Dn2_3subRL_perm3);
		fcn2_3subRL_perm3Ntrks1bin[2][fBin]->Fill(NtrksAfter, NtrksTrigger, v42Re_3subRL_perm3, Dn2_3subRL_perm3);

	}

	if(NtrksAfter3subM > 0 && NtrksAfter3subR > 0 && Dn2_3subRM_perm3 != 0)
	{//..middle+right
		TComplex v22_3subRM_perm3 = Two_3SubRM_perm3(2, -2);
		double v22Re_3subRM_perm3 = v22_3subRM_perm3.Re()/Dn2_3subRM_perm3;
		fChcn2_3subRM_perm3Ntrks1bin[0]->Fill(NtrksAfter, NtrksTrigger, v22Re_3subRM_perm3, Dn2_3subRM_perm3);
		fcn2_3subRM_perm3Ntrks1bin[0][fBin]->Fill(NtrksAfter, NtrksTrigger, v22Re_3subRM_perm3, Dn2_3subRM_perm3);

		TComplex v32_3subRM_perm3 = Two_3SubRM_perm3(3, -3);
		double v32Re_3subRM_perm3 = v32_3subRM_perm3.Re()/Dn2_3subRM_perm3;
		fChcn2_3subRM_perm3Ntrks1bin[1]->Fill(NtrksAfter, NtrksTrigger, v32Re_3subRM_perm3, Dn2_3subRM_perm3);
		fcn2_3subRM_perm3Ntrks1bin[1][fBin]->Fill(NtrksAfter, NtrksTrigger, v32Re_3subRM_perm3, Dn2_3subRM_perm3);

		TComplex v42_3subRM_perm3 = Two_3SubRM_perm3(4, -4);
		double v42Re_3subRM_perm3 = v42_3subRM_perm3.Re()/Dn2_3subRM_perm3;
		fChcn2_3subRM_perm3Ntrks1bin[2]->Fill(NtrksAfter, NtrksTrigger, v42Re_3subRM_perm3, Dn2_3subRM_perm3);
		fcn2_3subRM_perm3Ntrks1bin[2][fBin]->Fill(NtrksAfter, NtrksTrigger, v42Re_3subRM_perm3, Dn2_3subRM_perm3);

	}


	//..calculate 4-particle correlations
	//................................
	double Dn4 = Four(0, 0, 0, 0).Re();
	double Dn4Gap0 = FourGap0(0, 0, 0, 0).Re();
	double Dn4_3sub = Four_3SubEvts(0, 0, 0, 0).Re();
	double Dn4_3sub_perm2 = Four_3SubEvts_perm2(0, 0, 0, 0).Re();
	double Dn4_3sub_perm3 = Four_3SubEvts_perm3(0, 0, 0, 0).Re();

	if(NtrksAfter > 3 && Dn4 != 0)
	{
		
		TComplex v24 = Four(2, 2, -2, -2);
		double v24Re = v24.Re()/Dn4;
		fChcn4Ntrks1bin[0]->Fill(NtrksAfter, NtrksTrigger, v24Re, Dn4);
		fcn4Ntrks1bin[0][fBin]->Fill(NtrksAfter, NtrksTrigger, v24Re, Dn4);

		TComplex v34 = Four(3, 3, -3, -3);
		double v34Re = v34.Re()/Dn4;
		fChcn4Ntrks1bin[1]->Fill(NtrksAfter, NtrksTrigger, v34Re, Dn4);
		fcn4Ntrks1bin[1][fBin]->Fill(NtrksAfter, NtrksTrigger, v34Re, Dn4);

		TComplex v44 = Four(4, 4, -4, -4);
		double v44Re = v44.Re()/Dn4;
		fChcn4Ntrks1bin[2]->Fill(NtrksAfter, NtrksTrigger, v44Re, Dn4);
		fcn4Ntrks1bin[2][fBin]->Fill(NtrksAfter, NtrksTrigger, v44Re, Dn4);

	}

	if(NtrksAfterGap0M > 1 && NtrksAfterGap0P > 1 && Dn4Gap0 !=0)
	{

		TComplex v24Gap0 = FourGap0(2, 2, -2, -2);
		double v24Gap0Re = v24Gap0.Re()/Dn4Gap0;
		fChcn4Gap0Ntrks1bin[0]->Fill(NtrksAfter, NtrksTrigger, v24Gap0Re, Dn4Gap0);
		fcn4Gap0Ntrks1bin[0][fBin]->Fill(NtrksAfter, NtrksTrigger, v24Gap0Re, Dn4Gap0);

		TComplex v34Gap0 = FourGap0(3, 3, -3, -3);
		double v34Gap0Re = v34Gap0.Re()/Dn4Gap0;
		fChcn4Gap0Ntrks1bin[1]->Fill(NtrksAfter, NtrksTrigger, v34Gap0Re, Dn4Gap0);
		fcn4Gap0Ntrks1bin[1][fBin]->Fill(NtrksAfter, NtrksTrigger, v34Gap0Re, Dn4Gap0);

		TComplex v44Gap0 = FourGap0(4, 4, -4, -4);
		double v44Gap0Re = v44Gap0.Re()/Dn4Gap0;
		fChcn4Gap0Ntrks1bin[2]->Fill(NtrksAfter, NtrksTrigger, v44Gap0Re, Dn4Gap0);
		fcn4Gap0Ntrks1bin[2][fBin]->Fill(NtrksAfter, NtrksTrigger, v44Gap0Re, Dn4Gap0);

	}

	//..3-subevent method
	if(NtrksAfter3subL > 0 && NtrksAfter3subM > 1 && NtrksAfter3subR > 0 && Dn4_3sub != 0)
	{
		TComplex v24_3sub = Four_3SubEvts(2, 2, -2, -2);
		double v24_3subRe = v24_3sub.Re()/Dn4_3sub;
		fChcn4_3subNtrks1bin[0]->Fill(NtrksAfter, NtrksTrigger, v24_3subRe, Dn4_3sub);
		fcn4_3subNtrks1bin[0][fBin]->Fill(NtrksAfter, NtrksTrigger, v24_3subRe, Dn4_3sub);

		TComplex v34_3sub = Four_3SubEvts(3, 3, -3, -3);
		double v34_3subRe = v34_3sub.Re()/Dn4_3sub;
		fChcn4_3subNtrks1bin[1]->Fill(NtrksAfter, NtrksTrigger, v34_3subRe, Dn4_3sub);
		fcn4_3subNtrks1bin[1][fBin]->Fill(NtrksAfter, NtrksTrigger, v34_3subRe, Dn4_3sub);

		TComplex v44_3sub = Four_3SubEvts(4, 4, -4, -4);
		double v44_3subRe = v44_3sub.Re()/Dn4_3sub;
		fChcn4_3subNtrks1bin[2]->Fill(NtrksAfter, NtrksTrigger, v44_3subRe, Dn4_3sub);
		fcn4_3subNtrks1bin[2][fBin]->Fill(NtrksAfter, NtrksTrigger, v44_3subRe, Dn4_3sub);

	}

	//	c24_3sub perm 2
	if(NtrksAfter3subL > 1 && NtrksAfter3subM > 0 && NtrksAfter3subR > 0 && Dn4_3sub_perm2 != 0)
	{
		TComplex v24_3sub_perm2 = Four_3SubEvts_perm2(2, 2, -2, -2);
		double v24_3sub_perm2Re = v24_3sub_perm2.Re()/Dn4_3sub_perm2;
		fChcn4_3sub_perm2Ntrks1bin[0]->Fill(NtrksAfter, NtrksTrigger, v24_3sub_perm2Re, Dn4_3sub_perm2);
		fcn4_3sub_perm2Ntrks1bin[0][fBin]->Fill(NtrksAfter, NtrksTrigger, v24_3sub_perm2Re, Dn4_3sub_perm2);

		TComplex v34_3sub_perm2 = Four_3SubEvts_perm2(3, 3, -3, -3);
		double v34_3sub_perm2Re = v34_3sub_perm2.Re()/Dn4_3sub_perm2;
		fChcn4_3sub_perm2Ntrks1bin[1]->Fill(NtrksAfter, NtrksTrigger, v34_3sub_perm2Re, Dn4_3sub_perm2);
		fcn4_3sub_perm2Ntrks1bin[1][fBin]->Fill(NtrksAfter, NtrksTrigger, v34_3sub_perm2Re, Dn4_3sub_perm2);

		TComplex v44_3sub_perm2 = Four_3SubEvts_perm2(4, 4, -4, -4);
		double v44_3sub_perm2Re = v44_3sub_perm2.Re()/Dn4_3sub_perm2;
		fChcn4_3sub_perm2Ntrks1bin[2]->Fill(NtrksAfter, NtrksTrigger, v44_3sub_perm2Re, Dn4_3sub_perm2);
		fcn4_3sub_perm2Ntrks1bin[2][fBin]->Fill(NtrksAfter, NtrksTrigger, v44_3sub_perm2Re, Dn4_3sub_perm2);

	}

	//	c24_3sub perm 2
	if(NtrksAfter3subL > 0 && NtrksAfter3subM > 0 && NtrksAfter3subR > 1 && Dn4_3sub_perm3 != 0)
	{
		TComplex v24_3sub_perm3 = Four_3SubEvts_perm3(2, 2, -2, -2);
		double v24_3sub_perm3Re = v24_3sub_perm3.Re()/Dn4_3sub_perm3;
		fChcn4_3sub_perm3Ntrks1bin[0]->Fill(NtrksAfter, NtrksTrigger, v24_3sub_perm3Re, Dn4_3sub_perm3);
		fcn4_3sub_perm3Ntrks1bin[0][fBin]->Fill(NtrksAfter, NtrksTrigger, v24_3sub_perm3Re, Dn4_3sub_perm3);

		TComplex v34_3sub_perm3 = Four_3SubEvts_perm3(3, 3, -3, -3);
		double v34_3sub_perm3Re = v34_3sub_perm3.Re()/Dn4_3sub_perm3;
		fChcn4_3sub_perm3Ntrks1bin[1]->Fill(NtrksAfter, NtrksTrigger, v34_3sub_perm3Re, Dn4_3sub_perm3);
		fcn4_3sub_perm3Ntrks1bin[1][fBin]->Fill(NtrksAfter, NtrksTrigger, v34_3sub_perm3Re, Dn4_3sub_perm3);

		TComplex v44_3sub_perm3 = Four_3SubEvts_perm3(4, 4, -4, -4);
		double v44_3sub_perm3Re = v44_3sub_perm3.Re()/Dn4_3sub_perm3;
		fChcn4_3sub_perm3Ntrks1bin[2]->Fill(NtrksAfter, NtrksTrigger, v44_3sub_perm3Re, Dn4_3sub_perm3);
		fcn4_3sub_perm3Ntrks1bin[2][fBin]->Fill(NtrksAfter, NtrksTrigger, v44_3sub_perm3Re, Dn4_3sub_perm3);

	}


	//..SC(n,m)
	//................................
	if(NtrksAfter > 3 && Dn4 != 0)
	{
		//..SC(4,2,-4,-2)	
		TComplex sc4242 = Four(4, 2, -4, -2);
		double sc4242Re = sc4242.Re()/Dn4;
		fChsc4242->Fill(NtrksAfter, NtrksTrigger, sc4242Re, Dn4);
		fsc4242[fBin]->Fill(NtrksAfter, NtrksTrigger, sc4242Re, Dn4);

		//..SC(3,2,-3,-2)	
		TComplex sc3232 = Four(3, 2, -3, -2);
		double sc3232Re = sc3232.Re()/Dn4;
		fChsc3232->Fill(NtrksAfter, NtrksTrigger, sc3232Re, Dn4);
		fsc3232[fBin]->Fill(NtrksAfter, NtrksTrigger, sc3232Re, Dn4);

	}	

	if(NtrksAfterGap0M > 1 && NtrksAfterGap0P > 1 && Dn4Gap0 != 0)
	{

		TComplex sc4242Gap0 = FourGap0(4, 2, -4, -2);
		double sc4242Gap0Re = sc4242Gap0.Re()/Dn4Gap0;
		fChsc4242Gap0->Fill(NtrksAfter, NtrksTrigger, sc4242Gap0Re, Dn4Gap0);
		fsc4242Gap0[fBin]->Fill(NtrksAfter, NtrksTrigger, sc4242Gap0Re, Dn4Gap0);

		TComplex sc3232Gap0 = FourGap0(3, 2, -3, -2);
		double sc3232Gap0Re = sc3232Gap0.Re()/Dn4Gap0;
		fChsc3232Gap0->Fill(NtrksAfter, NtrksTrigger, sc3232Gap0Re, Dn4Gap0);
		fsc3232Gap0[fBin]->Fill(NtrksAfter, NtrksTrigger, sc3232Gap0Re, Dn4Gap0);

	}

	if(NtrksAfter3subL > 0 && NtrksAfter3subM > 1 && NtrksAfter3subR > 0 && Dn4_3sub != 0)
	{

		//..variant 1 (4, 2, -4, -2)
		TComplex sc4242_3sub = Four_3SubEvts(4, 2, -4, -2);
		double sc4242_3subRe = sc4242_3sub.Re()/Dn4_3sub;
		fChsc4242_3sub->Fill(NtrksAfter, NtrksTrigger, sc4242_3subRe, Dn4_3sub);
		fsc4242_3sub[fBin]->Fill(NtrksAfter, NtrksTrigger, sc4242_3subRe, Dn4_3sub);

		TComplex sc3232_3sub = Four_3SubEvts(3, 2, -3, -2);
		double sc3232_3subRe = sc3232_3sub.Re()/Dn4_3sub;
		fChsc3232_3sub->Fill(NtrksAfter, NtrksTrigger, sc3232_3subRe, Dn4_3sub);
		fsc3232_3sub[fBin]->Fill(NtrksAfter, NtrksTrigger, sc3232_3subRe, Dn4_3sub);

		//..variant 2 (4, 2, -2, -4)
		TComplex sc4224_3sub = Four_3SubEvts(4, 2, -2, -4);
		double sc4224_3subRe = sc4224_3sub.Re()/Dn4_3sub;
		fChsc4224_3sub->Fill(NtrksAfter, NtrksTrigger, sc4224_3subRe, Dn4_3sub);
		fsc4224_3sub[fBin]->Fill(NtrksAfter, NtrksTrigger, sc4224_3subRe, Dn4_3sub);

		TComplex sc3223_3sub = Four_3SubEvts(3, 2, -2, -3);
		double sc3223_3subRe = sc3223_3sub.Re()/Dn4_3sub;
		fChsc3223_3sub->Fill(NtrksAfter, NtrksTrigger, sc3223_3subRe, Dn4_3sub);
		fsc3223_3sub[fBin]->Fill(NtrksAfter, NtrksTrigger, sc3223_3subRe, Dn4_3sub);

	}

	//	perm 2
	if(NtrksAfter3subL > 1 && NtrksAfter3subM > 0 && NtrksAfter3subR > 0 && Dn4_3sub_perm2 != 0)
	{

		//..variant 1 (4, 2, -4, -2)
		TComplex sc4242_3sub_perm2 = Four_3SubEvts_perm2(4, 2, -4, -2);
		double sc4242_3subRe_perm2 = sc4242_3sub_perm2.Re()/Dn4_3sub_perm2;
		fChsc4242_3sub_perm2->Fill(NtrksAfter, NtrksTrigger, sc4242_3subRe_perm2, Dn4_3sub_perm2);
		fsc4242_3sub_perm2[fBin]->Fill(NtrksAfter, NtrksTrigger, sc4242_3subRe_perm2, Dn4_3sub_perm2);

		TComplex sc3232_3sub_perm2 = Four_3SubEvts_perm2(3, 2, -3, -2);
		double sc3232_3subRe_perm2 = sc3232_3sub_perm2.Re()/Dn4_3sub_perm2;
		fChsc3232_3sub_perm2->Fill(NtrksAfter, NtrksTrigger, sc3232_3subRe_perm2, Dn4_3sub_perm2);
		fsc3232_3sub_perm2[fBin]->Fill(NtrksAfter, NtrksTrigger, sc3232_3subRe_perm2, Dn4_3sub_perm2);

		//..variant 2 (4, 2, -2, -4)
		TComplex sc4224_3sub_perm2 = Four_3SubEvts_perm2(4, 2, -2, -4);
		double sc4224_3subRe_perm2 = sc4224_3sub_perm2.Re()/Dn4_3sub_perm2;
		fChsc4224_3sub_perm2->Fill(NtrksAfter, NtrksTrigger, sc4224_3subRe_perm2, Dn4_3sub_perm2);
		fsc4224_3sub_perm2[fBin]->Fill(NtrksAfter, NtrksTrigger, sc4224_3subRe_perm2, Dn4_3sub_perm2);

		TComplex sc3223_3sub_perm2 = Four_3SubEvts_perm2(3, 2, -2, -3);
		double sc3223_3subRe_perm2 = sc3223_3sub_perm2.Re()/Dn4_3sub_perm2;
		fChsc3223_3sub_perm2->Fill(NtrksAfter, NtrksTrigger, sc3223_3subRe_perm2, Dn4_3sub_perm2);
		fsc3223_3sub_perm2[fBin]->Fill(NtrksAfter, NtrksTrigger, sc3223_3subRe_perm2, Dn4_3sub_perm2);

	}

	//	perm 3
	if(NtrksAfter3subL > 0 && NtrksAfter3subM > 0 && NtrksAfter3subR > 1 && Dn4_3sub_perm3 != 0)
	{

		//..variant 1 (4, 2, -4, -2)
		TComplex sc4242_3sub_perm3 = Four_3SubEvts_perm3(4, 2, -4, -2);
		double sc4242_3subRe_perm3 = sc4242_3sub_perm3.Re()/Dn4_3sub_perm3;
		fChsc4242_3sub_perm3->Fill(NtrksAfter, NtrksTrigger, sc4242_3subRe_perm3, Dn4_3sub_perm3);
		fsc4242_3sub_perm3[fBin]->Fill(NtrksAfter, NtrksTrigger, sc4242_3subRe_perm3, Dn4_3sub_perm3);

		TComplex sc3232_3sub_perm3 = Four_3SubEvts_perm3(3, 2, -3, -2);
		double sc3232_3subRe_perm3 = sc3232_3sub_perm3.Re()/Dn4_3sub_perm3;
		fChsc3232_3sub_perm3->Fill(NtrksAfter, NtrksTrigger, sc3232_3subRe_perm3, Dn4_3sub_perm3);
		fsc3232_3sub_perm3[fBin]->Fill(NtrksAfter, NtrksTrigger, sc3232_3subRe_perm3, Dn4_3sub_perm3);

		//..variant 2 (4, 2, -2, -4)
		TComplex sc4224_3sub_perm3 = Four_3SubEvts_perm3(4, 2, -2, -4);
		double sc4224_3subRe_perm3 = sc4224_3sub_perm3.Re()/Dn4_3sub_perm3;
		fChsc4224_3sub_perm3->Fill(NtrksAfter, NtrksTrigger, sc4224_3subRe_perm3, Dn4_3sub_perm3);
		fsc4224_3sub_perm3[fBin]->Fill(NtrksAfter, NtrksTrigger, sc4224_3subRe_perm3, Dn4_3sub_perm3);

		TComplex sc3223_3sub_perm3 = Four_3SubEvts_perm3(3, 2, -2, -3);
		double sc3223_3subRe_perm3 = sc3223_3sub_perm3.Re()/Dn4_3sub_perm3;
		fChsc3223_3sub_perm3->Fill(NtrksAfter, NtrksTrigger, sc3223_3subRe_perm3, Dn4_3sub_perm3);
		fsc3223_3sub_perm3[fBin]->Fill(NtrksAfter, NtrksTrigger, sc3223_3subRe_perm3, Dn4_3sub_perm3);

	}

	//..calculate 6-particle correlations	
	//...................................
	double Dn6 = Six(0, 0, 0, 0, 0, 0).Re();
	double Dn6Gap0 = SixGap0(0, 0, 0, 0, 0, 0).Re();

	if(NtrksAfter > 5 && Dn6 != 0)
	{

		TComplex v26 = Six(2, 2, 2, -2, -2, -2);	
		double v26Re = v26.Re()/Dn6;
		fChcn6Ntrks1bin[0]->Fill(NtrksAfter, NtrksTrigger, v26Re, Dn6);
		fcn6Ntrks1bin[0][fBin]->Fill(NtrksAfter, NtrksTrigger, v26Re, Dn6);

		TComplex v36 = Six(3, 3, 3, -3, -3, -3);
		double v36Re = v36.Re()/Dn6;
		fChcn6Ntrks1bin[1]->Fill(NtrksAfter, NtrksTrigger, v36Re, Dn6);
		fcn6Ntrks1bin[1][fBin]->Fill(NtrksAfter, NtrksTrigger, v36Re, Dn6);

		TComplex v46 = Six(4, 4, 4, -4, -4, -4);
		double v46Re = v46.Re()/Dn6;
		fChcn6Ntrks1bin[2]->Fill(NtrksAfter, NtrksTrigger, v46Re, Dn6);
		fcn6Ntrks1bin[2][fBin]->Fill(NtrksAfter, NtrksTrigger, v46Re, Dn6);

	}

	if(NtrksAfterGap0M > 2 && NtrksAfterGap0P > 2 && Dn6Gap0 != 0)
	{

		TComplex v26Gap0 = SixGap0(2, 2, 2, -2, -2, -2);	
		double v26Gap0Re = v26Gap0.Re()/Dn6Gap0;
		fChcn6Gap0Ntrks1bin[0]->Fill(NtrksAfter, NtrksTrigger, v26Gap0Re, Dn6Gap0);
		fcn6Gap0Ntrks1bin[0][fBin]->Fill(NtrksAfter, NtrksTrigger, v26Gap0Re, Dn6Gap0);

		TComplex v36Gap0 = SixGap0(3, 3, 3, -3, -3, -3);
		double v36Gap0Re = v36Gap0.Re()/Dn6Gap0;
		fChcn6Gap0Ntrks1bin[1]->Fill(NtrksAfter, NtrksTrigger, v36Gap0Re, Dn6Gap0);
		fcn6Gap0Ntrks1bin[1][fBin]->Fill(NtrksAfter, NtrksTrigger, v36Gap0Re, Dn6Gap0);

		TComplex v46Gap0 = SixGap0(4, 4, 4, -4, -4, -4);
		double v46Gap0Re = v46Gap0.Re()/Dn6Gap0;
		fChcn6Gap0Ntrks1bin[2]->Fill(NtrksAfter, NtrksTrigger, v46Gap0Re, Dn6Gap0);
		fcn6Gap0Ntrks1bin[2][fBin]->Fill(NtrksAfter, NtrksTrigger, v46Gap0Re, Dn6Gap0);

	}

  //..calculate 8-particle correlations 
  //...................................
	double Dn8 = Eight(0, 0, 0, 0, 0, 0, 0, 0).Re();
  double Dn8Gap0 = EightGap0(0, 0, 0, 0, 0, 0, 0, 0).Re();

  if(NtrksAfter > 7 && Dn8 != 0)                                                                                                                                                                    
  {

    TComplex v28 = Eight(2, 2, 2, 2, -2, -2, -2, -2); 
    double v28Re = v28.Re()/Dn8;
    fChcn8Ntrks1bin[0]->Fill(NtrksAfter, NtrksTrigger, v28Re, Dn8);
    fcn8Ntrks1bin[0][fBin]->Fill(NtrksAfter, NtrksTrigger, v28Re, Dn8);

    TComplex v38 = Eight(3, 3, 3, 3, -3, -3, -3, -3);
    double v38Re = v38.Re()/Dn8;
    fChcn8Ntrks1bin[1]->Fill(NtrksAfter, NtrksTrigger, v38Re, Dn8);
    fcn8Ntrks1bin[1][fBin]->Fill(NtrksAfter, NtrksTrigger, v38Re, Dn8);

    TComplex v48 = Eight(4, 4, 4, 4, -4, -4, -4, -4);
    double v48Re = v48.Re()/Dn8;
    fChcn8Ntrks1bin[2]->Fill(NtrksAfter, NtrksTrigger, v48Re, Dn8);
    fcn8Ntrks1bin[2][fBin]->Fill(NtrksAfter, NtrksTrigger, v48Re, Dn8);

  } 

  if(NtrksAfterGap0M > 3 && NtrksAfterGap0P > 3 && Dn8Gap0 != 0)                                                                                                                                                                    
  {

    TComplex v28Gap0 = EightGap0(2, 2, 2, 2, -2, -2, -2, -2); 
    double v28Gap0Re = v28Gap0.Re()/Dn8Gap0;
    fChcn8Gap0Ntrks1bin[0]->Fill(NtrksAfter, NtrksTrigger, v28Gap0Re, Dn8Gap0);
    fcn8Gap0Ntrks1bin[0][fBin]->Fill(NtrksAfter, NtrksTrigger, v28Gap0Re, Dn8Gap0);

    TComplex v38Gap0 = EightGap0(3, 3, 3, 3, -3, -3, -3, -3);
    double v38Gap0Re = v38Gap0.Re()/Dn8Gap0;
    fChcn8Gap0Ntrks1bin[1]->Fill(NtrksAfter, NtrksTrigger, v38Gap0Re, Dn8Gap0);
    fcn8Gap0Ntrks1bin[1][fBin]->Fill(NtrksAfter, NtrksTrigger, v38Gap0Re, Dn8Gap0);

    TComplex v48Gap0 = EightGap0(4, 4, 4, 4, -4, -4, -4, -4);
    double v48Gap0Re = v48Gap0.Re()/Dn8Gap0;
    fChcn8Gap0Ntrks1bin[2]->Fill(NtrksAfter, NtrksTrigger, v48Gap0Re, Dn8Gap0);
    fcn8Gap0Ntrks1bin[2][fBin]->Fill(NtrksAfter, NtrksTrigger, v48Gap0Re, Dn8Gap0);

  } 


}
//_____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::Q(int n, int p)
{

	if(n>=0) return Qvector[n][p];
  else return TComplex::Conjugate(Qvector[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::QGapM(int n, int p)
{

	if(n>=0) return QvectorM[n][p];
  else return TComplex::Conjugate(QvectorM[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::QGapP(int n, int p)
{

	if(n>=0) return QvectorP[n][p];
  else return TComplex::Conjugate(QvectorP[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::QGap0M(int n, int p)
{

	if(n>=0) return Qvector0M[n][p];
  else return TComplex::Conjugate(Qvector0M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::QGap0P(int n, int p)
{

	if(n>=0) return Qvector0P[n][p];
  else return TComplex::Conjugate(Qvector0P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::QGap2M(int n, int p)
{

	if(n>=0) return Qvector2M[n][p];
  else return TComplex::Conjugate(Qvector2M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::QGap2P(int n, int p)
{

	if(n>=0) return Qvector2P[n][p];
  else return TComplex::Conjugate(Qvector2P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::QGap4M(int n, int p)
{

	if(n>=0) return Qvector4M[n][p];
  else return TComplex::Conjugate(Qvector4M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::QGap4P(int n, int p)
{

	if(n>=0) return Qvector4P[n][p];
  else return TComplex::Conjugate(Qvector4P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::QGap8M(int n, int p)
{

	if(n>=0) return Qvector8M[n][p];
  else return TComplex::Conjugate(Qvector8M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::QGap8P(int n, int p)
{

	if(n>=0) return Qvector8P[n][p];
  else return TComplex::Conjugate(Qvector8P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::QGap14M(int n, int p)
{

	if(n>=0) return Qvector14M[n][p];
  else return TComplex::Conjugate(Qvector14M[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::QGap14P(int n, int p)
{

	if(n>=0) return Qvector14P[n][p];
  else return TComplex::Conjugate(Qvector14P[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::QsubLeft(int n, int p)
{

	if(n>=0) return QvectorSubLeft[n][p];
  else return TComplex::Conjugate(QvectorSubLeft[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::QsubRight(int n, int p)
{

	if(n>=0) return QvectorSubRight[n][p];
  else return TComplex::Conjugate(QvectorSubRight[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::QsubMiddle(int n, int p)
{

	if(n>=0) return QvectorSubMiddle[n][p];
  else return TComplex::Conjugate(QvectorSubMiddle[-n][p]);

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::Two(int n1, int n2)
{

	TComplex formula = Q(n1,1)*Q(n2,1) - Q(n1+n2,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::TwoGap0(int n1, int n2)
{

	TComplex formula = QGap0M(n1,1)*QGap0P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::TwoGap2(int n1, int n2)
{

	TComplex formula = QGap2M(n1,1)*QGap2P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::TwoGap4(int n1, int n2)
{

	TComplex formula = QGap4M(n1,1)*QGap4P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::TwoGap8(int n1, int n2)
{

	TComplex formula = QGap8M(n1,1)*QGap8P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::TwoGap(int n1, int n2)
{

	TComplex formula = QGapM(n1,1)*QGapP(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::TwoGap14(int n1, int n2)
{

	TComplex formula = QGap14M(n1,1)*QGap14P(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::Two_3SubLM(int n1, int n2)
{

	TComplex formula = QsubLeft(n1,1)*QsubMiddle(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::Two_3SubRM(int n1, int n2)
{

	TComplex formula = QsubMiddle(n1,1)*QsubRight(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::Two_3SubLM_perm2(int n1, int n2)
{

	TComplex formula = QsubLeft(n1,1)*QsubMiddle(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::Two_3SubLR_perm2(int n1, int n2)
{

	TComplex formula = QsubLeft(n1,1)*QsubRight(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::Two_3SubRL_perm3(int n1, int n2)
{

	TComplex formula = QsubRight(n1,1)*QsubLeft(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::Two_3SubRM_perm3(int n1, int n2)
{

	TComplex formula = QsubRight(n1,1)*QsubMiddle(n2,1);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::Three(int n1, int n2, int n3)
{

	TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)-Q(n1+n2,2)*Q(n3,1)-Q(n2,1)*Q(n1+n3,2)
 		                 - Q(n1,1)*Q(n2+n3,2)+2.*Q(n1+n2+n3,3);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::ThreeGap0M(int n1, int n2, int n3)                                                                                                                                 
{

  TComplex formula = QGap0M(n1,1)*QGap0M(n2,1)*QGap0M(n3,1)-QGap0M(n1+n2,2)*QGap0M(n3,1)-QGap0M(n2,1)*QGap0M(n1+n3,2)
                     - QGap0M(n1,1)*QGap0M(n2+n3,2)+2.*QGap0M(n1+n2+n3,3);
  return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::ThreeGap0P(int n1, int n2, int n3)
{

  TComplex formula = QGap0P(n1,1)*QGap0P(n2,1)*QGap0P(n3,1)-QGap0P(n1+n2,2)*QGap0P(n3,1)-QGap0P(n2,1)*QGap0P(n1+n3,2)
                     - QGap0P(n1,1)*QGap0P(n2+n3,2)+2.*QGap0P(n1+n2+n3,3);
  return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::Four(int n1, int n2, int n3, int n4)
{

	TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
                 		- Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
                 		+ Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
                 		+ 2.*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
                 		+ 2.*Q(n2,1)*Q(n1+n3+n4,3)+2.*Q(n1,1)*Q(n2+n3+n4,3)-6.*Q(n1+n2+n3+n4,4);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::FourGap0(int n1, int n2, int n3, int n4)
{
	
	TComplex formula = QGap0P(n1,1)*QGap0P(n2,1)*QGap0M(n3,1)*QGap0M(n4,1)-QGap0P(n1+n2,2)*QGap0M(n3,1)*QGap0M(n4,1)
                    -QGap0P(n1,1)*QGap0P(n2,1)*QGap0M(n3+n4,2)+QGap0P(n1+n2,2)*QGap0M(n3+n4,2);
	return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::FourGap0M(int n1, int n2, int n3, int n4)                                                                                                                          
{

  TComplex formula = QGap0M(n1,1)*QGap0M(n2,1)*QGap0M(n3,1)*QGap0M(n4,1)-QGap0M(n1+n2,2)*QGap0M(n3,1)*QGap0M(n4,1)-QGap0M(n2,1)*QGap0M(n1+n3,2)*QGap0M(n4,1)
                    - QGap0M(n1,1)*QGap0M(n2+n3,2)*QGap0M(n4,1)+2.*QGap0M(n1+n2+n3,3)*QGap0M(n4,1)-QGap0M(n2,1)*QGap0M(n3,1)*QGap0M(n1+n4,2)
                    + QGap0M(n2+n3,2)*QGap0M(n1+n4,2)-QGap0M(n1,1)*QGap0M(n3,1)*QGap0M(n2+n4,2)+QGap0M(n1+n3,2)*QGap0M(n2+n4,2)
                    + 2.*QGap0M(n3,1)*QGap0M(n1+n2+n4,3)-QGap0M(n1,1)*QGap0M(n2,1)*QGap0M(n3+n4,2)+QGap0M(n1+n2,2)*QGap0M(n3+n4,2)
                    + 2.*QGap0M(n2,1)*QGap0M(n1+n3+n4,3)+2.*QGap0M(n1,1)*QGap0M(n2+n3+n4,3)-6.*QGap0M(n1+n2+n3+n4,4);
  return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::FourGap0P(int n1, int n2, int n3, int n4)
{

  TComplex formula = QGap0P(n1,1)*QGap0P(n2,1)*QGap0P(n3,1)*QGap0P(n4,1)-QGap0P(n1+n2,2)*QGap0P(n3,1)*QGap0P(n4,1)-QGap0P(n2,1)*QGap0P(n1+n3,2)*QGap0P(n4,1)
                    - QGap0P(n1,1)*QGap0P(n2+n3,2)*QGap0P(n4,1)+2.*QGap0P(n1+n2+n3,3)*QGap0P(n4,1)-QGap0P(n2,1)*QGap0P(n3,1)*QGap0P(n1+n4,2)
                    + QGap0P(n2+n3,2)*QGap0P(n1+n4,2)-QGap0P(n1,1)*QGap0P(n3,1)*QGap0P(n2+n4,2)+QGap0P(n1+n3,2)*QGap0P(n2+n4,2)
                    + 2.*QGap0P(n3,1)*QGap0P(n1+n2+n4,3)-QGap0P(n1,1)*QGap0P(n2,1)*QGap0P(n3+n4,2)+QGap0P(n1+n2,2)*QGap0P(n3+n4,2)
                    + 2.*QGap0P(n2,1)*QGap0P(n1+n3+n4,3)+2.*QGap0P(n1,1)*QGap0P(n2+n3+n4,3)-6.*QGap0P(n1+n2+n3+n4,4);
  return formula;

}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::Four_3SubEvts(int n1, int n2, int n3, int n4) 
{
    TComplex formula = QsubMiddle(n1,1)*QsubMiddle(n2,1)*QsubLeft(n3,1)*QsubRight(n4,1)-QsubMiddle(n1+n2,2)*QsubLeft(n3,1)*QsubRight(n4,1);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::Four_3SubEvts_perm2(int n1, int n2, int n3, int n4) 
{
    TComplex formula = QsubLeft(n1,1)*QsubLeft(n2,1)*QsubMiddle(n3,1)*QsubRight(n4,1)-QsubLeft(n1+n2,2)*QsubMiddle(n3,1)*QsubRight(n4,1);
    return formula;
}
//____________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::Four_3SubEvts_perm3(int n1, int n2, int n3, int n4) 
{
    TComplex formula = QsubRight(n1,1)*QsubRight(n2,1)*QsubMiddle(n3,1)*QsubLeft(n4,1)-QsubRight(n1+n2,2)*QsubMiddle(n3,1)*QsubLeft(n4,1);
    return formula;
}
//___________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::Five(int n1, int n2, int n3, int n4, int n5)
{

	TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)
                 - Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)
                 + 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)
                 + Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)
                 + Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)
                 - Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)
                 + 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)
                 - 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)
                 + Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)
                 + Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)
                 - Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)
                 + Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)
                 - 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)
                 - 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)
                 + Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)
                 + Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)
                 + 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)
                 + 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)
                 - 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)
                 + Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)
                 + Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)
                 + 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)
                 + 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)
                 - 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)
                 - 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)
                 - 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)+24.*Q(n1+n2+n3+n4+n5,5); 
	return formula;

}
//___________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::Six(int n1, int n2, int n3, int n4, int n5, int n6)
{

	
 TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)
              - Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)
              + 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)
              + Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)
              + Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)*Q(n6,1)
              - Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)
              + 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)*Q(n6,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)*Q(n6,1)
              - 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)
              + Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)*Q(n6,1)
              + Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)*Q(n6,1)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)*Q(n6,1)
              - Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)
              + Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)*Q(n6,1)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)*Q(n6,1)
              - 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)*Q(n6,1)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)*Q(n6,1)
              - 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)
              + Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)*Q(n6,1)
              + Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)*Q(n6,1)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)*Q(n6,1)
              + 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)*Q(n6,1)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)*Q(n6,1)
              + 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)*Q(n6,1)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)*Q(n6,1)
              - 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)
              + Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)*Q(n6,1)
              + Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)*Q(n6,1)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)*Q(n6,1)
              + 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)*Q(n6,1)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)*Q(n6,1)
              + 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)*Q(n6,1)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)*Q(n6,1)
              - 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)*Q(n6,1)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)*Q(n6,1)
              - 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)*Q(n6,1)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)*Q(n6,1)
              - 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)*Q(n6,1)+24.*Q(n1+n2+n3+n4+n5,5)*Q(n6,1)
              - Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)+Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)
              + Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n1+n6,2)+Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n1+n6,2)
              - 2.*Q(n2+n3+n4,3)*Q(n5,1)*Q(n1+n6,2)+Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n1+n6,2)
              - Q(n3+n4,2)*Q(n2+n5,2)*Q(n1+n6,2)+Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n1+n6,2)
              - Q(n2+n4,2)*Q(n3+n5,2)*Q(n1+n6,2)-2.*Q(n4,1)*Q(n2+n3+n5,3)*Q(n1+n6,2)
              + Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n1+n6,2)-Q(n2+n3,2)*Q(n4+n5,2)*Q(n1+n6,2)
              - 2.*Q(n3,1)*Q(n2+n4+n5,3)*Q(n1+n6,2)-2.*Q(n2,1)*Q(n3+n4+n5,3)*Q(n1+n6,2)
              + 6.*Q(n2+n3+n4+n5,4)*Q(n1+n6,2)-Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)
              + Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)+Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n2+n6,2)
              + Q(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n2+n6,2)-2.*Q(n1+n3+n4,3)*Q(n5,1)*Q(n2+n6,2)
              + Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n2+n6,2)-Q(n3+n4,2)*Q(n1+n5,2)*Q(n2+n6,2)
              + Q(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n2+n6,2)-Q(n1+n4,2)*Q(n3+n5,2)*Q(n2+n6,2)
              - 2.*Q(n4,1)*Q(n1+n3+n5,3)*Q(n2+n6,2)+Q(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n2+n6,2)
              - Q(n1+n3,2)*Q(n4+n5,2)*Q(n2+n6,2)-2.*Q(n3,1)*Q(n1+n4+n5,3)*Q(n2+n6,2)
              - 2.*Q(n1,1)*Q(n3+n4+n5,3)*Q(n2+n6,2)+6.*Q(n1+n3+n4+n5,4)*Q(n2+n6,2)
              + 2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n2+n6,3)-2.*Q(n3+n4,2)*Q(n5,1)*Q(n1+n2+n6,3)
              - 2.*Q(n4,1)*Q(n3+n5,2)*Q(n1+n2+n6,3)-2.*Q(n3,1)*Q(n4+n5,2)*Q(n1+n2+n6,3)
              + 4.*Q(n3+n4+n5,3)*Q(n1+n2+n6,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)
              + Q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)+Q(n2,1)*Q(n1+n4,2)*Q(n5,1)*Q(n3+n6,2)
              + Q(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n3+n6,2)-2.*Q(n1+n2+n4,3)*Q(n5,1)*Q(n3+n6,2)
              + Q(n2,1)*Q(n4,1)*Q(n1+n5,2)*Q(n3+n6,2)-Q(n2+n4,2)*Q(n1+n5,2)*Q(n3+n6,2)
              + Q(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n3+n6,2)-Q(n1+n4,2)*Q(n2+n5,2)*Q(n3+n6,2)
              - 2.*Q(n4,1)*Q(n1+n2+n5,3)*Q(n3+n6,2)+Q(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n3+n6,2)
              - Q(n1+n2,2)*Q(n4+n5,2)*Q(n3+n6,2)-2.*Q(n2,1)*Q(n1+n4+n5,3)*Q(n3+n6,2)
              - 2.*Q(n1,1)*Q(n2+n4+n5,3)*Q(n3+n6,2)+6.*Q(n1+n2+n4+n5,4)*Q(n3+n6,2)
              + 2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n1+n3+n6,3)-2.*Q(n2+n4,2)*Q(n5,1)*Q(n1+n3+n6,3)
              - 2.*Q(n4,1)*Q(n2+n5,2)*Q(n1+n3+n6,3)-2.*Q(n2,1)*Q(n4+n5,2)*Q(n1+n3+n6,3)
              + 4.*Q(n2+n4+n5,3)*Q(n1+n3+n6,3)+2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n3+n6,3)
              - 2.*Q(n1+n4,2)*Q(n5,1)*Q(n2+n3+n6,3)-2.*Q(n4,1)*Q(n1+n5,2)*Q(n2+n3+n6,3)
              - 2.*Q(n1,1)*Q(n4+n5,2)*Q(n2+n3+n6,3)+4.*Q(n1+n4+n5,3)*Q(n2+n3+n6,3)
              - 6.*Q(n4,1)*Q(n5,1)*Q(n1+n2+n3+n6,4)+6.*Q(n4+n5,2)*Q(n1+n2+n3+n6,4)
              - Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)+Q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)
              + Q(n2,1)*Q(n1+n3,2)*Q(n5,1)*Q(n4+n6,2)+Q(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n4+n6,2)
              - 2.*Q(n1+n2+n3,3)*Q(n5,1)*Q(n4+n6,2)+Q(n2,1)*Q(n3,1)*Q(n1+n5,2)*Q(n4+n6,2)
              - Q(n2+n3,2)*Q(n1+n5,2)*Q(n4+n6,2)+Q(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n4+n6,2)
              - Q(n1+n3,2)*Q(n2+n5,2)*Q(n4+n6,2)-2.*Q(n3,1)*Q(n1+n2+n5,3)*Q(n4+n6,2)
              + Q(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n4+n6,2)-Q(n1+n2,2)*Q(n3+n5,2)*Q(n4+n6,2)
              - 2.*Q(n2,1)*Q(n1+n3+n5,3)*Q(n4+n6,2)-2.*Q(n1,1)*Q(n2+n3+n5,3)*Q(n4+n6,2)
              + 6.*Q(n1+n2+n3+n5,4)*Q(n4+n6,2)+2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n1+n4+n6,3)
              - 2.*Q(n2+n3,2)*Q(n5,1)*Q(n1+n4+n6,3)-2.*Q(n3,1)*Q(n2+n5,2)*Q(n1+n4+n6,3)
              - 2.*Q(n2,1)*Q(n3+n5,2)*Q(n1+n4+n6,3)+4.*Q(n2+n3+n5,3)*Q(n1+n4+n6,3)
              + 2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n4+n6,3)-2.*Q(n1+n3,2)*Q(n5,1)*Q(n2+n4+n6,3)
              - 2.*Q(n3,1)*Q(n1+n5,2)*Q(n2+n4+n6,3)-2.*Q(n1,1)*Q(n3+n5,2)*Q(n2+n4+n6,3)
              + 4.*Q(n1+n3+n5,3)*Q(n2+n4+n6,3)-6.*Q(n3,1)*Q(n5,1)*Q(n1+n2+n4+n6,4)
              + 6.*Q(n3+n5,2)*Q(n1+n2+n4+n6,4)+2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n4+n6,3)
              - 2.*Q(n1+n2,2)*Q(n5,1)*Q(n3+n4+n6,3)-2.*Q(n2,1)*Q(n1+n5,2)*Q(n3+n4+n6,3)
              - 2.*Q(n1,1)*Q(n2+n5,2)*Q(n3+n4+n6,3)+4.*Q(n1+n2+n5,3)*Q(n3+n4+n6,3)
              - 6.*Q(n2,1)*Q(n5,1)*Q(n1+n3+n4+n6,4)+6.*Q(n2+n5,2)*Q(n1+n3+n4+n6,4)
              - 6.*Q(n1,1)*Q(n5,1)*Q(n2+n3+n4+n6,4)+6.*Q(n1+n5,2)*Q(n2+n3+n4+n6,4)
              + 24.*Q(n5,1)*Q(n1+n2+n3+n4+n6,5)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)
              + Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5+n6,2)
              + Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5+n6,2)-2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5+n6,2)
              + Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5+n6,2)-Q(n2+n3,2)*Q(n1+n4,2)*Q(n5+n6,2)
              + Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5+n6,2)-Q(n1+n3,2)*Q(n2+n4,2)*Q(n5+n6,2)
              - 2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5+n6,2)+Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5+n6,2)
              - Q(n1+n2,2)*Q(n3+n4,2)*Q(n5+n6,2)-2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5+n6,2)
              - 2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5+n6,2)+6.*Q(n1+n2+n3+n4,4)*Q(n5+n6,2)
              + 2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5+n6,3)-2.*Q(n2+n3,2)*Q(n4,1)*Q(n1+n5+n6,3)
              - 2.*Q(n3,1)*Q(n2+n4,2)*Q(n1+n5+n6,3)-2.*Q(n2,1)*Q(n3+n4,2)*Q(n1+n5+n6,3)
              + 4.*Q(n2+n3+n4,3)*Q(n1+n5+n6,3)+2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5+n6,3)
              - 2.*Q(n1+n3,2)*Q(n4,1)*Q(n2+n5+n6,3)-2.*Q(n3,1)*Q(n1+n4,2)*Q(n2+n5+n6,3)
              - 2.*Q(n1,1)*Q(n3+n4,2)*Q(n2+n5+n6,3)+4.*Q(n1+n3+n4,3)*Q(n2+n5+n6,3)
              - 6.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5+n6,4)+6.*Q(n3+n4,2)*Q(n1+n2+n5+n6,4)
              + 2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5+n6,3)-2.*Q(n1+n2,2)*Q(n4,1)*Q(n3+n5+n6,3)
              - 2.*Q(n2,1)*Q(n1+n4,2)*Q(n3+n5+n6,3)-2.*Q(n1,1)*Q(n2+n4,2)*Q(n3+n5+n6,3)
              + 4.*Q(n1+n2+n4,3)*Q(n3+n5+n6,3)-6.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5+n6,4)
              + 6.*Q(n2+n4,2)*Q(n1+n3+n5+n6,4)-6.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5+n6,4)
              + 6.*Q(n1+n4,2)*Q(n2+n3+n5+n6,4)+24.*Q(n4,1)*Q(n1+n2+n3+n5+n6,5)
              + 2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5+n6,3)-2.*Q(n1+n2,2)*Q(n3,1)*Q(n4+n5+n6,3)
              - 2.*Q(n2,1)*Q(n1+n3,2)*Q(n4+n5+n6,3)-2.*Q(n1,1)*Q(n2+n3,2)*Q(n4+n5+n6,3)
              + 4.*Q(n1+n2+n3,3)*Q(n4+n5+n6,3)-6.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5+n6,4)
              + 6.*Q(n2+n3,2)*Q(n1+n4+n5+n6,4)-6.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5+n6,4)
              + 6.*Q(n1+n3,2)*Q(n2+n4+n5+n6,4)+24.*Q(n3,1)*Q(n1+n2+n4+n5+n6,5)
              - 6.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5+n6,4)+6.*Q(n1+n2,2)*Q(n3+n4+n5+n6,4)
              + 24.*Q(n2,1)*Q(n1+n3+n4+n5+n6,5)+24.*Q(n1,1)*Q(n2+n3+n4+n5+n6,5)
              - 120.*Q(n1+n2+n3+n4+n5+n6,6);
	return formula;

}
//_____________________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::SixGap0(int n1, int n2, int n3, int n4, int n5, int n6)
{                                                                                                                                                                                                       

  TComplex formula = ThreeGap0M(n1, n2, n3)*ThreeGap0P(n4, n5, n6); 
  return formula;

}
//_________________________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::Seven(int n1, int n2, int n3, int n4, int n5, int n6, int n7)
{

	TComplex Correlation = {0, 0};
	int Narray[] = {n1, n2, n3, n4, n5, n6};

	for(int k=7; k-->0; )
	{// backward loop of k from m-1 until 0, where m is the m-particle correlation, in this case m=4

		int array[6] = {0,1,2,3,4,5};
		int iPerm = 0;
		int argument = 0;
		int count = 0;

		// k==6: there is just one combination, we can add it manually
		if(k==6){
			Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
										Six(n1, n2, n3, n4, n5, n6)*Q(n7, 7-k);
		}// k==6

		else if(k==5){
			do{
				iPerm += 1;
				if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4]){
					count += 1;
					Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
												Five(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
														Narray[int(array[3])], Narray[int(array[4])])*
												Q(Narray[int(array[5])]+n7, 7-k);
				}
			}while(std::next_permutation(array, array+6));
		}// k==5

		else if(k==4){
			do{
				iPerm += 1;
				if(iPerm%2 == 1){
					if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3]){
						Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
													Four(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
															Narray[int(array[3])])*
													Q(Narray[int(array[4])]+Narray[int(array[5])]+n7, 7-k);
					}
				}
			}while(std::next_permutation(array, array+6));
		}// k==4
		
		else if(k==3){
			do{
				iPerm += 1;
				if(iPerm%6 == 1){
					if(array[0] < array[1] && array[1] < array[2]){
						Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
													Three(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])])*
													Q(Narray[int(array[3])]+Narray[int(array[4])]+Narray[int(array[5])]+n7, 7-k);
					}
				}
			}while(std::next_permutation(array, array+6));
		}// k==3
	
		else if(k==2){
			do{
				iPerm += 1;
				if(iPerm%24 == 1){
					if(array[0] < array[1]){
						Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
													Two(Narray[int(array[0])], Narray[int(array[1])])*
													Q(Narray[int(array[2])]+Narray[int(array[3])]+Narray[int(array[4])]
														+Narray[int(array[5])]+n7, 7-k);
					}
				}
			}while(std::next_permutation(array, array+6));
		}// k==2

		else if(k == 1){
			Correlation = Correlation 
									+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n1, 1)*Q(n2+n3+n4+n5+n6+n7, 7-k)
									+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n2, 1)*Q(n1+n3+n4+n5+n6+n7, 7-k)
									+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n3, 1)*Q(n1+n2+n4+n5+n6+n7, 7-k)
									+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n4, 1)*Q(n1+n2+n3+n5+n6+n7, 7-k)
									+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n5, 1)*Q(n1+n2+n3+n4+n6+n7, 7-k)
									+ TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n6, 1)*Q(n1+n2+n3+n4+n5+n7, 7-k);
		}// k==1

		else if(k == 0){
				Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n1+n2+n3+n4+n5+n6+n7, 7-k);
		}// k==0

		else{
			cout<<"invalid range of k"<<endl;
			return {0,0};
		}

	}// loop over k

	return Correlation;

} 
//_____________________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::Eight(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8)
{

	TComplex Correlation = {0, 0};
	int Narray[] = {n1, n2, n3, n4, n5, n6, n7};

	for(int k=8; k-->0; )
	{// backward loop of k from m-1 until 0, where m is the m-particle correlation, in this case m=4

		int array[7] = {0,1,2,3,4,5,6};
		int iPerm = 0;
		int argument = 0;
		int count = 0;

		// k==7: there is just one combination, we can add it manually
		if(k==7){
			Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
										Seven(n1, n2, n3, n4, n5, n6, n7)*Q(n8, 8-k);
		}// k==7

		else if(k==6){
			do{
				iPerm += 1;
				if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4] && array[4] < array[5]){
					count += 1;
					Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
												Six(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
														Narray[int(array[3])], Narray[int(array[4])], Narray[int(array[5])])*
												Q(Narray[int(array[6])]+n8, 8-k);
				}
			}while(std::next_permutation(array, array+7));
		}// k==6

		else if(k==5){
			do{
				iPerm += 1;
				if(iPerm%2 == 1){
					if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4]){
						Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
													Five(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])],
															Narray[int(array[3])], Narray[int(array[4])])*
													Q(Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
					}
				}
			}while(std::next_permutation(array, array+7));
		}// k==5
		
		else if(k==4){
			do{
				iPerm += 1;
				if(iPerm%6 == 1){
					if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3]){
						Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
													Four(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])], Narray[int(array[3])])*
													Q(Narray[int(array[4])]+Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
					}
				}
			}while(std::next_permutation(array, array+7));
		}// k==4
		
		else if(k==3){
			do{
				iPerm += 1;
				if(iPerm%24 == 1){
					if(array[0] < array[1] && array[1] < array[2]){
						Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
													Three(Narray[int(array[0])], Narray[int(array[1])], Narray[int(array[2])])*
													Q(Narray[int(array[3])]+Narray[int(array[4])]+Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
					}
				}
			}while(std::next_permutation(array, array+7));
		}// k==3
	
		else if(k==2){
			do{
				iPerm += 1;
				if(iPerm%120 == 1){
					if(array[0] < array[1]){
						Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
													Two(Narray[int(array[0])], Narray[int(array[1])])*
													Q(Narray[int(array[2])]+Narray[int(array[3])]+Narray[int(array[4])]
														+Narray[int(array[5])]+Narray[int(array[6])]+n8, 8-k);
					}
				}
			}while(std::next_permutation(array, array+7));
		}// k==2

		else if(k == 1){
			Correlation = Correlation 
									+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n1, 1)*Q(n2+n3+n4+n5+n6+n7+n8, 8-k)
									+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n2, 1)*Q(n1+n3+n4+n5+n6+n7+n8, 8-k)
									+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n3, 1)*Q(n1+n2+n4+n5+n6+n7+n8, 8-k)
									+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n4, 1)*Q(n1+n2+n3+n5+n6+n7+n8, 8-k)
									+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n5, 1)*Q(n1+n2+n3+n4+n6+n7+n8, 8-k)
									+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n6, 1)*Q(n1+n2+n3+n4+n5+n7+n8, 8-k)
									+ TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n7, 1)*Q(n1+n2+n3+n4+n5+n6+n8, 8-k);
		}// k==1

		else if(k == 0){
				Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n1+n2+n3+n4+n5+n6+n7+n8, 8-k);
		}// k==0

		else{
			cout<<"invalid range of k"<<endl;
			return {0,0};
		}

	}// loop over k

	return Correlation;

}
//_____________________________________________________________________________
TComplex AliAnalysisTaskChargedFlow::EightGap0(int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8)                                                                                  
{

  TComplex formula = FourGap0M(n1, n2, n3, n4)*FourGap0P(n5, n6, n7, n8);
  return formula;

}
//_____________________________________________________________________________
void AliAnalysisTaskChargedFlow::Terminate(const Option_t *)
{ 
  // Terminate loop
  Printf("Terminate()");
}
