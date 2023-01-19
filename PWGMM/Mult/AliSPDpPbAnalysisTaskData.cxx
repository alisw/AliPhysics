#ifndef AliSPDpPbAnalysisTaskData_cxx
#define AliSPDpPbAnalysisTaskData_cxx

#include "AliSPDpPbAnalysisTaskData.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TArrayF.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TParticle.h"
#include "TObjArray.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliHeader.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliESDVZERO.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliPPVsMultUtils.h"
#include "AliAnalysisUtils.h"
#include "AliVEvent.h"


ClassImp(AliSPDpPbAnalysisTaskData)      // classimp: necessary for root

//..................................................//

AliSPDpPbAnalysisTaskData::AliSPDpPbAnalysisTaskData() : AliAnalysisTaskSE(), 
    fESD(0),
    fEsdV0(0),
    fOutputList(0),
    fCentEstimator("V0A"),
    fTrigSel("kINT7"), 
    fHistVtxZ(0),
    fHistTrksVtxZ(0),
    fHistVtxXY(0),
    fHistVtxZEta(0),
    fHistEtaVtxZ(0),
    fHistPhiEta(0),
    fHistEtaPhi(0),
    fHistEta(0),
    fHistMult(0),
    fHistTrks05(0),
    fHistTrks10(0),
    fHistTrks15(0),
    fNTracklets(0),
    fHistTotEvent(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}

AliSPDpPbAnalysisTaskData::AliSPDpPbAnalysisTaskData(const char *name) : AliAnalysisTaskSE(name),
    fESD(0),
    fEsdV0(0),
    fOutputList(0),
    fCentEstimator("V0A"),
    fTrigSel("kINT7"), 
    fHistVtxZ(0),
    fHistTrksVtxZ(0),
    fHistVtxXY(0),
    fHistVtxZEta(0),
    fHistEtaVtxZ(0),
    fHistPhiEta(0),
    fHistEtaPhi(0),
    fHistEta(0),
    fHistMult(0),
    fHistTrks05(0),
    fHistTrks10(0),
    fHistTrks15(0),
    fNTracklets(0),
    fHistTotEvent(0)
{
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

AliSPDpPbAnalysisTaskData::~AliSPDpPbAnalysisTaskData()
{
    //destructor
    if(fOutputList) {
        delete fOutputList;     // list is deleted fromm memory at the end of the task
    }
}

void AliSPDpPbAnalysisTaskData::UserCreateOutputObjects()
{
    // this function is called once at the start of the analysis
    // histograms are added to the list

    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);

    fHistTotEvent = new TH1F("fHistTotEvent","fHistTotEvent",5,0,5);
    fHistTotEvent->SetTitle("");
    fHistTotEvent->GetXaxis()->SetBinLabel(1, "Processed events");
    fHistTotEvent->GetXaxis()->SetBinLabel(2, (fTrigSel+" trigger").Data());
    fHistTotEvent->GetXaxis()->SetBinLabel(3, "Pileup");                                                                      
    fHistTotEvent->GetXaxis()->SetBinLabel(4, "NContributors > 2");
    fHistTotEvent->GetXaxis()->SetBinLabel(5, "|VtxZ|<10:Selected for Analysis");

    fHistVtxZ = new TH1F("fHistVtxZ","fHistVtxZ", 101, -15, 15);
    fHistVtxZ->SetXTitle("Vertex - Z");
    fHistVtxZ->SetYTitle("# of events");
    fHistVtxZ->SetTitle("");
    fHistVtxZ->GetXaxis()->SetLabelSize(0.035);
    fHistVtxZ->GetYaxis()->SetLabelSize(0.035);
    fHistVtxZ->GetXaxis()->SetTitleSize(0.050);
    fHistVtxZ->GetYaxis()->SetTitleSize(0.050);
    //fHistVtxZ->SetStats(0);

    fNTracklets = new TH1F("fNTracklets", "fNTracklets", 200, 0, 200);
    fNTracklets->SetTitle("");
    fNTracklets->SetXTitle("Number of Tracklets");
    fNTracklets->SetYTitle("# of events");
    fNTracklets->GetXaxis()->SetLabelSize(0.035);
    fNTracklets->GetYaxis()->SetLabelSize(0.035);
    fNTracklets->GetXaxis()->SetTitleSize(0.050);
    fNTracklets->GetYaxis()->SetTitleSize(0.050);
    //fNTracklets->SetStats(0);

    fHistTrksVtxZ = new TH2F("fHistTrksVtxZ", "fHistTrksVtxZ", 100, -15, 15, 50, 0, 70);
    fHistTrksVtxZ->SetTitle("");
    fHistTrksVtxZ->SetXTitle("Vertex - Z");
    fHistTrksVtxZ->SetYTitle("Number of Tracklets");
    fHistTrksVtxZ->GetXaxis()->SetLabelSize(0.035);
    fHistTrksVtxZ->GetYaxis()->SetLabelSize(0.035);
    fHistTrksVtxZ->GetXaxis()->SetTitleSize(0.050);
    fHistTrksVtxZ->GetYaxis()->SetTitleSize(0.050);
    fHistTrksVtxZ->SetStats(0);

    fHistVtxXY = new TH2F("fHistVtxXY", "fHistVtxXY", 100, -0.5, 0.5, 100, -0.5, 0.5);
    fHistVtxXY->SetTitle("");
    fHistVtxXY->SetXTitle("V_{X}");
    fHistVtxXY->SetYTitle("V_{Y}");
    fHistVtxXY->GetXaxis()->SetLabelSize(0.035);
    fHistVtxXY->GetYaxis()->SetLabelSize(0.035);
    fHistVtxXY->GetXaxis()->SetTitleSize(0.050);
    fHistVtxXY->GetYaxis()->SetTitleSize(0.050);
    fHistVtxXY->SetStats(0);

    fHistVtxZEta = new TH2F("fHistVtxZEta","fHistVtxZEta", 100, -10, 10, 100, -2, 2);
    fHistVtxZEta->SetTitle("");
    fHistVtxZEta->SetXTitle("Vertex - Z");
    fHistVtxZEta->SetYTitle("#eta");
    fHistVtxZEta->GetXaxis()->SetLabelSize(0.035);
    fHistVtxZEta->GetYaxis()->SetLabelSize(0.035);
    fHistVtxZEta->GetXaxis()->SetTitleSize(0.050);
    fHistVtxZEta->GetYaxis()->SetTitleSize(0.050);
    fHistVtxZEta->SetStats(0);

    fHistEtaVtxZ = new TH2F("fHistEtaVtxZ","fHistEtaVtxZ", 100, -2, 2, 100, -10, 10);
    fHistEtaVtxZ->SetTitle("");
    fHistEtaVtxZ->SetXTitle("#eta");
    fHistEtaVtxZ->SetYTitle("Vertex - Z");
    fHistEtaVtxZ->GetXaxis()->SetLabelSize(0.035);
    fHistEtaVtxZ->GetYaxis()->SetLabelSize(0.035);
    fHistEtaVtxZ->GetXaxis()->SetTitleSize(0.050);
    fHistEtaVtxZ->GetYaxis()->SetTitleSize(0.050);
    fHistEtaVtxZ->SetStats(0);


    fHistPhiEta = new TH2F("fHistPhiEta", "fHistPhiEta", 100, 0, 6, 100, -2, 2);
    fHistPhiEta->SetTitle("");
    fHistPhiEta->SetXTitle("#phi");
    fHistPhiEta->SetYTitle("#eta");
    fHistPhiEta->GetXaxis()->SetLabelSize(0.035);
    fHistPhiEta->GetYaxis()->SetLabelSize(0.035);
    fHistPhiEta->GetXaxis()->SetTitleSize(0.050);
    fHistPhiEta->GetYaxis()->SetTitleSize(0.050);
    fHistPhiEta->SetStats(0);

    fHistEta = new TH1F("fHistEta", "fHistEta", 100, -2, 2);
    fHistEta->SetTitle("");
    fHistEta->SetXTitle("#eta");
    fHistEta->SetYTitle("# of events");
    fHistEta->GetXaxis()->SetLabelSize(0.035);
    fHistEta->GetYaxis()->SetLabelSize(0.035);
    fHistEta->GetXaxis()->SetTitleSize(0.050);
    fHistEta->GetYaxis()->SetTitleSize(0.050);
    //fHistEta->SetStats(0);

    fHistMult = new TH1F("fHistMult", "fHistMult", 200, 0, 200);
    fHistMult->SetTitle("");

    fHistTrks05 = new TH1F("fHistTrks05", "fHistTrks05", 200, 0, 200);
    fHistTrks05->SetTitle("");
    fHistTrks05->SetXTitle("Number of Tracklets |#eta|<0.5");
    fHistTrks05->SetYTitle("# of events");
    fHistTrks05->GetXaxis()->SetLabelSize(0.035);
    fHistTrks05->GetYaxis()->SetLabelSize(0.035);
    fHistTrks05->GetXaxis()->SetTitleSize(0.050);
    fHistTrks05->GetYaxis()->SetTitleSize(0.050);
    //fHistTrks05->SetStats(0);

    fHistTrks10 = new TH1F("fHistTrks10", "fHistTrks10", 200, 0, 200);
    fHistTrks10->SetTitle("");
    fHistTrks10->SetXTitle("Number of Tracklets |#eta|<1.0");
    fHistTrks10->SetYTitle("# of events");
    fHistTrks10->GetXaxis()->SetLabelSize(0.035);
    fHistTrks10->GetYaxis()->SetLabelSize(0.035);
    fHistTrks10->GetXaxis()->SetTitleSize(0.050);
    fHistTrks10->GetYaxis()->SetTitleSize(0.050);
    //fHistTrks10->SetStats(0);

    fHistTrks15 = new TH1F("fHistTrks15", "fHistTrks15", 200, 0, 200);
    fHistTrks15->SetTitle("");
    fHistTrks15->SetXTitle("Number of Tracklets |#eta|<1.5");
    fHistTrks15->SetYTitle("# of events");
    fHistTrks15->GetXaxis()->SetLabelSize(0.035);
    fHistTrks15->GetYaxis()->SetLabelSize(0.035);
    fHistTrks15->GetXaxis()->SetTitleSize(0.050);
    fHistTrks15->GetYaxis()->SetTitleSize(0.050);
    //fHistTrks15->SetStats(0);


    fOutputList->Add(fHistTotEvent);
    fOutputList->Add(fHistVtxZ);
    fOutputList->Add(fNTracklets);
    fOutputList->Add(fHistTrksVtxZ);
    fOutputList->Add(fHistVtxXY);
    fOutputList->Add(fHistVtxZEta);
    fOutputList->Add(fHistEtaVtxZ);
    fOutputList->Add(fHistPhiEta);
    fOutputList->Add(fHistEta);
    fOutputList->Add(fHistMult);
    fOutputList->Add(fHistTrks05);
    fOutputList->Add(fHistTrks10);
    fOutputList->Add(fHistTrks15);

    PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the 

}

void AliSPDpPbAnalysisTaskData::UserExec(Option_t *)
{
    // main loop called for each event

    // Pointer to a event----------------------------------------------------
    fESD = dynamic_cast<AliESDEvent*> (InputEvent());
    if (!fESD) {
        printf("ERROR: fESD not available\n");
        return;
    }

    fHistTotEvent->Fill(0.5);

    // Event selection-------------------------------------------------------
    Bool_t isSelected = 0;
    if(fTrigSel == "kINT7")
        isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7);
    if (!isSelected) return;

    fHistTotEvent->Fill(1.5);

    // Pile-up rejection-----------------------------------------------------    
    AliAnalysisUtils util;                                                                                                                      
    util.SetMinPlpContribMV(5);
    util.SetMaxPlpChi2MV(5);
    util.SetMinWDistMV(15);
    util.SetCheckPlpFromDifferentBCMV(kFALSE);
    Bool_t IsPileUpMV = util.IsPileUpMV(fESD);                                                                                                
    if(IsPileUpMV)return;                                                                                                                       
    fHistTotEvent->Fill(2.5);

    // Get primary vertex---------------------------------------------------- 
    const AliESDVertex *vertex = fESD->GetPrimaryVertex();
    float VtxZ = vertex->GetZ();
    float VtxX = vertex->GetX();
    float VtxY = vertex->GetY();

    // Get primary vertex SPD------------------------------------------------
    const AliESDVertex *spdVtx     = fESD->GetPrimaryVertexSPD();
    if (!spdVtx) return;
    float spdVtxZ = spdVtx->GetZ();
    float spdVtxX = spdVtx->GetX();
    float spdVtxY = spdVtx->GetY();

    // Vertex selection: If I have SPD vertex, with >2 contributors, then get only vertex with 10cm from IP    
    if(!spdVtx) return;
    if(spdVtx->GetNContributors() <= 2) return;
    fHistTotEvent->Fill(3.5);
    if(TMath::Abs(spdVtxZ) > 10) return;
    
    fHistTotEvent->Fill(4.5);

    fHistVtxZ->Fill(spdVtxZ);
    fHistVtxXY->Fill(spdVtxX, spdVtxY);

    // Centrality selection---------------------------------------------------
    // doing nothing with this yet
    MultSelection = (AliMultSelection *) fESD->FindListObject("MultSelection");
	if(!MultSelection) {
        AliWarning("AliMultSelection object not found!");
	}
	else{
	  float nCentrality = MultSelection->GetMultiplicityPercentile(fCentEstimator);
	}

    // Multiplicities---------------------------------------------------------
    // doing nothing with this yet
	fEsdV0 = fESD->GetVZEROData();
	float fV0Amult = fEsdV0->GetMTotV0A(); //returns total multiplicity in V0A 
	float fV0Cmult = fEsdV0->GetMTotV0C(); //returns total multiplicity in V0C
	float fV0mult  = fV0Amult + fV0Cmult;   //returns total multiplicity in V0A+V0C
	int GlobalTracks  = fESD->GetNumberOfTracks();  
    
    // Get number os tracklets for each event---------------------------------
    fMultiplicity = fESD -> GetMultiplicity();
    Int_t ntrks = fMultiplicity->GetNumberOfTracklets();
    fNTracklets->Fill(ntrks); 
    // histogram of number of tracklets of the event Vs the vertex Z-position
    fHistTrksVtxZ->Fill(spdVtxZ, ntrks);

	// Reconstructed number of tracklets
    Long_t lTrksEta15 = 0;
    Long_t lTrksEta10 = 0;
    Long_t lTrksEta05 = 0;
    Long_t lTrksAll = 0;

	for (auto it = 0; it<ntrks; it++) {
		Double_t eta = fMultiplicity->GetEta(it);
		Double_t phi = fMultiplicity->GetPhi(it);
        fHistPhiEta->Fill(phi, eta);
        fHistEta->Fill(eta);
        fHistVtxZEta->Fill(spdVtxZ, eta);
        fHistEtaVtxZ->Fill(eta, spdVtxZ);

        if( TMath::Abs(eta) < 0.5 ) lTrksEta05++;
        if( TMath::Abs(eta) < 1.0 ) lTrksEta10++;
        if( TMath::Abs(eta) < 1.5 ) lTrksEta15++;
        lTrksAll++;
    }
    fHistMult->Fill(fV0Amult);
    if(lTrksEta05!=0) {
        fHistTrks05->Fill(lTrksEta05);        
    }
    if(lTrksEta10!=0) {
        fHistTrks10->Fill(lTrksEta10);        
    }
    if(lTrksEta15!=0) {
        fHistTrks15->Fill(lTrksEta15);        
    }

} 

void AliSPDpPbAnalysisTaskData::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query
  
  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}


#endif