#ifndef AliSPDpPbAnalysisTask_cxx
#define AliSPDpPbAnalysisTask_cxx

#include "AliSPDpPbAnalysisTask.h"
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
#include "AliGenEventHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliPWG0Helper.h"

using namespace std;

ClassImp(AliSPDpPbAnalysisTask)      // classimp: necessary for root

//..................................................//

AliSPDpPbAnalysisTask::AliSPDpPbAnalysisTask() : AliAnalysisTaskSE(), 
    fESD(0),
    fEsdV0(0),
    fMCEvent(0),
    fMCStack(0),
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
    fMCHistVtxZ(0),
    fMCHistXY(0),
    fMCHistPhiEta(0),
    fMCHistPhiEta10(0),
    fMCHistEtaVxtZ(0),
    fMCHistPt(0),
    fMCHistEta(0),
    fMCHistEta10(0),
    fMCHistPhi(0),
    fMCHistNch(0),
    fMCHistNch10(0),
    fHistTotEvent(0),
    fResponseMatrix(0),
    fResponseMatrix10(0),
    fResponseMatrix101(0),
    fMCHistNch05(0),
    fMCHistNch15(0),
    fResponseMatrix05(0),
    fResponseMatrix051(0),
    fResponseMatrix15(0),
    fResponseMatrix151(0),
    fEffVsNch(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}

AliSPDpPbAnalysisTask::AliSPDpPbAnalysisTask(const char *name) : AliAnalysisTaskSE(name),
    fESD(0),
    fEsdV0(0),
    fMCEvent(0),
    fMCStack(0),
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
    fMCHistVtxZ(0),
    fMCHistXY(0),
    fMCHistPhiEta(0),
    fMCHistPhiEta10(0),
    fMCHistEtaVxtZ(0),
    fMCHistPt(0),
    fMCHistEta(0),
    fMCHistEta10(0),
    fMCHistPhi(0),
    fMCHistNch(0),
    fMCHistNch10(0),
    fHistTotEvent(0),
    fResponseMatrix(0),
    fResponseMatrix10(0),
    fResponseMatrix101(0),
    fMCHistNch05(0),
    fMCHistNch15(0),
    fResponseMatrix05(0),
    fResponseMatrix051(0),
    fResponseMatrix15(0),
    fResponseMatrix151(0),
    fEffVsNch(0)
{
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

AliSPDpPbAnalysisTask::~AliSPDpPbAnalysisTask()
{
    //destructor
    if(fOutputList) {
        delete fOutputList;     // list is deleted fromm memory at the end of the task
    }
}

void AliSPDpPbAnalysisTask::UserCreateOutputObjects()
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

    
    fMCHistVtxZ = new TH1F("fMCHistVtxZ","fMCHistVtxZ", 101, -15, 15);
    fMCHistVtxZ->SetXTitle("MC Vertex - Z");
    fMCHistVtxZ->SetYTitle("# of events");
    fMCHistVtxZ->GetXaxis()->SetLabelSize(0.035);
    fMCHistVtxZ->GetYaxis()->SetLabelSize(0.035);
    fMCHistVtxZ->GetXaxis()->SetTitleSize(0.050);
    fMCHistVtxZ->GetYaxis()->SetTitleSize(0.050);

    fMCHistXY = new TH2F("fMCHistXY","fMCHistXY", 100, -0.5, 0.5, 100, -0.5, 0.5);
    fMCHistXY->SetXTitle("V_{X}");
    fMCHistXY->SetYTitle("V_{Y}");
    fMCHistXY->GetXaxis()->SetLabelSize(0.035);
    fMCHistXY->GetYaxis()->SetLabelSize(0.035);
    fMCHistXY->GetXaxis()->SetTitleSize(0.050);
    fMCHistXY->GetYaxis()->SetTitleSize(0.050);

    fMCHistPhiEta = new TH2F("fMCHistPhiEta", "fMCHistPhiEta", 100, 0, 6, 100, -2, 2);
    fMCHistPhiEta->SetXTitle("#phi");
    fMCHistPhiEta->SetYTitle("#eta");
    fMCHistPhiEta->GetXaxis()->SetLabelSize(0.035);
    fMCHistPhiEta->GetYaxis()->SetLabelSize(0.035);
    fMCHistPhiEta->GetXaxis()->SetTitleSize(0.050);
    fMCHistPhiEta->GetYaxis()->SetTitleSize(0.050);

    fMCHistPhiEta10 = new TH2F("fMCHistPhiEta10", "fMCHistPhiEta10", 100, 0, 6, 100, -2, 2);
    fMCHistPhiEta10->SetXTitle("#phi");
    fMCHistPhiEta10->SetYTitle("#eta");
    fMCHistPhiEta10->GetXaxis()->SetLabelSize(0.035);
    fMCHistPhiEta10->GetYaxis()->SetLabelSize(0.035);
    fMCHistPhiEta10->GetXaxis()->SetTitleSize(0.050);
    fMCHistPhiEta10->GetYaxis()->SetTitleSize(0.050);

    fMCHistEtaVxtZ = new TH2F("fMCHistEtaVxtZ", "fMCHistEtaVxtZ", 100, -2, 2, 100, -10, 10);
    fMCHistEtaVxtZ->SetXTitle("#eta");
    fMCHistEtaVxtZ->SetYTitle("Vertez-Z");
    fMCHistEtaVxtZ->GetXaxis()->SetLabelSize(0.035);
    fMCHistEtaVxtZ->GetYaxis()->SetLabelSize(0.035);
    fMCHistEtaVxtZ->GetXaxis()->SetTitleSize(0.050);
    fMCHistEtaVxtZ->GetYaxis()->SetTitleSize(0.050);


    fMCHistEta = new TH1F("fMCHistEta", "fMCHistEta", 100, -2, 2);
    fMCHistEta->SetXTitle("#eta");
    fMCHistEta->SetYTitle("# of events (MC)");
    fMCHistEta->GetXaxis()->SetLabelSize(0.035);
    fMCHistEta->GetYaxis()->SetLabelSize(0.035);
    fMCHistEta->GetXaxis()->SetTitleSize(0.050);
    fMCHistEta->GetYaxis()->SetTitleSize(0.050);

    fMCHistEta10 = new TH1F("fMCHistEta10","fMCHistEta10",100,-2,2);
    fMCHistEta10->SetXTitle("#eta (|#eta|<1.0)");
    fMCHistEta10->SetYTitle("# of events (MC)");
    fMCHistEta10->GetXaxis()->SetLabelSize(0.035);
    fMCHistEta10->GetYaxis()->SetLabelSize(0.035);
    fMCHistEta10->GetXaxis()->SetTitleSize(0.050);
    fMCHistEta10->GetYaxis()->SetTitleSize(0.050);

    fMCHistPhi = new TH1F("fMCHistPhi", "fMCHistPhi", 100, 0, 6);
    fMCHistPhi->SetXTitle("#phi");
    fMCHistPhi->SetYTitle("# of events (MC)");
    fMCHistPhi->GetXaxis()->SetLabelSize(0.035);
    fMCHistPhi->GetYaxis()->SetLabelSize(0.035);
    fMCHistPhi->GetXaxis()->SetTitleSize(0.050);
    fMCHistPhi->GetYaxis()->SetTitleSize(0.050);

    fMCHistPt = new TH1F("fMCHistPt", "fMCHistPt", 100, 0, 7);
    fMCHistPt->SetXTitle("p_{t}");
    fMCHistPt->SetYTitle("# of events (MC)");
    fMCHistPt->GetXaxis()->SetLabelSize(0.035);
    fMCHistPt->GetYaxis()->SetLabelSize(0.035);
    fMCHistPt->GetXaxis()->SetTitleSize(0.050);
    fMCHistPt->GetYaxis()->SetTitleSize(0.050);

    fMCHistNch = new TH1F("fMCHistNch", "fMCHistNch", 200, 0, 200);
    fMCHistNch->SetXTitle("N_{ch}");
    fMCHistNch->SetYTitle("# of events (MC)");
    fMCHistNch->GetXaxis()->SetLabelSize(0.035);
    fMCHistNch->GetYaxis()->SetLabelSize(0.035);
    fMCHistNch->GetXaxis()->SetTitleSize(0.050);
    fMCHistNch->GetYaxis()->SetTitleSize(0.050);

    fMCHistNch10 = new TH1F("fMCHistNch10", "fMCHistNch10", 200, 0, 200);
    fMCHistNch10->SetXTitle("N_{ch}");
    fMCHistNch10->SetYTitle("# of events (MC)");
    fMCHistNch10->GetXaxis()->SetLabelSize(0.035);
    fMCHistNch10->GetYaxis()->SetLabelSize(0.035);
    fMCHistNch10->GetXaxis()->SetTitleSize(0.050);
    fMCHistNch10->GetYaxis()->SetTitleSize(0.050);

    fMCHistNch05 = new TH1F("fMCHistNch05", "fMCHistNch05", 200, 0, 200);
    fMCHistNch05->SetXTitle("N_{ch}");
    fMCHistNch05->SetYTitle("# of events (MC)");
    fMCHistNch05->GetXaxis()->SetLabelSize(0.035);
    fMCHistNch05->GetYaxis()->SetLabelSize(0.035);
    fMCHistNch05->GetXaxis()->SetTitleSize(0.050);
    fMCHistNch05->GetYaxis()->SetTitleSize(0.050);

    fMCHistNch15 = new TH1F("fMCHistNch15", "fMCHistNch15", 200, 0, 200);
    fMCHistNch15->SetXTitle("N_{ch}");
    fMCHistNch15->SetYTitle("# of events (MC)");
    fMCHistNch15->GetXaxis()->SetLabelSize(0.035);
    fMCHistNch15->GetYaxis()->SetLabelSize(0.035);
    fMCHistNch15->GetXaxis()->SetTitleSize(0.050);
    fMCHistNch15->GetYaxis()->SetTitleSize(0.050);

    fResponseMatrix = new TH2F("fResponseMatrix", "fResponseMatrix", 200, 0, 200, 200, 0, 200);
    fResponseMatrix->SetXTitle("MC_{truth}");
    fResponseMatrix->SetYTitle("MC_{reconstructed}");
    fResponseMatrix->GetXaxis()->SetLabelSize(0.035);
    fResponseMatrix->GetYaxis()->SetLabelSize(0.035);
    fResponseMatrix->GetXaxis()->SetTitleSize(0.050);
    fResponseMatrix->GetYaxis()->SetTitleSize(0.050);

    fResponseMatrix10 = new TH2F("fResponseMatrix10", "fResponseMatrix10", 200, 0, 200, 200, 0, 200);
    fResponseMatrix10->SetXTitle("MC_{truth}");
    fResponseMatrix10->SetYTitle("MC_{reconstructed}");
    fResponseMatrix10->GetXaxis()->SetLabelSize(0.035);
    fResponseMatrix10->GetYaxis()->SetLabelSize(0.035);
    fResponseMatrix10->GetXaxis()->SetTitleSize(0.050);
    fResponseMatrix10->GetYaxis()->SetTitleSize(0.050);

    fResponseMatrix101 = new TH2F("fResponseMatrix101", "fResponseMatrix101", 200, 0, 200, 200, 0, 200);
    fResponseMatrix101->SetYTitle("MC_{truth}");
    fResponseMatrix101->SetXTitle("MC_{reconstructed}");
    fResponseMatrix101->GetXaxis()->SetLabelSize(0.035);
    fResponseMatrix101->GetYaxis()->SetLabelSize(0.035);
    fResponseMatrix101->GetXaxis()->SetTitleSize(0.050);
    fResponseMatrix101->GetYaxis()->SetTitleSize(0.050);

    fResponseMatrix05 = new TH2F("fResponseMatrix05", "fResponseMatrix05", 200, 0, 200, 200, 0, 200);
    fResponseMatrix05->SetXTitle("MC_{truth}");
    fResponseMatrix05->SetYTitle("MC_{reconstructed}");
    fResponseMatrix05->GetXaxis()->SetLabelSize(0.035);
    fResponseMatrix05->GetYaxis()->SetLabelSize(0.035);
    fResponseMatrix05->GetXaxis()->SetTitleSize(0.050);
    fResponseMatrix05->GetYaxis()->SetTitleSize(0.050);

    fResponseMatrix051 = new TH2F("fResponseMatrix051", "fResponseMatrix051", 200, 0, 200, 200, 0, 200);
    fResponseMatrix051->SetYTitle("MC_{truth}");
    fResponseMatrix051->SetXTitle("MC_{reconstructed}");
    fResponseMatrix051->GetXaxis()->SetLabelSize(0.035);
    fResponseMatrix051->GetYaxis()->SetLabelSize(0.035);
    fResponseMatrix051->GetXaxis()->SetTitleSize(0.050);
    fResponseMatrix051->GetYaxis()->SetTitleSize(0.050);

    fResponseMatrix15 = new TH2F("fResponseMatrix15", "fResponseMatrix15", 200, 0, 200, 200, 0, 200);
    fResponseMatrix15->SetXTitle("MC_{truth}");
    fResponseMatrix15->SetYTitle("MC_{reconstructed}");
    fResponseMatrix15->GetXaxis()->SetLabelSize(0.035);
    fResponseMatrix15->GetYaxis()->SetLabelSize(0.035);
    fResponseMatrix15->GetXaxis()->SetTitleSize(0.050);
    fResponseMatrix15->GetYaxis()->SetTitleSize(0.050);

    fResponseMatrix151 = new TH2F("fResponseMatrix151", "fResponseMatrix151", 200, 0, 200, 200, 0, 200);
    fResponseMatrix151->SetYTitle("MC_{truth}");
    fResponseMatrix151->SetXTitle("MC_{reconstructed}");
    fResponseMatrix151->GetXaxis()->SetLabelSize(0.035);
    fResponseMatrix151->GetYaxis()->SetLabelSize(0.035);
    fResponseMatrix151->GetXaxis()->SetTitleSize(0.050);
    fResponseMatrix151->GetYaxis()->SetTitleSize(0.050);
    
    fEffVsNch = new TH2F("fEffVsNch", "fEffVsNch", 201, 0, 200, 101, 0 ,1);
    fEffVsNch->SetYTitle("Eff");
    fEffVsNch->SetXTitle("N_{ch}");
    fEffVsNch->GetXaxis()->SetLabelSize(0.035);
    fEffVsNch->GetYaxis()->SetLabelSize(0.035);
    fEffVsNch->GetXaxis()->SetTitleSize(0.050);
    fEffVsNch->GetYaxis()->SetTitleSize(0.050);

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

    fOutputList->Add(fMCHistVtxZ);
    fOutputList->Add(fMCHistXY);
    fOutputList->Add(fMCHistPhiEta);
    fOutputList->Add(fMCHistPhiEta10);
    fOutputList->Add(fMCHistEtaVxtZ);
    fOutputList->Add(fMCHistPt);
    fOutputList->Add(fMCHistEta);
    fOutputList->Add(fMCHistEta10);
    fOutputList->Add(fMCHistPhi);
    fOutputList->Add(fMCHistNch);
    fOutputList->Add(fMCHistNch10);
    fOutputList->Add(fResponseMatrix);
    fOutputList->Add(fResponseMatrix10);
    fOutputList->Add(fResponseMatrix101);
    fOutputList->Add(fMCHistNch05);
    fOutputList->Add(fMCHistNch15);
    fOutputList->Add(fResponseMatrix05);
    fOutputList->Add(fResponseMatrix051);
    fOutputList->Add(fResponseMatrix15);
    fOutputList->Add(fResponseMatrix151);
 
    fOutputList->Add(fEffVsNch);
 

    PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the 

}

void AliSPDpPbAnalysisTask::UserExec(Option_t *)
{
    // main loop called for each event

    // Pointer to a event----------------------------------------------------
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
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
    Long_t lTrksEta05 = 0;
    Long_t lTrksEta10 = 0;
    Long_t lTrksEta15 = 0;
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
    if(lTrksEta05!=0) {
        fHistTrks05->Fill(lTrksEta05);        
    }
    if(lTrksEta10!=0) {
        fHistTrks10->Fill(lTrksEta10);        
    }
    if(lTrksEta15!=0) {
        fHistTrks15->Fill(lTrksEta15);        
    }
    
    // Pointer to a MC event-------------------------------------------------	    
	AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
	
	if (!eventHandler) {
	  Printf("ERROR: Could not retrieve MC event handler");
	  return;
	}//eventHandler
	
	AliMCEvent *fMCEvent = eventHandler->MCEvent();
	if (!fMCEvent) {
	  Printf("ERROR: Could not retrieve MC event \n");
	  return;
	}//mcEvent
    
    // Access MC generated event----------------------------------------------
    TArrayF MC_Vtx_true_XYZ(3);
    fMCEvent->GenEventHeader()->PrimaryVertex(MC_Vtx_true_XYZ);
    fMCHistVtxZ->Fill(MC_Vtx_true_XYZ[2]);
    fMCHistXY->Fill(MC_Vtx_true_XYZ[0],MC_Vtx_true_XYZ[1]);

    // Generated number of charged particles
    Long_t lNchEta10 = 0;
    Long_t lNchEta05 = 0;
    Long_t lNchEta15 = 0;
    Long_t lNchAll = 0;

    // Loop on mc Tracks Starts Here 
    for (Int_t mcTrack = 0;  mcTrack < (fMCEvent->GetNumberOfTracks()); mcTrack++){   
        AliMCParticle *particle = (AliMCParticle*)fMCEvent->GetTrack(mcTrack);
	    if (!particle) {
	        Printf("ERROR: Could not receive track %d", mcTrack);
	        continue;
	    }//!track

        Double_t MCeta  = particle->Eta();
	    Double_t MCpt  = particle->Pt();
	    Double_t MCphi  = particle->Phi();

        Bool_t TrackIsPrim = particle->IsPhysicalPrimary();
        Bool_t TrackCharge = (particle->Charge())!=0;
        Bool_t TrackEtaMax10 = TMath::Abs(MCeta)<1.0;
        Bool_t TrackEtaMax05 = TMath::Abs(MCeta)<0.5;
        Bool_t TrackEtaMax15 = TMath::Abs(MCeta)<1.5;

        fMCHistPhiEta->Fill(MCphi, MCeta);
        fMCHistEtaVxtZ->Fill(MCeta, MC_Vtx_true_XYZ[2]);
        fMCHistPt->Fill(MCpt);
        fMCHistEta->Fill(MCeta);
        fMCHistPhi->Fill(MCphi);

        if (TrackIsPrim && TrackEtaMax05 && TrackCharge) {
            lNchEta05++;
        }
        if (TrackIsPrim && TrackEtaMax10 && TrackCharge) {
            lNchEta10++;
            fMCHistEta10->Fill(MCeta);
            fMCHistPhiEta10->Fill(MCphi, MCeta);
        }
        if (TrackIsPrim && TrackEtaMax15 && TrackCharge) {
            lNchEta15++;
        }
        if (TrackIsPrim && TrackCharge) {
            lNchAll++;
        }
	}//mcTrack loop

    // Fill true distributions
    fHistMult->Fill(fV0Amult);
    fMCHistNch->Fill(lNchAll);
    if (lNchEta05 != 0){
        fMCHistNch05->Fill(lNchEta05);
    } 
    if (lNchEta10 != 0){
        fMCHistNch10->Fill(lNchEta10);
    } 
    if (lNchEta15 != 0){
        fMCHistNch15->Fill(lNchEta15);
    }
    // Fill response matrices and efficiency hist
    if (lNchEta05 != 0 && lTrksEta05 != 0){
        fResponseMatrix05->Fill(lNchEta05, lTrksEta05);
        fResponseMatrix051->Fill(lTrksEta05, lNchEta05);
    } 
    if (lNchEta10 != 0 && lTrksEta10 != 0){
        fResponseMatrix10->Fill(lNchEta10, lTrksEta10);
        fResponseMatrix101->Fill(lTrksEta10, lNchEta10);
        fEffVsNch->Fill(lNchEta10, float(lTrksEta10)/float(lNchEta10) );
        // cout << "gen: " << lNchEta10 << " rec: " << lTrksEta10 << " eff: "<< float(lTrksEta10)/float(lNchEta10) << endl;
    } 
    if (lNchEta15 != 0 && lTrksEta15 != 0){
        fResponseMatrix15->Fill(lNchEta15, lTrksEta15);
        fResponseMatrix151->Fill(lTrksEta15, lNchEta15);
    }  



    

} 

void AliSPDpPbAnalysisTask::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query
  
  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}


#endif