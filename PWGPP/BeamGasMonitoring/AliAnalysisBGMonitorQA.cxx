//#include <Riostream.h>
#include <iostream>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH1D.h>
#include "TCanvas.h"
#include "TArrayI.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliVEvent.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisBGMonitorQA.h"
#include "AliLog.h"
#include "AliAnalysisFilter.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliTriggerAnalysis.h"
#include "AliAnalysisCuts.h"
#include "AliMultiplicity.h"
#include "AliESDVZERO.h"
#include "AliESDVZEROfriend.h"
#include "AliESDTZERO.h"
#include "AliAnalysisUtils.h"
#include "AliESDAD.h"
#include "AliESDADfriend.h"
using namespace std;
ClassImp(AliAnalysisBGMonitorQA)

Bool_t IsItBGSPDClusterVsTracklet(AliVEvent *event, float slope=4.); // add function info and initial condition (blim)
//________________________________________________________________________
AliAnalysisBGMonitorQA::AliAnalysisBGMonitorQA(const char *name) :
AliAnalysisTaskSE(name),
fESD(0x0),
fESDfriend(0x0),
fTreeTrack(0),
fList(0),
fUseTree(kFALSE)

{
    
    // Constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
   if(fUseTree==kTRUE) DefineOutput(2, TTree::Class());

    
}

//________________________________________________________________________
void AliAnalysisBGMonitorQA::ConnectInputData(Option_t *)
{
    
    TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
    if (!tree) {
        Printf("ERROR: Could not read chain from input slot 0");
    } else {
        
        AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
        
        if (esdH) {
            fESD = esdH->GetEvent();
            if(fESD) {
                fESDfriend = (AliESDfriend*)fESD->FindListObject("AliESDfriend");
                if (!fESDfriend){
                    AliError("No friend found");
                }
            }
        } else {
            Printf("ERROR: Could not get ESDInputHandler");
        }
        
    }
}

//________________________________________________________________________
void AliAnalysisBGMonitorQA::CreateOutputObjects()
{
    // Called once
    
    fTreeTrack = new TTree("TreeTrack","Track Properties");
    
    fTreeTrack->Branch("ntr",&ntr,"ntr/s"); // number of trigger classes
    fTreeTrack->Branch("nbunch",&nbunch,"nbunch/s"); // number of bunch
    fTreeTrack->Branch("runNumber",&runNumber,"runNumber/I"); //run number
    fTreeTrack->Branch("trigger",&ftrigger,"trigger[ntr]/I"); //trigger classes
    fTreeTrack->Branch("vertX",&fvertX,"vertX/D"); // primary x-vertex
    fTreeTrack->Branch("vertY",&fvertY,"vertY/D"); // primary y-vertex
    fTreeTrack->Branch("vertZ",&fvertZ,"vertZ/D"); // primary z-vertex
    fTreeTrack->Branch("vertTPCZ",&fvertTPCZ,"vertTPCZ/D"); // primary z-vertex
    fTreeTrack->Branch("vertTPCX",&fvertTPCX,"vertTPCX/D"); // primary z-vertex
    fTreeTrack->Branch("vertTPCY",&fvertTPCY,"vertTPCY/D"); // primary z-vertex
    fTreeTrack->Branch("v0a",&fv0a,"v0a/F"); // v0a mean time
    fTreeTrack->Branch("v0c",&fv0c,"v0c/F"); // v0c mean time
    fTreeTrack->Branch("fad0a",&fad0a,"fad0a/F"); // ad0a mean time
    fTreeTrack->Branch("fad0c",&fad0c,"fad0c/F"); // ad0c mean time
    fTreeTrack->Branch("multa",&fMulta,"multa/F"); // v0a multiplicity
    fTreeTrack->Branch("multc",&fMultc,"multc/F"); // v0c multiplicity
    fTreeTrack->Branch("triCha",&fTriCha,"triCha/F"); // v0a trigger charges
    fTreeTrack->Branch("triChc",&fTriChc,"triChc/F"); // v0c trigger charges
    fTreeTrack->Branch("bx",&fbx,"bx/I"); // BXid
    fTreeTrack->Branch("time",&ftime,"time/I"); // timestamp
    fTreeTrack->Branch("SpdC1",&fSpdC1,"SpdC1/I"); // SPD clusters layer 1
    fTreeTrack->Branch("SpdC2",&fSpdC2,"SpdC2/I"); // SPD clusters layer 1
    fTreeTrack->Branch("SpdT",&fSpdT,"SpdT/I"); // SPD tracklets
    fTreeTrack->Branch("tracks",&ntracks,"tracks/I"); // number of tracks
    fTreeTrack->Branch("v0abb",&V0A,"v0abb/I"); //
    fTreeTrack->Branch("v0cbb",&V0C,"v0cbb/I"); //
    fTreeTrack->Branch("v0abg",&V0ABG,"v0abg/I"); //
    fTreeTrack->Branch("v0cbg",&V0CBG,"v0cbg/I"); //
    fTreeTrack->Branch("n0abb",&nV0A,"n0abb/I"); //
    fTreeTrack->Branch("n0cbb",&nV0C,"n0cbb/I"); //
    fTreeTrack->Branch("n0abg",&nV0ABG,"n0abg/I"); //
    fTreeTrack->Branch("n0cbg",&nV0CBG,"n0cbg/I"); //
    fTreeTrack->Branch("vba",&VBA,"vba/I"); //
    fTreeTrack->Branch("vbc",&VBC,"vbc/I"); //
    fTreeTrack->Branch("vga",&VGA,"vga/I"); //
    fTreeTrack->Branch("vgc",&VGC,"vgc/I"); //
    fTreeTrack->Branch("vtx",&VTX,"vtx/I"); //
 //   fTreeTrack->Branch("fastOR",&fastORHW,"fastOR/I"); //
  //  fTreeTrack->Branch("chips1",&SPD1,"chips1/I"); //
  //  fTreeTrack->Branch("chips2",&SPD2,"chips2/I"); //
  //  fTreeTrack->Branch("chips3",&SPDHw1,"chips3/I"); //
  //  fTreeTrack->Branch("chips4",&SPDHw2,"chips4/I"); //
    fTreeTrack->Branch("bgID",&bgID,"bgID/I"); //
    fTreeTrack->Branch("bgID2",&bgID2,"bgID2/I"); //
  //  fTreeTrack->Branch("t0PU",&t0PileUp,"t0PU/I"); // pile-up from T0
    fTreeTrack->Branch("spdPU",&spdPileUp,"spdPU/I"); // pile-up from SPD
    fTreeTrack->Branch("spdPUOOB",&spdPileUpOutOfBunch,"spdPUOOB/I"); // out-of-bunch pile-up from SPD
    fTreeTrack->Branch("BGFlagA",&BGFlagA,"BGFlagA[nbunch]/I"); // V0A BG flag for PF protection
    fTreeTrack->Branch("BGFlagC",&BGFlagC,"BGFlagC[nbunch]/I"); // V0C BG flag for PF protection
    fTreeTrack->Branch("BBFlagA",&BBFlagA,"BBFlagA[nbunch]/I"); // V0A BG flag for PF protection
    fTreeTrack->Branch("BBFlagC",&BBFlagC,"BBFlagC[nbunch]/I"); // V0C BG flag for PF protection
    if(fUseTree==kTRUE)PostData(2, fTreeTrack);
    
    if(fList != NULL){
        delete fList;
        fList = NULL;
    }
    if(fList == NULL){
        fList = new TList();
        fList->SetOwner(kTRUE);
    }
    
    
   
    
    
    TH1F *hNumEffPurityBC[3][3][3];
    TH1F *hDenomEffBC[3][3][3];
    TH1F *hDenomPurityBC[3][3][3];
    TH1F *hDenomRejecEffBC[3][3][3];
    TH1F *hNumRejecEffBC[3][3][3];
    TH1F *hSPDNumBC[3][3][3];
    TH1F *hSPDDenomBC[3][3][3];
    TH2F *hNumTrkVsClsSPID[3][3][3];
    TH2F *hDenomTrkVsClsSPID[3][3][3];
    TH2F *hNumV0[3][3][3];
    TH2F *hDenomV0[3][3][3];
    
    TH1F *hADNumEffPurityBC[3][3][3];
    TH1F *hADDenomPurityBC[3][3][3];
    TH1F *hADNumRejecEffBC[3][3][3];
    TH2F *hADNumTrkVsClsSPD[3][3][3];
    TH2F *hNumAD[3][3][3];
    TH2F *hDenomAD[3][3][3];
    
    TH2F *hTotalTrkVsClsSPID = new TH2F("hTotalTrkVsClsSPID","; Spd : total",140,0,140,500,0,500);
    hTotalTrkVsClsSPID->GetXaxis()->SetTitle("Tracklet");
    hTotalTrkVsClsSPID->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList->Add(hTotalTrkVsClsSPID);
    
    TH2F *hTotalV0 = new TH2F("hTotalV0","; V0 : total",600,-300,300,2000,-1000,1000);
    hTotalV0->GetXaxis()->SetTitle("V0A-V0C");
    hTotalV0->GetYaxis()->SetTitle("V0A+V0C");
    fList->Add(hTotalV0);
    
    TH2F *hTotalAD = new TH2F("hTotalAD","; AD : total",400,-200,200,2000,-1000,1000);
    hTotalAD->GetXaxis()->SetTitle("ADA-ADC");
    hTotalAD->GetYaxis()->SetTitle("ADA+ADC");
    fList->Add(hTotalAD);
    
    //histogram for event list(blim)
    TH1F *hNumEvents  = new TH1F("hNumEvents","total event",10,0,10);
    fList->Add(hNumEvents);
    
    //histogram for modified Cut(blim)
    
    TH2F *hTrkVsClsSPIDSlope3 = new TH2F("hTrkVsClsSPIDSlope3","; Spd : total",140,0,140,500,0,500);
    hTrkVsClsSPIDSlope3->GetXaxis()->SetTitle("Tracklet");
    hTrkVsClsSPIDSlope3->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList->Add(hTrkVsClsSPIDSlope3);

    TH2F *hTrkVsClsSPIDSlope4 = new TH2F("hTrkVsClsSPIDSlope4","; Spd : total",140,0,140,500,0,500);
    hTrkVsClsSPIDSlope4->GetXaxis()->SetTitle("Tracklet");
    hTrkVsClsSPIDSlope4->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList->Add(hTrkVsClsSPIDSlope4);
    
    TH2F *hTrkVsClsSPIDSlope5 = new TH2F("hTrkVsClsSPIDSlope5","; Spd : total",140,0,140,500,0,500);
    hTrkVsClsSPIDSlope5->GetXaxis()->SetTitle("Tracklet");
    hTrkVsClsSPIDSlope5->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList->Add(hTrkVsClsSPIDSlope5);
    
    TH2F *hTrkVsClsSPIDSlope0 = new TH2F("hTrkVsClsSPIDSlope0","; Spd : total",140,0,140,500,0,500);
    hTrkVsClsSPIDSlope0->GetXaxis()->SetTitle("Tracklet");
    hTrkVsClsSPIDSlope0->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList->Add(hTrkVsClsSPIDSlope0);
    
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            for(int k=0; k<3; k++){

        hNumEffPurityBC[i][j][k] = new TH1F(Form("hNumEffPurityBC%d_V0%d_Flag%d",i,j,k),"; #V0flags in PF", 35, 0, 35);
        hDenomEffBC[i][j][k] = new TH1F(Form("hDenomEffBC%d_V0%d_Flag%d",i,j,k),"; #V0flags in PF", 35, 0, 35);
        hDenomPurityBC[i][j][k] = new TH1F(Form("hDenomPurityBC%d_V0%d_Flag%d",i,j,k),"; #V0flags in PF", 35, 0, 35);
        hDenomRejecEffBC[i][j][k] = new TH1F(Form("hDenomRejecEffBC%d_V0%d_Flag%d",i,j,k),"; #V0flags in PF", 35, 0, 35);
        hNumRejecEffBC[i][j][k] = new TH1F(Form("hNumRejecEffBC%d_V0%d_Flag%d",i,j,k),"; #V0flags in PF", 35, 0, 35);

        hSPDNumBC[i][j][k] = new TH1F(Form("hSPDNumBC%d_V0%d_Flag%d",i,j,k),"; Spd tracklet", 200, 0, 200);
        hSPDDenomBC[i][j][k] = new TH1F(Form("hSPDDenomBC%d_V0%d_Flag%d",i,j,k),"; Spd tracklet", 200, 0, 200);

        hNumTrkVsClsSPID[i][j][k] = new TH2F(Form("hNumTrkVsClsSPID%d_V0%d_Flag%d",i,j,k),"; Spd : !BGid & GoodEvent",140,0,140,500,0,500);
        hNumTrkVsClsSPID[i][j][k]->GetXaxis()->SetTitle("Tracklet");
        hNumTrkVsClsSPID[i][j][k]->GetYaxis()->SetTitle("Cluster (fspdC1)");
        hDenomTrkVsClsSPID[i][j][k]= new TH2F(Form("hDenomTrkVsClsSPID%d_V0%d_Flag%d",i,j,k),"; Spd : !BGid",140,0,140,500,0,500);
        hDenomTrkVsClsSPID[i][j][k]->GetXaxis()->SetTitle("Tracklet");
        hDenomTrkVsClsSPID[i][j][k]->GetYaxis()->SetTitle("Cluster (fspdC1)");
                
       hNumV0[i][j][k] = new TH2F(Form("hNumV0%d_V0%d_Flag%d",i,j,k),"; V0 : !BGid & GoodEvent",600,-300,300,2000,-1000,1000);
       hNumV0[i][j][k]->GetXaxis()->SetTitle("V0A-V0C");
       hNumV0[i][j][k]->GetYaxis()->SetTitle("V0A+V0C");
       hDenomV0[i][j][k]= new TH2F(Form("hDenomV0%d_V0%d_Flag%d",i,j,k),"; V0 : !BGid",600,-300,300,2000,-1000,1000);
       hDenomV0[i][j][k]->GetXaxis()->SetTitle("V0A-V0C");
       hDenomV0[i][j][k]->GetYaxis()->SetTitle("V0A+V0C");

                
                hADNumEffPurityBC[i][j][k] = new TH1F(Form("hADNumEffPurityBC%d_AD%d_Flag%d",i,j,k),"; #ADflags in PF", 10, 0, 10);
                hADDenomPurityBC[i][j][k] = new TH1F(Form("hADDenomPurityBC%d_AD%d_Flag%d",i,j,k),"; #ADflags in PF", 10, 0, 10);
                hADNumRejecEffBC[i][j][k] = new TH1F(Form("hADNumRejecEffBC%d_AD%d_Flag%d",i,j,k),"; #ADflags in PF", 10, 0, 10);
                hADNumTrkVsClsSPD[i][j][k] = new TH2F(Form("hADNumTrkVsClsSPD%d_AD%d_Flag%d",i,j,k),"; Spd : !BGid & GoodEvent",140,0,140,500,0,500);
                hADNumTrkVsClsSPD[i][j][k]->GetXaxis()->SetTitle("Tracklet");
                hADNumTrkVsClsSPD[i][j][k]->GetYaxis()->SetTitle("Cluster (fspdC1)");
                
                hNumAD[i][j][k] = new TH2F(Form("hNumAD%d_AD%d_Flag%d",i,j,k),"; AD : !BGid & GoodEvent",400,-200,200,2000,-1000,1000);
                hNumAD[i][j][k]->GetXaxis()->SetTitle("ADA-ADC");
                hNumAD[i][j][k]->GetYaxis()->SetTitle("ADA+ADC");
                hDenomAD[i][j][k]= new TH2F(Form("hDenomAD%d_AD%d_Flag%d",i,j,k),"; AD : !BGid",400,-200,200,2000,-1000,1000);
                hDenomAD[i][j][k]->GetXaxis()->SetTitle("ADA-ADC");
                hDenomAD[i][j][k]->GetYaxis()->SetTitle("ADA+ADC");
        
        
        fList->Add(hNumEffPurityBC[i][j][k]);
        fList->Add(hDenomEffBC[i][j][k]);
        fList->Add(hDenomPurityBC[i][j][k]);
        fList->Add(hDenomRejecEffBC[i][j][k]);
        fList->Add(hNumRejecEffBC[i][j][k]);
        fList->Add(hSPDNumBC[i][j][k]);
        fList->Add(hSPDDenomBC[i][j][k]);
        fList->Add(hNumTrkVsClsSPID[i][j][k]);
        fList->Add(hDenomTrkVsClsSPID[i][j][k]);
        fList->Add(hNumV0[i][j][k]);
        fList->Add(hDenomV0[i][j][k]);
       
        fList->Add(hADNumEffPurityBC[i][j][k]);
        fList->Add(hADDenomPurityBC[i][j][k]);
        fList->Add(hADNumRejecEffBC[i][j][k]);
        fList->Add(hADNumTrkVsClsSPD[i][j][k]);
        fList->Add(hNumAD[i][j][k]);
        fList->Add(hDenomAD[i][j][k]);

                
                
            }
        }
        
    }
    
    
    
    PostData(1, fList);
}

//________________________________________________________________________
void AliAnalysisBGMonitorQA::Exec(Option_t *)
{
    // Called for each event
    
    if (!fESD) {
        Printf("ERROR: fESD not available");
        return;
    }
    
    Int_t iEv= 0;
    iEv = fESD->GetEventNumberInFile();
    runNumber = fESD->GetRunNumber();
    UInt_t timeGDC=fESD->GetTimeStamp();
    ftime=timeGDC;
    Int_t timeStampBX = fESD->GetBunchCrossNumber();
    fbx=timeStampBX;
    ntr = 10;
    nbunch = 21;
    ofstream ftxt;
    
    static AliTriggerAnalysis * triggerAnalysis = new AliTriggerAnalysis();
    
    V0A = 0;
    V0C = 0;
    V0ABG = 0;
    V0CBG = 0;
    bgID = 0;
    
    // additional value initialize (blim)
    bgID2=0;
    bgID3=0;
    bgID4=0;
    bgID5=0;
    
    VBA = 0;
    VBC = 0;
    VGA = 0;
    VTX = 0;
   fastORHW = 0;
   SPD1 = 0;
   SPD2 = 0;
   SPDHw1 = 0;
   SPDHw2 = 0;
 
    fastORHW = triggerAnalysis->EvaluateTrigger(fESD, AliTriggerAnalysis::kSPDGFO); // SPD number of chips from trigger bits (!)
    SPD1 = triggerAnalysis->SPDFiredChips(fESD,0,kFALSE,1);  //SPD Fired Chips in layer 1 (from cluster)
    SPD2 = triggerAnalysis->SPDFiredChips(fESD,0,kFALSE,2);  //SPD Fired Chips in layer 2 (from cluster)
    SPDHw1 = triggerAnalysis->SPDFiredChips(fESD,1,kFALSE,1);  //SPD Fired Chips in layer 1 (from hardware bit)
    SPDHw2 = triggerAnalysis->SPDFiredChips(fESD,1,kFALSE,2);  //SPD Fired Chips in layer 2 (from hardware bit)
    t0PileUp = triggerAnalysis->EvaluateTrigger(fESD, (AliTriggerAnalysis::Trigger) (AliTriggerAnalysis::kOfflineFlag | AliTriggerAnalysis::kT0Pileup)); //T0 pile-up
    
    AliVVZERO *vzero = fESD->GetVZEROData();
    V0A   = (vzero->GetV0ADecision()==AliVVZERO::kV0BB);
    V0ABG = (vzero->GetV0ADecision()==AliVVZERO::kV0BG);
    V0C   = (vzero->GetV0CDecision()==AliVVZERO::kV0BB);
    V0CBG = (vzero->GetV0CDecision()==AliVVZERO::kV0BG);
    
    AliAnalysisUtils *utils = new AliAnalysisUtils();
    bgID = utils->IsSPDClusterVsTrackletBG(fESD);
    
    // modified slope cut. the function is in below of this source(blim)
    bgID2 = IsItBGSPDClusterVsTracklet(fESD,0.); // bgID2->slope0
    bgID3 = IsItBGSPDClusterVsTracklet(fESD,3.); // bgID3->slope3
    bgID4 = IsItBGSPDClusterVsTracklet(fESD,4.); // bgID4->slope4
    bgID5 = IsItBGSPDClusterVsTracklet(fESD,5.); // bgID5->slope5
    
    spdPileUp = utils->IsPileUpSPD(fESD);
    spdPileUpOutOfBunch = utils->IsOutOfBunchPileUp(fESD);
    
    
    //CTP inputs
    VTX = fESD->GetHeader()->IsTriggerInputFired("0TVX");
    VGA = fESD->GetHeader()->IsTriggerInputFired("0VGA");
    VGC = fESD->GetHeader()->IsTriggerInputFired("0VGC");
    VBA = fESD->GetHeader()->IsTriggerInputFired("0VBA");
    VBC = fESD->GetHeader()->IsTriggerInputFired("0VBC");
    triMask = fESD->GetHeader()->GetTriggerMask();
    
    //--- vertex
    const AliESDVertex *vertSPD=fESD->GetPrimaryVertexSPD();
    if(vertSPD->GetNContributors()>0){
        fvertZ=vertSPD->GetZ();
        fvertX=vertSPD->GetX();
        fvertY=vertSPD->GetY();
    }
    else{
        fvertZ=-99999;
        fvertX=-99999;
        fvertY=-99999;
    }
    
    const AliESDVertex *vertTPC=fESD->GetPrimaryVertexTracks();
    if(vertTPC->GetNContributors()>0){
        fvertTPCZ=vertTPC->GetZ();
        fvertTPCX=vertTPC->GetX();
        fvertTPCY=vertTPC->GetY();
    }
    else{
        fvertTPCZ=-99999;
        fvertTPCX=-99999;
        fvertTPCY=-99999;
    }
    
    //--- SPD cluster and tracklets
    const AliMultiplicity* mult = fESD->GetMultiplicity();
    
    fSpdC1 = 0;
    fSpdC2 = 0;
    //for(Int_t ilayer = 0; ilayer < 2; ilayer++){
    //  fSpdC += mult->GetNumberOfITSClusters(ilayer);
    //}
    fSpdC1 = mult->GetNumberOfITSClusters(0);
    fSpdC2 = mult->GetNumberOfITSClusters(1);
    
    fSpdT = mult->GetNumberOfTracklets();
    
    //--- V0 data
    //AliESDVZERO* vzero = fESD->GetVZEROData();
    fv0a = vzero->GetV0ATime();  //V0A time
    fv0c = vzero->GetV0CTime();  //V0C time
    fMulta = vzero->GetMTotV0A();  //V0A multiplicity
    fMultc = vzero->GetMTotV0C();  //V0C multiplicity
    fTriCha = vzero->GetTriggerChargeA();  //Sum of the trigger (clock=10) charge on A side (Ring 0 excluded)
    fTriChc = vzero->GetTriggerChargeC();  //Sum of the trigger (clock=10) charge on A side
    
    
    //-- AD data
    AliESDAD* adzero = fESD->GetADData();
    fad0a = adzero->GetADATime();
    fad0c = adzero->GetADCTime();
    
    
    
    // "online" V0 flags
    nV0A = 0;
    nV0ABG = 0;
    for (Int_t i = 32; i < 64; ++i) {
        if (vzero->GetBBFlag(i)) nV0A++;
        if (vzero->GetBGFlag(i)) nV0ABG++;
    }
    nV0C = 0;
    nV0CBG = 0;
    for (Int_t i = 0; i < 32; ++i) {
        if (vzero->GetBBFlag(i)) nV0C++;
        if (vzero->GetBGFlag(i)) nV0CBG++;
    }
    
    memset(BGFlagA, 0, sizeof(Float_t)*nbunch);
    memset(BBFlagA, 0, sizeof(Float_t)*nbunch);
    memset(BGFlagC, 0, sizeof(Float_t)*nbunch);
    memset(BBFlagC, 0, sizeof(Float_t)*nbunch);
    
    memset(ADBGFlagA, 0, sizeof(Float_t)*nbunch);
    memset(ADBBFlagA, 0, sizeof(Float_t)*nbunch);
    memset(ADBGFlagC, 0, sizeof(Float_t)*nbunch);
    memset(ADBBFlagC, 0, sizeof(Float_t)*nbunch);
    
    AliESDVZEROfriend *esdV0friend = fESDfriend->GetVZEROfriend();
    if(esdV0friend) {
        for(Int_t j = 0; j < 20; j++){
            //V0 --- infor
            for (Int_t i = 32; i < 64; ++i) {
                //BBFlagA[j] |= esdV0friend->GetBBFlag(i,j);
                //BGFlagA[j] |= esdV0friend->GetBGFlag(i,j);
                if(esdV0friend->GetBBFlag(i,j)) BBFlagA[j]++;
                if(esdV0friend->GetBGFlag(i,j)) BGFlagA[j]++;
            }
            for (Int_t i = 0; i < 32; ++i) {
                //BBFlagC[j] |= esdV0friend->GetBBFlag(i,j);
                //BGFlagC[j] |= esdV0friend->GetBGFlag(i,j);
                if(esdV0friend->GetBBFlag(i,j)) BBFlagC[j]++;
                if(esdV0friend->GetBGFlag(i,j)) BGFlagC[j]++;
            }
            //AD --- infor
            for (Int_t i = 8; i < 16; ++i) {
                if(esdV0friend->GetBBFlag(i,j)) ADBBFlagA[j]++;
                if(esdV0friend->GetBGFlag(i,j)) ADBGFlagA[j]++;
            }
            for (Int_t i = 0; i < 8; ++i) {
                if(esdV0friend->GetBBFlag(i,j)) ADBBFlagC[j]++;
                if(esdV0friend->GetBGFlag(i,j)) ADBGFlagC[j]++;
            }
            
            
            
        }
    } else {
        Printf("No esdV0friend available");
        return;
    }
    
      ntracks = fESD->GetNumberOfTracks(); // number of tracks (no quality cuts)
    
    //--- Trigger classes --//
    memset(ftrigger, 0, sizeof(Float_t)*ntr);
    
    //Minimum Bias
    if(fESD->IsTriggerClassFired("CINT7-B-NOPF-ALLNOTRD") || fESD->IsTriggerClassFired("CINT7-S-NOPF-ALLNOTRD") || fESD->IsTriggerClassFired("CINT1-B-NOPF-ALLNOTRD") || fESD->IsTriggerClassFired("CINT1-S-NOPF-ALLNOTRD")) ftrigger[0] = 1;
    if(fESD->IsTriggerClassFired("CINT7-AC-NOPF-ALLNOTRD") || fESD->IsTriggerClassFired("CINT7-ACE-NOPF-ALLNOTRD") || fESD->IsTriggerClassFired("CINT1-AC-NOPF-ALLNOTRD")) ftrigger[1] = 1;
    if(fESD->IsTriggerClassFired("C0VGA-B-NOPF-ALLNOTRD") || fESD->IsTriggerClassFired("C0VGA-AC-NOPF-ALLNOTRD") || fESD->IsTriggerClassFired("C0VGA-ABCE-NOPF-ALLNOTRD")) ftrigger[2] = 1;
    if(fESD->IsTriggerClassFired("C0VGC-B-NOPF-ALLNOTRD") || fESD->IsTriggerClassFired("C0VGC-AC-NOPF-ALLNOTRD") || fESD->IsTriggerClassFired("C0VGC-ABCE-NOPF-ALLNOTRD")) ftrigger[3] = 1;
    if(fESD->IsTriggerClassFired("CVGO-ABCE-NOPF-ALLNOTRD")) ftrigger[4] = 1;
    //Zero Bias
    if(fESD->IsTriggerClassFired("CBEAMB-B-NOPF-ALLNOTRD") || fESD->IsTriggerClassFired("CTRUE-S-NOPF-ALLNOTRD")) ftrigger[5] = 1;
    //T0 triggers
    if(fESD->IsTriggerClassFired("CINT8-S-NOPF-ALLNOTRD")) ftrigger[6] = 1;
    //Power trigger
    if(fESD->IsTriggerClassFired("CSPI7-B-NOPF-ALLNOTRD") || fESD->IsTriggerClassFired("CSPI7-S-NOPF-ALLNOTRD")) ftrigger[7] = 1;
    //High-multiplicity triggers
    if(fESD->IsTriggerClassFired("CSHM8-S-NOPF-ALLNOTRD")) ftrigger[8] = 1;
    //    if(fESD->IsTriggerClassFired("CSHM8-ACE-NOPF-ALLNOTRD")) ftrigger[9] = 1; // for LHC11h

    if(fESD->IsTriggerClassFired("CTEST52-B-NOPF-ALLNOTRD")) ftrigger[9] = 1; // for LHC15f, run 297884

    // count total event number (blim)
    if(ftrigger[0]==1)((TH1F*)fList->FindObject("hNumEvents"))->Fill(1);
    if(ftrigger[9]==1)((TH1F*)fList->FindObject("hNumEvents"))->Fill(3);

    
    int nAfterBunch = 3;
    int nV0 = 3;
    int nFlag = 3;
    
   // Bool_t goodEvent;
    Bool_t SelGoodEvent[nAfterBunch][nV0][nFlag];
    Bool_t SelGoodEventAD[nAfterBunch][nV0][nFlag];


    if(ftrigger[9]) {  // trigger class //0 for MB, 9 for HM
        
        ((TH1F*)fList->FindObject("hTotalTrkVsClsSPID"))->Fill(fSpdT, fSpdC1+fSpdC2);
        ((TH1F*)fList->FindObject("hTotalV0"))->Fill(fv0a-fv0c, fv0a+fv0c);
        ((TH1F*)fList->FindObject("hTotalAD"))->Fill(fad0a-fad0c, fad0a+fad0c);
        
        
        // Modified Cut result, added by blim
        if (!bgID2) {
            ((TH1F*)fList->FindObject("hTrkVsClsSPIDSlope0"))->Fill(fSpdT, fSpdC1+fSpdC2); // bgID2->slope0
        }
        if (!bgID3) {
            ((TH1F*)fList->FindObject("hTrkVsClsSPIDSlope3"))->Fill(fSpdT, fSpdC1+fSpdC2); // bgID3->slope3
        }
        if (!bgID4) {
            ((TH1F*)fList->FindObject("hTrkVsClsSPIDSlope4"))->Fill(fSpdT, fSpdC1+fSpdC2); // bgID4->slope4
        }
        if (!bgID5) {
            ((TH1F*)fList->FindObject("hTrkVsClsSPIDSlope5"))->Fill(fSpdT, fSpdC1+fSpdC2); // bgID5->slope5
        }
        
        
        for(Int_t ii=1; ii<33; ii++){
            
            //***** Bunch range 3-7 and 11-17 ****/
            //V0 variation 0// Bunch range 3-7, V0 : V0AC, Flag : BB & BG
            SelGoodEvent[0][0][0] = BGFlagA[11]<ii & BGFlagA[12]<ii & BGFlagA[13]<ii & BGFlagA[14]<ii & BGFlagA[15]<ii & BGFlagA[16]<ii & BGFlagA[17]<ii;
            SelGoodEvent[0][0][0] &= BGFlagC[11]<ii & BGFlagC[12]<ii & BGFlagC[13]<ii & BGFlagC[14]<ii & BGFlagC[15]<ii & BGFlagC[16]<ii & BGFlagC[17]<ii;
            SelGoodEvent[0][0][0] &= BGFlagA[7]<ii  & BGFlagA[6]<ii  & BGFlagA[5]<ii  & BGFlagA[4]<ii  & BGFlagA[3]<ii;
            SelGoodEvent[0][0][0] &= BGFlagC[7]<ii  & BGFlagC[6]<ii  & BGFlagC[5]<ii  & BGFlagC[4]<ii  & BGFlagC[3]<ii;
            SelGoodEvent[0][0][0] &= BBFlagA[11]<ii & BBFlagA[12]<ii & BBFlagA[13]<ii & BBFlagA[14]<ii & BBFlagA[15]<ii & BBFlagA[16]<ii & BBFlagA[17]<ii;
            SelGoodEvent[0][0][0] &= BBFlagC[11]<ii & BBFlagC[12]<ii & BBFlagC[13]<ii & BBFlagC[14]<ii & BBFlagC[15]<ii & BBFlagC[16]<ii & BBFlagC[17]<ii;
            SelGoodEvent[0][0][0] &= BBFlagA[7]<ii  & BBFlagA[6]<ii  & BBFlagA[5]<ii  & BBFlagA[4]<ii  & BBFlagA[3]<ii;
            SelGoodEvent[0][0][0] &= BBFlagC[7]<ii  & BBFlagC[6]<ii  & BBFlagC[5]<ii  & BBFlagC[4]<ii  & BBFlagC[3]<ii;
            //V0 variation 1 // Bunch range 3-7, V0 : V0A, Flag : BB & BG
            SelGoodEvent[0][1][0] = BGFlagA[11]<ii & BGFlagA[12]<ii & BGFlagA[13]<ii & BGFlagA[14]<ii & BGFlagA[15]<ii & BGFlagA[16]<ii & BGFlagA[17]<ii;
            SelGoodEvent[0][1][0] &= BGFlagA[7]<ii  & BGFlagA[6]<ii  & BGFlagA[5]<ii  & BGFlagA[4]<ii  & BGFlagA[3]<ii;
            SelGoodEvent[0][1][0] &= BBFlagA[11]<ii & BBFlagA[12]<ii & BBFlagA[13]<ii & BBFlagA[14]<ii & BBFlagA[15]<ii & BBFlagA[16]<ii & BBFlagA[17]<ii;
            SelGoodEvent[0][1][0] &= BBFlagA[7]<ii  & BBFlagA[6]<ii  & BBFlagA[5]<ii  & BBFlagA[4]<ii  & BBFlagA[3]<ii;
            //V0 variation 2 // Bunch range 3-7, V0 : V0C, Flag : BB & BG
            SelGoodEvent[0][2][0] = BGFlagC[11]<ii & BGFlagC[12]<ii & BGFlagC[13]<ii & BGFlagC[14]<ii & BGFlagC[15]<ii & BGFlagC[16]<ii & BGFlagC[17]<ii;
            SelGoodEvent[0][2][0] &= BGFlagC[7]<ii  & BGFlagC[6]<ii  & BGFlagC[5]<ii  & BGFlagC[4]<ii  & BGFlagC[3]<ii;
            SelGoodEvent[0][2][0] &= BBFlagC[11]<ii & BBFlagC[12]<ii & BBFlagC[13]<ii & BBFlagC[14]<ii & BBFlagC[15]<ii & BBFlagC[16]<ii & BBFlagC[17]<ii;
            SelGoodEvent[0][2][0] &= BBFlagC[7]<ii  & BBFlagC[6]<ii  & BBFlagC[5]<ii  & BBFlagC[4]<ii  & BBFlagC[3]<ii;
            //Flag variation : BB// Bunch range 3-7, V0 : V0AC, Flag : BB
            SelGoodEvent[0][0][1] = BBFlagA[11]<ii & BBFlagA[12]<ii & BBFlagA[13]<ii & BBFlagA[14]<ii & BBFlagA[15]<ii & BBFlagA[16]<ii & BBFlagA[17]<ii;
            SelGoodEvent[0][0][1] &= BBFlagC[11]<ii & BBFlagC[12]<ii & BBFlagC[13]<ii & BBFlagC[14]<ii & BBFlagC[15]<ii & BBFlagC[16]<ii & BBFlagC[17]<ii;
            SelGoodEvent[0][0][1] &= BBFlagA[7]<ii  & BBFlagA[6]<ii  & BBFlagA[5]<ii  & BBFlagA[4]<ii  & BBFlagA[3]<ii;
            SelGoodEvent[0][0][1] &= BBFlagC[7]<ii  & BBFlagC[6]<ii  & BBFlagC[5]<ii  & BBFlagC[4]<ii  & BBFlagC[3]<ii;
            // Bunch range 3-7, V0 : V0A, Flag : BB
            SelGoodEvent[0][1][1] = BBFlagA[11]<ii & BBFlagA[12]<ii & BBFlagA[13]<ii & BBFlagA[14]<ii & BBFlagA[15]<ii & BBFlagA[16]<ii & BBFlagA[17]<ii;
            SelGoodEvent[0][1][1] &= BBFlagA[7]<ii  & BBFlagA[6]<ii  & BBFlagA[5]<ii  & BBFlagA[4]<ii  & BBFlagA[3]<ii;
            // Bunch range 3-7, V0 : V0C, Flag : BB
            SelGoodEvent[0][2][1] = BBFlagC[11]<ii & BBFlagC[12]<ii & BBFlagC[13]<ii & BBFlagC[14]<ii & BBFlagC[15]<ii & BBFlagC[16]<ii & BBFlagC[17]<ii;
            SelGoodEvent[0][2][1] &= BBFlagC[7]<ii  & BBFlagC[6]<ii  & BBFlagC[5]<ii  & BBFlagC[4]<ii  & BBFlagC[3]<ii;
            //Flag variation : BG // Bunch range 3-7, V0 : V0AC, Flag : BG
            SelGoodEvent[0][0][2] = BGFlagA[11]<ii & BGFlagA[12]<ii & BGFlagA[13]<ii & BGFlagA[14]<ii & BGFlagA[15]<ii & BGFlagA[16]<ii & BGFlagA[17]<ii;
            SelGoodEvent[0][0][2] &= BGFlagC[11]<ii & BGFlagC[12]<ii & BGFlagC[13]<ii & BGFlagC[14]<ii & BGFlagC[15]<ii & BGFlagC[16]<ii & BGFlagC[17]<ii;
            SelGoodEvent[0][0][2] &= BGFlagA[7]<ii  & BGFlagA[6]<ii  & BGFlagA[5]<ii  & BGFlagA[4]<ii  & BGFlagA[3]<ii;
            SelGoodEvent[0][0][2] &= BGFlagC[7]<ii  & BGFlagC[6]<ii  & BGFlagC[5]<ii  & BGFlagC[4]<ii  & BGFlagC[3]<ii;
            //V0 variation 1 // Bunch range 3-7, V0 : V0A, Flag : BG
            SelGoodEvent[0][1][2] = BGFlagA[11]<ii & BGFlagA[12]<ii & BGFlagA[13]<ii & BGFlagA[14]<ii & BGFlagA[15]<ii & BGFlagA[16]<ii & BGFlagA[17]<ii;
            SelGoodEvent[0][1][2] &= BGFlagA[7]<ii  & BGFlagA[6]<ii  & BGFlagA[5]<ii  & BGFlagA[4]<ii  & BGFlagA[3]<ii;
             //V0 variation 2 // Bunch range 3-7, V0 : V0C, Flag : BG
            SelGoodEvent[0][2][2] = BGFlagC[11]<ii & BGFlagC[12]<ii & BGFlagC[13]<ii & BGFlagC[14]<ii & BGFlagC[15]<ii & BGFlagC[16]<ii & BGFlagC[17]<ii;
            SelGoodEvent[0][2][2] &= BGFlagC[7]<ii  & BGFlagC[6]<ii  & BGFlagC[5]<ii  & BGFlagC[4]<ii  & BGFlagC[3]<ii;

            
            //***** Bunch range 3-8 and 11-17 ****/
            //V0 variation 0// Bunch range 3-8, V0 : V0AC, Flag : BB & BG
            SelGoodEvent[1][0][0] = BGFlagA[11]<ii & BGFlagA[12]<ii & BGFlagA[13]<ii & BGFlagA[14]<ii & BGFlagA[15]<ii & BGFlagA[16]<ii & BGFlagA[17]<ii;
            SelGoodEvent[1][0][0] &= BGFlagC[11]<ii & BGFlagC[12]<ii & BGFlagC[13]<ii & BGFlagC[14]<ii & BGFlagC[15]<ii & BGFlagC[16]<ii & BGFlagC[17]<ii;
            SelGoodEvent[1][0][0] &= BGFlagA[8]<ii  & BGFlagA[7]<ii  & BGFlagA[6]<ii  & BGFlagA[5]<ii  & BGFlagA[4]<ii  & BGFlagA[3]<ii;
            SelGoodEvent[1][0][0] &= BGFlagC[8]<ii  & BGFlagC[7]<ii  & BGFlagC[6]<ii  & BGFlagC[5]<ii  & BGFlagC[4]<ii  & BGFlagC[3]<ii;
            SelGoodEvent[1][0][0] &= BBFlagA[11]<ii & BBFlagA[12]<ii & BBFlagA[13]<ii & BBFlagA[14]<ii & BBFlagA[15]<ii & BBFlagA[16]<ii & BBFlagA[17]<ii;
            SelGoodEvent[1][0][0] &= BBFlagC[11]<ii & BBFlagC[12]<ii & BBFlagC[13]<ii & BBFlagC[14]<ii & BBFlagC[15]<ii & BBFlagC[16]<ii & BBFlagC[17]<ii;
            SelGoodEvent[1][0][0] &= BBFlagA[8]<ii  & BBFlagA[7]<ii  & BBFlagA[6]<ii  & BBFlagA[5]<ii  & BBFlagA[4]<ii  & BBFlagA[3]<ii;
            SelGoodEvent[1][0][0] &= BBFlagC[8]<ii  & BBFlagC[7]<ii  & BBFlagC[6]<ii  & BBFlagC[5]<ii  & BBFlagC[4]<ii  & BBFlagC[3]<ii;
            //V0 variation 1 // Bunch range 3-8, V0 : V0A, Flag : BB & BG
            SelGoodEvent[1][1][0] = BGFlagA[11]<ii & BGFlagA[12]<ii & BGFlagA[13]<ii & BGFlagA[14]<ii & BGFlagA[15]<ii & BGFlagA[16]<ii & BGFlagA[17]<ii;
            SelGoodEvent[1][1][0] &= BGFlagA[8]<ii  & BGFlagA[7]<ii  & BGFlagA[6]<ii  & BGFlagA[5]<ii  & BGFlagA[4]<ii  & BGFlagA[3]<ii;
            SelGoodEvent[1][1][0] &= BBFlagA[11]<ii & BBFlagA[12]<ii & BBFlagA[13]<ii & BBFlagA[14]<ii & BBFlagA[15]<ii & BBFlagA[16]<ii & BBFlagA[17]<ii;
            SelGoodEvent[1][1][0] &= BBFlagA[8]<ii  & BBFlagA[7]<ii  & BBFlagA[6]<ii  & BBFlagA[5]<ii  & BBFlagA[4]<ii  & BBFlagA[3]<ii;
            //V0 variation 2 // Bunch range 3-8, V0 : V0C, Flag : BB & BG
            SelGoodEvent[1][2][0] = BGFlagC[11]<ii & BGFlagC[12]<ii & BGFlagC[13]<ii & BGFlagC[14]<ii & BGFlagC[15]<ii & BGFlagC[16]<ii & BGFlagC[17]<ii;
            SelGoodEvent[1][2][0] &= BGFlagC[8]<ii  & BGFlagC[7]<ii  & BGFlagC[6]<ii  & BGFlagC[5]<ii  & BGFlagC[4]<ii  & BGFlagC[3]<ii;
            SelGoodEvent[1][2][0] &= BBFlagC[11]<ii & BBFlagC[12]<ii & BBFlagC[13]<ii & BBFlagC[14]<ii & BBFlagC[15]<ii & BBFlagC[16]<ii & BBFlagC[17]<ii;
            SelGoodEvent[1][2][0] &= BBFlagC[8]<ii  & BBFlagC[7]<ii  & BBFlagC[6]<ii  & BBFlagC[5]<ii  & BBFlagC[4]<ii  & BBFlagC[3]<ii;
            //Flag variation : BB// Bunch range 3-8, V0 : V0AC, Flag : BB
            SelGoodEvent[1][0][1] = BBFlagA[11]<ii & BBFlagA[12]<ii & BBFlagA[13]<ii & BBFlagA[14]<ii & BBFlagA[15]<ii & BBFlagA[16]<ii & BBFlagA[17]<ii;
            SelGoodEvent[1][0][1] &= BBFlagC[11]<ii & BBFlagC[12]<ii & BBFlagC[13]<ii & BBFlagC[14]<ii & BBFlagC[15]<ii & BBFlagC[16]<ii & BBFlagC[17]<ii;
            SelGoodEvent[1][0][1] &= BBFlagA[8]<ii  & BBFlagA[7]<ii  & BBFlagA[6]<ii  & BBFlagA[5]<ii  & BBFlagA[4]<ii  & BBFlagA[3]<ii;
            SelGoodEvent[1][0][1] &= BBFlagC[8]<ii  & BBFlagC[7]<ii  & BBFlagC[6]<ii  & BBFlagC[5]<ii  & BBFlagC[4]<ii  & BBFlagC[3]<ii;
            // Bunch range 3-8, V0 : V0A, Flag : BB
            SelGoodEvent[1][1][1] = BBFlagA[11]<ii & BBFlagA[12]<ii & BBFlagA[13]<ii & BBFlagA[14]<ii & BBFlagA[15]<ii & BBFlagA[16]<ii & BBFlagA[17]<ii;
            SelGoodEvent[1][1][1] &= BBFlagA[8]<ii  & BBFlagA[7]<ii  & BBFlagA[6]<ii  & BBFlagA[5]<ii  & BBFlagA[4]<ii  & BBFlagA[3]<ii;
            // Bunch range 3-8, V0 : V0C, Flag : BB
            SelGoodEvent[1][2][1] = BBFlagC[11]<ii & BBFlagC[12]<ii & BBFlagC[13]<ii & BBFlagC[14]<ii & BBFlagC[15]<ii & BBFlagC[16]<ii & BBFlagC[17]<ii;
            SelGoodEvent[1][2][1] &= BBFlagC[8]<ii  & BBFlagC[7]<ii  & BBFlagC[6]<ii  & BBFlagC[5]<ii  & BBFlagC[4]<ii  & BBFlagC[3]<ii;
            //Flag variation : BG // Bunch range 3-8, V0 : V0AC, Flag : BG
            SelGoodEvent[1][0][2] = BGFlagA[11]<ii & BGFlagA[12]<ii & BGFlagA[13]<ii & BGFlagA[14]<ii & BGFlagA[15]<ii & BGFlagA[16]<ii & BGFlagA[17]<ii;
            SelGoodEvent[1][0][2] &= BGFlagC[11]<ii & BGFlagC[12]<ii & BGFlagC[13]<ii & BGFlagC[14]<ii & BGFlagC[15]<ii & BGFlagC[16]<ii & BGFlagC[17]<ii;
            SelGoodEvent[1][0][2] &= BGFlagA[8]<ii  & BGFlagA[7]<ii  & BGFlagA[6]<ii  & BGFlagA[5]<ii  & BGFlagA[4]<ii  & BGFlagA[3]<ii;
            SelGoodEvent[1][0][2] &= BGFlagC[8]<ii  & BGFlagC[7]<ii  & BGFlagC[6]<ii  & BGFlagC[5]<ii  & BGFlagC[4]<ii  & BGFlagC[3]<ii;
            //V0 variation 1 // Bunch range 3-8, V0 : V0A, Flag : BG
            SelGoodEvent[1][1][2] = BGFlagA[11]<ii & BGFlagA[12]<ii & BGFlagA[13]<ii & BGFlagA[14]<ii & BGFlagA[15]<ii & BGFlagA[16]<ii & BGFlagA[17]<ii;
            SelGoodEvent[1][1][2] &= BGFlagA[8]<ii  & BGFlagA[7]<ii  & BGFlagA[6]<ii  & BGFlagA[5]<ii  & BGFlagA[4]<ii  & BGFlagA[3]<ii;
            //V0 variation 2 // Bunch range 3-8, V0 : V0C, Flag : BG
            SelGoodEvent[1][2][2] = BGFlagC[11]<ii & BGFlagC[12]<ii & BGFlagC[13]<ii & BGFlagC[14]<ii & BGFlagC[15]<ii & BGFlagC[16]<ii & BGFlagC[17]<ii;
            SelGoodEvent[1][2][2] &= BGFlagC[8]<ii  & BGFlagC[7]<ii  & BGFlagC[6]<ii  & BGFlagC[5]<ii  & BGFlagC[4]<ii  & BGFlagC[3]<ii;
            
            //***** Bunch range 3-9 and 11-17 ****/
            //V0 variation 0// Bunch range 3-9, V0 : V0AC, Flag : BB & BG
            SelGoodEvent[2][0][0] = BGFlagA[11]<ii & BGFlagA[12]<ii & BGFlagA[13]<ii & BGFlagA[14]<ii & BGFlagA[15]<ii & BGFlagA[16]<ii & BGFlagA[17]<ii;
            SelGoodEvent[2][0][0] &= BGFlagC[11]<ii & BGFlagC[12]<ii & BGFlagC[13]<ii & BGFlagC[14]<ii & BGFlagC[15]<ii & BGFlagC[16]<ii & BGFlagC[17]<ii;
            SelGoodEvent[2][0][0] &= BGFlagA[9]<ii  &  BGFlagA[8]<ii  & BGFlagA[7]<ii  & BGFlagA[6]<ii  & BGFlagA[5]<ii  & BGFlagA[4]<ii  & BGFlagA[3]<ii;
            SelGoodEvent[2][0][0] &= BGFlagA[9]<ii  &  BGFlagC[8]<ii  & BGFlagC[7]<ii  & BGFlagC[6]<ii  & BGFlagC[5]<ii  & BGFlagC[4]<ii  & BGFlagC[3]<ii;
            SelGoodEvent[2][0][0] &= BBFlagA[11]<ii & BBFlagA[12]<ii & BBFlagA[13]<ii & BBFlagA[14]<ii & BBFlagA[15]<ii & BBFlagA[16]<ii & BBFlagA[17]<ii;
            SelGoodEvent[2][0][0] &= BBFlagC[11]<ii & BBFlagC[12]<ii & BBFlagC[13]<ii & BBFlagC[14]<ii & BBFlagC[15]<ii & BBFlagC[16]<ii & BBFlagC[17]<ii;
            SelGoodEvent[2][0][0] &= BBFlagA[9]<ii  &  BBFlagA[8]<ii  & BBFlagA[7]<ii  & BBFlagA[6]<ii  & BBFlagA[5]<ii  & BBFlagA[4]<ii  & BBFlagA[3]<ii;
            SelGoodEvent[2][0][0] &= BBFlagC[9]<ii  &  BBFlagC[8]<ii  & BBFlagC[7]<ii  & BBFlagC[6]<ii  & BBFlagC[5]<ii  & BBFlagC[4]<ii  & BBFlagC[3]<ii;
            //V0 variation 1 // Bunch range 3-9, V0 : V0A, Flag : BB & BG
            SelGoodEvent[2][1][0] = BGFlagA[11]<ii & BGFlagA[12]<ii & BGFlagA[13]<ii & BGFlagA[14]<ii & BGFlagA[15]<ii & BGFlagA[16]<ii & BGFlagA[17]<ii;
            SelGoodEvent[2][1][0] &= BGFlagA[9]<ii  &  BGFlagA[8]<ii  & BGFlagA[7]<ii  & BGFlagA[6]<ii  & BGFlagA[5]<ii  & BGFlagA[4]<ii  & BGFlagA[3]<ii;
            SelGoodEvent[2][1][0] &= BBFlagA[11]<ii & BBFlagA[12]<ii & BBFlagA[13]<ii & BBFlagA[14]<ii & BBFlagA[15]<ii & BBFlagA[16]<ii & BBFlagA[17]<ii;
            SelGoodEvent[2][1][0] &= BBFlagA[9]<ii  &  BBFlagA[8]<ii  & BBFlagA[7]<ii  & BBFlagA[6]<ii  & BBFlagA[5]<ii  & BBFlagA[4]<ii  & BBFlagA[3]<ii;
            //V0 variation 2 // Bunch range 3-9, V0 : V0C, Flag : BB & BG
            SelGoodEvent[2][2][0] = BGFlagC[11]<ii & BGFlagC[12]<ii & BGFlagC[13]<ii & BGFlagC[14]<ii & BGFlagC[15]<ii & BGFlagC[16]<ii & BGFlagC[17]<ii;
            SelGoodEvent[2][2][0] &= BGFlagC[9]<ii  &  BGFlagC[8]<ii  & BGFlagC[7]<ii  & BGFlagC[6]<ii  & BGFlagC[5]<ii  & BGFlagC[4]<ii  & BGFlagC[3]<ii;
            SelGoodEvent[2][2][0] &= BBFlagC[11]<ii & BBFlagC[12]<ii & BBFlagC[13]<ii & BBFlagC[14]<ii & BBFlagC[15]<ii & BBFlagC[16]<ii & BBFlagC[17]<ii;
            SelGoodEvent[2][2][0] &= BBFlagC[9]<ii  &  BBFlagC[8]<ii  & BBFlagC[7]<ii  & BBFlagC[6]<ii  & BBFlagC[5]<ii  & BBFlagC[4]<ii  & BBFlagC[3]<ii;
            //Flag variation : BB// Bunch range 3-9, V0 : V0AC, Flag : BB
            SelGoodEvent[2][0][1] = BBFlagA[11]<ii & BBFlagA[12]<ii & BBFlagA[13]<ii & BBFlagA[14]<ii & BBFlagA[15]<ii & BBFlagA[16]<ii & BBFlagA[17]<ii;
            SelGoodEvent[2][0][1] &= BBFlagC[11]<ii & BBFlagC[12]<ii & BBFlagC[13]<ii & BBFlagC[14]<ii & BBFlagC[15]<ii & BBFlagC[16]<ii & BBFlagC[17]<ii;
            SelGoodEvent[2][0][1] &= BBFlagA[9]<ii  &  BBFlagA[8]<ii  & BBFlagA[7]<ii  & BBFlagA[6]<ii  & BBFlagA[5]<ii  & BBFlagA[4]<ii  & BBFlagA[3]<ii;
            SelGoodEvent[2][0][1] &= BBFlagC[9]<ii  &  BBFlagC[8]<ii  & BBFlagC[7]<ii  & BBFlagC[6]<ii  & BBFlagC[5]<ii  & BBFlagC[4]<ii  & BBFlagC[3]<ii;
            // Bunch range 3-9, V0 : V0A, Flag : BB
            SelGoodEvent[2][1][1] = BBFlagA[11]<ii & BBFlagA[12]<ii & BBFlagA[13]<ii & BBFlagA[14]<ii & BBFlagA[15]<ii & BBFlagA[16]<ii & BBFlagA[17]<ii;
            SelGoodEvent[2][1][1] &= BBFlagA[9]<ii  &  BBFlagA[8]<ii  & BBFlagA[7]<ii  & BBFlagA[6]<ii  & BBFlagA[5]<ii  & BBFlagA[4]<ii  & BBFlagA[3]<ii;
            // Bunch range 3-9, V0 : V0C, Flag : BB
            SelGoodEvent[2][2][1] = BBFlagC[11]<ii & BBFlagC[12]<ii & BBFlagC[13]<ii & BBFlagC[14]<ii & BBFlagC[15]<ii & BBFlagC[16]<ii & BBFlagC[17]<ii;
            SelGoodEvent[2][2][1] &= BBFlagC[9]<ii  &  BBFlagC[8]<ii  & BBFlagC[7]<ii  & BBFlagC[6]<ii  & BBFlagC[5]<ii  & BBFlagC[4]<ii  & BBFlagC[3]<ii;
            //Flag variation : BG // Bunch range 3-9, V0 : V0AC, Flag : BG
            SelGoodEvent[2][0][2] = BGFlagA[11]<ii & BGFlagA[12]<ii & BGFlagA[13]<ii & BGFlagA[14]<ii & BGFlagA[15]<ii & BGFlagA[16]<ii & BGFlagA[17]<ii;
            SelGoodEvent[2][0][2] &= BGFlagC[11]<ii & BGFlagC[12]<ii & BGFlagC[13]<ii & BGFlagC[14]<ii & BGFlagC[15]<ii & BGFlagC[16]<ii & BGFlagC[17]<ii;
            SelGoodEvent[2][0][2] &= BGFlagA[9]<ii  &  BGFlagA[8]<ii  & BGFlagA[7]<ii  & BGFlagA[6]<ii  & BGFlagA[5]<ii  & BGFlagA[4]<ii  & BGFlagA[3]<ii;
            SelGoodEvent[2][0][2] &= BGFlagC[9]<ii  &  BGFlagC[8]<ii  & BGFlagC[7]<ii  & BGFlagC[6]<ii  & BGFlagC[5]<ii  & BGFlagC[4]<ii  & BGFlagC[3]<ii;
            //V0 variation 1 // Bunch range 3-9, V0 : V0A, Flag : BG
            SelGoodEvent[2][1][2] = BGFlagA[11]<ii & BGFlagA[12]<ii & BGFlagA[13]<ii & BGFlagA[14]<ii & BGFlagA[15]<ii & BGFlagA[16]<ii & BGFlagA[17]<ii;
            SelGoodEvent[2][1][2] &= BGFlagA[9]<ii  &  BGFlagA[8]<ii  & BGFlagA[7]<ii  & BGFlagA[6]<ii  & BGFlagA[5]<ii  & BGFlagA[4]<ii  & BGFlagA[3]<ii;
            //V0 variation 2 // Bunch range 3-9, V0 : V0C, Flag : BG
            SelGoodEvent[2][2][2] = BGFlagC[11]<ii & BGFlagC[12]<ii & BGFlagC[13]<ii & BGFlagC[14]<ii & BGFlagC[15]<ii & BGFlagC[16]<ii & BGFlagC[17]<ii;
            SelGoodEvent[2][2][2] &= BGFlagC[9]<ii  &  BGFlagC[8]<ii  & BGFlagC[7]<ii  & BGFlagC[6]<ii  & BGFlagC[5]<ii  & BGFlagC[4]<ii  & BGFlagC[3]<ii;
            
            
            double check1=0;
            double check2=0;
            double check3=0;
            if(SelGoodEvent[2][2][0]) check1++;
            if(SelGoodEvent[2][2][1]) check2++;
            if(SelGoodEvent[2][2][2]) check3++;

            
            cout<< "3-9, v0c , bb+bg = " <<SelGoodEvent[2][2][0]<<",  num 1 = "<<check1<<endl;
            cout<< "3-9, v0c , bb = " <<SelGoodEvent[2][2][1]<<",  num 2 = "<<check2<<endl;
            cout<< "3-9, v0c , bg = " <<SelGoodEvent[2][2][2]<<",  num 3 = "<<check3<<endl;

            
            for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                    for(int k=0; k<3; k++){
                        
                        if(SelGoodEvent[i][j][k]) {
                            ((TH1F*)fList->FindObject(Form("hDenomPurityBC%d_V0%d_Flag%d",i,j,k)))->Fill(ii-1);
                        }
                        
                        if(!bgID) {
                            ((TH1F*)fList->FindObject(Form("hDenomEffBC%d_V0%d_Flag%d",i,j,k)))->Fill(ii-1);
                            ((TH1F*)fList->FindObject(Form("hSPDDenomBC%d_V0%d_Flag%d",i,j,k)))->Fill(fSpdT);
                            ((TH1F*)fList->FindObject(Form("hDenomTrkVsClsSPID%d_V0%d_Flag%d",i,j,k)))->Fill(fSpdT, fSpdC1+fSpdC2);
                            ((TH1F*)fList->FindObject(Form("hDenomV0%d_V0%d_Flag%d",i,j,k)))->Fill(fv0a-fv0c, fv0a+fv0c);

                            
                            if(SelGoodEvent[i][j][k]){
                                ((TH1F*)fList->FindObject(Form("hNumEffPurityBC%d_V0%d_Flag%d",i,j,k)))->Fill(ii-1);
                                ((TH1F*)fList->FindObject(Form("hSPDNumBC%d_V0%d_Flag%d",i,j,k)))->Fill(fSpdT);
                                ((TH1F*)fList->FindObject(Form("hNumTrkVsClsSPID%d_V0%d_Flag%d",i,j,k)))->Fill(fSpdT, fSpdC1+fSpdC2);
                                ((TH1F*)fList->FindObject(Form("hNumV0%d_V0%d_Flag%d",i,j,k)))->Fill(fv0a-fv0c, fv0a+fv0c);

                            }
                        }
                        if(bgID){
                            ((TH1F*)fList->FindObject(Form("hDenomRejecEffBC%d_V0%d_Flag%d",i,j,k)))->Fill(ii-1);
                            if(!SelGoodEvent[i][j][k]){
                                ((TH1F*)fList->FindObject(Form("hNumRejecEffBC%d_V0%d_Flag%d",i,j,k)))->Fill(ii-1);
                                
                            }
                        }
                    }
                }
            } // end of fill histograms
                          } // end of V0 flag loop
        
        for(Int_t ii=1; ii<9; ii++){
            
            //***** Bunch range 3-7 and 11-17 ****/
            //V0 variation 0// Bunch range 3-7, V0 : V0AC, Flag : BB & BG
            SelGoodEventAD[0][0][0] = ADBGFlagA[11]<ii & ADBGFlagA[12]<ii & ADBGFlagA[13]<ii & ADBGFlagA[14]<ii & ADBGFlagA[15]<ii & ADBGFlagA[16]<ii & ADBGFlagA[17]<ii;
            SelGoodEventAD[0][0][0] &= ADBGFlagC[11]<ii & ADBGFlagC[12]<ii & ADBGFlagC[13]<ii & ADBGFlagC[14]<ii & ADBGFlagC[15]<ii & ADBGFlagC[16]<ii & ADBGFlagC[17]<ii;
            SelGoodEventAD[0][0][0] &= ADBGFlagA[7]<ii  & ADBGFlagA[6]<ii  & ADBGFlagA[5]<ii  & ADBGFlagA[4]<ii  & ADBGFlagA[3]<ii;
            SelGoodEventAD[0][0][0] &= ADBGFlagC[7]<ii  & ADBGFlagC[6]<ii  & ADBGFlagC[5]<ii  & ADBGFlagC[4]<ii  & ADBGFlagC[3]<ii;
            SelGoodEventAD[0][0][0] &= ADBBFlagA[11]<ii & ADBBFlagA[12]<ii & ADBBFlagA[13]<ii & ADBBFlagA[14]<ii & ADBBFlagA[15]<ii & ADBBFlagA[16]<ii & ADBBFlagA[17]<ii;
            SelGoodEventAD[0][0][0] &= ADBBFlagC[11]<ii & ADBBFlagC[12]<ii & ADBBFlagC[13]<ii & ADBBFlagC[14]<ii & ADBBFlagC[15]<ii & ADBBFlagC[16]<ii & ADBBFlagC[17]<ii;
            SelGoodEventAD[0][0][0] &= ADBBFlagA[7]<ii  & ADBBFlagA[6]<ii  & ADBBFlagA[5]<ii  & ADBBFlagA[4]<ii  & ADBBFlagA[3]<ii;
            SelGoodEventAD[0][0][0] &= ADBBFlagC[7]<ii  & ADBBFlagC[6]<ii  & ADBBFlagC[5]<ii  & ADBBFlagC[4]<ii  & ADBBFlagC[3]<ii;
            //V0 variation 1 // Bunch range 3-7, V0 : V0A, Flag : BB & ADBG
            SelGoodEventAD[0][1][0] = ADBGFlagA[11]<ii & ADBGFlagA[12]<ii & ADBGFlagA[13]<ii & ADBGFlagA[14]<ii & ADBGFlagA[15]<ii & ADBGFlagA[16]<ii & ADBGFlagA[17]<ii;
            SelGoodEventAD[0][1][0] &= ADBGFlagA[7]<ii  & ADBGFlagA[6]<ii  & ADBGFlagA[5]<ii  & ADBGFlagA[4]<ii  & ADBGFlagA[3]<ii;
            SelGoodEventAD[0][1][0] &= ADBBFlagA[11]<ii & ADBBFlagA[12]<ii & ADBBFlagA[13]<ii & ADBBFlagA[14]<ii & ADBBFlagA[15]<ii & ADBBFlagA[16]<ii & ADBBFlagA[17]<ii;
            SelGoodEventAD[0][1][0] &= ADBBFlagA[7]<ii  & ADBBFlagA[6]<ii  & ADBBFlagA[5]<ii  & ADBBFlagA[4]<ii  & ADBBFlagA[3]<ii;
            //V0 variation 2 // Bunch range 3-7, V0 : V0C, Flag : BB & ADBG
            SelGoodEventAD[0][2][0] = ADBGFlagC[11]<ii & ADBGFlagC[12]<ii & ADBGFlagC[13]<ii & ADBGFlagC[14]<ii & ADBGFlagC[15]<ii & ADBGFlagC[16]<ii & ADBGFlagC[17]<ii;
            SelGoodEventAD[0][2][0] &= ADBGFlagC[7]<ii  & ADBGFlagC[6]<ii  & ADBGFlagC[5]<ii  & ADBGFlagC[4]<ii  & ADBGFlagC[3]<ii;
            SelGoodEventAD[0][2][0] &= ADBBFlagC[11]<ii & ADBBFlagC[12]<ii & ADBBFlagC[13]<ii & ADBBFlagC[14]<ii & ADBBFlagC[15]<ii & ADBBFlagC[16]<ii & ADBBFlagC[17]<ii;
            SelGoodEventAD[0][2][0] &= ADBBFlagC[7]<ii  & ADBBFlagC[6]<ii  & ADBBFlagC[5]<ii  & ADBBFlagC[4]<ii  & ADBBFlagC[3]<ii;
            //Flag variation : BB// Bunch range 3-7, V0 : V0AC, Flag : BB
            SelGoodEventAD[0][0][1] = ADBBFlagA[11]<ii & ADBBFlagA[12]<ii & ADBBFlagA[13]<ii & ADBBFlagA[14]<ii & ADBBFlagA[15]<ii & ADBBFlagA[16]<ii & ADBBFlagA[17]<ii;
            SelGoodEventAD[0][0][1] &= ADBBFlagC[11]<ii & ADBBFlagC[12]<ii & ADBBFlagC[13]<ii & ADBBFlagC[14]<ii & ADBBFlagC[15]<ii & ADBBFlagC[16]<ii & ADBBFlagC[17]<ii;
            SelGoodEventAD[0][0][1] &= ADBBFlagA[7]<ii  & ADBBFlagA[6]<ii  & ADBBFlagA[5]<ii  & ADBBFlagA[4]<ii  & ADBBFlagA[3]<ii;
            SelGoodEventAD[0][0][1] &= ADBBFlagC[7]<ii  & ADBBFlagC[6]<ii  & ADBBFlagC[5]<ii  & ADBBFlagC[4]<ii  & ADBBFlagC[3]<ii;
            // Bunch range 3-7, V0 : V0A, Flag : BB
            SelGoodEventAD[0][1][1] = ADBBFlagA[11]<ii & ADBBFlagA[12]<ii & ADBBFlagA[13]<ii & ADBBFlagA[14]<ii & ADBBFlagA[15]<ii & ADBBFlagA[16]<ii & ADBBFlagA[17]<ii;
            SelGoodEventAD[0][1][1] &= ADBBFlagA[7]<ii  & ADBBFlagA[6]<ii  & ADBBFlagA[5]<ii  & ADBBFlagA[4]<ii  & ADBBFlagA[3]<ii;
            // Bunch range 3-7, V0 : V0C, Flag : BB
            SelGoodEventAD[0][2][1] = ADBBFlagC[11]<ii & ADBBFlagC[12]<ii & ADBBFlagC[13]<ii & ADBBFlagC[14]<ii & ADBBFlagC[15]<ii & ADBBFlagC[16]<ii & ADBBFlagC[17]<ii;
            SelGoodEventAD[0][2][1] &= ADBBFlagC[7]<ii  & ADBBFlagC[6]<ii  & ADBBFlagC[5]<ii  & ADBBFlagC[4]<ii  & ADBBFlagC[3]<ii;
            //Flag variation : ADBG // Bunch range 3-7, V0 : V0AC, Flag : ADBG
            SelGoodEventAD[0][0][2] = ADBGFlagA[11]<ii & ADBGFlagA[12]<ii & ADBGFlagA[13]<ii & ADBGFlagA[14]<ii & ADBGFlagA[15]<ii & ADBGFlagA[16]<ii & ADBGFlagA[17]<ii;
            SelGoodEventAD[0][0][2] &= ADBGFlagC[11]<ii & ADBGFlagC[12]<ii & ADBGFlagC[13]<ii & ADBGFlagC[14]<ii & ADBGFlagC[15]<ii & ADBGFlagC[16]<ii & ADBGFlagC[17]<ii;
            SelGoodEventAD[0][0][2] &= ADBGFlagA[7]<ii  & ADBGFlagA[6]<ii  & ADBGFlagA[5]<ii  & ADBGFlagA[4]<ii  & ADBGFlagA[3]<ii;
            SelGoodEventAD[0][0][2] &= ADBGFlagC[7]<ii  & ADBGFlagC[6]<ii  & ADBGFlagC[5]<ii  & ADBGFlagC[4]<ii  & ADBGFlagC[3]<ii;
            //V0 variation 1 // Bunch range 3-7, V0 : V0A, Flag : ADBG
            SelGoodEventAD[0][1][2] = ADBGFlagA[11]<ii & ADBGFlagA[12]<ii & ADBGFlagA[13]<ii & ADBGFlagA[14]<ii & ADBGFlagA[15]<ii & ADBGFlagA[16]<ii & ADBGFlagA[17]<ii;
            SelGoodEventAD[0][1][2] &= ADBGFlagA[7]<ii  & ADBGFlagA[6]<ii  & ADBGFlagA[5]<ii  & ADBGFlagA[4]<ii  & ADBGFlagA[3]<ii;
            //V0 variation 2 // Bunch range 3-7, V0 : V0C, Flag : ADBG
            SelGoodEventAD[0][2][2] = ADBGFlagC[11]<ii & ADBGFlagC[12]<ii & ADBGFlagC[13]<ii & ADBGFlagC[14]<ii & ADBGFlagC[15]<ii & ADBGFlagC[16]<ii & ADBGFlagC[17]<ii;
            SelGoodEventAD[0][2][2] &= ADBGFlagC[7]<ii  & ADBGFlagC[6]<ii  & ADBGFlagC[5]<ii  & ADBGFlagC[4]<ii  & ADBGFlagC[3]<ii;
            
            
            //***** Bunch range 3-8 and 11-17 ****/
            //V0 variation 0// Bunch range 3-8, V0 : V0AC, Flag : BB & ADBG
            SelGoodEventAD[1][0][0] = ADBGFlagA[11]<ii & ADBGFlagA[12]<ii & ADBGFlagA[13]<ii & ADBGFlagA[14]<ii & ADBGFlagA[15]<ii & ADBGFlagA[16]<ii & ADBGFlagA[17]<ii;
            SelGoodEventAD[1][0][0] &= ADBGFlagC[11]<ii & ADBGFlagC[12]<ii & ADBGFlagC[13]<ii & ADBGFlagC[14]<ii & ADBGFlagC[15]<ii & ADBGFlagC[16]<ii & ADBGFlagC[17]<ii;
            SelGoodEventAD[1][0][0] &= ADBGFlagA[8]<ii  & ADBGFlagA[7]<ii  & ADBGFlagA[6]<ii  & ADBGFlagA[5]<ii  & ADBGFlagA[4]<ii  & ADBGFlagA[3]<ii;
            SelGoodEventAD[1][0][0] &= ADBGFlagC[8]<ii  & ADBGFlagC[7]<ii  & ADBGFlagC[6]<ii  & ADBGFlagC[5]<ii  & ADBGFlagC[4]<ii  & ADBGFlagC[3]<ii;
            SelGoodEventAD[1][0][0] &= ADBBFlagA[11]<ii & ADBBFlagA[12]<ii & ADBBFlagA[13]<ii & ADBBFlagA[14]<ii & ADBBFlagA[15]<ii & ADBBFlagA[16]<ii & ADBBFlagA[17]<ii;
            SelGoodEventAD[1][0][0] &= ADBBFlagC[11]<ii & ADBBFlagC[12]<ii & ADBBFlagC[13]<ii & ADBBFlagC[14]<ii & ADBBFlagC[15]<ii & ADBBFlagC[16]<ii & ADBBFlagC[17]<ii;
            SelGoodEventAD[1][0][0] &= ADBBFlagA[8]<ii  & ADBBFlagA[7]<ii  & ADBBFlagA[6]<ii  & ADBBFlagA[5]<ii  & ADBBFlagA[4]<ii  & ADBBFlagA[3]<ii;
            SelGoodEventAD[1][0][0] &= ADBBFlagC[8]<ii  & ADBBFlagC[7]<ii  & ADBBFlagC[6]<ii  & ADBBFlagC[5]<ii  & ADBBFlagC[4]<ii  & ADBBFlagC[3]<ii;
            //V0 variation 1 // Bunch range 3-8, V0 : V0A, Flag : BB & ADBG
            SelGoodEventAD[1][1][0] = ADBGFlagA[11]<ii & ADBGFlagA[12]<ii & ADBGFlagA[13]<ii & ADBGFlagA[14]<ii & ADBGFlagA[15]<ii & ADBGFlagA[16]<ii & ADBGFlagA[17]<ii;
            SelGoodEventAD[1][1][0] &= ADBGFlagA[8]<ii  & ADBGFlagA[7]<ii  & ADBGFlagA[6]<ii  & ADBGFlagA[5]<ii  & ADBGFlagA[4]<ii  & ADBGFlagA[3]<ii;
            SelGoodEventAD[1][1][0] &= ADBBFlagA[11]<ii & ADBBFlagA[12]<ii & ADBBFlagA[13]<ii & ADBBFlagA[14]<ii & ADBBFlagA[15]<ii & ADBBFlagA[16]<ii & ADBBFlagA[17]<ii;
            SelGoodEventAD[1][1][0] &= ADBBFlagA[8]<ii  & ADBBFlagA[7]<ii  & ADBBFlagA[6]<ii  & ADBBFlagA[5]<ii  & ADBBFlagA[4]<ii  & ADBBFlagA[3]<ii;
            //V0 variation 2 // Bunch range 3-8, V0 : V0C, Flag : BB & ADBG
            SelGoodEventAD[1][2][0] = ADBGFlagC[11]<ii & ADBGFlagC[12]<ii & ADBGFlagC[13]<ii & ADBGFlagC[14]<ii & ADBGFlagC[15]<ii & ADBGFlagC[16]<ii & ADBGFlagC[17]<ii;
            SelGoodEventAD[1][2][0] &= ADBGFlagC[8]<ii  & ADBGFlagC[7]<ii  & ADBGFlagC[6]<ii  & ADBGFlagC[5]<ii  & ADBGFlagC[4]<ii  & ADBGFlagC[3]<ii;
            SelGoodEventAD[1][2][0] &= ADBBFlagC[11]<ii & ADBBFlagC[12]<ii & ADBBFlagC[13]<ii & ADBBFlagC[14]<ii & ADBBFlagC[15]<ii & ADBBFlagC[16]<ii & ADBBFlagC[17]<ii;
            SelGoodEventAD[1][2][0] &= ADBBFlagC[8]<ii  & ADBBFlagC[7]<ii  & ADBBFlagC[6]<ii  & ADBBFlagC[5]<ii  & ADBBFlagC[4]<ii  & ADBBFlagC[3]<ii;
            //Flag variation : BB// Bunch range 3-8, V0 : V0AC, Flag : BB
            SelGoodEventAD[1][0][1] = ADBBFlagA[11]<ii & ADBBFlagA[12]<ii & ADBBFlagA[13]<ii & ADBBFlagA[14]<ii & ADBBFlagA[15]<ii & ADBBFlagA[16]<ii & ADBBFlagA[17]<ii;
            SelGoodEventAD[1][0][1] &= ADBBFlagC[11]<ii & ADBBFlagC[12]<ii & ADBBFlagC[13]<ii & ADBBFlagC[14]<ii & ADBBFlagC[15]<ii & ADBBFlagC[16]<ii & ADBBFlagC[17]<ii;
            SelGoodEventAD[1][0][1] &= ADBBFlagA[8]<ii  & ADBBFlagA[7]<ii  & ADBBFlagA[6]<ii  & ADBBFlagA[5]<ii  & ADBBFlagA[4]<ii  & ADBBFlagA[3]<ii;
            SelGoodEventAD[1][0][1] &= ADBBFlagC[8]<ii  & ADBBFlagC[7]<ii  & ADBBFlagC[6]<ii  & ADBBFlagC[5]<ii  & ADBBFlagC[4]<ii  & ADBBFlagC[3]<ii;
            // Bunch range 3-8, V0 : V0A, Flag : BB
            SelGoodEventAD[1][1][1] = ADBBFlagA[11]<ii & ADBBFlagA[12]<ii & ADBBFlagA[13]<ii & ADBBFlagA[14]<ii & ADBBFlagA[15]<ii & ADBBFlagA[16]<ii & ADBBFlagA[17]<ii;
            SelGoodEventAD[1][1][1] &= ADBBFlagA[8]<ii  & ADBBFlagA[7]<ii  & ADBBFlagA[6]<ii  & ADBBFlagA[5]<ii  & ADBBFlagA[4]<ii  & ADBBFlagA[3]<ii;
            // Bunch range 3-8, V0 : V0C, Flag : BB
            SelGoodEventAD[1][2][1] = ADBBFlagC[11]<ii & ADBBFlagC[12]<ii & ADBBFlagC[13]<ii & ADBBFlagC[14]<ii & ADBBFlagC[15]<ii & ADBBFlagC[16]<ii & ADBBFlagC[17]<ii;
            SelGoodEventAD[1][2][1] &= ADBBFlagC[8]<ii  & ADBBFlagC[7]<ii  & ADBBFlagC[6]<ii  & ADBBFlagC[5]<ii  & ADBBFlagC[4]<ii  & ADBBFlagC[3]<ii;
            //Flag variation : ADBG // Bunch range 3-8, V0 : V0AC, Flag : ADBG
            SelGoodEventAD[1][0][2] = ADBGFlagA[11]<ii & ADBGFlagA[12]<ii & ADBGFlagA[13]<ii & ADBGFlagA[14]<ii & ADBGFlagA[15]<ii & ADBGFlagA[16]<ii & ADBGFlagA[17]<ii;
            SelGoodEventAD[1][0][2] &= ADBGFlagC[11]<ii & ADBGFlagC[12]<ii & ADBGFlagC[13]<ii & ADBGFlagC[14]<ii & ADBGFlagC[15]<ii & ADBGFlagC[16]<ii & ADBGFlagC[17]<ii;
            SelGoodEventAD[1][0][2] &= ADBGFlagA[8]<ii  & ADBGFlagA[7]<ii  & ADBGFlagA[6]<ii  & ADBGFlagA[5]<ii  & ADBGFlagA[4]<ii  & ADBGFlagA[3]<ii;
            SelGoodEventAD[1][0][2] &= ADBGFlagC[8]<ii  & ADBGFlagC[7]<ii  & ADBGFlagC[6]<ii  & ADBGFlagC[5]<ii  & ADBGFlagC[4]<ii  & ADBGFlagC[3]<ii;
            //V0 variation 1 // Bunch range 3-8, V0 : V0A, Flag : ADBG
            SelGoodEventAD[1][1][2] = ADBGFlagA[11]<ii & ADBGFlagA[12]<ii & ADBGFlagA[13]<ii & ADBGFlagA[14]<ii & ADBGFlagA[15]<ii & ADBGFlagA[16]<ii & ADBGFlagA[17]<ii;
            SelGoodEventAD[1][1][2] &= ADBGFlagA[8]<ii  & ADBGFlagA[7]<ii  & ADBGFlagA[6]<ii  & ADBGFlagA[5]<ii  & ADBGFlagA[4]<ii  & ADBGFlagA[3]<ii;
            //V0 variation 2 // Bunch range 3-8, V0 : V0C, Flag : ADBG
            SelGoodEventAD[1][2][2] = ADBGFlagC[11]<ii & ADBGFlagC[12]<ii & ADBGFlagC[13]<ii & ADBGFlagC[14]<ii & ADBGFlagC[15]<ii & ADBGFlagC[16]<ii & ADBGFlagC[17]<ii;
            SelGoodEventAD[1][2][2] &= ADBGFlagC[8]<ii  & ADBGFlagC[7]<ii  & ADBGFlagC[6]<ii  & ADBGFlagC[5]<ii  & ADBGFlagC[4]<ii  & ADBGFlagC[3]<ii;
            
            //***** Bunch range 3-9 and 11-17 ****/
            //V0 variation 0// Bunch range 3-9, V0 : V0AC, Flag : BB & ADBG
            SelGoodEventAD[2][0][0] = ADBGFlagA[11]<ii & ADBGFlagA[12]<ii & ADBGFlagA[13]<ii & ADBGFlagA[14]<ii & ADBGFlagA[15]<ii & ADBGFlagA[16]<ii & ADBGFlagA[17]<ii;
            SelGoodEventAD[2][0][0] &= ADBGFlagC[11]<ii & ADBGFlagC[12]<ii & ADBGFlagC[13]<ii & ADBGFlagC[14]<ii & ADBGFlagC[15]<ii & ADBGFlagC[16]<ii & ADBGFlagC[17]<ii;
            SelGoodEventAD[2][0][0] &= ADBGFlagA[9]<ii  &  ADBGFlagA[8]<ii  & ADBGFlagA[7]<ii  & ADBGFlagA[6]<ii  & ADBGFlagA[5]<ii  & ADBGFlagA[4]<ii  & ADBGFlagA[3]<ii;
            SelGoodEventAD[2][0][0] &= ADBGFlagA[9]<ii  &  ADBGFlagC[8]<ii  & ADBGFlagC[7]<ii  & ADBGFlagC[6]<ii  & ADBGFlagC[5]<ii  & ADBGFlagC[4]<ii  & ADBGFlagC[3]<ii;
            SelGoodEventAD[2][0][0] &= ADBBFlagA[11]<ii & ADBBFlagA[12]<ii & ADBBFlagA[13]<ii & ADBBFlagA[14]<ii & ADBBFlagA[15]<ii & ADBBFlagA[16]<ii & ADBBFlagA[17]<ii;
            SelGoodEventAD[2][0][0] &= ADBBFlagC[11]<ii & ADBBFlagC[12]<ii & ADBBFlagC[13]<ii & ADBBFlagC[14]<ii & ADBBFlagC[15]<ii & ADBBFlagC[16]<ii & ADBBFlagC[17]<ii;
            SelGoodEventAD[2][0][0] &= ADBBFlagA[9]<ii  &  ADBBFlagA[8]<ii  & ADBBFlagA[7]<ii  & ADBBFlagA[6]<ii  & ADBBFlagA[5]<ii  & ADBBFlagA[4]<ii  & ADBBFlagA[3]<ii;
            SelGoodEventAD[2][0][0] &= ADBBFlagC[9]<ii  &  ADBBFlagC[8]<ii  & ADBBFlagC[7]<ii  & ADBBFlagC[6]<ii  & ADBBFlagC[5]<ii  & ADBBFlagC[4]<ii  & ADBBFlagC[3]<ii;
            //V0 variation 1 // Bunch range 3-9, V0 : V0A, Flag : BB & ADBG
            SelGoodEventAD[2][1][0] = ADBGFlagA[11]<ii & ADBGFlagA[12]<ii & ADBGFlagA[13]<ii & ADBGFlagA[14]<ii & ADBGFlagA[15]<ii & ADBGFlagA[16]<ii & ADBGFlagA[17]<ii;
            SelGoodEventAD[2][1][0] &= ADBGFlagA[9]<ii  &  ADBGFlagA[8]<ii  & ADBGFlagA[7]<ii  & ADBGFlagA[6]<ii  & ADBGFlagA[5]<ii  & ADBGFlagA[4]<ii  & ADBGFlagA[3]<ii;
            SelGoodEventAD[2][1][0] &= ADBBFlagA[11]<ii & ADBBFlagA[12]<ii & ADBBFlagA[13]<ii & ADBBFlagA[14]<ii & ADBBFlagA[15]<ii & ADBBFlagA[16]<ii & ADBBFlagA[17]<ii;
            SelGoodEventAD[2][1][0] &= ADBBFlagA[9]<ii  &  ADBBFlagA[8]<ii  & ADBBFlagA[7]<ii  & ADBBFlagA[6]<ii  & ADBBFlagA[5]<ii  & ADBBFlagA[4]<ii  & ADBBFlagA[3]<ii;
            //V0 variation 2 // Bunch range 3-9, V0 : V0C, Flag : BB & ADBG
            SelGoodEventAD[2][2][0] = ADBGFlagC[11]<ii & ADBGFlagC[12]<ii & ADBGFlagC[13]<ii & ADBGFlagC[14]<ii & ADBGFlagC[15]<ii & ADBGFlagC[16]<ii & ADBGFlagC[17]<ii;
            SelGoodEventAD[2][2][0] &= ADBGFlagC[9]<ii  &  ADBGFlagC[8]<ii  & ADBGFlagC[7]<ii  & ADBGFlagC[6]<ii  & ADBGFlagC[5]<ii  & ADBGFlagC[4]<ii  & ADBGFlagC[3]<ii;
            SelGoodEventAD[2][2][0] &= ADBBFlagC[11]<ii & ADBBFlagC[12]<ii & ADBBFlagC[13]<ii & ADBBFlagC[14]<ii & ADBBFlagC[15]<ii & ADBBFlagC[16]<ii & ADBBFlagC[17]<ii;
            SelGoodEventAD[2][2][0] &= ADBBFlagC[9]<ii  &  ADBBFlagC[8]<ii  & ADBBFlagC[7]<ii  & ADBBFlagC[6]<ii  & ADBBFlagC[5]<ii  & ADBBFlagC[4]<ii  & ADBBFlagC[3]<ii;
            //Flag variation : BB// Bunch range 3-9, V0 : V0AC, Flag : BB
            SelGoodEventAD[2][0][1] = ADBBFlagA[11]<ii & ADBBFlagA[12]<ii & ADBBFlagA[13]<ii & ADBBFlagA[14]<ii & ADBBFlagA[15]<ii & ADBBFlagA[16]<ii & ADBBFlagA[17]<ii;
            SelGoodEventAD[2][0][1] &= ADBBFlagC[11]<ii & ADBBFlagC[12]<ii & ADBBFlagC[13]<ii & ADBBFlagC[14]<ii & ADBBFlagC[15]<ii & ADBBFlagC[16]<ii & ADBBFlagC[17]<ii;
            SelGoodEventAD[2][0][1] &= ADBBFlagA[9]<ii  &  ADBBFlagA[8]<ii  & ADBBFlagA[7]<ii  & ADBBFlagA[6]<ii  & ADBBFlagA[5]<ii  & ADBBFlagA[4]<ii  & ADBBFlagA[3]<ii;
            SelGoodEventAD[2][0][1] &= ADBBFlagC[9]<ii  &  ADBBFlagC[8]<ii  & ADBBFlagC[7]<ii  & ADBBFlagC[6]<ii  & ADBBFlagC[5]<ii  & ADBBFlagC[4]<ii  & ADBBFlagC[3]<ii;
            // Bunch range 3-9, V0 : V0A, Flag : BB
            SelGoodEventAD[2][1][1] = ADBBFlagA[11]<ii & ADBBFlagA[12]<ii & ADBBFlagA[13]<ii & ADBBFlagA[14]<ii & ADBBFlagA[15]<ii & ADBBFlagA[16]<ii & ADBBFlagA[17]<ii;
            SelGoodEventAD[2][1][1] &= ADBBFlagA[9]<ii  &  ADBBFlagA[8]<ii  & ADBBFlagA[7]<ii  & ADBBFlagA[6]<ii  & ADBBFlagA[5]<ii  & ADBBFlagA[4]<ii  & ADBBFlagA[3]<ii;
            // Bunch range 3-9, V0 : V0C, Flag : BB
            SelGoodEventAD[2][2][1] = ADBBFlagC[11]<ii & ADBBFlagC[12]<ii & ADBBFlagC[13]<ii & ADBBFlagC[14]<ii & ADBBFlagC[15]<ii & ADBBFlagC[16]<ii & ADBBFlagC[17]<ii;
            SelGoodEventAD[2][2][1] &= ADBBFlagC[9]<ii  &  ADBBFlagC[8]<ii  & ADBBFlagC[7]<ii  & ADBBFlagC[6]<ii  & ADBBFlagC[5]<ii  & ADBBFlagC[4]<ii  & ADBBFlagC[3]<ii;
            //Flag variation : ADBG // Bunch range 3-9, V0 : V0AC, Flag : ADBG
            SelGoodEventAD[2][0][2] = ADBGFlagA[11]<ii & ADBGFlagA[12]<ii & ADBGFlagA[13]<ii & ADBGFlagA[14]<ii & ADBGFlagA[15]<ii & ADBGFlagA[16]<ii & ADBGFlagA[17]<ii;
            SelGoodEventAD[2][0][2] &= ADBGFlagC[11]<ii & ADBGFlagC[12]<ii & ADBGFlagC[13]<ii & ADBGFlagC[14]<ii & ADBGFlagC[15]<ii & ADBGFlagC[16]<ii & ADBGFlagC[17]<ii;
            SelGoodEventAD[2][0][2] &= ADBGFlagA[9]<ii  &  ADBGFlagA[8]<ii  & ADBGFlagA[7]<ii  & ADBGFlagA[6]<ii  & ADBGFlagA[5]<ii  & ADBGFlagA[4]<ii  & ADBGFlagA[3]<ii;
            SelGoodEventAD[2][0][2] &= ADBGFlagC[9]<ii  &  ADBGFlagC[8]<ii  & ADBGFlagC[7]<ii  & ADBGFlagC[6]<ii  & ADBGFlagC[5]<ii  & ADBGFlagC[4]<ii  & ADBGFlagC[3]<ii;
            //V0 variation 1 // Bunch range 3-9, V0 : V0A, Flag : ADBG
            SelGoodEventAD[2][1][2] = ADBGFlagA[11]<ii & ADBGFlagA[12]<ii & ADBGFlagA[13]<ii & ADBGFlagA[14]<ii & ADBGFlagA[15]<ii & ADBGFlagA[16]<ii & ADBGFlagA[17]<ii;
            SelGoodEventAD[2][1][2] &= ADBGFlagA[9]<ii  &  ADBGFlagA[8]<ii  & ADBGFlagA[7]<ii  & ADBGFlagA[6]<ii  & ADBGFlagA[5]<ii  & ADBGFlagA[4]<ii  & ADBGFlagA[3]<ii;
            //V0 variation 2 // Bunch range 3-9, V0 : V0C, Flag : ADBG
            SelGoodEventAD[2][2][2] = ADBGFlagC[11]<ii & ADBGFlagC[12]<ii & ADBGFlagC[13]<ii & ADBGFlagC[14]<ii & ADBGFlagC[15]<ii & ADBGFlagC[16]<ii & ADBGFlagC[17]<ii;
            SelGoodEventAD[2][2][2] &= ADBGFlagC[9]<ii  &  ADBGFlagC[8]<ii  & ADBGFlagC[7]<ii  & ADBGFlagC[6]<ii  & ADBGFlagC[5]<ii  & ADBGFlagC[4]<ii  & ADBGFlagC[3]<ii;
            
            
            double adcheck1=0;
            double adcheck2=0;
            double adcheck3=0;
            if(SelGoodEventAD[2][2][0]) adcheck1++;
            if(SelGoodEventAD[2][2][1]) adcheck2++;
            if(SelGoodEventAD[2][2][2]) adcheck3++;
            
            
            cout<< "AD ===3-9, v0c , bb+bg = " <<SelGoodEventAD[2][2][0]<<",  num 1 = "<<adcheck1<<endl;
            cout<< "AD ===3-9, v0c , bb = " <<SelGoodEventAD[2][2][1]<<",  num 2 = "<<adcheck2<<endl;
            cout<< "AD ===3-9, v0c , bg = " <<SelGoodEventAD[2][2][2]<<",  num 3 = "<<adcheck3<<endl;
            
          
            for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                    for(int k=0; k<3; k++){
                        
                        cout<< "AD i  = " <<i<<", AD j   = "<<j<< ",   AD k  = " <<k<<endl;

                        
                        if(SelGoodEventAD[i][j][k]) {
                            ((TH1F*)fList->FindObject(Form("hADDenomPurityBC%d_AD%d_Flag%d",i,j,k)))->Fill(ii-1);
                        }
                        
                        
                        if(!bgID) {
                               ((TH1F*)fList->FindObject(Form("hDenomAD%d_AD%d_Flag%d",i,j,k)))->Fill(fad0a-fad0c, fad0a+fad0c);

                            if(SelGoodEventAD[i][j][k]){
                                ((TH1F*)fList->FindObject(Form("hADNumEffPurityBC%d_AD%d_Flag%d",i,j,k)))->Fill(ii-1);
                                ((TH1F*)fList->FindObject(Form("hADNumTrkVsClsSPD%d_AD%d_Flag%d",i,j,k)))->Fill(fSpdT, fSpdC1+fSpdC2);
                                ((TH1F*)fList->FindObject(Form("hNumAD%d_AD%d_Flag%d",i,j,k)))->Fill(fad0a-fad0c, fad0a+fad0c);
                            }
                        }
                        
                        if(bgID){
                        if(!SelGoodEventAD[i][j][k]){
                                ((TH1F*)fList->FindObject(Form("hADNumRejecEffBC%d_AD%d_Flag%d",i,j,k)))->Fill(ii-1);
                            }
                        }
                        
                    
                    }
                }
            } // end of fill histograms
            
           
            
        } // end of AD flag loop

    } // end of events in trigger loop
    
    if(fUseTree==kTRUE)fTreeTrack->Fill();
    PostData(1, fList);
    if(fUseTree==kTRUE)PostData(2, fTreeTrack);
}


//________________________________________________________________________
void AliAnalysisBGMonitorQA::Terminate(Option_t *)
{
 //   fList = dynamic_cast<TList*> (GetOutputData(1));
 //   if(!fList)    Printf("ERROR: fList is not available");
}



//______________________________________________________________________ Modified cut function(blim)
Bool_t IsItBGSPDClusterVsTracklet(AliVEvent *event,float slope)
{
    Int_t nClustersLayer0 = event->GetNumberOfITSClusters(0);
    Int_t nClustersLayer1 = event->GetNumberOfITSClusters(1);
    Int_t nTracklets      = event->GetMultiplicity()->GetNumberOfTracklets();
    if (nClustersLayer0 + nClustersLayer1 > 65. + nTracklets*slope) return kTRUE;
    return kFALSE;
}

                