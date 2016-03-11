// Analysis Task for the Quality Assurence of Beam Gas Monitoring
//
// This code will check the each event for several parameters,
// and check if it is Background or Signal with below function.
// after that,
//
// Authors
// Alexander Borissov <aborisso@mail.cern.ch>
// Bong-Hwi Lim <bong-hwi.lim@cern.ch>
// Jihye Song <Jihye.Song@cern.ch>
//
// If you have any comment or question of this code,
// Please send a mail to Bong-Hwi or Jihye
//
// Last update: 2016.03.02 (blim)
//
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

Bool_t IsItBGSPDClusterVsTracklet(AliVEvent *event); // add function info and initial condition (blim)
Bool_t IsItBGSPDClusterVsTracklet2(AliVEvent *event); // add function to check bg in small tracklet 2015.09.14. (blim)

Bool_t SelGoodEvent[3][3][3];
Bool_t SelGoodEventAD[3][3][3];
Int_t BGFlagA[20];
Int_t BGFlagC[20];
Int_t BBFlagA[20];
Int_t BBFlagC[20];
Int_t ADBGFlagA[20];
Int_t ADBGFlagC[20];
Int_t ADBBFlagA[20];
Int_t ADBBFlagC[20];
Int_t bunchinputarray[7] = {201};  // the output file which we interested in. 2015.08.20. (blim)
//________________________________________________________________________
AliAnalysisBGMonitorQA::AliAnalysisBGMonitorQA(const char *name) :
AliAnalysisTaskSE(name),
fESD(0x0),
fESDfriend(0x0),
fTreeTrack(0),
fTreeTrack2(0),
fList(0),
fList2(0), //add new List for both result 2015.08.20. (blim)
fUseTree(kFALSE)

{
    
    // Constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class()); //add new line for both result 2015.08.20. (blim)
    DefineOutput(0, TTree::Class()); //add new line for Tree result 2016.01.29. (blim)
    
    if(fUseTree==kTRUE) DefineOutput(3, TTree::Class());
    
    
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
    //  fTreeTrack->Branch("t0PU",&t0PileUp,"t0PU/I"); // pile-up from T0
    fTreeTrack->Branch("spdPU",&spdPileUp,"spdPU/I"); // pile-up from SPD
    fTreeTrack->Branch("spdPUOOB",&spdPileUpOutOfBunch,"spdPUOOB/I"); // out-of-bunch pile-up from SPD
    fTreeTrack->Branch("BGFlagA",&BGFlagA,"BGFlagA[nbunch]/I"); // V0A BG flag for PF protection
    fTreeTrack->Branch("BGFlagC",&BGFlagC,"BGFlagC[nbunch]/I"); // V0C BG flag for PF protection
    fTreeTrack->Branch("BBFlagA",&BBFlagA,"BBFlagA[nbunch]/I"); // V0A BG flag for PF protection
    fTreeTrack->Branch("BBFlagC",&BBFlagC,"BBFlagC[nbunch]/I"); // V0C BG flag for PF protection
    if(fUseTree==kTRUE)PostData(3, fTreeTrack);
    
    if(fList != NULL){
        delete fList;
        fList = NULL;
    }
    if(fList2 != NULL){
        delete fList2;
        fList2 = NULL;
    }
    if(fList == NULL){
        fList = new TList();
        fList->SetOwner(kTRUE);
    }
    if(fList2 == NULL){
        fList2 = new TList();
        fList2->SetOwner(kTRUE); //add new List for both result 2015.08.20. (blim)
    }
    
    fTreeTrack2 = new TTree("TreeTrack","Track Properties2");
    fTreeTrack2->Branch("runNumber",&runNumber,"runNumber/I"); //run number
    PostData(0, fTreeTrack2);
    
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
    
    //______________________________add new List for both result 2015.08.20. (blim)
    TH1F *hNumEffPurityBC_HM[3][3][3];
    TH1F *hDenomEffBC_HM[3][3][3];
    TH1F *hDenomPurityBC_HM[3][3][3];
    TH1F *hDenomRejecEffBC_HM[3][3][3];
    TH1F *hNumRejecEffBC_HM[3][3][3];
    TH1F *hSPDNumBC_HM[3][3][3];
    TH1F *hSPDDenomBC_HM[3][3][3];
    TH2F *hNumTrkVsClsSPID_HM[3][3][3];
    TH2F *hDenomTrkVsClsSPID_HM[3][3][3];
    TH2F *hNumV0_HM[3][3][3];
    TH2F *hDenomV0_HM[3][3][3];
    
    TH1F *hADNumEffPurityBC_HM[3][3][3];
    TH1F *hADDenomPurityBC_HM[3][3][3];
    TH1F *hADNumRejecEffBC_HM[3][3][3];
    TH2F *hADNumTrkVsClsSPD_HM[3][3][3];
    TH2F *hNumAD_HM[3][3][3];
    TH2F *hDenomAD_HM[3][3][3];
    //______________________________
    
    
    
    TH1F *runNumber_hist;
    runNumber_hist = new TH1F("runNumber_hist","runNum", 1, 0, 1);
    fList->Add(runNumber_hist);
    TH1F *runNumber_hist_HM;
    runNumber_hist_HM = new TH1F("runNumber_hist_HM","runNum", 1, 0, 1);
    fList2->Add(runNumber_hist_HM);
    //______________________________
    
    TH2F *hTotalTrkVsClsSPID = new TH2F("hTotalTrkVsClsSPID","; Spd : total",140,0,140,500,0,500);
    hTotalTrkVsClsSPID->GetXaxis()->SetTitle("Tracklet");
    hTotalTrkVsClsSPID->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList->Add(hTotalTrkVsClsSPID);
    
    TH2F *hTotalTrkVsClsSPID_HM = new TH2F("hTotalTrkVsClsSPID_HM","; Spd : total",140,0,140,500,0,500);
    hTotalTrkVsClsSPID_HM->GetXaxis()->SetTitle("Tracklet");
    hTotalTrkVsClsSPID_HM->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList2->Add(hTotalTrkVsClsSPID_HM); //add new List for both result 2015.08.20. (blim)
    
    TH2F *hTotalV0 = new TH2F("hTotalV0","; V0 : total",600,-300,300,2000,-1000,1000);
    hTotalV0->GetXaxis()->SetTitle("V0A-V0C");
    hTotalV0->GetYaxis()->SetTitle("V0A+V0C");
    fList->Add(hTotalV0);
    
    TH2F *hTotalV0_HM = new TH2F("hTotalV0_HM","; V0 : total",600,-300,300,2000,-1000,1000);
    hTotalV0_HM->GetXaxis()->SetTitle("V0A-V0C");
    hTotalV0_HM->GetYaxis()->SetTitle("V0A+V0C");
    fList2->Add(hTotalV0_HM); //add new List for both result 2015.08.20. (blim)
    
    TH2F *hTotalAD = new TH2F("hTotalAD","; AD : total",400,-200,200,2000,-1000,1000);
    hTotalAD->GetXaxis()->SetTitle("ADA-ADC");
    hTotalAD->GetYaxis()->SetTitle("ADA+ADC");
    fList->Add(hTotalAD);
    
    TH2F *hTotalAD_HM = new TH2F("hTotalAD_HM","; AD : total",400,-200,200,2000,-1000,1000);
    hTotalAD_HM->GetXaxis()->SetTitle("ADA-ADC");
    hTotalAD_HM->GetYaxis()->SetTitle("ADA+ADC");
    fList2->Add(hTotalAD_HM); //add new List for both result 2015.08.20. (blim)
    
    //histogram for event list(blim)
    TH1F *hNumEvents  = new TH1F("hNumEvents","total event",10,0,10);
    fList->Add(hNumEvents);
    
    //histogram for modified Cut(blim)
    
    TH2F *hTrkVsClsSPIDSlopeM = new TH2F("hTrkVsClsSPIDSlopeM","; Spd : total",140,0,140,500,0,500);
    hTrkVsClsSPIDSlopeM->GetXaxis()->SetTitle("Tracklet");
    hTrkVsClsSPIDSlopeM->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList->Add(hTrkVsClsSPIDSlopeM);
    
    TH2F *hTrkVsClsSPIDSlopeM_HM = new TH2F("hTrkVsClsSPIDSlopeM_HM","; Spd : total",140,0,140,500,0,500);
    hTrkVsClsSPIDSlopeM_HM->GetXaxis()->SetTitle("Tracklet");
    hTrkVsClsSPIDSlopeM_HM->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList2->Add(hTrkVsClsSPIDSlopeM_HM); //add new List for both result 2015.08.20. (blim)

    
    TH2F *hTrkVsClsSPIDSlopeM2 = new TH2F("hTrkVsClsSPIDSlopeM2","; Spd : total",140,0,140,500,0,500);
    hTrkVsClsSPIDSlopeM2->GetXaxis()->SetTitle("Tracklet");
    hTrkVsClsSPIDSlopeM2->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList->Add(hTrkVsClsSPIDSlopeM2);
    
    TH2F *hTrkVsClsSPIDSlopeM_HM2 = new TH2F("hTrkVsClsSPIDSlopeM_HM2","; Spd : total",140,0,140,500,0,500);
    hTrkVsClsSPIDSlopeM_HM2->GetXaxis()->SetTitle("Tracklet");
    hTrkVsClsSPIDSlopeM_HM2->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList2->Add(hTrkVsClsSPIDSlopeM_HM2); //add new List for both result 2015.11.09. (blim)
    
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            for(int k=0; k<3; k++){
                Int_t check000 = i*100+j*10+k;
                for ( int l=0; l<8;l++){
                    if(check000==bunchinputarray[l]){
                        
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
                        
                        
                        //_______________________________________add new List for both result 2015.08.20. (blim)
                        
                        hNumEffPurityBC_HM[i][j][k] = new TH1F(Form("hNumEffPurityBC_HM%d_V0%d_Flag%d",i,j,k),"; #V0flags in PF", 35, 0, 35);
                        hDenomEffBC_HM[i][j][k] = new TH1F(Form("hDenomEffBC_HM%d_V0%d_Flag%d",i,j,k),"; #V0flags in PF", 35, 0, 35);
                        hDenomPurityBC_HM[i][j][k] = new TH1F(Form("hDenomPurityBC_HM%d_V0%d_Flag%d",i,j,k),"; #V0flags in PF", 35, 0, 35);
                        hDenomRejecEffBC_HM[i][j][k] = new TH1F(Form("hDenomRejecEffBC_HM%d_V0%d_Flag%d",i,j,k),"; #V0flags in PF", 35, 0, 35);
                        hNumRejecEffBC_HM[i][j][k] = new TH1F(Form("hNumRejecEffBC_HM%d_V0%d_Flag%d",i,j,k),"; #V0flags in PF", 35, 0, 35);
                        
                        hSPDNumBC_HM[i][j][k] = new TH1F(Form("hSPDNumBC_HM%d_V0%d_Flag%d",i,j,k),"; Spd tracklet", 200, 0, 200);
                        hSPDDenomBC_HM[i][j][k] = new TH1F(Form("hSPDDenomBC_HM%d_V0%d_Flag%d",i,j,k),"; Spd tracklet", 200, 0, 200);
                        
                        hNumTrkVsClsSPID_HM[i][j][k] = new TH2F(Form("hNumTrkVsClsSPID_HM%d_V0%d_Flag%d",i,j,k),"; Spd : !BGid & GoodEvent",140,0,140,500,0,500);
                        hNumTrkVsClsSPID_HM[i][j][k]->GetXaxis()->SetTitle("Tracklet");
                        hNumTrkVsClsSPID_HM[i][j][k]->GetYaxis()->SetTitle("Cluster (fspdC1)");
                        hDenomTrkVsClsSPID_HM[i][j][k]= new TH2F(Form("hDenomTrkVsClsSPID_HM%d_V0%d_Flag%d",i,j,k),"; Spd : !BGid",140,0,140,500,0,500);
                        hDenomTrkVsClsSPID_HM[i][j][k]->GetXaxis()->SetTitle("Tracklet");
                        hDenomTrkVsClsSPID_HM[i][j][k]->GetYaxis()->SetTitle("Cluster (fspdC1)");
                        
                        hNumV0_HM[i][j][k] = new TH2F(Form("hNumV0_HM%d_V0%d_Flag%d",i,j,k),"; V0 : !BGid & GoodEvent",600,-300,300,2000,-1000,1000);
                        hNumV0_HM[i][j][k]->GetXaxis()->SetTitle("V0A-V0C");
                        hNumV0_HM[i][j][k]->GetYaxis()->SetTitle("V0A+V0C");
                        hDenomV0_HM[i][j][k]= new TH2F(Form("hDenomV0_HM%d_V0%d_Flag%d",i,j,k),"; V0 : !BGid",600,-300,300,2000,-1000,1000);
                        hDenomV0_HM[i][j][k]->GetXaxis()->SetTitle("V0A-V0C");
                        hDenomV0_HM[i][j][k]->GetYaxis()->SetTitle("V0A+V0C");
                        
                        
                        hADNumEffPurityBC_HM[i][j][k] = new TH1F(Form("hADNumEffPurityBC_HM%d_AD%d_Flag%d",i,j,k),"; #ADflags in PF", 10, 0, 10);
                        hADDenomPurityBC_HM[i][j][k] = new TH1F(Form("hADDenomPurityBC_HM%d_AD%d_Flag%d",i,j,k),"; #ADflags in PF", 10, 0, 10);
                        hADNumRejecEffBC_HM[i][j][k] = new TH1F(Form("hADNumRejecEffBC_HM%d_AD%d_Flag%d",i,j,k),"; #ADflags in PF", 10, 0, 10);
                        hADNumTrkVsClsSPD_HM[i][j][k] = new TH2F(Form("hADNumTrkVsClsSPD_HM%d_AD%d_Flag%d",i,j,k),"; Spd : !BGid & GoodEvent",140,0,140,500,0,500);
                        hADNumTrkVsClsSPD_HM[i][j][k]->GetXaxis()->SetTitle("Tracklet");
                        hADNumTrkVsClsSPD_HM[i][j][k]->GetYaxis()->SetTitle("Cluster (fspdC1)");
                        
                        hNumAD_HM[i][j][k] = new TH2F(Form("hNumAD_HM%d_AD%d_Flag%d",i,j,k),"; AD : !BGid & GoodEvent",400,-200,200,2000,-1000,1000);
                        hNumAD_HM[i][j][k]->GetXaxis()->SetTitle("ADA-ADC");
                        hNumAD_HM[i][j][k]->GetYaxis()->SetTitle("ADA+ADC");
                        hDenomAD_HM[i][j][k]= new TH2F(Form("hDenomAD_HM%d_AD%d_Flag%d",i,j,k),"; AD : !BGid",400,-200,200,2000,-1000,1000);
                        hDenomAD_HM[i][j][k]->GetXaxis()->SetTitle("ADA-ADC");
                        hDenomAD_HM[i][j][k]->GetYaxis()->SetTitle("ADA+ADC");
                        
                        fList2->Add(hNumEffPurityBC_HM[i][j][k]);
                        fList2->Add(hDenomEffBC_HM[i][j][k]);
                        fList2->Add(hDenomPurityBC_HM[i][j][k]);
                        fList2->Add(hDenomRejecEffBC_HM[i][j][k]);
                        fList2->Add(hNumRejecEffBC_HM[i][j][k]);
                        fList2->Add(hSPDNumBC_HM[i][j][k]);
                        fList2->Add(hSPDDenomBC_HM[i][j][k]);
                        fList2->Add(hNumTrkVsClsSPID_HM[i][j][k]);
                        fList2->Add(hDenomTrkVsClsSPID_HM[i][j][k]);
                        fList2->Add(hNumV0_HM[i][j][k]);
                        fList2->Add(hDenomV0_HM[i][j][k]);
                        
                        fList2->Add(hADNumEffPurityBC_HM[i][j][k]);
                        fList2->Add(hADDenomPurityBC_HM[i][j][k]);
                        fList2->Add(hADNumRejecEffBC_HM[i][j][k]);
                        fList2->Add(hADNumTrkVsClsSPD_HM[i][j][k]);
                        fList2->Add(hNumAD_HM[i][j][k]);
                        fList2->Add(hDenomAD_HM[i][j][k]);
                        //_______________________________________
                        
                        
                    }
                }
            }
        }
        
    }
    
    
    
    PostData(1, fList);
    PostData(2, fList2);
    
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

    ((TH1F*)fList->FindObject("runNumber_hist"))->SetBinContent(1,runNumber);
    ((TH1F*)fList2->FindObject("runNumber_hist_HM"))->SetBinContent(1,runNumber);
    
    UInt_t timeGDC=fESD->GetTimeStamp();
    ftime=timeGDC;
    Int_t timeStampBX = fESD->GetBunchCrossNumber();
    fbx=timeStampBX;
    ntr = 10;
    nbunch = 21;
   // ofstream ftxt;
    
    static AliTriggerAnalysis * triggerAnalysis = new AliTriggerAnalysis();
    
    V0A = 0;
    V0C = 0;
    V0ABG = 0;
    V0CBG = 0;
    bgID = 0;
    
    // additional value initialize (blim)
    bgID2=0;
    
    VBA = 0;
    VBC = 0;
    VGA = 0;
    VTX = 0;
    fastORHW = 0;
    SPD1 = 0;
    SPD2 = 0;
    SPDHw1 = 0;
    SPDHw2 = 0;
    
    for (Int_t i=0; i<3; i++) {
        for (Int_t j=0; j<3; j++) {
            for (Int_t k=0; k<3; k++) {
                SelGoodEvent[i][j][k] = 0;
                //cout << "SelGoodEvent " << i << " " << j << " " << k << " is" << SelGoodEvent[i][j][k] << endl;
                //cout << BGFlagA[i] << " " << BGFlagC[i] << " " << BBFlagA[i] << " " << BBFlagC[i] << " " << endl;
            }
        }
    }
    // initialize 2015.08.17.(blim)
    
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
    //    bgID = utils->IsSPDClusterVsTrackletBG(fESD);
    
    // modified slope cut. the function is in below of this source(blim)
    bgID = IsItBGSPDClusterVsTracklet(fESD); // original modified function
    bgID2 = IsItBGSPDClusterVsTracklet2(fESD); // modified modified function
    
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
    /*
     memset(BGFlagA, 0, sizeof(Float_t)*nbunch);
     memset(BBFlagA, 0, sizeof(Float_t)*nbunch);
     memset(BGFlagC, 0, sizeof(Float_t)*nbunch);
     memset(BBFlagC, 0, sizeof(Float_t)*nbunch);
     
     memset(ADBGFlagA, 0, sizeof(Float_t)*nbunch);
     memset(ADBBFlagA, 0, sizeof(Float_t)*nbunch);
     memset(ADBGFlagC, 0, sizeof(Float_t)*nbunch);
     memset(ADBBFlagC, 0, sizeof(Float_t)*nbunch);
     */
    for (Int_t i=0; i<20; i++) {
        BBFlagA[i]=0;
        BGFlagA[i]=0;
        BBFlagC[i]=0;
        BGFlagC[i]=0;
        ADBBFlagA[i]=0;
        ADBBFlagC[i]=0;
        ADBGFlagA[i]=0;
        ADBGFlagC[i]=0;
    } // initialize 2015.08.17.(blim)
    
    
    
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
        Printf("No esdV0friend available : test Jihye");
        return;
    }
    
    /*
     for (Int_t i=0; i<20; i++) {
     cout << "BBflagA[" << i << "] = " << BBFlagA[i] << endl;
     cout << "BGflagA[" << i << "] = " << BGFlagA[i] << endl;
     cout << "BBflagC[" << i << "] = " << BBFlagC[i] << endl;
     cout << "BGflagC[" << i << "] = " << BGFlagC[i] << endl;
     cout << "ADBBflagA[" << i << "] = " << ADBBFlagA[i] << endl;
     cout << "ADBGflagA[" << i << "] = " << ADBGFlagA[i] << endl;
     cout << "ADBBflagC[" << i << "] = " << ADBBFlagC[i] << endl;
     cout << "ADBGflagC[" << i << "] = " << ADBGFlagC[i] << endl;
     } // show Flags 2015.08.17.(blim)
     */
    
    
    ntracks = fESD->GetNumberOfTracks(); // number of tracks (no quality cuts)
    
    //--- Trigger classes --//
    memset(ftrigger, 0, sizeof(Float_t)*ntr);
    
    //Minimum Bias
    if(fESD->IsTriggerClassFired("CINT7-B-NOPF-ALLNOTRD") || fESD->IsTriggerClassFired("CINT7-S-NOPF-ALLNOTRD") || fESD->IsTriggerClassFired("CINT1-B-NOPF-ALLNOTRD") || fESD->IsTriggerClassFired("CINT1-S-NOPF-ALLNOTRD") || fESD->IsTriggerClassFired("CINT7-A-NOPF-CENT") || fESD->IsTriggerClassFired("CINT7-B-NOPF-CENT") || fESD->IsTriggerClassFired("CINT7-C-NOPF-CENT") || fESD->IsTriggerClassFired("CINT7-E-NOPF-CENT")) ftrigger[0] = 1;
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
    
    if(fESD->IsTriggerClassFired("CTEST52-B-NOPF-ALLNOTRD")|| fESD->IsTriggerClassFired("CTEST52-B-NOPF-ALLNOTRD")|| fESD->IsTriggerClassFired("CTEST52-C-NOPF-ALLNOTRD")|| fESD->IsTriggerClassFired("CTEST52-E-NOPF-ALLNOTRD") || fESD->IsTriggerClassFired("CTEST52-B-NOPF-CENT") || fESD->IsTriggerClassFired("CTEST52-C-NOPF-CENT") || fESD->IsTriggerClassFired("CTEST52-A-NOPF-CENT") || fESD->IsTriggerClassFired("CTEST52-E-NOPF-CENT") || fESD->IsTriggerClassFired("CTEST52-B-NOPF-CENTNOTRD") || fESD->IsTriggerClassFired("CTEST52-A-NOPF-CENTNOTRD") || fESD->IsTriggerClassFired("CTEST52-C-NOPF-CENTNOTRD") || fESD->IsTriggerClassFired("CTEST52-E-NOPF-CENTNOTRD") || fESD->IsTriggerClassFired("CVHMV0M-A-NOPF-CENT") || fESD->IsTriggerClassFired("CVHMV0M-B-NOPF-CENT") || fESD->IsTriggerClassFired("CVHMV0M-C-NOPF-CENT") || fESD->IsTriggerClassFired("CVHMV0M-E-NOPF-CENT") || fESD->IsTriggerClassFired("CVHMV0M-B-spdl-CENT") || fESD->IsTriggerClassFired("CVHMSH2-A-NOPF-CENT") || fESD->IsTriggerClassFired("CVHMSH2-B-NOPF-CENT") || fESD->IsTriggerClassFired("CVHMSH2-C-NOPF-CENT") || fESD->IsTriggerClassFired("CVHMSH2-E-NOPF-CENT") || fESD->IsTriggerClassFired("CVHMV0M-A-NOPF-CENTNOTRD") || fESD->IsTriggerClassFired("CVHMV0M-B-NOPF-CENTNOTRD") || fESD->IsTriggerClassFired("CVHMV0M-C-NOPF-CENTNOTRD") || fESD->IsTriggerClassFired("CVHMV0M-E-NOPF-CENTNOTRD") || fESD->IsTriggerClassFired("CVHMSH2-A-NOPF-CENTNOTRD") || fESD->IsTriggerClassFired("CVHMSH2-B-NOPF-CENTNOTRD") || fESD->IsTriggerClassFired("CVHMSH2-C-NOPF-CENTNOTRD") || fESD->IsTriggerClassFired("CVHMSH2-E-NOPF-CENTNOTRD")) ftrigger[9] = 1; // for after LHC15h, run 263000 updated by blim(20160223)
    
    // count total event number (blim)
    if(ftrigger[0]==1)((TH1F*)fList->FindObject("hNumEvents"))->Fill(1);
    if(ftrigger[9]==1)((TH1F*)fList->FindObject("hNumEvents"))->Fill(3);
    
    
    int nAfterBunch = 3;
    int nV0 = 3;
    int nFlag = 3;
    
    // Bool_t goodEvent;
    //    static Bool_t SelGoodEvent[nAfterBunch][nV0][nFlag];
    //    static Bool_t SelGoodEventAD[nAfterBunch][nV0][nFlag];
    
    
    
    if(ftrigger[0]) {  // trigger class for MB
        
        ((TH1F*)fList->FindObject("hTotalTrkVsClsSPID"))->Fill(fSpdT, fSpdC1+fSpdC2);
        ((TH1F*)fList->FindObject("hTotalV0"))->Fill(fv0a-fv0c, fv0a+fv0c);
        ((TH1F*)fList->FindObject("hTotalAD"))->Fill(fad0a-fad0c, fad0a+fad0c);
        
        
        // Modified Cut result, added by blim
        if (!bgID) {
            ((TH1F*)fList->FindObject("hTrkVsClsSPIDSlopeM"))->Fill(fSpdT, fSpdC1+fSpdC2); // bgID
        }
        if (!bgID2) {
            ((TH1F*)fList->FindObject("hTrkVsClsSPIDSlopeM2"))->Fill(fSpdT, fSpdC1+fSpdC2); // bgID2
        }
        
        for(Int_t ii=1; ii<33; ii++){
            
            
            //___________ modified method 2015.08.12.(blim)
            for (Int_t windowvar=0; windowvar<3; windowvar++) {
                for (Int_t flagvar=0; flagvar<3; flagvar++) {
                    for (Int_t v0var=0; v0var<3; v0var++) {
                        SelectGoodEventWithV0Variation(windowvar,v0var,flagvar,ii);
                    }
                }
            }
            
            //___________
            
            double check1=0;
            double check2=0;
            double check3=0;
            if(SelGoodEvent[2][0][1]) check1++;
            if(SelGoodEvent[2][2][1]) check2++;
            if(SelGoodEvent[2][2][2]) check3++;
            
            /*
             cout<< "3-9, v0a+c , bb = " <<SelGoodEvent[2][0][1]<<",  num 1 = "<<check1<<endl;
             cout<< "3-9, v0c , bb = " <<SelGoodEvent[2][2][1]<<",  num 2 = "<<check2<<endl;
             cout<< "3-9, v0c , bg = " <<SelGoodEvent[2][2][2]<<",  num 3 = "<<check3<<endl;
             */
            
            for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                    for(int k=0; k<3; k++){
                        Int_t check000 = i*100+j*10+k;
                        for ( int l=0; l<8;l++){
                            if(check000==bunchinputarray[l]){
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
                    }
                }
            } // end of fill histograms
        } // end of V0 flag loop
        /*
         for (Int_t i=0; i<20; i++) {
         if(BBFlagA[i]) cout << "BBflagA[" << i << "] in main = " << BBFlagA[i] << endl;
         if(BBFlagC[i]) cout << "BBflagC[" << i << "] in main = " << BBFlagC[i] << endl;
         if(BGFlagA[i]) cout << "BGflagA[" << i << "] in main = " << BGFlagA[i] << endl;
         if(BGFlagC[i]) cout << "BGflagC[" << i << "] in main = " << BGFlagC[i] << endl;
         if(ADBBFlagA[i]) cout << "ADBBflagA[" << i << "] in main = " << ADBBFlagA[i] << endl;
         if(ADBBFlagC[i]) cout << "ADBBflagC[" << i << "] in main = " << ADBBFlagC[i] << endl;
         if(ADBGFlagA[i]) cout << "ADBGflagA[" << i << "] in main = " << ADBGFlagA[i] << endl;
         if(ADBGFlagC[i]) cout << "ADBGflagC[" << i << "] in main = " << ADBGFlagC[i] << endl;
         } // show Flags 2015.08.17.(blim)
         */
        for(Int_t ii=1; ii<9; ii++){
            
            //___________ modified method 2015.08.12.(blim)
            for (Int_t windowvar=0; windowvar<3; windowvar++) {
                for (Int_t flagvar=0; flagvar<3; flagvar++) {
                    for (Int_t v0var=0; v0var<3; v0var++) {
                        SelectADGoodEventWithV0Variation(windowvar,v0var,flagvar,ii);
                    }
                }
            }
            //___________
            
            
            double adcheck1=0;
            double adcheck2=0;
            double adcheck3=0;
            if(SelGoodEventAD[2][2][0]) adcheck1++;
            if(SelGoodEventAD[2][2][1]) adcheck2++;
            if(SelGoodEventAD[2][2][2]) adcheck3++;
            
            /*
             cout<< "AD ===3-9, v0c , bb+bg = " <<SelGoodEventAD[2][2][0]<<",  num 1 = "<<adcheck1<<endl;
             cout<< "AD ===3-9, v0c , bb = " <<SelGoodEventAD[2][2][1]<<",  num 2 = "<<adcheck2<<endl;
             cout<< "AD ===3-9, v0c , bg = " <<SelGoodEventAD[2][2][2]<<",  num 3 = "<<adcheck3<<endl;
             */
            
            for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                    for(int k=0; k<3; k++){
                        Int_t check000 = i*100+j*10+k;
                        for ( int l=0; l<8;l++){
                            if(check000==bunchinputarray[l]){
                                
                                //cout<< "AD i  = " <<i<<", AD j   = "<<j<< ",   AD k  = " <<k<<endl;
                                
                                
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
                        
                    }
                }
            } // end of fill histograms
            
        } // end of AD flag loop
        
    } // end of events in trigger loop
    
    for (Int_t i=0; i<3; i++) {
        for (Int_t j=0; j<3; j++) {
            for (Int_t k=0; k<3; k++) {
                SelGoodEvent[i][j][k] = 0;
                //cout << "SelGoodEvent " << i << " " << j << " " << k << " is" << SelGoodEvent[i][j][k] << endl;
                //cout << BGFlagA[i] << " " << BGFlagC[i] << " " << BBFlagA[i] << " " << BBFlagC[i] << " " << endl;
            }
        }
    }
    // initialize 2015.08.17.(blim)
    
    if(ftrigger[9]) {  // trigger class for HM // add new List for both result 2015.08.20. (blim)
        
        ((TH1F*)fList2->FindObject("hTotalTrkVsClsSPID_HM"))->Fill(fSpdT, fSpdC1+fSpdC2);
        ((TH1F*)fList2->FindObject("hTotalV0_HM"))->Fill(fv0a-fv0c, fv0a+fv0c);
        ((TH1F*)fList2->FindObject("hTotalAD_HM"))->Fill(fad0a-fad0c, fad0a+fad0c);
        
        
        // Modified Cut result, added by blim
        if (!bgID) {
            ((TH1F*)fList2->FindObject("hTrkVsClsSPIDSlopeM_HM"))->Fill(fSpdT, fSpdC1+fSpdC2); // bgID3->slope3
        }
        if (!bgID2) {
            ((TH1F*)fList2->FindObject("hTrkVsClsSPIDSlopeM_HM2"))->Fill(fSpdT, fSpdC1+fSpdC2); // bgID2
        }
        
        
        for(Int_t ii=1; ii<33; ii++){
            
            
            //___________ modified method 2015.08.12.(blim)
            for (Int_t windowvar=0; windowvar<3; windowvar++) {
                for (Int_t flagvar=0; flagvar<3; flagvar++) {
                    for (Int_t v0var=0; v0var<3; v0var++) {
                        SelectGoodEventWithV0Variation(windowvar,v0var,flagvar,ii);
                    }
                }
            }
            
            //___________
            /*
             double check1=0;
             double check2=0;
             double check3=0;
             if(SelGoodEvent[2][0][1]) check1++;
             if(SelGoodEvent[2][2][1]) check2++;
             if(SelGoodEvent[2][2][2]) check3++;
             
             
             cout<< "3-9, v0a+c , bb = " <<SelGoodEvent[2][0][1]<<",  num 1 = "<<check1<<endl;
             cout<< "3-9, v0c , bb = " <<SelGoodEvent[2][2][1]<<",  num 2 = "<<check2<<endl;
             cout<< "3-9, v0c , bg = " <<SelGoodEvent[2][2][2]<<",  num 3 = "<<check3<<endl;
             
             */
            for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                    for(int k=0; k<3; k++){
                        Int_t check000 = i*100+j*10+k;
                        for ( int l=0; l<8;l++){
                            if(check000==bunchinputarray[l]){
                                if(SelGoodEvent[i][j][k]) {
                                    ((TH1F*)fList2->FindObject(Form("hDenomPurityBC_HM%d_V0%d_Flag%d",i,j,k)))->Fill(ii-1);
                                }
                                
                                if(!bgID) {
                                    ((TH1F*)fList2->FindObject(Form("hDenomEffBC_HM%d_V0%d_Flag%d",i,j,k)))->Fill(ii-1);
                                    ((TH1F*)fList2->FindObject(Form("hSPDDenomBC_HM%d_V0%d_Flag%d",i,j,k)))->Fill(fSpdT);
                                    ((TH1F*)fList2->FindObject(Form("hDenomTrkVsClsSPID_HM%d_V0%d_Flag%d",i,j,k)))->Fill(fSpdT, fSpdC1+fSpdC2);
                                    ((TH1F*)fList2->FindObject(Form("hDenomV0_HM%d_V0%d_Flag%d",i,j,k)))->Fill(fv0a-fv0c, fv0a+fv0c);
                                    
                                    
                                    if(SelGoodEvent[i][j][k]){
                                        ((TH1F*)fList2->FindObject(Form("hNumEffPurityBC_HM%d_V0%d_Flag%d",i,j,k)))->Fill(ii-1);
                                        ((TH1F*)fList2->FindObject(Form("hSPDNumBC_HM%d_V0%d_Flag%d",i,j,k)))->Fill(fSpdT);
                                        ((TH1F*)fList2->FindObject(Form("hNumTrkVsClsSPID_HM%d_V0%d_Flag%d",i,j,k)))->Fill(fSpdT, fSpdC1+fSpdC2);
                                        ((TH1F*)fList2->FindObject(Form("hNumV0_HM%d_V0%d_Flag%d",i,j,k)))->Fill(fv0a-fv0c, fv0a+fv0c);
                                        
                                    }
                                }
                                if(bgID){
                                    ((TH1F*)fList2->FindObject(Form("hDenomRejecEffBC_HM%d_V0%d_Flag%d",i,j,k)))->Fill(ii-1);
                                    if(!SelGoodEvent[i][j][k]){
                                        ((TH1F*)fList2->FindObject(Form("hNumRejecEffBC_HM%d_V0%d_Flag%d",i,j,k)))->Fill(ii-1);
                                    }
                                }
                            }
                        }
                    }
                }
            } // end of fill histograms
        } // end of V0 flag loop
        /*
         for (Int_t i=0; i<20; i++) {
         if(BBFlagA[i]) cout << "BBflagA[" << i << "] in main = " << BBFlagA[i] << endl;
         if(BBFlagC[i]) cout << "BBflagC[" << i << "] in main = " << BBFlagC[i] << endl;
         if(BGFlagA[i]) cout << "BGflagA[" << i << "] in main = " << BGFlagA[i] << endl;
         if(BGFlagC[i]) cout << "BGflagC[" << i << "] in main = " << BGFlagC[i] << endl;
         if(ADBBFlagA[i]) cout << "ADBBflagA[" << i << "] in main = " << ADBBFlagA[i] << endl;
         if(ADBBFlagC[i]) cout << "ADBBflagC[" << i << "] in main = " << ADBBFlagC[i] << endl;
         if(ADBGFlagA[i]) cout << "ADBGflagA[" << i << "] in main = " << ADBGFlagA[i] << endl;
         if(ADBGFlagC[i]) cout << "ADBGflagC[" << i << "] in main = " << ADBGFlagC[i] << endl;
         } // show Flags 2015.08.17.(blim)
         */
        for(Int_t ii=1; ii<9; ii++){
            
            //___________ modified method 2015.08.12.(blim)
            for (Int_t windowvar=0; windowvar<3; windowvar++) {
                for (Int_t flagvar=0; flagvar<3; flagvar++) {
                    for (Int_t v0var=0; v0var<3; v0var++) {
                        SelectADGoodEventWithV0Variation(windowvar,v0var,flagvar,ii);
                    }
                }
            }
            //___________
            
            /*
             double adcheck1=0;
             double adcheck2=0;
             double adcheck3=0;
             if(SelGoodEventAD[2][2][0]) adcheck1++;
             if(SelGoodEventAD[2][2][1]) adcheck2++;
             if(SelGoodEventAD[2][2][2]) adcheck3++;
             
             
             cout<< "AD ===3-9, v0c , bb+bg = " <<SelGoodEventAD[2][2][0]<<",  num 1 = "<<adcheck1<<endl;
             cout<< "AD ===3-9, v0c , bb = " <<SelGoodEventAD[2][2][1]<<",  num 2 = "<<adcheck2<<endl;
             cout<< "AD ===3-9, v0c , bg = " <<SelGoodEventAD[2][2][2]<<",  num 3 = "<<adcheck3<<endl;
             */
            
            for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                    for(int k=0; k<3; k++){
                        Int_t check000 = i*100+j*10+k;
                        for ( int l=0; l<8;l++){
                            if(check000==bunchinputarray[l]){
                                
                                //cout<< "AD i  = " <<i<<", AD j   = "<<j<< ",   AD k  = " <<k<<endl;
                                
                                
                                if(SelGoodEventAD[i][j][k]) {
                                    ((TH1F*)fList2->FindObject(Form("hADDenomPurityBC_HM%d_AD%d_Flag%d",i,j,k)))->Fill(ii-1);
                                }
                                
                                
                                if(!bgID) {
                                    ((TH1F*)fList2->FindObject(Form("hDenomAD_HM%d_AD%d_Flag%d",i,j,k)))->Fill(fad0a-fad0c, fad0a+fad0c);
                                    
                                    if(SelGoodEventAD[i][j][k]){
                                        ((TH1F*)fList2->FindObject(Form("hADNumEffPurityBC_HM%d_AD%d_Flag%d",i,j,k)))->Fill(ii-1);
                                        ((TH1F*)fList2->FindObject(Form("hADNumTrkVsClsSPD_HM%d_AD%d_Flag%d",i,j,k)))->Fill(fSpdT, fSpdC1+fSpdC2);
                                        ((TH1F*)fList2->FindObject(Form("hNumAD_HM%d_AD%d_Flag%d",i,j,k)))->Fill(fad0a-fad0c, fad0a+fad0c);
                                    }
                                }
                                
                                if(bgID){
                                    if(!SelGoodEventAD[i][j][k]){
                                        ((TH1F*)fList2->FindObject(Form("hADNumRejecEffBC_HM%d_AD%d_Flag%d",i,j,k)))->Fill(ii-1);
                                    }
                                }
                            }
                        }
                        
                    }
                }
            } // end of fill histograms
            
        } // end of AD flag loop
        
    } // end of events in trigger loop
    
    if(fUseTree==kTRUE)fTreeTrack->Fill();
    PostData(1, fList);
    PostData(2, fList2);
    if(fUseTree==kTRUE)PostData(3, fTreeTrack);
    fTreeTrack2->Fill();
    PostData(0, fTreeTrack2);
}


//________________________________________________________________________
void AliAnalysisBGMonitorQA::Terminate(Option_t *)
{
    //   fList = dynamic_cast<TList*> (GetOutputData(1));
    //   if(!fList)    Printf("ERROR: fList is not available");
}



//______________________________________________________________________ Modified cut function(blim)
Bool_t IsItBGSPDClusterVsTracklet(AliVEvent *event)
{
    /*
     Int_t nClustersLayer0 = event->GetNumberOfITSClusters(0);
     Int_t nClustersLayer1 = event->GetNumberOfITSClusters(1);
     Int_t nTracklets      = event->GetMultiplicity()->GetNumberOfTracklets();
     if (nClustersLayer0 + nClustersLayer1 > 65. + nTracklets*slope) return kTRUE;
     return kFALSE;
     */
    Double_t nClustersLayer0 = event->GetNumberOfITSClusters(0);
    Double_t nClustersLayer1 = event->GetNumberOfITSClusters(1);
    Double_t trk = event->GetMultiplicity()->GetNumberOfTracklets();
    Double_t cls = nClustersLayer0 + nClustersLayer1;
    Bool_t spdBg = kFALSE;
    
    if (trk < 1.5) {
        spdBg = kTRUE;
    }
    else if (trk >= 1.5 && trk < 26 && cls > 20.0 + (378-20)/(26-1.5)*(trk-1.5)) {
        spdBg = kTRUE;
    }
    else if (trk >= 26 && trk < 40 && cls > 378.0 + (505-378)/(40-26.)*(trk-26.)) {
        spdBg = kTRUE;
    }
    else if (trk >= 40 && cls > 505.0 + (1770-505)/(300-40.)*(trk-40.)) {
        spdBg = kTRUE;
    }
    
    return spdBg;
}
//______________________________________________________________________ Modified cut function(blim)
Bool_t IsItBGSPDClusterVsTracklet2(AliVEvent *event)
{
    /*
     Int_t nClustersLayer0 = event->GetNumberOfITSClusters(0);
     Int_t nClustersLayer1 = event->GetNumberOfITSClusters(1);
     Int_t nTracklets      = event->GetMultiplicity()->GetNumberOfTracklets();
     if (nClustersLayer0 + nClustersLayer1 > 65. + nTracklets*slope) return kTRUE;
     return kFALSE;
     */
    Double_t nClustersLayer0 = event->GetNumberOfITSClusters(0);
    Double_t nClustersLayer1 = event->GetNumberOfITSClusters(1);
    Double_t trk = event->GetMultiplicity()->GetNumberOfTracklets();
    Double_t cls = nClustersLayer0 + nClustersLayer1;
    Bool_t spdBg = kFALSE;
    
    if (trk < 5.944 && cls > 65 + 4.*trk) {
        spdBg = kTRUE;
    }
    else if (trk >= 5.944 && trk < 26 && cls > 20.0 + (378-20)/(26-1.5)*(trk-1.5)) {
        spdBg = kTRUE;
    }
    else if (trk >= 26 && trk < 40 && cls > 378.0 + (505-378)/(40-26.)*(trk-26.)) {
        spdBg = kTRUE;
    }
    else if (trk >= 40 && cls > 505.0 + (1770-505)/(300-40.)*(trk-40.)) {
        spdBg = kTRUE;
    }
    
    return spdBg;
}




//______________________

void AliAnalysisBGMonitorQA::SelectGoodEventWithV0Variation(Int_t bunchrange, Int_t v0variation ,Int_t flagvariation, Int_t ii){
    
    
    //Bunch range 0-9 and 11-21
    if (bunchrange == 0) {
        if (v0variation == 0 || v0variation == 1) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEvent[bunchrange][v0variation][flagvariation] = BGFlagA[11]<ii & BGFlagA[12]<ii & BGFlagA[13]<ii & BGFlagA[14]<ii & BGFlagA[15]<ii & BGFlagA[16]<ii & BGFlagA[17]<ii & BGFlagA[18]<ii & BGFlagA[19]<ii & BGFlagA[20]<ii & BGFlagA[21]<ii;//BG-A 11-21
        }
        if (v0variation == 0) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEvent[bunchrange][v0variation][flagvariation] &= BGFlagC[11]<ii & BGFlagC[12]<ii & BGFlagC[13]<ii & BGFlagC[14]<ii & BGFlagC[15]<ii & BGFlagC[16]<ii & BGFlagC[17]<ii & BGFlagC[18]<ii & BGFlagC[19]<ii & BGFlagC[20]<ii & BGFlagC[21]<ii; //BG-C 11-21
        }
        if (v0variation == 2) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEvent[bunchrange][v0variation][flagvariation] = BGFlagC[11]<ii & BGFlagC[12]<ii & BGFlagC[13]<ii & BGFlagC[14]<ii & BGFlagC[15]<ii & BGFlagC[16]<ii & BGFlagC[17]<ii & BGFlagC[18]<ii & BGFlagC[19]<ii & BGFlagC[20]<ii & BGFlagC[21]<ii; //BG-C 11-21 for the first loop
        }
        if (v0variation == 0 ||v0variation == 1) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEvent[bunchrange][v0variation][flagvariation] &= BGFlagA[9]<ii & BGFlagA[8]<ii & BGFlagA[7]<ii  & BGFlagA[6]<ii  & BGFlagA[5]<ii  & BGFlagA[4]<ii & BGFlagA[3]<ii & BGFlagA[2]<ii & BGFlagA[1]<ii & BGFlagA[0]<ii; //BG-A 0-9
        }
        if (v0variation == 0 ||v0variation == 2) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEvent[bunchrange][v0variation][flagvariation] &= BGFlagC[9]<ii & BGFlagC[8]<ii & BGFlagC[7]<ii  & BGFlagC[6]<ii  & BGFlagC[5]<ii  & BGFlagC[4]<ii & BGFlagC[3]<ii & BGFlagC[2]<ii & BGFlagC[1]<ii & BGFlagC[0]<ii; //BG-C 0-9
        }
        if (v0variation == 0 ||v0variation == 1) if (flagvariation == 0){
            SelGoodEvent[bunchrange][v0variation][flagvariation] &= BBFlagA[11]<ii & BBFlagA[12]<ii & BBFlagA[13]<ii & BBFlagA[14]<ii & BBFlagA[15]<ii & BBFlagA[16]<ii & BBFlagA[17]<ii & BBFlagA[18]<ii & BBFlagA[19]<ii & BBFlagA[20]<ii & BBFlagA[21]<ii; //BB-A 11-21
        }
        if (v0variation == 0 ||v0variation == 1) if (flagvariation == 1){
            SelGoodEvent[bunchrange][v0variation][flagvariation] = BBFlagA[11]<ii & BBFlagA[12]<ii & BBFlagA[13]<ii & BBFlagA[14]<ii & BBFlagA[15]<ii & BBFlagA[16]<ii & BBFlagA[17]<ii & BBFlagA[18]<ii & BBFlagA[19]<ii & BBFlagA[20]<ii & BBFlagA[21]<ii; //BB-A 11-21 for the first loop
        }
        if (v0variation == 0 ) if (flagvariation == 0 || flagvariation == 1){
            SelGoodEvent[bunchrange][v0variation][flagvariation] &= BBFlagC[11]<ii & BBFlagC[12]<ii & BBFlagC[13]<ii & BBFlagC[14]<ii & BBFlagC[15]<ii & BBFlagC[16]<ii & BBFlagC[17]<ii & BBFlagC[18]<ii & BBFlagC[19]<ii & BBFlagC[20]<ii & BBFlagC[21]<ii; //BB-C 11-21
        }
        if (v0variation == 2) if (flagvariation == 0){
            SelGoodEvent[bunchrange][v0variation][flagvariation] &= BBFlagC[11]<ii & BBFlagC[12]<ii & BBFlagC[13]<ii & BBFlagC[14]<ii & BBFlagC[15]<ii & BBFlagC[16]<ii & BBFlagC[17]<ii & BBFlagC[18]<ii & BBFlagC[19]<ii & BBFlagC[20]<ii & BBFlagC[21]<ii; //BB-C 11-21
        }
        if (v0variation == 2) if (flagvariation == 1){
            SelGoodEvent[bunchrange][v0variation][flagvariation] = BBFlagC[11]<ii & BBFlagC[12]<ii & BBFlagC[13]<ii & BBFlagC[14]<ii & BBFlagC[15]<ii & BBFlagC[16]<ii & BBFlagC[17]<ii & BBFlagC[18]<ii & BBFlagC[19]<ii& BBFlagC[20]<ii& BBFlagC[21]<ii; //BB-C 11-21
        }
        if (v0variation == 0 ||v0variation == 1) if (flagvariation == 0 || flagvariation == 1){
            SelGoodEvent[bunchrange][v0variation][flagvariation] &= BBFlagA[9]<ii & BBFlagA[8]<ii & BBFlagA[7]<ii  & BBFlagA[6]<ii  & BBFlagA[5]<ii  & BBFlagA[4]<ii  & BBFlagA[3]<ii & BBFlagA[2]<ii & BBFlagA[1]<ii & BBFlagA[0]<ii; //BB-A 0-9
        }
        if (v0variation == 0 ||v0variation == 2) if (flagvariation == 0 || flagvariation == 1){
            SelGoodEvent[bunchrange][v0variation][flagvariation] &= BBFlagC[9]<ii & BBFlagC[8]<ii & BBFlagC[7]<ii  & BBFlagC[6]<ii  & BBFlagC[5]<ii  & BBFlagC[4]<ii & BBFlagC[3]<ii & BBFlagC[2]<ii & BBFlagC[1]<ii & BBFlagC[0]<ii; //BB-C 0-9
        }
    }
    
    //Bunch range 1-9 and 11-19
    else if (bunchrange == 1){
        if (v0variation == 0 || v0variation == 1) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEvent[bunchrange][v0variation][flagvariation] = BGFlagA[11]<ii & BGFlagA[12]<ii & BGFlagA[13]<ii & BGFlagA[14]<ii & BGFlagA[15]<ii & BGFlagA[16]<ii & BGFlagA[17]<ii & BGFlagA[18]<ii & BGFlagA[19]<ii;//BG-A 11-19
        }
        if (v0variation == 0) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEvent[bunchrange][v0variation][flagvariation] &= BGFlagC[11]<ii & BGFlagC[12]<ii & BGFlagC[13]<ii & BGFlagC[14]<ii & BGFlagC[15]<ii & BGFlagC[16]<ii & BGFlagC[17]<ii & BGFlagC[18]<ii & BGFlagC[19]<ii; //BG-C 11-19
        }
        if (v0variation == 2) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEvent[bunchrange][v0variation][flagvariation] = BGFlagC[11]<ii & BGFlagC[12]<ii & BGFlagC[13]<ii & BGFlagC[14]<ii & BGFlagC[15]<ii & BGFlagC[16]<ii & BGFlagC[17]<ii & BGFlagC[18]<ii & BGFlagC[19]<ii; //BG-C 11-19 for the first loop
        }
        if (v0variation == 0 ||v0variation == 1) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEvent[bunchrange][v0variation][flagvariation] &= BGFlagA[9]<ii & BGFlagA[8]<ii & BGFlagA[7]<ii  & BGFlagA[6]<ii  & BGFlagA[5]<ii  & BGFlagA[4]<ii & BGFlagA[3]<ii & BGFlagA[2]<ii & BGFlagA[1]<ii; //BG-A 1-9
        }
        if (v0variation == 0 ||v0variation == 2) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEvent[bunchrange][v0variation][flagvariation] &= BGFlagC[9]<ii & BGFlagC[8]<ii & BGFlagC[7]<ii  & BGFlagC[6]<ii  & BGFlagC[5]<ii  & BGFlagC[4]<ii & BGFlagC[3]<ii & BGFlagC[2]<ii & BGFlagC[1]<ii; //BG-C 1-9
        }
        if (v0variation == 0 ||v0variation == 1) if (flagvariation == 0){
            SelGoodEvent[bunchrange][v0variation][flagvariation] &= BBFlagA[11]<ii & BBFlagA[12]<ii & BBFlagA[13]<ii & BBFlagA[14]<ii & BBFlagA[15]<ii & BBFlagA[16]<ii & BBFlagA[17]<ii & BBFlagA[18]<ii & BBFlagA[19]<ii; //BB-A 11-19
        }
        if (v0variation == 0 ||v0variation == 1) if (flagvariation == 1){
            SelGoodEvent[bunchrange][v0variation][flagvariation] = BBFlagA[11]<ii & BBFlagA[12]<ii & BBFlagA[13]<ii & BBFlagA[14]<ii & BBFlagA[15]<ii & BBFlagA[16]<ii & BBFlagA[17]<ii & BBFlagA[18]<ii & BBFlagA[19]<ii; //BB-A 11-19 for the first loop
        }
        if (v0variation == 0 ) if (flagvariation == 0 || flagvariation == 1){
            SelGoodEvent[bunchrange][v0variation][flagvariation] &= BBFlagC[11]<ii & BBFlagC[12]<ii & BBFlagC[13]<ii & BBFlagC[14]<ii & BBFlagC[15]<ii & BBFlagC[16]<ii & BBFlagC[17]<ii & BBFlagC[18]<ii & BBFlagC[19]<ii; //BB-C 11-19
        }
        if (v0variation == 2) if (flagvariation == 0){
            SelGoodEvent[bunchrange][v0variation][flagvariation] &= BBFlagC[11]<ii & BBFlagC[12]<ii & BBFlagC[13]<ii & BBFlagC[14]<ii & BBFlagC[15]<ii & BBFlagC[16]<ii & BBFlagC[17]<ii & BBFlagC[18]<ii & BBFlagC[19]<ii; //BB-C 11-19
        }
        if (v0variation == 2) if (flagvariation == 1){
            SelGoodEvent[bunchrange][v0variation][flagvariation] = BBFlagC[11]<ii & BBFlagC[12]<ii & BBFlagC[13]<ii & BBFlagC[14]<ii & BBFlagC[15]<ii & BBFlagC[16]<ii & BBFlagC[17]<ii & BBFlagC[18]<ii & BBFlagC[19]<ii; //BB-C 11-19
        }
        if (v0variation == 0 ||v0variation == 1) if (flagvariation == 0 || flagvariation == 1){
            SelGoodEvent[bunchrange][v0variation][flagvariation] &= BBFlagA[9]<ii & BBFlagA[8]<ii & BBFlagA[7]<ii  & BBFlagA[6]<ii  & BBFlagA[5]<ii  & BBFlagA[4]<ii  & BBFlagA[3]<ii & BBFlagA[2]<ii & BBFlagA[1]<ii; //BB-A 1-9
        }
        if (v0variation == 0 ||v0variation == 2) if (flagvariation == 0 || flagvariation == 1){
            SelGoodEvent[bunchrange][v0variation][flagvariation] &= BBFlagC[9]<ii & BBFlagC[8]<ii & BBFlagC[7]<ii  & BBFlagC[6]<ii  & BBFlagC[5]<ii  & BBFlagC[4]<ii & BBFlagC[3]<ii & BBFlagC[2]<ii & BBFlagC[1]<ii; //BB-C 1-9
        }
    }
    
    //Bunch range 3-9 and 11-17
    else if (bunchrange == 2){
        if (v0variation == 0 || v0variation == 1) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEvent[bunchrange][v0variation][flagvariation] = BGFlagA[11]<ii & BGFlagA[12]<ii & BGFlagA[13]<ii & BGFlagA[14]<ii & BGFlagA[15]<ii & BGFlagA[16]<ii & BGFlagA[17]<ii;//BG-A 11-17
        }
        if (v0variation == 0) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEvent[bunchrange][v0variation][flagvariation] &= BGFlagC[11]<ii & BGFlagC[12]<ii & BGFlagC[13]<ii & BGFlagC[14]<ii & BGFlagC[15]<ii & BGFlagC[16]<ii & BGFlagC[17]<ii; //BG-C 11-17
        }
        if (v0variation == 2) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEvent[bunchrange][v0variation][flagvariation] = BGFlagC[11]<ii & BGFlagC[12]<ii & BGFlagC[13]<ii & BGFlagC[14]<ii & BGFlagC[15]<ii & BGFlagC[16]<ii & BGFlagC[17]<ii; //BG-C 11-17 for the first loop
        }
        if (v0variation == 0 ||v0variation == 1) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEvent[bunchrange][v0variation][flagvariation] &= BGFlagA[9]<ii & BGFlagA[8]<ii & BGFlagA[7]<ii  & BGFlagA[6]<ii  & BGFlagA[5]<ii  & BGFlagA[4]<ii  & BGFlagA[3]<ii; //BG-A 3-9
        }
        if (v0variation == 0 ||v0variation == 2) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEvent[bunchrange][v0variation][flagvariation] &= BGFlagC[9]<ii & BGFlagC[8]<ii & BGFlagC[7]<ii  & BGFlagC[6]<ii  & BGFlagC[5]<ii  & BGFlagC[4]<ii  & BGFlagC[3]<ii; //BG-C 3-9
        }
        if (v0variation == 0 ||v0variation == 1) if (flagvariation == 0){
            SelGoodEvent[bunchrange][v0variation][flagvariation] &= BBFlagA[11]<ii & BBFlagA[12]<ii & BBFlagA[13]<ii & BBFlagA[14]<ii & BBFlagA[15]<ii & BBFlagA[16]<ii & BBFlagA[17]<ii; //BB-A 11-17
        }
        if (v0variation == 0 ||v0variation == 1) if (flagvariation == 1){
            SelGoodEvent[bunchrange][v0variation][flagvariation] = BBFlagA[11]<ii  &  BBFlagA[12]<ii  &  BBFlagA[13]<ii  &  BBFlagA[14]<ii  &  BBFlagA[15]<ii  &  BBFlagA[16]<ii  & BBFlagA[17]<ii; //BB-A 11-17 for the first loop
        }
        if (v0variation == 0 ) if (flagvariation == 0 || flagvariation == 1){
            SelGoodEvent[bunchrange][v0variation][flagvariation] &= BBFlagC[11]<ii  &  BBFlagC[12]<ii  &  BBFlagC[13]<ii  &  BBFlagC[14]<ii  &  BBFlagC[15]<ii  &  BBFlagC[16]<ii  &  BBFlagC[17]<ii; //BB-C 11-17
        }
        if (v0variation == 2) if (flagvariation == 0){
            SelGoodEvent[bunchrange][v0variation][flagvariation] &= BBFlagC[11]<ii & BBFlagC[12]<ii & BBFlagC[13]<ii & BBFlagC[14]<ii & BBFlagC[15]<ii & BBFlagC[16]<ii & BBFlagC[17]<ii; //BB-C 11-17
        }
        if (v0variation == 2) if (flagvariation == 1){
            SelGoodEvent[bunchrange][v0variation][flagvariation] = BBFlagC[11]<ii & BBFlagC[12]<ii & BBFlagC[13]<ii & BBFlagC[14]<ii & BBFlagC[15]<ii & BBFlagC[16]<ii & BBFlagC[17]<ii; //BB-C 11-17
        }
        if (v0variation == 0 ||v0variation == 1) if (flagvariation == 0 || flagvariation == 1){
            SelGoodEvent[bunchrange][v0variation][flagvariation] &= BBFlagA[9]<ii  &  BBFlagA[8]<ii  &  BBFlagA[7]<ii  & BBFlagA[6]<ii  & BBFlagA[5]<ii  & BBFlagA[4]<ii  & BBFlagA[3]<ii; //BB-A 3-9
        }
        if (v0variation == 0 ||v0variation == 2) if (flagvariation == 0 || flagvariation == 1){
            SelGoodEvent[bunchrange][v0variation][flagvariation] &= BBFlagC[9]<ii  &  BBFlagC[8]<ii  &  BBFlagC[7]<ii  &  BBFlagC[6]<ii  &  BBFlagC[5]<ii  &  BBFlagC[4]<ii  & BBFlagC[3]<ii; //BB-C 3-9
        }
    }
    //    cout << "for " << ii << "loop, SelGoodEvent[" << bunchrange << "][" << v0variation << "][" << flagvariation << "] = " << SelGoodEvent[bunchrange][v0variation][flagvariation] << endl;
}


void AliAnalysisBGMonitorQA::SelectADGoodEventWithV0Variation(Int_t bunchrange, Int_t v0variation ,Int_t flagvariation, Int_t ii){
    
    
    
    //Bunch range 3-7 and 11-17
    if (bunchrange == 0) {
        if (v0variation == 0 || v0variation == 1) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] = ADBGFlagA[11]<ii & ADBGFlagA[12]<ii & ADBGFlagA[13]<ii & ADBGFlagA[14]<ii & ADBGFlagA[15]<ii & ADBGFlagA[16]<ii & ADBGFlagA[17]<ii;//BG-A 11-17
        }
        if (v0variation == 0) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] &= ADBGFlagC[11]<ii & ADBGFlagC[12]<ii & ADBGFlagC[13]<ii & ADBGFlagC[14]<ii & ADBGFlagC[15]<ii & ADBGFlagC[16]<ii & ADBGFlagC[17]<ii; //BG-C 11-17
        }
        if (v0variation == 2) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] = ADBGFlagC[11]<ii & ADBGFlagC[12]<ii & ADBGFlagC[13]<ii & ADBGFlagC[14]<ii & ADBGFlagC[15]<ii & ADBGFlagC[16]<ii & ADBGFlagC[17]<ii; //BG-C 11-17 for the first loop
        }
        if (v0variation == 0 ||v0variation == 1) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] &= ADBGFlagA[7]<ii  & ADBGFlagA[6]<ii  & ADBGFlagA[5]<ii  & ADBGFlagA[4]<ii  & ADBGFlagA[3]<ii; //BG-A 3-7
        }
        if (v0variation == 0 ||v0variation == 2) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] &= ADBGFlagC[7]<ii  & ADBGFlagC[6]<ii  & ADBGFlagC[5]<ii  & ADBGFlagC[4]<ii  & ADBGFlagC[3]<ii; //BG-C 3-7
        }
        if (v0variation == 0 ||v0variation == 1) if (flagvariation == 0){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] &= ADBBFlagA[11]<ii & ADBBFlagA[12]<ii & ADBBFlagA[13]<ii & ADBBFlagA[14]<ii & ADBBFlagA[15]<ii & ADBBFlagA[16]<ii & ADBBFlagA[17]<ii; //BB-A 11-17
        }
        if (v0variation == 0 ||v0variation == 1) if (flagvariation == 1){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] = ADBBFlagA[11]<ii & ADBBFlagA[12]<ii & ADBBFlagA[13]<ii & ADBBFlagA[14]<ii & ADBBFlagA[15]<ii & ADBBFlagA[16]<ii & ADBBFlagA[17]<ii; //BB-A 11-17 for the first loop
        }
        if (v0variation == 0 ) if (flagvariation == 0 || flagvariation == 1){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] &= ADBBFlagC[11]<ii & ADBBFlagC[12]<ii & ADBBFlagC[13]<ii & ADBBFlagC[14]<ii & ADBBFlagC[15]<ii & ADBBFlagC[16]<ii & ADBBFlagC[17]<ii; //BB-C 11-17
        }
        if (v0variation == 2) if (flagvariation == 0){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] &= ADBBFlagC[11]<ii & ADBBFlagC[12]<ii & ADBBFlagC[13]<ii & ADBBFlagC[14]<ii & ADBBFlagC[15]<ii & ADBBFlagC[16]<ii & ADBBFlagC[17]<ii; //BB-C 11-17
        }
        if (v0variation == 2) if (flagvariation == 1){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] = ADBBFlagC[11]<ii & ADBBFlagC[12]<ii & ADBBFlagC[13]<ii & ADBBFlagC[14]<ii & ADBBFlagC[15]<ii & ADBBFlagC[16]<ii & ADBBFlagC[17]<ii; //BB-C 11-17
        }
        if (v0variation == 0 ||v0variation == 1) if (flagvariation == 0 || flagvariation == 1){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] &= ADBBFlagA[7]<ii  & ADBBFlagA[6]<ii  & ADBBFlagA[5]<ii  & ADBBFlagA[4]<ii  & ADBBFlagA[3]<ii; //BB-A 3-7
        }
        if (v0variation == 0 ||v0variation == 2) if (flagvariation == 0 || flagvariation == 1){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] &= ADBBFlagC[7]<ii  & ADBBFlagC[6]<ii  & ADBBFlagC[5]<ii  & ADBBFlagC[4]<ii  & ADBBFlagC[3]<ii; //BB-C 3-7
        }
    }
    
    //Bunch range 3-8 and 11-17
    else if (bunchrange == 1){
        if (v0variation == 0 || v0variation == 1) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] = ADBGFlagA[11]<ii & ADBGFlagA[12]<ii & ADBGFlagA[13]<ii & ADBGFlagA[14]<ii & ADBGFlagA[15]<ii & ADBGFlagA[16]<ii & ADBGFlagA[17]<ii;//BG-A 11-17
        }
        if (v0variation == 0) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] &= ADBGFlagC[11]<ii & ADBGFlagC[12]<ii & ADBGFlagC[13]<ii & ADBGFlagC[14]<ii & ADBGFlagC[15]<ii & ADBGFlagC[16]<ii & ADBGFlagC[17]<ii; //BG-C 11-17
        }
        if (v0variation == 2) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] = ADBGFlagC[11]<ii & ADBGFlagC[12]<ii & ADBGFlagC[13]<ii & ADBGFlagC[14]<ii & ADBGFlagC[15]<ii & ADBGFlagC[16]<ii & ADBGFlagC[17]<ii; //BG-C 11-17 for the first loop
        }
        if (v0variation == 0 ||v0variation == 1) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] &= ADBGFlagA[8]<ii & ADBGFlagA[7]<ii  & ADBGFlagA[6]<ii  & ADBGFlagA[5]<ii  & ADBGFlagA[4]<ii  & ADBGFlagA[3]<ii; //BG-A 3-8
        }
        if (v0variation == 0 ||v0variation == 2) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] &= ADBGFlagC[8]<ii & ADBGFlagC[7]<ii  & ADBGFlagC[6]<ii  & ADBGFlagC[5]<ii  & ADBGFlagC[4]<ii  & ADBGFlagC[3]<ii; //BG-C 3-8
        }
        if (v0variation == 0 ||v0variation == 1) if (flagvariation == 0){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] &= ADBBFlagA[11]<ii & ADBBFlagA[12]<ii & ADBBFlagA[13]<ii & ADBBFlagA[14]<ii & ADBBFlagA[15]<ii & ADBBFlagA[16]<ii & ADBBFlagA[17]<ii; //BB-A 11-17
        }
        if (v0variation == 0 ||v0variation == 1) if (flagvariation == 1){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] = ADBBFlagA[11]<ii & ADBBFlagA[12]<ii & ADBBFlagA[13]<ii & ADBBFlagA[14]<ii & ADBBFlagA[15]<ii & ADBBFlagA[16]<ii & ADBBFlagA[17]<ii; //BB-A 11-17 for the first loop
        }
        if (v0variation == 0 ) if (flagvariation == 0 || flagvariation == 1){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] &= ADBBFlagC[11]<ii & ADBBFlagC[12]<ii & ADBBFlagC[13]<ii & ADBBFlagC[14]<ii & ADBBFlagC[15]<ii & ADBBFlagC[16]<ii & ADBBFlagC[17]<ii; //BB-C 11-17
        }
        if (v0variation == 2) if (flagvariation == 0){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] &= ADBBFlagC[11]<ii & ADBBFlagC[12]<ii & ADBBFlagC[13]<ii & ADBBFlagC[14]<ii & ADBBFlagC[15]<ii & ADBBFlagC[16]<ii & ADBBFlagC[17]<ii; //BB-C 11-17
        }
        if (v0variation == 2) if (flagvariation == 1){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] = ADBBFlagC[11]<ii & ADBBFlagC[12]<ii & ADBBFlagC[13]<ii & ADBBFlagC[14]<ii & ADBBFlagC[15]<ii & ADBBFlagC[16]<ii & ADBBFlagC[17]<ii; //BB-C 11-17
        }
        if (v0variation == 0 ||v0variation == 1) if (flagvariation == 0 || flagvariation == 1){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] &= ADBBFlagA[8]<ii & ADBBFlagA[7]<ii  & ADBBFlagA[6]<ii  & ADBBFlagA[5]<ii  & ADBBFlagA[4]<ii  & ADBBFlagA[3]<ii; //BB-A 3-8
        }
        if (v0variation == 0 ||v0variation == 2) if (flagvariation == 0 || flagvariation == 1){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] &= ADBBFlagC[8]<ii & ADBBFlagC[7]<ii  & ADBBFlagC[6]<ii  & ADBBFlagC[5]<ii  & ADBBFlagC[4]<ii  & ADBBFlagC[3]<ii; //BB-C 3-8
        }
    }
    
    //Bunch range 3-9 and 11-17
    else if (bunchrange == 2){
        if (v0variation == 0 || v0variation == 1) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] = ADBGFlagA[11]<ii & ADBGFlagA[12]<ii & ADBGFlagA[13]<ii & ADBGFlagA[14]<ii & ADBGFlagA[15]<ii & ADBGFlagA[16]<ii & ADBGFlagA[17]<ii;//BG-A 11-17
        }
        if (v0variation == 0) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] &= ADBGFlagC[11]<ii & ADBGFlagC[12]<ii & ADBGFlagC[13]<ii & ADBGFlagC[14]<ii & ADBGFlagC[15]<ii & ADBGFlagC[16]<ii & ADBGFlagC[17]<ii; //BG-C 11-17
        }
        if (v0variation == 2) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] = ADBGFlagC[11]<ii & ADBGFlagC[12]<ii & ADBGFlagC[13]<ii & ADBGFlagC[14]<ii & ADBGFlagC[15]<ii & ADBGFlagC[16]<ii & ADBGFlagC[17]<ii; //BG-C 11-17 for the first loop
        }
        if (v0variation == 0 ||v0variation == 1) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] &= ADBGFlagA[9]<ii & ADBGFlagA[8]<ii & ADBGFlagA[7]<ii  & ADBGFlagA[6]<ii  & ADBGFlagA[5]<ii  & ADBGFlagA[4]<ii  & ADBGFlagA[3]<ii; //BG-A 3-9
        }
        if (v0variation == 0 ||v0variation == 2) if (flagvariation == 0 || flagvariation == 2){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] &= ADBGFlagC[9]<ii & ADBGFlagC[8]<ii & ADBGFlagC[7]<ii  & ADBGFlagC[6]<ii  & ADBGFlagC[5]<ii  & ADBGFlagC[4]<ii  & ADBGFlagC[3]<ii; //BG-C 3-9
        }
        if (v0variation == 0 ||v0variation == 1) if (flagvariation == 0){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] &= ADBBFlagA[11]<ii & ADBBFlagA[12]<ii & ADBBFlagA[13]<ii & ADBBFlagA[14]<ii & ADBBFlagA[15]<ii & ADBBFlagA[16]<ii & ADBBFlagA[17]<ii; //BB-A 11-17
        }
        if (v0variation == 0 ||v0variation == 1) if (flagvariation == 1){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] = ADBBFlagA[11]<ii & ADBBFlagA[12]<ii & ADBBFlagA[13]<ii & ADBBFlagA[14]<ii & ADBBFlagA[15]<ii & ADBBFlagA[16]<ii & ADBBFlagA[17]<ii; //BB-A 11-17 for the first loop
        }
        if (v0variation == 0 ) if (flagvariation == 0 || flagvariation == 1){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] &= ADBBFlagC[11]<ii & ADBBFlagC[12]<ii & ADBBFlagC[13]<ii & ADBBFlagC[14]<ii & ADBBFlagC[15]<ii & ADBBFlagC[16]<ii & ADBBFlagC[17]<ii; //BB-C 11-17
        }
        if (v0variation == 2) if (flagvariation == 0){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] &= ADBBFlagC[11]<ii & ADBBFlagC[12]<ii & ADBBFlagC[13]<ii & ADBBFlagC[14]<ii & ADBBFlagC[15]<ii & ADBBFlagC[16]<ii & ADBBFlagC[17]<ii; //BB-C 11-17
        }
        if (v0variation == 2) if (flagvariation == 1){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] = ADBBFlagC[11]<ii & ADBBFlagC[12]<ii & ADBBFlagC[13]<ii & ADBBFlagC[14]<ii & ADBBFlagC[15]<ii & ADBBFlagC[16]<ii & ADBBFlagC[17]<ii; //BB-C 11-17
        }
        if (v0variation == 0 ||v0variation == 1) if (flagvariation == 0 || flagvariation == 1){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] &= ADBBFlagA[9]<ii & ADBBFlagA[8]<ii & ADBBFlagA[7]<ii  & ADBBFlagA[6]<ii  & ADBBFlagA[5]<ii  & ADBBFlagA[4]<ii  & ADBBFlagA[3]<ii; //BB-A 3-9
        }
        if (v0variation == 0 ||v0variation == 2) if (flagvariation == 0 || flagvariation == 1){
            SelGoodEventAD[bunchrange][v0variation][flagvariation] &= ADBBFlagC[9]<ii & ADBBFlagC[8]<ii & ADBBFlagC[7]<ii  & ADBBFlagC[6]<ii  & ADBBFlagC[5]<ii  & ADBBFlagC[4]<ii  & ADBBFlagC[3]<ii; //BB-C 3-9
        }
    }
    //    cout << "for " << ii << "loop, SelGoodEventAD[" << bunchrange << "][" << v0variation << "][" << flagvariation << "] = " << SelGoodEventAD[bunchrange][v0variation][flagvariation] << endl;
}