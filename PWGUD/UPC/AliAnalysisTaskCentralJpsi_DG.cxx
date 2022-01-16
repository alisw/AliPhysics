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

// david.grund(at)cern.ch

// C++ headers:
#include <iostream>
#include <string.h>

// Root headers
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TClonesArray.h"
#include "TRandom.h"
#include "TRandom3.h"

// AliRoot headers:
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliPIDResponse.h"
#include "AliTimeRangeCut.h"
#include "AliDataFile.h"
#include "AliOADBContainer.h"
#include "AliVVertex.h"
#include "AliVVZERO.h"
#include "AliVAD.h"
#include "AliESDZDC.h"
#include "AliTOFTriggerMask.h"

// Other AliRoot headers
#include "AliESDEvent.h" 
#include "AliESDtrack.h" 
#include "AliESDtrackCuts.h"
#include "AliMCEvent.h"

// My headers:
#include "AliAnalysisTaskCentralJpsi_DG.h"

class AliAnalysisTaskCentralJpsi_DG; // analysis class

ClassImp(AliAnalysisTaskCentralJpsi_DG) // classimp: necessary for root

AliAnalysisTaskCentralJpsi_DG::AliAnalysisTaskCentralJpsi_DG() : // initializer list
    AliAnalysisTaskSE(),
    fPIDResponse(0),
    fTimeRangeCut(),
    fTrackCutsBit4(0),
    isMC(kFALSE),
    fEvent(0),
    fOutputList(0),
    fTreeJpsi(0),
    fTreeJpsiMCGen(0),
    fRunNumber(0),
    fTriggerName(0),
    // Histograms:
    hCounterCuts(0),
    hCounterTrigger(0),
    hVertexContrib(0),
    hVertexZ(0),
    hADdecision(0),
    hV0decision(0),
    hTPCdEdx(0),
    hTPCdEdxMuon(0),
    hTPCdEdxElectron(0),
    hPtRecGen(0),
    // PID, sigmas:
    fTrk1dEdx(0),
    fTrk2dEdx(0),
    fTrk1SigIfMu(0),
    fTrk1SigIfEl(0),
    fTrk2SigIfMu(0),
    fTrk2SigIfEl(0),
    // Kinematics:
    fPt(0), fPhi(0), fY(0), fM(0),
    // Two tracks:
    fPt1(0), fPt2(0), fEta1(0), fEta2(0), fPhi1(0), fPhi2(0), fQ1(0), fQ2(0),
    // Vertex info:
    fVertexZ(0), fVertexContrib(0),
    // Info from the detectors:
    // ZDC
    fZNA_energy(0), fZNC_energy(0),
    // V0
    fV0A_dec(0), fV0C_dec(0),
    // AD
    fADA_dec(0), fADC_dec(0),
    // Matching SPD clusters with FOhits
    fMatchingSPD(0),
    fFOCrossFiredChips(),
    // Trigger inputs for MC data
    fSPDfile(0), fTOFfile(0), fLoadedRun(-1), hTOFeff(0), hSPDeff(0), fTOFmask(0),
    // MC kinematics on generator level
    fPtGen(0), fMGen(0), fYGen(0), fPhiGen(0)
{
    // default constructor
}
//_____________________________________________________________________________
AliAnalysisTaskCentralJpsi_DG::AliAnalysisTaskCentralJpsi_DG(const char* name) : // initializer list
    AliAnalysisTaskSE(name),
    fPIDResponse(0),
    fTimeRangeCut(),
    fTrackCutsBit4(0),
    isMC(kFALSE),
    fEvent(0),
    fOutputList(0),
    fTreeJpsi(0),
    fTreeJpsiMCGen(0),
    fRunNumber(0),
    fTriggerName(0),
    // Histograms:
    hCounterCuts(0),
    hCounterTrigger(0),
    hVertexContrib(0),
    hVertexZ(0),
    hADdecision(0),
    hV0decision(0),
    hTPCdEdx(0),
    hTPCdEdxMuon(0),
    hTPCdEdxElectron(0),
    hPtRecGen(0),
    // PID, sigmas:
    fTrk1dEdx(0),
    fTrk2dEdx(0),
    fTrk1SigIfMu(0),
    fTrk1SigIfEl(0),
    fTrk2SigIfMu(0),
    fTrk2SigIfEl(0),
    // Kinematics:
    fPt(0), fPhi(0), fY(0), fM(0),
    // Two tracks:
    fPt1(0), fPt2(0), fEta1(0), fEta2(0), fPhi1(0), fPhi2(0), fQ1(0), fQ2(0),
    // Vertex info:
    fVertexZ(0), fVertexContrib(0),
    // Info from the detectors:
    // ZDC
    fZNA_energy(0), fZNC_energy(0),
    // V0
    fV0A_dec(0), fV0C_dec(0),
    // AD
    fADA_dec(0), fADC_dec(0),
    // Matching SPD clusters with FOhits
    fMatchingSPD(0),
    fFOCrossFiredChips(),
    // Trigger inputs for MC data
    fSPDfile(0), fTOFfile(0), fLoadedRun(-1), hTOFeff(0), hSPDeff(0), fTOFmask(0),
    // MC kinematics on generator level
    fPtGen(0), fMGen(0), fYGen(0), fPhiGen(0)
{
    // constructor
    for(Int_t i = 0; i < 11; i++) fTriggerInputsMC[i] = kFALSE;

    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                        // it does its work automatically
    DefineOutput(1, TTree::Class());
    DefineOutput(2, TTree::Class());
    DefineOutput(3, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTaskCentralJpsi_DG::~AliAnalysisTaskCentralJpsi_DG()
{
    // destructor

    if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis){
        delete fOutputList;
        fOutputList = 0;
        delete fTreeJpsi;
        fTreeJpsi = 0;
        delete fTreeJpsiMCGen;
        fTreeJpsiMCGen = 0;
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskCentralJpsi_DG::UserCreateOutputObjects()
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();

    fTrackCutsBit4 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
    
    // ##########################################################
    // OUTPUT TREE:

    fTreeJpsi = new TTree("fTreeJpsi", "fTreeJpsi");
    // Basic things:
    fTreeJpsi->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");
    fTreeJpsi->Branch("fTriggerName", &fTriggerName);
    // PID, sigmas:
    fTreeJpsi->Branch("fTrk1dEdx", &fTrk1dEdx, "fTrk1dEdx/D");
    fTreeJpsi->Branch("fTrk2dEdx", &fTrk2dEdx, "fTrk2dEdx/D");
    fTreeJpsi->Branch("fTrk1SigIfMu", &fTrk1SigIfMu, "fTrk1SigIfMu/D");
    fTreeJpsi->Branch("fTrk1SigIfEl", &fTrk1SigIfEl, "fTrk1SigIfEl/D");
    fTreeJpsi->Branch("fTrk2SigIfMu", &fTrk2SigIfMu, "fTrk2SigIfMu/D");
    fTreeJpsi->Branch("fTrk2SigIfEl", &fTrk2SigIfEl, "fTrk2SigIfEl/D");
    // Kinematics:
    fTreeJpsi->Branch("fPt", &fPt, "fPt/D");
    fTreeJpsi->Branch("fPhi", &fPhi, "fPhi/D");
    fTreeJpsi->Branch("fY", &fY, "fY/D");
    fTreeJpsi->Branch("fM", &fM, "fM/D");
    // Two tracks:
    fTreeJpsi->Branch("fPt1", &fPt1, "fPt1/D");
    fTreeJpsi->Branch("fPt2", &fPt2, "fPt2/D");
    fTreeJpsi->Branch("fEta1", &fEta1, "fEta1/D");
    fTreeJpsi->Branch("fEta2", &fEta2, "fEta2/D");
    fTreeJpsi->Branch("fPhi1", &fPhi1, "fPhi1/D");
    fTreeJpsi->Branch("fPhi2", &fPhi2, "fPhi2/D");
    fTreeJpsi->Branch("fQ1", &fQ1, "fQ1/D");
    fTreeJpsi->Branch("fQ2", &fQ2, "fQ2/D");
    // Vertex info:
    fTreeJpsi->Branch("fVertexZ", &fVertexZ, "fVertexZ/D");
    fTreeJpsi->Branch("fVertexContrib", &fVertexContrib, "fVertexContrib/I");    
    // Info from the detectors:
    // ZDC:
    fTreeJpsi->Branch("fZNA_energy", &fZNA_energy, "fZNA_energy/D");
    fTreeJpsi->Branch("fZNC_energy", &fZNC_energy, "fZNC_energy/D");
    fTreeJpsi->Branch("fZNA_time", &fZNA_time[0], "fZNA_time[4]/D");
    fTreeJpsi->Branch("fZNC_time", &fZNC_time[0], "fZNC_time[4]/D");
    // V0:
    fTreeJpsi->Branch("fV0A_dec", &fV0A_dec, "fV0A_dec/I");
    fTreeJpsi->Branch("fV0C_dec", &fV0C_dec, "fV0C_dec/I");
    // AD:
    fTreeJpsi->Branch("fADA_dec", &fADA_dec, "fADA_dec/I");
    fTreeJpsi->Branch("fADC_dec", &fADC_dec, "fADC_dec/I");
    // Matching SPD clusters with FOhits:
    fTreeJpsi->Branch("fMatchingSPD", &fMatchingSPD, "fMatchingSPD/O");
    if(isMC){
        // Replayed trigger inputs:
        fTreeJpsi->Branch("fTriggerInputsMC", &fTriggerInputsMC[0], "fTriggerInputsMC[11]/O"); 
        // Kinematics, MC gen:
        fTreeJpsi->Branch("fPtGen", &fPtGen, "fPtGen/D");
        fTreeJpsi->Branch("fMGen", &fMGen, "fMGen/D");
        fTreeJpsi->Branch("fYGen", &fYGen, "fYGen/D");
        fTreeJpsi->Branch("fPhiGen", &fPhiGen, "fPhiGen/D");
    }    
    
    PostData(1, fTreeJpsi);

    // ##########################################################
    // OUTPUT TREE MC GEN:

    fTreeJpsiMCGen = new TTree("fTreeJpsiMCGen", "fTreeJpsiMCGen");
    if(isMC){
        // Run number:
        fTreeJpsiMCGen->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");
        // Kinematics:
        fTreeJpsiMCGen->Branch("fPtGen", &fPtGen, "fPtGen/D");
        fTreeJpsiMCGen->Branch("fMGen", &fMGen, "fMGen/D");
        fTreeJpsiMCGen->Branch("fYGen", &fYGen, "fYGen/D");
        fTreeJpsiMCGen->Branch("fPhiGen", &fPhiGen, "fPhiGen/D");
    }

    PostData(2, fTreeJpsiMCGen);

    // ##########################################################
    // OUTPUT LIST:

    fOutputList = new TList();       
    fOutputList->SetOwner(kTRUE); 

    // Counter for events passing each cut
    hCounterCuts = new TH1D("hCounterCuts", "# of events passing each cut", 5, -0.5, 4.5);
    hCounterCuts->GetXaxis()->SetBinLabel(1,"0: non-empty ev");
    //hCounterCuts->GetXaxis()->SetBinLabel(2,"1: vrtx contrib");
    //hCounterCuts->GetXaxis()->SetBinLabel(3,"2: vrtx Z dist");
    hCounterCuts->GetXaxis()->SetBinLabel(2,"1: two good trks");
    if(!isMC) hCounterCuts->GetXaxis()->SetBinLabel(3,"2: CCUP31 trigg");

    fOutputList->Add(hCounterCuts);

    // Number of triggered events per each run
    // For 2018q: first run = 295585, last run = 296623
    // For 2018r: first run = 296690, last run = 297595

    Int_t nFirstRun = 295585;
    Int_t nLastRun = 297595;
    Int_t nRuns = nLastRun - nFirstRun + 1;

    hCounterTrigger = new TH1D("hCounterTrigger", "# of events per run passing central triggers", nRuns, nFirstRun-0.5, nLastRun+0.5);
    fOutputList->Add(hCounterTrigger);

    hVertexZ = new TH1D("hVertexZ","hVertexZ",600,-30,30);
    fOutputList->Add(hVertexZ);

    hVertexContrib = new TH1D("hVertexContrib","hVertexContrib",103,-2,100);
    fOutputList->Add(hVertexContrib);

    hADdecision = new TH2I("hADdecision","hADdecision",7,-2,5,7,-2,5);
    fOutputList->Add(hADdecision);

    hV0decision = new TH2I("hV0decision","hV0decision",7,-2,5,7,-2,5);
    fOutputList->Add(hV0decision);

    hTPCdEdx = new TH2D("hTPCdEdx"," ",140,0,140.,140,0,140.);
    // lepton with a negative charge on a horizontal axis
    hTPCdEdx->GetXaxis()->SetTitle("dE/dx^{TPC}(l^{-}) (a.u.)");
    hTPCdEdx->GetYaxis()->SetTitle("dE/dx^{TPC}(l^{+}) (a.u.)");
    fOutputList->Add(hTPCdEdx);

    hTPCdEdxMuon = new TH2D("hTPCdEdxMuon"," ",140,0,140.,140,0,140.);
    // muon with a negative charge on a horizontal axis
    hTPCdEdxMuon->GetXaxis()->SetTitle("dE/dx^{TPC}(#mu^{-}) (a.u.)");
    hTPCdEdxMuon->GetYaxis()->SetTitle("dE/dx^{TPC}(#mu^{+}) (a.u.)");
    fOutputList->Add(hTPCdEdxMuon);

    hTPCdEdxElectron = new TH2D("hTPCdEdxElectron"," ",140,0,140.,140,0,140.);
    // electron with a negative charge on a horizontal axis
    hTPCdEdxElectron->GetXaxis()->SetTitle("dE/dx^{TPC}(e^{-}) (a.u.)");
    hTPCdEdxElectron->GetYaxis()->SetTitle("dE/dx^{TPC}(e^{+}) (a.u.)");
    fOutputList->Add(hTPCdEdxElectron);

    if(isMC){
        Int_t n_bins = 2000;
        // x axis = pT from the generator, y axis = pT reconstructed
        hPtRecGen = new TH2D("hPtRecGen", "pT rec vs pT gen", n_bins, 0., 2., n_bins, 0., 2.);
        fOutputList->Add(hPtRecGen);
    }

    // https://github.com/alisw/AliPhysics/blob/6015b235c21bc8b9d00af9a764be1fb58f7bb32d/PWGCF/Correlations/DPhi/AliPhiCorrelationsQATask.cxx
    fTrackCutsBit4->DefineHistograms(kBlue);
    fTrackCutsBit4->SetName("track_cuts");
    fOutputList->Add(fTrackCutsBit4);

    PostData(3, fOutputList); 

}
//_____________________________________________________________________________
void AliAnalysisTaskCentralJpsi_DG::UserExec(Option_t *)
{
    // Analysis of the 2018q, 2018r datasets in central rapidity
    // Counter of selection criteria
    Int_t iSelectionCounter = 0;

    // New event
    fEvent = InputEvent();

    // ##########################################################
        // CUT 0 
        // If the event is empty, it is skipped
        if(!fEvent)
        {                                          
            PostData(1, fTreeJpsi);
            PostData(2, fOutputList);
            return;
        }                               
        hCounterCuts->Fill(iSelectionCounter);
        iSelectionCounter++;
    // ##########################################################

    // Get the run number associated with the ESD event
    fRunNumber = fEvent->GetRunNumber(); 

    // Time range cuts: https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGRunList18r 
    fTimeRangeCut.InitFromEvent(InputEvent());
    if(fTimeRangeCut.CutEvent(InputEvent())) return;

    // Check the central UPC trigger CCUP31
    // Fill the hCounterTrigger (to calculate the integrated lumi of the analysed sample)
    // Cut on trigger will be performed later (# 4)
    fTriggerName = fEvent->GetFiredTriggerClasses();
    Bool_t triggered = kFALSE;
    if(fRunNumber < 295881){
        if(fTriggerName.Contains("CCUP31-B-NOPF-CENTNOTRD")){
            triggered = kTRUE;
            hCounterTrigger->Fill(fRunNumber);
        } 
    }
    if(fRunNumber >= 295881){
        if(fTriggerName.Contains("CCUP31-B-SPD2-CENTNOTRD")){
            triggered = kTRUE;
            hCounterTrigger->Fill(fRunNumber);
        }
    }

    // if MC data
    if(isMC){
        // Replay triggers
        if(fRunNumber != fLoadedRun){
            if(!fSPDfile) fSPDfile = AliDataFile::OpenOADB("PWGUD/UPC/SPDEfficiency18qr.root");
            if(!fTOFfile) fTOFfile = AliDataFile::OpenOADB("PWGUD/TOFTriggerEfficiency.root");
            AliOADBContainer* fTOFcont = (AliOADBContainer*)fTOFfile->Get("TOFTriggerEfficiency");
            hTOFeff  = (TH2F*)fTOFcont->GetObject(fRunNumber,"Default");
            
            if(fSPDfile->Get(Form("eff%i",fRunNumber))) hSPDeff = (TH1D*) fSPDfile->Get(Form("eff%i",fRunNumber));
            Int_t tempRun = fRunNumber;
            while(!hSPDeff){
                tempRun--;
                hSPDeff = (TH1D*) fSPDfile->Get(Form("eff%i",tempRun));
            }
            fLoadedRun = fRunNumber;
        }
        ReplayTriggersMC(fEvent);   
        // Run analysis on the generator level
        RunMCGenLevel();     
    }

    // Forward detectors
    // Save info from these detectors to the analysis tree
    // and fill the 2D histograms with V0 and AD responses
    
    // AD
    AliVAD *fADdata = fEvent->GetADData();
    fADA_dec = fADdata->GetADADecision();
    fADC_dec = fADdata->GetADCDecision();
    hADdecision->Fill(fADA_dec,fADC_dec);

    // V0
    AliVVZERO *fV0data = fEvent->GetVZEROData();
    fV0A_dec = fV0data->GetV0ADecision();
    fV0C_dec = fV0data->GetV0CDecision();
    hV0decision->Fill(fV0A_dec,fV0C_dec);

    // Get info about primary vertex
    const AliVVertex *fVertex = fEvent->GetPrimaryVertex();
    // Fill the trees and histograms
    fVertexZ = fVertex->GetZ();
    hVertexZ->Fill(fVertexZ);
    fVertexContrib = fVertex->GetNContributors();
    hVertexContrib->Fill(fVertexContrib);
    // Cuts on fVertexZ (must be lower than 15 cm) and fVertexContrib (at least two tracks associated with the vertex)
    // => will be done offline

    // ##########################################################
        // CUT 1
        // Select events with two good central tracks
        Int_t nTrks = fEvent->GetNumberOfTracks();
        Int_t nGoodTracksTPC = 0;
        Int_t nGoodTracksSPD = 0;
        Int_t *fIndicesOfGoodTracks = new Int_t[nTrks];
        for(Int_t iTrk(0); iTrk < nTrks; iTrk++){
            AliESDtrack *trk = dynamic_cast<AliESDtrack*>(fEvent->GetTrack(iTrk));

            // If the track is empty
            if(!trk) continue;

            // Count good TPC tracks:
            Bool_t goodTrackTPC = fTrackCutsBit4->AcceptTrack(trk);
            if(goodTrackTPC) nGoodTracksTPC++;
            
            // Count good SPD tracks:
            if(goodTrackTPC && trk->HasPointOnITSLayer(0) && trk->HasPointOnITSLayer(1)){
                // If the track satisfies both the SPD and TPC criterion, add it to the list:
                fIndicesOfGoodTracks[nGoodTracksSPD] = iTrk;
                nGoodTracksSPD++;  
            }
        }
        // Continue only if two good central tracks are found
        if(!(nGoodTracksSPD == 2 && nGoodTracksTPC == 2)){                                          
            PostData(1, fTreeJpsi);
            PostData(2, fTreeJpsiMCGen);
            PostData(3, fOutputList);
            delete [] fIndicesOfGoodTracks; 
            return;
        }
        hCounterCuts->Fill(iSelectionCounter);
        iSelectionCounter++;
    // ##########################################################
        // CUT 2
        // Central UPC trigger CCUP31
        if(!isMC){ // skipped for MC
            if(!triggered)
            {
                PostData(1, fTreeJpsi);
                PostData(2, fTreeJpsiMCGen);
                PostData(3, fOutputList);
                delete [] fIndicesOfGoodTracks; 
                return;        
            }
            hCounterCuts->Fill(iSelectionCounter);
            iSelectionCounter++;
        }
    // ##########################################################
    // Two track loop

    // PID and kinematics 
    Double_t massMuon = 0.105658;       // GeV/c^2
    Double_t massElectron = 0.000511;   // GeV/c^2
    
    // Get the tracks
    AliESDtrack *trk1 = dynamic_cast<AliESDtrack*>(fEvent->GetTrack(fIndicesOfGoodTracks[0]));
    AliESDtrack *trk2 = dynamic_cast<AliESDtrack*>(fEvent->GetTrack(fIndicesOfGoodTracks[1]));
    // Get sigmas if electrons/muons
    fTrk1SigIfMu = fPIDResponse->NumberOfSigmasTPC(trk1, AliPID::kMuon);
    fTrk1SigIfEl = fPIDResponse->NumberOfSigmasTPC(trk1, AliPID::kElectron);
    fTrk2SigIfMu = fPIDResponse->NumberOfSigmasTPC(trk2, AliPID::kMuon);
    fTrk2SigIfEl = fPIDResponse->NumberOfSigmasTPC(trk2, AliPID::kElectron);
    // Decide whether it is a muon/electron pair according to lower sum of squares of the sigmas
    Double_t isMuonPair = fTrk1SigIfMu*fTrk1SigIfMu + fTrk2SigIfMu*fTrk2SigIfMu;
    Double_t isElectronPair = fTrk1SigIfEl*fTrk1SigIfEl + fTrk2SigIfEl*fTrk2SigIfEl;

    // Decide if muons/electrons, then assign a proper mass
    Double_t massTracks = -1;
    if(isMuonPair < isElectronPair) massTracks = massMuon;
    else massTracks = massElectron;

    // Track-track kinematics
    // Fill the 4-vector of the first track
    TLorentzVector vTrk1;
    vTrk1.SetPtEtaPhiM(trk1->Pt(), trk1->Eta(), trk1->Phi(), massTracks);
    // Fill the 4-vector of the second track
    TLorentzVector vTrk2;
    vTrk2.SetPtEtaPhiM(trk2->Pt(), trk2->Eta(), trk2->Phi(), massTracks);
    // Vector of trk+trk: add up the two
    TLorentzVector vTrkTrk = vTrk1 + vTrk2;

    // Set tree variables
    fPt = vTrkTrk.Pt(); 
    fPhi = vTrkTrk.Phi();
    fY = vTrkTrk.Rapidity(); 
    fM = vTrkTrk.M();
    fPt1 = trk1->Pt(); 
    fPt2 = trk2->Pt();
    fEta1 = trk1->Eta(); 
    fEta2 = trk2->Eta();
    fPhi1 = trk1->Phi();
    fPhi2 = trk2->Phi();
    fQ1 = trk1->Charge(); 
    fQ2 = trk2->Charge();

    // Get dE/dx TPC for both tracks
    fTrk1dEdx = trk1->GetTPCsignal();
    fTrk2dEdx = trk2->GetTPCsignal();

    // If opposite sign tracks and mass of a J/psi candidate within (2.2, 5.0) GeV
    // => fill the histograms with dE/dx
    if(vTrkTrk.M() > 2.2 && vTrkTrk.M() < 5.0 && (fQ1 * fQ2 < 0)){
        if(fQ1 < 0) hTPCdEdx->Fill(fTrk1dEdx, fTrk2dEdx);
        else        hTPCdEdx->Fill(fTrk2dEdx, fTrk1dEdx);
        // if considered muons
        if(isMuonPair < isElectronPair){ 
            if(fQ1 < 0) hTPCdEdxMuon->Fill(fTrk1dEdx, fTrk2dEdx);
            else        hTPCdEdxMuon->Fill(fTrk2dEdx, fTrk1dEdx);
        // if considered electrons
        } else {
            if(fQ1 < 0) hTPCdEdxElectron->Fill(fTrk1dEdx, fTrk2dEdx);
            else        hTPCdEdxElectron->Fill(fTrk2dEdx, fTrk1dEdx);
        }
    }

    // Fill the 2D histogram of pT gen and pT rec
    if(isMC) hPtRecGen->Fill(fPtGen,fPt);

    // Check if SPD clusters match FOhits (according to macro by MB)
    Int_t crossedFO[4];
    TBits fFOCrossedChips(1200); 
    const AliVMultiplicity *mult = fEvent->GetMultiplicity();
    TBits fFOFiredChips = mult->GetFastOrFiredChips();

    fFOCrossedChips.ResetAllBits(kFALSE);

    for(Int_t iTrack = 0; iTrack < 2; iTrack++){
        AliESDtrack *trk = dynamic_cast<AliESDtrack*>(fEvent->GetTrack(fIndicesOfGoodTracks[iTrack]));
	    if(!trk) continue;
        crossedFO[0] = trk->GetITSModuleIndex(0);
        crossedFO[1] = trk->GetITSModuleIndex(1);
        crossedFO[2] = trk->GetITSModuleIndex(6);
        crossedFO[3] = trk->GetITSModuleIndex(7);
        SetCrossed(crossedFO, fFOCrossedChips);
    }    

    fFOCrossFiredChips = fFOCrossedChips & fFOFiredChips;
    fMatchingSPD = IsSTGFired(fFOCrossFiredChips,fRunNumber >= 295753 ? 9 : 3);

    // Clean up
    delete [] fIndicesOfGoodTracks;
    // ##########################################################

    // Data from the ZDC
  	AliESDZDC *fZDCdata = (AliESDZDC*)fEvent->GetZDCData();
  	fZNA_energy = fZDCdata->GetZNATowerEnergy()[0];
  	fZNC_energy = fZDCdata->GetZNCTowerEnergy()[0];
	Int_t detChZNA  = fZDCdata->GetZNATDCChannel();
    Int_t detChZNC  = fZDCdata->GetZNCTDCChannel();
	if(fEvent->GetRunNumber() >= 245726 && fEvent->GetRunNumber() <= 245793) detChZNA = 10;
  	for(Int_t i = 0; i < 4; i++){ 
  		fZNA_time[i] = fZDCdata->GetZDCTDCCorrected(detChZNA,i);
  		fZNC_time[i] = fZDCdata->GetZDCTDCCorrected(detChZNC,i);
	}

    // Fill the analysis tree
    fTreeJpsi->Fill();

    // Finally post the data
    PostData(1, fTreeJpsi);
    PostData(2, fTreeJpsiMCGen);
    PostData(3, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskCentralJpsi_DG::ReplayTriggersMC(AliVEvent *fEvent)
{
    // First set all triggers to kFALSE state, then we check them
    for(Int_t i = 0; i < 11; i++) fTriggerInputsMC[i] = kFALSE;
    UShort_t fTriggerAD = fEvent->GetADData()->GetTriggerBits();
    UShort_t fTriggerVZERO = fEvent->GetVZEROData()->GetTriggerBits();
    UInt_t fL0inputs = fEvent->GetHeader()->GetL0TriggerInputs();

    fTriggerInputsMC[0] = fTriggerVZERO & (1 << 12);// 0VBA VZERO A
    fTriggerInputsMC[1] = fTriggerVZERO & (1 << 13);// 0VBC VZERO C
    fTriggerInputsMC[2] = fTriggerAD & (1 << 12);   // 0UBA ADA
    fTriggerInputsMC[3] = fTriggerAD & (1 << 13);   // 0UBC ADC

    // ------------------------
    // TOF trigger input
    // ------------------------
    const AliTOFHeader *tofH = fEvent->GetTOFHeader();
    fTOFmask = tofH->GetTriggerMask();
  
    Bool_t firedMaxiPhi[36] = {0};
    Int_t NfiredMaxiPads = 0;
 
    for(Int_t ltm = 0; ltm < 72; ltm++){
        Int_t ip = ltm % 36;
        for(Int_t cttm = 0; cttm < 23; cttm++){
            if(fTOFmask->IsON(ltm,cttm) && gRandom->Rndm(1.0) < hTOFeff->GetBinContent(ltm+1,cttm+1)){
                firedMaxiPhi[ip] = kTRUE;
                NfiredMaxiPads++;
            }
        }
    }
    Bool_t offlineOMU = kFALSE;
    for(Int_t ip = 0;ip < 36; ip++){
        if(!firedMaxiPhi[ip]) continue;
        for(Int_t jp = ip + 1;jp < 36; jp++){
            if(!firedMaxiPhi[jp]) continue;
            Int_t DeSlots = jp - ip;
            Int_t AntiDeSlots = 36 - DeSlots;
            if(DeSlots >= 15 && DeSlots <= 18) offlineOMU = kTRUE;
            else if(AntiDeSlots >= 15 && AntiDeSlots <= 18) offlineOMU = kTRUE;
        }
    }
    if(NfiredMaxiPads>6)offlineOMU = kFALSE;
    if(offlineOMU)fTriggerInputsMC[4] = kTRUE;          //0OMU TOF two hits with topology
    if(NfiredMaxiPads >= 2)fTriggerInputsMC[5] = kTRUE;	//0OM2 TOF two hits	

    // ------------------------
    // STG trigger input
    // ------------------------
    //SPD inputs
    const AliVMultiplicity *mult = fEvent->GetMultiplicity();
    Bool_t vPhiInner[20]; for (Int_t i=0; i<20; ++i) vPhiInner[i]=kFALSE;
    Bool_t vPhiOuter[40]; for (Int_t i=0; i<40; ++i) vPhiOuter[i]=kFALSE;

    Int_t nInner(0), nOuter(0);
    for (Int_t i(0); i<1200; ++i){
        Bool_t isFired = (mult->TestFastOrFiredChips(i)) && (gRandom->Rndm(1.0) < hSPDeff->GetBinContent(i+1));
        if (i<400){
            if(isFired)vPhiInner[i/20] = kTRUE;
            nInner += isFired;
        } else {
            if(isFired)vPhiOuter[(i-400)/20] = kTRUE;
            nOuter += isFired;
        }
    }

    // STG input
    Int_t dphiMax = 10; 
    Int_t dphiMin = 9;
    Bool_t tolerance = 1;
    Bool_t firedSTG = 0;
    Bool_t phi[20] = { 0 };
    if(fRunNumber < 295753){ dphiMin = 3; }

    for (Int_t i = 0; i < 20; ++i) {
        if (tolerance) phi[i] = vPhiInner[i] & (vPhiOuter[(2*i)%40] | vPhiOuter[(2*i+1)%40] | vPhiOuter[(2*i+2)%40] | vPhiOuter[(2*i+39)%40]);
        else           phi[i] = vPhiInner[i] & (vPhiOuter[(2*i)%40] | vPhiOuter[(2*i+1)%40]);
    }
    for (Int_t dphi=dphiMin;dphi<=dphiMax;dphi++)
    for (Int_t i=0; i<20; ++i) firedSTG |= phi[i] & phi[(i+dphi)%20];

    // STP input
    Int_t fired = 0;
    for (Int_t i(0); i<10; ++i) {
        for (Int_t j(0); j<2; ++j) {
        const Int_t k(2*i+j);
        fired += (( vPhiOuter[k]    || vPhiOuter[k+1]       ||
                    vPhiOuter[k+2]      )
                && (vPhiOuter[k+20] || vPhiOuter[(k+21)%40] ||
                    vPhiOuter[(k+22)%40])
                && (vPhiInner[i]    || vPhiInner[i+1]       )
                && (vPhiInner[i+10] || vPhiInner[(i+11)%20]));
        }
    }
    //0SMB - At least one hit in SPD
    if (nOuter > 0 || nInner > 0) fTriggerInputsMC[6] = kTRUE;
    //0SM2 - Two hits on outer layer
    if (nOuter > 1) fTriggerInputsMC[7] = kTRUE;
    //0STP - Topological SPD trigger (two pairs)
    if (fired != 0) fTriggerInputsMC[8] = kTRUE;
    //0SH1 - More then 6 hits on outer layer
    if (nOuter >= 7) fTriggerInputsMC[9] = kTRUE;
    if (firedSTG) fTriggerInputsMC[10] = kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskCentralJpsi_DG::RunMCGenLevel()
{
    TLorentzVector vGen, vDecayProduct;
    TDatabasePDG *pdgdat = TDatabasePDG::Instance();

    vGen.SetPtEtaPhiM(0.,0.,0.,0.);

    AliMCEvent *mc = MCEvent();
    if(!mc){
        // Printf("Not found");
        return;
    } 

    for(Int_t imc = 0; imc < mc->GetNumberOfTracks(); imc++) {
        AliMCParticle *mcPart = (AliMCParticle*) mc->GetTrack(imc);
        if(!mcPart) continue;
    
        // kCohJpsiToMu, kIncohJpsiToMu, kTwoGammaToMuMedium:
        // if mu+ or mu- without a mother assigned (STARlight) 
        if(TMath::Abs(mcPart->PdgCode()) == 13 && mcPart->GetMother() == -1){
            // add its 4-vector to vGen
            TParticlePDG *partGen = pdgdat->GetParticle(mcPart->PdgCode());
            vDecayProduct.SetXYZM(mcPart->Px(),mcPart->Py(), mcPart->Pz(),partGen->Mass());
            vGen += vDecayProduct;
        }
        // kCohPsi2sToMuPi, kIncohPsi2sToMuPi:
        // if J/psi with a mother Psi(2s) (STARlight + EvtGen)
        if(TMath::Abs(mcPart->PdgCode()) == 443){
            // get its mother
            AliMCParticle *mcMother = (AliMCParticle*) mc->GetTrack(mcPart->GetMother());
            // if its mother is Psi(2s)
            if(TMath::Abs(mcMother->PdgCode()) == 100443){
                // add its 4-vector to vGen
                TParticlePDG *partGen = pdgdat->GetParticle(mcPart->PdgCode());
                vDecayProduct.SetXYZM(mcPart->Px(),mcPart->Py(), mcPart->Pz(),partGen->Mass());
                vGen += vDecayProduct;
            }
        }
    } // loop over mc particles

    FillMCGenTree(vGen);
}
//_____________________________________________________________________________
void AliAnalysisTaskCentralJpsi_DG::FillMCGenTree(TLorentzVector v)
{
    fPtGen      = v.Pt();
    if(v.E() != v.Pz()) fYGen = v.Rapidity();
    else fYGen  = -999; // when E = Pz, rapidity goes to infty
    fMGen       = v.M();
    fPhiGen     = v.Phi();

    fTreeJpsiMCGen->Fill();
}
//_____________________________________________________________________________
void AliAnalysisTaskCentralJpsi_DG::SetCrossed(Int_t spd[4], TBits &crossed){ 
    // from the macro by MB
    Int_t chipId2;
    for(Int_t iLayer = 0; iLayer < 4; iLayer++)
    if(spd[iLayer] > 0){
        crossed.SetBitNumber(GetChipId(spd[iLayer],chipId2)); 
        crossed.SetBitNumber(chipId2); 
    }
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskCentralJpsi_DG::GetChipId(Int_t index, Int_t &chipId2, Bool_t debug){
    // from the macro by MB
    Int_t status   = (index%1000000)/100000;
    Int_t iModule  = index/1000000;           // 0 - 239
    Int_t iPhi     = iModule/4;               // 0-19 - inner, 20-59 outer
    Int_t iModuleZ = iModule%4;               // 0-3
    Int_t iSign    = (index%100000)/10000;    // 1-4
    Int_t signZ    = iPhi<20 ? (iSign%2==1 ? 1 : -1) : (iSign%2==0 ? 1 : -1); // 1 or -1
    Int_t iX       = (index%10000)/100;       // ??
    Int_t iZ       = index%100;               // 0-36 [mm]
    Int_t signZiZ  = (36-signZ*iZ);
    Int_t chipId   = iModule*5+signZiZ*5/72;
    if (chipId<0) return 1200;
    if (chipId>=1200) return 1201;
    if (signZiZ<0) return 1202;
    if (signZiZ>72) return 1203;
    if (signZiZ==72 && chipId%20==0 && chipId>=400) return 1204;
    chipId2=chipId;

    if (signZiZ==0  && chipId%20!=0)  chipId2=chipId-1;
    if (signZiZ==72 && chipId%20!=19) chipId2=chipId+1;
    if (signZiZ==13)  chipId2=chipId+1;
    if (signZiZ==14)  chipId2=chipId+1;
    if (signZiZ==15)  chipId2=chipId-1;
    if (signZiZ==16)  chipId2=chipId-1;
    if (signZiZ==27)  chipId2=chipId+1;
    if (signZiZ==28)  chipId2=chipId+1;
    if (signZiZ==29)  chipId2=chipId-1;
    if (signZiZ==30)  chipId2=chipId-1;
    if (signZiZ==42)  chipId2=chipId+1;
    if (signZiZ==43)  chipId2=chipId+1;
    if (signZiZ==44)  chipId2=chipId-1;
    if (signZiZ==45)  chipId2=chipId-1;
    if (signZiZ==56)  chipId2=chipId+1;
    if (signZiZ==57)  chipId2=chipId+1;
    if (signZiZ==58)  chipId2=chipId-1;
    if (signZiZ==59)  chipId2=chipId-1;
    if (debug) printf("%4i %4i %3i %3i %3i\n",chipId,chipId2,iX,signZiZ,iSign);
    return chipId;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCentralJpsi_DG::IsSTGFired(TBits bits, Int_t dphiMin, Int_t dphiMax, Bool_t tolerance){
    // from the macro by MB
    Int_t n1 = bits.CountBits(400);
    Int_t n0 = bits.CountBits() - n1;
    //cout<<n0<<" "<<n1<<endl;
    if (n0<1 || n1<1) return 0;
    Bool_t stg = 0;
    Bool_t l0[20]={0};
    Bool_t l1[40]={0};
    Bool_t phi[20]={0};
    for (Int_t i=0;   i< 400; ++i) if (bits.TestBitNumber(i)) l0[      i/20] = 1;
    for (Int_t i=400; i<1200; ++i) if (bits.TestBitNumber(i)) l1[(i-400)/20] = 1;
    for (Int_t i=0; i<20; ++i){
        if (tolerance) phi[i] = l0[i] & (l1[(2*i)%40] | l1[(2*i+1)%40] | l1[(2*i+2)%40] | l1[(2*i+39)%40]);
        else           phi[i] = l0[i] & (l1[(2*i)%40] | l1[(2*i+1)%40]);
    }
    for (Int_t dphi=dphiMin;dphi<=dphiMax;dphi++) for (Int_t i=0; i<20; ++i) stg |= phi[i] & phi[(i+dphi)%20];
    return stg;
}
//_____________________________________________________________________________
void AliAnalysisTaskCentralJpsi_DG::Terminate(Option_t *)
{
    // the end

    fOutputList = dynamic_cast<TList*> (GetOutputData(3));
    if(!fOutputList)
    {
        Printf("ERROR: fOutputList not available");
        return;
    }

    fTrackCutsBit4 = dynamic_cast<AliESDtrackCuts*> (fOutputList->FindObject("track_cuts"));
    if(!fTrackCutsBit4)
    {
        Printf("ERROR: fTrackCutsBit4 not available");
        return;
    }

    TFile* file = TFile::Open("track_cuts.root", "RECREATE");
    fTrackCutsBit4->SaveHistograms();
    file->Write();
    file->Close();
}
//_____________________________________________________________________________