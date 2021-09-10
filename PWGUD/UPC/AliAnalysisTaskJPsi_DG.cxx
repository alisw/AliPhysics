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

/* AliAnaysisTaskMyTask
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */

// C++ headers:
#include <iostream>
#include <fstream>

// Root headers
#include <TMath.h>
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include <TFile.h>
#include <TF2.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TTree.h>
#include <TGraph2D.h>
#include <TStopwatch.h>
#include <TMatrixDSym.h>
#include <TFitResult.h>
#include <TLatex.h>
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TObjString.h"
#include "TList.h"
#include "TChain.h"
#include "TRandom3.h"

// AliRoot headers:
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"
#include "AliTimeRangeCut.h"
#include "AliDataFile.h"
#include "AliOADBContainer.h"

#include "AliESDEvent.h" 
#include "AliESDtrack.h" 
#include "AliESDtrackCuts.h"

// My headers:
#include "AliAnalysisTaskJPsi_DG.h"

class AliAnalysisTaskJPsi_DG; // your analysis class

//using namespace std; // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskJPsi_DG) // classimp: necessary for root

AliAnalysisTaskJPsi_DG::AliAnalysisTaskJPsi_DG() : // initializer list
    AliAnalysisTaskSE(),
    fPIDResponse(0),
    fTrackCutsBit4(0),
    fEvent(0),
    fOutputList(0),
    fTreeJPsi(0),
    fRunNumber(0),
    fTriggerName(0),
    // Histograms:
    hCounterCuts(0),
    hCounterTrigger(0),
    hVertexContrib(0),
    hVertexZ(0),
    hADdecision(0),
    hV0decision(0),
    // PID, sigmas:
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
    fV0A_dec(0), fV0C_dec(0), fV0A_time(0), fV0C_time(0),
    // AD
    fADA_dec(0), fADC_dec(0), fADA_time(0), fADC_time(0),
    // Matching SPD clusters with FOhits
    fMatchingSPD(0)
{
    // default constructor
}
//_____________________________________________________________________________
AliAnalysisTaskJPsi_DG::AliAnalysisTaskJPsi_DG(const char* name) : // initializer list
    AliAnalysisTaskSE(name),
    fPIDResponse(0),
    fTrackCutsBit4(0),
    fEvent(0),
    fOutputList(0),
    fTreeJPsi(0),
    fRunNumber(0),
    fTriggerName(0),
    // Histograms:
    hCounterCuts(0),
    hCounterTrigger(0),
    hVertexContrib(0),
    hVertexZ(0),
    hADdecision(0),
    hV0decision(0),
    // PID, sigmas:
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
    fV0A_dec(0), fV0C_dec(0), fV0A_time(0), fV0C_time(0),
    // AD
    fADA_dec(0), fADC_dec(0), fADA_time(0), fADC_time(0),
    // Matching SPD clusters with FOhits
    fMatchingSPD(0)
{ 
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                        // it does its work automatically
    DefineOutput(1, TTree::Class());
    DefineOutput(2, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
}
//_____________________________________________________________________________
AliAnalysisTaskJPsi_DG::~AliAnalysisTaskJPsi_DG()
{
    // destructor

    if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis){
        delete fOutputList;
        fOutputList = 0x0;
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskJPsi_DG::UserCreateOutputObjects()
{
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();

    fTrackCutsBit4 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
    
    // ##########################################################
    // OUTPUT TREE:

    fTreeJPsi = new TTree("fTreeJPsi", "fTreeJPsi");
    // Basic things:
    fTreeJPsi->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");
    fTreeJPsi->Branch("fTriggerName", &fTriggerName);
    // PID, sigmas:
    fTreeJPsi->Branch("fTrk1SigIfMu", &fTrk1SigIfMu, "fTrk1SigIfMu/D");
    fTreeJPsi->Branch("fTrk1SigIfEl", &fTrk1SigIfEl, "fTrk1SigIfEl/D");
    fTreeJPsi->Branch("fTrk2SigIfMu", &fTrk2SigIfMu, "fTrk2SigIfMu/D");
    fTreeJPsi->Branch("fTrk2SigIfEl", &fTrk2SigIfEl, "fTrk2SigIfEl/D");
    // Kinematics:
    fTreeJPsi->Branch("fPt", &fPt, "fPt/D");
    fTreeJPsi->Branch("fPhi", &fPhi, "fPhi/D");
    fTreeJPsi->Branch("fY", &fY, "fY/D");
    fTreeJPsi->Branch("fM", &fM, "fM/D");
    // Two tracks:
    fTreeJPsi->Branch("fPt1", &fPt1, "fPt1/D");
    fTreeJPsi->Branch("fPt2", &fPt2, "fPt2/D");
    fTreeJPsi->Branch("fEta1", &fEta1, "fEta1/D");
    fTreeJPsi->Branch("fEta2", &fEta2, "fEta2/D");
    fTreeJPsi->Branch("fPhi1", &fPhi1, "fPhi1/D");
    fTreeJPsi->Branch("fPhi2", &fPhi2, "fPhi2/D");
    fTreeJPsi->Branch("fQ1", &fQ1, "fQ1/D");
    fTreeJPsi->Branch("fQ2", &fQ2, "fQ2/D");
    // Vertex info:
    fTreeJPsi->Branch("fVertexZ", &fVertexZ, "fVertexZ/D");
    fTreeJPsi->Branch("fVertexContrib", &fVertexContrib, "fVertexContrib/I");    
    // Info from the detectors:
    // ZDC:
    fTreeJPsi->Branch("fZNA_energy", &fZNA_energy, "fZNA_energy/D");
    fTreeJPsi->Branch("fZNC_energy", &fZNC_energy, "fZNC_energy/D");
    fTreeJPsi->Branch("fZNA_time", &fZNA_time[0], "fZNA_time[4]/D");
    fTreeJPsi->Branch("fZNC_time", &fZNC_time[0], "fZNC_time[4]/D");
    // V0:
    fTreeJPsi->Branch("fV0A_dec", &fV0A_dec, "fV0A_dec/I");
    fTreeJPsi->Branch("fV0C_dec", &fV0C_dec, "fV0C_dec/I");
    fTreeJPsi->Branch("fV0A_time", &fV0A_time, "fV0A_time/D");
    fTreeJPsi->Branch("fV0C_time", &fV0C_time, "fV0C_time/D");
    // AD:
    fTreeJPsi->Branch("fADA_dec", &fADA_dec, "fADA_dec/I");
    fTreeJPsi->Branch("fADC_dec", &fADC_dec, "fADC_dec/I");
    fTreeJPsi->Branch("fADA_time", &fADA_time, "fADA_time/D");
    fTreeJPsi->Branch("fADC_time", &fADC_time, "fADC_time/D");
    // Matching SPD clusters with FOhits:
    fTreeJPsi->Branch("fMatchingSPD", &fMatchingSPD, "fMatchingSPD/O");
    
    PostData(1, fTreeJPsi);

    // ##########################################################
    // OUTPUT LIST:

    fOutputList = new TList();       
    fOutputList->SetOwner(kTRUE); 

    // Counter for events passing each cut
    hCounterCuts = new TH1F("hCounterCuts", "# of events passing each cut", 10, -0.5, 9.5);
    hCounterCuts->GetXaxis()->SetBinLabel(1,"0: non-empty ev");
    hCounterCuts->GetXaxis()->SetBinLabel(2,"1: vrtx contrib");
    hCounterCuts->GetXaxis()->SetBinLabel(3,"2: vrtx Z dist");
    hCounterCuts->GetXaxis()->SetBinLabel(4,"3: two good trks");
    hCounterCuts->GetXaxis()->SetBinLabel(5,"4: CCUP31 trigg");

    fOutputList->Add(hCounterCuts);

    // Number of triggered events per each run
    // For 2018q: first run = 295585, last run = 296623
    // For 2018r: first run = 296690, last run = 297595

    Int_t nFirstRun = 295585;
    Int_t nLastRun = 297595;
    Int_t nRuns = nLastRun - nFirstRun + 1;

    hCounterTrigger = new TH1F("hCounterTrigger", "# of events per run passing central triggers", nRuns, nFirstRun-0.5, nLastRun+0.5);
    fOutputList->Add(hCounterTrigger);

    hVertexZ = new TH1F("hVertexZ","hVertexZ",600,-30,30);
    fOutputList->Add(hVertexZ);

    hVertexContrib = new TH1F("hVertexContrib","hVertexContrib",103,-2,100);
    fOutputList->Add(hVertexContrib);

    hADdecision = new TH2I("hADdecision","hADdecision",7,-2,5,7,-2,5);
    fOutputList->Add(hADdecision);

    hV0decision = new TH2I("hV0decision","hV0decision",7,-2,5,7,-2,5);
    fOutputList->Add(hV0decision);

    PostData(2, fOutputList); 

}
//_____________________________________________________________________________
void AliAnalysisTaskJPsi_DG::TrkTrkKinematics(Int_t *fIndicesOfGoodTrks, Double_t fTrkMass)
{
  // Get the first track
  TLorentzVector fTrk1LorVec;
  AliESDtrack *trk1 = dynamic_cast<AliESDtrack*>(fEvent->GetTrack(fIndicesOfGoodTrks[0]));
  // Fill its 4-vector in the form: pt, eta, phi, mass
  fTrk1LorVec.SetPtEtaPhiM(trk1->Pt(), trk1->Eta(), trk1->Phi(), fTrkMass);
  // Get the second track
  TLorentzVector fTrk2LorVec;
  AliESDtrack *trk2 = dynamic_cast<AliESDtrack*>(fEvent->GetTrack(fIndicesOfGoodTrks[1]));
  // Fill its 4-vector
  fTrk2LorVec.SetPtEtaPhiM(trk2->Pt(), trk2->Eta(), trk2->Phi(), fTrkMass);
  // Vector of Trk+Trk: we add up the two
  TLorentzVector TrkTrk = fTrk1LorVec + fTrk2LorVec;

  // Set tree variables
  fPt = TrkTrk.Pt(); 
  fPhi = TrkTrk.Phi();
  fY = TrkTrk.Rapidity(); 
  fM = TrkTrk.M();
  fPt1 = trk1->Pt(); 
  fPt2 = trk2->Pt();
  fEta1 = trk1->Eta(); 
  fEta2 = trk2->Eta();
  fPhi1 = trk1->Phi();
  fPhi2 = trk2->Phi();
  fQ1 = trk1->Charge(); 
  fQ2 = trk2->Charge();
}
//_____________________________________________________________________________
void AliAnalysisTaskJPsi_DG::UserExec(Option_t *)
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
            PostData(1, fTreeJPsi);
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

    // Forward detectors
    // Store the data to the analysis tree
    // and fill the 2D histograms with V0 and AD responses
    
    // AD
    AliVAD *fADdata = fEvent->GetADData();

    fADA_dec = fADdata->GetADADecision();
    fADC_dec = fADdata->GetADCDecision();
    fADA_time = fADdata->GetADATime();
    fADC_time = fADdata->GetADCTime();
    hADdecision->Fill(fADA_dec,fADC_dec);

    // V0
    AliVVZERO *fV0data = fEvent->GetVZEROData();

    fV0A_dec = fV0data->GetV0ADecision();
    fV0C_dec = fV0data->GetV0CDecision();
    fV0A_time = fV0data->GetV0ATime();
    fV0C_time = fV0data->GetV0CTime();
    hV0decision->Fill(fV0A_dec,fV0C_dec);
    
    // ##########################################################
        // CUT 1 & 2
        // Check if each event has at maximum 1 vertex within 15 cm from the IP in beam direction
        const AliVVertex *fVertex = fEvent->GetPrimaryVertex();

        // Fill the trees and histograms
        fVertexZ = fVertex->GetZ();
        hVertexZ->Fill(fVertexZ);
        fVertexContrib = fVertex->GetNContributors();
        hVertexContrib->Fill(fVertexContrib);

        // At least two tracks associated with the vertex
        if(fVertex->GetNContributors()<2)
        {                                          
            PostData(1, fTreeJPsi);
            PostData(2, fOutputList);
            return;
        }
        hCounterCuts->Fill(iSelectionCounter);
        iSelectionCounter++;

        // Distance from the IP lower than 15 cm
        if(TMath::Abs(fVertex->GetZ())>15)
        {                                          
            PostData(1, fTreeJPsi);
            PostData(2, fOutputList);
            return;
        }
        hCounterCuts->Fill(iSelectionCounter);
        iSelectionCounter++;
    // ##########################################################
        // CUT 3 
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
            if(fTrackCutsBit4->AcceptTrack(trk)){
                nGoodTracksTPC++;
            }
            
            // Count good SPD tracks:
            if(fTrackCutsBit4->AcceptTrack(trk) && trk->HasPointOnITSLayer(0) && trk->HasPointOnITSLayer(1)){
                // If the track satisfies both the SPD and TPC criterion, add it to the list:
                fIndicesOfGoodTracks[nGoodTracksSPD] = iTrk;
                nGoodTracksSPD++;  
            }
        }
        // Continue only if two good central tracks are found
        if(!(nGoodTracksSPD == 2 && nGoodTracksTPC == 2)){                                          
            PostData(1, fTreeJPsi);
            PostData(2, fOutputList);
            delete [] fIndicesOfGoodTracks; 
            return;
        }
        hCounterCuts->Fill(iSelectionCounter);
        iSelectionCounter++;
    // ##########################################################
        // CUT 4
        // Central UPC trigger CCUP31
        if(!triggered)
        {
            PostData(1, fTreeJPsi);
            PostData(2, fOutputList);
            delete [] fIndicesOfGoodTracks; 
            return;        
        }
        hCounterCuts->Fill(iSelectionCounter);
        iSelectionCounter++;
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

    if(isMuonPair < isElectronPair) TrkTrkKinematics(fIndicesOfGoodTracks, massMuon);
    else TrkTrkKinematics(fIndicesOfGoodTracks, massElectron);
    
    // Check if SPD cluster matches FOhits (according to MB's macro)
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
    // Store the data to the analysis tree
  	fZNA_energy = fZDCdata->GetZNATowerEnergy()[0];
  	fZNC_energy = fZDCdata->GetZNCTowerEnergy()[0];

	Int_t detChZNA  = fZDCdata->GetZNATDCChannel();
    Int_t detChZNC  = fZDCdata->GetZNCTDCChannel();
	if (fEvent->GetRunNumber() >= 245726 && fEvent->GetRunNumber() <= 245793) detChZNA = 10;
  	for (Int_t i=0;i<4;i++){ 
  		fZNA_time[i] = fZDCdata->GetZDCTDCCorrected(detChZNA,i);
  		fZNC_time[i] = fZDCdata->GetZDCTDCCorrected(detChZNC,i);
	}

    // Fill the analysis tree
    fTreeJPsi->Fill();

    // Finally post the data
    PostData(1, fTreeJPsi);
    PostData(2, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskJPsi_DG::SetCrossed(Int_t spd[4], TBits &crossed){ 
    // from the macro by MB
    Int_t chipId2;
    for(Int_t iLayer = 0; iLayer < 4; iLayer++)
    if(spd[iLayer] > 0){
        crossed.SetBitNumber(GetChipId(spd[iLayer],chipId2)); 
        crossed.SetBitNumber(chipId2); 
    }
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskJPsi_DG::GetChipId(Int_t index, Int_t &chipId2, Bool_t debug){
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
Bool_t AliAnalysisTaskJPsi_DG::IsSTGFired(TBits bits, Int_t dphiMin, Int_t dphiMax, Bool_t tolerance){
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
void AliAnalysisTaskJPsi_DG::Terminate(Option_t *)
{
    // the end
}
//_____________________________________________________________________________
