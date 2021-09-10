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
#include <string.h>

// Root headers
#include "TROOT.h"
#include "TH1I.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TColor.h"
#include "TRandom.h"
#include "TList.h"
#include "TObjArray.h"
#include "TString.h"
#include "TChain.h"
#include "TMath.h"
#include "TParticle.h"

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
#include "AliAnalysisTaskJPsiMC_DG.h"

class AliAnalysisTaskJPsiMC_DG; // your analysis class

//using namespace std; // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskJPsiMC_DG) // classimp: necessary for root

AliAnalysisTaskJPsiMC_DG::AliAnalysisTaskJPsiMC_DG() : // initializer list
    AliAnalysisTaskSE(),
    fPIDResponse(0),
    fTrackCutsBit4(0),
    isNeutralPions(kFALSE),
    fEvent(0),
    fOutputList(0),
    fTreeJPsiMCRec(0),
    fTreeJPsiMCGen(0),
    fRunNumber(0),
    // Histograms:
    hCounterCuts(0),
    hPtRecGen(0),
    // PID, sigmas:
    fTrk1SigIfMu(0),
    fTrk1SigIfEl(0),
    fTrk2SigIfMu(0),
    fTrk2SigIfEl(0),
    // Kinematics:
    fPt(0), fPhi(0), fY(0), fM(0),
    // Two tracks:
    fPt1(0), fPt2(0), fEta1(0), fEta2(0), fPhi1(0), fPhi2(0), fQ1(0), fQ2(0),
    // Info from the detectors:
    // ZDC
    fZNA_energy(0), fZNC_energy(0),
    // V0
    fV0A_dec(0), fV0C_dec(0), fV0A_time(0), fV0C_time(0),
    // AD
    fADA_dec(0), fADC_dec(0), fADA_time(0), fADC_time(0),
    // Matching SPD clusters with FOhits
    fMatchingSPD(0),
    // Trigger inputs for MC data
    fSPDfile(0), fTOFfile(0), fLoadedRun(-1), hTOFeff(0), hSPDeff(0), fTOFmask(0),
    // MC kinematics on generated level
    fPtGen(0), fMGen(0), fYGen(0), fPhiGen(0)
{
    // default constructor
}
//_____________________________________________________________________________
AliAnalysisTaskJPsiMC_DG::AliAnalysisTaskJPsiMC_DG(const char* name) : // initializer list
    AliAnalysisTaskSE(name),
    fPIDResponse(0),
    fTrackCutsBit4(0),
    isNeutralPions(kFALSE),
    fEvent(0),
    fOutputList(0),
    fTreeJPsiMCRec(0),
    fTreeJPsiMCGen(0),
    fRunNumber(0),
    // Histograms:
    hCounterCuts(0),
    hPtRecGen(0),
    // PID, sigmas:
    fTrk1SigIfMu(0),
    fTrk1SigIfEl(0),
    fTrk2SigIfMu(0),
    fTrk2SigIfEl(0),
    // Kinematics:
    fPt(0), fPhi(0), fY(0), fM(0),
    // Two tracks:
    fPt1(0), fPt2(0), fEta1(0), fEta2(0), fPhi1(0), fPhi2(0), fQ1(0), fQ2(0),
    // Info from the detectors:
    // ZDC
    fZNA_energy(0), fZNC_energy(0),
    // V0
    fV0A_dec(0), fV0C_dec(0), fV0A_time(0), fV0C_time(0),
    // AD
    fADA_dec(0), fADC_dec(0), fADA_time(0), fADC_time(0),
    // Matching SPD clusters with FOhits
    fMatchingSPD(0),
    // Trigger inputs for MC data
    fSPDfile(0), fTOFfile(0), fLoadedRun(-1), hTOFeff(0), hSPDeff(0), fTOFmask(0),
    // MC kinematics on generated level
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
AliAnalysisTaskJPsiMC_DG::~AliAnalysisTaskJPsiMC_DG()
{
    // destructor

    if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis){
        delete fOutputList;
        fOutputList = 0x0;
    }
    if(fOutputList) {delete fOutputList;}
    if(fTreeJPsiMCRec) {delete fTreeJPsiMCRec;}
    if(fTreeJPsiMCGen) {delete fTreeJPsiMCGen;} 
    if(fPIDResponse) {delete fPIDResponse;} 
    // delete other things as well ??
}
//_____________________________________________________________________________
void AliAnalysisTaskJPsiMC_DG::UserCreateOutputObjects()
{
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();

    fTrackCutsBit4 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
    
    // ##########################################################
    // OUTPUT TREE MC REC:

    fTreeJPsiMCRec = new TTree("fTreeJPsiMCRec", "fTreeJPsiMCRec");
    // Basic things:
    fTreeJPsiMCRec->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");
    // PID, sigmas:
    fTreeJPsiMCRec->Branch("fTrk1SigIfMu", &fTrk1SigIfMu, "fTrk1SigIfMu/D");
    fTreeJPsiMCRec->Branch("fTrk1SigIfEl", &fTrk1SigIfEl, "fTrk1SigIfEl/D");
    fTreeJPsiMCRec->Branch("fTrk2SigIfMu", &fTrk2SigIfMu, "fTrk2SigIfMu/D");
    fTreeJPsiMCRec->Branch("fTrk2SigIfEl", &fTrk2SigIfEl, "fTrk2SigIfEl/D");
    // Kinematics:
    fTreeJPsiMCRec->Branch("fPt", &fPt, "fPt/D");
    fTreeJPsiMCRec->Branch("fM", &fM, "fM/D");
    fTreeJPsiMCRec->Branch("fY", &fY, "fY/D");
    fTreeJPsiMCRec->Branch("fPhi", &fPhi, "fPhi/D");
    // Two tracks:
    fTreeJPsiMCRec->Branch("fPt1", &fPt1, "fPt1/D");
    fTreeJPsiMCRec->Branch("fPt2", &fPt2, "fPt2/D");
    fTreeJPsiMCRec->Branch("fEta1", &fEta1, "fEta1/D");
    fTreeJPsiMCRec->Branch("fEta2", &fEta2, "fEta2/D");
    fTreeJPsiMCRec->Branch("fPhi1", &fPhi1, "fPhi1/D");
    fTreeJPsiMCRec->Branch("fPhi2", &fPhi2, "fPhi2/D");
    fTreeJPsiMCRec->Branch("fQ1", &fQ1, "fQ1/D");
    fTreeJPsiMCRec->Branch("fQ2", &fQ2, "fQ2/D");
    // Info from the detectors:
    // ZDC:
    fTreeJPsiMCRec->Branch("fZNA_energy", &fZNA_energy, "fZNA_energy/D");
    fTreeJPsiMCRec->Branch("fZNC_energy", &fZNC_energy, "fZNC_energy/D");
    fTreeJPsiMCRec->Branch("fZNA_time", &fZNA_time[0], "fZNA_time[4]/D");
    fTreeJPsiMCRec->Branch("fZNC_time", &fZNC_time[0], "fZNC_time[4]/D");
    // V0:
    fTreeJPsiMCRec->Branch("fV0A_dec", &fV0A_dec, "fV0A_dec/I");
    fTreeJPsiMCRec->Branch("fV0C_dec", &fV0C_dec, "fV0C_dec/I");
    fTreeJPsiMCRec->Branch("fV0A_time", &fV0A_time, "fV0A_time/D");
    fTreeJPsiMCRec->Branch("fV0C_time", &fV0C_time, "fV0C_time/D");
    // AD:
    fTreeJPsiMCRec->Branch("fADA_dec", &fADA_dec, "fADA_dec/I");
    fTreeJPsiMCRec->Branch("fADC_dec", &fADC_dec, "fADC_dec/I");
    fTreeJPsiMCRec->Branch("fADA_time", &fADA_time, "fADA_time/D");
    fTreeJPsiMCRec->Branch("fADC_time", &fADC_time, "fADC_time/D");
    // Matching SPD clusters with FOhits:
    fTreeJPsiMCRec->Branch("fMatchingSPD", &fMatchingSPD, "fMatchingSPD/O"); // O is for bool
    // Replayed trigger inputs:
    fTreeJPsiMCRec->Branch("fTriggerInputsMC", &fTriggerInputsMC[0], "fTriggerInputsMC[11]/O"); 
    // Kinematics, MC gen:
    fTreeJPsiMCRec->Branch("fPtGen", &fPtGen, "fPtGen/D");
    fTreeJPsiMCRec->Branch("fMGen", &fMGen, "fMGen/D");
    fTreeJPsiMCRec->Branch("fYGen", &fYGen, "fYGen/D");
    fTreeJPsiMCRec->Branch("fPhiGen", &fPhiGen, "fPhiGen/D");    

    PostData(1, fTreeJPsiMCRec);

    // ##########################################################
    // OUTPUT TREE MC GEN:

    fTreeJPsiMCGen = new TTree("fTreeJPsiMCGen", "fTreeJPsiMCGen");
    // Run number:
    fTreeJPsiMCGen->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");
    // Kinematics:
    fTreeJPsiMCGen->Branch("fPtGen", &fPtGen, "fPtGen/D");
    fTreeJPsiMCGen->Branch("fMGen", &fMGen, "fMGen/D");
    fTreeJPsiMCGen->Branch("fYGen", &fYGen, "fYGen/D");
    fTreeJPsiMCGen->Branch("fPhiGen", &fPhiGen, "fPhiGen/D");

    PostData(2, fTreeJPsiMCGen);

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
    //hCounterCuts->GetXaxis()->SetBinLabel(5,"4: CCUP31 trigg");

    fOutputList->Add(hCounterCuts);

    Int_t n_bins = 2000;
    // x axis = pt generated, y axis = pt reconstructed
    hPtRecGen = new TH2F("hPtRecGen", "pt rec vs pt gen", n_bins, 0., 2., n_bins, 0., 2.);
    fOutputList->Add(hPtRecGen);

    PostData(3, fOutputList); 

}
//_____________________________________________________________________________
void AliAnalysisTaskJPsiMC_DG::SetNeutralPions(Bool_t Neutral)
{
    isNeutralPions = Neutral;
}
//_____________________________________________________________________________
void AliAnalysisTaskJPsiMC_DG::TrkTrkKinematics(Int_t *fIndicesOfGoodTrks, Double_t fTrkMass)
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
void AliAnalysisTaskJPsiMC_DG::UserExec(Option_t *)
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
            PostData(1, fTreeJPsiMCRec);
            PostData(2, fTreeJPsiMCGen);
            PostData(3, fOutputList);
            return;
        }                               
        hCounterCuts->Fill(iSelectionCounter);
        iSelectionCounter++;
    // ##########################################################

    // Get the run number associated to the ESD event
    fRunNumber = fEvent->GetRunNumber(); 

    // Time range cuts: https://twiki.cern.ch/twiki/bin/view/ALICE/AliDPGRunList18r 
    fTimeRangeCut.InitFromEvent(InputEvent());
    if(fTimeRangeCut.CutEvent(InputEvent())) return;

    // Replay triggers
    if(fRunNumber != fLoadedRun){
      if(!fSPDfile) fSPDfile = AliDataFile::OpenOADB("PWGUD/UPC/SPDEfficiency18qr.root");
      if(!fTOFfile) fTOFfile = AliDataFile::OpenOADB("PWGUD/TOFTriggerEfficiency.root");
      AliOADBContainer* fTOFcont = (AliOADBContainer*)fTOFfile->Get("TOFTriggerEfficiency");
      hTOFeff  = (TH2F*)fTOFcont->GetObject(fRunNumber,"Default");
      
      if(fSPDfile->Get(Form("eff%i",fRunNumber))) hSPDeff  = (TH1D*) fSPDfile->Get(Form("eff%i",fRunNumber));
      Int_t tempRun = fRunNumber;
      while(!hSPDeff){
        tempRun--;
        hSPDeff = (TH1D*) fSPDfile->Get(Form("eff%i",tempRun));
      }
      fLoadedRun = fRunNumber;
      }
    ReplayTriggersMC(fEvent); 

    // Run analysis on the generated level
    RunMCGenerated();

    // ##########################################################
        // CUT 1 & 2
        // Check if each event has at maximum 1 vertex within 15 cm from the IP in beam direction
        const AliVVertex *fVertex = fEvent->GetPrimaryVertex();
        // At least two tracks associated with the vertex
        if(fVertex->GetNContributors()<2)
        {                                          
            PostData(1, fTreeJPsiMCRec);
            PostData(2, fTreeJPsiMCGen);
            PostData(3, fOutputList);
            return;
        }
        hCounterCuts->Fill(iSelectionCounter);
        iSelectionCounter++;
        // Distance from the IP lower than 15 cm
        if(TMath::Abs(fVertex->GetZ())>15)
        {                                          
            PostData(1, fTreeJPsiMCRec);
            PostData(2, fTreeJPsiMCGen);
            PostData(3, fOutputList);
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

            // *********************************************************************
            // FD correction for neutral pions: skip pion tracks to simulate FD with neutral pions
            if(isNeutralPions){
                if(trk->GetLabel() >= 0){
                    TDatabasePDG *pdgdat = TDatabasePDG::Instance();
                    AliMCEvent *mc = MCEvent();
                    if(!mc) return;
                    AliMCParticle *mcPart = (AliMCParticle*) mc->GetTrack(trk->GetLabel());
                    if(TMath::Abs(mcPart->PdgCode()) == 211) continue;
                }
            }
            // *********************************************************************

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
            PostData(1, fTreeJPsiMCRec);
            PostData(2, fTreeJPsiMCGen);
            PostData(3, fOutputList);
            delete [] fIndicesOfGoodTracks; 
            return;
        }
        hCounterCuts->Fill(iSelectionCounter);
        iSelectionCounter++;
    // ##########################################################
        // CUT 4
        // Central UPC trigger CCUP31
        // skipped for MC
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

    // Fill the 2D histogram of pt gen and pt rec
    hPtRecGen->Fill(fPtGen,fPt);

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
    // Data from the other detectors

    // AD
    AliVAD *fADdata = fEvent->GetADData();
    // Store the data to the analysis tree
    fADA_dec = fADdata->GetADADecision();
    fADC_dec = fADdata->GetADCDecision();
    fADA_time = fADdata->GetADATime();
    fADC_time = fADdata->GetADCTime();

    // V0
    AliVVZERO *fV0data = fEvent->GetVZEROData();
    // Store the data to the analysis tree
    fV0A_dec = fV0data->GetV0ADecision();
    fV0C_dec = fV0data->GetV0CDecision();
    fV0A_time = fV0data->GetV0ATime();
    fV0C_time = fV0data->GetV0CTime();

    // ZDC
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
    fTreeJPsiMCRec->Fill();

    // Finally post the data
    PostData(1, fTreeJPsiMCRec);
    PostData(2, fTreeJPsiMCGen);
    PostData(3, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskJPsiMC_DG::ReplayTriggersMC(AliVEvent *fEvent)
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
void AliAnalysisTaskJPsiMC_DG::RunMCGenerated()
{
    TLorentzVector vGenerated, vDecayProduct;
    TDatabasePDG *pdgdat = TDatabasePDG::Instance();

    vGenerated.SetPtEtaPhiM(0.,0.,0.,0.);

    AliMCEvent *mc = MCEvent();
    if(!mc){
        //Printf("Not found");
        return;
    } 

    for(Int_t imc = 0; imc < mc->GetNumberOfTracks(); imc++) {
        AliMCParticle *mcPart = (AliMCParticle*) mc->GetTrack(imc);
        if(!mcPart) continue;
    
        if(TMath::Abs(mcPart->PdgCode()) == 13){ 
            // if mu+ or mu-

            if(mcPart->GetMother() == -1){
                // without mother particle
                TParticlePDG *partGen = pdgdat->GetParticle(mcPart->PdgCode());
                vDecayProduct.SetXYZM(mcPart->Px(),mcPart->Py(), mcPart->Pz(),partGen->Mass());
                vGenerated += vDecayProduct;
            } else { // this branch not needed for kTwoGammaToMuMedium
                // with J/psi mother particle
                AliMCParticle *mcMother = (AliMCParticle*) mc->GetTrack(mcPart->GetMother());
                // Original code (manually selected):
                    // if(TMath::Abs(mcMother->PdgCode()) != 443) continue;     // for kCohJpsiToMu and kIncohJpsiToMu
                    // if(TMath::Abs(mcMother->PdgCode()) != 100443) continue;  // for kCohPsi2sToMuPi and kIncohPsi2sToMuPi
                // Previous two conditions merged (they cannot happen at the same time):
                    if(TMath::Abs(mcMother->PdgCode()) != 443 && TMath::Abs(mcMother->PdgCode()) != 100443) continue;
                TParticlePDG *partGen = pdgdat->GetParticle(mcPart->PdgCode());
                vDecayProduct.SetXYZM(mcPart->Px(),mcPart->Py(), mcPart->Pz(),partGen->Mass());
                vGenerated += vDecayProduct;
            }
        }
    } // loop over mc particles

    FillMCGenTree(vGenerated);
}
//_____________________________________________________________________________
void AliAnalysisTaskJPsiMC_DG::FillMCGenTree(TLorentzVector v)
{
    fPtGen      = v.Pt();
    if(v.E() != v.Pz()) fYGen = v.Rapidity();
    else fYGen  = -999; // when E = Pz, rapidity goes to infty
    fMGen       = v.M();
    fPhiGen     = v.Phi();

    fTreeJPsiMCGen->Fill();
}
//_____________________________________________________________________________
void AliAnalysisTaskJPsiMC_DG::SetCrossed(Int_t spd[4], TBits &crossed)
{ 
    // from the macro by MB
    Int_t chipId2;
    for(Int_t iLayer = 0; iLayer < 4; iLayer++)
    if(spd[iLayer] > 0){
        crossed.SetBitNumber(GetChipId(spd[iLayer],chipId2)); 
        crossed.SetBitNumber(chipId2); 
    }
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskJPsiMC_DG::GetChipId(Int_t index, Int_t &chipId2, Bool_t debug)
{
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
Bool_t AliAnalysisTaskJPsiMC_DG::IsSTGFired(TBits bits, Int_t dphiMin, Int_t dphiMax, Bool_t tolerance)
{
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
void AliAnalysisTaskJPsiMC_DG::Terminate(Option_t *)
{
    // the end
}
//_____________________________________________________________________________
