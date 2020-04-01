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

#include "AliAnalysisTaskDeuteronAbsorption.h"

#include <iostream>

#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>

#include <AliAnalysisManager.h>
#include <AliESDEvent.h>
#include <AliESDInputHandler.h>
#include <AliESDtrack.h>
#include <AliESDtrackCuts.h>
#include <AliMCEvent.h>
#include <AliMCEventHandler.h>
#include <AliPIDResponse.h>

class AliAnalysisTaskDeuteronAbsorption;

const AliPID::EParticleType AliAnalysisTaskDeuteronAbsorption::fgkSpecies[kNabsSpecies] = {AliPID::kKaon, AliPID::kProton, AliPID::kDeuteron, AliPID::kTriton, AliPID::kHe3};
const std::string AliAnalysisTaskDeuteronAbsorption::fgkParticleNames[kNabsSpecies] = {"Kaon", "Proton", "Deuteron", "Triton","3He"};
const double AliAnalysisTaskDeuteronAbsorption::fgkPhiParamPos[4][4] = {
      {1.38984e+00, -2.10187e+01, 5.81724e-02, 1.91938e+01},
      {2.02372e+00, -2.44456e+00, 8.99000e-01, 9.22399e-01},
      {4.21954e+00, -2.56555e+01, 4.17557e-02, 2.40301e+01},
      {5.17499e+00, -2.69241e+00, 6.97167e-01, 1.25974e+00}};
const double AliAnalysisTaskDeuteronAbsorption::fgkPhiParamNeg[4][4] = {
      {2.81984e+00, -1.81497e-01, -2.03494e+00, 2.64148e-01},
      {5.79322e+00, -5.44966e-02, -1.10803e+00, 1.29737e+00},
      {5.60000e+00, -2.06000e-01, -1.97130e+00, 2.67181e-01},
      {9.72180e+00, -4.35801e-02, -1.14550e+00, 1.49160e+00}};

ClassImp(AliAnalysisTaskDeuteronAbsorption); // classimp: necessary for root

AliAnalysisTaskDeuteronAbsorption::AliAnalysisTaskDeuteronAbsorption(const char *name) : AliAnalysisTaskSE(name),
                                                                                         fUseTRDboundariesCut{true},
                                                                                         fNtpcSigmas{5.},
                                                                                         fEventCuts{},
                                                                                         fMindEdx{100.},
                                                                                         fMinTPCsignalN{50},
                                                                                         fPIDResponse{nullptr},
                                                                                         fESDtrackCuts{*AliESDtrackCuts::GetStandardITSTPCTrackCuts2011()},
                                                                                         fOutputList{nullptr},
                                                                                         fTreeTrack{nullptr},
                                                                                         tPt{-999.},
                                                                                         tEta{-999.},
                                                                                         tPhi{-999.},
                                                                                         tnsigTPC{-999.},
                                                                                         tnsigTOF{-999.},
                                                                                         tmass2{-999.},
                                                                                         tnPIDclsTPC{0},
                                                                                         tTOFsigDx{-999.},
                                                                                         tTOFsigDz{-999.},
                                                                                         tTOFclsN{0},
                                                                                         tID{0},
                                                                                         fHistZv{nullptr},
                                                                                         fHist3TPCpid{nullptr},
                                                                                         fHist3TPCpidAll{nullptr},
                                                                                         fHist3TOFpid{nullptr},
                                                                                         fHist3TOFpidAll{nullptr},
                                                                                         fHist3TOFmass{nullptr},
                                                                                         fHist3TOFmassAll{nullptr},
                                                                                         fHist1AcceptanceAll{nullptr},
                                                                                         fHist2Matching{nullptr},
                                                                                         fHist2Phi{nullptr},
                                                                                         fHist2TPCnSigma{nullptr},
                                                                                         fHist2MatchingMC{nullptr},
                                                                                         fTRDboundariesPos{nullptr},
                                                                                         fTRDboundariesNeg{nullptr}
{
  fESDtrackCuts.SetEtaRange(-0.8, 0.8);

  // constructor
  DefineInput(0, TChain::Class()); // define the input of the analysis: in this case we take a 'chain' of events
                                   // this chain is created by the analysis manager, so no need to worry about it,
                                   // it does its work automatically
  DefineOutput(1, TList::Class()); // define the ouptut of the analysis: in this case it's a list of histograms
                                   // you can add more output objects by calling DefineOutput(2, classname::Class())
                                   // if you add more output objects, make sure to call PostData for all of them, and to
//  if (fTreemode)
  DefineOutput(2, TTree::Class());
}

AliAnalysisTaskDeuteronAbsorption::~AliAnalysisTaskDeuteronAbsorption()
{
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  for (int iFunction = 0; iFunction < 4; ++iFunction)
  {
    if (fTRDboundariesPos[iFunction])
      delete fTRDboundariesPos[iFunction];
    if (fTRDboundariesNeg[iFunction])
      delete fTRDboundariesNeg[iFunction];
  }

  if (fOutputList)
    delete fOutputList; // at the end of your task, it is deleted from memory by calling this function

  if (fTreeTrack)
    delete fTreeTrack;
}

void AliAnalysisTaskDeuteronAbsorption::UserCreateOutputObjects()
{
  // create output objects
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if (man)
  {
    AliInputEventHandler *inputHandler = (AliInputEventHandler *)(man->GetInputEventHandler());
    if (inputHandler)
      fPIDResponse = inputHandler->GetPIDResponse();
  }

  // histograms used in the analysis
  // to an output file
  fOutputList = new TList();    // this is a list which will contain all of your histograms
  fOutputList->SetOwner(kTRUE); // memory stuff: the list is owner of all objects it contains and will delete them
  //
  fHistZv = new TH1F("fHistZv", "fHistZv", 200, -40, 40); // histogram to monitor z-position of the primary vertex -- quality assurance
  fOutputList->Add(fHistZv);
  //

  std::string wTRD[2] = {"woTRD", "wTRD"};
  std::string wTOF[2] = {"woTOF", "wTOF"};
  std::string pos_neg[2] = {"neg", "pos"};

  for (int iSpecies = 0; iSpecies < kNabsSpecies; ++iSpecies)
  {
    fHist3TPCpid[iSpecies] = new TH3F(Form("fHist3TPCpid%s", fgkParticleNames[iSpecies].data()), Form("%s; #it{p}/#it{z} (Gev/#it{c}); d#it{E}/d#it{x} (arb. units); #Phi (rad)", fgkParticleNames[iSpecies].data()), 400, -10, 10, 400, 0, 1000, 18, 0, TMath::TwoPi());
    fHist3TOFpid[iSpecies] = new TH3F(Form("fHist3TOFpid%s", fgkParticleNames[iSpecies].data()), Form("%s; #it{p}/#it{z} (Gev/#it{c}); #beta; #Phi (rad)", fgkParticleNames[iSpecies].data()), 400, -10, 10, 300, 0, 1.2, 18, 0, 2 * TMath::Pi());
    fHist3TOFmass[iSpecies] = new TH3F(Form("fHist3TOFmass%s", fgkParticleNames[iSpecies].data()), Form("%s; #it{p}/#it{z} (Gev/#it{c}); TOF m^{2} (GeV/#it{c}^{2})^{2}; #Phi (rad)", fgkParticleNames[iSpecies].data()), 400, -10, 10, 160, 0, 6.5, 18, 0, 2 * TMath::Pi());
    fOutputList->Add(fHist3TPCpid[iSpecies]);
    fOutputList->Add(fHist3TOFpid[iSpecies]);
    fOutputList->Add(fHist3TOFmass[iSpecies]);
  }
  fHist3TPCpidAll = new TH3F("fHist3TPCpidAll", "; #it{p}/#it{z} (Gev/#it{c}); d#it{E}/d#it{x} (arb. units); #Phi (rad)", 400, -10, 10, 1000, 0, 1000, 18, 0, TMath::TwoPi());
  fHist3TOFpidAll = new TH3F("fHist3TOFpidAll", "; #it{p}/#it{z} (Gev/#it{c}); #beta; #Phi (rad)", 400, -10, 10, 300, 0, 1.2, 18, 0, 2 * TMath::Pi());
  fHist3TOFmassAll = new TH3F("fHist3TOFmassAll", "; #it{p}/#it{z} (Gev/#it{c}); TOF m^{2} (GeV/#it{c}^{2})^{2}; #Phi (rad)", 400, -10, 10, 160, 0, 6.5, 18, 0, 2 * TMath::Pi());
  fOutputList->Add(fHist3TPCpidAll);
  fOutputList->Add(fHist3TOFpidAll);
  fOutputList->Add(fHist3TOFmassAll);

  for (int iCharge = 0; iCharge < 2; ++iCharge)
  {
    for (int iTRD = 0; iTRD < 2; ++iTRD)
    {
      fHist2Phi[iCharge][iTRD] = new TH2F(Form("fHist2Phi_%s_%s", pos_neg[iCharge].data(), wTRD[iTRD].data()), Form("%s %s; #Phi (rad) ;#it{p}_{T} (Gev/#it{c});", pos_neg[iCharge].data(), wTRD[iTRD].data()), 100, 0, 2 * TMath::Pi(), 100, 0, 7);
      fOutputList->Add(fHist2Phi[iCharge][iTRD]);
      for (int iTOF = 0; iTOF < 2; ++iTOF) {
        fHist1AcceptanceAll[iCharge][iTRD][iTOF] = new TH1F(Form("fHist1AcceptanceAll%s_%s_%s", pos_neg[iCharge].data(), wTRD[iTRD].data(), wTOF[iTOF].data()), Form("%s %s %s; #it{p}/#it{z} (Gev/#it{c});", pos_neg[iCharge].data(), wTRD[iTRD].data(), wTOF[iTOF].data()), 200, 0, 10);
        fOutputList->Add(fHist1AcceptanceAll[iCharge][iTRD][iTOF]);
      }
      for (int iSpecies = 0; iSpecies < kNabsSpecies; ++iSpecies)
      {
        fHist2Matching[iSpecies][iCharge][iTRD] = new TH2F(Form("fHist2Matching%s_%s_%s", fgkParticleNames[iSpecies].data(), pos_neg[iCharge].data(), wTRD[iTRD].data()), Form("%s %s %s; #it{p}/#it{z} (Gev/#it{c}); TOF m^{2} (GeV/#it{c}^{2})^{2}", fgkParticleNames[iSpecies].data(), pos_neg[iCharge].data(), wTRD[iTRD].data()), 200, 0, 10, 160, 0, 6.5);
        fOutputList->Add(fHist2Matching[iSpecies][iCharge][iTRD]);
        fHist2MatchingMC[iSpecies][iCharge][iTRD] = new TH2F(Form("fHist2MatchingMC%s_%s_%s", fgkParticleNames[iSpecies].data(), pos_neg[iCharge].data(), wTRD[iTRD].data()), Form("%s %s %s; #it{p}/#it{z} (Gev/#it{c}); TOF m^{2} (GeV/#it{c}^{2})^{2}", fgkParticleNames[iSpecies].data(), pos_neg[iCharge].data(), wTRD[iTRD].data()), 200, 0, 10, 160, 0, 6.5);
        fOutputList->Add(fHist2MatchingMC[iSpecies][iCharge][iTRD]);
        fHist2TPCnSigma[iSpecies][iCharge][iTRD] = new TH2F(Form("fHist2TPCnSigma%s_%s_%s", fgkParticleNames[iSpecies].data(), pos_neg[iCharge].data(), wTRD[iTRD].data()), Form("%s %s %s; #it{p}/#it{z} (Gev/#it{c}); TPC n_{#sigma}", fgkParticleNames[iSpecies].data(), pos_neg[iCharge].data(), wTRD[iTRD].data()), 200, 0, 10, 40, -5, 5);
        fOutputList->Add(fHist2TPCnSigma[iSpecies][iCharge][iTRD]);
      }
    }
  }

  // Tree
  if (fTreemode)
  {
    OpenFile(2);
    fTreeTrack = new TTree("fTreeTrack", "Track Parameters");
    //fTreeTrack->Branch("tP", &tP, "tP/D");
    fTreeTrack->Branch("tPt", &tPt, "tPt/D");
    fTreeTrack->Branch("tEta", &tEta, "tEta/D");
    fTreeTrack->Branch("tPhi", &tPhi, "tPhi/D");
    fTreeTrack->Branch("tnsigTPC", &tnsigTPC, "tnsigTPC/D");
    fTreeTrack->Branch("tnsigTOF", &tnsigTOF, "tnsigTOF/D");
    fTreeTrack->Branch("tmass2", &tmass2, "tmass2/D");
    fTreeTrack->Branch("tnPIDclsTPC", &tnPIDclsTPC, "tnPIDclsTPC/I");
    fTreeTrack->Branch("tTOFsigDx", &tTOFsigDx, "tTOFsigDx/D");
    fTreeTrack->Branch("tTOFsigDz", &tTOFsigDz, "tTOFsigDz/D");
    fTreeTrack->Branch("tTOFclsN", &tTOFclsN, "tTOFclsN/I");
    fTreeTrack->Branch("tID", &tID, "tID/I");
  }
  fEventCuts.AddQAplotsToList(fOutputList);

  PostData(1, fOutputList); // postdata will notify the analysis manager of changes / updates to the

  if (fTreemode)
    PostData(2, fTreeTrack);

  for (int iFunction = 0; iFunction < 4; ++iFunction)
  {
    fTRDboundariesNeg[iFunction] = new TF1(Form("fNeg%i", iFunction), "[0]-exp([1]*pow(x,[2])+[3])", 0.2, 10);
    fTRDboundariesPos[iFunction] = new TF1(Form("fPos%i", iFunction), "[0]-exp([1]*pow(x,[2])+[3])", 0.2, 10);
    for (int iParam = 0; iParam < 4; ++iParam)
    {
      fTRDboundariesNeg[iFunction]->SetParameter(iParam, fgkPhiParamNeg[iFunction][iParam]);
      fTRDboundariesPos[iFunction]->SetParameter(iParam, fgkPhiParamPos[iFunction][iParam]);
    }
  }
}

void AliAnalysisTaskDeuteronAbsorption::UserExec(Option_t *)
{
  // main loop over events
  AliESDEvent* esdEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!esdEvent)
    ::Fatal("AliAnalysisTaskDeuteronAbsorption::UserExec","No ESD event found.");  // if the pointer to the event is empty (getting it failed) skip this event
  Int_t nTracks = esdEvent->GetNumberOfTracks(); // see how many tracks there are in the event

  Bool_t isMC = false;
  AliMCEvent *mcEvent = nullptr;

  AliMCEventHandler * eventHandlerMC = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (eventHandlerMC)
  {
    mcEvent = eventHandlerMC->MCEvent();
    isMC = (mcEvent != nullptr);
  }

  if (!fEventCuts.AcceptEvent(esdEvent))
    return;

  // check for a proper primary vertex and monitor
  const AliESDVertex *vertex = esdEvent->GetPrimaryVertexTracks();
  if (vertex->GetNContributors() < 1)
  {
    // SPD vertex
    vertex = esdEvent->GetPrimaryVertexSPD();
    if (vertex->GetNContributors() < 1)
      vertex = nullptr;
  }
  if (vertex == nullptr) {
    PostData(1, fOutputList);
    return;
  }

  fHistZv->Fill(vertex->GetZ());
  if (TMath::Abs(vertex->GetZ()) > 10.0) {
    PostData(1, fOutputList);
    return; // remove events with a vertex which is more than 10cm away
  }

  // track loop
  for (Int_t i = 0; i < nTracks; i++)
  {                                                                     // loop ove rall these tracks
    AliESDtrack *track = static_cast<AliESDtrack *>(esdEvent->GetTrack(i)); // get a track (type AliESDDTrack) from the event
    if (!track)
      continue;
    if (!fESDtrackCuts.AcceptTrack(track))
      continue; // check if track passes the cuts
    if (!track->GetInnerParam())
      continue;                                     // check if track is a proper TPC track
    if (track->GetTPCsignalN() < fMinTPCsignalN)
      continue;

    // Process TOF information
    ULong_t status = (ULong_t)track->GetStatus();
    bool hasTOFout  = status & AliVTrack::kTOFout;
    bool hasTOFtime = status & AliVTrack::kTIME;
    const float length = track->GetIntegratedLength();
    bool hasTOF = hasTOFout && hasTOFtime && (length > 350.);
    //
    double ptot = track->GetTPCmomentum(); // momentum for dEdx determination
    double tof = track->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetStartTime(track->P());
    //
    double beta = -1.;
    double mass2 = -1;
    //
    if (hasTOF)
    {
      beta = length / (TMath::C() * 1.e-10 * tof);
      if ((1 - beta * beta) > 0)
        mass2 = ptot * ptot * (1. / (beta * beta) - 1.);
    }

    if (fTreemode && track->GetTPCsignal() > fMindEdx && std::abs(fPIDResponse->NumberOfSigmasTPC(track, fgkSpecies[4])) < 6)
    {
      //tP = track->GetInnerParam()->GetP();
      tPt = track->GetInnerParam()->GetSignedPt();
      tEta = track->GetInnerParam()->Eta();
      tPhi = track->GetInnerParam()->Phi();
      tmass2 = mass2;
      tnPIDclsTPC = track->GetTPCsignalN();
      tTOFsigDx = track->GetTOFsignalDx();
      tTOFsigDz = track->GetTOFsignalDz();
      tTOFclsN = track->GetTOFclusterN();
      tID = track->GetID();
      tnsigTPC = fPIDResponse->NumberOfSigmasTPC(track, fgkSpecies[4]);
      tnsigTOF = fPIDResponse->NumberOfSigmasTOF(track, fgkSpecies[4]);
      fTreeTrack->Fill();
    }

    //
    double sign = track->GetSign();
    // fill QA histograms
    double tpcNsigmas[kNabsSpecies];
    for (int iSpecies = 0; iSpecies < kNabsSpecies; ++iSpecies) tpcNsigmas[iSpecies] = 999.;
    //
    fHist3TPCpidAll->Fill(ptot * sign, track->GetTPCsignal(), track->Phi());
    if (hasTOF && track->GetTPCsignal() > fMindEdx)
    {
      fHist3TOFpidAll->Fill(ptot * sign, beta, track->Phi());
      fHist3TOFmassAll->Fill(ptot * sign, mass2, track->Phi());
    }
    for (int iSpecies = 0; iSpecies < kNabsSpecies; ++iSpecies)
    {
      tpcNsigmas[iSpecies] = fPIDResponse->NumberOfSigmasTPC(track, fgkSpecies[iSpecies]);
      if (std::abs(tpcNsigmas[iSpecies]) < fNtpcSigmas)
      {
        fHist3TPCpid[iSpecies]->Fill(ptot * sign, track->GetTPCsignal(), track->Phi());
        if (hasTOF && track->GetTPCsignal() > fMindEdx)
        {
          fHist3TOFpid[iSpecies]->Fill(ptot * sign, beta, track->Phi());
          fHist3TOFmass[iSpecies]->Fill(ptot * sign, mass2, track->Phi());
        }
      }
    }

    // study using TRDin
    bool hasTRDin = bool(status & AliESDtrack::kTRDin); // 2D phi pt for TRD

    double pt = track->Pt();
    double phi = track->Phi();
    while (phi < 0)
      phi += TMath::TwoPi();
    while (phi > TMath::TwoPi())
      phi -= TMath::TwoPi();
    bool withTRD[2]{
        fUseTRDboundariesCut ? true :
        phi < fTRDboundariesNeg[0]->Eval(pt) ||
            (phi > fTRDboundariesNeg[1]->Eval(pt) && phi < fTRDboundariesNeg[2]->Eval(pt)) ||
            phi > fTRDboundariesNeg[3]->Eval(pt),
        fUseTRDboundariesCut ? true :
        phi < fTRDboundariesPos[0]->Eval(pt) ||
            (phi > fTRDboundariesPos[1]->Eval(pt) && phi < fTRDboundariesPos[2]->Eval(pt)) ||
            phi > fTRDboundariesPos[3]->Eval(pt)};
    bool positive = sign > 0;
    fHist1AcceptanceAll[positive][withTRD[positive]][hasTOF]->Fill(ptot);
    for (int iSpecies = 0; iSpecies < kNabsSpecies; ++iSpecies)
    {
      if (track->GetTPCsignal() > fMindEdx) {
        fHist2TPCnSigma[iSpecies][positive][withTRD[positive]]->Fill(ptot, tpcNsigmas[iSpecies]);
      }
      if (std::abs(tpcNsigmas[iSpecies]) < 3.)
      {
        if (track->GetTPCsignal() > fMindEdx) {
          fHist2Matching[iSpecies][positive][withTRD[positive]]->Fill(ptot, mass2);
        }
        if (isMC)
        {
          int tofLabel[3]{-1,-1,-1};
          track->GetTOFLabel(tofLabel);
          int label = TMath::Abs(track->GetLabel());
          bool trueMatch = false;
          for (int iLabel = 0; iLabel < 3; ++iLabel)
          {
            if (tofLabel[iLabel] == label)
              trueMatch = true;
          }
          double mcMass = trueMatch ? AliPID::ParticleMassZ(fgkSpecies[iSpecies]) : 0;
          AliVParticle *mcpart = mcEvent->GetTrack(TMath::Abs(track->GetLabel()));
          if (mcpart)
          {
            int pdg = TMath::Abs(mcpart->PdgCode());
            if (pdg == AliPID::ParticleCode(fgkSpecies[iSpecies]))
            {
              fHist2MatchingMC[iSpecies][positive][withTRD[positive]]->Fill(ptot, mcMass * mcMass);
            }
          }
        }
      }
    }
    fHist2Phi[positive][hasTRDin]->Fill(phi, pt);

  } // end the track loop

  // post the data
  PostData(1, fOutputList);
  PostData(2, fTreeTrack);
} // end the UserExec
