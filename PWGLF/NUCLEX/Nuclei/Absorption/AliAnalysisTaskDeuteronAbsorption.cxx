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

/* AliAnaysisTaskDeuteronAbsorption
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */

#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include <iostream>
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliPIDResponse.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisTaskDeuteronAbsorption.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"
#include "AliTRDCalDCS.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"

class AliAnalysisTaskDeuteronAbsorption;

const AliPID::EParticleType AliAnalysisTaskDeuteronAbsorption::fgkSpecies[4] = {AliPID::kKaon, AliPID::kProton, AliPID::kDeuteron, AliPID::kTriton};
const std::string AliAnalysisTaskDeuteronAbsorption::fgkParticleNames[4] = {"Kaon", "Proton", "Deuteron", "Triton"};

using namespace std; // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskDeuteronAbsorption); // classimp: necessary for root

//_____________________________________________________________________________
AliAnalysisTaskDeuteronAbsorption::AliAnalysisTaskDeuteronAbsorption(const char *name) : AliAnalysisTaskSE(name),
                                                                                         fMindEdx{100.},
                                                                                         fESD(0), fPIDResponse(0), fESDtrackCuts(0), fOutputList(0),
                                                                                         fHistZv{nullptr}, //
                                                                                         fHist3TPCpid{nullptr},
                                                                                         fHist3TPCpidAll{nullptr},
                                                                                         fHist3TOFpid{nullptr},
                                                                                         fHist3TOFpidAll{nullptr},
                                                                                         fHist3TOFmass{nullptr},
                                                                                         fHist3TOFmassAll{nullptr},
                                                                                         fHist2Matching{nullptr},
                                                                                         fHist2Phi{nullptr},
                                                                                         fHist2MatchingMC{nullptr},
                                                                                         fTRDboundariesPos{nullptr},
                                                                                         fTRDboundariesNeg{nullptr}
{

  //
  // constructor
  //
  DefineInput(0, TChain::Class()); // define the input of the analysis: in this case we take a 'chain' of events
                                   // this chain is created by the analysis manager, so no need to worry about it,
                                   // it does its work automatically
  DefineOutput(1, TList::Class()); // define the ouptut of the analysis: in this case it's a list of histograms
                                   // you can add more output objects by calling DefineOutput(2, classname::Class())
                                   // if you add more output objects, make sure to call PostData for all of them, and to
}

//_____________________________________________________________________________
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
  {
    delete fOutputList; // at the end of your task, it is deleted from memory by calling this function
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskDeuteronAbsorption::UserCreateOutputObjects()
{
  //
  // create output objects
  //
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if (man)
  {
    AliInputEventHandler *inputHandler = (AliInputEventHandler *)(man->GetInputEventHandler());
    if (inputHandler)
      fPIDResponse = inputHandler->GetPIDResponse();
  }
  //
  // histograms used in the analysis
  // to an output file
  //
  fOutputList = new TList();    // this is a list which will contain all of your histograms
  fOutputList->SetOwner(kTRUE); // memory stuff: the list is owner of all objects it contains and will delete them
  //
  fHistZv = new TH1F("fHistZv", "fHistZv", 200, -40, 40); // histogram to monitor z-position of the primary vertex -- quality assurance
  fOutputList->Add(fHistZv);
  //

  std::string TRDio[2] = {"woTRD", "wTRD"};
  std::string pos_neg[2] = {"neg", "pos"};

  for (int iSpecies = 0; iSpecies < 4; ++iSpecies)
  {
    fHist3TPCpid[iSpecies] = new TH3F(Form("fHist3TPCpid%s", fgkParticleNames[iSpecies].data()), Form("%s; #it{p}/#it{z} (Gev/#it{c}); d#it{E}/d#it{x} (arb. units); #Phi (rad)", fgkParticleNames[iSpecies].data()), 400, -10, 10, 400, 0, 1000, 18, 0, TMath::TwoPi());
    fHist3TOFpid[iSpecies] = new TH3F(Form("fHist3TOFpid%s", fgkParticleNames[iSpecies].data()), Form("%s; #it{p}/#it{z} (Gev/#it{c}); #beta; #Phi (rad)", fgkParticleNames[iSpecies].data()), 400, -10, 10, 300, 0, 1.2, 18, 0, 2 * TMath::Pi());
    fHist3TOFmass[iSpecies] = new TH3F(Form("fHist3TOFmass%s", fgkParticleNames[iSpecies].data()), Form("%s; #it{p}/#it{z} (Gev/#it{c}); TOF mass (GeV/#it{c}^{2}); #Phi (rad)", fgkParticleNames[iSpecies].data()), 400, -10, 10, 160, 0, 6.5, 18, 0, 2 * TMath::Pi());
    fOutputList->Add(fHist3TPCpid[iSpecies]);
    fOutputList->Add(fHist3TOFpid[iSpecies]);
    fOutputList->Add(fHist3TOFmass[iSpecies]);
  }
  fHist3TPCpidAll = new TH3F("fHist3TPCpidAll", "; #it{p}/#it{z} (Gev/#it{c}); d#it{E}/d#it{x} (arb. units); #Phi (rad)", 400, -10, 10, 1000, 0, 1000, 18, 0, TMath::TwoPi());
  fHist3TOFpidAll = new TH3F("fHist3TOFpidAll", "; #it{p}/#it{z} (Gev/#it{c}); #beta; #Phi (rad)", 400, -10, 10, 300, 0, 1.2, 18, 0, 2 * TMath::Pi());
  fHist3TOFmassAll = new TH3F("fHist3TOFmassAll", "; #it{p}/#it{z} (Gev/#it{c}); TOF mass (GeV/#it{c}^{2}); #Phi (rad)", 400, -10, 10, 160, 0, 6.5, 18, 0, 2 * TMath::Pi());
  fOutputList->Add(fHist3TPCpidAll);
  fOutputList->Add(fHist3TOFpidAll);
  fOutputList->Add(fHist3TOFmassAll);

  for (int iCharge = 0; iCharge < 2; ++iCharge)
  {
    for (int iTRD = 0; iTRD < 2; ++iTRD)
    {
      fHist2Phi[iCharge][iTRD] = new TH2F(Form("fHist2Phi_%s_%s", pos_neg[iCharge].data(), TRDio[iTRD].data()), Form("%s %s; #Phi (rad) ;#it{p}_{T} (Gev/#it{c});", pos_neg[iCharge].data(), TRDio[iTRD].data()), 100, 0, 2 * TMath::Pi(), 100, 0, 7);
      fOutputList->Add(fHist2Phi[iCharge][iTRD]);
      for (int iSpecies = 0; iSpecies < 4; ++iSpecies)
      {
        fHist2Matching[iSpecies][iCharge][iTRD] = new TH2F(Form("fHist2Matching%s_%s_%s", fgkParticleNames[iSpecies].data(), pos_neg[iCharge].data(), TRDio[iTRD].data()), Form("%s %s %s; #it{p}/#it{z} (Gev/#it{c}); TOF mass (GeV/#it{c}^{2})", fgkParticleNames[iSpecies].data(), pos_neg[iCharge].data(), TRDio[iTRD].data()), 400, -10, 10, 160, 0, 6.5);
        fOutputList->Add(fHist2Matching[iSpecies][iCharge][iTRD]);
        fHist2MatchingMC[iSpecies][iCharge][iTRD] = new TH2F(Form("fHist2MatchingMC%s_%s_%s", fgkParticleNames[iSpecies].data(), pos_neg[iCharge].data(), TRDio[iTRD].data()), Form("%s %s %s; #it{p}/#it{z} (Gev/#it{c}); TOF mass (GeV/#it{c}^{2})", fgkParticleNames[iSpecies].data(), pos_neg[iCharge].data(), TRDio[iTRD].data()), 400, -10, 10, 160, 0, 6.5);
        fOutputList->Add(fHist2MatchingMC[iSpecies][iCharge][iTRD]);
      }
    }
  }

  //
  // create track cuts object
  //
  fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE, kTRUE);
  // fESDtrackCuts->SetMaxDCAToVertexXY(3);
  // fESDtrackCuts->SetMaxDCAToVertexZ(2);
  fESDtrackCuts->SetEtaRange(-0.8, 0.8);

  PostData(1, fOutputList); // postdata will notify the analysis manager of changes / updates to the

  double paramsPos[4][4] = {
      {1.38984e+00, -2.10187e+01, 5.81724e-02, 1.91938e+01},
      {2.02372e+00, -2.44456e+00, 8.99000e-01, 9.22399e-01},
      {4.21954e+00, -2.56555e+01, 4.17557e-02, 2.40301e+01},
      {5.17499e+00, -2.69241e+00, 6.97167e-01, 1.25974e+00}};
  for (int iFunction = 0; iFunction < 4; ++iFunction)
  {
    fTRDboundariesPos[iFunction] = new TF1(Form("f%i", iFunction), "[0]-exp([1]*pow(x,[2])+[3])", 0.2, 10);
    for (int iParam = 0; iParam < 4; ++iParam)
    {
      fTRDboundariesPos[iFunction]->SetParameter(iParam, paramsPos[iFunction][iParam]);
    }
  }

  double paramsNeg[4][4] = {
      {2.81984e+00, -1.81497e-01, -2.03494e+00, 2.64148e-01},
      {5.79322e+00, -5.44966e-02, -1.10803e+00, 1.29737e+00},
      {5.60000e+00, -2.06000e-01, -1.97130e+00, 2.67181e-01},
      {9.72180e+00, -4.35801e-02, -1.14550e+00, 1.49160e+00}};
  for (int iFunction = 0; iFunction < 4; ++iFunction)
  {
    fTRDboundariesNeg[iFunction] = new TF1(Form("f%i", iFunction), "[0]-exp([1]*pow(x,[2])+[3])", 0.2, 10);
    for (int iParam = 0; iParam < 4; ++iParam)
    {
      fTRDboundariesNeg[iFunction]->SetParameter(iParam, paramsNeg[iFunction][iParam]);
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskDeuteronAbsorption::UserExec(Option_t *)
{
  //
  // main loop over events
  //
  fESD = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!fESD)
    ::Fatal("AliAnalysisTaskDeuteronAbsorption::UserExec","No ESD event found.");  // if the pointer to the event is empty (getting it failed) skip this event
  Int_t nTracks = fESD->GetNumberOfTracks(); // see how many tracks there are in the event

  Bool_t isMC = false;
  AliMCEvent *mcEvent = 0x0;

  AliMCEventHandler * eventHandlerMC = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (eventHandlerMC)
  {
    mcEvent = eventHandlerMC->MCEvent();
    isMC = true;
  }
  //
  // check for a proper primary vertex and monitor
  //
  const AliESDVertex *vertex = fESD->GetPrimaryVertexTracks();
  if (vertex->GetNContributors() < 1)
  {
    // SPD vertex
    vertex = fESD->GetPrimaryVertexSPD();
    if (vertex->GetNContributors() < 1)
      vertex = 0x0;
  }
  if (!vertex) {
    PostData(1, fOutputList);
    return;
  }

  fHistZv->Fill(vertex->GetZ());
  if (TMath::Abs(vertex->GetZ()) > 10.0) {
    PostData(1, fOutputList);
    return; // remove events with a vertex which is more than 10cm away
  }
  //
  // track loop
  //
  for (Int_t i = 0; i < nTracks; i++)
  {                                                                     // loop ove rall these tracks
    AliESDtrack *track = static_cast<AliESDtrack *>(fESD->GetTrack(i)); // get a track (type AliESDDTrack) from the event
    if (!track)
      continue;
    //
    if (!fESDtrackCuts->AcceptTrack(track))
      continue; // check if track passes the cuts
    if (!track->GetInnerParam())
      continue;                                     // check if track is a proper TPC track
    Double_t ptot = track->GetInnerParam()->GetP(); // momentum for dEdx determination
    if (track->GetTPCsignalN() < 50)
      continue;

    Double_t sign = track->GetSign();

    // Process TOF information
    ULong_t status = (ULong_t)track->GetStatus();
    Bool_t hasTOF = kFALSE;
    Bool_t hasTOFout = status & AliESDtrack::kTOFout;
    if (hasTOFout)
      hasTOF = kTRUE;
    Float_t length = track->GetIntegratedLength();
    if (length < 350.)
      hasTOF = kFALSE;
    //
    Double_t tof = track->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetStartTime(ptot);
    //
    Float_t beta = -1.;
    Float_t gamma = 0;
    Float_t mass = -99;
    //
    if (hasTOF)
    {
      beta = length / (TMath::C() * 1.e-10 * tof);
      if ((1 - beta * beta) > 0)
      {
        gamma = 1 / TMath::Sqrt(1 - beta * beta);
        mass = ptot / TMath::Sqrt(gamma * gamma - 1); // using inner TPC mom. as approx.}
      }
    }

    // fill QA histograms
    double tpcNsigmas[4]{999., 999., 999., 999.};
    fHist3TPCpidAll->Fill(ptot * sign, track->GetTPCsignal(), track->Phi());
    if (hasTOF && track->GetTPCsignal() > fMindEdx)
    {
      fHist3TOFpidAll->Fill(ptot * sign, beta, track->Phi());
      fHist3TOFmassAll->Fill(ptot * sign, mass * mass, track->Phi());
    }
    for (int iSpecies = 0; iSpecies < 4; ++iSpecies)
    {
      tpcNsigmas[iSpecies] = fPIDResponse->NumberOfSigmasTPC(track, fgkSpecies[iSpecies]);
      if (std::abs(tpcNsigmas[iSpecies]) < 3.)
      {
        fHist3TPCpid[iSpecies]->Fill(ptot * sign, track->GetTPCsignal(), track->Phi());
        if (hasTOF && track->GetTPCsignal() > fMindEdx)
        {
          fHist3TOFpid[iSpecies]->Fill(ptot * sign, beta, track->Phi());
          fHist3TOFmass[iSpecies]->Fill(ptot * sign, mass * mass, track->Phi());
        }
      }
    }

    // study using TRDin
    //
    bool hasTRDin = bool(status & AliESDtrack::kTRDin); // 2D phi pt for TRD

    float pt = track->Pt();
    float phi = track->Phi();
    while (phi < 0)
      phi += TMath::TwoPi();
    while (phi > TMath::TwoPi())
      phi -= TMath::TwoPi();
    bool withTRD[2]{
        phi < fTRDboundariesNeg[0]->Eval(pt) ||
            (phi > fTRDboundariesNeg[1]->Eval(pt) && phi < fTRDboundariesNeg[2]->Eval(pt)) ||
            phi > fTRDboundariesNeg[3]->Eval(pt),
        phi < fTRDboundariesPos[0]->Eval(pt) ||
            (phi > fTRDboundariesPos[1]->Eval(pt) && phi < fTRDboundariesPos[2]->Eval(pt)) ||
            phi > fTRDboundariesPos[3]->Eval(pt)};
    bool positive = sign > 0;
    for (int iSpecies = 0; iSpecies < 4; ++iSpecies)
    {
      if (std::abs(tpcNsigmas[iSpecies]) < 3.)
      {
        if (track->GetTPCsignal() > fMindEdx)
          fHist2Matching[iSpecies][positive][withTRD[positive]]->Fill(ptot, mass * mass);
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
          double mcMass = trueMatch ? AliPID::ParticleMass(fgkSpecies[iSpecies]) : 0;
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

  //
  // post the data
  //
  PostData(1, fOutputList);
} // end the UserExec
