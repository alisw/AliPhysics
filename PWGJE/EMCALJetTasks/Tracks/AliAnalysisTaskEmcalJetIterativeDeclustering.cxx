/************************************************************************************
 * Copyright (C) 2019, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#include <sstream>
#include <THistManager.h>

#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskEmcalJetIterativeDeclustering.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliEmcalDownscaleFactorsOCDB.h"
#include "AliEmcalTriggerStringDecoder.h"
#include "AliInputEventHandler.h"
#include "AliLundPlaneHelper.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetIterativeDeclustering)

    using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalJetIterativeDeclustering::AliAnalysisTaskEmcalJetIterativeDeclustering() : AliAnalysisTaskEmcalJet(),
                                                                                               fDecluster(nullptr),
                                                                                               fHistos(nullptr),
                                                                                               fUseNeutralConstituents(true),
                                                                                               fTriggerBits(0),
                                                                                               fTriggerString(),
                                                                                               fJetPtMin(0.),
                                                                                               fJetPtMax(1000.),
                                                                                               fHardCutoff(-1.)
{
}

AliAnalysisTaskEmcalJetIterativeDeclustering::AliAnalysisTaskEmcalJetIterativeDeclustering(EMCAL_STRINGVIEW name) : AliAnalysisTaskEmcalJet(name.data(), true),
                                                                                                                    fDecluster(nullptr),
                                                                                                                    fHistos(nullptr),
                                                                                                                    fUseDownscaleWeight(false),
                                                                                                                    fUseChargedConstituents(true),
                                                                                                                    fUseNeutralConstituents(true),
                                                                                                                    fTriggerBits(0),
                                                                                                                    fTriggerString(),
                                                                                                                    fJetPtMin(0.),
                                                                                                                    fJetPtMax(1000.),
                                                                                                                    fHardCutoff(-1.)
{
}

AliAnalysisTaskEmcalJetIterativeDeclustering::~AliAnalysisTaskEmcalJetIterativeDeclustering()
{
}

void AliAnalysisTaskEmcalJetIterativeDeclustering::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  fDecluster = new PWGJE::EMCALJetTasks::AliLundPlaneHelper();
  if (fHardCutoff >= 0.)
    fDecluster->SetHardCutoff(fHardCutoff);

  TAxis fJetPtBinning(300, 0., 300.),
      fLnPtRelBinning(25, -10., 2.),
      fLnDeltaRBinning(25, 0., 5.),
      fNSplittingsBinning(100, 0., 100);

  const TAxis *binnings[5] = {&fJetPtBinning, &fJetPtBinning, &fLnDeltaRBinning, &fLnPtRelBinning, &fNSplittingsBinning};

  fHistos = new THistManager(Form("Histos%s", GetName()));
  fHistos->CreateTH1("hEventCounter", "Number of events", 1, 0.5, 1.5);
  fHistos->CreateTH1("hEventCounterWeighted", "Number of events", 1, 0.5, 1.5);
  fHistos->CreateTH2("hNSplittings", "Number of splittings", 300, 0., 300., 100, 0., 100.);
  fHistos->CreateTHnSparse("hSplittings", "Splittings", 5, binnings);

  for (auto h : *fHistos->GetListOfHistograms()) fOutput->Add(h);
  PostData(1, fOutput);
}

Bool_t AliAnalysisTaskEmcalJetIterativeDeclustering::Run()
{
  auto jets = GetJetContainer("datajets");
  if (!jets)
  {
    AliErrorStream() << "Jet container not found" << std::endl;
    return false;
  }
  auto clusters = GetClusterContainer(0);
  if (fUseNeutralConstituents && !clusters)
  {
    AliErrorStream() << "Cluster container not found, but neutral constituents requested" << std::endl;
  }
  auto tracks = GetTrackContainer(0);
  if (fUseChargedConstituents && !tracks)
  {
    AliErrorStream() << "Track container not found, but charged constituent requested." << std::endl;
    return false;
  }

  Double_t weight = fUseDownscaleWeight ? GetDownscaleWeight() : 1.;
  fHistos->FillTH1("hEventCounter", 1.);
  if (fUseDownscaleWeight)
    fHistos->FillTH1("hEventCounterWeighted", 1., weight);

  for (auto jet : jets->accepted())
  {
    if (jet->Pt() < fJetPtMin || jet->Pt() > fJetPtMax)
      continue;
    AliDebugStream(1) << "Jet found, pt = " << jet->Pt() << " GeV/c" <<  std::endl;
    auto lundplane = fDecluster->Evaluate(*jet, tracks, clusters, fVertex);
    auto splittings = lundplane.GetSplittings();
    fHistos->FillTH2("hNSplittings", jet->Pt(), static_cast<double>(splittings.size()));
    for (auto splitting : splittings)
    {
      AliDebugStream(2) << "Next splitting: " << splitting.GetNSplittings() << ", angle " << splitting.GetLnDeltaR() << ", rel kt " << splitting.GetLnPtrel() << ", pt lower " << splitting.GetPtLower() << std::endl;
      Double_t point[] = {jet->Pt(), splitting.GetPtLower(), splitting.GetLnDeltaR(), splitting.GetLnPtrel(), static_cast<double>(splitting.GetNSplittings())};
      fHistos->FillTHnSparse("hSplittings", point, weight);
    }
    AliDebugStream(2) << "Jet done" << std::endl;
  }

  return true;
}

bool AliAnalysisTaskEmcalJetIterativeDeclustering::IsTriggerSelected()
{
  if (!(fInputHandler->IsEventSelected() & fTriggerBits))
    return false;
  if (fTriggerString.length())
  {
    if (!fInputEvent->GetFiredTriggerClasses().Contains(fTriggerString))
      return false;
  }
  return true;
}
void AliAnalysisTaskEmcalJetIterativeDeclustering::RunChanged(Int_t newrun)
{
  if (fUseDownscaleWeight)
    PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance()->SetRun(newrun);
}

Double_t AliAnalysisTaskEmcalJetIterativeDeclustering::GetDownscaleWeight() const
{
  Double_t weight = 1.;
  TString triggerclass;
  if (fTriggerString == "INT7")
    triggerclass = "CINT7-B-NOPF-CENT";
  else if (fTriggerString == "EJ1")
    triggerclass = "CEMC7EJ1-B-NOPF-CENTNOTRD";
  else if (fTriggerString == "EJ2")
    triggerclass = "CEMC7EJ1-B-NOPF-CENT";
  if (triggerclass.Length())
    weight = PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance()->GetDownscaleFactorForTriggerClass(triggerclass);
  return weight;
}

AliAnalysisTaskEmcalJetIterativeDeclustering *AliAnalysisTaskEmcalJetIterativeDeclustering::AddTaskEmcalJetIterativeDeclustering(Double_t jetradius, AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recombinationScheme, EMCAL_STRINGVIEW trigger)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  Bool_t isAOD(kFALSE);
  AliInputEventHandler *inputhandler = static_cast<AliInputEventHandler *>(mgr->GetInputEventHandler());
  if (inputhandler)
  {
    if (inputhandler->IsA() == AliAODInputHandler::Class())
    {
      std::cout << "Analysing AOD events\n";
      isAOD = kTRUE;
    }
    else
    {
      std::cout << "Analysing ESD events\n";
    }
  }

  std::stringstream taskname;
  taskname << "SoftdropDataMaker_R" << std::setw(2) << std::setfill('0') << int(jetradius * 10) << trigger;
  AliAnalysisTaskEmcalJetIterativeDeclustering *datamaker = new AliAnalysisTaskEmcalJetIterativeDeclustering(taskname.str());
  datamaker->SelectCollisionCandidates(AliVEvent::kINT7);
  mgr->AddTask(datamaker);

  ULong_t triggerbits = 0;
  if(trigger == "INT7") triggerbits = AliVEvent::kINT7;
  else if(trigger == "EJ1" || trigger == "EJ2") triggerbits = AliVEvent::kEMCEJE;
  datamaker->SetSelectTrigger(triggerbits, trigger.data());

  AliTrackContainer *tracks(nullptr);
  if ((jettype == AliJetContainer::kChargedJet) || (jettype == AliJetContainer::kFullJet))
  {
    tracks = datamaker->AddTrackContainer(AliEmcalAnalysisFactory::TrackContainerNameFactory(isAOD));
    std::cout << "Track container name: " << tracks->GetName() << std::endl;
    tracks->SetMinPt(0.15);
  }
  AliClusterContainer *clusters(nullptr);
  if ((jettype == AliJetContainer::kFullJet) || (jettype == AliJetContainer::kNeutralJet))
  {
    std::cout << "Using full or neutral jets ..." << std::endl;
    clusters = datamaker->AddClusterContainer(AliEmcalAnalysisFactory::ClusterContainerNameFactory(isAOD));
    std::cout << "Cluster container name: " << clusters->GetName() << std::endl;
    clusters->SetClusHadCorrEnergyCut(0.3); // 300 MeV E-cut
    clusters->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }
  else
  {
    std::cout << "Using charged jets ... " << std::endl;
  }

  AliJetContainer *datajets = datamaker->AddJetContainer(
      jettype,
      AliJetContainer::antikt_algorithm,
      recombinationScheme,
      jetradius,
      ((jettype == AliJetContainer::kFullJet) || (jettype == AliJetContainer::kNeutralJet)) ? AliEmcalJet::kEMCALfid : AliEmcalJet::kTPCfid,
      tracks, clusters);
  datajets->SetName("datajets");
  datajets->SetJetPtCut(0.);
  datajets->SetMaxTrackPt(1000.);

  std::string jettypestring;
  switch (jettype)
  {
  case AliJetContainer::kFullJet:
    jettypestring = "FullJets";
    break;
  case AliJetContainer::kChargedJet:
    jettypestring = "ChargedJets";
    break;
  case AliJetContainer::kNeutralJet:
    jettypestring = "NeutralJets";
    break;
  default:
    jettypestring = "Undef";
  };

  // Connecting containers
  std::stringstream outputfile, histname;
  outputfile << mgr->GetCommonFileName() << ":IterativeDeclustering" << jettypestring << "_R" << std::setw(2) << std::setfill('0') << int(jetradius * 10.) << "_" << trigger;
  histname << "IterativeDeclusteringHistos_" << jettypestring << "_R" << std::setw(2) << std::setfill('0') << int(jetradius * 10.) << "_" << trigger;
  mgr->ConnectInput(datamaker, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(datamaker, 1, mgr->CreateContainer(histname.str().data(), AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, outputfile.str().data()));

  return datamaker;
}
