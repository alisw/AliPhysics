/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
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
#include <algorithm>
#include <array>
#include <sstream>
#include <string>
#include <vector>

#include <THistManager.h>
#include <TLinearBinning.h>
#include <TVariableBinning.h>

#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEmcalJetEnergySpectrum.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliEmcalDownscaleFactorsOCDB.h"
#include "AliEmcalJet.h"
#include "AliEmcalTriggerDecisionContainer.h"
#include "AliEmcalTriggerStringDecoder.h"
#include "AliEventCuts.h"
#include "AliInputEventHandler.h"
#include "AliJetContainer.h"
#include "AliLog.h"
#include "AliMultSelection.h"
#include "AliVEvent.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergySpectrum);

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalJetEnergySpectrum::AliAnalysisTaskEmcalJetEnergySpectrum():
  AliAnalysisTaskEmcalJet(),
  fHistos(nullptr),
  fIsMC(false),
  fFillHSparse(false),
	fTriggerSelectionBits(AliVEvent::kAny),
  fTriggerSelectionString(""),
  fRequireSubsetMB(false),
  fMinBiasTrigger(AliVEvent::kAny),
  fNameTriggerDecisionContainer("EmcalTriggerDecision"),
  fUseTriggerSelectionForData(false),
  fUseDownscaleWeight(false),
  fNameJetContainer("datajets"),
  fRequestTriggerClusters(true),
  fRequestCentrality(false),
  fUseAliEventCuts(false),
  fUseSumw2(false),
  fUseMuonCalo(false),
  fScaleShift(0.),
  fCentralityEstimator("V0M"),
  fUserPtBinning()
{
  SetUseAliAnaUtils(true);
}

AliAnalysisTaskEmcalJetEnergySpectrum::AliAnalysisTaskEmcalJetEnergySpectrum(EMCAL_STRINGVIEW name):
  AliAnalysisTaskEmcalJet(name.data(), true),
  fHistos(nullptr),
  fIsMC(false),
  fFillHSparse(false),
	fTriggerSelectionBits(AliVEvent::kAny),
  fTriggerSelectionString(""),
  fRequireSubsetMB(false),
  fMinBiasTrigger(AliVEvent::kAny),
  fNameTriggerDecisionContainer("EmcalTriggerDecision"),
  fUseTriggerSelectionForData(false),
  fUseDownscaleWeight(false),
  fNameJetContainer("datajets"),
  fRequestTriggerClusters(true),
  fRequestCentrality(false),
  fUseAliEventCuts(false),
  fUseSumw2(false),
  fUseMuonCalo(false),
  fScaleShift(0.),
  fCentralityEstimator("V0M"),
  fUserPtBinning()
{
  SetUseAliAnaUtils(true);
}

AliAnalysisTaskEmcalJetEnergySpectrum::~AliAnalysisTaskEmcalJetEnergySpectrum(){
  if(fHistos) delete fHistos;
}

void AliAnalysisTaskEmcalJetEnergySpectrum::UserCreateOutputObjects(){
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  if(!fUserPtBinning.GetSize()) {
    // binning not set. apply default binning
    AliInfoStream() << "Using default pt binning";
    fUserPtBinning.Set(301);
    double current(0.);
    for(int istep = 0; istep < 301; istep++) {
      fUserPtBinning[istep] = current;
      current += 1; 
    }
  }

  fHistos = new THistManager(Form("Histos_%s", GetName()));
  fHistos->CreateTH1("hEventCounter", "Event counter histogram", 1, 0.5, 1.5);
  fHistos->CreateTH1("hEventCounterAbs", "Event counter histogram absolute", 1, 0.5, 1.5);
  fHistos->CreateTH1("hEventCentrality", "Event centrality", 100., 0., 100.);
  fHistos->CreateTH1("hEventCentralityAbs", "Event centrality absolute", 100., 0., 100.);
  fHistos->CreateTH1("hClusterCounter", "Event counter histogram", kTrgClusterN, -0.5, kTrgClusterN - 0.5);
  fHistos->CreateTH1("hClusterCounterAbs", "Event counter histogram absolute", kTrgClusterN, -0.5, kTrgClusterN - 0.5);
  fHistos->CreateTH2("hJetSpectrum", "Jet pt spectrum", kTrgClusterN, -0.5, kTrgClusterN - 0.5, 350., 0., 350., "s");
  fHistos->CreateTH2("hJetSpectrumMax", "Max jet pt spectrum", kTrgClusterN, -0.5, kTrgClusterN - 0.5, 350., 0., 350., "s");
  if(fFillHSparse) {
    TLinearBinning centralitybinning(100, 0., 100.), etabinning(100, -1., 1.), phibinning(100., 0., 7.), nefbinning(100, 0., 1.), trgclusterbinning(kTrgClusterN + 1, -0.5, kTrgClusterN -0.5);
    TVariableBinning jetptbinning(fUserPtBinning);
    const TBinning *binnings[6] = {&centralitybinning, &jetptbinning, &etabinning, &phibinning, &nefbinning, &trgclusterbinning};
    fHistos->CreateTHnSparse("hJetTHnSparse", "jet thnsparse", 6, binnings, fUseSumw2 ? "s" : "");
    fHistos->CreateTHnSparse("hMaxJetTHnSparse", "jet thnsparse", 6, binnings, fUseSumw2 ? "s" : "");
  }

  for(auto h : *fHistos->GetListOfHistograms()) fOutput->Add(h);
  PostData(1, fOutput);
}

Bool_t AliAnalysisTaskEmcalJetEnergySpectrum::CheckMCOutliers() {
  if(!fMCRejectFilter) return true;
  if(!(fIsPythia || fIsHerwig)) return true;    // Only relevant for pt-hard production
  AliDebugStream(1) << "Using custom MC outlier rejection" << std::endl;
  auto partjets = GetJetContainer("partjets");
  if(!partjets) return true;

  // Check whether there is at least one particle level jet with pt above n * event pt-hard
  auto jetiter = partjets->accepted();
  auto max = std::max_element(jetiter.begin(), jetiter.end(), [](const AliEmcalJet *lhs, const AliEmcalJet *rhs ) { return lhs->Pt() < rhs->Pt(); });
  if(max != jetiter.end())  {
    // At least one jet found with pt > n * pt-hard
    AliDebugStream(1) << "Found max jet with pt " << (*max)->Pt() << " GeV/c" << std::endl;
    if((*max)->Pt() > fPtHardAndJetPtFactor * fPtHard) return false;
  }
  return true;
}

bool AliAnalysisTaskEmcalJetEnergySpectrum::Run(){
  auto datajets = this->GetJetContainer(fNameJetContainer);
  if(!datajets) {
    AliErrorStream() << "Jet container " << fNameJetContainer << " not found" << std::endl;
    return false;
  }

  double eventCentrality = 99;   // without centrality put everything in the peripheral bin
  if(fRequestCentrality){
    AliMultSelection *mult = dynamic_cast<AliMultSelection *>(InputEvent()->FindListObject("MultSelection"));
    if(!mult){
      AliErrorStream() << GetName() << ": Centrality selection enabled but no centrality estimator found" << std::endl;
      return false;
    }
    if(mult->IsEventSelected()) return false;
    eventCentrality = mult->GetEstimator(fCentralityEstimator)->GetPercentile();
    AliDebugStream(1) << GetName() << ": Centrality " <<  eventCentrality << std::endl;
  } else {
    AliDebugStream(1) << GetName() << ": No centrality selection applied" << std::endl;
  }


  auto trgclusters = GetTriggerClustersANY();
  if(!fIsMC && fRequestTriggerClusters) trgclusters = GetTriggerClusterIndices(fInputEvent->GetFiredTriggerClasses().Data());
  Double_t weight = 1.;
  if(fUseDownscaleWeight) {
    weight = 1./PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance()->GetDownscaleFactorForTriggerClass(MatchTrigger(fInputEvent->GetFiredTriggerClasses().Data(), fTriggerSelectionString.Data(), fUseMuonCalo));
  }
  fHistos->FillTH1("hEventCounterAbs", 1.);
  fHistos->FillTH1("hEventCounter", weight);
  fHistos->FillTH1("hEventCentralityAbs", eventCentrality);
  fHistos->FillTH1("hEventCentrality", eventCentrality, weight);
  AliEmcalJet *maxjet(nullptr);
  for(auto t : trgclusters) {
    fHistos->FillTH1("hClusterCounterAbs", t);
    fHistos->FillTH1("hClusterCounter", t, weight);
  }
  for(auto j : datajets->accepted()){
    if(!maxjet || (j->E() > maxjet->E())) maxjet = j;
    Double_t ptjet = j->Pt();
    if(TMath::Abs(fScaleShift) > DBL_EPSILON){
      // Apply artificial (fixed) shift of the jet energy scale to det. level jets
      ptjet += fScaleShift * ptjet; 
    }
    double datapoint[6] = {eventCentrality, ptjet, j->Eta(), j->Phi(), j->NEF(), 0.};
    for(auto t : trgclusters){
      fHistos->FillTH2("hJetSpectrum", static_cast<double>(t), ptjet, weight);
      if(fFillHSparse) {
        datapoint[5] = static_cast<double>(t);
        fHistos->FillTHnSparse("hJetTHnSparse", datapoint, weight);
      }
    }
  }

  double maxdata[6];
  memset(maxdata, 0., sizeof(double) * 5);
  maxdata[0] = eventCentrality;
  if(maxjet){
    maxdata[1] = maxjet->Pt();
    if(TMath::Abs(fScaleShift) > DBL_EPSILON) {
      maxdata[1] += fScaleShift * maxdata[1];
    }
    maxdata[2] = maxjet->Eta();
    maxdata[3] = maxjet->Phi();
    maxdata[4] = maxjet->NEF();
  }
  for(auto t : trgclusters){
    fHistos->FillTH2("hJetSpectrumMax", t, maxdata[1], weight);
    if(fFillHSparse){
      maxdata[5] = static_cast<double>(t);
      fHistos->FillTHnSparse("hMaxJetTHnSparse", maxdata, weight);
    }
  }
  return true;
}

void AliAnalysisTaskEmcalJetEnergySpectrum::RunChanged(Int_t newrun){
  if(fUseDownscaleWeight) {
    auto downscalehandler = PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance();
    if(downscalehandler->GetCurrentRun() != newrun){
      downscalehandler->SetRun(newrun);
    }
  }
}

bool AliAnalysisTaskEmcalJetEnergySpectrum::IsTriggerSelected() {
  if(!fIsMC){
    // Pure data - do EMCAL trigger selection from selection string
    UInt_t triggerbits = fTriggerSelectionBits;
    if(fUseMuonCalo) fTriggerSelectionBits = AliVEvent::kMuonCalo;  // in case of the muon-calo / calo(fast) cluster all data is in the 
    if(!(fInputHandler->IsEventSelected() & triggerbits)) return false;
    if(fTriggerSelectionString.Length()) {
      if(!fInputEvent->GetFiredTriggerClasses().Contains(fTriggerSelectionString)) return false;
      if(fRequireSubsetMB && !(fInputHandler->IsEventSelected() & fMinBiasTrigger)) return false;   // Require EMCAL trigger to be subset of the min. bias trigger (for efficiency studies)
      if((fTriggerSelectionString.Contains("EJ") || fTriggerSelectionString.Contains("EG") || fTriggerSelectionString.Contains("DJ") || fTriggerSelectionString.Contains("DG")) && fUseTriggerSelectionForData) {
        auto trgselresult = static_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer *>(fInputEvent->FindListObject(fNameTriggerDecisionContainer));
        AliDebugStream(1) << "Found trigger decision object: " << (trgselresult ? "yes" : "no") << std::endl;
        if(!trgselresult){
          AliErrorStream() <<  "Trigger decision container with name " << fNameTriggerDecisionContainer << " not found in event - not possible to select EMCAL triggers" << std::endl;
          return false;
        }
        if(!trgselresult->IsEventSelected(fTriggerSelectionString)) return false;
      }
    }
  } else {
    if(!(fInputHandler->IsEventSelected() & AliVEvent::kINT7)) return false;
    if(IsSelectEmcalTriggers(fTriggerSelectionString.Data())){
      // Simulation - do EMCAL trigger selection from trigger selection object
      auto mctrigger = static_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer *>(fInputEvent->FindListObject(fNameTriggerDecisionContainer));
      AliDebugStream(1) << "Found trigger decision object: " << (mctrigger ? "yes" : "no") << std::endl;
      if(!mctrigger){
        AliErrorStream() <<  "Trigger decision container with name " << fNameTriggerDecisionContainer << " not found in event - not possible to select EMCAL triggers" << std::endl;
        return false;
      }
      if(!mctrigger->IsEventSelected(fTriggerSelectionString)) return false;
    }
  }
  return true;
}

AliAnalysisTaskEmcalJetEnergySpectrum *AliAnalysisTaskEmcalJetEnergySpectrum::AddTaskJetEnergySpectrum(Bool_t isMC, AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recoscheme, AliVCluster::VCluUserDefEnergy_t energydef, double radius, EMCAL_STRINGVIEW namepartcont, EMCAL_STRINGVIEW trigger, EMCAL_STRINGVIEW suffix){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) {
    std::cerr << "Analysis manager not initialized" << std::endl;
    return nullptr;
  }

  Bool_t isAOD(kFALSE);
  AliInputEventHandler *inputhandler = static_cast<AliInputEventHandler *>(mgr->GetInputEventHandler());
  if(inputhandler) {
    if(inputhandler->IsA() == AliAODInputHandler::Class()){
      std::cout << "Analysing AOD events\n";
      isAOD = kTRUE;
    } else {
      std::cout << "Analysing ESD events\n";
    }
  }

  std::string jettypestring;
  UInt_t acctype(AliJetContainer::kTPCfid);
  switch(jettype){
    case AliJetContainer::kChargedJet:  jettypestring = "ChargedJets"; acctype = AliJetContainer::kTPCfid; break;
    case AliJetContainer::kFullJet:     jettypestring = "FullJets";    acctype = AliJetContainer::kEMCALfid; break;
    case AliJetContainer::kNeutralJet:  jettypestring = "NeutralJets"; acctype = AliJetContainer::kEMCALfid; break;
    case AliJetContainer::kUndefinedJetType: break;
  };

  std::stringstream tag, outfilename;
  tag << jettypestring << "_R" << std::setw(2) << std::setfill('0') << int(radius * 10.) << "_" << trigger;
  if(suffix.length()) {
    tag << "_" << suffix;
  }
  auto task = new AliAnalysisTaskEmcalJetEnergySpectrum(Form("JetEnergySpectrum_%s", tag.str().data()));
  task->SetIsMC(isMC);
  mgr->AddTask(task);

  auto contains = [](EMCAL_STRINGVIEW str, EMCAL_STRINGVIEW test) {
    return str.find(test) != std::string::npos;
  };

  std::string trgstr(trigger);
  if(contains(trgstr, "INT7")) task->SetTriggerSelection(AliVEvent::kINT7, "INT7");
  else if(contains(trgstr, "EJ1")) task->SetTriggerSelection(AliVEvent::kEMCEJE, "EJ1");
  else if(contains(trgstr, "EJ2")) task->SetTriggerSelection(AliVEvent::kEMCEJE, "EJ2");
  else if(contains(trgstr, "EG1")) task->SetTriggerSelection(AliVEvent::kEMCEGA, "EG1");
  else if(contains(trgstr, "EG2")) task->SetTriggerSelection(AliVEvent::kEMCEGA, "EG2");

  // Connect particle and cluster container
  AliTrackContainer *tracks(nullptr);
  AliClusterContainer *clusters(nullptr);
  if(jettype == AliJetContainer::kChargedJet || jettype == AliJetContainer::kFullJet) {
    tracks = task->AddTrackContainer(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::TrackContainerNameFactory(isAOD));
    tracks->SetMinPt(0.15);
  }
  if(jettype == AliJetContainer::kNeutralJet || jettype == AliJetContainer::kFullJet){
    clusters = task->AddClusterContainer(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::ClusterContainerNameFactory(isAOD));
    clusters->SetDefaultClusterEnergy(energydef);
    clusters->SetClusUserDefEnergyCut(energydef, 0.3);
  }


  // Create proper jet container
  auto jetcont = task->AddJetContainer(jettype, AliJetContainer::antikt_algorithm, recoscheme, radius, acctype, tracks, clusters);
  jetcont->SetName("datajets");
  task->SetNameJetContainer("datajets");
  std::cout << "Adding jet container with underlying array:" << jetcont->GetArrayName() << std::endl;

  if(isMC){
    // Create also particle and particle level jet container for outlier rejection
    TString partcontname = namepartcont;
    if(partcontname == "usedefault") partcontname = "mcparticles";
    auto partcont = task->AddMCParticleContainer(partcontname.Data());
    partcont->SetMinPt(0.);
    
    //AliJetContainer::EJetType_t mcjettype = (jettype == AliJetContainer::kNeutralJet) ? AliJetContainer::kFullJet : jettype;
    AliJetContainer::EJetType_t mcjettype = AliJetContainer::kFullJet;
    auto pjcont = task->AddJetContainer(mcjettype, AliJetContainer::antikt_algorithm, recoscheme, radius, AliJetContainer::kTPCfid, partcont, nullptr);
    pjcont->SetName("partjets");
    pjcont->SetMinPt(0);
    pjcont->SetMaxTrackPt(1000.);
  }

  // Link input and output container
  outfilename << mgr->GetCommonFileName() << ":JetSpectrum_" << tag.str().data();
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(Form("JetSpectrum_%s", tag.str().data()), TList::Class(), AliAnalysisManager::kOutputContainer, outfilename.str().data()));

  return task;
}
