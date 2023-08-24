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
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <THistManager.h>
#include <TLinearBinning.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEmcalJetConstituentQA.h"
#include "AliAODInputHandler.h"
#include "AliClusterContainer.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliEmcalTriggerDecisionContainer.h"
#include "AliEmcalJet.h"
#include "AliInputEventHandler.h"
#include "AliJetContainer.h"
#include "AliLog.h"
#include "AliTrackContainer.h"
#include "AliVCluster.h"
#include "AliVEvent.h"
#include "AliVTrack.h"

/// \cond CLASSIMP
ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetConstituentQA);
/// \endcond

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalJetConstituentQA::AliAnalysisTaskEmcalJetConstituentQA():
  AliAnalysisTaskEmcalJet(),
  fHistos(nullptr),
  fJetType(AliJetContainer::kFullJet),
  fNameTrackContainer(""),
  fNameClusterContainer(""),
  fTriggerSelectionString(""),
  fUseTriggerSelection(kFALSE),
  fNameTriggerDecisionContainer("EmcalTriggerDecision"),
  fDoHighZClusters(kTRUE),
  fDoCharged(kTRUE),
  fDoNeutral(kTRUE),
  fDoLeadingCharged(kTRUE),
  fDoLeadingNeutral(kTRUE),
  fDoLeadingCell(kTRUE)
{
}

AliAnalysisTaskEmcalJetConstituentQA::AliAnalysisTaskEmcalJetConstituentQA(const char *name):
  AliAnalysisTaskEmcalJet(name, true),
  fHistos(nullptr),
  fJetType(AliJetContainer::kFullJet),
  fNameTrackContainer(""),
  fNameClusterContainer(""),
  fTriggerSelectionString(""),
  fUseTriggerSelection(kFALSE),
  fNameTriggerDecisionContainer("EmcalTriggerDecision"),
  fDoHighZClusters(kTRUE),
  fDoCharged(kTRUE),
  fDoNeutral(kTRUE),
  fDoLeadingCharged(kTRUE),
  fDoLeadingNeutral(kTRUE),
  fDoLeadingCell(kTRUE)
{
  this->SetMakeGeneralHistograms(true);
}

AliAnalysisTaskEmcalJetConstituentQA::~AliAnalysisTaskEmcalJetConstituentQA(){
  if(fHistos) delete fHistos;
}

void AliAnalysisTaskEmcalJetConstituentQA::UserCreateOutputObjects(){
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  TLinearBinning binningz(50, 0., 1), multbinning(51, -0.5, 50.5), binningnef(50, 0., 1.), binningR(60, 0., 0.6), binningptconst(300, 0., 300.), binningptjet(30, 0., 300.),
                 binningNCell(101, -0.5, 100.5), binningFracCellLeading(100, 0., 1.), binningM02(100, 0., 1.), etabinning(100, -0.8, 0.8), phibinning(100, 0., TMath::TwoPi());

  const TBinning *jetbinning[4] = {&binningptjet, &binningnef, &multbinning, &multbinning},
                 *chargedbinning[7] = {&binningptjet, &binningnef, &multbinning, &multbinning, &binningptconst, &binningz, &binningR},
                 *neutralbinning[9] = {&binningptjet, &binningnef, &multbinning, &multbinning, &binningptconst, &binningptconst, &binningz, &binningR, &binningptconst},
                 *binningHighZClusters[7] = {&binningptjet, &binningnef, &binningptconst, &binningz, &binningNCell, &binningFracCellLeading, &binningM02},
                 *leadingchargedbinning[5] = {&binningptjet, &binningnef, &binningptconst, &binningz, &binningR},
                 *leadingneutralbinning[7] = {&binningptjet, &binningnef, &binningptconst, &binningz, &binningR, &binningptconst, &binningFracCellLeading},
                 *leadingcellbinning[5] = {&binningptjet, &binningnef, &binningptconst, &binningptconst, &binningFracCellLeading},
                 *leadingjetvecbinning[4] = {&binningptjet, &etabinning, &phibinning, &binningptjet};

  fHistos = new THistManager(Form("histos_%s", GetName()));
  for(auto c : fNamesJetContainers){
    auto contname = dynamic_cast<TObjString *>(c);
    if(!contname) continue;
    fHistos->CreateTHnSparse(Form("hJetCounter%s", contname->String().Data()), Form("jet counter for jets %s", contname->String().Data()), 4, jetbinning);
    fHistos->CreateTHnSparse(Form("hPtEtaPhiELeadingJet%s", contname->String().Data()), Form("Momemtum vector of leading jets %s", contname->String().Data()), 4, leadingjetvecbinning);
    if(fDoCharged && (fJetType == AliJetContainer::kFullJet || fJetType == AliJetContainer::kChargedJet)){
      fHistos->CreateTHnSparse(Form("hChargedConstituents%s", contname->String().Data()), Form("charged constituents in jets %s", contname->String().Data()), 7, chargedbinning);
      if(fDoLeadingCharged){
        fHistos->CreateTHnSparse(Form("hLeadingTrack%s", contname->String().Data()), Form("leading charged constituent in jets %s", contname->String().Data()), 5, leadingchargedbinning);
        fHistos->CreateTHnSparse(Form("hLeadingJetLeadingTrack%s", contname->String().Data()), Form("leading charged constituent in jets %s", contname->String().Data()), 5, leadingchargedbinning);
      }
    }
    if(fDoNeutral && (fJetType == AliJetContainer::kFullJet || fJetType == AliJetContainer::kNeutralJet)){
      fHistos->CreateTHnSparse(Form("hNeutralConstituents%s", contname->String().Data()), Form("neutral constituents in jets %s", contname->String().Data()), 9, neutralbinning);
      if(fDoHighZClusters) fHistos->CreateTHnSparse(Form("hHighZClusters%s", contname->String().Data()), "Properties of high-z clusters", 7, binningHighZClusters);
      if(fDoLeadingCell) fHistos->CreateTHnSparse(Form("hLeadingCell%s", contname->String().Data()), Form("Leading cell in jets %s", contname->String().Data()), 5, leadingcellbinning);
      if(fDoLeadingNeutral){
        fHistos->CreateTHnSparse(Form("hLeadingCluster%s", contname->String().Data()), Form("leading neutral constituent in jets %s", contname->String().Data()), 7, leadingneutralbinning);
        fHistos->CreateTHnSparse(Form("hLeadingJetLeadingCluster%s", contname->String().Data()), Form("leading neutral constituent in jets %s", contname->String().Data()), 7, leadingneutralbinning);
      }
    }
  }

  for(auto h : *(fHistos->GetListOfHistograms())) fOutput->Add(h);  
  PostData(1, fOutput);
}

bool AliAnalysisTaskEmcalJetConstituentQA::Run(){
  AliParticleContainer * tracks = GetTrackContainer(fNameTrackContainer);
  if(!tracks) tracks = GetParticleContainer(fNameTrackContainer);
  const auto clusters = GetClusterContainer(fNameClusterContainer);
  if(fNameTrackContainer.Length() && !tracks){
    AliErrorStream() << "Track container " << fNameTrackContainer << " required but missing ..." << std::endl;
    return kFALSE;
  }
  if(fNameClusterContainer.Length() && !clusters){
    AliErrorStream() << "Cluster container " << fNameClusterContainer << " required but missing ..." << std::endl;
    return kFALSE;
  }

  for(auto jc : fNamesJetContainers){
    auto contname = dynamic_cast<TObjString *>(jc);
    if(!contname) {
      AliErrorStream() << "Non-string object in the list of jet container names" << std::endl;
      continue;
    } 
    const auto jetcont = GetJetContainer(contname->String().Data());
    if(!jetcont){
      AliErrorStream() << "Jet container with name " << contname->String() << " not found in the list of jet containers" << std::endl;
      continue;
    } 
    AliDebugStream(2) << "Reading " << jetcont->GetArray()->GetName() << std::endl;

    AliEmcalJet *leadingjet(nullptr);
    for(auto jet : jetcont->accepted()){
      if(!leadingjet || jet->Pt() > leadingjet->Pt()) leadingjet = jet;
      AliDebugStream(3) << "Next accepted jet, found " << jet->GetNumberOfTracks() << " tracks and " << jet->GetNumberOfClusters() << " clusters." << std::endl;
      Double_t pointjet[4] = {std::abs(jet->Pt()), jet->NEF(), static_cast<double>(jet->GetNumberOfTracks()), static_cast<double>(jet->GetNumberOfClusters())}, 
               pointcharged[7] = {std::abs(jet->Pt()), jet->NEF(), static_cast<double>(jet->GetNumberOfTracks()), static_cast<double>(jet->GetNumberOfClusters()), -1., 1., -1.}, 
               pointneutral[9] = {std::abs(jet->Pt()), jet->NEF(), static_cast<double>(jet->GetNumberOfTracks()), static_cast<double>(jet->GetNumberOfClusters()), -1., 1., -1., -1.},
               pointHighZCluster[7] = {std::abs(jet->Pt()), jet->NEF(), -1., -1., -1., -1., -1.};
      fHistos->FillTHnSparse(Form("hJetCounter%s", contname->String().Data()), pointjet);
      TVector3 jetvec{jet->Px(), jet->Py(), jet->Pz()};
      if(fDoCharged && (fJetType == AliJetContainer::kFullJet || fJetType == AliJetContainer::kChargedJet) && tracks){
        for(decltype(jet->GetNumberOfTracks()) itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++){
          const auto trk = jet->TrackAt(itrk, tracks->GetArray());
          if(!trk) continue;
          if(trk->Charge()){
            pointcharged[4] = std::abs(trk->Pt());
            pointcharged[5] = std::abs(jet->GetZ(trk));
            pointcharged[6] = jet->DeltaR(trk);
            fHistos->FillTHnSparse(Form("hChargedConstituents%s", contname->String().Data()), pointcharged);
          } else {
            // particle level jets
            pointneutral[4] = pointneutral[5] = std::abs(trk->E());
            pointneutral[6] = std::abs(jet->GetZ(trk));
            pointneutral[7] = jet->DeltaR(trk);
            fHistos->FillTHnSparse(Form("hNeutralConstituents%s", contname->String().Data()), pointneutral);
          }
        }
        if(fDoLeadingCharged){
          // Leading track
          auto leadingtrack = jet->GetLeadingTrack(tracks->GetArray());
          if(leadingtrack) {
            double ltrackpoint[5] = {std::abs(jet->Pt()), jet->NEF(), std::abs(leadingtrack->Pt()), jet->GetZ(leadingtrack), jet->DeltaR(leadingtrack)};
            fHistos->FillTHnSparse(Form("hLeadingTrack%s", contname->String().Data()), ltrackpoint);
          }
        }
      }
      if(fDoNeutral && (fJetType == AliJetContainer::kFullJet || fJetType == AliJetContainer::kNeutralJet) && clusters){
        for(decltype(jet->GetNumberOfClusters()) icl = 0; icl < jet->GetNumberOfClusters(); icl++){
          const auto clust = jet->ClusterAt(icl, clusters->GetArray());
          std::vector<double> fracamp(clust->GetNCells());
          memcpy(fracamp.data(), clust->GetCellsAmplitudeFraction(), sizeof(double) * clust->GetNCells());
          if(!clust) continue; 
          TLorentzVector ptvec;
          double maxEcell = 0.;
          for(auto icell = 0; icell < clust->GetNCells(); icell++) {
            auto ecell =  fInputEvent->GetEMCALCells()->GetCellAmplitude(clust->GetCellAbsId(icell));
            if(ecell > maxEcell) maxEcell = ecell;
          }
          double fracmax = maxEcell / clust->E();
          clust->GetMomentum(ptvec, this->fVertex, AliVCluster::kHadCorr);
          pointneutral[4] = std::abs(clust->GetHadCorrEnergy());
          pointneutral[5] = std::abs(clust->GetNonLinCorrEnergy());
          pointneutral[6] = jet->GetZ(ptvec.Px(), ptvec.Py(), ptvec.Pz());
          pointneutral[7] = jetvec.DeltaR(ptvec.Vect());
          pointneutral[8] = maxEcell;
          fHistos->FillTHnSparse(Form("hNeutralConstituents%s", contname->String().Data()), pointneutral);
          if(fDoHighZClusters && (pointneutral[6] > 0.95)) {
            pointHighZCluster[2] = pointneutral[4];
            pointHighZCluster[3] = pointneutral[6];
            pointHighZCluster[4] = clust->GetNCells();
            pointHighZCluster[5] = fracmax;
            pointHighZCluster[6] = clust->GetM02();
            fHistos->FillTHnSparse(Form("hHighZClusters%s", contname->String().Data()), pointHighZCluster);
          }
        }
        // Leading cluster
        auto leadingcluster = jet->GetLeadingCluster(clusters->GetArray());
        if(leadingcluster){
          TLorentzVector pvect;
          leadingcluster->GetMomentum(pvect, fVertex);
          double maxEcell = 0.;
          for(auto icell = 0; icell < leadingcluster->GetNCells(); icell++) {
            auto ecell =  fInputEvent->GetEMCALCells()->GetCellAmplitude(leadingcluster->GetCellAbsId(icell));
            AliDebugStream(3) << icell << " Pos "  << leadingcluster->GetCellAbsId(icell) <<  " Ecell " << ecell << std::endl;
            if(ecell > maxEcell) maxEcell = ecell;
          }
          double fracmax = maxEcell / leadingcluster->E();
          AliDebugStream(3) << "leading cluster Max E: " << maxEcell << ", frac " << fracmax << std::endl;
          double lclusterpoint[7] = {std::abs(jet->Pt()), jet->NEF(), std::abs(pvect.Pt()), jet->GetZ(pvect.Px(), pvect.Py(), pvect.Pz()), jetvec.DeltaR(pvect.Vect()), std::abs(maxEcell), fracmax};
          fHistos->FillTHnSparse(Form("hLeadingCluster%s", contname->String().Data()), lclusterpoint);
          double leadingCellE = 0.;

          double cellpoint[5] = {std::abs(jet->Pt()), jet->NEF(), std::abs(pvect.Pt()), leadingCellE, fracmax};
          fHistos->FillTHnSparse(Form("hLeadingCell%s", contname->String().Data()), cellpoint);
        }
      }
    }

    // Look at leading particles in the leading jet
    double leadingvec[4] = {0., 0., 0., 0.};
    if(leadingjet) {
      TVector3 leadingjetvec{leadingjet->Px(), leadingjet->Py(), leadingjet->Pz()};
      leadingvec[0] = std::abs(leadingjet->Pt());
      leadingvec[1] = leadingjet->Eta();
      leadingvec[2] = leadingjet->Phi();
      if(leadingvec[2] < 0) leadingvec[2] += TMath::TwoPi();
      leadingvec[3] = leadingjet->E();
      if(fDoCharged && fDoLeadingCharged && (fJetType == AliJetContainer::kFullJet || fJetType == AliJetContainer::kChargedJet) && tracks){
        auto leadingtrack = leadingjet->GetLeadingTrack(tracks->GetArray());
        if(leadingtrack) {
          double ltrackpoint[5] = {std::abs(leadingjet->Pt()), leadingjet->NEF(), std::abs(leadingtrack->Pt()), leadingjet->GetZ(leadingtrack), leadingjet->DeltaR(leadingtrack)};
          fHistos->FillTHnSparse(Form("hLeadingTrack%s", contname->String().Data()), ltrackpoint);
        }
      }
      if(fDoNeutral && fDoLeadingNeutral && (fJetType == AliJetContainer::kFullJet || fJetType == AliJetContainer::kNeutralJet) && clusters){
        auto leadingcluster = leadingjet->GetLeadingCluster(clusters->GetArray());
        if(leadingcluster){
          TLorentzVector pvect;
          leadingcluster->GetMomentum(pvect, fVertex);
          double maxEcell = 0.;
          for(auto icell = 0; icell < leadingcluster->GetNCells(); icell++) {
            auto ecell =  fInputEvent->GetEMCALCells()->GetCellAmplitude(leadingcluster->GetCellAbsId(icell));
            if(ecell > maxEcell) maxEcell = ecell;
          }
          double fracmax = maxEcell / leadingcluster->E();
          double lclusterpoint[7] = {std::abs(leadingjet->Pt()), leadingjet->NEF(), std::abs(pvect.Pt()), leadingjet->GetZ(pvect.Px(), pvect.Py(), pvect.Pz()), leadingjetvec.DeltaR(pvect.Vect()), maxEcell, fracmax};
          fHistos->FillTHnSparse(Form("hLeadingJetLeadingCluster%s", contname->String().Data()), lclusterpoint);
        }
      }
    }
    fHistos->FillTHnSparse(Form("hPtEtaPhiELeadingJet%s", contname->String().Data()), leadingvec);
  }
  
  return kTRUE;
}

Bool_t AliAnalysisTaskEmcalJetConstituentQA::IsTriggerSelected(){
  // Event selection
  AliDebugStream(1) << "Trigger selection string: " << fTriggerSelectionString << ", fired trigger classes: " << fInputEvent->GetFiredTriggerClasses() << std::endl;
  if(fTriggerSelectionString.Contains("INT7")){
    // INT7 trigger
    if(!(fInputHandler->IsEventSelected() & AliVEvent::kINT7)) return false;
  } else if(fTriggerSelectionString.Contains("EJ")){
    auto triggerclass = fTriggerSelectionString(fTriggerSelectionString.Index("EJ"),3); // cleanup trigger string from possible tags
    AliDebugStream(1) << "Inspecting trigger class " << triggerclass << std::endl;
    // EMCAL JET trigger
    if(!fMCEvent){ // Running on Data
      if(!(fInputHandler->IsEventSelected() & AliVEvent::kEMCEJE)) return false;
      if(!fInputEvent->GetFiredTriggerClasses().Contains(triggerclass)) return false;
    }

    if(fUseTriggerSelection) {
      auto triggerdecisions = dynamic_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer *>(fInputEvent->FindListObject(fNameTriggerDecisionContainer.Data()));
      if(!triggerdecisions) {
        AliErrorStream() << "No offline trigger selection available" << std::endl;
        return false;
      }
      else if(!triggerdecisions->IsEventSelected(triggerclass.Data())) return false;
    }
  } else return false;

  AliDebugStream(1) << "Event is selected" << std::endl;
  return true;
}

AliAnalysisTaskEmcalJetConstituentQA *AliAnalysisTaskEmcalJetConstituentQA::AddTaskEmcalJetConstituentQA(const char *trigger, AliJetContainer::EJetType_t jettype, bool partmode){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr) {
    std::cerr << "[AliAnalysisTaskJetConstituentQA::AddTaskEmcalJetConstituentQA(EE)] No analysis manager provided ..." << std::endl;
    return nullptr;
  }

  std::string jettypestring;
  switch(jettype){
    case AliJetContainer::kFullJet:    jettypestring = "fulljets"; break;
    case AliJetContainer::kChargedJet: jettypestring = "chargedjets"; break;
    case AliJetContainer::kNeutralJet: jettypestring = "neutraljets"; break;
    case AliJetContainer::kUndefinedJetType: break;
  };
  
  std::stringstream taskname;
  taskname << "constituentQA_" << jettypestring << "_" << trigger;
  auto task = new AliAnalysisTaskEmcalJetConstituentQA(taskname.str().data());
  task->SetTriggerSelection(trigger);
  task->SetJetType(jettype);
  mgr->AddTask(task);

  auto inputhandler = mgr->GetInputEventHandler();
  auto isAOD = false;
  if(inputhandler->IsA() == AliAODInputHandler::Class()) isAOD = true;

  TString tracksname, clustername;
  AliParticleContainer *tracks(nullptr);
  AliClusterContainer *clusters(nullptr);
  if(!partmode) {
    if(jettype ==  AliJetContainer::kChargedJet || jettype == AliJetContainer::kFullJet){
      tracksname = AliEmcalAnalysisFactory::TrackContainerNameFactory(isAOD);
      tracks = task->AddTrackContainer(tracksname);
      task->SetNameTrackContainer(tracksname);
      tracks->SetMinPt(0.15);
    }

    if(jettype == AliJetContainer::kNeutralJet || jettype == AliJetContainer::kFullJet){
      clustername = AliEmcalAnalysisFactory::ClusterContainerNameFactory(isAOD);
      clusters = task->AddClusterContainer(clustername);
      task->SetNameClusterContainer(clustername);
      clusters->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
      clusters->SetClusHadCorrEnergyCut(0.3);
    }
  } else {
    tracksname = "mcparticles";
    tracks = task->AddParticleContainer(tracksname);
    task->SetNameTrackContainer(tracksname);
    tracks->SetMinPt(0.);
  }

  // create jet containers for R02 and R04 jets
  std::array<double, 5> jetradii = {{0.2, 0.3, 0.4, 0.5, 0.6}};
  for(auto r : jetradii) {
    std::stringstream contname;
    contname << jettypestring << "_R" << std::setw(2) << std::setfill('0') << int(r*10.);
    auto jcont = task->AddJetContainer(jettype, AliJetContainer::antikt_algorithm, AliJetContainer::E_scheme, r, AliJetContainer::kEMCALfid, tracks, clusters, "Jet");
    jcont->SetName(contname.str().data());
    task->AddNameJetContainer(contname.str().data());
    jcont->SetMinPt(20.);
  }

  std::stringstream contname, outfilename;
  contname << "JetConstituentQA_" << jettypestring << "_" << trigger;
  outfilename << mgr->GetCommonFileName() << ":JetConstituentQA_" << trigger;
  if(partmode) {
    contname << "_part";
    outfilename << "_part";
  }
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, mgr->CreateContainer(contname.str().data(), AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, outfilename.str().data()));

  return task;
}
