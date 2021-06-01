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
#include <bitset>
#include <iostream>
#include <string>
#include <set>
#include <sstream>
#include <vector>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/contrib/Nsubjettiness.hh>
#include <fastjet/contrib/SoftDrop.hh>
#include <fastjet/config.h>
#if FASJET_VERSION_NUMBER >= 30302
#include <fastjet/tools/Recluster.hh>
#else 
#include <fastjet/contrib/Recluster.hh>
#endif

#include <THistManager.h>
#include <TLinearBinning.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjString.h>
#include <TString.h>
#include <TVector3.h>

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisTaskEmcalJetSubstructureTree.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliClusterContainer.h"
#include "AliJetContainer.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliEmcalDownscaleFactorsOCDB.h"
#include "AliEmcalJet.h"
#include "AliEmcalList.h"
#include "AliEmcalTriggerDecisionContainer.h"
#include "AliEmcalTriggerStringDecoder.h"
#include "AliLog.h"
#include "AliParticleContainer.h"
#include "AliRhoParameter.h"
#include "AliTrackContainer.h"
#include "AliTriggerCluster.h"
#include "AliTriggerConfiguration.h"
#include "AliVCluster.h"
#include "AliVParticle.h"

#ifdef EXPERIMENTAL_JETCONSTITUENTS
#include "AliEmcalClusterJetConstituent.h"
#include "AliEmcalParticleJetConstituent.h"
#endif

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetSubstructureTree);

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalJetSubstructureTree::AliAnalysisTaskEmcalJetSubstructureTree() :
    AliAnalysisTaskEmcalJet(),
    fJetSubstructureTree(nullptr),
    fGlobalTreeParams(nullptr),
    fSoftDropMeasured(nullptr),
    fSoftDropTrue(nullptr),
    fNSubMeasured(nullptr),
    fNSubTrue(nullptr),
    fKineRec(nullptr),
    fKineSim(nullptr),
    fJetStructureMeasured(nullptr),
    fJetStructureTrue(nullptr),
    fQAHistos(nullptr),
    fLumiMonitor(nullptr),
    fSDZCut(0.1),
    fSDBetaCut(0),
    fReclusterizer(kCAAlgo),
    fHasRecEvent(false),
    fHasTrueEvent(false),
    fTriggerSelectionBits(AliVEvent::kAny),
    fTriggerSelectionString(""),
    fNameTriggerDecisionContainer("EmcalTriggerDecision"),
    fUseTriggerSelectionForData(false),
    fUseDownscaleWeight(false),
    fUseChargedConstituents(true),
    fUseNeutralConstituents(true),
    fFillPart(true),
    fFillRho(true),
    fFillSoftDrop(true),
    fFillNSub(true),
    fFillStructGlob(true)
{
}

AliAnalysisTaskEmcalJetSubstructureTree::AliAnalysisTaskEmcalJetSubstructureTree(const char *name) :
    AliAnalysisTaskEmcalJet(name, kTRUE),
    fJetSubstructureTree(nullptr),
    fGlobalTreeParams(nullptr),
    fSoftDropMeasured(nullptr),
    fSoftDropTrue(nullptr),
    fNSubMeasured(nullptr),
    fNSubTrue(nullptr),
    fKineRec(nullptr),
    fKineSim(nullptr),
    fJetStructureMeasured(nullptr),
    fJetStructureTrue(nullptr),
    fQAHistos(nullptr),
    fLumiMonitor(nullptr),
    fSDZCut(0.1),
    fSDBetaCut(0),
    fReclusterizer(kCAAlgo),
    fHasRecEvent(false),
    fHasTrueEvent(false),
    fTriggerSelectionBits(AliVEvent::kAny),
    fTriggerSelectionString(""),
    fNameTriggerDecisionContainer("EmcalTriggerDecision"),
    fUseTriggerSelectionForData(false),
    fUseDownscaleWeight(false),
    fUseChargedConstituents(true),
    fUseNeutralConstituents(true),
    fFillPart(true),
    fFillRho(true),
    fFillSoftDrop(true),
    fFillNSub(true),
    fFillStructGlob(true)
{
  DefineOutput(2, TTree::Class());
  SetUseAliAnaUtils(true);
}

AliAnalysisTaskEmcalJetSubstructureTree::~AliAnalysisTaskEmcalJetSubstructureTree() { 
  if(fGlobalTreeParams) delete fGlobalTreeParams;
  if(fSoftDropMeasured) delete fSoftDropMeasured;
  if(fSoftDropTrue) delete fSoftDropTrue;
  if(fNSubMeasured) delete fNSubMeasured;
  if(fNSubTrue) delete fNSubTrue;
  if(fKineRec) delete fKineRec;
  if(fKineSim) delete fKineSim;
  if(fJetStructureMeasured) delete fJetStructureMeasured;
  if(fJetStructureTrue) delete fJetStructureTrue;
}

void AliAnalysisTaskEmcalJetSubstructureTree::UserCreateOutputObjects() {
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  // Make QA for constituent clusters
  TLinearBinning jetptbinning(9, 20, 200),
                 clusterenergybinning(200, 0., 200),
                 cellenergybinning(1000, 0., 100),
                 timebinning(1000, -500., 500.),
                 m02binning(100, 0., 1.),
                 ncellbinning(101, -0.5, 100.5),
                 exoticsbinning(2, -0.5, 1.5);
  fQAHistos = new THistManager("QAhistos");
  fQAHistos->CreateTH1("hEventCounter", "Event counter", 1, 0.5, 1.5);
  fQAHistos->CreateTH1("hTriggerClusterCounter", "Event counter separating into trigger clusters", 7, -1.5, 5.5);
  fQAHistos->CreateTH2("hClusterConstE", "EMCAL cluster energy vs jet pt; p_{t, jet} (GeV/c); E_{cl} (GeV)", jetptbinning, clusterenergybinning);
  fQAHistos->CreateTH2("hClusterConstTime", "EMCAL cluster time vs. jet pt; p_{t, jet} (GeV/c); t_{cl} (ns)", jetptbinning, timebinning);
  fQAHistos->CreateTH2("hClusterConstM02", "EMCAL cluster M02 vs. jet pt; p{t, jet} (GeV/c); M02", jetptbinning, m02binning);
  fQAHistos->CreateTH2("hClusterConstNcell", "EMCAL cluster ncell vs. jet pt; p{t, jet} (GeV/c); Number of cells", jetptbinning, ncellbinning);
  fQAHistos->CreateTH2("hClusterConstExotics", "EMCAL cluster exotics cut vs jet pt; p{t, jet} (GeV/c); Cluster exotics", jetptbinning, exoticsbinning);
  fQAHistos->CreateTH2("hClusterConstMinCellEnergy", "EMCAL Cluster const min cell energy; p{t, jet} (GeV/c); E_{cell} (GeV/c)", jetptbinning, cellenergybinning);
  fQAHistos->CreateTH2("hClusterConstMaxCellEnergy", "EMCAL Cluster const max (seed) cell energy; p{t, jet} (GeV/c); E_{cell} (GeV/c)", jetptbinning, cellenergybinning);

  // Test of constituent QA
#ifdef EXPERIMENTAL_JETCONSTITUENTS
  fQAHistos->CreateTH2("hChargedConstituentPt", "charged constituent pt vs jet pt (via constituent map); p_{t,jet} (GeV/c); p_{t,ch} (GeV/c)", jetptbinning, clusterenergybinning);
  fQAHistos->CreateTH2("hChargedIndexPt", "charged constituent pt vs jet pt (via index map); p_{t, jet} (GeV/c); p_{t, ch} (GeV/c)", jetptbinning, clusterenergybinning);

 fQAHistos->CreateTH2("hClusterConstituentEDefault", "cluster constituent default energy vs. jet pt (va constituent map); p_{t, jet} (GeV/c); E_{cl} (GeV)", jetptbinning, clusterenergybinning);
 fQAHistos->CreateTH2("hClusterConstituentENLC", "cluster constituent non-linearity-corrected energy vs. jet pt (va constituent map); p_{t, jet} (GeV/c); E_{cl} (GeV)", jetptbinning, clusterenergybinning);
 fQAHistos->CreateTH2("hClusterConstituentEHC", "cluster constituent hadronic-corrected energy vs. jet pt (va constituent map); p_{t, jet} (GeV/c); E_{cl} (GeV)", jetptbinning, clusterenergybinning);
 fQAHistos->CreateTH2("hClusterIndexENLC", "cluster constituent non-linearity-corrected energy vs. jet pt (via index map); p_{t, jet} (GeV/c); E_{cl} (GeV)", jetptbinning, clusterenergybinning);
 fQAHistos->CreateTH2("hClusterIndexEHC", "cluster constituent hadronic-corrected energy vs. jet pt (via index map); p_{t, jet} (GeV/c); E_{cl} (GeV)", jetptbinning, clusterenergybinning);
 fQAHistos->CreateTH2("hLeadingChargedConstituentPt", "Pt of the leading charged constituent in jet; p_{t,jet} (GeV/c); p_{t,ch} (GeV/c)", jetptbinning, clusterenergybinning);
 fQAHistos->CreateTH2("hLeadingClusterConstituentPt", "Pt of the leading cluster constituent in jet; p_{t,jet} (GeV/c); p_{t,ch} (GeV/c)", jetptbinning, clusterenergybinning);
#endif
  for(auto h : *(fQAHistos->GetListOfHistograms())) fOutput->Add(h);

  OpenFile(2);
  std::string treename = this->GetOutputSlot(2)->GetContainer()->GetName();
  fJetSubstructureTree = new TTree(treename.data(), "Tree with jet substructure information");

  fGlobalTreeParams = new AliJetTreeGlobalParameters;
  fGlobalTreeParams->LinkJetTreeBranches(fJetSubstructureTree, fFillRho);
  fKineRec = new AliJetKineParameters;
  fKineRec->LinkJetTreeBranches(fJetSubstructureTree, "Rec");
  if(fFillPart) {
    fKineSim = new AliJetKineParameters;
    fKineSim->LinkJetTreeBranches(fJetSubstructureTree, "Sim");
  }
  if(fFillSoftDrop) {
    fSoftDropMeasured = new AliSoftDropParameters;
    fSoftDropMeasured->LinkJetTreeBranches(fJetSubstructureTree, "Measured");
    if(fFillPart) {
      fSoftDropTrue = new AliSoftDropParameters;
      fSoftDropTrue->LinkJetTreeBranches(fJetSubstructureTree, "True");
    }
  }
  if(fFillNSub) {
    fNSubMeasured = new AliNSubjettinessParameters;
    fNSubMeasured->LinkJetTreeBranches(fJetSubstructureTree, "Measured");
    if(fFillPart) {
      fNSubTrue = new AliNSubjettinessParameters;
      fNSubTrue->LinkJetTreeBranches(fJetSubstructureTree, "True");
    }
  }

  if(fFillStructGlob){
    fJetStructureMeasured = new AliJetStructureParameters;
    fJetStructureMeasured->LinkJetTreeBranches(fJetSubstructureTree, "Measured");
    if(fFillPart){
      fJetStructureTrue = new AliJetStructureParameters;
      fJetStructureTrue->LinkJetTreeBranches(fJetSubstructureTree, "True");
    } 
  }

  PostData(1, fOutput);
  PostData(2, fJetSubstructureTree);
}

void AliAnalysisTaskEmcalJetSubstructureTree::RunChanged(Int_t newrun) {
  if(fUseDownscaleWeight){
    PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance()->SetRun(newrun);
  }
}

bool AliAnalysisTaskEmcalJetSubstructureTree::Run(){
  AliClusterContainer *clusters = GetClusterContainer(AliEmcalAnalysisFactory::ClusterContainerNameFactory(fInputEvent->IsA() == AliAODEvent::Class()));
  AliTrackContainer *tracks = GetTrackContainer(AliEmcalAnalysisFactory::TrackContainerNameFactory(fInputEvent->IsA() == AliAODEvent::Class()));
  AliParticleContainer *particles = GetParticleContainer("mcparticles");

  AliJetContainer *mcjets = GetJetContainer("mcjets");
  AliJetContainer *datajets = GetJetContainer("datajets");

  FillLuminosity(); // Makes only sense in data

  // for(auto e : *(fInputEvent->GetList())) std::cout << e->GetName() << std::endl;

  std::stringstream rhoTagData, rhoTagMC;
  if(datajets) rhoTagData << "R" << std::setw(2) << std::setfill('0') << static_cast<Int_t>(datajets->GetJetRadius() * 10.);
  if(mcjets) rhoTagMC << "R" << std::setw(2) << std::setfill('0') << static_cast<Int_t>(mcjets->GetJetRadius() * 10.);

  if(fFillRho){
    std::string rhoSparseData = "RhoSparse_Full_" + rhoTagData.str(), rhoSparseMC = "RhoSparse_Full_" + rhoTagMC.str(), 
                rhoMassData = "RhoMassSparse_Full_" + rhoTagData.str(), rhoMassMC = "RhoMassSparse_Full_" + rhoTagMC.str();
    AliRhoParameter *rhoPtRec = GetRhoFromEvent(rhoSparseData.data()),
                    *rhoMassRec = GetRhoFromEvent(rhoMassData.data()),
                    *rhoPtSim = GetRhoFromEvent(rhoSparseMC.data()),
                    *rhoMassSim = GetRhoFromEvent(rhoMassMC.data());
    AliDebugStream(2) << "Found rho parameter for reconstructed pt:    " << (rhoPtRec ? "yes" : "no") << ", value: " << (rhoPtRec ? rhoPtRec->GetVal() : 0.) << std::endl;
    AliDebugStream(2) << "Found rho parameter for sim pt:              " << (rhoPtSim ? "yes" : "no") << ", value: " << (rhoPtSim ? rhoPtSim->GetVal() : 0.) << std::endl;
    AliDebugStream(2) << "Found rho parameter for reconstructed Mass:  " << (rhoMassRec ? "yes" : "no") << ", value: " << (rhoMassRec ? rhoMassRec->GetVal() : 0.) << std::endl;
    AliDebugStream(2) << "Found rho parameter for sim Mass:            " << (rhoMassSim ? "yes" : "no") << ", value: " << (rhoMassSim ? rhoMassSim->GetVal() : 0.) << std::endl;
    Double_t rhopars[4] = {
                            rhoPtRec ? rhoPtRec->GetVal() : 0., 
                            rhoPtSim ? rhoPtSim->GetVal() : 0., 
                            rhoMassRec ? rhoMassRec->GetVal() : 0., 
                            rhoMassSim ? rhoMassSim->GetVal() : 0.
                          };
    memcpy(this->fGlobalTreeParams->fRhoParamters, rhopars, sizeof(Double_t) * 4);
  }

  AliDebugStream(1) << "Inspecting jet radius " << (datajets ? datajets->GetJetRadius() : mcjets->GetJetRadius()) << std::endl;
  this->fGlobalTreeParams->fJetRadius = (datajets ? datajets->GetJetRadius() : mcjets->GetJetRadius());
  fGlobalTreeParams->fTriggerClusterIndex = -1;       // Reset trigger cluster index

  if(datajets && !mcjets){
    // decode trigger string in order to determine the trigger clusters
    std::vector<std::string> clusternames;
    auto triggerinfos = PWG::EMCAL::Triggerinfo::DecodeTriggerString(fInputEvent->GetFiredTriggerClasses().Data());
    for(auto t : triggerinfos) {
      if(std::find(clusternames.begin(), clusternames.end(), t.Triggercluster()) == clusternames.end()) clusternames.emplace_back(t.Triggercluster());
    }
    bool isCENT = (std::find(clusternames.begin(), clusternames.end(), "CENT") != clusternames.end()),
         isCENTNOTRD = (std::find(clusternames.begin(), clusternames.end(), "CENTNOTRD") != clusternames.end()),
         isCALO = (std::find(clusternames.begin(), clusternames.end(), "CALO") != clusternames.end()),
         isCALOFAST = (std::find(clusternames.begin(), clusternames.end(), "CALOFAST") != clusternames.end());
    if(isCENT || isCENTNOTRD) {
      if(isCENT && isCENTNOTRD) fGlobalTreeParams->fTriggerClusterIndex = 0;    // CENTBOTH
      else if(isCENT) fGlobalTreeParams->fTriggerClusterIndex = 1;              // OnlyCENT
      else fGlobalTreeParams->fTriggerClusterIndex = 2;                         // OnlyCENTNOTRD
    }
    if(isCALO || isCALOFAST){
      // CALO(FAST) and CENT(NOTRD) clusters disjunct - CALO cluster not read out in parallel with CENT cluster
      // Therefore mixing between CALO and CENT triggers doesn't need to be handled.
      if(isCALO && isCALOFAST) fGlobalTreeParams->fTriggerClusterIndex = 3;
      else if(isCALO) fGlobalTreeParams->fTriggerClusterIndex = 4;
      else if(isCALOFAST) fGlobalTreeParams->fTriggerClusterIndex = 5;
    }
  }

  double weight = 1.;
  if(fUseDownscaleWeight){
    AliDebugStream(2) << "Trigger selection string: " << fTriggerSelectionString << std::endl;
    TString selectionString = (fTriggerSelectionBits & AliVEvent::kINT7) ? "INT7" : fTriggerSelectionString;
    auto triggerstring = MatchTrigger(selectionString.Data(), fTriggerSelectionString.Data());
    AliDebugStream(2) << "Getting downscale correction factor for trigger string " << triggerstring << std::endl;
    weight = 1./PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance()->GetDownscaleFactorForTriggerClass(triggerstring);
  }
  AliDebugStream(1) << "Using downscale weight " << weight << std::endl;
  this->fGlobalTreeParams->fEventWeight = weight;


  // Count events (for spectrum analysis)
  fQAHistos->FillTH1("hEventCounter", 1);
  fQAHistos->FillTH1("hTriggerClusterCounter", fGlobalTreeParams->fTriggerClusterIndex);

  AliSoftdropDefinition softdropSettings;
  softdropSettings.fBeta = fSDBetaCut;
  softdropSettings.fZ = fSDZCut;
  softdropSettings.fR0 = datajets->GetJetRadius();
  switch(fReclusterizer) {
  case kCAAlgo: softdropSettings.fRecluserAlgo = fastjet::cambridge_aachen_algorithm; break;
  case kKTAlgo: softdropSettings.fRecluserAlgo = fastjet::kt_algorithm; break;
  case kAKTAlgo: softdropSettings.fRecluserAlgo = fastjet::antikt_algorithm; break;
  };

  AliNSubjettinessDefinition nsubjettinessSettings;
  nsubjettinessSettings.fBeta = 1.;
  nsubjettinessSettings.fRadius = 0.4;

  if(datajets) {
    AliDebugStream(1) << "In data jets branch: found " <<  datajets->GetNJets() << " jets, " << datajets->GetNAcceptedJets() << " were accepted\n";
    AliDebugStream(1) << "Having MC information: " << (mcjets ? TString::Format("yes, with %d jets", mcjets->GetNJets()) : "no") << std::endl; 
    if(mcjets) {
      AliDebugStream(1) << "In MC jets branch: found " << mcjets->GetNJets() << " jets, " << mcjets->GetNAcceptedJets() << " were accepted\n";
    }
    for(auto jet : datajets->accepted()) {
      double pt = jet->Pt(), pz = jet->Pz(), E = jet->E(), M = TMath::Sqrt(E*E - pt*pt - pz*pz);
      AliDebugStream(2) << "Next jet: pt:" << jet->Pt() << ", E: " << E << ", pz: " << pz << ", M(self): " << M << "M(fj)" << jet->M() << std::endl;
      AliEmcalJet *associatedJet = jet->ClosestJet();

      if(mcjets) {
        if(!associatedJet) {
          AliDebugStream(2) << "Not found associated jet" << std::endl;
          continue;
        }
        if(!(SelectJet(*jet, tracks) && SelectJet(*associatedJet, particles))) continue;
        try {
          DoConstituentQA(jet, tracks, clusters);
          AliJetSubstructureData structureData =  MakeJetSubstructure(*jet, datajets->GetJetRadius() * 2., tracks, clusters, {softdropSettings, nsubjettinessSettings}),
                                 structureMC = fFillPart ? MakeJetSubstructure(*associatedJet, mcjets->GetJetRadius() * 2, particles, nullptr, {softdropSettings, nsubjettinessSettings}) : AliJetSubstructureData();
          if(fKineRec) *fKineRec = MakeJetKineParameters(*jet, kDetLevel, tracks, clusters);
          if(fKineSim) *fKineSim = MakeJetKineParameters(*associatedJet, kPartLevel, particles, nullptr);
          if(fSoftDropMeasured) *fSoftDropMeasured = structureData.fSoftDrop;
          if(fSoftDropTrue) *fSoftDropTrue = structureMC.fSoftDrop;
          if(fNSubMeasured) *fNSubMeasured = structureData.fNsubjettiness;
          if(fNSubTrue) *fNSubTrue = structureMC.fNsubjettiness;
          if(fJetStructureMeasured) *fJetStructureMeasured = {MakeAngularity(*jet, tracks, clusters), MakePtD(*jet, tracks, clusters)};
          if(fJetStructureTrue) *fJetStructureTrue = {MakeAngularity(*associatedJet, particles, nullptr), MakePtD(*associatedJet, particles, nullptr)};
          fJetSubstructureTree->Fill();
        } catch(ReclusterizerException &e) {
          AliErrorStream() << "Error in reclusterization - skipping jet" << std::endl;
        } catch(SubstructureException &e) {
          AliErrorStream() << "Error in substructure observable - skipping jet" << std::endl;
        }
      } else {
        if(!SelectJet(*jet, tracks)) continue;
        try {
          DoConstituentQA(jet, tracks, clusters);
          AliJetSubstructureData structure = MakeJetSubstructure(*jet, 0.4, tracks, clusters, {softdropSettings, nsubjettinessSettings});
          if(fKineRec) *fKineRec = MakeJetKineParameters(*jet, kDetLevel, tracks, clusters);
          if(fSoftDropMeasured) *fSoftDropMeasured = structure.fSoftDrop;
          if(fNSubMeasured) *fNSubMeasured = structure.fNsubjettiness;
          if(fJetStructureMeasured) *fJetStructureMeasured = {MakeAngularity(*jet, tracks, clusters), MakePtD(*jet, tracks, clusters)};
          fJetSubstructureTree->Fill();
        } catch(ReclusterizerException &e) {
          AliErrorStream() << "Error in reclusterization - skipping jet" << std::endl;
        } catch(SubstructureException &e) {
          AliErrorStream() << "Error in substructure observable - skipping jet" << std::endl;
        }
      }
    }
  } else {
    if(mcjets) {
      // for MCgen train
      AliDebugStream(1) << "In MC pure jet branch: found " << mcjets->GetNJets() << " jets, " << mcjets->GetNAcceptedJets() << " were accepted\n";
      for(auto j : mcjets->accepted()){
        AliEmcalJet *mcjet = static_cast<AliEmcalJet *>(j);
        try {
          AliJetSubstructureData structure = MakeJetSubstructure(*mcjet, mcjets->GetJetRadius() * 2., particles, nullptr,{softdropSettings, nsubjettinessSettings});
          if(this->fKineSim) *fKineSim = MakeJetKineParameters(*mcjet, kPartLevel, particles, nullptr);
          if(fSoftDropTrue) *fSoftDropTrue = structure.fSoftDrop;
          if(fNSubTrue) *fNSubTrue = structure.fNsubjettiness;
          if(fJetStructureTrue) *fJetStructureTrue = {MakeAngularity(*mcjet, particles, nullptr), MakePtD(*mcjet, particles, nullptr)};
          fJetSubstructureTree->Fill();
        } catch (ReclusterizerException &e) {
          AliErrorStream() << "Error in reclusterization - skipping jet" << std::endl;
        } catch (SubstructureException &e) {
          AliErrorStream() << "Error in substructure observable - skipping jet" << std::endl;
        }
      }
    }
  }

  return true;
}

Bool_t AliAnalysisTaskEmcalJetSubstructureTree::IsTriggerSelected(){
  AliDebugStream(1) << "Trigger selection called\n";
  if(!(fHasRecEvent || fHasTrueEvent)){
    AliErrorStream() << "Impossible combination: Neither rec nor true event available. Rejecting ..." << std::endl;
    return false;
  }
  // Run trigger selection (not on pure MCgen train - pure MCgen train has no rec event, ESD event is fake there)
  if(fHasRecEvent){
    if(!fHasTrueEvent){
      // Pure data - do EMCAL trigger selection from selection string
      AliDebugStream(1) << "Applying trigger selection for trigger bits " << std::bitset<sizeof(decltype(fTriggerSelectionBits)) * 8>(fTriggerSelectionBits) << "and trigger selection string " << fTriggerSelectionString << std::endl;
      if(!(fInputHandler->IsEventSelected() & fTriggerSelectionBits)) return false;
      AliDebugStream(1) << "Passed trigger bit selection" << std::endl;
      if(fTriggerSelectionString.Length()) {
        if(!fInputEvent->GetFiredTriggerClasses().Contains(fTriggerSelectionString)) return false;
        AliDebugStream(1) << "Passed trigger string section" << std::endl;
        if(fTriggerSelectionString.Contains("EJ") && fUseTriggerSelectionForData) {
          auto trgselresult = static_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer *>(fInputEvent->FindListObject(fNameTriggerDecisionContainer));
          AliDebugStream(1) << "Found trigger decision object: " << (trgselresult ? "yes" : "no") << std::endl;
          if(!trgselresult){
            AliErrorStream() <<  "Trigger decision container with name " << fNameTriggerDecisionContainer << " not found in event - not possible to select EMCAL triggers" << std::endl;
            return false;
          }
          if(!trgselresult->IsEventSelected(fTriggerSelectionString)) return false;
          AliDebugStream(1) << "Data event selected" << std::endl;
        }
      }
    } else {
      // Simulation - do EMCAL trigger selection from trigger selection object
      if(!(fInputHandler->IsEventSelected() & AliVEvent::kINT7)) return false;        // Require INT7 trigger - EMCAL triggers will be a subset
      if(IsSelectEmcalTriggers(fTriggerSelectionString.Data())){
        auto mctrigger = static_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer *>(fInputEvent->FindListObject(fNameTriggerDecisionContainer));
        AliDebugStream(1) << "Found trigger decision object: " << (mctrigger ? "yes" : "no") << std::endl;
        if(!mctrigger){
          AliErrorStream() <<  "Trigger decision container with name " << fNameTriggerDecisionContainer << " not found in event - not possible to select EMCAL triggers" << std::endl;
          return false;
        }
        if(!mctrigger->IsEventSelected(fTriggerSelectionString)) return false;
      }
    }
  }
  return true;      // trigger selected or pure MCgen information
}

void AliAnalysisTaskEmcalJetSubstructureTree::UserExecOnce() {
  AliCDBManager * cdb = AliCDBManager::Instance();
  if(!fMCEvent && cdb){
    // Get List of trigger clusters
    AliCDBEntry *en = cdb->Get("GRP/CTP/Config");
    AliTriggerConfiguration *trg = static_cast<AliTriggerConfiguration *>(en->GetObject());
    std::vector<std::string> clusternames;
    for(auto c : trg->GetClusters()) {
      AliTriggerCluster *clust = static_cast<AliTriggerCluster *>(c);
      std::string clustname = clust->GetName();
      auto iscent = clustname.find("CENT") != std::string::npos, iscalo = clustname.find("CALO") != std::string::npos; 
      if(!(iscalo || iscent)) continue;
      AliInfoStream() << "Adding trigger cluster " << clustname << " to cluster lumi monitor" << std::endl;
      clusternames.emplace_back(clustname);
   }

    // Set the x-axis of the luminosity monitor histogram
    fLumiMonitor = new TH1F("hLumiMonitor", "Luminosity monitor", clusternames.size(), 0, clusternames.size());
    int currentbin(1);
    for(auto c : clusternames) {
      fLumiMonitor->GetXaxis()->SetBinLabel(currentbin++, c.data());
    }
    fOutput->Add(fLumiMonitor);  
  }
}

void AliAnalysisTaskEmcalJetSubstructureTree::FillLuminosity() {
  if(fLumiMonitor && fUseDownscaleWeight){
    auto downscalefactors = PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance();
    if(fInputEvent->GetFiredTriggerClasses().Contains("INT7")) {
      for(auto trigger : PWG::EMCAL::Triggerinfo::DecodeTriggerString(fInputEvent->GetFiredTriggerClasses().Data())){
        auto int7trigger = trigger.IsTriggerClass("INT7");
        auto bunchcrossing = trigger.BunchCrossing() == "B";
        auto nopf = trigger.PastFutureProtection() == "NOPF";
        bool centcalo = (trigger.Triggercluster().find("CENT") != std::string::npos) || (trigger.Triggercluster().find("CALO") != std::string::npos);
        AliDebugStream(4) << "Full name: " << trigger.ExpandClassName() << ", INT7 trigger:  " << (int7trigger ? "Yes" : "No") << ", bunch crossing: " << (bunchcrossing ? "Yes" : "No") << ", no past-future protection: " << (nopf ? "Yes" : "No")  << ", Cluster: " << trigger.Triggercluster() << std::endl;
        if(int7trigger && bunchcrossing && nopf && centcalo) {
          double downscale = downscalefactors->GetDownscaleFactorForTriggerClass(trigger.ExpandClassName());
          AliDebugStream(5) << "Using downscale " << downscale << std::endl;
          fLumiMonitor->Fill(trigger.Triggercluster().data(), 1./downscale);
        }
      }
    }
  }
}

AliJetKineParameters AliAnalysisTaskEmcalJetSubstructureTree::MakeJetKineParameters(const AliEmcalJet &jet, JetRecType_t rectype, const AliParticleContainer *const tracks, const AliClusterContainer *const clusters) const {
  AliJetKineParameters result;
  result.fPt = TMath::Abs(jet.Pt());
  result.fE = jet.E();
  result.fEta = jet.Eta();
  result.fPhi = jet.Phi();
  result.fArea = jet.Area();
  result.fMass = jet.M();
  result.fNEF = jet.NEF();
  result.fNCharged = jet.GetNumberOfTracks();
  result.fNNeutral = jet.GetNumberOfClusters();
  std::vector<double> zcharged, zneutral;
  if(tracks) {
    // Find the leading track
    for(auto icharged = 0; icharged < jet.GetNumberOfTracks(); icharged++){
      auto trk = jet.TrackAt(icharged, tracks->GetArray());
      bool charged = true;
      if(rectype == kPartLevel){
        if(!trk->Charge()) charged = false;
      }
      auto z = jet.GetZ(trk);
      if(charged) zcharged.push_back(z);
      else zneutral.push_back(z);
    }
  } 
  if(clusters) {
    for(auto iclust = 0; iclust < jet.GetNumberOfClusters(); iclust++){
      auto clust = jet.ClusterAt(iclust, clusters->GetArray());
      TLorentzVector clustervec;
      clust->GetMomentum(clustervec, fVertex, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy());
      auto z = jet.GetZ(clustervec.Px(), clustervec.Py(), clustervec.Pz());
      zneutral.push_back(z);
    }
  }
  result.fZLeading = -1.;
  result.fZLeadingCharged = -1.;
  result.fZLeadingNeutral = -1;
  if(zcharged.size()) {
    std::sort(zcharged.begin(), zcharged.end(), std::greater<double>());
    result.fZLeadingCharged = zcharged[0];
    result.fZLeading = result.fZLeadingCharged;
  }
  if(zneutral.size()){
    std::sort(zneutral.begin(), zneutral.end(), std::greater<double>());
    result.fZLeadingNeutral = zneutral[0];
    if(result.fZLeadingNeutral > result.fZLeading) result.fZLeading = result.fZLeadingNeutral;
  }
  return result;
}

AliJetSubstructureData AliAnalysisTaskEmcalJetSubstructureTree::MakeJetSubstructure(const AliEmcalJet &jet, double jetradius, const AliParticleContainer *tracks, const AliClusterContainer *clusters, const AliJetSubstructureSettings &settings) const {
  const int kClusterOffset = 30000; // In order to handle tracks and clusters in the same index space the cluster index needs and offset, large enough so that there is no overlap with track indices
  std::vector<fastjet::PseudoJet> constituents;
  bool isMC = dynamic_cast<const AliMCParticleContainer *>(tracks);
  AliDebugStream(2) << "Make new jet substrucutre for " << (isMC ? "MC" : "data") << " jet: Number of tracks " << jet.GetNumberOfTracks() << ", clusters " << jet.GetNumberOfClusters() << std::endl;
  if(tracks && (fUseChargedConstituents || isMC)){                    // Neutral particles part of particle container in case of MC
    AliDebugStream(1) << "Jet substructure: Using charged constituents" << std::endl;
    for(int itrk = 0; itrk < jet.GetNumberOfTracks(); itrk++){
      auto track = jet.TrackAt(itrk, tracks->GetArray());
      if(!track->Charge() && !fUseNeutralConstituents) continue;      // Reject neutral constituents in case of using only charged consituents
      if(track->Charge() && !fUseChargedConstituents) continue;       // Reject charged constituents in case of using only neutral consituents
      fastjet::PseudoJet constituentTrack(track->Px(), track->Py(), track->Pz(), track->E());
      constituentTrack.set_user_index(jet.TrackAt(itrk));
      constituents.push_back(constituentTrack);
    }
  }

  if(clusters && fUseNeutralConstituents){
    AliDebugStream(1) << "Jet substructure: Using neutral constituents" << std::endl;
    for(int icl = 0; icl < jet.GetNumberOfClusters(); icl++) {
      auto cluster = jet.ClusterAt(icl, clusters->GetArray());
      TLorentzVector clustervec;
      cluster->GetMomentum(clustervec, fVertex, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy());
      fastjet::PseudoJet constituentCluster(clustervec.Px(), clustervec.Py(), clustervec.Pz(), cluster->GetHadCorrEnergy());
      constituentCluster.set_user_index(jet.ClusterAt(icl) + kClusterOffset);
      constituents.push_back(constituentCluster);
    }
  }

  AliDebugStream(3) << "Found " << constituents.size() << " constituents for jet with pt=" << jet.Pt() << " GeV/c" << std::endl;
  if(!constituents.size()){
    AliErrorStream() << "Jet has 0 constituents." << std::endl;
    throw ReclusterizerException();
  }
  // Redo jet finding on constituents with a
  fastjet::JetDefinition jetdef(fastjet::antikt_algorithm, jetradius*2, static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30 );
  std::vector<fastjet::PseudoJet> outputjets;
  try {
    fastjet::ClusterSequence jetfinder(constituents, jetdef);
    outputjets = jetfinder.inclusive_jets(0);
    AliJetSubstructureData result({fFillSoftDrop ? MakeSoftDropParameters(outputjets[0], settings.fSoftdropSettings) : AliSoftDropParameters(), fFillNSub ? MakeNsubjettinessParameters(outputjets[0], settings.fSubjettinessSettings): AliNSubjettinessParameters()});
    return result;
  } catch (fastjet::Error &e) {
    AliErrorStream() << " FJ Exception caught: " << e.message() << std::endl;
    throw ReclusterizerException();
  } catch (SoftDropException &e) {
    AliErrorStream() << "Softdrop exception caught: " << e.what() << std::endl;
    throw ReclusterizerException();
  }
}

AliSoftDropParameters AliAnalysisTaskEmcalJetSubstructureTree::MakeSoftDropParameters(const fastjet::PseudoJet &jet, const AliSoftdropDefinition &cutparameters) const {
  fastjet::contrib::SoftDrop softdropAlgorithm(cutparameters.fBeta, cutparameters.fZ, cutparameters.fR0);
  softdropAlgorithm.set_verbose_structure(kTRUE);
#if FASTJET_VERSION_NUMBER >= 30302
  fastjet::Recluster reclusterizer(cutparameters.fRecluserAlgo, 1, fastjet::Recluster::keep_only_hardest);
#else
  fastjet::contrib::Recluster reclusterizer(cutparameters.fRecluserAlgo, 1, true);
#endif
  softdropAlgorithm.set_reclustering(kTRUE, &reclusterizer);
  AliDebugStream(4) << "Jet has " << jet.constituents().size() << " constituents" << std::endl;
  auto groomed = softdropAlgorithm(jet);
  try {
    auto softdropstruct = groomed.structure_of<fastjet::contrib::SoftDrop>();

    AliSoftDropParameters result({softdropstruct.symmetry(),
                                  groomed.m(),
                                  softdropstruct.delta_R(),
                                  groomed.perp(),
                                  softdropstruct.delta_R(),
                                  softdropstruct.mu(),
                                  softdropstruct.dropped_count()});
    return result;
  } catch(std::bad_cast &e) {
    throw SoftDropException(); 
  }
}

AliNSubjettinessParameters AliAnalysisTaskEmcalJetSubstructureTree::MakeNsubjettinessParameters(const fastjet::PseudoJet &jet, const AliNSubjettinessDefinition &cut) const {
  AliNSubjettinessParameters result({
    fastjet::contrib::Nsubjettiness (1,fastjet::contrib::KT_Axes(),fastjet::contrib::NormalizedCutoffMeasure(cut.fBeta, cut.fRadius, 1e100)).result(jet),
    fastjet::contrib::Nsubjettiness (2,fastjet::contrib::KT_Axes(),fastjet::contrib::NormalizedCutoffMeasure(cut.fBeta, cut.fRadius, 1e100)).result(jet)
  });
  return result;
}

Double_t AliAnalysisTaskEmcalJetSubstructureTree::MakeAngularity(const AliEmcalJet &jet, const AliParticleContainer *tracks, const AliClusterContainer *clusters) const {
  if(!(jet.GetNumberOfTracks() || jet.GetNumberOfClusters()))
    throw SubstructureException();
  TVector3 jetvec(jet.Px(), jet.Py(), jet.Pz());
  Double_t den(0.), num(0.);
  bool isMC = dynamic_cast<const AliMCParticleContainer *>(tracks);
  if(tracks && (fUseChargedConstituents || isMC)){
    AliDebugStream(1) << "Angularity: Using charged constituents" << std::endl;
    for(UInt_t itrk = 0; itrk < jet.GetNumberOfTracks(); itrk++) {
      auto track = jet.TrackAt(itrk, tracks->GetArray());
      if(!track){
        AliErrorStream() << "Associated constituent particle / track not found\n";
        continue;
      }
      if(!track->Charge() && !fUseNeutralConstituents) continue;      // Reject neutral constituents in case of using only charged consituents
      if(track->Charge() && !fUseChargedConstituents) continue;       // Reject charged constituents in case of using only neutral consituents
      TVector3 trackvec(track->Px(), track->Py(), track->Pz());

      num +=  track->Pt() * trackvec.DrEtaPhi(jetvec);
      den += +track->Pt();
    }
  }
  if(clusters && fUseNeutralConstituents) {
    AliDebugStream(1) << "Using neutral constituents" << std::endl;
    for(UInt_t icl = 0; icl < jet.GetNumberOfClusters(); icl++){
      auto clust = jet.ClusterAt(icl, clusters->GetArray());
      if(!clust) {
        AliErrorStream() << "Associated constituent cluster not found\n";
        continue;
      }
      TLorentzVector clusterp;
      clust->GetMomentum(clusterp, fVertex, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy());

      num += clusterp.Pt() * clusterp.Vect().DrEtaPhi(jetvec);
      den += clusterp.Pt();
    }
  }
  return num/den;
}

Double_t AliAnalysisTaskEmcalJetSubstructureTree::MakePtD(const AliEmcalJet &jet, const AliParticleContainer *const particles, const AliClusterContainer *const clusters) const {
  if (!(jet.GetNumberOfTracks() || jet.GetNumberOfClusters()))
    throw SubstructureException();
  Double_t den(0.), num(0.);
  bool isMC = dynamic_cast<const AliMCParticleContainer *>(particles);
  if(particles && (fUseChargedConstituents || isMC)){
    AliDebugStream(1) << "Using charged constituents" << std::endl;
    for(UInt_t itrk = 0; itrk < jet.GetNumberOfTracks(); itrk++) {
      auto trk = jet.TrackAt(itrk, particles->GetArray());
      if(!trk){
        AliErrorStream() << "Associated constituent particle / track not found\n";
        continue;
      }
      if(!trk->Charge() && !fUseNeutralConstituents) continue;      // Reject neutral constituents in case of using only charged consituents
      if(trk->Charge() && !fUseChargedConstituents) continue;       // Reject charged constituents in case of using only neutral consituents
      num += trk->Pt() * trk->Pt();
      den += trk->Pt();
    }
  }
  if(clusters && fUseNeutralConstituents){
    AliDebugStream(1) << "Using neutral constituents" << std::endl;
    for(UInt_t icl = 0; icl < jet.GetNumberOfClusters(); icl++){
      auto clust = jet.ClusterAt(icl, clusters->GetArray());
      if(!clust) {
        AliErrorStream() << "Associated constituent cluster not found\n";
        continue;
      }
      TLorentzVector clusterp;
      clust->GetMomentum(clusterp, fVertex, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy());
      num += clusterp.Pt() * clusterp.Pt();
      den += clusterp.Pt();
    }
  }
  return TMath::Sqrt(num)/den;
}

void AliAnalysisTaskEmcalJetSubstructureTree::DoConstituentQA(const AliEmcalJet *jet, const AliParticleContainer *cont, const AliClusterContainer *clusters){
  for(int icl = 0; icl < jet->GetNumberOfClusters(); icl++){
    auto clust = jet->ClusterAt(icl, clusters->GetArray());
    AliDebugStream(3) << "cluster time " << clust->GetTOF() << std::endl;
    fQAHistos->FillTH2("hClusterConstE", jet->Pt(),clust->GetUserDefEnergy(clusters->GetDefaultClusterEnergy()));
    fQAHistos->FillTH2("hClusterConstTime", jet->Pt(), clust->GetTOF()*1e9);    // convert to nanoseconds
    fQAHistos->FillTH2("hClusterConstM02", jet->Pt(), clust->GetM02());
    fQAHistos->FillTH2("hClusterConstNcell", jet->Pt(), clust->GetNCells());
    fQAHistos->FillTH2("hClusterConstExotics", jet->Pt(), clust->GetIsExotic() ? 1. : 0.);

    double mincell(100000.), maxcell(0.);
    for(int icell = 0; icell < clust->GetNCells(); icell++){
      double ecell = clust->E() * clust->GetCellAmplitudeFraction(icell);
      if(ecell < mincell) mincell = ecell;
      if(ecell > maxcell) maxcell = ecell;
    }
    fQAHistos->FillTH2("hClusterConstMinCellEnergy", jet->Pt(), mincell);
    fQAHistos->FillTH2("hClusterConstMaxCellEnergy", jet->Pt(), maxcell);

#ifdef EXPERIMENTAL_JETCONSTITUENTS
    fQAHistos->FillTH2("hClusterIndexENLC", jet->Pt(), clust->GetNonLinCorrEnergy());
    fQAHistos->FillTH2("hClusterIndexEHC", jet->Pt(), clust->GetHadCorrEnergy());
#endif
  }

#ifdef EXPERIMENTAL_JETCONSTITUENTS
  // Loop over charged particles - fill test histogram
  for(int itrk = 0; itrk < jet->GetNumberOfTracks(); itrk++){
    auto part = jet->TrackAt(itrk, cont->GetArray());
    fQAHistos->FillTH2("hChargedIndexPt", jet->Pt(), part->Pt());
  }

  // Look over charged constituents
  AliDebugStream(2) << "Jet: Number of particle constituents: " << jet->GetParticleConstituents().size() << std::endl;
  for(auto part : jet->GetParticleConstituents()) {
    //auto part = static_cast<PWG::JETFW::AliEmcalParticleJetConstituent *>(pconst);
    AliDebugStream(3) << "Found particle constituent with pt " << part.Pt() << ", from VParticle " << part.GetParticle()->Pt() << std::endl;
    fQAHistos->FillTH2("hChargedConstituentPt", jet->Pt(), part.Pt());
  }

  // Loop over neutral constituents
  AliDebugStream(2) << "Jet: Number of cluster constituents: " << jet->GetClusterConstituents().size() << std::endl;
  for(auto clust : jet->GetClusterConstituents()){
    //auto clust = static_cast<PWG::JETFW::AliEmcalClusterJetConstituent *>(cconst);
    AliDebugStream(3) << "Found cluster constituent with energy " << clust.E() << " using energy definition " << static_cast<int>(clust.GetDefaultEnergyType()) << std::endl;
    fQAHistos->FillTH2("hClusterConstituentEDefault", jet->Pt(), clust.E());
    fQAHistos->FillTH2("hClusterConstituentENLC", jet->Pt(), clust.GetCluster()->GetNonLinCorrEnergy());
    fQAHistos->FillTH2("hClusterConstituentEHC", jet->Pt(), clust.GetCluster()->GetHadCorrEnergy());
  }

  // Fill global observables: Leading charged and cluster constituents
  auto leadingcharged = jet->GetLeadingParticleConstituent();
  auto leadingcluster = jet->GetLeadingClusterConstituent();
  if(leadingcluster){
    fQAHistos->FillTH1("hLeadingClusterConstituentPt", jet->Pt(), leadingcluster->GetCluster()->GetHadCorrEnergy());
  }
  if(leadingcharged) {
    fQAHistos->FillTH1("hLeadingChargedConstituentPt", jet->Pt(), leadingcharged->GetParticle()->Pt());
  }
#endif
}

bool AliAnalysisTaskEmcalJetSubstructureTree::SelectJet(const AliEmcalJet &jet, const AliParticleContainer *particles) const {
  int ncharged = 0, nneutral = jet.GetNumberOfClusters();
  if(particles) {
    for(decltype(jet.GetNumberOfTracks()) ipart = 0; ipart < jet.GetNumberOfTracks(); ipart++){
      auto part = jet.TrackAt(ipart, particles->GetArray());
      if(!part) continue;
      if(part->Charge()) ncharged++;
      else nneutral++;
   }
  }
  // check if the jet has at least one consituent for jet substructure
  int nallowed = 0;
  nallowed += fUseChargedConstituents ? ncharged : 0;
  nallowed += fUseNeutralConstituents ? nneutral : 0;
  return nallowed > 0;
}

AliAnalysisTaskEmcalJetSubstructureTree *AliAnalysisTaskEmcalJetSubstructureTree::AddEmcalJetSubstructureTreeMaker(Bool_t isMC, Bool_t isData, Double_t jetradius, AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recombinationScheme, Bool_t useDCAL, const char *trigger){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

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

  std::stringstream taskname;
  taskname << "JetSubstructureTreemaker_R" << std::setw(2) << std::setfill('0') << int(jetradius*10) << trigger;  
  AliAnalysisTaskEmcalJetSubstructureTree *treemaker = new AliAnalysisTaskEmcalJetSubstructureTree(taskname.str().data());
  mgr->AddTask(treemaker);
  treemaker->SetMakeGeneralHistograms(kTRUE);
  if(isMC) treemaker->SetHasTrueEvent(true);
  if(isData) treemaker->SetHasRecEvent(true);

  // Adding containers
  if(isMC) {
    AliParticleContainer *particles = treemaker->AddMCParticleContainer("mcparticles");
    particles->SetMinPt(0.);

    AliJetContainer *mcjets = treemaker->AddJetContainer(
                              jettype,
                              AliJetContainer::antikt_algorithm,
                              recombinationScheme,
                              jetradius,
                              (isData && ((jettype == AliJetContainer::kFullJet) || (jettype == AliJetContainer::kNeutralJet))) ? (useDCAL ? AliEmcalJet::kDCALfid : AliEmcalJet::kEMCALfid) : AliEmcalJet::kTPC,
                              particles, nullptr);
    mcjets->SetName("mcjets");
    mcjets->SetJetPtCut(20.);
  }

  if(isData) {
    AliTrackContainer *tracks(nullptr);
    if((jettype == AliJetContainer::kChargedJet) || (jettype == AliJetContainer::kFullJet)){
      tracks = treemaker->AddTrackContainer(AliEmcalAnalysisFactory::TrackContainerNameFactory(isAOD));
      std::cout << "Track container name: " << tracks->GetName() << std::endl;
      tracks->SetMinPt(0.15);
    }
    AliClusterContainer *clusters(nullptr);
    if((jettype == AliJetContainer::kFullJet) || (jettype == AliJetContainer::kNeutralJet)){
      std::cout << "Using full or neutral jets ..." << std::endl;
      clusters = treemaker->AddClusterContainer(AliEmcalAnalysisFactory::ClusterContainerNameFactory(isAOD));
      std::cout << "Cluster container name: " << clusters->GetName() << std::endl;
      clusters->SetClusHadCorrEnergyCut(0.3); // 300 MeV E-cut
      clusters->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
    } else {
      std::cout << "Using charged jets ... " << std::endl;
    }

    AliJetContainer *datajets = treemaker->AddJetContainer(
                              jettype,
                              AliJetContainer::antikt_algorithm,
                              recombinationScheme,
                              jetradius,
                              ((jettype == AliJetContainer::kFullJet) || (jettype == AliJetContainer::kNeutralJet)) ? (useDCAL ? AliEmcalJet::kDCALfid : AliEmcalJet::kEMCALfid): AliEmcalJet::kTPCfid,
                              tracks, clusters);
    datajets->SetName("datajets");
    datajets->SetJetPtCut(20.);

    treemaker->SetUseAliAnaUtils(true, true);
    treemaker->SetVzRange(-10., 10);

    // configure trigger selection
    std::string triggerstring(trigger);
    if(triggerstring.find("INT7") != std::string::npos) {
      treemaker->SetTriggerBits(AliVEvent::kINT7);
    } else if(triggerstring.find("EJ1") != std::string::npos) {
      treemaker->SetTriggerBits(AliVEvent::kEMCEJE);
      treemaker->SetTriggerString("EJ1");
    } else if(triggerstring.find("EJ2") != std::string::npos) {
      treemaker->SetTriggerBits(AliVEvent::kEMCEJE);
      treemaker->SetTriggerString("EJ2");
    } else if(triggerstring.find("EG1") != std::string::npos) {
      treemaker->SetTriggerBits(AliVEvent::kEMCEGA);
      treemaker->SetTriggerString("EG1");
    } else if(triggerstring.find("EG2") != std::string::npos) {
      treemaker->SetTriggerBits(AliVEvent::kEMCEGA);
      treemaker->SetTriggerString("EG2");
    }
  }
  
  std::string jettypestring;
  switch(jettype) {
    case AliJetContainer::kFullJet: jettypestring = "FullJets"; break;
    case AliJetContainer::kChargedJet: jettypestring = "ChargedJets"; break;
    case AliJetContainer::kNeutralJet: jettypestring = "NeutralJets"; break;
    default: jettypestring = "Undef";
  };

  // Connecting containers
  std::stringstream outputfile, histname, treename;
  outputfile << mgr->GetCommonFileName() << ":JetSubstructure_" << jettypestring << "_R" << std::setw(2) << std::setfill('0') << int(jetradius * 10.) << "_" << trigger;
  histname << "JetSubstructureHistos_" << jettypestring << "_R" << std::setw(2) << std::setfill('0') << int(jetradius * 10.) << "_" << trigger;
  treename << "JetSubstructureTree_" << jettypestring << "_R" << std::setw(2) << std::setfill('0') << int(jetradius * 10.) << "_" << trigger;
  mgr->ConnectInput(treemaker, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(treemaker, 1, mgr->CreateContainer(histname.str().data(), AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, outputfile.str().data()));
  mgr->ConnectOutput(treemaker, 2, mgr->CreateContainer(treename.str().data(), TTree::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName()));

  return treemaker;
}

void AliSoftDropParameters::LinkJetTreeBranches(TTree *jettree, const char *tag) {
  LinkBranch(jettree, &fZg, Form("Zg%s", tag), "D");
  LinkBranch(jettree, &fRg, Form("Rg%s", tag), "D");
  LinkBranch(jettree, &fMg, Form("Mg%s", tag), "D");
  LinkBranch(jettree, &fPtg, Form("Ptg%s", tag), "D");
  LinkBranch(jettree, &fMug, Form("Mug%s", tag), "D");
  LinkBranch(jettree, &fDeltaR, Form("DeltaRg%s", tag), "D");
  LinkBranch(jettree, &fNDropped, Form("NDropped%s", tag), "I");
};

void AliNSubjettinessParameters::LinkJetTreeBranches(TTree *jettree, const char *tag) {
  LinkBranch(jettree, &fOneSubjettiness, Form("OneSubjettiness%s", tag), "D");
  LinkBranch(jettree, &fTwoSubjettiness, Form("TwoSubjettiness%s", tag), "D");
}

void AliJetStructureParameters::LinkJetTreeBranches(TTree *jettree, const char *tag){
  LinkBranch(jettree, &fAngularity, Form("Angularity%s", tag), "D");
  LinkBranch(jettree, &fPtD, Form("PtD%s", tag), "D");
}

void AliJetKineParameters::LinkJetTreeBranches(TTree *jettree, const char *tag){
  LinkBranch(jettree, &fPt, Form("PtJet%s", tag), "D");
  LinkBranch(jettree, &fE, Form("EJet%s", tag), "D");
  LinkBranch(jettree, &fEta, Form("Eta%s", tag), "D");
  LinkBranch(jettree, &fPhi, Form("Phi%s", tag), "D");
  LinkBranch(jettree, &fArea, Form("Area%s", tag), "D");
  LinkBranch(jettree, &fMass, Form("Mass%s", tag), "D");
  LinkBranch(jettree, &fNEF, Form("NEF%s", tag), "D");
  LinkBranch(jettree, &fNCharged, Form("NCharged%s", tag), "I");
  LinkBranch(jettree, &fNNeutral, Form("NNeutral%s", tag), "I");
  LinkBranch(jettree, &fZLeading, Form("ZLeading%s", tag), "D");
  LinkBranch(jettree, &fZLeadingCharged, Form("ZLeadingCharged%s", tag), "D");
  LinkBranch(jettree, &fZLeadingNeutral, Form("ZLeadingNeutral%s", tag), "D");
}

void AliJetTreeGlobalParameters::LinkJetTreeBranches(TTree *jettree, bool fillRho) {
  LinkBranch(jettree, &fJetRadius, "Radius", "D");
  LinkBranch(jettree, &fEventWeight, "EventWeight", "D");
  LinkBranch(jettree, &fTriggerClusterIndex, "TriggerClusterIndex", "I");
  if(fillRho) {
    std::string varnames[] = {"RhoPtRec", "RhoPtSim", "RhoMassRec", "RhoMassSim"};
    for(int i = 0; i < 4; i++){
      LinkBranch(jettree, fRhoParamters + i, varnames[i].data(), "D");
    }
  }
}

void PWGJE::EMCALJetTasks::LinkBranch(TTree *jettree, void *data, const char *branchname, const char *type) {
  jettree->Branch(branchname, data, Form("%s/%s", branchname, type));
}
