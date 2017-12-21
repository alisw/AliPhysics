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
#include <iostream>
#include <string>
#include <set>
#include <sstream>
#include <vector>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/contrib/Nsubjettiness.hh>
#include <fastjet/contrib/SoftDrop.hh>

#include <THistManager.h>
#include <TLinearBinning.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TString.h>
#include <TVector3.h>

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisTaskEmcalJetSubstructureTree.h"
#include "AliClusterContainer.h"
#include "AliJetContainer.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliEmcalDownscaleFactorsOCDB.h"
#include "AliEmcalJet.h"
#include "AliEmcalList.h"
#include "AliEmcalTriggerDecisionContainer.h"
#include "AliLog.h"
#include "AliParticleContainer.h"
#include "AliTrackContainer.h"
#include "AliRhoParameter.h"
#include "AliVCluster.h"
#include "AliVParticle.h"

#ifdef EXPERIMENTAL_JETCONSTITUENTS
#include "AliEmcalClusterJetConstituent.h"
#include "AliEmcalParticleJetConstituent.h"
#endif


/// \cond CLASSIMP
ClassImp(EmcalTriggerJets::AliAnalysisTaskEmcalJetSubstructureTree);
/// \endcond

namespace EmcalTriggerJets {

AliAnalysisTaskEmcalJetSubstructureTree::AliAnalysisTaskEmcalJetSubstructureTree() :
    AliAnalysisTaskEmcalJet(),
    fJetSubstructureTree(nullptr),
    fQAHistos(nullptr),
    fSDZCut(0.1),
    fSDBetaCut(0),
    fReclusterizer(kCAAlgo),
    fTriggerSelectionBits(AliVEvent::kAny),
    fTriggerSelectionString(""),
    fUseDownscaleWeight(false),
    fFillPart(true),
    fFillAcceptance(true),
    fFillRho(true),
    fFillMass(true),
    fFillSoftDrop(true),
    fFillNSub(true),
    fFillStructGlob(true)
{
  memset(fJetTreeData, 0, sizeof(Double_t) * kTNVar);
}

AliAnalysisTaskEmcalJetSubstructureTree::AliAnalysisTaskEmcalJetSubstructureTree(const char *name) :
    AliAnalysisTaskEmcalJet(name, kTRUE),
    fJetSubstructureTree(nullptr),
    fQAHistos(nullptr),
    fSDZCut(0.1),
    fSDBetaCut(0),
    fReclusterizer(kCAAlgo),
    fTriggerSelectionBits(AliVEvent::kAny),
    fTriggerSelectionString(""),
    fUseDownscaleWeight(false),
    fFillPart(true),
    fFillAcceptance(true),
    fFillRho(true),
    fFillMass(true),
    fFillSoftDrop(true),
    fFillNSub(true),
    fFillStructGlob(true)
{
  memset(fJetTreeData, 0, sizeof(Double_t) * kTNVar);
  DefineOutput(2, TTree::Class());
}

AliAnalysisTaskEmcalJetSubstructureTree::~AliAnalysisTaskEmcalJetSubstructureTree() {

}

void AliAnalysisTaskEmcalJetSubstructureTree::UserCreateOutputObjects() {
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  // Make QA for constituent clusters
  TLinearBinning jetptbinning(9, 20, 200),
                 clusterenergybinning(200, 0., 200),
                 timebinning(1000, -500., 500.),
                 m02binning(100, 0., 1.),
                 ncellbinning(101, -0.5, 100.5);
  fQAHistos = new THistManager("QAhistos");
  fQAHistos->CreateTH2("hClusterConstE", "EMCAL cluster energy vs jet pt; p_{t, jet} (GeV/c); E_{cl} (GeV)", jetptbinning, clusterenergybinning);
  fQAHistos->CreateTH2("hClusterConstTime", "EMCAL cluster time vs. jet pt; p_{t, jet} (GeV/c); t_{cl} (ns)", jetptbinning, timebinning);
  fQAHistos->CreateTH2("hClusterConstM02", "EMCAL cluster M02 vs. jet pt; p{t, jet} (GeV/c); M02", jetptbinning, m02binning);
  fQAHistos->CreateTH2("hClusterConstNcell", "EMCAL cluster ncell vs. jet pt; p{t, jet} (GeV/c); Number of cells", jetptbinning, ncellbinning);

  // Test of constituent QA
#ifdef EXPERIMENTAL_JETCONSTITUENTS
  fQAHistos->CreateTH2("hChargedConstituentPt", "charged constituent pt vs jet pt (via constituent map); p_{t,jet} (GeV/c); p_{t,ch} (GeV/c)", jetptbinning, clusterenergybinning);
  fQAHistos->CreateTH2("hChargedIndexPt", "charged constituent pt vs jet pt (via index map); p_{t, jet} (GeV/c); p_{t, ch} (GeV/c)", jetptbinning, clusterenergybinning);

 fQAHistos->CreateTH2("hClusterConstituentEDefault", "cluster constituent default energy vs. jet pt (va constituent map); p_{t, jet} (GeV/c); E_{cl} (GeV)", jetptbinning, clusterenergybinning);
 fQAHistos->CreateTH2("hClusterConstituentENLC", "cluster constituent non-linearity-corrected energy vs. jet pt (va constituent map); p_{t, jet} (GeV/c); E_{cl} (GeV)", jetptbinning, clusterenergybinning);
 fQAHistos->CreateTH2("hClusterConstituentEHC", "cluster constituent hadronic-corrected energy vs. jet pt (va constituent map); p_{t, jet} (GeV/c); E_{cl} (GeV)", jetptbinning, clusterenergybinning);
 fQAHistos->CreateTH2("hClusterIndexENLC", "cluster constituent non-linearity-corrected energy vs. jet pt (via index map); p_{t, jet} (GeV/c); E_{cl} (GeV)", jetptbinning, clusterenergybinning);
 fQAHistos->CreateTH2("hClusterIndexEHC", "cluster constituent hadronic-corrected energy vs. jet pt (via index map); p_{t, jet} (GeV/c); E_{cl} (GeV)", jetptbinning, clusterenergybinning);
#endif
  for(auto h : *(fQAHistos->GetListOfHistograms())) fOutput->Add(h);

  OpenFile(2);
  TString treename = this->GetOutputSlot(2)->GetContainer()->GetName();
  fJetSubstructureTree = new TTree(treename, "Tree with jet substructure information");
  TString varnames[kTNVar];
  varnames[0] = "Radius";
  varnames[1] = "EventWeight";
  varnames[2] = "PtJetRec";
  varnames[3] = "PtJetSim";
  varnames[4] = "EJetRec";
  varnames[5] = "EJetSim";
  varnames[6] = "EtaRec";
  varnames[7] = "EtaSim";
  varnames[8] = "PhiRec";
  varnames[9] = "PhiSim";
  varnames[10] = "RhoPtRec";
  varnames[11] = "RhoPtSim";
  varnames[12] = "RhoMassRec";
  varnames[13] = "RhoMassSim";
  varnames[14] = "AreaRec";
  varnames[15] = "AreaSim";
  varnames[16] = "NEFRec";
  varnames[17] = "NEFSim";
  varnames[18] = "MassRec";
  varnames[19] = "MassSim";
  varnames[20] = "ZgMeasured";
  varnames[21] = "ZgTrue";
  varnames[22] = "RgMeasured";
  varnames[23] = "RgTrue";
  varnames[24] = "MgMeasured";
  varnames[25] = "MgTrue";
  varnames[26] = "PtgMeasured";
  varnames[27] = "PtgTrue";
  varnames[28] = "MugMeasured";
  varnames[29] = "MugTrue";
  varnames[30] = "OneSubjettinessMeasured";
  varnames[31] = "OneSubjettinessTrue";
  varnames[32] = "TwoSubjettinessMeasured";
  varnames[33] = "TwoSubjettinessTrue";
  varnames[34] = "AngularityMeasured";
  varnames[35] = "AngularityTrue";
  varnames[36] = "PtDMeasured";
  varnames[37] = "PtDTrue";
  varnames[38] = "NCharged";
  varnames[39] = "NNeutral";
  varnames[40] = "NConstTrue";
  varnames[41] = "NDroppedMeasured";
  varnames[42] = "NDroppedTrue";

  for(int ib = 0; ib < kTNVar; ib++){
    LinkOutputBranch(varnames[ib], fJetTreeData + ib);
  }
  PostData(1, fOutput);
  PostData(2, fJetSubstructureTree);
}

void AliAnalysisTaskEmcalJetSubstructureTree::LinkOutputBranch(const TString &branchname, Double_t *datalocation) {
  // Check whether branch is rejected
  if(!fFillPart && IsPartBranch(branchname)) return;
  if(!fFillAcceptance && IsAcceptanceBranch(branchname)) return;
  if(!fFillRho && IsRhoBranch(branchname)) return;
  if(!fFillMass && IsMassBranch(branchname)) return;
  if(!fFillSoftDrop && IsSoftdropBranch(branchname)) return;
  if(!fFillNSub && IsNSubjettinessBranch(branchname)) return;
  if(!fFillStructGlob && IsStructbranch(branchname)) return;

  std::cout << "Adding branch " << branchname << std::endl;
  fJetSubstructureTree->Branch(branchname, datalocation, Form("%s/D", branchname.Data()));  
}

void AliAnalysisTaskEmcalJetSubstructureTree::RunChanged(Int_t newrun) {
  if(fUseDownscaleWeight){
    AliEmcalDownscaleFactorsOCDB::Instance()->SetRun(newrun);
  }
}

bool AliAnalysisTaskEmcalJetSubstructureTree::Run(){
  AliClusterContainer *clusters = GetClusterContainer(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::ClusterContainerNameFactory(fInputEvent->IsA() == AliAODEvent::Class()));
  AliTrackContainer *tracks = GetTrackContainer(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::TrackContainerNameFactory(fInputEvent->IsA() == AliAODEvent::Class()));
  AliParticleContainer *particles = GetParticleContainer("mcparticles");

  AliJetContainer *mcjets = GetJetContainer("mcjets");
  AliJetContainer *datajets = GetJetContainer("datajets");

  // for(auto e : *(fInputEvent->GetList())) std::cout << e->GetName() << std::endl;

  TString rhoTagData = datajets ? TString::Format("R%02d", static_cast<Int_t>(datajets->GetJetRadius() * 10.)) : "",
          rhoTagMC = mcjets ? TString::Format("R%02d", static_cast<Int_t>(mcjets->GetJetRadius() * 10.)) : "";

  AliRhoParameter *rhoPtRec = GetRhoFromEvent("RhoSparse_Full_" + rhoTagData),
                  *rhoMassRec = GetRhoFromEvent("RhoMassSparse_Full_" + rhoTagData),
                  *rhoPtSim = GetRhoFromEvent("RhoSparse_Full_" + rhoTagMC),
                  *rhoMassSim = GetRhoFromEvent("RhoMassSparse_Full_" + rhoTagMC);
  AliDebugStream(2) << "Found rho parameter for reconstructed pt:    " << (rhoPtRec ? "yes" : "no") << ", value: " << (rhoPtRec ? rhoPtRec->GetVal() : 0.) << std::endl;
  AliDebugStream(2) << "Found rho parameter for sim pt:              " << (rhoPtSim ? "yes" : "no") << ", value: " << (rhoPtSim ? rhoPtSim->GetVal() : 0.) << std::endl;
  AliDebugStream(2) << "Found rho parameter for reconstructed Mass:  " << (rhoMassRec ? "yes" : "no") << ", value: " << (rhoMassRec ? rhoMassRec->GetVal() : 0.) << std::endl;
  AliDebugStream(2) << "Found rho parameter for sim Mass:            " << (rhoMassSim ? "yes" : "no") << ", value: " << (rhoMassSim ? rhoMassSim->GetVal() : 0.) << std::endl;

  AliDebugStream(1) << "Inspecting jet radius " << (datajets ? datajets->GetJetRadius() : mcjets->GetJetRadius()) << std::endl;

  double weight = 1.;
  if(fUseDownscaleWeight){
    weight = AliEmcalDownscaleFactorsOCDB::Instance()->GetDownscaleFactorForTriggerClass(this->fTriggerSelectionString);
  }

  // Run trigger selection (not on pure MCgen train)
  if(datajets){
    if(!(fInputHandler->IsEventSelected() & fTriggerSelectionBits)) return false;
    if(!mcjets){
      // Pure data - do EMCAL trigger selection from selection string
      if(fTriggerSelectionString.Length()) {
        if(!fInputEvent->GetFiredTriggerClasses().Contains(fTriggerSelectionString)) return false;
      }
    } else {
      // Simulation - do EMCAL trigger selection from trigger selection object
      PWG::EMCAL::AliEmcalTriggerDecisionContainer *mctrigger = static_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer *>(fInputEvent->FindListObject("EmcalTriggerDecision"));
      AliDebugStream(1) << "Found trigger decision object: " << (mctrigger ? "yes" : "no") << std::endl;
      if(fTriggerSelectionString.Length()){
        if(!mctrigger){
          AliErrorStream() <<  "Trigger decision container not found in event - not possible to select EMCAL triggers" << std::endl;
          return false;
        }
        if(!mctrigger->IsEventSelected(fTriggerSelectionString)) return false;
      }
    }
  }

  Double_t rhoparameters[4]; memset(rhoparameters, 0, sizeof(Double_t) * 4);
  if(rhoPtRec) rhoparameters[0] = rhoPtRec->GetVal();
  if(rhoPtSim) rhoparameters[1] = rhoPtSim->GetVal();
  if(rhoMassRec) rhoparameters[2] = rhoMassRec->GetVal();
  if(rhoMassSim) rhoparameters[3] = rhoMassSim->GetVal();

  AliSoftdropDefinition softdropSettings;
  softdropSettings.fBeta = fSDBetaCut;
  softdropSettings.fZ = fSDZCut;
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
        try {
          DoConstituentQA(jet, tracks, clusters);
          AliJetSubstructureData structureData =  MakeJetSubstructure(*jet, datajets->GetJetRadius() * 2., tracks, clusters, {softdropSettings, nsubjettinessSettings}),
                                 structureMC = fFillPart ? MakeJetSubstructure(*associatedJet, mcjets->GetJetRadius() * 2, particles, nullptr, {softdropSettings, nsubjettinessSettings}) : AliJetSubstructureData();
          Double_t angularity[2] = {fFillStructGlob ? MakeAngularity(*jet, tracks, clusters) : 0., (fFillStructGlob && fFillPart) ? MakeAngularity(*associatedJet, particles, nullptr) : 0.},
                   ptd[2] = {fFillStructGlob ? MakePtD(*jet, tracks, clusters) : 0., (fFillStructGlob && fFillPart) ? MakePtD(*associatedJet, particles, nullptr) : 0};
          FillTree(datajets->GetJetRadius(), weight, jet, associatedJet, &(structureData.fSoftDrop), &(structureMC.fSoftDrop), &(structureData.fNsubjettiness), &(structureMC.fNsubjettiness), angularity, ptd, rhoparameters);
        } catch(ReclusterizerException &e) {
          AliErrorStream() << "Error in reclusterization - skipping jet" << std::endl;
        } catch(SubstructureException &e) {
          AliErrorStream() << "Error in substructure observable - skipping jet" << std::endl;
        }
      } else {
        try {
          DoConstituentQA(jet, tracks, clusters);
          AliJetSubstructureData structure = MakeJetSubstructure(*jet, 0.4, tracks, clusters, {softdropSettings, nsubjettinessSettings});
          Double_t angularity[2] = {fFillStructGlob ? MakeAngularity(*jet, tracks, clusters): 0., 0.},
                   ptd[2] = {fFillStructGlob ? MakePtD(*jet, tracks, clusters) : 0., 0.};
          FillTree(datajets->GetJetRadius(), weight, jet, nullptr, &(structure.fSoftDrop), nullptr, &(structure.fNsubjettiness), nullptr, angularity, ptd, rhoparameters);
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
          Double_t angularity[2] = {0., MakeAngularity(*mcjet, particles, nullptr)},
                   ptd[2] = {0., MakePtD(*mcjet, particles, nullptr)};
          FillTree(mcjets->GetJetRadius(), weight, nullptr, mcjet, nullptr, &(structure.fSoftDrop), nullptr, &(structure.fNsubjettiness), angularity, ptd, rhoparameters);
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

void AliAnalysisTaskEmcalJetSubstructureTree::FillTree(double r, double weight,
                                                       const AliEmcalJet *datajet, const AliEmcalJet *mcjet,
                                                       AliSoftDropParameters *dataSoftdrop, AliSoftDropParameters *mcSoftdrop,
                                                       AliNSubjettinessParameters *dataSubjettiness, AliNSubjettinessParameters *mcSubjettiness,
                                                       Double_t *angularity, Double_t *ptd, Double_t *rhoparameters){

  memset(fJetTreeData,0, sizeof(Double_t) * kTNVar);
  fJetTreeData[kTRadius] = r;
  fJetTreeData[kTWeight] = weight;
  if(fFillRho){
    fJetTreeData[kTRhoPtRec] = rhoparameters[0]; 
    if(fFillMass) fJetTreeData[kTRhoMassRec] = rhoparameters[2]; 
    if(fFillPart){
      fJetTreeData[kTRhoPtSim] = rhoparameters[1]; 
      if(fFillMass) fJetTreeData[kTRhoMassSim] = rhoparameters[3]; 
    }
  }
  if(datajet) {
    fJetTreeData[kTPtJetRec] = TMath::Abs(datajet->Pt());
    fJetTreeData[kTNCharged] = datajet->GetNumberOfTracks();
    fJetTreeData[kTNNeutral] = datajet->GetNumberOfClusters();
    fJetTreeData[kTAreaRec] = datajet->Area();
    fJetTreeData[kTNEFRec] = datajet->NEF();
    if(fFillMass) fJetTreeData[kTMassRec] = datajet->M();
    fJetTreeData[kTEJetRec] = datajet->E();
    if(fFillAcceptance) {
      fJetTreeData[kTEtaRec] = datajet->Eta();
      fJetTreeData[kTPhiRec] = datajet->Phi();
    }
  }

  if(fFillPart && mcjet){
    fJetTreeData[kTPtJetSim] = TMath::Abs(mcjet->Pt());
    fJetTreeData[kTNConstTrue] = mcjet->GetNumberOfConstituents();
    fJetTreeData[kTAreaSim] = mcjet->Area();
    fJetTreeData[kTNEFSim] = mcjet->NEF();
    if(fFillMass) fJetTreeData[kTMassSim] = mcjet->M();
    fJetTreeData[kTEJetSim] = mcjet->E();
    if(fFillAcceptance){
      fJetTreeData[kTEtaSim] = mcjet->Eta();
      fJetTreeData[kTPhiSim] = mcjet->Phi();
    }
  }

  if(fFillSoftDrop){
    if(dataSoftdrop) {
      fJetTreeData[kTZgMeasured] = dataSoftdrop->fZg;
      fJetTreeData[kTRgMeasured] = dataSoftdrop->fRg;
      fJetTreeData[kTMgMeasured] = dataSoftdrop->fMg;
      fJetTreeData[kTPtgMeasured] = dataSoftdrop->fPtg;
      fJetTreeData[kTMugMeasured] = dataSoftdrop->fMug;
      fJetTreeData[kTNDroppedMeasured] = dataSoftdrop->fNDropped;
    }

    if(fFillPart && mcSoftdrop) {
      fJetTreeData[kTZgTrue] = mcSoftdrop->fZg;
      fJetTreeData[kTRgTrue] = mcSoftdrop->fRg;
      fJetTreeData[kTMgTrue] = mcSoftdrop->fMg;
      fJetTreeData[kTPtgTrue] = mcSoftdrop->fPtg;
      fJetTreeData[kTMugTrue] = mcSoftdrop->fMug;
      fJetTreeData[kTNDroppedTrue] = mcSoftdrop->fNDropped;
    }
  }

  if(fFillNSub){
    if(dataSubjettiness) {
      fJetTreeData[kTOneNSubjettinessMeasured] = dataSubjettiness->fOneSubjettiness;
      fJetTreeData[kTTwoNSubjettinessMeasured] = dataSubjettiness->fTwoSubjettiness;
    }

    if(fFillPart && mcSubjettiness) {
      fJetTreeData[kTOneNSubjettinessTrue] = mcSubjettiness->fOneSubjettiness;
      fJetTreeData[kTTwoNSubjettinessTrue] = mcSubjettiness->fTwoSubjettiness;
    } 
  }

  if(fFillStructGlob){
    fJetTreeData[kTAngularityMeasured] = angularity[0];
    fJetTreeData[kTPtDMeasured] = ptd[0];

    if(fFillPart){
      fJetTreeData[kTAngularityTrue] = angularity[1];
      fJetTreeData[kTPtDTrue] = ptd[1];
    }
  }

  fJetSubstructureTree->Fill();
}


AliJetSubstructureData AliAnalysisTaskEmcalJetSubstructureTree::MakeJetSubstructure(const AliEmcalJet &jet, double jetradius, const AliParticleContainer *tracks, const AliClusterContainer *clusters, const AliJetSubstructureSettings &settings) const {
  const int kClusterOffset = 30000; // In order to handle tracks and clusters in the same index space the cluster index needs and offset, large enough so that there is no overlap with track indices
  std::vector<fastjet::PseudoJet> constituents;
  bool isMC = dynamic_cast<const AliTrackContainer *>(tracks);
  AliDebugStream(2) << "Make new jet substrucutre for " << (isMC ? "MC" : "data") << " jet: Number of tracks " << jet.GetNumberOfTracks() << ", clusters " << jet.GetNumberOfClusters() << std::endl;
  for(int itrk = 0; itrk < jet.GetNumberOfTracks(); itrk++){
    auto track = jet.TrackAt(itrk, tracks->GetArray());
    fastjet::PseudoJet constituentTrack(track->Px(), track->Py(), track->Pz(), track->E());
    constituentTrack.set_user_index(jet.TrackAt(itrk));
    constituents.push_back(constituentTrack);
  }

  if(clusters){
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
  if(!constituents.size())
    throw ReclusterizerException();
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
  fastjet::contrib::SoftDrop softdropAlgorithm(cutparameters.fBeta, cutparameters.fZ);
  softdropAlgorithm.set_verbose_structure(kTRUE);
  std::unique_ptr<fastjet::contrib::Recluster> reclusterizer(new fastjet::contrib::Recluster(cutparameters.fRecluserAlgo, 1, true));
  softdropAlgorithm.set_reclustering(kTRUE, reclusterizer.get());
  AliDebugStream(4) << "Jet has " << jet.constituents().size() << " constituents" << std::endl;
  auto groomed = softdropAlgorithm(jet);
  try {
    auto softdropstruct = groomed.structure_of<fastjet::contrib::SoftDrop>();

    AliSoftDropParameters result({softdropstruct.symmetry(),
                                  groomed.m(),
                                  softdropstruct.delta_R(),
                                  groomed.perp(),
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
  if(tracks){
    for(UInt_t itrk = 0; itrk < jet.GetNumberOfTracks(); itrk++) {
      auto track = jet.TrackAt(itrk, tracks->GetArray());
      if(!track){
        AliErrorStream() << "Associated constituent particle / track not found\n";
        continue;
      }
      TVector3 trackvec(track->Px(), track->Py(), track->Pz());

      num +=  track->Pt() * trackvec.DrEtaPhi(jetvec);
      den += +track->Pt();
    }
  }
  if(clusters) {
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
  if(particles){
    for(UInt_t itrk = 0; itrk < jet.GetNumberOfTracks(); itrk++) {
      auto trk = jet.TrackAt(itrk, particles->GetArray());
      if(!trk){
        AliErrorStream() << "Associated constituent particle / track not found\n";
        continue;
      }
      num += trk->Pt() * trk->Pt();
      den += trk->Pt();
    }
  }
  if(clusters){
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
    fQAHistos->FillTH2("hClusterConstTime", jet->Pt(), clust->GetTOF());
    fQAHistos->FillTH2("hClusterConstM02", jet->Pt(), clust->GetM02());
    fQAHistos->FillTH2("hClusterConstNcell", jet->Pt(), clust->GetNCells());

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
  AliDebugStream(2) << "Jet: Number of particle constituents: " << jet->GetParticleConstituents().GetEntriesFast() << std::endl;
  for(auto pconst : jet->GetParticleConstituents()) {
    auto part = static_cast<PWG::JETFW::AliEmcalParticleJetConstituent *>(pconst);
    AliDebugStream(3) << "Found particle constituent with pt " << part->Pt() << ", from VParticle " << part->GetParticle()->Pt() << std::endl;
    fQAHistos->FillTH2("hChargedConstituentPt", jet->Pt(), part->Pt());
  }

  // Loop over neutral constituents
  AliDebugStream(2) << "Jet: Number of cluster constituents: " << jet->GetClusterConstituents().GetEntriesFast() << std::endl;
  for(auto cconst : jet->GetClusterConstituents()){
    auto clust = static_cast<PWG::JETFW::AliEmcalClusterJetConstituent *>(cconst);
    AliDebugStream(3) << "Found cluster constituent with energy " << clust->E() << " using energy definition " << static_cast<int>(clust->GetDefaultEnergyType()) << std::endl;
    fQAHistos->FillTH2("hClusterConstituentEDefault", jet->Pt(), clust->E());
    fQAHistos->FillTH2("hClusterConstituentENLC", jet->Pt(), clust->GetCluster()->GetNonLinCorrEnergy());
    fQAHistos->FillTH2("hClusterConstituentEHC", jet->Pt(), clust->GetCluster()->GetHadCorrEnergy());
  }
#endif
}

bool AliAnalysisTaskEmcalJetSubstructureTree::IsPartBranch(const TString &branchname) const{
  return branchname.Contains("Sim") || branchname.Contains("True");
}

bool AliAnalysisTaskEmcalJetSubstructureTree::IsAcceptanceBranch(const TString &branchname) const {
  return branchname.Contains("Eta") || branchname.Contains("Phi");
}

bool AliAnalysisTaskEmcalJetSubstructureTree::IsRhoBranch(const TString &branchname) const{
  return branchname.Contains("Rho");
}

bool AliAnalysisTaskEmcalJetSubstructureTree::IsMassBranch(const TString &branchname) const{
  return branchname.Contains("Mass");     // also disable rho mass branch
}

bool AliAnalysisTaskEmcalJetSubstructureTree::IsSoftdropBranch(const TString &branchname) const{
  return branchname.Contains("gMeasured") || branchname.Contains("gTrue") || branchname.Contains("NDropped");
}

bool AliAnalysisTaskEmcalJetSubstructureTree::IsNSubjettinessBranch(const TString &branchname) const{
  return branchname.Contains("Subjettiness");
}

bool AliAnalysisTaskEmcalJetSubstructureTree::IsStructbranch(const TString &branchname) const{
  return branchname.Contains("Angularity") || branchname.Contains("PtD");
}

AliAnalysisTaskEmcalJetSubstructureTree *AliAnalysisTaskEmcalJetSubstructureTree::AddEmcalJetSubstructureTreeMaker(Bool_t isMC, Bool_t isData, Double_t jetradius, AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recombinationScheme, const char *trigger){
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

  AliAnalysisTaskEmcalJetSubstructureTree *treemaker = new AliAnalysisTaskEmcalJetSubstructureTree("JetSubstructureTreemaker_" + TString::Format("R%02d_", int(jetradius * 10.)) + trigger);
  mgr->AddTask(treemaker);
  treemaker->SetMakeGeneralHistograms(kTRUE);

  // Adding containers
  if(isMC) {
    AliParticleContainer *particles = treemaker->AddMCParticleContainer("mcparticles");
    particles->SetMinPt(0.);

    AliJetContainer *mcjets = treemaker->AddJetContainer(
                              jettype,
                              AliJetContainer::antikt_algorithm,
                              recombinationScheme,
                              jetradius,
                              (isData && ((jettype == AliJetContainer::kFullJet) || (jettype == AliJetContainer::kNeutralJet))) ? AliEmcalJet::kEMCALfid : AliEmcalJet::kTPC,
                              particles, nullptr);
    mcjets->SetName("mcjets");
    mcjets->SetJetPtCut(20.);
  }

  if(isData) {
    AliTrackContainer *tracks(nullptr);
    if((jettype == AliJetContainer::kChargedJet) || (jettype == AliJetContainer::kFullJet)){
      tracks = treemaker->AddTrackContainer(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::TrackContainerNameFactory(isAOD));
      std::cout << "Track container name: " << tracks->GetName() << std::endl;
      tracks->SetMinPt(0.15);
    }
    AliClusterContainer *clusters(nullptr);
    if((jettype == AliJetContainer::kFullJet) || (AliJetContainer::kNeutralJet)){
      std::cout << "Using full or neutral jets ..." << std::endl;
      clusters = treemaker->AddClusterContainer(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::ClusterContainerNameFactory(isAOD));
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
                              ((jettype == AliJetContainer::kFullJet) || (jettype == AliJetContainer::kNeutralJet)) ? AliEmcalJet::kEMCALfid : AliEmcalJet::kTPCfid,
                              tracks, clusters);
    datajets->SetName("datajets");
    datajets->SetJetPtCut(20.);

    treemaker->SetUseAliAnaUtils(true, true);
    treemaker->SetVzRange(-10., 10);

    // configure trigger selection
    TString triggerstring(trigger);
    if(triggerstring.Contains("INT7")) {
      treemaker->SetTriggerBits(AliVEvent::kINT7);
    } else if(triggerstring.Contains("EJ1")) {
      treemaker->SetTriggerBits(AliVEvent::kEMCEJE);
      treemaker->SetTriggerString("EJ1");
    } else if(triggerstring.Contains("EJ2")) {
      treemaker->SetTriggerBits(AliVEvent::kEMCEJE);
      treemaker->SetTriggerString("EJ2");
    }
  }

  // Connecting containers
  TString outputfile = mgr->GetCommonFileName();
  outputfile += TString::Format(":JetSubstructure_R%02d_%s", int(jetradius * 10.), trigger);
  mgr->ConnectInput(treemaker, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(treemaker, 1, mgr->CreateContainer("JetSubstructureHistos_" + TString::Format("R%02d_", int(jetradius * 10.)) + trigger, AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, outputfile));
  mgr->ConnectOutput(treemaker, 2, mgr->CreateContainer("JetSubstructureTree_" + TString::Format("R%02d_", int(jetradius * 10.)) + trigger, TTree::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName()));

  return treemaker;
}

} /* namespace EmcalTriggerJets */
