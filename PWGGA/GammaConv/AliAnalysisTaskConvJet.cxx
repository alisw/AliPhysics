/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

#include <cstddef>
#include <cstring>

#include <TObjString.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>

#include <AliAnalysisManager.h>
#include <AliVEventHandler.h>
#include <AliVCluster.h>
#include <AliVParticle.h>
#include <AliLog.h>

#include "AliTLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

#include "AliMCEvent.h"

#include "AliAnalysisTaskConvJet.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskConvJet);
/// \endcond

AliAnalysisTaskConvJet::AliAnalysisTaskConvJet() : AliAnalysisTaskEmcalJet(),
                                                   fNJets(0),
                                                   fVectorJetPt(0),
                                                   fVectorJetPx(0),
                                                   fVectorJetPy(0),
                                                   fVectorJetPz(0),
                                                   fVectorJetEta(0),
                                                   fVectorJetPhi(0),
                                                   fVectorJetR(0),
                                                   fVectorJetNEF({}),
                                                   fVectorJetNClus({}),
                                                   fVectorJetNCh({}),
                                                   fTrueNJets(0),
                                                   fTrueVectorJetPt(0),
                                                   fTrueVectorJetPx(0),
                                                   fTrueVectorJetPy(0),
                                                   fTrueVectorJetPz(0),
                                                   fTrueVectorJetEta(0),
                                                   fTrueVectorJetPhi(0),
                                                   fTrueVectorJetR(0),
                                                   fTrueVectorJetNPart(0),
                                                   fTrueVectorJetParton(0),
                                                   fTrueVectorJetPartonPt(0),
                                                   fTrueVectorJetPartonPx(0),
                                                   fTrueVectorJetPartonPy(0),
                                                   fTrueVectorJetPartonPz(0),
                                                   fVecJetClusters({}),
                                                   fVecJetTracks({}),
                                                   fVecTrueJetParticles({}),
                                                   fAccType(0),
                                                   fAccTypeMC(0),
                                                   fDistToEMCBorder(0),
                                                   fEMCSMEdgesMode(0),
                                                   fDistEMCSMEdge(0),
                                                   fApplyEnergyWeight(false),
                                                   funcEnergyWeights(nullptr),
                                                   fVecMeasurable({22, 211, 321, 2212, 11, 13}),
                                                   funcJES(nullptr),
                                                   funcJER(nullptr),
                                                   funcJERCut(nullptr)
{
}

AliAnalysisTaskConvJet::AliAnalysisTaskConvJet(const char* name) : AliAnalysisTaskEmcalJet(name, kTRUE),
                                                                   fNJets(0),
                                                                   fVectorJetPt(0),
                                                                   fVectorJetPx(0),
                                                                   fVectorJetPy(0),
                                                                   fVectorJetPz(0),
                                                                   fVectorJetEta(0),
                                                                   fVectorJetPhi(0),
                                                                   fVectorJetR(0),
                                                                   fVectorJetNEF({}),
                                                                   fVectorJetNClus({}),
                                                                   fVectorJetNCh({}),
                                                                   fTrueVectorJetPt(0),
                                                                   fTrueVectorJetPx(0),
                                                                   fTrueVectorJetPy(0),
                                                                   fTrueVectorJetPz(0),
                                                                   fTrueVectorJetEta(0),
                                                                   fTrueVectorJetPhi(0),
                                                                   fTrueVectorJetR(0),
                                                                   fTrueVectorJetNPart(0),
                                                                   fTrueVectorJetParton(0),
                                                                   fTrueVectorJetPartonPt(0),
                                                                   fTrueVectorJetPartonPx(0),
                                                                   fTrueVectorJetPartonPy(0),
                                                                   fTrueVectorJetPartonPz(0),
                                                                   fVecJetClusters({}),
                                                                   fVecJetTracks({}),
                                                                   fVecTrueJetParticles({}),
                                                                   fAccType(0),
                                                                   fAccTypeMC(0),
                                                                   fDistToEMCBorder(0),
                                                                   fEMCSMEdgesMode(0),
                                                                   fDistEMCSMEdge(0),
                                                                   fApplyEnergyWeight(false),
                                                                   funcEnergyWeights(nullptr),
                                                                   fVecMeasurable({22, 211, 321, 2212, 11, 13}),
                                                                   funcJES(nullptr),
                                                                   funcJER(nullptr),
                                                                   funcJERCut(nullptr)
{
  SetMakeGeneralHistograms(kTRUE);
}

AliAnalysisTaskConvJet::~AliAnalysisTaskConvJet()
{
}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskConvJet::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

/**
 * The body of this function should contain instructions to fill the vectors.
 * This function is called inside the event loop, after the function Run() has been
 * executed successfully (i.e. it returned kTRUE).
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskConvJet::FillHistograms()
{
  DoJetLoop();

  return kTRUE;
}

/**
 * This function performs a loop over the reconstructed jets
 * in the current event and fills the vectors.
 */
void AliAnalysisTaskConvJet::DoJetLoop()
{
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    TString JetName = jetCont->GetTitle();
    TObjArray* arr = JetName.Tokenize("__");
    TObjString* testObjString = (TObjString*)arr->At(2);
    if (!(static_cast<TString>(testObjString->GetString())).Contains("mcparticles") || (static_cast<TString>(testObjString->GetString()).Contains("mcparticles") && static_cast<TString>(testObjString->GetString()).Contains("Rec"))) {
      fAccType = jetCont->GetAcceptanceType();
      UInt_t count = 0;
      fNJets = 0;
      fVectorJetPt.clear();
      fVectorJetPx.clear();
      fVectorJetPy.clear();
      fVectorJetPz.clear();
      fVectorJetEta.clear();
      fVectorJetPhi.clear();
      fVectorJetR.clear();
      fVectorJetNEF.clear();
      fVectorJetNClus.clear();
      fVectorJetNCh.clear();
      fVecJetClusters.clear();
      fVecJetTracks.clear();

      for (auto const& jet : jetCont->accepted()) {
        if (!jet)
          continue;
        if (!IsJetAccepted(jet))
          continue;

        std::vector<AliVCluster*> vecTmpClus;
        for(size_t cl = 0; cl < jet->GetNumberOfClusters(); ++cl){
          vecTmpClus.push_back(jet->Cluster(cl));
        }
        fVecJetClusters.push_back(vecTmpClus);

        std::vector<AliVParticle*> vecTmpTracks;
        for(size_t tr = 0; tr < jet->GetNumberOfTracks(); ++tr){
          vecTmpTracks.push_back(jet->Track(tr));
        }
        fVecJetTracks.push_back(vecTmpTracks);

        double jetEnergyWeight = 1.;
        if(fApplyEnergyWeight == 1){
          jetEnergyWeight = funcEnergyWeights->Eval(jet->Pt());
        } 
        count++;
        fVectorJetPt.push_back(jet->Pt()*jetEnergyWeight);
        fVectorJetPx.push_back(jet->Px()*jetEnergyWeight);
        fVectorJetPy.push_back(jet->Py()*jetEnergyWeight);
        fVectorJetPz.push_back(jet->Pz()*jetEnergyWeight);
        fVectorJetEta.push_back(jet->Eta());
        fVectorJetPhi.push_back(jet->Phi());
        fVectorJetR.push_back(jet->Area());
        fVectorJetNEF.push_back(jet->NEF());
        fVectorJetNClus.push_back(jet->Nn());
        fVectorJetNCh.push_back(jet->Nch());
      }
      fNJets = count;
    } else {
      fAccTypeMC = jetCont->GetAcceptanceType();
      UInt_t count = 0;
      fTrueNJets = 0;
      fTrueVectorJetPt.clear();
      fTrueVectorJetPx.clear();
      fTrueVectorJetPy.clear();
      fTrueVectorJetPz.clear();
      fTrueVectorJetEta.clear();
      fTrueVectorJetPhi.clear();
      fTrueVectorJetR.clear();
      fTrueVectorJetNPart.clear();
      fVecTrueJetMaxPartPt.clear();
      fVecTrueJetMaxPartPDG.clear();
      for(auto & vec: fVecTrueJetParticles){
        vec.clear();
      }
      fVecTrueJetParticles.clear();
      for (auto const& jet : jetCont->accepted()) {
        if (!jet)
          continue;
        count++;
        fTrueVectorJetPt.push_back(jet->Pt());
        fTrueVectorJetPx.push_back(jet->Px());
        fTrueVectorJetPy.push_back(jet->Py());
        fTrueVectorJetPz.push_back(jet->Pz());
        fTrueVectorJetEta.push_back(jet->Eta());
        fTrueVectorJetPhi.push_back(jet->Phi());
        fTrueVectorJetR.push_back(jet->Area());
        fTrueVectorJetNPart.push_back(jet->N());

        std::vector<AliVParticle*> vecTmpPart;
        for(size_t tr = 0; tr < jet->GetNumberOfTracks(); ++tr){
          vecTmpPart.push_back(jet->Track(tr));
        }
        fVecTrueJetParticles.push_back(vecTmpPart);
        auto [ptLead, pdgcodeLead] = GetLeadingPartPt(jet, true);
        fVecTrueJetMaxPartPt.push_back(ptLead);
        fVecTrueJetMaxPartPDG.push_back(pdgcodeLead);
      }
      fTrueNJets = count;
    }
  }
  if(fApplyEnergyWeight == 2){
    // match true and rec jets and assign weights to rec. jet accoring to particle composition in true jet
    for(size_t iRec = 0; iRec < fVectorJetPt.size(); ++iRec){
      double etaRec = fVectorJetEta[iRec];
      double phiRec = fVectorJetPhi[iRec];
      double Rmatch = 100.;
      size_t Imatch = 0;
      for(size_t iTrue = 0; iTrue < fTrueVectorJetPt.size(); ++iTrue){
        double etaTrue = fTrueVectorJetEta[iTrue];
        double phiTrue = fTrueVectorJetPhi[iTrue];
        double deltaEta = etaRec - etaTrue;
        double deltaPhi = phiRec - phiTrue;
        if (deltaPhi > M_PI) {
          deltaPhi = 2 * M_PI - deltaPhi;
        }
        double R_jetjet = sqrt(pow((deltaEta), 2) + pow((deltaPhi), 2));
        if (R_jetjet < Rmatch) {
          Rmatch = R_jetjet;
          Imatch = iTrue;
        }
      }
      if(Rmatch < 0.3){ // hardcoded min distance requirement
        double ptNonMeas = 0;
        for(const auto & tr : fVecTrueJetParticles[Imatch]){
          if(!tr) {
            AliError("track from vecTmpPart is null");
            continue;
          }
          int mcLabel = std::abs(tr->GetLabel());
          auto mcParticle = dynamic_cast<AliAODMCParticle*>(fMCEvent->GetTrack(mcLabel));
          if(!mcParticle) {
            AliError("mcParticle is null");
            continue;
          }
          if(IsNonMeasurable(std::abs(mcParticle->GetPdgCode()), tr->Charge())){
            ptNonMeas+=tr->Pt();
          }
        }
        double jetEnergyWeight = funcEnergyWeights->Eval(fTrueVectorJetPt[Imatch]) - 1.;
        if(fTrueVectorJetPt[Imatch]>0){
          jetEnergyWeight*= ptNonMeas/fTrueVectorJetPt[Imatch];
        }
        jetEnergyWeight += 1.;
        // now scale it
        fVectorJetPt[iRec]*=jetEnergyWeight;
        fVectorJetPx[iRec]*=jetEnergyWeight;
        fVectorJetPy[iRec]*=jetEnergyWeight;
        fVectorJetPz[iRec]*=jetEnergyWeight;
      }
    }
  }
}

/**
 * This function finds the leading /hardest parton in each jet
 * The id of the mc particles is written to fTrueVectorJetParton
 */
void AliAnalysisTaskConvJet::FindPartonsJet(TClonesArray* arrMCPart)
{
  // Loop over all primary MC particle
  fTrueVectorJetParton.resize(fTrueNJets);
  fTrueVectorJetPartonPt.resize(fTrueNJets);
  fTrueVectorJetPartonPx.resize(fTrueNJets);
  fTrueVectorJetPartonPy.resize(fTrueNJets);
  fTrueVectorJetPartonPz.resize(fTrueNJets);
  std::vector<double> partonEnergy(fTrueNJets, -1);
  double JetR2 = Get_Jet_Radius() * Get_Jet_Radius();

  for (Long_t i = 0; i < arrMCPart->GetEntriesFast(); i++) {
    AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(arrMCPart->At(i));
    if (!particle)
      continue;

    // particle has to be quark or gluon
    if (std::abs(particle->GetPdgCode()) == 21 || (std::abs(particle->GetPdgCode()) < 9 && std::abs(particle->GetPdgCode()) > 0)) {

      int indexNearestJet = -1;
      double deltaR2 = 100000;
      // loop over all jets and find the closest jet to the parton inside the jet radius
      for (int j = 0; j < fTrueNJets; ++j) {
        // check if mc particle and jet coincide
        double dEta = particle->Eta() - fTrueVectorJetEta[j];
        double dPhi = particle->Phi() - fTrueVectorJetPhi[j];
        double deltaR2tmp = dEta * dEta + dPhi * dPhi;
        if (deltaR2tmp < JetR2 && deltaR2tmp < deltaR2) {
          indexNearestJet = j;
          deltaR2 = deltaR2tmp;
        }
      }
      if (indexNearestJet == -1) {
        continue;
      }

      // assign a parton to the jet if the parton has an energy thats larger than the previously assigned value
      if (partonEnergy[indexNearestJet] < particle->Pt()) {
        partonEnergy[indexNearestJet] = particle->Pt();
        fTrueVectorJetParton[indexNearestJet] = i;
        fTrueVectorJetPartonPt[indexNearestJet] = particle->Pt();
        fTrueVectorJetPartonPx[indexNearestJet] = particle->Px();
        fTrueVectorJetPartonPy[indexNearestJet] = particle->Py();
        fTrueVectorJetPartonPz[indexNearestJet] = particle->Pz();
      }
    }
  }
}

/**
 * This function checks if the jet is inside the desired acceptance
 * if a distance to the EMCal border is required
 * For example, a distance of 0.4 to the border ensures, that every 0.4 jet is fully on the EMCal
 */
bool AliAnalysisTaskConvJet::IsJetAccepted(const AliEmcalJet* jet)
{
  double jetPhi = jet->Phi();
  if(jetPhi < 0) jetPhi += 2*TMath::Pi();
  bool accept = true;
  if(fDistToEMCBorder > 0){
    // geometry values from https://arxiv.org/pdf/2209.04216.pdf (page 10)
    if (jetPhi < 1.40 + fDistToEMCBorder ||
        jetPhi > 5.70 - fDistToEMCBorder ||
        (jetPhi > 3.26 - fDistToEMCBorder && jetPhi < 4.54 + fDistToEMCBorder)) {
      accept = false;
    }
    if (std::abs(jet->Eta()) > 0.7 - fDistToEMCBorder ||
        (std::abs(jet->Eta()) < 0.23 + fDistToEMCBorder && jetPhi > 4.54 && jetPhi < 5.58)) {
      accept = false;
    }
  }
  if(fEMCSMEdgesMode > 0){
    // first, check in steps of 20 degrees
    const double angle = TMath::Pi()/9.; // twenty degree
    // const double angleDiff = TMath::Pi()/60.; // 3 degrees
    // Calculate the remainder when dividing by 20
    double remainder = fmod(jetPhi, angle);

    // Check if it's within the range ±angleDiff
    bool isAtBorder = false;
    if (remainder <= fDistEMCSMEdge || (angle - remainder) <= fDistEMCSMEdge) {
      isAtBorder = true;
    }
    if(fEMCSMEdgesMode == 1){
      accept = !isAtBorder;
    } else if(fEMCSMEdgesMode == 2){
      accept = isAtBorder;
    }
  }
  return accept;
}

/**
 * This function is used to get the leading particle pt and MC
 */
std::tuple<double, int> AliAnalysisTaskConvJet::GetLeadingPartPt(AliEmcalJet * jet, const bool isTrueJet){
  double maxTrackPt = 0;
  double maxClusterPt = 0;

  auto particle = jet->GetLeadingParticleConstituent();
  if (particle) {
    maxTrackPt = particle->Pt();
  }
  auto cluster = jet->GetLeadingClusterConstituent();
  if (cluster) {
    // Uses the energy definition that was used when the constituent was created
    // to calculate the Pt(). Usually, this would be the hadronic corrected energy
    maxClusterPt = cluster->Pt();
  }
  int pdgcode = 0;
  if(isTrueJet){
    auto mcpart =particle->GetParticle();
    pdgcode = mcpart->PdgCode();
  }
  
  double maxPt = (maxTrackPt > maxClusterPt) ? maxTrackPt : maxClusterPt;
  return std::make_tuple(maxPt, pdgcode);  
}

/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskConvJet::ExecOnce()
{
  AliAnalysisTaskEmcalJet::ExecOnce();
}

/**
 * Run analysis code here, if needed.
 * It will be executed before FillHistograms().
 * If this function return kFALSE, FillHistograms() will *not*
 * be executed for the current event
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskConvJet::Run()
{
  return kTRUE;
}

/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskConvJet::Terminate(Option_t*)
{
}

/**
 * This function adds the task to the analysis manager. Often, this function is called
 * by an AddTask C macro. However, by compiling the code, it ensures that we do not
 * have to deal with difficulties caused by CINT.
 */
AliAnalysisTaskConvJet* AliAnalysisTaskConvJet::AddTask_GammaConvJet(
  const char* ntracks,
  const char* nclusters,
  const char* ncells,
  const char* suffix,
  const double distToEMCBorder,
  const double distToSMEdges
  )
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTask_GammaConvJet", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler) {
    ::Error("AddTask_GammaConvJet", "This task requires an input event handler");
    return 0;
  }

  enum EDataType_t {
    kUnknown,
    kESD,
    kAOD
  };

  EDataType_t dataType = kUnknown;

  if (handler->InheritsFrom("AliESDInputHandler")) {
    dataType = kESD;
  } else if (handler->InheritsFrom("AliAODInputHandler")) {
    dataType = kAOD;
  }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString trackName(ntracks);
  TString clusName(nclusters);
  TString cellName(ncells);

  if (trackName == "usedefault") {
    if (dataType == kESD) {
      trackName = "Tracks";
    } else if (dataType == kAOD) {
      trackName = "tracks";
    } else {
      trackName = "";
    }
  }

  if (clusName == "usedefault") {
    if (dataType == kESD) {
      clusName = "CaloClusters";
    } else if (dataType == kAOD) {
      clusName = "caloClusters";
    } else {
      clusName = "";
    }
  }

  if (cellName == "usedefault") {
    if (dataType == kESD) {
      cellName = "EMCALCells";
    } else if (dataType == kAOD) {
      cellName = "emcalCells";
    } else {
      cellName = "";
    }
  }

  TString name(Form("AliAnalysisTaskConvJet%s", strlen(suffix) == 0 ? "" : Form("_%s", suffix)));

  AliAnalysisTaskConvJet* sampleTask = new AliAnalysisTaskConvJet(name);
  sampleTask->SetCaloCellsName(cellName);
  sampleTask->SetVzRange(-10, 10);
  sampleTask->SetDistToEMCBorder(distToEMCBorder);
  if(distToSMEdges != 0) sampleTask->SetDistToEMCSMEdge(std::abs(distToSMEdges), 1 + (distToSMEdges < 0));

  if (trackName == "mcparticles") {
    sampleTask->AddMCParticleContainer(trackName);
  } else if (trackName == "tracks" || trackName == "Tracks") {
    sampleTask->AddTrackContainer(trackName);
  } else if (!trackName.IsNull()) {
    sampleTask->AddParticleContainer(trackName);
  }
  sampleTask->AddClusterContainer(clusName);

  sampleTask->GetClusterContainer(0)->SetClusECut(0.);
  sampleTask->GetClusterContainer(0)->SetClusPtCut(0.);
  sampleTask->GetClusterContainer(0)->SetClusNonLinCorrEnergyCut(0.);
  sampleTask->GetClusterContainer(0)->SetClusHadCorrEnergyCut(0.30);
  sampleTask->GetClusterContainer(0)->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  sampleTask->GetParticleContainer(0)->SetParticlePtCut(0.15);
  sampleTask->GetParticleContainer(0)->SetParticleEtaLimits(-0.8, 0.8);

  if (trackName != "mcparticles") {
    sampleTask->GetTrackContainer(0)->SetFilterHybridTracks(kTRUE);
    sampleTask->GetTrackContainer(0)->SetParticlePtCut(0.15);
    sampleTask->GetTrackContainer(0)->SetParticleEtaLimits(-0.8, 0.8);
  }

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  mgr->AddTask(sampleTask);

  // Create containers for input/output
  AliAnalysisDataContainer* cinput1 = mgr->GetCommonInputContainer();
  TString contname(trackName);
  contname += "_histos";
  contname += strlen(suffix) == 0 ? "" : Form("_%s", suffix);
  printf("contname: %s", contname.Data());
  AliAnalysisDataContainer* coutput1 = mgr->CreateContainer(contname.Data(),
                                                            TList::Class(), AliAnalysisManager::kOutputContainer,
                                                            Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput(sampleTask, 0, cinput1);
  mgr->ConnectOutput(sampleTask, 1, coutput1);

  return sampleTask;
}


void AliAnalysisTaskConvJet::setWeightEnergyJets(const char * formula, const int mode){
  printf("setWeightEnergyJets %i\n", mode);
  fApplyEnergyWeight = mode;
  funcEnergyWeights = new TF1("funcEnergyWeightsJets", formula, 0, 10000);
}

void AliAnalysisTaskConvJet::SetMeasurablePart(TString str){
  fVecMeasurable.clear();  // Ensure the vector is empty

  // Split the input by semicolons
  TObjArray* tokens = str.Tokenize(";");

  if (tokens) {
      for (int i = 0; i < tokens->GetEntries(); ++i) {
          TObjString* objStr = (TObjString*)tokens->At(i);
          if (objStr) {
              TString token = objStr->GetString();
              if (!token.IsWhitespace()) {  // Skip empty tokens
                  fVecMeasurable.push_back(token.Atoi());  // Convert to int and add to vector
              }
          }
      }
      delete tokens; // Clean up to avoid memory leaks
    }
}

bool AliAnalysisTaskConvJet::IsNonMeasurable(const int pdg, const int charge){
  if(std::find(fVecMeasurable.begin(), fVecMeasurable.end(), pdg) != fVecMeasurable.end()) {
    return false;
  }
  return true;
}
