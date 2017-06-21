/*************************************************************************
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

#include "AliAnalysisManager.h"
#include "AliVEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliHFTrackContainer.h"
#include "AliHFAODMCParticleContainer.h"

#include "AliAnalysisTaskDmesonJetsDetectorResponse.h"

// Definitions of class AliAnalysisTaskDmesonJetsDetectorResponse::AliDmesonMatchInfoSummary

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJetsDetectorResponse::AliDmesonMatchInfoSummary);
/// \endcond

// Definitions of class AliAnalysisTaskDmesonJetsDetectorResponse::AliD0MatchInfoSummary

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJetsDetectorResponse::AliD0MatchInfoSummary);
/// \endcond

/// Reset the object
void AliAnalysisTaskDmesonJetsDetectorResponse::AliD0MatchInfoSummary::Reset()
{
  fGenerated.Reset();
  fReconstructed.Reset();
}

/// Set the current object using instances of AliDmesonJetInfo as its source
///
/// \param reco A const reference to a AliDmesonJetInfo object representing the reconstructed D meson
void AliAnalysisTaskDmesonJetsDetectorResponse::AliD0MatchInfoSummary::SetReconstructed(const AliDmesonJetInfo& reco)
{
  fReconstructed.Set(reco);
}

/// Set the current object using instances of AliDmesonJetInfo as its source
///
/// \param truth A const reference to a AliDmesonJetInfo object representing the generated D meson
void AliAnalysisTaskDmesonJetsDetectorResponse::AliD0MatchInfoSummary::SetGenerated(const AliDmesonJetInfo& truth)
{
  fGenerated.Set(truth);
}

// Definitions of class AliAnalysisTaskDmesonJetsDetectorResponse::AliDStarMatchInfoSummary

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJetsDetectorResponse::AliDStarMatchInfoSummary);
/// \endcond

/// Reset the object
void AliAnalysisTaskDmesonJetsDetectorResponse::AliDStarMatchInfoSummary::Reset()
{
  fGenerated.Reset();
  fReconstructed.Reset();
}

/// Set the current object using instances of AliDmesonJetInfo as its source
///
/// \param reco A const reference to a AliDmesonJetInfo object representing the reconstructed D meson
void AliAnalysisTaskDmesonJetsDetectorResponse::AliDStarMatchInfoSummary::SetReconstructed(const AliDmesonJetInfo& reco)
{
  fReconstructed.Set(reco);
}

/// Set the current object using instances of AliDmesonJetInfo as its source
///
/// \param truth A const reference to a AliDmesonJetInfo object representing the generated D meson
void AliAnalysisTaskDmesonJetsDetectorResponse::AliDStarMatchInfoSummary::SetGenerated(const AliDmesonJetInfo& truth)
{
  fGenerated.Set(truth);
}

// Definitions of class AliAnalysisTaskDmesonJetsDetectorResponse::ResponseEngine

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJetsDetectorResponse::ResponseEngine);
/// \endcond

/// Default constructor, for ROOT I/O
AliAnalysisTaskDmesonJetsDetectorResponse::ResponseEngine::ResponseEngine() :
  TObject(),
  fCandidateType(kD0toKpi),
  fInhibit(kFALSE),
  fMaxJetDmesonDistance(0),
  fName(),
  fTree(0),
  fCurrentDmeson(0),
  fCurrentJetInfoReco(0),
  fCurrentJetInfoTruth(0),
  fDataSlotNumber(-1),
  fRecontructed(0),
  fGenerated(0)
{

}

/// Standard constructor with candidate type
///
/// \param type D meson candidate type (D0, D*, ...)
AliAnalysisTaskDmesonJetsDetectorResponse::ResponseEngine::ResponseEngine(ECandidateType_t type) :
  TObject(),
  fCandidateType(type),
  fInhibit(kFALSE),
  fMaxJetDmesonDistance(1),
  fName(),
  fTree(0),
  fCurrentDmeson(0),
  fCurrentJetInfoReco(0),
  fCurrentJetInfoTruth(0),
  fDataSlotNumber(-1),
  fRecontructed(0),
  fGenerated(0)
{

}

/// Copy constructor
///
/// \param source Reference to a valid AnalysisEngine to copy from.
AliAnalysisTaskDmesonJetsDetectorResponse::ResponseEngine::ResponseEngine(const AliAnalysisTaskDmesonJetsDetectorResponse::ResponseEngine &source) :
  TObject(source),
  fCandidateType(source.fCandidateType),
  fInhibit(source.fInhibit),
  fMaxJetDmesonDistance(1),
  fName(source.fName),
  fTree(0),
  fCurrentDmeson(0),
  fCurrentJetInfoReco(0),
  fCurrentJetInfoTruth(0),
  fDataSlotNumber(source.fDataSlotNumber),
  fRecontructed(source.fRecontructed),
  fGenerated(source.fGenerated)
{
}

/// Assignement operator
///
/// \param source Reference to a valid AnalysisEngine to copy from.
AliAnalysisTaskDmesonJetsDetectorResponse::ResponseEngine& AliAnalysisTaskDmesonJetsDetectorResponse::ResponseEngine::operator=(const ResponseEngine& source)
{
  new (this) ResponseEngine(source);
  return *this;
}

/// Checks if the response engine is properly set up. Otherwise inhibit the engine.
///
/// \param kTRUE if successful
Bool_t AliAnalysisTaskDmesonJetsDetectorResponse::ResponseEngine::CheckInit()
{
  fInhibit = (fRecontructed == 0 || fGenerated == 0);
  fName = fRecontructed->GetCandidateName();
  return !fInhibit;
}

/// Run the requested analysis for the current event
void AliAnalysisTaskDmesonJetsDetectorResponse::ResponseEngine::RunAnalysis()
{

}

/// Builds the tree where the output will be posted. Called from UserCreateOutputObject
///
/// \param taskName Name of the underlying task
/// \return Newly created tree
TTree* AliAnalysisTaskDmesonJetsDetectorResponse::ResponseEngine::BuildTree(const char* taskName)
{
  TString classname;
  switch (fCandidateType) {
  case kD0toKpi:
  case kD0toKpiLikeSign:
    classname = "AliAnalysisTaskDmesonJetsDetectorResponse::AliD0MatchInfoSummary";
    fCurrentDmeson = new AliD0MatchInfoSummary();
    break;
  case kDstartoKpipi:
    classname = "AliAnalysisTaskDmesonJetsDetectorResponse::AliDStarMatchInfoSummary";
    fCurrentDmeson = new AliDStarMatchInfoSummary();
    break;
  }
  TString treeName = TString::Format("%s_%s", taskName, GetName());
  fTree = new TTree(treeName, treeName);
  fTree->Branch("DmesonJet", classname, &fCurrentDmeson);

  fCurrentJetInfoTruth = new AliJetInfoSummary*[fGenerated->GetJetDefinitions().size()];
  for (Int_t i = 0; i < fGenerated->GetJetDefinitions().size(); i++) {
    if (fGenerated->GetJetDefinitions()[i].GetRhoName().IsNull()) {
      fCurrentJetInfoTruth[i] = new AliJetInfoSummary();
      TString bname = TString::Format("%s_truth", fGenerated->GetJetDefinitions()[i].GetName());
      fTree->Branch(bname, "AliAnalysisTaskDmesonJets::AliJetInfoSummary", &fCurrentJetInfoTruth[i]);
    }
    else {
      fCurrentJetInfoTruth[i] = new AliJetInfoPbPbSummary();
      TString bname = TString::Format("%s_truth", fGenerated->GetJetDefinitions()[i].GetName());
      fTree->Branch(bname, "AliAnalysisTaskDmesonJets::AliJetInfoPbPbSummary", &fCurrentJetInfoTruth[i]);
    }
  }

  fCurrentJetInfoReco = new AliJetInfoSummary*[fRecontructed->GetJetDefinitions().size()];
  for (Int_t i = 0; i < fRecontructed->GetJetDefinitions().size(); i++) {
    if (fGenerated->GetJetDefinitions()[i].GetRhoName().IsNull()) {
      fCurrentJetInfoReco[i] = new AliJetInfoSummary();
      TString bname = TString::Format("%s_reco", fRecontructed->GetJetDefinitions()[i].GetName());
      fTree->Branch(bname, "AliAnalysisTaskDmesonJets::AliJetInfoSummary", &fCurrentJetInfoReco[i]);
    }
    else {
      fCurrentJetInfoReco[i] = new AliJetInfoPbPbSummary();
      TString bname = TString::Format("%s_reco", fRecontructed->GetJetDefinitions()[i].GetName());
      fTree->Branch(bname, "AliAnalysisTaskDmesonJets::AliJetInfoPbPbSummary", &fCurrentJetInfoReco[i]);
    }
  }

  return fTree;
}

/// Loops over the analysis output and posts the output in the tree. Apply kinematic cuts if requested.
///
/// \param applyKinCuts Whether kinematic cuts should be applied
/// \return Always kTRUE.
Bool_t AliAnalysisTaskDmesonJetsDetectorResponse::ResponseEngine::FillTree(Bool_t applyKinCuts)
{
  fRecontructed->FillQA(applyKinCuts);
  fGenerated->FillQA(applyKinCuts);

  std::map<int, AliDmesonJetInfo>& recoDmesons = fRecontructed->GetDmesons();
  std::map<int, AliDmesonJetInfo>& truthDmesons = fGenerated->GetDmesons();

  // Loop over reconstructed D meson and look for their generated counterparts
  for (auto& dmeson_reco : recoDmesons) {
    // Reset D meson and jet tree branches
    fCurrentDmeson->Reset();
    for (UInt_t ij = 0; ij < fRecontructed->GetJetDefinitions().size(); ij++) {
      fCurrentJetInfoReco[ij]->Reset();
    }
    for (UInt_t ij = 0; ij < fGenerated->GetJetDefinitions().size(); ij++) {
      fCurrentJetInfoTruth[ij]->Reset();
    }

    // Copy the reconstructed D meson information in the tree branch
    fCurrentDmeson->SetReconstructed(dmeson_reco.second);

    Int_t accRecJets = 0;
    Int_t accGenJets = 0;

    // Copy the reconstructed jet information in the tree branch
    for (UInt_t ij = 0; ij < fRecontructed->GetJetDefinitions().size(); ij++) {
      AliJetInfo* jet = dmeson_reco.second.GetJet(fRecontructed->GetJetDefinitions()[ij].GetName());
      if (!jet) continue;
      if (applyKinCuts && !fRecontructed->GetJetDefinitions()[ij].IsJetInAcceptance(*jet)) continue;
      fCurrentJetInfoReco[ij]->Set(dmeson_reco.second, fRecontructed->GetJetDefinitions()[ij].GetName());
      accRecJets++;
    }

    // Look for the generated D meson counterpart
    if (dmeson_reco.second.fMCLabel >= 0) {
      std::map<int, AliDmesonJetInfo>::iterator it = truthDmesons.find(dmeson_reco.second.fMCLabel);
      if (it != truthDmesons.end()) {
        std::pair<const int, AliDmesonJetInfo>& dmeson_truth = (*it);
        // Flag the generated D meson as reconstructed
        dmeson_truth.second.fReconstructed = kTRUE;
        // Copy the generated D meson information in the tree branch
        fCurrentDmeson->SetGenerated((*it).second);

        // Copy the reconstructed jet information in the tree branch
        for (UInt_t ij = 0; ij < fGenerated->GetJetDefinitions().size(); ij++) {
          AliJetInfo* jet = dmeson_truth.second.GetJet(fGenerated->GetJetDefinitions()[ij].GetName());
          if (!jet) continue;
          if (!applyKinCuts || fGenerated->GetJetDefinitions()[ij].IsJetInAcceptance(*jet)) accGenJets++;
          fCurrentJetInfoTruth[ij]->Set(dmeson_truth.second, fGenerated->GetJetDefinitions()[ij].GetName());
        }
      }
    }

    // Fill the tree with the current D meson
    if (accRecJets > 0 || accGenJets > 0) fTree->Fill();
  }

  // Reset jet tree branches
  for (UInt_t ij = 0; ij < fRecontructed->GetJetDefinitions().size(); ij++) fCurrentJetInfoReco[ij]->Reset();
  for (UInt_t ij = 0; ij < fGenerated->GetJetDefinitions().size(); ij++) fCurrentJetInfoTruth[ij]->Reset();

  // Loop over generated D meson that were not reconstructed
  for (auto& dmeson_truth : truthDmesons) {
    // Reset D meson tree branch
    fCurrentDmeson->Reset();
    // Skip if the generated D meson was reconstructed
    if (dmeson_truth.second.fReconstructed) continue;
    // Copy the generated D meson information in the tree branch
    fCurrentDmeson->SetGenerated(dmeson_truth.second);

    Int_t accRecJets = 0;
    Int_t accGenJets = 0;

    // Copy the reconstructed jet information in the tree branch
    for (UInt_t ij = 0; ij < fGenerated->GetJetDefinitions().size(); ij++) {
      fCurrentJetInfoTruth[ij]->Reset();
      AliJetInfo* jet = dmeson_truth.second.GetJet(fGenerated->GetJetDefinitions()[ij].GetName());
      if (!jet) continue;
      if (applyKinCuts && !fGenerated->GetJetDefinitions()[ij].IsJetInAcceptance(*jet)) continue;
      fCurrentJetInfoTruth[ij]->Set(dmeson_truth.second, fGenerated->GetJetDefinitions()[ij].GetName());
      accGenJets++;
    }

    // Fill the tree with the current D meson
    if (accRecJets > 0 || accGenJets > 0) fTree->Fill();
  }

  return kTRUE;
}

// Definitions of class AliAnalysisTaskDmesonJetsDetectorResponse

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskDmesonJetsDetectorResponse);
/// \endcond

/// This is the default constructor, used for ROOT I/O purposes.
AliAnalysisTaskDmesonJetsDetectorResponse::AliAnalysisTaskDmesonJetsDetectorResponse() :
  AliAnalysisTaskDmesonJets(),
  fResponseEngines()
{
  fOutputType = kNoOutput;
}

/// This is the standard named constructor.
///
/// \param name Name of the task
AliAnalysisTaskDmesonJetsDetectorResponse::AliAnalysisTaskDmesonJetsDetectorResponse(const char* name, Int_t nOutputTrees) :
  AliAnalysisTaskDmesonJets(name, nOutputTrees),
  fResponseEngines()
{
  fOutputType = kNoOutput;
}

/// Creates the output containers.
void AliAnalysisTaskDmesonJetsDetectorResponse::UserCreateOutputObjects()
{
  ::Info("UserCreateOutputObjects", "CreateOutputObjects of task %s", GetName());

  AliAnalysisTaskDmesonJets::UserCreateOutputObjects();

  for (auto &param : fAnalysisEngines) {
    if (param.IsInhibit()) continue;
    if (param.GetMCMode() != kSignalOnly && param.GetMCMode() != kMCTruth) continue;

    std::map<ECandidateType_t, ResponseEngine>::iterator it = fResponseEngines.find(param.GetCandidateType());

    if (it == fResponseEngines.end()) {
      it = (fResponseEngines.insert(std::pair<const ECandidateType_t, ResponseEngine>(param.GetCandidateType(), ResponseEngine(param.GetCandidateType())))).first;
    }

    if (param.GetMCMode() == kMCTruth) {
      (*it).second.SetGeneratedAnalysisEngine(&param);
    }
    else {
      (*it).second.SetReconstructedAnalysisEngine(&param);
    }
  }

  Int_t treeSlot = 0;

  for (auto &resp : fResponseEngines) {
    if (!resp.second.CheckInit()) continue;

    resp.second.BuildTree(GetName());
    if (treeSlot < fNOutputTrees) {
      resp.second.AssignDataSlot(treeSlot+2);
      treeSlot++;
      PostDataFromResponseEngine(resp.second);
    }
    else {
      AliError(Form("Number of data output slots %d not sufficient. Tree of response engine %s will not be posted!", fNOutputTrees, resp.second.GetName()));
    }
  }
}

/// Does some specific initializations for the analysis engines,
/// then calls the base class ExecOnce() method.
void AliAnalysisTaskDmesonJetsDetectorResponse::ExecOnce()
{
  AliAnalysisTaskDmesonJets::ExecOnce();
}

/// Run the analysis
///
/// \return kTRUE on success
Bool_t AliAnalysisTaskDmesonJetsDetectorResponse::Run()
{
  return AliAnalysisTaskDmesonJets::Run();
}

/// Fill the histograms.
///
/// \return Always kTRUE
Bool_t AliAnalysisTaskDmesonJetsDetectorResponse::FillHistograms()
{
  for (auto &resp : fResponseEngines) {
    if (resp.second.IsInhibit()) continue;

    resp.second.FillTree(fApplyKinematicCuts);

    PostDataFromResponseEngine(resp.second);
  }

  if (fMCContainer) FillPartonLevelHistograms();
  return kTRUE;

}

/// This method overrides the base class method and forbids changing the output type
///
/// \param b Output type (none, tree, thn)
void AliAnalysisTaskDmesonJetsDetectorResponse::SetOutputTypeInternal(EOutputType_t b)
{
  if (fOutputType != kNoOutput && fOutputType != kTreeOutput) {
    AliWarning("This class only provides a tree output.");
  }

  // Always set it to kNoOutput: base class does not generate any output (other than QA histograms), all the output comes from the derived class
  fOutputType = kNoOutput;
}

/// Post the tree of an response engine in the data slot (if the tree exists and the data slot has been assigned)
///
/// \param eng Constant reference to a response engine
///
/// \return -1 if unsuccessful, an integer number corresponding to the data slot if successful
Int_t AliAnalysisTaskDmesonJetsDetectorResponse::PostDataFromResponseEngine(const ResponseEngine& eng)
{
  if (eng.GetDataSlotNumber() >= 0 && eng.GetTree()) {
    PostData(eng.GetDataSlotNumber(), eng.GetTree());
    return eng.GetDataSlotNumber();
  }
  else {
    return -1;
  }
}

/// Create an instance of this class and add it to the analysis manager
///
/// \param ntracks name of the track collection
/// \param nclusters name of the calorimeter cluster collection
/// \param nMCpart name of the MC particle collection
/// \param nMaxTrees number of output trees
/// \param suffix additional suffix that can be added at the end of the task name
/// \return pointer to the new AliAnalysisTaskDmesonJets task
AliAnalysisTaskDmesonJetsDetectorResponse* AliAnalysisTaskDmesonJetsDetectorResponse::AddTaskDmesonJetsDetectorResponse(TString trackName, TString clusName, TString mcPartName, Int_t nMaxTrees, TString suffix)
{
  // Get the pointer to the existing analysis manager via the static access method.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskDmesonJetsDetectorResponse", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler) {
    ::Error("AddTaskDmesonJetsDetectorResponse", "This task requires an input event handler");
    return 0;
  }

  EDataType_t dataType = kUnknownDataType;

  if (handler->InheritsFrom("AliESDInputHandler")) {
    dataType = kESD;
  }
  else if (handler->InheritsFrom("AliAODInputHandler")) {
    dataType = kAOD;
  }

  // Init the task and do settings
  if (trackName == "usedefault") {
    if (dataType == kESD) {
      trackName = "Tracks";
    }
    else if (dataType == kAOD) {
      trackName = "tracks";
    }
    else {
      trackName = "";
    }
  }

  if (clusName == "usedefault") {
    if (dataType == kESD) {
      clusName = "CaloClusters";
    }
    else if (dataType == kAOD) {
      clusName = "caloClusters";
    }
    else {
      clusName = "";
    }
  }

  if (mcPartName == "usedefault") {
    mcPartName = "mcparticles"; // Always needs AliAODMCParticle objects
  }

  TString name("AliAnalysisTaskDmesonJetsDetectorResponse");
  if (!suffix.IsNull()) {
    name += TString::Format("_%s", suffix.Data());
  }

  AliAnalysisTaskDmesonJetsDetectorResponse* jetTask = new AliAnalysisTaskDmesonJetsDetectorResponse(name, nMaxTrees);

  if (!trackName.IsNull()) {
    AliHFTrackContainer* trackCont = new AliHFTrackContainer(trackName);
    jetTask->AdoptParticleContainer(trackCont);
  }

  if (!mcPartName.IsNull()) {
    AliMCParticleContainer* partCont = new AliHFAODMCParticleContainer(mcPartName);
    partCont->SetEtaLimits(-1.5, 1.5);
    partCont->SetPtLimits(0, 1000);
    jetTask->AdoptParticleContainer(partCont);
  }

  jetTask->AddClusterContainer(clusName.Data());

  // Final settings, pass to manager and set the containers
  mgr->AddTask(jetTask);

  // Create containers for input/output
  AliAnalysisDataContainer* cinput1 = mgr->GetCommonInputContainer();
  TString contname1(name);
  contname1 += "_histos";
  AliAnalysisDataContainer* coutput1 = mgr->CreateContainer(contname1.Data(),
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s", AliAnalysisManager::GetCommonFileName()));

  mgr->ConnectInput(jetTask, 0, cinput1);
  mgr->ConnectOutput(jetTask, 1, coutput1);

  for (Int_t i = 0; i < nMaxTrees; i++) {
    TString contname = TString::Format("%s_tree_%d", name.Data(), i);
    AliAnalysisDataContainer *coutput = mgr->CreateContainer(contname.Data(),
        TTree::Class(),AliAnalysisManager::kOutputContainer,
        Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(jetTask, 2+i, coutput);
  }
  return jetTask;
}

