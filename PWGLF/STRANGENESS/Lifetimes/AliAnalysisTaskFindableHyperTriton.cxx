#include "AliAnalysisTaskFindableHyperTriton.h"

#include <array>
#include <algorithm>
#include <cmath>
#include <unordered_map>

#include <TChain.h>
#include <TFile.h>
#include <TParticle.h>

#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliPDG.h"
#include "AliPID.h"
#include "AliVVertex.h"

ClassImp(AliAnalysisTaskFindableHyperTriton);

namespace
{
  constexpr double Sq(double x) { return x * x; }
  constexpr float kEps = 1.e-6;

  template<int N> int FillVector(std::vector<FindableHyperTriton<N> >& vec, TParticle* mom, TParticle* prong) {
    FindableHyperTriton<N> hyper;
    std::array<double, 3> decayVtx = {prong->Vx(), prong->Vy(), prong->Vz()};
    std::array<double, 3> momMom = {mom->Px(), mom->Py(), mom->Pz()};
    std::copy(decayVtx.begin(), decayVtx.end(), hyper.fDecayVertex);
    std::copy(momMom.begin(), momMom.end(), hyper.fMomentum);
    hyper.fDeltaT = prong->T() - mom->T();
    hyper.fFoundTracks = 0;
    hyper.SetSign(mom->GetPdgCode() > 0);
    vec.push_back(hyper);
    return vec.size() - 1;
  }

  template<int N> void UpdateElement(FindableHyperTriton<N> & hyper, TParticle* mom, TParticle* prong, AliESDtrack* track) {
    int index = std::abs(track->GetLabel()) - mom->GetFirstDaughter();
    hyper.fTracks[index] = new AliESDtrack(*track);
    hyper.fPDG[index] = prong->GetPdgCode();
    hyper.fFoundTracks |= BIT(index);
  }
} // namespace

AliAnalysisTaskFindableHyperTriton::AliAnalysisTaskFindableHyperTriton(std::string name)
  : AliAnalysisTaskSE(name.data()),
  fEventCuts{},
  fEventSummary{},
  f2Body{},
  f3Body{}
{
  // Standard output
  DefineInput(0, TChain::Class());
  DefineOutput(1, TTree::Class()); // Tree output
}

AliAnalysisTaskFindableHyperTriton::~AliAnalysisTaskFindableHyperTriton()
{
  if (fTree)
    delete fTree;
}

void AliAnalysisTaskFindableHyperTriton::UserCreateOutputObjects()
{
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler = (AliInputEventHandler *)(man->GetInputEventHandler());
  inputHandler->SetNeedField();

  fTree = new TTree("fTree", "Findable hypertritons");
  fTree->Branch("Event", &fEventSummary);
  fTree->Branch("HyperTriton2body", &f2Body);
  fTree->Branch("HyperTriton3body", &f3Body);

  PostData(1, fTree);

  AliPDG::AddParticlesToPdgDataBase();
} // end UserCreateOutputObjects

void AliAnalysisTaskFindableHyperTriton::UserExec(Option_t *)
{
  AliESDEvent *esdEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!esdEvent)
    ::Fatal("AliAnalysisTaskFindableHyperTriton::UserExec",
        "AliESDEvent not found.");

  AliMCEvent *mcEvent = MCEvent();
  if (!mcEvent)
    ::Fatal("AliAnalysisTaskFindableHyperTriton::UserExec", "Could not retrieve MC event");

  if (!fEventCuts.AcceptEvent(esdEvent))
    return;

  std::array<double, 3> primaryVertex;
  fEventSummary.fMultiplicity = fEventCuts.GetCentrality();
  fEventCuts.GetPrimaryVertex()->GetXYZ(primaryVertex.data());

  std::copy(primaryVertex.begin(), primaryVertex.end(), fEventSummary.fRecVertex);

  const AliVVertex *mcV = mcEvent->GetPrimaryVertex();
  mcV->GetXYZ(primaryVertex.data());
  std::copy(primaryVertex.begin(), primaryVertex.end(), fEventSummary.fTrueVertex);
  fEventSummary.fMagField = esdEvent->GetMagneticField();

  std::unordered_map<int, int> threeBodyMom;
  std::unordered_map<int, int> twoBodyMom;
  f2Body.clear();
  f3Body.clear();
  for (int ilab = 0;  ilab < mcEvent->GetNumberOfTracks(); ilab++) {   // This is the begining of the loop on tracks
    TParticle* part = mcEvent->Particle( ilab );
    if(!part) {
      ::Warning("AliAnalysisTaskFindableHyperTriton::UserExec","Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skipping.", ilab );
      continue;
    }
    if (!mcEvent->IsPhysicalPrimary(ilab))
      continue;
    if (std::abs(part->GetPdgCode()) != 1010010030)
      continue;
    TParticle* prong = nullptr;
    int nProngs = 0;
    for (int iProng = part->GetFirstDaughter(); iProng < part->GetLastDaughter(); ++iProng) {
      if (!mcEvent->IsSecondaryFromWeakDecay(iProng))
        continue;
      prong = mcEvent->Particle(iProng);
      if (std::abs(prong->GetPdgCode()) == 11)
        continue;
      nProngs++;
    }
    if (nProngs > 3 && fBreakOnMultiBody)
      ::Fatal("AliAnalysisTaskFindableHyperTriton::UserExec","Multi body found... exiting");
    if (nProngs == 3)
      threeBodyMom[ilab] = FillVector(f3Body, part, prong);
    else
      twoBodyMom[ilab] = FillVector(f2Body, part, prong);
  }

  for (int iEv = 0; iEv < esdEvent->GetNumberOfTracks(); ++iEv) {
    AliESDtrack *track = esdEvent->GetTrack(iEv);
    int label = std::abs(track->GetLabel());

    if (mcEvent->IsPhysicalPrimary(label) ||
        mcEvent->IsSecondaryFromMaterial(label))
      continue;
    TParticle *prong = mcEvent->Particle(label);
    int momLabel = prong->GetFirstMother();
    TParticle *mom = mcEvent->Particle(momLabel);
    if (std::abs(mom->GetPdgCode()) != 1010010030)
      continue;

    bool newElement2B = twoBodyMom.find(momLabel) == twoBodyMom.end();
    bool newElement3B = threeBodyMom.find(momLabel) == threeBodyMom.end();
    if ((newElement2B && newElement3B) || (!newElement2B && !newElement3B))
      ::Fatal("AliAnalysisTaskFindableHyperTriton::UserExec","Unexpected additional MC hypertriton %d", newElement2B && newElement3B);
    else {
      auto element = newElement2B ? threeBodyMom[momLabel] : twoBodyMom[momLabel];
      if (newElement2B)
        UpdateElement(f3Body[element], mom, prong, track);
      else
        UpdateElement(f2Body[element], mom, prong, track);
    }
  }

  fTree->Fill();
  PostData(1, fTree);
}

AliAnalysisTaskFindableHyperTriton* AliAnalysisTaskFindableHyperTriton::AddTask(std::string tskname, std::string suffix) {
  // Get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Fatal("AliAnalysisTaskFindableHyperTriton::AddTaskFindableHyperTriton", "No analysis manager found.");
    return nullptr;
  }

  // Check the analysis type using the event handlers connected to the analysis
  // manager.
  if (!mgr->GetInputEventHandler()) {
    ::Fatal("AliAnalysisTaskFindableHyperTriton::AddTaskFindableHyperTriton",
            "This task requires an input event handler");
    return nullptr;
  }

  tskname.append(suffix.data());
  AliAnalysisTaskFindableHyperTriton *task =
      new AliAnalysisTaskFindableHyperTriton(tskname.data());

  AliAnalysisDataContainer *coutput1 =
      mgr->CreateContainer(Form("FindableTrees%s", suffix.data()), TTree::Class(),
                           AliAnalysisManager::kOutputContainer, "FindableTrees.root");
  coutput1->SetSpecialOutput();

  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task, 1, coutput1);
  return task;
}
