#include "AliFemtoDreamBaseDump.h"
#include "AliFemtoDreamHigherPairMath.h"
#include "AliAnalysisManager.h"
#include "TFile.h"

// ------------------------------------------------------------------------------------------------
ClassImp(AliFemtoDreamPairDump)

AliFemtoDreamPairDump::AliFemtoDreamPairDump()
    : fkstar(0),
      fPart1(),
      fPart2() {
}

AliFemtoDreamPairDump::AliFemtoDreamPairDump(AliFemtoDreamBasePart part1,
                                             AliFemtoDreamBasePart part2,
                                             float kstar)
    : fkstar(kstar),
      fPart1(),
      fPart2() {
  part1.KillGlobalTrackArray();  // set the global track array to nullptr
  fPart1 = part1;
  part2.KillGlobalTrackArray();  // set the global track array to nullptr
  fPart2 = part2;
}

AliFemtoDreamPairDump::AliFemtoDreamPairDump(const AliFemtoDreamPairDump& pair)
    : fkstar(pair.fkstar),
      fPart1(pair.fPart1),
      fPart2(pair.fPart2) {
}

AliFemtoDreamPairDump &AliFemtoDreamPairDump::operator=(
    const AliFemtoDreamPairDump &obj) {
  if (this == &obj) {
    return *this;
  }
  fkstar = obj.fkstar;
  fPart1 = obj.fPart1;
  fPart2 = obj.fPart2;
  return (*this);
}

// ------------------------------------------------------------------------------------------------
ClassImp(AliFemtoDreamEventDump)

AliFemtoDreamEventDump::AliFemtoDreamEventDump()
    : fFilePath(),
      fMultiplicity(0),
      fVertexZ(0),
      fPairDump() {
}

void AliFemtoDreamEventDump::Reset() {
  fFilePath = "";
  fMultiplicity = 0;
  fVertexZ = 0;
  fPairDump.resize(0);
}

void AliFemtoDreamEventDump::SetEventProperties(float mult, float zvertex,
                                                TString filePath) {
  fFilePath = filePath;
  fMultiplicity = mult;
  fVertexZ = zvertex;
}

// ------------------------------------------------------------------------------------------------
ClassImp(AliFemtoDreamDump)

AliFemtoDreamDump::AliFemtoDreamDump()
    : fkstarThreshold(0.4),
      fOutputTree(nullptr),
      fEventDump(nullptr) {
}

AliFemtoDreamDump::AliFemtoDreamDump(const char *name)
    : fkstarThreshold(0.4),
      fOutputTree(nullptr),
      fEventDump(nullptr) {
  fOutputTree = new TTree(name, name);
  fOutputTree->SetDirectory(0);

  fEventDump = new AliFemtoDreamEventDump();
  fOutputTree->Branch("EventDump", &fEventDump);
}

AliFemtoDreamDump::~AliFemtoDreamDump() {
  delete fOutputTree;
}

void AliFemtoDreamDump::SetEvent(std::vector<AliFemtoDreamBasePart> &vec1,
                                 std::vector<AliFemtoDreamBasePart> &vec2,
                                 AliFemtoDreamEvent *evt, const int pdg1,
                                 const int pdg2) {
  if (vec1.size() == 0 && vec2.size() == 0) {
    return;
  }

  TString filePath;
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if (man) {
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man
        ->GetInputEventHandler());
    if (inputHandler) {
      filePath = inputHandler->GetTree()->GetCurrentFile()->GetName();
    }
  }

  // Reset the event and set global event properties
  fEventDump->Reset();
  fEventDump->SetEventProperties(evt->GetMultiplicity(), evt->GetZVertex(),
                                 filePath);

  for (auto part1 : vec1) {
    if (!part1.UseParticle()) {
      continue;
    }
    for (auto part2 : vec2) {
      if (!part2.UseParticle()) {
        continue;
      }
      const float kstar = AliFemtoDreamHigherPairMath::RelativePairMomentum(
          &part1, pdg1, &part2, pdg2);
      if (kstar < fkstarThreshold) {
        fEventDump->AddPair(part1, part2, kstar);
      }
    }
  }

  // Fill the tree when we found some particles
  if (fEventDump->GetNPairs() > 0) {
    fOutputTree->Fill();
  }
}

void AliFemtoDreamDump::SetEvent(std::vector<AliFemtoDreamBasePart> &vec,
                                 AliFemtoDreamEvent *evt, const int pdg) {
  if (vec.size() == 0) {
    return;
  }

  TString filePath;
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  if (man) {
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man
        ->GetInputEventHandler());
    if (inputHandler) {
      filePath = inputHandler->GetTree()->GetCurrentFile()->GetName();
    }
  }

  // Reset the event and set global event properties
  fEventDump->Reset();
  fEventDump->SetEventProperties(evt->GetMultiplicity(), evt->GetZVertex(),
                                 filePath);

  for (auto iterPart1 = vec.begin(); iterPart1 < vec.end(); ++iterPart1) {
    auto part1 = *iterPart1;
    if (!part1.UseParticle()) {
      continue;
    }
    for (auto iterPart2 = iterPart1 + 1; iterPart2 < vec.end(); ++iterPart2) {
      auto part2 = *iterPart2;
      if (!part2.UseParticle()) {
        continue;
      }
      const float kstar = AliFemtoDreamHigherPairMath::RelativePairMomentum(
          &part1, pdg, &part2, pdg);
      if (kstar < fkstarThreshold) {
        fEventDump->AddPair(part1, part2, kstar);
      }
    }
  }

  // Fill the tree
  if (fEventDump->GetNPairs() > 0) {
    fOutputTree->Fill();
  }
}
