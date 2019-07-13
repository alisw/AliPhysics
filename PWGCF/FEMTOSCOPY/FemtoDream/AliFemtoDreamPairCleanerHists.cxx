/*
 * AliFemtoPairCleanerHists.cxx
 *
 *  Created on: 22 Dec 2017
 *      Author: bernhardhohlweger
 */

#include "AliFemtoDreamPairCleanerHists.h"
ClassImp(AliFemtoDreamPairCleanerHists)
AliFemtoDreamPairCleanerHists::AliFemtoDreamPairCleanerHists()
    : fTrackDecays(nullptr),
      fDecayDecays(nullptr),
      fPairInvMass(nullptr),
      fOutput(0) {
}

AliFemtoDreamPairCleanerHists::AliFemtoDreamPairCleanerHists(
    const AliFemtoDreamPairCleanerHists& hists)
    : fTrackDecays(hists.fTrackDecays),
      fDecayDecays(hists.fDecayDecays),
      fPairInvMass(hists.fPairInvMass),
      fOutput(hists.fOutput) {
}

AliFemtoDreamPairCleanerHists::AliFemtoDreamPairCleanerHists(
    int nTrackDecays, int nDecayDecays)
    : fTrackDecays(nullptr),
      fDecayDecays(nullptr),
      fPairInvMass(nullptr),
      fOutput(0)
{
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("PairCleaner");

  fTrackDecays = new TH1F*[nTrackDecays];
  for (int i = 0; i < nTrackDecays; ++i) {
    TString histName = Form("DaugthersSharedTracks_%d", i);
    fTrackDecays[i] = new TH1F(histName.Data(), histName.Data(), 20, 0, 20);
    fOutput->Add(fTrackDecays[i]);
  }
  if (nDecayDecays > 0)
    fDecayDecays = new TH1F*[nDecayDecays];
  for (int i = 0; i < nDecayDecays; ++i) {
    TString histName = Form("DaugthersSharedDaughters_%d", i);
    fDecayDecays[i] = new TH1F(histName.Data(), histName.Data(), 20, 0, 20);
    fOutput->Add(fDecayDecays[i]);
  }
}

AliFemtoDreamPairCleanerHists& AliFemtoDreamPairCleanerHists::operator=(
    const AliFemtoDreamPairCleanerHists& hists) {
  if (this != &hists) {
    this->fTrackDecays = hists.fTrackDecays;
    this->fDecayDecays = hists.fDecayDecays;
    this->fPairInvMass = hists.fPairInvMass;
    this->fOutput = hists.fOutput;
  }
  return *this;
}

AliFemtoDreamPairCleanerHists::~AliFemtoDreamPairCleanerHists() {
  // TODO Auto-generated destructor stub
}

