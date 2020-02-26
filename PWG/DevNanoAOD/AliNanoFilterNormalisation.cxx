#include "AliNanoFilterNormalisation.h"
#include <TCollection.h>
#include <TH2D.h>

#include <string>

ClassImp(AliNanoFilterNormalisation);

AliNanoFilterNormalisation::AliNanoFilterNormalisation(TString name, TString title, int nMultBins, float multBegin, float multEnd) : TNamed(name, title) {
  fCandidateEvents = new TH2D("fCandidateEvents", ";Multiplicity estimator;", nMultBins, multBegin, multEnd, 5, -0.5, 4.5);
  fSelectedEvents  = (TH2D*) fCandidateEvents->Clone("fSelectedEvents");

  std::string labels[5]{"Input events", "Triggered events", "Triggered + Quality cuts", "Triggered + QC + Reco vertex", "Analysis events"};
  for (int iType{0}; iType < 5; ++iType) {
    fCandidateEvents->GetYaxis()->SetBinLabel(iType + 1, labels[iType].data());
    fSelectedEvents->GetYaxis()->SetBinLabel(iType + 1, labels[iType].data());
  }
}

AliNanoFilterNormalisation::AliNanoFilterNormalisation(TString name, TString title, int nMultBins, float* mBins) : TNamed(name, title) {
  float selBins[5]{-0.5,0.5,1.5,2.5,3.5,4.5};
  fCandidateEvents = new TH2D("fCandidateEvents", ";Multiplicity estimator;", nMultBins, mBins, 5, selBins);
  fSelectedEvents  = (TH2D*) fCandidateEvents->Clone("fSelectedEvents");

  std::string labels[5]{"Input events", "Triggered events", "Triggered + Quality cuts", "Triggered + QC + Reco vertex", "Analysis events"};
  for (int iType{0}; iType < 5; ++iType) {
    fCandidateEvents->GetYaxis()->SetBinLabel(iType + 1, labels[iType].data());
    fSelectedEvents->GetYaxis()->SetBinLabel(iType + 1, labels[iType].data());
  }
}

AliNanoFilterNormalisation::~AliNanoFilterNormalisation() {
  delete fCandidateEvents;
  delete fSelectedEvents;
}

Long64_t AliNanoFilterNormalisation::Merge(TCollection* hlist) {
  if (hlist) {
    TIter nxh(hlist);
    AliNanoFilterNormalisation* xh;
    while ((xh = (AliNanoFilterNormalisation *)nxh())) {
      fCandidateEvents->Add(xh->GetCandidateEventsHistogram());
      fSelectedEvents->Add(xh->GetSelectedEventsHistogram());
    }
    return hlist->GetEntries();
  }
  return 0ll;
}

double AliNanoFilterNormalisation::GetScalingFactor(NormBin ybin, float mult) {
  int centBin = fCandidateEvents->GetXaxis()->FindBin(mult);
  double cand = fCandidateEvents->GetBinContent(centBin, ybin);
  double sel = fSelectedEvents->GetBinContent(centBin, ybin);
  if (cand < 1)
    return 0;
  else
    return sel / cand;
}

void AliNanoFilterNormalisation::FillCandidate(bool triggered, bool nonVertexRelatedSel, bool recoVertex, bool allCuts, float mult) {
  fCandidateEvents->Fill(mult, kAnyEvent);
  if (triggered) fCandidateEvents->Fill(mult, kTriggeredEvent);
  if (nonVertexRelatedSel) fCandidateEvents->Fill(mult, kTriggeredEventWithQualityCuts);
  if (recoVertex) fCandidateEvents->Fill(mult, kTriggeredEventWithQualityCutsAndRecoVertex);
  if (allCuts) fCandidateEvents->Fill(mult, kAnalysisEvent);
}

void AliNanoFilterNormalisation::FillSelected(bool triggered, bool nonVertexRelatedSel, bool recoVertex, bool allCuts, float mult) {
  fSelectedEvents->Fill(mult, kAnyEvent);
  if (triggered) fSelectedEvents->Fill(mult, kTriggeredEvent);
  if (nonVertexRelatedSel) fSelectedEvents->Fill(mult, kTriggeredEventWithQualityCuts);
  if (recoVertex) fSelectedEvents->Fill(mult, kTriggeredEventWithQualityCutsAndRecoVertex);
  if (allCuts) fSelectedEvents->Fill(mult, kAnalysisEvent);
}

double AliNanoFilterNormalisation::GetNcanditateEvents(NormBin ybin, float mult) {
  int centBin = fCandidateEvents->GetXaxis()->FindBin(mult);
  return fCandidateEvents->GetBinContent(centBin, ybin);
}

double AliNanoFilterNormalisation::GetNselectedEvents(NormBin ybin, float mult) {
  int centBin = fSelectedEvents->GetXaxis()->FindBin(mult);
  return fSelectedEvents->GetBinContent(centBin, ybin);
}
