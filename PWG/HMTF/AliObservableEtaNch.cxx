#include <iostream>

#include "TString.h"

#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"

#include "AliObservableBase.h"
#include "AliObservableEtaNch.h"
#include "AliEventClassifierBase.h"

#include "AliIsPi0PhysicalPrimary.h"

using namespace std;

ClassImp(AliObservableEtaNch)

AliObservableEtaNch::AliObservableEtaNch()
  : AliObservableBase()
{
}


AliObservableEtaNch::AliObservableEtaNch(AliEventClassifierBase *classifier)
  :  AliObservableBase("feta_Nch", "#eta vs. classifier vs. Nch")
{
  fclassifier = classifier;
  const Float_t eta_max = 10;
  const Int_t eta_bins  = 200;
  const Int_t classifier_bins  = 250;

  fhistogram = new TH2F("eta_classifier_" + TString(classifier->GetName()),
			"#eta vs. classifier vs. N_{ch}",
			eta_bins, -eta_max, eta_max,
			classifier_bins, fclassifier->GetExpectedMinValue(), fclassifier->GetExpectedMaxValue());
  fhistogram->GetXaxis()->SetTitle("#eta");
  fhistogram->GetYaxis()->SetTitle(classifier->GetTitle());
  fhistogram->GetZaxis()->SetTitle("N_{ch}");
  fhistogram->Sumw2();
  fhistogram->SetDirectory(0);

  // Add this histogram to the "folder" of the classifier in which it is binned
  classifier->GetClassifierOutputList()->Add(fhistogram);
}

void AliObservableEtaNch::Fill(AliMCEvent *event, AliStack *stack) {
  Double_t classifier_value = fclassifier->GetClassifierValue(event, stack);
  Double_t event_weight = event->GenEventHeader()->EventWeight();

  for (Int_t iTrack = 0; iTrack < event->GetNumberOfTracks(); iTrack++) {
    AliMCParticle *track = static_cast<AliMCParticle*>(event->GetTrack(iTrack));
    // load track
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTrack);
      continue;
    }
    // discard unphysical particles from some generators
    if (track->Pt() == 0 || track->E() <= 0)
      continue;

    // is it a primary particle or a pi0? Else, skip it.
    if (!stack->IsPhysicalPrimary(iTrack)) continue;

    // we want only charged particles here!
    if (track->Charge() == 0) continue;

    // Ok, lets fill!
    fhistogram->Fill(track->Eta(), classifier_value, event_weight);
  }
}
