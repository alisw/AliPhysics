#include <iostream>

#include "TString.h"
#include "TH3F.h"

#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"

#include "AliObservableBase.h"
#include "AliObservableClassifierpTPID.h"
#include "AliEventClassifierBase.h"

#include "AliIsPi0PhysicalPrimary.h"

using namespace std;

ClassImp(AliObservableClassifierpTPID)

AliObservableClassifierpTPID::AliObservableClassifierpTPID()
  : AliObservableBase()
{
}


AliObservableClassifierpTPID::AliObservableClassifierpTPID(AliEventClassifierBase *classifier)
  :  AliObservableBase("fClassifier_pT_PID", "classifier vs. p_{T} vs. PID")
{
  fclassifier = classifier;
  const Int_t classifier_min   = fclassifier->GetExpectedMinValue();
  const Int_t classifier_max   = fclassifier->GetExpectedMaxValue();
  const Int_t classifier_bins  = 250;
  const Float_t step_size = Float_t(classifier_max - classifier_min) / Float_t(classifier_bins);
  
  const Float_t pt_bin_edges[] =
    {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8, 0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,  // 0.1 * 20
     2.2,2.4,2.6,2.8,3.0,                                                                // 0.2 * 5
     3.3,3.6,3.9,4.2,                                                                    // 0.3 * 4
     4.6,5,5.4,                                                                          // 0.4 * 3
     5.9,
     6.5,7,7.5,8,8.5,
     9.2,
     10,11,12,
     13.5,15,
     17,20
    };

  // Construct array for classifier values
  Float_t classifier_bin_edges[classifier_bins + 1];
  Float_t current_edge = classifier_min;
  for (Int_t idx = classifier_min; idx < classifier_bins+1; idx++, current_edge += step_size){
    classifier_bin_edges[idx] = current_edge;
  }

  // Construct the pid bin edges (each particle type has its own bin)
  Float_t pid_bin_edges[kNPID + 1];
  for (Int_t idx = 0; idx < kNPID +1; idx++) {
    pid_bin_edges[idx] = Float_t(idx) - 0.5;  // shift bins to center around int values
  }
  
  fhistogram = new TH3F("classifier_pT_PID_" + TString(classifier->GetName()),
			"classifier vs. p_{T} vs. PID",
			classifier_bins, classifier_bin_edges,
			sizeof(pt_bin_edges)/sizeof(*pt_bin_edges) - 1, pt_bin_edges,
			kNPID, pid_bin_edges);
  fhistogram->GetXaxis()->SetTitle("Classifier value");
  fhistogram->GetYaxis()->SetTitle("p_{T} (GeV)");
  fhistogram->GetZaxis()->SetTitle("PID");

  // Name bins with pdg code:
  for (Int_t ipid = 0; ipid < kNPID; ipid++) {
    fhistogram->GetZaxis()->SetBinLabel(ipid + 1, Form("%d", Pid_enum_to_pdg(ipid)));
  }
  fhistogram->Sumw2();
  fhistogram->SetDirectory(0);

  // Add this histogram to the "folder" of the classifier in which it is binned
  classifier->GetClassifierOutputList()->Add(fhistogram);
}


void AliObservableClassifierpTPID::Fill(AliMCEvent *event, AliStack *stack) {
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

    // We only want primaries and pi0's!
    if (!(stack->IsPhysicalPrimary(iTrack) ||
	  AliIsPi0PhysicalPrimary(iTrack, stack))) {
      continue;
    }

    // Ok, lets fill!
    Int_t pdgCode = track->PdgCode();
    if(((TMath::Abs(track->Y()) < 0.5))) {  //y since this is for identified particles. Region is TPC+ITS
      for (Int_t ipid = 0; ipid < kNPID; ipid++) {
	if (pdgCode == this->Pid_enum_to_pdg(ipid)){
	  fhistogram->Fill(classifier_value,
			   track->Pt(),
			   ipid,
			   event_weight);
	  break;
	}
      }
      // Fill the "allcharged" bin of the 3d histogram. Note that this is still restricted to the y<.5 region
      if (track->Charge() != 0) {
      	fhistogram->Fill(classifier_value,
			 track->Pt(),
			 this->kALLCHARGED,
			 event_weight);
      }
    }
  }
}

Int_t AliObservableClassifierpTPID::Pid_enum_to_pdg(Int_t pid_enum) {
  if (pid_enum == kPROTON) return 2212;
  else if (pid_enum == kANTIPROTON) return -2212;
  else if (pid_enum == kLAMBDA) return 3122;
  else if (pid_enum == kANTILAMBDA) return -3122;
  else if (pid_enum == kK0S) return 310;
  else if (pid_enum == kKPLUS) return 321;
  else if (pid_enum == kKMINUS) return -321;
  else if (pid_enum == kPIPLUS) return 211;
  else if (pid_enum == kPIMINUS) return -211;
  else if (pid_enum == kPI0) return 111;
  else if (pid_enum == kXI) return 3312;
  else if (pid_enum == kANTIXI) return -3312;
  else if (pid_enum == kOMEGAMINUS) return 3334;
  else if (pid_enum == kOMEGAPLUS) return -3334;
  else if (pid_enum == kLAMBDA0B) return 5122;
  else if (pid_enum == kANITLAMBDA0B) return -5122;
  else return 99999;
}

