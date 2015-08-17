#include <iostream>

#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"

#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenDPMjetEventHeader.h"

#include "AliMultiplicityEstimator.h"

using namespace std;

ClassImp(AliMultiplicityEstimator)

AliMultiplicityEstimator::AliMultiplicityEstimator()
  : TNamed(), feta_Nch(0), fNch_pT_pid(0), fNchTuple(0),
    fuseWeights(kTRUE)
{
}

AliMultiplicityEstimator::AliMultiplicityEstimator(const char* name, const char* title,
							   Float_t feta_min_backwards, Float_t feta_max_backwards,
							   Float_t feta_min_forwards, Float_t feta_max_forwards)
  : TNamed(name, title), feta_Nch(0), fNch_pT_pid(0), fNchTuple(0),
    fuseWeights(kTRUE),
    feta_min_backwards(feta_min_backwards), feta_max_backwards(feta_max_backwards),
    feta_min_forwards(feta_min_forwards), feta_max_forwards(feta_max_forwards),
    fbypass_eta_selection(false), festimator_bins(200)
{
}

// Constructor for bypassing the eta selection to get the full range
AliMultiplicityEstimator::AliMultiplicityEstimator(const char* name, const char* title)
  : TNamed(name, title), feta_Nch(0), fNch_pT_pid(0), fNchTuple(0),
    fuseWeights(kTRUE),
    feta_min_backwards(0), feta_max_backwards(0),
    feta_min_forwards(0), feta_max_forwards(0),
    fbypass_eta_selection(true), festimator_bins(200)
{
}

Int_t AliMultiplicityEstimator::pid_enum_to_pdg(Int_t pid_enum) {
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
  else return 99999;
}

void AliMultiplicityEstimator::RegisterHistograms(TList *outputList){
  std::cout << "Registering Histogram: " << GetName() << std::endl;
  // Put all histograms of one estimator in their own sub-list
  TList *curr_est = new TList();
  curr_est->SetName(GetName());
  outputList->Add(curr_est);

  /////////////////////////////////////////////
  // Histograms to be filled in a track loop //
  /////////////////////////////////////////////
  TString postfix = fuseWeights?"":"_unweighted";
  feta_Nch = new TH2F("feta_Nch" + postfix ,
		     "N_{ch} vs. #eta, " + GetTitlePostfix(),
		     200, -10.0, 10.0,
		     festimator_bins, 0, festimator_bins);
  feta_Nch->GetXaxis()->SetTitle("#eta");
  feta_Nch->GetYaxis()->SetTitle("N_{ch} in " + GetTitlePostfix());
  feta_Nch->GetZaxis()->SetTitle("N_{ch} per #eta bin");
  feta_Nch->SetMarkerStyle(kFullCircle);
  feta_Nch->Sumw2();
  feta_Nch->SetDirectory(0);
  curr_est->Add(feta_Nch);

  // 3D: Mult, pT, PID
  Float_t mult_bin_edges[festimator_bins + 1];
  for (Int_t idx = 0; idx < festimator_bins+1; idx++){
    mult_bin_edges[idx] = idx;
  }
  Float_t pt_bin_edges[] =
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

  Float_t pid_bin_edges[kNPID + 1];
  for (Int_t idx = 0; idx < kNPID +1; idx++) {
    pid_bin_edges[idx] = Float_t(idx) - 0.5;  // shift bins to center around int values
  }

  fNch_pT_pid = new TH3F("fNch_pT_pid" + postfix,
			  Form("Event class vs. p_{T} vs. pid, %s", GetTitlePostfix().Data()),
			  festimator_bins, mult_bin_edges,
			  sizeof(pt_bin_edges)/sizeof(*pt_bin_edges) - 1, pt_bin_edges,
			  kNPID, pid_bin_edges);
  fNch_pT_pid->GetXaxis()->SetTitle("Multiplicity");
  fNch_pT_pid->GetYaxis()->SetTitle("p_{T} [GeV]");
  fNch_pT_pid->GetZaxis()->SetTitle("PID");
  // Name bins with pdg code:
  for (Int_t ipid = 0; ipid < kNPID; ipid++) {
    fNch_pT_pid->GetZaxis()->SetBinLabel(ipid + 1, Form("%d",pid_enum_to_pdg(ipid)));
  }
  fNch_pT_pid->Sumw2();
  fNch_pT_pid->SetDirectory(0);
  curr_est->Add(fNch_pT_pid);

  //////////////////////////////////////////////////////////////////////
  // Objects to be filled on the event loop (ie not in the track loop //
  //////////////////////////////////////////////////////////////////////
  fNchTuple = new TNtuple("fEventTuple", "N_{ch}^{est}", "nch");
  curr_est->Add(fNchTuple);
}

void AliMultiplicityEstimator::PreEvent(Float_t ev_weight){
  // Clear counters and chaches for the following event:
  fnch_in_estimator_region = 0;
  feventWeight = ev_weight;
}

Bool_t AliMultiplicityEstimator::TrackSelection(AliMCParticle* track) {
  if (track->Charge() != 0){
    Double_t eta = track->Eta();
    if(fbypass_eta_selection ||
       ((eta >= feta_min_backwards &&
       eta <= feta_max_backwards) ||
       (eta >= feta_min_forwards &&
	eta <= feta_max_forwards)))
      {
	return true;
      }
    }
  return false;
}

/*
  Increment counters based on whether the track meets track selection
  criteria.
*/
void AliMultiplicityEstimator::ProcessTrackForMultiplicityEstimation(AliMCParticle *track){
  if (TrackSelection(track)){
	fnch_in_estimator_region++;
  }
}

// loop over tracks again, now that mult is known:
void AliMultiplicityEstimator::ProcessTrackWithKnownMultiplicity(AliMCParticle *track){
  if (track->Charge() != 0){
    // only enforce charged tracks for dN/deta!
    feta_Nch->Fill(track->Eta(), fnch_in_estimator_region, fuseWeights?feventWeight:1);
  }
  // z axis are the different particles defined as enum in the header file.
  // ipid is not the histogram bin number! The bins of the histogram are chosen to consume the enum
  // value for the pid.
  Int_t pdgCode = track->PdgCode();
  if(((TMath::Abs(track->Y()) < 0.5))) {  //y since this is for identified particles. Region is TPC+ITS
    for (Int_t ipid = 0; ipid < kNPID; ipid++) {
      if (pdgCode == pid_enum_to_pdg(ipid)){
	fNch_pT_pid->Fill(fnch_in_estimator_region,
			   track->Pt(),
			   ipid,
			   fuseWeights?feventWeight:1);
	break;
      }
    }
  }
}

void AliMultiplicityEstimator::PostEvent(){
  // Fill event counters
  fNchTuple->Fill(fnch_in_estimator_region);
}

void AliMultiplicityEstimator::Terminate(TList* outputlist){
  // recover pointers to histograms since they are null on master
  std::cout << "Terminate " << fName << std::endl;
  TList *curr_est = static_cast<TList*>(outputlist->FindObject(GetName()));
  feta_Nch = static_cast<TH2F*>(curr_est->FindObject("feta_Nch" ));
}


