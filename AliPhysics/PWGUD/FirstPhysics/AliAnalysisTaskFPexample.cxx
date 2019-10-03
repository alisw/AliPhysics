
#include "AliAnalysisTaskFPexample.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TList.h"
#include "TString.h"

#include "AliAnalysisManager.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliMCEvent.h"
#include "AliMultiplicity.h"
#include "AliTriggerAnalysis.h"

ClassImp(AliAnalysisTaskFPexample)

//________________________________________________________________________
AliAnalysisTaskFPexample::AliAnalysisTaskFPexample(const char *name) :
  AliAnalysisTaskFirstPhysics(name),
  fhTrackPt(0),
  fh2TrackPhiEta(0),
  fhMulITSTPC(0),
  fhMulITSSA(0),
  fhMulSPD(0),
  fh2TrackletsPhiEta(0),
  fh2TracksPhiTPCchi2(0)
{
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskFPexample::~AliAnalysisTaskFPexample()
{
  // Destructor. Clean-up the output list, but not the histograms that are put inside
  // (the list is owner and will clean-up these histograms). Protect in PROOF case.
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutput;
  }
  for (Int_t i = 0; i < knTrackCuts; i ++) {
    delete fTrackCuts[i];
    fTrackCuts[i] = 0;
  }
}

//________________________________________________________________________
void AliAnalysisTaskFPexample::UserCreateOutputObjects()
{
  PrepareOutputList();
  // Define cuts
  PrepareDefaultTrackCuts();

  // Create histograms
  const Int_t ptbins = 15;
  const Double_t ptlow = 0.1, ptup = 3.1;
  const Int_t etabins = 40;
  const Double_t etalow = -2.0, etaup = 2.0;
  const Int_t mulbins = 200;
  const Int_t phibins = 120;

  fhTrackPt = UserHisto1d("fhTrackPt", "p_{T} distribution for ESD",
			  "P_{T} (GeV/c)", ptbins, ptlow, ptup);
  fh2TrackPhiEta = UserHisto2d("fh2TrackPhiEta", "ESD tracks",
			       "#Phi", phibins, 0, 2 * TMath::Pi(),
			       "#eta", etabins, etalow, etaup);
  fhMulITSTPC = UserHisto1d("fhMulITSTPC", "N_{CH} distribution for ITS+TPC tracks",
			    "N_{CH}", mulbins, 0, mulbins);
  fhMulITSSA = UserHisto1d("fhMulITSSA", "N_{CH} distribution for ITS SA tracks",
			   "N_{CH}", mulbins, 0, mulbins);
  fhMulSPD = UserHisto1d("fhMulSPD", "N_{CH} distribution for SPD tracklets",
			 "N_{CH}", mulbins, 0, mulbins);
  fh2TrackletsPhiEta = UserHisto2d("fh2TrackletsPhiEta", "Tracklets",
				   "#Phi", phibins, 0, 2 * TMath::Pi(),
				   "#eta", etabins, etalow, etaup);
  fh2TracksPhiTPCchi2 = UserHisto2d("fh2TracksPhiTPCchi2", "ESD tracks #Phi vs. TPC #chi^{2}",
				    "#Phi", phibins, 0, 2 * TMath::Pi(),
				    "#chi^{2}", 500, 0, 500);

  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliAnalysisTaskFPexample::UserExec(Option_t *)
{
  if (!GetESDEvent()) {
    AliError("Unable to read the ESD");
    return;
  }
  if (!CheckVertex()) {
    return;
  }

  const Int_t nESDTracks = fESD->GetNumberOfTracks();

  Int_t iHighestID = 0;
  for (Int_t iTrack = 0; iTrack < nESDTracks; iTrack++) {
    if (fESD->GetTrack(iTrack)->GetLabel() > iHighestID) {
      iHighestID = fESD->GetTrack(iTrack)->GetLabel();
    }
  }
  const Int_t nMaxID = iHighestID + 1;
  bool aGlobalBits[nMaxID], aPureITSBits[nMaxID];
  for (Int_t iFlag = 0; iFlag < nMaxID; iFlag++) {
    aGlobalBits[iFlag] = false;
    aPureITSBits[iFlag] = false;
  }

  // flags for secondary and rejected tracks
  const int kRejectBit = BIT(15); // set this bit in ESD tracks if it is rejected by a cut
  const int kSecondaryBit = BIT(16); // set this bit in ESD tracks if it is secondary according to a cut

  Int_t nTracksITSTPC = 0; // multiplicity counters
  Int_t nTracksITSSA = 0;
  Int_t nTrackletsSPD = 0;

  for(Int_t iTrack = 0; iTrack < nESDTracks; iTrack++){
    AliESDtrack* esdTrack = fESD->GetTrack(iTrack);
    // if track is a secondary from a V0, flag as a secondary
    if (esdTrack->IsOn(AliESDtrack::kMultInV0)) {
      esdTrack->SetBit(kSecondaryBit);
      continue;
    }
    AliESDtrack* tpcTrack = AliESDtrackCuts::GetTPCOnlyTrack(fESD, esdTrack->GetID());
    // check tracks with ITS part
    if (esdTrack->IsOn(AliESDtrack::kITSin)) {
      if (!esdTrack->IsOn(AliESDtrack::kITSpureSA)) { // track has ITS part but is not an ITS_SA
	// TPC+ITS
	if (esdTrack->IsOn(AliESDtrack::kTPCin)) {  // Global track, has ITS and TPC contributions
	  if (fTrackCuts[kTrackCutQGlo]->AcceptTrack(esdTrack)) { // good ITS+TPC track
	    if (fTrackCuts[kTrackCutDCAwSPD]->AcceptTrack(esdTrack)
		|| fTrackCuts[kTrackCutDCAwoSPD]->AcceptTrack(esdTrack)) {
	      nTracksITSTPC++;
	      fhTrackPt->Fill(esdTrack->Pt());
	      fh2TrackPhiEta->Fill(esdTrack->Phi(), esdTrack->Eta());
	      if (tpcTrack) {
		fh2TracksPhiTPCchi2->Fill(esdTrack->Phi(), tpcTrack->GetTPCchi2());
	      }
	    } else {
	      esdTrack->SetBit(kSecondaryBit); // large DCA -> secondary, don't count either track not associated tracklet
	    }
	  } else {
	    esdTrack->SetBit(kRejectBit); // bad quality, don't count the track, but may count tracklet if associated
	  }
	} else if (fTrackCuts[kTrackCutQITS]->AcceptTrack(esdTrack)) { // good ITS complementary track
	  if (fTrackCuts[kTrackCutDCAwSPD]->AcceptTrack(esdTrack)
	      || fTrackCuts[kTrackCutDCAwoSPD]->AcceptTrack(esdTrack)) {
	    nTracksITSTPC++;
	    fhTrackPt->Fill(esdTrack->Pt());
	    fh2TrackPhiEta->Fill(esdTrack->Phi(), esdTrack->Eta());
	    if (tpcTrack) {
	      fh2TracksPhiTPCchi2->Fill(esdTrack->Phi(), tpcTrack->GetTPCchi2());
	    }
	  } else {
	    esdTrack->SetBit(kSecondaryBit); // large DCA -> secondary, don't count either track not associated tracklet
	  }
	} else {
	  esdTrack->SetBit(kRejectBit); // bad quality, don't count the track, but may count tracklet if associated
	}
      } else { // pure ITS SA tracks
	if (fTrackCuts[kTrackCutQITS]->AcceptTrack(esdTrack)) { // good ITSSA track
	  if (fTrackCuts[kTrackCutDCAwSPD]->AcceptTrack(esdTrack)
	      || fTrackCuts[kTrackCutDCAwoSPD]->AcceptTrack(esdTrack)) {
	    nTracksITSSA++;
	  } else {
	    esdTrack->SetBit(kRejectBit);
	  }
	} else {
	  esdTrack->SetBit(kRejectBit);
	}
      }
    }
    if (tpcTrack) {
      delete tpcTrack;
    }
  }

  // get multiplicity from ITS tracklets to complement TPC+ITS, and ITS pure SA
  const AliMultiplicity* multiplicitySPD = fESD->GetMultiplicity();    // spd multiplicity object
  Int_t id1, id2, id3, id4;
  AliESDtrack *tr1 = 0, *tr3 = 0;
  for (Int_t iTracklet = 0; iTracklet < multiplicitySPD->GetNumberOfTracklets(); iTracklet++) {
    if (TMath::Abs(multiplicitySPD->GetEta(iTracklet)) > GetCutEta()) {
      continue; // eta selection for tracklets
    }
    nTrackletsSPD++;
    fh2TrackletsPhiEta->Fill(multiplicitySPD->GetPhi(iTracklet), multiplicitySPD->GetEta(iTracklet));
    // if counting tracks+tracklets, check if clusters were already used in tracks
    // and get the id of the tracks in which they were used

    // references for eventual Global/ITS_SA tracks
    multiplicitySPD->GetTrackletTrackIDs(iTracklet, 0, id1, id2);
    tr1 = id1 >= 0 ? fESD->GetTrack(id1) : 0;

    // references for eventual ITS_SA_pure tracks
    multiplicitySPD->GetTrackletTrackIDs(iTracklet, 1, id3, id4);
    tr3 = id3 >= 0 ? fESD->GetTrack(id3) : 0;

    // are both clusters from the same tracks? If not, skip the
    // tracklet (shouldn't change things much)
    if (id1 != id2 || id3 != id4) {
      continue;
    }

    // has associated global track been associated to a previous tracklet?
    bool bUsedInGlobal = (id1 != -1) ? aGlobalBits[id1] : false;
    // has associated pure ITS track been associated to a previous tracklet?
    bool bUsedInPureITS = (id3 != -1) ? aPureITSBits[id3] : false;

    // counting tracklet as global+complementary track
    if ((tr1 && !tr1->TestBit(kSecondaryBit)) // reject as secondary
	&& (tr1 && tr1->TestBit(kRejectBit))) { // already accounted
      if (!bUsedInGlobal) {
	nTracksITSTPC++;
	if (id1 > 0) {
	  // mark global track linked to this tracklet as associated
	  aGlobalBits[id1] = true;
        }
      }
    } else if (id1 < 0) {
      nTracksITSTPC++;
    }

    // counting tracklet as ITS SA pure track
    if ((tr3 && tr3->TestBit(kSecondaryBit))
	&& (tr3 && !tr3->TestBit(kRejectBit))) {
      if (!bUsedInPureITS) {
	nTracksITSSA++;
	if (id3 > 0) {
	  // mark global track linked to this tracklet as associated
	  aPureITSBits[id3] = true;
	}
      }
    } else if (id3 < 0) {
      nTracksITSSA++;
    }

  }
  fhMulITSTPC->Fill(nTracksITSTPC);
  fhMulITSSA->Fill(nTracksITSSA);
  fhMulSPD->Fill(nTrackletsSPD);

  PostData(1, fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskFPexample::Terminate(Option_t *)
{
  // Draw result to screen, or perform fitting, normalizations
  // don't get too fancy, keep your histos raw and massage them with macros

  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if(!fOutput) {
    Printf("ERROR: could not retrieve TList fOutput");
    return;
  }

  if (!(GetHisto1FromOutput("fhTrackPt", fhTrackPt) &&
	GetHisto2FromOutput("fh2TrackPhiEta", fh2TrackPhiEta) &&
	GetHisto1FromOutput("fhMulITSTPC", fhMulITSTPC) &&
	GetHisto1FromOutput("fhMulITSSA", fhMulITSSA) &&
	GetHisto1FromOutput("fhMulSPD", fhMulSPD) &&
	GetHisto2FromOutput("fh2TrackletsPhiEta", fh2TrackletsPhiEta) &&
	GetHisto2FromOutput("fh2TracksPhiTPCchi2", fh2TracksPhiTPCchi2))) {
    AliError("Couldn't load every histogram from output.");
    return;
  }

  TCanvas *c = new TCanvas("AliAnalysisTaskFPexample", "Data Quality Quick Overview");
  c->Divide(2, 2);
  c->cd(1)->SetLogy();
  fhTrackPt->DrawCopy("E");
  c->cd(2);
  fh2TrackPhiEta->DrawCopy("");
  c->cd(3)->SetLogy();
  fhMulITSTPC->DrawCopy("E");
  c->cd(4)->SetLogy();
  fhMulITSSA->DrawCopy("");
}
