#define dStarAnalysisSelector_cxx


#include "dStarAnalysisSelector.h"
#include <TH2.h>
#include <TAxis.h>
#include <TStyle.h>

using namespace std;


dStarAnalysisSelector::dStarAnalysisSelector(TTree *) :
fChain{nullptr},
fTOFpid{nullptr},
fPool{nullptr},
fNSigmaDe{3.}, fNSigmaPi{3.},
fZVertex{nullptr},
fDedxPiDeu{ nullptr }, fDedx{ nullptr }, fBeta{nullptr},
fDeuPT{nullptr}, fPiPlusPT{nullptr}, fPiMinusPT{nullptr}, fMultiplicity{nullptr},
fMInvBlind{nullptr},
fMEMInvBackPPSame{nullptr},
fLSPlusMInvBackground{nullptr}, fLSMinusMInvBackground{nullptr},
fOutputFileName{"MInvPureEMfast.root"} {}

//______________________________________________________________________________
//______________________________________________________________________________



void dStarAnalysisSelector::Begin(TTree * /*tree*/) {
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
}

//______________________________________________________________________________
//______________________________________________________________________________


void dStarAnalysisSelector::SlaveBegin(TTree * /*tree*/) {
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  /// event mixing pool instantiation and setting
  fPool = new AliAnalysisCODEX::EventMixingPool<1, 10, 3>(100, 10., 5);
  fPool->SetPartMassI(0, AliAnalysisCODEX::kDeuteronMass);
  fPool->SetPartMassI(1, AliAnalysisCODEX::kPionMass); // piplus
  fPool->SetPartMassI(2, AliAnalysisCODEX::kPionMass); // piminus

  /// output histogram instantiation
  fZVertex      = new TH1F("z_vertex", "", 22, -11, 11);
  fPiPlusPT     = new TH1F("pion_plus_pT", "", 200, 0., 10.);
  fPiMinusPT    = new TH1F("pion_minus_pT", "", 200, 0., 10.);
  fDeuPT        = new TH1F("deuteron_pT", "", 200, 0., 10.);
  fMultiplicity = new TH1F("event_multiplicity", "", 15, 0, 15);

  /// energy loss in TPC for Pions
  fDedxPiDeu  = new TH2D("fDedxPiDeu", "dE/dX vs p; p (GeV/c); TPC dE/dx (a.u.)", 1000, 0.1, 15., 1000, 0., 1000.);
  fDedx       = new TH2D("fDedx", "dE/dX vs p; p (GeV/c); TPC dE/dx (a.u.)", 1000, 0.1, 15., 1000, 0., 1000.);
  fBeta       = new TH2D("fBeta", "Beta vs p; p (GeV/c); Beta ", 10000, 0.2, 15., 1000, 0., 1.2);

  /// invariant mass distribution with blinded region for different cuts on π+ π- invariant mass.
  // for (int i = 0; i < 10; i++) {
  //   fMInvBlind[i]    = new TH2D(Form("minv_blind_%.2fcut", 0.4 + (0.02*i)), "", 840, 2.120, 8.000, 10, 0, 10);
  // }
  /// the last one with no cuts
  fMInvBlind[10]    = new TH2D("minv_blind_nocut", "", 840, 2.120, 8.000, 10, 0, 10);

  /// event mixing background with cuts on pi pi invariant mass
  // for (int i = 0; i < 10; i++) {
  //   fMEMInvBackPPSame[i] = new TH2D(Form("background_em_%.2fcut", 0.4 + (0.02*i)), "", 840, 2.120, 8.000, 10, 0, 10);
  // }
  /// last oune with no cuts
  fMEMInvBackPPSame[10] = new TH2D("em_nocut_ppsame", "", 840, 2.120, 8.000, 10, 0, 10);
  fMEMInvBack3Event[10] = new TH2D("em_nocut_3event", "", 840, 2.120, 8.000, 10, 0, 10);

  /// like-sign background
  // for (int i = 0; i < 10; i++) {
  //   fLSPlusMInvBackground[i]  = new TH2D(Form("background_lsPlus_%.2fcut", 0.4 + (0.02*i)), "", 840, 2.120, 8.000, 10, 0, 10);
  //   fLSMinusMInvBackground[i] = new TH2D(Form("background_lsMinus_%.2fcut", 0.4 + (0.02*i)), "", 840, 2.120, 8.000, 10, 0, 10);
  // }
  /// last one with no cuts
  fLSPlusMInvBackground[10]  = new TH2D("background_lsPlus_nocut", "", 840, 2.120, 8.000, 10, 0, 10);
  fLSMinusMInvBackground[10] = new TH2D("background_lsMinus_nocut", "", 840, 2.120, 8.000, 10, 0, 10);

  GetOutputList()->Add(fZVertex);
  GetOutputList()->Add(fPiPlusPT);
  GetOutputList()->Add(fPiMinusPT);
  GetOutputList()->Add(fDeuPT);
  GetOutputList()->Add(fMultiplicity);
  GetOutputList()->Add(fDedxPiDeu);
  GetOutputList()->Add(fDedx);
  GetOutputList()->Add(fBeta);
  for (int i = 0; i < 10; i++) {
    GetOutputList()->Add(fMInvBlind[i]);
    GetOutputList()->Add(fMEMInvBackPPSame[i]);
    GetOutputList()->Add(fLSPlusMInvBackground[i]);
    GetOutputList()->Add(fLSMinusMInvBackground[i]);
  }
  GetOutputList()->Add(fMInvBlind[10]);
  GetOutputList()->Add(fMEMInvBackPPSame[10]);
  GetOutputList()->Add(fMEMInvBack3Event[10]);
  // GetOutputList()->Add(fLSPlusMInvBackground[10]);
  // GetOutputList()->Add(fLSMinusMInvBackground[10]);

  BinLogAxis(fDedxPiDeu);
  BinLogAxis(fDedx);

}

//______________________________________________________________________________
//______________________________________________________________________________


Bool_t dStarAnalysisSelector::Process(Long64_t entry) {
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  bool eventmixing3 = true;
  bool eventmixing = false;
  bool likesign = false;

  fReader.SetEntry(entry);


  fDeuteron.clear();
  fPiPlus.clear();
  fPiMinus.clear();

  int multiplicity = std::distance(fTracks.begin(), fTracks.end())/10 < 10 ? std::distance(fTracks.begin(), fTracks.end())/10 : 10;
  fMultiplicity->Fill(multiplicity);

  Short_t zvtx = ((short)*fZvert.Get())/1000.;


  /// PID cuts on deuterons and pions (definition below)
  /// and conversion from AliAnalysisCODEX::Track to FourVector_t
  for (const auto track : fTracks) {
    // fill PID histograms
    // fDedx->Fill(track.GetTPCmomentum(), track.GetTPCsignal());
    // if (track.HasTOF()) fBeta->Fill(track.GetTPCmomentum(), track.TOFbeta());

    // PID checks
    if (DeuteronCutGood(track, fNSigmaDe)) {
      fDeuteron.push_back(TrackToFourVector(track, 4));
      // fDedxPiDeu->Fill(track.GetTPCsignal(), track.GetTPCsignal());
      // fZVertex->Fill(zvtx);
    }
    if (PionCutGood(track, fNSigmaPi)) {
      track.Charge() > 0 ? fPiPlus.push_back(TrackToFourVector(track, 1)) : fPiMinus.push_back(TrackToFourVector(track, 1));
      // fDedxPiDeu->Fill(track.GetTPCmomentum(), track.GetTPCsignal());
    }
  }

  /// deuterons pT for check
  for (const auto &deu : fDeuteron) {
    fDeuPT->Fill(deu.Pt());
  }
  /// pions pT for check
  for (const auto &pip : fPiPlus) {
    fPiPlusPT->Fill(pip.Pt());
  }
  /// like above
  for (const auto &pim : fPiMinus) {
    fPiMinusPT->Fill(pim.Pt());
  }

  //______________________________________________________________________________

  /// invarian mass distribution for π+ π- d  with blinded region
  for (const auto &deu : fDeuteron) {
    for (const auto &pip : fPiPlus) {
      for (const auto &pim : fPiMinus) {
        FourVector_t vpp = pip + pim;
        FourVector_t v = vpp + deu;
        /// cuts on π+ π- invariant mass
        // for (int i = 0; i < 10; i++) {
        //   if (vpp.M() < (0.4+(0.02*i))) {
        //     if (v.M() < 2.480 && v.M() > 2.280) continue;
        //     fMInvBlind[i]->Fill(v.M(), v.Pt());
        //   }
        // }
        if (v.M() < 2.477 && v.M() > 2.281) continue;
        fMInvBlind[10]->Fill(v.M(), v.Pt());
      }
    }
  }

  //______________________________________________________________________________

  if (eventmixing3) {
    /// Event Mixing (da valutare le varie possibilita di mixing)
    if (!fDeuteron.empty()) {
      for (int i = 0; i < fPool->GetNLevelOccupied(multiplicity, zvtx, 0); i++ ) {
        auto deuvec = fPool->GetVectorFV(multiplicity, zvtx, 0, i);
        for (const auto &deu : deuvec) {
          for (int j = 0; j < fPool->GetNLevelOccupied(multiplicity, zvtx, 1); j++) {
            auto pipvec = fPool->GetVectorFV(multiplicity, zvtx, 1, j);
            for (const auto &pip : pipvec) {
              for (const auto &pim : fPiMinus) {
                FourVector_t vpp = pip + pim;
                FourVector_t v = vpp + deu;
                fMEMInvBack3Event[10]->Fill(v.M(), v.Pt());
              }
            }
          }
        }
      }
      fPool->FillEvent(fDeuteron, multiplicity, zvtx, 0);
    }

    for (int i = 0; i < fPool->GetNLevelOccupied(multiplicity, zvtx, 1); i++) {
      auto pipvec = fPool->GetVectorFV(multiplicity, zvtx, 1, i);
      for (const auto &pip : pipvec) {
        for (int j = 0; j < fPool->GetNLevelOccupied(multiplicity, zvtx, 2); j++) {
          auto pimvec = fPool->GetVectorFV(multiplicity, zvtx, 2, j);
          for (const auto &pim : pimvec) {
            for (const auto &deu : fDeuteron) {
              FourVector_t vpp = pip + pim;
              FourVector_t v = vpp + deu;
              fMEMInvBack3Event[10]->Fill(v.M(), v.Pt());
            }
          }
        }
      }
    }

    for (int i = 0; i < fPool->GetNLevelOccupied(multiplicity, zvtx, 0); i++) {
      auto deuvec = fPool->GetVectorFV(multiplicity, zvtx, 0, i);
      for (const auto &deu : deuvec) {
        for (int j = 0; j < fPool->GetNLevelOccupied(multiplicity, zvtx, 2); j++) {
          auto pimvec = fPool->GetVectorFV(multiplicity, zvtx, 2, j);
          for (const auto &pim : pimvec) {
            for (const auto &pip : fPiPlus) {
              FourVector_t vpp = pip + pim;
              FourVector_t v = vpp + deu;
              fMEMInvBack3Event[10]->Fill(v.M(), v.Pt());
            }
          }
        }
      }
    }

    fPool->FillEvent(fPiPlus, multiplicity, zvtx, 1);
    fPool->FillEvent(fPiMinus, multiplicity, zvtx, 2);
  }

  //______________________________________________________________________________

  if (eventmixing) {
    /// Event Mixing (da valutare le varie possibilita di mixing)
    if (!fDeuteron.empty()) {
      for (int i = 0; i < fPool->GetNLevelOccupied('a', zvtx, 0); i++ ) {
        auto deu = fPool->GetVectorFV('a', zvtx, 0, i);
        for (const auto &deuvec : deu) {
          for (const auto &pip : fPiPlus) {
            for (const auto &pim : fPiMinus) {
              FourVector_t vpp = pip + pim;
              FourVector_t v = vpp + deuvec;
              // for (int j = 0; j < 10; j++) {
              //   if (vpp.M() < (0.4+(0.02*j))) {
              //     fMEMInvBackPPSame[j]->Fill(v.M(), v.Pt());
              //   }
              // }
              fMEMInvBackPPSame[10]->Fill(v.M(), v.Pt());
            }
          }
        }
      }
      fPool->FillEvent(fDeuteron, 'a', zvtx, 0);
    }
  }


  if (likesign) {
    /// like-sign background model
    if (!fDeuteron.empty()) {
      /// like sign with pi plus
      for (const auto &deu : fDeuteron) {
        for (const auto &pip1 : fPiPlus) {
          for (const auto &pip2 : fPiPlus) {
            FourVector_t vpp = pip1 + pip2;
            FourVector_t v = vpp + deu;
            // for (int j = 0; j < 10; j++) {
            //   if (vpp.M() < (0.4+(0.02*j))) {
            //     fLSPlusMInvBackground[j]->Fill(v.M(), v.Pt());
            //   }
            // }
            fLSPlusMInvBackground[10]->Fill(v.M(), v.Pt());
          }
        }
      }
      /// like sign with pi minus
      for (const auto &deu : fDeuteron) {
        for (const auto &pim1 : fPiMinus) {
          for (const auto &pim2 : fPiMinus) {
            FourVector_t vpp = pim1 + pim2;
            FourVector_t v = vpp + deu;
            for (int j = 0; j < 10; j++) {
              //   if (vpp.M() < (0.4+(0.02*j))) {
              //     fLSMinusMInvBackground[j]->Fill(v.M(), v.Pt());
              //   }
              // }
              fLSMinusMInvBackground[10]->Fill(v.M(), v.Pt());
            }
          }
        }
      }
    }
  }

  return kTRUE;
}

//______________________________________________________________________________
//______________________________________________________________________________


void dStarAnalysisSelector::SlaveTerminate() {
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

}

//______________________________________________________________________________
//______________________________________________________________________________


void dStarAnalysisSelector::Terminate() {
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  TFile output(Form("~/results/%s", fOutputFileName.data()),"RECREATE");
  GetOutputList()->Write();
  output.Close();

}

//______________________________________________________________________________
//_______________________________________________________________________________


void dStarAnalysisSelector::Init(TTree *tree) {
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the reader is initialized.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  fReader.SetTree(tree);

  /// set TOFpidLite parmas based on TOF config in analyzed runs
  float TOFparams[4] = {0.004, 0., 0.009, 15.};
  fTOFpid.SetParams(TOFparams);
  fTOFpid.SetTOFresolution(80.);
  fTOFpid.SetTOFtail(0.95);

  /// setting header branch
  header = new Header();
  header_b = tree->GetBranch("Header");
  header_b->SetAddress(&header);

  /// connect header to TOFpidLite
  fTOFpid.ConnectHeader(header);

}

//______________________________________________________________________________
//_______________________________________________________________________________


Bool_t dStarAnalysisSelector::Notify() {
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

//______________________________________________________________________________
//______________________________________________________________________________


bool dStarAnalysisSelector::PionCutGood(const Track& tr, float nSigmaTPC) {
  //
  // Method that provides true if the input track passes all standard cuts for Pions
  //
  if ( tr.GetNumberOfSigmasTPC(1) < nSigmaTPC && tr.GetNumberOfSigmasTPC(1) > -nSigmaTPC &&
  tr.GetTPCNcls() > 60 &&
  tr.GetTPCChi2NDF() < 4 &&
  !tr.Is(AliAnalysisCODEX::kIsKink) && tr.TPCrefit() && tr.ITSrefit() &&
  (tr.HasPointOnITSLayer(0) || tr.HasPointOnITSLayer(1)) &&
  fabs(tr.GetDCAxy()) < (0.0105 + 0.0350 / TMath::Power(tr.Pt(), 1.1)) &&
  fabs(tr.GetDCAz()) < 2 && (tr.Pt() > 0.001 && tr.Pt() < 10.0) &&
  (TMath::Abs(tr.eta) < 0.8)) {
    return true;
  } else {
    return false;
  }

}

//_______________________________________________________________________________
//_______________________________________________________________________________


bool dStarAnalysisSelector::DeuteronCutGood(const Track& tr, float nSigmaTPC) {
  //
  // Method that provides true if the input track passes all cuts for Deuterons
  //
  if ( tr.GetNumberOfSigmasTPC(4) < nSigmaTPC && tr.GetNumberOfSigmasTPC(4) > -nSigmaTPC &&
  tr.GetTPCNcls() > 70 &&
  tr.GetTPCChi2NDF() < 4 &&
  !tr.Is(AliAnalysisCODEX::kIsKink) && tr.TPCrefit() && tr.ITSrefit() &&
  (tr.HasPointOnITSLayer(0) || tr.HasPointOnITSLayer(1)) &&
  fabs(tr.GetDCAxy()) < (0.0105 + 0.0350 / TMath::Power(tr.Pt(), 1.1)) &&
  fabs(tr.GetDCAz()) < 2 && (tr.Pt() > 0.001 && tr.Pt() < 10.0) &&
  (TMath::Abs(tr.eta) < 0.9)) {
    if (tr.P() > 1.5) {
      return (tr.HasTOF() && fTOFpid.GetNumberOfSigmas(tr, AliAnalysisCODEX::kDeuteronMass) < 3. && fTOFpid.GetNumberOfSigmas(tr, AliAnalysisCODEX::kDeuteronMass) > -3.) ? true : false;
    } else {
      return true;
    }
  } else {
    return false;
  }

}

//_______________________________________________________________________________
//_______________________________________________________________________________


FourVector_t dStarAnalysisSelector::TrackToFourVector(const Track &tr, int species) {
  //
  // Method that convert AliAnalysisCODEX::Track in FourVector_t for using much less memory
  //
  FourVector_t v = {(float)tr.Pt(), (float)tr.Eta(), (float)tr.Phi(), (float)AliAnalysisCODEX::kMasses[species]};
  return v;
}

//_______________________________________________________________________________
//_______________________________________________________________________________

void dStarAnalysisSelector::BinLogAxis(TH2 *h) {
  //
  // Method for the correct logarithmic binning of histograms
  //
  TAxis* axis = h->GetXaxis();
  int bins = axis->GetNbins();
  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins + 1];

  newBins[0] = from;
  Double_t factor = pow(to / from, 1. / bins);

  for (int i = 1; i <= bins; i++) {
    newBins[i] = factor * newBins[i - 1];
  }
  axis->Set(bins, newBins);
  delete[] newBins;
}
