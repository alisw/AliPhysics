////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnDEtaDPhi - A correlation function that analyzes            //
// two particle correlations with respect to the azimuthal angle (phi)        //
// and pseudorapidity (eta) difference                                        //
//                                                                            //
// Authors: Adam Kisiel Adam.Kisiel@cern.ch                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoCorrFctnDEtaDPhiSimple.h"
#include "AliFemtoModelHiddenInfo.h"
//#include "AliFemtoHisto.hh"
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
ClassImp(AliFemtoCorrFctnDEtaDPhiSimple)
#endif

#define PIH 1.57079632679489656
#define PIT 6.28318530717958623

AliFemtoCorrFctnDEtaDPhiSimple
  ::AliFemtoCorrFctnDEtaDPhiSimple(const char *title,
                                   const int aPhiBins,
                                   const int aEtaBins)
  : AliFemtoCorrFctn()
  , fDPhiDEtaNumerator(nullptr)
  , fDPhiDEtaDenominator(nullptr)
  , fDPhiDEtaHiddenNumerator(nullptr)
  , fDPhiDEtaHiddenDenominator(nullptr)
  , fDPhiDEtaHiddenPrimaryNumerator(nullptr)
  , fDPhiDEtaHiddenPrimaryDenominator(nullptr)
  , fDPhiDEtaHiddenSecWeakNumerator(nullptr)
  , fDPhiDEtaHiddenSecWeakDenominator(nullptr)
  , fDPhiDEtaHiddenSecMatNumerator(nullptr)
  , fDPhiDEtaHiddenSecMatDenominator(nullptr)
  , fDPhiDEtaHiddenPrimaryNumeratorData(nullptr)
  , fDPhiDEtaHiddenPrimaryDenominatorData(nullptr)
  , fDPhiDEtaHiddenSecWeakNumeratorData(nullptr)
  , fDPhiDEtaHiddenSecWeakDenominatorData(nullptr)
  , fDPhiDEtaHiddenSecMatNumeratorData(nullptr)
  , fDPhiDEtaHiddenSecMatDenominatorData(nullptr)
  , fphiL(0.0)
  , fphiT(0.0)
  , fEtaBins(aEtaBins)
  , fPhiBins(aPhiBins)
  , fTitle(title)
  , fReadHiddenInfo(false)
{
  fphiL = (-(int)(aPhiBins / 4) + 0.5) * 2. * TMath::Pi() / aPhiBins;
  fphiT = 2 * TMath::Pi() + (-(int)(aPhiBins / 4) + 0.5) * 2. * TMath::Pi() / aPhiBins;

  fDPhiDEtaNumerator = new TH2D("NumDPhiDEta" + fTitle,
                                "Numerator; #Delta#phi; #Delta#eta;",
                                aPhiBins, fphiL, fphiT,
                                aEtaBins, -2.0, 2.0);

  fDPhiDEtaDenominator = new TH2D("DenDPhiDEta" + fTitle,
                                  "Denominator; #Delta #phi; #Delta#eta",
                                  aPhiBins, fphiL, fphiT,
                                  aEtaBins, -2.0, 2.0);

  fDPhiDEtaNumerator->Sumw2();
}

//____________________________
AliFemtoCorrFctnDEtaDPhiSimple
  ::AliFemtoCorrFctnDEtaDPhiSimple(const AliFemtoCorrFctnDEtaDPhiSimple &aCorrFctn)
  : AliFemtoCorrFctn(aCorrFctn)
  , fDPhiDEtaNumerator(new TH2D(*aCorrFctn.fDPhiDEtaNumerator))
  , fDPhiDEtaDenominator(new TH2D(*aCorrFctn.fDPhiDEtaDenominator))
  , fDPhiDEtaHiddenNumerator(nullptr)
  , fDPhiDEtaHiddenDenominator(nullptr)
  , fDPhiDEtaHiddenPrimaryNumerator(nullptr)
  , fDPhiDEtaHiddenPrimaryDenominator(nullptr)
  , fDPhiDEtaHiddenSecWeakNumerator(nullptr)
  , fDPhiDEtaHiddenSecWeakDenominator(nullptr)
  , fDPhiDEtaHiddenSecMatNumerator(nullptr)
  , fDPhiDEtaHiddenSecMatDenominator(nullptr)
  , fDPhiDEtaHiddenPrimaryNumeratorData(nullptr)
  , fDPhiDEtaHiddenPrimaryDenominatorData(nullptr)
  , fDPhiDEtaHiddenSecWeakNumeratorData(nullptr)
  , fDPhiDEtaHiddenSecWeakDenominatorData(nullptr)
  , fDPhiDEtaHiddenSecMatNumeratorData(nullptr)
  , fDPhiDEtaHiddenSecMatDenominatorData(nullptr)
  , fphiL(aCorrFctn.fphiL)
  , fphiT(aCorrFctn.fphiT)
  , fEtaBins(aCorrFctn.fEtaBins)
  , fPhiBins(aCorrFctn.fPhiBins)
  , fTitle(aCorrFctn.fTitle)
  , fReadHiddenInfo(aCorrFctn.fReadHiddenInfo)
{
  if (fReadHiddenInfo) {
    fDPhiDEtaHiddenNumerator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenNumerator);
    fDPhiDEtaHiddenDenominator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenDenominator);

    fDPhiDEtaHiddenPrimaryNumerator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenPrimaryNumerator);
    fDPhiDEtaHiddenPrimaryDenominator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenPrimaryDenominator);

    fDPhiDEtaHiddenSecWeakNumerator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecWeakNumerator);
    fDPhiDEtaHiddenSecWeakDenominator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecWeakDenominator);

    fDPhiDEtaHiddenSecMatNumerator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecMatNumerator);
    fDPhiDEtaHiddenSecMatDenominator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecMatDenominator);

    fDPhiDEtaHiddenPrimaryNumeratorData = new TH2D(*aCorrFctn.fDPhiDEtaHiddenPrimaryNumeratorData);
    fDPhiDEtaHiddenPrimaryDenominatorData = new TH2D(*aCorrFctn.fDPhiDEtaHiddenPrimaryDenominatorData);

    fDPhiDEtaHiddenSecWeakNumeratorData = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecWeakNumeratorData);
    fDPhiDEtaHiddenSecWeakDenominatorData = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecWeakDenominatorData);

    fDPhiDEtaHiddenSecMatNumeratorData = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecMatNumeratorData);
    fDPhiDEtaHiddenSecMatDenominatorData = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecMatDenominatorData);
  }
}

AliFemtoCorrFctnDEtaDPhiSimple::~AliFemtoCorrFctnDEtaDPhiSimple()
{
  // destructor

  delete fDPhiDEtaNumerator;
  delete fDPhiDEtaDenominator;

  if (fReadHiddenInfo) {
    delete fDPhiDEtaHiddenNumerator;
    delete fDPhiDEtaHiddenDenominator;
    delete fDPhiDEtaHiddenPrimaryNumerator;
    delete fDPhiDEtaHiddenPrimaryDenominator;
    delete fDPhiDEtaHiddenSecWeakNumerator;
    delete fDPhiDEtaHiddenSecWeakDenominator;
    delete fDPhiDEtaHiddenSecMatNumerator;
    delete fDPhiDEtaHiddenSecMatDenominator;
    delete fDPhiDEtaHiddenPrimaryNumeratorData;
    delete fDPhiDEtaHiddenPrimaryDenominatorData;
    delete fDPhiDEtaHiddenSecWeakNumeratorData;
    delete fDPhiDEtaHiddenSecWeakDenominatorData;
    delete fDPhiDEtaHiddenSecMatNumeratorData;
    delete fDPhiDEtaHiddenSecMatDenominatorData;
  }
}

AliFemtoCorrFctnDEtaDPhiSimple&
AliFemtoCorrFctnDEtaDPhiSimple::operator=(const AliFemtoCorrFctnDEtaDPhiSimple &aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn) {
    return *this;
  }

  AliFemtoCorrFctn::operator=(aCorrFctn);

  fEtaBins = aCorrFctn.fEtaBins;
  fPhiBins = aCorrFctn.fPhiBins;

  fTitle = aCorrFctn.fTitle;

  *fDPhiDEtaNumerator = *aCorrFctn.fDPhiDEtaNumerator;
  *fDPhiDEtaDenominator = *aCorrFctn.fDPhiDEtaDenominator;

  if (aCorrFctn.fReadHiddenInfo) {
    if (fReadHiddenInfo) {
      *fDPhiDEtaHiddenNumerator = *aCorrFctn.fDPhiDEtaHiddenNumerator;
      *fDPhiDEtaHiddenDenominator = *aCorrFctn.fDPhiDEtaHiddenDenominator;
      *fDPhiDEtaHiddenPrimaryNumerator = *aCorrFctn.fDPhiDEtaHiddenPrimaryNumerator;
      *fDPhiDEtaHiddenPrimaryDenominator = *aCorrFctn.fDPhiDEtaHiddenPrimaryDenominator;
      *fDPhiDEtaHiddenSecWeakNumerator = *aCorrFctn.fDPhiDEtaHiddenSecWeakNumerator;
      *fDPhiDEtaHiddenSecWeakDenominator = *aCorrFctn.fDPhiDEtaHiddenSecWeakDenominator;
      *fDPhiDEtaHiddenSecMatNumerator = *aCorrFctn.fDPhiDEtaHiddenSecMatNumerator;
      *fDPhiDEtaHiddenSecMatDenominator = *aCorrFctn.fDPhiDEtaHiddenSecMatDenominator;
      *fDPhiDEtaHiddenPrimaryNumeratorData = *aCorrFctn.fDPhiDEtaHiddenPrimaryNumeratorData;
      *fDPhiDEtaHiddenPrimaryDenominatorData = *aCorrFctn.fDPhiDEtaHiddenPrimaryDenominatorData;
      *fDPhiDEtaHiddenSecWeakNumeratorData = *aCorrFctn.fDPhiDEtaHiddenSecWeakNumeratorData;
      *fDPhiDEtaHiddenSecWeakDenominatorData = *aCorrFctn.fDPhiDEtaHiddenSecWeakDenominatorData;
      *fDPhiDEtaHiddenSecMatNumeratorData = *aCorrFctn.fDPhiDEtaHiddenSecMatNumeratorData;
      *fDPhiDEtaHiddenSecMatDenominatorData = *aCorrFctn.fDPhiDEtaHiddenSecMatDenominatorData;
    } else {
      fDPhiDEtaHiddenNumerator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenNumerator);
      fDPhiDEtaHiddenDenominator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenDenominator);
      fDPhiDEtaHiddenPrimaryNumerator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenPrimaryNumerator);
      fDPhiDEtaHiddenPrimaryDenominator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenPrimaryDenominator);
      fDPhiDEtaHiddenSecWeakNumerator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecWeakNumerator);
      fDPhiDEtaHiddenSecWeakDenominator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecWeakDenominator);
      fDPhiDEtaHiddenSecMatNumerator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecMatNumerator);
      fDPhiDEtaHiddenSecMatDenominator = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecMatDenominator);
      fDPhiDEtaHiddenPrimaryNumeratorData = new TH2D(*aCorrFctn.fDPhiDEtaHiddenPrimaryNumeratorData);
      fDPhiDEtaHiddenPrimaryDenominatorData = new TH2D(*aCorrFctn.fDPhiDEtaHiddenPrimaryDenominatorData);
      fDPhiDEtaHiddenSecWeakNumeratorData = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecWeakNumeratorData);
      fDPhiDEtaHiddenSecWeakDenominatorData = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecWeakDenominatorData);
      fDPhiDEtaHiddenSecMatNumeratorData = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecMatNumeratorData);
      fDPhiDEtaHiddenSecMatDenominatorData = new TH2D(*aCorrFctn.fDPhiDEtaHiddenSecMatDenominatorData);
    }
  } else if (fReadHiddenInfo) {
    delete fDPhiDEtaHiddenNumerator;
    delete fDPhiDEtaHiddenDenominator;
    delete fDPhiDEtaHiddenPrimaryNumerator;
    delete fDPhiDEtaHiddenPrimaryDenominator;
    delete fDPhiDEtaHiddenSecWeakNumerator;
    delete fDPhiDEtaHiddenSecWeakDenominator;
    delete fDPhiDEtaHiddenSecMatNumerator;
    delete fDPhiDEtaHiddenSecMatDenominator;
    delete fDPhiDEtaHiddenPrimaryNumeratorData;
    delete fDPhiDEtaHiddenPrimaryDenominatorData;
    delete fDPhiDEtaHiddenSecWeakNumeratorData;
    delete fDPhiDEtaHiddenSecWeakDenominatorData;
    delete fDPhiDEtaHiddenSecMatNumeratorData;
    delete fDPhiDEtaHiddenSecMatDenominatorData;

    fDPhiDEtaHiddenNumerator = nullptr;
    fDPhiDEtaHiddenDenominator = nullptr;
    fDPhiDEtaHiddenPrimaryNumerator = nullptr;
    fDPhiDEtaHiddenPrimaryDenominator = nullptr;
    fDPhiDEtaHiddenSecWeakNumerator = nullptr;
    fDPhiDEtaHiddenSecWeakDenominator = nullptr;
    fDPhiDEtaHiddenSecMatNumerator = nullptr;
    fDPhiDEtaHiddenSecMatDenominator = nullptr;
    fDPhiDEtaHiddenPrimaryNumeratorData = nullptr;
    fDPhiDEtaHiddenPrimaryDenominatorData = nullptr;
    fDPhiDEtaHiddenSecWeakNumeratorData = nullptr;
    fDPhiDEtaHiddenSecWeakDenominatorData = nullptr;
    fDPhiDEtaHiddenSecMatNumeratorData = nullptr;
    fDPhiDEtaHiddenSecMatDenominatorData = nullptr;
  }

  fReadHiddenInfo = aCorrFctn.fReadHiddenInfo;

  return *this;
}
//_________________________
void AliFemtoCorrFctnDEtaDPhiSimple::Finish()
{
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  //mShareDenominator->Draw();
  //mRatio->Draw();
}

//____________________________
AliFemtoString AliFemtoCorrFctnDEtaDPhiSimple::Report()
{
  // create report
  AliFemtoString report = "TPC Ncls Correlation Function Report:\n";
  report += Form("Number of entries in numerator:\t%E\n", fDPhiDEtaNumerator->GetEntries());
  report += Form("Number of entries in denominator:\t%E\n", fDPhiDEtaDenominator->GetEntries());

  if (fReadHiddenInfo)
  {
    report += Form("Number of entries in hidden numerator:\t%E\n", fDPhiDEtaHiddenNumerator->GetEntries());
    report += Form("Number of entries in hidden denominator:\t%E\n", fDPhiDEtaHiddenDenominator->GetEntries());

    report += Form("Number of entries in hidden primary numerator:\t%E\n", fDPhiDEtaHiddenPrimaryNumerator->GetEntries());
    report += Form("Number of entries in hidden primary denominator:\t%E\n", fDPhiDEtaHiddenPrimaryDenominator->GetEntries());

    report += Form("Number of entries in hidden second. weak numerator:\t%E\n", fDPhiDEtaHiddenSecWeakNumerator->GetEntries());
    report += Form("Number of entries in hidden second. weak denominator:\t%E\n", fDPhiDEtaHiddenSecWeakDenominator->GetEntries());

    report += Form("Number of entries in hidden second. material numerator:\t%E\n", fDPhiDEtaHiddenSecMatNumerator->GetEntries());
    report += Form("Number of entries in hidden second. material denominator:\t%E\n", fDPhiDEtaHiddenSecMatDenominator->GetEntries());

    report += Form("Number of entries in hidden primary numerator data:\t%E\n", fDPhiDEtaHiddenPrimaryNumeratorData->GetEntries());
    report += Form("Number of entries in hidden primary denominator data:\t%E\n", fDPhiDEtaHiddenPrimaryDenominatorData->GetEntries());

    report += Form("Number of entries in hidden second. weak numerator data:\t%E\n", fDPhiDEtaHiddenSecWeakNumeratorData->GetEntries());
    report += Form("Number of entries in hidden second. weak denominator data:\t%E\n", fDPhiDEtaHiddenSecWeakDenominatorData->GetEntries());

    report += Form("Number of entries in hidden second. material numerator data:\t%E\n", fDPhiDEtaHiddenSecMatNumeratorData->GetEntries());
    report += Form("Number of entries in hidden second. material denominator data:\t%E\n", fDPhiDEtaHiddenSecMatDenominatorData->GetEntries());
  }

  return report;
}

inline
void AddPair(const AliFemtoPair &pair,
             bool read_hidden_info,
             TH2D &data_hist,
             TH2D *mc_hist,
             TH2D *primary_hist,
             TH2D *primary_data_hist,
             TH2D *secweak_hist,
             TH2D *secweak_data_hist,
             TH2D *secmat_hist,
             TH2D *secmat_data_hist)
{
  const auto &track1 = *pair.Track1(),
             &track2 = *pair.Track2();

  const double PHI_LO = data_hist.GetXaxis()->GetXmin(),
               PHI_HI = data_hist.GetXaxis()->GetXmax();

  // Calculates the difference in phi and shifts to within the histogram axis
  auto calc_delta_phi = [&] (const AliFemtoThreeVector &p1, const AliFemtoThreeVector &p2) {
    double delta_phi = p1.Phi() - p2.Phi();
    while (delta_phi < PHI_LO) {
      delta_phi += PIT;
    }
    while (delta_phi > PHI_HI) {
      delta_phi -= PIT;
    }
    return delta_phi;
  };


  const auto &P1 = track1.FourMomentum(),
             &P2 = track2.FourMomentum();

  const double dphi = calc_delta_phi(P1.vect(), P2.vect()),
               deta = P1.PseudoRapidity() - P2.PseudoRapidity();

  data_hist.Fill(dphi, deta);

  if (read_hidden_info) {

    // tries to get the
    auto get_track_info = [] (const AliFemtoParticle &particle) {
      AliFemtoModelHiddenInfo *result = nullptr;

      if (auto *track = particle.Track()) {
        result = static_cast<AliFemtoModelHiddenInfo *>(track->GetHiddenInfo());
      } else if (auto *v0 = particle.V0()) {
        result = static_cast<AliFemtoModelHiddenInfo*>(v0->GetHiddenInfo());
      }

      return result;
    };

    AliFemtoModelHiddenInfo const *hInfo1 = get_track_info(track1),
                                  *hInfo2 = get_track_info(track2);

    if (hInfo1 && hInfo2) {

      const AliFemtoThreeVector
        &v1 = *hInfo1->GetTrueMomentum(),
        &v2 = *hInfo2->GetTrueMomentum();

      const double
        dhphi = calc_delta_phi(v1, v2),
        dheta = v1.PseudoRapidity() - v2.PseudoRapidity();

      const int
        origin_1 = hInfo1->GetOrigin(),
        origin_2 = hInfo2->GetOrigin();

      mc_hist->Fill(dhphi, dheta);
      if (origin_1 == 0 && origin_2 == 0) {
        primary_hist->Fill(dhphi, dheta);
        primary_data_hist->Fill(dphi, deta);
      }
      else if (origin_1 == 1 || origin_2 == 1) {
        secweak_hist->Fill(dhphi, dheta);
        secweak_data_hist->Fill(dphi, deta);
      }
      else if (origin_1 == 2 || origin_2 == 2) {
        secmat_hist->Fill(dhphi, dheta);
        secmat_data_hist->Fill(dphi, deta);
      }
    }
  }
}

void AliFemtoCorrFctnDEtaDPhiSimple::AddRealPair(AliFemtoPair *pair)
{
  // add real (effect) pair
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  // call with numerator targets
  AddPair(*pair,
          fReadHiddenInfo,
          *fDPhiDEtaNumerator,
          fDPhiDEtaHiddenNumerator,
          fDPhiDEtaHiddenPrimaryNumerator,
          fDPhiDEtaHiddenPrimaryNumeratorData,
          fDPhiDEtaHiddenSecWeakNumerator,
          fDPhiDEtaHiddenSecWeakNumeratorData,
          fDPhiDEtaHiddenSecMatNumerator,
          fDPhiDEtaHiddenSecMatNumeratorData);
}

void AliFemtoCorrFctnDEtaDPhiSimple::AddMixedPair(AliFemtoPair *pair)
{
  // add mixed (background) pair
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  // call with denominator targets
  AddPair(*pair,
          fReadHiddenInfo,
          *fDPhiDEtaDenominator,
          fDPhiDEtaHiddenDenominator,
          fDPhiDEtaHiddenPrimaryDenominator,
          fDPhiDEtaHiddenPrimaryDenominatorData,
          fDPhiDEtaHiddenSecWeakDenominator,
          fDPhiDEtaHiddenSecWeakDenominatorData,
          fDPhiDEtaHiddenSecMatDenominator,
          fDPhiDEtaHiddenSecMatDenominatorData);
}

void AliFemtoCorrFctnDEtaDPhiSimple::WriteHistos()
{
  // Write out result histograms
  fDPhiDEtaNumerator->Write();
  fDPhiDEtaDenominator->Write();
  if (fReadHiddenInfo) {
    fDPhiDEtaHiddenNumerator->Write();
    fDPhiDEtaHiddenDenominator->Write();

    fDPhiDEtaHiddenPrimaryNumerator->Write();
    fDPhiDEtaHiddenPrimaryDenominator->Write();

    fDPhiDEtaHiddenSecWeakNumerator->Write();
    fDPhiDEtaHiddenSecWeakDenominator->Write();

    fDPhiDEtaHiddenSecMatNumerator->Write();
    fDPhiDEtaHiddenSecMatDenominator->Write();

    fDPhiDEtaHiddenPrimaryNumeratorData->Write();
    fDPhiDEtaHiddenPrimaryDenominatorData->Write();

    fDPhiDEtaHiddenSecWeakNumeratorData->Write();
    fDPhiDEtaHiddenSecWeakDenominatorData->Write();

    fDPhiDEtaHiddenSecMatNumeratorData->Write();
    fDPhiDEtaHiddenSecMatDenominatorData->Write();
  }
}

TList *AliFemtoCorrFctnDEtaDPhiSimple::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fDPhiDEtaNumerator);
  tOutputList->Add(fDPhiDEtaDenominator);
  if (fReadHiddenInfo) {
    tOutputList->Add(fDPhiDEtaHiddenNumerator);
    tOutputList->Add(fDPhiDEtaHiddenDenominator);

    tOutputList->Add(fDPhiDEtaHiddenPrimaryNumerator);
    tOutputList->Add(fDPhiDEtaHiddenPrimaryDenominator);

    tOutputList->Add(fDPhiDEtaHiddenSecWeakNumerator);
    tOutputList->Add(fDPhiDEtaHiddenSecWeakDenominator);

    tOutputList->Add(fDPhiDEtaHiddenSecMatNumerator);
    tOutputList->Add(fDPhiDEtaHiddenSecMatDenominator);

    tOutputList->Add(fDPhiDEtaHiddenPrimaryNumeratorData);
    tOutputList->Add(fDPhiDEtaHiddenPrimaryDenominatorData);

    tOutputList->Add(fDPhiDEtaHiddenSecWeakNumeratorData);
    tOutputList->Add(fDPhiDEtaHiddenSecWeakDenominatorData);

    tOutputList->Add(fDPhiDEtaHiddenSecMatNumeratorData);
    tOutputList->Add(fDPhiDEtaHiddenSecMatDenominatorData);
  }

  return tOutputList;
}

void AliFemtoCorrFctnDEtaDPhiSimple::SetParticleTypes(ParticleType partType1, ParticleType partType2)
{
  part1 = partType1;
  part2 = partType2;
}

void AliFemtoCorrFctnDEtaDPhiSimple::SetParticle1Type(ParticleType partType)
{
  part1 = partType;
}

void AliFemtoCorrFctnDEtaDPhiSimple::SetParticle2Type(ParticleType partType)
{
  part2 = partType;
}

void AliFemtoCorrFctnDEtaDPhiSimple::SetReadHiddenInfo(bool read)
{
  if (read == fReadHiddenInfo) {
    return;
  }
  fReadHiddenInfo = read;

  // we had hidden-info, now we don't - delete the
  if (!fReadHiddenInfo) {
    delete fDPhiDEtaHiddenNumerator;
    delete fDPhiDEtaHiddenDenominator;
    delete fDPhiDEtaHiddenPrimaryNumerator;
    delete fDPhiDEtaHiddenPrimaryDenominator;
    delete fDPhiDEtaHiddenSecWeakNumerator;
    delete fDPhiDEtaHiddenSecWeakDenominator;
    delete fDPhiDEtaHiddenSecMatNumerator;
    delete fDPhiDEtaHiddenSecMatDenominator;
    delete fDPhiDEtaHiddenPrimaryNumeratorData;
    delete fDPhiDEtaHiddenPrimaryDenominatorData;
    delete fDPhiDEtaHiddenSecWeakNumeratorData;
    delete fDPhiDEtaHiddenSecWeakDenominatorData;
    delete fDPhiDEtaHiddenSecMatNumeratorData;
    delete fDPhiDEtaHiddenSecMatDenominatorData;

    fDPhiDEtaHiddenNumerator = nullptr;
    fDPhiDEtaHiddenDenominator = nullptr;
    fDPhiDEtaHiddenPrimaryNumerator = nullptr;
    fDPhiDEtaHiddenPrimaryDenominator = nullptr;
    fDPhiDEtaHiddenSecWeakNumerator = nullptr;
    fDPhiDEtaHiddenSecWeakDenominator = nullptr;
    fDPhiDEtaHiddenSecMatNumerator = nullptr;
    fDPhiDEtaHiddenSecMatDenominator = nullptr;
    fDPhiDEtaHiddenPrimaryNumeratorData = nullptr;
    fDPhiDEtaHiddenPrimaryDenominatorData = nullptr;
    fDPhiDEtaHiddenSecWeakNumeratorData = nullptr;
    fDPhiDEtaHiddenSecWeakDenominatorData = nullptr;
    fDPhiDEtaHiddenSecMatNumeratorData = nullptr;
    fDPhiDEtaHiddenSecMatDenominatorData = nullptr;

    return;
  }

  int aEtaBins = fEtaBins;
  int aPhiBins = fPhiBins;

  TString ax_title = "; #Delta#Phi; #Delta#eta";

  fDPhiDEtaHiddenNumerator = new TH2D("NumDPhiDEtaHidden" + fTitle,
                                      "Numerator - All pairs (MonteCarlo Momentum)" + ax_title,
                                      aPhiBins, fphiL, fphiT,
                                      aEtaBins, -2.0, 2.0);
  fDPhiDEtaHiddenDenominator = new TH2D("DenDPhiDEtaHidden" + fTitle,
                                        "Denominator - All pairs (MonteCarlo Momentum)" + ax_title,
                                        aPhiBins, fphiL, fphiT,
                                        aEtaBins, -2.0, 2.0);

  fDPhiDEtaHiddenPrimaryNumerator = new TH2D("NumDPhiDEtaHiddenPrimary" + fTitle,
                                             "Numerator - Pairs of Primary Particles (MonteCarlo Momentum)" + ax_title,
                                             aPhiBins, fphiL, fphiT,
                                             aEtaBins, -2.0, 2.0);

  fDPhiDEtaHiddenPrimaryDenominator = new TH2D("DenDPhiDEtaHiddenPrimary" + fTitle,
                                               "Denominator - Pairs of Primary Particles (MonteCarlo Momentum)" + ax_title,
                                               aPhiBins, fphiL, fphiT,
                                               aEtaBins, -2.0, 2.0);

  fDPhiDEtaHiddenSecWeakNumerator = new TH2D("NumDPhiDEtaHiddenSecWeak" + fTitle,
                                             "Numerator - Pairs With a Secondary (weak) Particle (MonteCarlo Momentum)" + ax_title,
                                             aPhiBins, fphiL, fphiT,
                                             aEtaBins, -2.0, 2.0);

  fDPhiDEtaHiddenSecWeakDenominator = new TH2D("DenDPhiDEtaHiddenSecWeak" + fTitle,
                                               "Denominator - Pairs With a Secondary (weak) Particle (MonteCarlo Momentum)" + ax_title,
                                               aPhiBins, fphiL, fphiT,
                                               aEtaBins, -2.0, 2.0);

  fDPhiDEtaHiddenSecMatNumerator = new TH2D("NumDPhiDEtaHiddenSecMat" + fTitle,
                                            "Numerator - Pairs With a Secondary (material) Particle (MonteCarlo momentum)" + ax_title,
                                            aPhiBins, fphiL, fphiT,
                                            aEtaBins, -2.0, 2.0);

  fDPhiDEtaHiddenSecMatDenominator = new TH2D("DenDPhiDEtaHiddenSecMat" + fTitle,
                                              "Denominator - Pairs With a Secondary (material) Particle (MonteCarlo momentum)" + ax_title,
                                              aPhiBins, fphiL, fphiT,
                                              aEtaBins, -2.0, 2.0);

  fDPhiDEtaHiddenPrimaryNumeratorData = new TH2D("NumDPhiDEtaHiddenPrimaryData" + fTitle,
                                                 "Numerator - Pairs of Primary Particles (Reconstructed Momentum)" + ax_title,
                                                  aPhiBins, fphiL, fphiT,
                                                  aEtaBins, -2.0, 2.0);

  fDPhiDEtaHiddenPrimaryDenominatorData = new TH2D("DenDPhiDEtaHiddenPrimaryData" + fTitle,
                                                   "Denominator - Pairs of Primary Particles (Reconstructed Momentum)" + ax_title,
                                                   aPhiBins, fphiL, fphiT,
                                                   aEtaBins, -2.0, 2.0);

  fDPhiDEtaHiddenSecWeakNumeratorData = new TH2D("NumDPhiDEtaHiddenSecWeakData" + fTitle,
                                                 "Numerator - Pairs With a Secondary (weak) Particle (Reconstructed Momentum)" + ax_title,
                                                 aPhiBins, fphiL, fphiT,
                                                 aEtaBins, -2.0, 2.0);

  fDPhiDEtaHiddenSecWeakDenominatorData = new TH2D("DenDPhiDEtaHiddenSecWeakData" + fTitle,
                                                   "Denominator - Pairs With a Secondary (weak) Particle (Reconstructed Momentum)" + ax_title,
                                                   aPhiBins, fphiL, fphiT,
                                                   aEtaBins, -2.0, 2.0);

  fDPhiDEtaHiddenSecMatNumeratorData = new TH2D("NumDPhiDEtaHiddenSecMatData" + fTitle,
                                                "Numerator - Pairs of Seconadry (material) Particle (Reconstructed momentum)" + ax_title,
                                                aPhiBins, fphiL, fphiT,
                                                aEtaBins, -2.0, 2.0);

  fDPhiDEtaHiddenSecMatDenominatorData = new TH2D("DenDPhiDEtaHiddenSecMatData" + fTitle,
                                                  "Denominator - Pairs With a Secondary (material) Particle (Reconstructed momentum)" + ax_title,
                                                  aPhiBins, fphiL, fphiT,
                                                  aEtaBins, -2.0, 2.0);

  fDPhiDEtaHiddenNumerator->Sumw2();
  fDPhiDEtaHiddenPrimaryNumerator->Sumw2();
  fDPhiDEtaHiddenSecWeakNumerator->Sumw2();
  fDPhiDEtaHiddenSecMatNumerator->Sumw2();
  fDPhiDEtaHiddenPrimaryNumeratorData->Sumw2();
  fDPhiDEtaHiddenSecWeakNumeratorData->Sumw2();
  fDPhiDEtaHiddenSecMatNumeratorData->Sumw2();
}
