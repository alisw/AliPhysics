///
/// \file AliFemtoModelCorrFctnDEtaDPhiStar.cxx
///

#include "AliFemtoModelCorrFctnDEtaDPhiStar.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"
#include "AliFemtoModelHiddenInfo.h"


#include <TObjArray.h>
#include <cmath>


AliFemtoModelCorrFctnDEtaDPhiStar::Parameters
AliFemtoModelCorrFctnDEtaDPhiStar::Parameters::Default = {
  TString("DetaDphistar"), // title
  {73, -0.1, 0.1}, // phi hist data
  {73, -0.1, 0.1}, // eta hist data
  1.2, // radius
  false // group_output
};

AliFemtoModelCorrFctnDEtaDPhiStar*
AliFemtoModelCorrFctnDEtaDPhiStar::Builder::Build() const
{
  return static_cast<AliFemtoModelCorrFctnDEtaDPhiStar*>(*this);
}


AliFemtoModelCorrFctnDEtaDPhiStar::AliFemtoModelCorrFctnDEtaDPhiStar(
  const AliFemtoModelCorrFctnDEtaDPhiStar::Parameters &p):
    fRadius(p.radius),
    fCurrentMagneticField(0.0),
    fGroupOutputList(p.group_output),
    fTitle(p.title)
{
  const TString xAxisTitle = "#Delta#eta",
                yAxisTitle = "#Delta#phi*";

  // const TString hist_title = p.title + ";" + xAxisTitle + ";" + yAxisTitle;
  const auto hist_title = [] (const TString &title) { return title + ";#Delta#eta;#Delta#phi*"; };

  const TString fix = (fGroupOutputList == true) ? "" : p.title + "_";

  auto hist_name = [&] (const TString &title) { return fix + title; };

  fDPhiStarDEtaNumeratorWeighted = new TH2F(hist_name("NumDPhiStarDEtaWeighted"),
                                            hist_title("Weighted Real Pairs"),
                                            p.eta.bin_count, p.eta.low, p.eta.high,
                                            p.phi.bin_count, p.phi.low, p.phi.high);
  fDPhiStarDEtaNumeratorUnweighted = new TH2F(hist_name("NumDPhiStarDEtaUnweighted"),
                                              hist_title("Unweighted Real Pairs"),
                                              p.eta.bin_count, p.eta.low, p.eta.high,
                                              p.phi.bin_count, p.phi.low, p.phi.high);

  fDPhiStarDEtaDenominator = new TH2F(hist_name("DenDPhiStarDEta"),
                                      hist_title("Mixed Pairs"),
                                      p.eta.bin_count, p.eta.low, p.eta.high,
                                      p.phi.bin_count, p.phi.low, p.phi.high);


  fDPhiStarDEtaNumeratorIdealWeighted = new TH2F(hist_name("NumDPhiStarDEtaIdealWeighted"),
                                            hist_title("(Ideal) Weighted Real Pairs"),
                                            p.eta.bin_count, p.eta.low, p.eta.high,
                                            p.phi.bin_count, p.phi.low, p.phi.high);
  fDPhiStarDEtaNumeratorIdealUnweighted = new TH2F(hist_name("NumDPhiStarDEtaIdealUnweighted"),
                                              hist_title("(Ideal) Unweighted Real Pairs"),
                                              p.eta.bin_count, p.eta.low, p.eta.high,
                                              p.phi.bin_count, p.phi.low, p.phi.high);

  fDPhiStarDEtaIdealDenominator = new TH2F(hist_name("DenDPhiStarDEtaIdeal"),
                                      hist_title("(Ideal) Mixed Pairs"),
                                      p.eta.bin_count, p.eta.low, p.eta.high,
                                      p.phi.bin_count, p.phi.low, p.phi.high);


  fDPhiStarDEtaColNumerator = new TH2F(hist_name("NumDPhiDEtaCol"),
                                       hist_title("Colinear #Delta#Phi* Numerator"),
                                       p.eta.bin_count, p.eta.low, p.eta.high,
                                       p.phi.bin_count, p.phi.low, p.phi.high);
  fDPhiStarDEtaColDenominator = new TH2F(hist_name("DenDPhiDEtaCol"),
                                         hist_title("Colinear #Delta#Phi* Denominator"),
                                         p.eta.bin_count, p.eta.low, p.eta.high,
                                         p.phi.bin_count, p.phi.low, p.phi.high);

  const UInt_t pt_bincount = 60,
               cosphi_bincount = p.phi.bin_count * 2;

  fDPhiStarPtNumerator = new TH2F(hist_name("NumDPhiPt"),
                                  "#Delta#phi* vs p_{T} Numinator; #Delta#phi*; p_{T} (GeV)",
                                   p.phi.bin_count, p.phi.low, p.phi.high,
                                   pt_bincount, 0.0, 3.0);
  fDPhiStarPtDenominator = new TH2F(hist_name("DenDPhiPt"),
                                    "#Delta#phi* vs p_{T} Denominator; #Delta#phi*; p_{T} (GeV)",
                                    p.phi.bin_count, p.phi.low, p.phi.high,
                                    pt_bincount, 0.0, 3.0);

  fDCosPtNumerator = new TH2F(hist_name("NumDCosPt"),
                              "Cos vs p_{T} Numinator; cos(#phi); p_T (GeV)",
                              cosphi_bincount, -1.0, 1.0,
                              pt_bincount, 0.0, 3.0);
  fDCosPtDenominator = new TH2F(hist_name("DenDCosPt"),
                                "Cos vs P_{T} Denominator; cos(#phi); p_T (GeV)",
                                cosphi_bincount, -1.0, 1.0,
                                pt_bincount, 0.0, 3.0);
}

AliFemtoModelCorrFctnDEtaDPhiStar::AliFemtoModelCorrFctnDEtaDPhiStar(const char* title, const UInt_t phi_bins, const UInt_t eta_bins):
  AliFemtoModelCorrFctnDEtaDPhiStar(Parameters::Default.WithTitle(title).WithPhiBins(phi_bins).WithEtaBins(eta_bins))
{
  // no-op
}

AliFemtoModelCorrFctnDEtaDPhiStar::~AliFemtoModelCorrFctnDEtaDPhiStar()
{
  delete fDPhiStarDEtaNumeratorWeighted;
  delete fDPhiStarDEtaNumeratorUnweighted;
  delete fDPhiStarDEtaDenominator;

  delete fDPhiStarDEtaColNumerator;
  delete fDPhiStarDEtaColDenominator;

  delete fDPhiStarPtNumerator;
  delete fDPhiStarPtDenominator;

  delete fDCosPtNumerator;
  delete fDCosPtDenominator;
}

TList*
AliFemtoModelCorrFctnDEtaDPhiStar::GetOutputList()
{
  return AppendOutputList(*(new TList()));
}

TList*
AliFemtoModelCorrFctnDEtaDPhiStar::AppendOutputList(TList &output_list)
{
  TCollection *output = &output_list;
  if (fGroupOutputList) {
    output = new TObjArray();
    output->SetName(fTitle);
    output_list.Add(output);
  }

  output->Add(fDPhiStarDEtaNumeratorWeighted);
  output->Add(fDPhiStarDEtaNumeratorUnweighted);
  output->Add(fDPhiStarDEtaDenominator);

  output->Add(fDPhiStarDEtaNumeratorIdealWeighted);
  output->Add(fDPhiStarDEtaNumeratorIdealUnweighted);
  output->Add(fDPhiStarDEtaIdealDenominator);

  output->Add(fDPhiStarDEtaColNumerator);
  output->Add(fDPhiStarDEtaColDenominator);

  output->Add(fDPhiStarPtNumerator);
  output->Add(fDPhiStarPtDenominator);

  output->Add(fDCosPtNumerator);
  output->Add(fDCosPtDenominator);

  return &output_list;
}

inline
static double GetCosPhi(const AliFemtoTrack& t1, const AliFemtoTrack& t2)
{
  const auto p1 = t1.P(),
             p2 = t2.P();

  const double cosphi = p1.Dot(p2) / std::sqrt(p1.Mag2() * p2.Mag2());
  return cosphi;
}


inline static
EtaPhiStar
ApplyIdealCalculations(const AliFemtoTrack &t1, const AliFemtoTrack &t2, float radius, float magfield)
{
  // ideal momentum
  AliFemtoModelHiddenInfo *mc1 = dynamic_cast<AliFemtoModelHiddenInfo*>(t1.GetHiddenInfo()),
                          *mc2 = dynamic_cast<AliFemtoModelHiddenInfo*>(t2.GetHiddenInfo());

  if (mc1 == nullptr || mc2 == nullptr) {
      return {NAN, NAN};
  }

  const AliFemtoThreeVector &p1 = *mc1->GetTrueMomentum(),
                            &p2 = *mc2->GetTrueMomentum();

  // skip no-momentum particles
  if (p1.Mag2() == 0 || p2.Mag2() == 0) {
    return {NAN, NAN};
  }

  const Short_t charge1 = t1.Charge(),
                charge2 = t2.Charge();

  const Float_t delta_eta = AliFemtoPairCutDetaDphi::CalculateDEta(p1, p2),
                delta_phi_star = AliFemtoPairCutDetaDphi::CalculateDPhiStar(p1, charge1,
                                                                            p2, charge2,
                                                                            radius,
                                                                            magfield);
  return {delta_eta, delta_phi_star};
}


void AliFemtoModelCorrFctnDEtaDPhiStar::AddRealPair(AliFemtoPair* pair)
{
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }
  const AliFemtoTrack &track1 = *pair->Track1()->Track(),
                      &track2 = *pair->Track2()->Track();

  auto delta = GetDeltaEtaDeltaPhiStar(track1, track2);

  const Double_t cosphi = GetCosPhi(track1, track2);
  const Double_t ptmin = std::min(track1.Pt(), track2.Pt());
  const Double_t weight = fManager->GetWeight(pair);

  fDPhiStarDEtaNumeratorWeighted->Fill(delta.eta, delta.phi_star, weight);
  fDPhiStarDEtaNumeratorUnweighted->Fill(delta.eta, delta.phi_star);

  // not sure what this is (See AliFemtoModelCorrFctnDEtaDPhi::AddMixedPair)
  auto alt_deta = cosphi > 0 ? -delta.eta : delta.eta - 2.0 * track2.P().PseudoRapidity();
  fDPhiStarDEtaColNumerator->Fill(alt_deta, delta.phi_star);

  fDPhiStarPtNumerator->Fill(delta.phi_star, ptmin);
  fDCosPtNumerator->Fill(cosphi, ptmin);

  auto ideal_delta = ApplyIdealCalculations(track1, track2, fRadius, fCurrentMagneticField);
  if (!std::isnan(ideal_delta.eta)) {
    fDPhiStarDEtaNumeratorIdealWeighted->Fill(ideal_delta.eta, ideal_delta.phi_star, weight);
    fDPhiStarDEtaNumeratorIdealUnweighted->Fill(ideal_delta.eta, ideal_delta.phi_star);
  }
}


void AliFemtoModelCorrFctnDEtaDPhiStar::AddMixedPair(AliFemtoPair* pair)
{
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  const AliFemtoTrack &track1 = *pair->Track1()->Track(),
                      &track2 = *pair->Track2()->Track();

  auto delta = GetDeltaEtaDeltaPhiStar(track1, track2);

  const Double_t cosphi = GetCosPhi(track1, track2);
  const Double_t ptmin = std::min(track1.Pt(), track2.Pt());

  fDPhiStarDEtaDenominator->Fill(delta.eta, delta.phi_star);

  auto alt_deta = cosphi > 0 ? -delta.eta : delta.eta - 2.0 * track2.P().PseudoRapidity();
  fDPhiStarDEtaColDenominator->Fill(alt_deta, delta.phi_star);

  fDPhiStarPtDenominator->Fill(delta.phi_star, ptmin);
  fDCosPtDenominator->Fill(cosphi, ptmin);

  auto ideal_delta = ApplyIdealCalculations(track1, track2, fRadius, fCurrentMagneticField);
  if (!std::isnan(ideal_delta.eta)) {
    fDPhiStarDEtaIdealDenominator->Fill(ideal_delta.eta, ideal_delta.phi_star);
  }
}

AliFemtoString
AliFemtoModelCorrFctnDEtaDPhiStar::Report()
{
  TString report;
  report += TString::Format("Number of entries in weighted numerator:\t%lu\n", fDPhiStarDEtaNumeratorWeighted->GetEntries());
  report += TString::Format("Number of entries in unweighted numerator:\t%lu\n", fDPhiStarDEtaNumeratorUnweighted->GetEntries());
  report += TString::Format("Number of entries in denominator:\t%lu\n", fDPhiStarDEtaDenominator->GetEntries());

  report += TString::Format("Number of entries in ideal weighted numerator:\t%lu\n", fDPhiStarDEtaNumeratorIdealWeighted->GetEntries());
  report += TString::Format("Number of entries in ideal unweighted numerator:\t%lu\n", fDPhiStarDEtaNumeratorIdealUnweighted->GetEntries());
  report += TString::Format("Number of entries in denominator:\t%lu\n", fDPhiStarDEtaIdealDenominator->GetEntries());

  return AliFemtoString(report);
}
