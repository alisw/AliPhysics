///
/// \file AliFemtoModelCorrFctnDEtaDPhiStar.cxx
///

#include "AliFemtoModelCorrFctnDEtaDPhiStar.h"
#include "AliFemtoPair.h"
#include "AliFemtoModelManager.h"

#include <TObjArray.h>


AliFemtoModelCorrFctnDEtaDPhiStar::Parameters
AliFemtoModelCorrFctnDEtaDPhiStar::Parameters::Default = {
  TString(""), // title
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

  const TString hist_title = p.title + ";" + xAxisTitle + ";" + yAxisTitle;

  // numerator
  fDPhiStarDEtaNumeratorTrue = new TH2F("NumDPhiStarDEtaTrue" + p.title, hist_title,
                                    p.eta.bin_count, p.eta.low, p.eta.high,
                                    p.phi.bin_count, p.phi.low, p.phi.high);
  fDPhiStarDEtaNumeratorFake = new TH2F("NumDPhiStarDEtaFake" + p.title, hist_title,
                                    p.eta.bin_count, p.eta.low, p.eta.high,
                                    p.phi.bin_count, p.phi.low, p.phi.high);

  fDPhiStarDEtaDenominator = new TH2F("DenDPhiStarDEta" + p.title, hist_title,
                                    p.eta.bin_count, p.eta.low, p.eta.high,
                                    p.phi.bin_count, p.phi.low, p.phi.high);
  fDPhiStarDEtaColNumerator = new TH2F("NumDPhiDEtaCol"+ p.title, hist_title,
                                    p.eta.bin_count, p.eta.low, p.eta.high,
                                    p.phi.bin_count, p.phi.low, p.phi.high);
  fDPhiStarDEtaColDenominator = new TH2F("DenDPhiDEtaCol" + p.title, hist_title,
                                         p.eta.bin_count, p.eta.low, p.eta.high,
                                         p.phi.bin_count, p.phi.low, p.phi.high);

  fDPhiStarNumeratorTrue = new TH1F("NumDPhiTrue" + p.title, "#Delta#phi* Numerator_{true}; #Delta#phi*",
                                    p.phi.bin_count, p.phi.low, p.phi.high);
  fDPhiStarNumeratorTrue = new TH1F("NumDPhiFake" + p.title,"#Delta#phi* Numerator_{fake}; #Delta#phi*",
                                    p.phi.bin_count, p.phi.low, p.phi.high);
  fDPhiStarDenominator = new TH1F("DenDPhi" + p.title, "#Delta#phi* Denominator; #Delta#phi*",
                                  p.phi.bin_count, p.phi.low, p.phi.high);

  fDCosNumeratorTrue = new TH1F("NumDCosTrue_" + p.title, "Cos Numerator_{true}; Cos(#phi)",
                                p.phi.bin_count * 2, -1.0, 1.0);
  fDCosNumeratorFake = new TH1F("NumDCosFake_" + p.title, "Cos Numerator_{fake}; Cos(#phi)",
                                p.phi.bin_count * 2, -1.0, 1.0);
  fDCosDenominator = new TH1F("DenDCos_" + p.title, "Cos() Denominator; Cos(#phi)",
                              p.phi.bin_count * 2, -1.0, 1.0);

  fDPhiStarPtNumerator = new TH2F("NumDPhiPt_" + p.title, "#Delta#phi* vs p_T Numinator; #Delta#phi*; p_T (GeV)",
                                  p.phi.bin_count * 2, p.phi.low, p.phi.high,
                                  30, 0.0, 3.0);
  fDPhiStarPtDenominator = new TH2F("DenDPhiPt_" + p.title, "#Delta#phi* vs p_T Numinator; #Delta#phi*; p_T (GeV)",
                                    p.phi.bin_count*2, p.phi.low, p.phi.high,
                                    30, 0.0, 3.0);

  fDCosPtNumerator = new TH2F("NumDCosPt_" + p.title, "Cos vs p_{T} Numinator; cos(#phi); p_T (GeV)",
                              p.phi.bin_count * 2, -1.0, 1.0,
                              30, 0.0, 3.0);
  fDCosPtDenominator = new TH2F("DenDCosPt_" + p.title, "Cos vs P_{T} Denominator; cos(#phi); p_T (GeV)",
                                p.phi.bin_count * 2, -1.0, 1.0,
                                30, 0.0, 3.0);
}

AliFemtoModelCorrFctnDEtaDPhiStar::AliFemtoModelCorrFctnDEtaDPhiStar(const char* title, const UInt_t phi_bins, const UInt_t eta_bins):
  AliFemtoModelCorrFctnDEtaDPhiStar(Parameters::Default.WithTitle(title).WithPhiBins(phi_bins).WithEtaBins(eta_bins))
{
  // no-op
}

AliFemtoModelCorrFctnDEtaDPhiStar::~AliFemtoModelCorrFctnDEtaDPhiStar()
{
  delete fDPhiStarDEtaNumeratorTrue;
  delete fDPhiStarDEtaNumeratorFake;
  delete fDPhiStarDEtaDenominator;
  delete fDPhiStarDEtaColNumerator;
  delete fDPhiStarDEtaColDenominator;
  delete fDPhiStarNumeratorTrue;
  delete fDPhiStarNumeratorFake;
  delete fDPhiStarDenominator;
  delete fDCosNumeratorTrue;
  delete fDCosNumeratorFake;
  delete fDCosDenominator;
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

  output->Add(fDPhiStarDEtaNumeratorTrue);
  output->Add(fDPhiStarDEtaNumeratorFake);
  output->Add(fDPhiStarDEtaDenominator);
  output->Add(fDPhiStarDEtaColNumerator);
  output->Add(fDPhiStarDEtaColDenominator);
  output->Add(fDPhiStarNumeratorTrue);
  output->Add(fDPhiStarNumeratorFake);
  output->Add(fDPhiStarDenominator);
  output->Add(fDCosNumeratorTrue);
  output->Add(fDCosNumeratorFake);
  output->Add(fDCosDenominator);
  output->Add(fDPhiStarPtNumerator);
  output->Add(fDPhiStarPtDenominator);
  output->Add(fDCosPtNumerator);
  output->Add(fDCosPtDenominator);

  return &output_list;
}

inline
static double GetCosPhi(const AliFemtoTrack& t1, const AliFemtoTrack& t2)
{
  const auto p1 = t1.P(), p2 = t2.P();

  const double cosphi = p1.Dot(p2) / std::sqrt(p1.Mag2() * p2.Mag2());
  return cosphi;
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

  fDPhiStarDEtaNumeratorTrue->Fill(delta.eta, delta.phi_star, weight);

  // not sure what this is (See AliFemtoModelCorrFctnDEtaDPhi::AddMixedPair)
  auto alt_deta = cosphi > 0 ? -delta.eta : delta.eta - 2.0 * track2.P().PseudoRapidity();
  fDPhiStarDEtaColNumerator->Fill(alt_deta, delta.phi_star);

  fDPhiStarNumeratorTrue->Fill(delta.phi_star, weight);
  fDCosNumeratorTrue->Fill(cosphi, weight);
  fDPhiStarPtNumerator->Fill(delta.phi_star, ptmin);
  fDCosPtNumerator->Fill(cosphi, ptmin);
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

  const Double_t weight = fManager->GetWeight(pair);

  fDPhiStarDEtaNumeratorFake->Fill(delta.eta, delta.phi_star, weight);
  fDPhiStarDEtaDenominator->Fill(delta.eta, delta.phi_star);

  auto alt_deta = cosphi > 0 ? -delta.eta : delta.eta - 2.0 * track2.P().PseudoRapidity();
  fDPhiStarDEtaColDenominator->Fill(alt_deta, delta.phi_star);

  fDPhiStarNumeratorFake->Fill(delta.phi_star, weight);
  fDCosNumeratorFake->Fill(cosphi, weight);

  fDPhiStarDenominator->Fill(delta.phi_star);
  fDCosDenominator->Fill(cosphi);

  fDPhiStarPtDenominator->Fill(delta.phi_star, ptmin);
  fDCosPtDenominator->Fill(cosphi, ptmin);

}

AliFemtoString
AliFemtoModelCorrFctnDEtaDPhiStar::Report()
{
  TString report;
  report += TString::Format("Number of entries in numerator true:\t%lu\n", fDPhiStarDEtaNumeratorTrue->GetEntries());
  report += TString::Format("Number of entries in numerator fake:\t%E\n", fDPhiStarDEtaNumeratorFake->GetEntries());
  report += TString::Format("Number of entries in denominator:\t%E\n", fDPhiStarDEtaDenominator->GetEntries());
  return AliFemtoString(report);
}
