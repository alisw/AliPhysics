///
/// \file AliFemtoUser/AliFemtoCorrFctn3DLCMSPosQuad.cxx
///

#include "AliFemtoCorrFctn3DLCMSPosQuad.h"

#include <TH3F.h>


AliFemtoCorrFctn3DLCMSPosQuad::AliFemtoCorrFctn3DLCMSPosQuad(const Parameters &p)
  : AliFemtoCorrFctn3DLCMSPosQuad(p.prefix, p.suffix,
                                  p.nbins_out, p.nbins_side, p.nbins_long,
                                  p.QoHi, p.QsHi, p.QlHi)
{
}

AliFemtoCorrFctn3DLCMSPosQuad::AliFemtoCorrFctn3DLCMSPosQuad(const TString &name, const int nbins, const float QHi)
  : AliFemtoCorrFctn3DLCMSPosQuad(name, "", nbins / 2 + 1, nbins, nbins, QHi, QHi, QHi)
{
}

AliFemtoCorrFctn3DLCMSPosQuad::AliFemtoCorrFctn3DLCMSPosQuad(const TString &prefix,
                                                             const TString &suffix,
                                                             const size_t nbins_out,
                                                             const size_t nbins_side,
                                                             const size_t nbins_long,
                                                             const float QoHi,
                                                             const float QsHi,
                                                             const float QlHi)
  : AliFemtoCorrFctn()
  , fNumerator(nullptr)
  , fDenominator(nullptr)
  , fQinvWeight(nullptr)
{
  auto hist_name = [&] (TString name)
    {
      return prefix + name + suffix;
    };

  auto new_hist = [&] (TString name, TString title)
    {
      return new TH3F(hist_name(name), title,
                      nbins_out, 0.0, QoHi,
                      nbins_side, -QsHi, QsHi,
                      nbins_long, -QlHi, QlHi);
    };

  fNumerator = new_hist("Num", "Numerator; q_{out}; q_{side}; q_{long}");
  fDenominator = new_hist("Den", "Denominator; q_{out}; q_{side}; q_{long}");
  fQinvWeight = new_hist("Qinv", "q_{inv} Weighted Denominator; q_{out}; q_{side}; q_{long}");
  fQinvWeight->Sumw2();
}


AliFemtoCorrFctn3DLCMSPosQuad
  ::AliFemtoCorrFctn3DLCMSPosQuad(const TString &prefix,
                                  const TString &suffix,
                                  const std::vector<double> &out_bins,
                                  const std::vector<double> &side_bins,
                                  const std::vector<double> &long_bins)
  : AliFemtoCorrFctn()
  , fNumerator(nullptr)
  , fDenominator(nullptr)
  , fQinvWeight(nullptr)
{
  auto hist_name = [&] (TString name)
    {
      return prefix + name + suffix;
    };

  auto new_hist = [&] (TString name, TString title)
    {
      return new TH3F(hist_name(name), title,
                      out_bins.size() - 1, out_bins.data(),
                      side_bins.size() - 1, side_bins.data(),
                      long_bins.size() - 1, long_bins.data());
    };

  fNumerator = new_hist("Num", "Numerator; q_{out}; q_{side}; q_{long}");
  fDenominator = new_hist("Den", "Denominator; q_{out}; q_{side}; q_{long}");
  fQinvWeight = new_hist("Qinv", "q_{inv} Weighted Denominator; q_{out}; q_{side}; q_{long}");
  fQinvWeight->Sumw2();
}

AliFemtoCorrFctn3DLCMSPosQuad::AliFemtoCorrFctn3DLCMSPosQuad(const AliFemtoCorrFctn3DLCMSPosQuad &orig)
  : AliFemtoCorrFctn(orig)
  , fNumerator(new TH3F(*orig.fNumerator))
  , fDenominator(new TH3F(*orig.fDenominator))
  , fQinvWeight(new TH3F(*orig.fQinvWeight))
{
}

AliFemtoCorrFctn3DLCMSPosQuad& AliFemtoCorrFctn3DLCMSPosQuad::operator=(const AliFemtoCorrFctn3DLCMSPosQuad &rhs)
{
  if (this != &rhs) {
    *fNumerator = *rhs.fNumerator;
    *fDenominator = *rhs.fDenominator;
    *fQinvWeight = *rhs.fQinvWeight;
  }
  return *this;
}

AliFemtoCorrFctn3DLCMSPosQuad::~AliFemtoCorrFctn3DLCMSPosQuad()
{
  delete fNumerator;
  delete fDenominator;
  delete fQinvWeight;
}

inline void fill_histograms(const AliFemtoPair &pair, TH3 &dest, TH3 *qinv_dest=nullptr)
{
  const double
    qo = pair.QOutCMS(),
    qs = pair.QSideCMS(),
    ql = pair.QLongCMS(),
    // flip about origin if qout negative
    f = std::copysign(1.0, qo),

    qOut = f * qo,
    qSide = f * qs,
    qLong = f * ql;

  dest.Fill(qOut, qSide, qLong);

  if (qinv_dest) {
    qinv_dest->Fill(qOut, qSide, qLong, fabs(pair.QInv()));
  }
}

AliFemtoString
AliFemtoCorrFctn3DLCMSPosQuad::Report()
{
  AliFemtoString report("AliFemtoCorrFctn3DLCMSPosQuad Report\n");
  report += Form("Number of entries in numerator:\t%E\n", fNumerator->GetEntries());
  report += Form("Number of entries in denominator:\t%E\n", fDenominator->GetEntries());

  return report;
}


void AliFemtoCorrFctn3DLCMSPosQuad::AddRealPair(const AliFemtoPair &pair)
{
  if (fPairCut && !fPairCut->Pass(&pair)) {
    return;
  }
  fill_histograms(pair, *fNumerator);
}


void AliFemtoCorrFctn3DLCMSPosQuad::AddMixedPair(const AliFemtoPair &pair)
{
  if (fPairCut && !fPairCut->Pass(&pair)) {
    return;
  }
  fill_histograms(pair, *fDenominator, fQinvWeight);
}


TList* AliFemtoCorrFctn3DLCMSPosQuad::GetOutputList()
{
  /// Prepare the list of objects to be written to the output

  TList *output_list = new TList();

  output_list->Add(fNumerator);
  output_list->Add(fDenominator);
  output_list->Add(fQinvWeight);

  return output_list;
}
