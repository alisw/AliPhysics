///
/// \file AliFemtoCorrFctnQLCMS.cxx
///


#include "AliFemtoCorrFctnQLCMS.h"

AliFemtoCorrFctnQLCMS::AliFemtoCorrFctnQLCMS(const Parameters &p)
  : AliFemtoCorrFctnQLCMS(p.prefix.Data(), p.suffix.Data(), p.nbins, p.qmin, p.qmax)
{
}

AliFemtoCorrFctnQLCMS::AliFemtoCorrFctnQLCMS(const char* prefix,
                                             const char* suffix,
                                             const int nbins,
                                             const float qhi)
  : AliFemtoCorrFctnQLCMS(prefix, suffix, nbins, -qhi, qhi)
{
}

AliFemtoCorrFctnQLCMS::AliFemtoCorrFctnQLCMS(const char* prefix,
                                             const char* suffix,
                                             const int nbins,
                                             float qlo,
                                             float qhi)
  : AliFemtoCorrFctn()
  , fNumerator(nullptr)
  , fDenominator(nullptr)
  , fQinv(nullptr)
{
  if (std::isnan(qlo)) {
    qlo = -qhi;
  }

  fNumerator = new TH1F(Form("%sNum%s", prefix, suffix),
                        "Numerator; q_{lcms} (GeV);",
                        nbins, qlo, qhi);

  fDenominator = new TH1F(Form("%sDen%s", prefix, suffix),
                          "Denominator; q_{lcms} (GeV);",
                          nbins, qlo, qhi);

  fQinv = new TProfile(Form("%sQinv%s", prefix, suffix),
                       "Qinv Weighted Profile (Num + Den); q_{lcms} (GeV); q_{inv} (GeV)",
                       nbins, qlo, qhi);

}

AliFemtoCorrFctnQLCMS::AliFemtoCorrFctnQLCMS(const AliFemtoCorrFctnQLCMS& orig)
  : AliFemtoCorrFctn(orig)
  , fNumerator(new TH1F(*orig.fNumerator))
  , fDenominator(new TH1F(*orig.fDenominator))
  , fQinv(static_cast<TProfile*>(orig.fQinv->Clone()))
{
}

AliFemtoCorrFctnQLCMS&
AliFemtoCorrFctnQLCMS::operator=(const AliFemtoCorrFctnQLCMS &rhs)
{
  if (this != &rhs) {
    AliFemtoCorrFctn::operator=(rhs);
    *fNumerator = *rhs.fNumerator;
    *fDenominator = *rhs.fDenominator;
    delete fQinv;
    fQinv = static_cast<TProfile*>(rhs.fQinv->Clone());
  }

  return *this;
}

AliFemtoCorrFctnQLCMS::~AliFemtoCorrFctnQLCMS()
{
  delete fNumerator;
  delete fDenominator;
  delete fQinv;
}

AliFemtoString
AliFemtoCorrFctnQLCMS::Report()
{
  // Construct the report
  AliFemtoString report = "=== Begin QLCMS Correlation Function";
  report += Form("Number of entries in numerator:\t%E\n", fNumerator->GetEntries());
  report += Form("Number of entries in denominator:\t%E\n", fDenominator->GetEntries());

  if (fPairCut) {
    report += "PairCut report:\n";
    report += fPairCut->Report();
  } else {
    report += "No PairCut specific to this CorrFctn\n";
  }
  report += "===\n";

  return report;
}
