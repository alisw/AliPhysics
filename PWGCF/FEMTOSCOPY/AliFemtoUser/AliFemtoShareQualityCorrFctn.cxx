////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoShareQualityCorrFctn - A correlation function that saves the     ///
/// amount of sharing and splitting hits per pair as a function of qinv      ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoShareQualityCorrFctn.h"
//#include "AliFemtoHisto.hh"
#include <tuple>

#ifdef __ROOT__
ClassImp(AliFemtoShareQualityCorrFctn)
#endif

//____________________________
AliFemtoShareQualityCorrFctn::AliFemtoShareQualityCorrFctn(const char* title, const int& nbins, const float& QinvLo, const float& QinvHi):
  AliFemtoCorrFctn(),
  fShareNumerator(0),
  fShareDenominator(0),
  fQualityNumerator(0),
  fQualityDenominator(0),
  fTPCSepNumerator(0),
  fTPCSepDenominator(0)
{

  const auto
     share_params = std::make_tuple("; q_{inv} (GeV); Share Fraction", 100, 0.0, 1.00001),
     qual_params = std::make_tuple("; q_{inv} (GeV); Share Quality", 150, -0.500001, 1.000001),
     sep_params = std::make_tuple("; q_{inv} (GeV); TPC Entrace Separation (cm)", 150, 0.0, 100.0);

  auto _new_hist = [=] (TString name, TString title, std::tuple<const char *, int, double, double> params)
    {
      return new TH2D(name, title + std::get<0>(params),
                      nbins, QinvLo, QinvHi,
                      std::get<1>(params), std::get<2>(params), std::get<3>(params));
    };

  fShareNumerator = _new_hist("NumShare", "Numerator Share Fraction vs Qinv", share_params);
  fShareDenominator = _new_hist("DenShare", "Denominator Share Fraction vs Qinv", share_params);

  fQualityNumerator = _new_hist("NumQuality", "Numerator Quality vs Qinv", qual_params);
  fQualityDenominator = _new_hist("DenQuality", "Denominator Quality vs Qinv", qual_params);

  fTPCSepNumerator = _new_hist("NumTPCSep", "Numerator TPC Separation vs Qinv", sep_params);
  fTPCSepDenominator = _new_hist("DenTPCSep", "Denominator TPC Separation vs Qinv", sep_params);

  // to enable error bar calculation...
  fShareNumerator->Sumw2();
  fShareDenominator->Sumw2();

  fQualityNumerator->Sumw2();
  fQualityDenominator->Sumw2();

  fTPCSepNumerator->Sumw2();
  fTPCSepDenominator->Sumw2();
}

//____________________________
AliFemtoShareQualityCorrFctn::AliFemtoShareQualityCorrFctn(const AliFemtoShareQualityCorrFctn& aCorrFctn) :
  AliFemtoCorrFctn(),
  fShareNumerator(0),
  fShareDenominator(0),
  fQualityNumerator(0),
  fQualityDenominator(0),
  fTPCSepNumerator(0),
  fTPCSepDenominator(0)
{
  // copy constructor
  fShareNumerator = new TH2D(*aCorrFctn.fShareNumerator);
  fShareDenominator = new TH2D(*aCorrFctn.fShareDenominator);
  fQualityNumerator = new TH2D(*aCorrFctn.fQualityNumerator);
  fQualityDenominator = new TH2D(*aCorrFctn.fQualityDenominator);
  fTPCSepNumerator = new TH2D(*aCorrFctn.fTPCSepNumerator);
  fTPCSepDenominator = new TH2D(*aCorrFctn.fTPCSepDenominator);
}
//____________________________
AliFemtoShareQualityCorrFctn::~AliFemtoShareQualityCorrFctn(){
  // destructor
  delete fShareNumerator;
  delete fShareDenominator;
  delete fQualityNumerator;
  delete fQualityDenominator;
  delete fTPCSepNumerator;
  delete fTPCSepDenominator;
}
//_________________________
AliFemtoShareQualityCorrFctn& AliFemtoShareQualityCorrFctn::operator=(const AliFemtoShareQualityCorrFctn& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  *fShareNumerator = *aCorrFctn.fShareNumerator;
  *fShareDenominator = *aCorrFctn.fShareDenominator;
  *fQualityNumerator = *aCorrFctn.fQualityNumerator;
  *fQualityDenominator = *aCorrFctn.fQualityDenominator;
  *fTPCSepNumerator = *aCorrFctn.fTPCSepNumerator;
  *fTPCSepDenominator = *aCorrFctn.fTPCSepDenominator;

  return *this;
}
//_________________________
void AliFemtoShareQualityCorrFctn::Finish()
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
AliFemtoString AliFemtoShareQualityCorrFctn::Report()
{
  // create report
  AliFemtoString report = "Qinv Correlation Function Report:\n";
  report += Form("Number of entries in numerator:\t%E\n", fShareNumerator->GetEntries());
  report += Form("Number of entries in denominator:\t%E\n", fShareDenominator->GetEntries());
  return report;
}

void AliFemtoShareQualityCorrFctn::AddPair(const AliFemtoPair* pair, TH2 &share, TH2 &qual, TH2 &sep)
{
  double tQinv = fabs(pair->QInv());   // note - qInv() will be negative for identical pairs...

  double share_fraction, share_quality;
  pair->CalcTrackShareQualFractions(share_fraction, share_quality);

  const auto
    &x1 = pair->Track1()->Track()->NominalTpcEntrancePoint(),
    &x2 = pair->Track2()->Track()->NominalTpcEntrancePoint();

  const double dist = (x1 - x2).Mag();

  share.Fill(tQinv, share_fraction);
  qual.Fill(tQinv, share_quality);
  sep.Fill(tQinv, dist);
}

void AliFemtoShareQualityCorrFctn::AddRealPair(AliFemtoPair *pair)
{
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  AddPair(pair, *fShareNumerator, *fQualityNumerator, *fTPCSepNumerator);
}

void AliFemtoShareQualityCorrFctn::AddMixedPair(AliFemtoPair *pair)
{
  if (fPairCut && !fPairCut->Pass(pair)) {
    return;
  }

  AddPair(pair, *fShareDenominator, *fQualityDenominator, *fTPCSepDenominator);
}

void AliFemtoShareQualityCorrFctn::WriteHistos()
{
  // Write out result histograms
  fShareNumerator->Write();
  fShareDenominator->Write();
  fQualityNumerator->Write();
  fQualityDenominator->Write();
  fTPCSepNumerator->Write();
  fTPCSepDenominator->Write();
}

TList* AliFemtoShareQualityCorrFctn::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fShareNumerator);
  tOutputList->Add(fShareDenominator);
  tOutputList->Add(fQualityNumerator);
  tOutputList->Add(fQualityDenominator);
  tOutputList->Add(fTPCSepNumerator);
  tOutputList->Add(fTPCSepDenominator);

  return tOutputList;
}
