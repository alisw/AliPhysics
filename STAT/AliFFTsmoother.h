#ifndef ALIFFTSMOOTHER_H
#define ALIFFTSMOOTHER_H

class TGraph;
class TH1;
class TTreeSRedirector;

namespace AliFFTsmoother
{
  TGraph  * ReplaceOutlierFrequenciesMedian(TH1 *hinput, Double_t outlierCut, Int_t medianRange, Float_t smoothSigma,  Int_t lowBand, TTreeSRedirector * pcstream);
  TGraph  * BandFilter(TH1 *hinput, Int_t lowBand, Int_t highBand, TTreeSRedirector * pcstream);
};


#endif

