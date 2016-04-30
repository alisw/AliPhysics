#ifndef ALIFFTSMOOTHER_H
#define ALIFFTSMOOTHER_H

class TGraph;
class TH1;
class TTreeSRedirector;

namespace AliFFTsmoother
{
  TGraph  * ReplaceOutlierFrequencies(TH1 *hinput, Int_t firstBin, Int_t lastBin, Double_t outlierCut, Int_t skipFreq, TTreeSRedirector * pcstream);
  TGraph  * SmoothFrequencies(TH1 *hinput, Int_t firstBin, Int_t lastBin, Double_t outlierCut, Int_t skipFreq, TTreeSRedirector * pcstream);

};


#endif

