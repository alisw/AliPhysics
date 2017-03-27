#ifndef ALIANALYSISMUMUSPECTRA_H
#define ALIANALYSISMUMUSPECTRA_H

/// @ingroup pwg_muondep_mumu
/// @class AliAnalysisMuMuSpectra
/// @brief Encapsulates results from MuMu analysis
/// @details AliAnalysisMuMuSpectra is a spectra with a binning (AliAnalysisMuMuBinning)
//  and the results per bin (AliAnalysisMuMuResult). It is mergeable, so
//  it can be put into an AliMergeableCollection. It can also be converted into
 // histograms.
/// @author Laurent Aphecetche (Subatech)
//

#include "TNamed.h"

#ifndef ALIANALYSISMUMUBINNING_H
#  include "AliAnalysisMuMuBinning.h"
#endif

class AliAnalysisMuMuResult;
class TCollection;
class TH1;
class TObjArray;

class AliAnalysisMuMuSpectra : public TNamed
{
public:
  AliAnalysisMuMuSpectra(const char* name="", const char* title="");
  AliAnalysisMuMuSpectra(const AliAnalysisMuMuSpectra& rhs);
  AliAnalysisMuMuSpectra& operator=(const AliAnalysisMuMuSpectra& rhs);

  virtual ~AliAnalysisMuMuSpectra();

  Bool_t AdoptResult(const AliAnalysisMuMuBinning::Range& bin, AliAnalysisMuMuResult* result);

  Bool_t IsEmpty() const;

  Long64_t Merge(TCollection* list);

  Long64_t Merge(AliAnalysisMuMuSpectra* spectraToAdd);

  TH1* Plot(const char* what="NofJpsi", const char* subresult="", Bool_t divideByBinWidth=kTRUE) const;

  void Print(Option_t* opt="") const;

  TObjArray* BinContentArray() const { return fBins; }

  AliAnalysisMuMuBinning* Binning() const { return fBinning; }

  Bool_t Correct(const AliAnalysisMuMuSpectra& accEff, const char* particle, const char* subResultName="");

  AliAnalysisMuMuResult* GetResultForBin(const AliAnalysisMuMuBinning::Range& bin) const;

  AliAnalysisMuMuResult* GetResultForBin(const char* binName) const;

  Bool_t HasValue(const char* what="NofJpsi") const;

  void Scale(Double_t value);

  void SetWeight(Double_t w);

  Double_t Weight() const { return fWeight; }

private:
  AliAnalysisMuMuBinning* fBinning; // internal binning
  TObjArray* fBins; // the results (bin by bin)
  Double_t fWeight; // weight of this spectra (assumed to be a normalized weight)

  ClassDef(AliAnalysisMuMuSpectra,2) // class to hold spectra (with its associated binning and errors)
};

#endif
