#ifndef ALIANALYSISMUMUSPECTRA_H
#define ALIANALYSISMUMUSPECTRA_H

///
/// AliAnalysisMuMuSpectra : a spectra is a binning (AliAnalysisMuMuBinning)
/// and the results per bin (AliAnalysisMuMuResult), and is mergeable, so
/// it can be put into an AliMergeableCollection
///
/// author : Laurent Aphecetche (Subatech)
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

  void AdoptResult(const AliAnalysisMuMuBinning::Range& bin, AliAnalysisMuMuResult* result);

  Bool_t IsEmpty() const;
  
  Long64_t Merge(TCollection* list);

  TH1* Plot(const char* what="NofJpsi", const char* subresult="", Bool_t divideByBinWidth=kTRUE) const;

  void Print(Option_t* opt="") const;

  TObjArray* Bins() const { return fBins; }
  
  AliAnalysisMuMuBinning* Binning() const { return fBinning; }
  
  Bool_t Correct(const AliAnalysisMuMuSpectra& accEff, const char* particle, const char* subResultName="");
  
  AliAnalysisMuMuResult* GetResultForBin(const AliAnalysisMuMuBinning::Range& bin) const;
  
  Bool_t HasValue(const char* what="NofJpsi") const;
  
private:
  AliAnalysisMuMuBinning* fBinning; // internal binning
  TObjArray* fBins; // the results (bin by bin)
  
  ClassDef(AliAnalysisMuMuSpectra,1) // class to hold spectra (with its associated binning and errors)
};

#endif
