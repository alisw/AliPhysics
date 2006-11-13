#ifndef AliTPCRawHistograms_H
#define AliTPCRawHistograms_H

/* $Id$ */

//
// This class contains a number of histograms for diagnostics of a TPC
// read out chamber from the raw data
//

#include <TNamed.h>

class TH3F;
class TH1F;
class TCanvas;
class TTree;
class TNtuple;

class AliTPCRawStream;

class AliTPCRawHistograms : public TNamed
{
public:
  AliTPCRawHistograms();
  AliTPCRawHistograms(Int_t detector, const Char_t* comment="", Int_t timeStart=-1, Int_t timeStop=-1);
  
  AliTPCRawHistograms(const AliTPCRawHistograms& c);
  virtual ~AliTPCRawHistograms();
  AliTPCRawHistograms& operator=(const AliTPCRawHistograms& corrMatrix);

  virtual Long64_t Merge(TCollection* list);

  virtual void SaveHistograms();

  void FillDigit(AliTPCRawStream* rawStream, Int_t time=-1);

  TCanvas* DrawHistograms(const Char_t* opt="");

protected:
  Int_t       fTimeStart;               // begin time of run(s)
  Int_t       fTimeStop;                // end time of runs(s)

  TH3F*       fhDigits;                  // cluster of all digits above threshold
  TH1F*       fhSignal;                  // signal distribution
  
  TNtuple*    fDigitTree;                 // row:pad:time:signal
  
  ClassDef(AliTPCRawHistograms,1)
};

#endif

