// -*- C++ -*-
#ifndef ALIADTRENDING_H
#define ALIADTRENDING_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved.
 *
 * See cxx source for full Copyright notice
 */


//
// Class AliADTrending
// ---------------------------
//
//  class used in QA to publish variables evolution versus time in AMORE.
//  These histo are the one which will be looked at by QA Shifter
//


#include <TH1.h>

class TGraph;
class TMultiGraph;

class AliADTrending  : public TH1 {
public:
  AliADTrending();
  AliADTrending(const char* name, const char* title);
  virtual ~AliADTrending();
  AliADTrending(const AliADTrending &trend);

  Double_t * GetTime(){return fTime;};
  Double_t * GetChannel(Int_t i){return fData[i];};
  Double_t  GetLastTime(){return fTime[fNEntries-1];};
  Double_t  GetLastChannel(Int_t i){return fData[i][fNEntries];};
  UInt_t GetNEntries(){return fNEntries;};
  void AddEntry(Double_t * data, UInt_t time);
  void PrintEntry(UInt_t entry);
  virtual void Draw(Option_t  *option="");

private:
  AliADTrending& operator= (const AliADTrending & /*trend*/); // Not implemented
  enum{kDataSize = 500};
  Double_t fData[4][kDataSize];
  Double_t fTime[kDataSize];
  UInt_t fNEntries;
  TMultiGraph *fMultiGraphs;
  TGraph * fGraphs[4];

  ClassDef( AliADTrending, 2 );
};

#endif // ALIADTRENDING_H

