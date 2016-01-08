#ifndef AliITSMFTSDigit_H
#define AliITSMFTSDigit_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

#include <TObject.h>

using std::istream;

class AliITSMFTSDigit: public TObject 
{
 public:
  enum {kBuffSize=10};
  //
  AliITSMFTSDigit();
  AliITSMFTSDigit(Int_t track,Int_t hit,UInt_t chip,UInt_t index,Double_t signal,Int_t roCycle=0);
  AliITSMFTSDigit(UInt_t chip,UInt_t index,Double_t noise,Int_t roCycle=0);
  AliITSMFTSDigit(const AliITSMFTSDigit &source);
  AliITSMFTSDigit& operator=(const AliITSMFTSDigit &source);
  virtual ~AliITSMFTSDigit() {}
  Double_t GetSignal(Int_t i)    const {return ( (i>=0&&i<kBuffSize) ? fSignal[i] : 0.0);}
  Double_t GetSignal()           const {return fTsignal;}
  Double_t GetSignalAfterElect() const {return fSignalAfterElect;}
  Double_t GetSumSignal()        const {return fTsignal+fNoise;}
  Double_t GetNoise()            const {return fNoise;}
  Int_t GetNsignals()            const {return kBuffSize;}
  void AddSignal(Int_t track,Int_t hit,Double_t signal);
  void AddSignalAfterElect(Double_t signal) {fSignalAfterElect += signal;}
  void AddNoise(Double_t noise)  {fNoise += noise;}
  void SetNoise(Double_t noise)  {fNoise = noise;}
  void SetROCycle(Int_t cl)      {fROCycle=cl;}
  //
  Int_t GetTrack(Int_t i)        const {return ((i>=0&&i<kBuffSize) ? fTrack[i] : 0);}
  Int_t GetHit(Int_t i)          const {return ((i>=0&&i<kBuffSize) ? fHits[i] : 0);}
  Int_t GetChip()              const {return fChip;}
  Int_t GetNTracks()             const {return fNTracks;}
  Int_t GetROCycle()             const {return fROCycle;}
  //
  void Add(const AliITSMFTSDigit *pl);
  void AddTo(Int_t fileIndex, const AliITSMFTSDigit *pl);
  void ShiftIndices(Int_t fileIndex);
  void Print(Option_t *option="")                 const;
  Int_t Read(const char *name)                          {return TObject::Read(name);}
  //
  virtual Bool_t IsSortable()                     const {return kTRUE;}
  virtual Bool_t IsEqual(const TObject* obj)      const {return GetUniqueID()==obj->GetUniqueID();}
  virtual Int_t  Compare(const TObject* obj)      const;
  //
  static Int_t GetBuffSize() {return kBuffSize;};
  //
 private:
  UShort_t fChip;            // chip number
  UShort_t fNTracks;           // number of tracks contributing
  Int_t    fROCycle;           // readOut cycle
  Int_t    fTrack[kBuffSize];  // track Number
  Int_t    fHits[kBuffSize];   // hit number
  Float_t  fSignal[kBuffSize]; // Signals
  Float_t  fTsignal;           // Total signal (no noise)
  Float_t  fNoise;             // Total noise, coupling, ...
  Float_t  fSignalAfterElect;  // Signal after electronics
  //
  ClassDef(AliITSMFTSDigit,2) // Item list of signals and track numbers
};	



#endif
