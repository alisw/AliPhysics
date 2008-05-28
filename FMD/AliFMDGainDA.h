#ifndef ALIFMDGAINDA_H
#define ALIFMDGAINDA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
// This class implements the Gain detector algorithm (DA) for the FMD. 
// The gain is the response of the VA chips to a known pulse and has to be 
// calculated strip-by-strip, pulse-by-pulse. 


#include "AliFMDBaseDA.h"
#include "TH1.h"
#include "TObjArray.h"
#include "TGraphErrors.h"

class AliFMDGainDA: public AliFMDBaseDA {
  
 public:
  
  AliFMDGainDA() ;
  AliFMDGainDA(const AliFMDGainDA & gainDA) ;
  //  AliFMDGainDA& operator = (const AliFMDGainDA & gainDA) ; 
  virtual ~AliFMDGainDA();
  void Init();
  // void SetPulseSize(Int_t pulseSize = 32) {fPulseSize = pulseSize; }
  void SetMaxPulse(Int_t highPulse = 256) {fHighPulse = highPulse; }
  //  void SetPulseLength(Int_t pulseLength = 100) {fPulseLength = pulseLength; }
  void SetNumberOfStrips(Int_t nStrips) {fNumberOfStripsPerChip = nStrips;}

 protected:
 
  void AddChannelContainer(TObjArray* sectorArray, UShort_t det, Char_t ring, UShort_t sec, UShort_t strip);
  void FillChannels(AliFMDDigit* digit);
  void Analyse(UShort_t det, Char_t ring, UShort_t sec, UShort_t strip);
  void WriteHeaderToFile();
  void UpdatePulseAndADC(UShort_t det, Char_t ring, UShort_t sec, UShort_t strip);
  void ResetPulseAndUpdateChannel();
  void FinishEvent();
  
 private:
  
  TH1S* GetChannelHistogram(UShort_t det, Char_t ring, UShort_t sec, UShort_t strip);
  TGraphErrors* GetChannel(UShort_t det, Char_t ring, UShort_t sec, UShort_t strip);
  TObjArray fGainArray;
  //  Int_t fPulseSize;
  Int_t fHighPulse;
  //Int_t fPulseLength;
  TArrayS fEventsPerChannel;
  TArrayS fCurrentPulse;
  TArrayS fCurrentChannel;
  Int_t fNumberOfStripsPerChip;
  
  ClassDef(AliFMDGainDA,0)

};
#endif
