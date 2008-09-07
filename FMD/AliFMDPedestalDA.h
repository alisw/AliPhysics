#ifndef ALIFMDPEDESTALDA_H
#define ALIFMDPEDESTALDA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
// This class implements the pedestal detector algorithm (DA) for the FMD.
// It uses 51200 TH1S histograms to store the data for each channel of the FMD.
// The mean and standard deviation of a histogram define the pedestal and 
// the noise for that channel.


#include "AliFMDBaseDA.h"
#include "TH1.h"
#include "TObjArray.h"

class AliFMDPedestalDA: public AliFMDBaseDA {
  
 public:
  AliFMDPedestalDA() ;
  AliFMDPedestalDA(const AliFMDPedestalDA & pedDA) ;
  //  AliFMDPedestalDA& operator = (const AliFMDPedestalDA & pedDA) ; 
  virtual ~AliFMDPedestalDA();
  void Init();
 
 protected:
 
  void AddChannelContainer(TObjArray* sectorArray, UShort_t det, Char_t ring, UShort_t sec, UShort_t strip);
  void FillChannels(AliFMDDigit* digit);
  void Analyse(UShort_t det, Char_t ring, UShort_t sec, UShort_t strip);
  void WriteHeaderToFile();
  void FinishEvent() {}
  void Terminate(TFile* );
 private:
  TH1S* GetChannel(UShort_t det, Char_t ring, UShort_t sec, UShort_t strip, UInt_t sample);
 
  Int_t fCurrentChannel;
  TH1F  fPedSummary;
  TH1F  fNoiseSummary;
  std::ofstream fZSfileFMD1;
  std::ofstream fZSfileFMD2;
  std::ofstream fZSfileFMD3;
  
  ClassDef(AliFMDPedestalDA,0)
    
};
#endif
