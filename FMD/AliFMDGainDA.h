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
#include "TH1F.h"
#include "TObjArray.h"
// #include "TGraphErrors.h"
class TGraphErrors;
class TH1S;

class AliFMDGainDA: public AliFMDBaseDA 
{
public:
  /** 
   * Constructor 
   * 
   */
  AliFMDGainDA() ;
  /** 
   * Copy constructor 
   * 
   * @param gainDA  Object to copy from
   */  
  AliFMDGainDA(const AliFMDGainDA & gainDA) ;
  /** 
   * Destructor 
   * 
   */
  virtual ~AliFMDGainDA();
  /** 
   * Initialize 
   * 
   */
  void Init();
  /** 
   * Set the maximum pulse size 
   * 
   * @param highPulse Maximum pulse size
   */
  void SetMaxPulse(Int_t highPulse = 256) {fHighPulse = highPulse; }
  /** 
   * Set the number of strips per input channel 
   * 
   * @param nStrips Number of strips per channel
   */
  void SetNumberOfStrips(Int_t nStrips) {fNumberOfStripsPerChip = nStrips;}

protected:
  /** 
   * Make a container for a channel 
   * 
   * @param sectorArray Sectors 
   * @param det         Detector # 
   * @param ring        Ring identifier 
   * @param sec         Sector number
   * @param strip       Strip number
   */ 
  void AddChannelContainer(TObjArray* sectorArray, 
			   UShort_t det, Char_t ring, 
			   UShort_t sec, UShort_t strip);
  /** 
   * Fill channel histogram 
   * 
   * @param digit Digit to fill from
   */
  void FillChannels(AliFMDDigit* digit);
  /** 
   * Analyse the result
   * 
   * @param det   Detector # 
   * @param ring  Ring identifier 
   * @param sec   Sector number
   * @param strip Strip number
   * @param h     Summary histogram with bins for sector and strip
   */
  void Analyse(UShort_t det, Char_t ring, UShort_t sec, UShort_t strip);
  /**
   * Write out the header 
   */
  void WriteHeaderToFile();
  /** 
   * Go to the next sample
   * 
   * @param det         Detector # 
   * @param ring        Ring identifier 
   * @param sec         Sector number
   * @param strip       Strip number
   */
  void UpdatePulseAndADC(UShort_t det, 
			 Char_t ring, 
			 UShort_t sec, 
			 UShort_t strip);
  /** 
   * Reset all 
   * 
   */
  void ResetPulseAndUpdateChannel();
  /** 
   * End of event 
   * 
   */
  void FinishEvent();
  /** 
   * End of job
   * 
   * @param dummy Not used
   */
  void Terminate(TFile* dummy);

private:
  /** 
   * Get the channel histogram
   * 
   * @param det         Detector # 
   * @param ring        Ring identifier 
   * @param sec         Sector number
   * @param strip       Strip number
   * 
   * @return Histogram
   */  
  TH1S* GetChannelHistogram(UShort_t det, Char_t ring, UShort_t sec, 
			    UShort_t strip);
  /** 
   * Get strip graph
   * 
   * @param det         Detector # 
   * @param ring        Ring identifier 
   * @param sec         Sector number
   * @param strip       Strip number
   * 
   * @return Graph
   */
  TGraphErrors* GetChannel(UShort_t det, Char_t ring, 
			   UShort_t sec, UShort_t strip);
  TObjArray fGainArray;          // Array of gains 
  Int_t     fHighPulse;          // Highest pulse
  TArrayS   fEventsPerChannel;   // # of events per pulse step
  TArrayS   fCurrentPulse;       // The current pulse size 
  TArrayS   fCurrentChannel;     // The current strip number
  Int_t     fNumberOfStripsPerChip; // Number of strips
  
  TH1F      fSummaryGains;         // Summary histogram 
  Int_t     fCurrentSummaryStrip;  // Current strip for summary

  void  MakeSummary(UShort_t det, Char_t ring);

  TH2* fGainFMD1i;
  TH2* fGainFMD2i;
  TH2* fGainFMD2o;
  TH2* fGainFMD3i;
  TH2* fGainFMD3o;
  TH2* fChi2FMD1i;
  TH2* fChi2FMD2i;
  TH2* fChi2FMD2o;
  TH2* fChi2FMD3i;
  TH2* fChi2FMD3o;
  
  ClassDef(AliFMDGainDA,0) // Detector algorithm for gain runs

};
#endif
// Local Variables: 
//  mode: C++ 
// End:
