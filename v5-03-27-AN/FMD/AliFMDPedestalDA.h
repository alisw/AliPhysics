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
class TH2;

class AliFMDPedestalDA: public AliFMDBaseDA 
{
public:
  /** 
   * Constructor.
   * 
   */  
  AliFMDPedestalDA();
  /** 
   * Copy constructor 
   * 
   * @param pedDA Object to copy from
   */  
  AliFMDPedestalDA(const AliFMDPedestalDA & pedDA);
  /** 
   * Assignment operator 
   * 
   * @param pedDA Object to assign from
   */  
  AliFMDPedestalDA& operator=(const AliFMDPedestalDA&) { return *this; }
  /** 
   * Destructor
   * 
   */
  virtual ~AliFMDPedestalDA();
  /** 
   * Initialiser
   * 
   */  
  void Init();
 
protected:
  /** 
   * Add a channel to the containers. 
   * 
   * @param sectorArray  Array of sectors
   * @param det          Detector 
   * @param ring         Ring
   * @param sec          Sector 
   * @param strip        Strip
   */
  void AddChannelContainer(TObjArray* sectorArray, UShort_t det, 
			   Char_t ring, UShort_t sec, UShort_t strip);
  /** 
   * Fill ADC values from a digit into the corresponding histogram.
   * 
   * @param digit Digit to fill ADC values for.
   */
  void FillChannels(AliFMDDigit* digit);
  /** 
   * Analyse a strip.  That is, compute the mean and spread of the ADC
   * spectra for all strips.  Also output on files the values. 
   * 
   * @param det   Detector
   * @param ring  Ring
   * @param sec   Sector 
   * @param strip Strip.
   * @param h     Summary histogram with bins for sector and strip
   */
  void Analyse(UShort_t det, Char_t ring, UShort_t sec, UShort_t strip);
  /** 
   * Write headers to files. 
   * 
   */  
  void WriteHeaderToFile();
  /** 
   * Called at the end of an event.
   * 
   */  
  void FinishEvent() {}
  /** 
   * Called at the end of a job.  Fills in missing time-bins and
   * closes output files  
   * 
   */  
  void Terminate(TFile* );
private:
  /** 
   * Get the histogram corresponding to a strip sample.
   * 
   * @param det    Detector
   * @param ring   Ring
   * @param sec    Sector
   * @param strip  Strip
   * @param sample Sample
   * 
   * @return ADC spectra of a strip.
   */
  TH1S* GetChannel(UShort_t det, Char_t ring, UShort_t sec, 
		   UShort_t strip, UInt_t sample);
  /** 
   * Calculate the hardware index
   * 
   * @param ddl    DDL number
   * @param board  Board number
   * @param altro  ALTRO number
   * @param chan   Channel number
   * 
   * @return Index into hardware cache.
   */
  Int_t HWIndex(UShort_t ddl, UShort_t board, UShort_t altro, 
		UShort_t chan) const;
  void FillinTimebins(std::ofstream& out, UShort_t ddl);
  /** Current strip */ 
  Int_t fCurrentChannel;                           //The current channel
  /** Pedestal summary */ 
  TH1F  fPedSummary;                               //Summary of pedestals
  /** Noise summary */
  TH1F  fNoiseSummary;                             //Summary of noises
  /** Output file for zero-suppression for FMD1 */
  std::ofstream fZSfileFMD1;                       //Stream for ZS FMD1
  /** Output file for zero-suppression for FMD2 */
  std::ofstream fZSfileFMD2;                       //Stream for ZS FMD2
  /** Output file for zero-suppression for FMD3 */
  std::ofstream fZSfileFMD3;                       //Stream for ZS FMD3 
  /** The minimum timebin seen for all channels */
  TArrayS fMinTimebin;                             //minimum timebin
  /** The maximum timebin seen for all channels */
  TArrayS fMaxTimebin;                             //maximum timebin
  
  void  MakeSummary(UShort_t det, Char_t ring);

  TH2* fSummaryFMD1i;                              //Summary of FMD1
  TH2* fSummaryFMD2i;                              //Summary of FMD2I 
  TH2* fSummaryFMD2o;                              //Summary of FMD2O
  TH2* fSummaryFMD3i;                              //Summary of FMD3I
  TH2* fSummaryFMD3o;                              //Summary of FMD3O
  
  ClassDef(AliFMDPedestalDA,0)
};

inline Int_t
AliFMDPedestalDA::HWIndex(UShort_t ddl, UShort_t b, 
			  UShort_t a, UShort_t c) const
{
  // Save some array entries 
  UShort_t lb = (b > 1 ? b-16+2 : b);
  const Int_t kNDDL     = 3;
  const Int_t kNBoard   = 4;
  const Int_t kNAltro   = 3;
  const Int_t kNChannel = 16;
  Int_t idx =  c + kNChannel * (a + kNAltro * (lb + kNBoard * ddl));
  if (idx > kNDDL * kNBoard * kNAltro * kNChannel) return -1;
  return idx;
}

#endif
//
// Local Variables:
//  mode: C++
// End:
//
