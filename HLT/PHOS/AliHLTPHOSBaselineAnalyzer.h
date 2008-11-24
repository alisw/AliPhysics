//-*- Mode: C++ -*-
// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIHLTPHOSBASELINEANALYZER_H
#define ALIHLTPHOSBASELINEANALYZER_H

/**
 * Measures the baselines and calculates relevant values
 *
 * @file   AliHLTPHOSBaselineAnalyzer.h
 * @author Oystein Djuvsland
 * @date
 * @brief  Baseline analyzer for PHOS HLT
 */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSConstants.h"
#include "AliHLTPHOSBase.h"

class AliHLTPHOSRcuCellEnergyDataStruct;
class AliHLTPHOSValidCellDataStruct;
class AliHLTPHOSSanityInspector;
class TTree;
class TClonesArray;
class TH1F;
class TH2F;

using namespace PhosHLTConst;


/**
 * @class AliHLTPHOSBaselineAnalyzer
 * Measures the baseline in a baseline run. It also calculates 
 * RMS of the signal noise in each channel. It takes as input raw signals,
 * and outputs an object of AliHLTPHOSBaseline class for each channel. 
 * It also creates several histograms for each channel.
 *
 * @ingroup alihlt_phos
 */

class AliHLTPHOSBaselineAnalyzer : public AliHLTPHOSBase
{
  
public:

  /** Constructor */
  AliHLTPHOSBaselineAnalyzer();
  
  /** Destructor */ 
  virtual ~AliHLTPHOSBaselineAnalyzer();
  
  /** Copy constructor */
  AliHLTPHOSBaselineAnalyzer(const AliHLTPHOSBaselineAnalyzer &) : 
    AliHLTPHOSBase(),
    fSampleStart(5),
    fMaxCrazyDifference(0),	
    fMaxSignal(0),  
    fChannelCount(0),
    fTreePtr(0),
    fBaselineArrayPtr(0),
    fRMSHistogramPtr(0),
    fRMSMapHighGainHistogramPtr(0),
    fRMSMapLowGainHistogramPtr(0),
    fFixedRMSHistogramPtr(0),
    fFixedRMSMapHighGainHistogramPtr(0),
    fFixedRMSMapLowGainHistogramPtr(0),
    fSanityInspector(0)
  {
    //Copy constructor not implemented
  }
  
  /** Assignment */
  AliHLTPHOSBaselineAnalyzer & operator = (const AliHLTPHOSBaselineAnalyzer)
    {
      //Assignment
      return *this; 
    }
    
  /** 
   * Calculate baselines for an RCU
   * @param rcuData is a pointer to energy and timing data from an RCU
   */
  void CalculateRcuBaselines(AliHLTPHOSRcuCellEnergyDataStruct* rcuData);

  /**
   * Calculate the baseline of channels
   * @param cellData is a pointer to a valid channel
   * @param rawDataPtr is a pointer to the raw data
   * @param xOff is the offset in the x position given by the RCU number
   * @param zOff is the offset in the z position given by the RCU number
   * @return the baseline
   */
  Float_t CalculateChannelBaseline(AliHLTPHOSValidCellDataStruct *cellData, Int_t* rawDataPtr, Int_t xOff, Int_t zOff);

  /** 
   * Calculate the accumulated baseline of a channel
   * @param x is the x coord of the channel
   * @param z is the z coord of the channel
   * @param gain is the gain of the channel
   * @param baseline is the most recent measured baseline
   */
  void CalculateAccumulatedChannelBaseline(Int_t x, Int_t z, Int_t gain, Float_t baseline);

  /**
   * Calculate the channels baseline RMS
   */
  void CalculateChannelsBaselineRMS();


  // void CalculateAccumulatedBaselines();

  /** Set the ROOT objects to be filled */
  void SetRootObjects(TTree *tree, TClonesArray *array);

  /** Fill the tree */
  void FillTree();

  /** Write the accumulated baselines */
  void WriteAccumulatedBaselines(const Char_t* filename);

  /** Write the channel histograms */
  void WriteChannelHistograms(const Char_t* filename);

  /** Write the RMS histograms */
  void WriteRMSHistogram(const Char_t* filename);

  /** Reset the baseline array */
  void ResetBaselines();

  /** Reset the channel count */
  void ResetChannelCount();

  /** Reset the accumulated baseline array */ 
  void ResetAccumulatedBaselines();

  /** Set the number of samples to be used for the baseline calculation */
  void SetNumberOfSamples(Int_t nSamples) { fNSamples = nSamples; }

  /** Set the max difference between 2 samples before flagging crazy  */
  void SetMaxCrazyDifference(Int_t diff);

  /** Set the maximum signal before flagging that a real signal is in the samples */
  void SetMaxSignal(Int_t max) { fMaxSignal = max; }
  
 // void SetChannelsHistograms(TH1F *channelLowGainHistArray[N_XCOLUMNS_MOD][N_ZROWS_MOD], TH1F *channelLowGainHistArray[N_XCOLUMNS_MOD][N_ZROWS_MOD]);
 
     
private:

  /** At which time index to start the measuring of the baseline */
  Int_t fSampleStart;  // Shutting up rule checker

  /** Not used anymore */
  Int_t fMaxCrazyDifference;  // Shutting up rule checker

  /** Maximum signal level, used to not count channels with signals for baseline measurement */
  Int_t fMaxSignal;  // Shutting up rule checker

  /** Count of channels, not used anymore */
  Int_t fChannelCount;   // Shutting up rule checker

  /** Array containing the baselines of all channels */
  Float_t fBaselines[N_XCOLUMNS_MOD][N_ZROWS_MOD][N_GAINS];      // Shutting up rule checker

  /** Array containing the accumulated baselines for all channels */
  Float_t fAccumulatedBaselines[N_XCOLUMNS_MOD][N_ZROWS_MOD][N_GAINS][2];    // Shutting up rule checker

  /** Pointer to a tree containing TClonesArray of AliHLTPHOSBaseline objects */ 
  TTree *fTreePtr;                                      //! transient

  /** Pointer to an array of AliHLTPHOSBaseline objects */
  TClonesArray *fBaselineArrayPtr;                      //! transient


 // TH1F *fChannelLowGainHistogramsPtr[N_XCOLUMNS_MOD][N_ZROWS_MOD]; //comment
 // TH1F *fChannelHighGainHistogramsPtr[N_XCOLUMNS_MOD][N_ZROWS_MOD]; //comment


  /** Pointer to an array of histograms containg the baselines */ 
  TH1F *fChannelHistogramsPtr[N_XCOLUMNS_MOD][N_ZROWS_MOD][N_GAINS];      //! transient

  /** Pointer to an array of histograms containg the "fixed" baselines */
  TH1F *fFixedChannelHistogramsPtr[N_XCOLUMNS_MOD][N_ZROWS_MOD][N_GAINS]; //! transient

  /** Pointer to a histograms containing the RMS values for all channels */ 
  TH1F *fRMSHistogramPtr;                             //! transient

  /** Pointer to a 2D histogram of the high gain RMS */ 
  TH2F *fRMSMapHighGainHistogramPtr;                         //! transient

  /** Pointer to a 2D histogram of the low gain RMS */ 
  TH2F *fRMSMapLowGainHistogramPtr;                          //! transient

  /** Pointer to a histogram containing the "fixed" RMS values for all channels */ 
  TH1F *fFixedRMSHistogramPtr;                               //! transient

  /** Pointer to a 2D histogram of the "fixed" high gain channels */ 
  TH2F *fFixedRMSMapHighGainHistogramPtr;                    //! transient

  /** Pointer to a 2D histogram of the "fixed" low gain channels */
  TH2F *fFixedRMSMapLowGainHistogramPtr;                     //! transient

  /** Pointer to a sanity inspector */
  AliHLTPHOSSanityInspector *fSanityInspector;               //! transient
 
  ClassDef(AliHLTPHOSBaselineAnalyzer, 1);
};

#endif
 
