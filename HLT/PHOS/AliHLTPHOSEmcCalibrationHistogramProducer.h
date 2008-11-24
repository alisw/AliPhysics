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

#ifndef ALIHLTPHOSEMCCALIBRATIONHISTOGRAMPRODUCER_H
#define ALIHLTPHOSEMCCALIBRATIONHISTOGRAMPRODUCER_H


/**
 *
 * @file   AliHLTPHOSEmcCalibrationHistogramProducer.h
 * @author Oystein Djuvsland
 * @date
 * @brief  EmcCalibration histogram producer for PHOS HLT
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


#include "AliHLTPHOSConstants.h"
#include "AliHLTPHOSBase.h"

class TH1F;
class TH2F;
class TH2I;

using namespace PhosHLTConst;

/** 
 * @class AliHLTPHOSEmcCalibrationHistogramProducer
 * 
 *
 * @ingroup alihlt_phos
 */
class AliHLTPHOSEmcCalibrationHistogramProducer: public AliHLTPHOSBase
{
  
 public:

  /** Constructor */  
  AliHLTPHOSEmcCalibrationHistogramProducer();

  /** Destructor */
  virtual ~AliHLTPHOSEmcCalibrationHistogramProducer();

    /** Copy constructor */  
  AliHLTPHOSEmcCalibrationHistogramProducer(const AliHLTPHOSEmcCalibrationHistogramProducer &) : AliHLTPHOSBase()
  {
    //Copy constructor not implemented
  }
  
  /** Assignment */
  AliHLTPHOSEmcCalibrationHistogramProducer & operator = (const AliHLTPHOSEmcCalibrationHistogramProducer)
  {
    //Assignment
    return *this; 
  }

  /** Reset the histograms */
  void Reset();

  /** 
   * Fill energy in channel histogram 
   * @param module is the module number [0 - 4]
   * @param xCol is the x-coordinate (column) [0 - 63]
   * @param zRow is the z-coordinate (row) [0 - 55]
   * @param gain is the gain [PhosHLTConst::HIGH_GAIN / PhosHLTConst::LOW_GAIN]
   */
  void FillChannelEnergy(Int_t module, Int_t xCol, Int_t zRow, Int_t gain);

  /**
   * Fill histograms from a histogram containing energies from one RCU from one RCU
   * @param rcuEnergyHistPtr is a pointer to the energy histogram
   * @param module is the module number
   * @param rcuX is the rcu "x coordinate"
   * @param rcuZ is the rcu "z coordinate"
   */
  void FillRCUEnergies(TH1F* rcuChannelEnergyHistPtr[][N_ZROWS_RCU][N_GAINS], Int_t module,Int_t rcuX, Int_t rcuZ);

  /**
   * Fill bad channel histogram from a histogram containing bad channels from one RCU 
   * @param rcuBadChannelHistPtr is a pointer to the bad channel histogram
   * @param module is the module number
   */
  void FillBadChannels(TH2F* rcuBadChannelHistPtr[], Int_t module);
  
  /** 
   * Get the energy distribution histogram for one channel 
   * @param module is the module number [0 - 4]
   * @param xCol is the x-coordinate (column) [0 - 63]
   * @param zRow is the z-coordinate (row) [0 - 55]
   * @param gain is the gain [PhosHLTConst::HIGH_GAIN / PhosHLTConst::LOW_GAIN]
   * @return a pointer to the histogram (TH1F*)
   */
  TH1F* GetChannelEnergyHistogram(Int_t module, Int_t xCol, Int_t zRow, Int_t gain);

  /**
   * Get mean energy 2D histogram for one module, one gain 
   * @param module is the module number [0 - 4]
   * @param gain is the gain [PhosHLTConst::HIGH_GAIN / PhosHLTConst::LOW_GAIN]
   * @return a pointer to the histogram (TH2F*);
   */
  TH2F* GetMeanEnergyHistogram(Int_t module, Int_t gain);

  /**
   * Get max signal 2D histogram for one module, one gain 
   * @param module is the module number [0 - 4]
   * @param gain is the gain [PhosHLTConst::HIGH_GAIN / PhosHLTConst::LOW_GAIN]
   * @return a pointer to the histogram (TH2F*);
   */
  TH2F* GetMaxSignalHistogram(Int_t module, Int_t gain);
  
  /**
   * Get number of hits 2D histogram for one module, one gain 
   * @param module is the module number [0 - 4]
   * @param gain is the gain [PhosHLTConst::HIGH_GAIN / PhosHLTConst::LOW_GAIN]
   * @return a pointer to the histogram (TH2F*);
   */  
  TH2I* GetNHitsHistogram(Int_t module, Int_t gain);
  
  /**
   * Get live channel 2D histogram for one module, one gain 
   * @param module is the module number [0 - 4]
   * @param gain is the gain [PhosHLTConst::HIGH_GAIN / PhosHLTConst::LOW_GAIN]
   * @return a pointer to the histogram (TH2F*);
   */  
  TH2I* GetLiveChannelHistogram(Int_t module, Int_t gain);

  /**
   * Get bad channel 2D histogram for one module, one gain 
   * @param module is the module number [0 - 4]
   * @param gain is the gain [PhosHLTConst::HIGH_GAIN / PhosHLTConst::LOW_GAIN]
   * @return a pointer to the histogram (TH2F*);
   */  
  TH2I* GetBadChannelHistogram(Int_t module, Int_t gain);
  
 protected:
  
  /** Array of pointers to the energy histograms */
  TH1F* fChannelEnergyHistogramPtr[N_MODULES][N_XCOLUMNS_MOD][N_ZROWS_MOD][N_GAINS];    //! transient

  /** Array of pointers to the mean energy histograms */
  TH2F* fMeanEnergyHistogramPtr[N_MODULES][N_GAINS];                                    //! transient

  /** Array of pointers to max signal histograms */
  TH2F* fMaxSignalHistogramPtr[N_MODULES][N_GAINS];                                     //! transient

  /** Array of pointers to hit count histograms */
  TH2I* fNHitsHistogramPtr;                                                             //! transient
  
  /** Array of pointers to live channel histograms */
  TH2I* fLiveChannelHistogramPtr[N_MODULES][N_GAINS];                                   //! transient
  
  /** Array of pointers to bad channel histograms */
  TH2I* fBadChannelHistogramPtr[N_MODULES][N_GAINS];                                    //! transient

  /** ADC -> energy conversion factors */
  Float_t fEnergyConversionFactors[N_MODULES][N_XCOLUMNS_MOD][N_ZROWS_MOD][N_GAINS];    //COMMENT
  
  // ClassDef(AliHLTPHOSEmcCalibrationHistogramProducer,1);
};

#endif
 
