
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
#ifndef ALIHLTPHOSDIGITMAKER_H
#define ALIHLTPHOSDIGITMAKER_H

/**
 * Class makes digits from information from raw data
 *
 * @file   AliHLTPHOSDigitMaker.h
 * @author Oystein Djuvsland
 * @date
 * @brief  Digit maker for PHOS HLT
 */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSBase.h"
#include "AliHLTPHOSConstants.h"

/**
 * @class AliHLTPHOSDigitMaker
 * Digit maker for PHOS HLT. Takes input from AliHLTPHOSRawAnalyzer, and 
 * outputs an AliHLTPHOSDigitDataStruct container, or a TClonesArray of
 * AliHLTPHOSDigit. Can do software zero suppression
 * @ingroup alihlt_phos
 */

class TClonesArray;
class TTree;
class TH2F;
class AliHLTPHOSValidCellDataStruct;
class AliHLTPHOSRcuCellEnergyDataStruct;
class AliHLTPHOSDigitContainerDataStruct;
class AliHLTPHOSDigitDataStruct;
class AliHLTPHOSSharedMemoryInterfacev2; // added by PTH
class AliHLTPHOSChannelDataHeaderStruct;
       
using namespace PhosHLTConst;

class AliHLTPHOSDigitMaker : public AliHLTPHOSBase
{
public:

  /** Constructor */
  AliHLTPHOSDigitMaker();

  /** Destructor */
  virtual ~AliHLTPHOSDigitMaker();

  /** Copy constructor */  
  AliHLTPHOSDigitMaker(const AliHLTPHOSDigitMaker &) : 
    AliHLTPHOSBase(),
    fCellDataPtr(0),
    fShmPtr(0),
    fDigitContainerStructPtr(0),
    fDigitArrayPtr(0),
    fDigitStructPtr(0),
    fDigitCount(0),
    fOrdered(true)
  {
    //Copy constructor not implemented
  }
  
  /** Assignment */
  AliHLTPHOSDigitMaker & operator = (const AliHLTPHOSDigitMaker)
  {
    //Assignment
    return *this; 
  }

 
  // void SetValidCellData(AliHLTPHOSValidCellDataStruct *data) { fCellDataPtr = data; }
  //  void SetDigitContainerStruct(AliHLTPHOSDigitContainerStruct *container) 
  //{ fDigitContainerStructPtr = container; }
   
  /**
   * Sets the AliHLTPHOSDigitDataStruct container 
   * @param container the digit container
   */
  void SetDigitContainerStruct(AliHLTPHOSDigitContainerDataStruct *container) 
  { fDigitContainerStructPtr = container; }

  /** 
   * Sets the TClonesArray of AliHLTPHOSDigit 
   * @param array the array
   */
  void SetDigitArray(TClonesArray *array) { fDigitArrayPtr = array; }

  /** 
   * Sets the digit thresholds
   * @param filepath is the path to a file containing RMS
   *                 histograms for both gains
   * @param nSigmas is the number of sigmas to put the threshold
   */
  void SetDigitThresholds(const char* filepath, Int_t nSigmas);

  /** 
   * Sets the digit thresholds
   * @param thresholdHG is the high gain threshold for digitisation
   * @param thresholdLG is the low gain threshold for digitisation
   */
  void SetDigitThresholds(const float thresholdHG, float thresholdLG);
  

  /**
   * Sets the number of pre samples
   * @param n the number of pre samples
   */
  //void SetNrPresamples(Int_t n) { fNrPresamples = n; }

  /**
   * Set the global high gain conversion factory 
   * @param factor is the conversion factor
   */
  void SetGlobalHighGainFactor(Float_t factor);

  /**
   * Set the global low gain conversion factory 
   * @param factor is the conversion factor
   */
  void SetGlobalLowGainFactor(Float_t factor);

  /**
   * Make the digits for one event.
   * @param channelDataHeader is the data header from the AliHLTPHOSRawAnalyzer
   * @return the number of digits found
   */
  Int_t MakeDigits(AliHLTPHOSChannelDataHeaderStruct* channelDataHeader);

  /**
   * Set the mask for dead channels
   * @param badChannelHGHist is a pointer to a high gain bad channel histogram
   * @param badChannelLGHist is a pointer to a low gain bad channel histogram
   * @param qCut is the cut 
   */
  void SetBadChannelMask(TH2F* badChannelHGHist, TH2F* badChannelLGHist, Float_t qCut);

  /**
   * Set ordering of gains or not
   */
  void SetOrdered(bool val) { fOrdered = val; }

  /** Reset the digit maker */
  void Reset();

private:

  /** Pointer to valid cell list */
  AliHLTPHOSValidCellDataStruct *fCellDataPtr;                   //! transient

  /** Pointer to shared memory interface */
  AliHLTPHOSSharedMemoryInterfacev2* fShmPtr;                    //! transient

  /** Pointer to the digit container */
  AliHLTPHOSDigitContainerDataStruct *fDigitContainerStructPtr;  //! transient

  /** Pointer to the digit TClonesArray */
  TClonesArray *fDigitArrayPtr;                                  //! transient

  /** Pointer to a AliHLTPHOSDigitDataStruct */
  AliHLTPHOSDigitDataStruct *fDigitStructPtr;                    //! transient

  /** Digit count */
  Int_t fDigitCount;                                             //COMMENT

  /** Is the gains ordered? */
  bool fOrdered;                                                 //COMMENT

  /** Array containing the energies of all RCU channels */
  Float_t fEnergyArray[N_XCOLUMNS_RCU][N_ZROWS_RCU][N_GAINS];    //COMMENT

  /** High gain energy conversion factors */
  Float_t fHighGainFactors[N_XCOLUMNS_MOD][N_ZROWS_MOD];         //COMMENT

  /** Low gain energy conversion factors */
  Float_t fLowGainFactors[N_XCOLUMNS_MOD][N_ZROWS_MOD];          //COMMENT

   /** Threshold for making digit ( zero suppression threshold) */
  Float_t fDigitThresholds[N_XCOLUMNS_MOD][N_ZROWS_MOD][N_GAINS]; //COMMENT

  /** Bad channel mask */
  Float_t fBadChannelMask[N_XCOLUMNS_MOD][N_ZROWS_MOD][N_GAINS]; //COMMENT

  ClassDef(AliHLTPHOSDigitMaker, 1); 
};


#endif
 
