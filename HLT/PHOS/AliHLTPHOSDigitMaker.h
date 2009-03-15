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
#include "AliHLTPHOSDigitDataStruct.h"
#include "AliHLTPHOSChannelDataStruct.h"

/**
 * @class AliHLTPHOSDigitMaker
 * Digit maker for PHOS HLT. Takes input from AliHLTPHOSRawAnalyzer, and 
 * outputs a block of AliHLTPHOSDigitDataStruct container
 * @ingroup alihlt_phos
 */

class TH2F;
class AliHLTPHOSSharedMemoryInterfacev2; // added by PTH
class AliHLTPHOSChannelDataHeaderStruct;
class AliHLTPHOSMapper;
       
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
    fShmPtr(0),
    fDigitStructPtr(0),
    fDigitCount(0),
    fOrdered(true),
    fMapperPtr(0)
  {
    //Copy constructor not implemented
  }
  
  /** Assignment */
  AliHLTPHOSDigitMaker & operator = (const AliHLTPHOSDigitMaker)
  {
    //Assignment
    return *this; 
  }

  /**
   * Sets the pointer to the output
   * @param the output pointer
   */
  void SetDigitDataPtr(AliHLTPHOSDigitDataStruct *digitDataPtr) 
  { fDigitStructPtr = digitDataPtr; }

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
  Int_t MakeDigits(AliHLTPHOSChannelDataHeaderStruct* channelDataHeader, AliHLTUInt32_t availableSize);


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
  
  void Reset() { fDigitCount = 0; }

private:

  /**
   * Add a new digit
   * @param channelData is the channel data
   * @param coordinates is the coordinates of the channel, including gain and module
   */
  void AddDigit(AliHLTPHOSChannelDataStruct* channelData, UShort_t* coordinates) 
  {
    fDigitStructPtr->fX = coordinates[0];
    fDigitStructPtr->fZ = coordinates[1];
    if(coordinates[2] == HIGHGAIN)
      fDigitStructPtr->fAmplitude = channelData->fEnergy*fHighGainFactors[coordinates[0]][coordinates[1]];
    else
      fDigitStructPtr->fAmplitude = channelData->fEnergy*fLowGainFactors[coordinates[0]][coordinates[1]];
    fDigitStructPtr->fTime = channelData->fTime * 0.0000001; //CRAP
    fDigitStructPtr->fCrazyness = channelData->fCrazyness;
    fDigitStructPtr->fModule = coordinates[3];
    fDigitStructPtr++;
  }

  /** Pointer to shared memory interface */
  AliHLTPHOSSharedMemoryInterfacev2* fShmPtr;                    //! transient

  /** Pointer to the AliHLTPHOSDigitDataStruct */
  AliHLTPHOSDigitDataStruct *fDigitStructPtr;                    //! transient

  /** Digit count */
  Int_t fDigitCount;                                             //COMMENT

  /** Are the gains ordered? */
  bool fOrdered;                                                 //COMMENT

  /** Mapper */
  AliHLTPHOSMapper* fMapperPtr;                                  //COMMENT

  /** High gain energy conversion factors */
  Float_t fHighGainFactors[NXCOLUMNSMOD][NZROWSMOD];         //COMMENT

  /** Low gain energy conversion factors */
  Float_t fLowGainFactors[NXCOLUMNSMOD][NZROWSMOD];          //COMMENT

  /** Bad channel mask */
  Float_t fBadChannelMask[NXCOLUMNSMOD][NZROWSMOD][NGAINS]; //COMMENT

  ClassDef(AliHLTPHOSDigitMaker, 1); 

};


#endif
 
