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

//#include "AliHLTPHOSBase.h"
#include "AliHLTPHOSConstants.h"
#include "AliHLTPHOSDigitDataStruct.h"
#include "AliHLTPHOSChannelDataStruct.h"
#include "AliHLTDataTypes.h"
#include "AliHLTLogging.h"

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

//class AliHLTPHOSDigitMaker : public AliHLTPHOSBase
class AliHLTPHOSDigitMaker : public AliHLTLogging
{
public:

  /** Constructor */
  AliHLTPHOSDigitMaker();

  /** Destructor */
  virtual ~AliHLTPHOSDigitMaker();

  /** Copy constructor */  
  AliHLTPHOSDigitMaker(const AliHLTPHOSDigitMaker &) : 
    AliHLTLogging(),
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
  void SetDigitHeaderPtr(AliHLTPHOSDigitHeaderStruct *digitHeaderPtr) 
  { 
    fDigitHeaderPtr = digitHeaderPtr;
    fDigitStructPtr = reinterpret_cast<AliHLTPHOSDigitDataStruct*>(reinterpret_cast<Long_t>(digitHeaderPtr) + sizeof(AliHLTPHOSDigitHeaderStruct));
  }

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

  /**
   * Reset the digit maker
   */
  void Reset() { fDigitCount = 0; }

  /**
   * Sort the digits and make internal links between them
   */
  void SortDigits();

  /** 
   * Compare two digits, used during the sorting
   */
  static Int_t CompareDigits(const void *dig0, const void *dig);

  
private:

  /**
   * Add a new digit
   * @param channelData is the channel data
   * @param coordinates is the coordinates of the channel, including gain and module
   * @return true if the digit is added correctly, false if out of buffer
   */
  bool AddDigit(AliHLTPHOSChannelDataStruct* channelData, UShort_t* channelCoordinates, Float_t* localCoordinates)
  {
    if(fAvailableSize < sizeof(AliHLTPHOSDigitDataStruct))
      {
	HLTError("Output buffer is full, stopping digit making");
	return false;
      }

    fAvailableSize -= sizeof(AliHLTPHOSDigitDataStruct);

    fDigitStructPtr->fX = channelCoordinates[0];
    fDigitStructPtr->fZ = channelCoordinates[1];

    fDigitStructPtr->fID = fDigitStructPtr->fZ * NXCOLUMNSRCU + fDigitStructPtr->fX;

    fDigitStructPtr->fLocX = localCoordinates[0];
    fDigitStructPtr->fLocZ = localCoordinates[1];

    if(channelCoordinates[2] == HIGHGAIN)
      {
	fDigitStructPtr->fEnergy = channelData->fEnergy*fHighGainFactors[channelCoordinates[0]][channelCoordinates[1]];
	//	HLTDebug("HG channel (x = %d, z = %d) with amplitude: %f --> Digit with energy: %f \n", channelCoordinates[0], channelCoordinates[1], channelData->fEnergy, fDigitStructPtr->fEnergy);
      }
    else
      {
	fDigitStructPtr->fEnergy = channelData->fEnergy*fLowGainFactors[channelCoordinates[0]][channelCoordinates[1]];
	//	HLTDebug("LG channel (x = %d, z = %d) with amplitude: %f --> Digit with energy: %f\n", channelCoordinates[0], channelCoordinates[1], channelData->fEnergy, fDigitStructPtr->fEnergy); 
      }
    fDigitStructPtr->fTime = channelData->fTime * 0.0000001; //TODO
    fDigitStructPtr->fCrazyness = channelData->fCrazyness;
    fDigitStructPtr->fModule = channelCoordinates[3];

    fDigitPtrArray[fDigitCount] = fDigitStructPtr;
    fDigitCount++;
    fDigitStructPtr++;
  }

  /** Pointer to shared memory interface */
  AliHLTPHOSSharedMemoryInterfacev2* fShmPtr;                    //! transient

  /** Pointer to the AliHLTPHOSDigitDataStruct */
  AliHLTPHOSDigitDataStruct *fDigitStructPtr;                    //! transient

  /** Pointer to the AliHLTPHOSDigitDataStruct */
  AliHLTPHOSDigitHeaderStruct *fDigitHeaderPtr;                  //! transient

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

  /** Array of digit pointers */
  AliHLTPHOSDigitDataStruct **fDigitPtrArray;               //COMMENT
  
  /** The available size of the output buffer */
  AliHLTUInt32_t fAvailableSize;                            //COMMENT

  ClassDef(AliHLTPHOSDigitMaker, 0); 

};


#endif
 
