//-*- Mode: C++ -*-
// $Id: AliHLTPHOSSanityInspector.h 35107 2009-09-30 01:45:06Z phille $

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

#ifndef ALIHLTCALOSANITYINSPECTOR_H
#define ALIHLTCALOSANITYINSPECTOR_H

/**
 * Class checks data for insanity
 * for use in HLT, but can also be used offline
 *
 * @file   AliHLTPHOSSanityInspector.h
 * @author Oystein Djuvsland
 * @date
 * @brief  Sanity inspector for PHOS HLT
 */

//#include "AliHLTPHOSBase.h"
#include "Rtypes.h"


class AliHLTAltroData;

/** 
 * @class AliHLTPHOSSanityInspector
 * Sanity inspector for PHOS HLT. It takes raw data as input and checks it for insanity
 * It will then flag it.
 *
 * @ingroup alihlt_phos
 */

//class AliHLTCaloSanityInspector : public AliHLTCaloBase
class AliHLTCaloSanityInspector
{
  
public:

  /** Constructor */
  AliHLTCaloSanityInspector();
  
  /* Destructor */
  virtual ~AliHLTCaloSanityInspector();
  
  //  Int_t CheckInsanity(UInt_t* data, Int_t nSamples);
  
  
  // Int_t CheckInsanity(Int_t* data, Int_t nSamples);
  
  
  
  /** Check for insanity */
  Int_t CheckInsanity(const UInt_t* data, const Int_t nSamples) const;
  
  /** Check for and heal insanity */
  Int_t CheckAndHealInsanity(UInt_t* data, Int_t nSamples);  //Not completely reliable

  /** Check for and heal insanity */
  Int_t CheckAndHealInsanity(Int_t* data, Int_t nSamples);  //Not completely reliable

  /** Set the max difference between 2 samples before flagging insanity */
  void SetMaxDifference(Int_t maxDiff) { fMaxDifference = maxDiff; }
    
private:
  /** Copy constructor */  
  AliHLTCaloSanityInspector(const AliHLTCaloSanityInspector &);
  /** Assignment */
  AliHLTCaloSanityInspector & operator = (const AliHLTCaloSanityInspector &);

  /** The max difference between 2 samples */
  Int_t fMaxDifference;           //COMMENT
 
  ClassDef(AliHLTCaloSanityInspector, 1);
};

#endif
