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

#ifndef ALIHLTPHOSSANITYINSPECTOR_H
#define ALIHLTPHOSSANITYINSPECTOR_H

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


/** 
 * @class AliHLTPHOSSanityInspector
 * Sanity inspector for PHOS HLT. It takes raw data as input and checks it for insanity
 * It will then flag it.
 *
 * @ingroup alihlt_phos
 */

//class AliHLTPHOSSanityInspector : public AliHLTPHOSBase
class AliHLTPHOSSanityInspector
{
  
public:

  /** Constructor */
  AliHLTPHOSSanityInspector();
  
  /* Destructor */
  virtual ~AliHLTPHOSSanityInspector();
  
  //  Int_t CheckInsanity(UInt_t* data, Int_t nSamples);
  
  
  // Int_t CheckInsanity(Int_t* data, Int_t nSamples);
  
  
  
  /** Check for insanity */
  Int_t CheckInsanity(const UShort_t* data, const Int_t nSamples) const;
  
  /** Check for and heal insanity */
  Int_t CheckAndHealInsanity(UShort_t* data, Int_t nSamples);  //Not completely reliable

  /** Set the max difference between 2 samples before flagging insanity */
  void SetMaxDifference(Int_t maxDiff) { fMaxDifference = maxDiff; }
    
private:
  /** Copy constructor */  
  AliHLTPHOSSanityInspector(const AliHLTPHOSSanityInspector &);
  /** Assignment */
  AliHLTPHOSSanityInspector & operator = (const AliHLTPHOSSanityInspector &);

  /** The max difference between 2 samples */
  Int_t fMaxDifference;           //COMMENT
 
  ClassDef(AliHLTPHOSSanityInspector, 1);
};

#endif
