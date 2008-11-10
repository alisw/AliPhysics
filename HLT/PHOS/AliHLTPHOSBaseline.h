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

#ifndef ALIHLTPHOSBASELINE_H
#define ALIHLTPHOSBASELINE_H

/**
 * Baseline class for PHOS HLT
 *
 * @file   AliHLTPHOSBaseline.h
 * @author Oystein Djuvsland
 * @date
 * @brief  Baseline class for PHOS HLT
 */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "TObject.h"
#include "AliHLTPHOSConstants.h"

using namespace PhosHLTConst;

/**
 * @class AliHLTPHOSBaseline
 * Baseline class for PHOS HLT, used for writing to ROOT files
 *
 * @ingroup alihlt_phos
 */
class AliHLTPHOSBaseline : public TObject
{
  
public: 
  /**
   * Default constructor
   */
  AliHLTPHOSBaseline();
  /**
   * Default destructor
   */
  virtual ~AliHLTPHOSBaseline();

  /** Copy constructor */  
  AliHLTPHOSBaseline(const AliHLTPHOSBaseline &) : 
    TObject(),
    fBaseline(-1),
    fX(-1),
    fZ(-1),
    fGain(-1),
    fEntries(0)
  {
    //Copy constructor not implemented
  }
  
  /** Assignment */
  AliHLTPHOSBaseline & operator = (const AliHLTPHOSBaseline)
    
    {
      //Assignment
      return *this; 
    }
  

  /**
   * Set the baseline
   * @param baseline is the baseline
   */
  void SetBaseline(Float_t baseline) { fBaseline = baseline; }

  /** 
   * Set the x coordinate
   * @param x is the x coordinate
   */
  void SetX(Int_t x) { fX = x; }

  /** 
   * Set the z coordinate
   * @param z is the z coordinate
   */
  void SetZ(Int_t z) { fZ = z; }

  /** 
   * Set the gain 
   * @param gain is the gain 
   */
  void SetGain(Int_t gain) { fGain = gain; }


  /** 
   * Set the number of entries 
   * @param entries is the entries
   */
  void SetEntries(Int_t entries) { fEntries = entries; }

  /**
   * Get the baseline
   * @return the baseline
   */
  Float_t GetBaseline() { return fBaseline; }  

  /**
   * Get the x
   * @return the x
   */
  Int_t GetX() { return fX; }

  /**
   * Get the z
   * @return the z
   */
  Int_t GetZ() { return fZ; }

  /**
   * Get the gain
   * @return the gain
   */
  Int_t GetGain() { return fGain; }

  /**
   * Get the number of entries
   * @return the number of entries
   */
  Int_t GetEntries() { return fEntries; }

private:

  /** The baseline */
  Float_t fBaseline;      //COMMENT

  /** The x coordinate */
  Int_t fX;               //COMMENT

  /** The z coordinate */
  Int_t fZ;               //COMMENT

  /** The gain */
  Int_t fGain;            //COMMENT

  /** The number of entries */ 
  Int_t fEntries;         //COMMENT
  
  ClassDef(AliHLTPHOSBaseline, 1);
  
};

#endif
