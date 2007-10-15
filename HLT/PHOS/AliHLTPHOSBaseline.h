
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

#include "TObject.h"
#include "AliHLTPHOSConstants.h"

using namespace PhosHLTConst;

/**
 @author Øystein Djuvsland <oystein.djuvsland@gmail.com>
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

  /**
   * 
   * @param baseline 
   */
  void SetBaseline(Float_t baseline) { fBaseline = baseline; }
  void SetX(Int_t x) { fX = x; }
  void SetZ(Int_t z) { fZ = z; }
  void SetGain(Int_t gain) { fGain = gain; }
  void SetEntries(Int_t entries) { fEntries = 0; }

  Float_t GetBaseline() { return fBaseline; }  
  Int_t GetX() { return fX; }
  Int_t GetZ() { return fZ; }
  Int_t GetGain() { return fGain; }
  Int_t GetEntries() { return fEntries; }

private:

  Float_t fBaseline;
  Int_t fX;
  Int_t fZ;
  Int_t fGain;
  Int_t fEntries;
  
  ClassDef(AliHLTPHOSBaseline, 1);
  
};

#endif
