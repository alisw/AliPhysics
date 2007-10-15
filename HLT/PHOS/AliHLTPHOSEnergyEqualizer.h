
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

#ifndef ALIHLTPHOSENERGYEQUALIZER_H
#define ALIHLTPHOSENERGYEQUALIZER_H

class AliHLTPHOSEnergyEqualizer : public AliHLTPHOSBase
{
  
public:
  AliHLTPHOSEnergyEqualizer();
  
  virtual ~AliHLTPHOSEnergyEqualizer();
  
  void SetGlobalHighGainConversionFactor(Float_t factor) { fGlobalHighGainFactor = factor; }
  void SetGlobalLowGainConversionFactor(Float_t factor) { fGlobalLowGainFactor = factor; }
  void SetGlobalConversion() { fGlobalConversion = true; }
  
  void SetDigitContainer(AliHLTPHOSDigitContainerStruct *digConPtr) { fDigitContainerPtr = digConPtr; }
  
  Int_t MakeDigits(AliHLTPHOSRcuCellEnergyDataStruct*);
  
private:
  AliHLTPHOSDigitContainerStruct *fDigitContainerPtr;
  
  Float_t fGlobalHighGainFactor;
  Float_t fGlobalLowGainFactor;
  Bool_t fGlobalConversion;

  ClassImp(AliHLTPHOSEnergyEqualizer, 1);

#endif
