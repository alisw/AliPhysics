
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

#include "AliHLTPHOSBase.h"
#include "AliHLTPHOSConstants.h"


class AliHLTPHOSDigit;
class TClonesArray;
class TTree;
class AliHLTPHOSValidCellDataStruct;
class AliHLTPHOSRcuCellEnergyDataStruct;
class AliHLTPHOSDigitContainerDataStruct;
class AliHLTPHOSDigitDataStruct;
       
using namespace PhosHLTConst;

class AliHLTPHOSDigitMaker : public AliHLTPHOSBase
{
public:
  AliHLTPHOSDigitMaker();
  virtual ~AliHLTPHOSDigitMaker();
 
  // void SetValidCellData(AliHLTPHOSValidCellDataStruct *data) { fCellDataPtr = data; }
  //  void SetDigitContainerStruct(AliHLTPHOSDigitContainerStruct *container) 
  //{ fDigitContainerStructPtr = container; }
                                
  void SetDigitContainerStruct(AliHLTPHOSDigitContainerDataStruct *container) 
  { fDigitContainerStructPtr = container; }
  
  void SetDigitArray(TClonesArray *array) { fDigitArrayPtr = array; }
  void ResetDigitCount() { fDigitCount = 0; }
  void SetDigitThreshold(Int_t threshold) { fDigitThreshold = threshold; }
  void SetNrPresamples(Int_t n) { fNrPresamples = n; }
  Int_t MakeDigits(AliHLTPHOSRcuCellEnergyDataStruct*);
  void Reset();

private:

  AliHLTPHOSValidCellDataStruct *fCellDataPtr; //comment
  AliHLTPHOSDigitContainerDataStruct *fDigitContainerStructPtr; //comment
  TClonesArray *fDigitArrayPtr; //comment
  AliHLTPHOSDigit *fDigitPtr; //comment
  //AliHLTPHOSDigitDataStruct *fDigitStructPtr; //comment
  AliHLTPHOSDigitDataStruct *fDigitStructPtr; //comment
  Int_t fDigitCount;  //comment
  Int_t fNrPresamples; //comment

  Float_t fDigitThreshold; //comment

  //  ClassDef(AliHLTPHOSDigitMaker, 1); 
};


#endif
 
