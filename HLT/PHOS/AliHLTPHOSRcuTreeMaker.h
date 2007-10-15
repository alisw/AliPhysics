 
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


#ifndef ALIHLTPHOSRCUTREEMAKER_H
#define ALIHLTPHOSRCUTREEMAKER_H

#include "AliHLTPHOSBase.h"
#include "TTree.h"
class AliHLTPHOSRcuDigitContainerDataStruct;

class TClonesArray;
//class TTree;

class AliHLTPHOSRcuTreeMaker : public AliHLTPHOSBase
{
public:

  AliHLTPHOSRcuTreeMaker();
  ~AliHLTPHOSRcuTreeMaker();

  Int_t MakeDigitArray(AliHLTPHOSRcuDigitContainerDataStruct* digitContainer, Int_t nDigits);

  void FillDigitTree();
  
  void ResetDigitTree() { fDigitTreePtr->Reset(); }
  
  void SetDigitTree(TTree* tree);
  
  TTree* GetDigitTree() { return fDigitTreePtr; }
  
private:
  TClonesArray *fDigitArrayPtr;
  TTree* fDigitTreePtr;
 
  ClassDef(AliHLTPHOSRcuTreeMaker, 1);

};


#endif
