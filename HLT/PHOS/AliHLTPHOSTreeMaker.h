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


#ifndef ALIHLTPHOSTREEMAKER_H
#define ALIHLTPHOSTREEMAKER_H

/**
 * Class makes tree of digits of the class AliHLTPHOSDigit 
 * stored in an TClonesArray
 *
 * @file   AliHLTPHOSTreeMaker.h
 * @author Oystein Djuvsland
 * @date
 * @brief  Tree maker  for PHOS HLT
 */

#include "AliHLTPHOSBase.h"
#include "TTree.h"
class AliHLTPHOSDigitContainerDataStruct;

class TClonesArray;
//class TTree;


/** 
 * @class AliHLTPHOSTreeMaker
 * Tree maker for PHOS HLT. Takes AliHLTPHOSDigitContainerDataStruct as input,
 * makes a TClonesArray of objecst of the class AliHLTPHOSDigit, and then makes
 * a tree out of it.
 *
 * @ingroup alihlt_phos
 */
class AliHLTPHOSTreeMaker : public AliHLTPHOSBase
{
public:

  /** Constructor */
  AliHLTPHOSTreeMaker();

  /** Destructor */ 
  ~AliHLTPHOSTreeMaker();

  /** 
   * Make a TClonesArray AliHLTPHOSDigits from a container of digit structs
   * @param digitContainer is the container of digit structs
   * @param nDigits is the number of digits in the container
   * @return the number of digits
   */
  Int_t MakeDigitArray(AliHLTPHOSDigitContainerDataStruct* digitContainer, Int_t nDigits);

  /** Fill the digit tree */
  void FillDigitTree();
  
  /** Reset the digit tree */
  void ResetDigitTree() { fDigitTreePtr->Reset(); }
  
  /** Set the digit tree */
  void SetDigitTree(TTree* tree);
  
  /** 
   * Get the digit tree 
   * @return a pointer to the tree
   */
  TTree* GetDigitTree() { return fDigitTreePtr; }
  
private:

  /** The array of digit object */
  TClonesArray *fDigitArrayPtr;     //COMMENT

  /** The digit tree */
  TTree* fDigitTreePtr;             //COMMENT
 
  ClassDef(AliHLTPHOSTreeMaker, 1);

};


#endif
