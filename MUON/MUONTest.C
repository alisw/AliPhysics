/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id$
//
// Macro for testing detection elements transformations
// and segmentations
// To be run from aliroot:
// .L MUONTest.C
// MUONTest(option, testNumber); > testN.out
//     option     = "./Config.C", 
//                  "default", "FactoryV2", "FactoryV3","FactoryV4"
//     testNumber = 1, 2, 3
//
//  Author: I. Hrivnacova, IPN Orsay

void MUONTest(const TString& option = "./Config.C", 
              Int_t testNumber = 1)
{
  AliMUONTest test(option);
  switch (testNumber) {
    case 1: test.DetElemTransforms();  break;
    case 2: test.ForWhole(kPrintPads); break; 
    case 3: test.ForWhole(kDrawPads);  break; 
    default: ;
  }    
}  
