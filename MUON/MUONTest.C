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
// MUONTest(testNumber); > testN.out
//     testNumber = 1  ...    
//     testNumber = 2  ...   
//     testNumber = 3  ...   
//
//  Author: I. Hrivnacova, IPN Orsay

void MUONTest(Int_t testNumber)
{
  gAlice->Init("./Config_MUON_test.C");
  cout << "Init done " << endl;

  AliMUONTest test("./Config_MUON_test.C");
  switch (testNumber) {
    case 1: test.DetElemTransforms();       break;
    case 2: test.PrintPadPositions1(); break; 
    case 3: test.PrintPadPositions2(); break;
    default: ;
  }    
}  
