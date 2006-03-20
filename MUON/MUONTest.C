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
//  MUONTest(option); > test.out
//     option     = "./Config.C", 
//                  "default", "FactoryV2", "FactoryV3","FactoryV4"
//
//
//  Author: I. Hrivnacova, IPN Orsay

void MUONTest(const TString& option = "./Config.C")
{
  // Load tests library
  gSystem->Load("libMUONtests");

  AliMUONTest test(option);
  
  // Print pads for all DEs
  // test.PrintPadsForAll();
  
  // Print pads for first chamber, first cathod
  // test.PrintPadsForSegmentation(0, 0);

  // Print pads for detElem 100, first cathod
  // test.PrintPadsForDetElement(100, 0);

  
  // Print pads for all DEs
  test.DrawPadsForAll();
  
  // Draw pads for first chamber, first cathod
  // test.DrawForSegmentation(0, 0);

  // Draw pads for detElem 100, first cathod
  // test.DrawForDetElement(100, 0);
}  
