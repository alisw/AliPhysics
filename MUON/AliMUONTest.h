/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$
//
// Class AliMUONTest
// -----------------
// Class with functions for testing
//
// Author: Ivana Hrivnacova, IPN Orsay

#ifndef ALI_MUON_TEST_H
#define ALI_MUON_TEST_H

#include <TObject.h>

#include "AliDetector.h"
#include "AliMUONData.h"
#include "AliMUONChamber.h"

class TVector;
class TFile;
class TTree;


class AliMUONTest : public  TObject 
{
  public:
    AliMUONTest(const TString& configMacro);
    AliMUONTest();
    virtual ~AliMUONTest();
   
    // tests
    void DetElemTransforms();
    void PrintPadPositions1();
    void PrintPadPositions2();

  protected:
    AliMUONTest(const AliMUONTest& rhs);
    AliMUONTest& operator = (const AliMUONTest& rhs);

    ClassDef(AliMUONTest,0)  // MUON class for tests
};

#endif //ALI_MUON_TEST_H

