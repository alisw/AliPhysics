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

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Pre-Trigger Control-Box A or C for simulation                         //
//                                                                        //
//  Authors: F. Reidt (Felix.Reidt@cern.ch)                               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <stdio.h>

#include "AliRun.h"
#include "AliRunLoader.h"

#include "AliLog.h"

#include "AliTRDptrgParam.h"
#include "AliTRDptrgFEB.h"
#include "AliTRDptrgLUT.h"

#include "AliTRDptrgCBAC.h"

ClassImp(AliTRDptrgCBAC)

//______________________________________________________________________________
AliTRDptrgCBAC::AliTRDptrgCBAC(AliRunLoader *rl) 
  : TObject(),
  fRunLoader(rl),
  fLUTArray(),
  fFEBArray(),
  fPosition(kUnknown),
  fOperatingMode(kDigits),
  fParam(0x0)
{
  // ctor
  AliError("default ctor - usage not recommended");
}

//______________________________________________________________________________
AliTRDptrgCBAC::AliTRDptrgCBAC(AliRunLoader *rl, 
                               AliTRDptrgFEBPosition_t position,
                               AliTRDptrgOperatingMode_t operatingMode,
                               AliTRDptrgParam *param) 
  : TObject(),
  fRunLoader(rl),
  fLUTArray(),
  fFEBArray(),
  fPosition(position),
  fOperatingMode(operatingMode),
  fParam(param)
{
  // ctor  
  this->LoadParams(); // load parameters
 
  // T0
  AliTRDptrgFEB *FEB = new AliTRDptrgFEB(this->fRunLoader, kTZERO, 
                                         this->fOperatingMode, this->fPosition,
                                         0, this->fParam);
  this->fFEBArray.AddLast(FEB);

  // V0-1
  FEB = new AliTRDptrgFEB(this->fRunLoader, kVZERO, this->fOperatingMode, 
                          this->fPosition, 1, this->fParam);
  this->fFEBArray.AddLast(FEB);

  // V0-2
  FEB = new AliTRDptrgFEB(this->fRunLoader, kVZERO, this->fOperatingMode, 
                          this->fPosition, 2, this->fParam);
  this->fFEBArray.AddLast(FEB);

  // V0-3
  FEB = new AliTRDptrgFEB(this->fRunLoader, kVZERO, this->fOperatingMode, 
                          this->fPosition, 3, this->fParam);
  this->fFEBArray.AddLast(FEB);

  // V0-4
  FEB = new AliTRDptrgFEB(this->fRunLoader, kVZERO, this->fOperatingMode, 
                          this->fPosition, 4, this->fParam);
  this->fFEBArray.AddLast(FEB);

}

//______________________________________________________________________________
Bool_t AliTRDptrgCBAC::LoadParams() 
{
  // load configuration parameters

  if (this->fParam != 0x0) {
    // read AliTRDptrgParam

    // get LUTs
    AliTRDptrgLUT* LUT = new AliTRDptrgLUT();
    // 0
    LUT = new AliTRDptrgLUT();
    LUT->InitTable(10, 10, this->fParam->GetCBLUT(this->fPosition, 0), kFALSE); 
    // do not copy table data 
    this->fLUTArray.AddLast(LUT);
    // 1
    LUT = new AliTRDptrgLUT();
    LUT->InitTable(10, 10, this->fParam->GetCBLUT(this->fPosition, 1), kFALSE); 
    // do not copy table data 
    this->fLUTArray.AddLast(LUT);
  }
  else {
    // load default parameters 
    AliTRDptrgLUT* LUT = new AliTRDptrgLUT();
    this->fLUTArray.AddLast(LUT);
    LUT = new AliTRDptrgLUT();
    this->fLUTArray.AddLast(LUT);
    // the following lines are only needed for test reasons
    LUT = dynamic_cast<AliTRDptrgLUT*>(this->fLUTArray.At(0));
    Int_t* initData = new Int_t[1024]; // 2^10
    for (Int_t i = 0; i < 1024; i++ ) {
      initData[i] = i;
    }
    LUT->InitTable(10, 10, initData, kTRUE); // copy initData
    LUT = dynamic_cast<AliTRDptrgLUT*>(this->fLUTArray.At(1));
    for (Int_t i = 1023; i >= 0; i--) {
      initData[31 - i] = i;  // inverse ramp
    }
    LUT->InitTable(10, 10, initData, kTRUE); // copy initData
  }  
  return false;
}

//______________________________________________________________________________
AliTRDptrgCBAC::~AliTRDptrgCBAC() 
{
  // destructor

  this->fLUTArray.Delete();
  this->fFEBArray.Delete();
}

//______________________________________________________________________________
Int_t* AliTRDptrgCBAC::Simulate()
{ 
  // Simulate the CBAC behavior of event
  Int_t nFEBs = this->fFEBArray.GetEntries();
  Int_t nLUTs = this->fLUTArray.GetEntries();

  Int_t inputVector = 0x0;

  Int_t** partResults = 0x0;  
  partResults = new Int_t* [nFEBs];


  for (Int_t iFEB = 0; iFEB < nFEBs; iFEB++) {
    partResults[iFEB] = 
      dynamic_cast<AliTRDptrgFEB*>(this->fFEBArray.At(iFEB))->Simulate();
  }
  
  // combine partResults and create inputVector  
  Int_t iBit = 0;
  Int_t mask = 0x1;
  for (Int_t iFEB = 0; iFEB < nFEBs ; iFEB++) {
    for (Int_t j = 1; j <= partResults[iFEB][0]; j++) {
      if ((iBit >
           dynamic_cast<AliTRDptrgLUT*>(this->fLUTArray[0])->GetInputWidth()) 
          || (iBit >
           dynamic_cast<AliTRDptrgLUT*>(this->fLUTArray[0])->GetInputWidth())) {
        AliError("FEB result width does not match CB-A/C input with!");
      }
      iBit++;
      if (partResults[iFEB][j] > 0) {
        inputVector |= mask; // Add bit to the corresponding inputVector
        mask <<= 1;
      } 
    }
  }

  AliDebug(5, Form("Inputvector: 0x%x", inputVector));
    
  // perform look up
  Int_t* result = new Int_t[nLUTs + 1]; // generate new return array
  result[0] = nLUTs; // storage array length in the first array value
  for (Int_t iLUT = 0; iLUT < nLUTs; iLUT++) { 
    // process the return value for each LUT and store the result in the array
    result[iLUT + 1] = 
      dynamic_cast<AliTRDptrgLUT*>(
        this->fLUTArray[iLUT])->LookUp(inputVector);
    AliDebug(4, Form("CBAC result[%d] = 0x%x",(iLUT + 1),result[iLUT + 1])); 
  }

  return result;
}



