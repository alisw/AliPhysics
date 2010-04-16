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
//  Pre-Trigger Control-Box bottom class                                  //
//                                                                        //
//  Authors: F. Reidt (Felix.Reidt@cern.ch)                               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


#include "TROOT.h"

#include "AliRun.h"
#include "AliRunLoader.h"

#include "AliTRDptrgLUT.h"
#include "AliTRDptrgParam.h"
#include "AliTRDptrgCBAC.h"
#include "AliTRDptrgTLMU.h"
#include "AliTRDptrgCBB.h"
ClassImp(AliTRDptrgCBB)

AliTRDptrgCBB::AliTRDptrgCBB(AliRunLoader *rl) 
  : TObject(),
  fRunLoader(rl),
  fParam(0),
  fOperatingMode(kDigits),
  fCBA(0),
  fCBC(0),
  fTLMU(0),
  fLUTArray(0),
  fPTmasks(0x0)
{
  // default ctor
  AliError("default ctor - usage not recommended\n");
}


AliTRDptrgCBB::AliTRDptrgCBB(AliRunLoader *rl, AliTRDptrgParam* param, 
                             AliTRDptrgOperatingMode_t operatingMode)
  : TObject(),
  fRunLoader(rl),
  fParam(param),
  fOperatingMode(operatingMode),
  fCBA(0),
  fCBC(0),
  fTLMU(0),
  fLUTArray(0),
  fPTmasks(0x0) 
{
  // recommended ctor
  this->fCBA = new AliTRDptrgCBAC(rl, kA, operatingMode, param);
  this->fCBC = new AliTRDptrgCBAC(rl, kC, operatingMode, param);
  this->fTLMU = new AliTRDptrgTLMU(rl, param, operatingMode);

  this->LoadParams();
}

AliTRDptrgCBB::~AliTRDptrgCBB() 
{
  // destructor

  if (this->fCBA != 0x0) {
    delete this->fCBA;
    this->fCBA = 0x0;
  }

  if (this->fCBC != 0x0) {
    delete this->fCBC;
    this->fCBC = 0x0;
  }

  if (this->fTLMU != 0x0) {
    delete this->fTLMU;
    this->fTLMU = 0x0;
  }

  this->fLUTArray.Delete();
}


//______________________________________________________________________________
Bool_t AliTRDptrgCBB::LoadParams() 
{
  // load configuration parameters

  if (this->fParam != 0x0) {
    // read AliTRDptrgParam
    
    // get LUTs
    // 0
    AliTRDptrgLUT* LUT = new AliTRDptrgLUT();
    LUT->InitTable(12, 12, this->fParam->GetCBLUT(0, 0), kFALSE);
    // get CB-B_0 and do not copy lut content
    this->fLUTArray.AddLast(LUT);
 
    // 1
    LUT = new AliTRDptrgLUT();
    LUT->InitTable(12, 12, this->fParam->GetCBLUT(0, 1), kFALSE);
    // get CB-B_0 and do not copy lut content
    this->fLUTArray.AddLast(LUT);

    // masks
    this->fPTmasks = this->fParam->GetPTmasks();
  }
  else {
    // load default parameters 
    // initialize LUTsoutputWidth=<value optimized out>
    AliTRDptrgLUT* LUT = new AliTRDptrgLUT();
    this->fLUTArray.AddLast(LUT);
    LUT = new AliTRDptrgLUT(); // this->fRunLoader
    this->fLUTArray.AddLast(LUT);
    // the following lines are only needed for test reasons
    LUT = dynamic_cast<AliTRDptrgLUT*>(this->fLUTArray.At(0));
    Int_t* initData = new Int_t[4096]; // 2^12
    for (Int_t i = 0; i < 4096; i++ ) {
      initData[i] = i;
    }
    LUT->InitTable(12, 12, initData, kTRUE); // make a copy
    LUT = dynamic_cast<AliTRDptrgLUT*>(this->fLUTArray.At(1));
    for (Int_t i = 4096; i >= 0; i--) {
      initData[4096 - i] = i;  // inverse ramp
    }
    LUT->InitTable(12, 12, initData, kTRUE); // make a copy 
  
    AliTRDptrgParam::AliTRDptrgPTmasks* masks = 
      new AliTRDptrgParam::AliTRDptrgPTmasks();  
    masks->fLUTs[0] = kTRUE;
    masks->fLUTs[1] = kTRUE;
    this->fPTmasks = masks;
  }  
  return false;
}

//______________________________________________________________________________
Int_t* AliTRDptrgCBB::Simulate()
{ 
  // Simulate the CBB behavior of event
  //
  // returns array containing:
  // 0: array element count
  // 1..count-2: LUT results
  // count-1: pretrigger decision

  Int_t nLUTs = this->fLUTArray.GetEntries();

  Int_t inputVector = 0x0;
  // initialize partResults
  Int_t** partResults = 0x0;  
  partResults = new Int_t* [3]; // CB-A, CB-C, TLMU
 
  // get partResults
  partResults[0] = this->fCBA->Simulate(); // CB-A
  partResults[1] = this->fCBC->Simulate(); // CB-C
  partResults[2] = this->fTLMU->Simulate(); // TLMU
  
  
  // combine partResults and create inputVectors  
  Int_t mask = 0x1;
  for (Int_t i = 0; i < 3 ; i++) {
    for (Int_t j = 1; j <= partResults[i][0]; j++) {
      if (partResults[i][j] > 0) {
        inputVector |= mask; // Add bit to the  inputVector
      }
      mask <<= 1;
    }
  }
  
  AliDebug(5, Form("Inputvectors: 0x%x", inputVector));
    
  // perform look up
  Int_t* result = new Int_t[nLUTs + 2]; // generate new return array
  result[0] = nLUTs + 1; // storage array length in the first array value
  for (Int_t iLUT = 0; iLUT < nLUTs; iLUT++) { 
    // process the return value for each LUT and store the result in the array
    result[iLUT + 1] = 
      dynamic_cast<AliTRDptrgLUT*>(this->fLUTArray[iLUT])->LookUp(inputVector);
    AliDebug(4, Form("CBB result[%d] = 0x%x\n",(iLUT + 1),result[iLUT + 1])); 
  }
  
  // evaluate PT decision
  // stored in result[nLUTs + 1]
  result[nLUTs + 1] = 0;

  for (Int_t i = 0; i < 2; i++) {
    // CB-A
    if (this->fPTmasks->fCBA[i] && partResults[0][i + 1]) {
      result[nLUTs + 1]++;
    }
    // CB-C
    if (this->fPTmasks->fCBC[i] && partResults[1][i + 1]) {
      result[nLUTs + 1]++;
    }
    // CB-B (own LUTs)
    if (this->fPTmasks->fLUTs[i] && result[i + 1]) {
      result[nLUTs + 1]++;
    }       
  }
  // TLMU
  for (Int_t i = 0; i < 8; i++) {
    if (this->fPTmasks->fTLMU[i] && partResults[2][i + 1]) {
      result[nLUTs + 1]++;
    }
  }
  AliDebug(4, Form("CBB TRD Wake up result = %d", result[nLUTs + 1]));
  return result;
}

//______________________________________________________________________________
Bool_t AliTRDptrgCBB::GetPT() {
  // evaluates the pre trigger decision

  Int_t* LUTresults = this->Simulate();
  if (LUTresults[(LUTresults[0] - 1)]) {
    delete[] LUTresults;
    return kTRUE;
  }
  else {
    delete[] LUTresults;
    return kFALSE;
  }
}
