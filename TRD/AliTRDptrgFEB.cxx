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
//                                                                        
//  Pre-Trigger simulation                                                
//                                                                        
//  Authors: F. Reidt (Felix.Reidt@cern.ch)                               
//		  
//  This class is used to simulate the front end box behavior of the 
//  pretrigger system. Digits of T0 and V0 are used as input. A threshold
//  discrimination, masking and first processing with look up tables  is 
//  done during the simulation process
//                                                                  
////////////////////////////////////////////////////////////////////////////

#include <TClonesArray.h> 

#include "AliRunLoader.h"
#include "AliLog.h"

#include "../VZERO/AliVZEROdigit.h" 
#include "../T0/AliT0digit.h"

#include "AliTRDptrgParam.h"
#include "AliTRDptrgLUT.h"
#include "AliTRDptrgFEB.h"

ClassImp(AliTRDptrgFEB)

//______________________________________________________________________________
AliTRDptrgFEB::AliTRDptrgFEB(AliRunLoader *rl) 
  : TObject(),
  fRunLoader(rl),
  fParam(0),
  fLUTArray(0),
  fType(kUndefined),
  fOperatingMode(kDigits),
  fInputChannelCount(0),
  fPosition(kUnknown),  
  fID(0),
  fThreshold(0)
{
  // default constructor
  AliError("default ctor - not recommended");
}

//______________________________________________________________________________
AliTRDptrgFEB::AliTRDptrgFEB(AliRunLoader *rl, AliTRDptrgFEBType_t febType, 
                             AliTRDptrgOperatingMode_t operatingMode,
                             AliTRDptrgFEBPosition_t position, Int_t id, 
                             AliTRDptrgParam *param)
  : TObject(),
  fRunLoader(rl),
  fParam(param),
  fLUTArray(0),
  fType(febType),
  fOperatingMode(operatingMode),
  fInputChannelCount(0),
  fPosition(position),
  fID(id),
  fThreshold(0x0) 
{
  // prefered constructor
  
  this->LoadParams(); // load configuration parameters

}

//______________________________________________________________________________
AliTRDptrgFEB::~AliTRDptrgFEB() 
{
  // destructor
  if (this->fParam == 0x0) {  
    if (this->fThreshold != 0x0) {
      delete[] this->fThreshold;
      this->fThreshold = 0x0;
   }
  }
  // delete LUTArray
  this->fLUTArray.Delete();
}

//______________________________________________________________________________
Int_t AliTRDptrgFEB::LoadDigits()
{
  // loads T0 or V0 digits and discriminates them automatically
 
  if (this->fType == kVZERO) {
    // load V0's digits --------------------------------------------------------
    
    // get V0 run loader
    AliLoader* loader = this->fRunLoader->GetLoader( "VZEROLoader" );

    if (!loader) {
      AliError("Can not get VZERO loader");
      return -1;
    }
    loader->LoadDigits("READ");
    TTree* vzeroDigitsTree = loader->TreeD();

    if (!vzeroDigitsTree) {
      AliError("Can not get the VZERO digit tree");
      return -1;
    }
    
		
    TClonesArray* vzeroDigits = NULL;
    TBranch* digitBranch = vzeroDigitsTree->GetBranch("VZERODigit");
    digitBranch->SetAddress(&vzeroDigits);
    vzeroDigitsTree->GetEvent(0);	
    
  
    Int_t nDigits = vzeroDigits->GetEntriesFast(); // get digit count
		
    AliDebug(5, Form("Found a whole of %d digits", nDigits));    
		
    Int_t inputVector = 0x0; // Vector which is feed into the LUT
		
    for (Int_t iDigit=0; iDigit<nDigits; iDigit++) {
      // loop over all digits
      AliDebug(5, "Looping over digit");
      AliVZEROdigit* digit = (AliVZEROdigit*)vzeroDigits->At(iDigit);
			      
      Int_t pmNumber   = digit->PMNumber();
      Int_t board   = pmNumber / 8;
      Int_t channel = pmNumber % 8;
      Int_t position = pmNumber / 32 + 1;

      // check whether the digits belongs to the current FEB, otherwise omit it
      if ((position == this->fPosition) && (board == this->fID)) {
        AliDebug(5, "Found an digit corresponding to the current FEB");
        Float_t value = digit->ADC();
        AliDebug(5, Form("ADC value: %f\n", value));
        Int_t channelBitMask = 0x01; 
        // channel0 => 0x01; channel1=> 0x02;  2^(channel number)
        channelBitMask <<= channel;
        if (value >= this->fThreshold[channel]) {
          inputVector |= channelBitMask;
          AliDebug(5,
            Form("Threshold exceeded in channel %d, new inputVector 0x%x", 
                 channel, inputVector));
        }
      }
    }

    AliDebug(5, Form("inputVector: 0x%x", inputVector));
    
    return inputVector;
  }
  else if (this->fType == kTZERO) {
    // load T0's digits --------------------------------------------------------
    AliLoader * fT0Loader = this->fRunLoader->GetLoader("T0Loader");
    //   AliT0digit *fDigits; 
    fT0Loader->LoadDigits("READ");
    // Creating T0 data container

    TTree* treeD = fT0Loader->TreeD();
    if (!treeD) {
      AliError("no digits tree");
      return -1;
    }
    AliT0digit* digits = new AliT0digit();
    TBranch *brDigits = treeD->GetBranch("T0");

    if (brDigits) {
      brDigits->SetAddress(&digits);
    }
    else {
      AliError("Branch T0 DIGIT not found");
      return -1;
    }     
    brDigits->GetEntry(0);		

    TArrayI qtc0(24); // Array must have 24 entries!
    TArrayI qtc1(24); // Array must have 24 entries!
    
    digits->GetQT0(qtc0); // baseline (reference level)
    digits->GetQT1(qtc1); // measurement value

    Int_t inputVector = 0x0; // vector to be fed into the look up table
    
    // PMT Positions
    // C: 0  to 11
    // A: 12 to 23
    // positions according to AliT0Digitizer.cxx Revision 37491
    Int_t nStart = 0;
    if (this->fPosition == kC) { // C
      nStart = 0;
    }
    else if (this->fPosition == kA) { // A
      nStart = 12;
    }

    Int_t channelBitMask = 0x01;
    for (Int_t i = 0 + nStart; i < nStart + 12; i++) {
      //Int_t channelBitMask = 0x01;
      AliDebug(5, Form("channel: %d", i));

      Int_t value = qtc1[i] - qtc0[i]; // calculate correct measurement value

      if (value > (Int_t)this->fThreshold[i - nStart]) {
        inputVector |= channelBitMask;    // Add bit
          
        AliDebug(5, Form("Threshold exceeded in channel %d,", i));
        AliDebug(5, Form("new inputVector 0x%x", inputVector));       
        AliDebug(5, Form("channelBitMask 0x%x", channelBitMask));
      }
      channelBitMask <<= 1; // go on to the next channel
    }
    
    delete digits;
    return inputVector;
  }
  return -1;
}

//______________________________________________________________________________
Int_t AliTRDptrgFEB::LoadAndProcessHits()
{
  // loads TO or VO hits and converts them to digits optimized for ptrg  
  // afterwards the digits will be discriminated
  AliError("LoadAndProcessHits() - not yet implemented!\n");
  if (this->fType == kVZERO) {		
    return 0;
  }
  else if (this->fType == kTZERO) {
    return 0;
  }
  return -1;
}

//______________________________________________________________________________
Bool_t AliTRDptrgFEB::LoadParams()
{
  // Load Parameters

  if (this->fParam == 0x0) {
    AliWarning("No paramater object specified - start loading defaults\n");
    if (this->fType == kVZERO) {		
      // initialize threshold
      this->fThreshold = new UInt_t[8]; 
      for (Int_t i = 0; i < 8; i++) {
        this->fThreshold[i] = 10; 
      }
      // initialize LUTsoutputWidth=<value optimized out>
      AliTRDptrgLUT* lut = new AliTRDptrgLUT();
      this->fLUTArray.AddLast(lut);
      lut = new AliTRDptrgLUT(); 
      this->fLUTArray.AddLast(lut);
			// the following lines are only needed for test reasons
      lut = dynamic_cast<AliTRDptrgLUT*>(this->fLUTArray.At(0));
      Int_t* initData = new Int_t[256]; // 2^8
      for (Int_t i = 0; i < 256; i++ ) {
        initData[i] = i;
      }
      lut->InitTable(8, 8, initData, kTRUE); // make copy of initData
      lut = dynamic_cast<AliTRDptrgLUT*>(this->fLUTArray.At(1));
      for (Int_t i = 255; i >= 0; i--) {
        initData[255 - i] = i;  // inverse ramp
      }
      lut->InitTable(8, 8, initData, kTRUE);
    }
    else {
      // initialize threshold
      this->fThreshold = new UInt_t[12];
      for (Int_t i = 0; i < 12; i++) {
        this->fThreshold[i] = 10; 
      }
      
      // initialize LUTsoutputWidth=<value optimized out>
      AliTRDptrgLUT* lut = new AliTRDptrgLUT();
      this->fLUTArray.AddLast(lut);
      lut = new AliTRDptrgLUT(); // this->fRunLoader
      this->fLUTArray.AddLast(lut);
      // the following lines are only needed for test reasons
      lut = dynamic_cast<AliTRDptrgLUT*>(this->fLUTArray.At(0));
      Int_t* initData = new Int_t[4096]; // 2^12
      for (Int_t i = 0; i < 4096; i++ ) {
        initData[i] = i;
      }
      lut->InitTable(12, 12, initData, kTRUE); // make a copy of the table
      lut = dynamic_cast<AliTRDptrgLUT*>(this->fLUTArray.At(1));
      for (Int_t i = 4095; i >= 0; i--) {
        initData[4096 - i] = i;  // inverse ramp
      }
      lut->InitTable(12, 12, initData, kTRUE); // make a copy of the table
      delete[] initData;    
    }
    return false;
  }
  else {
    // load parameters from object
    if (this->fType == kVZERO) {		
      // threshold
      this->fThreshold = 
        this->fParam->GetFEBV0Thresholds(this->fPosition, (this->fID - 1));

      // look up tables
      // 1
      AliTRDptrgLUT* LUT = new AliTRDptrgLUT();
      LUT->InitTable(8, 8, this->fParam->GetFEBV0LUT(this->fPosition, 
                                                     (this->fID - 1),
                                                     0), kFALSE);
      // do not make a copy of the table due to performance reasons
      this->fLUTArray.AddLast(LUT);
 
      // 2
      LUT = new AliTRDptrgLUT(); 
      LUT->InitTable(8, 8, this->fParam->GetFEBV0LUT(this->fPosition, 
                                                     (this->fID - 1),
                                                     1), kFALSE);
      // do not make a copy of the table due to performance reasons
      this->fLUTArray.AddLast(LUT);
    }
    else {		
      // threshold
      this->fThreshold =
        this->fParam->GetFEBT0Thresholds(this->fPosition);

      // look up tables
      // 1
      AliTRDptrgLUT* LUT = new AliTRDptrgLUT();
      LUT->InitTable(12, 12, fParam->GetFEBT0LUT(this->fPosition, 0), kFALSE); 
      // do not make a copy of the table due to performance reasosn
      this->fLUTArray.AddLast(LUT);

      // 2
      LUT = new AliTRDptrgLUT(); 
      LUT->InitTable(12, 12, fParam->GetFEBT0LUT(this->fPosition, 1), kFALSE); 
      // do not make a copy of the table due to performance reasosn      
      this->fLUTArray.AddLast(LUT);
    }
    return true;
  }
  
  return false;
}

//______________________________________________________________________________
Int_t* AliTRDptrgFEB::Simulate()
{
  // simulates the FEB behavior and returns a 2 bit ouput 
  // (least significant bits)
  
  Int_t *result = new Int_t;
  (*result) = -1; 
  if (this->fOperatingMode == kDigits) {
    Int_t inputVector = this->LoadDigits();
    delete result; // delete error return value

    // perform look up
    Int_t nLUTs = this->fLUTArray.GetEntriesFast();  // get LUT count
    result = new Int_t[nLUTs + 1]; // generate new return array
    result[0] = nLUTs; // storage array length in the first array value
    for (Int_t iLUT = 0; iLUT < nLUTs; iLUT++) { 
      // process the return value for each LUT and store the result in the array
      AliDebug(4, Form("FEB: (pos=%d,id=%d,lut=%d,vector=0x%x)", 
                       this->fPosition, this->fID, iLUT, inputVector));

      result[iLUT + 1] = 
       dynamic_cast<AliTRDptrgLUT*>(this->fLUTArray[iLUT])->LookUp(inputVector);
      AliDebug(4, Form("FEB result[%d] = 0x%x",(iLUT + 1),result[iLUT + 1])); 
    }
  }
  else if (this->fOperatingMode == kHits) {
    return result;
  }
  return result;
}
