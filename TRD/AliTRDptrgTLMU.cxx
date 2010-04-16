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
//  Pre-Trigger simulation                                                //
//                                                                        //
//  Authors: F. Reidt (Felix.Reidt@cern.ch)                               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TROOT.h"

#include "TClonesArray.h"

#include "AliRun.h"
#include "AliRunLoader.h"

#include "../TOF/AliTOFdigit.h"
#include "../RAW/AliRawReader.h" // needed for AliTOFTrigger's raw digit support
#include "../TOF/AliTOFTrigger.h"

#include "AliTRDptrgParam.h"

#include "AliTRDptrgTLMU.h"

ClassImp(AliTRDptrgTLMU)

//______________________________________________________________________________
AliTRDptrgTLMU::AliTRDptrgTLMU(AliRunLoader *rl) 
  : TObject(),
  fRunLoader(rl),
  fParam(0x0),
  fOperatingMode(kDigits),
  fInputMask(0x0),
  fCMatrices(0x0),
  fMultiplicity(0x0),
  fOutput(0x0) 
{
  // default ctor
 
  for (Int_t i = 0; i < 18; i++) {
    this->fTOFinputBits[i] = 0;
  }
  
  this->LoadParams();
}

//______________________________________________________________________________
AliTRDptrgTLMU::AliTRDptrgTLMU(AliRunLoader *rl,  AliTRDptrgParam *param,
                               AliTRDptrgOperatingMode_t operatingMode)
  : TObject(),
  fRunLoader(rl),
  fParam(param),
  fOperatingMode(operatingMode),
  fInputMask(0x0),
  fCMatrices(0x0),
  fMultiplicity(0x0),
  fOutput(0x0) 
{
  // recommended ctor
  
  for (Int_t i = 0; i < 18; i++) {
    this->fTOFinputBits[i] = 0;
  }
 
  this->LoadParams();
}

//______________________________________________________________________________
AliTRDptrgTLMU::~AliTRDptrgTLMU() 
{
  // destructor
  
  this->fCMatrices = 0x0;
  this->fMultiplicity = 0x0;
  this->fOutput = 0x0;
}

//______________________________________________________________________________
Int_t* AliTRDptrgTLMU::Simulate() 
{
  // starts a simulation
	
  if (this->fOperatingMode == kDigits) {
    this->LoadDigits();
  }	
  else if (this->fOperatingMode == kHits) {
    return 0x0; // TODO
  }

  // which super modules are active? - needed for fast coincidence processing:  
  // M(sm | coincidence condition)>=2
  UInt_t temp = 0x1;
  UInt_t sm = 0x0; 
  for (Int_t iSM = 0; iSM < 18; iSM++) {
    if (this->Or(iSM)) {
      sm |= temp;
    }
    temp <<= 1;
  }
  AliDebug(4, Form("Active supermodules: 0x%x", sm));

  // get multiplicity
  UInt_t multiplicity = this->GetMultiplicitySum();
  AliDebug(4, Form("TOF/TLMU multiplicity: %d", multiplicity)); 	

  Int_t* result = new Int_t[9]; 
  result[0] = 8;
  for (Int_t iResult = 0; iResult < 8; iResult++) {
    result[iResult + 1] = 0;

    // coincidence matrix
    if (this->fOutput[iResult][0] != -1) {
      for (Int_t iLine = 0; iLine < 18; iLine++) {
        AliDebug(5, Form("Entry: %d, matrix: %d, line: %d, output: 0x%x", 
                         iResult, this->fOutput[iResult][0], iLine, sm & 
                         this->fCMatrices[this->fOutput[iResult][0]][iLine]));
        if (this->GetBitVectorMultiplicity(
	     sm & this->fCMatrices[this->fOutput[iResult][0]][iLine]) > 1) {
	  result[iResult + 1] = 1;
          break;
	}        
      }
    }
    
    // multiplicity conditions
    if (this->fOutput[iResult][1] != -1) {
      AliDebug(5, Form("Entry: %d, slice: %d", iResult, 
                       this->fOutput[iResult][1]));
      if ((this->fMultiplicity[this->fOutput[iResult][1]][0] < multiplicity) &&
          (this->fMultiplicity[this->fOutput[iResult][1]][1] >= multiplicity)) {
        result[iResult + 1] = 1;
      } 
    }
    AliDebug(4, Form("TLMU result[%d] = %d", iResult, result[iResult + 1]));
  }

  return result;
}

//______________________________________________________________________________
Int_t AliTRDptrgTLMU::LoadDigits()
{
  // loads Digits (for usage with aquired data)
  this->GetInputBits(); // get bits from AliTOFTrigger
  return 0;
}


//______________________________________________________________________________
Bool_t AliTRDptrgTLMU::LoadParams()
{
  // load AliTRDprtgParam content
	
  if (this->fParam == 0x0) {
    // no parameter object assigned
    AliWarning("no parameter object assigned - using default settings!");

    UInt_t* imask = 0x0;
    imask = new UInt_t[18];
    for (Int_t i = 0; i < 18; i++) {
      imask[i] = 0xFFFFFFFF;
    }

    this->fInputMask = imask;
    
    // TLMU Coincidence Matrices
    this->fCMatrices = new UInt_t*[3];
    this->fCMatrices[0] = new UInt_t[18];
    this->fCMatrices[1] = new UInt_t[18];
    this->fCMatrices[2] = new UInt_t[18];    

    // Matrix 0: Back-To-Back
    // Matrix 1: Back-To-Back +/-1
    // Matrix 2: Back-To-Back +/-2
    for (UInt_t iMatrix = 0; iMatrix < 3; iMatrix++) {
      for (UInt_t iSlice = 0; iSlice < 18; iSlice++) {
        if (iMatrix == 0) {
          if (iSlice < 9) {
            this->fCMatrices[iMatrix][iSlice] = 0x201 << iSlice; 
            // Back-To-Back 
            AliDebug(5, Form("fCMatrices[%d][%d]=0x%x",iMatrix,iSlice,
                             this->fCMatrices[iMatrix][iSlice]));
	  }
          // because of symmetrie the other slices are not necessary
        } 
        else if (iMatrix == 1)  {
          // Back-To-Back +/- 1
          if (iSlice < 8) {
            this->fCMatrices[iMatrix][iSlice] = 0x381 << iSlice;
          }
          else if (iSlice == 8) {
            this->fCMatrices[iMatrix][iSlice] = 0x30101;
          }
          else if (iSlice == 9) {
            this->fCMatrices[iMatrix][iSlice] = 0x20203;
          }
          else {
            this->fCMatrices[iMatrix][iSlice] = 0x407 << (iSlice - 10);
          } 
          AliDebug(5, Form("fCMatrices[%d][%d]=0x%x",iMatrix,iSlice,
                           this->fCMatrices[iMatrix][iSlice])); 
        }
        else if (iMatrix == 2) {
          // Back-To-Back +/-2
          if (iSlice < 7 ) {
            this->fCMatrices[iMatrix][iSlice] = 0xF81 << iSlice;
          }
          else if (iSlice == 7) {
            this->fCMatrices[iMatrix][iSlice] = 0x3C081;
          }
          else if (iSlice == 8) {
            this->fCMatrices[iMatrix][iSlice] = 0x38103;
          }
          else if (iSlice == 9) {
            this->fCMatrices[iMatrix][iSlice] = 0x30207;
          }
          else if (iSlice == 10) {
            this->fCMatrices[iMatrix][iSlice] = 0x2040F;
          }
          else {
            this->fCMatrices[iMatrix][iSlice] = 0x81F << (iSlice - 11);
          } 
          AliDebug(5, Form("fCMatrices[%d][%d]=0x%x",iMatrix,iSlice,
                           this->fCMatrices[iMatrix][iSlice]));     
        }
      } 
    }
 
    // Mulitplicity
    this->fMultiplicity = new UInt_t*[9];
    for (Int_t i = 0; i < 9; i++) {
      this->fMultiplicity[i] = new UInt_t[2];
    }
    this->fMultiplicity[0][0] = 0;
    this->fMultiplicity[0][1] = 10;
    this->fMultiplicity[1][0] = 10;
    this->fMultiplicity[1][1] = 25;
    this->fMultiplicity[2][0] = 25;
    this->fMultiplicity[2][1] = 50;
    this->fMultiplicity[3][0] = 50;
    this->fMultiplicity[3][1] = 100;
    this->fMultiplicity[4][0] = 100;
    this->fMultiplicity[4][1] = 200;
    this->fMultiplicity[5][0] = 200;
    this->fMultiplicity[5][1] = 350;
    this->fMultiplicity[6][0] = 350;
    this->fMultiplicity[6][1] = 400;
    this->fMultiplicity[7][0] = 400;
    this->fMultiplicity[7][1] = 576;
    this->fMultiplicity[8][0] = 100;
    this->fMultiplicity[8][1] = 576;
 
    // TLMU output
    this->fOutput = new Int_t*[8];
    for (Int_t i = 0; i < 9; i++) {
      this->fOutput[i] = new Int_t[2];
      this->fOutput[i][0] = -1;
      this->fOutput[i][1] = -1;
    }
    this->fOutput[0][0] = 0;
    this->fOutput[1][0] = 1;
    this->fOutput[2][0] = 2;
    this->fOutput[3][1] = 0;
    this->fOutput[4][1] = 1;
    this->fOutput[5][1] = 2;
    this->fOutput[6][1] = 3;
    this->fOutput[7][1] = 8;


  }
  else {
    // parameter object assigned
    
    this->fInputMask = this->fParam->GetTLMUInputMask(); 
    // input mask for TOF-bits (18x32=576)
    
    this->fCMatrices = this->fParam->GetTLMUcmatrices();
    // get coincidence matrices
 
    this->fMultiplicity = this->fParam->GetTLMUmultiplicity();
    // get multiplicity slices
  
    this->fOutput = this->fParam->GetTLMUoutput();
    // get output signal assignment
  }

  return false;
}

//______________________________________________________________________________
void AliTRDptrgTLMU::GetInputBits() {
  // Gets TOF-to-TRD input bits from AliTOFTrigger as Bool_t array

  AliTOFTrigger *toftrig = new AliTOFTrigger(); // create AliTOFTrigger 
  toftrig->CreateLTMMatrixFromDigits(); // Generate LTMMatrix from AliTOFdigits
  
  // prepare map  
  Bool_t** map = 0x0;
  map = new Bool_t*[72];
  for (Int_t i=0; i < 72; i++) 
    map[i] = new Bool_t[8];
  
  // initialise map
  for (Int_t i=0; i < 72; i++)
    for (Int_t j=0; j < 8; j++)
      map[i][j] = kFALSE;

  // get 576 TOF-to-TRD bits
  toftrig->GetTRDmap(map);

  //* DEBUG output 
  // used to determine the correct bit assignment
  AliDebug(5, "AliTOFTrigger->GetTRDmap(map):");
  for (Int_t i=0; i < 72; i++) {
    AliDebug(5, Form("%d %d%d%d%d%d%d%d%d", i, map[i][0], map[i][1], map[i][2],
                      map[i][3], map[i][4], map[i][5], map[i][6], map[i][7]));
  }
  //*/ // end of DEBUG output

  // initialise fTOFinputBits
  for (Int_t i=0; i < 18; i++) {
    this->fTOFinputBits[i] = 0;
  }

 
  // transform Bool_t array to UInt_t bitvectors according to
  // http://www.physi.uni-heidelberg.de/~schicker/cbtof/cbtof_docu.pdf
  // chapter 1.4 and 2.1 to a supermodule based structured 
  // integer (bit) array
  Int_t supermodule = -1;
  UInt_t tempA = 0x00000001;
  UInt_t tempC = 0x00010000;
  for (Int_t iLTM = 0; iLTM < (kNLTM / 2); iLTM++) {
    if (!(iLTM % 2)) { // renew temp vectors, update supermodule
      tempA = 0x00000001;
      tempC = 0x00010000;
      supermodule++;
    }
    // AliDebug(5, Form("(%2d,0x%8x,0x%8x)", iLTM, tempA, tempC));
    for (Int_t iLTMchan = 0; iLTMchan < 8; iLTMchan++) {
      // A-side
      if (map[iLTM][iLTMchan]) {
        this->fTOFinputBits[supermodule] |= tempA;        
      }
      // C-side
      if (map[iLTM + 36][iLTMchan]) {
        this->fTOFinputBits[supermodule] |= tempC;
      }
      // change temp vectors
      tempA <<= 1;
      tempC <<= 1;
    }
  }

  // handle input mask
  for (Int_t iSM = 0; iSM < 18; iSM++) {
    AliDebug(5, Form("fInputTOFinputBits[%d]=0x%x", iSM, 
             this->fTOFinputBits[iSM]));
    this->fTOFinputBits[iSM] &= this->fInputMask[iSM];
  }
}


//______________________________________________________________________________
Int_t AliTRDptrgTLMU::BackToBack(Int_t iSM, Int_t range) {
  // Check whether there is an back-to-back particle trace

  // get the counterpart of supermodule iSM
  Int_t counterPart = -1;
  if (iSM >= 9) { 
    counterPart = iSM - 9;
  }
  else  {
    counterPart = iSM + 9;
  }

  if (this->Or(iSM)) { // is there are active bits in supermodule iSM
    Int_t result = 0;
    for (Int_t i = counterPart - range; i <= counterPart + range; i++) {
      // check whether there are active bits in supermodule i (counterParts)
      if ((i >= 0) && (i < 18)) {
        if (Or(i)) {
          result++;
        }
      }
      else {
	if (i < 0) {
          if (Or(17 - i)) {
            result++;
          }
        }
        if (i > 17) {
          if (Or(i - 18)) {
            result++;
          }
        }
      }
    }
    AliDebug(5, Form("BackToBack of %d and %d+-%d\n: %d", iSM, counterPart, 
                     range, result));
    return result; // return whether there was a possible back-to-back trace
  }    
  else {
    AliDebug(5, Form("BackToBack unsuccessful, not hit in sm%d", iSM));
    return 0; // iSM did not recognize anything
  }
}

//______________________________________________________________________________
Int_t AliTRDptrgTLMU::Coincidence(Int_t iSM1, Int_t iSM2) {
  // checks whether there is an coincidence in iSM1 and iSM2

  if (this->Or(iSM1) && this->Or(iSM2)) {
    return 1;
  }
  else
    return 0;
}

//______________________________________________________________________________
inline Int_t AliTRDptrgTLMU::Or(Int_t  iSM) {
  // returns 1 if one or more bits are active in supermodule iSM
 
  if ((iSM >= 0) && (iSM < 18)) {
    if (this->fTOFinputBits[iSM] > 0)
      return 1;
    else
      return 0; 
  }
  else {
    return -1;
  }
}

//______________________________________________________________________________
Int_t AliTRDptrgTLMU::GetMultiplicity(Int_t iSM) {
  // counts how many bits equal one are in class member fTOFinputBits[iSM]
  // (32bits from TOF to TRD of supermodule iSM)

  UInt_t temp = this->fTOFinputBits[iSM];
  UInt_t mask = 0x01;
  Int_t multiplicity = 0;  
	
  for (int iBit = 0; iBit < 32; iBit++) {
    if ((mask & temp) != 0x0) { // is the bit equal one?
      multiplicity++;
    }
    mask <<= 1; // rotate mask to the left after each iteration
  }
  AliDebug(5, Form("Multiplicity of supermodule %d: %d", iSM, multiplicity));
  return multiplicity;
}

//______________________________________________________________________________
Int_t AliTRDptrgTLMU::GetMultiplicitySum() {
  // returns the multiplicity of the whole detector (all 576bit TOF to TRD bits)
  Int_t sum = 0;
  
  for (Int_t i = 0; i < 18; i++) {
    sum += this->GetMultiplicity(i);
  }
  AliDebug(5, Form("Whole multiplicity: %d", sum));
  return sum;
}

//______________________________________________________________________________
UInt_t AliTRDptrgTLMU::GetBitVectorMultiplicity(UInt_t BitVector) {
  // returns the multiplicity of a given bit vector
  
  UInt_t result = 0;
  UInt_t temp = 0x01;
  for (UInt_t iBit = 0; iBit < 32; iBit++) {
    if (BitVector & temp) {
      result++;
    }
    temp <<= 1;
  }

  return result;
}
