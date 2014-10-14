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

////////////////////////////////////////////////////////////////////////////////
//                                                                        
//  Pre-Trigger simulation                                                
//                                                                        
//  Authors: F. Reidt (Felix.Reidt@cern.ch)                               
//                                                                        	
//                                                                        
//  Limitations: input/output width: 32 bits (UInt_t)                      
//                                                                        
//  Annotation: That LUT is usually used to provide a single output bit   
//              In that case every output value bigger 0 means true       
//
////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <fstream>
#include <string>
#include <math.h>

#include "TFile.h"
#include "TROOT.h"

#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliLog.h"

#include "AliTRDptrgLUT.h"

ClassImp(AliTRDptrgLUT)
//_____________________________________________________________________________
AliTRDptrgLUT::AliTRDptrgLUT() 
  : TObject(),
  fLUTData(0),
  fInputWidth(0),
  fOutputWidth(0),
  fTableEntryCount(0),
  fCopiedTable(kFALSE)
{
  // ctor
}

//_____________________________________________________________________________
AliTRDptrgLUT::~AliTRDptrgLUT() 
{
  // destructor
  if (this->fCopiedTable) {
    AliDebug(5, "Deleted LUT data");
    if (this->fLUTData != 0x0) {
      delete[] this->fLUTData;
    }
    this->fLUTData = 0x0;
  }
}

//_____________________________________________________________________________
Int_t AliTRDptrgLUT::LookUp(UInt_t input)
{
  // perform a look up 
  
  if (input > (UInt_t)this->fTableEntryCount) {
    // check whether the input value is out of bounds
    AliWarning("incorrect LUT input value");
    return -1;
  }
  return this->fLUTData[input]; // do look up and output 
}

//______________________________________________________________________________
Int_t AliTRDptrgLUT::InitTable(Int_t inputWidth, Int_t outputWidth, 
                               Int_t *tableData, Bool_t copy)
{
  // load and initialize the look up table
  
  // assign width
  this->fInputWidth = inputWidth;
  this->fOutputWidth = outputWidth;
  
  // calculated table entry count 
  this->fTableEntryCount = 0x1;
  this->fTableEntryCount <<= inputWidth; 
  AliDebug(5,Form("fTableEntryCount=%d", this->fTableEntryCount));
 
  this->fCopiedTable = copy; 
 
  if (copy) {
    this->fLUTData = new Int_t[this->fTableEntryCount]; // allocate data table
    for (Int_t i=0; i < this->fTableEntryCount; i++) {
      this->fLUTData[i] = tableData[i];
    }
  }
  else { // do not copy (due to performace reasons)
    this->fLUTData = tableData;
  }
  return 0;
}

