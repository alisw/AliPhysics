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

/*
 *  Class for reading the Trigger Data Stream from Raw.
 *   please read documentation for AliPHOSTriggerRawReader::ReadFromStream .
 */


#include "AliPHOSTriggerRawReader.h"

#include "AliCaloRawStreamV3.h"
#include "AliPHOSTRURawReader.h"

ClassImp(AliPHOSTriggerRawReader)


//________________________________________________________________
AliPHOSTriggerRawReader::AliPHOSTriggerRawReader()
  : TObject(),
    fTRUs()
{
  // default constructor.

  // Initialise fTRUs
  for(Int_t mod=0; mod<kNMods; mod++)
    for(Int_t row=0; row<kNTRURows; ++row)
      for(Int_t branch=0; branch<kNBranches; ++branch)
	fTRUs[mod][row][branch] = 0;
}

    
//________________________________________________________________
AliPHOSTriggerRawReader::~AliPHOSTriggerRawReader()
{
  // destructor
  
  for(Int_t mod = 0; mod < kNMods; ++mod)
    for(Int_t row = 0; row < kNTRURows; ++row)
      for(Int_t branch = 0; branch < kNBranches; branch++)
	delete fTRUs[mod][row][branch];
}
 

//________________________________________________________________
AliPHOSTRURawReader* AliPHOSTriggerRawReader::GetTRU(Int_t mod, Int_t truRow, Int_t branch)
{
  // Get TRU Raw Reader.

  if (mod<0 || mod>=kNMods) return 0x0;
  if (truRow<0 || truRow>=kNTRURows) return 0x0;
  if (branch<0 || branch>=kNBranches) return 0x0;

  if( ! fTRUs[mod][truRow][branch] )
    fTRUs[mod][truRow][branch] = new AliPHOSTRURawReader();
  return fTRUs[mod][truRow][branch];
}
 

//________________________________________________________________
void AliPHOSTriggerRawReader::ReadFromStream(AliCaloRawStreamV3* rawStream)
{
  // Give a AliCaloRawStreamV3* to an instance of this class.
  // It will read from the stream. The stream is passed to 'AliPHOSTRURawReader's
  // which are accesible through 'AliPHOSTriggerRawReader::GetTRU'.
  // note that @param rawStream will not be left in the same state in terms of
  // bunch position, i.e. rawStream->NextBunch() will be called.
  //
  // It is up to the user to check that
  // the is at a channel which is tru data, i.e.:
  //   while (rawStream->NextDDL()) {
  //     while (rawStream->NextChannel()) {
  //       if ( rawStream->IsTRUData() ) 
  // 	  triggerReader->ReadFromStream(rawStream);
  //       else
  // 	  // do something else
  //     }
  //   } 
  // . Other uses will result in undefined behaviour!

  while (rawStream->NextBunch()) {
    Int_t module = rawStream->GetModule();
    Int_t rcuRow = rawStream->GetRow();
    Int_t branch = 1 - rawStream->GetBranch(); // !!! Found this to be necessary, -Henrik Qvigstad <henrik.qvigstad@cern.ch>
    
    AliPHOSTRURawReader* reader = GetTRU(module, rcuRow, branch);
    if (reader) reader->ReadFromStream(rawStream);
  } // end while 
}

//________________________________________________________________
void AliPHOSTriggerRawReader::Reset()
{
  // Reset

  for(Int_t mod = 0; mod < kNMods; ++mod)
    for(Int_t truRow = 0; truRow < kNTRURows; ++truRow)
      for(Int_t branch = 0; branch < kNBranches; branch++)
	if( fTRUs[mod][truRow][branch] )
	  fTRUs[mod][truRow][branch]->Reset();
}
