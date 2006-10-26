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

#include "AliMUONDDLTracker.h"
#include "AliMUONBusStruct.h"
#include "AliMUONDspHeader.h"
#include "AliMUONBlockHeader.h"

///
/// \class AliMUONDDLTracker
///
/// A wrapper object for 1 DDL of the MUON tracking chambers.
///
/// \author C. Finck

/// \cond CLASSIMP
ClassImp(AliMUONDDLTracker)
/// \endcond

//___________________________________________
AliMUONDDLTracker::AliMUONDDLTracker()
  :  TObject(),
     fBlkHeaderArray(new TClonesArray("AliMUONBlockHeader", 2))
{
  ///
  ///ctor
  ///

}

//___________________________________________
AliMUONDDLTracker::~AliMUONDDLTracker()
{
  ///
  ///dtor
  ///
  fBlkHeaderArray->Delete();
  delete fBlkHeaderArray;

}

//___________________________________________
void AliMUONDDLTracker::AddBusPatch(const AliMUONBusStruct& busPatch, Int_t iBlock, Int_t iDsp )
{
  /// adding bus patch informations
  /// for a given block & Dsp structure
  /// using TClonesArrays

  AliMUONBlockHeader* blockHeader = (AliMUONBlockHeader*)fBlkHeaderArray->At(iBlock);
  AliMUONDspHeader* dspHeader     = (AliMUONDspHeader*)blockHeader->GetDspHeaderEntry(iDsp);

  TClonesArray* busPatchArray = (TClonesArray*)dspHeader->GetBusPatchArray();

  TClonesArray &eventArray = *busPatchArray;
  new(eventArray[eventArray.GetEntriesFast()]) AliMUONBusStruct(busPatch);
}

//___________________________________________
void AliMUONDDLTracker::AddDspHeader(const AliMUONDspHeader& dspHeader, Int_t iBlock)
{
  /// adding DspHeader informations
  /// for a given block structure
  /// using TClonesArrays

  AliMUONBlockHeader* blockHeader = (AliMUONBlockHeader*)fBlkHeaderArray->At(iBlock);

  TClonesArray* dspHeaderArray = (TClonesArray*)blockHeader->GetDspHeaderArray();

  TClonesArray &dspArray = *dspHeaderArray;
  new(dspArray[dspArray.GetEntriesFast()]) AliMUONDspHeader(dspHeader);
}

//___________________________________________
void AliMUONDDLTracker::AddBlkHeader(const AliMUONBlockHeader& blkHeader)
{
  /// adding Block header informations
  /// for a given block structure
  /// using TClonesArrays

  TClonesArray &blkArray = *fBlkHeaderArray;
  new(blkArray[blkArray.GetEntriesFast()]) AliMUONBlockHeader(blkHeader);
}

//___________________________________________
void AliMUONDDLTracker::Clear(Option_t* )
{
  /// Clear TClones arrays
  /// instead of deleting
  ///
  fBlkHeaderArray->Clear("C");

}
