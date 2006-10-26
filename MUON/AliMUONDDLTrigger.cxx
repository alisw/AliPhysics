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
 
#include <TClonesArray.h>

#include "AliMUONDDLTrigger.h"
#include "AliMUONLocalStruct.h"
#include "AliMUONRegHeader.h"
#include "AliMUONDarcHeader.h"

///
/// \class AliMUONDDLTrigger
///
/// One DDL for trigger part of the MUON spectrometer.
/// 
/// \author C. Finck
///

/// \cond CLASSIMP
ClassImp(AliMUONDDLTrigger)
/// \endcond

//___________________________________________
AliMUONDDLTrigger::AliMUONDDLTrigger()
  :  TObject(),
     fDarcHeader(new AliMUONDarcHeader())
{
 
  ///ctor
 
}

//___________________________________________
AliMUONDDLTrigger::~AliMUONDDLTrigger()
{
  ///dtor
  delete fDarcHeader;

}

//___________________________________________
void AliMUONDDLTrigger::AddRegHeader(const AliMUONRegHeader& regHeader)
{
  /// adding regional header informations 
  /// in the given TClonesArray
  ///
  TClonesArray* regHeaderArray = (TClonesArray*)fDarcHeader->GetRegHeaderArray();

  TClonesArray &eventArray = *regHeaderArray;
  new(eventArray[eventArray.GetEntriesFast()]) AliMUONRegHeader(regHeader);
}

//___________________________________________
void AliMUONDDLTrigger::AddLocStruct(const AliMUONLocalStruct& loc, Int_t iReg)
{
  /// adding local card informations 
  /// for a given regional structure
  /// using TClonesArrays

  AliMUONRegHeader* regHeader = (AliMUONRegHeader*)fDarcHeader->GetRegHeaderEntry(iReg);

  TClonesArray* localArray = (TClonesArray*)regHeader->GetLocalArray();

  TClonesArray &locArray = *localArray;
  new(locArray[locArray.GetEntriesFast()]) AliMUONLocalStruct(loc);
}
