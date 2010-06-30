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

//-----------------------------------------------------------------------------
/// \class AliMUONRawStreamTrigger
/// This class provides access to MUON digits in raw data.
///
/// It loops over all MUON digits in the raw data given by the AliRawReader.
/// The Next method goes to the next local response. If there are no local response left
/// it returns kFALSE.
/// It can loop also over DDL and store the decoded rawdata in TClonesArrays
/// in payload class.
/// 
/// Version implement for Trigger
/// \author Christian Finck
//-----------------------------------------------------------------------------

#include <TArrayS.h>

#include "AliMUONRawStreamTrigger.h"
#include "AliMUONDarcHeader.h"
#include "AliMUONRegHeader.h"
#include "AliMUONLocalStruct.h"
#include "AliMUONDDLTrigger.h"
#include "AliMUONLogger.h"

#include "AliRawReader.h"
#include "AliRawDataHeader.h"
#include "AliDAQ.h"
#include "AliLog.h"

#include <cassert>

/// \cond CLASSIMP
ClassImp(AliMUONRawStreamTrigger)
/// \endcond

const Int_t AliMUONRawStreamTrigger::fgkMaxDDL = 2;

//___________________________________________
AliMUONRawStreamTrigger::AliMUONRawStreamTrigger()
:   AliMUONVRawStreamTrigger(),
    fPayload(new AliMUONPayloadTrigger()),
    fCurrentDDL(0x0),
    fCurrentDDLIndex(fgkMaxDDL),
    fCurrentDarcHeader(0x0),
    fCurrentRegHeader(0x0),
    fCurrentRegHeaderIndex(0),
    fCurrentLocalStruct(0x0),
    fCurrentLocalStructIndex(0),
    fLocalStructRead(kFALSE),
    fDDL(0),
    fNextDDL(kFALSE)
{
  ///
  /// create an object to read MUON raw digits
  /// Default ctor for monitoring purposes
  ///


}

//_________________________________________________________________
AliMUONRawStreamTrigger::AliMUONRawStreamTrigger(AliRawReader* rawReader)
  : AliMUONVRawStreamTrigger(rawReader),
    fPayload(new AliMUONPayloadTrigger()),
    fCurrentDDL(0x0),
    fCurrentDDLIndex(fgkMaxDDL),
    fCurrentDarcHeader(0x0),
    fCurrentRegHeader(0x0),
    fCurrentRegHeaderIndex(0),
    fCurrentLocalStruct(0x0),
    fCurrentLocalStructIndex(0),
    fLocalStructRead(kFALSE),
    fDDL(0),
    fNextDDL(kFALSE)
{
  ///
  /// ctor with AliRawReader as argument
  /// for reconstruction purpose
  ///

}

//___________________________________
AliMUONRawStreamTrigger::~AliMUONRawStreamTrigger()
{
  ///
  /// clean up
  ///
  delete fPayload;
}

//_____________________________________________________________
Bool_t AliMUONRawStreamTrigger::Next(UChar_t& id,   UChar_t& dec,      Bool_t& trigY, 
				     UChar_t& yPos, UChar_t& sXDev,    UChar_t& xDev,
				     UChar_t& xPos, Bool_t& triggerY,  Bool_t& triggerX,
				     TArrayS& xPattern, TArrayS& yPattern)
{
  ///
  /// read the next raw digit (local structure)
  /// returns kFALSE if there is no digit left
  /// Should call First() before this method to start the iteration.
  ///
  
  if ( IsDone() ) return kFALSE;
  
  if ( fLocalStructRead ) {

    Bool_t ok = GetNextLocalStruct();
    if (!ok)
    {
      // this is the end
      return kFALSE;
    } 
  }

  fLocalStructRead = kTRUE;

  id    = fCurrentLocalStruct->GetId();
  dec   = fCurrentLocalStruct->GetDec(); 
  trigY = fCurrentLocalStruct->GetTrigY();
  yPos  = fCurrentLocalStruct->GetYPos();
  sXDev = fCurrentLocalStruct->GetSXDev();
  xDev  = fCurrentLocalStruct->GetXDev();
  xPos  = fCurrentLocalStruct->GetXPos();

  triggerX = fCurrentLocalStruct->GetTriggerX();
  triggerY = fCurrentLocalStruct->GetTriggerY();

  fCurrentLocalStruct->GetXPattern(xPattern);
  fCurrentLocalStruct->GetYPattern(yPattern);

  return kTRUE;
}

//______________________________________________________
Bool_t AliMUONRawStreamTrigger::IsDone() const
{
  /// Whether the iteration is finished or not
  return (fCurrentLocalStruct==0);
}

//______________________________________________________
void AliMUONRawStreamTrigger::First()
{
  /// Initialize the iteration process.
  
  fCurrentDDLIndex = -1;
  // Must reset all the pointers because if we return before calling
  // GetNextLocalStruct() the user might call CurrentDDL(), CurrentBlockHeader(),
  // CurrentRegHeader() or CurrentLocalStruct() which should return reasonable
  // results in that case.
  fCurrentDDL         = 0;
  fCurrentDarcHeader  = 0;
  fCurrentRegHeader   = 0;
  fCurrentLocalStruct = 0;
  
  // Find the first non-empty structure
  if (not GetNextDDL()) return;
  if (not GetNextRegHeader()) return;
  GetNextLocalStruct();
}

//______________________________________________________
Bool_t AliMUONRawStreamTrigger::GetNextDDL()
{
  /// Returns the next DDL present
  
  assert( GetReader() != 0 );

  
  Bool_t kFound(kFALSE);
  
  while ( fCurrentDDLIndex < fgkMaxDDL-1 && !kFound ) 
  {
    ++fCurrentDDLIndex;
    GetReader()->Reset();
    GetReader()->Select("MUONTRG",fCurrentDDLIndex,fCurrentDDLIndex);
    if ( GetReader()->ReadHeader() ) 
    {
      kFound = kTRUE;
    }
  }
  
  if ( !kFound ) 
  {
    // fCurrentDDLIndex is set to fgkMaxDDL so that we exit the above loop immediately
    // for a subsequent call to this method, unless NextEvent is called in between.
    fCurrentDDLIndex = fgkMaxDDL;
    // We have not actually been able to complete the loading of the new DDL so
    // we are still on the old one. In this case we do not need to reset fCurrentDDL.
    //fCurrentDDL = 0;
    if (IsErrorLogger()) AddErrorMessage();
    return kFALSE;
  }
  
  Int_t totalDataWord  = GetReader()->GetDataSize(); // in bytes
  
  Bool_t scalerEvent =  GetReader()->GetDataHeader()->GetL1TriggerMessage() & 0x1;

  AliDebug(3, Form("DDL Number %d totalDataWord %d\n", fCurrentDDLIndex,
                   totalDataWord));

  UInt_t *buffer = new UInt_t[totalDataWord/4];

  if ( !GetReader()->ReadNext((UChar_t*)buffer, totalDataWord) )
  {
    // We have not actually been able to complete the loading of the new DDL so
    // we are still on the old one. In this case we do not need to reset fCurrentDDL.
    //fCurrentDDL = 0;
    delete [] buffer;
    return kFALSE;
  }

#ifndef R__BYTESWAP  
  Swap(buffer, totalDataWord / sizeof(UInt_t)); // swap needed for mac power pc
#endif

  fPayload->ResetDDL();
  


  Bool_t ok = fPayload->Decode(buffer, scalerEvent);

  delete[] buffer;
  
  fCurrentDDL = fPayload->GetDDLTrigger();
  
  fCurrentDarcHeader = fCurrentDDL->GetDarcHeader();
  
  fCurrentRegHeaderIndex = -1;


  return ok;
}


//______________________________________________________
Bool_t AliMUONRawStreamTrigger::GetNextRegHeader()
{
  /// Returns the next Reg Header present

  assert( fCurrentDarcHeader != 0 );
  assert( fCurrentDDL != 0 );

  fCurrentRegHeader = 0;
  
  Int_t i = fCurrentRegHeaderIndex;
  
  while ( fCurrentRegHeader == 0 && i < fCurrentDarcHeader->GetRegHeaderEntries()-1 )
  {
    ++i;
    fCurrentRegHeader = fCurrentDarcHeader->GetRegHeaderEntry(i);
  }
     
  if ( !fCurrentRegHeader ) 
  {
    Bool_t ok = GetNextDDL();
    if (!ok) 
    {
      return kFALSE;
    }
    else
    {
      return GetNextRegHeader();
    }
  }
  
  fCurrentRegHeaderIndex = i;
  
  fCurrentLocalStructIndex = -1;
  
  return kTRUE;
}

//______________________________________________________
Bool_t AliMUONRawStreamTrigger::GetNextLocalStruct()
{
  /// Find the next non-empty local structure
  
  assert( fCurrentRegHeader != 0 );
  
  fCurrentLocalStruct = 0;

  Int_t i = fCurrentLocalStructIndex;
  
  while ( fCurrentLocalStruct == 0 && i < fCurrentRegHeader->GetLocalEntries()-1 ) 
  {
    ++i;
    fCurrentLocalStruct = fCurrentRegHeader->GetLocalEntry(i);
  }
    
  if ( !fCurrentLocalStruct ) 
  {
    Bool_t ok = GetNextRegHeader();
    if (!ok)
    {
      return kFALSE;
    }
    else
    {
      return GetNextLocalStruct();
    }
  }
  
  fCurrentLocalStructIndex = i;
  
  fLocalStructRead = kFALSE;
  
  return kTRUE;
}

//______________________________________________________
Bool_t AliMUONRawStreamTrigger::NextDDL()
{
  /// reading tracker DDL
  /// store local info into Array
  /// store only non-empty structures

  // reset TClones
  fPayload->ResetDDL();


  // loop over the two ddl's

  while ( fDDL < fgkMaxDDL ) {
    GetReader()->Reset();
    GetReader()->Select("MUONTRG", fDDL, fDDL);  //Select the DDL file to be read  
    if (GetReader()->ReadHeader()) break;
    AliDebug(3,Form("Skipping DDL %d which does not seem to be there",fDDL));
    ++fDDL;
  }

  if (fDDL >= fgkMaxDDL) {
    fDDL = 0;
    if (IsErrorLogger()) AddErrorMessage();
    return kFALSE;
  }

  AliDebug(3, Form("DDL Number %d\n", fDDL ));

  Int_t totalDataWord = GetReader()->GetDataSize(); // in bytes

  Bool_t scalerEvent =  GetReader()->GetDataHeader() && GetReader()->GetDataHeader()->GetL1TriggerMessage() & 0x1;


  UInt_t *buffer = new UInt_t[totalDataWord/4];

  // check not necessary yet, but for future developments
  if (!GetReader()->ReadNext((UChar_t*)buffer, totalDataWord)) return kFALSE; 
  
#ifndef R__BYTESWAP
  Swap(buffer, totalDataWord / sizeof(UInt_t)); // swap needed for mac power pc
#endif

  fPayload->Decode(buffer, scalerEvent);


  fDDL++;

  delete [] buffer;


  return kTRUE;
}

// //______________________________________________________
// void AliMUONRawStreamTrigger::SetMaxReg(Int_t reg) 
// {
//   /// set regional card number
//   fPayload->SetMaxReg(reg);
// }

//______________________________________________________
void AliMUONRawStreamTrigger::SetMaxLoc(Int_t loc) 
{
  /// set local card number
  fPayload->SetMaxLoc(loc);
}

//______________________________________________________
void AliMUONRawStreamTrigger::AddErrorMessage()
{
/// add message into logger of AliRawReader per event

  TString msg;
  Int_t occurance = 0;
  AliMUONLogger* log = fPayload->GetErrorLogger();
  
  log->ResetItr();
  while(log->Next(msg, occurance))
  { 
    if (msg.Contains("Darc"))
      GetReader()->AddMajorErrorLog(kDarcEoWErr, msg.Data());

    if (msg.Contains("Global"))
      GetReader()->AddMajorErrorLog(kGlobalEoWErr, msg.Data());

    if (msg.Contains("Regional"))
      GetReader()->AddMajorErrorLog(kRegEoWErr, msg.Data());

    if (msg.Contains("Local"))
      GetReader()->AddMajorErrorLog(kLocalEoWErr, msg.Data());
  }
  
  log->Clear(); // clear after each event
}
