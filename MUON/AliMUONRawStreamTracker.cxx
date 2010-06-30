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
/// \class AliMUONRawStreamTracker
/// This class provides access to MUON digits in raw data.
///
/// It loops over all MUON digits in the raw data given by the AliRawReader.
/// The Next method goes to the next digit. If there are no digits left
/// it returns kFALSE
/// It can loop also over DDL and store the decoded rawdata in TClonesArray
/// in Payload class.
/// 
/// Implement for Tracker
///
/// \author Christian Finck & Laurent Aphecetche
//-----------------------------------------------------------------------------

#include "AliMUONRawStreamTracker.h"

#include "AliMUONLogger.h"
#include "AliRawReader.h"
#include "AliRawDataHeader.h"
#include "AliDAQ.h"
#include "AliLog.h"
#include "AliMUONPayloadTracker.h"
#include "AliMUONBlockHeader.h"
#include "AliMUONDspHeader.h"
#include "AliMUONBusStruct.h"
#include "AliMUONDDLTracker.h"
#include "Riostream.h"
#include <cassert>

/// \cond CLASSIMP
ClassImp(AliMUONRawStreamTracker)
/// \endcond

//___________________________________________
AliMUONRawStreamTracker::AliMUONRawStreamTracker()
 : AliMUONVRawStreamTracker(),
   fPayload(new AliMUONPayloadTracker()),
   fCurrentDDL(0),
   fCurrentDDLIndex(fgkMaxDDL),
   fCurrentBlockHeader(0),
   fCurrentBlockHeaderIndex(0),
   fCurrentDspHeader(0),
   fCurrentDspHeaderIndex(0),
   fCurrentBusStruct(0),
   fCurrentBusStructIndex(0),
   fCurrentDataIndex(0),
   fDDL(0)
{
  ///
  /// create an object to read MUON raw digits
  /// Default ctor for monitoring purposes
  ///
  
  
}

//_________________________________________________________________
AliMUONRawStreamTracker::AliMUONRawStreamTracker(AliRawReader* rawReader)
: AliMUONVRawStreamTracker(rawReader),
  fPayload(new AliMUONPayloadTracker()),
  fCurrentDDL(0),
  fCurrentDDLIndex(fgkMaxDDL),
  fCurrentBlockHeader(0),
  fCurrentBlockHeaderIndex(0),
  fCurrentDspHeader(0),
  fCurrentDspHeaderIndex(0),
  fCurrentBusStruct(0),
  fCurrentBusStructIndex(0),
  fCurrentDataIndex(0),
  fDDL(0)
{
  ///
  /// ctor with AliRawReader as argument
  /// for reconstruction purpose
  ///
  
  
}

//___________________________________
AliMUONRawStreamTracker::~AliMUONRawStreamTracker()
{
  ///
  /// clean up
  ///
  delete fPayload;
}

//_____________________________________________________________
Bool_t 
AliMUONRawStreamTracker::Next(Int_t& busPatchId,
                              UShort_t& manuId, UChar_t& manuChannel,
                              UShort_t& adc)
{
  ///
  /// read the next raw digit (buspatch structure)
  /// returns kFALSE if there is no digit left
  /// Should call First() before this method to start the iteration.
  ///
  
  if ( IsDone() ) return kFALSE;
  
  if ( fCurrentDataIndex >= fCurrentBusStruct->GetLength()-1 )
  {
    Bool_t ok = GetNextBusStruct();
    if (!ok)
    {
      // this is the end
      return kFALSE;
    } 
  }

  ++fCurrentDataIndex;

  busPatchId = fCurrentBusStruct->GetBusPatchId();
  manuId = fCurrentBusStruct->GetManuId(fCurrentDataIndex);
  manuChannel = fCurrentBusStruct->GetChannelId(fCurrentDataIndex);
  adc = fCurrentBusStruct->GetCharge(fCurrentDataIndex);

  return kTRUE;
}

//______________________________________________________
Bool_t
AliMUONRawStreamTracker::IsDone() const
{
  /// Whether the iteration is finished or not
  return (fCurrentBusStruct==0);
}

//______________________________________________________
void
AliMUONRawStreamTracker::First()
{
  /// Initialize the iteration process.
  
  fCurrentDDLIndex = -1;
//  fCurrentDspHeaderIndex = -1; // Not necessary since this gets reset in the GetNextXXX methods.
//  fCurrentBusStructIndex = -1;

  // Must reset all the pointers because if we return before calling
  // GetNextBusStruct() the user might call CurrentDDL(), CurrentBlockHeader(),
  // CurrentDspHeader() or CurrentBusStruct() which should return reasonable
  // results in that case.
  fCurrentDDL = 0;
  fCurrentBlockHeader = 0;
  fCurrentDspHeader = 0;
  fCurrentBusStruct = 0;
  
  // Find the first non-empty structure
  if (not GetNextDDL()) return;
  if (not GetNextBlockHeader()) return;
  if (not GetNextDspHeader()) return;
  GetNextBusStruct();
}

//______________________________________________________
Bool_t
AliMUONRawStreamTracker::GetNextDDL()
{
  /// Returns the next DDL present
  
  assert( GetReader() != 0 );
  
  Bool_t kFound(kFALSE);
  
  while ( fCurrentDDLIndex < fgkMaxDDL-1 && !kFound ) 
  {
    ++fCurrentDDLIndex;
    GetReader()->Reset();
    GetReader()->Select("MUONTRK",fCurrentDDLIndex,fCurrentDDLIndex);
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
  fPayload->ResetDDL();
  
#ifndef R__BYTESWAP  
  Swap(buffer, totalDataWord / sizeof(UInt_t)); // swap needed for mac power pc
#endif

  Bool_t ok = fPayload->Decode(buffer, totalDataWord/4);
  
  delete[] buffer;
  
  fCurrentDDL = fPayload->GetDDLTracker();
  
  fCurrentBlockHeaderIndex = -1;
  
  return ok;
}

//______________________________________________________
Bool_t
AliMUONRawStreamTracker::GetNextBlockHeader()
{
  /// Returns the next block Header present
  
  assert( fCurrentDDL != 0 );

  fCurrentBlockHeader = 0;

  Int_t i(fCurrentBlockHeaderIndex);
  
  while ( fCurrentBlockHeader == 0 && i < fCurrentDDL->GetBlkHeaderEntries()-1 ) 
  {
    ++i;
    fCurrentBlockHeader = fCurrentDDL->GetBlkHeaderEntry(i);
  }
  
  if ( !fCurrentBlockHeader ) 
  {
    Bool_t ok = GetNextDDL();
    if (!ok) 
    {
      return kFALSE;
    }
    else
    {
      return GetNextBlockHeader();
    }
  }
  
  fCurrentBlockHeaderIndex = i;
  
  fCurrentDspHeaderIndex = -1;
  
  return kTRUE;
}

//______________________________________________________
Bool_t
AliMUONRawStreamTracker::GetNextDspHeader()
{
  /// Returns the next Dsp Header present

  assert( fCurrentBlockHeader != 0 );
  
  fCurrentDspHeader = 0;
  
  Int_t i(fCurrentDspHeaderIndex);
  
  while ( fCurrentDspHeader == 0 && i < fCurrentBlockHeader->GetDspHeaderEntries()-1 )
  {
    ++i;
    fCurrentDspHeader = fCurrentBlockHeader->GetDspHeaderEntry(i);
  }
  
  if ( !fCurrentDspHeader ) 
  {
    Bool_t ok = GetNextBlockHeader();
    if (!ok) 
    {
      return kFALSE;
    }
    else
    {
      return GetNextDspHeader();
    }
  }
  
  fCurrentDspHeaderIndex = i;
  
  fCurrentBusStructIndex = -1;
  
  return kTRUE;
}

//______________________________________________________
Bool_t
AliMUONRawStreamTracker::GetNextBusStruct()
{
  /// Find the next non-empty busPatch structure
  
  assert( fCurrentDspHeader != 0 );
  
  fCurrentBusStruct = 0;

  Int_t i(fCurrentBusStructIndex);
  
  while ( fCurrentBusStruct == 0 && i < fCurrentDspHeader->GetBusPatchEntries()-1 ) 
  {
    ++i;
    fCurrentBusStruct = fCurrentDspHeader->GetBusPatchEntry(i);
  }
    
  if ( !fCurrentBusStruct ) 
  {
    Bool_t ok = GetNextDspHeader();
    if (!ok)
    {
      return kFALSE;
    }
    else
    {
      return GetNextBusStruct();
    }
  }
  
  fCurrentBusStructIndex = i;
  
  fCurrentDataIndex = -1;
  
  return kTRUE;
}

//______________________________________________________
Bool_t AliMUONRawStreamTracker::NextDDL()
{
  /// reading tracker DDL
  
  assert( GetReader() != 0 );
  
  fPayload->ResetDDL();
  
  while ( fDDL < fgkMaxDDL ) 
  {
    GetReader()->Reset();
    GetReader()->Select("MUONTRK", fDDL, fDDL);  //Select the DDL file to be read  
    if (GetReader()->ReadHeader()) break;
    AliDebug(3,Form("Skipping DDL %d which does not seem to be there",fDDL));
    ++fDDL;
  }
  
  if ( fDDL == fgkMaxDDL ) 
  {
    fDDL = 0;
    if ( IsErrorLogger()) AddErrorMessage();
    return kFALSE;
  }
  
  AliDebug(3, Form("DDL Number %d\n", fDDL ));
  
  Int_t totalDataWord  = GetReader()->GetDataSize(); // in bytes
  
  UInt_t *buffer = new UInt_t[totalDataWord/4];

  if(!GetReader()->ReadNext((UChar_t*)buffer, totalDataWord))
  {
    delete[] buffer;
    return kFALSE;
  }

#ifndef R__BYTESWAP  
  Swap(buffer, totalDataWord / sizeof(UInt_t)); // swap needed for mac power pc
#endif
  
  Bool_t ok = fPayload->Decode(buffer, totalDataWord/4);

  delete[] buffer;
  
  fDDL++;
  
  return ok;
}

//______________________________________________________
void AliMUONRawStreamTracker::SetMaxBlock(Int_t blk) 
{
  /// set regional card number
  fPayload->SetMaxBlock(blk);
}

//______________________________________________________
void AliMUONRawStreamTracker::AddErrorMessage()
{
  /// add message into logger of AliRawReader per event

    assert( GetReader() != 0 );
    TString msg;
    Int_t occurance = 0;
    AliMUONLogger* log = fPayload->GetErrorLogger();
    
    log->ResetItr();
    while(log->Next(msg, occurance))
    { 
      if (msg.Contains("Parity"))
         GetReader()->AddMinorErrorLog(kParityErr, msg.Data());

      if (msg.Contains("Glitch"))
        GetReader()->AddMajorErrorLog(kGlitchErr, msg.Data());

      if (msg.Contains("Padding"))
        GetReader()->AddMinorErrorLog(kPaddingWordErr, msg.Data());
    }
    
    log->Clear(); // clear logger after each event
}

//______________________________________________________
Bool_t AliMUONRawStreamTracker::IsErrorMessage() const
{
  /// true if there is any error/warning 
  if (GetPayLoad()->GetParityErrors() || 
        GetPayLoad()->GetGlitchErrors() || 
        GetPayLoad()->GetPaddingErrors())
    return kTRUE;
  
  return kFALSE;
}  
    
    
    
    

