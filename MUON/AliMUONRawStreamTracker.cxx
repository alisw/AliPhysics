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

/* $Id $ */


///////////////////////////////////////////////////////////////////////////////
///
/// \class AliMUONRawStreamTracker
/// This class provides access to MUON digits in raw data.
///
/// It loops over all MUON digits in the raw data given by the AliRawReader.
/// The Next method goes to the next digit. If there are no digits left
/// it returns kFALSE (under develpment)
/// It can loop also over DDL and store the decoded rawdata in TClonesArray
/// in Payload class.
/// 
/// Version implement for Tracker
///
/// \author Christian Finck & Laurent Aphecetche
///////////////////////////////////////////////////////////////////////////////

#include "AliMUONRawStreamTracker.h"

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

/// \cond CLASSIMP
ClassImp(AliMUONRawStreamTracker)
/// \endcond

AliMUONRawStreamTracker::AliMUONRawStreamTracker()
: TObject(),
  fRawReader(0x0),
  fDDL(0),
  fMaxDDL(20),
  fPayload(new AliMUONPayloadTracker()),
  fCurrentDDL(0),
  fCurrentDDLIndex(fMaxDDL),
  fCurrentBlockHeader(0),
  fCurrentBlockHeaderIndex(0),
  fCurrentDspHeader(0),
  fCurrentDspHeaderIndex(0),
  fCurrentBusStruct(0),
  fCurrentBusStructIndex(0),
  fCurrentDataIndex(0),
  fEnableErrorLogger(kFALSE)
{
  ///
  /// create an object to read MUON raw digits
  /// Default ctor for monitoring purposes
  ///
  
  
}

//_________________________________________________________________
AliMUONRawStreamTracker::AliMUONRawStreamTracker(AliRawReader* rawReader)
: TObject(),
  fRawReader(rawReader),
  fDDL(0),
  fMaxDDL(20),
  fPayload(new AliMUONPayloadTracker()),
  fCurrentDDL(0L),
  fCurrentDDLIndex(fMaxDDL),
  fCurrentBlockHeader(0),
  fCurrentBlockHeaderIndex(0),
  fCurrentDspHeader(0),
  fCurrentDspHeaderIndex(0),
  fCurrentBusStruct(0),
  fCurrentBusStructIndex(0),
  fCurrentDataIndex(0),
  fEnableErrorLogger(kFALSE)
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
  /// Initialize the iteration process
  
  fCurrentDDLIndex = -1;
  fCurrentDspHeaderIndex = -1;
  fCurrentBusStructIndex = -1;

  fCurrentDspHeader = 0;
  fCurrentBusStruct = 0;
  
  // Find the first non-empty structure
  GetNextDDL();
  GetNextBlockHeader();
  GetNextDspHeader();
  GetNextBusStruct();
}

//______________________________________________________
Bool_t
AliMUONRawStreamTracker::GetNextDDL()
{
  /// Returns the next DDL present
  
  Bool_t kFound(kFALSE);
  
  while ( fCurrentDDLIndex < fMaxDDL-1 && !kFound ) 
  {
    ++fCurrentDDLIndex;
    fRawReader->Reset();
    fRawReader->Select("MUONTRK",fCurrentDDLIndex,fCurrentDDLIndex);
    if ( fRawReader->ReadHeader() ) 
    {
      kFound = kTRUE;
    }
  }
  
  if ( !kFound ) 
  {
    fCurrentDDLIndex = fMaxDDL;
    return kFALSE;
  }
  
  Int_t totalDataWord  = fRawReader->GetDataSize(); // in bytes
  
  AliDebug(3, Form("DDL Number %d totalDataWord %d\n", fCurrentDDLIndex,
                   totalDataWord));
  
  UInt_t *buffer = new UInt_t[totalDataWord/4];
  
  if ( !fRawReader->ReadNext((UChar_t*)buffer, totalDataWord) )
  {
    fCurrentDDL = 0;
    return kFALSE;
  }
  fPayload->ResetDDL();
  
  Bool_t ok = fPayload->Decode(buffer, totalDataWord/4);
  
  if (fEnableErrorLogger) AddErrorMessage();

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
  
  fPayload->ResetDDL();
  
  while ( fDDL < 20 ) 
  {
    fRawReader->Reset();
    fRawReader->Select("MUONTRK", fDDL, fDDL);  //Select the DDL file to be read  
    if (fRawReader->ReadHeader()) break;
    AliDebug(3,Form("Skipping DDL %d which does not seem to be there",fDDL));
    ++fDDL;
  }
  
  if ( fDDL == 20 ) 
  {
    fDDL = 0;
    return kFALSE;
  }
  
  AliDebug(3, Form("DDL Number %d\n", fDDL ));
  
  Int_t totalDataWord  = fRawReader->GetDataSize(); // in bytes
  
  UInt_t *buffer = new UInt_t[totalDataWord/4];
  
  if(!fRawReader->ReadNext((UChar_t*)buffer, totalDataWord))
  {
    delete[] buffer;
    return kFALSE;
  }
  
  Bool_t ok = fPayload->Decode(buffer, totalDataWord/4);

  if (fEnableErrorLogger) AddErrorMessage();

  delete[] buffer;
  
  fDDL++;
  
  return ok;
}

//______________________________________________________
void AliMUONRawStreamTracker::SetMaxDDL(Int_t ddl) 
{
  /// set DDL number
  if (ddl > 20) ddl = 20;
  fMaxDDL = ddl;
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

    for (Int_t i = 0; i < fPayload->GetParityErrors(); ++i)
	fRawReader->AddMinorErrorLog(kParityErr, Form("Parity error for buspatch %s",  
						      fPayload->GetParityErrBus()[i]));

    for (Int_t i = 0; i < fPayload->GetGlitchErrors(); ++i)
	fRawReader->AddMajorErrorLog(kGlitchErr, "Glitch error occurs skip event");

    for (Int_t i = 0; i < fPayload->GetPaddingErrors(); ++i)
	fRawReader->AddMinorErrorLog(kPaddingWordErr, "Padding word error");

}
