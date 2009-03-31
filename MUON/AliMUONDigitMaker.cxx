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

// $Id$

//-----------------------------------------------------------------------------
/// \class AliMUONDigitMaker
/// MUON Digit maker from rawdata.
///
/// Raw2Digits:
/// Using real mapping  for tracker
/// Indranil Das (Adapted for runloader: Ch. Finck) july 05
///
/// Implemented non-constant buspatch numbers for tracking
/// with correct DDL id.
/// (Ch. Finck, dec 05)
///
/// Add reader for scaler trigger events
/// Use memcpy instead of assignment elt by elt
/// (Ch. Finck, Jan 06)
///
/// Using new interface with AliMUONRawStreamTracker(Trigger)
/// (New interface of AliMUONRawReader class)
/// (further details could be found in Alice-note)
/// (Ch. Finck, March 06)
///
/// Add (S)Digit maker tracker (for free)
/// and for trigger. Create trigger inverse mapping.
///
/// \author Ch. Finck, oct 06 
//-----------------------------------------------------------------------------

#include "AliMUONDigitMaker.h"

#include "AliMUONDDLTrigger.h"
#include "AliMUONDarcHeader.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONLocalStruct.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONLogger.h"
#include "AliMUONRawStreamTracker.h"
#include "AliMUONRawStreamTrackerHP.h"
#include "AliMUONRawStreamTrigger.h"
#include "AliMUONRawStreamTriggerHP.h"
#include "AliMUONRegHeader.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONVTriggerStore.h"
#include "AliMpDetElement.h"
#include "AliMpTriggerCrate.h"
#include "AliMpLocalBoard.h"
#include "AliMpCathodType.h"
#include "AliMpDDLStore.h"
#include "AliMpDEManager.h"
#include "AliMpPad.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliCodeTimer.h"
#include "AliLog.h"
#include "AliRawReader.h"
#include <TArrayS.h>

/// \cond CLASSIMP
ClassImp(AliMUONDigitMaker) // Class implementation in ROOT context
/// \endcond

//__________________________________________________________________________
AliMUONDigitMaker::AliMUONDigitMaker(
      Bool_t enableErrorLogger,
      Bool_t useFastTrackerDecoder, Bool_t useFastTriggerDecoder
  ) :
    TObject(),
    fScalerEvent(kFALSE),
    fMakeTriggerDigits(kFALSE),
    fRawStreamTracker(NULL),
    fRawStreamTrigger(NULL),
    fDigitStore(0x0),
    fTriggerStore(0x0),
  fLogger(new AliMUONLogger(10000))
{
  /// ctor 

  AliDebug(1,"");
  
  CreateRawStreamTracker(useFastTrackerDecoder);
  CreateRawStreamTrigger(useFastTriggerDecoder);

  // Standard Constructor
  if (enableErrorLogger) {
    fRawStreamTracker->EnabbleErrorLogger();
    fRawStreamTrigger->EnabbleErrorLogger();
  }

  SetMakeTriggerDigits();

}

//__________________________________________________________________________
AliMUONDigitMaker::~AliMUONDigitMaker()
{
  /// clean up
  /// and time processing measure

  delete fRawStreamTracker;
  delete fRawStreamTrigger;
  delete fLogger;
}

//__________________________________________________________________________
void AliMUONDigitMaker::CreateRawStreamTracker(Bool_t useFastDecoder)
{
/// Create raw stream tracker according to the passed option

  if (useFastDecoder)
  {
    fRawStreamTracker = new AliMUONRawStreamTrackerHP();
  }
  else {
    AliInfo("Using non-high performance tracker decoder.");
    fRawStreamTracker = new AliMUONRawStreamTracker();
  }  
}

//__________________________________________________________________________
void AliMUONDigitMaker::CreateRawStreamTrigger(Bool_t useFastDecoder)
{
/// Create raw stream trigger according to the passed option

  if (useFastDecoder)
  {
    fRawStreamTrigger = new AliMUONRawStreamTriggerHP();
  }
  else {
    AliInfo("Using non-high performance tracker decoder.");
    fRawStreamTrigger = new AliMUONRawStreamTrigger();
  }  
}

//____________________________________________________________________
void
AliMUONDigitMaker::Print(Option_t*) const
{
  /// Printout

  cout << "RawStreamerTracker class=" << fRawStreamTracker->ClassName()
       << " MakeTriggerDigits=" << fMakeTriggerDigits
       << " ScalerEvent=" << fScalerEvent
       << " DigitStore=" << fDigitStore
       << " TriggerStore=" << fTriggerStore << endl;

  if ( fLogger ) fLogger->Print();
}

//____________________________________________________________________
Int_t
AliMUONDigitMaker::Raw2Digits(AliRawReader* rawReader, 
                                        AliMUONVDigitStore* digitStore,
                                        AliMUONVTriggerStore* triggerStore)
{
  /// Main method to creates digit
  /// for tracker 
  /// and trigger

  AliDebug(1,Form("rawReader=%p digitStore=%p triggerStore=%p",
                  rawReader,digitStore,triggerStore));
  
  fDigitStore = digitStore;
  fTriggerStore = triggerStore;
  
  if (!fDigitStore && !fTriggerStore)
  {
    fLogger->Log("No digit or trigger store given. Nothing to do...");
    return kTriggerBAD & kTrackerBAD;
  }
  
  Int_t tracker(kOK);
  Int_t trigger(kOK);
  
  if ( fDigitStore ) 
  {
    fDigitStore->Clear(); // insure we start with an empty container
    tracker = ReadTrackerDDL(rawReader);
  }
  
  if ( fTriggerStore || fMakeTriggerDigits ) 
  {
    if ( fTriggerStore ) fTriggerStore->Clear();
    if ( fMakeTriggerDigits && !fDigitStore ) 
    {
      fLogger->Log("Asking for trigger digits but digitStore is null");
    }
    else
    {
      trigger = ReadTriggerDDL(rawReader);
    }
  }
  
  return tracker | trigger;
}

//____________________________________________________________________
Int_t
AliMUONDigitMaker::ReadTrackerDDL(AliRawReader* rawReader)
{
  /// Reading tracker DDL
  /// filling the fDigitStore container, which must not be null

  AliDebug(1,"");
  
  AliCodeTimerAuto("");

  // elex info
  Int_t    buspatchId;
  UChar_t  channelId;
  UShort_t manuId;
  UShort_t charge; 

  fRawStreamTracker->SetReader(rawReader);
  fRawStreamTracker->First();
  
  while ( fRawStreamTracker->Next(buspatchId,manuId,channelId,charge) )
  {    
    // getting DE from buspatch
    Int_t detElemId = AliMpDDLStore::Instance()->GetDEfromBus(buspatchId);

    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);

    if (!de)
      {
	fLogger->Log(Form("DE %04d does not exist !"));
	continue;
      }

    if (!de->IsConnectedChannel(manuId,channelId))
      {
	// non connected pad, do nothing (this is not an error !)
	continue;
      }

    const AliMpVSegmentation* seg 
      = AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId, 
                                                                      manuId);  

    if (!seg)
    {
      fLogger->Log(Form("(DE,MANUID)=(%04d,%04d) is not valid",detElemId,manuId));
      continue;
    }
    
    AliMp::CathodType cathodeType = de->GetCathodType(seg->PlaneType());

    AliMpPad pad = seg->PadByLocation(manuId,channelId,kFALSE);

    if (!pad.IsValid())
    {
      fLogger->Log(Form("No pad for detElemId: %d, manuId: %d, channelId: %d",
                    detElemId, manuId, channelId));
      continue;
    } 

    AliMUONVDigit* digit = fDigitStore->Add(detElemId,manuId,channelId,cathodeType,
                                            AliMUONVDigitStore::kDeny);

    if (!digit)
    {
      fLogger->Log(Form("Digit DE %04d Manu %04d Channel %02d could not be added",
                    detElemId, manuId, channelId));
      continue;
    }
    
    digit->SetPadXY(pad.GetIx(),pad.GetIy());
    
    digit->SetADC(charge);

  }
  
  if ( fRawStreamTracker->IsErrorMessage() ) 
  {
    return kTrackerBAD;
  }
  
  return kOK;
}

//____________________________________________________________________
Int_t
AliMUONDigitMaker::ReadTriggerDDL(AliRawReader* rawReader)
{
  /// reading tracker DDL
  /// filling the fTriggerStore container, which must not be null

  AliDebug(1,"");
  
  AliMUONDDLTrigger*       ddlTrigger      = 0x0;
  AliMUONDarcHeader*       darcHeader      = 0x0;
  AliMUONRegHeader*        regHeader       = 0x0;
  AliMUONLocalStruct*      localStruct     = 0x0;

  Int_t loCircuit;

  AliCodeTimerAuto("");
  
  if (UsingFastTriggerDecoder()) return ReadTriggerDDLFast(rawReader);

  fRawStreamTrigger->SetReader(rawReader);

  while (fRawStreamTrigger->NextDDL()) 
  {
    ddlTrigger = fRawStreamTrigger->GetDDLTrigger();
    darcHeader = ddlTrigger->GetDarcHeader();
    
    // fill global trigger information
    if (fTriggerStore) 
    {
      if (darcHeader->GetGlobalFlag()) 
      {
          AliMUONGlobalTrigger globalTrigger;
          globalTrigger.SetFromGlobalResponse(darcHeader->GetGlobalOutput());
          fTriggerStore->SetGlobal(globalTrigger);
      }
    }
    
    Int_t nReg = darcHeader->GetRegHeaderEntries();
    
    for(Int_t iReg = 0; iReg < nReg ;iReg++)
    {   //reg loop
      

      // crate info  
      AliMpTriggerCrate* crate = AliMpDDLStore::Instance()->
                                GetTriggerCrate(fRawStreamTrigger->GetDDL(), iReg);
      
      if (!crate) 
        fLogger->Log(Form("Missing crate number %d in DDL %d\n", iReg, fRawStreamTrigger->GetDDL()));
     
      
      regHeader =  darcHeader->GetRegHeaderEntry(iReg);
      
      Int_t nLocal = regHeader->GetLocalEntries();
      for(Int_t iLocal = 0; iLocal < nLocal; iLocal++) 
      {  
        
        localStruct = regHeader->GetLocalEntry(iLocal);
        
        // if card exist
        if (localStruct) {
          
   	  loCircuit = crate->GetLocalBoardId(localStruct->GetId());

	  if ( !loCircuit ) continue; // empty slot

	  AliMpLocalBoard* localBoard = AliMpDDLStore::Instance()->GetLocalBoard(loCircuit, false);

	  // skip copy cards
	  if( !localBoard->IsNotified()) 
	     continue;
          
          if (fTriggerStore) 
          {
            // fill local trigger
            AliMUONLocalTrigger localTrigger;
            localTrigger.SetLocalStruct(loCircuit, *localStruct);
            fTriggerStore->Add(localTrigger);
          } 
          
          if ( fMakeTriggerDigits )
          {
            //FIXEME should find something better than a TArray
            TArrayS xyPattern[2];
            
	    localStruct->GetXPattern(xyPattern[0]);
	    localStruct->GetYPattern(xyPattern[1]);
            
            TriggerDigits(loCircuit, xyPattern, *fDigitStore);
          }          
        } // if triggerY
      } // iLocal
    } // iReg
  } // NextDDL
  
  return kOK;
}

//____________________________________________________________________
Int_t
AliMUONDigitMaker::ReadTriggerDDLFast(AliRawReader* rawReader)
{
  /// reading tracker DDL like ReadTriggerDDL but with fast decoder interface.
  /// filling the fTriggerStore container, which must not be null

  const AliMUONRawStreamTriggerHP::AliHeader*          darcHeader  = 0x0;
  const AliMUONRawStreamTriggerHP::AliRegionalHeader*  regHeader   = 0x0;
  const AliMUONRawStreamTriggerHP::AliLocalStruct*     localStruct = 0x0;

  Int_t loCircuit;

  fRawStreamTrigger->SetReader(rawReader);
  AliMUONRawStreamTriggerHP* rawStreamTrigger =
    dynamic_cast<AliMUONRawStreamTriggerHP*>(fRawStreamTrigger);

  while (fRawStreamTrigger->NextDDL())
  {
    darcHeader = rawStreamTrigger->GetHeaders();
    
    // fill global trigger information
    if (fTriggerStore) 
    {
      if (darcHeader->GetGlobalFlag()) 
      {
          AliMUONGlobalTrigger globalTrigger;
          globalTrigger.SetFromGlobalResponse(darcHeader->GetGlobalOutput());
          fTriggerStore->SetGlobal(globalTrigger);
      }
    }
    
    Int_t nReg = rawStreamTrigger->GetRegionalHeaderCount();
    
    for(Int_t iReg = 0; iReg < nReg ;iReg++)
    {   //reg loop
      

      // crate info  
      AliMpTriggerCrate* crate = AliMpDDLStore::Instance()->
                                GetTriggerCrate(fRawStreamTrigger->GetDDL(), iReg);
      
      if (!crate) 
        fLogger->Log(Form("Missing crate number %d in DDL %d\n", iReg, fRawStreamTrigger->GetDDL()));
     
      
      regHeader =  rawStreamTrigger->GetRegionalHeader(iReg);
      
      Int_t nLocal = regHeader->GetLocalStructCount();
      for(Int_t iLocal = 0; iLocal < nLocal; iLocal++) 
      {
        
        localStruct = regHeader->GetLocalStruct(iLocal);
        
        // if card exist
        if (localStruct) {
          
   	  loCircuit = crate->GetLocalBoardId(localStruct->GetId());

	  if ( !loCircuit ) continue; // empty slot

	  AliMpLocalBoard* localBoard = AliMpDDLStore::Instance()->GetLocalBoard(loCircuit, kTRUE);

	  // skip copy cards
	  if( !localBoard->IsNotified()) 
	     continue;
          
          if (fTriggerStore)
          {
            // fill local trigger
            AliMUONLocalTrigger localTrigger;
            localTrigger.SetLoCircuit(loCircuit);
            localTrigger.SetLoStripX((Int_t)localStruct->GetXPos());
            localTrigger.SetLoStripY((Int_t)localStruct->GetYPos());
            localTrigger.SetLoDev((Int_t)localStruct->GetXDev());
            localTrigger.SetLoSdev((Int_t)localStruct->GetSXDev());
            localTrigger.SetLoTrigY((Int_t)localStruct->GetTrigY());
            localTrigger.SetLoLpt(localStruct->GetLpt());
            localTrigger.SetLoHpt(localStruct->GetHpt());
            localTrigger.SetX1Pattern(localStruct->GetX1());
            localTrigger.SetX2Pattern(localStruct->GetX2());
            localTrigger.SetX3Pattern(localStruct->GetX3());
            localTrigger.SetX4Pattern(localStruct->GetX4());
            localTrigger.SetY1Pattern(localStruct->GetY1());
            localTrigger.SetY2Pattern(localStruct->GetY2());
            localTrigger.SetY3Pattern(localStruct->GetY3());
            localTrigger.SetY4Pattern(localStruct->GetY4());
            fTriggerStore->Add(localTrigger);
          }
          
          if ( fMakeTriggerDigits )
          {
            //FIXEME should find something better than a TArray
            TArrayS xyPattern[2];
            
	    localStruct->GetXPattern(xyPattern[0]);
	    localStruct->GetYPattern(xyPattern[1]);
            
            TriggerDigits(loCircuit, xyPattern, *fDigitStore);
          }
        } // if triggerY
      } // iLocal
    } // iReg
  } // NextDDL
  
  return kOK;
}

//____________________________________________________________________
Int_t AliMUONDigitMaker::TriggerDigits(Int_t nBoard, 
                                       TArrayS* xyPattern,
                                       AliMUONVDigitStore& digitStore) const
{
  /// make digits for trigger from pattern, and add them to digitStore

  AliCodeTimerAuto("");
  
  Int_t detElemId;

  AliMpLocalBoard* localBoard = AliMpDDLStore::Instance()->GetLocalBoard(nBoard);

  Int_t n,b;

  // loop over x1-4 and y1-4
  for (Int_t iChamber = 0; iChamber < 4; ++iChamber)
  {
    for (Int_t iCath = 0; iCath < 2; ++iCath)
    {
      Int_t pattern = (Int_t)xyPattern[iCath].At(iChamber); 
      if (!pattern) continue;
      
      // get detElemId
      detElemId = AliMpDDLStore::Instance()->GetDEfromLocalBoard(nBoard, iChamber);
        
        const AliMpVSegmentation* seg 
          = AliMpSegmentation::Instance()
          ->GetMpSegmentation(detElemId, AliMp::GetCathodType(iCath));  
        
        // loop over the 16 bits of pattern
        for (Int_t ibitxy = 0; ibitxy < 16; ++ibitxy) 
        {
          if ((pattern >> ibitxy) & 0x1) 
          {            
            // not quite sure about this
            Int_t offset = 0;
            if (iCath && localBoard->GetSwitch(6)) offset = -8;
            
            AliMpPad pad = seg->PadByLocation(nBoard,ibitxy+offset,kTRUE);
                        
            if (!pad.IsValid()) 
            {
              fLogger->Log(Form("No pad for detElemId: %d, nboard %d, ibitxy: %d\n",
                              detElemId, nBoard, ibitxy));
              continue;
            }

            n = pad.GetLocalBoardId(0); // always take first location so that digits are not inserted several times
	    b = pad.GetLocalBoardChannel(0);

	    AliDebug(1,Form("Using localBoard %d ixy %d instead of %d,%d",
			    n,b,nBoard,ibitxy));

	    AliMUONVDigit* digit = digitStore.Add(detElemId,n,b,iCath,AliMUONVDigitStore::kDeny);
            
            if (!digit)
            {
		AliDebug(1, Form("Digit DE %04d LocalBoard %03d ibitxy %02d cath %d already in store",
				 detElemId,nBoard,ibitxy,iCath));
		continue;
            }
            
            Int_t padX = pad.GetIx();
            Int_t padY = pad.GetIy();
            
            // fill digit
            digit->SetPadXY(padX,padY);
            digit->SetCharge(1.);
          }// xyPattern
        }// ibitxy
    }// cath
  } // ichamber
  
  return kTRUE;
}

//______________________________________________________________________________
Bool_t 
AliMUONDigitMaker::TriggerToDigitsStore(const AliMUONVTriggerStore& triggerStore,
					AliMUONVDigitStore& digitStore) const
{
  //
  /// make (S)Digit for trigger
  //
  
  digitStore.Clear();
  
  AliMUONLocalTrigger* locTrg;
  TIter next(triggerStore.CreateLocalIterator());
  
  while ( ( locTrg = static_cast<AliMUONLocalTrigger*>(next()) ) ) 
  {
    if (locTrg->IsNull()) continue;
   
    TArrayS xyPattern[2];
    locTrg->GetXPattern(xyPattern[0]);
    locTrg->GetYPattern(xyPattern[1]);

    Int_t nBoard = locTrg->LoCircuit();
    TriggerDigits(nBoard, xyPattern, digitStore);
  }
  return kTRUE;
}

//____________________________________________________________________
Bool_t AliMUONDigitMaker::UsingFastTrackerDecoder() const
{
/// Returns kTRUE if the digit maker is using the high performance decoder for
/// tracker DDL stream decoding.

  return (fRawStreamTracker->IsA() == AliMUONRawStreamTrackerHP::Class());
}

//____________________________________________________________________
Bool_t AliMUONDigitMaker::UsingFastTriggerDecoder() const
{
/// Returns kTRUE if the digit maker is using the high performance decoder for
/// trigger DDL stream decoding.

  return (fRawStreamTrigger->IsA() == AliMUONRawStreamTriggerHP::Class());
}

//____________________________________________________________________
void  AliMUONDigitMaker::SetFastTrackerDecoder(Bool_t useFastDecoder)
{
/// Set fast raw data decoder

  delete fRawStreamTracker;
  CreateRawStreamTracker(useFastDecoder);
}

//____________________________________________________________________
void  AliMUONDigitMaker::SetFastTriggerDecoder(Bool_t useFastDecoder)
{
/// Set fast raw data decoder

  delete fRawStreamTrigger;
  CreateRawStreamTrigger(useFastDecoder);
}

