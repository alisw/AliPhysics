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
#include "AliMUONRawStreamTrackerHP.h"
#include "AliMUONRawStreamTriggerHP.h"
#include "AliMUONRegHeader.h"
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
AliMUONDigitMaker::AliMUONDigitMaker(Bool_t enableErrorLogger, Bool_t a, Bool_t b) :
TObject(),
fScalerEvent(kFALSE),
fMakeTriggerDigits(kFALSE),
fRawStreamTracker(new AliMUONRawStreamTrackerHP),
fRawStreamTrigger(new AliMUONRawStreamTriggerHP),
fDigitStore(0x0),
fTriggerStore(0x0),
fLogger(new AliMUONLogger(10000)){
  /// ctor 
  
  if  ( !a || !b ) AliFatal("no longer supported");
  
  AliDebug(1,"");
  
  // Standard Constructor
  if (enableErrorLogger) 
  {
    fRawStreamTracker->EnabbleErrorLogger();
    fRawStreamTrigger->EnabbleErrorLogger();
  }
  else
  {
    fRawStreamTracker->DisableWarnings();
  }
  
  SetMakeTriggerDigits();
  
}

//__________________________________________________________________________
AliMUONDigitMaker::AliMUONDigitMaker(Bool_t enableErrorLogger) :
TObject(),
    fScalerEvent(kFALSE),
    fMakeTriggerDigits(kFALSE),
    fRawStreamTracker(new AliMUONRawStreamTrackerHP),
    fRawStreamTrigger(new AliMUONRawStreamTriggerHP),
    fDigitStore(0x0),
    fTriggerStore(0x0),
  fLogger(new AliMUONLogger(10000))
{
  /// ctor 

  AliDebug(1,"");
  
  // Standard Constructor
  if (enableErrorLogger) 
  {
    fRawStreamTracker->EnabbleErrorLogger();
    fRawStreamTrigger->EnabbleErrorLogger();
  }
  else
  {
    fRawStreamTracker->DisableWarnings();
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
  
  AliCodeTimerAuto("",0);

  // elex info
  Int_t    buspatchId;
  UChar_t  channelId;
  UShort_t manuId;
  UShort_t charge; 

  fRawStreamTracker->SetReader(rawReader);
  fRawStreamTracker->First();
  
  while ( fRawStreamTracker->Next(buspatchId,manuId,channelId,charge,kTRUE) )
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
  /// reading tracker DDL like ReadTriggerDDL but with fast decoder interface.
  /// filling the fTriggerStore container, which must not be null

  const AliMUONRawStreamTriggerHP::AliHeader*          darcHeader  = 0x0;
  const AliMUONRawStreamTriggerHP::AliRegionalHeader*  regHeader   = 0x0;
  const AliMUONRawStreamTriggerHP::AliLocalStruct*     localStruct = 0x0;

  Int_t loCircuit;

  fRawStreamTrigger->SetReader(rawReader);

  while (fRawStreamTrigger->NextDDL())
  {
    darcHeader = fRawStreamTrigger->GetHeaders();
    
    // fill global trigger information
    if (fTriggerStore) 
    {
      if (darcHeader->GetGlobalFlag()) 
      {
          AliMUONGlobalTrigger globalTrigger;
          globalTrigger.SetFromGlobalResponse(darcHeader->GetGlobalOutput());
          globalTrigger.SetFromGlobalInput(darcHeader->GetGlobalHeader()->fInput);
          fTriggerStore->SetGlobal(globalTrigger);
      }
    }
    
    Int_t nReg = fRawStreamTrigger->GetRegionalHeaderCount();
    
    for(Int_t iReg = 0; iReg < nReg ;iReg++)
    {   //reg loop
      

      // crate info  
      AliMpTriggerCrate* crate = AliMpDDLStore::Instance()->
                                GetTriggerCrate(fRawStreamTrigger->GetDDL(), iReg);
      
      if (!crate) 
        fLogger->Log(Form("Missing crate number %d in DDL %d\n", iReg, fRawStreamTrigger->GetDDL()));
     
      
      regHeader =  fRawStreamTrigger->GetRegionalHeader(iReg);
      
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
Int_t AliMUONDigitMaker::TriggerDigits(Int_t nBoard, 
                                       const TArrayS* xyPattern,
                                       AliMUONVDigitStore& digitStore) const
{
  /// make digits for trigger from pattern, and add them to digitStore

  AliCodeTimerAuto("",0);
  
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
            if (iCath && localBoard->GetSwitch(AliMpLocalBoard::kZeroAllYLSB)) offset = -8;
            
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

//______________________________________________________________________________
void 
AliMUONDigitMaker::SetTryRecover(Bool_t flag) 
{
  /// Instruct the decoder to try to recover corrupted raw data.
  /// Only use for specific cases for which you know it will work...
  fRawStreamTracker->TryRecover(flag);
}
