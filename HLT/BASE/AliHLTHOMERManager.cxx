//-*- Mode: C++ -*-
// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTHOMERManager.cxx
    @author Jochen Thaeder
    @date
    @brief  Manger for HOMER in aliroot
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#if __GNUC__>= 3
   using namespace std;
#endif

#define EVE_DEBUG 0

#include "AliHLTHOMERManager.h"
// -- -- -- -- -- -- -- 
#include "AliHLTHOMERLibManager.h"
#include "AliHLTHOMERSourceDesc.h"
#include "AliHLTHOMERBlockDesc.h"
// -- -- -- -- -- -- -- 
#include "AliHLTGlobalTriggerDecision.h"
#include "AliHLTTriggerDecision.h"
//---------------------------

ClassImp(AliHLTHOMERManager)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
  AliHLTHOMERManager::AliHLTHOMERManager() :
  fLibManager(new AliHLTHOMERLibManager),
  fStateHasChanged(kTRUE),
  fProxyHandler(NULL),
  fCurrentReader(NULL),
  fReaderList(NULL),
  fSourceList(NULL),
  fNBlks(0),
  fEventID(),
  fEventId(-1),
  fCurrentBlk(0),
  fAsyncBlockList(NULL),
  fBlockList(NULL),
   fEventBuffer(NULL),
  fBufferTopIdx(-1),
  fBufferLowIdx(-1),
  fCurrentBufferIdx(-1),
  fNavigateBufferIdx(-1),
  fNEventsAvailable(0),
  fConnected(kFALSE), 
  fTriggerString("ALL"), 
  fNEventsNotTriggered(0),
  fRetryNextEvent(kFALSE),
  fIsBlockOwner(kTRUE)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

//##################################################################################
AliHLTHOMERManager::~AliHLTHOMERManager() {
  // see header file for class documentation

  if ( fLibManager ) {

    if ( fReaderList ) {
      TIter next(fReaderList);
      TObject * object = NULL;
      while ( ( object = next()) )
	fLibManager->DeleteReader(static_cast<AliHLTHOMERReader*>(object) );
      
      fReaderList->Clear();
      delete fReaderList;
    }
    fReaderList = NULL;   
    
    delete fLibManager;
  } 
  fLibManager = NULL;

  if ( fProxyHandler != NULL )
    delete fProxyHandler;
  fProxyHandler = NULL;

  if ( fSourceList != NULL )
    delete fSourceList;
  fSourceList = NULL;

  if ( fEventBuffer ) {
    fEventBuffer->Clear();
    delete fEventBuffer;
  }
  fEventBuffer = NULL;

  if(fBlockList) {
    fBlockList->Clear();
    delete fBlockList;
  }
  fBlockList = NULL;

  if ( fAsyncBlockList ) {
    fAsyncBlockList->Delete();
    delete fAsyncBlockList;
  }
  fAsyncBlockList = NULL;
}

//##################################################################################
Int_t AliHLTHOMERManager::Initialize() {
  // see header file for class documentation

  HLTInfo("Initializing");

  Int_t iResult = 0;

  // -- Initialize ProxyHandler
  if ( !fProxyHandler )
    fProxyHandler = new AliHLTHOMERProxyHandler();
  
  if ( fProxyHandler ) {
    iResult = fProxyHandler->Initialize();
    if (iResult)
      HLTError(Form("Initialize of ProxyHandler failed."));
  
  } else {
    iResult = -1;
    HLTError(Form("Creating of ProxyHandler failed."));
  }
 
  // -- Initialize ReaderList
  //    List ist not owner, as reader have to be created/deleted by the LibManager
  if( !fReaderList )
    fReaderList = new TList();
  
  // -- Initialize asynchronous BlockList
  if( !fAsyncBlockList ) {
    fAsyncBlockList = new TList();
    fAsyncBlockList->SetOwner(kTRUE);
  }

  //initialize normal block list
  if( !fBlockList ) {
    fBlockList = new TList();
    fBlockList->SetOwner(kFALSE);
  }

  // -- Initialize Event Buffer and EventID array
  if ( !fEventBuffer ) {
    fEventBuffer = new TClonesArray( "TList", BUFFERSIZE );
  }

  for ( Int_t idx = 0; idx < BUFFERSIZE; ++idx ) {
    new ((*fEventBuffer)[idx]) TList( );
    (reinterpret_cast<TList*>((*fEventBuffer)[idx]))->SetOwner(kTRUE);
    
    fEventID[idx] = 0;
  }

  return iResult;
}

/*
 * ---------------------------------------------------------------------------------
 *                                 Source Handling
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Int_t AliHLTHOMERManager::CreateSourcesList() {
  // see header file for class documentation

  Int_t iResult = 0;

  if ( fSourceList != NULL )
    delete fSourceList;
  fSourceList = NULL;

  fSourceList = new TList();
  fSourceList->SetOwner( kTRUE );

  iResult = fProxyHandler->FillSourceList( fSourceList );
  if ( iResult < 0 ) {
    HLTWarning(Form("There have been errors, while creating the sources list."));
  }
  else if ( iResult > 0 ) {
    HLTWarning(Form("No active services found."));
  }
  else if ( fSourceList->IsEmpty() ) {
    HLTWarning(Form("No active services in the list."));
    iResult = 2;
  }
  else {
     HLTInfo(Form("New sources list created."));

    // -- New SourceList has been created 
    // --> All Sources are new --> State has changed
    fStateHasChanged = kTRUE;
  }

  return iResult;
}

//##################################################################################
void AliHLTHOMERManager::SetSourceState( AliHLTHOMERSourceDesc * source, Bool_t state ) {
  // see header file for class documentation

  if ( source->IsSelected() != state ) {
    source->SetState( state );
    fStateHasChanged = kTRUE;
  }

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                         Connection Handling - public
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Int_t AliHLTHOMERManager::ConnectHOMER( TString detector ){
  // see header file for class documentation

  Int_t iResult = 0;

  // HAck Jochen
  //----
  detector="ALL";

  // -- Check if LibManager is present
  if ( ! fLibManager ) {
    HLTError(Form("No LibManager present."));
    return -1;
  }
  
  // -- Check if already connected and state has not changed
  if ( fStateHasChanged == kFALSE && IsConnected() ) {
    HLTInfo(Form("No need for reconnection."));
    return 0;
  }
  
  // -- If already connected, disconnect before connect
  //    or if ReaderList already filled
  if ( IsConnected() || fReaderList->GetSize() != 0 ) {
    HLTInfo(Form("IsConnected: %d      fReaderList.Size:   %d", IsConnected(), fReaderList->GetSize()));
    DisconnectHOMER();
  }
  // -- Create the Readoutlist
  UShort_t* sourcePorts = new UShort_t [fSourceList->GetEntries()];
  const Char_t ** sourceHostnames = new const Char_t* [fSourceList->GetEntries()];
  for(Int_t i = 0; i < fSourceList->GetEntries(); i++) {
   sourceHostnames[i] = "";
  }
  UInt_t sourceCount = 0;

  CreateReadoutList( sourceHostnames, sourcePorts, sourceCount, detector );
  if ( sourceCount == 0 ) {
    HLTError(Form("No sources selected, aborting."));
    delete [] sourcePorts;
    delete [] sourceHostnames;
    return -2;
  }

  // ***
  // *** Connect to data sources
  // ***
  
  for (UInt_t idx = 0; idx < sourceCount; idx++) {
    
    if (sourcePorts[idx] > 60000)
      continue;

    HLTInfo(Form("Adding source %d as %s : %d", idx, sourceHostnames[idx], sourcePorts[idx]));
    
    fReaderList->Add(dynamic_cast<TObject*>(fLibManager->OpenReader(sourceHostnames[idx], sourcePorts[idx])));
    AliHLTHOMERReader *reader = static_cast<AliHLTHOMERReader*>(fReaderList->Last());
    if ( !reader ) {
      HLTError(Form("Adding reader failed, aborting"));
      delete [] sourcePorts;
      delete [] sourceHostnames;
      return -3;
    }

    if ( (iResult = reader->GetConnectionStatus()) )  {

      // -- Connection to source failed
      
      HLTError(Form("Error establishing connection to TCP source %s:%hu: %s (%d)",
		    sourceHostnames[idx], sourcePorts[idx], strerror(iResult), iResult));

      if( !(TString(sourceHostnames[idx]).CompareTo("localhost")) ) {
	HLTInfo("The failed connection is on localhost. is SSH tunnel up????? ");
	HLTInfo(Form("Do: 'ssh -L %d:alihlt-vobox0.cern.ch:%d cernUser@lxplus.cern.ch -fN'",
		     sourcePorts[idx], sourcePorts[idx]));
      }
      
      // -- Remove reader
      fReaderList->RemoveLast();

      if ( reader )
	fLibManager->DeleteReader( reader );
      reader = NULL;
      
      HLTInfo(Form("Removed source %d,  %s : %d from sourceList", idx, sourceHostnames[idx], sourcePorts[idx]));
      
    } 
    else {
      // -- Connection succeded
      fConnected = kTRUE;

      HLTInfo(Form("Connection established to source %s on port %d", sourceHostnames[idx], sourcePorts[idx]));
    }
    
  } // for (Int_t idx = 0; idx < sourceCount; idx++) {
  
  delete[] sourceHostnames;
  delete[] sourcePorts;

  return iResult;

}

//##################################################################################
void AliHLTHOMERManager::DisconnectHOMER(){
  // see header file for class documentation

  HLTInfo("Disconnecting");

  if ( fReaderList && fLibManager ) {
    HLTInfo("Deleting readerlist and libmanager");
    TIter next(fReaderList);
    TObject * object = NULL;
    while ( ( object = next()) ) 
      fLibManager->DeleteReader(static_cast<AliHLTHOMERReader*>(object) );
      

    HLTInfo(Form("fReaderList size %d", fReaderList->GetSize()));
    fReaderList->Clear();
    HLTInfo(Form("fReaderList size %d", fReaderList->GetSize()));
    delete fReaderList;
    fReaderList = new TList ();
    HLTInfo(Form("fReaderList size %d", fReaderList->GetSize()));
  }
  
  fStateHasChanged = kTRUE;
  fConnected = kFALSE;
  
  HLTInfo(Form("Connection closed."));

  return;
}

//##################################################################################
Int_t AliHLTHOMERManager::ReconnectHOMER( TString detector="" ){
  // see header file for class documentation
  
  Int_t iResult = 0;

  if ( IsConnected() )
    DisconnectHOMER();

  iResult = ConnectHOMER(detector);
  if ( iResult ) {
    HLTError(Form("Error reconnecting."));
  }

  return iResult;
}

/*
 * ---------------------------------------------------------------------------------
 *                            Event Handling - public
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Int_t AliHLTHOMERManager::NextEvent(){
 
  // see header file for class documentation
  
  Int_t iResult = 0;
  Int_t iRetryCount = 0;
  
  if ( !IsConnected() || fStateHasChanged ) {
    HLTInfo("Not connected or state has changed, returning to AliEveHOMERManager, which will deal with this situation");
    //    cout << "connectecd  " << IsConnected()  << "haschanged  "<<fStateHasChanged << endl;
    return 55;//ConnectHOMER();
  }
  if ( !IsConnected() ) {
    HLTWarning(Form( "Not connected yet." ));
    return -1;
  }

  // -- Reset asyncronous BlockList
  fAsyncBlockList->Delete();

  // ***
  // *** Loop over all readers and get new event data
  // ***
  
  TIter next(fReaderList);
  TObject * object = NULL;
  
  while( (object = next()) ) {
    
    fCurrentReader = static_cast<AliHLTHOMERReader*>(object);
    
    // -- Read next event data and error handling for HOMER (error codes and empty blocks)
    while ( 1 ) {
      
      iResult = fCurrentReader->ReadNextEvent( 40000000 /*timeout in us*/);
      
      if ( iResult == 111 || iResult == 32 || iResult == 6 ) {
	HLTError(Form("No connection to source %d: %s (%d)", 
		      fCurrentReader->GetErrorConnectionNdx(), strerror(iResult), iResult));
	break;
      } 
      else if ( iResult == 110 ) {
	HLTError(Form("Timeout occured, reading event from source %d: %s (%d)", 
		      fCurrentReader->GetErrorConnectionNdx(), strerror(iResult), iResult));
	break;
      } 
      else if ( iResult == 56 ) {
	++iRetryCount;
      
	if ( iRetryCount >= 20 ) {
	  HLTError(Form("Retry Failed: Error reading event from source %d: %s (%d), returning", 
			fCurrentReader->GetErrorConnectionNdx(), strerror(iResult), iResult));
	  break;
	} 
	else {
	  HLTError(Form("Retry: Error reading event from source %d: %s (%d), making another attempt (no %d out of 20)", 
			fCurrentReader->GetErrorConnectionNdx(), strerror(iResult), iResult, iRetryCount));
	  //break;
	  continue;
	}
      }
      else if ( iResult ) {
	HLTError(Form("General Error reading event from source %d: %s (%d), giving up", 
		      fCurrentReader->GetErrorConnectionNdx(), strerror(iResult), iResult));
	fConnected = kFALSE;
	break;
      } 
      else {
	HLTDebug("Successfully read out event from source");
	break;
      }

    } // while( 1 ) {
    
    // -- Check if event could be read
    if ( iResult ) {
      HLTInfo("Reading event from source failed");
      continue;
    }

    // -- Handle Blocks from current reader
    iResult = HandleBlocks();
    if ( iResult ) {
      HLTError(Form("Handling of blocks failed."));
    }

  } // while( (object = next()) ) {

  return iResult;  
}

/* ---------------------------------------------------------------------------------
 *                           Buffer Handling - public
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Int_t AliHLTHOMERManager::NavigateEventBufferBack() { 
  // see header file for class documentation

  // -- reached the end of the buffer
  if ( fNavigateBufferIdx == fBufferLowIdx )
    return -1;

  Int_t newIdx = fNavigateBufferIdx - 1;
  if ( newIdx == -1 )
    newIdx = BUFFERSIZE-1;

  fCurrentBufferIdx = fNavigateBufferIdx = newIdx;

  return newIdx;
}

//##################################################################################
Int_t AliHLTHOMERManager::NavigateEventBufferFwd() {
  // see header file for class documentation

  HLTInfo(Form("fNavigateBufferIdx: %d, fCurrentBufferIdx %d, fBufferTopIdx %d", fNavigateBufferIdx, fCurrentBufferIdx, fBufferTopIdx));

  // -- reached the top of the buffer
  if ( fNavigateBufferIdx == fBufferTopIdx )
    return -1;

  Int_t newIdx = fNavigateBufferIdx + 1;
  if ( newIdx == BUFFERSIZE )
    newIdx = 0;
  
  fCurrentBufferIdx = fNavigateBufferIdx = newIdx;
  fNEventsAvailable -= 1;

  HLTInfo(Form("fNavigateBufferIdx: %d, fCurrentBufferIdx %d, fBufferTopIdx %d", fNavigateBufferIdx, fCurrentBufferIdx, fBufferTopIdx));


  return 0;
}

 ///////////////////////////////////////////////////////////////////////////////////

/*
 * ---------------------------------------------------------------------------------
 *                            Connection Handling - private
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
void AliHLTHOMERManager::CreateReadoutList( const char** sourceHostnames, UShort_t *sourcePorts, 
					    UInt_t &sourceCount, TString detector ){
  // see header file for class documentation

  AliHLTHOMERSourceDesc * source= NULL;

  // -- Read all sources and check if they should be read out
  TIter next( fSourceList );
  while ( ( source = dynamic_cast<AliHLTHOMERSourceDesc*>(next()) ) ) {

    ///Don't use sources from dev cluster
    if(source->GetPort() > 60000) continue;

    // -- If detector NO detector name given
    if ( ! detector.CompareTo("ALL") ) {
      // -- Continue if source is not selected
      // HACK Jochen
      //if ( ! source->IsSelected() )
      //	continue;
    }
    // -- DetectorName given
    else {
      // -- Continue if detector name doesn't match
      if ( detector.CompareTo(source->GetDetector()) )
	continue;
      else
	source->Select();
    }
    
    Bool_t exists = kFALSE;
    
    // -- Loop over existing entries and check if entry is already in readout list
    for ( UInt_t ii = 0; ii < sourceCount; ii++ ){
      if ( !strcmp( sourceHostnames[ii], source->GetHostname().Data() ) 
	   && sourcePorts[ii] == source->GetPort() ) {
	exists = kTRUE;
	break;
      }
    }

    // -- Add new entires to readout list
    if ( ! exists ) {
      sourcePorts[sourceCount] = source->GetPort();
      sourceHostnames[sourceCount] = source->GetHostname().Data();
      sourceCount++;
    }

  } // while ( ( source = (AliHLTHOMERSourceDesc*)next() ) ) {

  fStateHasChanged = kFALSE;

  return;
}

/*
 * ---------------------------------------------------------------------------------
 *                          Buffer Handling - private
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
void AliHLTHOMERManager::AddBlockListToBuffer() {
  // see header file for class documentation
  // -- Check if event is already in buffer
  ULong_t eventID = static_cast<ULong64_t>(fCurrentReader->GetEventID());  
  
  if ( fEventID[fBufferTopIdx] == eventID ) {
    HLTInfo(Form("Event 0x%016lX (%lu) already in buffer.", eventID, eventID));
    return;
  }

  // -- Check if event should be selected on basis of trigger string
  if( fTriggerString.CompareTo("ALL") ){
    if ( !CheckTriggerDecision() ) {
      HLTInfo("Event not triggered");
      return;
    } else {
      HLTInfo("Event triggered");
    }
  }
  else {
    HLTDebug("No trigger selection.");
  }

  // -- Set Top mark 
  ++fBufferTopIdx;
  if ( fBufferTopIdx == BUFFERSIZE )
    fBufferTopIdx = 0;

  // -- Change the low mark if necessary
  if ( fBufferLowIdx == -1 )
    fBufferLowIdx = 0;
  else if ( fBufferTopIdx == fBufferLowIdx ) {
    ++fBufferLowIdx;
    if ( fBufferLowIdx == BUFFERSIZE )
      fBufferLowIdx = 0;
  }


  // -- Fill EventID
  fEventID[fBufferTopIdx] = eventID;

  // -- Clear Buffer slot
  (reinterpret_cast<TList*>((*fEventBuffer)[fBufferTopIdx]))->Clear();
  if(fBlockList->IsOwner()) HLTWarning("block list is owner!!");
  HLTInfo(Form("fBlockList size %d", fBlockList->GetSize()));
  //fBlockList->Clear();
  fBlockList = new TList();
  HLTInfo(Form("fBlockList size %d", fBlockList->GetSize()));

  GetFirstBlk();

  // -- Fill block list
  do {

    // -- Create new block
    AliHLTHOMERBlockDesc * block = new AliHLTHOMERBlockDesc();
    block->SetBlock( GetBlk(), GetBlkSize(), GetBlkOrigin(),
		     GetBlkType(), GetBlkSpecification() );
    
    // -- Check sources list if block is requested
    if ( CheckIfRequested( block ) ) {
      (reinterpret_cast<TList*>((*fEventBuffer)[fBufferTopIdx]))->Add( block );
      fBlockList->Add(block);
    }
    else {
      // XXX HACK Jochen
      (reinterpret_cast<TList*>((*fEventBuffer)[fBufferTopIdx]))->Add( block );
      fBlockList->Add(block);
      // delete block;
      // block = NULL;
    }
 
  } while( GetNextBlk() );

  //We have one more event available
  fNEventsAvailable++;
  HLTInfo(Form("fNEventsAvailable %d", fNEventsAvailable));
  return;
}

//##################################################################################
void AliHLTHOMERManager::AddToAsyncBlockList() {
  // see header file for class documentation

  HLTInfo("Adding blocks to the asynchroneous block list");

  GetFirstBlk();

  do {
    
    AliHLTHOMERBlockDesc * block = new AliHLTHOMERBlockDesc();
    block->SetBlock( GetBlk(), GetBlkSize(), GetBlkOrigin(),
		     GetBlkType(), GetBlkSpecification() );
    

    fAsyncBlockList->Add( block );
 
  } while( GetNextBlk() );

  return;
}
//__________________________________________________________________________________
void AliHLTHOMERManager::AddToBlockList() {
  // see header file for class documentation
  HLTInfo("Adding blocks to the synchroneous block list");

  ULong_t eventID = static_cast<ULong64_t>(fCurrentReader->GetEventID());  
  
  if ( fEventId == eventID ) {
    HLTInfo(Form("Event 0x%016lX (%lu) already in buffer.", eventID, eventID));
    return;
  }

  fEventId = eventID;

  GetFirstBlk();
  do {

    AliHLTHOMERBlockDesc * block = new AliHLTHOMERBlockDesc();
    block->SetBlock( GetBlk(), GetBlkSize(), GetBlkOrigin(),
		     GetBlkType(), GetBlkSpecification() );
    fBlockList->Add( block );
  
  } while( GetNextBlk() );  
}

//__________________________________________________________________________________
TList* AliHLTHOMERManager::GetBlockListEventBuffer() {
  // see header file for class documentation

  if(fBlockList)
    return fBlockList;
  else 
    return NULL;


}


//__________________________________________________________________________________
Int_t AliHLTHOMERManager::HandleBlocks() {
  // see header file for class documentation
  
  Int_t iResult = 0;

  // -- Get blockCnt and eventID
  fNBlks = static_cast<ULong_t>(fCurrentReader->GetBlockCnt());
  ULong_t eventID = static_cast<ULong64_t>(fCurrentReader->GetEventID());  
  fCurrentBlk = 0;

  // -- Check if blocks present
  if ( fNBlks ==  0 ) {
    HLTWarning(Form("Event 0x%016lX (%lu) with no blocks", eventID, eventID));
    return -1;
  }

  HLTInfo(Form("Event 0x%016lX (%lu) with %lu blocks", eventID, eventID, fNBlks));
    
  if ( IsSyncBlocks() ) {
    //AddBlockListToBuffer();
    fBlockList->Clear();
    AddToBlockList();
  } else {
    AddToAsyncBlockList();
  }

  return iResult;
}

//__________________________________________________________________________________
Bool_t AliHLTHOMERManager::IsSyncBlocks() {
  // see header file for class documentation
  
  Bool_t bResult = kFALSE;

  GetFirstBlk();
  
  do {          
    
    if ( !GetBlkType().CompareTo("ALIESDV0")) {
      bResult = kTRUE;
      break;
    }

    if ( !GetBlkType().CompareTo("GLOBTRIG")) {
      bResult = kTRUE;
      break;
    }
    
    if ( !GetBlkType().CompareTo("ROOTTOBJ") ) {
      AliHLTHOMERBlockDesc blockDesc;
      
      blockDesc.SetBlock( GetBlk(), GetBlkSize(), GetBlkOrigin(),
			  GetBlkType(), GetBlkSpecification() );
      if ( !blockDesc.GetClassName().CompareTo("AliHLTGlobalTriggerDecision") ) {

	bResult = kTRUE;
	break;
      }
    }

  } while( GetNextBlk() );


  return bResult;
}

//##################################################################################
void* AliHLTHOMERManager::GetBlk( Int_t ndx ) {
  // see header file for class documentation
  // Get pointer to current block in current event
   
  if ( !fCurrentReader || !IsConnected() ) {
    HLTError(Form("Not connected yet."));
    return NULL;
  }
  if ( ndx < static_cast<Int_t>(fNBlks) )
    return  const_cast<void*> (fCurrentReader->GetBlockData(ndx));
  else
    return NULL;
}

//##################################################################################
ULong_t AliHLTHOMERManager::GetBlkSize( Int_t ndx ) {
  // see header file for class documentation
   
  if ( !fCurrentReader || !IsConnected() ) {
    HLTError(Form("Not connected yet."));
    return 0;
  }
  
  if ( ndx < static_cast<Int_t>(fNBlks) )
    return static_cast<ULong_t> (fCurrentReader->GetBlockDataLength(ndx));
  else
    return 0;
}

//##################################################################################
TString AliHLTHOMERManager::GetBlkOrigin( Int_t ndx ) {
  // see header file for class documentation

  TString origin = "";

  // -- Check for Connection
  if ( !fCurrentReader || ! IsConnected() ) {
    HLTError(Form("Not connected yet."));
    return origin;
  }

  // -- Check block index
  if ( ndx >= static_cast<Int_t>(fNBlks) ) {
    HLTError(Form("Block index %d out of range.", ndx ));
    return origin;
  }

  // -- Get origin
  union{
    UInt_t data;
    Char_t array[4];
  } reverseOrigin;

  reverseOrigin.data = static_cast<UInt_t>(fCurrentReader->GetBlockDataOrigin(ndx));

  // -- Reverse the order
  for (Int_t ii = 3; ii >= 0; ii-- )
    if ( reverseOrigin.array[ii] != ' ')
      origin.Append( reverseOrigin.array[ii] );

  origin.Remove( TString::kTrailing, ' ' );

  return origin;
}

//##################################################################################
TString AliHLTHOMERManager::GetBlkType( Int_t ndx ) {
  // see header file for class documentation

  TString type = "";

  // -- Check for Connection
  if ( !fCurrentReader || ! IsConnected() ) {
    HLTError(Form("Not connected yet."));
    return type;
  }

  // -- Check block index
  if ( ndx >= static_cast<Int_t>(fNBlks) ) {
    HLTError(Form("Block index %d out of range.", ndx ));
    return type;
  }

  // -- Get type
  union{
    ULong64_t data;
    Char_t array[8];
  } reverseType;

  reverseType.data = static_cast<ULong64_t> (fCurrentReader->GetBlockDataType(ndx));

  // -- Reverse the order
  for (Int_t ii = 7; ii >= 0; ii-- )
    if ( reverseType.array[ii] != ' ')
      type.Append( reverseType.array[ii] );
  
  type.Remove( TString::kTrailing, ' ' );

  return type;
}

//##################################################################################
ULong_t AliHLTHOMERManager::GetBlkSpecification( Int_t ndx ) {
  // see header file for class documentation

  // -- Check for Connection
  if ( !fCurrentReader || ! IsConnected() ) {
    HLTError(Form("Not connected yet."));
    return 0;
  }

  // -- Check block index
  if ( ndx >= static_cast<Int_t>(fNBlks) ) {
    HLTError(Form("Block index %d out of range.", ndx ));
    return 0;
  }

  return static_cast<ULong_t>(fCurrentReader->GetBlockDataSpec(ndx));
}

//##################################################################################
Bool_t AliHLTHOMERManager::CheckIfRequested( AliHLTHOMERBlockDesc * block ) {
  // see header file for class documentation

  Bool_t requested = kFALSE;

  AliHLTHOMERSourceDesc * source= NULL;

  // -- Read all sources and check if they should be read out
  TIter next( fSourceList );
  while ( ( source = dynamic_cast<AliHLTHOMERSourceDesc*>(next()) ) ) {
    
    // -- Check if source is selected
    if ( ! source->IsSelected() )
      continue;
    
    // -- Check if detector matches
    if ( source->GetSourceName().CompareTo( block->GetBlockName() ) )
      continue;

    requested = kTRUE;
    break;

  } // while ( ( source = dynamic_cast<AliHLTHOMERSourceDesc*>(next()) ) ) {
  
#if EVE_DEBUG
  if ( requested ) {
    HLTInfo(Form("Block requested : %s", block->GetBlockName().Data())); 
  }
  else {
    HLTInfo(Form("Block NOT requested : %s", block->GetBlockName().Data())); 
  }
#endif

  return requested;
}

/* ---------------------------------------------------------------------------------
 *                          Trigger Handling - private
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Bool_t AliHLTHOMERManager::CheckTriggerDecision() {
  // see header file for class documentation

  Bool_t triggered = kFALSE;

  if ( !fCurrentReader || !IsConnected() ) {
    HLTError(Form("Not connected yet."));
    return kFALSE;
  }

  AliHLTHOMERBlockDesc blockDesc;

  GetFirstBlk();
  
  // -- Fill block list
  Bool_t foundTriggerBlock = kFALSE;
  
  do {
    if ( (GetBlkType().CompareTo("ROOTTOBJ") == 0) ) {
      blockDesc.SetBlock( GetBlk(), GetBlkSize(), GetBlkOrigin(),
			  GetBlkType(), GetBlkSpecification() );

      if ( ! blockDesc.GetClassName().CompareTo("AliHLTGlobalTriggerDecision") ) {

	foundTriggerBlock = kTRUE;
	break;
      }
      
    }
  } while( GetNextBlk() );
  
  if ( !foundTriggerBlock ) {
    HLTError(Form("No trigger decision object found"));
    return kFALSE;
  }

  // -- Get the global decision object
  AliHLTGlobalTriggerDecision* globalDecision = 
    static_cast<AliHLTGlobalTriggerDecision*>(blockDesc.GetTObject());

  if ( fTriggerString.CompareTo("HLTGlobalTrigger") == 0 ) {
    triggered = globalDecision->EventTriggered();
  } 
  else {
    
    for (Int_t idx = 0; idx < globalDecision->NumberOfInputObjects(); idx++) {
       
      const AliHLTTriggerDecision* triggerDecision = 
	reinterpret_cast<const AliHLTTriggerDecision*>(globalDecision->InputObject(idx));
    
      if ( !(fTriggerString.CompareTo(triggerDecision->Description())) ) {
	triggered = triggerDecision->EventTriggered();
	break;
      }
    } // for (Int_t idx = 0; idx < globalDecision->NumberOfInputObjects(); idx++) {
  }



  if ( triggered ) {
    fRetryNextEvent = kFALSE;
    fNEventsNotTriggered = 0;
  }
  else {
    fRetryNextEvent = kTRUE;
    ++fNEventsNotTriggered;
  }

  return triggered;
}
