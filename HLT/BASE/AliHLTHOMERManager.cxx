//-*- Mode: C++ -*-
// $Id: AliHLTHOMERManager.cxx  $
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
  fReader(NULL),
  fSourceList(NULL),
  fNBlks(0),
  fEventID(),
  fCurrentBlk(0),
  fEventBuffer(NULL),
  fBufferTopIdx(-1),
  fBufferLowIdx(-1),
  fCurrentBufferIdx(-1),
  fNavigateBufferIdx(-1),
  fConnected(kFALSE) {
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
    if ( fReader )
      fLibManager->DeleteReader(fReader);
    delete fLibManager;
    fLibManager = NULL;
    fReader = NULL;
  }

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

}

//##################################################################################
Int_t AliHLTHOMERManager::Initialize() {
  // see header file for class documentation

  Int_t iResult = 0;

  if ( !fProxyHandler )
    fProxyHandler = new AliHLTHOMERProxyHandler();
  
  if ( fProxyHandler ) {
    iResult = fProxyHandler->Initialize();
    if (iResult)
      HLTError(Form("Initialize of ProxyHandler failed."));
  }
  else {
    iResult = -1;
    HLTError(Form("Creating of ProxyHandler failed."));
  }
 
  // -- Initialize Event Buffer
  if ( !fEventBuffer ) {
    fEventBuffer = new TClonesArray( "TList", BUFFERSIZE );
  }

  for ( Int_t idx = 0; idx < BUFFERSIZE; ++idx ) {
    new ((*fEventBuffer)[idx]) TList();
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

  // -- Check if LibManager is present
  if ( ! fLibManager ) {
    HLTError(Form("No LibManager present."));
    return -1;
  }
  
  // -- Check if already connected and state has not changed
  if ( fStateHasChanged == kFALSE && IsConnected() ) {
    HLTInfo(Form("No need for reconnection."));
    return iResult;
  }
  
  // -- If already connected, disconnect before connect
  if ( IsConnected() )
    DisconnectHOMER();
  
  // -- Create the Readoutlist
  UShort_t* sourcePorts = new UShort_t [fSourceList->GetEntries()];
  const char ** sourceHostnames = new const char* [fSourceList->GetEntries()];
  UInt_t sourceCount = 0;
  
  CreateReadoutList( sourceHostnames, sourcePorts, sourceCount, detector );
  if ( sourceCount == 0 ) {
    HLTError(Form("No sources selected, aborting."));
    return -2;
  }

  // *** Connect to data sources
  if ( !fReader )
    fReader = fLibManager->OpenReader( sourceCount, sourceHostnames, sourcePorts );
  
  iResult = fReader->GetConnectionStatus();
  if ( iResult ) {
    // -- Connection failed

    UInt_t ndx = fReader->GetErrorConnectionNdx();

    if ( ndx < sourceCount ) {
      HLTError(Form("Error establishing connection to TCP source %s:%hu: %s (%d)",
		    sourceHostnames[ndx], sourcePorts[ndx], strerror(iResult), iResult));
    }
    else {
      HLTError(Form("Error establishing connection to unknown source with index %d: %s (%d)",
		    ndx, strerror(iResult), iResult));
    }

    if ( fReader )
      fLibManager->DeleteReader( fReader );
    fReader = NULL;
  }
  else {
    // -- Connection ok - set reader
    fConnected = kTRUE;

    HLTInfo(Form("Connection established."));
  }

  delete[] sourceHostnames;
  delete[] sourcePorts;

  return iResult;
}

//##################################################################################
void AliHLTHOMERManager::DisconnectHOMER(){
  // see header file for class documentation

  if ( ! IsConnected() )
    return;

  if ( fReader )
    fLibManager->DeleteReader( fReader );
  fReader = NULL;

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
  
  if ( !IsConnected() || fStateHasChanged )
    ConnectHOMER();
  
  if ( !fReader || !IsConnected() ) {
    HLTWarning(Form( "Not connected yet." ));
    return -1;
  }

  //  fReader->SetEventRequestAdvanceTime( 20000000 /*timeout in us*/ );

  // -- Read next event data and error handling for HOMER (error codes and empty blocks)
  while( 1 ) {
    
    iResult = fReader->ReadNextEvent( 40000000 /*timeout in us*/);

    if ( iResult == 111 || iResult == 32 || iResult == 6 ) {
      HLTError(Form("No Connection to source %d: %s (%d)", 
		    fReader->GetErrorConnectionNdx(), strerror(iResult), iResult));
      return -iResult;
    }
    else if ( iResult == 110 ) {
      HLTError(Form("Timeout occured, reading event from source %d: %s (%d)", 
		    fReader->GetErrorConnectionNdx(), strerror(iResult), iResult));
      return -iResult;
    }
    else if ( iResult == 56) {
      ++iRetryCount;

      if ( iRetryCount >= 20 ) {
	HLTError(Form("Retry Failed: Error reading event from source %d: %s (%d)", 
		      fReader->GetErrorConnectionNdx(), strerror(iResult), iResult));
	return -iResult;
      }
      else {
	HLTError(Form("Retry: Error reading event from source %d: %s (%d)", 
		      fReader->GetErrorConnectionNdx(), strerror(iResult), iResult));
	continue;
      }
    }
    else if ( iResult ) {
      HLTError(Form("General Error reading event from source %d: %s (%d)", 
		    fReader->GetErrorConnectionNdx(), strerror(iResult), iResult));
      fConnected = kFALSE;
      return -iResult;
    }
    else {
      break;
    }
  } // while( 1 ) {

  // -- Get blockCnt and eventID
  fNBlks = static_cast<ULong_t>(fReader->GetBlockCnt());
  ULong_t eventID = static_cast<ULong64_t>(fReader->GetEventID());  
  fCurrentBlk = 0;

  HLTInfo(Form("Event 0x%016LX (%Lu) with %lu blocks", eventID,eventID, fNBlks));

#if EVE_DEBUG
  // Loop for Debug only
  for ( ULong_t ii = 0; ii < fNBlks; ii++ ) {
    Char_t tmp1[9], tmp2[5];
    memset( tmp1, 0, 9 );
    memset( tmp2, 0, 5 );
    void *tmp11 = tmp1;
    ULong64_t* tmp12 = static_cast<ULong64_t*>(tmp11);
    *tmp12 = fReader->GetBlockDataType(ii);
    void *tmp21 = tmp2;
    ULong_t* tmp22 = static_cast<ULong_t*>(tmp21);
    *tmp22 = fReader->GetBlockDataOrigin(ii);
    HLTInfo(Form( "Block %lu length: %lu - type: %s - origin: %s - spec 0x%08X",
		  ii, fReader->GetBlockDataLength(ii), tmp1, tmp2, fReader->GetBlockDataSpec(ii) ));
  } // end for ( ULong_t ii = 0; ii < fNBlks; ii++ ) {
#endif

  // -- Create BlockList
  if ( fNBlks > 0 ) {
    HLTInfo(Form("Add Block List to buffer"));
    AddBlockListToBuffer();
  }
  else {
    HLTWarning(Form("Event 0x%016LX (%Lu) with %lu blocks", eventID, eventID, fNBlks));
  }
    
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

  // -- reached the top of the buffer
  if ( fNavigateBufferIdx == fBufferTopIdx )
    return -1;

  Int_t newIdx = fNavigateBufferIdx + 1;
  if ( newIdx == BUFFERSIZE )
    newIdx = 0;
  
  fCurrentBufferIdx = fNavigateBufferIdx = newIdx;

  return newIdx;
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

    // -- If detector NO detector name given
    if ( ! detector.CompareTo("ALL") ) {
      // -- Continue if source is not selected
      if ( ! source->IsSelected() )
	continue;
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

  fNavigateBufferIdx = fCurrentBufferIdx = fBufferTopIdx;    

  // -- Fill EventID
  fEventID[fBufferTopIdx] = static_cast<ULong64_t>(fReader->GetEventID());

  // -- Clear Buffer slot
  (reinterpret_cast<TList*>((*fEventBuffer)[fBufferTopIdx]))->Clear();

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
    }
    else {
      // XXX HACK Jochen
      (reinterpret_cast<TList*>((*fEventBuffer)[fBufferTopIdx]))->Add( block );
      //      delete block;
      //      block = NULL;
    }
 
  } while( GetNextBlk() );

  return;
}

//##################################################################################
TList* AliHLTHOMERManager::GetBlockListEventBuffer( Int_t idx ) {
  // see header file for class documentation
  
  if ( idx == -1 )
    return NULL;

  return reinterpret_cast<TList*>((*fEventBuffer)[idx]);
}

/*
 * ---------------------------------------------------------------------------------
 *                          Block Handling - private
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
void* AliHLTHOMERManager::GetBlk( Int_t ndx ) {
  // see header file for class documentation
  // Get pointer to current block in current event
   
  if ( !fReader || !IsConnected() ) {
    HLTError(Form("Not connected yet."));
    return NULL;
  }
  if ( ndx < static_cast<Int_t>(fNBlks) )
    return  const_cast<void*> (fReader->GetBlockData(ndx));
  else
    return NULL;
}

//##################################################################################
ULong_t AliHLTHOMERManager::GetBlkSize( Int_t ndx ) {
  // see header file for class documentation
   
  if ( !fReader || !IsConnected() ) {
    HLTError(Form("Not connected yet."));
    return 0;
  }
  
  if ( ndx < static_cast<Int_t>(fNBlks) )
    return static_cast<ULong_t> (fReader->GetBlockDataLength(ndx));
  else
    return 0;
}

//##################################################################################
TString AliHLTHOMERManager::GetBlkOrigin( Int_t ndx ) {
  // see header file for class documentation

  TString origin = "";

  // -- Check for Connection
  if ( !fReader || ! IsConnected() ) {
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

  reverseOrigin.data = static_cast<UInt_t>(fReader->GetBlockDataOrigin(ndx));

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
  if ( !fReader || ! IsConnected() ) {
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

  reverseType.data = static_cast<ULong64_t> (fReader->GetBlockDataType(ndx));

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
  if ( !fReader || ! IsConnected() ) {
    HLTError(Form("Not connected yet."));
    return 0;
  }

  // -- Check block index
  if ( ndx >= static_cast<Int_t>(fNBlks) ) {
    HLTError(Form("Block index %d out of range.", ndx ));
    return 0;
  }

  return static_cast<ULong_t>(fReader->GetBlockDataSpec(ndx));
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

