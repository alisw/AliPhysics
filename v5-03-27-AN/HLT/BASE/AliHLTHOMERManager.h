//-*- Mode: C++ -*-

// $Id$

#ifndef ALIHLTHOMERMANAGER_H
#define ALIHLTHOMERMANAGER_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice     
 */

/** @file   AliHLTHOMERManager.h
    @author Jochen Thaeder
    @author Svein Lindal <slindal@fys.uio.no>
    @date   October 2010
    @brief  Manager for HOMER in aliroot
*/


#include "TClonesArray.h"
#include "TString.h"
#include "TList.h"

#include "AliHLTHOMERSourceDesc.h"
#include "AliHLTHOMERBlockDesc.h"
#include "AliHLTHOMERReader.h"
#include "AliHLTHOMERProxyHandler.h"

#include "AliHLTLoggingVariadicFree.h"

#define BUFFERSIZE 15

class AliHLTHOMERLibManager;

/**
 * @class AliHLTHOMERManager
 * This Class should handle the communication
 * from the HLT to aliroot. The HLT sends data via 
 * the HOMER interface on several TCP ports of nodes 
 * in the CERN GPN and DCS network.
 * All this communication is hidden from the user.
 * 
 * Right now, a xml file ( SCC1 ) is used to get the
 * configuration, this will/ has to change to a proxy
 * running on dedicated nodes.
 *
 * @ingroup alihlt_homer
 */

class AliHLTHOMERManager : public AliHLTLogging 
{
public:
  
  /** default constructor */
  AliHLTHOMERManager();

  /** destructor */
  virtual ~AliHLTHOMERManager();

  /** Initialize */
  Int_t Initialize();

  /** Create Sources List from HOMER-Proxy */
  virtual Int_t CreateSourcesList();

  /** Set state of a source */
  void   SetSourceState( AliHLTHOMERSourceDesc* source, Bool_t state);

  /** Get pointer to source List */
  TList* GetSourceList() { return fSourceList; }

  /** Connect to HOMER sources, of a certain detector. */
  Int_t ConnectHOMER( TString detector="ALL" );

  /** Disconnect from HOMER sources */
  void  DisconnectHOMER();

  /** Reconnect from HOMER sources */
  Int_t ReconnectHOMER( TString detector);


  /** Loads the next Event, after being connected */
  virtual Int_t NextEvent();

  /** Loads the next Cycle, after being connected */
  virtual Int_t NextCycle() { return NextEvent(); }

  /** Get event ID */
  ULong_t GetEventID() { return fEventId; }

  Int_t GetNAvailableEvents() { return fNEventsAvailable;}
  
  /** Get pointer to last requested BlockList */
  TList* GetBlockList() { return fBlockList; }
  TList* GetAsyncBlockList() { return fAsyncBlockList; }

  /** Navigate backwards in event buffer */
  Int_t  NavigateEventBufferBack();

  /** Navigate forwards in event buffer */
  Int_t  NavigateEventBufferFwd();

  /** Set and get the string used to select triggers */
  void SetTriggerString ( TString triggerString ) { fTriggerString = triggerString; }

  /** Get TriggerString */
  TString GetTriggerString () { return fTriggerString; }

  void SetBlockOwner(Bool_t owner) { fBlockList->SetOwner(owner); }
  Bool_t GetBlockOwner() const { return fBlockList->IsOwner(); }

protected:

  /** Dynamic loader manager for the HOMER library */
  AliHLTHOMERLibManager* fLibManager;             //! transient

  /** Indicates, if a sources have changes,  so that one has to reconnect. */
  Bool_t    fStateHasChanged;                     //  see above

  Bool_t Connected() const { return fConnected; }

private:

  /** copy constructor prohibited */
  AliHLTHOMERManager(const AliHLTHOMERManager&);

  /** assignment operator prohibited */
  AliHLTHOMERManager& operator=(const AliHLTHOMERManager&);

  //==============Connection to homer ==========================

  /** Create a readout list for Hostname and ports */
  void CreateReadoutList( const char** sourceHostnames, UShort_t* sourcePorts, 
			  UInt_t &sourceCount, TString detector );

  /** Checks if already connected to HOMER sources */
  Bool_t IsConnected() { return fConnected; }  

  /** Create and add Block List to Buffer */
  void AddBlockListToBuffer();

  /** Add bocks to asynchronous BlockList */
  void AddToAsyncBlockList();
  void AddToBlockList();


  //============ Block Handling ====================

  /** Get pointer to block list in event buffer */
  TList* GetBlockListEventBuffer( );
    
  /** Get Number of blocks in current event */
  ULong_t GetNBlks() { return fNBlks; }

  /** Handle Blocks and fill them in event buffer or asyncronous BlockList */
  Int_t HandleBlocks();

  /** Check is block are from syncronous source */
  Bool_t IsSyncBlocks();

  /** Get pointer to block ndx in current event */
  void* GetBlk( Int_t ndx );

  /** Get pointer to current block in current event */
  void* GetBlk() { return GetBlk(fCurrentBlk); }

  /** Get first block in current event */
  void* GetFirstBlk() { fCurrentBlk=0; return GetBlk(0); }

  /** Get next block in current event */
  void* GetNextBlk() { return GetBlk(++fCurrentBlk); }

  /** Get size of block ndx */
  ULong_t GetBlkSize( Int_t ndx );

  /** Get size of current block */ 
  ULong_t GetBlkSize() { return GetBlkSize( fCurrentBlk ); }

  /** Get origin of block ndx  */
  TString GetBlkOrigin( Int_t ndx );

  /** Get origin of current block */
  TString GetBlkOrigin(){ return GetBlkOrigin( fCurrentBlk ); }

  /** Get type of block ndx */
  TString GetBlkType( Int_t ndx ); 

  /** Get type of current block */
  TString GetBlkType() { return GetBlkType( fCurrentBlk ); } 
  
  //Get specification of block at ndx in bufferindex
  ULong_t GetBlkSpecification( Int_t ndx );

  /** Get specification of current block */
  ULong_t GetBlkSpecification() { return GetBlkSpecification( fCurrentBlk ); } 

  //Check if requested in eve
  Bool_t CheckIfRequested( AliHLTHOMERBlockDesc* block );
    
  //Check trigger decision
  Bool_t CheckTriggerDecision();
  
  AliHLTHOMERProxyHandler * fProxyHandler;  /** Proxy Handler to get the list of sources */  //! transient 
  AliHLTHOMERReader* fCurrentReader;   /** Pointer to current HOMER reader */ //! transient 
  TList* fReaderList;                 /** List to pointer of HOMER readers */

  // == sources ==
  TList* fSourceList;                /** List to HOMER sources */
  ULong_t fNBlks;                    /** Number of blockes in current event */
  ULong64_t fEventID[BUFFERSIZE];    /** EventID of current event */
  ULong64_t fEventId;
  ULong_t fCurrentBlk;               /** Current block in current event */
  TList* fAsyncBlockList;            /** List containing asychronous blocks */
  TList* fBlockList;            /** List containing asychronous blocks */

  // == event buffer ==
  TClonesArray     * fEventBuffer;  /** Event Buffer */
  Int_t  fBufferTopIdx;             /** Buffer index to last received event */
  Int_t  fBufferLowIdx;             /** Buffer index to last received event */
  Int_t  fCurrentBufferIdx;         /** Buffer index to current event */
  Int_t  fNavigateBufferIdx;        //  Navigate index through event buffer */
  Int_t  fNEventsAvailable;         //Number of available events
  
  Bool_t fConnected;                /** Shows connection status */
  TString fTriggerString;           /** String indicating which trigger should be used to select events */
  Int_t  fNEventsNotTriggered;      /** Number Events not triggered, before next triggered event is found */
  
  Bool_t fRetryNextEvent;           /** Retry reading next event */

  Bool_t fIsBlockOwner;

  ClassDef(AliHLTHOMERManager, 1); // Manage connections to HLT data-sources.
};

#endif
