//-*- Mode: C++ -*-

// $Id: AliHLTHOMERManager.h $

#ifndef ALIHLTHOMERMANAGER_H
#define ALIHLTHOMERMANAGER_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice     
 */

/** @file   AliHLTHOMERManager.h
    @author Jochen Thaeder
    @date
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

#define BUFFERSIZE 10

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
  
  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** default constructor */
  AliHLTHOMERManager();

  /** destructor */
  virtual ~AliHLTHOMERManager();

  /** Initialize 
   *  @return 0 on success, <0 for failure
   */
  Int_t Initialize();

  /*
   * ---------------------------------------------------------------------------------
   *                            Source Handling - public
   * ---------------------------------------------------------------------------------
   */

  /** Create Sources List from HOMER-Proxy 
   *  @return 0 on success, <0 for failure, 1 for no active service
   */
  virtual Int_t CreateSourcesList();

  /** Set state of a source 
   *  @param source      Pointer to AliHLTHOMERSourceDesc object.
   *  @param state       New (selected/not selected) state.
   */
  void   SetSourceState( AliHLTHOMERSourceDesc* source, Bool_t state);

  /** Get pointer to source List */
  TList* GetSourceList() { return fSourceList; }

  /*
   * ---------------------------------------------------------------------------------
   *                            Connection Handling - public
   * ---------------------------------------------------------------------------------
   */

  /** Connect to HOMER sources, of a certain detector.
   *  which gets created when state has changed 
   *  @param detector    Detector to be connected to
   *  @return            0 on success, <0 for failure
   */
  Int_t ConnectHOMER( TString detector="ALL" );

  /** Disconnect from HOMER sources */
  void  DisconnectHOMER();

  /** Reconnect from HOMER sources 
   *  @param detector    Detector to be connected to
   *  @return            0 on success, <0 for failure
   */
  Int_t ReconnectHOMER( TString detector);

  /*
   * ---------------------------------------------------------------------------------
   *                            Event Handling - public
   * ---------------------------------------------------------------------------------
   */

  /** Loads the next Event, after being connected 
   *  @return 0 on success, <0 for failure
   */
  virtual Int_t NextEvent();

  /** Loads the next Cycle, after being connected 
   *  @return 0 on success, <0 for failure
   */
  virtual Int_t NextCycle() { return NextEvent(); }

  /** Get event ID */
  ULong_t GetEventID() { return fEventID[fCurrentBufferIdx]; }

  /* ---------------------------------------------------------------------------------
   *                           Buffer Handling - public
   * ---------------------------------------------------------------------------------
   */

  /** Get pointer to last requested BlockList 
   *  @return     ptr to buffer, NULL if buffer boundary reached                
   */
  TList* GetBlockList() { return GetBlockListEventBuffer(fCurrentBufferIdx); }

  /** Navigate backwards in event buffer 
   *  @return      index in buffer, -1 if boundary reached                
   */
  Int_t NavigateEventBufferBack();

  /** Navigate forwards in event buffer 
   *  @return      index in buffer, -1 if boundary reached                
   */
  Int_t NavigateEventBufferFwd();

  ///////////////////////////////////////////////////////////////////////////////////

protected:

  /** Dynamic loader manager for the HOMER library */
  AliHLTHOMERLibManager* fLibManager;             //! transient

  /** Indicates, if a sources have changes, 
   *  so that one has to reconnect. */
  Bool_t    fStateHasChanged;                     //  see above

  ///////////////////////////////////////////////////////////////////////////////////

private:

  /** copy constructor prohibited */
  AliHLTHOMERManager(const AliHLTHOMERManager&);

  /** assignment operator prohibited */
  AliHLTHOMERManager& operator=(const AliHLTHOMERManager&);

  /*
   * ---------------------------------------------------------------------------------
   *                            Connection Handling - private
   * ---------------------------------------------------------------------------------
   */

  /** Create a readout list for Hostname and ports 
   *  @param socurceHostnames   Array of selected hostnames
   *  @param socurcePorts       Array of selected ports
   *  @param socurceCount       Number of selected hostname:port
   *  @param detector           detector to be selected
   */
  void CreateReadoutList( const char** sourceHostnames, UShort_t* sourcePorts, 
			  UInt_t &sourceCount, TString detector );

  /** Checks if already connected to HOMER sources */
  Bool_t IsConnected() { return fConnected; }  
  
  /* ---------------------------------------------------------------------------------
   *                           Buffer Handling - private
   * ---------------------------------------------------------------------------------
   */

  /** Create and add Block List to Buffer */
  void AddBlockListToBuffer();

  /** Get pointer to block list in event buffer 
   *  @return     ptr to buffer, NULL if not present
   */
  TList* GetBlockListEventBuffer( Int_t idx );
    
  /*
   * ---------------------------------------------------------------------------------
   *                            Block Handling - private
   * ---------------------------------------------------------------------------------
   */

  /** Get Number of blocks in current event */
  ULong_t GetNBlks() { return fNBlks; }

  // ----------------------------------------------------

  /** Get pointer to block ndx in current event 
   *  @param ndx        Block index
   *  @return           returns pointer to blk, NULL if no block present
   */
  void* GetBlk( Int_t ndx );

  /** Get pointer to current block in current event */
  void* GetBlk() { return GetBlk(fCurrentBlk); }

  /** Get first block in current event */
  void* GetFirstBlk() { fCurrentBlk=0; return GetBlk(0); }

  /** Get next block in current event */
  void* GetNextBlk() { return GetBlk(++fCurrentBlk); }

  // ----------------------------------------------------

  /** Get size of block ndx 
   *  @param ndx        Block index
   *  @return           returns size blk, 0 otherwise
   */
  ULong_t GetBlkSize( Int_t ndx );

  /** Get size of current block */ 
  ULong_t GetBlkSize() { return GetBlkSize( fCurrentBlk ); }

  // ---------------------------------------------------- 

  /** Get origin of block ndx 
   *  @param ndx        Block index
   *  @return           origin of block
   */
  TString GetBlkOrigin( Int_t ndx );

  /** Get origin of current block */
  TString GetBlkOrigin(){ return GetBlkOrigin( fCurrentBlk ); }

  // ----------------------------------------------------

  /** Get type of block ndx 
   *  @param ndx        Block index
   *  @return           type of block
   */
  TString GetBlkType( Int_t ndx ); 

  /** Get type of current block */
  TString GetBlkType() { return GetBlkType( fCurrentBlk ); } 
  
  // ----------------------------------------------------
  
  /** Get specification of block ndx 
   *  @param ndx        Block index
   *  @return           specification of block
   */
  ULong_t GetBlkSpecification( Int_t ndx );

  /** Get specification of current block */
  ULong_t GetBlkSpecification() { return GetBlkSpecification( fCurrentBlk ); } 

  // ----------------------------------------------------

  /** Checks if current Block should was requested 
   *  @return           returns kTRUE, if block should was requested
   */
  Bool_t CheckIfRequested( AliHLTHOMERBlockDesc* block );

  /*
   * ---------------------------------------------------------------------------------
   *                            Members - private
   * ---------------------------------------------------------------------------------
   */

  /** Proxy Handler to get the list of sources */
  AliHLTHOMERProxyHandler *fProxyHandler;               //! transient 

  // == connection ==

  /** Pointer to HOMER reader */
  AliHLTHOMERReader       *fReader;                     //! transient 

  // == sources ==

  /** List to HOMER sources */
  TList                   *fSourceList;                 //! transient

  // == events ==

  /** Number of blockes in current event */
  ULong_t                  fNBlks;                      //  see above

  /** EventID of current event */
  ULong64_t                fEventID[BUFFERSIZE];                    //  see above

  /** Current block in current event */
  ULong_t                  fCurrentBlk;                 //  see above

  // == event buffer ==

  /** Event Buffer */
  TClonesArray            *fEventBuffer;                //  see above

  /** Buffer index to last received event */
  Int_t                    fBufferTopIdx;               //  see above

  /** Buffer index to last received event */
  Int_t                    fBufferLowIdx;               //  see above

  /** Buffer index to current event */
  Int_t                    fCurrentBufferIdx;           //  see above

  /** Navigate index through event buffer */
  Int_t                    fNavigateBufferIdx;          //  see above
  
  // == states ==
  
  /** Shows connection status */
  Bool_t        fConnected;                              //  see above


  ClassDef(AliHLTHOMERManager, 1); // Manage connections to HLT data-sources.
};

#endif
