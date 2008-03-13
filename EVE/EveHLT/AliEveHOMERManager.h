// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

//-*- Mode: C++ -*-
#ifndef ALIEVEHOMERMANGER_H
#define ALIEVEHOMERMANGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliEveHOMERManager.h
    @author Jochen Thaeder
    @date
    @brief  Manager for HOMER in offline
*/

#include <TEveElement.h>

#include "TString.h"
#include "TDOMParser.h"
#include "TXMLNode.h"
#include "TList.h"

#include "AliHLTHOMERSourceDesc.h"
#include "AliHLTHOMERBlockDesc.h"
#include "AliHLTHOMERReader.h"

#include "AliTPCPreprocessorOnline.h"


class AliHLTHOMERLibManager;

/**
 * @class AliEveHOMERManager
 *
 * This class is the main class of the AliEveHOMERManager
 * ... more to come
 *
 * @ingroup alihlt_homer
 */

class AliEveHOMERManager : public TEveElementList
{
private:
  AliEveHOMERManager(const AliEveHOMERManager&);            // Not implemented.
  AliEveHOMERManager& operator=(const AliEveHOMERManager&); // Not implemented.

public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** constructor
   * @param argc    Number of command line arguments.
   * @param argv    Array of command line arguments.
   */
  AliEveHOMERManager(TString xmlFile="" );
  /** destructor */
  virtual ~AliEveHOMERManager();

  /*
   * ---------------------------------------------------------------------------------
   *                            Source Handling - public
   * ---------------------------------------------------------------------------------
   */

  /** Create Sources List from HOMER-Proxy */
  Int_t CreateHOMERSourcesList();

  /** Set state of a source
   * @param source      Pointer to AliHLTHOMERSourceDesc object.
   * @param state       New (selected/not selected) state.
   */
  void SetSourceState( AliHLTHOMERSourceDesc* source, Bool_t state);

  /** Get pointer to source List
   * @return            returns pointer to TList of sources
   */
  TList* GetSourceList() { return fSourceList; }

  /*
   * ---------------------------------------------------------------------------------
   *                            Connection Handling - public
   * ---------------------------------------------------------------------------------
   */

  /** Connect to HOMER sources, out of Readout List, which gets created when state has changed
   * @return            0 on sucess, "HOMER" errors on error
   */
  Int_t ConnectHOMER();

  /** Disconnect from HOMER sources */
  void DisconnectHOMER();

  /** Reconnect from HOMER sources
   * @return            0 on sucess, "ConnectHOMER()" errors on error
   */
  Int_t ReconnectHOMER();

  /*
   * ---------------------------------------------------------------------------------
   *                            Event Handling - public
   * ---------------------------------------------------------------------------------
   */

  /** Loads the next Event, after being connected
   * @return            0 on sucess, "HOMER" errors on error
   */
  Int_t NextEvent();

  /** Get event ID
   * @return            Returns eventID
   */
  ULong_t GetEventID() { return fEventID; }

  /** Get pointer to block List
   * @return            returns pointer to TList of blocks
   */
  TList* GetBlockList() { return fBlockList; }


  ///////////////////////////////////////////////////////////////////////////////////

  void SelectRawTPC();
  void SelectClusterTPC();
  void SelectESDTPC();

  void TestSelect();
  void TestSelectClass( TString objectName );

  void DumpTPCCalib(TString objectName, Bool_t dumpToFile);

  ///////////////////////////////////////////////////////////////////////////////////

protected:

  /** Dynamic loader manager for the HOMER library */
  AliHLTHOMERLibManager* fLibManager;             //! transient


private:

  /*
   * ---------------------------------------------------------------------------------
   *                            Source Handling - private
   * ---------------------------------------------------------------------------------
   */

  /** Get Information out of a TDS process in XML file
   * @param xmlNode   Pointer to childs of TDS node
   * @return          0 on sucess, > 0 on error
   */
  Int_t GetTDSAttributes( TXMLNode * xmlNode );

  /** Resolve Information of hostname and port for source which has to be used by HOMER
   * ( due to port mapping inside the HLT )
   * @param xmlHostname  Hostname out of the XML
   * @param xmlPort      Port out of the XML
   * @param hostname     Return of the hostname
   * @param port         Return of the port
   * @return             0 on sucess, 1 if hostname couldn't be resolved, 2 if port couldn't be resolved,
   */
  Int_t ResolveHostPortInformation( TString xmlHostname, TString xmlPort, TString &hostname, Int_t &port );

  /** Resolve information of source
   * @param xmlParent   ParentString out of the XML
   * @param source      Return the filled AliHLTHOMERSourceDesc object
   * @return            0 on sucess, 1 on error
   */
  Int_t ResolveSourceInformation( TString xmlParent, AliHLTHOMERSourceDesc * source );

  /*
   * ---------------------------------------------------------------------------------
   *                            Connection Handling - private
   * ---------------------------------------------------------------------------------
   */

  /** Create a readout list for Hostname and ports
   * @param socurceHostnames   Array of selected hostnames
   * @param socurcePorts       Array of selected ports
   * @param socurceCount       Number of selected hostname:port
   */
  void CreateReadoutList( const char** socurceHostnames, UShort_t* sourcePorts, UInt_t &sourceCount);

  /** Checks if already connected to HOMER sources
   * @return             kTRUE or kFALSE, depending on connection state
   */
  Bool_t IsConnected() { return fConnected; }

  /* ---------------------------------------------------------------------------------
   *                            Eve AliEveHOMERManager::foo(nt Handling - private
   * ---------------------------------------------------------------------------------
   */

  /** Create a TList of blocks, which have been readout */
  Int_t CreateBlockList();

  /*
   * ---------------------------------------------------------------------------------
   *                            Block Handling - private
   * ---------------------------------------------------------------------------------
   */

  /** Get Number of blocks in current event
   * @return           returns number of blocks in current event
   */
  ULong_t GetNBlks() { return fNBlks; }

  /** Get pointer to block ndx in current event
   * @param ndx        Block index
   * @return           returns pointer to blk, NULL if no block present
   */
  void* GetBlk( Int_t ndx );

  /** Get pointer to current block in current event
   * @param ndx        Block index
   * @return           returns pointer to blk, NULL if no block present
   */
  void* GetBlk() { return GetBlk( fCurrentBlk ); }

  /** Get first block in current event
   * @return           returns pointer to blk, NULL if no block present
   */
  void* GetFirstBlk() { return GetBlk( 0 ); }

  /** Get next block in current event
   * @return           returns pointer to blk, NULL if no block present
   */
  void* GetNextBlk() { return GetBlk( ++fCurrentBlk ); }

  /** Get size of block ndx
   * @param ndx        Block index
   * @return           returns pointer to blk, 0 if no block present
   */
  ULong_t GetBlkSize( Int_t ndx );

  /** Get size of current block
   * @param ndx        Block index
   * @return           returns pointer to blk, 0 if no block present
   */
  ULong_t GetBlkSize() { return GetBlkSize( fCurrentBlk ); }

  /** Get origin of block ndx
   * @param ndx        Block index
   * @return           origin of block
   */
  TString GetBlkOrigin( Int_t ndx );

  /** Get origin of current block
   * @param ndx        Block index
   * @return           origin of block
   */
  TString GetBlkOrigin(){ return GetBlkOrigin( fCurrentBlk ); }

  /** Get type of block ndx
   * @param ndx        Block index
   * @return           type of block
   */
  TString GetBlkType( Int_t ndx );

  /** Get type of current block
   * @param ndx        Block index
   * @return           type of block
   */
  TString GetBlkType() { return GetBlkType( fCurrentBlk ); }

  /** Get specification of block ndx
   * @param ndx        Block index
   * @return           specification of block
   */
  ULong_t GetBlkSpecification( Int_t ndx );

  /** Get specification of current block
   * @param ndx        Block index
   * @return           specification of block
   */
  ULong_t GetBlkSpecification(){ return GetBlkSpecification( fCurrentBlk ); }

  /** Checks if current Block should was requested
   * @return           returns kTRUE, if block should was requested
   */
  Bool_t CheckIfRequested( AliHLTHOMERBlockDesc* block );

  /*
   * ---------------------------------------------------------------------------------
   *                            Members - private
   * ---------------------------------------------------------------------------------
   */

  // == XML parser ==

  /** 	 */
  TString fXMLFile;                               // XML input file

  /**  */
  TDOMParser* fXMLParser;                         //! XML parser into DOM model

  /**  */
  TXMLNode * fRootNode;                           //! Root node of parsed config file

  // == sources ==

  /**  */
  TList * fSourceList;                            //! List to HOMER sources

  // == connection ==


  /**  */
  AliHLTHOMERReader* fReader;                     //! Pointer to HOMER reader

  // == blocks ==

  /**  */
  TList * fBlockList;                             //! List to HOMER blocks

  // == events ==

  /**  */
  ULong_t fNBlks;                                 // Number of blockes in current event

  /**  */
  ULong64_t fEventID;                             // EventID of current event

  /**  */
  ULong_t fCurrentBlk;                            // Current block in current event

  // == states ==

  /**  */
  Bool_t fConnected;                              // Shows connection status

  /**  .
   */
  Bool_t fStateHasChanged;                        // Indicates, if a sources have changes, so that one has to reconnect.


  //----

  AliTPCPreprocessorOnline* fTPCPre;              // Preprocessor for TPC calibration.

  ClassDef(AliEveHOMERManager, 0); // Manage connections to HLT data-sources.
};

#endif
