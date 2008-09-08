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

#include "AliEveHOMERSourceList.h"
#include "AliHLTHOMERSourceDesc.h"
#include "AliHLTHOMERBlockDesc.h"
#include "AliHLTHOMERReader.h"
#include "AliEveHOMERXMLHandler.h"

#include "AliTPCPreprocessorOnline.h"


class AliHLTHOMERLibManager;

class AliEveHOMERManager : public TEveElementList
{
public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** constructor */
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

  /** Set state of a source */
  void SetSourceState( AliHLTHOMERSourceDesc* source, Bool_t state);

  /** Get pointer to source List */
  TList* GetSourceList() { return fSourceList; } // Get pointer to source List

  /*
   * ---------------------------------------------------------------------------------
   *                            Connection Handling - public
   * ---------------------------------------------------------------------------------
   */

  /** Connect to HOMER sources, out of Readout List, which gets created when state has changed */
  Int_t ConnectHOMER();

  /** Disconnect from HOMER sources */
  void DisconnectHOMER();

  /** Reconnect from HOMER sources */
  Int_t ReconnectHOMER();

  /*
   * ---------------------------------------------------------------------------------
   *                            Event Handling - public
   * ---------------------------------------------------------------------------------
   */

  /** Loads the next Event, after being connected */
  Int_t NextEvent();

  /** Get event ID */
  ULong_t GetEventID() { return fEventID; }    // Get event ID

  /** Get pointer to block List */
  TList* GetBlockList() { return fBlockList; } // Get pointer to block List
  
  /*
   * ---------------------------------------------------------------------------------
   *                            Test Realm ....
   * ---------------------------------------------------------------------------------
   */

  /** Still under testing ... */
  void DumpTPCCalib(TString objectName, Bool_t dumpToFile);

  ///////////////////////////////////////////////////////////////////////////////////

protected:

  AliHLTHOMERLibManager* fLibManager;             //! Dynamic loader manager for the HOMER library

  ///////////////////////////////////////////////////////////////////////////////////

private:

  AliEveHOMERManager(const AliEveHOMERManager&);            // Not implemented.
  AliEveHOMERManager& operator=(const AliEveHOMERManager&); // Not implemented.

  /*
   * ---------------------------------------------------------------------------------
   *                            Connection Handling - private
   * ---------------------------------------------------------------------------------
   */

  /** Create a readout list for Hostname and ports */
  void CreateReadoutList( const char** sourceHostnames, UShort_t* sourcePorts, UInt_t &sourceCount);

  /** Checks if already connected to HOMER sources */
  Bool_t IsConnected() { return fConnected; }  // Checks if already connected to HOMER sources

  /** Sets realm ( which can be ACR, GPN, HLT, KIP ) */ 
  void SetRealm( TString s ) { fXMLHandler->SetRealm(s); } // Sets realm ( which can be ACR, GPN, HLT, KIP )

  /* ---------------------------------------------------------------------------------
   *                            Event Handling - private
   * ---------------------------------------------------------------------------------
   */

  /** Create a TList of blocks, which have been readout */
  Int_t CreateBlockList();

  /*
   * ---------------------------------------------------------------------------------
   *                            Block Handling - private
   * ---------------------------------------------------------------------------------
   */

  /** Get Number of blocks in current event */
  ULong_t GetNBlks() { return fNBlks; }                                        // Get Number of blocks in current event

  /** Get pointer to block ndx in current event */
  void* GetBlk( Int_t ndx );

  /** Get pointer to current block in current event */
  void* GetBlk() { return GetBlk( fCurrentBlk ); }                             // Get pointer to current block in current event

  /** Get first block in current event */
  void* GetFirstBlk() { return GetBlk( 0 ); }                                  // Get first block in current event

  /** Get next block in current event */
  void* GetNextBlk() { return GetBlk( ++fCurrentBlk ); }                       // Get next block in current event

  /** Get size of block ndx */
  ULong_t GetBlkSize( Int_t ndx );

  /** Get size of current block */ 
  ULong_t GetBlkSize() { return GetBlkSize( fCurrentBlk ); }                   // Get size of current block 
 
  /** Get origin of block ndx */
  TString GetBlkOrigin( Int_t ndx );

  /** Get origin of current block */
  TString GetBlkOrigin(){ return GetBlkOrigin( fCurrentBlk ); }                // Get origin of current block

  /** Get type of block ndx */
  TString GetBlkType( Int_t ndx ); 

  /** Get type of current block */
  TString GetBlkType() { return GetBlkType( fCurrentBlk ); }                   // Get type of current block

  /** Get specification of block ndx */
  ULong_t GetBlkSpecification( Int_t ndx );

  /** Get specification of current block */
  ULong_t GetBlkSpecification() { return GetBlkSpecification( fCurrentBlk ); } // Get specification of current block

  /** Checks if current Block should was requested */
  Bool_t CheckIfRequested( AliHLTHOMERBlockDesc* block );

  /*
   * ---------------------------------------------------------------------------------
   *                            Members - private
   * ---------------------------------------------------------------------------------
   */

  // == XML handler ==
  AliEveHOMERXMLHandler* fXMLHandler;             //! Handles HLT XML Config Files

  // == sources ==
  TList * fSourceList;                            //! List to HOMER sources

  // == connection ==
  AliHLTHOMERReader* fReader;                     //! Pointer to HOMER reader

  // == blocks ==
  TList * fBlockList;                             //! List to HOMER blocks

  // == events ==
  ULong_t   fNBlks;                               // Number of blockes in current event
  ULong64_t fEventID;                             // EventID of current event
  ULong_t   fCurrentBlk;                          // Current block in current event

  // == states ==
  Bool_t fConnected;                              // Shows connection status
  Bool_t fStateHasChanged;                        // Indicates, if a sources have changes, so that one has to reconnect.

  // == sources ==
  AliEveHOMERSourceList* fSrcList;                // List of HOMER Sources

  //-----------------------------------------------------------------------------------------
  AliTPCPreprocessorOnline* fTPCPre;              // Preprocessor for TPC calibration.

  ClassDef(AliEveHOMERManager, 0); // Manage connections to HLT data-sources.
};

#endif
