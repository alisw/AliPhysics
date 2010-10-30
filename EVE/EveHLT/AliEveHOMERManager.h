//-*- Mode: C++ -*-

// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVEHOMERMANAGER_H
#define ALIEVEHOMERMANAGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliEveHOMERManager.h
    @author Jochen Thaeder, Svein Lindal
    @date
    @brief  Manager for HOMER in HLT . Inherits some functionalitye
    from AliHLTHOMERManager, mainly from TEveEventManager
*/
#include "AliHLTHOMERManager.h"
#include <TEveEventManager.h>
#include <TGeoManager.h>

class AliEveHOMERSourceList;
class TString;
class TTimer;

class AliEveHOMERManager : public TEveElementList, public AliHLTHOMERManager {

public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** constructor */
  AliEveHOMERManager();

  /** destructor */
  virtual ~AliEveHOMERManager();

  /** Connect to avahi and get the list for sources */
  Int_t CreateEveSourcesList();

  /** Keep on looking for sources until some are found */
  Int_t CreateEveSourcesListLoop();
  
  /** Connect to the sources found */
  Int_t ConnectEVEtoHOMER(TString detector="ALL");

  /** Delete current connections to sources and reconnect */
  Int_t ReConnectHOMER( TString detector="" );

  /** Set the retry count for source list loop */
  void SetRetryCount(Int_t count, Int_t sleeptime) { fRetryCount = count; fRetrySleeptime = sleeptime; }
  
  /** Get next event from the readers */
  TList * NextHOMEREvent();

  void StartEveSourceListLoop();
  void StopEveSourceListLoop();

 
private:

  /** copy constructor prohibited */
  AliEveHOMERManager(const AliEveHOMERManager&);

  /** assignment operator prohibited */
  AliEveHOMERManager& operator=(const AliEveHOMERManager&);

  // == sources ==
  AliEveHOMERSourceList* fSrcList;        // List of Connected HOMER Sources


  Int_t fRetryCount;                     //How many times to retry creating source list before giving up
  Int_t fRetrySleeptime;                 //Sleep time between attempt at craeting source list

  
  TTimer * fSourceListTimer;             //Timer to attempt source list creation!
    
  ClassDef(AliEveHOMERManager, 0); // Manage connections to HLT data-sources.
  

};
#endif
