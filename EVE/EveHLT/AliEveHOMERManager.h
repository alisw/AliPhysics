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
    @author Jochen Thaeder
    @date
    @brief  Manager for HOMER in offline. Inherits most functionalitye
    from AliHLTHOMERManager, with small additions for alieve interface
*/
#include <AliHLTHOMERManager.h>
#include <TEveEventManager.h>
#include "AliEveHOMERSourceList.h"


class AliEveHOMERManager : public AliHLTHOMERManager, public TEveElementList
{
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


  /*
   * ---------------------------------------------------------------------------------
   *                            Source Handling - public
   * ---------------------------------------------------------------------------------
   */

  /** Create Sources List from HOMER-Proxy 
   *  @return 0 on success, <0 for failure, 1 for no active service
   */
  Int_t CreateEveSourcesList();

  Int_t CreateEveSourcesListLoop();
  
  Int_t ConnectEVEtoHOMER(TString detector="ALL");

  void SetRetryCount(Int_t count, Int_t sleeptime) { fRetryCount = count; fRetrySleeptime = sleeptime; }

  ///////////////////////////////////////////////////////////////////////////////////

private:

  /** copy constructor prohibited */
  AliEveHOMERManager(const AliEveHOMERManager&);

  /** assignment operator prohibited */
  AliEveHOMERManager& operator=(const AliEveHOMERManager&);

  /*
   * ---------------------------------------------------------------------------------
   *                            Members - private
   * ---------------------------------------------------------------------------------
   */

  // == sources ==
  AliEveHOMERSourceList* fSrcList;                // List of Connected HOMER Sources

  Int_t fRetryCount;

  Int_t fRetrySleeptime;

  ClassDef(AliEveHOMERManager, 0); // Manage connections to HLT data-sources.

};



#endif
