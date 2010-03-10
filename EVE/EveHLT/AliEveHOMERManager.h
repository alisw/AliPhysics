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
class TEveManager;
class TEveScene;
class TEveProjectionManager;
class TTimer;

class AliHLTEvePhos;
class AliHLTEveEmcal;
class AliHLTEveTPC;
class AliHLTEveHLT;
class AliHLTEveITS;
class AliHLTEveISSD;
class AliHLTEveISDD;
class AliHLTEveISPD;
class AliHLTEveTRD;
class AliHLTEveAny;
class AliHLTEveMuon;

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
  
  /** Get next event and process it */
  Int_t NextHOMEREvent();

  /** Set flag for event loop */
  void SetEventLoopStarted (Bool_t started) {fEventLoopStarted = started;}

  /** Set flag for showing barrel */
  void SetBarrelFlag(Bool_t flag) { fShowBarrel = flag;}
  /** Set flag for showing muon arm */
  void SetMuonFlag(Bool_t flag) { fShowMuon = flag;}

  /**Set and get the global instance of the Eve manager */
  void SetEveManager(TEveManager * manager) {fEveManager = manager;}
  TEveManager * GetEveManager() const {return fEveManager;}

  /**Set and get the global instance of TGeoManager */
  void SetGeoManager(TGeoManager * manager) {fGeoManager = manager;}
  TGeoManager * GetGeoManager() const {return fGeoManager;}

  /** Set the projection scenes and their managers */
  void SetRPhiManager (TEveProjectionManager * mgr) {fRPhiManager = mgr;}
  void SetRPhiEventScene (TEveScene * scene ) {fRPhiEventScene = scene;}
  void SetRhoZManager(TEveProjectionManager * mgr) {fRhoZManager = mgr;}
  void SetRhoZEventScene(TEveScene * scene ) {fRhoZEventScene = scene;}
  

  /** Start and stop the automatic event loop */
  void StartLoop();
  void StopLoop();
 
private:

  /** copy constructor prohibited */
  AliEveHOMERManager(const AliEveHOMERManager&);

  /** assignment operator prohibited */
  AliEveHOMERManager& operator=(const AliEveHOMERManager&);

  void DestroyDetectorElements();
  
  /** Process block */
  void ProcessBlock(AliHLTHOMERBlockDesc * block);  //Process block
  /** Reset the elements in the display */
  void ResetDisplay();  
  /** Update the display  */
  void UpdateDisplay(); 

  // == sources ==
  AliEveHOMERSourceList* fSrcList;        // List of Connected HOMER Sources

  Int_t fRetryCount;                     //How many times to retry creating source list before giving up
  Int_t fRetrySleeptime;                 //Sleep time between attempt at craeting source list

  TGeoManager * fGeoManager;              //The global TGeoManager instance
  TEveManager * fEveManager;              //The global TEveManager instance
  TEveProjectionManager * fRPhiManager;   //The R - Phi projection scene manager
  TEveProjectionManager * fRhoZManager;   //The Rho- Z projection sene manager
  TEveScene * fRPhiEventScene;            //The R - Phi projection scene
  TEveScene * fRhoZEventScene;            //The Rho - Z projection sene


  TTimer * fTimer;                   //Timer for event loop
  //TTimer * fSourceListTimer;       //Timer for source list loop
 
  AliHLTEvePhos  * fPhosElement;     //Phos eve processor
  AliHLTEveEmcal * fEmcalElement;    //Emcal eve processor
  AliHLTEveTPC   * fTPCElement;      //TPC eve processor
  AliHLTEveHLT   * fHLTElement;      //HLT
  AliHLTEveITS   * fITSElement;      //ITS
  AliHLTEveISPD  * fISPDElement;     //ISPD
  AliHLTEveISSD  * fISSDElement;     //ISSD
  AliHLTEveISDD  * fISDDElement;     //ISDD
  AliHLTEveTRD   * fTRDElement;      //TRD
  AliHLTEveMuon  * fMuonElement;     //MUON
  AliHLTEveAny   * fAnyElement;      //Catch all

  Bool_t fEventLoopStarted;                    // Flag indicating whether the loop is running
  Bool_t fCenterProjectionsAtPrimaryVertex;    // Flag indicating whether to center the projection scenes at primary vertex (as opposed to 0, 0, 0)
  Bool_t fShowBarrel;                               // Display barrel detectors ?
  Bool_t fShowMuon;                                 // Display Muon arm ?
   
  ClassDef(AliEveHOMERManager, 0); // Manage connections to HLT data-sources.

};
#endif
