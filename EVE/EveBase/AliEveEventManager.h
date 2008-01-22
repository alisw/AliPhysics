// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVE_EventAlieve_H
#define ALIEVE_EventAlieve_H

#include <TEveEventManager.h>

class AliRunLoader;
class AliESDEvent;
class AliESDfriend;

class AliMagF;

class TFile;
class TTree;
class TGeoManager;


class AliEveEventManager : public TEveEventManager
{
private:
  AliEveEventManager(const AliEveEventManager&);            // Not implemented
  AliEveEventManager& operator=(const AliEveEventManager&); // Not implemented

protected:
  TString       fPath;			// URL to event-data.
  Int_t         fEventId;		// Id of current event.

  AliRunLoader* fRunLoader;		// Run loader.

  TFile*        fESDFile;		// ESD file.
  TTree*        fESDTree;		// ESD tree.
  AliESDEvent*  fESD;			// ESDEvent object.
  AliESDfriend* fESDfriend;		// ESDfriend object.
  Bool_t        fESDfriendExists;	// Flag specifying if ESDfriend was found during opening of the event-data.

  static TString  fgCdbUri;		// Global URI to CDB.
  static Bool_t   fgAssertRunLoader;	// Global flag specifying if AliRunLoader must be asserted during opening of the event-data.
  static Bool_t   fgAssertESD;		// Global flag specifying if ESDEvent must be asserted during opening of the event-data.

  static AliMagF* fgMagField;		// Global pointer to magneti field.

public:
  static void SetCdbUri(const Text_t* cdb) { if (cdb) fgCdbUri = cdb; }
  static void SetAssertElements(Bool_t assertRunloader, Bool_t assertEsd)
  { fgAssertRunLoader = assertRunloader; fgAssertESD = assertEsd; }

  AliEveEventManager();
  AliEveEventManager(TString path, Int_t ev=0);

  virtual void Open();
  virtual void GotoEvent(Int_t event);
  virtual void NextEvent() { GotoEvent(fEventId + 1); }
  virtual void PrevEvent() { GotoEvent(fEventId - 1); }
  virtual void Close();

  Int_t         GetEventId()   const { return fEventId; }
  AliRunLoader* GetRunLoader() const { return fRunLoader; }
  TTree*        GetESDTree()   const { return fESDTree; }
  AliESDEvent*  GetESD()       const { return fESD; }
  AliESDfriend* GetESDfriend()       const { return fESDfriend; }
  Bool_t        GetESDfriendExists() const { return fESDfriendExists; }
  virtual const Text_t* GetTitle()   const { return fPath.Data(); }

  static AliRunLoader* AssertRunLoader();
  static AliESDEvent*  AssertESD();
  static AliESDfriend* AssertESDfriend();

  static AliMagF*      AssertMagField();

  static TGeoManager*  AssertGeometry();

  ClassDef(AliEveEventManager, 1);
}; // endclass AliEveEventManager

extern AliEveEventManager* gEvent;

#endif
