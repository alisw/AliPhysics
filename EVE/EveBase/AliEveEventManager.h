// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveEventManager_H
#define AliEveEventManager_H

#include <TEveEventManager.h>
#include <TQObject.h>
#include <TObjArray.h>

class AliEveMacroExecutor;

class AliRunLoader;
class AliESDEvent;
class AliESDfriend;
class AliAODEvent;
class AliRawReader;

class AliMagF;

class TFile;
class TTree;
class TGeoManager;

//==============================================================================
//
// AliEveEventManager
//
// Interface to ALICE event-data (RunLoader, ESD), magnetic field and
// geometry.
//


class AliEveEventManager : public TEveEventManager,
                           public TQObject
{
public:
  static void SetESDFileName(const TString& esd);
  static void SetAODFileName(const TString& aod);
  static void AddAODfriend  (const TString& friendFileName);
  static void SetRawFileName(const TString& raw);
  static void SetCdbUri     (const TString& cdb);
  static void SetAssertElements(Bool_t assertRunloader, Bool_t assertEsd,
				Bool_t assertAod, Bool_t assertRaw);

  AliEveEventManager(const TString& name="Event");
  AliEveEventManager(const TString& name, const TString& path, Int_t ev=0);
  virtual ~AliEveEventManager();

  virtual void  Open();
  void          SetEvent(AliRunLoader *runLoader, AliRawReader *rawReader, AliESDEvent *esd, AliESDfriend *esdf);
  virtual Int_t GetMaxEventId(Bool_t refreshESD=kFALSE) const;
  virtual void  GotoEvent(Int_t event);
  virtual void  NextEvent();
  virtual void  PrevEvent();
  virtual void  Close();
  Bool_t        FindNextByTrigger(Int_t& i);
  Bool_t        FindPrevByTrigger(Int_t& i);


  Int_t         GetEventId()         const { return fEventId; }
  AliRunLoader* GetRunLoader()       const { return fRunLoader; }
  TFile*        GetESDFile()         const { return fESDFile; }
  TTree*        GetESDTree()         const { return fESDTree; }
  AliESDEvent*  GetESD()             const { return fESD;     }
  AliESDfriend* GetESDfriend()       const { return fESDfriend; }
  Bool_t        GetESDfriendExists() const { return fESDfriendExists; }
  TFile*        GetAODFile()         const { return fAODFile; }
  TTree*        GetAODTree()         const { return fAODTree; }
  AliAODEvent*  GetAOD()             const { return fAOD;     }
  virtual const Text_t* GetTitle()   const { return fPath.Data(); }
  TString       GetEventInfoHorizontal() const;
  TString       GetEventInfoVertical()   const;

  static Bool_t HasRunLoader();
  static Bool_t HasESD();
  static Bool_t HasESDfriend();
  static Bool_t HasAOD();
  static Bool_t HasRawReader();

  static AliRunLoader* AssertRunLoader();
  static AliESDEvent*  AssertESD();
  static AliESDfriend* AssertESDfriend();
  static AliAODEvent*  AssertAOD();
  static AliRawReader* AssertRawReader();

  static TGeoManager*  AssertGeometry();

  static AliEveEventManager* AddDependentManager(const TString& name, const TString& path);
  static AliEveEventManager* GetDependentManager(const TString& name);

  static AliEveEventManager* GetMaster();
  static AliEveEventManager* GetCurrent();

  Double_t      GetAutoLoadTime()        const { return fAutoLoadTime; }
  Bool_t        GetAutoLoad()            const { return fAutoLoad;     }
  void          SetAutoLoadTime(Float_t time);
  void          SetAutoLoad(Bool_t autoLoad);
  void          AutoLoadNextEvent();

  Bool_t        GetSelectOnTriggerType()     const { return fSelectOnTriggerType; }
  TString       GetTriggerType()             const { return fTriggerType; }
  void          SetTriggerType(const TString& triggertype) { fTriggerType = triggertype; }
  void          SetSelectOnTriggerType(Bool_t sel)         { fSelectOnTriggerType = sel; }

  Bool_t        AreEventFilesOpened()    const { return fIsOpen;       }
  Bool_t        IsEventAvailable()       const { return fHasEvent;     }
  Bool_t        IsUnderExternalControl() const { return fExternalCtrl; }

  virtual void  AfterNewEventLoaded();
  void          NewEventLoaded();      // *SIGNAL*

  AliEveMacroExecutor* GetExecutor() const { return fExecutor; }

protected:
  TString       fPath;			// URL to event-data.
  Int_t         fEventId;		// Id of current event.

  AliRunLoader* fRunLoader;		// Run loader.

  TFile        *fESDFile;		// ESD file.
  TTree        *fESDTree;		// ESD tree.
  AliESDEvent  *fESD;			// ESDEvent object.
  AliESDfriend *fESDfriend;		// ESDfriend object.
  Bool_t        fESDfriendExists;	// Flag specifying if ESDfriend was found during opening of the event-data.
  TFile        *fAODFile;		// AOD file.
  TTree        *fAODTree;		// AOD tree.
  AliAODEvent  *fAOD;			// AODEvent object.

  AliRawReader *fRawReader;             // Raw-data reader.

  Bool_t        fAutoLoad;              // Automatic loading of events (online)
  Float_t       fAutoLoadTime;          // Auto-load time in seconds
  TTimer       *fAutoLoadTimer;         // Timer for automatic event loading

  Bool_t        fIsOpen;                // Are event-files opened.
  Bool_t        fHasEvent;              // Is an event available.
  Bool_t        fExternalCtrl;          // Are we under external event-loop.

  Bool_t        fSelectOnTriggerType;   // Whether to select on trigger-type.
  TString       fTriggerType;           // Trigger-type to select on.

  AliEveMacroExecutor *fExecutor;       // Executor for std macros

  TList        *fSubManagers;           // Dependent event-managers, used for event embedding.

  static TString  fgESDFileName;        // Name by which to open ESD.
  static TString  fgAODFileName;        // Name by which to open AOD.
  static TString  fgRawFileName;        // Name by which to open raw-data file.
  static TString  fgCdbUri;		// Global URI to CDB.
  static Bool_t   fgAssertRunLoader;	// Global flag specifying if AliRunLoader must be asserted during opening of the event-data.
  static Bool_t   fgAssertESD;		// Global flag specifying if ESDEvent must be asserted during opening of the event-data.
  static Bool_t   fgAssertAOD;		// Global flag specifying if AODEvent must be asserted during opening of the event-data.
  static Bool_t   fgAssertRaw;		// Global flag specifying if raw-data presence must be asserted during opening of the event-data.

  static TList   *fgAODfriends;         // Global list of AOD friend names to be attached during opening of the event-data (empty by default).

private:
  AliEveEventManager(const AliEveEventManager&);            // Not implemented
  AliEveEventManager& operator=(const AliEveEventManager&); // Not implemented

  void InitInternals();
  void StartAutoLoadTimer();
  void StopAutoLoadTimer();

  Bool_t fAutoLoadTimerRunning; // State of auto-load timer.

  static AliEveEventManager* fgMaster;
  static AliEveEventManager* fgCurrent;

  ClassDef(AliEveEventManager, 0); // Interface for getting all event components in a uniform way.
};

#endif
