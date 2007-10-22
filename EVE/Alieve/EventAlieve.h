// $Header$

#ifndef ALIEVE_Event_H
#define ALIEVE_Event_H

#include <Reve/EventBase.h>

class AliRunLoader;
class AliESDEvent;
class AliESDfriend;

class AliMagF;

class TFile;
class TTree;
class TGeoManager;

namespace Alieve {

class Event : public Reve::EventBase
{
private:
  Event(const Event&);            // Not implemented
  Event& operator=(const Event&); // Not implemented

protected:
  TString       fPath;
  Int_t         fEventId;

  AliRunLoader* fRunLoader;

  TFile*        fESDFile;
  TTree*        fESDTree;
  AliESDEvent*  fESD;
  AliESDfriend* fESDfriend;
  Bool_t        fESDfriendExists;

  static TString  fgCdbUri;
  static Bool_t   fgAssertRunLoader;
  static Bool_t   fgAssertESDTree;

  static AliMagF* fgMagField;

public:
  static void SetCdbUri(const Text_t* cdb) { if (cdb) fgCdbUri = cdb; }
  static void SetAssertElements(Bool_t assertRunloader, Bool_t assertEsd)
  { fgAssertRunLoader = assertRunloader; fgAssertESDTree = assertEsd; }

  Event();
  Event(TString path, Int_t ev=0);

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

  ClassDef(Event, 1);
}; // endclass Event

extern Event* gEvent;

}

#endif
