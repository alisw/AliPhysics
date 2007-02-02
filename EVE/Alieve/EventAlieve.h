// $Header$

#ifndef ALIEVE_Event_H
#define ALIEVE_Event_H

#include <TList.h>

#include <Reve/EventBase.h>

class AliRunLoader;
class AliESD;
class AliESDfriend;

class TFile;
class TTree;

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
  AliESD*       fESD;
  AliESDfriend* fESDfriend;
  Bool_t        fESDfriendExists;

  TList         fNewEventCommands;

  static Bool_t fgUseRunLoader;
  static Bool_t fgUseESDTree;
  static Bool_t fgAvoidExcOnOpen;

public:
  static void Initialize(Bool_t use_runloader=kTRUE, Bool_t use_esd=kTRUE,
			 Bool_t avoid_exc_on_open=kTRUE);

  Event();
  Event(TString path, Int_t ev=0);

  void Open();
  void GotoEvent(Int_t event);
  void NextEvent() { GotoEvent(fEventId + 1); }
  void PrevEvent() { GotoEvent(fEventId - 1); }
  void Close();

  virtual void  AfterNewEventLoaded();

  TList& GetNewEventCommands() { return fNewEventCommands; }
  void   AddNewEventCommand(const Text_t* cmd);

  Int_t         GetEventId()   const { return fEventId; }
  AliRunLoader* GetRunLoader() const { return fRunLoader; }
  TTree*        GetESDTree()   const { return fESDTree; }
  AliESD*       GetESD()       const { return fESD; }
  AliESDfriend* GetESDfriend()       const { return fESDfriend; }
  Bool_t        GetESDfriendExists() const { return fESDfriendExists; }
  virtual const Text_t* GetTitle()   const { return fPath.Data(); }

  static AliRunLoader* AssertRunLoader();
  static AliESD*       AssertESD();
  static AliESDfriend* AssertESDfriend();

  ClassDef(Event, 1);
}; // endclass Event

extern Event* gEvent;

}

#endif
