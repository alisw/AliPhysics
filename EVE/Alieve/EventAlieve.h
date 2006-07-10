// $Header$

#ifndef ALIEVE_Event_H
#define ALIEVE_Event_H

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
  void Init();

protected:
  TString       fPath;
  Int_t         fEventId;

  AliRunLoader* fRunLoader;

  TFile*        fESDFile;
  TTree*        fESDTree;
  AliESD*       fESD;
  TFile*        fESDfriendFile;
  TTree*        fESDfriendTree;
  AliESDfriend* fESDfriend;

  static Bool_t fgUseRunLoader;
  static Bool_t fgUseESDTree;

public:
  static void Initialize(Bool_t use_runloader=true, Bool_t use_esd=true);

  Event();
  Event(TString path, Int_t ev=0);

  void Open();
  void Close();

  Int_t         GetEventId()   const { return fEventId; }
  AliRunLoader* GetRunLoader() const { return fRunLoader; }
  TTree*        GetESDTree()   const { return fESDTree; }
  AliESD*       GetESD()       const { return fESD; }
  TTree*        GetESDfriendTree() const { return fESDfriendTree; }
  AliESDfriend* GetESDfriend()     const { return fESDfriend; }

  virtual const Text_t* GetTitle() const { return fPath.Data(); }

  static AliRunLoader* AssertRunLoader();
  static AliESD*       AssertESD();
  static AliESDfriend* AssertESDfriend();

  ClassDef(Event, 1);
}; // endclass Event

extern Event* gEvent;

}

#endif
