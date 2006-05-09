// $Header$

#ifndef REVE_Event_H
#define REVE_Event_H

#include <Reve/Reve.h>
#include <TDirectory.h>
#include <map>

#include <TTree.h>
#include <TFolder.h>
class TBrowser;

namespace Reve {

class VSD;

class EvTree : public TTree
{
public:
  TFolder fFolder;

  EvTree() : TTree() {}
  EvTree(const char* name, const char* title, Int_t splitlevel = 99) :
    TTree(name, title, splitlevel), fFolder("Folder", "Additional event data") {}
  virtual ~EvTree() {}

  virtual void Browse(TBrowser* b);

  ClassDef(EvTree, 1);
};

class Event : public TNamed
{
private:
  void Init();

protected:
  Int_t        fRun;
  Int_t        fEvent;
  TString      fUrl;

  std::map<TString, TString> fTags;

  TDirectory*  fDirectory; //!

  VSD*         fCreator;   //!
  VSD*         fSelector;  //!

public:
  Event();
  Event(Int_t run, Int_t evt, const TString& url=".");
  Event(const TString& url);
  // Static ctors
  static Event* OpenDirectory(const TString& dir_name);
  static Event* OpenDirectory(TDirectory* dir);

  TDirectory* GetDirectory() { return fDirectory; }
  void        SetDirectory(TDirectory* dir);
  TDirectory* MakeDirectory(const Text_t* name, const Text_t* title="");

  TObject*    Get(const Text_t* obj_name) { return fDirectory->Get(obj_name); }

  virtual void Print(Option_t* opt="") const;

  static const TString sVSDHeaderName;

  ClassDef(Event, 1);
}; // endclass Event

}

#endif
