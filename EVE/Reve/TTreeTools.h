// $Header

#ifndef REVE_TTreeTools_H
#define REVE_TTreeTools_H

#include <TSelectorDraw.h>
#include <TEventList.h>

/**************************************************************************/
// TSelectorToEventList
/**************************************************************************/

class TSelectorToEventList : public TSelectorDraw
{
  TSelectorToEventList(const TSelectorToEventList&);            // Not implemented
  TSelectorToEventList& operator=(const TSelectorToEventList&); // Not implemented

protected:
  TEventList* fEvList;
  TList       fInput;
public:
  TSelectorToEventList(TEventList* evl, const Text_t* sel);

  virtual Int_t  Version() const { return 1; }
  virtual Bool_t Process(Long64_t entry);

  ClassDef(TSelectorToEventList, 1)
};

/**************************************************************************/
// TTreeQuery
/**************************************************************************/

class TTreeQuery : public TEventList
{
public:
  TTreeQuery() : TEventList() {}

  Int_t Select(TTree* t, const Text_t* selection);

  ClassDef(TTreeQuery, 1)
};

/**************************************************************************/
// TPointSelectorConsumer, TPointSelector
/**************************************************************************/

class TPointSelector;

class TPointSelectorConsumer
{
public:
  enum TreeVarType_e { TVT_XYZ, TVT_RPhiZ };

protected:
  TreeVarType_e fSourceCS; // Coordinate-System of the source tree variables

public:
  TPointSelectorConsumer(TreeVarType_e cs=TVT_XYZ) :fSourceCS(cs) {}
  virtual ~TPointSelectorConsumer() {}

  virtual void TakeAction(TSelectorDraw*) = 0;

  TreeVarType_e GetSourceCS() const  { return fSourceCS; }
  void SetSourceCS(TreeVarType_e cs) { fSourceCS = cs; }

  ClassDef(TPointSelectorConsumer, 1);
};

class TPointSelector : public TSelectorDraw
{
  TPointSelector(const TPointSelector&);            // Not implemented
  TPointSelector& operator=(const TPointSelector&); // Not implemented

protected:
  TTree                  *fTree;
  TPointSelectorConsumer *fConsumer;

  TString                 fVarexp;
  TString                 fSelection;

  TList                   fInput;

public:
  TPointSelector(TTree* t=0, TPointSelectorConsumer* c=0,
		 const Text_t* vexp="", const Text_t* sel="");
  virtual ~TPointSelector() {}

  virtual Long64_t Select(const Text_t* selection=0);
  virtual Long64_t Select(TTree* t, const Text_t* selection=0);
  virtual void  TakeAction();


  TTree* GetTree() const   { return fTree; }
  void   SetTree(TTree* t) { fTree = t; }

  TPointSelectorConsumer* GetConsumer() const { return fConsumer; }
  void SetConsumer(TPointSelectorConsumer* c) { fConsumer = c; }

  const Text_t* GetVarexp() const { return fVarexp; }
  void SetVarexp(const Text_t* v) { fVarexp = v; }

  const Text_t* GetSelection() const { return fSelection; }
  void SetSelection(const Text_t* s) { fSelection = s; }

  ClassDef(TPointSelector, 1);
};

/**************************************************************************/
/**************************************************************************/

#endif
