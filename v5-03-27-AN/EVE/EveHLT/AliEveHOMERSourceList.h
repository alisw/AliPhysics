//-*- Mode: C++ -*-

// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveAliEVEHOMERSourceList_H
#define AliEveAliEVEHOMERSourceList_H

#include <TEveElement.h>

#include <TObject.h>

#include <map>

class AliEveHOMERManager;
class AliEveHOMERSourceMap;

class AliEveHOMERSourceList : public TEveElementList
{
public:
  AliEveHOMERSourceList(const Text_t* n="HOMER Source List", const Text_t* t="");
  virtual ~AliEveHOMERSourceList();

  // void InitMap(TList* srcHandles, ESourceGrouping_e

  AliEveHOMERManager* GetManager() const { return fManager; }
  void SetManager(AliEveHOMERManager* m) { fManager = m; }

  Bool_t GetSelectedSources();

  void CreateByDet();  // *MENU*
  void CreateByType(); // *MENU*

  void RebuildSourceReps();

  //void SelectAll();   // *MENU*
  //void DeselectAll(); // *MENU*

protected:
  //SourceMap_t       fByType;
  //ESourceGrouping_e fView;
  //Bool_t            fDefaultState;

  AliEveHOMERManager   *fManager;
  AliEveHOMERSourceMap *fSrcMap;

private:
  AliEveHOMERSourceList(const AliEveHOMERSourceList&);            // Not implemented
  AliEveHOMERSourceList& operator=(const AliEveHOMERSourceList&); // Not implemented

  ClassDef(AliEveHOMERSourceList, 0); // Interface to a list of HOMER sourcces.
};

#endif
