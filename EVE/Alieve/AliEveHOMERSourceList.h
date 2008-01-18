// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVE_AliEVEHOMERSourceList_H
#define ALIEVE_AliEVEHOMERSourceList_H

#include <TEveElement.h>

#include <TObject.h>

class AliEveHOMERSourceList : public TEveElementList
{
private:
  AliEveHOMERSourceList(const AliEveHOMERSourceList&);            // Not implemented
  AliEveHOMERSourceList& operator=(const AliEveHOMERSourceList&); // Not implemented

protected:

public:
  AliEveHOMERSourceList(const Text_t* n="HOMER Source List", const Text_t* t="");
  virtual ~AliEveHOMERSourceList() {}

  void SelectAll();   // *MENU*
  void DeselectAll(); // *MENU*

  ClassDef(AliEveHOMERSourceList, 1);
}; // endclass AliEveHOMERSourceList

#endif
