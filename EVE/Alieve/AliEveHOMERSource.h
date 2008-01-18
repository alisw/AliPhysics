// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 * 
 **************************************************************************/

#ifndef ALIEVE_AliEVEHOMERSource_H
#define ALIEVE_AliEVEHOMERSource_H

#include <TEveElement.h>

#include <TNamed.h>

class AliHLTHOMERSourceDesc;

class AliEveHOMERSource : public TEveElement,
  public TNamed
{
private:
  AliEveHOMERSource(const AliEveHOMERSource&);            // Not implemented
  AliEveHOMERSource& operator=(const AliEveHOMERSource&); // Not implemented

protected:
  AliHLTHOMERSourceDesc *fSource;

public:
  AliEveHOMERSource(const Text_t* n="HOMER Source", const Text_t* t="");
  AliEveHOMERSource(AliHLTHOMERSourceDesc* src, const Text_t* n="HOMER Source", const Text_t* t="");
  virtual ~AliEveHOMERSource() {}

  AliHLTHOMERSourceDesc* GetSource() const   { return fSource; }
  void SetSource(AliHLTHOMERSourceDesc* src) { fSource = src; }

  ClassDef(AliEveHOMERSource, 1);
}; // endclass AliEveHOMERSource

#endif
