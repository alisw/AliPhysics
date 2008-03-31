// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveKineTools_H
#define AliEveKineTools_H

#include <TEveUtil.h>
#include <TObject.h>

class TTree;
class AliStack;

class AliEveKineTools
{
public:
  AliEveKineTools() {}
  virtual ~AliEveKineTools() {}

  void SetDaughterPathMarks(TEveElement* cont, AliStack* stack, Bool_t recurse=kFALSE);
  void SetTrackReferences  (TEveElement* cont, TTree* treeTR=0, Bool_t recurse=kFALSE);
  void SortPathMarks       (TEveElement* cont, Bool_t recurse=kFALSE);

private:
  AliEveKineTools(const AliEveKineTools&);            // Not implemented
  AliEveKineTools& operator=(const AliEveKineTools&); // Not implemented

  ClassDef(AliEveKineTools, 0); // Tools for import of kinematics.
};

#endif
