// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 * 
 **************************************************************************/

#include "AliEveHOMERSource.h"

//______________________________________________________________________
// AliEveHOMERSource
//

ClassImp(AliEveHOMERSource)

AliEveHOMERSource::AliEveHOMERSource(const Text_t* n, const Text_t* t) :
  TEveElement(),
  TNamed(n, t),
  fSource(0)
{}

AliEveHOMERSource::AliEveHOMERSource(AliHLTHOMERSourceDesc* src, const Text_t* n, const Text_t* t) :
  TEveElement(),
  TNamed(n, t),
  fSource(src)
{}

