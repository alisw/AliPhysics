// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TStyle.h>

#include <AliRawReader.h>
#include <AliEveEventManager.h>
#include <AliEveT0Module.h>
#endif

void t0_raw()
{
    printf("*** RAW T0 ***");
    
  AliRawReader *reader = AliEveEventManager::AssertRawReader();

  gStyle->SetPalette(1, 0);

  AliEveT0Module::LoadRaw(reader);
}
