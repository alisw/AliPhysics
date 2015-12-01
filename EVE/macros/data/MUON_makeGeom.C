// $Id$

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

/// \ingroup evemacros
/// \file MUON_makeGeom.C
///
/// \author B. Vulpescu, LPC
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TGeoManager.h>

#include <AliMpCDB.h>
#include <AliRun.h>
#endif

void MUON_makeGeom()
{
  AliMpCDB::LoadMpSegmentation2(); 
  gAlice->SetConfigFunction("$ALICE_ROOT/MUON/Config.C");
  //gAlice->Init("Config.C");

  gGeoManager->Export("geometry.root");

}
