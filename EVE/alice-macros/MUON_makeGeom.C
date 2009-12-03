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

{
  AliMpCDB::LoadMpSegmentation2(); 
  gAlice->Init("$ALICE_ROOT/MUON/Config.C");
  //gAlice->Init("Config.C");

  gGeoManager->Export("geometry.root");

}
