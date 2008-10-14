// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

// MT comment, 14.10.2008 - This is wrong and not in sync with AliEveT0Module.
// It can not work.

void t0_raw()
{
  AliRunLoader* rl =  AliEveEventManager::AssertRunLoader();

  gStyle->SetPalette(1, 0);

  AliEveT0Module::LoadRaw("raw.root",ievt);
}
