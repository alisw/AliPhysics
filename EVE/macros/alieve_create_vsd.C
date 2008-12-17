// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

void alieve_create_vsd()
{
  // Invoke from a running alieve, current event will be dumped.

  TEveVSD::DisableTObjectStreamersForVSDStruct();

  AliEveVSDCreator vc;
  vc.SetDebugLevel(2);
  vc.CreateVSD("AliVSD.root");
}
