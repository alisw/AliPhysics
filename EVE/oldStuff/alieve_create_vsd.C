// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

void alieve_create_vsd(Int_t nEvents=1, Int_t minTracks=0)
{
  // Invoke from a running alieve.
  // nEvents will be domped starting from current one.
  // If minTracks is set at least that many ESD tracks must exist.

  TEveVSD::DisableTObjectStreamersForVSDStruct();

  AliEveVSDCreator vc;
  vc.SetDebugLevel(2);

  Int_t nDone = 0;
  while (nDone < nEvents)
  {
    if (minTracks)
    {
      AliESDEvent* esd = AliEveEventManager::Instance()->AssertESD();
      while (esd->GetNumberOfTracks() < minTracks)
      {
	AliEveEventManager::Instance()->NextEvent();
	esd = AliEveEventManager::Instance()->AssertESD();
      }
    }
    vc.CreateVSD("AliVSD.root");
    ++nDone;

    AliEveEventManager::Instance()->NextEvent();
  }
}
