// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
void alieve_create_vsd()
{
  // Invoke as: aliroot alieve_create_vsd.C

  gSystem->Load("libPhysics");
  gSystem->Load("libEG");
  gSystem->Load("libTreePlayer");
  gSystem->Load("libGed");
  gSystem->Load("libRGL");

  gSystem->Load("libReve");
  gSystem->Load("libAlieve");

  DisablePODTObjectStreamers();

  TGeoManager::Import("geometry.root");

  AliEveVSDCreator vc;
  vc.SetDebugLevel(2);
  vc.CreateVSD(".", 0, "AliVSD.root");
}
