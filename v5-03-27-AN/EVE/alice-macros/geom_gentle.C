// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TEveManager.h>
#include <TEveGeoNode.h>
#include <TEveElement.h>
#include <TEveGeoShape.h>
#include <TEveGeoShapeExtract.h>
#endif

TEveGeoShape* geom_gentle(Bool_t register_as_global=kTRUE)
{
  TFile f("$ALICE_ROOT/EVE/alice-data/gentle_geo.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
  TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
  f.Close();

  TEveElement* elPHOS = gsre->FindChild("PHOS");
  elPHOS->SetRnrState(kTRUE);
  elPHOS->FindChild("PHOS_4")->SetRnrState(kFALSE);
  elPHOS->FindChild("PHOS_5")->SetRnrState(kFALSE);

  if (register_as_global)
  {
    gEve->AddGlobalElement(gsre);
  }

  return gsre;
}

TEveGeoShape* geom_gentle_rphi()
{
  // The resulting geometry is NOT added into the global scene!

  TFile f("$ALICE_ROOT/EVE/alice-data/gentle_rphi_geo.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
  TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
  f.Close();

  TEveElement* elPHOS = gsre->FindChild("PHOS");
  elPHOS->SetRnrState(kTRUE);
  elPHOS->FindChild("PHOS_4")->SetRnrState(kFALSE);
  elPHOS->FindChild("PHOS_5")->SetRnrState(kFALSE);

  return gsre;
}

TEveGeoShape* geom_gentle_rhoz()
{
  // The resulting geometry is NOT added into the global scene!

  TFile f("$ALICE_ROOT/EVE/alice-data/gentle_rhoz_geo.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
  TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
  f.Close();

  return gsre;
}
