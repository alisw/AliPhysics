// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

TEveGeoShape* geom_gentle(Bool_t register_as_global=kTRUE)
{
  TFile f("$ALICE_ROOT/EVE/alice-data/gentle_geo.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
  TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
  f.Close();

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
