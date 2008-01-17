// $Id$

TEveGeoShape* geom_gentle()
{
  TFile f("$REVESYS/alice-data/gentle_geo.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
  TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse, 0);
  f.Close();

  return gsre;
}
