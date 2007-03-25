// $Id$

void geom_gentle()
{
  TFile f("$REVESYS/alice-data/gentle_geo.root");
  TGeoShapeExtract* gse = (TGeoShapeExtract*) f.Get("Gentle");
  Reve::GeoShapeRnrEl::ImportShapeExtract(gse, 0);
  f.Close();
}
