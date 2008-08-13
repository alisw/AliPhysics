void VizDB_scan()
{
  TEvePointSet* ps;

  // Clusters

  ps = new TEvePointSet();
  ps->SetMarkerColor(5);
  ps->SetMarkerSize(0.2);
  ps->SetMarkerStyle(2);
  gEve->InsertVizDBEntry("ITS Clusters", ps);
  
  ps = new TEvePointSet();
  ps->SetMarkerColor(4);
  ps->SetMarkerSize(0.2);
  ps->SetMarkerStyle(2);
  gEve->InsertVizDBEntry("TPC Clusters", ps);

  ps = new TEvePointSet();
  ps->SetMarkerColor(7);
  ps->SetMarkerSize(0.5);
  ps->SetMarkerStyle(4);
  gEve->InsertVizDBEntry("TRD Clusters", ps);

  ps = new TEvePointSet();
  ps->SetMarkerColor(kOrange);
  ps->SetMarkerSize(0.5);
  ps->SetMarkerStyle(4);
  gEve->InsertVizDBEntry("TOF Clusters", ps);
}
