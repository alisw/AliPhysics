void VizDB_scan()
{
  TEvePointSet        *ps = 0;
  TEveStraightLineSet *ls = 0;


  //============================================================================
  // Clusters
  //============================================================================

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

  //============================================================================
  // Primary vertex
  //============================================================================

  // Combined vertex

  ls = new TEveStraightLineSet;
  ls->SetMarkerStyle(2);
  ls->SetMarkerColor(7);
  ls->SetLineColor(7);
  ls->SetLineWidth(3);
  gEve->InsertVizDBEntry("PVTX", ls);

  ls = new TEveStraightLineSet;
  ls->SetMarkerStyle(2);
  ls->SetMarkerColor(7);
  ls->SetLineColor(7);
  ls->SetLineWidth(1);
  gEve->InsertVizDBEntry("PVTX Ellipse", ls);

  ls = new TEveStraightLineSet;
  ls->SetMarkerStyle(2);
  ls->SetMarkerColor(7);
  ls->SetLineColor(7);
  ls->SetLineWidth(1);
  gEve->InsertVizDBEntry("PVTX Box", ls);

  // SPD vertex

  ls = new TEveStraightLineSet;
  ls->SetMarkerStyle(2);
  ls->SetMarkerColor(6);
  ls->SetLineColor(6);
  ls->SetLineWidth(3);
  gEve->InsertVizDBEntry("PVTX SPD", ls);

  ls = new TEveStraightLineSet;
  ls->SetMarkerStyle(2);
  ls->SetMarkerColor(6);
  ls->SetLineColor(6);
  ls->SetLineWidth(1);
  gEve->InsertVizDBEntry("PVTX Ellipse SPD", ls);

  ls = new TEveStraightLineSet;
  ls->SetMarkerStyle(2);
  ls->SetMarkerColor(6);
  ls->SetLineColor(6);
  ls->SetLineWidth(1);
  gEve->InsertVizDBEntry("PVTX Box SPD", ls);

  // TPC vertex

  ls = new TEveStraightLineSet;
  ls->SetMarkerStyle(2);
  ls->SetMarkerColor(5);
  ls->SetLineColor(5);
  ls->SetLineWidth(3);
  gEve->InsertVizDBEntry("PVTX TPC", ls);

  ls = new TEveStraightLineSet;
  ls->SetMarkerStyle(2);
  ls->SetMarkerColor(5);
  ls->SetLineColor(5);
  ls->SetLineWidth(1);
  gEve->InsertVizDBEntry("PVTX Ellipse TPC", ls);

  ls = new TEveStraightLineSet;
  ls->SetMarkerStyle(2);
  ls->SetMarkerColor(5);
  ls->SetLineColor(5);
  ls->SetLineWidth(1);
  gEve->InsertVizDBEntry("PVTX Box TPC", ls);
}
