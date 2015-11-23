void VizDB_scan()
{

  TEvePointSet        *ps = 0;
  TEveStraightLineSet *ls = 0;
  TEveTrackList       *tl = 0;

  //============================================================================
  // Hits
  //============================================================================

  ps = new TEvePointSet();
  ps->SetMarkerColor(2);
  ps->SetMarkerSize(0.5);
  ps->SetMarkerStyle(2);
  gEve->InsertVizDBEntry("Hits", ps);

  ps = new TEvePointSet();
  ps->SetMarkerColor(2);
  ps->SetMarkerSize(0.5);
  ps->SetMarkerStyle(2);
  gEve->InsertVizDBEntry("SIM Hits ITS", ps);

  ps = new TEvePointSet();
  ps->SetMarkerColor(3);
  ps->SetMarkerSize(0.5);
  ps->SetMarkerStyle(2);
  gEve->InsertVizDBEntry("SIM Hits TPC", ps);

  ps = new TEvePointSet();
  ps->SetMarkerColor(3);
  ps->SetMarkerSize(0.5);
  ps->SetMarkerStyle(4);
  gEve->InsertVizDBEntry("SIM Hits T0", ps);

  ps = new TEvePointSet();
  ps->SetMarkerColor(2);
  ps->SetMarkerSize(0.5);
  ps->SetMarkerStyle(4);
  gEve->InsertVizDBEntry("SIM Hits FMD", ps);

  ps = new TEvePointSet();
  ps->SetMarkerColor(2);
  ps->SetMarkerSize(.5);
  ps->SetMarkerStyle(4);
  gEve->InsertVizDBEntry("SIM Hits ACORDE", ps);

  ps = new TEvePointSet();
  ps->SetMarkerColor(2);
  ps->SetMarkerSize(.5);
  ps->SetMarkerStyle(4);
  gEve->InsertVizDBEntry("SIM Hits EMCAL", ps);

  ps = new TEvePointSet();
  ps->SetMarkerColor(2);
  ps->SetMarkerSize(.5);
  ps->SetMarkerStyle(4);
  gEve->InsertVizDBEntry("SIM Hits PMD", ps);

  ps = new TEvePointSet();
  ps->SetMarkerColor(2);
  ps->SetMarkerSize(.5);
  ps->SetMarkerStyle(4);
  gEve->InsertVizDBEntry("SIM Hits TOF", ps);

  ps = new TEvePointSet();
  ps->SetMarkerColor(7);
  ps->SetMarkerSize(.5);
  ps->SetMarkerStyle(4);
  gEve->InsertVizDBEntry("SIM Hits TRD", ps);

  ps = new TEvePointSet();
  ps->SetMarkerColor(2);
  ps->SetMarkerSize(.5);
  ps->SetMarkerStyle(4);
  gEve->InsertVizDBEntry("SIM Hits VZERO", ps);

  //============================================================================
  // Clusters
  //============================================================================

  ps = new TEvePointSet();
  ps->SetMarkerColor(2);
  ps->SetMarkerSize(0.5);
  ps->SetMarkerStyle(2);
  gEve->InsertVizDBEntry("Clusters", ps);

  ps = new TEvePointSet();
  ps->SetMarkerColor(kBlue);
  ps->SetMarkerSize(0.2);
  ps->SetMarkerStyle(2);
  gEve->InsertVizDBEntry("REC Clusters ITS", ps);

  ps = new TEvePointSet();
  ps->SetMarkerColor(kBlue);
  ps->SetMarkerSize(0.2);
  ps->SetMarkerStyle(2);
  gEve->InsertVizDBEntry("REC Clusters TPC", ps);

  ps = new TEvePointSet();
  ps->SetMarkerColor(7);
  ps->SetMarkerSize(0.5);
  ps->SetMarkerStyle(4);
  gEve->InsertVizDBEntry("REC Clusters TRD", ps);

  ps = new TEvePointSet();
  ps->SetMarkerColor(kOrange+9);
  ps->SetMarkerSize(0.5);
  ps->SetMarkerStyle(4);
  gEve->InsertVizDBEntry("REC Clusters TOF", ps);

  ps = new TEvePointSet();
  ps->SetMarkerColor(4);
  ps->SetMarkerSize(0.2);
  ps->SetMarkerStyle(2);
  gEve->InsertVizDBEntry("REC Clusters HMPID", ps);

  ps = new TEvePointSet();
  ps->SetMarkerColor(4);
  ps->SetMarkerSize(0.5);
  ps->SetMarkerStyle(2);
  gEve->InsertVizDBEntry("REC Clusters PHOS", ps);

  //============================================================================
  // Primary vertex
  //============================================================================

  // Combined vertex

  ls = new TEveStraightLineSet;
  ls->SetMarkerStyle(2);
  ls->SetMarkerColor(7);
  ls->SetLineColor(7);
  ls->SetLineWidth(3);
  gEve->InsertVizDBEntry("REC PVTX", ls);

  ls = new TEveStraightLineSet;
  ls->SetMarkerStyle(2);
  ls->SetMarkerColor(7);
  ls->SetLineColor(7);
  ls->SetLineWidth(1);
  gEve->InsertVizDBEntry("REC PVTX Ellipse", ls);

  ls = new TEveStraightLineSet;
  ls->SetMarkerStyle(2);
  ls->SetMarkerColor(7);
  ls->SetLineColor(7);
  ls->SetLineWidth(1);
  gEve->InsertVizDBEntry("REC PVTX Box", ls);

  // SPD vertex

  ls = new TEveStraightLineSet;
  ls->SetMarkerStyle(2);
  ls->SetMarkerColor(6);
  ls->SetLineColor(6);
  ls->SetLineWidth(3);
  gEve->InsertVizDBEntry("REC PVTX SPD", ls);

  ls = new TEveStraightLineSet;
  ls->SetMarkerStyle(2);
  ls->SetMarkerColor(6);
  ls->SetLineColor(6);
  ls->SetLineWidth(1);
  gEve->InsertVizDBEntry("REC PVTX Ellipse SPD", ls);

  ls = new TEveStraightLineSet;
  ls->SetMarkerStyle(2);
  ls->SetMarkerColor(6);
  ls->SetLineColor(6);
  ls->SetLineWidth(1);
  gEve->InsertVizDBEntry("REC PVTX Box SPD", ls);

  // TPC vertex

  ls = new TEveStraightLineSet;
  ls->SetMarkerStyle(2);
  ls->SetMarkerColor(5);
  ls->SetLineColor(5);
  ls->SetLineWidth(3);
  gEve->InsertVizDBEntry("REC PVTX TPC", ls);

  ls = new TEveStraightLineSet;
  ls->SetMarkerStyle(2);
  ls->SetMarkerColor(5);
  ls->SetLineColor(5);
  ls->SetLineWidth(1);
  gEve->InsertVizDBEntry("REC PVTX Ellipse TPC", ls);

  ls = new TEveStraightLineSet;
  ls->SetMarkerStyle(2);
  ls->SetMarkerColor(5);
  ls->SetLineColor(5);
  ls->SetLineWidth(1);
  gEve->InsertVizDBEntry("REC PVTX Box TPC", ls);


  //============================================================================
  // Tracks
  //============================================================================

  tl = new TEveTrackList("ESD Tracks");
  tl->SetLineStyle(6);
  tl->SetMainColor(1);
  tl->SetLineWidth(1);
  gEve->InsertVizDBEntry("REC Tracks", tl);

  tl = new TEveTrackList("ESD Tracks MI");
  tl->SetLineStyle(6);
  tl->SetMainColor(1);
  tl->SetLineWidth(1);
  gEve->InsertVizDBEntry("REC Tracks MI", tl);


  // esd_tracks_by_category()

  tl = new TEveTrackList("Sigma < 3");
  tl->SetLineStyle(6);
  tl->SetLineColor(1);
  tl->SetLineWidth(1);
  gEve->InsertVizDBEntry("REC Tracks ByCat 0", tl);

  tl = new TEveTrackList("3 < Sigma < 5");
  tl->SetLineStyle(6);
  tl->SetLineColor(1);
  tl->SetLineWidth(1);
  gEve->InsertVizDBEntry("REC Tracks ByCat 1", tl);

  tl = new TEveTrackList("5 < Sigma");
  tl->SetLineStyle(6);
  tl->SetLineColor(1);
  tl->SetLineWidth(1);
  gEve->InsertVizDBEntry("REC Tracks ByCat 2", tl);

  tl = new TEveTrackList("no ITS refit; Sigma < 5");
  tl->SetLineStyle(6);
  tl->SetLineColor(1);
  tl->SetLineWidth(1);
  gEve->InsertVizDBEntry("REC Tracks ByCat 3", tl);

  tl = new TEveTrackList("no ITS refit; Sigma > 5");
  tl->SetLineStyle(6);
  tl->SetLineColor(1);
  tl->SetLineWidth(1);
  gEve->InsertVizDBEntry("REC Tracks ByCat 4", tl);

  tl = new TEveTrackList("no TPC refit");
  tl->SetLineStyle(6);
  tl->SetLineColor(1);
  tl->SetLineWidth(1);
  gEve->InsertVizDBEntry("REC Tracks ByCat 5", tl);

  tl = new TEveTrackList("ITS stand-alone");
  tl->SetLineStyle(6);
  tl->SetLineColor(1);
  tl->SetLineWidth(1);
  gEve->InsertVizDBEntry("REC Tracks ByCat 6", tl);


  // esd_tracks_by_anal_cuts()

  tl = new TEveTrackList("Passed");
  tl->SetLineStyle(6);
  tl->SetMainColor(1);
  tl->SetLineWidth(1);
  gEve->InsertVizDBEntry("REC Tracks ByAnCuts Passed", tl);

  tl = new TEveTrackList("Rejected");
  tl->SetLineStyle(6);
  tl->SetMainColor(1);
  tl->SetLineWidth(1);
  gEve->InsertVizDBEntry("REC Tracks ByAnCuts Rejected", tl);


  //============================================================================
  // SPD tracklets
  //============================================================================

  tl = new TEveTrackList("Good");
  tl->SetLineStyle(6);
  tl->SetMainColor(1);
  tl->SetLineWidth(1);
  gEve->InsertVizDBEntry("REC Tracklet Good", tl);

  tl = new TEveTrackList("Bad");
  tl->SetLineStyle(6);
  tl->SetMainColor(1);
  tl->SetLineWidth(1);
  gEve->InsertVizDBEntry("REC Tracklet Bad", tl);
}
