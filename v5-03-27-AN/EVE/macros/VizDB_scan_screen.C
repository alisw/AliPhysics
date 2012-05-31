void VizDB_scan_screen()
{

  TEvePointSet        *ps = 0;
  TEveStraightLineSet *ls = 0;

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
  ps->SetMarkerColor(5);
  ps->SetMarkerSize(0.2);
  ps->SetMarkerStyle(2);
  gEve->InsertVizDBEntry("REC Clusters ITS", ps);

  ps = new TEvePointSet();
  ps->SetMarkerColor(4);
  ps->SetMarkerSize(0.2);
  ps->SetMarkerStyle(2);
  gEve->InsertVizDBEntry("REC Clusters TPC", ps);

  ps = new TEvePointSet();
  ps->SetMarkerColor(7);
  ps->SetMarkerSize(0.5);
  ps->SetMarkerStyle(4);
  gEve->InsertVizDBEntry("REC Clusters TRD", ps);

  ps = new TEvePointSet();
  ps->SetMarkerColor(kOrange);
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

  //Tracks

  tl = new TEveTrackList("ESD Tracks");
  tl->SetLineStyle(6);
  tl->SetMainColor(1);
  tl->SetLineWidth(1);
  gEve->InsertVizDBEntry("REC Tracks",tl);

  tl = new TEveTrackList("ESD Tracks MI");
  tl->SetLineStyle(6);
  tl->SetMainColor(1);
  tl->SetLineWidth(1);
  gEve->InsertVizDBEntry("REC Tracks MI",tl);

  TEveElementList* el = new TEveElementList("ESD Tracks by category");
  TEveTrackList *tltemp[7];
  tltemp[0] = new TEveTrackList("Sigma < 3");
  tltemp[0]->SetLineStyle(6);
  tltemp[0]->SetLineColor(1);
  tltemp[0]->SetLineWidth(1);
  el->AddElement(tltemp[0]);

  tltemp[1] = new TEveTrackList("3 < Sigma < 5");
  tltemp[1]->SetLineStyle(6);
  tltemp[1]->SetLineColor(1);
  tltemp[1]->SetLineWidth(1);
  el->AddElement(tltemp[1]);

  tltemp[2] = new TEveTrackList("5 < Sigma");
  tltemp[2]->SetLineStyle(6);
  tltemp[2]->SetLineColor(1);
  tltemp[2]->SetLineWidth(1);
  el->AddElement(tltemp[2]);

  tltemp[3] = new TEveTrackList("no ITS refit; Sigma < 5");
  tltemp[3]->SetLineStyle(6);
  tltemp[3]->SetLineColor(1);
  tltemp[3]->SetLineWidth(1);
  el->AddElement(tltemp[3]);

  tltemp[4] = new TEveTrackList("no ITS refit; Sigma > 5");
  tltemp[4]->SetLineStyle(6);
  tltemp[4]->SetLineColor(1);
  tltemp[4]->SetLineWidth(1);
  el->AddElement(tltemp[4]);

  tltemp[5] = new TEveTrackList("no TPC refit");
  tltemp[5]->SetLineStyle(6);
  tltemp[5]->SetLineColor(1);
  tltemp[5]->SetLineWidth(1);
  el->AddElement(tltemp[5]);

  tltemp[6] = new TEveTrackList("ITS stand-alone");
  tltemp[6]->SetLineStyle(6);
  tltemp[6]->SetLineColor(1);
  tltemp[6]->SetLineWidth(1);
  el->AddElement(tltemp[6]);

  el->SetVizTag("ESD Tracks by category");
  gEve->AddElement(el);

  TEveElementList* el = new TEveElementList("ESD Tracks by anal cuts");
  TEveTrackList *tlac[2];
  tlac[0] = new TEveTrackList("Passed");
  tlac[0]->SetLineStyle(6);
  tlac[0]->SetMainColor(1);
  tlac[0]->SetLineWidth(1);
  el->AddElement(tlac[0]);

  tlac[1] = new TEveTrackList("Rejected");
  tlac[1]->SetLineStyle(6);
  tlac[1]->SetMainColor(1);
  tlac[1]->SetLineWidth(1);
  el->AddElement(tlac[1]);

  el->SetVizTag("ESD Tracks by anal cut");
  gEve->AddElement(el);

  TEveElementList* el = new TEveElementList("ESD Tracklets SPD");
  TEveTrackList *tlac[2];
  tlac[0] = new TEveTrackList("Good");
  tlac[0]->SetLineStyle(6);
  tlac[0]->SetMainColor(1);
  tlac[0]->SetLineWidth(1);
  el->AddElement(tlac[0]);

  tlac[1] = new TEveTrackList("Bad");
  tlac[1]->SetLineStyle(6);
  tlac[1]->SetMainColor(1);
  tlac[1]->SetLineWidth(1);
  el->AddElement(tlac[1]);

  el->SetVizTag("ESD Tracklets SPD");
  gEve->AddElement(el);

}
