void DataSelection_init(){

  AliEveMacroExecutor *exec = AliEveEventManager::Instance()->GetExecutor();
  exec->RemoveMacros();
  TEveBrowser *browser = gEve->GetBrowser();
  browser->ShowCloseTab(kFALSE);

  exec->AddMacro(new AliEveMacro(1, "SIM Track", "kine_tracks.C+", "kine_tracks", "", 0));

  exec->AddMacro(new AliEveMacro(1, "SIM Hits ITS", "its_hits.C", "its_hits", "", 0));

  exec->AddMacro(new AliEveMacro(1, "SIM Hits TPC", "tpc_hits.C", "tpc_hits", "", 0));

  exec->AddMacro(new AliEveMacro(1, "SIM Hits T0", "t0_hits.C", "t0_hits", "", 0));

  exec->AddMacro(new AliEveMacro(1, "SIM Hits FMD", "fmd_hits.C", "fmd_hits", "", 0));

  exec->AddMacro(new AliEveMacro(1, "DIG ITS", "its_digits.C", "its_digits", "", 0));

  exec->AddMacro(new AliEveMacro(1, "DIG TPC", "tpc_digits.C", "tpc_digits", "", 0));

  exec->AddMacro(new AliEveMacro(1, "DIG TOF", "tof_digits.C", "tof_digits", "", 0));

  exec->AddMacro(new AliEveMacro(1, "DIG HMPID", "hmpid_digits.C", "hmpid_digits", "", 0));

  exec->AddMacro(new AliEveMacro(1, "DIG FMD", "fmd_digits.C", "fmd_digits", "", 0));

  exec->AddMacro(new AliEveMacro(8, "RAW ITS", "its_raw.C", "its_raw", "", 0));

  exec->AddMacro(new AliEveMacro(8, "RAW TPC", "tpc_raw.C", "tpc_raw", "", 0));

  exec->AddMacro(new AliEveMacro(8, "RAW TOF", "tof_raw.C", "tof_raw", "", 0));

  exec->AddMacro(new AliEveMacro(8, "RAW HMPID", "hmpid_raw.C", "hmpid_raw", "", 0));

  exec->AddMacro(new AliEveMacro(8, "RAW T0", "t0_raw.C", "t0_raw", "", 0));

  exec->AddMacro(new AliEveMacro(8, "RAW FMD", "fmd_raw.C", "fmd_raw", "", 0));

  exec->AddMacro(new AliEveMacro(8, "RAW VZERO", "vzero_raw.C", "vzero_raw", "", 0));

  exec->AddMacro(new AliEveMacro(8, "RAW ACORDE", "acorde_raw.C", "acorde_raw", "", 0));

  exec->AddMacro(new AliEveMacro(2, "REC PVTX", "primary_vertex.C+", "primary_vertex", "", 1));

  exec->AddMacro(new AliEveMacro(2, "REC PVTX Ellipse", "primary_vertex.C+", "primary_vertex_ellipse", "", 1));

  exec->AddMacro(new AliEveMacro(2, "REC PVTX Box", "primary_vertex.C+", "primary_vertex_box", "kFALSE, 3, 3, 3", 0));

  exec->AddMacro(new AliEveMacro(2, "REC PVTX SPD", "primary_vertex.C+", "primary_vertex_spd", "", 1));

  exec->AddMacro(new AliEveMacro(2, "REC PVTX Ellipse SPD", "primary_vertex.C+", "primary_vertex_ellipse_spd", "", 1));

  exec->AddMacro(new AliEveMacro(2, "REC PVTX Box SPD", "primary_vertex.C+", "primary_vertex_box_spd", "kFALSE, 3, 3, 3", 0));

  exec->AddMacro(new AliEveMacro(2, "REC PVTX TPC", "primary_vertex.C+", "primary_vertex_tpc", "", 0));

  exec->AddMacro(new AliEveMacro(2, "REC PVTX Ellipse TPC", "primary_vertex.C+", "primary_vertex_ellipse_tpc", "", 0));

  exec->AddMacro(new AliEveMacro(2, "REC PVTX Box TPC", "primary_vertex.C+", "primary_vertex_box_tpc", "kFALSE, 3, 3, 3", 0));

  exec->AddMacro(new AliEveMacro(2, "REC V0", "esd_V0_points.C", "esd_V0_points_onfly", "", 1));

  exec->AddMacro(new AliEveMacro(2, "REC V0", "esd_V0_points.C", "esd_V0_points_offline", "", 1));

  exec->AddMacro(new AliEveMacro(2, "REC V0", "esd_V0.C", "esd_V0", "", 1));

  exec->AddMacro(new AliEveMacro(2, "REC CSCD", "esd_cascade_points.C", "esd_cascade_points", "", 1));

  exec->AddMacro(new AliEveMacro(2, "REC CSCD", "esd_cascade.C", "esd_cascade", "", 1));

  exec->AddMacro(new AliEveMacro(2, "REC KINK", "esd_kink_points.C", "esd_kink_points", "", 1));

  exec->AddMacro(new AliEveMacro(2, "REC KINK", "esd_kink.C", "esd_kink", "", 1));

  exec->AddMacro(new AliEveMacro(2, "REC Tracks", "esd_tracks.C+", "esd_tracks", "", 0));

  exec->AddMacro(new AliEveMacro(2, "REC Tracks MI", "esd_tracks.C+", "esd_tracks_MI", "", 0));

  exec->AddMacro(new AliEveMacro(2, "REC Tracks by category", "esd_tracks.C+", "esd_tracks_by_category", "", 1));

  exec->AddMacro(new AliEveMacro(2, "REC Tracks by anal cuts", "esd_tracks.C+", "esd_tracks_by_anal_cuts", "", 0));

  exec->AddMacro(new AliEveMacro(2, "REC Tracklets SPD", "esd_spd_tracklets.C", "esd_spd_tracklets", "", 1));

  exec->AddMacro(new AliEveMacro(2, "REC ZDC", "esd_zdc.C", "esd_zdc", "", 0));

  exec->AddMacro(new AliEveMacro(1, "REC Clusters", "clusters.C+", "clusters", "", 0));

  exec->AddMacro(new AliEveMacro(1, "REC Clusters ITS", "its_clusters.C+", "its_clusters", "", 1));

  exec->AddMacro(new AliEveMacro(1, "REC Clusters TPC", "tpc_clusters.C+", "tpc_clusters", "", 1));

  exec->AddMacro(new AliEveMacro(1, "REC Clusters TRD", "trd_clusters.C+", "trd_clusters", "", 1));

  exec->AddMacro(new AliEveMacro(1, "REC Clusters TOF", "tof_clusters.C+", "tof_clusters", "", 1));

  exec->AddMacro(new AliEveMacro(1, "REC Clusters TPC", "vplot_tpc.C+", "vplot_tpc", "", 0));

  exec->AddMacro(new AliEveMacro(16, "ANA HF", "aod_HF.C", "aod_HF", "", 0));

  exec->AddMacro(new AliEveMacro(16, "ANA Jets", "jetplane.C", "jetplane", "", 0));

  exec->AddMacro(new AliEveMacro(2, "DUMP VZERO", "vzero_dump.C", "vzero_dump", "", 0));

  TEveWindowSlot *slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  slot->StartEmbedding();
  AliEveMacroExecutorWindow* exewin = new AliEveMacroExecutorWindow(exec);
  slot->StopEmbedding("DataSelection");
  exewin->PopulateMacros();

}
