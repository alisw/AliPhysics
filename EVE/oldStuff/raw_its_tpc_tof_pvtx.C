//Macro generated automatically by AliEveMacroExecutor
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TEveManager.h>
#include <TEveBrowser.h>
#include <TEveWindow.h>
#include <TEvePointSet.h>

#include <AliEveEventManager.h>
#include <AliEveMacro.h>
#include <AliEveMacroExecutor.h>
#include <AliEveMacroExecutorWindow.h>
#include <AliEveMultiView.h>
#endif

void raw_its_tpc_tof_pvtx(){

  AliEveMacroExecutor *exec = AliEveEventManager::Instance()->GetExecutor();
  exec->RemoveMacros();
  TEveBrowser *browser = gEve->GetBrowser();
  browser->ShowCloseTab(kFALSE);
  exec->AddMacro(new AliEveMacro(2, "REC PVTX", "primary_vertex.C", "primary_vertex", "", 1));

  exec->AddMacro(new AliEveMacro(2, "REC PVTX Ellipse", "primary_vertex.C", "primary_vertex_ellipse", "", 0));

  exec->AddMacro(new AliEveMacro(2, "REC PVTX Box", "primary_vertex.C", "primary_vertex_box", "kFALSE, 3, 3, 3", 0));

  exec->AddMacro(new AliEveMacro(2, "REC PVTX", "primary_vertex.C", "primary_vertex_spd", "", 0));

  exec->AddMacro(new AliEveMacro(2, "REC PVTX Ellipse", "primary_vertex.C", "primary_vertex_ellipse_spd", "", 0));

  exec->AddMacro(new AliEveMacro(2, "REC PVTX Box", "primary_vertex.C", "primary_vertex_box_spd", "kFALSE, 3, 3, 3", 0));

  exec->AddMacro(new AliEveMacro(2, "REC PVTX", "primary_vertex.C", "primary_vertex_tpc", "", 0));

  exec->AddMacro(new AliEveMacro(2, "REC PVTX Ellipse", "primary_vertex.C", "primary_vertex_ellipse_tpc", "", 0));

  exec->AddMacro(new AliEveMacro(2, "REC PVTX Box", "primary_vertex.C", "primary_vertex_box_tpc", "kFALSE, 3, 3, 3", 0));

  exec->AddMacro(new AliEveMacro(1, "REC Clus ITS", "its_clusters.C++", "its_clusters", "", 0));

  exec->AddMacro(new AliEveMacro(1, "REC Clus TPC", "tpc_clusters.C++", "tpc_clusters", "", 0));

  exec->AddMacro(new AliEveMacro(1, "REC Clus TRD", "trd_clusters.C++", "trd_clusters", "", 0));

  exec->AddMacro(new AliEveMacro(1, "REC Clus TOF", "tof_clusters.C++", "tof_clusters", "", 0));

  exec->AddMacro(new AliEveMacro(1, "REC Clus HMPID", "hmpid_clusters.C++", "hmpid_clusters", "", 0));

  exec->AddMacro(new AliEveMacro(1, "REC Clus MUON", "muon_clusters.C++", "muon_clusters", "", 0));

  exec->AddMacro(new AliEveMacro(1, "REC Clus TOF", "emcal_digits.C++", "emcal_digits", "", 0));

  exec->AddMacro(new AliEveMacro(8, "RAW ITS", "its_raw.C", "its_raw", "", 1));

  exec->AddMacro(new AliEveMacro(8, "RAW TPC", "tpc_raw.C", "tpc_raw", "", 1));

  exec->AddMacro(new AliEveMacro(8, "RAW TOF", "tof_raw.C", "tof_raw", "", 1));

  exec->AddMacro(new AliEveMacro(8, "RAW VZERO", "vzero_raw.C", "vzero_raw", "", 0));

  exec->AddMacro(new AliEveMacro(8, "RAW ACORDE", "acorde_raw.C", "acorde_raw", "", 0));

  exec->AddMacro(new AliEveMacro(8, "RAW MUON", "muon_raw.C++", "muon_raw", "", 0));

  exec->AddMacro(new AliEveMacro(8, "RAW FMD", "fmd_raw.C", "fmd_raw", "", 0));

  exec->AddMacro(new AliEveMacro(2, "REC Track", "esd_tracks.C", "esd_tracks", "", 0));

  exec->AddMacro(new AliEveMacro(2, "REC Track", "esd_tracks.C", "esd_tracks_MI", "", 0));

  exec->AddMacro(new AliEveMacro(2, "REC Track", "esd_tracks.C", "esd_tracks_by_category", "", 0));

  exec->AddMacro(new AliEveMacro(2, "REC Track MUON", "esd_muon_tracks.C++", "esd_muon_tracks", "kTRUE,kFALSE", 0));

  exec->AddMacro(new AliEveMacro(2, "REC FMD", "fmd_esd.C", "fmd_esd", "", 0));

  exec->AddMacro(new AliEveMacro(2, "REC TRD", "trd_detectors.C++", "trd_detectors", "", 0));

  TEveWindowSlot *slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  slot->StartEmbedding();
  AliEveMacroExecutorWindow* exewin = new AliEveMacroExecutorWindow(exec);
  slot->StopEmbedding("DataSelection");
  exewin->PopulateMacros();


}
