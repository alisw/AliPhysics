// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <AliQAHistViewer.h>

#include <TString.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TInterpreter.h>
#include <TMath.h>
#include <TGListTree.h>
#include <TEveVSDStructs.h>
#include <TEveManager.h>
#include <TEveTrackPropagator.h>
#include <TEnv.h>
#include <TEveWindowManager.h>
#include <TGTab.h>
#include <TTimeStamp.h>

#include <AliESDtrackCuts.h>
//#include <AliPWG0Helper.h>
#include <AliESDEvent.h>
#include <AliESDfriend.h>
#include <AliESDtrack.h>
#include <AliESDfriendTrack.h>
#include <AliExternalTrackParam.h>
#include <AliEveTrack.h>
#include <AliEveTrackCounter.h>
#include <AliEveMagField.h>
#include <AliEveEventManager.h>
#include <AliEveEventManagerEditor.h>
#include <AliEveMultiView.h>
#include <AliEveMacroExecutor.h>
#include <AliEveMacro.h>
#include <AliEveMacroExecutorWindow.h>
#include <AliEveEventSelectorWindow.h>
#include <AliCDBManager.h>

#include <AliEveDataSourceOffline.h>

#endif

class AliEveMacroExecutor;
class TEveProjectionManager;
class TEveGeoShape;
class TEveUtil;
class TSystem;
class TInterpreter;


Bool_t gShowMuonRPhi = kFALSE;
Bool_t gShowMuonRhoZ = kTRUE;

Bool_t gCenterProjectionsAtPrimaryVertex = kFALSE;


void visscan_init(const TString& cdburi = "",
                  const TString& path   = ".",
                  Bool_t showHLTESDTree=kFALSE,
                  Bool_t showMuon = kTRUE,
                  Bool_t showTrd = kFALSE)
{
    AliEveEventManager *man = new AliEveEventManager(AliEveEventManager::kSourceOffline);
    
    if (showMuon)
    {
        if (gSystem->Getenv("ALICE_ROOT") != 0)
        {
            gInterpreter->AddIncludePath(Form("%s/MUON", gSystem->Getenv("ALICE_ROOT")));
            gInterpreter->AddIncludePath(Form("%s/MUON/mapping", gSystem->Getenv("ALICE_ROOT")));
        }
    }
    else
    {
        gShowMuonRPhi = gShowMuonRhoZ = kFALSE;
    }
    
    if(cdburi.EqualTo(""))
    {
        cdburi = Form("%s/../src/OCDB/", gSystem->Getenv("ALICE_ROOT"));
        cout<<"\n\nsetting cdb uri:"<<cdburi<<endl;
    }
    if (cdburi.IsNull() && ! AliCDBManager::Instance()->IsDefaultStorageSet())
    {
        gEnv->SetValue("Root.Stacktrace", "no");
        Fatal("visscan_init.C", "OCDB path MUST be specified as the first argument.");
    }
    
    AliEveDataSourceOffline *dataSource = (AliEveDataSourceOffline*)man->GetDataSourceOffline();
    
    gInterpreter->AddIncludePath(Form("%s/../src/MONITOR/MONITOR", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(Form("%s/../src/STEER/", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(Form("%s/../src/", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(Form("%s/../src/ANALYSIS/ANALYSISalice", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(Form("%s/../src/STEER/STEERbase", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(Form("%s/../src/STEER/ESD", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(Form("%s/../src/EVE/EveBase", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(Form("%s/../src/STEER/STEER", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(Form("%s/../src/STEER/AOD", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(Form("%s/../src/RAW/RAWdatarec", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(Form("%s/../src/RAW/RAWDatabase", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(Form("%s/../src/STEER/CDB", gSystem->Getenv("ALICE_ROOT")));
    
    TEveUtil::LoadMacro("alieve_init.C+");
    alieve_init(cdburi, path, -1, showHLTESDTree);
    
    man->Open();
    
    // TEveLine::SetDefaultSmooth(1);
    
    TEveUtil::AssertMacro("VizDB_scan.C");
    
    AliEveMacroExecutor *exec    = man->GetExecutor();
    TEveBrowser         *browser = gEve->GetBrowser();
    browser->ShowCloseTab(kFALSE);
    
    
    //==============================================================================
    // Geometry, scenes, projections and viewers
    //==============================================================================
    
    AliEveMultiView *mv = new AliEveMultiView;
    
    mv->SetDepth(-10);
    
    TEveUtil::LoadMacro("geom_gentle.C");
    mv->InitGeomGentle(geom_gentle(), geom_gentle_rphi(), geom_gentle_rhoz(), 0);
    
    if (showTrd) {
        TEveUtil::LoadMacro("geom_gentle_trd.C+");
        mv->InitGeomGentleTrd(geom_gentle_trd());
    }
    
    if (gShowMuonRPhi || gShowMuonRhoZ) {
        TEveUtil::LoadMacro("geom_gentle_muon.C+");
        mv->InitGeomGentleMuon(geom_gentle_muon(kFALSE), gShowMuonRPhi, gShowMuonRhoZ, kFALSE);
    }
    
    mv->SetDepth(0);
    
    //==============================================================================
    // Registration of per-event macros
    //==============================================================================
    
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Track",   "kine_tracks.C", "kine_tracks", "", kTRUE));
    
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hits ITS", "its_hits.C",    "its_hits",    "", kTRUE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hits TPC", "tpc_hits.C",    "tpc_hits",    "", kTRUE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hits T0",  "t0_hits.C",     "t0_hits",     "", kTRUE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hits FMD", "fmd_hits.C",    "fmd_hits",    "", kTRUE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hits ACORDE", "acorde_hits.C",    "acorde_hits",    "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hits EMCAL", "emcal_hits.C",    "emcal_hits",    "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hits TOF",  "tof_hits.C",     "tof_hits",     "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hits TRD", "trd_hits.C",    "trd_hits",    "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hits VZERO", "vzero_hits.C",    "vzero_hits",    "", kFALSE));
    
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "DIG ITS",     "its_digits.C",  "its_digits",  "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "DIG TPC",     "tpc_digits.C",  "tpc_digits",  "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "DIG TOF",     "tof_digits.C",  "tof_digits",  "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "DIG HMPID",   "hmpid_digits.C","hmpid_digits","", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "DIG FMD",     "fmd_digits.C",  "fmd_digits",  "", kFALSE));
    
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW ITS",     "its_raw.C",     "its_raw",     "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW TPC",     "tpc_raw.C",     "tpc_raw",     "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW TOF",     "tof_raw.C",     "tof_raw",     "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW HMPID",   "hmpid_raw.C",   "hmpid_raw",   "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW T0",      "t0_raw.C",      "t0_raw",      "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW FMD",     "fmd_raw.C",     "fmd_raw",     "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW VZERO",   "vzero_raw.C",   "vzero_raw",   "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW ACORDE",  "acorde_raw.C",  "acorde_raw",  "", kFALSE));
    /*
     exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX",             "primary_vertex.C", "primary_vertex",             "",                kTRUE));
     exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Ellipse",     "primary_vertex.C", "primary_vertex_ellipse",     "",                kTRUE));
     exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Box",         "primary_vertex.C", "primary_vertex_box",         "kFALSE, 3, 3, 3", kFALSE));
     exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX SPD",         "primary_vertex.C", "primary_vertex_spd",         "",                kTRUE));
     exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Ellipse SPD", "primary_vertex.C", "primary_vertex_ellipse_spd", "",                kTRUE));
     exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Box SPD",     "primary_vertex.C", "primary_vertex_box_spd",     "kFALSE, 3, 3, 3", kFALSE));
     exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX TPC",         "primary_vertex.C", "primary_vertex_tpc",         "",                kFALSE));
     exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Ellipse TPC", "primary_vertex.C", "primary_vertex_ellipse_tpc", "",                kFALSE));
     exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Box TPC",     "primary_vertex.C", "primary_vertex_box_tpc",     "kFALSE, 3, 3, 3", kFALSE));
     */
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC V0",   "esd_V0_points.C",       "esd_V0_points_onfly"));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC V0",   "esd_V0_points.C",       "esd_V0_points_offline"));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC V0",   "esd_V0.C",              "esd_V0"));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC CSCD", "esd_cascade_points.C",  "esd_cascade_points"));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC CSCD", "esd_cascade.C",         "esd_cascade"));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC KINK", "esd_kink_points.C",     "esd_kink_points"));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC KINK", "esd_kink.C",            "esd_kink"));
    
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Tracks",              "esd_tracks.C", "esd_tracks",              "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Tracks ITS standalone",          "esd_tracks.C", "esd_tracks_ITS_standalone",              "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Tracks ITS",          "esd_tracks.C", "esd_tracks_ITS",              "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Tracks TPC",           "esd_tracks.C", "esd_tracks_TPC",              "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Tracks MI",           "esd_tracks.C", "esd_tracks_MI",           "", kFALSE));
    
    // default appearance:
    //  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Tracks by category",  "esd_tracks.C", "esd_tracks_by_category",  "", kTRUE));
    
    // preset for cosmics:
    //  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Tracks by category",  "esd_tracks.C", "esd_tracks_by_category",  "kGreen,kGreen,kGreen,kGreen,kGreen,kGreen,kGreen,kGreen,kGreen,kFALSE", kTRUE));
    
    
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Tracks by anal cuts", "esd_tracks.C", "esd_tracks_by_anal_cuts", "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Tracks Lego", "lego.C", "lego", "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Tracks Beams Info", "beams_info.C", "beams_info", "", kFALSE));
    
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Tracklets SPD", "esd_spd_tracklets.C", "esd_spd_tracklets", "", kTRUE));
    
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC ZDC",      "esd_zdc.C", "esd_zdc", "", kFALSE));
    
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clusters",     "clusters.C",     "clusters", "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clusters ITS", "its_clusters.C", "its_clusters"));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clusters TPC", "tpc_clusters.C", "tpc_clusters"));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clusters TRD", "trd_clusters.C", "trd_clusters"));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clusters TOF", "tof_clusters.C", "tof_clusters"));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clusters HMPID", "hmpid_clusters.C", "hmpid_clusters"));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clusters PHOS", "phos_clusters.C", "phos_clusters"));
    
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clusters TPC", "vplot_tpc.C",    "vplot_tpc", "", kFALSE));
    
    exec->AddMacro(new AliEveMacro(AliEveMacro::kAOD, "ANA HF",   "aod_HF.C",   "aod_HF",   "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kAOD, "ANA Jets", "jetplane.C", "jetplane", "", kFALSE));
    
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "DUMP VZERO",   "vzero_dump.C",   "vzero_dump",   "", kFALSE));
    
    if (showMuon)
    {
        exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM TrackRef MUON", "muon_trackRefs.C", "muon_trackRefs", "kTRUE", kFALSE));
        exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW MUON", "muon_raw.C", "muon_raw", "", kFALSE));
        exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "DIG MUON", "muon_digits.C", "muon_digits", "", kFALSE));
        exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clusters MUON", "muon_clusters.C", "muon_clusters", "", kTRUE));
        exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Tracks MUON", "esd_muon_tracks.C", "esd_muon_tracks", "kTRUE,kFALSE", kTRUE));
    }
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "ESD AD", "ad_esd.C", "ad_esd", "", kTRUE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "ESD EMCal", "emcal_esdclustercells.C", "emcal_esdclustercells", "", kTRUE));
    
    
    //==============================================================================
    // Additional GUI components
    //==============================================================================
    
    // Macro / data selection
    TEveWindowSlot *slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
    slot->StartEmbedding();
    AliEveMacroExecutorWindow* exewin = new AliEveMacroExecutorWindow(exec);
    slot->StopEmbedding("DataSelection");
    exewin->PopulateMacros();
    
    // Event selection tab
    slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
    slot->StartEmbedding();
    new AliEveEventSelectorWindow(gClient->GetRoot(), 600, 400, man->GetEventSelector());
    slot->StopEmbedding("Selections");
    
    // QA viewer
    /*
     slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
     slot->StartEmbedding();
     new AliQAHistViewer(gClient->GetRoot(), 600, 400, kTRUE);
     slot->StopEmbedding("QA histograms");
     
     browser->GetTabRight()->SetTab(1);
     */
    browser->StartEmbedding(TRootBrowser::kBottom);
    new AliEveEventManagerWindow(man);
    browser->StopEmbedding("EventCtrl");
    
    slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
    TEveWindowTab *store_tab = slot->MakeTab();
    store_tab->SetElementNameTitle("WindowStore",
                                   "Undocked windows whose previous container is not known\n"
                                   "are placed here when the main-frame is closed.");
    gEve->GetWindowManager()->SetDefaultContainer(store_tab);
    
    
    //==============================================================================
    // AliEve objects - global tools
    //==============================================================================
    
    AliEveTrackCounter* g_trkcnt = new AliEveTrackCounter("Primary Counter");
    gEve->AddToListTree(g_trkcnt, kFALSE);
    
    
    //==============================================================================
    // Final stuff
    //==============================================================================
    
    // A refresh to show proper window.
    //gEve->GetViewers()->SwitchColorSet();
    browser->MoveResize(0, 0, gClient->GetDisplayWidth(),gClient->GetDisplayHeight() - 32);
    gEve->Redraw3D(kTRUE);
    gSystem->ProcessEvents();
    
    // Register command to call on each event.
    man->AddNewEventCommand("on_new_event();");
    man->GotoEvent(0);
    
    gEve->EditElement(g_trkcnt);
    gEve->Redraw3D(kTRUE);
    
    //move multiview to the front
    browser->GetTabRight()->SetTab(1);
    TGLViewer *glv1 = mv->Get3DView()->GetGLViewer();
    glv1->CurrentCamera().RotateRad(-0.4, 0.6);
    gEve->FullRedraw3D();
    gSystem->ProcessEvents();
    gEve->Redraw3D(kTRUE);
    // set autoload by default
    
    
    
    Color_t colors[9] = {kCyan,kCyan,kCyan,kRed,kRed,kRed,kGreen,kGreen,kGreen};
    
    man->SetESDcolorsByCategory(colors);
    man->SetESDwidth(2);
    man->SetESDdashNoRefit(true);
    man->SetESDdrawNoRefit(true);
    
    man->SetESDtracksByCategory(false);
    man->SetESDtracksByType(true);
    
    man->SetAutoLoad(false);
}

/******************************************************************************/

void on_new_event()
{
    Double_t x[3] = { 0, 0, 0 };
    
    if (AliEveEventManager::HasESD())
    {
        AliESDEvent* esd = AliEveEventManager::Instance()->AssertESD();
        esd->GetPrimaryVertex()->GetXYZ(x);
        
        TTimeStamp ts(esd->GetTimeStamp());
        TString win_title("Eve Main Window -- Timestamp: ");
        win_title += ts.AsString("s");
        win_title += "; Event # in ESD file: ";
        win_title += esd->GetEventNumberInFile();
        gEve->GetBrowser()->SetWindowName(win_title);
    }
    
    TEveElement* top = gEve->GetCurrentEvent();
    
    AliEveMultiView *mv = AliEveMultiView::Instance();
    
    //mv->DestroyEventRPhi();
    if (gCenterProjectionsAtPrimaryVertex)
        mv->SetCenterRPhi(x[0], x[1], x[2]);
    mv->ImportEventRPhi(top);
    
    //mv->DestroyEventRhoZ();
    if (gCenterProjectionsAtPrimaryVertex)
        mv->SetCenterRhoZ(x[0], x[1], x[2]);
    mv->ImportEventRhoZ(top);
}
