#include <TG3DLine.h>
#include <TGButton.h>
#include <TGMenu.h>
#include <TGPicture.h>
#include <TGSplitter.h>
#include <TGToolBar.h>
#include <TGMsgBox.h>

#include <TGrid.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TROOT.h>

#include <TEveManager.h>
#include <TEveSelection.h>

#include <AliEveEventManager.h>
#include <AliEveMultiView.h>
#include <AliEveMacro.h>
#include <AliEveMacroExecutor.h>
#include <AliEveTrackFitter.h>
#include <AliEveTrackCounter.h>
#include <AliEveDataSourceOffline.h>

#include "AliEveMainWindow.h"
#include "AliEveUtil.h"
#include "AliEveFileDialog.h"

AliEveMainWindow::AliEveMainWindow(const char* title, UInt_t width, UInt_t height)
    : TGMainFrame(gClient->GetRoot(), width, height),
      fMenuBar(0),
      fMenuFile(0),
      fMenuEdit(0),
      fMenuView(0),
      fMenuViewToolbars(0),
      fMenuViewSidebars(0),
      fMenuGo(0),
      fMenuTools(0),
      fMenuHelp(0),
      fToolBar(0),
      fPicturePool(0),
//      fEve(0),
      fFileDialog(0)
{
	static const TEveException kEH("AliEveMainWindow");
	Info(kEH.Data(),"Constructor called");
	
    AliEveUtil::Init();
    fPicturePool = AliEveUtil::GetPicturePool();
   
//    fEve = TEveManager::Create(kFALSE, "VVV");
    gEve->GetDefaultViewer()->SetElementName("3D View");
    gEve->GetSelection()->SetPickToSelect(TEveSelection::kPS_PableCompound);
    gEve->GetHighlight()->SetPickToSelect(TEveSelection::kPS_PableCompound);

    TString evedir(Form("%s/EVE", gSystem->Getenv("ALICE_ROOT")));
    gEve->RegisterGeometryAlias("Default", Form("%s/resources/geometry/default_geo.root", evedir.Data()));
    
    setupMenus();
    setupToolbars();

    TGHorizontalFrame* hf = new TGHorizontalFrame(this,200,200);
    
    // 3D View Frame
    TGFrame* towerViewFrame = gEve->GetDefaultViewer()->GetGUIFrame();
    towerViewFrame->MapWindow();
    towerViewFrame->ReparentWindow(hf);
    
    hf->AddFrame(towerViewFrame, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY , 3, 3, 3, 3));
    
    AddFrame(hf, new TGLayoutHints(kLHintsNormal | kLHintsExpandX | kLHintsExpandY , 3, 3, 3, 3));

    SetWindowName(title);
    Resize(width,height);
    MapSubwindows();
    MapWindow();
    Layout();

    gEve->Redraw3D(kTRUE);
    gSystem->ProcessEvents();
}

AliEveMainWindow::~AliEveMainWindow()
{
//    if(fEve) delete fEve;
}

void AliEveMainWindow::onMenuFileItem(UInt_t id)
{
    switch(id){
    case MENU_FILE_OPEN:
    {
        /*
        if(!fFileDialog) fFileDialog = new AliEveFileDialog( gClient->GetRoot(), this, kAliEveFDLocal);

        fFileDialog->setMode(kAliEveFDLocal);
        fFileDialog->MapWindow();
        
        if(fFileDialog->accepted()) {
            
            AliEveEventManager* evMan = AliEveEventManager::Instance();
            AliEveDataSourceOffline *dataSource = (AliEveDataSourceOffline*)evMan->GetDataSourceOffline();
            
            if(dataSource)
            {
             dataSource->SetESDFileName(fFileDialog->GetPathESD());
             dataSource->SetESDfriendFileName(fFileDialog->GetPathESDfriend());
             dataSource->SetAODFileName(fFileDialog->GetPathAOD());
             dataSource->AddAODfriend(fFileDialog->GetPathAODfriend());
             dataSource->SetRawFileName(fFileDialog->GetPathRaw());
             evMan->SetCdbUri(fFileDialog->GetCDBStoragePath());
             loadFiles();
            }
        }*/
        break;
    }
    case MENU_FILE_OPEN_URL:
    {
        if(!fFileDialog) fFileDialog = new AliEveFileDialog( gClient->GetRoot(), this, kAliEveFDRemote);

        fFileDialog->setMode(kAliEveFDRemote);
        fFileDialog->MapWindow();
        if(fFileDialog->accepted()) {
            AliEveEventManager* evMan = AliEveEventManager::Instance();
            AliEveDataSourceOffline *dataSource = (AliEveDataSourceOffline*)evMan->GetDataSourceOffline();
            
            if(dataSource)
            {
                 dataSource->SetFilesPath(fFileDialog->GetUrl());
                 evMan->SetCdbUri(fFileDialog->GetCDBStoragePath());
            }

    // Open event
    if (fFileDialog->GetUrl().BeginsWith("alien:"))
    {
        if (gGrid != 0)
        {
            Info("AliEveMainWindow::openFile", "TGrid already initializied. Skiping checks and initialization.");
        }
        else
        {
            Info("AliEveMainWindow::openFile", "AliEn requested - connecting.");
            if (gSystem->Getenv("GSHELL_ROOT") == 0)
            {
                Error("AliEveMainWindow::openFile", "AliEn environment not initialized. Aborting.");
                new TGMsgBox(gClient->GetRoot(), this, "AliEve", "AliEn environment not initialized. Aborting.", kMBIconStop);
                return;
            }
            if (TGrid::Connect("alien") == 0)
            {
                Error("AliEveMainWindow::openFile", "TGrid::Connect() failed. Aborting.");
                new TGMsgBox(gClient->GetRoot(), this, "AliEve", "TGrid::Connect() failed. Aborting.", kMBIconStop);
                return;
            }
        }
    }
         
         loadFiles();        
        }

        break;
    }
    default:
    {
        break;
    }
    }
}

void AliEveMainWindow::onMenuEditItem(UInt_t /*id*/)
{

}

void AliEveMainWindow::onMenuViewItem(UInt_t /*id*/)
{

}

void AliEveMainWindow::onMenuGoItem(UInt_t id)
{
    switch(id){
    case MENU_GO_NEXT_EVENT:
    {
        AliEveEventManager::Instance()->NextEvent();
        break;
    }
    case MENU_GO_PREV_EVENT:
    {
        AliEveEventManager::Instance()->PrevEvent();
        break;
    }
    default:
    {
        break;
    }
    }

    TEveElement* top = gEve->GetCurrentEvent();

    AliEveMultiView *mv = AliEveMultiView::Instance();

    mv->ImportEventRPhi(top);
    mv->ImportEventRhoZ(top);

    gEve->Redraw3D(kTRUE);
}

void AliEveMainWindow::setupMenus()
{

    fMenuBar = new TGMenuBar(this);

    // File Menu
    fMenuFile = new TGPopupMenu(gClient->GetRoot());
    fMenuFile->AddEntry("&Open...", MENU_FILE_OPEN, 0, fPicturePool->GetPicture("menu/document-open.png"));
    fMenuFile->AddEntry("&Open URL...", MENU_FILE_OPEN_URL, 0, fPicturePool->GetPicture("menu/document-open-remote.png"));
    fMenuFile->AddSeparator();
    fMenuFile->AddEntry("&Connect To Server...", MENU_FILE_OPEN_CONNECTION, 0, fPicturePool->GetPicture("menu/network-connect.png"));
    fMenuFile->AddSeparator();
    fMenuFile->AddEntry("Export View(s)...", MENU_FILE_EXPORT_VIEWS, 0, fPicturePool->GetPicture("menu/document-export.png"));
    fMenuFile->AddSeparator();
    fMenuFile->AddEntry("E&xit", MENU_FILE_EXIT, 0, fPicturePool->GetPicture("menu/application-exit.png"));
    // --

    // Edit Menu
    fMenuEdit = new TGPopupMenu(gClient->GetRoot());
    fMenuEdit->AddEntry("&Undo",   MENU_EDIT_UNDO, 0, fPicturePool->GetPicture("menu/edit-undo.png"));
    fMenuEdit->AddEntry("&Redo",   MENU_EDIT_REDO, 0, fPicturePool->GetPicture("menu/edit-redo.png"));
    fMenuEdit->AddSeparator();
    fMenuEdit->AddEntry("&Cut",   MENU_EDIT_CUT, 0, fPicturePool->GetPicture("menu/edit-cut.png"));
    fMenuEdit->AddEntry("C&opy",  MENU_EDIT_COPY, 0, fPicturePool->GetPicture("menu/edit-copy.png"));
    fMenuEdit->AddEntry("&Paste", MENU_EDIT_PASTE, 0, fPicturePool->GetPicture("menu/edit-paste.png"));
    fMenuEdit->AddEntry("&Delete",MENU_EDIT_DELETE, 0, fPicturePool->GetPicture("menu/edit-delete.png"));
    fMenuEdit->AddSeparator();
    fMenuEdit->AddEntry("P&references", MENU_EDIT_PROP, 0, fPicturePool->GetPicture("menu/document-properties.png"));
    // --

    // View Menu
    fMenuView = new TGPopupMenu(gClient->GetRoot());

    fMenuViewToolbars = new TGPopupMenu(gClient->GetRoot());
    fMenuViewToolbars->AddEntry("&Main Toolbar", MENU_VIEW_TOOLBAR_MAIN);
    fMenuViewToolbars->AddEntry("&Navigation Toolbar", MENU_VIEW_TOOLBAR_NAV);
    fMenuView->AddPopup("Toolbars", fMenuViewToolbars);

    fMenuViewSidebars = new TGPopupMenu(gClient->GetRoot());
    fMenuViewSidebars->AddEntry("Hi&story", MENU_VIEW_TOOLBAR_HIST);
    fMenuViewSidebars->AddEntry("&Properties", MENU_VIEW_TOOLBAR_NAV);
    fMenuView->AddPopup("Sidebars", fMenuViewSidebars);

    fMenuView->AddSeparator();
    fMenuView->AddEntry("&Reload", MENU_VIEW_RELOAD, 0, fPicturePool->GetPicture("menu/view-refresh.png"));
    fMenuView->AddSeparator();
    fMenuView->AddEntry("Zoom &In", MENU_VIEW_ZOOM_IN, 0, fPicturePool->GetPicture("menu/zoom-in.png"));
    fMenuView->AddEntry("Zoom &Out",MENU_VIEW_ZOOM_OUT, 0, fPicturePool->GetPicture("menu/zoom-out.png"));
    fMenuView->AddEntry("Zoom &Reset",MENU_VIEW_ZOOM_RESET, 0, fPicturePool->GetPicture("menu/zoom-original.png"));
    // --

    // Go Menu
    fMenuGo = new TGPopupMenu(gClient->GetRoot());
    fMenuGo->AddEntry("&Next Event",   MENU_GO_NEXT_EVENT, 0, fPicturePool->GetPicture("navigation/media-seek-forward.png"));
    fMenuGo->AddEntry("P&revious Event",   MENU_GO_PREV_EVENT, 0, fPicturePool->GetPicture("navigation/media-seek-backward.png"));
    fMenuGo->AddSeparator();
    fMenuGo->AddEntry("&First Event", MENU_GO_FIRST_EVENT,0, fPicturePool->GetPicture("navigation/media-skip-backward.png"));
    fMenuGo->AddEntry("&Last Event",  MENU_GO_LAST_EVENT, 0, fPicturePool->GetPicture("navigation/media-skip-forward.png"));
    fMenuGo->AddSeparator();
    fMenuGo->AddEntry("&Play", MENU_GO_PLAY, 0, fPicturePool->GetPicture("navigation/media-playback-start.png"));
    // --

    // Tools Menu
    fMenuTools = new TGPopupMenu(gClient->GetRoot());
    fMenuTools->AddEntry("&QA Histograms", MENU_TOOLS_QA);
    fMenuTools->AddEntry("&Macros", MENU_TOOLS_MACROS);
    // --

    // Help Menu
    fMenuHelp = new TGPopupMenu(gClient->GetRoot());
    fMenuHelp->AddEntry("&Contents", MENU_HELP_CONTENTS, 0, fPicturePool->GetPicture("menu/help-contents.png"));
    fMenuHelp->AddEntry("&About", MENU_HELP_ABOUT, 0, fPicturePool->GetPicture("menu/help-about.png"));
    // --

    // Add popupmenus to MenuBar
    fMenuBar->AddPopup("&File", fMenuFile, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));
    fMenuBar->AddPopup("&Edit", fMenuEdit, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));
    fMenuBar->AddPopup("&View", fMenuView, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));
    fMenuBar->AddPopup("&Go", fMenuGo, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));
    fMenuBar->AddPopup("&Tools", fMenuTools, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));
    fMenuBar->AddPopup("&Help", fMenuHelp, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));

    // MenuBar to the window
    AddFrame(fMenuBar,  new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX));

    // Menu signals
    fMenuFile->Connect("Activated(Int_t)", "AliEveMainWindow", this, "onMenuFileItem(Int_t)");
    fMenuEdit->Connect("Activated(Int_t)", "AliEveMainWindow", this, "onMenuEditItem(Int_t)");
    fMenuView->Connect("Activated(Int_t)", "AliEveMainWindow", this, "onMenuViewItem(Int_t)");
    fMenuGo->Connect("Activated(Int_t)", "AliEveMainWindow", this, "onMenuGoItem(Int_t)");
}

void AliEveMainWindow::setupToolbars()
{
    TGMenuEntry* tmpMenuEntry;

    fToolBar = new TGToolBar(this);

    tmpMenuEntry = fMenuFile->GetEntry("Open...");
    fToolBar->AddButton(this, new TGPictureButton(fToolBar, tmpMenuEntry->GetPic(), tmpMenuEntry->GetEntryId() ));
    tmpMenuEntry = fMenuFile->GetEntry("Open URL...");
    fToolBar->AddButton(this, new TGPictureButton(fToolBar, tmpMenuEntry->GetPic(), tmpMenuEntry->GetEntryId() ));
    tmpMenuEntry = fMenuFile->GetEntry("Connect To Server...");
    fToolBar->AddButton(this, new TGPictureButton(fToolBar, tmpMenuEntry->GetPic(), tmpMenuEntry->GetEntryId() ));
    tmpMenuEntry = fMenuFile->GetEntry("Export View(s)...");
    fToolBar->AddButton(this, new TGPictureButton(fToolBar, tmpMenuEntry->GetPic(), tmpMenuEntry->GetEntryId() ));

    fToolBar->AddFrame(new TGVertical3DLine(fToolBar),  new TGLayoutHints(kLHintsExpandY));

    tmpMenuEntry = fMenuView->GetEntry("Reload");
    fToolBar->AddButton(this, new TGPictureButton(fToolBar, tmpMenuEntry->GetPic(), tmpMenuEntry->GetEntryId() ));
    tmpMenuEntry = fMenuView->GetEntry("Zoom In");
    fToolBar->AddButton(this, new TGPictureButton(fToolBar, tmpMenuEntry->GetPic(), tmpMenuEntry->GetEntryId() ));
    tmpMenuEntry = fMenuView->GetEntry("Zoom Out");
    fToolBar->AddButton(this, new TGPictureButton(fToolBar, tmpMenuEntry->GetPic(), tmpMenuEntry->GetEntryId() ));
    tmpMenuEntry = fMenuView->GetEntry("Zoom Reset");
    fToolBar->AddButton(this, new TGPictureButton(fToolBar, tmpMenuEntry->GetPic(), tmpMenuEntry->GetEntryId() ));




    AddFrame(new TGHorizontal3DLine(this), new TGLayoutHints(kLHintsExpandX));
    AddFrame(fToolBar, new TGLayoutHints(kLHintsNormal));
    AddFrame(new TGHorizontal3DLine(this), new TGLayoutHints(kLHintsExpandX));


    //fToolBar->Connect("Clicked(Int_t)", "RCMainWindow", this, "openFile()");
}

void AliEveMainWindow::loadFiles()
{
    TString name("Event"); // CINT has trouble with direct "Event".
//    new AliEveEventManager(name, 0);
    gEve->AddEvent(AliEveEventManager::Instance());

    TEveUtil::AssertMacro("VizDB_scan.C");

    AliEveMacroExecutor *exec    = AliEveEventManager::Instance()->GetExecutor();
    //==============================================================================
    // Geometry, scenes, projections and viewers
    //==============================================================================

    AliEveMultiView *mv = new AliEveMultiView;

    mv->SetDepth(-10);

    TEveUtil::LoadMacro("geom_gentle.C");
    mv->InitGeomGentle(geom_gentle(), geom_gentle_rphi(), geom_gentle_rhoz(), 0);


    TEveUtil::LoadMacro("geom_gentle_trd.C");
    mv->InitGeomGentleTrd(geom_gentle_trd());

    TEveUtil::LoadMacro("geom_gentle_muon.C");
    mv->InitGeomGentleMuon(geom_gentle_muon(kFALSE), kTRUE, kTRUE, kFALSE);

    mv->SetDepth(0);

    //==============================================================================
    // Registration of per-event macros
    //==============================================================================

    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Track",   "kine_tracks.C", "kine_tracks", "", kFALSE));

    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hits ITS", "its_hits.C",    "its_hits",    "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hits TPC", "tpc_hits.C",    "tpc_hits",    "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hits T0",  "t0_hits.C",     "t0_hits",     "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hits FMD", "fmd_hits.C",    "fmd_hits",    "", kFALSE));
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
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Tracks by category",  "esd_tracks.C", "esd_tracks_by_category",  "", kTRUE));
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


    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM TrackRef MUON", "muon_trackRefs.C", "muon_trackRefs", "kTRUE", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW MUON", "muon_raw.C", "muon_raw", "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "DIG MUON", "muon_digits.C", "muon_digits", "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clusters MUON", "muon_clusters.C", "muon_clusters", "", kTRUE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Tracks MUON", "esd_muon_tracks.C", "esd_muon_tracks", "kTRUE,kFALSE", kTRUE));


    //==============================================================================
    // AliEve objects - global tools
    //==============================================================================

    AliEveTrackFitter* fitter = new AliEveTrackFitter();
    gEve->AddToListTree(fitter, 1);
    gEve->AddElement(fitter, gEve->GetEventScene());

    AliEveTrackCounter* g_trkcnt = new AliEveTrackCounter("Primary Counter");
    gEve->AddToListTree(g_trkcnt, kFALSE);

    // A refresh to show proper window.
    //gEve->GetViewers()->SwitchColorSet();
    gEve->Redraw3D(kTRUE);
    gSystem->ProcessEvents();

    // Register command to call on each event.
    // AliEveEventManager::Instance()->AddNewEventCommand("on_new_event();");
    AliEveEventManager::Instance()->GotoEvent(0);

    gEve->Redraw3D(kTRUE);
}

