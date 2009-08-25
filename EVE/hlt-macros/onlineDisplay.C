//-*- Mode: C++ -*-

// ** USED macros :
// ***************************************************
// - hlt_alieve_init.C
// - VizDB_scan.C
// - geom_gentle_hlt.C
// - geom_gentle_muon.C
// ***************************************************

#include "unistd.h"
#include <TEvePointSet.h>
#include "EveBase/AliEveEventManager.h"
#include <AliCluster.h>
#include <TPC/AliTPCClustersRow.h>

class TEveProjectionManager;
class TEveGeoShape;
class TEveUtil;
class AliEveMacroExecutor;
class TEveScene;
class TEveElement;
class TEveText;
class AliHLTTriggerDecision;
class TEvePointSet;

// -----------------------------------------------------------------
// --                       Geometry / Scenes                     --
// -----------------------------------------------------------------

TEveProjectionManager *gRPhiMgr = 0;
TEveProjectionManager *gRhoZMgr = 0;

// -----------------------------------------------------------------
// --                Geometry / Scenes Parameters                 --
// -----------------------------------------------------------------

// -- Parameters to show different geometries
Bool_t gShowMUON     = kTRUE;
Bool_t gShowMUONRPhi = kFALSE;
Bool_t gShowMUONRhoZ = kTRUE;

// -----------------------------------------------------------------
// --                         Members                            --
// -----------------------------------------------------------------

// -- Timer for automatic event loop
TTimer                                    eventTimer;

// -- HOMERManager
AliEveHOMERManager*                       gHomerManager    = 0;

// -- Cluster members
TEvePointSet*                             gPHOSClusters    = 0;
TEvePointSet*                             gTPCClusters     = 0;
TEvePointSet*                             gSPDClusters     = 0;

// -- Tracks members
TEveTrackList*                            gTPCTrack        = 0;

// -----------------------------------------------------------------
// --                          Methods                            --
// -----------------------------------------------------------------

Int_t initializeEveViewer( Bool_t showExtraGeo );

Int_t nextEvent();

Int_t processPHOSClusters( AliHLTHOMERBlockDesc* block);

Int_t processEsdTracks( AliHLTHOMERBlockDesc* block, TEveTrackList* cont );

Int_t processHLTRDLST( AliHLTHOMERBlockDesc* block );

Int_t processROOTTOBJ( AliHLTHOMERBlockDesc* block );

Int_t processTPCClusters (AliHLTHOMERBlockDesc * block, TEvePointSet cont );

// #################################################################
// #################################################################
// #################################################################

// -----------------------------------------------------------------
void onlineDisplay(Bool_t showMuonGeo=kFALSE) {
  
  // -- Loading Geometry
  // ---------------------
  Int_t run = 67179;
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetRun(run);
  AliGeomManager::LoadGeometry();

  // -- Create new hM object
  // -------------------------
  gHomerManager = new AliEveHOMERManager();
  gHomerManager->SetRetryCount(50,5);

  Int_t iResult = gHomerManager->Initialize();
  if (iResult) { 
    printf("Error Initializing AliHLTHOMERManager, quitting");
    return; 
  }

  // -- Add hM to EveTree
  // ----------------------
  gEve->AddToListTree(gHomerManager, kTRUE);

  // -- Create SourceList
  // ----------------------
  iResult = gHomerManager->CreateEveSourcesListLoop();
  if (iResult) {
    printf ("Couldn't find active services. returning\n");
    return;
  } 

  // -- Initialize pointsets and add macros
  // ----------------------------------------
  TEveUtil::LoadMacro("hlt_alieve_init.C");
  hlt_alieve_init(".", -1);

  // -- Initialize Eve
  // -------------------
  initializeEveViewer( showMuonGeo );

  // -- Finalize Eve
  // -----------------
  gSystem->ProcessEvents();
  gEve->Redraw3D(kTRUE);
}

// -----------------------------------------------------------------
Int_t initializeEveViewer( Bool_t showMuonGeo ) {
  
  //==============================================================================
  // Geometry, scenes, projections and viewers
  //==============================================================================

  TEveGeoShape *geomGentle     = 0;
  TEveGeoShape *geomGentleRPhi = 0;
  TEveGeoShape *geomGentleRhoZ = 0;
  TEveGeoShape *geomGentleTRD  = 0;
  TEveGeoShape *geomGentleMUON = 0;

  // -- Disable extra geometry
  // ---------------------------
  if ( ! showMuonGeo ) {
    gShowMUON = gShowMUONRPhi = gShowMUONRhoZ = kFALSE;
  }
  
  // -- Load Geometry
  // ------------------
  TEveUtil::LoadMacro("geom_gentle_hlt.C");
  geomGentle = geom_gentle_hlt();
  geomGentleRPhi = geom_gentle_rphi(); geomGentleRPhi->IncDenyDestroy();
  geomGentleRhoZ = geom_gentle_rhoz(); geomGentleRhoZ->IncDenyDestroy();
  geomGentleTRD  = geom_gentle_trd();
  
  if (gShowMUON) {
    TEveUtil::LoadMacro("geom_gentle_muon.C");
    geomGentleMUON = geom_gentle_muon(kFALSE);
  }

  // -- Scenes
  // -----------

  TEveScene *rPhiGeomScene  = gEve->SpawnNewScene("RPhi Geometry",
		                "Scene holding projected geometry for the RPhi view.");
  TEveScene *rhoZGeomScene  = gEve->SpawnNewScene("RhoZ Geometry",
				"Scene holding projected geometry for the RhoZ view.");
  TEveScene *rPhiEventScene = gEve->SpawnNewScene("RPhi Event Data",
		                "Scene holding projected geometry for the RPhi view.");
  TEveScene *rhoZEventScene = gEve->SpawnNewScene("RhoZ Event Data",
				"Scene holding projected geometry for the RhoZ view.");

  // -- Projection managers
  // ------------------------

  // -- R-Phi Projection
  gRPhiMgr = new TEveProjectionManager();
  gRPhiMgr->SetProjection(TEveProjection::kPT_RPhi);
  gEve->AddToListTree(gRPhiMgr, kFALSE);
  {
    TEveProjectionAxes* a = new TEveProjectionAxes(gRPhiMgr);
    a->SetMainColor(kWhite);
    a->SetTitle("R-Phi");
    a->SetTitleSize(0.05);
    a->SetTitleFont(102);
    a->SetLabelSize(0.025);
    a->SetLabelFont(102);
    rPhiGeomScene->AddElement(a);
  }

  gRPhiMgr->SetCurrentDepth(-10);
  gRPhiMgr->ImportElements(geomGentleRPhi, rPhiGeomScene);
  gRPhiMgr->SetCurrentDepth(0);
  gRPhiMgr->ImportElements(geomGentleTRD, rPhiGeomScene);
  if (gShowMUONRPhi) gRPhiMgr->ImportElements(geomGentleMUON, rPhiGeomScene);

  // -- Rho-Z Projection
  gRhoZMgr = new TEveProjectionManager();
  gRhoZMgr->SetProjection(TEveProjection::kPT_RhoZ);
  gEve->AddToListTree(gRhoZMgr, kFALSE);
  {
    TEveProjectionAxes* a = new TEveProjectionAxes(gRhoZMgr);
    a->SetMainColor(kWhite);
    a->SetTitle("Rho-Z");
    a->SetTitleSize(0.05);
    a->SetTitleFont(102);
    a->SetLabelSize(0.025);
    a->SetLabelFont(102);
    rhoZGeomScene->AddElement(a);
  }
  gRhoZMgr->SetCurrentDepth(-10);
  gRhoZMgr->ImportElements(geomGentleRhoZ, rhoZGeomScene);
  gRhoZMgr->SetCurrentDepth(0);
  gRhoZMgr->ImportElements(geomGentleTRD, rhoZGeomScene);
  if (gShowMUONRhoZ) gRhoZMgr->ImportElements(geomGentleMUON, rhoZGeomScene);

  // -- Viewers
  // ------------
  TEveBrowser *browser = gEve->GetBrowser();
  browser->ShowCloseTab(kFALSE);

  TEveViewer *threeDView  = 0;
  TEveViewer *rPhiView = 0;
  TEveViewer *rhoZView = 0;
  
  TEveWindowSlot *slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  TEveWindowPack *pack = slot->MakePack();
  pack->SetElementName("Multi View");
  pack->SetHorizontal();
  pack->SetShowTitleBar(kFALSE);
  pack->NewSlot()->MakeCurrent();
  threeDView = gEve->SpawnNewViewer("3D View", "");
  threeDView->AddScene(gEve->GetGlobalScene());
  threeDView->AddScene(gEve->GetEventScene());
    
  pack = pack->NewSlot()->MakePack();
  pack->SetShowTitleBar(kFALSE);
  pack->NewSlot()->MakeCurrent();
  rPhiView = gEve->SpawnNewViewer("RPhi View", "");
  rPhiView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  rPhiView->AddScene(rPhiGeomScene);
  rPhiView->AddScene(rPhiEventScene);

  pack->NewSlot()->MakeCurrent();
  rhoZView = gEve->SpawnNewViewer("RhoZ View", "");
  rhoZView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  rhoZView->AddScene(rhoZGeomScene);
  rhoZView->AddScene(rhoZEventScene);

  TEveViewerList *viewerlist = new TEveViewerList();
  viewerlist->AddElement(gEve->GetDefaultViewer());

  viewerlist->AddElement(threeDView);
  viewerlist->AddElement(rhoZView);
  viewerlist->AddElement(rPhiView);
  viewerlist->AddElement(threeDView);
  viewerlist->SwitchColorSet();

  //==============================================================================
  // Macros / QA histograms
  //==============================================================================

  // -- Registration of per-event macros
  // -------------------------------------

  AliEveMacroExecutor *exec    = new AliEveMacroExecutor();
#if 0
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Track",   "kine_tracks.C", "kine_tracks", "", kFALSE));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hit ITS", "its_hits.C",    "its_hits",    "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hit TPC", "tpc_hits.C",    "tpc_hits",    "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hit T0",  "t0_hits.C",     "t0_hits",     "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "SIM Hit FMD", "fmd_hits.C",    "fmd_hits",    "", kFALSE));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "DIG FMD",     "fmd_digits.C",  "fmd_digits",  "", kFALSE));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW TPC",     "tpc_raw.C",     "tpc_raw",     "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW T0",      "t0_raw.C",      "t0_raw",      "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW FMD",     "fmd_raw.C",     "fmd_raw",     "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW VZERO",   "vzero_raw.C",   "vzero_raw",   "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW ACORDE",  "acorde_raw.C",  "acorde_raw",  "", kFALSE));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX",         "primary_vertex.C", "primary_vertex",             "",                kTRUE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Ellipse", "primary_vertex.C", "primary_vertex_ellipse",     "",                kTRUE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Box",     "primary_vertex.C", "primary_vertex_box",         "kFALSE, 3, 3, 3", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX",         "primary_vertex.C", "primary_vertex_spd",         "",                kTRUE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Ellipse", "primary_vertex.C", "primary_vertex_ellipse_spd", "",                kTRUE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Box",     "primary_vertex.C", "primary_vertex_box_spd",     "kFALSE, 3, 3, 3", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX",         "primary_vertex.C", "primary_vertex_tpc",         "",                kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Ellipse", "primary_vertex.C", "primary_vertex_ellipse_tpc", "",                kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Box",     "primary_vertex.C", "primary_vertex_box_tpc",     "kFALSE, 3, 3, 3", kFALSE));
#endif
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC V0",   "esd_V0_points.C",       "esd_V0_points_onfly"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC V0",   "esd_V0_points.C",       "esd_V0_points_offline"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC V0",   "esd_V0.C",              "esd_V0"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC CSCD", "esd_cascade_points.C",  "esd_cascade_points"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC CSCD", "esd_cascade.C",         "esd_cascade"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC KINK", "esd_kink_points.C",     "esd_kink_points"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC KINK", "esd_kink.C",            "esd_kink"));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track", "esd_tracks.C", "esd_tracks",             "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track", "esd_tracks.C", "esd_tracks_MI",          "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track", "esd_tracks.C", "esd_tracks_by_category", "", kTRUE));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Tracklet", "esd_spd_tracklets.C", "esd_spd_tracklets", "", kFALSE));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC ZDC",      "esd_zdc.C", "esd_zdc", "", kFALSE));
#if 0
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus",     "clusters.C+",     "clusters", "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus ITS", "its_clusters.C+", "its_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus TPC", "tpc_clusters.C+", "tpc_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus TRD", "trd_clusters.C+", "trd_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus TOF", "tof_clusters.C+", "tof_clusters"));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus TPC", "vplot_tpc.C+",    "vplot_tpc", "", kFALSE));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kAOD, "ANA HF",   "aod_HF.C",   "aod_HF",   "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kAOD, "ANA Jets", "jetplane.C", "jetplane", "", kFALSE));

  // -- QA Viewer
  // --------------

  slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  slot->StartEmbedding();
  new AliQAHistViewer(gClient->GetRoot(), 600, 400, kTRUE);
  slot->StopEmbedding("QA histograms");

#endif

  //==============================================================================
  // Additional GUI components
  //==============================================================================
  
  slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  TEveWindowTab *store_tab = slot->MakeTab();
  store_tab->SetElementNameTitle("WindowStore",
   				 "Undocked windows whose previous container is not known\n"
 				 "are placed here when the main-frame is closed.");
  gEve->GetWindowManager()->SetDefaultContainer(store_tab);
  
  return 0;
}

// -----------------------------------------------------------------
Int_t nextEvent() {

  Int_t iResult = 0;

  gStyle->SetPalette(1, 0);
  gEve->DisableRedraw();

  // -- Get Next Event from HOMER
  // ------------------------------
  if ( ( iResult = gHomerManager->NextEvent()) ){
    return iResult;
  }

  // -- Reset
  // ----------
  if ( gTPCClusters ) gTPCClusters->Reset();
  if ( gPHOSClusters ) gPHOSClusters->Reset();
  if ( gTPCTrack )    gTPCTrack->DestroyElements();

  if (gHomerManager->GetBlockList()->IsEmpty() ) {
    printf ("No Blocks in list ... ");
    return;
  }

  TIter next(gHomerManager->GetBlockList());
  AliHLTHOMERBlockDesc* block = 0;

  // -- Iterate over blocks in the block list
  // ------------------------------------------
  while ((block = (AliHLTHOMERBlockDesc*)next())) {
    
    printf( "------------------- xxxxxxxxxxxxxxx ----------------------\n");
    printf( "Detector           : %s\n", block->GetDetector().Data() );
    printf( "Datatype           : %s\n", block->GetDataType().Data() );
    if (block->IsTObject() )
      printf( "Is TObject of class: %s\n", block->GetClassName().Data() );
    printf( "------------------- xxxxxxxxxxxxxxx ----------------------\n");

    // -- CHECK SOURCE
    // -----------------------------------------------------
    
    // ++ HLT BLOCK
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if ( ! block->GetDetector().CompareTo("HLT") ) {

      // -- ESDTREE
      if ( ! block->GetDataType().CompareTo("ALIESDV0") ) {
	cout << "ALIESDV0 ------- ALIESDV0 ------ ALIESDV0" << endl;

	if(!gTPCTrack){
	  gTPCTrack = new TEveTrackList("ESD Tracks");
	  gTPCTrack->SetMainColor(6);
	  gEve->AddElement(gTPCTrack);
	}

	iResult = processEsdTracks(block, gTPCTrack);
      }

      // -- Process ROOTObj
      else if ( ! block->GetDataType().CompareTo("ROOTTOBJ") ) {
       	processROOTTOBJ( block );
      } 

      // -- Process HLT RDLST
      else if ( ! block->GetDataType().CompareTo("HLTRDLST") ) {
 	processHLTRDLST( block );
      } 
    } // if ( ! block->GetDetector().CompareTo("HLT") ) {

    // ++ TPC BLOCK
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    else if ( ! block->GetDetector().CompareTo("TPC") ) {
      
      // -- ESDTREE
      if ( ! block->GetDataType().CompareTo("ALIESDV0") ) {

      }
    } // else if ( ! block->GetDetector().CompareTo("HLT") ) {

    // ++ ITS - SPD  BLOCK
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    else if ( ! block->GetDetector().CompareTo("ISPD") ){
      if ( block->GetDataType().CompareTo("CLUSTERS") == 0 ) {
 	//processISPDClusters( block );
      } 
    } // else if ( ! block->GetDetector().CompareTo("ISPD") ){

  } // while ((block = (AliHLTHOMERBlockDesc*)next())) {

  if ( gTPCClusters ) gTPCClusters->ResetBBox();
  if ( gPHOSClusters ) gPHOSClusters->ResetBBox();
  if ( gSPDClusters ) gSPDClusters->ResetBBox();
  if ( gTPCTrack ) gTPCTrack->ElementChanged();

  
#if 0
  TTimeStamp ts(esd->GetTimeStamp());
  TString win_title("Eve Main Window -- Timestamp: ");
  win_title += ts.AsString("s");
  win_title += "; Event # in ESD file: ";
  win_title += esd->GetEventNumberInFile();
  gEve->GetBrowser()->SetWindowName(win_title);
#endif

  // -- Set Projections
  // --------------------------------------------
  TEveElement* top = gEve->GetCurrentEvent();

  // XXX Primary vertex ... to be retrieved from the ESD
  Double_t x[3] = { 0, 0, 0 };

  if (gRPhiMgr && top) {
    gRPhiMgr->DestroyElements();
    gRPhiMgr->SetCenter(x[0], x[1], x[2]);
    gRPhiMgr->ImportElements(geomGentleRPhi);
    gRPhiMgr->ImportElements(top);
  }
  if (gRhoZMgr && top) {
    gRhoZMgr->DestroyElements();
    gRhoZMgr->SetCenter(x[0], x[1], x[2]);
    gRhoZMgr->ImportElements(geomGentleRhoZ);
    gRhoZMgr->ImportElements(geomMuon);
    gRhoZMgr->ImportElements(top);
  }
  // --------------------------------------------

  gEve->Redraw3D(0,1); // (0, 1)
  gEve->EnableRedraw(); 
  
  return iResult;
}

// -----------------------------------------------------------------
void loopEvent() {
  eventTimer.SetCommand("nextEvent()");
  eventTimer.Start(6000);
}

// -----------------------------------------------------------------
void stopLoopEvent() {
  eventTimer.Stop();
}

// -----------------------------------------------------------------
Int_t processHLTRDLST(AliHLTHOMERBlockDesc* block) {

  return 0;
}

// -----------------------------------------------------------------
Int_t processISPDClusters(AliHLTHOMERBlockDesc* block) {
  cout<<"ISPD dump:"<<endl;
  TObject ob = block->GetTObject();
  ob.Dump();

  return 0;
}

// -----------------------------------------------------------------
Int_t processROOTTOBJ(AliHLTHOMERBlockDesc* block) {
  
  // -- AliHLTGlobalTriggerDecision
  if ( ! block->GetClassName().CompareTo("AliHLTGlobalTriggerDecision") ) {

    AliHLTGlobalTriggerDecision *trig = dynamic_cast<AliHLTGlobalTriggerDecision*> block->GetTObject();
    trig->Print(); 

    TEveText* tt = new TEveText("Trigger: Class is known ;-) ");
    gEve->AddElement(tt);

  }
  else {
    printf(" Unknown root object %s",block->GetClassName().Data() );
  }

  return 0;
}

// -----------------------------------------------------------------
Int_t processEsdTracks( AliHLTHOMERBlockDesc* block, TEveTrackList* cont ) {

  AliESDEvent* esd = (AliESDEvent *) (block->GetTObject());
  esd->GetStdContent();

  esd_track_propagator_setup(cont->GetPropagator(),0.1*esd->GetMagneticField(), 520);

  printf( "Number of ESD Tracks : %d \n", esd->GetNumberOfTracks());

  for (Int_t iter = 0; iter < esd->GetNumberOfTracks(); ++iter) {
    AliEveTrack* track = esd_make_track(esd->GetTrack(iter), cont);
    cont->AddElement(track);
  }
  
  cont->SetTitle(Form("N=%d", esd->GetNumberOfTracks()) );
  cont->MakeTracks();

  return 0;
}


// Int_t tpc_clusters(TEveElement* cont=0, Float_t maxR=270)
// {


//   AliTPCClustersRow *clrow = new AliTPCClustersRow();
//   clrow->SetClass("AliTPCclusterMI");
//   clrow->SetArray(kMaxCl);
//   cTree->SetBranchAddress("Segment", &clrow);

//   tTPCClusters->SetOwnIds(kTRUE);


//   Float_t maxRsqr = maxR*maxR;
//   Int_t nentr=(Int_t)cTree->GetEntries();
//   for (Int_t i=0; i<nentr; i++)
//   {
//     if (!cTree->GetEvent(i)) continue;

//     TClonesArray *cl = clrow->GetArray();
//     Int_t ncl = cl->GetEntriesFast();

//     while (ncl--)
//     {
//       AliCluster *c = (AliCluster*) cl->UncheckedAt(ncl);
//       Float_t g[3]; //global coordinates
//       c->GetGlobalXYZ(g);
//       if (g[0]*g[0]+g[1]*g[1] < maxRsqr)
//       {
// 	clusters->SetNextPoint(g[0], g[1], g[2]);
// 	AliCluster *atp = new AliCluster(*c);
// 	clusters->SetPointId(atp);
//       }
//     }
//     cl->Clear();
//   }

//   delete clrow;

//   if (clusters->Size() == 0 && gEve->GetKeepEmptyCont() == kFALSE)
//   {
//     Warning("tpc_clusters.C", "No TPC clusters");
//     delete clusters;
//     return 1;
//   }

//   char form[1000];
//   sprintf(form,"TPC Clusters");
//   gTPCClusters->SetName(form);

//   char tip[1000];
//   sprintf(tip,"N=%d", gTPCClusters->Size());
//   gTPCClusters->SetTitle(tip);

//   const TString viz_tag("TPC Clusters");
//   gTPCClusters->ApplyVizTag(viz_tag, "Clusters");

//   return 0;
// }



 Int_t processTPCClusters(AliHLTHOMERBlockDesc* block, TEvePointSet *  cont) {
  Int_t iResult = 0;
  
  Int_t   slice = block->GetSubDetector().Atoi();
  Int_t   patch = block->GetSubSubDetector().Atoi();
  Float_t phi   = ( slice + 0.5 ) * TMath::Pi() / 9.0;  
  Float_t cos   = TMath::Cos( phi );
  Float_t sin   = TMath::Sin( phi );
    
  AliHLTTPCClusterData *cd = (AliHLTTPCClusterData*) block->GetData();
  UChar_t *data            = (UChar_t*) cd->fSpacePoints;

  if ( cd->fSpacePointCnt == 0 ) {
    printf ("No Clusters found in sector %d patch %d.\n", slice, patch );
    iResult = -1;
  } 
  else {
    
    for (Int_t ii = 0; ii < cd->fSpacePointCnt; ++ii, data += sizeof(AliHLTTPCSpacePointData)) {
      AliHLTTPCSpacePointData *sp = (AliHLTTPCSpacePointData *) data;
      
      cont->SetNextPoint(cos*sp->fX - sin*sp->fY, sin*sp->fX + cos*sp->fY, sp->fZ);
    }
  }
  
  return iResult;
}
 
