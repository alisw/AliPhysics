//-*- Mode: C++ -*-

// ** USED macros :
// ***************************************************
// - hlt_alieve_init.C
// - VizDB_scan.C
// - geom_gentle_hlt.C
// - geom_gentle_muon.C
// ***************************************************

#ifndef __CINT__
#include "unistd.h"
#include <TEvePointSet.h>
#include "EveBase/AliEveEventManager.h"
#include <AliCluster.h>
#include <TPC/AliTPCClustersRow.h>
#include <EveHLT/AliEveHOMERManager.h>
#include <TTimer.h>
#endif

class TEveTrackList;
class TEveProjectionManager;
class TEveGeoShape;
class TEveUtil;
class AliEveMacroExecutor;
class TEveScene;
class TEveElement;
class TEveText;
class AliHLTTriggerDecision;
class TEvePointSet;
class AliHLTHOMERBlockDesc;

class TEveViewer;
// -----------------------------------------------------------------
// --                       Geometry / Scenes                     --
// -----------------------------------------------------------------

TEveGeoShape *gGeomGentle     = 0;
TEveGeoShape *gGeomGentleRPhi = 0;
TEveGeoShape *gGeomGentleRhoZ = 0;
TEveGeoShape *gGeomGentleTRD  = 0;
TEveGeoShape *gGeomGentleMUON = 0;

TEveScene *gRPhiGeomScene  = 0;
TEveScene *gRhoZGeomScene  = 0;
TEveScene *gRPhiEventScene = 0;
TEveScene *gRhoZEventScene = 0;

TEveProjectionManager *gRPhiMgr = 0;
TEveProjectionManager *gRhoZMgr = 0;

TEveViewer *g3DView   = 0;
TEveViewer *gRPhiView = 0;
TEveViewer *gRhoZView = 0;

// -----------------------------------------------------------------
// --                Geometry / Scenes Parameters                 --
// -----------------------------------------------------------------

// -- Parameters to show different geometries
Bool_t gShowMUON     = kTRUE;
Bool_t gShowMUONRPhi = kFALSE;
Bool_t gShowMUONRhoZ = kTRUE;
Bool_t gShowTRD      = kFALSE;

Bool_t gCenterProjectionsAtPrimaryVertex = kFALSE;

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

// -- Text output members
TEveText*                                 gHLTText         = 0;

// -- Tracks members
TEveTrackList*                            gTPCTrack        = 0;

// -- Canvas for histos
TCanvas*                                  gCanvas          = 0;

// -----------------------------------------------------------------
// --                          Methods                            --
// -----------------------------------------------------------------

Int_t initializeEveViewer( Bool_t showExtraGeo );

Int_t nextEvent();

Int_t processEvent();

//Int_t processPHOSClusters( AliHLTHOMERBlockDesc* block);

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

  gHomerManager->ConnectEVEtoHOMER("TPC" );
}

// -----------------------------------------------------------------
Int_t initializeEveViewer( Bool_t showMuonGeo ) {
  
  //==============================================================================
  // Geometry, scenes, projections and viewers
  //==============================================================================

  TEveBrowser         *browser = gEve->GetBrowser();
  browser->ShowCloseTab(kFALSE);
  
  // -- Disable extra geometry
  // ---------------------------
  if ( ! showMuonGeo ) {
    gShowMUON = gShowMUONRPhi = gShowMUONRhoZ = kFALSE;
  }
  
  // -- Load Geometry
  // ------------------
  TEveUtil::LoadMacro("geom_gentle_hlt.C");
  gGeomGentle = geom_gentle_hlt();
  gGeomGentleRPhi = geom_gentle_rphi(); gGeomGentleRPhi->IncDenyDestroy();
  gGeomGentleRhoZ = geom_gentle_rhoz(); gGeomGentleRhoZ->IncDenyDestroy();
  gGeomGentleTRD  = geom_gentle_trd();

  if (gShowMUON) {
    TEveUtil::LoadMacro("geom_gentle_muon.C");
    gGeomGentleMUON = geom_gentle_muon(kFALSE);
  }

  // -- Scenes
  // -----------

  gRPhiGeomScene  = gEve->SpawnNewScene("RPhi Geometry",
                    "Scene holding projected geometry for the RPhi view.");
  gRhoZGeomScene  = gEve->SpawnNewScene("RhoZ Geometry",
		    "Scene holding projected geometry for the RhoZ view.");
  gRPhiEventScene = gEve->SpawnNewScene("RPhi Event Data",
		    "Scene holding projected geometry for the RPhi view.");
  gRhoZEventScene = gEve->SpawnNewScene("RhoZ Event Data",
		    "Scene holding projected geometry for the RhoZ view.");

  // -- Projection managers
  // ------------------------

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
    gRPhiGeomScene->AddElement(a);
  }
  gRPhiMgr->SetCurrentDepth(-10);
  gRPhiMgr->ImportElements(gGeomGentleRPhi, gRPhiGeomScene);
  gRPhiMgr->SetCurrentDepth(0);
  if (gShowTRD)      gRPhiMgr->ImportElements(gGeomGentleTRD, gRPhiGeomScene);
  if (gShowMUONRPhi) gRPhiMgr->ImportElements(gGeomGentleMUON, gRPhiGeomScene);

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
    gRhoZGeomScene->AddElement(a);
  }
  gRhoZMgr->SetCurrentDepth(-10);
  gRhoZMgr->ImportElements(gGeomGentleRhoZ, gRhoZGeomScene);
  gRhoZMgr->SetCurrentDepth(0);
  if (gShowTRD)      gRhoZMgr->ImportElements(gGeomGentleTRD, gRhoZGeomScene);
  if (gShowMUONRhoZ) gRhoZMgr->ImportElements(gGeomGentleMUON, gRhoZGeomScene);

  // -- Viewers
  // ------------

  TEveWindowSlot *slot = 0;
  TEveWindowPack *pack = 0;

  slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  pack = slot->MakePack();
  pack->SetElementName("Multi View");
  pack->SetHorizontal();
  pack->SetShowTitleBar(kFALSE);
  pack->NewSlot()->MakeCurrent();
  g3DView = gEve->SpawnNewViewer("3D View", "");
  g3DView->AddScene(gEve->GetGlobalScene());
  g3DView->AddScene(gEve->GetEventScene());

  pack = pack->NewSlot()->MakePack();
  pack->SetShowTitleBar(kFALSE);
  pack->NewSlot()->MakeCurrent();
  gRPhiView = gEve->SpawnNewViewer("RPhi View", "");
  gRPhiView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  gRPhiView->AddScene(gRPhiGeomScene);
  gRPhiView->AddScene(gRPhiEventScene);

  pack->NewSlot()->MakeCurrent();
  gRhoZView = gEve->SpawnNewViewer("RhoZ View", "");
  gRhoZView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
  gRhoZView->AddScene(gRhoZGeomScene);
  gRhoZView->AddScene(gRhoZEventScene);

  // -- List of Viewers
  // --------------------

  TEveViewerList *viewerlist = new TEveViewerList();
  viewerlist->AddElement(gEve->GetDefaultViewer());

  viewerlist->AddElement(g3DView);
  viewerlist->AddElement(gRhoZView);
  viewerlist->AddElement(gRPhiView);
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
  
//  gCanvas = new TCanvas("c1","c1", 600, 400);
  slot->StopEmbedding("QA histograms");
#endif


  //==============================================================================
  // Additional GUI components
  //==============================================================================
  
  slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  TEveWindowTab *storeTab = slot->MakeTab();
  store_tab->SetElementNameTitle("WindowStore",
   				 "Undocked windows whose previous container is not known\n"
 				 "are placed here when the main-frame is closed.");
  gEve->GetWindowManager()->SetDefaultContainer(storeTab);
  
  return 0;
}

// -----------------------------------------------------------------
void nextEvent() {

  if ( gHomerManager->NextEvent() )
    return;
  
  processEvent();
}

// -----------------------------------------------------------------
Int_t processEvent() {

  Int_t iResult = 0;

  gStyle->SetPalette(1, 0);
  gEve->DisableRedraw();

  // -- Reset
  // ----------
  if ( gTPCClusters ) gTPCClusters->Reset();
  if ( gPHOSClusters ) gPHOSClusters->Reset();
  if ( gTPCTrack )    gTPCTrack->DestroyElements();

  if ( gHomerManager == NULL) {
    printf ("No BlockList ... ");
    return -1;
  }
  if (gHomerManager->GetBlockList()->IsEmpty() ) {
    printf ("No Blocks in list ... ");
    return -2;
  }

  TIter next(gHomerManager->GetBlockList());
  AliHLTHOMERBlockDesc* block = 0;

  // -- Iterate over blocks in the block list
  // ------------------------------------------
  while ((block = (AliHLTHOMERBlockDesc*)next())) {
        
#if 0
    printf( "------------------- xxxxxxxxxxxxxxx ----------------------\n");
    printf( "Detector           : %s\n", block->GetDetector().Data() );
    printf( "Datatype           : %s\n", block->GetDataType().Data() );
    if (block->IsTObject() )
      printf( "Is TObject of class: %s\n", block->GetClassName().Data() );
    printf( "------------------- xxxxxxxxxxxxxxx ----------------------\n");
#endif

    // -- CHECK SOURCE
    // -----------------------------------------------------
    
    // ++ HLT BLOCK
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if ( ! block->GetDetector().CompareTo("HLT") ) {

      // -- ESDTREE
      if ( ! block->GetDataType().CompareTo("ALIESDV0") ) {
	if(!gTPCTrack){
	  gTPCTrack = new TEveTrackList("ESD Tracks");
	  gTPCTrack->SetMainColor(6);
	  gEve->AddElement(gTPCTrack);
	}
	iResult = processEsdTracks(block, gTPCTrack);
      }
      
      // -- Process ROOTObj
      else if ( ! block->GetDataType().CompareTo("ROOTTOBJ") ) {
	if(!gHLTText){
	  gHLTText = new TEveText();
	  //gHLTText->BBoxZero(5, -5, -5, 0);
	  //gHLTText->SetExtrude(25);
	  //gHLTText->AssertBBoxExtents(25,25,25);
	  gEve->AddElement(gHLTText);
	} 
 	processROOTTOBJ( block, gHLTText );
      } 

      // -- Process HLT RDLST
      else if ( ! block->GetDataType().CompareTo("HLTRDLST") ) {
 	;
	//cout<<"Readlist"<<endl;
	//processHLTRDLST( block );
      }
    } // if ( ! block->GetDetector().CompareTo("HLT") ) {

    // ++ TPC BLOCK
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    else if ( ! block->GetDetector().CompareTo("TPC") ) {
      
      // -- Process TPC Clusters
      if ( ! block->GetDataType().CompareTo("CLUSTERS") ) {
	if(!gTPCClusters){
	  gTPCClusters = new TEvePointSet("TPC Clusters");
	  gTPCClusters->SetMainColor(kRed);
	  gTPCClusters->SetMarkerStyle((Style_t)kFullDotSmall);
	  gEve->AddElement(gTPCClusters);
	} 
	iResult = processTPCClusters( block , gTPCClusters);
      }
    } // else if ( ! block->GetDetector().CompareTo("TPC") ) {


    // ++ SPD  BLOCK
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    else if ( ! block->GetDetector().CompareTo("ISPD") ){
      if ( block->GetDataType().CompareTo("CLUSTERS") == 0 ) {
 	//processISPDClusters( block );
      } 
    } // else if ( ! block->GetDetector().CompareTo("ISPD") ){


    // -- ITS
   else if ( ! block->GetDetector().CompareTo("ITS") ){
      if ( block->GetDataType().CompareTo("ROOTHIST") == 0 ) {
	iResult = 0;
	//iResult = processITSHist( block );
      } 
    } // else if ( ! block->GetDetector().CompareTo("ISPD") ){


    // ---------------------------------------------------------
  } // while ((block = (AliHLTHOMERBlockDesc*)next())) {


  if ( gTPCClusters ) gTPCClusters->ResetBBox();
  if ( gPHOSClusters ) gPHOSClusters->ResetBBox();
  if ( gSPDClusters ) gSPDClusters->ResetBBox();
  if ( gTPCTrack ) gTPCTrack->ElementChanged();


  // -- Set EventID in Window Title  
  // --------------------------------------------
  TString winTitle("Eve Main Window -- Event ID : ");
  winTitle += Form("0x%016X ", gHomerManager->GetEventID() );
  gEve->GetBrowser()->SetWindowName(winTitle);

  // -- Set Projections
  // --------------------------------------------

  // XXX Primary vertex ... to be retrieved from the ESD
  // Double_t x[3] = { 0, 0, 0 };

  TEveElement* top = gEve->GetCurrentEvent();

  if (gRPhiMgr && top) {
    gRPhiEventScene->DestroyElements();
    if (gCenterProjectionsAtPrimaryVertex)
      gRPhiMgr->SetCenter(x[0], x[1], x[2]);
    gRPhiMgr->ImportElements(top, gRPhiEventScene);
  }
  
  if (gRhoZMgr && top) {
    gRhoZEventScene->DestroyElements();
    if (gCenterProjectionsAtPrimaryVertex)
      gRhoZMgr->SetCenter(x[0], x[1], x[2]);
    gRhoZMgr->ImportElements(top, gRhoZEventScene);
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
Int_t processITSHist(AliHLTHOMERBlockDesc* block) {
  TH2F* hist = block->GetTObject();
  
  gCanvas->cd();
  hist->Draw();
  return 0;
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
Int_t processROOTTOBJ(AliHLTHOMERBlockDesc* block, TEveText* et) {
  
  // -- AliHLTGlobalTriggerDecision
  if ( ! block->GetClassName().CompareTo("AliHLTGlobalTriggerDecision") ) {

    AliHLTGlobalTriggerDecision *trig = dynamic_cast<AliHLTGlobalTriggerDecision*> block->GetTObject();
    trig->Print(); 
    
    // et->SetText("balle");;

    // TEveText* tt = new TEveText("Trigger: Class is known ;-) ");
    // gEve->AddElement(tt);

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

// -----------------------------------------------------------------
Int_t processTPCClusters(AliHLTHOMERBlockDesc* block, TEvePointSet* cont) {
  
  Int_t   slice = block->GetSubDetector();
  Int_t   patch = block->GetSubSubDetector();
  Float_t phi   = ( slice + 0.5 ) * TMath::Pi() / 9.0;  
  Float_t cos   = TMath::Cos( phi );
  Float_t sin   = TMath::Sin( phi );
  
  AliHLTTPCClusterData *cd = reinterpret_cast<AliHLTTPCClusterData*> (block->GetData());
  UChar_t *data            = reinterpret_cast<UChar_t*> (cd->fSpacePoints);

  if ( cd->fSpacePointCnt != 0 ) {
    for (Int_t iter = 0; iter < cd->fSpacePointCnt; ++iter, data += sizeof(AliHLTTPCSpacePointData)) {
      AliHLTTPCSpacePointData *sp = reinterpret_cast<AliHLTTPCSpacePointData*> (data);
      cont->SetNextPoint(cos*sp->fX - sin*sp->fY, sin*sp->fX + cos*sp->fY, sp->fZ);
    }
  }
  
  return 0;
}
 
