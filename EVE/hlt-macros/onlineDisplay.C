//-*- Mode: C++ -*-

// ** USED macros :
// ***************************************************
// - hlt_alieve_init.C
// - VizDB_scan.C
// - geom_gentle_hlt.C
// - geom_gentle_muon.C
// ***************************************************

#if !defined(__CINT__) || defined(__MAKECINT__)

//****************** ROOT ******************************************
#include "TTimer.h"
#include "TRandom.h"
#include "TVirtualPad.h"
#include "TGLViewer.h"
#include "TThread.h"
#include "TGFileBrowser.h"
#include "TStyle.h"
#include "TList.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TVector3.h"

//****************** ROOT/EVE **************************************
#include "TEveManager.h"
#include "TEvePointSet.h"
#include "TEveTrack.h"
#include "TEveVSDStructs.h"
#include "TEveTrackPropagator.h"
#include "TEveScene.h"
#include "TEveElement.h"
#include "TEveUtil.h"
#include "TEveEventManager.h"
#include "TEveProjectionAxes.h"
#include "TEveWindowManager.h"
#include "TEveViewer.h"
#include "TEveText.h"
#include "TEveProjectionManager.h"
#include "TEveGeoShape.h"
#include "TEveBoxSet.h"
#include "TEveTrans.h"
#include "TEveRGBAPalette.h"

//****************** AliRoot ***************************************
#include "AliESDEvent.h"
#include "AliCDBManager.h"
#include "AliRawReaderMemory.h"
#include "AliTPCRawStream.h"
#include "AliGeomManager.h"

//****************** AliRoot/EVE ***********************************
#include "AliHLTHOMERManager.h"
#include "AliEveHOMERManager.h"
#include "AliEveTPCLoader.h" 
#include "AliEveTPCData.h"
#include "AliEveITSDigitsInfo.h"
#include "AliEveITSModule.h"
#include "AliEveMacroExecutor.h"
#include "AliEveMacro.h"
#include "AliEveTrack.h"

//****************** AliRoot/HLT ***********************************
#include "AliHLTHOMERBlockDesc.h"
#include "AliHLTHOMERReader.h"
#include "AliHLTMUONUtils.h"
#include "AliHLTMUONDataBlockReader.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTGlobalTriggerDecision.h"
#include "AliHLTTPCCATrackParam.h"
#include "AliHLTCaloClusterDataStruct.h"
#include "AliHLTCaloClusterReader.h"

//****************** AliRoot/MUON **********************************
#include "AliMUONCalibrationData.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryDetElement.h"

#include "AliMpCDB.h"
#include "AliMpDDLStore.h"
#include "AliMpDetElement.h"
#include "AliMpDEIterator.h"
#include "AliMpVSegmentation.h"
#include "AliMpSegmentation.h"
#include "AliMpTriggerCrate.h"
#include "AliMpLocalBoard.h"

// ***************** AliRoot/ITS **********************************
#include "AliITSRecPoint.h"

//****************** AliRoot/TRD ***********************************
#include "AliHLTTRDCluster.h"
#include "AliTRDcluster.h"

//#####################AliRoot PHOS ##################################3
#include "AliPHOSGeometry.h"

//****************** Macros ****************************************
#include "hlt_structs.C"
#include "hlt_alieve_init.C"
#include "geom_gentle_hlt.C"
#include "alice-macros/esd_tracks.C"

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
class TEvePointSetArray;
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
TTimer                                    eventTimerFast;

// -- HOMERManager
AliEveHOMERManager*                       gHomerManager      = 0;

// -- Geometry Manager 
TGeoManager*                              gGeoManager        = 0;
AliPHOSGeometry*                          gPHOSGeom          = 0;

// -- Cluster members
TEvePointSet*                             gSPDClusters       = 0;
TEvePointSet*                             gSSDClusters       = 0;
TEvePointSet*                             gSDDClusters       = 0;
TEvePointSet*                             gTRDClusters       = 0;
TEvePointSetArray*                        gTRDColClusters    = 0;
TEvePointSet*                             gTPCClusters       = 0;
TEvePointSetArray*                        gTPCColClusters    = 0;
TEveBoxSet*                               gPHOSBoxSet[5]     = {0, 0, 0, 0, 0}; 
TEvePointSet*                             gMUONClusters      = 0;

// -- Text output members
TEveText*                                 gHLTText           = 0;

// -- Tracks members
TEveTrackList*                            gTPCTrack          = 0;

// -- Canvas for histograms
TCanvas*                                  gTRDCanvas         = 0;
TCanvas*                                  gTPCCanvas         = 0;

// -- TRD --
Int_t                                     gTRDHistoCount     = 0;
Int_t                                     gTRDEvents         = 0;
Int_t                                     gTRDBins           = 12;

// -- TPC --
Int_t                                     gTPCBins           = 15;
TH1F*                                     gTPCCharge         = 0;
TH1F*                                     gTPCQMax           = 0;
TH1F*                                     gTPCQMaxOverCharge = 0;

// -- PHOS --
TEveElementList*                          gPHOSElementList   = 0;

// --- Flag if eventloop is running
Bool_t                                    gEventLoopStarted = kFALSE;

// -----------------------------------------------------------------
// --                          Methods                            --
// -----------------------------------------------------------------

Int_t initializeEveViewer( Bool_t TPCMode, Bool_t MUONMode, Bool_t TRDMode );

void nextEvent();

Int_t processEvent();

Int_t processPHOSClusters( AliHLTHOMERBlockDesc* block);

Int_t processEsdTracks( AliHLTHOMERBlockDesc* block, TEveTrackList* cont );

Int_t processHLTRDLST( AliHLTHOMERBlockDesc* block );

Int_t processROOTTOBJ( AliHLTHOMERBlockDesc* block, TEveText* text );

Int_t processTPCClusters( AliHLTHOMERBlockDesc * block, TEvePointSet *cont, TEvePointSetArray *contCol );

Int_t processTRDClusters( AliHLTHOMERBlockDesc * block, TEvePointSet *cont, TEvePointSetArray *contCol);

Int_t processTRDHistograms (AliHLTHOMERBlockDesc * block, TCanvas * canvas );

Int_t processMUONClusters( AliHLTHOMERBlockDesc* block);

Int_t processITSClusters(AliHLTHOMERBlockDesc* block, TEvePointSet* cont);

Int_t processITSHist(AliHLTHOMERBlockDesc* block);

// #################################################################
// #################################################################
// #################################################################

// -----------------------------------------------------------------
void onlineDisplay(Bool_t TPCMode = kTRUE, Bool_t MUONMode = kFALSE, Bool_t TRDMode = kFALSE) {

  // -- Loading Geometry
  // ---------------------
  Int_t run = 67179;
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetRun(run);
  AliGeomManager::LoadGeometry();

  // Get the pointer to gGeoManager before it's broken (bug in alieve)
  TGeoManager *gm = gGeoManager;

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
  //TEveUtil::LoadMacro("hlt_alieve_init.C");
  //hlt_alieve_init(".", -1);

  // -- Initialize Eve
  // -------------------
  initializeEveViewer( TPCMode, MUONMode, TRDMode);

  // -- Reset gGeoManager to the original pointer
  // ----------------------------------------------
  gGeoManager = gm;
  gPHOSGeom = AliPHOSGeometry::GetInstance("IHEP", "IHEP");

  // -- Finalize Eve
  // -----------------
  gSystem->ProcessEvents();
  gEve->Redraw3D(kTRUE);

  if ( TPCMode ) {
    gHomerManager->ConnectEVEtoHOMER("TPC" );
  } else if ( MUONMode ) {
    gHomerManager->ConnectEVEtoHOMER("MUON");
  } else if( TRDMode ) {
    gHomerManager->ConnectEVEtoHOMER("TRD");  
  } else {
    cout<<" No detectors selected, nothing will be displayed"<<endl;
  }	

  g_esd_tracks_true_field = kFALSE;
}

// -------------------------------------------------------------------------
Int_t initializeEveViewer( Bool_t TPCMode, Bool_t MUONMode, Bool_t TRDMode) {
  
  //==============================================================================
  // -- Geometry, scenes, projections and viewers
  //==============================================================================

  TEveBrowser         *browser = gEve->GetBrowser();
  browser->ShowCloseTab(kFALSE);
  
  // -- Disable extra geometry
  // ---------------------------
  if (!MUONMode)
    gShowMUON = gShowMUONRPhi = gShowMUONRhoZ = kFALSE;
  
  // -- Load Geometry
  // ------------------
  TEveUtil::LoadMacro("geom_gentle_hlt.C");
  gGeomGentle = geom_gentle_hlt();
  gGeomGentleRPhi = geom_gentle_rphi(); gGeomGentleRPhi->IncDenyDestroy();
  gGeomGentleRhoZ = geom_gentle_rhoz(); gGeomGentleRhoZ->IncDenyDestroy();
  gGeomGentleTRD  = geom_gentle_trd();

  if (gShowMUON) 
    gGeomGentleMUON = geom_gentle_muon(kFALSE);
  
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
  gRPhiMgr->ImportElements(gGeomGentleTRD, gRPhiGeomScene);
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

  TEveWindowSlot *slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  TEveWindowPack *pack = slot->MakePack();
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
  // -- Macros / QA histograms
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
# endif


  //==============================================================================
  // -- Histograms
  //==============================================================================

  if(TRDMode) {
    slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
    slot->StartEmbedding();
    
    gTRDCanvas = new TCanvas("canvasTRD","canvasTRD", 600, 400);
    gTRDCanvas->Divide(3,2);
    slot->StopEmbedding("TRD histograms");
  }
  else if(TPCMode){
    slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
    slot->StartEmbedding();
    
    gTPCCanvas = new TCanvas("canvasTPC","canvasTPC", 600, 400);
    gTPCCharge = new TH1F("ClusterCharge","ClusterCharge",100,0,500);
    gTPCQMax = new TH1F("QMax","QMax",50,0,250);
    gTPCQMaxOverCharge = new TH1F("QMaxOverCharge","QMaxOverCharge",50,0,1);
    slot->StopEmbedding("TPC histograms");
  }

  //==============================================================================
  // -- Additional GUI components
  //==============================================================================
  
  slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  TEveWindowTab *storeTab = slot->MakeTab();
  storeTab->SetElementNameTitle("WindowStore",
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

  //==============================================================================
  // -- Reset
  //==============================================================================

  if ( gTRDCanvas ) {
    gTRDCanvas->Clear();
    gTRDCanvas->Divide(3,2);
  }

  if ( gTPCTrack )     gTPCTrack->DestroyElements();

  if ( gSPDClusters )  gSPDClusters->Reset();
  if ( gSSDClusters )  gSSDClusters->Reset();
  if ( gSDDClusters )  gSSDClusters->Reset();
  if ( gTPCClusters )  gTPCClusters->Reset();
  if ( gTRDClusters )  gTRDClusters->Reset();
  if ( gMUONClusters ) gMUONClusters->Reset();

  if ( gPHOSBoxSet[1] )
    for(int im = 0; im < 5; im++)
      gPHOSBoxSet[im]->Reset();   

  if ( gTPCColClusters )
    for (Int_t ii = 0; ii <= gTPCBins+1; ++ii) 
      gTPCColClusters->GetBin(ii)->Reset();
  
  if ( gTRDColClusters )
    for (Int_t ii = 0; ii <= gTRDBins+1; ++ii) 
      gTRDColClusters->GetBin(ii)->Reset();

  gTRDHistoCount = 0;

  //==============================================================================
  // -- Process Blocks
  //==============================================================================

  if ( gHomerManager->GetBlockList() == NULL) {
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
        
#if DEBUG
    printf( "------------------- xxxxxxxxxxxxxxx ----------------------\n");
    printf( "Detector           : %s\n", block->GetDetector().Data() );
    printf( "Datatype           : %s\n", block->GetDataType().Data() );
    if (block->IsTObject() )
      printf( "Is TObject of class: %s\n", block->GetClassName().Data() );
    printf( "------------------- xxxxxxxxxxxxxxx ----------------------\n");
#endif

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
	gTPCTrack->ElementChanged();
      } // if ( ! block->GetDataType().CompareTo("ALIESDV0") ) {
      
      // -- Process ROOTObj
      else if ( ! block->GetDataType().CompareTo("ROOTTOBJ") ) {
	if(!gHLTText){
	  //gHLTText = new TEveText();
	  //gHLTText->BBoxZero(5, -5, -5, 0);
	  //gHLTText->SetExtrude(25);
	  //gHLTText->AssertBBoxExtents(25,25,25);
	  //gEve->AddElement(gHLTText);
	} 
 	processROOTTOBJ( block, gHLTText );
      } // else if ( ! block->GetDataType().CompareTo("ROOTTOBJ") ) {

      // -- Process HLT RDLST
      else if ( ! block->GetDataType().CompareTo("HLTRDLST") ) {
	processHLTRDLST( block );
      } // else if ( ! block->GetDataType().CompareTo("HLTRDLST") ) {

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

	if(!gTPCColClusters){
	  gTPCColClusters = new TEvePointSetArray("TPC Clusters Colorized");
	  gTPCColClusters->SetMainColor(kRed);
	  gTPCColClusters->SetMarkerStyle(4); // antialiased circle
	  gTPCColClusters->SetMarkerSize(0.8);
	  gTPCColClusters->InitBins("Cluster Charge", gTPCBins, 0., gTPCBins*20.);

	  const Int_t nCol = TColor::GetNumberOfColors();
	  for (Int_t ii = 0; ii < gTPCBins+1; ++ii)
	    gTPCColClusters->GetBin(ii)->SetMainColor(TColor::GetColorPalette(ii * nCol / (gTPCBins+2)));
	  
	  gEve->AddElement(gTPCColClusters);
	} 

	iResult = processTPCClusters(block, gTPCClusters, gTPCColClusters);
	gTPCClusters->ElementChanged();
	gTPCColClusters->ElementChanged();
      }
      
    } // else if ( ! block->GetDetector().CompareTo("TPC") ) {

    // ++ TRD BLOCK
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    else if ( ! block->GetDetector().CompareTo("TRD") ) {
       
      // -- Process TRD Clusters
      if ( ! block->GetDataType().CompareTo("CLUSTERS") ) {

	if(!gTRDClusters){
	  gTRDClusters = new TEvePointSet("TRD Clusters");
	  gTRDClusters->SetMainColor(kBlue);
	  gTRDClusters->SetMarkerStyle((Style_t)kFullDotSmall);
	  gEve->AddElement(gTRDClusters);
	} 

	if(!gTRDColClusters){
	  gTRDColClusters = new TEvePointSetArray("TRD Clusters Colorized");
	  gTRDColClusters->SetMainColor(kRed);
	  gTRDColClusters->SetMarkerStyle(4); // antialiased circle
	  //	  gTRDColClusters->SetMarkerStyle((Style_t)kFullDotSmall);
	  gTRDColClusters->SetMarkerSize(0.8);
	  gTRDColClusters->InitBins("Cluster Charge", gTRDBins, 0., gTRDBins*100.);

	  //TColor::SetPalette(1, 0); // Spectrum palette
	  const Int_t nCol = TColor::GetNumberOfColors();
	  for (Int_t ii = 0; ii < gTRDBins+1; ++ii)
	    gTRDColClusters->GetBin(ii)->SetMainColor(TColor::GetColorPalette(ii * nCol / (gTRDBins+2)));
	  
	  gEve->AddElement(gTRDColClusters);
	} 

	iResult = processTRDClusters( block, gTRDClusters, gTRDColClusters );
	gTRDClusters->ElementChanged();
	gTRDColClusters->ElementChanged();
      }

      // -- Process TRD Histograms
      else if ( block->GetDataType().CompareTo("ROOTHIST") == 0 ) {
	iResult = processTRDHistograms( block, gTRDCanvas );     
	if ( gTRDCanvas) gTRDCanvas->Update();
      }

    } // else if ( ! block->GetDetector().CompareTo("TRD") ) {
    
    // ++ MUON BLOCK
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    else if ( ! block->GetDetector().CompareTo("MUON") && gShowMUON ) {

      // -- Process MUON Clusters
      if ( (block->GetDataType().CompareTo("RECHITS") == 0) || (block->GetDataType().CompareTo("TRIGRECS") == 0) ) {
	
	if ( !gMUONClusters ) {
	  gMUONClusters = new TEvePointSet("MUON RecHits");
	  gMUONClusters->SetMainColor(kBlue);
	  gMUONClusters->SetMarkerStyle(20);
	  gEve->AddElement(gMUONClusters);
	}
	
 	processMUONClusters( block );
	gMUONClusters->ElementChanged();
	
      } 
    } // else if ( ! block->GetDetector().CompareTo("MUON") && gShowMUON ) {

    // ++ ISPD BLOCK
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    else if ( ! block->GetDetector().CompareTo("ISPD") ){
      if ( block->GetDataType().CompareTo("CLUSTERS") == 0 ) {
	
	if(!gSPDClusters){
	  gSPDClusters = new TEvePointSet("SPD Clusters");
	  gSPDClusters->SetMainColor(kBlack);
	  gSPDClusters->SetMarkerStyle((Style_t)kFullDotLarge);
	  gEve->AddElement(gSPDClusters);
	} 
 	
	processITSClusters( block , gSPDClusters);
	gSPDClusters->ElementChanged();
      } 
    } // else if ( ! block->GetDetector().CompareTo("ISPD") ){

    // ++ ISDD BLOCK
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    else if ( ! block->GetDetector().CompareTo("ISDD") ){
      if ( block->GetDataType().CompareTo("CLUSTERS") == 0 ) {
	
	if(!gSDDClusters){
	  gSDDClusters = new TEvePointSet("SDD Clusters");
	  gSDDClusters->SetMainColor(kPink);
	  gSDDClusters->SetMarkerStyle((Style_t)kFullDotLarge);
	  gEve->AddElement(gSDDClusters);
	} 
 	
	processITSClusters( block , gSDDClusters);
	gSDDClusters->ElementChanged();
      } 
    } // else if ( ! block->GetDetector().CompareTo("ISDD") ){

    // ++ ISSD BLOCK
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    else if ( ! block->GetDetector().CompareTo("ISSD") ){
      if ( block->GetDataType().CompareTo("CLUSTERS") == 0 ) {
	
	if(!gSSDClusters){
	  gSSDClusters = new TEvePointSet("SSD Clusters");
	  gSSDClusters->SetMainColor(kPink);
	  gSSDClusters->SetMarkerStyle((Style_t)kFullDotLarge);
	  gEve->AddElement(gSSDClusters);
	} 
 	
	processITSClusters( block , gSSDClusters);
	gSSDClusters->ElementChanged();
      } 
    } // else if ( ! block->GetDetector().CompareTo("ISSD") ){
        
    // ++ ITS BLOCK
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    else if ( ! block->GetDetector().CompareTo("ITS") ){
      if ( block->GetDataType().CompareTo("ROOTHIST") == 0 ) {
	iResult = processITSHist( block );
      } 
    } // else if ( ! block->GetDetector().CompareTo("ITS") ){
    
    // ++ PHOS BLOCK
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    else if ( ! block->GetDetector().CompareTo("PHOS") ) {

      // -- Process Digits
      if ( block->GetDataType().CompareTo("DIGITTYP") == 0 ) {
	
	if( !gPHOSElementList ){
	  gPHOSElementList = new TEveElementList("PHOS Cells");
	  gEve->AddElement(gPHOSElementList);
	  
	  TVector3 center;
	  Float_t angle;
	  
	  // -- Create boxsets
	  for(int im = 0; im < 5; im++) {
	  
	    TEveRGBAPalette* pal = new TEveRGBAPalette(40,120);
	    gPHOSBoxSet[im] = new TEveBoxSet(Form("Cells Module %d" , im));
	    gPHOSBoxSet[im]->SetPalette(pal);
	    gPHOSBoxSet[im]->Reset(TEveBoxSet::kBT_AABox, kFALSE, 64);
	    gPHOSBoxSet[im]->SetOwnIds(kTRUE);
	    
	    gPHOSGeom->GetModuleCenter(center, "CPV", im+1);
	    angle = gPHOSGeom->GetPHOSAngle(im+1)*TMath::Pi()/180;
	  
	    gPHOSBoxSet[im]->RefitPlex();
	    TEveTrans& t = gPHOSBoxSet[im]->RefMainTrans();
	    t.SetupRotation(1, 2, angle );
	    t.SetPos(center.X(), center.Y(), center.Z());
	    
	    gPHOSElementList->AddElement(gPHOSBoxSet[im]);
	  }
	} // for(int im = 0; im < 5; im++) {
	
	iResult = processPHOSClusters( block );
	
	for(int im = 0; im < 5; im++)
	  gPHOSBoxSet[im]->ElementChanged();

      } // if ( block->GetDataType().CompareTo("DIGITTYP") == 0 ) {

    } // else if ( ! block->GetDetector().CompareTo("PHOS") ){

    // ---------------------------------------------------------
  } // while ((block = (AliHLTHOMERBlockDesc*)next())) {

  //==============================================================================
  // -- Update Objects
  //==============================================================================

  // -- TPC Histograms
  if ( gTPCCanvas && gTPCCharge && gTPCQMax) {
    gTPCCanvas->Clear();    
    gTPCCanvas->Divide(1,3);

    gTPCCanvas->cd(1);
    gTPCCharge->Draw();

    gTPCCanvas->cd(2);
    gTPCQMax->Draw();

    gTPCCanvas->cd(3);
    gTPCQMaxOverCharge->Draw();

    gTPCCanvas->Update();
  }

  if ( gTPCClusters ) gTPCClusters->ResetBBox();
  if ( gTRDClusters ) gTRDClusters->ResetBBox();
  if ( gSPDClusters ) gSPDClusters->ResetBBox();
  if ( gSDDClusters ) gSDDClusters->ResetBBox();
  if ( gSSDClusters ) gSSDClusters->ResetBBox();
  if ( gMUONClusters ) gMUONClusters->ResetBBox();

  if ( gPHOSBoxSet[1] )
    for(int im = 0; im < 5; im++)
      gPHOSBoxSet[im]->ResetBBox();      

  //==============================================================================
  // -- Set EventID in Window Title  
  // -- Update Objects
  //==============================================================================

  TString winTitle("Eve Main Window -- Event ID : ");
  winTitle += Form("0x%016X ", gHomerManager->GetEventID() );
  gEve->GetBrowser()->SetWindowName(winTitle);

  //==============================================================================
  // -- Set Projections
  //==============================================================================

  // XXX Primary vertex ... to be retrieved from the ESD
  Double_t x[3] = { 0, 0, 0 };
  
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

  //==============================================================================

  gEve->Redraw3D(0,1); // (0, 1)
  gEve->EnableRedraw(); 

  return iResult;
}

// -----------------------------------------------------------------
Int_t processITSHist(AliHLTHOMERBlockDesc* /*block*/) {
  return 0;
}

// -----------------------------------------------------------------
Int_t processHLTRDLST(AliHLTHOMERBlockDesc* /*block*/) {
  return 0;
}

// -----------------------------------------------------------------
Int_t processROOTTOBJ(AliHLTHOMERBlockDesc* block, TEveText* /*et*/) {
  
  // -- AliHLTGlobalTriggerDecision
  if ( ! block->GetClassName().CompareTo("AliHLTGlobalTriggerDecision") ) {

    AliHLTGlobalTriggerDecision *trig = dynamic_cast<AliHLTGlobalTriggerDecision*>( block->GetTObject());
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
    AliEveTrack* track = dynamic_cast<AliEveTrack*>(esd_make_track(esd->GetTrack(iter), cont));
    cont->AddElement(track);
  }
  
  cont->SetTitle(Form("N=%d", esd->GetNumberOfTracks()) );
  cont->MakeTracks();

  return 0;
}


// -----------------------------------------------------------------
// Int_t processPHOSClusters(AliHLTHOMERBlockDesc* block) {

//    AliHLTCaloClusterHeaderStruct *hd = reinterpret_cast<AliHLTCaloClusterHeaderStruct*> (block->GetData());
//    AliHLTCaloClusterReader * cReader = new AliHLTCaloClusterReader();
//    cReader->SetMemory(hd);
   
//    AliHLTCaloClusterDataStruct* cd = 0;
   
//    while( (cd = cReader->NextCluster()) ){

//      cout << cd->fEnergy << endl;
     
//      Float_t e = cd->fEnergy;
//      Float_t x = cd->fGlobalPos[0];
//      Float_t y = cd->fGlobalPos[1];
//      Float_t z = cd->fGlobalPos[2];

     
//      gPHOSBoxSet->AddBox(x, y, z, 2.2, e*10, 2.2);
     
//    }

//   return 0;
// }



// -----------------------------------------------------------------
Int_t processPHOSClusters(AliHLTHOMERBlockDesc* block) {

  AliHLTPHOSDigitDataStruct *ds = reinterpret_cast<AliHLTPHOSDigitDataStruct*> (block->GetData());
    UInt_t nClusters = block->GetSize()/sizeof(AliHLTPHOSDigitDataStruct);

  for(UInt_t i = 0; i < nClusters; i++, ds++) {
    gPHOSBoxSet[ds->fModule]->AddBox(ds->fX, 0, ds->fZ, 2.2, ds->fEnergy*20, 2.2);
    gPHOSBoxSet[ds->fModule]->DigitValue(static_cast<Int_t>(ds->fEnergy*100));
  }

  return 0;
}

// -----------------------------------------------------------------
Int_t processITSClusters(AliHLTHOMERBlockDesc* block, TEvePointSet* cont) {

  AliHLTITSClusterData *cd = reinterpret_cast<AliHLTITSClusterData*> (block->GetData());
  UChar_t *data            = reinterpret_cast<UChar_t*> (cd->fSpacePoints);
  
  if ( cd->fSpacePointCnt != 0 ) {
    for (Int_t iter = 0; iter < cd->fSpacePointCnt; ++iter, data += sizeof(AliHLTITSSpacePointData)) {
      AliHLTITSSpacePointData *sp = reinterpret_cast<AliHLTITSSpacePointData*> (data);
  
      Int_t lab[4]   = {0,0,0,0};
      Float_t hit[6] = {0,0,0,0,0,0};
      Int_t info[3]  = {0,0,0};
 				 
      lab[0]  = sp->fTracks[0];
      lab[1]  = sp->fTracks[1];
      lab[2]  = sp->fTracks[2];
      lab[3]  = sp->fIndex;
      hit[0]  = sp->fY;
      hit[1]  = sp->fZ;
      hit[2]  = sp->fSigmaY2;
      hit[3]  = sp->fSigmaZ2;
      hit[4]  = sp->fQ;
      hit[5]  = sp->fSigmaYZ;
      info[0] = sp->fNy;
      info[1] = sp->fNz;
      info[2] = sp->fLayer;
      
      Float_t xyz[3];
      AliITSRecPoint recpoint(lab,hit,info);
      recpoint.GetGlobalXYZ(xyz);

      cont->SetNextPoint(xyz[0], xyz[1], xyz[2]);
    }
  }
  return 0;
}

// -----------------------------------------------------------------
Int_t processTPCClusters(AliHLTHOMERBlockDesc* block, TEvePointSet* cont, TEvePointSetArray *contCol ) {

  Int_t   slice = block->GetSubDetector();
  Float_t phi   = ( slice + 0.5 ) * TMath::Pi() / 9.0;  
  Float_t cos   = TMath::Cos( phi );
  Float_t sin   = TMath::Sin( phi );
  
  AliHLTTPCClusterData *cd = reinterpret_cast<AliHLTTPCClusterData*> (block->GetData());
  UChar_t *data            = reinterpret_cast<UChar_t*> (cd->fSpacePoints);

  if ( cd->fSpacePointCnt != 0 ) {
    for (Int_t iter = 0; iter < cd->fSpacePointCnt; ++iter, data += sizeof(AliHLTTPCSpacePointData)) {
      AliHLTTPCSpacePointData *sp = reinterpret_cast<AliHLTTPCSpacePointData*> (data);
      cont->SetNextPoint(cos*sp->fX - sin*sp->fY, sin*sp->fX + cos*sp->fY, sp->fZ);
      contCol->Fill(cos*sp->fX - sin*sp->fY, sin*sp->fX + cos*sp->fY, sp->fZ, sp->fCharge);

      gTPCCharge->Fill(sp->fCharge);
      gTPCQMax->Fill(sp->fQMax);
      gTPCQMaxOverCharge->Fill(((Float_t)sp->fQMax)/((Float_t)sp->fCharge));
    }
  }
  
  return 0;
}

// -----------------------------------------------------------------
Int_t processMUONClusters(AliHLTHOMERBlockDesc* block) {
  
  Int_t iResult = 0;
  
  unsigned long size = block->GetSize();
  Int_t * buffer ;

  buffer = (Int_t *)block->GetData();
//   cout<<"block size : "<<size<<", buffer : "<<buffer<<", DataType : "<<block->GetDataType()<<endl;

// //   for(int idata=0;idata<int(size);idata++)
// //     printf("\tbuffer[%d] : %d\n",idata,buffer[idata]);
  
  
  
  if(block->GetDataType().CompareTo("RECHITS") == 0){

    AliHLTMUONRecHitsBlockReader trackblock((char*)buffer, size);
    const AliHLTMUONRecHitStruct* hit = trackblock.GetArray();
    
    for(AliHLTUInt32_t ientry = 0; ientry < trackblock.Nentries(); ientry++){
//       cout << setw(13) << left << hit->fX << setw(0);
//       cout << setw(13) << left << hit->fY << setw(0);
//       cout << hit->fZ << setw(0) << endl;
      if(hit->fX!=0.0 && hit->fY!=0.0 && hit->fZ!=0.0)
	gMUONClusters->SetNextPoint(hit->fX,hit->fY,hit->fZ);
      hit++;
      
    }// track hit loop
  }

  else{// if rechits
    //     if(!strcmp((BlockType(ULong64_t(reader->GetBlockDataType(i)))).Data(),"TRIGRECS")){
  
    AliHLTMUONTriggerRecordsBlockReader trigblock(buffer, size);
    const AliHLTMUONTriggerRecordStruct* trigrec = trigblock.GetArray();
    for(AliHLTUInt32_t ientry = 0; ientry < trigblock.Nentries(); ientry++){
      
      const AliHLTMUONRecHitStruct* hit = &trigrec->fHit[0];
      for(AliHLTUInt32_t ch = 0; ch < 4; ch++)
	{
	  cout << setw(10) << left << ch + 11 << setw(0);
	  cout << setw(13) << left << hit->fX << setw(0);
	  cout << setw(13) << left << hit->fY << setw(0);
	  cout << hit->fZ << setw(0) << endl;
	  if(hit->fX!=0.0 && hit->fY!=0.0 && hit->fZ!=0.0)
	    gMUONClusters->SetNextPoint(hit->fX,hit->fY,hit->fZ);
	  hit++;
	}// trig chamber loop
      trigrec++;
    }//trig hit loop
  }//else trigger

  return iResult;
}

 
// -----------------------------------------------------------------
Int_t processTRDClusters(AliHLTHOMERBlockDesc* block, TEvePointSet *cont, TEvePointSetArray *contCol) {
  
  Int_t iResult = 0;

  Int_t sm = block->GetSubDetector();
  if ( sm == 6 ) sm = 7;
  
  Float_t phi   = ( sm + 0.5 ) * TMath::Pi() / 9.0;  
  Float_t cos   = TMath::Cos( phi );
  Float_t sin   = TMath::Sin( phi );
  
  Byte_t* ptrData = reinterpret_cast<Byte_t*>(block->GetData());
  UInt_t ptrSize = block->GetSize();

  for (UInt_t size = 0; size+sizeof(AliHLTTRDCluster) <= ptrSize; size+=sizeof(AliHLTTRDCluster) ) {
    AliHLTTRDCluster *cluster = reinterpret_cast<AliHLTTRDCluster*>(&(ptrData[size]));
   
    AliTRDcluster *trdCluster = new AliTRDcluster;
    cluster->ExportTRDCluster( trdCluster );
   
    contCol->Fill(cos*trdCluster->GetX() - sin*trdCluster->GetY(), 
		   sin*trdCluster->GetX() + cos*trdCluster->GetY(), 
		   trdCluster->GetZ(),
		   trdCluster->GetQ() );    
     
    cont->SetNextPoint(cos*trdCluster->GetX() - sin*trdCluster->GetY(), 
    		       sin*trdCluster->GetX() + cos*trdCluster->GetY(), trdCluster->GetZ());
  }
  
  return iResult;
}

// -----------------------------------------------------------------
Int_t processTRDHistograms(AliHLTHOMERBlockDesc* block, TCanvas * canvas) {

  Int_t iResult = 0;

  if ( ! block->GetClassName().CompareTo("TH1D")) {
    TH1D* histo = reinterpret_cast<TH1D*>(block->GetTObject());
    ++gTRDHistoCount;
  
    TVirtualPad* pad = canvas->cd(gTRDHistoCount);
    histo->Draw();
    pad->SetGridy();
    pad->SetGridx();

    if ( ! strcmp(histo->GetName(), "nscls") ) {
      gTRDEvents = static_cast<Int_t>(histo->GetEntries());
	 histo->GetXaxis()->SetRangeUser(0.,15.);
    }

    if ( ! strcmp(histo->GetName(),"sclsdist") ||
	 ! strcmp(histo->GetName(),"qClsCand") )
      pad->SetLogy();
  }
  else if ( ! block->GetClassName().CompareTo("TH2F")) {
    TH2F *hista = reinterpret_cast<TH2F*>(block->GetTObject());
    ++gTRDHistoCount;
    
    TVirtualPad* pad = canvas->cd(gTRDHistoCount);

    if (gTRDEvents > 0)
      hista->Scale(1./gTRDEvents);

    hista->Draw("COLZ");
    pad->SetLogz();
    pad->SetGridy();
    pad->SetGridx();
  }

  return iResult;
}

//****************************************************************************

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
void loopEventFast() {
  eventTimerFast.SetCommand("nextEvent()");
  eventTimerFast.Start(500);
}

// -----------------------------------------------------------------
void stopLoopEventFast() {
  eventTimerFast.Stop();
}

// -----------------------------------------------------------------
void EventLoopFast() {
  
  // Start/stop event loop
  if ( gEventLoopStarted ) {
    loopEventFast();
    gEventLoopStarted = kTRUE;
  }
  else {
    stopLoopEventFast();
    gEventLoopStarted = kFALSE;
  }
}
