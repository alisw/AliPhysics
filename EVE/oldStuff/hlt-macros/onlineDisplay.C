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
#include "TLine.h"
#include "TEveStraightLineSet.h"
#include "TGeoNode.h"
#include "TGeoManager.h"
#include "TGLOverlayButton.h"

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

//****************** AliRoot/MUON **********************************
#include "AliMUONCalibrationData.h"
#include "AliMUONVCalibParam.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMUONConstants.h"

#include "AliMpCDB.h"
#include "AliMpDDLStore.h"
#include "AliMpDetElement.h"
#include "AliMpDEIterator.h"
#include "AliMpVSegmentation.h"
#include "AliMpSegmentation.h"
#include "AliMpTriggerCrate.h"
#include "AliMpLocalBoard.h"

//****************** AliRoot/ITS ***********************************
#include "AliITSRecPoint.h"

//****************** AliRoot/TRD ***********************************
#include "AliHLTTRDCluster.h"
#include "AliTRDcluster.h"
#include "AliTRDCalibraVdriftLinearFit.h"

//****************** AliRoot/EMCAL *********************************
#include "AliHLTCaloClusterDataStruct.h"
#include "AliHLTCaloClusterReader.h"

#include "HLT/CALO/AliHLTCaloChannelDataHeaderStruct.h"
#include "HLT/CALO/AliHLTCaloChannelDataStruct.h"

//****************** AliRoot/PHOS **********************************
#include "AliPHOSGeometry.h"
#include "HLT/PHOS/AliHLTPHOSDigitDataStruct.h"
#include  "AliHLTPHOSChannelDataHeaderStruct.h"
#include  "AliHLTPHOSChannelDataStruct.h"
#include  "HLT/CALO/AliHLTCaloChannelDataHeaderStruct.h"
#include  "HLT/CALO/AliHLTCaloChannelDataStruct.h"
#include "AliHLTCaloDigitDataStruct.h"
//****************** Macros ****************************************
#include "hlt_structs.C"
#include "hlt_alieve_init.C"
#include "geom_gentle_hlt.C"
#include "macros/esd_tracks.C"
#include "hlt_esd_tracks.C"

//#include "alieve_vizdb.C"

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
TEvePointSet*                             gTPCTestClusters       = 0;
TEvePointSetArray*                        gTPCColClusters    = 0;
TEveBoxSet*                               gPHOSBoxSet[5]     = {0, 0, 0, 0, 0}; 
TEveBoxSet*                               gEMCALBoxSet[13]   = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
TEvePointSet*                             gMUONClusters      = 0;
TEveStraightLineSet*                      gMUONTracks        = 0;

// -- Text output members
TEveText*                                 gHLTText           = 0;

// -- Tracks members
TEveTrackList*                            gTPCTrack          = 0;

// -- Canvas for histograms
TCanvas*                                  gTRDCanvas         = 0;
TCanvas*                                  gTPCCanvas         = 0;
TCanvas*                                  gTPCClustCanvas          = 0;
TCanvas*                                  gTRDCalibCanvas    = 0;
TCanvas*                                  gTRDEORCanvas      = 0;
TCanvas*                                  gPrimVertexCanvas  = 0;
TCanvas*                                  gSPDVertexCanvas   = 0;
TCanvas*                                  gITSCanvas         = 0;
TCanvas*                                  gSSDCanvas0        = 0;
TCanvas*                                  gSSDCanvas1        = 0;
TCanvas*                                  gV0Canvas          = 0;
TCanvas*                                  gPHOSCanvas          = NULL;
TCanvas*                                  gEMCALCanvas          = 0;

// -- vertex --
Int_t                                     gSPDVertexHistoCount  = 0;


// -- ITS --
Int_t                                     gITSHistoCount     = 0;
Int_t                                     gSSDHistoCount0    = 0;
Int_t                                     gSSDHistoCount1    = 0;

// -- TRD --
Int_t                                     gTRDHistoCount     = 0;
Int_t                                     gTRDEvents         = 0;
Int_t                                     gTRDBins           = 12;

// -- TPC --
Int_t                                     gTPCBins           = 15;
TH1F*                                     gTPCCharge         = 0;
TH1F*                                     gTPCQMax           = 0;
TH1F*                                     gTPCQMaxOverCharge = 0;

TH1F*                                     gTPCPt        = 0; // KK
TH1F*                                     gTPCEta       = 0; 
TH1F*                                     gTPCPsi       = 0; 
TH1F*                                     gTPCnClusters = 0; 
TH1F*                                     gTPCMult      = 0;

// -- PHOS --
TEveElementList*                          gPHOSElementList   = 0;
Int_t                                     gPHOSHistoCount    =0;
// -- EMCAL
TEveElementList*                          gEMCALElementList  = 0;
TGeoNode*                                 gEMCALNode         = 0;
Int_t                                     gEMCALHistoCount    =0;

// --- Flag if eventloop is running
Bool_t                                    gEventLoopStarted = kFALSE;


//Container for gGeoManager till it is broken
TGeoManager *fGeoManager = 0;
// -----------------------------------------------------------------
// --                          Methods                            --
// -----------------------------------------------------------------

Int_t initializeEveViewer( Bool_t TPCMode, Bool_t MUONMode, Bool_t TRDMode );

void nextEvent();

Int_t processEvent();

Int_t processPHOSClusters( AliHLTHOMERBlockDesc* block);
Int_t processPHOSDigits(AliHLTHOMERBlockDesc* block);

Int_t processEMCALClusters( AliHLTHOMERBlockDesc* block);

Int_t processEsdTracks( AliHLTHOMERBlockDesc* block, TEveTrackList* cont );

Int_t processHLTRDLST( AliHLTHOMERBlockDesc* block );

Int_t processROOTTOBJ( AliHLTHOMERBlockDesc* block, TEveText* text );

Int_t processTPCClusters( AliHLTHOMERBlockDesc * block, TEvePointSet *cont, TEvePointSetArray *contCol = NULL );

Int_t processTRDClusters( AliHLTHOMERBlockDesc * block, TEvePointSet *cont, TEvePointSetArray *contCol);

Int_t processTRDHistograms( AliHLTHOMERBlockDesc * block, TCanvas * canvas );

Int_t processITSHistograms( AliHLTHOMERBlockDesc * block, TCanvas * canvas );

Int_t processSSDHistograms( AliHLTHOMERBlockDesc * block, TCanvas * canvas, TCanvas * canvas1 );

Int_t processPrimVertexHistograms( AliHLTHOMERBlockDesc * block, TCanvas * canvas );

Int_t processSPDVertexHistograms( AliHLTHOMERBlockDesc * block, TCanvas * canvas );

Int_t processV0Histograms( AliHLTHOMERBlockDesc * block, TCanvas * canvas );

Int_t processTRDCalibHistograms( AliHLTHOMERBlockDesc *block, TCanvas *canvas );

Int_t processMUONClusters( AliHLTHOMERBlockDesc* block );

Int_t processMUONTracks( AliHLTHOMERBlockDesc* block );

Int_t processITSClusters( AliHLTHOMERBlockDesc* block, TEvePointSet* cont );

Int_t processITSHist( AliHLTHOMERBlockDesc* block );

//Global functions

Int_t processBlock( AliHLTHOMERBlockDesc* block );

void resetDisplay();

//Detector specific block handlers
Int_t processTRDBlock(AliHLTHOMERBlockDesc * block);
void processISPDBlock(AliHLTHOMERBlockDesc * block);
void processISDDBlock(AliHLTHOMERBlockDesc * block);
void processISSDBlock(AliHLTHOMERBlockDesc * block);
void processMUONBlock(AliHLTHOMERBlockDesc * block);
Int_t processPHOSBlock (AliHLTHOMERBlockDesc * block);
TEveElementList * createPHOSElementList();
void processEMCALBlock(AliHLTHOMERBlockDesc * block);
TEveElementList * createEMCALElementList();
void processHLTBlock(AliHLTHOMERBlockDesc * block);
void processTPCBlock(AliHLTHOMERBlockDesc * block);
void processITSBlock(AliHLTHOMERBlockDesc * block);
//Generic functions

TCanvas * CreateCanvas(TString  tabTitle, TString  canvasTitle );

Int_t updateDisplay();

Int_t AddHistogramsToCanvas(AliHLTHOMERBlockDesc * block, TCanvas * canvas, Int_t &cdCount );

void writeToFile();


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
  fGeoManager = gGeoManager;

  // -- Create new hM object
  // -------------------------
  gHomerManager = new AliEveHOMERManager();
  gHomerManager->SetRetryCount(100,10);

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
  gGeoManager = fGeoManager;
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
  
  //=============================================================================
  // Visualization database
  //============================================================================

  TEveUtil::AssertMacro("VizDB_scan.C");
  
  //  alieve_vizdb();

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

  gGeoManager = fGeoManager;

  gEMCALNode = gGeoManager->GetTopVolume()->FindNode("XEN1_1");

  TEveGeoTopNode* emcal_re = new TEveGeoTopNode(gGeoManager, gEMCALNode);
  gEve->AddGlobalElement(emcal_re);
  gEve->Redraw3D();

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


   
  //Add HLT Text to windows
 
  TGLOverlayButton *ob1 = new TGLOverlayButton(g3DView->GetGLViewer(),  "HLT", 0, 20, 110, 60);
  ob1->SetAlphaValues(0.8, 0.8);
  //  cout << "color" << ob1->GetBackColor() << endl;
  //ob1->SetBackColor(8421631);
  ob1->SetBackColor(10492431);
  TGLOverlayButton *ob2 = new TGLOverlayButton(g3DView->GetGLViewer(),  "ALICE", 0, 0, 110, 20);
  ob2->SetAlphaValues(0.8, 0.8);
  ob1->SetBackColor(0.2);
  TGLOverlayButton *ob3 = new TGLOverlayButton(gEve->GetDefaultGLViewer(),  "HLT", 0, 20, 110, 60);
  ob3->SetAlphaValues(0.8, 0.8);
  TGLOverlayButton *ob4 = new TGLOverlayButton(gEve->GetDefaultGLViewer(),  "ALICE", 0, 0, 110, 20);
  ob4->SetAlphaValues(0.8, 0.8);


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


  slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  slot->StartEmbedding();  
  gPrimVertexCanvas = new TCanvas("canvasPrimVertex","canvasPrimVertex", 600, 400);
  slot->StopEmbedding("Primary Vertex");

  slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  slot->StartEmbedding();  
  gSPDVertexCanvas = new TCanvas("canvasSPDVertex","canvasSPDVertex", 600, 400);
  slot->StopEmbedding("SPD Vertex");

  slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  slot->StartEmbedding();  
  gV0Canvas = new TCanvas("canvasV0","canvasV0", 600, 400);
  slot->StopEmbedding("V0");


  slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  slot->StartEmbedding();
  gSSDCanvas0 = new TCanvas("canvasSSD0","canvasSSD0", 600, 400);
  slot->StopEmbedding("SSD QA 0");

  slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  slot->StartEmbedding();
  gSSDCanvas1 = new TCanvas("canvasSSD1","canvasSSD1", 600, 400);
  slot->StopEmbedding("SSD QA 1");

  slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  slot->StartEmbedding();
  gTPCClustCanvas = new TCanvas("canvasTPCClust","canvasTPCClust", 600, 400);
  gTPCClustCanvas->Divide(3, 2);
  slot->StopEmbedding("TPC Cluster QA Histograms ");


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

  if ( gHomerManager->NextEvent() ) {
    if (gEventLoopStarted) {
      cout << "HomerManager failed getting next event, trying to reconnect" << endl;

      gHomerManager->DisconnectHOMER();
      gHomerManager->ConnectEVEtoHOMER();
      nextEvent();
   
    } else {
      return;
    }
  }

  processEvent();
}

TCanvas * CreateCanvas(TString  tabTitle, TString  canvasTitle ) {

  

//   TEveWindowSlot *slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
//   slot->StartEmbedding();
//   TCanvas * canvas = new TCanvas(canvasTitle.Data(),canvasTitle.Data(), 600, 400);
//   slot->StopEmbedding(tabTitle.Data());
 
  cout << "Create"<<endl;

  TCanvas * canvas = gEve->AddCanvasTab(canvasTitle.Data());
  return canvas;
}


// -----------------------------------------------------------------
Int_t processEvent() {

  Int_t iResult = 0;

  gStyle->SetPalette(1, 0);
  gEve->DisableRedraw();

  resetDisplay();
  
  //==============================================================================
  // -- Process Blocks
  //==============================================================================

  AliHLTHOMERBlockDesc* block = 0;

  if ( gHomerManager->GetBlockList() == NULL) {
    printf ("onlineDisplay:   No regular BlockList ... \n");
    cout << endl;
    //return -1;
    
  } else {
    
    if (gHomerManager->GetBlockList()->IsEmpty() ) {
      printf ("onlineDisplay:   No Sync Blocks in list ... \n");
      cout<<endl;
      return -2;
    }  
    
    
    TIter next(gHomerManager->GetBlockList());
    
    while ((block = (AliHLTHOMERBlockDesc*)next())) {
      iResult = processBlock(block);
    } 
  }

  //Read out histograms and elements from detectors outside physics 1 partition
  TIter anext(gHomerManager->GetAsyncBlockList());
  
  while ( (block = (AliHLTHOMERBlockDesc*)anext()) ) {
    iResult = processBlock(block);
  }
  
  updateDisplay();

  return iResult;
}


void resetDisplay() {


  //==============================================================================
  // -- Reset
  //==============================================================================
  
  if ( gTRDCanvas ) {
    gTRDCanvas->Clear();
    gTRDCanvas->Divide(3,2);
  }
  if(gTRDCalibCanvas){
    gTRDCalibCanvas->Clear();
    gTRDCalibCanvas->Divide(2,2);
  }
  if(gTRDEORCanvas){
    gTRDEORCanvas->Clear();
    gTRDEORCanvas->Divide(3,2);
  }
  
  if(gSSDCanvas0) {
    gSSDCanvas0->Clear();
    gSSDCanvas0->Divide(6,6);
  }

  if(gSSDCanvas1) {
    gSSDCanvas1->Clear();
    gSSDCanvas1->Divide(2,3);
  }
  
  if(gITSCanvas) {
    gITSCanvas->Clear();
    gITSCanvas->Divide(3,3);
  }

  if(gPrimVertexCanvas) {
    gPrimVertexCanvas->Clear();
    gPrimVertexCanvas->Divide(2,2);
  }

  if(gSPDVertexCanvas) {
    gSPDVertexCanvas->Clear();
    gSPDVertexCanvas->Divide(2,2);
  }

  if(gV0Canvas) {
    gV0Canvas->Clear();
    gV0Canvas->Divide(2,2);
  }

  if(gPHOSCanvas) {
    gPHOSCanvas->Clear();
    gPHOSCanvas->Divide(3, 3);
  }

  if(gEMCALCanvas) {
    gEMCALCanvas->Clear();
    gEMCALCanvas->Divide(3, 3);
  }

  if ( gTPCTrack )     gTPCTrack->DestroyElements();

  if ( gSPDClusters )  gSPDClusters->Reset();
  if ( gSSDClusters )  gSSDClusters->Reset();
  if ( gSDDClusters )  gSDDClusters->Reset();
  if ( gTPCClusters )  gTPCClusters->Reset();
  if ( gTPCTestClusters )  gTPCTestClusters->Reset();
  if ( gTRDClusters )  gTRDClusters->Reset();
  if ( gMUONClusters ) gMUONClusters->Reset();
  if ( gMUONTracks ){
    gMUONTracks->Destroy();
    gMUONTracks = 0x0;
  }

  if ( gPHOSBoxSet[1] )
    for(int im = 0; im < 5; im++)
      gPHOSBoxSet[im]->Reset();   

  if ( gEMCALElementList)
    for(int i = 0; i < 12; i++) {
      gEMCALBoxSet[i]->Reset();
    }
  
  if ( gTPCColClusters )
    for (Int_t ii = 0; ii <= gTPCBins+1; ++ii) 
      gTPCColClusters->GetBin(ii)->Reset();
  
  if ( gTRDColClusters )
    for (Int_t ii = 0; ii <= gTRDBins+1; ++ii) 
      gTRDColClusters->GetBin(ii)->Reset();

  gTRDHistoCount = 0;
  gSPDVertexHistoCount = 0;
  gSSDHistoCount0 = 0;
  gSSDHistoCount1 = 0;
  gITSHistoCount = 0;
  gPHOSHistoCount = 0;
  gEMCALHistoCount = 0;

  
}


Int_t updateDisplay() {
  Int_t iResult = 0;

  
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
  if ( gTPCTestClusters ) gTPCTestClusters->ResetBBox();
  if ( gTRDClusters ) gTRDClusters->ResetBBox();
  if ( gSPDClusters ) gSPDClusters->ResetBBox();
  if ( gSDDClusters ) gSDDClusters->ResetBBox();
  if ( gSSDClusters ) gSSDClusters->ResetBBox();
  if ( gMUONClusters ) gMUONClusters->ResetBBox();

  if ( gPHOSBoxSet[1] )
    for(int im = 0; im < 5; im++)
      gPHOSBoxSet[im]->ResetBBox();      

  if ( gEMCALElementList )
    for(int sm = 0; sm < 12; sm++) 
      gEMCALBoxSet[sm]->ResetBBox();
    
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
  
  esd_track_propagator_setup(cont->GetPropagator(),-0.1*esd->GetMagneticField(), 520);

  //  printf( "Number of ESD Tracks : %d \n", esd->GetNumberOfTracks());

  for (Int_t iter = 0; iter < esd->GetNumberOfTracks(); ++iter) {
    //AliEveTrack* track = dynamic_cast<AliEveTrack*>(esd_make_track(esd->GetTrack(iter), cont));
    AliEveTrack* track = dynamic_cast<AliEveTrack*>(hlt_esd_make_track(esd->GetTrack(iter), cont));
    cont->AddElement(track);
    
//     gTPCPt->Fill(esd->GetTrack(iter)->GetSignedPt()); // KK
//     gTPCEta->Fill(esd->GetTrack(iter)->GetSnp());
//     gTPCPsi->Fill(esd->GetTrack(iter)->GetTgl());
//     gTPCnClusters->Fill(esd->GetTrack(iter)->GetTPCNcls());  
  }
  
//   gTPCMult->Fill(esd->GetNumberOfTracks()); // KK
  
//   Int_t icd = 0;
//   gTPCClustCanvas->Clear();
//   gTPCClustCanvas->Divide(2, 2);
//   gTPCClustCanvas->cd(icd++);
//   gTPCPt->Draw();
//   gTPCClustCanvas->cd(icd++);
//   gTPCEta->Draw();
//   gTPCClustCanvas->cd(icd++);
//   gTPCPsi->Draw();
//   gTPCClustCanvas->cd(icd++);
//   gTPCnClusters->Draw();
//   gTPCClustCanvas->cd(icd++);
//   gTPCMult->Draw();
//   gTPCClustCanvas->Update();

  
  cont->SetTitle(Form("N=%d", esd->GetNumberOfTracks()) );
  cont->MakeTracks();
  
  return 0;
}

// // -----------------------------------------------------------------
// Int_t processPHOSClusters(AliHLTHOMERBlockDesc* block) {


//   cout <<"Processgin phos clusters"<<endl;
   
//   AliHLTCaloChannelDataHeaderStruct *chh = reinterpret_cast<AliHLTCaloChannelDataHeaderStruct*> (block->GetData());
  
//   AliHLTCaloChannelDataStruct *chd = reinterpret_cast<AliHLTCaloChannelDataStruct*>(chh+1);

//   for(Int_t i = 0; i < chh->fNChannels; i++, chd++) {
    
//     Int_t gain = (chd->fChannelID >> 12)&0x1;    
//     Int_t module = (chd->fChannelID >> 13)&0x1f;
//     module = 4 -module;
   

//     if(gain == 1)
      
//       {
// 	Float_t x = (static_cast<Float_t>(chd->fChannelID&0x3f) - 32)* 2.2;
// 	Float_t z = (static_cast<Float_t>((chd->fChannelID >> 6)&0x3f) - 28) * 2.2;
// 	//	gPHOSBoxSet[ds->fModule]->AddBox(ds->fLocX, 0, ds->fLocZ, 2.2, ds->fEnergy*20, 2.2);
// 	gPHOSBoxSet[module]->AddBox(x, 0, z, 2.2, chd->fEnergy/1, 2.2);
// 	gPHOSBoxSet[module]->DigitValue(static_cast<Int_t>(chd->fEnergy));
//       }



//   }
// cout <<"done with phos"<<endl;


//   return 0;
// }



// // -----------------------------------------------------------------
Int_t processPHOSDigits(AliHLTHOMERBlockDesc* block) {

  AliHLTCaloDigitDataStruct *ds = reinterpret_cast<AliHLTCaloDigitDataStruct*> (block->GetData());
  UInt_t nDigits = block->GetSize()/sizeof(AliHLTCaloDigitDataStruct);
    

  for(UInt_t i = 0; i < nDigits; i++, ds++) {

    Float_t x = (ds->fX - 32)* 2.2;
      Float_t z = (ds->fZ - 28) * 2.2;

      cout << "MODULE DIGITTYP  :" << ds->fModule;

    gPHOSBoxSet[4-ds->fModule]->AddBox(x, 0, z, 2.2, ds->fEnergy*20, 2.2);
    gPHOSBoxSet[4-ds->fModule]->DigitValue(static_cast<Int_t>(ds->fEnergy*10));
  }

  return 0;
}


// -----------------------------------------------------------------
Int_t processPHOSClusters(AliHLTHOMERBlockDesc* block) {

  AliHLTCaloClusterHeaderStruct *dh = reinterpret_cast<AliHLTCaloClusterHeaderStruct*> (block->GetData());
  AliHLTCaloClusterReader * clusterReader = new AliHLTCaloClusterReader();
  clusterReader->SetMemory(dh);  

  AliHLTCaloClusterDataStruct * ds;

  Int_t nc = 0;

  
  while( (ds = clusterReader->NextCluster()) ){
    nc++;
    cout << "Cluster nr " << nc  << "   Energy  " << ds->fEnergy<< "      " << ds->fGlobalPos[0] << " " << ds->fGlobalPos[1] << " " << ds->fGlobalPos[2] << endl;


    gPHOSBoxSet[2]->AddBox(ds->fGlobalPos[0], 0, ds->fGlobalPos[2], 2.2, -ds->fEnergy*20, 2.2);
    gPHOSBoxSet[2]->DigitValue(static_cast<Int_t>(ds->fEnergy*10));
  }

  return 0;
}



// -----------------------------------------------------------------
Int_t processEMCALClusters(AliHLTHOMERBlockDesc* block) {

  cout << "EMCAL"<<endl;

  AliHLTCaloChannelDataHeaderStruct *dhs = reinterpret_cast<AliHLTCaloChannelDataHeaderStruct*> (block->GetData());
  Short_t nC = dhs->fNChannels;
  AliHLTUInt8_t *ui = reinterpret_cast<AliHLTUInt8_t*>(dhs) + sizeof(AliHLTCaloChannelDataHeaderStruct);
  AliHLTCaloChannelDataStruct *ds = reinterpret_cast<AliHLTCaloChannelDataStruct*>(ui);
  
  UShort_t fX =0;
  UShort_t fZ =0;
  UShort_t fGain =0;
  UShort_t fModuleId =0;
  


  for(Short_t s = 0; s < nC; s++ ) {
      
    fX = ds->fChannelID&0x3f;
    fZ = (ds->fChannelID >> 6)&0x3f;
    fGain = (ds->fChannelID >> 12)&0x1;
    fModuleId  = (ds->fChannelID >> 13)&0x1f;
    
    //    cout << fX << " " << fZ << " " << fGain << " " << fModuleId <<endl;
   
    if ( (  fModuleId < 12 ) ) {
      gEMCALBoxSet[fModuleId]->AddBox(10, fX*6-12*6, fZ*6-24*6, ds->fEnergy/4, 6, 6);
      gEMCALBoxSet[fModuleId]->DigitValue(static_cast<Int_t>(ds->fEnergy));
    } 
    ds++;
    
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


  if(!gTPCCanvas){
    gTPCCanvas = gEve->AddCanvasTab("TPC QA");
    gTPCCharge = new TH1F("ClusterCharge","ClusterCharge",100,0,500);
    gTPCQMax = new TH1F("QMax","QMax",50,0,250);
    gTPCQMaxOverCharge = new TH1F("QMaxOverCharge","QMaxOverCharge",50,0,1);
  }



  





  if ( cd->fSpacePointCnt != 0 ) {
    for (Int_t iter = 0; iter < cd->fSpacePointCnt; ++iter, data += sizeof(AliHLTTPCSpacePointData)) {
      AliHLTTPCSpacePointData *sp = reinterpret_cast<AliHLTTPCSpacePointData*> (data);
      cont->SetNextPoint(cos*sp->fX - sin*sp->fY, sin*sp->fX + cos*sp->fY, sp->fZ);
      if (contCol)
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
// 	  cout << setw(10) << left << ch + 11 << setw(0);
// 	  cout << setw(13) << left << hit->fX << setw(0);
// 	  cout << setw(13) << left << hit->fY << setw(0);
// 	  cout << hit->fZ << setw(0) << endl;
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
Int_t processMUONTracks(AliHLTHOMERBlockDesc* block) {
  
  Int_t iResult = 0;
  
  unsigned long size = block->GetSize();
  Int_t * buffer = (Int_t *)block->GetData();
  AliHLTMUONRecHitStruct hit1,hit2;
  hit1.fX = hit1.fY = hit1.fZ = hit2.fX = hit2.fY = hit2.fZ = 0;
  Int_t ch1=0, ch2=0;
  Float_t x0=0.0,y0=0.0,z0=0.0;
  Float_t x3=0.0,y3=0.0,z3=0.0;
  if(block->GetDataType().CompareTo("MANTRACK") == 0){  
    AliHLTMUONMansoTracksBlockReader mantrackblock(buffer, size);
    const AliHLTMUONMansoTrackStruct* mtrack = mantrackblock.GetArray();
    for(AliHLTUInt32_t ientry = 0; ientry < mantrackblock.Nentries(); ientry++){
      const AliHLTMUONRecHitStruct* hit = &mtrack->fHit[0];
      for(AliHLTUInt32_t ch = 0; ch < 4; ch++){
	// cout << setw(10) << left << ch + 7 << setw(0);
	// cout << setw(13) << left << hit->fX << setw(0);
	// cout << setw(13) << left << hit->fY << setw(0);
	// cout << hit->fZ << setw(0) << endl;
	if(hit->fZ != 0.0){
	  if(ch==0 || ch==1){
	    hit1 = *hit; ch1 = ch+6;
	  }else{
	    hit2 = *hit; ch2 = ch+6;
	  }
	}
	hit++;
      }// trig chamber loop
      // printf("ch : %d, (X,Y,Z) : (%f,%f,%f)\n",ch1,hit1.fX,hit1.fY,hit1.fZ);
      // printf("ch : %d, (X,Y,Z) : (%f,%f,%f)\n",ch2,hit2.fX,hit2.fY,hit2.fZ);
      // meminfo();
      z3 = AliMUONConstants::DefaultChamberZ(ch2+4);
      y3 =  hit1.fY - (hit1.fZ-z3)*(hit1.fY - hit2.fY)/(hit1.fZ - hit2.fZ) ;
      x3 =  hit1.fX - (hit1.fZ-z3)*(hit1.fX - hit2.fX)/(hit1.fZ - hit2.fZ) ;

      z0 = AliMUONConstants::DefaultChamberZ(ch1);
      y0 =  hit1.fY - (hit1.fZ-z0)*(hit1.fY - hit2.fY)/(hit1.fZ - hit2.fZ) ;
      x0 =  hit1.fX - (hit1.fZ-z0)*(hit1.fX - hit2.fX)/(hit1.fZ - hit2.fZ) ;
      

      gMUONTracks->AddLine(x0,y0,z0,x3,y3,z3);
      mtrack++;
    }
    cout<<"NofManso Tracks : "<<mantrackblock.Nentries()<<endl;
  }
  
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

  if ( ! block->GetClassName().CompareTo("TH1F")) {
    TH1F* histo = reinterpret_cast<TH1F*>(block->GetTObject());
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
	 ! strcmp(histo->GetName(),"evSize") )
      pad->SetLogy();
  }

  gTRDCanvas->Update();

  return iResult;
}

// -----------------------------------------------------------------
Int_t processPrimVertexHistograms(AliHLTHOMERBlockDesc* block, TCanvas * canvas) {

  Int_t iResult = 0;

  if ( ! block->GetClassName().CompareTo("TH1F")) {
    TH1F* histo = reinterpret_cast<TH1F*>(block->GetTObject());
    if( histo ){
      TString name(histo->GetName());
      if( !name.CompareTo("primVertexZ") ){
	canvas->cd(2);
	histo->Draw();
      }else if( !name.CompareTo("primVertexX") ){
	canvas->cd(3);
	histo->Draw();
      }else if( !name.CompareTo("primVertexY") ){
	canvas->cd(4);
	histo->Draw();
      }
    }
  }  else if ( ! block->GetClassName().CompareTo("TH2F")) {
    TH2F *hista = reinterpret_cast<TH2F*>(block->GetTObject());
    if (hista ){
       TString name(hista->GetName());
       if( !name.CompareTo("primVertexXY")) {      
	 canvas->cd(1);
	 hista->Draw();
       }
    }
  }
  canvas->cd();

  return iResult;
}

// -----------------------------------------------------------------
Int_t processSPDVertexHistograms(AliHLTHOMERBlockDesc* block, TCanvas * canvas) {

  Int_t iResult = 0;

  if ( ! block->GetClassName().CompareTo("TH1F")) {
    TH1F* histo = reinterpret_cast<TH1F*>(block->GetTObject());
    ++gSPDVertexHistoCount;
  
    canvas->cd(gSPDVertexHistoCount);
    histo->Draw();

  }  
  else if ( ! block->GetClassName().CompareTo("TH2F")) {
    TH2F *hista = reinterpret_cast<TH2F*>(block->GetTObject());
    if (hista) {
      ++gSPDVertexHistoCount;
  
      canvas->cd(gSPDVertexHistoCount);
      hista->Draw();
    }
  }
  canvas->cd();

  return iResult;
}

// -----------------------------------------------------------------
Int_t processV0Histograms(AliHLTHOMERBlockDesc* block, TCanvas * canvas) {

  cout << "Processing to see if it's V0 histogram, !!!!!!!!!"<<endl;

  Int_t iResult = 0;
  bool update = 0;
  if ( ! block->GetClassName().CompareTo("TH1F")) {
    TH1F* histo = reinterpret_cast<TH1F*>(block->GetTObject());
    if( histo ){
      TString name(histo->GetName());
      if( !name.CompareTo("hKShort") ){
	canvas->cd(1);
	histo->Draw();
	update = 1;
      }else if( !name.CompareTo("hLambda") ){
	canvas->cd(3);
	histo->Draw();
	update = 1;
      }
    }
  }  else if ( ! block->GetClassName().CompareTo("TH2F")) {
    TH2F *hista = reinterpret_cast<TH2F*>(block->GetTObject());
    if (hista ){
       TString name(hista->GetName());
       if( !name.CompareTo("hAP")) {      
	 canvas->cd(2);
	 hista->Draw();
	 update = 1;
       }
       else if( !name.CompareTo("hGammaXY")) {      
	 canvas->cd(4);
	 hista->Draw();
	 update = 1;
       }
    }
  }
  if( update ){
    canvas->cd();
    canvas->Update();
  }
  return iResult;
}


// -----------------------------------------------------------------
Int_t processITSHistograms(AliHLTHOMERBlockDesc* block, TCanvas * canvas) {

  Int_t iResult = 0;

  if ( ! block->GetClassName().CompareTo("TH1F")) {
    TH1F* histo = reinterpret_cast<TH1F*>(block->GetTObject());
    ++gITSHistoCount;
  
    canvas->cd(gITSHistoCount);
    histo->Draw();

  } 
  else if ( ! block->GetClassName().CompareTo("TH2F")) {
    TH2F *hista = reinterpret_cast<TH2F*>(block->GetTObject());
    if (hista) {
      ++gITSHistoCount;
      
      canvas->cd(gITSHistoCount);
      hista->Draw("COLZ");
    }
  }
  canvas->cd();

  return iResult;
}


Int_t AddHistogramsToCanvas(AliHLTHOMERBlockDesc * block, TCanvas * canvas, Int_t &cdCount ) {

  if (canvas == NULL) 
    cout <<"fucked canvasn"<<endl;
    

  Int_t iResult = 0;

  if ( ! block->GetClassName().CompareTo("TObjArray")) {
    TIter next((TObjArray*)(block->GetTObject()));
    TObject *object;
    while (( object = (TObject*) next())) {
      TH1* histo = dynamic_cast<TH1*>(object);
      cout <<" cd "<< cdCount << endl;
      canvas->cd(++cdCount);
      cout <<"draw"<<endl;
      histo->Draw();
      cout <<"drew"<<endl;
    }
  }

  //  canvas->cd();
  canvas->Update();

  return iResult;
}


// -----------------------------------------------------------------
Int_t processSSDHistograms(AliHLTHOMERBlockDesc* block, TCanvas *canvas0, TCanvas *canvas1) {

  Int_t iResult = 0;

  if ( ! block->GetClassName().CompareTo("TH1F")) {
    TH1F* histo = reinterpret_cast<TH1F*>(block->GetTObject());
    ++gSSDHistoCount0;
  
    canvas0->cd(gSSDHistoCount0);
    histo->Draw();
  }  
  else if ( ! block->GetClassName().CompareTo("TH2F")) {
    TH2F *hista = reinterpret_cast<TH2F*>(block->GetTObject());
    if (hista) {
      ++gSSDHistoCount1;
  
      canvas1->cd(gSSDHistoCount1);
      hista->Draw("COLZ");
    }
  }  else if ( ! block->GetClassName().CompareTo("TObjArray")) {
    TIter next((TObjArray*)(block->GetTObject()));
    TObject *object;
    while (( object = (TObject*) next())) {
      TString string;
      string = "TH1F";
      TString string2;
      string2 = "TH2F";

      if ( !(string.CompareTo(object->ClassName())) ) {
	TH1F* histo = reinterpret_cast<TH1F*>(object);
	++gSSDHistoCount0;
	
	canvas0->cd(gSSDHistoCount0);
	histo->Draw();
	
      } 
      else if ( !(string2.CompareTo(object->ClassName()) ) ) {
	TH2F* histo = reinterpret_cast<TH2F*>(object);
	++gSSDHistoCount1;
	
	canvas1->cd(gSSDHistoCount1);
	histo->Draw("COLZ");
      }
    }
  }

  canvas0->cd();  canvas1->cd();
  
  return iResult;
}


//*-------------------------------------------------------------------------------------- 
Int_t processTRDCalibHistograms(AliHLTHOMERBlockDesc* block, TCanvas * canvas) {
  Int_t iResult = 0;

  TObjArray *HistArray=(TObjArray*)block->GetTObject();
  Int_t nCalibHistos=HistArray->GetEntriesFast();
  for(Int_t CalibHistoCount=0;CalibHistoCount<nCalibHistos;CalibHistoCount++){
    canvas->cd(CalibHistoCount+1);
    
    if(HistArray->At(CalibHistoCount)->InheritsFrom("TH2S")){
	 TH2S *histCalib=(TH2S*)(HistArray->At(CalibHistoCount));
	 histCalib->Draw("colz");
       }
     else if(HistArray->At(CalibHistoCount)->InheritsFrom("TH2")){
      //TH2D *histCalib=dynamic_cast<TH2D*>(HistArray->At(CalibHistoCount));
      TH2D *histCalib=(TH2D*)(HistArray->At(CalibHistoCount));
      histCalib->Draw("lego2");
    }
    else if(HistArray->At(CalibHistoCount)->InheritsFrom("TH1")){
      //TH1D *histCalib=dynamic_cast<TH1D*>(HistArray->At(CalibHistoCount));
      TH1D *histCalib=(TH1D*)(HistArray->At(CalibHistoCount));
      histCalib->Draw();
    }
    else if(HistArray->At(CalibHistoCount)->InheritsFrom("AliTRDCalibraVdriftLinearFit")){
      //TH2S *histCalib = ((dynamic_cast<AliTRDCalibraVdriftLinearFit*>(HistArray->At(CalibHistoCount)))->GetLinearFitterHisto(10,kTRUE));
      TH2S *histCalib =(TH2S*)(((AliTRDCalibraVdriftLinearFit*)HistArray->At(CalibHistoCount))->GetLinearFitterHisto(10,kTRUE));

      histCalib->Draw();
    }
    
   
  }
  
  gTRDCalibCanvas->Update();

 return iResult;
}
//****************************************************************************



void writeToFile(){

  TList * bList = gHomerManager->GetBlockList();
  TFile * file = TFile::Open(Form("Event_0x%016X_ITS.root", gHomerManager->GetEventID()), "RECREATE"); 
  bList->Write("blockList", TObject::kSingleKey);
  file->Close();
  
  bList = gHomerManager->GetAsyncBlockList();
  TFile * afile = TFile::Open(Form("Event_0x%016X_Async.root", gHomerManager->GetEventID()), "RECREATE"); 
  bList->Write("blockList", TObject::kSingleKey);
  afile->Close();


//   TIter next(bList);
  
//   AliHLTHOMERBlockDesc* block = 0;
  
//   // -- Iterate over blocks in the block list
//   // ------------------------------------------
//   while ((block = (AliHLTHOMERBlockDesc*)next())) {
//     cout << "balle " << endl;
//   }
  
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
  
  } else {  
    stopLoopEventFast();
    gEventLoopStarted = kFALSE;
  
  }
}



Int_t processBlock (AliHLTHOMERBlockDesc * block ){

  
#if 1//DEBUG
   printf( "------------------- xxxxxxxxxxxxxxx ----------------------\n");
   printf( "Detector           : %s\n", block->GetDetector().Data() );
   printf( "Datatype           : %s\n", block->GetDataType().Data() );
   if (block->IsTObject() )
     printf( "Is TObject of class: %s\n", block->GetClassName().Data() );
   printf( "------------------- xxxxxxxxxxxxxxx ----------------------\n");
#endif

  
  Int_t iResult = 0;

  if ( ! block->GetDetector().CompareTo("HLT") ) 
    processHLTBlock(block);
  
  else if ( ! block->GetDetector().CompareTo("TPC") ) 
    processTPCBlock(block);
  
  else if ( ! block->GetDetector().CompareTo("TRD") ) 
    iResult = processTRDBlock(block);

  else if ( ! block->GetDetector().CompareTo("MUON") && gShowMUON ) 
    processMUONBlock(block);

  else if ( ! block->GetDetector().CompareTo("ISPD") )
    processISPDBlock(block);

  else if ( ! block->GetDetector().CompareTo("ISDD") )
    processISDDBlock(block);

  else if ( ! block->GetDetector().CompareTo("ISSD") )
    processISSDBlock(block);

  else if ( ! block->GetDetector().CompareTo("PHOS") ) 
    processPHOSBlock(block);

  else if ( ! block->GetDetector().CompareTo("EMCA") )
    processEMCALBlock(block);

  else if ( ! block->GetDetector().CompareTo("ITS") )
    processITSBlock(block);
 
  

  if ( block->GetDataType().CompareTo("ROOTHIST") == 0 )     
    processV0Histograms( block , gV0Canvas);
   
  return iResult;
} 


void processITSBlock(AliHLTHOMERBlockDesc * block) {
  
  if ( block->GetDataType().CompareTo("ROOTHIST") == 0 ) {
    if(!gITSCanvas){
      gITSCanvas = gEve->AddCanvasTab("ITS QA");
    }
    processITSHistograms( block , gITSCanvas);
    gITSCanvas->Update();
  } 
}


 Int_t processTRDBlock (AliHLTHOMERBlockDesc * block) {

   Int_t iResult = 0;

  if ( ! block->GetDataType().CompareTo("CLUSTERS") ) {
     
    if(!gTRDClusters){
      gTRDClusters = new TEvePointSet("TRD Clusters");
      gTRDClusters->SetMainColor(kBlue);
      gTRDClusters->SetMarkerStyle((Style_t)kFullDotSmall);
      //gEve->AddElement(gTRDClusters);
    } 

    if(!gTRDColClusters){
      gTRDColClusters = new TEvePointSetArray("TRD Clusters Colorized");
      gTRDColClusters->SetMainColor(kRed);
      gTRDColClusters->SetMarkerStyle(4); // antialiased circle
      //	  gTRDColClusters->SetMarkerStyle((Style_t)kFullDotSmall);
      gTRDColClusters->SetMarkerSize(0.4);
      gTRDColClusters->InitBins("Cluster Charge", gTRDBins, 0., gTRDBins*100.);

      //TColor::SetPalette(1, 0); // Spectrum palette
      const Int_t nCol = TColor::GetNumberOfColors();
      for (Int_t ii = 0; ii < gTRDBins+1; ++ii)
	gTRDColClusters->GetBin(ii)->SetMainColor(TColor::GetColorPalette(ii * nCol / (gTRDBins+2)));
	  
      gEve->AddElement(gTRDColClusters);
    } 

    iResult = processTRDClusters( block, gTRDClusters, gTRDColClusters );
    //gTRDClusters->ElementChanged();
    gTRDColClusters->ElementChanged();
  }

  // -- Process TRD Histograms
  else if ( block->GetDataType().CompareTo("ROOTHIST") == 0 ) {
    if(!gTRDCanvas) {
      gTRDCanvas = gEve->AddCanvasTab("TRD");
      gTRDCanvas->Divide(3,2);
    }
    iResult = processTRDHistograms( block, gTRDCanvas );     
  }

  else if(block->GetDataType().CompareTo("CALIBRAH")==0){
     
    if(!gTRDCalibCanvas){
      gTRDCalibCanvas = gEve->AddCanvasTab("TRD Calib");
      gTRDCalibCanvas->Divide(2,2);
    }
     
    iResult=processTRDCalibHistograms(block,gTRDCalibCanvas);
  }

  else if(block->GetDataType().CompareTo("CALIBEOR")==0){
     
    if(!gTRDEORCanvas){
      gTRDEORCanvas = gEve->AddCanvasTab("TRD QA");
      gTRDEORCanvas->Divide(3,2);       
    }
  
    iResult=processTRDCalibHistograms(block,gTRDEORCanvas);
  }
  return iResult;
}

void processMUONBlock(AliHLTHOMERBlockDesc * block){

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
	
   }else if(block->GetDataType().CompareTo("MANTRACK") == 0){
	
     if ( !gMUONTracks ) {
       gMUONTracks = new TEveStraightLineSet("MUON Tracks");
       gMUONTracks->SetMainColor(kRed);
       gMUONTracks->SetLineWidth(3);
       gEve->AddElement(gMUONTracks);
     }

     processMUONTracks( block );
     gMUONTracks->ElementChanged();


   } 
}


void processISPDBlock( AliHLTHOMERBlockDesc * block) {

  if ( block->GetDataType().CompareTo("CLUSTERS") == 0 ) {
    
    if(!gSPDClusters){
      gSPDClusters = new TEvePointSet("SPD Clusters");
      gSPDClusters->SetMainColor(kBlack);
      gSPDClusters->SetMarkerStyle((Style_t)kFullDotMedium);
      gEve->AddElement(gSPDClusters);
    } 
    
    processITSClusters( block , gSPDClusters);
    gSPDClusters->ElementChanged();

  } else if ( block->GetDataType().CompareTo("ROOTHIST") == 0 ) {
    processSPDVertexHistograms( block , gSPDVertexCanvas);
    gSPDVertexCanvas->Update();
  } 
} 

void processISDDBlock(AliHLTHOMERBlockDesc * block){
  
  if ( block->GetDataType().CompareTo("CLUSTERS") == 0 ) {
    
    if(!gSDDClusters){
      gSDDClusters = new TEvePointSet("SDD Clusters");
      gSDDClusters->SetMainColor(kPink);
      gSDDClusters->SetMarkerStyle((Style_t)kFullDotMedium);
      gEve->AddElement(gSDDClusters);
    } 
    
    processITSClusters( block , gSDDClusters);
    gSDDClusters->ElementChanged();
  } 
}

void processISSDBlock(AliHLTHOMERBlockDesc * block){
  if ( block->GetDataType().CompareTo("CLUSTERS") == 0 ) {
    
    if(!gSSDClusters){
      gSSDClusters = new TEvePointSet("SSD Clusters");
      gSSDClusters->SetMainColor(kBlue);
      gSSDClusters->SetMarkerStyle((Style_t)kFullDotMedium);
      gEve->AddElement(gSSDClusters);
    } 
    
    processITSClusters( block , gSSDClusters);
    gSSDClusters->ElementChanged();
  } 
  
  else if ( block->GetDataType().CompareTo("ROOTOBAR") == 0 ) {
    processSSDHistograms( block , gSSDCanvas0, gSSDCanvas1);
    gSSDCanvas0->Update();
    gSSDCanvas1->Update();
  } 
}
 
        

Int_t processPHOSBlock(AliHLTHOMERBlockDesc * block) { 
  
  Int_t iResult = 0;

  if ( block->GetDataType().CompareTo("ROOTHIST") == 0 ) { 
    
    if(!gPHOSCanvas) {
      gPHOSCanvas = CreateCanvas("PHOS QA", "PHOS QA");
      gPHOSCanvas->Divide(3, 2);
    }
    
    AddHistogramsToCanvas(block, gPHOSCanvas, gPHOSHistoCount);
   
 
  } else {
    
    cout <<"PHOS block not histo"<<endl;

    if( !gPHOSElementList )
      gPHOSElementList = createPHOSElementList();
        
    if ( block->GetDataType().CompareTo("CALOCLUS") == 0 )
      iResult = processPHOSClusters( block );
      
    else if ( block->GetDataType().CompareTo("DIGITTYP") == 0 )
      iResult = processPHOSDigits( block);
    
    else if ( block->GetDataType().CompareTo("CHANNELT") == 0 ) 
      iResult = processPHOSClusters( block );
  
  
    //Update the eve elements
    for(int im = 0; im < 5; im++)
      gPHOSBoxSet[im]->ElementChanged();
    
  }
  
  return iResult;
} 


TEveElementList * createPHOSElementList() {
  
  TEveElementList * elementList  = new TEveElementList("PHOS Cells");
  
  gEve->AddElement(elementList);
  
  TVector3 center;
  Float_t angle;
  
  // -- Create boxsets
  for(int im = 0; im < 5; im++) {
    
    TEveRGBAPalette* pal = new TEveRGBAPalette(0,120);
    pal->SetLimits(-0.1, 1024);
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
    
    elementList->AddElement(gPHOSBoxSet[im]);
  }

  return elementList;
}


void processEMCALBlock(AliHLTHOMERBlockDesc * block) {	

  Int_t iResult = 0;
  
  if ( block->GetDataType().CompareTo("ROOTHIST") == 0 ) { 
    
    if(!gEMCALCanvas) {
      gEMCALCanvas = CreateCanvas("EMCAL QA", "EMCAL QA");
      gEMCALCanvas->Divide(3, 2);
    }
    AddHistogramsToCanvas(block, gEMCALCanvas, gEMCALHistoCount);
   
 
  } else {

    if( !gEMCALElementList ){
      gEMCALElementList = createEMCALElementList();
    }

    if ( block->GetDataType().CompareTo("CALOCLUS") == 0 ) {
      iResult = processEMCALClusters( block );
      for(int sm = 0; sm < 12; sm++) {
	gEMCALBoxSet[sm]->ElementChanged();
      }
    }
  } 
}

TEveElementList * createEMCALElementList() { 
  TEveElementList * elementList = new TEveElementList("EMCAL Cells");
  elementList->SetTitle("Tooltip");
  gEve->AddElement(elementList);
  
  gStyle->SetPalette(1, 0);
  TEveRGBAPalette* pal = new TEveRGBAPalette(0, 512);
  pal->SetLimits(0, 1024);
  
  
  for (Int_t sm=0; sm<12; ++sm) {
    
    TEveBoxSet* q = new TEveBoxSet(Form("SM %d", sm+1));
    q->SetOwnIds(kTRUE);
    
    q->Reset(TEveBoxSet::kBT_AABox, kFALSE, 64);
    q->RefMainTrans().SetFrom(*gEMCALNode->GetDaughter(sm)->GetMatrix());
    q->SetPalette(pal);
    
    gEve->AddElement(q, elementList);
    gEMCALBoxSet[sm] = q;
  }

  return elementList;
}

void processHLTBlock(AliHLTHOMERBlockDesc * block) {

  if ( ! block->GetDataType().CompareTo("ALIESDV0") ) {
    if(!gTPCTrack){
      gTPCTrack = new TEveTrackList("ESD Tracks");
      gTPCTrack->SetMainColor(6);
      gEve->AddElement(gTPCTrack);
    }
    
    processEsdTracks(block, gTPCTrack);
    gTPCTrack->ElementChanged();
  } 
  
  else if ( ! block->GetDataType().CompareTo("ROOTTOBJ") ) {
    processROOTTOBJ( block, gHLTText );
  } 

  else if ( ! block->GetDataType().CompareTo("HLTRDLST") ) {
    processHLTRDLST( block );
  } 

  else if ( !block->GetDataType().CompareTo("ROOTHIST") ) {      
    
    processPrimVertexHistograms( block , gPrimVertexCanvas);
    gPrimVertexCanvas->Update();    
  }

} 


void processTPCBlock(AliHLTHOMERBlockDesc * block) {

  cout <<"Processing tpc clusters"<<endl;

  if ( ! block->GetDataType().CompareTo("CLUSTERS") ) {
    if(!gTPCClusters){	  
      gTPCClusters = new TEvePointSet("TPC Clusters");
      //gTPCClusters->ApplyVizTag("TPC Clusters");
      gTPCClusters->SetMainColor(kRed);
      gTPCClusters->SetMarkerStyle((Style_t)kFullDotSmall);
      //gEve->AddElement(gTPCClusters);
    } 
    
    if(!gTPCColClusters){
      gTPCColClusters = new TEvePointSetArray("TPC Clusters Colorized");
      gTPCColClusters->SetMainColor(kRed);
      gTPCColClusters->SetMarkerStyle(4); // antialiased circle
      gTPCColClusters->SetMarkerSize(0.4);
      gTPCColClusters->InitBins("Cluster Charge", gTPCBins, 0., gTPCBins*20.);
      
      const Int_t nCol = TColor::GetNumberOfColors();
      
      for (Int_t ii = 0; ii < gTPCBins+1; ++ii)
	gTPCColClusters->GetBin(ii)->SetMainColor(TColor::GetColorPalette(ii * nCol / (gTPCBins+2)));
      
      gEve->AddElement(gTPCColClusters);
    } 
    
    processTPCClusters(block, gTPCClusters, gTPCColClusters);
    //gTPCClusters->ElementChanged();
    gTPCColClusters->ElementChanged();
  }
      
  else if ( ! block->GetDataType().CompareTo("HWCL_ALT") ) {
    if(!gTPCTestClusters){	  
      
      gTPCTestClusters = new TEvePointSet("TPC Clusters Test");
      //ggTPCTestClusters->ApplyVizTag("TPC Clusters");
      gTPCTestClusters->SetMainColor(kBlue);
      gTPCTestClusters->SetMarkerStyle((Style_t)kFullDotSmall);
      gEve->AddElement(gTPCTestClusters);
    }
    
    processTPCClusters(block, gTPCTestClusters);
    gTPCTestClusters->ElementChanged();
  }
} 
