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

#include "AliEveHOMERManager.h"
#include "AliEveHLTEventManager.h"
#include "geom_gentle_hlt.C"

//***************************************************************
#include "HLT/rec/AliHLTReconstructor.h"



#endif



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


// -- HOMERManager
AliEveHOMERManager*                       gHomerManager      = 0;
AliEveHLTEventManager*                       geventManager      = 0;

// -- Geometry Manager 
TGeoManager*                              gGeoManager        = 0;


//Container for gGeoManager till it is broken
TGeoManager *fGeoManager = 0;
// -----------------------------------------------------------------
// --                          Methods                            --
// -----------------------------------------------------------------

Int_t initializeEveViewer( Bool_t showBarrel, Bool_t showMuon );

void writeToFile();


// #################################################################
// #################################################################
// #################################################################

// -----------------------------------------------------------------
void od ( Bool_t showBarrel = kTRUE, Bool_t showMuon = kFALSE ) {

  // -- Loading Geometry
  // ---------------------
  Int_t run = 0;
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetRun(run);
  AliGeomManager::LoadGeometry();
  // The default in the simulation is the following line
  // TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1, AliMagF::k5kG));
  // However for the current setting of +ve L3 and +ve Dipole magnetic field
  // the following setting creates the field close to real field with currect polarity
  if(showMuon)
    TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1, AliMagF::k5kG));

  // Get the pointer to gGeoManager before it's broken (bug in alieve)
  fGeoManager = gGeoManager;

    // -- Initialize pointsets and add macros
  // ----------------------------------------
  //TEveUtil::LoadMacro("hlt_alieve_init.C");
  //hlt_alieve_init(".", -1);

  // -- Initialize Eve
  // -------------------
  cout << "Initializing the EVE viewer"<<endl;
  initializeEveViewer( showBarrel, showMuon, fGeoManager );

  // -- Reset gGeoManager to the original pointer
  // ----------------------------------------------

  // -- Finalize Eve
  // -----------------
  gSystem->ProcessEvents();
  gEve->Redraw3D(kTRUE);

  // -- Create new hM object
  // -------------------------

  cout << "Creating the Event Manager"<<endl;
  gEventManager = new AliEveHLTEventManagerHomer();
  gEventManager->SetShowMuon(showMuon);
  gEventManager->SetEveManager(gEve);
  gEventManager->SetGeoManager(gGeoManager);
  gEventManager->SetRPhiManager(gRPhiMgr);
  gEventManager->SetRPhiEventScene(gRPhiEventScene);
  gEventManager->SetRPhiViewer(gRPhiView);
  gEventManager->SetRhoZManager(gRhoZMgr);
  gEventManager->SetRhoZEventScene(gRhoZEventScene);
  gEventManager->SetRhoZViewer(gRhoZView);
 
  //gEventManager->SetBarrelFlag(showBarrel);
  //gEventManager->SetMuonFlag(showMuon);

  // Int_t iResult = gHomerManager->Initialize();
  // if (iResult) { 
  //   printf("Error Initializing AliHLTHOMERManager, quitting");
  //   return; 
  // }
  
  // gEventManager->SetHomerManager(gHomerManager);

  // -- Add hM to EveTree
  // ----------------------
  //gEve->AddToListTree(gHomerManager, kTRUE);
  gEve->AddToListTree(gEventManager, kTRUE);

  // -- Create SourceList
  // ----------------------
  // iResult = gHomerManager->CreateEveSourcesListLoop();
  // if (iResult) {
  //   printf ("Couldn't find active services. Giving up. \n");
  //   return;
  // } 


  // if ( showBarrel ) {
  //   gHomerManager->ConnectEVEtoHOMER("TPC" );
  // } else if ( MUONMode ) {
  //   gHomerManager->ConnectEVEtoHOMER("MUON");
  // } else if( TRDMode ) {
  //   gHomerManager->ConnectEVEtoHOMER("TRD");  
  // } else {
  //   cout<<" No detectors selected, nothing will be displayed"<<endl;
  // }	
  // THIS LINE DOES NOT WORK. PLEASE USE ANOTHER LINE!


  gGeoManager = fGeoManager;
  

}

// -------------------------------------------------------------------------
Int_t initializeEveViewer( Bool_t showBarrel, Bool_t showMuon, TGeoManager * manager ) {
  
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
  if (!showMuon)
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
  emcal_re->SetVisLevel(1);
  //  emcal_re->FirstChild()->Dump();

  // for(Int_t i = 4; i < 11; i++) {
  //   emcal_re->FindChild(Form("SMOD_%d", i))->SetRnrState(kFALSE);
  // }
  // emcal_re->FindChild("SM10_1")->SetRnrState(kFALSE);
  // emcal_re->FindChild("SM10_2")->SetRnrState(kFALSE);



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
  gRhoZMgr->ImportElements(gGeomGentleTRD, gRhoZGeomScene);
  
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
  //ob1->SetBackColor(10492431);
  TGLOverlayButton *ob2 = new TGLOverlayButton(g3DView->GetGLViewer(),  "ALICE", 0, 0, 110, 20);
  ob2->SetAlphaValues(0.8, 0.8);
  //ob2->SetBackColor(0.2);
  TGLOverlayButton *ob3 = new TGLOverlayButton(gEve->GetDefaultGLViewer(),  "HLT", 0, 20, 110, 60);
  ob3->SetAlphaValues(0.8, 0.8);
  TGLOverlayButton *ob4 = new TGLOverlayButton(gEve->GetDefaultGLViewer(),  "ALICE", 0, 0, 110, 20);
  ob4->SetAlphaValues(0.8, 0.8);


  TGLOverlayButton *ne = new TGLOverlayButton(gEve->GetDefaultGLViewer(),  "Next Event", 110, 0, 210, 20);
  ne->SetAlphaValues(0.0, 0.8);

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



  gStyle->SetPalette(1, 0);
  gStyle->SetOptFit(1);


  
  return 0;
}

//****************************************************************************
void writeToFile(){

  TList * bList = gHomerManager->GetBlockList();
  if(bList){
    TFile * file = TFile::Open(Form("Event_0x%016X_ITS.root", gHomerManager->GetEventID()), "RECREATE"); 
    bList->Write("blockList", TObject::kSingleKey);
    file->Close();
  }
  
  bList = gHomerManager->GetAsyncBlockList();
  if(bList){
    TFile * afile = TFile::Open(Form("Event_0x%016X_Async.root", gHomerManager->GetEventID()), "RECREATE"); 
    bList->Write("blockList", TObject::kSingleKey);
    afile->Close();
  }
}



        
