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

#include "geom_gentle_hlt.C"

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

Int_t initializeEveViewer( Bool_t showBarrel, Bool_t showMuon );

void writeToFile();


// #################################################################
// #################################################################
// #################################################################

// -----------------------------------------------------------------
void od ( Bool_t showBarrel = kTRUE, Bool_t showMuon = kFALSE ) {

  // -- Loading Geometry
  // ---------------------
  Int_t run = 67179;
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetRun(run);
  AliGeomManager::LoadGeometry();


  // Get the pointer to gGeoManager before it's broken (bug in alieve)
  fGeoManager = gGeoManager;

    // -- Initialize pointsets and add macros
  // ----------------------------------------
  //TEveUtil::LoadMacro("hlt_alieve_init.C");
  //hlt_alieve_init(".", -1);

  // -- Initialize Eve
  // -------------------
  initializeEveViewer( showBarrel, showMuon );

  // -- Reset gGeoManager to the original pointer
  // ----------------------------------------------

  // -- Finalize Eve
  // -----------------
  gSystem->ProcessEvents();
  gEve->Redraw3D(kTRUE);

  // -- Create new hM object
  // -------------------------
  gHomerManager = new AliEveHOMERManager();
  gHomerManager->SetRetryCount(1000,15);
  gHomerManager->SetEveManager(gEve);
  gHomerManager->SetGeoManager(gGeoManager);
  gHomerManager->SetRPhiManager(gRPhiMgr);
  gHomerManager->SetRPhiEventScene(gRPhiEventScene);
  gHomerManager->SetRPhiViewer(gRPhiView);
  gHomerManager->SetRhoZManager(gRhoZMgr);
  gHomerManager->SetRhoZEventScene(gRhoZEventScene);
  gHomerManager->SetRhoZViewer(gRhoZView);
  gHomerManager->SetBarrelFlag(showBarrel);
  gHomerManager->SetMuonFlag(showMuon);

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
    printf ("Couldn't find active services. Giving up. \n");
    return;
  } 


  if ( showBarrel ) {
    gHomerManager->ConnectEVEtoHOMER("TPC" );
  } else if ( MUONMode ) {
    gHomerManager->ConnectEVEtoHOMER("MUON");
  } else if( TRDMode ) {
    gHomerManager->ConnectEVEtoHOMER("TRD");  
  } else {
    cout<<" No detectors selected, nothing will be displayed"<<endl;
  }	

  gGeoManager = fGeoManager;
  

}

// -------------------------------------------------------------------------
Int_t initializeEveViewer( Bool_t showBarrel, Bool_t showMuon ) {
  
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

  //  processEvent();
}



// -----------------------------------------------------------------


Int_t updateDisplay() {

  Int_t iResult = 0;

  
    
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



//*****************************************************************************
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


// -----------------------------------------------------------------
void loopEvent() {
  eventTimer.SetCommand("nextEvent()");
  eventTimer.Start(3000);
}

// -----------------------------------------------------------------
void stopLoopEvent() {
  eventTimer.Stop();
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
      gTRDCanvas = createCanvas("TRD", "TRD");
      gTRDCanvas->Divide(3,2);
    }
    iResult = processTRDHistograms( block, gTRDCanvas );     
  }

  else if(block->GetDataType().CompareTo("CALIBRAH")==0){
     
    if(!gTRDCalibCanvas){
      gTRDCalibCanvas = createCanvas("TRD Calib", "TRD Calib");
      gTRDCalibCanvas->Divide(2,2);
    }
     
    iResult=processTRDCalibHistograms(block,gTRDCalibCanvas);
  }

  else if(block->GetDataType().CompareTo("CALIBEOR")==0){
     
    if(!gTRDEORCanvas){
      gTRDEORCanvas = createCanvas("TRD QA", "TRD QA");
      gTRDEORCanvas->Divide(3,2);       
    }
  
    iResult=processTRDCalibHistograms(block,gTRDEORCanvas);
  }
  return iResult;
}



        

