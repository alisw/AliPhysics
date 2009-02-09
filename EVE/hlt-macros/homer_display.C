// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

// Functions to read data from HOMER.
//
// Setup: edit location of HLT configuration in first line of
// homer_display(). This is a temporary solution.
//
// Run as: alieve command_queue.C+ hlt_structs.C+ homer_display.C
//
// nextEvent() will get next event from HOMER.


class AliRawReaderMemory;
class AliEveHOMERManager;
class AliHLTHOMERBlockDesc;

class TEvePointSet;
class TEveTrackList;
class TEveTrack;

class AliEveTPCLoader;
class AliEveTPCData;
class AliEveTPCSector2D;
class AliEveTPCSector3D;
class AliEveITSDigitsInfo;

//***********************************************************
#include "TTimer.h"
#include "TRandom.h"
#include "TVirtualPad.h"
#include "TGLViewer.h"
#include "TThread.h"
#include "TGFileBrowser.h"
#include "TStyle.h"
#include "TList.h"
#include "TDirectory.h"
//***********************************************************
#include "TEveManager.h"
#include "TEvePointSet.h"
#include "TEveTrack.h"
#include "TEveVSDStructs.h"
#include "TEveTrackPropagator.h"
#include "TEvePointSet.h"
#include "TEveScene.h"
#include "TEveElement.h"
// #include "TEveElementList.h"
#include "TEveEventManager.h"
//***********************************************************
#include "AliESDEvent.h"
#include "AliCDBManager.h"
#include "AliRawReaderMemory.h"
#include "AliTPCRawStream.h"
#include "AliGeomManager.h"
//***********************************************************
#include "AliEveHOMERManager.h"
#include "AliEveTPCLoader.h" 
#include "AliEveTPCData.h"
#include "AliEveITSDigitsInfo.h"
#include "AliEveITSModule.h"
//***********************************************************
#include "AliHLTHOMERBlockDesc.h"
#include "AliHLTHOMERReader.h"
//***********************************************************
#include "hlt_structs.C"
#include "TFile.h"
//***********************************************************
#include <AliHLTMUONUtils.h>
#include "AliHLTMUONDataBlockReader.h"
#include "tracking-ca/AliHLTTPCCATrackParam.h"

// -- globals --
 
AliEveTPCLoader*                          gTPCLoader    = 0;
AliEveTPCData*                            gTPCData      = 0;
TEvePointSet*                             gTPCClusters  = 0;
TEveTrackList*                            gTPCTrack     = 0;

AliEveITSDigitsInfo*                      gITSDigits    = 0; 

AliRawReaderMemory*                       gMemReader    = 0;
AliEveHOMERManager*                       gHomerManager = 0;

TEvePointSet*                             gMUONClusters  = 0;
Double_t                                  gSolenoidField = 5;
//***********************************************************

Int_t globMaxPoint = 0 ;

//***********************************************************

// -- needed below ??

Int_t    event  = -1;

TTimer   timer;
TTimer   event_timer;

TThread* ldthread = 0;

TRandom  rnd(0);

Bool_t vC = kFALSE;
TGFileBrowser *g_hlt_browser = 0;
TCanvas       *g_hlt_canvas  = 0;

TGLViewer::ECameraType camera = TGLViewer::kCameraPerspXOZ;

//****************************************************************************
Int_t nextEvent();

//****************************************************************************
Int_t processSPDRawData( AliHLTHOMERBlockDesc* block );
Int_t drawSPDRawData();

//****************************************************************************
Int_t processTPCRawData( AliHLTHOMERBlockDesc* block );
Int_t processTPCClusters( AliHLTHOMERBlockDesc* block );
Int_t processTPCTracks( AliHLTHOMERBlockDesc* block );

Int_t processMUONClusters( AliHLTHOMERBlockDesc* block );

//****************************************************************************
TEveTrack* makeESDTrack( TEveTrackPropagator*   rnrStyle,
			 Int_t                  idx,
			 AliESDtrack*           at,
			 AliExternalTrackParam* tp = 0 );

//****************************************************************************
void homer_display( Int_t run = 0) {

  AliCDBManager::Instance()->SetRun(run);
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  
  gMemReader = new AliRawReaderMemory(0, 0);

  gStyle->SetPalette(1, 0);
  gEve->DisableRedraw();

  // -- Create HOMER Manager
  gHomerManager = new AliEveHOMERManager("/local/home/hlt/AliEVE-Config.xml");

  // -- Set Realm  ( can be "GPN","ACR","HLT","KIP" )
  gHomerManager->SetRealm("ACR"); 
  
  gEve->AddToListTree(gHomerManager, kTRUE);

  // -- Create list of HOMER sources
  gHomerManager->CreateHOMERSourcesList();

  // -- TPC Loader
  gTPCLoader = new AliEveTPCLoader;
  gTPCLoader->SetDoubleSR(kTRUE);
  gTPCLoader->SetInitParams(40, 900, 2, 100);   // Sector params (mint, maxt, thr, maxval)
  
  // -- TPC Data
  gTPCData = gTPCLoader->GetData();
  gTPCData->SetLoadPedestal(0);
  gTPCData->SetLoadThreshold(0);
  gTPCData->SetAutoPedestal(kFALSE);             // For zero suppressed data.

  gEve->AddElement(gTPCLoader);

  gEve->Redraw3D(0,1); // (0, 1)
  gEve->EnableRedraw();  
}

//****************************************************************************
Int_t nextEvent() {

  Int_t iResult = 0;

  gStyle->SetPalette(1, 0);
  gEve->DisableRedraw();

  // ** Get Next Event from HOMER
  if ( gHomerManager->NextEvent() ) 
    return ++iResult;

  // ** Reset
  if ( gTPCClusters ) gTPCClusters->Reset();
  if ( gTPCTrack )    gTPCTrack->DestroyElements();
  if ( gTPCData )     gTPCData->DropAllSectors();

  if ( gITSDigits ) {
    delete gITSDigits;
    gITSDigits = NULL;
  }
  
  if ( gMUONClusters ) gMUONClusters->Reset();


  // ----------------------------------- foo A
  vC = kFALSE;

  // TList* hListKR = new TList;
  // TList* hListCF = new TList;

  // ----------------------------------- foo A

  TIter next(gHomerManager->GetBlockList());
  AliHLTHOMERBlockDesc* block = 0;

  while ((block = (AliHLTHOMERBlockDesc*)next())) {
   
//     printf ( "Det : %s\n" ,block->GetDetector().Data() );
//     printf ( "Datatype : %s\n" ,block->GetDataType().Data() );

    
    
    // -- TPC
    // -----------------------------------------------------
    
    if ( ! block->GetDetector().CompareTo("TPC") ){
      

    // +++ CLUSTERS BLOCK
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      if ( block->GetDataType().CompareTo("CLUSTERS") == 0 ) {
	if ( !gTPCClusters ) {
	  gTPCClusters = new TEvePointSet("TPC Clusters");
	  gTPCClusters->SetMainColor(kRed);
	  gTPCClusters->SetMarkerStyle((Style_t)kFullDotSmall);
	  gEve->AddElement(gTPCClusters);
	}
	
	// ** Process Clusters
	processTPCClusters( block );
	
	gTPCClusters->ElementChanged();
	
      }else if ( block->GetDataType().CompareTo("ESD_TREE") == 0 ) {
      // +++ ESD BLOCK
      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++

	// ** Initialize TPC Tracks
	if ( !gTPCTrack ) {
	  gTPCTrack = new TEveTrackList("TPC Tracks");
	  gTPCTrack->SetMainColor(kBlue);
	  
	  gEve->AddElement(gTPCTrack);
	
	  TEveTrackPropagator* rnrStyle = gTPCTrack->GetPropagator();
	  rnrStyle->SetMagField( 0 );
	  rnrStyle->SetFitDecay( 1 );
	}
	cout<<"SIZE : "<<block->GetSize()<<endl;
	// ** Proces Tracks
	processTPCTracks( block );
      }else if ( block->GetDataType().CompareTo("DDL_RAW") == 0 ) {
	processTPCRawData( block );
      }
      
    } else if ( ! block->GetDetector().CompareTo("MUON") ) {

      // -- MUON
      //-----------------------------------------------------
      if ( (block->GetDataType().CompareTo("RECHITS") == 0) || (block->GetDataType().CompareTo("TRIGRECS") == 0) ) {

//  	printf ( "Inside : Datatype : %s\n" ,block->GetDataType().Data() );
// 	printf ( "Inside : DataSize : %d\n" ,block->GetSize() );
	if ( !gMUONClusters ) {
	  gMUONClusters = new TEvePointSet("MUON RecHits");
	  gMUONClusters->SetMainColor(kBlue);
	  gMUONClusters->SetMarkerStyle(20);
	  gEve->AddElement(gMUONClusters);
	}
	
	// ** Process Clusters
 	processMUONClusters( block );
	
	gMUONClusters->ElementChanged();
	
      }//MUON Clusters


    }  else if ( ! block->GetDetector().CompareTo("SPD") ) {

      // -- SPD 
      // -----------------------------------------------------

      if ( block->GetDataType().CompareTo("DDL_RAW") == 0 ) {
	// ** Initialize SPD Digits
	if ( !gITSDigits ) {
	  gITSDigits = new AliEveITSDigitsInfo();
	}
	
	processSPDRawData( block );		
      }

    }else{

      printf ("Detector \"%s\" has not been recognized",block->GetDetector().Data());
    }
    
  }
  
  if ( gTPCLoader )   gTPCLoader->UpdateSectors( kTRUE );
  if ( gTPCClusters ) gTPCClusters->ResetBBox();
  if ( gTPCTrack )    gTPCTrack->MakeTracks();
  if ( gITSDigits )   drawSPDRawData();
  
  if ( gMUONClusters ) gMUONClusters->ResetBBox();
  
  gEve->Redraw3D(0,1); // (0, 1)
  gEve->EnableRedraw(); 
    
  return iResult;
}

//****************************************************************************
//****************************************************************************
//****************************************************************************
void loopEvent() {
  event_timer.SetCommand("nextEvent()");
  event_timer.Start(6000);
}

//****************************************************************************
void stopLoopEvent() {
  event_timer.Stop();
}

//****************************************************************************
Int_t processTPCRawData(AliHLTHOMERBlockDesc* block) {

  Int_t iResult = 0;

  Int_t sector = block->GetSubDetector().Atoi();
  Int_t patch  = block->GetSubSubDetector().Atoi();
  Int_t eqId   = 768 + patch;
  
  if ( patch >= 2) eqId += 4 * sector + 70;
  else   	   eqId += 2 * sector;

  printf("sector : %d %d %d -- %p %lu\n", sector, patch, eqId, block->GetData(), block->GetSize());

  gMemReader->SetMemory( reinterpret_cast<UChar_t*> ( block->GetData() ), block->GetSize() );
  gMemReader->SetEquipmentID( eqId );
  gMemReader->Reset();

  AliTPCRawStream tpcStream( gMemReader );
  gMemReader->Select("TPC"); 
  
  gTPCData->LoadRaw( tpcStream, kTRUE, kTRUE );

  return iResult;
}

//****************************************************************************
Int_t processTPCClusters(AliHLTHOMERBlockDesc* block) {
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

      gTPCClusters->SetNextPoint(cos*sp->fX - sin*sp->fY, sin*sp->fX + cos*sp->fY, sp->fZ);
    }
  }

  return iResult;
}

//****************************************************************************
TEveTrack* makeESDTrack( TEveTrackPropagator*   rnrStyle,
			 Int_t                  idx,
			 AliESDtrack*           esdTrack,
			 AliExternalTrackParam* trackParam  ) {
  // Helper function
  
  Double_t     pbuf[3], vbuf[3];
  TEveRecTrack rt;
  TEvePathMark startPoint(TEvePathMark::kReference);
  TEvePathMark midPoint(TEvePathMark::kReference);
  TEvePathMark mid1Point(TEvePathMark::kReference);
  TEvePathMark endPoint(TEvePathMark::kReference);
  TEvePathMark decPoint(TEvePathMark::kDecay);

  /*printf("ESD track: %f, %f, %f, %f, %f, %f, %f",
	 esdTrack->GetAlpha(),
	 esdTrack->GetX(), 
	 esdTrack->GetY(),
	 esdTrack->GetZ(),    
	 esdTrack->GetSnp(),  
	 esdTrack->GetTgl(),  
	 esdTrack->GetSigned1Pt()
	 );
  */
  cout<<"TPCPoints::"<<esdTrack->GetTPCPoints(0)<<" "<<esdTrack->GetTPCPoints(1)<<" "<<esdTrack->GetTPCPoints(2)<<" "<<esdTrack->GetTPCPoints(3)<<endl;

  if ( trackParam == 0 ) 
    trackParam = esdTrack;

  rt.fLabel  = esdTrack->GetLabel();
  rt.fIndex  = idx;
  rt.fStatus = (Int_t) esdTrack->GetStatus();
  rt.fSign   = (Int_t) trackParam->GetSign();

  Double_t x0 = trackParam->GetX();
  Double_t dx = esdTrack->GetTPCPoints(2) - x0;

  for( Double_t x1=x0; x1<x0+dx; x1+=(dx)*.1 ){//SG
    AliExternalTrackParam startParam = *trackParam;
    AliHLTTPCCATrackParam t;
    t.SetExtParam(startParam, gSolenoidField );
    if( !t.TransportToX(x1) ) continue; 
    t.GetExtParam( startParam, startParam.GetAlpha(), gSolenoidField );
    if( TMath::Abs(startParam.GetSnp())>.99 ) continue;

    startParam.GetXYZ(vbuf);
    
    { // get momentum manually because trackParam->GetPxPyPz doesn't works for straight lines

      Double_t pt= TMath::Abs(startParam.GetSigned1Pt());
      pt = (pt>kAlmost0) ?1./pt :100.;
      
      Double_t cA=TMath::Cos(startParam.GetAlpha()), sA=TMath::Sin(startParam.GetAlpha());
    
      Double_t sT=startParam.GetSnp();
      if( sT>kAlmost1 ){ sT = kAlmost1; }
      else if( sT<-kAlmost1 ){ sT = -kAlmost1; }
      Double_t cT = TMath::Sqrt(TMath::Abs(1 - sT*sT));
      
      pbuf[0] = pt*(cT*cA - sT*sA); 
      pbuf[1] = pt*(sT*cA + cT*sA); 
      pbuf[2] = pt*startParam.GetTgl();
    }
  
    break;
  }

  rt.fV.Set(vbuf);
  rt.fP.Set(pbuf);
  startPoint.fV.Set(vbuf);
  startPoint.fP.Set(pbuf);
  
  Double_t ep = esdTrack->GetP(), mc = esdTrack->GetMass();
  rt.fBeta = ep/TMath::Sqrt(ep*ep + mc*mc);

  TEveTrack* track = new TEveTrack(&rt, rnrStyle);

  cout<<"startPoint = "<<vbuf[0]<<" "<<vbuf[1]<<" "<<vbuf[2]<<" "<<pbuf[0]<<" "<<pbuf[1]<<" "<<pbuf[2]<<endl;

 
  for( ; TMath::Abs(dx)>=1.; dx*=.9 ){
    AliExternalTrackParam endParam = *trackParam;
    //if( !endParam.PropagateTo(x0+dx, gSolenoidField) ) continue; 
    AliHLTTPCCATrackParam t;
    t.SetExtParam(endParam, gSolenoidField );
    if( !t.TransportToX(x0+dx) ) continue; 
    t.GetExtParam( endParam, endParam.GetAlpha(), gSolenoidField );

    if( TMath::Abs(endParam.GetSnp())>.99 ) continue;
    
    { // get momentum manually because trackParam->GetPxPyPz doesn't works for straight lines

      Double_t pt= TMath::Abs(endParam.GetSigned1Pt());
      pt = (pt>kAlmost0) ?1./pt :100.;
      
      Double_t cA=TMath::Cos(endParam.GetAlpha()), sA=TMath::Sin(endParam.GetAlpha());
      
      Double_t sT=endParam.GetSnp();
      if( sT>=kAlmost1 ){ sT = kAlmost1; }
      else if( sT<-kAlmost1 ){ sT = -kAlmost1; }
      Double_t cT = TMath::Sqrt(TMath::Abs(1 - sT*sT));
      
      endParam.GetXYZ(vbuf);
      pbuf[0] = pt*(cT*cA - sT*sA);
      pbuf[1] = pt*(sT*cA + cT*sA);
      pbuf[2] = pt*endParam.GetTgl();
    }
    break;
  }
  endPoint.fV.Set(vbuf);
  endPoint.fP.Set(pbuf);
  decPoint.fV.Set(vbuf);
  decPoint.fP.Set(pbuf);

  cout<<"endPoint = "<<vbuf[0]<<" "<<vbuf[1]<<" "<<vbuf[2]<<" "<<pbuf[0]<<" "<<pbuf[1]<<" "<<pbuf[2]<<endl;

  dx*=.6;
  for( ; TMath::Abs(dx)>=.5; dx*=.8 ){
    AliExternalTrackParam endParam = *trackParam;
    //if( !endParam.PropagateTo(x0+dx, gSolenoidField) ) continue; 
    AliHLTTPCCATrackParam t;
    t.SetExtParam(endParam, gSolenoidField );
    if( !t.TransportToX(x0+dx) ) continue; 
    t.GetExtParam( endParam, endParam.GetAlpha(), gSolenoidField );
    if( TMath::Abs(endParam.GetSnp())>.99 ) continue;
    
    { // get momentum manually because trackParam->GetPxPyPz doesn't works for straight lines

      Double_t pt= TMath::Abs(endParam.GetSigned1Pt());
      pt = (pt>kAlmost0) ?1./pt :100.;
      
      Double_t cA=TMath::Cos(endParam.GetAlpha()), sA=TMath::Sin(endParam.GetAlpha());
      
      Double_t sT=endParam.GetSnp();
      if( sT>=kAlmost1 ){ sT = kAlmost1; }
      else if( sT<-kAlmost1 ){ sT = -kAlmost1; }
      Double_t cT = TMath::Sqrt(TMath::Abs(1 - sT*sT));
      
      endParam.GetXYZ(vbuf);
      pbuf[0] = pt*(cT*cA - sT*sA);
      pbuf[1] = pt*(sT*cA + cT*sA);
      pbuf[2] = pt*endParam.GetTgl();
    }
    break;
  }

  mid1Point.fV.Set(vbuf);
  mid1Point.fP.Set(pbuf);

  //cout<<"midPoint = "<<vbuf[0]<<" "<<vbuf[1]<<" "<<vbuf[2]<<" "<<pbuf[0]<<" "<<pbuf[1]<<" "<<pbuf[2]<<endl;
 
  dx*=.5;
  for( ; TMath::Abs(dx)>=.5; dx*=.8 ){
    AliExternalTrackParam endParam = *trackParam;
    //if( !endParam.PropagateTo(x0+dx, gSolenoidField) ) continue; 
    AliHLTTPCCATrackParam t;
    t.SetExtParam(endParam, gSolenoidField );
    if( !t.TransportToX(x0+dx) ) continue; 
    t.GetExtParam( endParam, endParam.GetAlpha(), gSolenoidField );
    if( TMath::Abs(endParam.GetSnp())>.99 ) continue;
    
    { // get momentum manually because trackParam->GetPxPyPz doesn't works for straight lines

      Double_t pt= TMath::Abs(endParam.GetSigned1Pt());
      pt = (pt>kAlmost0) ?1./pt :100.;
      
      Double_t cA=TMath::Cos(endParam.GetAlpha()), sA=TMath::Sin(endParam.GetAlpha());
      
      Double_t sT=endParam.GetSnp();
      if( sT>=kAlmost1 ){ sT = kAlmost1; }
      else if( sT<-kAlmost1 ){ sT = -kAlmost1; }
      Double_t cT = TMath::Sqrt(TMath::Abs(1 - sT*sT));
      
      endParam.GetXYZ(vbuf);
      pbuf[0] = pt*(cT*cA - sT*sA);
      pbuf[1] = pt*(sT*cA + cT*sA);
      pbuf[2] = pt*endParam.GetTgl();
    }
    break;
  }

  midPoint.fV.Set(vbuf);
  midPoint.fP.Set(pbuf);

  track->AddPathMark( startPoint );
  track->AddPathMark( midPoint );
  track->AddPathMark( mid1Point );
  track->AddPathMark( endPoint );
  track->AddPathMark( decPoint );
  
  //PH The line below is replaced waiting for a fix in Root
  //PH which permits to use variable siza arguments in CINT
  //PH on some platforms (alphalinuxgcc, solariscc5, etc.)
  //PH  track->SetName(Form("ESDTrack %d", rt.label));
  //PH  track->SetTitle(Form("pT=%.3f, pZ=%.3f; V=(%.3f, %.3f, %.3f)",
  //PH		       rt.sign*TMath::Hypot(rt.P.x, rt.P.y), rt.P.z,
  //PH		       rt.V.x, rt.V.y, rt.V.z));
  char form[1000];
  sprintf(form,"TEveTrack %d", rt.fIndex);
  track->SetName(form);
  track->SetStdTitle();
  return track;
}

//****************************************************************************
Int_t processTPCTracks(AliHLTHOMERBlockDesc* block) {

  TTree *tr = (TTree*) block->GetTObject();

//   ofstream fout("ESD_TPC.dat",ios::binary);
//   fout.write((char*)block->GetData(),block->GetSize());
//   fout.close();
  TFile f("ESD_TPC.root","recreate");
  tr->Write();
  f.Close();

  TFile* esdFile = TFile::Open("ESD_TPC.root");
  
  AliESDEvent* esd = new AliESDEvent();
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  esd->ReadFromTree(tree);
  //tr->SetBranchAddress("ESD", &esd);
  //   if(tr->GetBranch("ESD"))
  //     return 0;

  tree->GetEntry(0);

  TEveTrackPropagator* rnrStyle = gTPCTrack->GetPropagator();
  rnrStyle->SetMagField( 0.1*esd->GetMagneticField() );
  gSolenoidField = esd->GetMagneticField();

  cout << "Number of tracks found :" << esd->GetNumberOfTracks() << endl;

  Double_t pin[3];
  
  for (Int_t ii=0; ii< esd->GetNumberOfTracks(); ii++) {

    AliESDtrack           *esdTrack = esd->GetTrack(ii);
    AliExternalTrackParam *trackParam = esdTrack;
    cout<<"\nESD track N"<<ii<<":"<<endl;
    trackParam->Print();

    TEveTrack* track = makeESDTrack( rnrStyle, ii, esdTrack, trackParam );
    esdTrack->GetPxPyPz(pin);

    //cout<<"pt : "<<sqrt(pin[0]*pin[0] + pin[1]*pin[1])<<endl;
    track->SetAttLineAttMarker(gTPCTrack);    
   gEve->AddElement(track, gTPCTrack);
 }
  
  delete esd;
  
  return 0;
}

//****************************************************************************
Int_t processSPDRawData(AliHLTHOMERBlockDesc* block) {
  Int_t iResult = 0;

  Int_t partition = block->GetSubDetector().Atoi();
  Int_t eqId      = partition;

  gMemReader->SetMemory( reinterpret_cast<UChar_t*> ( block->GetData() ), block->GetSize() );
  gMemReader->SetEquipmentID( eqId );
  gMemReader->Reset();

  gITSDigits->ReadRaw( gMemReader, 3);

  return iResult;
}

//****************************************************************************
Int_t drawSPDRawData() {

  Int_t iResult = 0;

  TString sSector;
  TString bsSector="Sector";
  TString sStave;
  TString bsStave="Stave";

  Int_t ndx=0;
  Int_t sector, stave, module;

  //gEve->DisableRedraw();

  // ** first layer **

  TEveElementList* layer1 = new TEveElementList( "SPD0" );
  layer1->SetTitle( "SPDs' first layer" );
  layer1->SetMainColor(2);
  gEve->AddElement( layer1 );
  
  for ( sector=0; sector<10; sector++ ) {
    sSector  = bsSector;
    sSector += sector;

    TEveElementList* relSector = new TEveElementList( sSector.Data() );
    relSector->SetMainColor(2);
    gEve->AddElement( relSector, layer1 );

    for ( stave=0; stave<2; stave++ ) {
      sStave  = bsStave;
      sStave += stave;
      
      TEveElementList* relStave = new TEveElementList( sStave.Data() );
      relStave->SetMainColor(2);
      gEve->AddElement( relStave, relSector );

      for ( module=0; module<4; module++ ) {
	
	if ( gITSDigits->GetDigits( ndx, 0 ) && 
	     gITSDigits->GetDigits( ndx, 0 )->GetEntriesFast() > 0) {
	  
	  AliEveITSModule* moduleITS = new AliEveITSModule( ndx, gITSDigits );
	  gEve->AddElement( moduleITS, relStave );
	}

	++ndx; 

      } // for ( module=0; module<4; module++ ) {
    } // for ( stave=0; stave<2; stave++ ) {
  } // for ( sector=0; sector<10; sector++ ) {

  // ** second layer **

  TEveElementList* layer2 = new TEveElementList( "SPD1" );
  layer2->SetTitle( "SPDs' second layer" );
  layer2->SetMainColor(2);
  gEve->AddElement(layer2);
  
  for ( sector=0; sector<10; sector++ ) {
    sSector  = bsSector;
    sSector += sector;
    
    TEveElementList* relSector = new TEveElementList( sSector.Data() );
    relSector->SetMainColor(2);
    gEve->AddElement(relSector, layer2 );
    
    for ( stave=0; stave<4; stave++ ) {
      sStave  = bsStave;
      sStave += stave;

      TEveElementList* relStave = new TEveElementList( sStave.Data() );
      relStave->SetMainColor(2);
      gEve->AddElement( relStave, relSector );

      for ( module=0; module<4; module++) {

	if ( gITSDigits->GetDigits( ndx, 0 ) && 
	     gITSDigits->GetDigits( ndx, 0 )->GetEntriesFast() > 0) {

	  AliEveITSModule* moduleITS = new AliEveITSModule( ndx, gITSDigits );
	  gEve->AddElement( moduleITS, relStave );
	}
	 
	++ndx;
      } // for ( module=0; module<4; module++) {
    } // for ( stave=0; stave<2; stave++ ) {
  } //for ( sector=0; sector<10; sector++ ) {

  gEve->EnableRedraw();
    
  return iResult;
}

//****************************************************************************
Int_t processMUONClusters(AliHLTHOMERBlockDesc* block) {
  Int_t iResult = 0;

      unsigned long size = block->GetSize();
      int * buffer ;

      //       for(int idata=0;idata<int(size);idata++)
      //  	printf("\tbuffer[%d] : %d\n",idata,buffer[idata]);
      
      buffer = (int *)block->GetData();
      
      if(block->GetDataType().CompareTo("RECHITS") == 0){
	
	AliHLTMUONRecHitsBlockReader trackblock((char*)buffer, size);
	const AliHLTMUONRecHitStruct* hit = trackblock.GetArray();
      
	for(AliHLTUInt32_t ientry = 0; ientry < trackblock.Nentries(); ientry++){
// 	  cout << setw(13) << left << hit->fX << setw(0);
// 	  cout << setw(13) << left << hit->fY << setw(0);
// 	  cout << hit->fZ << setw(0) << endl;
	  gMUONClusters->SetNextPoint(hit->fX,hit->fY,hit->fZ);
	  hit++;
	  
	}// track hit loop
      }else{// if rechits
      
      
//       if(!strcmp((BlockType(ULong64_t(reader->GetBlockDataType(i)))).Data(),"TRIGRECS")){

	AliHLTMUONTriggerRecordsBlockReader trigblock(buffer, size);
	const AliHLTMUONTriggerRecordStruct* trigrec = trigblock.GetArray();
	for(AliHLTUInt32_t ientry = 0; ientry < trigblock.Nentries(); ientry++){

	  const AliHLTMUONRecHitStruct* hit = &trigrec->fHit[0];
	  for(AliHLTUInt32_t ch = 0; ch < 4; ch++)
	    {
// 	      cout << setw(10) << left << ch + 11 << setw(0);
// 	      cout << setw(13) << left << hit->fX << setw(0);
// 	      cout << setw(13) << left << hit->fY << setw(0);
// 	      cout << hit->fZ << setw(0) << endl;
	      gMUONClusters->SetNextPoint(hit->fX,hit->fY,hit->fZ);
	      hit++;
	    }// trig chamber loop
	  
	}//trig hit loop
      }// if trigrecs
      
      
//       //delete[] buffer;
//     }//nof Block received
//   }// if any event found
  
  
//   delete reader;


  return iResult;
}


//****************************************************************************
