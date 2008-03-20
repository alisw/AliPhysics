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

#include "TTimer.h"
#include "TRandom.h"
#include "TVirtualPad.h"

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

#include "TGLViewer.h"
#include "TThread.h"
#include "TGFileBrowser.h"
#include "TStyle.h"
#include "TList.h"
#include "TDirectory.h"

#include "TEveManager.h"
#include "TEvePointSet.h"
#include "TEveTrack.h"
#include "TEveVSDStructs.h"
#include "TEveTrackPropagator.h"
#include "TEveElement.h"

#include "AliESDEvent.h"
#include "AliCDBManager.h"
#include "AliRawReaderMemory.h"

#include "AliEveHOMERManager.h"
#include "AliEveTPCLoader.h"
#include "AliEveTPCData.h"
#include "AliEveITSDigitsInfo.h"
#include "AliEveITSModule.h"

#include "AliHLTHOMERBlockDesc.h"
#include "AliTPCRawStream.h"

//***********************************************************
#include "DIMUONRawReader.C"

#include "TEvePointSet.h"
#include "TEveScene.h"
#include "AliEveMUONData.h"
#include "AliEveMUONChamber.h"
#include "AliEveMUONChamberData.h"

#include "AliGeomManager.h"

#include "AliMpCDB.h"
#include "AliMpDDLStore.h"
#include "AliMpDetElement.h"
#include "AliMpDEIterator.h"
#include "AliMpVSegmentation.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMpSegmentation.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONVCalibParam.h"
#include "TEveEventManager.h"
#include "AliMpTriggerCrate.h"
#include "AliMpLocalBoard.h"
#include "AliMUONGeometryDetElement.h"

//***********************************************************

// -- globals --

AliEveTPCLoader*                          gTPCLoader    = 0;
AliEveTPCData*                            gTPCData      = 0;
TEvePointSet*                             gTPCClusters  = 0;
TEveTrackList*                            gTPCTrack     = 0;

AliEveITSDigitsInfo*                      gITSDigits    = 0; 

AliRawReaderMemory*                       gMemReader    = 0;
AliEveHOMERManager*                       gHomerManager = 0;

//***********************************************************

map <UInt_t,AliHLTMUONTrackerMappingData> gDimuTrackerMapping;
vector<UInt_t>                            gDimuTrackerDataList;
vector<AliHLTMUONTriggerPointData>        gDimuTriggerDataList;
AliHLTMUONTriggerMappingData              gDimuTriggerMapping[2];
Bool_t                                    gMUONRawData = false;


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
Int_t nextEvent(Int_t ADCCut = 6);

//****************************************************************************
Int_t processSPDRawData( AliHLTHOMERBlockDesc* block );
Int_t drawSPDRawData();

//****************************************************************************
Int_t processTPCRawData( AliHLTHOMERBlockDesc* block );
Int_t processTPCClusters( AliHLTHOMERBlockDesc* block );
Int_t processTPCTracks( AliHLTHOMERBlockDesc* block );

//****************************************************************************
Int_t initDiMuonMapping();
Int_t initTrackerMapping() ;
Int_t processDiMuonRawData( AliHLTHOMERBlockDesc* block );
Int_t drawDiMuonRawData(Int_t);

//****************************************************************************
TEveTrack* makeESDTrack( TEveTrackPropagator*   rnrStyle,
			 Int_t                  idx,
			 AliESDtrack*           at,
			 AliExternalTrackParam* tp = 0 );

//****************************************************************************
void homer_display(Int_t run = 0) {

  AliCDBManager::Instance()->SetRun(run);
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  
  gMemReader = new AliRawReaderMemory(0, 0);

  gStyle->SetPalette(1, 0);

  // -- Create HOMER Manager
  gHomerManager = new AliEveHOMERManager("/local/home/hlt/AliEVE-Config.xml");

  // -- Set Realm  ( can be "GPN","ACR","HLT" )
  gHomerManager->SetRealm("ACR"); 
  
  gEve->AddToListTree(gHomerManager, kTRUE);

  // -- Create list of HOMER sources
  gHomerManager->CreateHOMERSourcesList();

  // -- TPC Loader
  gTPCLoader = new AliEveTPCLoader;
  gTPCLoader->SetDoubleSR(kTRUE);
  gTPCLoader->SetInitParams(40, 900, 10, 100);   // Sector params (mint, maxt, thr, maxval)
  
  // -- TPC Data
  gTPCData = gTPCLoader->GetData();
  gTPCData->SetLoadPedestal(0);
  gTPCData->SetLoadThreshold(0);
  gTPCData->SetAutoPedestal(kFALSE);             // For zero suppressed data.

  // -- Load MUON mapping
  initDiMuonMapping();  

  gEve->AddElement(gTPCLoader);
}

//****************************************************************************
Int_t nextEvent(Int_t ADCCut) {

  Int_t iResult = 0;

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
  
  if(gDimuTrackerDataList.size()>0)
    gDimuTrackerDataList.clear();

  if(gDimuTriggerDataList.size()>0)
    gDimuTriggerDataList.clear();

  for(Int_t ipoint=0;ipoint<globMaxPoint;ipoint++)
    point3d[ipoint].SetPoint(0,0.0,0.0,0.0);


  // ----------------------------------- foo A
  vC = kFALSE;

  // TList* hListKR = new TList;
  // TList* hListCF = new TList;

  gDirectory = 0;
  // ----------------------------------- foo A

  TIter next(gHomerManager->GetBlockList());
  AliHLTHOMERBlockDesc* block = 0;

  while ((block = (AliHLTHOMERBlockDesc*)next())) {

    printf ( "Det : %s\n" ,block->GetDetector().Data() );
    printf ( "Datatype : %s\n" ,block->GetDataType().Data() );
    
    // +++ CLUSTERS BLOCK
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if ( block->GetDataType().CompareTo("CLUSTERS") == 0 ) {

      // ** Initialize TPC Clusters
      if ( !gTPCClusters ) {
	gTPCClusters = new TEvePointSet("TPC Clusters");
	gTPCClusters->SetMainColor((Color_t)kRed);
	gTPCClusters->SetMarkerStyle((Style_t)kFullDotSmall);
	gEve->AddElement(gTPCClusters);
      }

      // ** Process Clusters
      processTPCClusters( block );
    
      gTPCClusters->ElementChanged();

      // +++ ESD BLOCK
      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    } else if ( block->GetDataType().CompareTo("ESD_TREE") == 0 ) {

      // ** Initialize TPC Tracks
      if ( !gTPCTrack ) {
	gTPCTrack = new TEveTrackList("TPC Tracks");
	gTPCTrack->SetMainColor((Color_t)kBlue);

	gEve->AddElement(gTPCTrack);
	
	TEveTrackPropagator* rnrStyle = gTPCTrack->GetPropagator();
	rnrStyle->SetMagField( 0 );
	rnrStyle->SetFitDecay( 1 );
      }

      // ** Proces Tracks
      processTPCTracks( block );
      
      // +++ RAW DATA BLOCK
      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    } else if ( block->GetDataType().CompareTo("DDL_RAW") == 0 ) {

      // -- TPC
      // -----------------------------------------------------
      if ( ! block->GetDetector().CompareTo("TPC") )
	processTPCRawData( block );

      // -- SPD 
      // -----------------------------------------------------
      else if ( ! block->GetDetector().CompareTo("SPD") ) {
	// !!!!!! TODO add also other ITS detectors

	// ** Initialize SPD Digits
	if ( !gITSDigits ) {
	  gITSDigits = new AliEveITSDigitsInfo();
	}
	  
	processSPDRawData( block );		
      }

      // -- MUON
      // -----------------------------------------------------
      else if ( ! block->GetDetector().CompareTo("MUON") ) {

	if( processDiMuonRawData ( block ) )
	  gMUONRawData = true;
      }
      
      // +++ KRYPTON HISTOGRAM BLOCK
      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    } else if (block->GetDataType().CompareTo("KRPTHIST") == 0) {
      
      


      // ++ else
      // +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    } else {


      printf ("This is nothing");
    }
    
  }

  if ( gTPCLoader )   gTPCLoader->UpdateSectors( kTRUE );
  if ( gTPCClusters ) gTPCClusters->ResetBBox();
  if ( gTPCTrack )    gTPCTrack->MakeTracks();
  if ( gITSDigits )   drawSPDRawData();
  if ( gMUONRawData ) drawDiMuonRawData ( ADCCut ) ;

  gEve->Redraw3D(0,1); // (0, 1)
  gMUONRawData = false;
  
  return iResult;
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

  gEve->DisableRedraw();

  // ** first layer **

  TEveElementList* layer1 = new TEveElementList( "SPD0" );
  layer1->SetTitle( "SPDs' first layer" );
  layer1->SetMainColor( (Color_t) 2 );
  gEve->AddElement( layer1 );
  
  for ( sector=0; sector<10; sector++ ) {
    sSector  = bsSector;
    sSector += sector;

    TEveElementList* relSector = new TEveElementList( sSector.Data() );
    relSector->SetMainColor( (Color_t) 2 );
    gEve->AddElement( relSector, layer1 );

    for ( stave=0; stave<2; stave++ ) {
      sStave  = bsStave;
      sStave += stave;
      
      TEveElementList* relStave = new TEveElementList( sStave.Data() );
      relStave->SetMainColor( (Color_t) 2 );
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
  layer2->SetMainColor( (Color_t) 2 );
  gEve->AddElement(layer2);
  
  for ( sector=0; sector<10; sector++ ) {
    sSector  = bsSector;
    sSector += sector;
    
    TEveElementList* relSector = new TEveElementList( sSector.Data() );
    relSector->SetMainColor( (Color_t) 2 );
    gEve->AddElement(relSector, layer2 );
    
    for ( stave=0; stave<4; stave++ ) {
      sStave  = bsStave;
      sStave += stave;

      TEveElementList* relStave = new TEveElementList( sStave.Data() );
      relStave->SetMainColor( (Color_t) 2 );
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
Int_t processTPCRawData(AliHLTHOMERBlockDesc* block) {

  Int_t iResult = 0;

  Int_t sector = block->GetSubDetector().Atoi();
  Int_t patch  = block->GetSubSubDetector().Atoi();
  Int_t eqId   = 768 + patch;
  
  if ( patch >= 2) eqId += 4 * sector + 70;
  else   	   eqId += 2 * sector;

  printf("%d %d %d -- %p %lu\n", sector, patch, eqId, block->GetData(), block->GetSize());

  gMemReader->SetMemory( reinterpret_cast<UChar_t*> ( block->GetData() ), block->GetSize() );
  gMemReader->SetEquipmentID( eqId );
  gMemReader->Reset();

  AliTPCRawStream tpcStream( gMemReader );
  tpcStream.SetOldRCUFormat(kTRUE);      
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

  if ( trackParam == 0 ) 
    trackParam = esdTrack;

  rt.fLabel  = esdTrack->GetLabel();
  rt.fIndex  = idx;
  rt.fStatus = (Int_t) esdTrack->GetStatus();
  rt.fSign   = (Int_t) trackParam->GetSign();
  trackParam->GetXYZ(vbuf);
  rt.fV.Set(vbuf);
  
  { // get momentum manually because trackParam->GetPxPyPz doesn't works for straight lines

    Double_t pt= TMath::Abs(trackParam->GetSigned1Pt());
    pt = (pt>kAlmost0) ?1./pt :100.;

    Double_t cA=TMath::Cos(trackParam->GetAlpha()), sA=TMath::Sin(trackParam->GetAlpha());

    Double_t sT=trackParam->GetSnp();
    if( sT>=kAlmost1 ){ sT = kAlmost1; }
    else if( sT<-kAlmost1 ){ sT = -kAlmost1; }
    Double_t cT = TMath::Sqrt(TMath::Abs(1 - sT*sT));

    pbuf[0] = pt*(cT*cA - sT*sA); 
    pbuf[1] = pt*(sT*cA + cT*sA); 
    pbuf[2] = pt*trackParam->GetTgl();
  }
  
  rt.fP.Set(pbuf);
  
  Double_t ep = esdTrack->GetP(), mc = esdTrack->GetMass();
  rt.fBeta = ep/TMath::Sqrt(ep*ep + mc*mc);

  TEveTrack* track = new TEveTrack(&rt, rnrStyle);

  AliExternalTrackParam endParam = *trackParam;
  if( endParam.PropagateTo(esdTrack->GetTPCPoints(2), 0.) ){ // 5 kG

    TEvePathMark *startPoint = new TEvePathMark(TEvePathMark::kReference);
    TEvePathMark *endPoint = new TEvePathMark(TEvePathMark::kDecay);
    startPoint->fV.Set(vbuf);
    startPoint->fP.Set(pbuf);

    endParam.GetXYZ(vbuf);
    endPoint->fV.Set(vbuf);
    cout<<"endPoint = "<<vbuf[0]<<" "<<vbuf[1]<<" "<<vbuf[2]<<endl;
   
    { // get momentum manually because trackParam->GetPxPyPz doesn't works for straight lines

      Double_t pt= TMath::Abs(endParam.GetSigned1Pt());
      pt = (pt>kAlmost0) ?1./pt :100.;
      
      Double_t cA=TMath::Cos(endParam.GetAlpha()), sA=TMath::Sin(endParam.GetAlpha());
      
      Double_t sT=endParam.GetSnp();
      if( sT>=kAlmost1 ){ sT = kAlmost1; }
      else if( sT<-kAlmost1 ){ sT = -kAlmost1; }
      Double_t cT = TMath::Sqrt(TMath::Abs(1 - sT*sT));
      
      pbuf[0] = pt*(cT*cA - sT*sA); 
      pbuf[1] = pt*(sT*cA + cT*sA); 
      pbuf[2] = pt*endParam.GetTgl();
    }

    endPoint->fP.Set(pbuf);
    
    track->AddPathMark( startPoint );
     track->AddPathMark( endPoint );
  }
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

  AliESDEvent* esd = new AliESDEvent;
  esd->ReadFromTree(tr);
  tr->GetEntry(0);

  TEveTrackPropagator* rnrStyle = gTPCTrack->GetPropagator();
  rnrStyle->SetMagField( 0.1*esd->GetMagneticField() );

  cout << "Number of tracks found :" << esd->GetNumberOfTracks() << endl;

  for (Int_t ii=0; ii< esd->GetNumberOfTracks(); ii++) {

    AliESDtrack           *esdTrack = esd->GetTrack(ii);
    AliExternalTrackParam *trackParam = esdTrack;
    
    TEveTrack* track = makeESDTrack( rnrStyle, ii, esdTrack, trackParam );
    
    track->SetAttLineAttMarker(gTPCTrack);
    
    gEve->AddElement(track, gTPCTrack);
  }
  
  //  delete esd;
  
  return 0;
}
//****************************************************************************
Int_t initTrackerMapping()
{
  gDimuTrackerMapping.clear();

  if (! AliMpCDB::LoadDDLStore(true)){
    cerr<<__FILE__<<": Failed to Load DDLStore specified for CDBPath "
	<<AliCDBManager::Instance()->GetDefaultStorage()<<endl;
    return kFALSE;
  }

  AliMpSegmentation *mpSegFactory = AliMpSegmentation::Instance(); 
  AliGeomManager::LoadGeometry();
  AliMUONGeometryTransformer* chamberGeometryTransformer = new AliMUONGeometryTransformer();
  
  if(! chamberGeometryTransformer->LoadGeometryData()){
    cerr<<__FILE__<<": Failed to Load Geomerty Data "<<endl;
    return kFALSE;
  }

  AliMUONCalibrationData cd((AliCDBManager::Instance())->GetRun());
  AliHLTMUONTrackerMappingData md;

  cout<<"Loading Mapping for Dimuon Tracking Chambers.....be patient..."<<endl;

  for(Int_t iCh = 0; iCh < 10; iCh++){ 
    
    AliMpDEIterator it;
    for ( it.First(iCh); ! it.IsDone(); it.Next() ) {
    
      Int_t detElemId = it.CurrentDEId();
      
      for(Int_t iCath = 0 ; iCath <= 1 ; iCath++){
	
	AliMp::CathodType cath;

	if(iCath == 0)
	  cath = AliMp::kCath0 ;
	else
	  cath = AliMp::kCath1 ;
	
	const AliMpVSegmentation* seg = mpSegFactory->GetMpSegmentation(detElemId, cath);
	AliMp::PlaneType plane = seg->PlaneType(); 
	Int_t maxIX = seg->MaxPadIndexX();  
	Int_t maxIY = seg->MaxPadIndexY(); 
	UInt_t idManuChannel;
	Int_t buspatchId, manuId, channelId;
	Double_t realX, realY, realZ;
	Double_t localX, localY, localZ;

	//Pad Info of a segment
	for(Int_t iX = 0; iX<= maxIX ; iX++){
	  for(Int_t iY = 0; iY<= maxIY ; iY++){
	    if(seg->HasPad(AliMpIntPair(iX,iY))){
	      AliMpPad pad = seg->PadByIndices(AliMpIntPair(iX,iY),kFALSE);
	      

	      // Getting Manu id
	      manuId = pad.GetLocation().GetFirst();
	      manuId &= 0x7FF; // 11 bits 

	      // Getting channel id
	      channelId =  pad.GetLocation().GetSecond();
	      channelId &= 0x3F; // 6 bits

	      buspatchId = AliMpDDLStore::Instance()->GetBusPatchId(detElemId,manuId);
	      
	      idManuChannel &= 0x0;
	      idManuChannel = (idManuChannel|buspatchId)<<11;  
	      idManuChannel = (idManuChannel|manuId)<<6 ;
	      idManuChannel |= channelId ;
	      
	      localX = pad.Position().X();
	      localY = pad.Position().Y();
	      localZ = 0.0;

 	      chamberGeometryTransformer->Local2Global(detElemId,localX,localY,localZ,
 						       realX,realY,realZ);
	      md.fDetElemId = detElemId;
	      md.fIX = iX ;
	      md.fIY = iY ;
	      md.fCath = cath ;
	      md.fRealX = realX ;
	      md.fRealY = realY;
	      md.fRealZ = realZ;
  	      md.fPed = (cd.Pedestals(detElemId,manuId))->ValueAsFloat(channelId);
	      
	      gDimuTrackerMapping[idManuChannel] = md;
	      
	    }// HasPad condn
	  }// iY loop
	}// iX loop
	
      }// iPlane

    } // detElemId loop

  }// ichamber loop

  return 0;
}
//****************************************************************************
Int_t initTrackerMappingPed()
{
  cout<<"Loading Pedestal for Dimuon..."<<endl;
  Int_t linenum = 0;
  
  char line[90];
  int buspatch,manu,channel;
  float mean,sigma;
  int idManuChannel;

  FILE *fin[2] ;
  fin[0] = fopen("/local/home/hlt/pedestal/MUONTRKda_ped_25948.St1.ped","r");
  fin[1] = fopen("/local/home/hlt/pedestal/MUONTRKda_ped_25948.St2.ped","r");

  while(feof(fin[0])==0){
    fgets(line,100,fin[0]);
    sscanf(line,"%d\t%d\t%d\t%f\t%f\n",&buspatch,&manu,&channel,&mean,&sigma);
    if(linenum>12){
//       printf("%d\t%d\t%d\t%f\t%f\n",buspatch,manu,channel,mean,sigma);

      idManuChannel &= 0x0;
      idManuChannel = (idManuChannel|buspatch)<<11;  
      idManuChannel = (idManuChannel|manu)<<6 ;
      idManuChannel |= channel ;

//       data.fBuspatchId = UInt_t(buspatch) ;
//       data.fManuId = UInt_t(manu) ;
//       data.fChannelId = UInt_t(channel) ;
//       data.fADC = UShort_t(mean);
      
//       lookuptable[idManuChannel] = data;

      gDimuTrackerMapping[idManuChannel].fPed = mean;
      gDimuTrackerMapping[idManuChannel].fSigma = sigma;

    }
    linenum++;
  }

  linenum = 0;

  while(feof(fin[1])==0){
    fgets(line,100,fin[1]);
    sscanf(line,"%d\t%d\t%d\t%f\t%f\n",&buspatch,&manu,&channel,&mean,&sigma);
    if(linenum>12){
//       printf("%d\t%d\t%d\t%f\t%f\n",buspatch,manu,channel,mean,sigma);

      idManuChannel &= 0x0;
      idManuChannel = (idManuChannel|buspatch)<<11;  
      idManuChannel = (idManuChannel|manu)<<6 ;
      idManuChannel |= channel ;

//       data.fBuspatchId = UInt_t(buspatch) ;
//       data.fManuId = UInt_t(manu) ;
//       data.fChannelId = UInt_t(channel) ;
//       data.fADC = UShort_t(mean);
      
//       lookuptable[idManuChannel] = data;

      gDimuTrackerMapping[idManuChannel].fPed = mean;
      gDimuTrackerMapping[idManuChannel].fSigma = sigma;
    }
    linenum++;

  }

  return 0;
}
//****************************************************************************
Int_t initTriggerMapping()
{

  for (Int_t d = 0; d < 2; d++)
    for (Int_t i = 0; i < 8; i++)
      for (Int_t j = 0; j < 16; j++)
	for (Int_t k = 0; k < 4; k++)
	  for (Int_t n = 0; n < 2; n++)
	    for (Int_t m = 0; m < 16; m++){
	      gDimuTriggerMapping[d].fLut[i][j][k][n][m].fDetElemId = 0;
	      gDimuTriggerMapping[d].fLut[i][j][k][n][m].fX = 0;
	      gDimuTriggerMapping[d].fLut[i][j][k][n][m].fY = 0;
	      gDimuTriggerMapping[d].fLut[i][j][k][n][m].fZ = 0;
	    }
  
  
  AliCDBManager* cdbManager = AliCDBManager::Instance();
  cdbManager->SetRun((AliCDBManager::Instance())->GetRun());
  
  AliGeomManager::LoadGeometry();
	
  AliMUONGeometryTransformer transformer;
  if (! transformer.LoadGeometryData()){
    cerr << "ERROR: Could not load geometry into transformer." << endl;
    return false;
  }
  
  if (! AliMpCDB::LoadDDLStore()){
      cerr << "ERROR: Could not load DDL mapping." << endl;
      return false;
  }
  
  AliMpSegmentation* segmentation = AliMpSegmentation::Instance();
  if (segmentation == NULL){
    cerr << "ERROR: AliMpSegmentation::Instance() was NULL." << endl;
    return false;
  }


  AliMpDDLStore* ddlStore = AliMpDDLStore::Instance();
  if (ddlStore == NULL){
    cerr << "ERROR: AliMpDDLStore::Instance() was NULL." << endl;
    return false;
  }
	
  cout << "Loading Mapping for Dimuon Trigger Chambers....." << endl;
	
  AliMpDEIterator detElemIter;
  for (Int_t iDDL = 20; iDDL <= 21; iDDL++){

    for (Int_t iReg = 0; iReg < 8; iReg++){

      AliMpTriggerCrate* crate = ddlStore->GetTriggerCrate(iDDL, iReg);
      
      if (crate == NULL){
	cerr << "ERROR: Could not get crate for regional header = " << iReg
	     << ", and DDL ID = " << iDDL << endl;
	continue;
      }
			
      for (Int_t iLocBoard = 0; iLocBoard < 16; iLocBoard++){
	Int_t boardId = crate->GetLocalBoardId(iLocBoard);
	if (boardId == 0) continue;
				
	AliMpLocalBoard* localBoard = ddlStore->GetLocalBoard(boardId);
	if (localBoard == NULL){
	  cerr << "ERROR: Could not get loacl board: " << boardId << endl;
	  continue;
	}

	// skip copy cards
	if (! localBoard->IsNotified()) continue;
	
	for (Int_t iChamber = 0; iChamber < 4; iChamber++){

	  Int_t detElemId = ddlStore->GetDEfromLocalBoard(boardId, iChamber);
	  
	  const AliMUONGeometryDetElement* detElemTransform = transformer.GetDetElement(detElemId);
	  if (detElemTransform == NULL){
	    cerr << "ERROR: Got NULL pointer for geometry transformer for detection element ID = "
		 << detElemId << endl;
	    continue;
	  }
					
	  for (Int_t iCathode = 0; iCathode <= 1; iCathode++){
	    const AliMpVSegmentation* seg = segmentation->GetMpSegmentation
	      (detElemId, AliMp::GetCathodType(iCathode));
	    
	    for (Int_t bitxy = 0; bitxy < 16; bitxy++){
	      Int_t offset = 0;
	      if (iCathode && localBoard->GetSwitch(6)) offset = -8;
	      
	      AliMpPad pad = seg->PadByLocation(AliMpIntPair(boardId, bitxy+offset), kFALSE);
	      
	      if (! pad.IsValid()){
		  // There is no pad associated with the given local board and bit pattern.
		  continue;
	      }
	      
	      // Get the global coodinates of the pad.
	      Float_t lx = pad.Position().X();
	      Float_t ly = pad.Position().Y();
	      Float_t gx, gy, gz;
	      detElemTransform->Local2Global(lx, ly, 0, gx, gy, gz);
							
	      // Fill the LUT
	      gDimuTriggerMapping[(iDDL%20)].fLut[iReg][iLocBoard][iChamber][iCathode][bitxy].fDetElemId = detElemId;
	      gDimuTriggerMapping[(iDDL%20)].fLut[iReg][iLocBoard][iChamber][iCathode][bitxy].fX = gx;
	      gDimuTriggerMapping[(iDDL%20)].fLut[iReg][iLocBoard][iChamber][iCathode][bitxy].fY = gy;
	      gDimuTriggerMapping[(iDDL%20)].fLut[iReg][iLocBoard][iChamber][iCathode][bitxy].fZ = gz;

	    }// ibitxy loop
	  }// cathode loop
	}// chamber loop
      }// local board loop
    }// regional card loop
  }// ddl loop
  
	
  return true;
}
//****************************************************************************
Int_t initDiMuonMapping()
{
  initTrackerMapping();
  initTrackerMappingPed();
  initTriggerMapping();
  
  return 0;
}
//****************************************************************************

int ReadRawData(int iEvent = 0)
{
  cout<<"Executing for Event : "<<iEvent<<endl;
  
  TString rawDataPath = "$HOME/MUON_RawData";

  bool MakeTrackDDLStream(int*& rawData, int* rawDataSize, TString rawDataPath, int iEvent, int iDDL);  
  
  int *buffer;
  int bufferSize;

  AliMUONTrackerDDLDecoder<DiMuonTrackerCustomHandler> trackerDDLDecoder;
  DiMuonTrackerCustomHandler& handler = 
    reinterpret_cast<DiMuonTrackerCustomHandler&>(trackerDDLDecoder.GetHandler());
  
  for(Int_t iDDL=1;iDDL<=20;iDDL++){

    if(!MakeTrackDDLStream(buffer,&bufferSize,rawDataPath,iEvent,iDDL)){
      printf("Cannot Read DDL Stream, Check the Event number in RawData Directory or DDL number \n");
      continue ;
    }
    
    //       for(int i=0;i<int(bufferSize/sizeof(int));i++)
    // 	printf("rawData[%d] : %-8x :: %15d\n",i,buffer[i],buffer[i]);
    
    handler.ResetDataCounter();
    if(!trackerDDLDecoder.Decode(buffer, UInt_t(bufferSize))){
      cerr<<"Cannot decode. "<<endl;
      continue ;

    }
    
    for(Int_t i=0;i<handler.GetDataSize();i++){
      if(handler.GetData(i).fADC > 5){
	gDimuTrackerMapping[handler.GetData(i).fDataId].fADC = 
	  Int_t(handler.GetData(i).fADC) ;//- 
// 	  Int_t(gDimuTrackerMapping[handler.GetData(i).fDataId].fPed + 4*gDimuTrackerMapping[handler.GetData(i).fDataId].fSigma);

// 	cout<<"detElem : "<<gDimuTrackerMapping[handler.GetData(i).fDataId].fDetElemId
// 	    <<", iX : "<<gDimuTrackerMapping[handler.GetData(i).fDataId].fIX
// 	    <<", iY : "<<gDimuTrackerMapping[handler.GetData(i).fDataId].fIY
// 	    <<", realX : "<<gDimuTrackerMapping[handler.GetData(i).fDataId].fRealX
// 	    <<", realY : "<<gDimuTrackerMapping[handler.GetData(i).fDataId].fRealY
// 	    <<", realZ : "<<gDimuTrackerMapping[handler.GetData(i).fDataId].fRealZ
// 	    <<", plane : "<<gDimuTrackerMapping[handler.GetData(i).fDataId].fCath
// 	    <<", ADC : "<<gDimuTrackerMapping[handler.GetData(i).fDataId].fADC
// 	    <<endl;

	gDimuTrackerDataList.push_back(handler.GetData(i).fDataId);
      }//fADC > refADC value
      
    }
    
    delete []buffer ;
  }// DDL loop

  cout<<"Tracker Data size : "<<gDimuTrackerDataList.size()<<endl;
  return kTRUE;
}

//****************************************************************************

bool MakeTrackDDLStream(int*& rawData, int* rawDataSize, TString rawDataPath, int iEvent, int iDDL)
{
  char rawDataFile[500];
  sprintf(rawDataFile,"%s/raw%d/MUONTRK_%d.ddl",gSystem->ExpandPathName(rawDataPath.Data()),iEvent,0xA00 + iDDL - 1);

  
  FILE *fin = fopen(rawDataFile,"r");
  if(fin == NULL){
    printf("Failed to open file %s\n",rawDataFile);
    return false;
  }else{
//     printf("Opening file %s\n",rawDataFile);
    
    int ddlHeader[8];
    int rawDDLSize;

    if(feof(fin)==0){
      fread((char *) & ddlHeader,(size_t)4*8,1,fin);
      rawDDLSize = ddlHeader[0]/sizeof(int) - 8;
//       rawDDLSize = ddlHeader[0]/sizeof(int) - 8 - 2 ; // temporary solution
      //rewind(fin);
      
//       cout<<"rawDDLSize : "<<rawDDLSize<<endl;
//        for(int i=0;i<8;i++)
// 	 printf("ddlHeader[%d] : %d\n",i,ddlHeader[i]);
      
      if(ddlHeader[0]>10){ // 10 is temporary inprinciple 8 should be enough
	int* buffer = new int[rawDDLSize];
	fread(buffer,sizeof(int), rawDDLSize,fin);
	if(UInt_t(buffer[0])==0xFC0000FC){
	  rawData = buffer;
	  *rawDataSize = rawDDLSize*sizeof(int) ;
	  //cout<<"buffer : "<<buffer<<endl;
	}else{
	  fclose(fin);
	  return false;
	}
      }else{
	fclose(fin);
	return false;
      }
    }
    
  }// 
  
  fclose(fin);
  return true;
}

//****************************************************************************
Int_t ReadPort( AliHLTHOMERBlockDesc* block )
{
  int *buffer;
  int bufferSize;

  AliMUONTrackerDDLDecoder<DiMuonTrackerCustomHandler> trackerDDLDecoder;
  DiMuonTrackerCustomHandler& handler = 
    reinterpret_cast<DiMuonTrackerCustomHandler&>(trackerDDLDecoder.GetHandler());
  
  DiMuonTriggerDDLDecoder triggerDecoder;

  unsigned long size = block->GetSize();
  bufferSize = UInt_t(int(size/sizeof(int))-10);
  buffer = reinterpret_cast<int *>(block->GetData());
  buffer += 8;

  //   if(buffer[0]==0xFC0000FC){
  //     for(int i=0;i<int(bufferSize);i++)
  //       printf("rawData[%d] : %-8x :: %15d\n",i,buffer[i],buffer[i]);
  //   }

  if(buffer[0]!=int(0xFC0000FC)){
    
//     for(int i=0;i<int(bufferSize);i++)
//       printf("rawData[%d] : %-8x :: %15d\n",i,buffer[i],buffer[i]);
    
    triggerDecoder.ResetDataCounter();
    triggerDecoder.SetTriggerMappingData(&gDimuTriggerMapping[0]);
    if(!triggerDecoder.Decode(buffer)){
      cerr<<"Cannot decode DiMuon Trigger RawData "<<endl;
      //return kFALSE;
    }
  
    for(Int_t i=0;i<triggerDecoder.GetDataSize();i++)
      gDimuTriggerDataList.push_back(triggerDecoder.GetData(i));
    
      
  }else{

    handler.ResetDataCounter();
    if(!trackerDDLDecoder.Decode(buffer,UInt_t(bufferSize*sizeof(int)))){
      cerr<<"Cannot decode DiMuon Tracker RawData "<<endl;
      //return kFALSE;
    }
  
    for(Int_t i=0;i<handler.GetDataSize();i++){
      if(handler.GetData(i).fADC > 0){
	
	gDimuTrackerMapping[handler.GetData(i).fDataId].fADC = 
	  Int_t(handler.GetData(i).fADC) - 
	  Int_t(gDimuTrackerMapping[handler.GetData(i).fDataId].fPed);// + 4*gDimuTrackerMapping[handler.GetData(i).fDataId].fSigma);
	
	if(gDimuTrackerMapping[handler.GetData(i).fDataId].fADC > 0)
	  gDimuTrackerDataList.push_back(handler.GetData(i).fDataId);
	
	
      }// adc cut
    }// tracker dataSize
  }// if trigger/tracker data
  
  
  return true;
}
//****************************************************************************
Int_t processDiMuonRawData( AliHLTHOMERBlockDesc* block )
{
  if(!ReadRawData(0))
    return false;

//   if(!ReadPort(block))
//     return false;

  return true ;
}
//****************************************************************************
Int_t drawDiMuonRawData(Int_t ADCCut) 
{

  cout<<"Trying to read RawData "<<endl;


  char ChamberName[50];
  TEveElement *eveMuonEventElmt = gEve->GetCurrentEvent()->FindChild("MUONChamberData");
  TEveElement *eveMuonGeomElmt =  gEve->GetGlobalScene()->FindChild("MUONChambers");
  
  eveMuonEventElmt->DestroyElements();

  AliEveMUONChamber* mch = 0;
  int ipoint = 0 ;


  for(Int_t ichamber = 0 ; ichamber < 14; ichamber++){

    if(ichamber<10)
      sprintf(ChamberName,"Chamber 0%d (trac)",ichamber);
    else
      sprintf(ChamberName,"Chamber %2d (trig)",ichamber);

    mch = (AliEveMUONChamber*)(eveMuonGeomElmt->FindChild(ChamberName));
			       
//     cout<<"ichamber : "<<ichamber<<", mch : "<<mch<<endl;
    mch->GetChamberData()->DropData();
  }
  
  
  // -- DrawTracker RawData
  Int_t dataId,ich, busPatchId, manuId ;
  
  cout<<"Total number of Tracker data : "<<Int_t(gDimuTrackerDataList.size())<<endl;

  for(Int_t idata = 0 ; idata < Int_t(gDimuTrackerDataList.size()); idata++){

    dataId = gDimuTrackerDataList.at(idata);
    cout<<"ADC : "<<gDimuTrackerMapping[dataId].fADC<<endl;
    if(gDimuTrackerMapping[dataId].fADC > ADCCut){

      if(Int_t(gDimuTrackerMapping[dataId].fDetElemId)<100 or Int_t(gDimuTrackerMapping[dataId].fDetElemId) > 1025 ) 
      continue ;
      
      ich = Int_t(gDimuTrackerMapping[dataId].fDetElemId/100) - 1;
      sprintf(ChamberName,"Chamber 0%d (trac)",ich);
      
      mch = (AliEveMUONChamber*)(eveMuonGeomElmt->FindChild(ChamberName));
      
      //if(ich==0) continue ;

      manuId =  (dataId>>6) & 0x7FF ;
      busPatchId =  (dataId>>17) & 0xFFF ;
      
      if(busPatchId==728 and manuId==1122) continue ;

      cout<<"ich : "<<ich<<", idata : "<<idata
	  <<", detElemId : "<<gDimuTrackerMapping[dataId].fDetElemId
	  <<" bus : "<<busPatchId
	  <<" manu : "<<manuId
	  <<" iX : "<<gDimuTrackerMapping[dataId].fIX
	  <<" iY : "<<gDimuTrackerMapping[dataId].fIY
	  <<" realX : "<<gDimuTrackerMapping[dataId].fRealX
	  <<" realY : "<<gDimuTrackerMapping[dataId].fRealY
	  <<" realZ : "<<gDimuTrackerMapping[dataId].fRealZ
	  <<" ADC : "<<gDimuTrackerMapping[dataId].fADC
	  <<endl;
      
      if(gDimuTrackerMapping[dataId].fCath == 0)
	mch->GetChamberData()->RegisterDigit(gDimuTrackerMapping[dataId].fDetElemId,
					     0,
					     gDimuTrackerMapping[dataId].fIX,
					     gDimuTrackerMapping[dataId].fIY,
					     gDimuTrackerMapping[dataId].fADC);
      
      mch->GetChamberData()->RegisterHit(gDimuTrackerMapping[dataId].fDetElemId,
					 gDimuTrackerMapping[dataId].fRealX,
					 gDimuTrackerMapping[dataId].fRealY,
					 gDimuTrackerMapping[dataId].fRealZ);

      TEvePointSet* ps = new TEvePointSet(Form("Point-%d-Ch-0%d",ipoint,ich),1);
      ps->SetPoint(0,
		   gDimuTrackerMapping[dataId].fRealX,
		   gDimuTrackerMapping[dataId].fRealY,
		   gDimuTrackerMapping[dataId].fRealZ);
      
      //point3d[ipoint].SetElementName();
      //point3d[ipoint].SetRnrState(kFALSE);
      //point3d[ipoint].SetMarkerSize(15);
      //point3d[ipoint].SetMarkerColor(4);
      //point3d[ipoint].SetMarkerStyle(4);
      
      //gEve->AddElement(&point3d[ipoint],eveMuonEventElmt);

      ps->SetMarkerSize(15);
      ps->SetMarkerColor(4);
      ps->SetMarkerStyle(4);
      
      eveMuonEventElmt->AddElement(ps);
      ipoint++;

    }// ADC Cut
  }//foor loop
  
  // -- DrawTrigger RawData
//   cout<<"Total number of Trigger data : "<<Int_t(gDimuTriggerDataList.size())<<endl;
//   AliHLTMUONTriggerPointData data;
//   for(Int_t idata = 0 ; idata < Int_t(gDimuTriggerDataList.size()); idata++){

//     data = gDimuTriggerDataList.at(idata) ;
    
//     if(Int_t(data.fDetElemId)<1100 or Int_t(data.fDetElemId) > 1417 ) 
//       continue ;

//     ich = Int_t(data.fDetElemId/100) - 1;
//     sprintf(ChamberName,"Chamber %2d (trig)",ich);

//     mch = (AliEveMUONChamber*)(eveMuonGeomElmt->FindChild(ChamberName));

      
//     cout<<"ich : "<<ich<<", idata : "<<idata
// 	<<", detElemId : "<<data.fDetElemId
// 	  <<" realX : "<<data.fX
// 	  <<" realY : "<<data.fY
// 	  <<" realZ : "<<data.fZ
// 	  <<endl;
    
//     mch->GetChamberData()->RegisterHit(data.fDetElemId,data.fX,data.fY,data.fZ);

    
//     point3d[ipoint].SetPoint(0,data.fX,data.fY,data.fZ);


//     sprintf(ChamberName,"Point-%d-Ch-0%d",ipoint,ich);
//     point3d[ipoint].SetName(ChamberName);
//     point3d[ipoint].SetRnrState(kFALSE);
//     point3d[ipoint].SetMarkerSize(15);
//     point3d[ipoint].SetMarkerColor(4);
//     point3d[ipoint].SetMarkerStyle(4);
    
//     gEve->AddElement(&point3d[ipoint],eveMuonEventElmt);
//     ipoint++;

//   }// data loop


  globMaxPoint = ipoint ;



  return 0;
}
//****************************************************************************
void loopEvent() {
  event_timer.SetCommand("nextEvent()");
  event_timer.Start(60);
}

//****************************************************************************
void stopLoopEvent() {
  event_timer.Stop();
}

