//-*- Mode: C++ -*-

// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007              *
// Author: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>                *
// Author: 2010 Svein Lindal <slindal@fys.uio.no>                        *
//         for The ALICE HLT Project.                                    *



/** @file   AliEveHOMERManager.cxx
    @author Jochen Thaeder,  Svein Lindal <slindal@fys.uio.no>
    @date
    @brief  Manager for HOMER online
*/

#if __GNUC__>= 3
   using namespace std;
#endif

#include "unistd.h"
#include "TEveManager.h"
#include "TTimer.h"
#include "TEveScene.h"
#include "TEveProjectionManager.h"
#include "TEveBrowser.h"
#include "TGLViewer.h"
#include "TEveViewer.h"

#include "AliEveHOMERManager.h"
#include "AliHLTHOMERBlockDesc.h"
#include "AliHLTHOMERManager.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTEvePhos.h"
#include "AliHLTEveEmcal.h"
#include "AliHLTEveTPC.h"
#include "AliHLTEveHLT.h"
#include "AliHLTEveITS.h"
#include "AliHLTEveISPD.h"
#include "AliHLTEveISSD.h"
#include "AliHLTEveISDD.h"
#include "AliHLTEveTRD.h"
#include "AliHLTEveMuon.h"
#include "AliHLTEveAny.h"
#include "AliEveHOMERSourceList.h"

//#include "TTimer.h"

ClassImp(AliEveHOMERManager)
  
/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */
  
//##################################################################################
AliEveHOMERManager::AliEveHOMERManager() :
TEveElementList("Homer Manager"),
  AliHLTHOMERManager(), 
  fSrcList(NULL),
  fRetryCount(1000),
  fRetrySleeptime(15),
  fGeoManager(NULL),
  fEveManager(NULL),
  fRPhiManager(NULL),
  fRhoZManager(NULL),
  fRPhiEventScene(NULL),
  fRhoZEventScene(NULL),
  fRhoZViewer(NULL),
  fRPhiViewer(NULL),
  fTimer(NULL),
//  fSourceListTimer(NULL),
  fPhosElement(NULL), 
  fEmcalElement(NULL), 
  fTPCElement(NULL),
  fHLTElement(NULL),
  fITSElement(NULL),
  fISPDElement(NULL),
  fISSDElement(NULL),
  fISDDElement(NULL),
  fTRDElement(NULL),
  fMuonElement(NULL),
  fAnyElement(NULL),
  fEventLoopStarted(kFALSE),
  fCenterProjectionsAtPrimaryVertex(kFALSE),
  fShowBarrel(kTRUE),
  fShowMuon(kFALSE) 
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
   
  fTimer = new TTimer();
  fTimer->Connect("Timeout()", "AliEveHOMERManager", this, "NextHOMEREvent()" );

  //fSourceListTimer = new TTimer();

  fPhosElement = new AliHLTEvePhos();
  fPhosElement->SetEventManager(this);
   
  fEmcalElement = new AliHLTEveEmcal();
  fEmcalElement->SetEventManager(this);
   
  fTPCElement = new AliHLTEveTPC();
  fTPCElement->SetEventManager(this);
   
  fHLTElement = new AliHLTEveHLT();
  fHLTElement->SetEventManager(this);
   
  fITSElement = new AliHLTEveITS();
  fITSElement->SetEventManager(this);
  
  fISPDElement = new AliHLTEveISPD();
  fISPDElement->SetEventManager(this);

  fISDDElement = new AliHLTEveISDD();
  fISDDElement->SetEventManager(this);

  fISSDElement = new AliHLTEveISSD();
  fISSDElement->SetEventManager(this);

  fTRDElement = new AliHLTEveTRD();
  fTRDElement->SetEventManager(this);

}

//##################################################################################
AliEveHOMERManager::~AliEveHOMERManager() {
 // see header file for class documentation 

  DestroyElements();
  DestroyDetectorElements();

}



void AliEveHOMERManager::DestroyDetectorElements(){
  //See header file for documentation
  if (fSrcList)
    delete fSrcList;
  fSrcList = NULL;

  if (fPhosElement)
    delete fPhosElement;
  fPhosElement = NULL;

  if(fEmcalElement)
    delete fEmcalElement;
  fEmcalElement = NULL;

  if(fTPCElement)
    delete fTPCElement;
  fTPCElement = NULL;

  if(fHLTElement)
    delete fHLTElement;
  fHLTElement = NULL;

  if(fITSElement)
    delete fITSElement;
  fITSElement = NULL;

  if(fISSDElement)
    delete fISSDElement;
  fISSDElement = NULL;

  if(fISDDElement)
    delete fISDDElement;
  fISDDElement = NULL;

  if(fISPDElement)
    delete fISPDElement;
  fISPDElement = NULL;

  if(fTRDElement)
    delete fTRDElement;
  fTRDElement = NULL;

  if(fMuonElement)
    delete fMuonElement;
  fMuonElement = NULL;
 
  if(fAnyElement)
    delete fAnyElement;
  fAnyElement = NULL;
  

}



/*
 * ---------------------------------------------------------------------------------
 *                                 Source Handling
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
Int_t AliEveHOMERManager::CreateEveSourcesList() {
  // see header file for class documentation

  DestroyElements();

  Int_t iResult = CreateSourcesList();

  //  fStateHasChanged = kTRUE;
  
  HLTDebug(Form("iResult XXX %d", iResult));
  if ( iResult )
    return iResult;


  HLTDebug(Form("iResult %d", iResult));
  if (fSrcList) {
    HLTInfo(Form("delete source list", iResult));
    DestroyElements();
    //delete fSrcList;
    fSrcList = NULL;
    //fSrcList->Clear();
    HLTInfo(Form("cleared source list", iResult));
  }


  // -- Create new AliEVE sources list 
  if(!fSrcList){
    HLTInfo(Form("no source list", iResult));
    fSrcList = new AliEveHOMERSourceList("HLT Sources");
    fSrcList->SetManager(this);
    
    AddElement(fSrcList);
  }
  
  //HLTInfo(Form("createbytype", iResult));
  fSrcList->CreateByDet();
  
  HLTDebug(Form("Done createing source list %d", iResult));    
  
  
  
  return iResult;
}

//##################################################################################
Int_t AliEveHOMERManager::CreateEveSourcesListLoop() {
  // see header file for class documentation

  Int_t iResult = 0;

  for ( Int_t retry = 0; retry < fRetryCount ; retry++ ) {

    iResult = CreateEveSourcesList();

    if (!iResult) {
      HLTInfo("Source list successfully created.");
      break;
    }
    else if (iResult == 1) {
      HLTWarning( Form("Couldn't find active services, sleeping %d s before making attempt %d out of %d", fRetrySleeptime, retry, fRetryCount) ) ;
    }   
    else if (iResult == 2) {
      HLTWarning( Form("Services List empty, sleeping %d s before making new attempt.", fRetrySleeptime) ) ;
    }
    else {
      HLTError( Form("Other problem ... \n") ); 
      return iResult;
    } 

    //fSourceListTimer->Start(fRetrySleeptime, kFALSE);
    sleep(fRetrySleeptime);

  }

  if ( iResult ) {
    HLTWarning( Form("Couldn't find active services.") );
    return iResult;
  } 
  
  return iResult;
}

//##################################################################################
Int_t AliEveHOMERManager::ConnectEVEtoHOMER( TString detector ) {
  // see header file for class documentation

  HLTInfo("");
  fStateHasChanged = fSrcList->GetSelectedSources();
  HLTInfo(Form("has state changed % d", fStateHasChanged));
  return ConnectHOMER(detector);
}


//##################################################################################
Int_t AliEveHOMERManager::ReConnectHOMER( TString /*detector*/ ){
  // see header file for class documentation
  
  Int_t iResult = 0;

  DisconnectHOMER();
  iResult = CreateEveSourcesListLoop();
  HLTInfo("Created new source list, reconnect to HOMER");
  iResult = ConnectEVEtoHOMER();
  if ( iResult ) {
    HLTError(Form("Error reconnecting."));
  }

  return iResult;
}


//_____________________________________________________________________________________
Int_t AliEveHOMERManager::NextHOMEREvent() {
  //See header file for documentation  
  Int_t iResult = 0;

  
 
  if ( NextEvent() ) {
    
    HLTInfo("Failed getting next event, trying to reconnect");
    iResult = ReConnectHOMER();
    return NextHOMEREvent();
  }
 
  return ProcessEvent();
 
}


//______________________________________________________________________________________________
Int_t AliEveHOMERManager::ProcessEvent() {

  //We have a new event, reset display items (need to check if there really is anything interesting in event before resetting. ie not just histos)
  ResetDisplay();
 

  AliHLTHOMERBlockDesc * block = NULL;


  //Process the SYNCED block list
  if ( GetBlockList() == NULL) {
    printf ("onlineDisplay:   No regular BlockList ... \n");
    cout << endl;
    //return -1;
   
  } else {
    
    if (GetBlockList()->IsEmpty() ) {
      printf ("onlineDisplay:   No Sync Blocks in list ... \n");
      cout<<endl;
      //return;
    }  
   
    
    TIter next(GetBlockList());
    while ((block = (AliHLTHOMERBlockDesc*)next())) {
      ProcessBlock(block);
     
    } 
  }


  //Read out histograms and elements from detectors outside physics 1 partition
  TIter anext(GetAsyncBlockList());

  while ( (block = (AliHLTHOMERBlockDesc*)anext()) ) {
    ProcessBlock(block);
  }
  
  UpdateDisplay();

  return 0;

}

void  AliEveHOMERManager::UpdateDisplay() {
  //See header file for documentation
  fPhosElement->UpdateElements();
  fEmcalElement->UpdateElements();
  fTPCElement->UpdateElements();
  fHLTElement->UpdateElements();
  fITSElement->UpdateElements();
  fISSDElement->UpdateElements();
  fISDDElement->UpdateElements();
  fISPDElement->UpdateElements();
  fTRDElement->UpdateElements();
  if(fAnyElement) fAnyElement->UpdateElements();
  if(fMuonElement) fMuonElement->UpdateElements();


  // -- Set EventID in Window Title  
  TString winTitle("Eve Main Window -- Event ID : ");
  winTitle += Form("0x%016X ", GetEventID() );
  GetEveManager()->GetBrowser()->SetWindowName(winTitle);

  //==============================================================================
  // -- Import global scene into projection scenes
  //==============================================================================

  // XXX Primary vertex ... to be retrieved from the ESD
  Double_t x[3] = { 0, 0, 0 };
  
  TEveElement* top = GetEveManager()->GetCurrentEvent();
  
  if (fRPhiManager && top) {
    fRPhiEventScene->DestroyElements();
    if (fCenterProjectionsAtPrimaryVertex)
      fRPhiManager->SetCenter(x[0], x[1], x[2]);
    fRPhiManager->ImportElements(top, fRPhiEventScene);
  }
  
  if (fRhoZManager && top) {
    fRhoZEventScene->DestroyElements();
    if (fCenterProjectionsAtPrimaryVertex)
      fRhoZManager->SetCenter(x[0], x[1], x[2]);
    fRhoZManager->ImportElements(top, fRhoZEventScene);
  }


  //Redraw the display
  GetEveManager()->Redraw3D(0,1); // (0, 1)
  GetEveManager()->EnableRedraw(); 

}

void AliEveHOMERManager::ProcessBlock(AliHLTHOMERBlockDesc * block) {
  //See header file for documentation
  
#if 0//DEBUG
  printf( "------------------- xxxxxxxxxxxxxxx ----------------------\n");
  printf( "Detector           : %s\n", block->GetDetector().Data() );
  printf( "Datatype           : %s\n", block->GetDataType().Data() );
  if (block->IsTObject() )
    printf( "Is TObject of class: %s\n", block->GetClassName().Data() );
  printf( "------------------- xxxxxxxxxxxxxxx ----------------------\n");
#endif


  if(fShowBarrel) {

    if ( ! block->GetDetector().CompareTo("PHOS") ) 
      fPhosElement->ProcessBlock(block);
    
    else if ( ! block->GetDetector().CompareTo("EMCA") ) 
      fEmcalElement->ProcessBlock(block);
  
    else if ( ! block->GetDetector().CompareTo("TPC") ) 
      fTPCElement->ProcessBlock(block);
  
    else if ( ! block->GetDetector().CompareTo("HLT") ) 
      fHLTElement->ProcessBlock(block);
  
    else if ( ! block->GetDetector().CompareTo("ITS") ) 
      fITSElement->ProcessBlock(block);
  
    else if ( ! block->GetDetector().CompareTo("ISDD") ) 
      fISDDElement->ProcessBlock(block);
  
    else if ( ! block->GetDetector().CompareTo("ISPD") ) 
      fISPDElement->ProcessBlock(block);
  
    else if ( ! block->GetDetector().CompareTo("ISSD") ) 
      fISSDElement->ProcessBlock(block);

    else if ( ! block->GetDetector().CompareTo("TRD") ) 
      fTRDElement->ProcessBlock(block);
   
    else if ( ! block->GetDetector().CompareTo("MUON") ) {
      //Do Nothing
       

    } else {
      if(!fAnyElement) {
	fAnyElement = new AliHLTEveAny();
	fAnyElement->SetEventManager(this);
      } 
      fAnyElement->ProcessBlock(block);
    }
     
  }

   
  if(fShowMuon) {
    if ( ! block->GetDetector().CompareTo("MUON") ) {
      if(!fMuonElement) {
	fMuonElement = new AliHLTEveMuon();
	fMuonElement->SetEventManager(this);
      }
      fMuonElement->ProcessBlock(block);
    }
  }

}

void AliEveHOMERManager::ResetDisplay () {
  //See header file for documentation

 if(fPhosElement)
   fPhosElement->ResetElements();
 
 if(fEmcalElement)
   fEmcalElement->ResetElements();
 
 if(fTPCElement)
   fTPCElement->ResetElements();
 
 if(fHLTElement)
   fHLTElement->ResetElements();
 
 if(fITSElement)
   fITSElement->ResetElements();
 
 if(fISPDElement)
   fISPDElement->ResetElements();
 
 if(fISDDElement)
   fISDDElement->ResetElements();
 
 if(fISSDElement)
   fISSDElement->ResetElements();

 if (fTRDElement)
   fTRDElement->ResetElements();

 if(fAnyElement)
   fAnyElement->ResetElements();

 if(fMuonElement)
   fMuonElement->ResetElements();

}


void AliEveHOMERManager::PrintScreens() {
  //See header file for documentation

  fEveManager->GetDefaultGLViewer()->SavePicture(Form("0x%016X_3D.gif", GetEventID()));
  fRhoZViewer->GetGLViewer()->SavePicture(Form("0x%016X_RhoZ.gif", GetEventID()));
  fRPhiViewer->GetGLViewer()->SavePicture(Form("0x%016X_RPhi.gif", GetEventID()));

}

void AliEveHOMERManager::StartLoop() {
  //See header file for documentation
  //fTimer->SetCommand("NextEvent()", "AliEveHOMERManager", this);
  SetEventLoopStarted(kTRUE);
  fTimer->Start(10000);
}

void AliEveHOMERManager::StopLoop() {
  //See header file for documentation
  fTimer->Stop();
  SetEventLoopStarted(kFALSE);
}

// void AliEveHOMERManager::TimeOut(Int_t sec) {
  
//   TTimer t(sec*1000, kFALSE);

//   cout << "Start timer for " << sec << " seconds" << endl;

//   while (!t.CheckTimer(gSystem->Now()))
//     gSystem->Sleep(100);  // 100 ms sleep

//   cout << "Timed out." << endl;

// }
