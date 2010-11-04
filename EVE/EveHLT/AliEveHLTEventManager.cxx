#include "TEveManager.h"
#include "TEveScene.h"
#include "TEveProjectionManager.h"
#include "TEveBrowser.h"
#include "TGLViewer.h"
#include "TEveViewer.h"
#include "TEveEventManager.h"

//#include "AliHLTTriggerDecision.h"
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

#include "AliEveHLTEventManager.h"
#include "AliEveHOMERManager.h"
#include "AliEveEventBuffer.h"

#include "TList.h"
#include "TTimer.h"

#include "TThread.h"

ClassImp(AliEveHLTEventManager);

AliEveHLTEventManager::AliEveHLTEventManager() : 
  TEveElementList("Event Manager"),
  fGeoManager(NULL),
  fEveManager(NULL),
  fRPhiManager(NULL),
  fRhoZManager(NULL),
  fRPhiEventScene(NULL),
  fRhoZEventScene(NULL),
  fRhoZViewer(NULL),
  fRPhiViewer(NULL),
  fTimer(NULL),
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
  fShowMuon(kFALSE), 
  fRunNumber(-1)  
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
   
  fTimer = new TTimer();
  fTimer->Connect("Timeout()", "AliEveHLTEventManager", this, "NextEvent()" );

}
 
AliEveHLTEventManager::~AliEveHLTEventManager() {
  
  //DestroyElements();
  //DestroyDetectorElements();  
  
}


void AliEveHLTEventManager::DestroyDetectorElements(){
  //See header file for documentation

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

///_______________________________________________________________________________________
void AliEveHLTEventManager::ConnectEventBuffer() {
  GetEventBuffer()->ConnectToSource();
}


///___________________________________________________________________________________________
void AliEveHLTEventManager::StartBufferMonitor() { 
  AliEveEventBuffer * buffer = GetEventBuffer();
  buffer->StartBufferMonitor();
}

//______________________________________________________________________________________________
Int_t AliEveHLTEventManager::ProcessEvent(AliESDEvent * event) {

  //We have a new event, reset display items (need to check if there really is anything interesting in event before resetting. ie not just histos)

  gEve->DisableRedraw();

  
  cout << "reset()"<<endl;
  
  ResetDisplay();
  
  cout << "process()"<<endl;
  if(!fHLTElement) {
    fHLTElement = new AliHLTEveHLT();
    fHLTElement->SetEventManager(this);
    gEve->AddElement(fHLTElement);
 }
  fHLTElement->ProcessEsdEvent(event);

  if(!fPhosElement) CreatePhosElement();
  fPhosElement->ProcessEvent(event);
  
  if(!fEmcalElement) CreateEmcalElement();
  fEmcalElement->ProcessEvent(event);
  
  cout << "update()"<<endl;
  UpdateDisplay();

  gEve->Redraw3D(0, 1);
  gEve->EnableRedraw();

  return 0;

}



//______________________________________________________________________________________________
Int_t AliEveHLTEventManager::ProcessEvent(TList * blockList) {

  //We have a new event, reset display items (need to check if there really is anything interesting in event before resetting. ie not just histos)
  
  if(!blockList) {
    cout << "Block list is NULL pointer, return " << endl;
    return -1;
  }
 
  cout << "reset ()"<<endl;
  ResetDisplay();
  cout << "done reset()"<<endl;
  AliHLTHOMERBlockDesc * block = NULL;
  TIter next(blockList);
  while ((block = (AliHLTHOMERBlockDesc*)next())) {
    cout <<"Process Block"<<endl;
    ProcessBlock(block);
  } 

  cout << "update() "<<endl;
  UpdateDisplay();
  
  return 0;

}
///___________________________________________________________________________________________

void AliEveHLTEventManager::ProcessBlock(AliHLTHOMERBlockDesc * block) {
  //See header file for documentation
  
#if 1//DEBUG
  printf( "------------------- xxxxxxxxxxxxxxx ----------------------\n");
  printf( "Detector           : %s\n", block->GetDetector().Data() );
  printf( "Datatype           : %s\n", block->GetDataType().Data() );
  if (block->IsTObject() )
    printf( "Is TObject of class: %s\n", block->GetClassName().Data() );
  printf( "------------------- xxxxxxxxxxxxxxx ----------------------\n");
#endif


  if(fShowBarrel) {

    if ( ! block->GetDetector().CompareTo("PHOS") ) {
      if(!fPhosElement) CreatePhosElement();
      fPhosElement->ProcessBlock(block);
    }
    
    else if ( ! block->GetDetector().CompareTo("EMCA") ) {
      if(!fEmcalElement) CreateEmcalElement();
      fEmcalElement->ProcessBlock(block);
    }  
    
    else if ( ! block->GetDetector().CompareTo("TPC") ) {
      if(!fTPCElement) CreateTPCElement();
      fTPCElement->ProcessBlock(block);
    }
    
    else if ( ! block->GetDetector().CompareTo("HLT") ) {
      if(!fHLTElement) CreateHLTElement();
      fHLTElement->ProcessBlock(block);
    }

    else if ( ! block->GetDetector().CompareTo("ITS") ) {
      if(!fITSElement) CreateITSElement();
      fITSElement->ProcessBlock(block);
    }
    
    else if ( ! block->GetDetector().CompareTo("ISDD") ) {
      if(!fISDDElement) CreateISDDElement();
      fISDDElement->ProcessBlock(block);
    }
    
    else if ( ! block->GetDetector().CompareTo("ISPD") ) {
      if(!fISPDElement) CreateISPDElement();
      fISPDElement->ProcessBlock(block);
    }
    
    else if ( ! block->GetDetector().CompareTo("ISSD") ) {
      if(!fISSDElement) CreateISSDElement();
      fISSDElement->ProcessBlock(block);
    }
    
    else if ( ! block->GetDetector().CompareTo("TRD") ) {
      if(!fTRDElement) CreateTRDElement();
      fTRDElement->ProcessBlock(block);
    }
    
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
	gEve->AddElement(fMuonElement);
      }
      fMuonElement->ProcessBlock(block);
    }
  }

}

void AliEveHLTEventManager::ResetDisplay () {
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


void AliEveHLTEventManager::PrintScreens() {
//   //See header file for documentation

//   fEveManager->GetDefaultGLViewer()->SavePicture(Form("%d_0x%lu_3D.gif", fRunNumber, GetEventID(fCurrentEventIdx)));
//   fRhoZViewer->GetGLViewer()->SavePicture(Form("%d_0x%lu_RhoZ.gif", fRunNumber, GetEventID(fCurrentEventIdx)));
//   fRPhiViewer->GetGLViewer()->SavePicture(Form("%d_0x%lu_RPhi.gif", fRunNumber, GetEventID(fCurrentEventIdx)));
  return;
}


void AliEveHLTEventManager::StartLoop() {
  //See header file for documentation
  //fTimer->SetCommand("NextEvent()", "AliEveHLTEventManager", this);
  SetEventLoopStarted(kTRUE);
  fTimer->Start(10000);
}

void AliEveHLTEventManager::StopLoop() {
  //See header file for documentation
  fTimer->Stop();
  SetEventLoopStarted(kFALSE);
}


// void AliEveHLTEventManager::NavigateBack() {
  
//   if (fHomerManager->NavigateEventBufferBack()) {
//     //return -1;
//   } else {
    
//     TList * blockList = fHomerManager->GetBlockList();
//     if(blockList) {      
//       ProcessEvent(blockList);
//     } else {
//       cout << "BALLE Error blocklist NULL pointer even though it's navigateable"<<endl;
//     }
//   }   

// }

// void AliEveHLTEventManager::NavigateFwd() {

//   if (fHomerManager->NavigateEventBufferFwd()) {
//     cout << "No event available" << endl;
//     return;
//     //return -1;
//   } else {
//     cout << "Getting block list" << endl;
//     TList * blockList = fHomerManager->GetBlockList();
//     if (blockList){
//       ProcessEvent(blockList);
//     } else {
//       cout << "blockList is NULL pointer"<<endl;
//     }
//   }

// }

void  AliEveHLTEventManager::UpdateDisplay() {
  //See header file for documentation
  if(fPhosElement) fPhosElement->UpdateElements();
  if(fEmcalElement) fEmcalElement->UpdateElements();
  if(fTPCElement) fTPCElement->UpdateElements();
  if(fHLTElement) fHLTElement->UpdateElements();
  if(fITSElement) fITSElement->UpdateElements();
  if(fISSDElement) fISSDElement->UpdateElements();
  if(fISDDElement) fISDDElement->UpdateElements();
  if(fISPDElement) fISPDElement->UpdateElements();
  if(fTRDElement) fTRDElement->UpdateElements();
  if(fAnyElement) fAnyElement->UpdateElements();
  if(fMuonElement) fMuonElement->UpdateElements();


  // -- Set EventID in Window Title  
  TString winTitle("Eve Main Window");
  winTitle += Form("-- Run Number: %d", fRunNumber); 
  //winTitle += Form("-- Event ID : 0x%lu ", GetEventID(fCurrentEventIdx) );
  GetEveManager()->GetBrowser()->SetWindowName(winTitle);

  //==============================================================================
  // -- Import global scene into projection scenes
  //==============================================================================

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

void AliEveHLTEventManager::SaveEveryThing() {

  GetEventBuffer()->WriteToFile();
  //Save everything to file
  //fEventBuffer->SaveBlockList();
  //fEventBuffer->SaveAsyncBlockList();


}



void AliEveHLTEventManager::CreatePhosElement() {
  fPhosElement = new AliHLTEvePhos();
  fPhosElement->SetEventManager(this);
  gEve->AddElement(fPhosElement);
}

void AliEveHLTEventManager::CreateEmcalElement() {
  fEmcalElement = new AliHLTEveEmcal();
  fEmcalElement->SetEventManager(this);
  gEve->AddElement(fEmcalElement);
}
void AliEveHLTEventManager::CreateTPCElement() {
  fTPCElement = new AliHLTEveTPC();
  fTPCElement->SetEventManager(this);
  gEve->AddElement(fTPCElement);
}
void AliEveHLTEventManager::CreateITSElement() {
  fITSElement = new AliHLTEveITS();
  fITSElement->SetEventManager(this);
  gEve->AddElement(fITSElement);
}
void AliEveHLTEventManager::CreateISPDElement() {
  fISPDElement = new AliHLTEveISPD();
  fISPDElement->SetEventManager(this);
  gEve->AddElement(fISPDElement);
}
void AliEveHLTEventManager::CreateISDDElement() {
  fISDDElement = new AliHLTEveISDD();
  fISDDElement->SetEventManager(this);
  gEve->AddElement(fISSDElement);
}
void AliEveHLTEventManager::CreateISSDElement() {
  fISSDElement = new AliHLTEveISSD();
  fISSDElement->SetEventManager(this);
  gEve->AddElement(fISSDElement);
}
void AliEveHLTEventManager::CreateTRDElement() {
  fTRDElement = new AliHLTEveTRD();
  fTRDElement->SetEventManager(this);
  gEve->AddElement(fTRDElement);
}
void AliEveHLTEventManager::CreateHLTElement() {
  fHLTElement = new AliHLTEveHLT();
  fHLTElement->SetEventManager(this);
  gEve->AddElement(fHLTElement);
}


