//-*- Mode: C++ -*-

// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007              *
// Author: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>                *
// Author: 2010 Svein Lindal <slindal@fys.uio.no>                        *
//         for The ALICE HLT Project.                                    *



/** @file   AliEveHOMERManager.cxx
    @author Jochen Thaeder,  
    @author Svein Lindal <slindal@fys.uio.no>
    @date
    @brief  Manager for HOMER events
*/

#if __GNUC__>= 3
   using namespace std;
#endif

#include "unistd.h"

#include "AliHLTTriggerDecision.h"
#include "AliEveHOMERManager.h"
#include "AliHLTHOMERBlockDesc.h"
#include "AliHLTHOMERManager.h"
#include "AliEveHOMERSourceList.h"

#include "TTimer.h"

ClassImp(AliEveHOMERManager)
  
//____________________________________________________________________________________
AliEveHOMERManager::AliEveHOMERManager() :
TEveElementList("Homer Manager"),
  AliHLTHOMERManager(), 
  fSrcList(NULL),
  fRetryCount(1000),
  fRetrySleeptime(15),
  fSourceListTimer(NULL)
{
  fSourceListTimer = new TTimer();
  fSourceListTimer->Connect("Timeout()", "AliEveHOMERManager", this, "CreateEveSourcesListLoop()");
}

//____________________________________________________________________________________
AliEveHOMERManager::~AliEveHOMERManager() {
  // see header file for class documentation 
  
  if (fSrcList)
    delete fSrcList;
  fSrcList = NULL;
  
  


}

 //____________________________________________________________________________________
Int_t AliEveHOMERManager::CreateEveSourcesList() {
  // see header file for class documentation

  DestroyElements();

  Int_t iResult = CreateSourcesList();

  fStateHasChanged = kTRUE;
  
  HLTDebug(Form("iResult XXX %d", iResult));
  if ( iResult )
    return iResult;

 

  HLTDebug(Form("iResult %d", iResult));
  if (fSrcList) {
    HLTInfo("delete source list");
    DestroyElements();
    //delete fSrcList;
    fSrcList = NULL;
    //fSrcList->Clear();
    HLTInfo("cleared source list");
  }


  // -- Create new AliEVE sources list 
  if(!fSrcList){
    HLTInfo("no source list");
    fSrcList = new AliEveHOMERSourceList("HLT Sources");
    fSrcList->SetManager(this);
    AddElement(fSrcList);
  }
  
  //HLTInfo(Form("createbytype", iResult));
  fSrcList->CreateByDet();
  
  HLTDebug(Form("Done creating source list %d", iResult));    
    
  return iResult;
}

///_______________________________________________________________
void AliEveHOMERManager::StartEveSourceListLoop() {
  HLTInfo("Starting source list timer");
  fSourceListTimer->Start(5000); 
}
///_______________________________________________________________
void AliEveHOMERManager::StopEveSourceListLoop() {
  fSourceListTimer->Stop(); 
}


//________________________________________________________________###
Int_t AliEveHOMERManager::CreateEveSourcesListLoop() {
  // see header file for class documentation

  HLTInfo("Attempting to create source list");
  Int_t iResult = CreateEveSourcesList();

  if (!iResult) {
    HLTInfo("Source list successfully created.");
    StopEveSourceListLoop();
    HLTInfo("Conneting to sources");
    return ConnectEVEtoHOMER();
  }
  

  else if (iResult == 1) {
    HLTWarning( Form("Couldn't find active services,"));
  
  } else if (iResult == 2) {
    HLTWarning( Form("Services List empty, sleeping %d s before making new attempt.", fRetrySleeptime) ) ;
  
  } else {
    HLTError( Form("Other problem ... \n") ); 
  } 
  
  return iResult;
}

//________________________________________________________________
Int_t AliEveHOMERManager::ConnectEVEtoHOMER( TString detector ) {
  // see header file for class documentation
  HLTInfo("");
  fStateHasChanged = fSrcList->GetSelectedSources();
  HLTInfo(Form("has state changed % d", fStateHasChanged));
  return ConnectHOMER(detector);
}


//________________________________________________________________
Int_t AliEveHOMERManager::ReConnectHOMER( TString /*detector*/ ){
  // see header file for class documentation
  Int_t iResult = 0;
  if (Connected()) DisconnectHOMER();
  CreateEveSourcesListLoop();
  return iResult;
}

//_____________________________________________________________________________________
TList * AliEveHOMERManager::NextHOMEREvent() {
  //See header file for documentation  
  
  if(!Connected()) {
    HLTInfo("Homer is not connected, trying to reconnect!");
    ReConnectHOMER();
    return NULL;
  }
  
  
  if ( NextEvent() ) {
    HLTInfo("Failed getting next event, trying to reconnect");
    ReConnectHOMER();
    return NULL;
  }
  
  return GetBlockList();
  
}
