//-*- Mode: C++ -*-

// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007
// Author: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>                *
//         for The ALICE HLT Project.                                    *



/** @file   AliEveHOMERManager.cxx
    @author Jochen Thaeder
    @date
    @brief  Manger for HOMER in offline
*/

#if __GNUC__>= 3
   using namespace std;
#endif

#include "unistd.h"

#include "AliEveHOMERManager.h"


ClassImp(AliEveHOMERManager)
  
/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */
  
//##################################################################################
AliEveHOMERManager::AliEveHOMERManager() :
  AliHLTHOMERManager(), 
  TEveElementList("Homer Manager"),
  fSrcList(NULL),
  fRetryCount(1),
  fRetrySleeptime(10) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

//##################################################################################
AliEveHOMERManager::~AliEveHOMERManager() {
 // see header file for class documentation 

  if (fSrcList)
    delete fSrcList;
  fSrcList = NULL;
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
  
  fStateHasChanged = kTRUE;
  
  if ( iResult )
    return iResult;

  if (fSrcList)
    delete fSrcList;
  fSrcList = NULL;

  // -- Create new AliEVE sources list 
  fSrcList = new AliEveHOMERSourceList("HLT Sources");
  fSrcList->SetManager(this);
    
  AddElement(fSrcList);
  fSrcList->CreateByType();
    
  return iResult;
}

//##################################################################################
Int_t AliEveHOMERManager::CreateEveSourcesListLoop() {
  // see header file for class documentation

  Int_t iResult = 0;

  for ( Int_t retry = 0; retry < fRetryCount ; retry++ ) {
  
    iResult = CreateEveSourcesList();
    if (!iResult) 
      break;
    
    else if (iResult == 1) {
      HLTWarning( Form("Couldn't find active services, sleeping %d s", fRetryCount) ) ;
    }   
    else if (iResult == 2) {
      HLTWarning( Form("Services List empty, sleeping %d s", fRetryCount) ) ;
    }
    else {
      HLTError( Form("Other problem ... \n") ); 
      return iResult;
    } 
    
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

  fStateHasChanged = fSrcList->GetSelectedSources();
  
  return ConnectHOMER(detector);
}
