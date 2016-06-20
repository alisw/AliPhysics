// Author: Stefano Carrazza 2010

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TEveElement.h>

#include <AliEveBeamsInfo.h>
#include <AliEveEventManager.h>
#endif

AliEveBeamsInfo* beams_info(){

  AliEveEventManager *mng = AliEveEventManager::Instance();
  AliEveBeamsInfo *beamsinfo = dynamic_cast<AliEveBeamsInfo *>(mng->FindGlobal("BeamsInfo"));

  if ( beamsinfo == 0) {
     beamsinfo = new AliEveBeamsInfo();
     mng->InsertGlobal("BeamsInfo", beamsinfo);
   } else {
     beamsinfo->Update();
   }

  return beamsinfo;

}
