// Author: Stefano Carrazza 2010

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
 
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TEveElement.h>

#include <AliEveEventManager.h>
#include <AliEveLego.h>
#endif

AliEveLego* lego(){

   AliEveEventManager *mng = AliEveEventManager::Instance();
   AliEveLego *lego = dynamic_cast<AliEveLego*>(mng->FindGlobal("LegoHisto2D"));

   if ( lego == 0) {
      lego = new AliEveLego();
      mng->InsertGlobal("LegoHisto2D",lego);
   } else {
      lego->Update();
   } 

   return lego;
}
