/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

// Macro to visualise rootified raw-data from VZERO.
//

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TStyle.h>
#include <TEveManager.h>

#include <AliRawReader.h>
#include <AliEveEventManager.h>
#include <AliEveVZEROModule.h>
#else
class AliRawReader;
class AliEveVZEROModule;
#endif

void vzero_raw()
{
    printf("*** RAW VZero ***");
    
  gStyle->SetPalette(1, 0);

  AliRawReader *reader = AliEveEventManager::AssertRawReader();
  reader->Reset();

  gEve->DisableRedraw();

  AliEveVZEROModule* rawA = new AliEveVZEROModule("VZERO_RAW_A", kTRUE);
  rawA->LoadRaw(reader);


  AliEveVZEROModule* rawC = new AliEveVZEROModule("VZERO_RAW_C", kFALSE);
  rawC->LoadRaw(reader);

  gEve->EnableRedraw();
}
