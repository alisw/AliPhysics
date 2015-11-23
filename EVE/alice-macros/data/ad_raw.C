/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

// Macro to visualise rootified raw-data from AD.
//

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TStyle.h>
#include <TEveManager.h>

#include <AliRawReader.h>
#include <AliEveEventManager.h>
#include <AliEveADModule.h>
#else
class AliRawReader;
class AliEveADModule;
#endif

void ad_raw(Int_t maxCharge = 1023, Bool_t showLegend = kFALSE)
{
  gStyle->SetPalette(1, 0);

  AliRawReader *reader = AliEveEventManager::AssertRawReader();
  reader->Reset();

  gEve->DisableRedraw();

  AliEveADModule* rawA = new AliEveADModule("AD_RAW_A", kTRUE, maxCharge, showLegend);
  rawA->LoadRaw(reader);


  AliEveADModule* rawC = new AliEveADModule("AD_RAW_C", kFALSE, maxCharge, showLegend);
  rawC->LoadRaw(reader);
  
  

  gEve->EnableRedraw();
}
