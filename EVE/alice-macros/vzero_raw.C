/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

// Macro to visualise rootified raw-data from VZERO.
//

class AliRawReader;

class AliEveVZEROModule;

void vzero_raw()
{
  gStyle->SetPalette(1, 0);

  gEve->DisableRedraw();

  AliRawReader *reader = AliEveEventManager::AssertRawReader();
  reader->Reset();

  AliEveVZEROModule* rawA = new AliEveVZEROModule("VZERO_RAW_A", kTRUE);
  rawA->LoadRaw(reader);


  AliEveVZEROModule* rawC = new AliEveVZEROModule("VZERO_RAW_C", kFALSE);
  rawC->LoadRaw(reader);

  gEve->EnableRedraw();
  gEve->Redraw3D();
}
