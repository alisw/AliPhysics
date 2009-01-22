// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

// Macro to visualise rootified raw-data from TPC.
//
// Use tpc_raw(Int_t mode) in order to run it
// Needs that alieve_init() is already called
// mode = 1 - show only 2D sectors
// mode = 2 - show only 3D sectors
// mode = 3 - show both 2D and 3D sectors

void tpc_raw(Int_t mode = 3)
{
  gStyle->SetPalette(1, 0);

  gEve->DisableRedraw();

  AliRawReader *reader = AliEveEventManager::AssertRawReader();
  reader->Reset();
  AliTPCRawStream input(reader);
  reader->Select("TPC"); // ("TPC", firstRCU, lastRCU);

  AliEveTPCData *x = new AliEveTPCData;
  // x->SetLoadPedestal(5);
  x->SetLoadThreshold(5);
  x->SetAutoPedestal(kTRUE);

  x->LoadRaw(input, kTRUE, kTRUE);

  TEveElementList* sec2d = new TEveElementList("TPC 2D");
  gEve->AddElement(sec2d);

  TEveElementList* sec3d = new TEveElementList("TPC 3D");
  gEve->AddElement(sec3d);

  for (Int_t i=0; i<=35; ++i) {
    if (mode & 1) {
      s = new AliEveTPCSector2D(Form("2D sector %d",i));
      s->SetSectorID(i);
      s->SetAutoTrans(kTRUE); // place on proper 3D coordinates
      s->SetDataSource(x);
      s->SetFrameColor(36);
      sec2d->AddElement(s);
      s->IncRTS();
    }
    if (mode & 2) {
      t = new AliEveTPCSector3D(Form("3D sector %d",i));
      t->SetSectorID(i);
      t->SetAutoTrans(kTRUE);
      t->SetDataSource(x);
      t->SetMinTime(40);
      t->SetMaxTime(980);
      t->SetDriftVel(2.273);
      sec3d->AddElement(t);
      t->IncRTS();
    }
  }

  gEve->EnableRedraw();
  gEve->Redraw3D();
}
