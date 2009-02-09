/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

/*
 *   Geometry as used for HLT 
 */

#include "AliEveMUONData.h"
#include "AliEveMUONChamber.h"
#include "AliCDBManager.h"

#include "TEveGeoShapeExtract.h"
#include "TEveGeoNode.h"
#include "TEveManager.h"
#include "TEveEventManager.h"
#include "TEveElement.h"

#include "TFile.h"
#include "TStyle.h"

TEveGeoShape* geom_hlt()
{
  TFile f("$ALICE_ROOT/EVE/alice-data/gentle_geo.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
  TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse, 0);
  gEve->AddGlobalElement(gsre);
  f.Close();

  TEveElement* elTRD = gsre->FindChild("TRD+TOF");
  elTRD->SetRnrState(kFALSE);

  TEveElement* elPHOS = gsre->FindChild("PHOS");
  elPHOS->SetRnrState(kFALSE);

  TEveElement* elHMPID = gsre->FindChild("HMPID");
  elHMPID->SetRnrState(kFALSE);

  Int_t MUON_geom();
  MUON_geom();

  return gsre;
}

//#if 0

Int_t MUON_geom()
{
  gStyle->SetPalette(1, 0);
  gEve->DisableRedraw();

  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);

  AliEveMUONData *g_muon_data = new AliEveMUONData;

  TEveElementList* l = new TEveElementList("MUONChambers");
  l->SetTitle("MUON chambers");
  l->SetMainColor(2);
  gEve->AddGlobalElement(l);

  for (Int_t ic = 0; ic < 14; ic++)
  {
    AliEveMUONChamber* mucha = new AliEveMUONChamber(ic);

    mucha->SetFrameColor(2);
    mucha->SetChamberID(ic);

    mucha->SetDataSource(g_muon_data);

    gEve->AddElement(mucha, l);
  }

  gEve->Redraw3D(kTRUE);
  gEve->EnableRedraw();



  
  return true;
}

//#endif
