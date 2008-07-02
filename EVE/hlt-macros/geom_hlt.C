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

Int_t MUON_geom()
{

  AliEveMUONChamber*   mucha[14];


  AliEveMUONData * g_muon_data = new AliEveMUONData;

  gStyle->SetPalette(1, 0);

  gEve->DisableRedraw();

  TEveElementList* mul = new TEveElementList("MUONChambers");
  TEveElementList* muChData = new TEveElementList("MUONChamberData");
  mul->SetTitle("MUON chambers");
  mul->SetMainColor(3);
  gEve->AddGlobalElement(mul);
  gEve->AddElement(muChData);
  
  for (Int_t ic = 0; ic < 14; ic++){
    
    mucha[ic] = new AliEveMUONChamber(ic);
    
    mucha[ic]->SetFrameColor(3);
    mucha[ic]->SetChamberID(ic);

    mucha[ic]->SetDataSource(g_muon_data);
    
    gEve->AddElement(mucha[ic], mul);
    
  }

  gEve->Redraw3D(kTRUE);
  gEve->EnableRedraw();

  
  return true;
}
