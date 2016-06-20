// $Id: geom_gentle.C 54257 2012-01-30 20:52:05Z quark $
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TEveManager.h>
#include <TEveGeoNode.h>
#include <TEveElement.h>
#include <TEveGeoShape.h>
#include <TEveGeoShapeExtract.h>
#include <AliITSUGeomTGeo.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TGeoTube.h>
#include <TSystem.h>
#endif

TEveGeoShape* GetItsUpgradeGeom(Bool_t bPrint=kTRUE);

TEveGeoShape* geom_gentle_itsU(Bool_t register_as_global=kTRUE)
{

  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeSim");

  TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_geo.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
  TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
  f.Close();

  TEveElement* elPHOS = gsre->FindChild("PHOS");
  elPHOS->SetRnrState(kTRUE);
  elPHOS->FindChild("PHOS_4")->SetRnrState(kFALSE);
  elPHOS->FindChild("PHOS_5")->SetRnrState(kFALSE);
 
  TEveElement* elITS = gsre->FindChild("ITS");
  elITS->SetRnrState(kFALSE);
  
    TEveElement* elTPC = gsre->FindChild("TPC");
    elTPC->SetRnrState(kFALSE);
    TEveElement* elPHOS = gsre->FindChild("PHOS");
    elPHOS->SetRnrState(kFALSE);
    TEveElement* elTRD = gsre->FindChild("TRD+TOF");
    elTRD->SetRnrState(kFALSE);
    TEveElement* elHMPID = gsre->FindChild("HMPID");
    elHMPID->SetRnrState(kFALSE);
  

  TEveGeoShape* elITSU = GetItsUpgradeGeom(); 
  elITSU->SetRnrState(kTRUE);
  gsre->AddElement(elITSU);

  if (register_as_global)
  {
    gEve->AddGlobalElement(gsre);
  }

  return gsre;
}

TEveGeoShape* geom_gentle_itsU_rphi()
{
  // The resulting geometry is NOT added into the global scene!

  TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_rphi_geo.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
  TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
  f.Close();

  TEveElement* elPHOS = gsre->FindChild("PHOS");
  elPHOS->SetRnrState(kTRUE);
  elPHOS->FindChild("PHOS_4")->SetRnrState(kFALSE);
  elPHOS->FindChild("PHOS_5")->SetRnrState(kFALSE);

  TEveElement* elITS = gsre->FindChild("ITS");
  elITS->SetRnrState(kFALSE);
 
  TEveGeoShape* elITSupgr = GetItsUpgradeGeom(kFALSE); 
  elITSupgr->SetRnrState(kTRUE);
  gsre->AddElement(elITSupgr);

  return gsre;
}

TEveGeoShape* geom_gentle_itsU_rhoz()
{
  // The resulting geometry is NOT added into the global scene!

  TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_rhoz_geo.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
  TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
  f.Close();
 
  TEveElement* elITS = gsre->FindChild("ITS");
  elITS->SetRnrState(kFALSE);
 
  TEveGeoShape* elITSupgr = GetItsUpgradeGeom(kFALSE); 
  elITSupgr->SetRnrState(kTRUE);
  gsre->AddElement(elITSupgr);

  return gsre;
}



TEveGeoShape* GetItsUpgradeGeom(Bool_t bPrint) {
  
  gGeoManager = gEve->GetGeometry("geometry.root");

  AliITSUGeomTGeo* gm = new AliITSUGeomTGeo(kTRUE);
  TGeoVolume *itsV = gGeoManager->GetVolume(gm->GetITSVolPattern());
  if (!itsV) printf("ITS volume %s is not in the geometry\n",gm->GetITSVolPattern());



  // Re-create the Volume
  // --------------------------------------------

  TEveGeoShape *itsU = new TEveGeoShape("ITSU");  
 
  // Loop on all ITSV nodes, count Layer volumes by checking names
  Int_t nLayers=gm->GetNLayers();
  for (Int_t l=0; l<nLayers; l++) {
      Int_t id=gm->GetFirstChipIndex(l);
      Double_t loc[3]={0.,0.,0.}, glo[3];
      gm->LocalToGlobal(id,loc,glo);
      Double_t rMid=TMath::Sqrt(glo[0]*glo[0] + glo[1]*glo[1]);
      Double_t dz=TMath::Abs(glo[2]); 
      
      // Add Element to EVE
      TGeoTube *layer = new TGeoTube(rMid-0.25,rMid+0.25,dz+1.5);
      TEveGeoShape *elLayer = 
      new TEveGeoShape(Form("%s%d",gm->GetITSLayerPattern(),l));
      elLayer->SetShape(layer);
      elLayer->SetMainColor(kGreen);
      elLayer->SetMainTransparency(70/*90*/);
      itsU->AddElement(elLayer);
  }

  itsU->SetMainColor(kGreen);
  
  return itsU;


}
