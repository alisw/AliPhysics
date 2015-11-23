// $Id$
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
#endif

TEveGeoShape* geom_gentle(Bool_t register_as_global=kTRUE)
{
  TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_geo.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
  TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
  f.Close();

    // set PHOS's colors
    TEveElement* elPHOS = gsre->FindChild("PHOS");
    elPHOS->SetRnrState(kTRUE);
    elPHOS->FindChild("PHOS_5")->SetRnrState(kFALSE);
    
    elPHOS->FindChild("PHOS_1")->SetMainColor(593);
    elPHOS->FindChild("PHOS_2")->SetMainColor(593);
    elPHOS->FindChild("PHOS_3")->SetMainColor(593);
    elPHOS->FindChild("PHOS_4")->SetMainColor(593);
    elPHOS->FindChild("PHOS_5")->SetMainColor(593);
    
    elPHOS->FindChild("PHOS_1")->SetMainTransparency(70);
    elPHOS->FindChild("PHOS_2")->SetMainTransparency(70);
    elPHOS->FindChild("PHOS_3")->SetMainTransparency(70);
    elPHOS->FindChild("PHOS_4")->SetMainTransparency(70);
    elPHOS->FindChild("PHOS_5")->SetMainTransparency(70);
    
    // set TPC's color
    TEveElement *elTPC = gsre->FindChild("TPC");
    TEveElement *elTPC_M_1 = elTPC->FindChild("TPC_M_1");
    TEveElement *elTPC_Drift_1 = elTPC_M_1->FindChild("TPC_Drift_1");
    elTPC_Drift_1->SetMainColor(3);
    
    // set HMPID's color
    TEveElement* elHMPID = gsre->FindChild("HMPID");
    
    elHMPID->FindChild("HMPID_0")->SetMainColor(5);
    elHMPID->FindChild("HMPID_1")->SetMainColor(5);
    elHMPID->FindChild("HMPID_2")->SetMainColor(5);
    elHMPID->FindChild("HMPID_3")->SetMainColor(5);
    elHMPID->FindChild("HMPID_4")->SetMainColor(5);
    elHMPID->FindChild("HMPID_5")->SetMainColor(5);
    elHMPID->FindChild("HMPID_6")->SetMainColor(5);
    
    elHMPID->FindChild("HMPID_0")->SetMainTransparency(80);
    elHMPID->FindChild("HMPID_1")->SetMainTransparency(80);
    elHMPID->FindChild("HMPID_2")->SetMainTransparency(80);
    elHMPID->FindChild("HMPID_3")->SetMainTransparency(80);
    elHMPID->FindChild("HMPID_4")->SetMainTransparency(80);
    elHMPID->FindChild("HMPID_5")->SetMainTransparency(80);
    elHMPID->FindChild("HMPID_6")->SetMainTransparency(80);


  if (register_as_global)
  {
    gEve->AddGlobalElement(gsre);
  }

  return gsre;
}

TEveGeoShape* geom_gentle_rphi()
{
  // The resulting geometry is NOT added into the global scene!

  TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_rphi_geo.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
  TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
  f.Close();

    // set PHOS's colors
    TEveElement* elPHOS = gsre->FindChild("PHOS");
    elPHOS->SetRnrState(kTRUE);
    elPHOS->FindChild("PHOS_5")->SetRnrState(kFALSE);
    
    elPHOS->FindChild("PHOS_1")->SetMainColor(593);
    elPHOS->FindChild("PHOS_2")->SetMainColor(593);
    elPHOS->FindChild("PHOS_3")->SetMainColor(593);
    elPHOS->FindChild("PHOS_4")->SetMainColor(593);
    elPHOS->FindChild("PHOS_5")->SetMainColor(593);
    
    elPHOS->FindChild("PHOS_1")->SetMainTransparency(70);
    elPHOS->FindChild("PHOS_2")->SetMainTransparency(70);
    elPHOS->FindChild("PHOS_3")->SetMainTransparency(70);
    elPHOS->FindChild("PHOS_4")->SetMainTransparency(70);
    elPHOS->FindChild("PHOS_5")->SetMainTransparency(70);
    
    // set TPC's color
    TEveElement *elTPC = gsre->FindChild("TPC");
    TEveElement *elTPC_M_1 = elTPC->FindChild("TPC_M_1");
    TEveElement *elTPC_SSWHEEL_1 = elTPC_M_1->FindChild("TPC_SSWHEEL_1");
    TEveElement *elTPC_SSWSEC[18];
    TEveElement *elTPC_SWSEG[18];
    TEveElement *elTPC_SWS1[18];
    
    for(int i=0;i<18;i++)
    {
        elTPC_SSWSEC[i] = elTPC_SSWHEEL_1->FindChild(Form("TPC_SSWSEC_%d",i+1));
        elTPC_SWSEG[i] = elTPC_SSWSEC[i]->FindChild("TPC_SWSEG_1");
        elTPC_SWS1[i] = elTPC_SWSEG[i]->FindChild("TPC_SWS1_1");
        elTPC_SWS1[i]->SetMainColor(3);
        elTPC_SWS1[i]->SetMainTransparency(70);
    }
    
    // set HMPID's color
    TEveElement* elHMPID = gsre->FindChild("HMPID");
    
    elHMPID->FindChild("HMPID_0")->SetMainColor(5);
    elHMPID->FindChild("HMPID_1")->SetMainColor(5);
    elHMPID->FindChild("HMPID_2")->SetMainColor(5);
    elHMPID->FindChild("HMPID_3")->SetMainColor(5);
    elHMPID->FindChild("HMPID_4")->SetMainColor(5);
    elHMPID->FindChild("HMPID_5")->SetMainColor(5);
    elHMPID->FindChild("HMPID_6")->SetMainColor(5);
    
    elHMPID->FindChild("HMPID_0")->SetMainTransparency(80);
    elHMPID->FindChild("HMPID_1")->SetMainTransparency(80);
    elHMPID->FindChild("HMPID_2")->SetMainTransparency(80);
    elHMPID->FindChild("HMPID_3")->SetMainTransparency(80);
    elHMPID->FindChild("HMPID_4")->SetMainTransparency(80);
    elHMPID->FindChild("HMPID_5")->SetMainTransparency(80);
    elHMPID->FindChild("HMPID_6")->SetMainTransparency(80);
  return gsre;
}

TEveGeoShape* geom_gentle_rhoz()
{
  // The resulting geometry is NOT added into the global scene!

  TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_rhoz_geo.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
  TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
  f.Close();

    TEveElement* elPHOS = gsre->FindChild("PHOS");
    elPHOS->FindChild("PHOS_5")->SetRnrState(kFALSE);
    
    
    elPHOS->FindChild("PHOS_1")->SetMainColor(593);
    elPHOS->FindChild("PHOS_2")->SetMainColor(593);
    elPHOS->FindChild("PHOS_3")->SetMainColor(593);
    elPHOS->FindChild("PHOS_4")->SetMainColor(593);
    elPHOS->FindChild("PHOS_5")->SetMainColor(593);
    
    elPHOS->FindChild("PHOS_1")->SetMainTransparency(70);
    elPHOS->FindChild("PHOS_2")->SetMainTransparency(70);
    elPHOS->FindChild("PHOS_3")->SetMainTransparency(70);
    elPHOS->FindChild("PHOS_4")->SetMainTransparency(70);
    elPHOS->FindChild("PHOS_5")->SetMainTransparency(70);
    
    // set TPC's color
    TEveElement *elTPC = gsre->FindChild("TPC");
    TEveElement *elTPC_M_1 = elTPC->FindChild("TPC_M_1");
    TEveElement *elTPC_Drift_1 = elTPC_M_1->FindChild("TPC_Drift_1");
    elTPC_Drift_1->SetMainColor(3);
    
    // set HMPID's color
    TEveElement* elHMPID = gsre->FindChild("HMPID");
    
    elHMPID->FindChild("HMPID_0")->SetMainColor(5);
    elHMPID->FindChild("HMPID_1")->SetMainColor(5);
    elHMPID->FindChild("HMPID_2")->SetMainColor(5);
    elHMPID->FindChild("HMPID_3")->SetMainColor(5);
    elHMPID->FindChild("HMPID_4")->SetMainColor(5);
    elHMPID->FindChild("HMPID_5")->SetMainColor(5);
    elHMPID->FindChild("HMPID_6")->SetMainColor(5);
    
    elHMPID->FindChild("HMPID_0")->SetMainTransparency(80);
    elHMPID->FindChild("HMPID_1")->SetMainTransparency(80);
    elHMPID->FindChild("HMPID_2")->SetMainTransparency(80);
    elHMPID->FindChild("HMPID_3")->SetMainTransparency(80);
    elHMPID->FindChild("HMPID_4")->SetMainTransparency(80);
    elHMPID->FindChild("HMPID_5")->SetMainTransparency(80);
    elHMPID->FindChild("HMPID_6")->SetMainTransparency(80);
    
  return gsre;
}
