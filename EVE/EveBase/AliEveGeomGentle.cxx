//
//  AliEveGeomGentle.cxx
//  xAliRoot
//
//  Created by Jeremi Niedziela on 11/05/15.
//
//

#include "AliEveGeomGentle.h"

#include <TFile.h>
#include <TEveManager.h>
#include <TEveGeoNode.h>
#include <TEveElement.h>
#include <TEveGeoShapeExtract.h>
#include <TGLViewer.h>

#include <iostream>

using namespace std;

AliEveGeomGentle::AliEveGeomGentle()
{
    
}
AliEveGeomGentle::~AliEveGeomGentle()
{
    
}

TEveGeoShape* AliEveGeomGentle::GetGeomGentle(bool register_as_global)
{
    TFile f("$ALICE_ROOT/EVE/alice-data/gentle_geo.root");
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
    
    // finish
    if (register_as_global){
        gEve->AddGlobalElement(gsre);
    }
    return gsre;
}
TEveGeoShape* AliEveGeomGentle::GetGeomGentleRphi()
{
    // The resulting geometry is NOT added into the global scene!
    
    TFile f("$ALICE_ROOT/EVE/alice-data/gentle_rphi_geo.root");
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
TEveGeoShape* AliEveGeomGentle::GetGeomGentleRhoz()
{
    // The resulting geometry is NOT added into the global scene!
    
    TFile f("$ALICE_ROOT/EVE/alice-data/gentle_rhoz_geo.root");
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

TEveGeoShape* AliEveGeomGentle::GetGeomGentleTRD(Color_t color)
{
    TFile f("$ALICE_ROOT/EVE/alice-data/gentle_geo_trd.root");
    TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle TRD");
    TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
    gEve->AddGlobalElement(gsre);
    f.Close();
    
    const Int_t smInstalled[]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};
    const Int_t nInstalled = static_cast<Int_t>(sizeof(smInstalled)/sizeof(Int_t));
    Int_t sm = 0;
    // Fix visibility, color and transparency
    gsre->SetRnrSelf(kFALSE);
    for (TEveElement::List_i i = gsre->BeginChildren(); i != gsre->EndChildren(); ++i)
    {
        TEveGeoShape* lvl1 = (TEveGeoShape*) *i;
        lvl1->SetRnrSelf(kFALSE);
        for (TEveElement::List_i j = lvl1->BeginChildren(); j != lvl1->EndChildren(); ++j)
        {
            TEveGeoShape* lvl2 = (TEveGeoShape*) *j;
            lvl2->SetRnrSelf(kFALSE);
            for(Int_t ism(nInstalled); ism--;){
                if ( sm == smInstalled[ism] ){
                    lvl2->SetRnrSelf(kTRUE);
                    break;
                }
            }
            lvl2->SetMainColor(color);
            lvl2->SetMainTransparency(80);
            
            ++sm;
        }
    }
    
    return gsre;
}

TEveGeoShape* AliEveGeomGentle::GetGeomGentleMUON(bool updateScene, Color_t color)
{
    TFile f("$ALICE_ROOT/EVE/alice-data/gentle_geo_muon.root");
    TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle MUON");
    TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
    gEve->AddGlobalElement(gsre);
    f.Close();
    
    DrawDeep(gsre,color);
    
    if ( updateScene ) {
        TGLViewer* v = gEve->GetDefaultGLViewer();
        v->UpdateScene();
    }
    
    return gsre;
}


void AliEveGeomGentle::DrawDeep(TEveGeoShape *gsre, Color_t color)
{
    if(gsre->HasChildren())
    {
        gsre->SetRnrSelf(kFALSE);
        for (TEveElement::List_i i = gsre->BeginChildren(); i != gsre->EndChildren(); ++i)
        {
            TEveGeoShape* lvl = (TEveGeoShape*) *i;
            DrawDeep(lvl,color);
        }
    }
    else
    {
        gsre->SetRnrSelf(kTRUE);
        gsre->SetMainColor(color);
        gsre->SetMainTransparency(80);
    }
}


