//
//  AliEveGeomGentle.cxx
//  xAliRoot
//
//  Created by Jeremi Niedziela on 11/05/15.
//
//

#include "AliEveGeomGentle.h"
#include "AliEveInit.h"

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
    AliEveInit::GetConfig(&fSettings);
}
AliEveGeomGentle::~AliEveGeomGentle()
{
    
}

TEveGeoShape* AliEveGeomGentle::GetGeomGentle(bool register_as_global)
{
    TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_geo.root");
    TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
    TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
    f.Close();
    
    // set PHOS colors
    TEveElement* elPHOS = gsre->FindChild("PHOS");
    elPHOS->SetRnrState(kTRUE);
    elPHOS->FindChild("PHOS_5")->SetRnrState(kFALSE);
    
    for (int i=0; i<5; i++) {
        elPHOS->FindChild(Form("PHOS_%i",i+1))->SetMainColor(fSettings.GetValue("PHOS.color",593));
        elPHOS->FindChild(Form("PHOS_%i",i+1))->SetMainTransparency(70);
    }
    
    // set TPC color
    TEveElement *elTPC = gsre->FindChild("TPC");
    TEveElement *elTPC_M_1 = elTPC->FindChild("TPC_M_1");
    TEveElement *elTPC_Drift_1 = elTPC_M_1->FindChild("TPC_Drift_1");
    elTPC_Drift_1->SetMainColor(fSettings.GetValue("TPC.color",3));
    
    // set HMPID color
    TEveElement* elHMPID = gsre->FindChild("HMPID");
    
    for (int i=0; i<7; i++) {
        elHMPID->FindChild(Form("HMPID_%i",i))->SetMainColor(fSettings.GetValue("HMPID.color",5));
        elHMPID->FindChild(Form("HMPID_%i",i))->SetMainTransparency(80);
    }
    
    // set ITS color
    TEveElement* elITS = gsre->FindChild("ITS");
    TEveElement* elITS_Dets = elITS->FindChild("ITS_Dets");
    TEveElement* elIT12_1 = elITS_Dets->FindChild("IT12_1");
    TEveElement* elIT34_1 = elITS_Dets->FindChild("IT34_1");
    TEveElement* elIT56_1 = elITS_Dets->FindChild("IT56_1");
    
    elIT12_1->SetMainColor(fSettings.GetValue("ITS.SPD.color",924));
    elIT34_1->SetMainColor(fSettings.GetValue("ITS.SDD.color",925));
    elIT56_1->SetMainColor(fSettings.GetValue("ITS.SSD.color",926));
    
    // set TOF color
    TEveElement* elTRD_TOF = gsre->FindChild("TRD+TOF");
    TEveElement* elB076_1 = elTRD_TOF->FindChild("B076_1");
    TEveElement* elBREF_1 = elTRD_TOF->FindChild("BREF_1");
    elB076_1->SetMainColor(fSettings.GetValue("TOF.color",930));
    elBREF_1->SetMainColor(fSettings.GetValue("TRD.ref.color",929));
    
    	
    // finish
    if (register_as_global){
        gEve->AddGlobalElement(gsre);
    }
    return gsre;
}
TEveGeoShape* AliEveGeomGentle::GetGeomGentleRphi()
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
    
    for (int i=0;i<5;i++) {
        elPHOS->FindChild(Form("PHOS_%i",i+1))->SetMainColor(fSettings.GetValue("PHOS.color",593));
        elPHOS->FindChild(Form("PHOS_%i",i+1))->SetMainTransparency(70);
    }

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
        elTPC_SWS1[i]->SetMainColor(fSettings.GetValue("TPC.color",3));
        elTPC_SWS1[i]->SetMainTransparency(70);
    }
    
    // set HMPID's color
    TEveElement* elHMPID = gsre->FindChild("HMPID");
    
    for (int i=0; i<7; i++) {
        elHMPID->FindChild(Form("HMPID_%i",i))->SetMainColor(fSettings.GetValue("HMPID.color",5));
        elHMPID->FindChild(Form("HMPID_%i",i))->SetMainTransparency(80);
    }

    // set ITS color
    TEveElement* elITS = gsre->FindChild("ITS");
    TEveElement* elITS_Dets = elITS->FindChild("ITS_Dets");
    TEveElement* elIT12_1 = elITS_Dets->FindChild("IT12_1");
    TEveElement* elIT34_1 = elITS_Dets->FindChild("IT34_1");
    TEveElement* elIT56_1 = elITS_Dets->FindChild("IT56_1");
    
    elIT12_1->SetMainColor(fSettings.GetValue("ITS.SPD.color",924));
    elIT34_1->SetMainColor(fSettings.GetValue("ITS.SDD.color",925));
    elIT56_1->SetMainColor(fSettings.GetValue("ITS.SSD.color",926));
    
    // set TOF color
    TEveElement* elTRD_TOF = gsre->FindChild("TRD+TOF");
    TEveElement* elB076_1 = elTRD_TOF->FindChild("B076_1");
    TEveElement* elBREF_1 = elTRD_TOF->FindChild("BREF_1");
    elB076_1->SetMainColor(fSettings.GetValue("TOF.color",930));
    elBREF_1->SetMainColor(fSettings.GetValue("TRD.ref.color",929));

    
    return gsre;
}
TEveGeoShape* AliEveGeomGentle::GetGeomGentleRhoz()
{
    // The resulting geometry is NOT added into the global scene!
    
    TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_rhoz_geo.root");
    TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
    TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
    f.Close();
    
    TEveElement* elPHOS = gsre->FindChild("PHOS");
    elPHOS->FindChild("PHOS_5")->SetRnrState(kFALSE);

    for (int i=0; i<5; i++) {
        elPHOS->FindChild(Form("PHOS_%i",i+1))->SetMainColor(fSettings.GetValue("PHOS.color",593));
        elPHOS->FindChild(Form("PHOS_%i",i+1))->SetMainTransparency(70);
    }
    
    // set TPC's color
    TEveElement *elTPC = gsre->FindChild("TPC");
    TEveElement *elTPC_M_1 = elTPC->FindChild("TPC_M_1");
    TEveElement *elTPC_Drift_1 = elTPC_M_1->FindChild("TPC_Drift_1");
    elTPC_Drift_1->SetMainColor(fSettings.GetValue("TPC.color",3));
    
    // set HMPID's color
    TEveElement* elHMPID = gsre->FindChild("HMPID");
    
    for(int i=0;i<7;i++){
        elHMPID->FindChild(Form("HMPID_%i",i))->SetMainColor(fSettings.GetValue("HMPID.color",5));
        elHMPID->FindChild(Form("HMPID_%i",i))->SetMainTransparency(80);
    }
    
    // set ITS color
    TEveElement* elITS = gsre->FindChild("ITS");
    TEveElement* elITS_Dets = elITS->FindChild("ITS_Dets");
    TEveElement* elIT12_1 = elITS_Dets->FindChild("IT12_1");
    TEveElement* elIT34_1 = elITS_Dets->FindChild("IT34_1");
    TEveElement* elIT56_1 = elITS_Dets->FindChild("IT56_1");
    
    elIT12_1->SetMainColor(fSettings.GetValue("ITS.SPD.color",924));
    elIT34_1->SetMainColor(fSettings.GetValue("ITS.SDD.color",925));
    elIT56_1->SetMainColor(fSettings.GetValue("ITS.SSD.color",926));
    
    // set TOF color
    TEveElement* elTRD_TOF = gsre->FindChild("TRD+TOF");
    TEveElement* elB076_1 = elTRD_TOF->FindChild("B076_1");
    TEveElement* elBREF_1 = elTRD_TOF->FindChild("BREF_1");
    elB076_1->SetMainColor(fSettings.GetValue("TOF.color",930));
    elBREF_1->SetMainColor(fSettings.GetValue("TRD.ref.color",929));
    
    return gsre;
}

TEveGeoShape* AliEveGeomGentle::GetGeomGentleTRD()
{
    TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_geo_trd.root");
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
            lvl2->SetMainColor(fSettings.GetValue("TRD.color", 920));
            lvl2->SetMainTransparency(80);
            
            ++sm;
        }
    }
    
    return gsre;
}

TEveGeoShape* AliEveGeomGentle::GetGeomGentleEMCAL()
{
    TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_geo_emcal.root");
    TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle EMCAL");
    TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
    gEve->AddGlobalElement(gsre);
    f.Close();
    
    // Fix visibility, color and transparency
    gsre->SetRnrSelf(kTRUE);
    gsre->SetMainColor(fSettings.GetValue("EMCAL.color", 953));
    gsre->SetMainTransparency(70);
    
    return gsre;
}

TEveGeoShape* AliEveGeomGentle::GetGeomGentleZDC()
{
    TFile f1("$ALICE_ROOT/EVE/resources/geometry/gentle_geom_zdca.root");
    TFile f2("$ALICE_ROOT/EVE/resources/geometry/gentle_geom_zdcc.root");
    TEveGeoShapeExtract* gse1 = (TEveGeoShapeExtract*) f1.Get("Gentle ZDCA");
    TEveGeoShapeExtract* gse2 = (TEveGeoShapeExtract*) f2.Get("Gentle ZDCC");
    TEveGeoShape* gsre1 = TEveGeoShape::ImportShapeExtract(gse1);
    TEveGeoShape* gsre2 = TEveGeoShape::ImportShapeExtract(gse2);
    gEve->AddGlobalElement(gsre1);
    gEve->AddGlobalElement(gsre2);
    f1.Close();
    f2.Close();
    
    // Fix visibility, color and transparency
    gsre1->SetRnrSelf(kTRUE);
    gsre1->SetMainColor(fSettings.GetValue("EMCAL.color", 953));
    gsre1->SetMainTransparency(70);
    
    gsre2->SetRnrSelf(kTRUE);
    gsre2->SetMainColor(fSettings.GetValue("EMCAL.color", 953));
    gsre2->SetMainTransparency(70);
    
    TEveGeoShape *zdc = new TEveGeoShape("ZDC");
    zdc->AddElement(gsre1);
    zdc->AddElement(gsre2);
    
    return zdc;
}

TEveGeoShape* AliEveGeomGentle::GetGeomGentleMUON(bool updateScene)
{
    TFile f("$ALICE_ROOT/EVE/resources/geometry/gentle_geo_muon.root");
    TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle MUON");
    TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
    gEve->AddGlobalElement(gsre);
    f.Close();
    
    DrawDeep(gsre);
    
    if ( updateScene ) {
        TGLViewer* v = gEve->GetDefaultGLViewer();
        v->UpdateScene();
    }
    
    return gsre;
}


void AliEveGeomGentle::DrawDeep(TEveGeoShape *gsre)
{
    if(gsre->HasChildren())
    {
        gsre->SetRnrSelf(kFALSE);
        for (TEveElement::List_i i = gsre->BeginChildren(); i != gsre->EndChildren(); ++i)
        {
            TEveGeoShape* lvl = (TEveGeoShape*) *i;
            DrawDeep(lvl);
        }
    }
    else
    {
        gsre->SetRnrSelf(kTRUE);
        gsre->SetMainColor(fSettings.GetValue("MUON.color",920));
        gsre->SetMainTransparency(80);
    }
}


