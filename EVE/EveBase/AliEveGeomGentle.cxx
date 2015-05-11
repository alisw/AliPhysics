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
    
    TEveElement* elPHOS = gsre->FindChild("PHOS");
    elPHOS->SetRnrState(kTRUE);
    //  elPHOS->FindChild("PHOS_4")->SetRnrState(kFALSE);
    elPHOS->FindChild("PHOS_5")->SetRnrState(kFALSE);
    
    if (register_as_global)
    {
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
    
    TEveElement* elPHOS = gsre->FindChild("PHOS");
    elPHOS->SetRnrState(kTRUE);
    //  elPHOS->FindChild("PHOS_4")->SetRnrState(kFALSE);
    elPHOS->FindChild("PHOS_5")->SetRnrState(kFALSE);
    
    return gsre;
}
TEveGeoShape* AliEveGeomGentle::GetGeomGentleRhoz()
{
    // The resulting geometry is NOT added into the global scene!
    
    TFile f("$ALICE_ROOT/EVE/alice-data/gentle_rhoz_geo.root");
    TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle");
    TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
    f.Close();
    
    return gsre;
}

TEveGeoShape* AliEveGeomGentle::GetGeomGentleTRD()
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
            lvl2->SetMainColor(3);
            lvl2->SetMainTransparency(80);
            
            ++sm;
        }
    }
    
    return gsre;
}

TEveGeoShape* AliEveGeomGentle::GetGeomGentleMUON(bool updateScene)
{
    TFile f("$ALICE_ROOT/EVE/alice-data/gentle_geo_muon.root");
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
    if (gsre->HasChildren()) {
        
        gsre->SetRnrSelf(kFALSE);
        for (TEveElement::List_i i = gsre->BeginChildren(); i != gsre->EndChildren(); ++i) {
            TEveGeoShape* lvl = (TEveGeoShape*) *i;
            DrawDeep(lvl);
        }
        
    } else {
        
        gsre->SetRnrSelf(kTRUE);
        gsre->SetMainColor(3);
        gsre->SetMainTransparency(80);
        
    }
    
}


