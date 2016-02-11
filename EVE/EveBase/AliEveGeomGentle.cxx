//
//  AliEveGeomGentle.cxx
//
//  Created by Jeremi Niedziela on 11/05/15.
//
//

#include "AliEveGeomGentle.h"
#include "AliEveInit.h"

#include <TFile.h>

TEveGeoShape* AliEveGeomGentle::GetSimpleGeom(char* detector)
{
    TFile f(Form("$ALICE_ROOT/EVE/resources/geometry/simple_geom_%s.root",detector));
    TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get(detector);
    TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
    f.Close();
    
    // tricks for different R-Phi geom of TPC:
    if(strcmp(detector,"RPH")!=0) // don't add RPhi TPC geom to the 3D view
    {
        gEve->AddGlobalElement(gsre);
    }
    else // but use all other parameters of regular TPC geom
    {
        detector = "TPC";
    }

    TEnv settings;
    AliEveInit::GetConfig(&settings);
    
    DrawDeep(gsre,
             settings.GetValue(Form("%s.color",detector),-1),
             settings.GetValue(Form("%s.trans",detector),-1),
             settings.GetValue(Form("%s.line.color",detector),-1));
    
    gEve->GetDefaultGLViewer()->UpdateScene();
    
    return gsre;
}


void AliEveGeomGentle::DrawDeep(TEveGeoShape *gsre,Color_t color, Char_t transparency, Color_t lineColor)
{
    if(gsre->HasChildren())
    {
        gsre->SetRnrSelf(kFALSE);
        
        if(strcmp(gsre->GetElementName(),"TPC_Drift_1")==0) // hack for TPC drift chamber
        {
            gsre->SetRnrSelf(kTRUE);
            if(color>=0) gsre->SetMainColor(color);
            if(lineColor>=0){
                gsre->SetLineColor(lineColor);
                gsre->SetLineWidth(0.1);
                gsre->SetDrawFrame(true);
            }
            else{
                gsre->SetDrawFrame(false);
            }
            if(transparency>=0) gsre->SetMainTransparency(transparency);
        }
        
        for (TEveElement::List_i i = gsre->BeginChildren(); i != gsre->EndChildren(); ++i)
        {
            TEveGeoShape* lvl = (TEveGeoShape*) *i;
            DrawDeep(lvl,color,transparency,lineColor);
        }
    }
    else
    {
        gsre->SetRnrSelf(kTRUE);
        if(color>=0) gsre->SetMainColor(color);
        if(lineColor>=0){
            gsre->SetLineColor(lineColor);
            gsre->SetLineWidth(0.1);
            gsre->SetDrawFrame(true);
        }
        else{
            gsre->SetDrawFrame(false);
        }
        if(transparency>=0) gsre->SetMainTransparency(transparency);
        
        if(strcmp(gsre->GetElementName(),"PHOS_5")==0) // hack for PHOS module which is not installed
        {
            gsre->SetRnrSelf(false);
        }
    }
}


