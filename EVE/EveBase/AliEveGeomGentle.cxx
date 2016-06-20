//
//  AliEveGeomGentle.cxx
//
//  Created by Jeremi Niedziela on 11/05/15.
//
//

#include "AliEveGeomGentle.h"
#include "AliEveInit.h"

#include <TFile.h>
#include <TSystem.h>
#include <iostream>

using namespace std;

TEveGeoShape* AliEveGeomGentle::GetSimpleGeom(char* detector)
{
    TEnv settings;
    AliEveInit::GetConfig(&settings);
    
    string geomPath = settings.GetValue("simple.geom.path","${ALICE_ROOT}/EVE/resources/geometry/run2/");
    string alirootBasePath = gSystem->Getenv("ALICE_ROOT");
    size_t alirootPos = geomPath.find("${ALICE_ROOT}");
    
    if(alirootPos != string::npos){
        geomPath.replace(alirootPos,alirootPos+13,alirootBasePath);
    }
    
    TFile *f = TFile::Open(Form("%s/simple_geom_%s.root",geomPath.c_str(),detector));
    if(!f){
        cout<<"AliEveGeomGentle::GetSimpleGeom -- no file with geometry found!"<<endl;
        return nullptr;
    }
    TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f->Get(detector);
    TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
    f->Close();
    
    // tricks for different R-Phi geom of TPC:
    if(strcmp(detector,"RPH")==0) // use all other parameters of regular TPC geom
    {
        detector = "TPC";
    }
    
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


