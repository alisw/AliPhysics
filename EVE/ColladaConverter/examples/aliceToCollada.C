//#include "AliColladaBuffer.h"
//
//#include "TGeoShape.h"
//#include "TEveGeoShapeExtract.h"
//#include "TEveGeoShape.h"
//#include "TEveManager.h"
//#include "TFile.h"
//
//#include <iostream>
//#include <string>
//
//using namespace std;

AliColladaBuffer *colladaBuffer;
int iter;

void ReadShapes(TEveGeoShapeExtract *shapeExtract,TString name)
{
    if(shapeExtract->HasElements())
    {
        TIter next(shapeExtract->GetElements());
        TEveGeoShapeExtract* nextShape;
        while ((nextShape = (TEveGeoShapeExtract*) next()))
        {
            if(name.EqualTo("MAIN/TPC/TPC_M_1/TPC_Drift_1"))
            {
                cout<<name<<endl;
                colladaBuffer->AddShape(name,shapeExtract,"lambert-TPC");
            }
            
            ReadShapes(nextShape,Form("%s/%s",name.Data(),nextShape->GetName()));
        }
    }
    else
    {
        TEveGeoShape *shape = TEveGeoShape::ImportShapeExtract(shapeExtract);
        string materialName;
        
        // select correct material for this detector
        if(strcmp(shape->GetShape()->ClassName(),"TGeoTrd1")==0){materialName = "lambert-MOUN-TRD";}
        else if(name.Contains("HMPID")){materialName = "lambert-HMPID";}
        else if(name.Contains("TPC")){materialName = "lambert-TPC";}
        else if(name.Contains("EMCAL")){materialName = "lambert-EMCAL";}
        else if(name.Contains("MUON")){materialName = "lambert-MUON";}
        else{materialName = "default_lambert";}
        
        // add shape to Collada buffer
        cout<<name<<endl;
        colladaBuffer->AddShape(name,shapeExtract,materialName.c_str());
    }
}

void aliceToCollada(const char* outFile="aliceGeom.dae")
{
    // prepare collada buffer
    colladaBuffer = new AliColladaBuffer();
    iter=0;
    
    // Create new materials
    TColor *emission = new TColor();
    TColor *ambient = new TColor();
    TColor *diffuse = new TColor();
    
    emission->SetAlpha(1.0);
    ambient->SetRGB(0.0,0.0,0.0);
    ambient->SetAlpha(1.0);
    diffuse->SetRGB(0.0,0.0,0.0);
    diffuse->SetAlpha(1.0);
    
    emission->SetRGB(0.,1.0,0.0);
    colladaBuffer->AddNewLambertMaterial("lambert-TPC",emission,ambient,diffuse,0.5);

    emission->SetRGB(0.5,0.5,0.5);
    colladaBuffer->AddNewLambertMaterial("lambert-MUON-TRD",emission,ambient,diffuse,0.5);

    emission->SetRGB(1.0,1.0,0.0);
    colladaBuffer->AddNewLambertMaterial("lambert-HMPID",emission,ambient,diffuse,0.5);

    emission->SetRGB(0.0,0.0,1.0);
    colladaBuffer->AddNewLambertMaterial("lambert-EMCAL",emission,ambient,diffuse,0.5);
    
    
    // Get TGeo shapes from simplified geometry files
    TEveManager::Create();
    TFile f("../../alice-data/gentle_geo_muon.root");
    TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle MUON");
    f.Close();
    ReadShapes(gse,"MUON");
    
    TFile f2("../../alice-data/gentle_geo.root");
    gse = (TEveGeoShapeExtract*) f2.Get("Gentle");
    f2.Close();
    ReadShapes(gse,"MAIN");
    
    cout<<"EMCAL"<<endl;
    
    TFile f3("../../alice-data/gentle_geo_emcal.root");
    gse = (TEveGeoShapeExtract*) f3.Get("Gentle EMCAL");
    f3.Close();
    ReadShapes(gse,"EMCAL");
    //
    
    // save collada buffer to file
    colladaBuffer->SaveAs(outFile);
    
    // clean up
    delete colladaBuffer;
    return;
}