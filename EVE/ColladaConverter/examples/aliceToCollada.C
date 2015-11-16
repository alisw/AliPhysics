//#include "AliColladaBuffer.h"
//
//#include "TGeoShape.h"
//#include "TEveGeoShapeExtract.h"
//#include "TEveGeoShape.h"
//#include "TEveManager.h"
//#include "TFile.h"
//#include "TString.h"
//
//#include <iostream>
//#include <string>
//
//using namespace std;

AliColladaBuffer *colladaBuffer;
TString removeSpaces(TString input);

void ReadShapes(TEveGeoShapeExtract *shapeExtract,TString parent)
{
    if(shapeExtract->HasElements())
    {
        TIter next(shapeExtract->GetElements());
        TEveGeoShapeExtract* nextShape;
        
        TString shapeName = shapeExtract->GetName();
        TString currentName = Form("%s_%s",parent.Data(),shapeName.Data());
        
        if(parent.Contains("Gentle Geometry")) parent = removeSpaces(parent);
        if(currentName.Contains("Gentle Geometry")) currentName = removeSpaces(currentName);
        
        colladaBuffer->CreateNode(currentName,parent);
        
        if(shapeName.EqualTo("TPC_Drift_1")){ // this one is both node and a shape
            colladaBuffer->AddShape(TString::Format("%s_geom",currentName.Data()),shapeExtract,parent,"lambert-TPC");
        }
        
        while ((nextShape = (TEveGeoShapeExtract*) next()))
        {
            ReadShapes(nextShape,currentName); // current name is a parent for next level shapes
        }
    }
    else
    {
        TEveGeoShape *shape = TEveGeoShape::ImportShapeExtract(shapeExtract);
        string materialName;
        
        TString shapeName = shapeExtract->GetName();
        TString name = Form("%s_%s",parent.Data(),shapeName.Data());
        
        if(name.Contains("Gentle Geometry")) name = removeSpaces(name);
        
        // select correct material for this detector
             if(name.Contains("BSEG")){materialName = "lambert-MUON-TRD";}
        else if(name.Contains("HMPID")){materialName = "lambert-HMPID";}
        else if(name.Contains("TPC")){materialName = "lambert-TPC";}
        else if(name.Contains("EMCAL")){materialName = "lambert-EMCAL";}
        else if(name.Contains("PHOS")){materialName = "lambert-EMCAL";}
        else if(name.Contains("IT34")){materialName = "lambert-TPC";}
        else if(name.Contains("IT56")){materialName = "lambert-EMCAL";}
        else if(name.Contains("MUON")){materialName = "lambert-MUON-TRD";}
        else if(name.Contains("B076")){materialName = "lambert-TOF";}
        else {materialName = "default_lambert";}
        
        // add shape to Collada buffer
        colladaBuffer->AddShape(name,shapeExtract,parent,materialName.c_str());
    }
}

void aliceToCollada(const char* outFile="aliceGeom.dae")
{
    // prepare collada buffer
    colladaBuffer = new AliColladaBuffer();
    
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

    emission->SetRGB(1.0,0.0,0.0);
    colladaBuffer->AddNewLambertMaterial("lambert-TOF",emission,ambient,diffuse,0.5);
    
    delete emission;delete ambient;delete diffuse;
    
    // Get TGeo shapes from simplified geometry files
    TEveManager::Create();
    TFile f("../../alice-data/gentle_geo_muon.root");
    TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle MUON");
    f.Close();
    colladaBuffer->CreateNode("MUON");
    ReadShapes(gse,"MUON");
    
    TFile f2("../../alice-data/gentle_geo.root");
    gse = (TEveGeoShapeExtract*) f2.Get("Gentle");
    f2.Close();
    colladaBuffer->CreateNode("MAIN");
    ReadShapes(gse,"MAIN");
    
    TFile f3("../../alice-data/gentle_geo_emcal.root");
    gse = (TEveGeoShapeExtract*) f3.Get("Gentle EMCAL");
    f3.Close();
    colladaBuffer->CreateNode("EMCAL");
    ReadShapes(gse,"EMCAL");
    //
    
    // save collada buffer to file
    colladaBuffer->SaveAs(outFile);
    
    // clean up
    delete colladaBuffer;
    return;
}

TString removeSpaces(TString input)
{
    TString output("");
    for(int i=0;i<input.Length();i++)
    {
        if(input[i]!=char(32))
        {
            output += input[i];
        }
    }
    return output;
}
