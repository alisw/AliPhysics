// Created by Jeremi Niedziela on 09/02/2016.

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TEveManager.h>
#include <TEveElement.h>
#include <TEveGeoNode.h>
#include <TSystem.h>
#include <TRegexp.h>

#include <iostream>
#include <string>
#include <vector>

using namespace std;

#endif

void AddNodes(TGeoNode *node, TEveGeoNode *parent, Int_t depth, Int_t depthmax,TObjArray *list)
{
    if (--depth <= 0) return;
    
    TObjArray *nlist = node->GetVolume()->GetNodes(); // all nodes in current level
    if (!nlist) return;
    
    TObjString *nname = (TObjString*)list->At(depthmax-depth); // name of required node in current level
    
    for (int i = 0; i < nlist->GetEntries(); i++)
    {   // loop over nodes in current level and find the one with matching name
        TGeoNode *node2 = (TGeoNode*) nlist->At(i);
        
        if (strcmp(node2->GetName(),nname->GetString().Data()) == 0)
        {
            TEveGeoNode *son = dynamic_cast<TEveGeoNode*>(parent->FindChild(nname->GetName()));
            if (!son)
            {
                son = new TEveGeoNode(node2);
                parent->AddElement(son);
            }
            AddNodes(node2,son, depth, depthmax, list);
        }
    }
}

void simple_geom_generate(char *detectorName="", int runNumber=0)
{
    if(strcmp(detectorName,"")==0){
        cout<<"Give name of the detector as a first argument!"<<endl;
        return;
    }

    
    // read all files with names matching "geom_list_XYZ.txt"
    vector<string> detectorsList;
    TSystemDirectory dir(".",".");
    TList *files = dir.GetListOfFiles();

    if (files)
    {
        TRegexp e("geom_list_[A-Z][A-Z][A-Z].txt");
        TRegexp e2("[A-Z][A-Z][A-Z]");
        
        TSystemFile *file;
        TString fname;
        TIter next(files);
        
        while ((file=(TSystemFile*)next()))
        {
            fname = file->GetName();
            if(fname.Contains(e))
            {
                TString detName = fname(e2);
                detName.Resize(3);
                detectorsList.push_back(detName.Data());
            }
        }
    }
    
    // load geometry library
    gSystem->Load("libGeom");
    
    // create visualisation manager
    TEveManager::Create();
    
    // load config file
    TEnv settings;
    AliEveInit::GetConfig(&settings);
    
    // set OCDB path from config and set run number for which we want to generate geometry
    AliCDBManager::Instance()->SetDefaultStorage(settings.GetValue("OCDB.default.path","local://$ALICE_ROOT/../src/OCDB"));
    AliCDBManager::Instance()->SetRun(runNumber);
    
    // load geometry from OCDB
    AliGeomManager::LoadGeometry();
    gGeoManager = AliGeomManager::GetGeometry();
    gGeoManager->DefaultColors();
    
    // find main node for our detector
    TGeoNode* tnode = gGeoManager->GetTopNode();
    tnode->SetVisibility(kFALSE);
    
    
    for(int i=0;i<detectorsList.size();i++)
    {
        TString path;
        TObjArray *list;
        Int_t depth;
        Char_t line[256];
        
        const char *currentDetector;
        if(strcmp(detectorName,"ALL")==0)
        {
            // if we update all detectors
            currentDetector = detectorsList[i].c_str();
        }
        else
        {
            // if we update specific detector, exit loop after one pass
            currentDetector = detectorName;
            i=detectorsList.size();
        }

        TEveGeoTopNode* eve_tnode = new TEveGeoTopNode(gGeoManager, tnode);
        eve_tnode->SetVisLevel(0);
        
        gEve->AddGlobalElement(eve_tnode);
        
        ifstream in(Form("geom_list_%s.txt",currentDetector), ios::in);
        cout<<"Adding shapes from file:"<<Form("geom_list_%s.txt",currentDetector)<<endl;
        
        int lineIter=0;

        while (true)
        {
            in >> line;
            if(in.eof())break;
            
            path = TString(line);

            if (!path.Contains("ALIC")) continue;
            
            list = path.Tokenize("/");
            depth = list->GetEntries();
            AddNodes(tnode,eve_tnode,depth,depth,list);
            lineIter++;
        }
        in.close();
        
        if(lineIter==0)
        {
            cout<<"File for "<<currentDetector<<" is empty. Skipping..."<<endl;
        }
        else
        {
            eve_tnode->SaveExtract(Form("../../resources/geometry/run3/simple_geom_%s.root",currentDetector),
                                   currentDetector, kTRUE);
        }
    }
}



