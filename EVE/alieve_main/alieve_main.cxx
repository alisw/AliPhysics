// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#include <AliEveConfigManager.h>

#include <TInterpreter.h>
#include <TRint.h>
#include <TROOT.h>
#include <TPRegexp.h>
#include <TSystem.h>
#include <TError.h>
#include <RVersion.h>

#include <AliLog.h>

#include <TEveUtil.h>
#include <TEveManager.h>
#include <TEveSelection.h>
#include <TEveBrowser.h>
#include <TEveViewer.h>

#include <Getline.h>

// DEV Version includes
#include <AliEveApplication.h>
#include <AliEveMainWindow.h>

#include "AliEveOnline.h"

#include <iostream>
using namespace std;

int main(int argc, char **argv)
{
    if (gSystem->Getenv("ALICE_ROOT") == 0){
        cout<<"ALICE_ROOT is not defined, aborting."<<endl;
        gSystem->Exit(1);
    }

    TString evedir(Form("%s/EVE", gSystem->Getenv("ALICE_ROOT")));

    if (gSystem->AccessPathName(evedir) == kTRUE){
        cout<<"Directory $ALICE_ROOT/EVE does not exist."<<endl;
        gSystem->Exit(1);
    }

    TString macPath(gROOT->GetMacroPath());
    macPath += Form(":%s/macros", evedir.Data());
    gInterpreter->AddIncludePath(evedir);
    macPath += Form(":%s/alice-macros", evedir.Data());
    gInterpreter->AddIncludePath(Form("%s/EVE", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(Form("%s/PWG0", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(Form("%s/include", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(gSystem->Getenv("ALICE_ROOT"));
    
    // TabCom fails on double-colon in macro-path.
    // I fixed this in ROOT sometime ago ... could be removed
    // when we go to 5.26.
//    TPMERegexp doubleColon(":{2,}", "og");
//    doubleColon.Substitute(macPath, ":");

    gROOT->SetMacroPath(macPath);

    bool onlineMode=false;
    
    if(argc==2)
    {
        if(strcmp(argv[1],"online")==0)
        {
            onlineMode = true;
            argc--;
            strcpy(argv[1],"");
        }
    }
    else if(argc==3)
    {
        if(strcmp(argv[2],"online")==0)
        {
            onlineMode = true;
            argc--;
            strcpy(argv[2],"");
        }
    }

    
    // make sure logger is instantiated
    AliLog::GetRootLogger();
    TRint *app = new TRint("App", &argc, argv);

    TEveManager::Create();
    gEve->GetDefaultViewer()->SetElementName("3D View");
    gEve->GetSelection()->SetPickToSelect(TEveSelection::kPS_PableCompound);
    gEve->GetHighlight()->SetPickToSelect(TEveSelection::kPS_PableCompound);
    gEve->RegisterGeometryAlias("Default", Form("%s/alice-data/default_geo.root", evedir.Data()));

    try {
        AliEveConfigManager::InitializeMaster();
    }
    catch (TEveException exc) {
        cout<<"\n\nException while initializing config manager"<<endl;
        AliErrorGeneral("alieve_main",exc.Data());
    }

    
    if(onlineMode)
    {
        AliEveOnline *online = new AliEveOnline();
    }
    
    app->Connect("TEveBrowser", "CloseWindow()", "TRint", app, "Terminate()");
    app->Run(kTRUE);

    if (gEve && gEve->GetBrowser())	gEve->GetBrowser()->UnmapWindow();
    TEveManager::Terminate();
    if(gEve) {delete gEve; gEve = 0;}
    
    app->Terminate(0);
    return 0;
}
