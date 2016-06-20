// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#include <AliEveConfigManager.h>
#include <AliEveInit.h>
//#include <AliEveMyoListener.h>

#include <AliLog.h>

#include <TEveUtil.h>
#include <TEveManager.h>
#include <TEveSelection.h>
#include <TEveBrowser.h>
#include <TEveViewer.h>

#include <TInterpreter.h>
#include <TRint.h>
#include <TROOT.h>
#include <TPRegexp.h>
#include <TSystem.h>
#include <TError.h>
//#include <TThread.h>

#include <RVersion.h>
#include <Getline.h>
#include <iostream>

using namespace std;

//void *MyoListenerThreadHandle(void*)
//{
//    cout<<"AliEVE -- Starting Myo Listener thread"<<endl;
//    AliEveMyoListener *listener = new AliEveMyoListener();
//    if(listener){delete listener;}
//    return 0;
//}

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
    macPath += Form(":%s/macros/data", evedir.Data());
    macPath += Form(":%s/macros/common", evedir.Data());
    gInterpreter->AddIncludePath(Form("%s/EVE", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(Form("%s/PWG0", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(Form("%s/include", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(gSystem->Getenv("ALICE_ROOT"));

    gROOT->SetMacroPath(macPath);

    bool storageManager=false;
    const char* cdbPath="local:///local/cdb";
    bool classesMode=false;
    AliEveEventManager::EDataSource dataSource=AliEveEventManager::kSourceOffline;
    
    for (int i=0; i<argc; i++)
    {
        if(strcmp(argv[i],"online")==0){    dataSource = AliEveEventManager::kSourceOnline;classesMode=true; }
        if(strcmp(argv[i],"hlt")==0){       dataSource = AliEveEventManager::kSourceHLT;classesMode=true; }
        if(strcmp(argv[i],"local") ==0){    cdbPath = "local:///local/cdb";classesMode=true; }
        if(strcmp(argv[i],"mc")==0){        cdbPath = "mcideal://"; }
        if(strcmp(argv[i],"raw")==0){       cdbPath = "raw://";}
        if(strcmp(argv[i],"mcfull")==0){    cdbPath = "mcfull://"; }
        if(strcmp(argv[i],"mcresidual")==0){cdbPath = "mcresidual://"; }

        if(strcmp(argv[i],"sm")==0){storageManager=true;}
    }
    
    // make sure logger is instantiated
    AliLog::GetRootLogger();
    TRint *app = new TRint("App", &argc, argv);

    TEveManager::Create();
    gEve->GetDefaultViewer()->SetElementName("Default View");
    gEve->GetSelection()->SetPickToSelect(TEveSelection::kPS_PableCompound);
    gEve->GetHighlight()->SetPickToSelect(TEveSelection::kPS_PableCompound);
    gEve->RegisterGeometryAlias("Default", Form("%s/resources/geometry/default_geo.root", evedir.Data()));

    try {
        AliEveConfigManager::Instance(storageManager);
    }
    catch (TEveException exc) {
        cout<<"\n\nException while initializing config manager"<<endl;
        AliErrorGeneral("alieve_main",exc.Data());
    }

    if(classesMode){
        AliEveInit *init = new AliEveInit(".",dataSource);
    }
//
//    TThread *listenerThread = new TThread("listenerThread", MyoListenerThreadHandle,NULL);
//    listenerThread->Run();
    
    app->Connect("TEveBrowser", "CloseWindow()", "TRint", app, "Terminate()");
    app->Run(kTRUE);

    if (gEve && gEve->GetBrowser())	gEve->GetBrowser()->UnmapWindow();
    TEveManager::Terminate();
    if(gEve) {delete gEve; gEve = 0;}
    
    app->Terminate(0);
//    listenerThread->Join();
    return 0;
}
