/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

class TEveProjectionManager;
class TEveGeoShape;
class TEveUtil;
class AliTriggerAnalysis;
class AliSysInfo;

TH2D* V0StateHistogram;
Bool_t gCenterProjectionsAtPrimaryVertex = kFALSE;

Int_t      g_pic_id  = 0;
Int_t      g_pic_max = 100;
TTimeStamp g_pic_prev(0, 0);

void alieve_online_new()
{
    // set OCDB path:
    //AliEveEventManager::SetCdbUri("local://$ALICE_ROOT/OCDB"); // default OCDB from aliroot
    //AliEveEventManager::SetCdbUri("local:///local/OCDB/2013"); // OCDB snapshot for particular run
    AliEveEventManager::SetCdbUri("local:///local/cdb");         // current OCDB snapshot
    //AliEveEventManager::SetCdbUri("raw://");                   // reading OCDB from alien
    
    
    Info("alieve_init", "Adding standard macros.");
    TString  hack = gSystem->pwd(); // Problem with TGFileBrowser cding
    alieve_init_import_macros();
    gSystem->cd(Form("%s/../src/",gSystem->Getenv("ALICE_ROOT")));
    gROOT->ProcessLine(".L saveViews.C++");
    gROOT->ProcessLine(".L geom_gentle_muon.C++");
    gROOT->ProcessLine(".L geom_gentle_trd.C++");
    TEveUtil::LoadMacro("saveViews.C");
    gSystem->cd(hack);
    cout<<"Standard macros added"<<endl;
    
    new AliEveEventManager("online", -1);
    gEve->AddEvent(AliEveEventManager::GetMaster());
    
    TEveUtil::AssertMacro("VizDB_scan.C");
    gSystem->ProcessEvents();
    
    TEveBrowser *browser = gEve->GetBrowser();
    browser->ShowCloseTab(kFALSE);
    

    cout<<"Creating multiview"<<endl;
    AliEveMultiView *multiView = new AliEveMultiView(kTRUE);
    TEveUtil::LoadMacro("geom_gentle.C");
    multiView->InitGeomGentle(geom_gentle(),
                              geom_gentle_rphi(),
                              geom_gentle_rhoz(),
                              geom_gentle_rhoz());
    
    TEveUtil::LoadMacro("geom_gentle_trd.C");
    multiView->InitGeomGentleTrd(geom_gentle_trd());

    TEveUtil::LoadMacro("geom_gentle_muon.C");
    multiView->InitGeomGentleMuon(geom_gentle_muon(), kFALSE, kFALSE, kTRUE);
 
    //============================================================================
    // Standard macros to execute -- not all are enabled by default.
    //============================================================================
    
    printf("============ Setting macro executor ============\n");
    
    AliEveMacroExecutor *exec = AliEveEventManager::GetMaster()->GetExecutor();
    
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX",         "primary_vertex.C", "primary_vertex",             "",                kTRUE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Ellipse", "primary_vertex.C", "primary_vertex_ellipse",     "",                kTRUE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Box",     "primary_vertex.C", "primary_vertex_box",         "kFALSE, 3, 3, 3", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX",         "primary_vertex.C", "primary_vertex_spd",         "",                kTRUE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Ellipse", "primary_vertex.C", "primary_vertex_ellipse_spd", "",                kTRUE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Box",     "primary_vertex.C", "primary_vertex_box_spd",     "kFALSE, 3, 3, 3", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX",         "primary_vertex.C", "primary_vertex_tpc",         "",                kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Ellipse", "primary_vertex.C", "primary_vertex_ellipse_tpc", "",                kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC PVTX Box",     "primary_vertex.C", "primary_vertex_box_tpc",     "kFALSE, 3, 3, 3", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track",      "esd_tracks.C",        "esd_tracks",             "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track",      "esd_tracks.C",        "esd_tracks_MI",          "", kFALSE));
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track MUON", "esd_muon_tracks.C", "esd_muon_tracks",        "kTRUE,kFALSE", kTRUE));
    
    // these macros were leaking:
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track",      "esd_tracks.C",        "esd_tracks_by_category", "", kTRUE));// just a little
    exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC FMD",        "fmd_esd.C",           "fmd_esd",                "", kTRUE));//huge leak
    //
    
    // ???
    // exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC TRD", "trd_detectors.C", "trd_detectors",         "", kFALSE));
    // trd_tracks disabled due to memory leaks
    
    //----------------------------------------------------------------------------
    
    slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
    slot->StartEmbedding();
    AliEveMacroExecutorWindow* exewin = new AliEveMacroExecutorWindow(exec);
    slot->StopEmbedding("DataSelection");
    exewin->PopulateMacros();
    
    //============================================================================
    // Final GUI setup
    //============================================================================
    
    browser->GetTabRight()->SetTab(1);
    browser->StartEmbedding(TRootBrowser::kBottom);
    new AliEveEventManagerWindow(AliEveEventManager::GetMaster());
    browser->StopEmbedding("EventCtrl");
 
    browser->MoveResize(0, 0, gClient->GetDisplayWidth(),gClient->GetDisplayHeight() - 32);
    
    gEve->FullRedraw3D(kTRUE);
    gSystem->ProcessEvents();
    
    // move and rotate sub-views
    TGLViewer *glv1 = multiView->Get3DView()->GetGLViewer();
    TGLViewer *glv2 = multiView->GetRPhiView()->GetGLViewer();
    TGLViewer *glv3 = multiView->GetRhoZView()->GetGLViewer();
    
    glv1->CurrentCamera().RotateRad(-0.4, 0.6);
    glv2->CurrentCamera().Dolly(90, kFALSE, kFALSE);
    glv3->CurrentCamera().Dolly(1700, kFALSE, kFALSE);
    
    AliEveEventManager::GetMaster()->AddNewEventCommand("alieve_online_on_new_event();");
    gEve->FullRedraw3D();
    gSystem->ProcessEvents();
    gEve->Redraw3D(kTRUE);
    
    // set autoload by default
    AliEveEventManager::GetMaster()->SetAutoLoad(true);
}

void alieve_online_on_new_event()
{
    AliSysInfo::AddStamp("on_new_event_start");
    
    Double_t x[3] = { 0, 0, 0 };
    
    if (AliEveEventManager::HasESD())
    {
        AliESDEvent* esd = AliEveEventManager::AssertESD();
        esd->GetPrimaryVertex()->GetXYZ(x);
        
        TTimeStamp ts(esd->GetTimeStamp());
        TString win_title("Eve Main Window -- Timestamp: ");
        win_title += ts.AsString("s");
        win_title += "; Event # in ESD file: ";
        win_title += esd->GetEventNumberInFile();
        gEve->GetBrowser()->SetWindowName(win_title);
    }
    
    TEveElement* top = gEve->GetCurrentEvent();
    
    AliEveMultiView *mv = AliEveMultiView::Instance();
    
    //mv->DestroyEventRPhi();
    if (gCenterProjectionsAtPrimaryVertex)
        mv->SetCenterRPhi(x[0], x[1], x[2]);
    mv->ImportEventRPhi(top);
    
    //mv->DestroyEventRhoZ();
    if (gCenterProjectionsAtPrimaryVertex)
        mv->SetCenterRhoZ(x[0], x[1], x[2]);
    mv->ImportEventRhoZ(top);
    
    if (gCenterProjectionsAtPrimaryVertex)
        mv->SetCenterMuon(x[0], x[1], x[2]);
    mv->ImportEventMuon(top);
    
    
     // Register image to amore.
     // const TString pichost("aldaqacrs3");
     const TString pichost(gEnv->GetValue("AliEve.imageDumpHost", "localhost"));
     TTimeStamp now;
     Double_t delta = now.AsDouble() - g_pic_prev.AsDouble();
     
     printf("Pre image dump: host='%s', delta=%f.\n",gSystem->HostName(), delta);
     
     AliSysInfo::AddStamp("on_new_event_pic");
     // if (pichost == gSystem->HostName() && delta >= 30)
     {
       TString id;      id.Form("online-viz-%03d", g_pic_id);
       TString pic(id); pic += ".png";
     
       printf("In image dump: file='%s'.\n", pic.Data());
     
       gEve->GetBrowser()->RaiseWindow();
       gEve->FullRedraw3D();
       gSystem->ProcessEvents();
     
       Int_t status;
     
       // create screenshots from OpenGL views
       saveViews(pic.Data());
     
       // send screenshot to AMORE
       cout<<"Sending:"<<TString::Format("SendImageToAmore %s %s %d",id.Data(), pic.Data(),AliEveEventManager::AssertESD()->GetRunNumber())<<endl;

       status = gSystem->Exec(TString::Format("SendImageToAmore %s %s %d",
					      id.Data(), pic.Data(),
					      AliEveEventManager::AssertESD()->GetRunNumber()));
     
       printf("Post AMORE reg -- status=%d, run=%d.\n", status,
	      AliEveEventManager::AssertESD()->GetRunNumber());
     
       if (++g_pic_id >= g_pic_max)
	 g_pic_id = 0;
       g_pic_prev.Set();
     }
    AliSysInfo::AddStamp("on_new_event_end");
}

void alieve_init_import_macros()
{
    // Put macros in the list of browsables, add a macro browser to
    // top-level GUI.
    
    if (gSystem->Getenv("ALICE_ROOT") != 0)
    {
        gInterpreter->AddIncludePath(Form("%s/MUON", gSystem->Getenv("ALICE_ROOT")));
        gInterpreter->AddIncludePath(Form("%s/MUON/mapping", gSystem->Getenv("ALICE_ROOT")));
        TString macdir("$(ALICE_ROOT)/EVE/alice-macros");
        gSystem->ExpandPathName(macdir);
    }
    
    
    TFolder* f = gEve->GetMacroFolder();
    void* dirhandle = gSystem->OpenDirectory(macdir.Data());
    if (dirhandle != 0)
    {
        char* filename;
        TPMERegexp re("\\.C$");
        TObjArray names;
        while ((filename = gSystem->GetDirEntry(dirhandle)) != 0)
        {
            if (re.Match(filename))
                names.AddLast(new TObjString(filename));
        }
        names.Sort();
        
        for (Int_t ii=0; ii<names.GetEntries(); ++ii)
        {
            TObjString * si = (TObjString*) names.At(ii);
            f->Add(new TEveMacro(Form("%s/%s", macdir.Data(), (si->GetString()).Data())));
        }
    }
    gSystem->FreeDirectory(dirhandle);
    
    gROOT->GetListOfBrowsables()->Add(new TSystemDirectory(macdir.Data(), macdir.Data()));
    
    {
        TEveBrowser   *br = gEve->GetBrowser();
        TGFileBrowser *fb = 0;
        fb = br->GetFileBrowser();
        fb->GotoDir(macdir);
        {
            br->StartEmbedding(0);
            fb = br->MakeFileBrowser();
            fb->BrowseObj(f);
            fb->Show();
            br->StopEmbedding();
            br->SetTabTitle("Macros", 0);
            br->SetTab(0, 0);
        }
    }
}
