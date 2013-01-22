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

void alieve_online_init()
{
  
  if (gSystem->Getenv("ALICE_ROOT") != 0)
  {
    gInterpreter->AddIncludePath(Form("%s/MUON", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(Form("%s/MUON/mapping", gSystem->Getenv("ALICE_ROOT")));
  }
  
  TEveUtil::AssertMacro("VizDB_scan.C");

  TEveBrowser *browser = gEve->GetBrowser();
  browser->ShowCloseTab(kFALSE);

  // Gentle-geom loading changes gGeoManager.
  TEveGeoManagerHolder mgrRestore;

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

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus ITS",   "its_clusters.C++",   "its_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus TPC",   "tpc_clusters.C++",   "tpc_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus TRD",   "trd_clusters.C++",   "trd_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus TOF",   "tof_clusters.C++",   "tof_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus HMPID", "hmpid_clusters.C++", "hmpid_clusters"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "REC Clus MUON",  "muon_clusters.C++",  "muon_clusters"));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRunLoader, "DIG EMCAL",   "emcal_digits.C++",   "emcal_digits"));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW ITS",     "its_raw.C",     "its_raw"));
  //  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW TPC",     "tpc_raw.C",     "tpc_raw"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW TOF",     "tof_raw.C",     "tof_raw"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW VZERO",   "vzero_raw.C",   "vzero_raw", "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW ACORDE",  "acorde_raw.C",  "acorde_raw", "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW MUON",    "muon_raw.C++",  "muon_raw"));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kRawReader, "RAW FMD",     "fmd_raw.C",     "fmd_raw"));

  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track",      "esd_tracks.C++",        "esd_tracks",             "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track",      "esd_tracks.C++",        "esd_tracks_MI",          "", kFALSE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track",      "esd_tracks.C++",        "esd_tracks_by_category", "", kTRUE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Track MUON", "esd_muon_tracks.C++", "esd_muon_tracks",        "kTRUE,kFALSE", kTRUE));
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC FMD",        "fmd_esd.C",           "fmd_esd",                "", kTRUE));

  // ???
  // exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC TRD", "trd_detectors.C++", "trd_detectors",         "", kFALSE));
  // trd_tracks disabled due to memory leaks

  //----------------------------------------------------------------------------

  slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  slot->StartEmbedding();
  AliEveMacroExecutorWindow* exewin = new AliEveMacroExecutorWindow(exec);
  slot->StopEmbedding("DataSelection");
  exewin->PopulateMacros();

  //============================================================================
  // VZERO state histogram
  //============================================================================

  slot = TEveWindow::CreateWindowInTab(browser->GetTabRight());
  slot->StartEmbedding();  
  TCanvas* pad = new TCanvas();

  gStyle->SetCanvasColor(0);

  V0StateHistogram = new TH2D("V0 Histogram","V0 Trigger State", 4, 0, 4, 4, 0, 4);
  V0StateHistogram->Draw("colz");

  V0StateHistogram->GetXaxis()->SetBinLabel(1,"V0A Invalid");
  V0StateHistogram->GetXaxis()->SetBinLabel(2,"V0A Empty");
  V0StateHistogram->GetXaxis()->SetBinLabel(3,"V0A BB");
  V0StateHistogram->GetXaxis()->SetBinLabel(4,"V0A BG");

  V0StateHistogram->GetYaxis()->SetBinLabel(1,"V0C Invalid");
  V0StateHistogram->GetYaxis()->SetBinLabel(2,"V0C Empty");
  V0StateHistogram->GetYaxis()->SetBinLabel(3,"V0C BB");
  V0StateHistogram->GetYaxis()->SetBinLabel(4,"V0C BG");

  slot->StopEmbedding("V0 Trigger State");

  //============================================================================
  // Final GUI setup
  //============================================================================

  browser->GetTabRight()->SetTab(1);

  browser->StartEmbedding(TRootBrowser::kBottom);
  new AliEveEventManagerWindow(AliEveEventManager::GetMaster());
  browser->StopEmbedding("EventCtrl");

  browser->MoveResize(0, 0, gClient->GetDisplayWidth(),
		      gClient->GetDisplayHeight() - 32);

  gEve->GetViewers()->SwitchColorSet();

  TString autoRun(gSystem->Getenv("ONLINERECO_AUTORUN"));
  if (autoRun == "1" || autoRun.CompareTo("true", TString::kIgnoreCase) == 0)
  {
    AliEveEventManager::GetMaster()->SetAutoLoad(kTRUE);
  }

  {
    TGTab *tab = gEve->GetBrowser()->GetTab(2);

    TGHorizontalFrame *hf = (TGHorizontalFrame*) tab->GetParent();
    TGVerticalFrame   *vf = (TGVerticalFrame*)   hf ->GetParent();

    hf->Resize(hf->GetWidth(), hf->GetHeight() + 80);
    vf->Layout();
  }

  gEve->GetWindowManager()->HideAllEveDecorations();

  if(gEve->GetScenes()->FindChild("Geometry scene")->FindChild("Gentle MUON"))
  {
    gEve->GetScenes()->FindChild("Geometry scene")->FindChild("Gentle MUON")->SetRnrSelf(kFALSE);
    gEve->GetScenes()->FindChild("Geometry scene")->FindChild("Gentle MUON")->SetRnrChildren(kFALSE);
  }

  if(gEve->GetScenes()->FindChild("Muon Geometry")->FindChild("Gentle MUON [P]"))
  {
    gEve->GetScenes()->FindChild("Muon Geometry")->FindChild("Gentle MUON [P]")->SetRnrSelf(kTRUE);
    gEve->GetScenes()->FindChild("Muon Geometry")->FindChild("Gentle MUON [P]")->SetRnrChildren(kTRUE);
  }

  gEve->FullRedraw3D(kTRUE);

  TGLViewer *glv1 = multiView->Get3DView()->GetGLViewer();
  TGLViewer *glv2 = multiView->GetRPhiView()->GetGLViewer();
  TGLViewer *glv3 = multiView->GetRhoZView()->GetGLViewer();

  glv1->CurrentCamera().RotateRad(-0.4, -1.8);
  glv2->CurrentCamera().Dolly(450, kFALSE, kFALSE);
  glv3->CurrentCamera().Dolly(1500, kFALSE, kFALSE);

  gEve->FullRedraw3D();

}

Int_t      g_pic_id  = 0;
Int_t      g_pic_max = 100;
TTimeStamp g_pic_prev(0, 0);

void alieve_online_on_new_event()
{
  AliSysInfo::AddStamp("on_new_event_start");

  AliTriggerAnalysis atr;

  AliESDEvent* esd = AliEveEventManager::AssertESD();
  Double_t x[3] = { 0, 0, 0 };
  esd->GetPrimaryVertex()->GetXYZ(x);

  TEveElement* top = gEve->GetCurrentEvent();

  AliEveMultiView *multiView = AliEveMultiView::Instance();

/*
  TGLViewer *glv = (dynamic_cast<TEveViewer*>(gEve->GetViewers()->FindChild("3D View")))->GetGLViewer();
  
  if(gEve->GetScenes()->FirstChild()->FindChild("Gentle MUON"))
  {
    if (esd->GetNumberOfMuonTracks() == 0 && !gEve->GetKeepEmptyCont())
    {
      gEve->GetScenes()->FirstChild()->FindChild("Gentle MUON")->SetRnrChildren(kFALSE);
      
      if(gEve->GetEventScene()->FirstChild()->FindChild("MUON Clusters"))
        gEve->GetEventScene()->FirstChild()->FindChild("MUON Clusters")->SetRnrSelf(kFALSE);
      if(gEve->GetEventScene()->FirstChild()->FindChild("MUON Raw digits"))
        gEve->GetEventScene()->FirstChild()->FindChild("MUON Raw digits")->SetRnrChildren(kFALSE);

      gEve->FullRedraw3D(kTRUE);
      glv->CurrentCamera().RotateRad(-0.4, -1.8);
    }
    else
    {
      gEve->GetScenes()->FirstChild()->FindChild("Gentle MUON")->SetRnrChildren(kTRUE);

      if(gEve->GetEventScene()->FirstChild()->FindChild("MUON Clusters"))
        gEve->GetEventScene()->FirstChild()->FindChild("MUON Clusters")->SetRnrSelf(kTRUE);
      if(gEve->GetEventScene()->FirstChild()->FindChild("MUON Raw digits"))
        gEve->GetEventScene()->FirstChild()->FindChild("MUON Raw digits")->SetRnrChildren(kTRUE);

      gEve->FullRedraw3D(kTRUE);
      glv->CurrentCamera().RotateRad(-0.4, 1);
    }
  }

  glv->DoDraw();
*/

  AliTriggerAnalysis::V0Decision decisionV0a = 
    atr.V0Trigger(esd, AliTriggerAnalysis::kASide, kFALSE);
  AliTriggerAnalysis::V0Decision decisionV0c = 
    atr.V0Trigger(esd, AliTriggerAnalysis::kCSide, kFALSE);

  Double_t a = 0;
  Double_t c = 0;

  if( decisionV0a == AliTriggerAnalysis::kV0Invalid ) a = 0.5;
  if( decisionV0a == AliTriggerAnalysis::kV0Empty ) a = 1.5;
  if( decisionV0a == AliTriggerAnalysis::kV0BB ) a = 2.5;
  if( decisionV0a == AliTriggerAnalysis::kV0BG ) a = 3.5;

  if( decisionV0c == AliTriggerAnalysis::kV0Invalid ) c = 0.5;
  if( decisionV0c == AliTriggerAnalysis::kV0Empty ) c = 1.5;
  if( decisionV0c == AliTriggerAnalysis::kV0BB ) c = 2.5;
  if( decisionV0c == AliTriggerAnalysis::kV0BG ) c = 3.5;

  V0StateHistogram->Fill(a,c);
  AliSysInfo::AddStamp("on_new_event_after_trig");

  TGLViewer *glv = multiView->Get3DView()->GetGLViewer();
  TGLViewer *glv1 = multiView->GetRPhiView()->GetGLViewer();
  TGLViewer *glv2 = multiView->GetRhoZView()->GetGLViewer();
 
  Double_t RPhiCameraFrustrumCenter = TMath::Sqrt(glv1->CurrentCamera().FrustumCenter().X()*glv1->CurrentCamera().FrustumCenter().X() + glv1->CurrentCamera().FrustumCenter().Y()*glv1->CurrentCamera().FrustumCenter().Y());

  Double_t RhoZCameraFrustrumCenter = TMath::Sqrt(glv2->CurrentCamera().FrustumCenter().X()*glv2->CurrentCamera().FrustumCenter().X() + glv2->CurrentCamera().FrustumCenter().Y()*glv2->CurrentCamera().FrustumCenter().Y());

  if(RPhiCameraFrustrumCenter > 500 || RhoZCameraFrustrumCenter > 500)
  {

    glv->ResetCurrentCamera();
    glv1->ResetCurrentCamera();
    glv2->ResetCurrentCamera();

    glv->CurrentCamera().RotateRad(-0.4, -1.8);
    glv1->CurrentCamera().Dolly(450, kFALSE, kFALSE);
    glv2->CurrentCamera().Dolly(1500, kFALSE, kFALSE);

    gEve->FullRedraw3D();

  }

  multiView->DestroyEventRPhi();
  if (gCenterProjectionsAtPrimaryVertex)
    multiView->SetCenterRPhi(x[0], x[1], x[2]);
  multiView->ImportEventRPhi(top);

  multiView->DestroyEventRhoZ();
  if (gCenterProjectionsAtPrimaryVertex)
    multiView->SetCenterRhoZ(x[0], x[1], x[2]);
  multiView->ImportEventRhoZ(top);

  AliSysInfo::AddStamp("on_new_event_after_rozphi");

  if(multiView->IsMuonView()) { multiView->DestroyEventMuon(); multiView->ImportEventMuon(top); }

  if(gEve->GetScenes()->FindChild("Event scene")->FindChild("Online Event")->FindChild("MUON Clusters"))
  {
    gEve->GetScenes()->FindChild("Event scene")->FindChild("Online Event")->FindChild("MUON Clusters")->SetRnrSelf(kFALSE);
    gEve->GetScenes()->FindChild("Event scene")->FindChild("Online Event")->FindChild("MUON Clusters")->SetRnrChildren(kFALSE);
  }

  if(gEve->GetScenes()->FindChild("Event scene")->FindChild("Online Event")->FindChild("ESD MUON Clusters"))
  {
    gEve->GetScenes()->FindChild("Event scene")->FindChild("Online Event")->FindChild("ESD MUON Clusters")->SetRnrSelf(kFALSE);
    gEve->GetScenes()->FindChild("Event scene")->FindChild("Online Event")->FindChild("ESD MUON Clusters")->SetRnrChildren(kFALSE);
  }

  if(gEve->GetScenes()->FindChild("Event scene")->FindChild("Online Event")->FindChild("ESD MUON Tracks"))
  {
    gEve->GetScenes()->FindChild("Event scene")->FindChild("Online Event")->FindChild("ESD MUON Tracks")->SetRnrSelf(kFALSE);
    gEve->GetScenes()->FindChild("Event scene")->FindChild("Online Event")->FindChild("ESD MUON Tracks")->SetRnrChildren(kFALSE);
  }

  if(gEve->GetScenes()->FindChild("Event scene")->FindChild("Online Event")->FindChild("MUON Raw digits"))
  {
    gEve->GetScenes()->FindChild("Event scene")->FindChild("Online Event")->FindChild("MUON Raw digits")->SetRnrSelf(kFALSE);
    gEve->GetScenes()->FindChild("Event scene")->FindChild("Online Event")->FindChild("MUON Raw digits")->SetRnrChildren(kFALSE);
  }

  if(gEve->GetScenes()->FindChild("Muon Event Data"))
  {
    if(gEve->GetScenes()->FindChild("Muon Event Data")->FindChild("Online Event [P]"))
    {
      if(gEve->GetScenes()->FindChild("Muon Event Data")->FindChild("Online Event [P]")->FindChild("MUON Clusters [P]"))
      {
	gEve->GetScenes()->FindChild("Muon Event Data")->FindChild("Online Event [P]")->FindChild("MUON Clusters [P]")->SetRnrSelf(kTRUE);
	gEve->GetScenes()->FindChild("Muon Event Data")->FindChild("Online Event [P]")->FindChild("MUON Clusters [P]")->SetRnrChildren(kTRUE);
      }

      if(gEve->GetScenes()->FindChild("Muon Event Data")->FindChild("Online Event [P]")->FindChild("ESD MUON Clusters [P]"))
      {
	gEve->GetScenes()->FindChild("Muon Event Data")->FindChild("Online Event [P]")->FindChild("ESD MUON Clusters [P]")->SetRnrSelf(kTRUE);
	gEve->GetScenes()->FindChild("Muon Event Data")->FindChild("Online Event [P]")->FindChild("ESD MUON Clusters [P]")->SetRnrChildren(kTRUE);
      }

      if(gEve->GetScenes()->FindChild("Muon Event Data")->FindChild("Online Event [P]")->FindChild("ESD MUON Tracks [P]"))
      {
	gEve->GetScenes()->FindChild("Muon Event Data")->FindChild("Online Event [P]")->FindChild("ESD MUON Tracks [P]")->SetRnrSelf(kTRUE);
	gEve->GetScenes()->FindChild("Muon Event Data")->FindChild("Online Event [P]")->FindChild("ESD MUON Tracks [P]")->SetRnrChildren(kTRUE);
      }
    }
  }

  if(gEve->GetScenes()->FindChild("RhoZ Event Data"))
  {
    if(gEve->GetScenes()->FindChild("RhoZ Event Data")->FindChild("Online Event [P]"))
    {
      if(gEve->GetScenes()->FindChild("RhoZ Event Data")->FindChild("Online Event [P]")->FindChild("FMD [P]"))
      {
	gEve->GetScenes()->FindChild("RhoZ Event Data")->FindChild("Online Event [P]")->FindChild("FMD [P]")->SetRnrSelf(kFALSE);
	gEve->GetScenes()->FindChild("RhoZ Event Data")->FindChild("Online Event [P]")->FindChild("FMD [P]")->SetRnrChildren(kFALSE);
      }
    }
  }

  if(gEve->GetScenes()->FindChild("RPhi Event Data"))
  {
    if(gEve->GetScenes()->FindChild("RPhi Event Data")->FindChild("Online Event [P]"))
    {
      if(gEve->GetScenes()->FindChild("RPhi Event Data")->FindChild("Online Event [P]")->FindChild("FMD [P]"))
      {
        gEve->GetScenes()->FindChild("RPhi Event Data")->FindChild("Online Event [P]")->FindChild("FMD [P]")->SetRnrSelf(kFALSE);
        gEve->GetScenes()->FindChild("RPhi Event Data")->FindChild("Online Event [P]")->FindChild("FMD [P]")->SetRnrChildren(kFALSE);
      }
    }
  }
  AliSysInfo::AddStamp("on_new_event_after_muon");

  gEve->FullRedraw3D();

  // Register image to amore.
  // const TString pichost("aldaqacrs3");
  const TString pichost(gEnv->GetValue("AliEve.imageDumpHost", "aldaqacrs3"));
  TTimeStamp now;
  Double_t delta = now.AsDouble() - g_pic_prev.AsDouble();

  printf("Pre image dump: host='%s', delta=%f.\n",
	 gSystem->HostName(), delta);

  AliSysInfo::AddStamp("on_new_event_pic");
  if (pichost == gSystem->HostName() && delta >= 30)
  {
    TString id;      id.Form("online-viz-%03d", g_pic_id);
    TString pic(id); pic += ".png";

    printf("In image dump: file='%s'.\n", pic.Data());

    gEve->GetBrowser()->RaiseWindow();
    gEve->FullRedraw3D();
    gSystem->ProcessEvents();

    Int_t status;

    // create screenshots from OpenGL views
    TEveUtil::LoadMacro("saveViews.C+");
    saveViews(pic.Data()); 

    // send screenshot to AMORE
    status = gSystem->Exec(TString::Format("SendImageToAmore %s %s %d",
		          id.Data(), pic.Data(),
		          AliEveEventManager::AssertRawReader()->GetRunNumber()));

    printf("Post AMORE reg -- status=%d, run=%d.\n", status,
	   AliEveEventManager::AssertRawReader()->GetRunNumber());

    if (++g_pic_id >= g_pic_max)
      g_pic_id = 0;
    g_pic_prev.Set();
  }
  AliSysInfo::AddStamp("on_new_event_end");

}
