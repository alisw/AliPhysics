//
//  AliEveHLT.cxx
//  xAliRoot
//
//  Created by Jeremi Niedziela on 11/05/15.
//
//

#include "AliEveHLT.h"
#include "AliEveGeomGentle.h"
#include "AliEveHLTZMQeventManager.h"
#include "AliEveEventManagerEditor.h"
#include "AliEveMultiView.h"
#include "AliEveMacroExecutor.h"
#include "AliEveMacro.h"

#include <TROOT.h>
#include <TEveManager.h>
#include <TEveBrowser.h>
#include <TBrowser.h>
#include <TEveMacro.h>
#include <TGTab.h>
#include <TGFileBrowser.h>
#include <TInterpreter.h>
#include <TPRegexp.h>
#include <TFolder.h>
#include <TSystemDirectory.h>
#include "AliCDBManager.h"
#include "TTimeStamp.h"

#include <iostream>

using namespace std;

AliEveHLT::AliEveHLT(bool storageManager)
{
  cout<<"Creating AliEveHLT"<<endl;

  //-----------------------------------------------------------------------------------------
  //  Set all preferences here
  //
  Color_t colorTRD = kGray;      // color of TRD modules
  Color_t colorMUON = kGray;     // color of MUON modules

  bool customPreset = true;            // should one of the following custom presets be used
  Color_t colors[9] = {kCyan,kCyan,kCyan,kCyan,kCyan,kCyan,kCyan,kCyan,kCyan};
  Width_t widths[9] = {3,3,3,3,3,3,3,3,3};
  bool dashBad = true;

  //    Color_t colors[9] = {kGreen,kGreen,kGreen,kGreen,kGreen,kGreen,kGreen,kGreen,kGreen}; // preset for cosmics

  bool saveViews = false;          // should screenshot be saved and sent to ALICE LIVE
  

  //
  //-----------------------------------------------------------------------------------------


  TString ocdbStorage;
  if (gSystem->Getenv("ocdbStorage"))
    ocdbStorage=gSystem->Getenv("ocdbStorage");
  AliEveEventManager::SetCdbUri(ocdbStorage);         // current OCDB snapshot
  AliCDBManager* cdbManager = AliCDBManager::Instance();
  cdbManager->SetDefaultStorage(ocdbStorage);

  AliEveEventManager* eventManager = new AliEveHLTZMQeventManager();
  gEve->AddEvent(eventManager);
  cout<<"Event manager created"<<endl;

  cout<<"Creating multiview...";
  AliEveMultiView *multiView = new AliEveMultiView(kTRUE);
  cout<<"created"<<endl;

  Info("alieve_init", "Adding standard macros.");
  InitImportMacros();
  gROOT->ProcessLine(".L geom_gentle.C+");
  gROOT->ProcessLine(".L geom_gentle_trd.C+");
  gROOT->ProcessLine(".L geom_gentle_muon.C+");
  TEveUtil::LoadMacro("saveViews.C");
  cout<<"Standard macros added"<<endl;

  TEveUtil::AssertMacro("VizDB_scan.C");
  gSystem->ProcessEvents();
  cout<<"VizDB_scan loaded"<<endl;
  TEveBrowser *browser = gEve->GetBrowser();
  browser->ShowCloseTab(kFALSE);
  cout<<"browser created"<<endl;

  TEveUtil::LoadMacro("geom_gentle.C");
  cout<<"geom gentle loaded"<<endl;

  //multiView->InitGeomGentle(geom_gentle(),
  //    geom_gentle_rphi(),
  //    geom_gentle_rhoz(),
  //    geom_gentle_rhoz());
  //cout<<"geom gentl inited"<<endl;
  //multiView->InitGeomGentleTrd(geom_gentle_trd());
  //multiView->InitGeomGentleMuon(geom_gentle_muon(), kFALSE, kFALSE, kTRUE);

  printf("============ Setting macro executor ============\n");

  AliEveMacroExecutor *exec = AliEveEventManager::GetMaster()->GetExecutor();
  printf("exec created\n");

  // default appearance:
  exec->AddMacro(new AliEveMacro(AliEveMacro::kESD, "REC Tracks by category",  "esd_tracks.C", "esd_tracks_by_category",  "", kTRUE));

  cout<<"macros added to exec"<<endl;

  //============================================================================
  // Final GUI setup
  //============================================================================

  browser->GetTabRight()->SetTab(1);
  browser->StartEmbedding(TRootBrowser::kBottom);
  new AliEveEventManagerWindow(AliEveEventManager::GetMaster());
  browser->StopEmbedding("EventCtrl");

  gEve->FullRedraw3D(kTRUE);
  gSystem->ProcessEvents();

  // move and rotate sub-views
  TGLViewer *glv1 = multiView->Get3DView()->GetGLViewer();
  TGLViewer *glv2 = multiView->GetRPhiView()->GetGLViewer();
  TGLViewer *glv3 = multiView->GetRhoZView()->GetGLViewer();
  TGLViewer *glv4 = multiView->GetMuonView()->GetGLViewer();

  glv1->CurrentCamera().RotateRad(-0.4, 0.6);
  glv2->CurrentCamera().Dolly(90, kFALSE, kFALSE);
  glv3->CurrentCamera().Dolly(2300, kFALSE, kFALSE);
  glv4->CurrentCamera().Dolly(1, kFALSE, kFALSE);

  AliEveEventManager::GetMaster()->AddNewEventCommand("AliEveHLT::OnNewEvent();");

  if(customPreset)
  {
    eventManager->SetESDcolors(colors);
    eventManager->SetESDwidths(widths);
    eventManager->SetESDdashBad(dashBad);
  }

  eventManager->SetSaveViews(saveViews);
  eventManager->SetESDtracksByCategory(true);

  gEve->FullRedraw3D();
  gSystem->ProcessEvents();
  gEve->Redraw3D(kTRUE);

}

AliEveHLT::~AliEveHLT()
{
}

void AliEveHLT::OnNewEvent()
{
  if (!AliEveEventManager::HasESD()) return;

  AliESDEvent* esd = AliEveEventManager::AssertESD();
	
  Double_t x[3] = { 0, 0, 0 };
  esd->GetPrimaryVertex()->GetXYZ(x);
  
  TTimeStamp ts(esd->GetTimeStamp());
  TString win_title("Eve Main Window -- Timestamp: ");
  win_title += ts.AsString("s");
  win_title += "; Event # in ESD file: ";
  win_title += esd->GetEventNumberInFile();
  win_title += " | run: ";
  win_title += esd->GetRunNumber();
  gEve->GetBrowser()->SetWindowName(win_title);
    
	TEveElement* top = gEve->GetCurrentEvent();
    
	AliEveMultiView *mv = AliEveMultiView::Instance();
    
	mv->DestroyEventRPhi();
}

void AliEveHLT::InitImportMacros()
{
    // Put macros in the list of browsables, add a macro browser to
  // top-level GUI.

  TString macdir("$(ALICE_ROOT)/EVE/alice-macros");

  if (gSystem->Getenv("ALICE_ROOT") != 0)
  {
    gInterpreter->AddIncludePath(Form("%s/MUON", gSystem->Getenv("ALICE_ROOT")));
    gInterpreter->AddIncludePath(Form("%s/MUON/mapping", gSystem->Getenv("ALICE_ROOT")));
    gSystem->ExpandPathName(macdir);
  }

  TFolder* f = gEve->GetMacroFolder();
  void* dirhandle = gSystem->OpenDirectory(macdir.Data());
  if (dirhandle != 0)
  {
    char* filename;
    TPMERegexp re("\\.C$");
    TObjArray names;
    while ((filename = (char*)(gSystem->GetDirEntry(dirhandle))) != 0)
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
