/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id$

#include "AliMUONMchViewApplication.h"

#include "AliCDBManager.h"
#include "AliCodeTimer.h"
#include "AliLog.h"
#include "AliMUONPainterDataSourceFrame.h"
#include "AliMUONPainterEnv.h"
#include "AliMUONPainterHelper.h"
#include "AliMUONPainterMasterFrame.h"
#include "AliMUONPainterRegistry.h"
#include "AliMUONTrackerDataCompareDialog.h"
#include "AliMUONTrackerDataWrapper.h"
#include "AliMUONVTrackerData.h"
#include "AliMUONVTrackerDataMaker.h"
#include <Riostream.h>
#include <TCanvas.h>
#include <TEnv.h>
#include <TFile.h>
#include <TGClient.h>
#include <TGFileDialog.h>
#include <TGMenu.h>
#include <TGTab.h>
#include <TGTextView.h>
#include <TGrid.h>
#include <TKey.h>
#include <TList.h>
#include <TRegexp.h>
#include <TString.h>
#include <TSystem.h>

/// \class AliMUONMchViewApplication
///
/// Main class for the mchview program
///
///\author Laurent Aphecetche, Subatech

/// \cond CLASSIMP
ClassImp(AliMUONMchViewApplication)
/// \endcond CLASSIMP

const Int_t AliMUONMchViewApplication::fgkFILESAVEAS(1);
const Int_t AliMUONMchViewApplication::fgkFILEOPEN(2);
const Int_t AliMUONMchViewApplication::fgkFILEEXIT(3);
const Int_t AliMUONMchViewApplication::fgkFILEPRINTAS(4);
const Int_t AliMUONMchViewApplication::fgkABOUT(5);
const Int_t AliMUONMchViewApplication::fgkCOMPAREDATA(6);

const char* AliMUONMchViewApplication::fgkFileTypes[] = { 
  "ROOT files",    "*.root", 
  "All files",     "*", 
  0,               0 }; 

//______________________________________________________________________________
AliMUONMchViewApplication::AliMUONMchViewApplication(const char* name,
                                                     int* argc, char** argv,
                                                     UInt_t w, UInt_t h,
                                                     UInt_t ox, UInt_t oy) 
: TRint(name,argc,argv),
  fMainFrame(0x0),
  fPainterMasterFrame(0x0)
{

  /// ctor
  /// (w,h) is the size in pixel (if 0,0 it will be computed as 70%,90% of display size)
  /// (ox,oy) is the offset from the top-left of the display

  if (!w | !h)
  {
    w = (UInt_t)(gClient->GetDisplayWidth()*0.7);
    h = (UInt_t)(gClient->GetDisplayHeight()*0.9); 
  }

  fMainFrame = new TGMainFrame(gClient->GetRoot(),w,h);
  
  CreateMenuBar(w);

  const Int_t kbs = 2;
  
//  h -= 60; // menubar
  
  TGTab* tabs = new TGTab(fMainFrame,w,h);
  
  TGCompositeFrame* t = tabs->AddTab("Painter Master Frame");

  fPainterMasterFrame =
    new AliMUONPainterMasterFrame(t,t->GetWidth()-kbs*2,t->GetHeight()-kbs*2);
  
  t->AddFrame(fPainterMasterFrame, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,kbs,kbs,kbs,kbs));

  t = tabs->AddTab("Data Sources");
  
  AliMUONPainterDataSourceFrame* dsf = 
    new AliMUONPainterDataSourceFrame(t,t->GetWidth()-kbs*2,t->GetHeight()-kbs*2);
  
  t->AddFrame(dsf,new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,kbs,kbs,kbs,kbs));
  
  fMainFrame->AddFrame(tabs,new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,0,0,0,0));

  fMainFrame->SetWindowName("mchview - Visualization of MUON Tracker detector");

  fMainFrame->MapSubwindows();
  fMainFrame->Resize();
  
  fPainterMasterFrame->Update();
  
  fMainFrame->MapWindow();
  
  fMainFrame->Connect("CloseWindow()","AliMUONMchViewApplication",this,"Terminate()");

  fMainFrame->MoveResize(ox,oy, w, h); 
  fMainFrame->SetWMPosition(ox, oy);
  fMainFrame->SetWMSizeHints(w,h,w,h,0,0);
  
  cout << "***************************************************" << endl;
  cout << "   Welcome to mchview" << endl;
  cout << "   " << FullVersion() << endl;
  cout << "***************************************************" << endl;

}

//______________________________________________________________________________
AliMUONMchViewApplication::~AliMUONMchViewApplication()
{
  /// dtor
}

//______________________________________________________________________________
void
AliMUONMchViewApplication::CompareData()
{
  /// Launch compare data dialog
  TGTransientFrame* t = new AliMUONTrackerDataCompareDialog(gClient->GetRoot(),
                                                            gClient->GetRoot(),
                                                            400,400);

  t->MapSubwindows();
  t->Resize();
  t->MapWindow();
  t->CenterOnParent();
  
  // set names
  
  t->SetWindowName("mchview compare data tool");
  t->SetIconName("mchview compare data tool");
  
  t->MapRaised();
  
}

//______________________________________________________________________________
void
AliMUONMchViewApplication::CreateMenuBar(UInt_t w)
{
  /// Create the application menu bar
  
  TGPopupMenu* file = new TGPopupMenu(gClient->GetRoot());
  
  file->AddEntry("&Open...",fgkFILEOPEN);
  file->AddEntry("&Save As...",fgkFILESAVEAS);
  file->AddEntry("&Print As...",fgkFILEPRINTAS);
  file->AddEntry("&Exit",fgkFILEEXIT);
  
  TGMenuBar* bar = new TGMenuBar(fMainFrame,w);
  
  TGPopupMenu* tools = new TGPopupMenu(gClient->GetRoot());
  tools->AddEntry("&Compare data",fgkCOMPAREDATA);
  
  TGPopupMenu* about = new TGPopupMenu(gClient->GetRoot());  
  about->AddEntry(FullVersion(),fgkABOUT);

  file->Connect("Activated(Int_t)","AliMUONMchViewApplication",this,"HandleMenu(Int_t)");
  about->Connect("Activated(Int_t)","AliMUONMchViewApplication",this,"HandleMenu(Int_t)");
  tools->Connect("Activated(Int_t)","AliMUONMchViewApplication",this,"HandleMenu(Int_t)");
  
  bar->AddPopup("&File",file,new TGLayoutHints(kLHintsLeft|kLHintsTop));
  bar->AddPopup("&Tools",tools,new TGLayoutHints(kLHintsLeft|kLHintsTop));
  bar->AddPopup("&About",about,new TGLayoutHints(kLHintsRight|kLHintsTop));
  
  fMainFrame->AddFrame(bar,new TGLayoutHints(kLHintsLeft|kLHintsExpandX));
  
  AliMUONPainterRegistry::Instance()->SetMenuBar(bar);
}

//______________________________________________________________________________
void
AliMUONMchViewApplication::HandleMenu(Int_t i)
{
  /// Handle the click of one menu item

  switch (i)
    {
    case fgkFILEEXIT:
      Terminate(1);
      break;
    case fgkFILEOPEN:
      Open();
      break;
    case fgkFILESAVEAS:
      Save();
      break;
    case fgkFILEPRINTAS:
      PrintAs();
      break;
    case fgkABOUT:
      ReleaseNotes();
      break;
    case fgkCOMPAREDATA:
      CompareData();
      break;
    default:
      break;
    }
}

//______________________________________________________________________________
void
AliMUONMchViewApplication::Open()
{
  /// Open file dialog
  
  TGFileInfo fileInfo;
  
  fileInfo.fFileTypes = fgkFileTypes;
  
  delete[] fileInfo.fIniDir;
  
  AliMUONPainterEnv* env = AliMUONPainterHelper::Instance()->Env();
  
  fileInfo.fIniDir = StrDup(env->String("LastOpenDir","."));
  
  new TGFileDialog(gClient->GetRoot(),gClient->GetRoot(),
                   kFDOpen,&fileInfo);

  env->Set("LastOpenDir",fileInfo.fIniDir);
  env->Save();  
    
  Open(gSystem->ExpandPathName(Form("%s",fileInfo.fFilename)));
}  

//______________________________________________________________________________
void
AliMUONMchViewApplication::Open(const char* filename)
{
  /// Open a given file containing saved VTrackerDataMaker objects
  
  TString sfilename(gSystem->ExpandPathName(filename));
  
  if ( sfilename.Contains(TRegexp("^alien")) )
  {
    // insure we've initialized the grid...
    if (!gGrid)
    {
      TGrid::Connect("alien://");
    }
  }
  
  TFile* f = TFile::Open(filename);
  
	ReadDir(*f);
	
	delete f;
}

//______________________________________________________________________________
void
AliMUONMchViewApplication::ReadDir(TDirectory& dir)
{
  /// Read the given directory and import VTrackerData objects found
  
  TList* keys = dir.GetListOfKeys();
  TIter next(keys);
  
  TKey* k;
  
  while ( ( k = static_cast<TKey*>(next()) ) )
  {
    TObject* object = k->ReadObj();

		if ( object->InheritsFrom("TDirectory") )
		{
			TDirectory* d = static_cast<TDirectory*>(object);
			ReadDir(*d);
			continue;
		}
		
    if ( object->InheritsFrom("AliMUONVTrackerDataMaker") )
    {
      AliMUONVTrackerDataMaker* maker = dynamic_cast<AliMUONVTrackerDataMaker*>(object);
      if ( maker ) 
      {
        AliMUONPainterRegistry::Instance()->Register(maker);
      }
    }
    
    if ( object->InheritsFrom("AliMUONVTrackerData") )
    {
      // this is for backward compatibility. Early versions of mchview 
      // wrote VTrackerData objects, and not VTrackerDataMaker ones.
      
      AliMUONVTrackerData* data = dynamic_cast<AliMUONVTrackerData*>(object);
      if ( data ) 
      {
        AliMUONVTrackerDataMaker* maker = new AliMUONTrackerDataWrapper(data);
        AliMUONPainterRegistry::Instance()->Register(maker);
      }
    }
  }
  
} 

//______________________________________________________________________________
void
AliMUONMchViewApplication::PrintAs()
{
  /// Print as...
  
  TGFileInfo fileInfo;
  
  new TGFileDialog(gClient->GetRoot(),gClient->GetRoot(),
                   kFDSave,&fileInfo);
  
  fPainterMasterFrame->SaveAs(gSystem->ExpandPathName(Form("%s",fileInfo.fFilename)));
}

//______________________________________________________________________________
void
AliMUONMchViewApplication::ReleaseNotes()
{
  /// Display release notes
  
  UInt_t width = 600;
  UInt_t height = 400;
  
  TGTransientFrame* t = new TGTransientFrame(gClient->GetRoot(),gClient->GetRoot(),width,height);
  
  TGTextView* rn = new TGTextView(t);

  rn->AddLine("0.9%");
  rn->AddLine("");
  rn->AddLine("New features");
  rn->AddLine("");
  rn->AddLine("- Can now read and display HV values from OCDB");
  rn->AddLine("- New program option --geometry to force geometry of the window");
  rn->AddLine("- Added possibility, in painters' context menu, to include or exclude part of the detector");
  rn->AddLine("  (which will be used later on to communicate with LC2 which parts should be read out or not)");
  rn->AddLine("");
  rn->AddLine("Improvement");
  rn->AddLine("");
  rn->AddLine("- When displaying Gains, the quality information is now decoded");
  rn->AddLine("");
  
  rn->AddLine("0.94");
  rn->AddLine("");
  rn->AddLine("New features");
  rn->AddLine("");
  rn->AddLine("Can now read ASCII calibration files produced by the DA");
  rn->AddLine("");
  
  rn->AddLine("0.93");
  rn->AddLine("");
  rn->AddLine("New features");
  rn->AddLine("");
  rn->AddLine("- Adding a Lock button under the color slider to lock the range shown");
  rn->AddLine("  when switching between views");
  rn->AddLine("- Default display now shows bending plane (instead of cathode 0 before)");
  rn->AddLine("- If pad is responder and there's some histo for that pad, ");
  rn->AddLine("  clicking on it will display an histo");
  rn->AddLine("- Right-click on a painter will now display several histogram options");
  rn->AddLine("  (e.g. raw charge as before, but also simple distributions of mean");
  rn->AddLine("  and sigma");
  rn->AddLine("- In the Data Sources Tab, each data source can now be removed and saved");
  rn->AddLine("- There's a new Tool menu which allow to produce a TrackerData from two others");
  rn->AddLine("  in order to compare data.");
  rn->AddLine("  - The --use option can now reference alien files");
  rn->AddLine("");    
  rn->AddLine("Bug fixes");
  rn->AddLine("");    
  rn->AddLine("- Can now read Capacitances from OCDB");
    
  rn->Resize(width,height);
  
  t->AddFrame(rn, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
  
  t->MapSubwindows();
  t->Resize();
  t->MapWindow();
  t->CenterOnParent();
    
  // set names
  
  t->SetWindowName("mchview release notes");
  t->SetIconName("mchview release notes");
  
//  t->SetMWMHints(kMWMDecorAll | kMWMDecorResizeH  | kMWMDecorMaximize |
//              kMWMDecorMinimize | kMWMDecorMenu,
//              kMWMFuncAll  | kMWMFuncResize    | kMWMFuncMaximize |
//              kMWMFuncMinimize,
//              kMWMInputModeless);
  
  t->MapRaised();
//  gClient->WaitFor(t);
}

//______________________________________________________________________________
void
AliMUONMchViewApplication::Save()
{
  /// Open "Save VTrackerData objects to file" dialog
  
  TGFileInfo fileInfo;
  
  new TGFileDialog(gClient->GetRoot(),gClient->GetRoot(),
                   kFDSave,&fileInfo);
  
  Save(gSystem->ExpandPathName(Form("%s",fileInfo.fFilename)));
}  

//______________________________________________________________________________
void
AliMUONMchViewApplication::Save(const char* filename)
{
  /// Save VTrackerDataMaker objects into file of given name
  
  AliMUONPainterRegistry* reg = AliMUONPainterRegistry::Instance();

  TFile f(filename,"RECREATE");
  
  for ( Int_t i = 0; i < reg->NumberOfDataMakers(); ++i )
  {
    AliMUONVTrackerDataMaker* maker = reg->DataMaker(i);
    maker->Write();
  }
  
  f.Close();
}  
