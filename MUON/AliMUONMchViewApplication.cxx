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
#include "AliMUONChamberPainter.h"
#include "AliMUONDEPainter.h"
#include "AliMUONPainterDataRegistry.h"
#include "AliMUONPainterDataSourceFrame.h"
#include "AliMUONPainterEnv.h"
#include "AliMUONPainterHelper.h"
#include "AliMUONPainterGroup.h"
#include "AliMUONPainterMasterFrame.h"
#include "AliMUONPainterMatrix.h"
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
  fPainterMasterFrameList(new TList),
  fTabs(0x0)
{

  /// ctor
  /// (w,h) is the size in pixel (if 0,0 it will be computed as 70%,90% of display size)
  /// (ox,oy) is the offset from the top-left of the display

  if (!w || !h)
  {
    w = (UInt_t)(gClient->GetDisplayWidth()*0.7);
    h = (UInt_t)(gClient->GetDisplayHeight()*0.9); 
  }

  fMainFrame = new TGMainFrame(gClient->GetRoot(),w,h);

  CreateMenuBar(w);

  const Int_t kbs = 2;
  
//  h -= 60; // menubar
  
  fTabs = new TGTab(fMainFrame,w,h);
  
  TGCompositeFrame* t = fTabs->AddTab("Painter Master Frame");

  fPainterMasterFrameList->SetOwner(kTRUE);
  
  
  AliMUONPainterMasterFrame* pmf = new AliMUONPainterMasterFrame(t,t->GetWidth()-kbs*2,t->GetHeight()-kbs*2,
                                                                 GenerateStartupMatrix());

  fPainterMasterFrameList->Add(pmf);
  
  t->AddFrame(pmf, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,kbs,kbs,kbs,kbs));

  t = fTabs->AddTab("Data Sources");
  
  AliMUONPainterDataSourceFrame* dsf = 
    new AliMUONPainterDataSourceFrame(t,t->GetWidth()-kbs*2,t->GetHeight()-kbs*2);
  
  t->AddFrame(dsf,new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,kbs,kbs,kbs,kbs));
  
  fMainFrame->AddFrame(fTabs,new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,0,0,0,0));

  fMainFrame->SetWindowName("mchview - Visualization of MUON Tracker detector");

  fMainFrame->MapSubwindows();
  fMainFrame->Resize();
  
  pmf->Update();
  
  fMainFrame->MapWindow();
  
  fMainFrame->Connect("CloseWindow()","AliMUONMchViewApplication",this,"Terminate()");

//  fMainFrame->MoveResize(ox,oy, w, h); 
  fMainFrame->SetWMPosition(ox, oy);
//  fMainFrame->SetWMSizeHints(w,h,w,h,0,0);
//  fMainFrame->SetWMSizeHints(w,h,w,h,10,10);
  
  cout << "***************************************************" << endl;
  cout << "   Welcome to mchview" << endl;
  cout << "   " << FullVersion() << endl;
  cout << "***************************************************" << endl;
  
  // Trying to see if we're requested to draw something specific instead
  // of the global view of all the chambers
  
  AliMUONVPainter* painter(0x0);
  TObjArray args;
  args.SetOwner(kTRUE);
  
  for ( int i = 1; i < argc[0]; ++i ) 
  {
    args.Add(new TObjString(argv[i]));
  }
  
  for ( Int_t i = 0; i <= args.GetLast(); ++i ) 
  {
    TString a(static_cast<TObjString*>(args.At(i))->String());

    AliMUONAttPainter att;
    
    att.SetPlane(kTRUE,kFALSE);
    att.SetCathode(kFALSE,kFALSE);
    att.SetViewPoint(kTRUE,kFALSE);
        
    if ( a == "--de" )
    {
      Int_t detElemId = static_cast<TObjString*>(args.At(i+1))->String().Atoi();
      
      painter = new AliMUONDEPainter(att,detElemId);
      
      painter->SetOutlined("*",kFALSE);      
      painter->SetOutlined("BUSPATCH",kTRUE);

      painter->SetLine(1,4,3);      
      ++i;
    }

    if ( a == "--chamber" )
    {
      Int_t chamberId = static_cast<TObjString*>(args.At(i+1))->String().Atoi();
      
      painter = new AliMUONChamberPainter(att,chamberId-1);
      
      painter->SetOutlined("*",kFALSE);      
      painter->SetOutlined("DE",kTRUE);
      
      painter->SetLine(1,4,3);      
      ++i;
    }
    
  }
  
  if ( painter ) 
  {
    pmf->ShiftClicked(painter,0x0);
    
    pmf->Update();
  }
      
}

//______________________________________________________________________________
AliMUONMchViewApplication::~AliMUONMchViewApplication()
{
  /// dtor
  delete fPainterMasterFrameList;
}

//_____________________________________________________________________________
AliMUONPainterMatrix*
AliMUONMchViewApplication::GenerateStartupMatrix()
{
  /// Kind of bootstrap method to trigger the generation of all contours
  
  AliCodeTimerAuto("",0);
  
  AliMUONAttPainter att;
  
  att.SetViewPoint(kTRUE,kFALSE);
  att.SetCathode(kFALSE,kFALSE);
  att.SetPlane(kTRUE,kFALSE);

  AliMUONPainterMatrix* matrix = new AliMUONPainterMatrix("Tracker",5,2);
    
  for ( Int_t i = 0; i < 10; ++i )
  {
    AliMUONVPainter* painter = new AliMUONChamberPainter(att,i);
    
    painter->SetResponder("Chamber");
    
    painter->SetOutlined("*",kFALSE);
    
    painter->SetOutlined("MANU",kTRUE);
    
    for ( Int_t j = 0; j < 3; ++j ) 
    {
      painter->SetLine(j,1,4-j);
    }
    
    matrix->Adopt(painter);    
  }
  AliMUONPainterRegistry::Instance()->Register(matrix);
  return matrix;
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
        AliMUONPainterDataRegistry::Instance()->Register(maker);
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
        AliMUONPainterDataRegistry::Instance()->Register(maker);
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
  
  TIter next(fPainterMasterFrameList);
  AliMUONPainterMasterFrame* pmf;
  Bool_t first(kTRUE);
  
  while ( ( pmf = static_cast<AliMUONPainterMasterFrame*>(next()) ) )
  {
    pmf->SaveAs(gSystem->ExpandPathName(Form("%s",fileInfo.fFilename)),
                first ? "RECREATE" : "UPDATE");
    first = kFALSE;
  }
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

  rn->AddLine("1.10");
  rn->AddLine("");
  rn->AddLine("Make the raw OCDB more obvious in the data source tab");
  rn->AddLine("");
  
  rn->AddLine("1.08");
  rn->AddLine("");
  rn->AddLine("Changed the default OCDB to 2011 version");
  rn->AddLine("");
  
  rn->AddLine("1.07");
  rn->AddLine("");
  rn->AddLine("Added the RejectList as a possible OCDB data source");
  rn->AddLine("");
  
  rn->AddLine("1.06");
  rn->AddLine("");
  rn->AddLine("Changed a bit the HV display. Now a trip is indicated with a value of -1");
  rn->AddLine("");
  
  rn->AddLine("1.05");
  rn->AddLine("");
  rn->AddLine("Added the possibility to select an event range when reading raw data");
  rn->AddLine("Usefull e.g. to look at a single suspect event...");
  rn->AddLine("");
  
  rn->AddLine("1.04");
  rn->AddLine("");
  rn->AddLine("Changed the default OCDB to 2010 version");
  rn->AddLine("");
  
  rn->AddLine("1.03");
  rn->AddLine("");
  rn->AddLine("Add Print buttons");
  rn->AddLine("Add the automatic creation of often used canvases when using pedestal source");
  // Internal reorganization to allow several independent tabs to be created to 
  // show different master frames (not used yet). Important for the moment
  // is the ability to create a PainterMatrix and pass it to the PainterMasterFrame
  rn->AddLine("");
  
  rn->AddLine("1.02");
  rn->AddLine("");
  rn->AddLine("Internal change (merging of AliMUONTrackerACFDataMaker and AliMUONTrackerOCDBDataMaker into AliMUONTrackerConditionDataMaker)");
  rn->AddLine("Added --ocdb option");
  rn->AddLine("Corrected the display of the configuration");
  rn->AddLine("Corrected the interpretation of the switches for the HV display");
  rn->AddLine("");
  
  rn->AddLine("1.01");
  rn->AddLine("");
  rn->AddLine("Added the configuration as a possible OCDB data source");
  rn->AddLine("");
  
  rn->AddLine("1.00");
  rn->AddLine("");
  rn->AddLine("Added the Status and StatusMap as a possible OCDB data source");
  rn->AddLine("");
  rn->AddLine("Added one (computed) dimension to the Gains data source = 1/a1/0.2 (mV/fC)");
  rn->AddLine("");
    
  rn->AddLine("0.99a");
  rn->AddLine("");
  rn->AddLine("Added the --de and --chamber options");
  rn->AddLine("");
  
  rn->AddLine("0.99");
  rn->AddLine("");
  rn->AddLine("The chamberid in the label (top right of panel) is now starting at 1 as in common usage");  
  rn->AddLine("");
  
  rn->AddLine("0.98");
  rn->AddLine("");
  rn->AddLine("Added --asciimapping option");
  rn->AddLine("");
  
  rn->AddLine("0.97");
  rn->AddLine("");
  rn->AddLine("Adding calibration option with Emelec (aka injection) gain");
  rn->AddLine("");
  
  rn->AddLine("0.96a");
  rn->AddLine("");
  rn->AddLine("Internal reorganization of the contour computations, that lead to improved performance. ");
  rn->AddLine("Improved enough to be able to remove completely the usage of the padstore.root file with precomputed contours.");
  rn->AddLine("");
  
  rn->AddLine("0.96");
  rn->AddLine("");
  rn->AddLine("New features");
  rn->AddLine("");
  rn->AddLine("- Can now read raw data from memory (using the mem://@gdc: syntax)");
  rn->AddLine("- Raw data decoder now automatically skips buspatches with parity errors");
  rn->AddLine("");
  
  rn->AddLine("0.95");
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
  
  AliMUONPainterDataRegistry* reg = AliMUONPainterDataRegistry::Instance();

  TFile f(filename,"RECREATE");
  
  for ( Int_t i = 0; i < reg->NumberOfDataMakers(); ++i )
  {
    AliMUONVTrackerDataMaker* maker = reg->DataMaker(i);
    maker->Write();
  }
  
  f.Close();
}  
