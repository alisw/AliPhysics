#include "AliMUONMchViewApplication.h"

#include "AliCDBManager.h"
#include "AliCodeTimer.h"
#include "AliLog.h"
#include "AliMUONPainterDataSourceFrame.h"
#include "AliMUONPainterHelper.h"
#include "AliMUONPainterMasterFrame.h"
#include "AliMUONPainterRegistry.h"
#include "AliMUONVTrackerData.h"
#include <Riostream.h>
#include <TCanvas.h>
#include <TEnv.h>
#include <TFile.h>
#include <TGClient.h>
#include <TGFileDialog.h>
#include <TGMenu.h>
#include <TGTab.h>
#include <TKey.h>
#include <TList.h>
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

//______________________________________________________________________________
AliMUONMchViewApplication::AliMUONMchViewApplication(const char* name,
                                                     int* argc, char** argv,
                                                     Float_t wfraction,
                                                     Float_t hfraction) 
: TRint(name,argc,argv),
  fMainFrame(0x0)
{

  /// ctor
  /// wfraction,hfraction are the fractions of display width and height
  /// we want to draw on
  
  UInt_t dw = gClient->GetDisplayWidth(); 
  UInt_t dh = gClient->GetDisplayHeight(); 
                   
  UInt_t w = (UInt_t)(wfraction*dw);
  UInt_t h = (UInt_t)(hfraction*dh);

  fMainFrame = new TGMainFrame(gClient->GetRoot(),w,h);
  
  CreateMenuBar(w);

  const Int_t kbs = 2;
  
//  h -= 60; // menubar
  
  TGTab* tabs = new TGTab(fMainFrame,w,h);
  
  TGCompositeFrame* t = tabs->AddTab("Painter Master Frame");

  AliMUONPainterMasterFrame* pf = 
    new AliMUONPainterMasterFrame(t,t->GetWidth()-kbs*2,t->GetHeight()-kbs*2);
  
  t->AddFrame(pf, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,kbs,kbs,kbs,kbs));

  t = tabs->AddTab("Data Sources");
  
  AliMUONPainterDataSourceFrame* dsf = 
    new AliMUONPainterDataSourceFrame(t,t->GetWidth()-kbs*2,t->GetHeight()-kbs*2);
  
  t->AddFrame(dsf,new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,kbs,kbs,kbs,kbs));
  
  fMainFrame->AddFrame(tabs,new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,0,0,0,0));

  fMainFrame->SetWindowName("mchview - Visualization of MUON Tracker detector");

  fMainFrame->MapSubwindows();
  fMainFrame->Resize();
  fMainFrame->MapWindow();
  
  fMainFrame->Connect("CloseWindow()","AliMUONMchViewApplication",this,"Terminate()");
  
  UInt_t x = dw/2 - w/2;
  UInt_t y = 0;
  
  fMainFrame->MoveResize(x, y, w, h); 
  fMainFrame->SetWMPosition(x, y);
  
  fMainFrame->SetWMSizeHints(w,h,w,h,0,0);
}

//______________________________________________________________________________
AliMUONMchViewApplication::~AliMUONMchViewApplication()
{
  /// dtor
}

//______________________________________________________________________________
void
AliMUONMchViewApplication::CreateMenuBar(UInt_t w)
{
  /// Create the application menu bar
  
  TGPopupMenu* file = new TGPopupMenu(gClient->GetRoot());
  
  file->AddEntry("&Open...",fgkFILEOPEN);
  file->AddEntry("&Save As...",fgkFILESAVEAS);
  file->AddEntry("&Exit",fgkFILEEXIT);
  
  file->Connect("Activated(Int_t)","AliMUONMchViewApplication",this,"HandleMenu(Int_t)");
  
  TGMenuBar* bar = new TGMenuBar(fMainFrame,w);
  
  bar->AddPopup("&File",file,new TGLayoutHints(kLHintsLeft|kLHintsTop));
  
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
  
  new TGFileDialog(gClient->GetRoot(),gClient->GetRoot(),
                   kFDOpen,&fileInfo);
  
  Open(gSystem->ExpandPathName(Form("%s",fileInfo.fFilename)));
}  

//______________________________________________________________________________
void
AliMUONMchViewApplication::Open(const char* filename)
{
  /// Open a given file containing saved VTrackerData objects
  
  TFile f(filename);
  
  TList* keys = f.GetListOfKeys();
  TIter next(keys);
  
  TKey* k;
  
  while ( ( k = static_cast<TKey*>(next()) ) )
  {
    TObject* object = k->ReadObj();
    
    if ( object->InheritsFrom("AliMUONVTrackerData") )
    {
      AliMUONVTrackerData* data = static_cast<AliMUONVTrackerData*>(object);
      if ( data ) 
      {
        AliMUONPainterRegistry::Instance()->Register(data);
      }
    }
  }
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
  /// Save VTrackerData objects into file of given name
  
  AliMUONPainterRegistry* reg = AliMUONPainterRegistry::Instance();

  TFile f(filename,"RECREATE");
  
  for ( Int_t i = 0; i < reg->NumberOfDataSources(); ++i )
  {
    AliMUONVTrackerData* data = reg->DataSource(i);
    data->Write();
  }
  
  f.Close();
}  
