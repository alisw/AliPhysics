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

#include "AliMUONPainterDataSourceItem.h"

#include "AliMUONPainterEnv.h"
#include "AliMUONPainterHelper.h"
#include "AliMUONPainterDataRegistry.h"
#include "AliMUONVTrackerDataMaker.h"
#include "AliMUONVTrackerData.h"
#include "AliLog.h"
#include <TFile.h>
#include <TGFileDialog.h>
#include <TGLabel.h>
#include <TGButton.h>
#include <TSystem.h>
#include <TThread.h>
#include <Riostream.h>

///\class AliMUONPainterDataSourceItem
///
/// Widget to show one data source, and allow to run/stop/rewind/remove it
///
/// WARNING : the thread business is not really working yet (AliRawReaders are
/// not really thread-safe for the moment). So please use a single raw data
/// source at a time, otherwise you'll get a crash !
///
///\author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONPainterDataSourceItem)
///\endcond

namespace
{
  void* RunFunction(void* args)
  {
    Long_t* params = (Long_t*)(args);
    
    AliMUONPainterDataSourceItem* calling = reinterpret_cast<AliMUONPainterDataSourceItem*>(params[0]);
    AliMUONVTrackerDataMaker* reader = reinterpret_cast<AliMUONVTrackerDataMaker*> (params[1]);
    
    Bool_t ok(kTRUE);
    
    while ( ok ) 
    {
      ok = reader->NextEvent();
      if ( reader->IsZombie() ) 
      {
        AliMUONPainterDataRegistry::Instance()->DeleteZombies();
        return 0x0;
      }
      if ( !reader->IsRunning() ) gSystem->Sleep(1000);
    }
    
    calling->Rewind();
    
    return 0x0;
  }
}

//_____________________________________________________________________________
AliMUONPainterDataSourceItem::AliMUONPainterDataSourceItem(const TGWindow* p,
                                                           UInt_t w, UInt_t h,
                                                           AliMUONVTrackerDataMaker* maker)
: TGCompositeFrame(p,w,h,kHorizontalFrame),
  fDataMaker(maker),
  fSourceName(new TGLabel(this,maker->Data()->Name())),
  fSource(new TGLabel(this,maker->Source().Data())),
  fNumberOfEvents(new TGLabel(this,Form("%10d",0))),
  fRun(0x0),
  fStop(0x0),
  fRewind(0x0),
  fRemove(new TGTextButton(this,"Remove")),
  fSave(new TGTextButton(this,"Save")),
  fSaveAs(new TGTextButton(this,"Save As...")),
  fThread(0x0),
  fShouldReset(kFALSE)
{
    /// ctor
 
    SetCleanup(kDeepCleanup);
    
    Update();
    
    AddFrame(fSourceName, new TGLayoutHints(kLHintsNormal | kLHintsCenterY,5,5,5,5));
    AddFrame(fSource,new TGLayoutHints(kLHintsExpandX | kLHintsCenterY,5,5,5,5));
    AddFrame(fNumberOfEvents,new TGLayoutHints(kLHintsNormal | kLHintsCenterY,5,5,5,5));

    if ( fDataMaker->IsRunnable() ) 
    {
      fRun = new TGTextButton(this,"Run");
      fStop = new TGTextButton(this,"Stop");
      fRewind = new TGTextButton(this,"Rewind");
      
      fRun->SetEnabled(!maker->Data()->IsSingleEvent());
      fRun->Connect("Clicked()",
                    "AliMUONPainterDataSourceItem",
                    this,
                    "Run()");
      
      fStop->SetEnabled(kFALSE);
      fStop->Connect("Clicked()",
                     "AliMUONPainterDataSourceItem",
                     this,
                     "Stop()");
      
      fRewind->SetEnabled(kFALSE);
      fRewind->Connect("Clicked()",
                       "AliMUONPainterDataSourceItem",
                       this,
                       "Rewind()");
      
      AddFrame(fRun,new TGLayoutHints(kLHintsCenterY | kLHintsCenterY,5,5,5,5));
      AddFrame(fStop,new TGLayoutHints(kLHintsCenterY | kLHintsCenterY,5,5,5,5));
      AddFrame(fRewind,new TGLayoutHints(kLHintsCenterY | kLHintsCenterY,5,5,5,5));    
    }

    AddFrame(fRemove,new TGLayoutHints(kLHintsCenterY | kLHintsCenterY,5,5,5,5));    

    AddFrame(fSave,new TGLayoutHints(kLHintsCenterY | kLHintsCenterY,5,5,5,5));    

    AddFrame(fSaveAs,new TGLayoutHints(kLHintsCenterY | kLHintsCenterY,5,5,5,5));    

    maker->Data()->Connect("NumberOfEventsChanged()",
                            "AliMUONPainterDataSourceItem",
                            this,
                            "Update()");
    
    fRemove->Connect("Clicked()",
                     "AliMUONPainterDataSourceItem",
                     this,
                     "Remove()");
    
    fSave->Connect("Clicked()",
                   "AliMUONPainterDataSourceItem",
                   this,
                   "Save()");

    fSaveAs->Connect("Clicked()",
                   "AliMUONPainterDataSourceItem",
                   this,
                   "SaveWithDialog()");
    
    Resize();
}

//_____________________________________________________________________________
AliMUONPainterDataSourceItem::~AliMUONPainterDataSourceItem()
{
  /// dtor
  TThread::Delete(fThread);
  delete fThread;
}


//_____________________________________________________________________________
void 
AliMUONPainterDataSourceItem::EnableRun() 
{ 
  /// Enable run button
  if ( fRun ) 
  {
    fRun->SetEnabled(kTRUE); 
  }
}
  
//_____________________________________________________________________________
void 
AliMUONPainterDataSourceItem::DisableRun() 
{ 
  /// Disable run button
  if ( fRun )
  {
    fRun->SetEnabled(kFALSE); 
  }
}

//_____________________________________________________________________________
void
AliMUONPainterDataSourceItem::Remove()
{
  /// Remove
  
  MakeZombie();
  AliMUONPainterDataRegistry::Instance()->Unregister(fDataMaker);
}

//_____________________________________________________________________________
void
AliMUONPainterDataSourceItem::Reset()
{
  /// Reset the data
  fDataMaker->Data()->Clear();
}

//_____________________________________________________________________________
void
AliMUONPainterDataSourceItem::Rewind()
{
  /// Rewind button was clicked
  
  fRewind->SetEnabled(kTRUE);
  
  Stop();
  
  TThread::Delete(fThread);
  delete fThread;
  fThread = 0x0;
  
  if ( fRun && fStop && fRewind ) 
  {
    fRun->SetEnabled(kTRUE);
    fStop->SetEnabled(kFALSE);
    fRewind->SetEnabled(kFALSE);
  }
  
  fDataMaker->Rewind();
  
  fShouldReset = kTRUE;
}

//_____________________________________________________________________________
void
AliMUONPainterDataSourceItem::Run()
{
  /// Run button was clicked
  
  StartRunning();
  
  if ( fShouldReset ) 
  {
    Reset();
    fShouldReset = kFALSE;
  }
  
  fRemove->SetEnabled(kFALSE);
  
  if (!fThread)
  {
    fParams[0] = (Long_t)(this);
    fParams[1] = (Long_t)(fDataMaker);
    fThread = new TThread(RunFunction,(void*)(&fParams[0]));
    fThread->Run();
  }
  
  fDataMaker->SetRunning(kTRUE);
  
  if ( fRun && fStop )
  {
    fRun->SetEnabled(kFALSE);
    fStop->SetEnabled(kTRUE);
  }
}

//_____________________________________________________________________________
void
AliMUONPainterDataSourceItem::Save(const char* filename)
{
  /// Save the data maker
  
  TFile* f = TFile::Open(filename,"RECREATE");
  
  fDataMaker->Write();
  
  f->Write();
  f->Close();
  
  delete f;
}

//_____________________________________________________________________________
void
AliMUONPainterDataSourceItem::Save()
{
  /// Save the data maker (filename is fixed)
  
  TString dname(fDataMaker->Data()->GetName());
  dname.ToLower();
  
  TString outputDir(AliMUONPainterHelper::Instance()->Env()->String("LastSaveDir","."));

  TString filename(Form("%s/mchview.%s.root",gSystem->ExpandPathName(outputDir.Data()),dname.Data()));
  
  Save(filename.Data());
}

//_____________________________________________________________________________
void
AliMUONPainterDataSourceItem::SaveWithDialog()
{
  /// Save the data maker (filename given by dialog)
  
  TGFileInfo fileInfo;
  
//  fileInfo.fFileTypes = fgkFileTypes;
  
  delete[] fileInfo.fIniDir;
  
  AliMUONPainterEnv* env = AliMUONPainterHelper::Instance()->Env();
  
  fileInfo.fIniDir = StrDup(env->String("LastSaveDir","."));
  
  new TGFileDialog(gClient->GetRoot(),gClient->GetRoot(),
                   kFDSave,&fileInfo);
  
  env->Set("LastSaveDir",fileInfo.fIniDir);
  env->Save();  
  
  Save(fileInfo.fFilename);  
}

//_____________________________________________________________________________
void
AliMUONPainterDataSourceItem::Stop()
{
  /// Stop button was clicked
  
  StopRunning();
  
  fDataMaker->SetRunning(kFALSE);
  
  if ( fStop && fRun ) 
  {
    fStop->SetEnabled(kFALSE);
    fRun->SetEnabled(kTRUE);
  }
  
  fRemove->SetEnabled(kTRUE);
}

//_____________________________________________________________________________
void
AliMUONPainterDataSourceItem::Update()
{
  /// Update ourselves
  
  fNumberOfEvents->SetText(Form("%10d",fDataMaker->Data()->NumberOfEvents(-1)));
}

//_____________________________________________________________________________
void
AliMUONPainterDataSourceItem::StartRunning()
{
  /// Signal we start to run
  Emit("StartRunning()");
}  

//_____________________________________________________________________________
void
AliMUONPainterDataSourceItem::StopRunning()
{
  /// Signal we stop to run
  Emit("StopRunning()");
}
