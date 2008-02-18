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

#include "AliMUONVTrackerDataMaker.h"
#include "AliMUONVTrackerData.h"
#include "AliLog.h"
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
      if ( !reader->IsRunning() ) gSystem->Sleep(1000);
    }
    
    calling->Rewind();
    
    return 0x0;
  }
}

//_____________________________________________________________________________
AliMUONPainterDataSourceItem::AliMUONPainterDataSourceItem(const TGWindow* p,
                                                           UInt_t w, UInt_t h,
                                                           AliMUONVTrackerDataMaker* reader)
: TGCompositeFrame(p,w,h,kHorizontalFrame),
  fDataReader(reader),
  fSourceName(new TGLabel(this,reader->Data()->Name())),
  fSource(new TGLabel(this,reader->Source().Data())),
  fNumberOfEvents(new TGLabel(this,Form("%10d",0))),
  fRun(new TGTextButton(this,"Run")),
  fStop(new TGTextButton(this,"Stop")),
  fRewind(new TGTextButton(this,"Rewind")),
  fRemove(0x0),//new TGTextButton(this,"Remove")),
  fThread(0x0),
  fShouldReset(kFALSE)
{
    /// ctor
    
    Update();
    
    fRun->SetEnabled(reader->Data()->IsRunnable());
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
    
//    fRemove->Connect("Clicked()",
//                     "AliMUONPainterDataSourceItem",
//                     this,
//                     "Remove()");
    
    AddFrame(fSourceName, new TGLayoutHints(kLHintsNormal | kLHintsCenterY,5,5,5,5));
    AddFrame(fSource,new TGLayoutHints(kLHintsExpandX | kLHintsCenterY,5,5,5,5));
    AddFrame(fNumberOfEvents,new TGLayoutHints(kLHintsNormal | kLHintsCenterY,5,5,5,5));
    AddFrame(fRun,new TGLayoutHints(kLHintsCenterY | kLHintsCenterY,5,5,5,5));
    AddFrame(fStop,new TGLayoutHints(kLHintsCenterY | kLHintsCenterY,5,5,5,5));
    AddFrame(fRewind,new TGLayoutHints(kLHintsCenterY | kLHintsCenterY,5,5,5,5));    
//    AddFrame(fRemove,new TGLayoutHints(kLHintsCenterY | kLHintsCenterY,5,5,5,5));    
    
    reader->Data()->Connect("NumberOfEventsChanged()",
                            "AliMUONPainterDataSourceItem",
                            this,
                            "Update()");
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
  fRun->SetEnabled(kTRUE); 
}
  
//_____________________________________________________________________________
void 
AliMUONPainterDataSourceItem::DisableRun() 
{ 
  /// Disable run button
  fRun->SetEnabled(kFALSE); 
}
  
//_____________________________________________________________________________
void
AliMUONPainterDataSourceItem::Reset()
{
  /// Reset the data
  fDataReader->Data()->Clear();
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
  
  fRun->SetEnabled(kTRUE);
  fStop->SetEnabled(kFALSE);
  fRewind->SetEnabled(kFALSE);
  
  fDataReader->Rewind();
  
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
  
//  fRemove->SetEnabled(kFALSE);
  
  if (!fThread)
  {
    fParams[0] = (Long_t)(this);
    fParams[1] = (Long_t)(fDataReader);
    fThread = new TThread(RunFunction,(void*)(&fParams[0]));
    fThread->Run();
  }
  
  fDataReader->SetRunning(kTRUE);
  
  fRun->SetEnabled(kFALSE);
  fStop->SetEnabled(kTRUE);
}

//_____________________________________________________________________________
void
AliMUONPainterDataSourceItem::Stop()
{
  /// Stop button was clicked
  
  StopRunning();
  
  fDataReader->SetRunning(kFALSE);
  
  fStop->SetEnabled(kFALSE);
  fRun->SetEnabled(kTRUE);

  //  fRemove->SetEnabled(kTRUE);
}

//_____________________________________________________________________________
void
AliMUONPainterDataSourceItem::Update()
{
  /// Update ourselves
  
  fNumberOfEvents->SetText(Form("%10d",fDataReader->Data()->NumberOfEvents()));
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
