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

#include "AliMUONPainterDataSourceFrame.h"

#include "AliMUONPainterDataSourceItem.h"
#include "AliMUONPainterEnv.h"
#include "AliMUONPainterHelper.h"
#include "AliMUONPainterRegistry.h"
#include "AliMUONTrackerOCDBDataMaker.h"
#include "AliMUONTrackerRawDataMaker.h"
#include "AliMUONVTrackerDataMaker.h"
#include "AliLog.h"
#include "AliRawReaderDate.h"
#include "AliRawReaderRoot.h"
#include <TGButton.h>
#include <TGComboBox.h>
#include <TGFileDialog.h>
#include <TGNumberEntry.h>
#include <TGTextEntry.h>
#include <TGrid.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TRegexp.h>
#include <TString.h>
#include <TSystem.h>

///\class AliMUONPainterDataSourceFrame
///
/// A complete frame to select and display various data sources to 
/// be displayed : either raw data or OCDB data. 
/// Later on we might add digits and clusters for instance.
///
///\author Laurent Aphecetche, Subatech

const char* AliMUONPainterDataSourceFrame::fgkNumberOfDataSourcesKey = "NumberOfDataSources";
const char* AliMUONPainterDataSourceFrame::fgkDataSourceURIKey = "DataSourceURI.%d";

///\cond CLASSIMP
ClassImp(AliMUONPainterDataSourceFrame)
///\endcond

//_____________________________________________________________________________
AliMUONPainterDataSourceFrame::AliMUONPainterDataSourceFrame(const TGWindow* p, UInt_t w, UInt_t h)
: TGCompositeFrame(p,w,h,kVerticalFrame),
  fRecentSourceSelector(new TGGroupFrame(p,"Recent sources",kHorizontalFrame)),
  fRawSelector(new TGGroupFrame(p,"Raw file URI",kHorizontalFrame)),
  fOCDBSelector(new TGGroupFrame(p,"OCDB Path",kHorizontalFrame)),
  fDataReaders(new TGGroupFrame(p,"Data sources")),
  fFilePath(new TGTextEntry(fRawSelector,"")),
  fRawOCDBPath(new TGTextEntry(fRawSelector,"")),
  fOCDBPath(new TGTextEntry(fOCDBSelector,"")),
  fRunSelector(new TGNumberEntry(fOCDBSelector,0)),
  fOCDBTypes(new TGComboBox(fOCDBSelector)),
  fRecentSources(new TGComboBox(fRecentSourceSelector)),
  fItems(new TObjArray)
{
  /// Ctor
  
    AliMUONPainterRegistry* reg = AliMUONPainterRegistry::Instance();
    
    reg->Connect("DataReaderWasRegistered(AliMUONVTrackerDataMaker*)",
                 "AliMUONPainterDataSourceFrame",
                 this,
                 "DataReaderWasRegistered(AliMUONVTrackerDataMaker*)");
    
    reg->Connect("DataReaderWasUnregistered(AliMUONVTrackerDataMaker*)",
                 "AliMUONPainterDataSourceFrame",
                 this,
                 "DataReaderWasUnregistered(AliMUONVTrackerDataMaker*)");
    
    fItems->SetOwner(kFALSE);
    
    /// Recent source selection
    
    AliMUONPainterEnv* env = AliMUONPainterHelper::Instance()->Env();
    
    Int_t nsources = env->Integer(fgkNumberOfDataSourcesKey);
    
    for ( Int_t i = 0; i < nsources; ++i )
    {
      AddRecentSource(env->String(Form(fgkDataSourceURIKey,i)));
    }

    fRecentSources->Resize(100,20);
    
    TGButton* createRecentButton = new TGTextButton(fRecentSourceSelector,"Create data source");
    createRecentButton->Connect("Clicked()",
                                "AliMUONPainterDataSourceFrame",
                                this,
                                "OpenRecentSource()");
    
    fRecentSourceSelector->AddFrame(fRecentSources,new TGLayoutHints(kLHintsExpandX | kLHintsTop,5,5,5,5));
    fRecentSourceSelector->AddFrame(createRecentButton,new TGLayoutHints(kLHintsTop,5,5,5,5));
                                    
    /// Raw file selection
    
    TGButton* openButton = new TGPictureButton(fRawSelector,
                                           gClient->GetPicture("fileopen.xpm"));
    openButton->SetToolTipText("Click to open file dialog");
                                        
    TGButton* createRawButton = new TGTextButton(fRawSelector,"Create data source");
    
    fRawSelector->AddFrame(fFilePath, new TGLayoutHints(kLHintsExpandX | kLHintsTop,5,5,5,5));
    fRawSelector->AddFrame(fRawOCDBPath, new TGLayoutHints(kLHintsTop,5,5,5,5));
    fRawSelector->AddFrame(openButton,new TGLayoutHints(kLHintsTop,5,5,5,5));
    fRawSelector->AddFrame(createRawButton,new TGLayoutHints(kLHintsTop,5,5,5,5));
    
    openButton->Connect("Clicked()",
                        "AliMUONPainterDataSourceFrame",
                        this,
                        "OpenFileDialog()");

    createRawButton->Connect("Clicked()",
                        "AliMUONPainterDataSourceFrame",
                        this,
                        "CreateRawDataSource()");
    
    /// OCDB selection
    
    fOCDBTypes->AddEntry("Pedestals",0);
    fOCDBTypes->AddEntry("Gains",1);
    fOCDBTypes->AddEntry("Capacitances",2);
    fOCDBTypes->Select(0);
    fOCDBTypes->Resize(100,20);
    
    TGButton* createOCDBButton = new TGTextButton(fOCDBSelector,"Create data source");
    createOCDBButton->Connect("Clicked()",
                             "AliMUONPainterDataSourceFrame",
                             this,
                             "CreateOCDBDataSource()");
    
    
    fOCDBSelector->AddFrame(fOCDBPath,new TGLayoutHints(kLHintsExpandX | kLHintsTop,5,5,5,5));    
    fOCDBSelector->AddFrame(fRunSelector,new TGLayoutHints(kLHintsTop,5,5,5,5));
    fOCDBSelector->AddFrame(fOCDBTypes,new TGLayoutHints(kLHintsExpandX | kLHintsTop,5,5,5,5));
    fOCDBSelector->AddFrame(createOCDBButton,new TGLayoutHints(kLHintsTop,5,5,5,5));

    AddFrame(fRecentSourceSelector,new TGLayoutHints(kLHintsExpandX,10,10,10,10));

    AddFrame(fRawSelector,new TGLayoutHints(kLHintsExpandX,10,10,10,10));

    AddFrame(fOCDBSelector,new TGLayoutHints(kLHintsExpandX,10,10,10,10));

    AddFrame(fDataReaders, new TGLayoutHints(kLHintsExpandX,10,10,10,10));
    
}

//_____________________________________________________________________________
AliMUONPainterDataSourceFrame::~AliMUONPainterDataSourceFrame()
{
  /// dtor
  
  delete fItems;
}

//_____________________________________________________________________________
void
AliMUONPainterDataSourceFrame::AddRecentSource(const char* name)
{  
  /// Add a source to the list of recently used sources
  
  TGListBox* lb = fRecentSources->GetListBox();
  
  for ( Int_t i = 0; i < lb->GetNumberOfEntries(); ++i ) 
  {
    TGTextLBEntry* t = (TGTextLBEntry*)lb->GetEntry(i);
    TString s(t->GetText()->GetString());
    if ( s == name ) 
    {
      return;
    }
  }
  
  fRecentSources->AddEntry(name,lb->GetNumberOfEntries());
  fRecentSources->MapSubwindows();
  fRecentSources->Layout();
}

//_____________________________________________________________________________
void
AliMUONPainterDataSourceFrame::CreateOCDBDataSource()
{
  /// Create an OCDB data source (using information from the widgets)
  
  TString cdbPath = fOCDBPath->GetText();
  Int_t runNumber = fRunSelector->GetIntNumber();
  TGTextLBEntry* t = static_cast<TGTextLBEntry*>(fOCDBTypes->GetSelectedEntry());
  TString type = t->GetText()->GetString();
  
  CreateOCDBDataSource(cdbPath,runNumber,type);
  
  fOCDBPath->SetText("");
  fRunSelector->SetNumber(0);  
}

//_____________________________________________________________________________
void
AliMUONPainterDataSourceFrame::CreateOCDBDataSource(const TString& uri)
{
  /// Create an OCDB data source, given it's URI
  
  TObjArray* a = uri.Tokenize(";");
  TString cdbPath = static_cast<TObjString*>(a->At(1))->String();
  TString srun = static_cast<TObjString*>(a->At(2))->String();
  TString type = static_cast<TObjString*>(a->At(3))->String();
  
  CreateOCDBDataSource(cdbPath,atoi(srun.Data()),type);
  
  delete a;
}

//_____________________________________________________________________________
void
AliMUONPainterDataSourceFrame::CreateOCDBDataSource(const TString& cdbPath,
                                                    Int_t runNumber,
                                                    const TString& type)
{
  /// Create an OCDB data source for a given (path,runnumber,type) triplet
  
  AliMUONVTrackerDataMaker* reader = new AliMUONTrackerOCDBDataMaker(cdbPath.Data(),
                                                                       runNumber,
                                                                       type.Data());
  
  if ( reader->IsValid() ) 
  {
    AliMUONPainterRegistry::Instance()->Register(reader);
    
    AliMUONPainterEnv* env = AliMUONPainterHelper::Instance()->Env();
    
    Int_t n = env->Integer(fgkNumberOfDataSourcesKey);
    
    env->Set(fgkNumberOfDataSourcesKey,n+1);
    
    TString ds(Form("OCDB;%s;%d;%s",cdbPath.Data(),runNumber,type.Data()));
    
    env->Set(Form(fgkDataSourceURIKey,n),ds.Data());
    
    env->Save();
    
    AddRecentSource(ds.Data());
  }
}

//_____________________________________________________________________________
void 
AliMUONPainterDataSourceFrame::CreateRawDataSource()
{
  /// Create a new raw data source (using info from the widgets)
  
  TString uri(gSystem->ExpandPathName(fFilePath->GetText()));
  
  if ( gSystem->AccessPathName(uri.Data()) )
  {
    AliError(Form("File %s does not exist",uri.Data()));
    fFilePath->SetText("");
    return;
  }

  uri = Form("RAW;%s;%s",uri.Data(),fRawOCDBPath->GetText());
  
  if ( CreateRawDataSource(uri) )
  {
    fFilePath->SetText("");
    fRawOCDBPath->SetText("");
  }
}

//_____________________________________________________________________________
Bool_t 
AliMUONPainterDataSourceFrame::CreateRawDataSource(const TString& uri)
{
  /// Create a new raw data source, given its URI
  
  TString filename;
  TString ocdbPath;
  
  TObjArray* a = uri.Tokenize(";");
  
  filename = static_cast<TObjString*>(a->At(1))->String();
  
  if ( a->GetLast() > 1 ) 
  {
    ocdbPath = static_cast<TObjString*>(a->At(2))->String();
  }
  
  AliRawReader* rawReader = 0x0;

  if ( filename.Contains(TRegexp(".root$")) ) 
  {
    AliDebug(1,"Using RawReaderRoot");
    if ( filename.Contains(TRegexp("^alien")) )
    {
      // insure we've initialized the grid...
      if (!gGrid)
      {
        TGrid::Connect("alien://");
      }
    }
         
    rawReader = new AliRawReaderRoot(filename.Data());
  }
  else if ( uri.Contains(TRegexp(".raw")) )
  {
    AliDebug(1,"Using RawReaderDate");
    rawReader = new AliRawReaderDate(filename.Data());
  }
  else
  {
    AliError(Form("Don't know how to open that file : %s\nURI=%s",filename.Data(),uri.Data()));
    return kFALSE;
  }
  
  AliMUONTrackerRawDataMaker* reader = new AliMUONTrackerRawDataMaker(rawReader,ocdbPath.Data());
  
  reader->SetSource(filename.Data());
  
  AliMUONPainterRegistry::Instance()->Register(reader);
  
  AliMUONPainterEnv* env = AliMUONPainterHelper::Instance()->Env();
  
  Int_t n = env->Integer(fgkNumberOfDataSourcesKey);
  
  env->Set(fgkNumberOfDataSourcesKey,n+1);
  
  TString ds(Form("RAW;%s;%s",filename.Data(),ocdbPath.Data()));
  
  env->Set(Form(fgkDataSourceURIKey,n),ds.Data());
  
  AddRecentSource(ds.Data());
  
  env->Save();

  return kTRUE;
}

//_____________________________________________________________________________
void 
AliMUONPainterDataSourceFrame::DataReaderWasRegistered(AliMUONVTrackerDataMaker* reader)
{
  /// Update ourselves as a new data reader was created
  
  AliInfo(Form("%s",reader->GetName()));
  
  AliMUONPainterDataSourceItem* item = new AliMUONPainterDataSourceItem(fDataReaders,100,20,reader);
      
  fDataReaders->AddFrame(item);
  
  fItems->Add(item);

  fDataReaders->MapSubwindows();
  fDataReaders->Resize();
}

//_____________________________________________________________________________
void 
AliMUONPainterDataSourceFrame::DataReaderWasUnregistered(AliMUONVTrackerDataMaker* reader)
{
  /// Update ourselves as a new data reader was deleted
  
  AliInfo(Form("%s",reader->GetName()));
          
}

//_____________________________________________________________________________
void
AliMUONPainterDataSourceFrame::OpenFileDialog()
{
  /// Open a file dialog to select a file to be read
  
  TGFileInfo fileInfo;
  
  new TGFileDialog(gClient->GetRoot(),gClient->GetRoot(),
                   kFDOpen,&fileInfo);
  
  
  fFilePath->SetText(gSystem->ExpandPathName(Form("%s",fileInfo.fFilename)));
}


//_____________________________________________________________________________
void
AliMUONPainterDataSourceFrame::OpenRecentSource()
{
  /// Open one source from the recently used ones
  
  TGTextLBEntry* t = (TGTextLBEntry*)fRecentSources->GetSelectedEntry();

  TString uri(t->GetText()->GetString());
  
  if ( uri.Contains(TRegexp("^RAW")) )
  {
    CreateRawDataSource(uri);
  }
  else if ( uri.Contains(TRegexp("^OCDB")) )
  {
    CreateOCDBDataSource(uri);
  }
  
  fRecentSources->Select(-1);
}

