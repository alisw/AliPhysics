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

#include "AliMUONPainterRegistry.h"

#include "AliMUONPainterMatrix.h"
#include "AliMUONVTrackerData.h"
#include "AliMUONVTrackerDataMaker.h"
#include "AliLog.h"
#include <TGMenu.h>
#include <TGWindow.h>
#include <THashList.h>
#include <TObjArray.h>
#include <TString.h>
#include <Riostream.h>

///\class AliMUONPainterRegistry
///
/// Registry for AliMUONVPainter related stuff : painter data sources
/// and painter matrices
///
///\author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONPainterRegistry)
///\endcond

AliMUONPainterRegistry* AliMUONPainterRegistry::fgInstance(0x0);

//_____________________________________________________________________________
AliMUONPainterRegistry::AliMUONPainterRegistry() : TObject(), TQObject(),
fPainterMatrices(new TObjArray),
fDataMakers(new TObjArray),
fHistoryMenu(0x0),
fMenuBar(0x0),
fHistoryCounter(0),
fZombies(new TObjArray)
{
  /// ctor
  fPainterMatrices->SetOwner(kTRUE);
  fDataMakers->SetOwner(kTRUE);
  fZombies->SetOwner(kTRUE);
}

//_____________________________________________________________________________
AliMUONPainterRegistry::~AliMUONPainterRegistry()
{
  /// dtor
  delete fPainterMatrices;
  delete fDataMakers;
}

//_____________________________________________________________________________
AliMUONVTrackerDataMaker* 
AliMUONPainterRegistry::DataMaker(Int_t i) const
{
  /// Get one data source
  if ( i >= 0 && i <= fDataMakers->GetLast() )
  {
    return static_cast<AliMUONVTrackerDataMaker*>(fDataMakers->At(i));
  }
  else
  {
    AliError(Form("Index out of bounds : %d / %d",i,fDataMakers->GetLast()+1));
    return 0x0;
  }
}

//_____________________________________________________________________________
AliMUONVTrackerData* 
AliMUONPainterRegistry::DataSource(Int_t i) const
{
  /// Get one data source
  
  AliMUONVTrackerDataMaker* maker = DataMaker(i);
  if ( maker ) return maker->Data();
  return 0x0;
}

//_____________________________________________________________________________
void 
AliMUONPainterRegistry::DataMakerWasRegistered(AliMUONVTrackerDataMaker* data)
{
  /// A new reader source was registered
  Long_t param[] = { (Long_t)data };
  
  Emit("DataMakerWasRegistered(AliMUONVTrackerDataMaker*)",param);
}

//_____________________________________________________________________________
void
AliMUONPainterRegistry::DataMakerWasUnregistered(AliMUONVTrackerDataMaker* data)
{
  /// A data reader was unregistered
  Long_t param[] = { (Long_t)data };
  
  Emit("DataMakerWasUnregistered(AliMUONVTrackerDataMaker*)",param);
  
}

//_____________________________________________________________________________
void 
AliMUONPainterRegistry::DataSourceWasRegistered(AliMUONVTrackerData* data)
{
  /// A new data source was registered
  Long_t param[] = { (Long_t)data };
  
  Emit("DataSourceWasRegistered(AliMUONVTrackerData*)",param);
}

//_____________________________________________________________________________
void
AliMUONPainterRegistry::DataSourceWasUnregistered(AliMUONVTrackerData* data)
{
  /// A data source was unregistered
  Long_t param[] = { (Long_t)data };
  
  Emit("DataSourceWasUnregistered(AliMUONVTrackerData*)",param);
  
}

//_____________________________________________________________________________
AliMUONVTrackerData*
AliMUONPainterRegistry::DataSource(const char* name) const
{
  /// Find a data source by name
  for ( Int_t i = 0; i < NumberOfDataMakers(); ++i )
  {
    AliMUONVTrackerData* data = DataMaker(i)->Data();
    if ( data ) 
    {
      TString dname(data->GetName());
      if ( dname == name ) return data;
    }
  }
  return 0x0;
}

//_____________________________________________________________________________
Int_t 
AliMUONPainterRegistry::FindIndexOf(AliMUONPainterMatrix* group) const
{
  /// Get the index of a given painterMatrix
  return fPainterMatrices->IndexOf(group);
}

//_____________________________________________________________________________
AliMUONPainterMatrix*
AliMUONPainterRegistry::PainterMatrix(const char* name) const
{
  /// Get a painterMatrix by name
  return static_cast<AliMUONPainterMatrix*>(fPainterMatrices->FindObject(name));
}

//_____________________________________________________________________________
void
AliMUONPainterRegistry::HistoryMenuActivated(Int_t i)
{
  /// A painterMatrix was chosen from the history menu
  
  AliDebug(1,Form("i=%d",i));
  
  TGMenuEntry* entry = fHistoryMenu->GetEntry(i);
  
  AliMUONPainterMatrix* group = reinterpret_cast<AliMUONPainterMatrix*>(entry->GetUserData());
  
  PainterMatrixWantToShow(group);
}

//_____________________________________________________________________________
AliMUONPainterRegistry*
AliMUONPainterRegistry::Instance()
{
  /// Get unique instance of this class
  if ( !fgInstance ) fgInstance = new AliMUONPainterRegistry;
  return fgInstance;
}

//_____________________________________________________________________________
AliMUONPainterMatrix* 
AliMUONPainterRegistry::PainterMatrix(Int_t i) const
{
  /// Get one painter matrix
  if ( i >= 0 && i <= fPainterMatrices->GetLast() )
  {
    return static_cast<AliMUONPainterMatrix*>(fPainterMatrices->At(i));
  }
  else
  {
    AliError(Form("Index out of bounds : %d / %d",i,fPainterMatrices->GetLast()+1));
    return 0x0;
  }
}

//_____________________________________________________________________________
void 
AliMUONPainterRegistry::PainterMatrixWantToShow(AliMUONPainterMatrix* group)
{
  /// A given paintermatrix want to appear on screen
  Long_t param[] = { (Long_t)group };
  
  Emit("PainterMatrixWantToShow(AliMUONPainterMatrix*)",param);  
}

//_____________________________________________________________________________
void
AliMUONPainterRegistry::AddToHistory(AliMUONPainterMatrix* group)
{
  /// Add a matrix to the history
  
  if ( !fHistoryMenu && fMenuBar ) 
  {
    fHistoryMenu =  new TGPopupMenu(gClient->GetRoot());
    TGPopupMenu* before = 0x0; //FIXME: could try to find a place where to put it (e.g. before Help ?)
    
    fMenuBar->AddPopup("&History",fHistoryMenu, new TGLayoutHints(kLHintsNormal),before);
    
    fHistoryMenu->Connect("Activated(Int_t)",
                          "AliMUONPainterRegistry",this,
                          "HistoryMenuActivated(Int_t)");
    
    AliDebug(1,Form("HistoryMenu create at %x",fHistoryMenu));
  }
  
  if ( fHistoryMenu ) 
  {
    TIter next(fHistoryMenu->GetListOfEntries());
    TGMenuEntry* e(0x0);
    
    while ( ( e = static_cast<TGMenuEntry*>(next()) ) )
    {
      if ( e->GetUserData() == group ) 
      {
        fHistoryMenu->DeleteEntry(e);
        break;
      }
    }
    
    e = static_cast<TGMenuEntry*>(fHistoryMenu->GetListOfEntries()->First());
    
    fHistoryMenu->AddEntry(group->GetName(),++fHistoryCounter,(void*)group,0x0,e);
  }
  else
  {
    AliError("fHistoryMenu is null. We probably did not find the relevant menu entry ?");
  }
}

//_____________________________________________________________________________
void 
AliMUONPainterRegistry::PainterMatrixWasRegistered(AliMUONPainterMatrix* group)
{
  /// A new painter matrix was registered
  Long_t param[] = { (Long_t)group };
  
  Emit("PainterMatrixWasRegistered(AliMUONPainterMatrix*)",param);
}

//_____________________________________________________________________________
void 
AliMUONPainterRegistry::PainterMatrixWasUnregistered(AliMUONPainterMatrix* group)
{
  /// A painter matrix was unregistered
  Long_t param[] = { (Long_t)group };
  
  Emit("PainterMatrixWasUnregistered(AliMUONPainterMatrix*)",param);
}

//_____________________________________________________________________________
void 
AliMUONPainterRegistry::Print(Option_t* opt) const
{
  /// Printout
  TString sopt(opt);
  sopt.ToUpper();
  
  cout << "Number of data readers = " << NumberOfDataMakers() << endl;
  cout << "Number of painter matrices = " << NumberOfPainterMatrices() << endl;
  
  if ( sopt.Contains("FULL") || sopt.Contains("READER") )
  {
    TIter next(fDataMakers);
    AliMUONVTrackerDataMaker* reader;
    
    while ( ( reader = static_cast<AliMUONVTrackerDataMaker*>(next()) ) )
    {
      reader->Print();
    }
  }
  
  if ( sopt.Contains("FULL") || sopt.Contains("DATA") )
  {
    TIter next(fDataMakers);
    AliMUONVTrackerDataMaker* reader;
    
    while ( ( reader = static_cast<AliMUONVTrackerDataMaker*>(next()) ) )
    {
      AliMUONVTrackerData* data = reader->Data();
      if ( data ) data->Print();
    }
  }

  if ( sopt.Contains("FULL") || sopt.Contains("MATRIX") )
  {
    TIter next(fPainterMatrices);
    AliMUONPainterMatrix* matrix;
    
    while ( ( matrix = static_cast<AliMUONPainterMatrix*>(next()) ) )
    {
      matrix->Print();
    }
  }
  
}

//_____________________________________________________________________________
Int_t
AliMUONPainterRegistry::Register(AliMUONPainterMatrix* group)
{
  /// group is adopted, i.e. the registry becomes the owner of it.
  fPainterMatrices->AddLast(group);
  
  PainterMatrixWasRegistered(group);
  
  return fPainterMatrices->IndexOf(group);
}

//_____________________________________________________________________________
void
AliMUONPainterRegistry::Register(AliMUONVTrackerDataMaker* reader)
{
  /// reader is adopted, i.e. the registry becomes the owner of it.
  fDataMakers->AddLast(reader);
  DataMakerWasRegistered(reader);
  if ( reader->Data() ) DataSourceWasRegistered(reader->Data());
}

//_____________________________________________________________________________
Int_t 
AliMUONPainterRegistry::NumberOfDataMakers() const
{
  /// The number of data readers we handle
  return fDataMakers->GetLast()+1;
}

//_____________________________________________________________________________
Int_t 
AliMUONPainterRegistry::NumberOfPainterMatrices() const
{
  /// The number of painter matrices we handle
  return fPainterMatrices->GetLast()+1;
}

//_____________________________________________________________________________
void 
AliMUONPainterRegistry::DeleteZombies()
{
  /// Delete zombies
  fZombies->Delete();
}

//_____________________________________________________________________________
Bool_t 
AliMUONPainterRegistry::Unregister(AliMUONVTrackerDataMaker* reader)
{
  /// Unregister some reader
  
  if (!reader) return kFALSE;
  
  if ( reader->Data() ) 
  {
    DataSourceWasUnregistered(reader->Data());
    reader->Data()->Destroyed(); // we pretend it's deleted now, even
    // if it will be only later on when zombie are killed, so that
    // for instance painters depending on it will no longer try to access it
  }

  DataMakerWasUnregistered(reader);
  
  TObject* o = fDataMakers->Remove(reader);
  
  fZombies->Add(o); // for later deletion
  
//  if ( o ) 
//  {
//    delete o;
//  }
//  else
//  {
//    AliError(Form("Could not unregister data named %s title %s",reader->GetName(),
//                  reader->GetTitle()));
//  }
  return ( o != 0x0 );
}

//_____________________________________________________________________________
Bool_t 
AliMUONPainterRegistry::Unregister(AliMUONPainterMatrix* group)
{
  /// Unregister some matrix
  
  if (!group) return kFALSE;
  
  PainterMatrixWasUnregistered(group);
  
  TObject* o = fPainterMatrices->Remove(group);
  if ( o ) 
  {
    delete o;
  }
  else
  {
    AliError(Form("Could not unregister group named %s",group->GetName()));
  }
  return ( o != 0x0 );
}
