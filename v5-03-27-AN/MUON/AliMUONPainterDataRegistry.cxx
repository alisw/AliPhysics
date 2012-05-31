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

// $Id: AliMUONPainterDataRegistry.cxx 26812 2008-06-20 15:22:59Z laphecet $

#include "AliMUONPainterDataRegistry.h"

#include "AliMpManuIterator.h"
#include "AliMUON2DMap.h"
#include "AliMUONCalibParamND.h"
#include "AliMUONTrackerData.h"
#include "AliMUONVTrackerDataMaker.h"
#include "AliLog.h"
#include <THashList.h>
#include <TObjArray.h>
#include <TString.h>
#include <Riostream.h>

///\class AliMUONPainterDataRegistry
///
/// Registry for AliMUONVPainter related stuff : painter data sources
/// and painter matrices
///
///\author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONPainterDataRegistry)
///\endcond

AliMUONPainterDataRegistry* AliMUONPainterDataRegistry::fgInstance(0x0);

//_____________________________________________________________________________
AliMUONPainterDataRegistry::AliMUONPainterDataRegistry() : TObject(), TQObject(),
fDataMakers(new TObjArray),
fZombies(new TObjArray),
fInteractiveReadOutConfig(0x0)
{
  /// ctor
  fDataMakers->SetOwner(kTRUE);
  fZombies->SetOwner(kTRUE);
}

//_____________________________________________________________________________
AliMUONPainterDataRegistry::~AliMUONPainterDataRegistry()
{
  /// dtor
  delete fDataMakers;
  delete fInteractiveReadOutConfig;
}

//_____________________________________________________________________________
void
AliMUONPainterDataRegistry::CreateInteractiveReadOutConfig() const
{
  /// Create a base config
  
  fInteractiveReadOutConfig = new AliMUONTrackerData("IROC","IROC",1);
  fInteractiveReadOutConfig->SetDimensionName(0,"Switch");
  fInteractiveReadOutConfig->DisableChannelLevel();
  AliMpManuIterator it;
  Int_t detElemId;
  Int_t manuId;
  AliMUON2DMap store(true);

  while ( it.Next(detElemId,manuId) )
  {
    AliMUONVCalibParam* param = new AliMUONCalibParamND(1,64,detElemId,manuId,1);
    store.Add(param);
  }
  fInteractiveReadOutConfig->Add(store);
}

//_____________________________________________________________________________
AliMUONVTrackerDataMaker* 
AliMUONPainterDataRegistry::DataMaker(Int_t i) const
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
AliMUONPainterDataRegistry::DataSource(Int_t i) const
{
  /// Get one data source
  
  AliMUONVTrackerDataMaker* maker = DataMaker(i);
  if ( maker ) return maker->Data();
  return 0x0;
}

//_____________________________________________________________________________
void 
AliMUONPainterDataRegistry::DataMakerWasRegistered(const AliMUONVTrackerDataMaker* data)
{
  /// A new reader source was registered
  Long_t param[] = { (Long_t)data };
  
  Emit("DataMakerWasRegistered(AliMUONVTrackerDataMaker*)",param);
}

//_____________________________________________________________________________
void
AliMUONPainterDataRegistry::DataMakerWasUnregistered(const AliMUONVTrackerDataMaker* data)
{
  /// A data reader was unregistered
  Long_t param[] = { (Long_t)data };
  
  Emit("DataMakerWasUnregistered(AliMUONVTrackerDataMaker*)",param);
  
}

//_____________________________________________________________________________
void 
AliMUONPainterDataRegistry::DataSourceWasRegistered(const AliMUONVTrackerData* data)
{
  /// A new data source was registered
  Long_t param[] = { (Long_t)data };
  
  Emit("DataSourceWasRegistered(AliMUONVTrackerData*)",param);
}

//_____________________________________________________________________________
void
AliMUONPainterDataRegistry::DataSourceWasUnregistered(const AliMUONVTrackerData* data)
{
  /// A data source was unregistered
  Long_t param[] = { (Long_t)data };
  
  Emit("DataSourceWasUnregistered(AliMUONVTrackerData*)",param);
  
}

//_____________________________________________________________________________
AliMUONVTrackerData*
AliMUONPainterDataRegistry::DataSource(const char* name) const
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
AliMUONPainterDataRegistry*
AliMUONPainterDataRegistry::Instance()
{
  /// Get unique instance of this class
  if ( !fgInstance ) fgInstance = new AliMUONPainterDataRegistry;
  return fgInstance;
}

//_____________________________________________________________________________
AliMUONVTrackerData*
AliMUONPainterDataRegistry::InteractiveReadOutConfig() const
{
  /// Return an object that contains the parts of the detector selected
  /// (using the mouse) to be part of the readout.

  if (!fInteractiveReadOutConfig) CreateInteractiveReadOutConfig();
  return fInteractiveReadOutConfig;
}

//_____________________________________________________________________________
void 
AliMUONPainterDataRegistry::Print(Option_t* opt) const
{
  /// Printout
  TString sopt(opt);
  sopt.ToUpper();
  
  cout << "Number of data readers = " << NumberOfDataMakers() << endl;
  
  if ( sopt.Contains("FULL") || sopt.Contains("READER") || sopt.Contains("DATA") )
  {
    TIter next(fDataMakers);
    AliMUONVTrackerDataMaker* reader;
    
    while ( ( reader = static_cast<AliMUONVTrackerDataMaker*>(next()) ) )
    {
      if ( sopt.Contains("DATA") ) 
      {
        AliMUONVTrackerData* data = reader->Data();
        if ( data ) data->Print();
      }
      else
      {
        reader->Print();
      }
    }
  }
}

//_____________________________________________________________________________
void
AliMUONPainterDataRegistry::Register(AliMUONVTrackerDataMaker* reader)
{
  /// reader is adopted, i.e. the registry becomes the owner of it.
  fDataMakers->AddLast(reader);
  DataMakerWasRegistered(reader);
  if ( reader->Data() ) DataSourceWasRegistered(reader->Data());
}

//_____________________________________________________________________________
Int_t 
AliMUONPainterDataRegistry::NumberOfDataMakers() const
{
  /// The number of data readers we handle
  return fDataMakers->GetLast()+1;
}

//_____________________________________________________________________________
void 
AliMUONPainterDataRegistry::DeleteZombies()
{
  /// Delete zombies
  fZombies->Delete();
}

//_____________________________________________________________________________
Bool_t 
AliMUONPainterDataRegistry::Unregister(AliMUONVTrackerDataMaker* reader)
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
