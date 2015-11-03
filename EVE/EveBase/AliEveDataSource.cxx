//
//  AliEveDataSource
//
//  blame: Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch
//
//

#include "AliEveDataSource.h"

AliEveDataSource::AliEveDataSource(bool storageManager)
  : TNamed("","")
  , fCurrentData()
  , fSourceURL()
    
{
}

AliEveDataSource::~AliEveDataSource()
{
}

void AliEveDataSource::Init()
{
}

void AliEveDataSource::GotoEvent(Int_t /*event*/)
{
}

void AliEveDataSource::NextEvent()
{
}

void AliEveDataSource::SetEventFromStorageManager(AliESDEvent *event)
{
    
}

void AliEveDataSource::StorageManagerOk(){}
void AliEveDataSource::StorageManagerDown(){}
