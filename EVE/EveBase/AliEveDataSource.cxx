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

void AliEveDataSource::StorageManagerOk(){}
void AliEveDataSource::StorageManagerDown(){}
void AliEveDataSource::EventServerOk(){}
void AliEveDataSource::EventServerDown(){}
