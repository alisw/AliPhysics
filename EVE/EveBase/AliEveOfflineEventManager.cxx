//
//  AliEveOfflineEventManager.cpp
//  xAliRoot
//
//  Created by Jeremi Niedziela on 01/06/15.
//
//

#include "AliEveOfflineEventManager.h"

ClassImp(AliEveOfflineEventManager);

AliEveOfflineEventManager::AliEveOfflineEventManager(bool storageManager) :
AliEveEventManager("offline")
{
    AliEveEventManager::SetMaster(this);
}

AliEveOfflineEventManager::~AliEveOfflineEventManager()
{
    
}