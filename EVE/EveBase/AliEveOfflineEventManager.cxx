//
//  AliEveOfflineEventManager.cpp
//  xAliRoot
//
//  Created by Jeremi Niedziela on 01/06/15.
//
//

#include "AliEveOfflineEventManager.h"

AliEveOfflineEventManager::AliEveOfflineEventManager(Int_t ev, bool storageManager) :
AliEveEventManager("offline",ev,storageManager)
{
    AliEveEventManager::SetMaster(this);
}

AliEveOfflineEventManager::~AliEveOfflineEventManager()
{
    
}