// $Id$
// Author: Matevz Tadel 2009

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveConfigManager_H
#define AliEveConfigManager_H

#include "TObject.h"
#include "TGMenu.h"

class AliEveConfigManager : public TObject
{
public:
    static AliEveConfigManager* Instance(bool storageManager=false);
    
    void AliEvePopupHandler(int id);
    void SetEventInEventManager();
    void StorageManagerChangedState(int state);

    void ConnectEventManagerSignals();
    
private:
    static AliEveConfigManager* fInstance;  // Main instance.
    AliEveConfigManager(bool storageManager=false);
    ~AliEveConfigManager(){}
    
    TGPopupMenu *fAliEvePopup;  // AliEve menu
    TGPopupMenu *fStoragePopup; // Storage Manager menu
    
    enum EAliEveMenu{ kStorageListEvents, kStorageMarkEvent, kPreferences };
    
    AliEveConfigManager(const AliEveConfigManager&);            // Not implemented
    AliEveConfigManager& operator=(const AliEveConfigManager&); // Not implemented
    
    ClassDef(AliEveConfigManager, 0); // Short description.
};

#endif
