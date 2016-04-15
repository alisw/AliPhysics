// $Id$
// Author: Matevz Tadel 2009

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveConfigManager.h"

#include "AliEveEventManager.h"
#include "AliEvePreferencesWindow.h"
#ifdef ZMQ
#include "AliStorageAdministratorPanelListEvents.h"
#include "AliStorageAdministratorPanelMarkEvent.h"
#endif

#include <TEveBrowser.h>

#include <iostream>

using namespace std;

ClassImp(AliEveConfigManager)

AliEveConfigManager* AliEveConfigManager::fInstance = 0;

AliEveConfigManager* AliEveConfigManager::Instance(bool storageManager)
{
    // Get main instance.
    if (fInstance){
        cout<<"AliEveConfigManger -- Master already initialized"<<endl;
    }
    else{
        cout<<"AliEveConfigManger -- Creating instance of AliEveConfigManager"<<endl;
        fInstance = new AliEveConfigManager(storageManager);
    }
    return fInstance;
}

AliEveConfigManager::AliEveConfigManager(bool storageManager) :
TObject(),
fAliEvePopup(0),
fStoragePopup(0)
{
    fAliEvePopup = new TGPopupMenu(gClient->GetRoot());
    fAliEvePopup->AddEntry("&Preferences", kPreferences);
    fAliEvePopup->Connect("Activated(Int_t)", "AliEveConfigManager",this, "AliEvePopupHandler(Int_t)");
    
    //Storage Manager:
#ifdef ZMQ
    if(storageManager)
    {
        gEve->GetBrowser()->StartEmbedding(0);
        AliStorageAdministratorPanelListEvents* fListEventsTab = AliStorageAdministratorPanelListEvents::GetInstance();
        gEve->GetBrowser()->StopEmbedding("List");
        
        fListEventsTab->Connect("SelectedEvent()","AliEveConfigManager",this,"SetEventInEventManager()");
    }
#endif
    
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,25,4)
    TGMenuBar *mBar = gEve->GetBrowser()->GetMenuBar();
    mBar->AddPopup("&AliEve", fAliEvePopup, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));
    ((TGCompositeFrame*)mBar->GetParent()->GetParent())->Layout();
    
    gEve->GetBrowser()->GetTopMenuFrame()->Layout();
#else
    // Uber hack as TRootBrowser does not provede manu-bar getter.
    TGFrameElement   *xxFE = (TGFrameElement*)   gEve->GetBrowser()->GetList()->First();
    TGCompositeFrame *xxCF = (TGCompositeFrame*) xxFE->fFrame;
    xxFE = (TGFrameElement*)   xxCF->GetList()->First();
    xxCF = (TGCompositeFrame*) xxFE->fFrame;
    xxFE = (TGFrameElement*)   xxCF->GetList()->First();
    xxCF = (TGCompositeFrame*) xxFE->fFrame;
    xxFE = (TGFrameElement*)   xxCF->GetList()->First();
    TGMenuBar *mBar = (TGMenuBar*) xxFE->fFrame;
    mBar->AddPopup("&AliEve", fAliEvePopup, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0));
    ((TGCompositeFrame*)mBar->GetParent()->GetParent())->Layout();
    mBar->AddPopup("&Storage Manager",fStoragePopup,new TGLayoutHints(kLHintsTop | kLHintsLeft,0,4,0,0));
    ((TGCompositeFrame*)mBar->GetParent()->GetParent())->Layout();
#endif
}

void AliEveConfigManager::ConnectEventManagerSignals()
{
    AliEveEventManager *manager = AliEveEventManager::Instance();
    manager->Connect("StorageManagerOk()","AliEveConfigManager",this,"StorageManagerChangedState(=1)");
    manager->Connect("StorageManagerDown()","AliEveConfigManager",this,"StorageManagerChangedState(=0)");
}

void AliEveConfigManager::AliEvePopupHandler(int id)
{
    // Handle user selections from AliEve popup.
    switch (id)
    {
        case kStorageMarkEvent:
        {
#ifdef ZMQ
            AliStorageAdministratorPanelMarkEvent *markEventWindow =
            AliStorageAdministratorPanelMarkEvent::GetInstance();
#endif
            break;
        }
        case kPreferences:
        {
            AliEvePreferencesWindow::Instance();
            
            break;
        }
            
        default:
        {
            cout<<"AliEveConfigManager -- Unknown menu entry."<<endl;
            break;
        }
    }
}

void AliEveConfigManager::SetEventInEventManager()
{
#ifdef ZMQ
    AliEveEventManager *manager = AliEveEventManager::Instance();
    AliStorageAdministratorPanelListEvents* fListEventsTab = AliStorageAdministratorPanelListEvents::GetInstance();
    AliESDEvent *event = fListEventsTab->GetSelectedEvent();
    
    if(event)
    {
        cout<<"AliEveConfigManager -- setting event from Storage Manager"<<endl;
        
        manager->SetAutoLoad(kFALSE);
        AliEveDataSource *dataSource = manager->GetCurrentDataSource();
        dataSource->SetEventFromStorageManager(event);
        dataSource->NextEvent();
    }
#endif
}

void AliEveConfigManager::StorageManagerChangedState(int state)
{
#ifdef ZMQ
    AliEveEventManager *manager = AliEveEventManager::Instance();
    AliStorageAdministratorPanelListEvents* listEventsTab = AliStorageAdministratorPanelListEvents::GetInstance();
    
    if (state == 0)// storage manager is down
    {
        listEventsTab->SetOfflineMode(kTRUE);
    }
    else if(state == 1)// storage manager is alive
    {
        listEventsTab->SetOfflineMode(kFALSE);
    }
#endif
}


