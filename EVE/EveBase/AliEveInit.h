//
//  AliEveInit.h
//  xAliRoot
//
//  Created by Jeremi Niedziela on 01/06/15.
//
//

#ifndef __AliEveInit__
#define __AliEveInit__

#include <AliEveEventManager.h>

#include <TString.h>
#include <TEnv.h>

class AliEveInit
{
public:
    AliEveInit(const TString& path = ".",
                  AliEveEventManager::EDataSource defaultDataSource= AliEveEventManager::kSourceOffline,
                  bool storageManager=false);
    ~AliEveInit(){};
    
    static void GetConfig(TEnv *settings);
    static void AddMacros();
private:
    const TString& fPath;
    
    bool fShowHLTESDtree;
    
    void Init();
    void ImportMacros();
};

#endif
