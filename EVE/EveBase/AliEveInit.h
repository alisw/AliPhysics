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

class AliEveInit
{
public:
    AliEveInit(const TString& path = ".",
                  const TString& cdbUri="local:///local/cdb",
                  AliEveEventManager::EDataSource defaultDataSource= AliEveEventManager::kSourceOffline,
                  bool storageManager=false);
    ~AliEveInit(){};
private:
    const TString& fCDBuri;
    const TString& fPath;
    
    bool fShowHLTESDtree;
    
    void Init();
    void ImportMacros();
};

#endif
