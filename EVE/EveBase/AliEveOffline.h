//
//  AliEveOffline.h
//  xAliRoot
//
//  Created by Jeremi Niedziela on 01/06/15.
//
//

#ifndef __AliEveOffline__
#define __AliEveOffline__

#include <TString.h>

class AliEveOffline
{
public:
    AliEveOffline(const TString& path = ".", const TString& cdbUri="local:///local/cdb");
    ~AliEveOffline();
private:
    const TString& fCDBuri;
    const TString& fPath;
    
    bool fShowHLTESDtree;
    
    void Init();
    void ImportMacros();
};

#endif
