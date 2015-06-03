//
//  AliEveHLT.h
//  xAliRoot
//
//  Created by Jeremi Niedziela on 11/05/15.
//
//

#ifndef __AliEveHLT__
#define __AliEveHLT__

class AliEveHLT {
public:
    AliEveHLT(bool storageManager=false);
    ~AliEveHLT();
    static void OnNewEvent();
    
private:
    void InitImportMacros();
};


#endif
