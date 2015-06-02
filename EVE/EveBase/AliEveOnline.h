//
//  AliEveOnline.h
//  xAliRoot
//
//  Created by Jeremi Niedziela on 11/05/15.
//
//

#ifndef __AliEveOnline__
#define __AliEveOnline__

class AliEveOnline {
public:
    AliEveOnline(bool storageManager=false);
    ~AliEveOnline();
    
private:
    void InitImportMacros();
};


#endif
