//
//  AliEveGeomGentle.h
//  xAliRoot
//
//  Created by Jeremi Niedziela on 11/05/15.
//
//

#ifndef __AliEveGeomGentle__
#define __AliEveGeomGentle__

#include <TEveGeoShape.h>
#include <TEnv.h>

class AliEveGeomGentle
{
public:
    AliEveGeomGentle();
    ~AliEveGeomGentle();
    
    TEveGeoShape* GetGeomGentle(bool register_as_global=kTRUE);
    TEveGeoShape* GetGeomGentleRphi();
    TEveGeoShape* GetGeomGentleRhoz();
    TEveGeoShape* GetGeomGentleTRD();
    TEveGeoShape* GetGeomGentleEMCAL();
    TEveGeoShape* GetGeomGentleZDC();
    TEveGeoShape* GetGeomGentleMUON(bool updateScene = kTRUE);

private:
    TEnv fSettings;
    void DrawDeep(TEveGeoShape *gsre);
};


#endif
