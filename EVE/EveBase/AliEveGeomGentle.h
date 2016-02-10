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
    
    TEveGeoShape* GetSimpleGeom(char* detector);
    
    TEveGeoShape* GetGeomGentle(bool register_as_global=kTRUE);
    TEveGeoShape* GetGeomGentleRphi();
    TEveGeoShape* GetGeomGentleRhoz();
    TEveGeoShape* GetGeomGentleTRD();

private:
    TEnv fSettings;
    void DrawDeep(TEveGeoShape *gsre, Color_t color, Char_t transparency, bool drawLine);
};


#endif
