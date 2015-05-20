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

class AliEveGeomGentle
{
public:
    AliEveGeomGentle();
    ~AliEveGeomGentle();
    
    TEveGeoShape* GetGeomGentle(bool register_as_global=kTRUE);
    TEveGeoShape* GetGeomGentleRphi();
    TEveGeoShape* GetGeomGentleRhoz();
    TEveGeoShape* GetGeomGentleTRD(Color_t color=3);
    TEveGeoShape* GetGeomGentleMUON(bool updateScene = kTRUE, Color_t color=3);

private:
    void DrawDeep(TEveGeoShape *gsre, Color_t color);
};


#endif
