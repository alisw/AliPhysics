#ifndef SPACEPOINTDATA_H
#define SPACEPOINTDATA_H

#include "AliL3RootTypes.h"
struct AliL3SpacePointData{
    Float_t fX;
    Float_t fY;
    Float_t fZ;
    UInt_t fID;  //contains slice patch and number
    UChar_t fPadRow;
    Float_t fXYErr;
    Float_t fZErr;
};
typedef struct AliL3SpacePointData AliL3SpacePointData;


#endif /* SPACEPOINTDATA_H */
