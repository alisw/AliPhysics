#ifndef VERTEXDATA_H
#define VERTEXDATA_H

#include "AliL3RootTypes.h"

struct AliL3VertexData{
    Double_t fX;
    Double_t fY;
    Double_t fZ;
    Double_t fXErr;
    Double_t fYErr;
    Double_t fZErr;
};
typedef struct AliL3VertexData AliL3VertexData;

#endif /* VERTEXDATA_H */
