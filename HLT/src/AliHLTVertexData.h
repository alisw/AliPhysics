// @(#) $Id$

#ifndef VERTEXDATA_H
#define VERTEXDATA_H

#include "AliHLTRootTypes.h"

struct AliHLTVertexData{
    Double_t fX;
    Double_t fY;
    Double_t fZ;
    Double_t fXErr;
    Double_t fYErr;
    Double_t fZErr;
};
typedef struct AliHLTVertexData AliHLTVertexData;

#endif /* VERTEXDATA_H */
