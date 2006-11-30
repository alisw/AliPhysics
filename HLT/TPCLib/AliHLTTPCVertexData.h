// @(#) $Id$
// Original: AliHLTVertexData.h,v 1.2 2003/07/27 21:02:09 loizides 

#ifndef VERTEXDATA_H
#define VERTEXDATA_H

#include "AliHLTTPCRootTypes.h"

struct AliHLTTPCVertexData{
    Double_t fX;
    Double_t fY;
    Double_t fZ;
    Double_t fXErr;
    Double_t fYErr;
    Double_t fZErr;
};
typedef struct AliHLTTPCVertexData AliHLTTPCVertexData;

#endif /* VERTEXDATA_H */
