#ifndef ALIGAKFVERTEX_H
#define ALIGAKFVERTEX_H

#include "AliVVertex.h"

//#include "AliKFVertex.h"
#include "KFVertex.h"

class AliGAKFVertex : public KFVertex {
public :
    AliGAKFVertex() : KFVertex() {}
    AliGAKFVertex(const AliVVertex&);

    ClassDef(AliGAKFVertex, 1)
};

#endif // ALIGAKFVERTEX_H