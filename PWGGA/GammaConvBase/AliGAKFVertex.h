#ifndef ALIGAKFVERTEX_H
#define ALIGAKFVERTEX_H

#ifndef PWGGAUSEKFPARTICLE
#include "AliKFVertex.h"
#define AliGAKFVertex AliKFVertex
#else
#include "AliVVertex.h"

//#include "AliKFVertex.h"
#include "KFVertex.h"
#include "KFPVertex.h"

class AliGAKFVertex : public KFVertex {
public :
    AliGAKFVertex() : KFVertex() {}
    AliGAKFVertex(const AliVVertex&);

    ClassDef(AliGAKFVertex, 1)
};
#endif // PWGGAUSEKFPARTICLE
#endif // ALIGAKFVERTEX_H
