#include "AliLog.h"
#include "AliGAKFVertex.h"
#ifdef PWGGAUSEKFPARTICLE
ClassImp(AliGAKFVertex)

AliGAKFVertex::AliGAKFVertex(const AliVVertex&) : KFVertex() {
    AliError("THIS CLASS MEMBER NEEDS TO BE IMPLEMENTED");
}
#endif // PWGGAUSEKFPARTICLE
