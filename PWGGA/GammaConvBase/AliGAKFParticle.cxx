#include "AliLog.h"
#include "AliGAKFParticle.h"

#ifdef PWGGAUSEKFPARTICLE
ClassImp(AliGAKFParticle);


AliGAKFParticle::AliGAKFParticle( const AliVTrack &track, int PID ) : KFParticle() {
    AliError("THIS CLASS MEMBER NEEDS TO BE IMPLEMENTED");
}

AliGAKFParticle::AliGAKFParticle( const AliExternalTrackParam &track, double Mass, int Charge ) : KFParticle() {
    AliError("THIS CLASS MEMBER NEEDS TO BE IMPLEMENTED");
}

void AliGAKFParticle::GetDStoParticle(AliGAKFParticle&, double &, double &) const {
    AliError("THIS CLASS MEMBER NEEDS TO BE IMPLEMENTED");
}

void AliGAKFParticle::TransportToPoint(double [3]) {
    AliError("THIS CLASS MEMBER NEEDS TO BE IMPLEMENTED");
}

void AliGAKFParticle::ConstructGamma(const KFParticle&, const KFParticle&) {
    AliError("THIS CLASS MEMBER NEEDS TO BE IMPLEMENTED");
}

void AliGAKFParticle::GetArmenterosPodolanski(KFParticle&, KFParticle&, double *&) {
    AliErrorClass("THIS CLASS MEMBER NEEDS TO BE IMPLEMENTED");
}

#endif // PWGGAUSEKFPARTICLE