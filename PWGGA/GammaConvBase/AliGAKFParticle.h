#ifndef ALIGAKFPARTICLE_H
#define ALIGAKFPARTICLE_H

#include "AliExternalTrackParam.h"

//#include "AliKFParticle.h"
// includes added to play with KFParticle
#ifndef HomogeneousField
#define HomogeneousField 
#endif
#include "KFPTrack.h"
#include "KFParticle.h"

class AliGAKFParticle : public KFParticle {
public :
    AliGAKFParticle() : KFParticle() {}
    AliGAKFParticle( const AliGAKFParticle &d1, const AliGAKFParticle &d2, Bool_t gamma = kFALSE ) : KFParticle(d1,d2) {}
    AliGAKFParticle( const AliGAKFParticle &d1, const AliGAKFParticle &d2, 
		             const AliGAKFParticle &d3 ) : KFParticle(d1,d2,d3) {}
    AliGAKFParticle( const AliGAKFParticle &d1, const AliGAKFParticle &d2, 
		             const AliGAKFParticle &d3, const AliGAKFParticle &d4 ) : KFParticle(d1,d2) {}
    AliGAKFParticle( const AliVTrack &track, int PID );
    AliGAKFParticle( const AliExternalTrackParam &track, double Mass, int Charge );

    void GetDStoParticle(AliGAKFParticle&, double &, double &) const;
    void TransportToPoint(double [3]);
    void ConstructGamma(const KFParticle&, const KFParticle&);
    static void GetArmenterosPodolanski(KFParticle&, KFParticle&, double *&);

    ClassDef(AliGAKFParticle, 1)
};

#endif // ALIGAKFPARTICLE_H