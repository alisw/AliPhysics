//
//  AliEveESDCascades.h
//
//  Created by Jeremi Niedziela on 3/12/15
//  Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007
//

#ifndef __AliEveESDCascades__
#define __AliEveESDCascades__

#include <AliEveCascade.h>

#include <AliESDcascade.h>
#include <AliESDtrack.h>
#include <AliExternalTrackParam.h>

class AliEveESDCascades
{
public:
    AliEveESDCascades(){}
    ~AliEveESDCascades(){}
    
    AliEveCascadeList* Draw();
    TEvePointSet* DrawPoints();
private:
    void InitRecTrack(TEveRecTrack& rt, const AliExternalTrackParam* tp);
    AliEveCascade* MakeCascade(TEveTrackPropagator* rnrStyleBac,TEveTrackPropagator* rnrStyleNeg,TEveTrackPropagator* rnrStylePos, AliESDVertex* primVtx,AliESDtrack* bac, AliESDcascade* cascade, Int_t i);
    
    AliEveESDCascades(const AliEveESDCascades&);
    AliEveESDCascades& operator=(const AliEveESDCascades&);
};

#endif
