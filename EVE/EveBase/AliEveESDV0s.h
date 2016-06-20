//
//  AliEveESDV0s.h
//
//  Created by Jeremi Niedziela on 3/12/15
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007
//

#ifndef __AliEveESDV0s__
#define __AliEveESDV0s__

#include <AliEveV0.h>

#include <AliESDv0.h>
#include <AliESDVertex.h>


class AliEveESDV0s
{
public:
    AliEveESDV0s(){}
    ~AliEveESDV0s(){}
    
    AliEveV0List* Draw(Bool_t onFly=kFALSE);
    TEvePointSet* DrawPointsOffline();
    TEvePointSet* DrawPointsOnfly();
private:
    void InitRecTracks(TEveRecTrack& rt, const AliExternalTrackParam* tp);
    AliEveV0* MakeV0(TEveTrackPropagator* rnrStyleNeg,TEveTrackPropagator* rnrStylePos, AliESDVertex* primVtx, AliESDtrack* neg, AliESDtrack* pos, AliESDv0* v0, Int_t i);
    void FillPointSet(TEvePointSet* ps, Bool_t onFly);
    
    AliEveESDV0s(const AliEveESDV0s&);
    AliEveESDV0s& operator=(const AliEveESDV0s&);
};

#endif
