//
//  AliEveESDKinks.h
//
//  Created by Jeremi Niedziela on 3/12/15
//  Main authors: Paraskevi Ganoti: 2009
//

#ifndef __AliEveESDKinks__
#define __AliEveESDKinks__

#include <AliEveKink.h>

#include <AliESDtrack.h>
#include <AliESDkink.h>
#include <AliExternalTrackParam.h>

class AliEveESDKinks
{
public:
    AliEveESDKinks(){}
    ~AliEveESDKinks(){}
    
    AliEveKinkList* Draw();
    TEvePointSet* DrawPoints();
private:
    void InitRecTrackDother(TEveRecTrack& rt, const AliExternalTrackParam* tp);
    void InitRecTrackDaughter(TEveRecTrack& rt, const AliExternalTrackParam* tp, TEveVector* svt,TEveVector* spt);
    AliEveKink* MakeKink(TEveTrackPropagator* rnrStyleMoth,TEveTrackPropagator* rnrStyleDaugh, AliESDtrack* moth, AliESDtrack* daug, AliESDkink* kink, Int_t i);
    
    
    AliEveESDKinks(const AliEveESDKinks&);
    AliEveESDKinks& operator=(const AliEveESDKinks&);
};

#endif
