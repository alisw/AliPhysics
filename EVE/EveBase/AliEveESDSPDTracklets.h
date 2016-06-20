//
//  AliEveESDSPDTracklets.h
//
//  Created by Jeremi Niedziela on 1/12/15
//  Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007
//

#ifndef __AliEveESDSPDTracklets__
#define __AliEveESDSPDTracklets__

#include <TEveElement.h>

class AliEveESDSPDTracklets
{
public:
    AliEveESDSPDTracklets(){}
    ~AliEveESDSPDTracklets(){}
    
    TEveElementList* Draw(Float_t radius=8, Width_t line_width=2,
                          Float_t dPhiWindow=0.080, Float_t dThetaWindow=0.025,
                          Float_t dPhiShift05T=0.0045);
    
private:
    AliEveESDSPDTracklets(const AliEveESDSPDTracklets&);
    AliEveESDSPDTracklets& operator=(const AliEveESDSPDTracklets&);
};

#endif
