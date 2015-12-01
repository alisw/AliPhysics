//
//  AliEvePrimaryVertex.h
//
//  Created by Jeremi Niedziela on 1/12/15
//  Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007
//

#ifndef __AliEvePrimaryVertex__
#define __AliEvePrimaryVertex__

#include <TEveStraightLineSet.h>

class AliEvePrimaryVertex
{
public:
    AliEvePrimaryVertex(){}
    ~AliEvePrimaryVertex(){}
    
    enum EVertexType { kGlobal, kSPD, kTPC };
    enum EVertexStyle{ kCross, kEllipse, kBox };
    
    void PrimaryVertex(EVertexType type = kGlobal,EVertexStyle style = kCross, Bool_t use_sigma=kTRUE, Float_t fx=30, Float_t fy=30, Float_t fz=10);
    
private:
    TEveStraightLineSet* MakeVertexBox(const AliESDVertex* v, Bool_t use_sigma, Float_t fx, Float_t fy, Float_t fz);
    TEveStraightLineSet* MakeVertexEllipse(const AliESDVertex* v, Bool_t use_sigma, Float_t fx, Float_t fy, Float_t fz);
    TEveStraightLineSet* MakeVertexCross(const AliESDVertex* v, Bool_t use_sigma, Float_t fx, Float_t fy, Float_t fz);
    
    AliEvePrimaryVertex(const AliEvePrimaryVertex&);
    AliEvePrimaryVertex& operator=(const AliEvePrimaryVertex&);
};

#endif
