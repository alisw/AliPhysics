#ifndef ALIHBTTRACKPOINTS_H
#define ALIHBTTRACKPOINTS_H

#include <TObject.h>

class AliTPCtrack;
class AliESDtrack;

class AliHBTTrackPoints: public TObject
{
  public:
    AliHBTTrackPoints();
    AliHBTTrackPoints(Int_t n, AliTPCtrack* track, Float_t dr=30,Float_t r0 = 84.1); //min TPC R  = 84.1; max TPC R =  246.6cm, 
    
    virtual ~AliHBTTrackPoints();
    
    Double_t AvarageDistance(const AliHBTTrackPoints& tr);
    void PositionAt(Int_t n, Float_t &x, Float_t &y, Float_t &z);
    Int_t GetDebug() const {return 10;}
    static void tp(Int_t entr);
  protected:
  private:
    Int_t    fN;//number of points
    Float_t* fX;//[fN]
    Float_t* fY;//[fN]
    Float_t* fZ;//[fN]
    Float_t* fR;//! [fN] radii
    
    ClassDef(AliHBTTrackPoints,1)
};
#endif
