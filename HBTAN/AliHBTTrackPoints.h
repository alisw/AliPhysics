#ifndef ALIHBTTRACKPOINTS_H
#define ALIHBTTRACKPOINTS_H
//_________________________________
////////////////////////////////////////////////////////////
//                                                        //
// class AliHBTTrackPoints                                //
//                                                        //
// used by Anti-Merging cut                               //
// contains set of poits the lay on track trajectory      //
// according to reconstructed track parameters -          //
// NOT CLUSTERS POSITIONS!!!                              //
// Anti-Merging cut is applied only on tracks coming from //
// different events (that are use to fill deniminators)   //
//                                                        //
////////////////////////////////////////////////////////////
#include <TObject.h>

class AliTPCtrack;
class AliESDtrack;

class AliHBTTrackPoints: public TObject
{
  public:
    enum ETypes{kITS = 1};

    AliHBTTrackPoints();
    AliHBTTrackPoints(Int_t n, AliTPCtrack* track, Float_t dr=30, Float_t r0 = 84.1); //min TPC R  = 84.1; max TPC R =  246.6cm, 
    AliHBTTrackPoints(Int_t n, AliESDtrack* track, Float_t mf, Float_t dr=30,Float_t r0 = 84.1); //min TPC R  = 84.1; max TPC R =  246.6cm, 
    AliHBTTrackPoints(AliHBTTrackPoints::ETypes type, AliESDtrack* track);
    AliHBTTrackPoints(const AliHBTTrackPoints& in);
    
    virtual ~AliHBTTrackPoints();
    AliHBTTrackPoints& operator=(const AliHBTTrackPoints& in);
    
    Double_t AvarageDistance(const AliHBTTrackPoints& tr);
    void PositionAt(Int_t n, Float_t &x, Float_t &y, Float_t &z);
    void Move(Float_t x, Float_t y, Float_t z);

    Int_t GetDebug() const {return fgDebug;}
    void  SetDebug(Int_t deblevel){fgDebug = deblevel;} 
    static void testtpc(Int_t entr);
    static void testesd(Int_t entr,const char* fname = "AliESDs.root");

  protected:
    void MakePoints( Float_t dr, Float_t r0, Double_t x, Double_t* par, Double_t c, Double_t alpha);
    void MakeITSPoints(AliESDtrack* track);
    
  private:
    Int_t    fN;//number of points
    Float_t* fX;//[fN]positions at x
    Float_t* fY;//[fN]positions at y
    Float_t* fZ;//[fN] positions at z
//    Float_t* fR;//! [fN] radii
    static Int_t fgDebug;//! debug level
    ClassDef(AliHBTTrackPoints,1)
};
#endif
