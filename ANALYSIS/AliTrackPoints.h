#ifndef AliTrackPoints_H
#define AliTrackPoints_H
//_________________________________
////////////////////////////////////////////////////////////
//                                                        //
// class AliTrackPoints                                //
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

class AliTrackPoints: public TObject
{
  public:
    typedef enum
      {//defines the type and the track points are calculated
        kITS = 1,//points calculated in ITS, made out of the vector at vertex
        kITSInnerFromVertexOuterFromTPC //for ITS, for 3 inner layers calculated out of the vector at vertex
                                        //         for 3 outer ---------------//------------------ at inner TPC
      } ETypes;

    AliTrackPoints();
    AliTrackPoints(AliTrackPoints::ETypes type, AliESDtrack* track, Float_t mf);
    AliTrackPoints(Int_t n, AliESDtrack* track, Float_t mf, Float_t dr=30,Float_t r0 = 84.1); //min TPC R  = 84.1; max TPC R =  246.6cm, 
    AliTrackPoints(Int_t n, AliTPCtrack* track, Float_t dr=30, Float_t r0 = 84.1); //min TPC R  = 84.1; max TPC R =  246.6cm, 
    
    virtual ~AliTrackPoints();
    
    Double_t AvarageDistance(const AliTrackPoints& tr);
    void PositionAt(Int_t n, Float_t &x, Float_t &y, Float_t &z);
    void Move(Float_t x, Float_t y, Float_t z);
    Int_t GetN() const {return fN;}

    Int_t GetDebug() const {return fgDebug;}
    void  SetDebug(Int_t deblevel){fgDebug = deblevel;} 
    static void Testtpc(Int_t entr);
    static void Testesd(Int_t entr,const char* fname = "AliESDs.root");
    void Print(Option_t* /*option*/ = "") const;
    
  protected:
    void MakePoints( Float_t dr, Float_t r0, Double_t x, Double_t* par, Double_t c, Double_t alpha);
    void MakeITSPoints(AliESDtrack* track);
    void MakeITSPointsInnerFromVertexOuterFromTPC(AliESDtrack* track,Float_t mf);
    
  private:
    AliTrackPoints(const AliTrackPoints& in);
    AliTrackPoints& operator=(const AliTrackPoints& in);

    Int_t    fN;//number of points
    Float_t* fX;//[fN]positions at x
    Float_t* fY;//[fN]positions at y
    Float_t* fZ;//[fN] positions at z
//    Float_t* fR;//! [fN] radii
    static Int_t fgDebug;//! debug level
    ClassDef(AliTrackPoints,1)
};
#endif
