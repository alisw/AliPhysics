#ifndef AliL3Model_Track
#define AliL3Model_Track

#include "AliL3Track.h"


struct ClusterComp {
  Bool_t fEmpty;
  Float_t fDTime;
  Float_t fDPad;
  Float_t fDCharge;
  Float_t fDShape;
};


class AliL3ModelTrack : public AliL3Track {

 private:
  
  Short_t fClusterCharge; //Average cluster charge
  ClusterComp *fClusters; //!
  Short_t fNClusters;
  Int_t fOverlap;

  //Crossing points with padrows
  Float_t *fPad; //!
  Float_t *fTime; //!
  
 public:
  AliL3ModelTrack();
  virtual ~AliL3ModelTrack();
  
  void Init(Int_t slice,Int_t patch);
  void SetCluster(Float_t dpad,Float_t dtime,Float_t charge,Float_t sigmaY2,Float_t sigmaZ2);
  
  void SetPadHit(Int_t row,Float_t f) {fPad[row]=f;}
  void SetTimeHit(Int_t row,Float_t f) {fTime[row]=f;}
  void SetOverlap(Int_t i) {fOverlap=i;}
  
  Int_t GetOverlap() {return fOverlap;}
  Float_t GetPadHit(Int_t row) {return fPad[row];}
  Float_t GetTimeHit(Int_t row) {return fTime[row];}

  ClassDef(AliL3ModelTrack,1)

};

#endif
