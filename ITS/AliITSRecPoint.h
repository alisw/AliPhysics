#ifndef ALIITSRECPOINT_H
#define ALIITSRECPOINT_H 

////////////////////////////////////////////////////
//  Reconstructed space point class for set:ITS   //
////////////////////////////////////////////////////

#include <TObject.h>


class AliITSRecPoint : public TObject {


 public:

  AliITSRecPoint();
  virtual ~AliITSRecPoint() {}; // distructor
  Bool_t IsSortable() const {return kTRUE;} // allows for sorting
  Int_t   GetLabel(Int_t i) {return fTracks[i];} // get track label
  Float_t GetX(){return fX;} // gets fX
  Float_t GetZ(){return fZ;} // gets fZ
  Float_t GetQ(){return fQ;} // gets fQ
  Float_t GetdEdX(){return fdEdX;} // gets fdEdX
  Float_t GetSigmaX2(){return fSigmaX2;} // gets fSigmaX2
  Float_t GetSigmaZ2(){return fSigmaZ2;} // gets fSigmaZ2
  void SetLabel(Int_t i, Int_t lab){fTracks[i]=lab;} // sets track label
  void SetX(Float_t x){fX=x;} // sets fX
  void SetZ(Float_t z){fZ=z;} // sets fZ
  void SetQ(Float_t q){fQ=q;} // sets fQ
  void SetdEdX(Float_t dedx){fdEdX=dedx;} // sets fdEdX
  void SetSigmaX2(Float_t sx2){fSigmaX2=sx2;} // sets fSigmaX2
  void SetSigmaZ2(Float_t sz2){fSigmaZ2=sz2;} // sets fSigmaZ2
  void  Use() {
    //if fQ<0 cluster is already associated with a track
    fQ=-fQ;
  }
  Int_t IsUsed() const {return (fQ<0) ? 1 : 0;} // checks Use condision
  Int_t Compare(TObject *o) {
    //to be defined
    return 0;
  } 

 public:
    Int_t     fTracks[3]; //labels of overlapped tracks
    Float_t   fX ;        //X of cluster
    Float_t   fZ ;        //Z of cluster
    Float_t   fQ ;        //Q of cluster (in ADC counts)
    Float_t   fdEdX;      //dE/dX inside this cluster
    Float_t   fSigmaX2;   //Sigma X square of cluster
    Float_t   fSigmaZ2;   //Sigma Z square of cluster

  ClassDef(AliITSRecPoint,1)  // AliITSRecPoint class
};

#endif



