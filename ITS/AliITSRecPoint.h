#ifndef ALIITSRECPOINT_H
#define ALIITSRECPOINT_H 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  $Id$
*/

///////////////////////////////////////////////////////////////////////////////
//  Reconstructed space point class for set:ITS   
//  Reconstructed points are expressed simultaneously in two different 
//  reference frames, both differing from the global system.
//  The first is referred to the sensor (see AliITSsegmentation for the
//  definition) and each point is represented by two coordinates: fXloc and
//  fZloc. This system in the code is referred to as "local"
//  The second is used for tracking (V2, SA and MI versions) and the X axis 
//  represents the radial coordinate (this system is, in the bending plane, 
//  a rotated system w.r.t. the global reference system). 
//  Each reaconstructed point is represented by two coordinates: fY and fZ, 
//  inherited from AliCluster. This system in the code is referred to as 
//  "trackingV2".
///////////////////////////////////////////////////////////////////////////////

#include <AliCluster.h>
#include <Riostream.h>
#include <AliLog.h>
#include <AliITSgeom.h>


class AliITSRecPoint : public AliCluster {
 public:
  AliITSRecPoint();
  AliITSRecPoint(AliITSgeom* geom);
  AliITSRecPoint(Int_t *lab,Float_t *hit, Int_t *info);
  AliITSRecPoint(Int_t module,AliITSgeom* geom,Int_t *lab,
		 Float_t *hit, Int_t *info);
  AliITSRecPoint(const AliITSRecPoint& pt);
  AliITSRecPoint& operator=(const AliITSRecPoint &source);
  
  virtual ~AliITSRecPoint() {}; // distructor
  Bool_t IsSortable() const {return kTRUE;} // allows for sorting
  Int_t  *GetTracks(){return fTracks;}// Returns pointer to track array
  Int_t   GetNTracks() const {return 3;} // returns track array size
  Float_t GetDetLocalX() const {return fXloc;} // gets fX
  Float_t GetDetLocalZ() const {return fZloc;} // gets fZ
  Float_t GetdEdX() const {return fdEdX;} // gets fdEdX
  Float_t GetSigmaDetLocX2() const {return fSigmaY2;} // gets fSigmaX2
  void SetdEdX(Float_t dedx){fdEdX=dedx;} // sets fdEdX
  void SetSigmaDetLocX2(Float_t sx2){fSigmaY2=sx2;} // sets fSigmaX2
  Int_t Compare(const TObject *) const {return 0;} //to be defined
  void Print(ostream *os); 
  // Reads in the content of this class in the format of Print
  void Read(istream *is);
  virtual void Print(Option_t *option="") const {TObject::Print(option);}
  virtual Int_t Read(const char *name) {return TObject::Read(name);}
  
  void SetITSgeom(AliITSgeom* geom) {fGeom=geom;}

  virtual void SetY(Float_t  y ){fY=y;
    AliError("For consistency, Use method SetYZ. Data members are only partially set\n");}
  virtual void SetZ(Float_t  z ){fZ=z;
    AliError("For consistency, Use method SetYZ. Data members are only partially set\n");}
  void SetYZ(Int_t module, Float_t y, Float_t z){
    fY=y;fZ=z;
    if(fGeom)fGeom->TrackingV2ToDetL(module,y,z,fXloc,fZloc);
    else AliError("Geometry not set. \n");
  }
  void SetXZ(Int_t module, Float_t x, Float_t z){
    fXloc=x;fZloc=z;
    if(fGeom)fGeom->DetLToTrackingV2(module,x,z,fY,fZ);
    else AliError("Geometry not set. Nothing done!!!!!\n");
  }
  void Use(Int_t = 0) {fQ=-fQ;}
  void UnUse() {fQ=TMath::Abs(fQ);}
  void SetQ(Float_t q) {fQ=q;}
  void SetDetectorIndex(Int_t i) { fIndex=i; }
  void SetLayer(Int_t layer) {fLayer=layer;}
  void SetNz(Int_t nz) {fNz =nz;}
  void SetNy(Int_t ny){fNy=ny;}
  void SetChargeRatio(Float_t ratio) { fChargeRatio = ratio;}
  void SetPhiR(Float_t y) { fChargeRatio=y; }
  void SetType(Int_t type){ fType=type;}
  void SetDeltaProbability(Float_t prob){fDeltaProb = prob;}
  
  Int_t IsUsed() const {return (fQ<0)?1:0;}
  Float_t GetQ() const {return TMath::Abs(fQ);}
  Int_t GetDetectorIndex() const { return 0x3FF&fIndex; }
  Int_t GetLayer() const {return fLayer;}
  Int_t GetNz() const {return fNz;}
  Int_t GetNy() const {return fNy;}
  Float_t GetChargeRatio() const {return fChargeRatio;}
  Float_t GetPhiR() const {return fChargeRatio;}
  Int_t GetPindex() const { return 0xFFF00000&fIndex; }  //SSD clusters only
  Int_t GetNindex() const { return 0xFFC00&fIndex; }  //SSD clusters only
  Int_t GetType() const {return fType;}  // type of the cluster
  Float_t GetDeltaProbability() const{return fDeltaProb;} //probability to belong to the delta ray
  
  // The following two methods deal with the misaligned x position
  // of a cluster in the tracking coordidate system 
  virtual Float_t GetX() const {return fX;}
  virtual void SetX(Float_t x) {fX=x;}

 protected:
  Float_t   fX;        // X coordinate in the misaligned tracking system 

  Float_t   fXloc ;        //X of cluster (local coordinates)
  Float_t   fZloc ;        //Z of cluster (local coordinates)
  Float_t   fdEdX;      //dE/dX inside this cluster

  Int_t    fIndex;    // detector index
  Float_t  fQ ;       // Q of cluster (in ADC counts)
  Char_t   fLayer;    // layer number
  Short_t  fNz;       //number of digits in Z direction
  Short_t  fNy;       //number of digits in y direction 
  Float_t  fChargeRatio; //charge ratio
  Int_t    fType;         //quality factor of the cluster
  Float_t  fDeltaProb;    // probability to be deleta electron
    
  AliITSgeom* fGeom;     //!pointer to ITS geometry

  ClassDef(AliITSRecPoint,3)  // AliITSRecPoint class
};
// Input and output function for standard C++ input/output.
ostream& operator<<(ostream &os,AliITSRecPoint &source);
istream& operator>>(istream &is,AliITSRecPoint &source);
#endif
