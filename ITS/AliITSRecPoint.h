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

class AliITSRecPoint : public AliCluster {
 public:
  AliITSRecPoint();
  AliITSRecPoint(Int_t *lab,Float_t *hit, Int_t *info, Bool_t local = kFALSE);
  AliITSRecPoint(const AliITSRecPoint& pt);
  AliITSRecPoint& operator=(const AliITSRecPoint &source);
  
  virtual ~AliITSRecPoint() {}; // distructor
  Bool_t  IsSortable() const {return kTRUE;} // allows for sorting
  Float_t GetDetLocalX() const {return fXloc;} // gets fX
  Float_t GetDetLocalZ() const {return fZloc;} // gets fZ
  Float_t GetdEdX() const {return fdEdX;} // gets fdEdX
  Float_t GetSigmaDetLocX2() const {return GetSigmaY2();} // gets fSigmaX2
  void    SetdEdX(Float_t dedx){fdEdX=dedx;} // sets fdEdX
  Int_t Compare(const TObject *) const {return 0;} //to be defined
  void Print(ostream *os); 
  // Reads in the content of this class in the format of Print
  void Read(istream *is);
  virtual void Print(Option_t *option="") const {TObject::Print(option);}
  virtual Int_t Read(const char *name) {return TObject::Read(name);}
  
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
  void SetDriftTime(Float_t tim) {fDriftTime=tim;}
  void SetDriftSide(Int_t sid) {fDriftSide=sid;}
 
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
  Int_t GetType() const {return fType;}  // type of the cluster (for SPD the number of pixels in the cluster)
  Float_t GetDeltaProbability() const{return fDeltaProb;} //probability to belong to the delta ray
  Float_t GetDriftTime() const{return  fDriftTime;}
  Int_t GetDriftSide() const {return fDriftSide;}
  Int_t GetNpixels() const; // for SPD returns fType, i.e. the number of pixels in the cluster (-1 for SDD and SSD)
  Int_t GetSPDclusterType() const; // for SPD returns cluster type according to conventional numbering (-1 for SDD and SSD)
  Int_t GetSDDclusterType() const; 
  Int_t GetSSDclusterType() const; 
  static void DecodeSDDclusterType(Int_t cluType, Int_t &cluSizAn, Int_t& cluSizTb, Int_t &drSide);

  Int_t GetClusterType() const {
    if(fLayer<=1) return GetSPDclusterType();
    if(fLayer==2 || fLayer==3) return GetSDDclusterType();
    return GetSSDclusterType();
  }
 protected:

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
  Float_t  fDeltaProb;    // probability to be delta electron
  Float_t  fDriftTime;    // drift time in SDD
  Char_t   fDriftSide;    // drift region in SDD (0=left=positive xlocal, 1=right)
    
  ClassDef(AliITSRecPoint,7)  // AliITSRecPoint class
};
// Input and output function for standard C++ input/output.
ostream& operator<<(ostream &os,AliITSRecPoint &source);
istream& operator>>(istream &is,AliITSRecPoint &source);
#endif
