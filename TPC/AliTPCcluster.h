#ifndef ALITPCCLUSTER_H
#define ALITPCCLUSTER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------
//                    TPC Cluster Class
//
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------

#include <TObject.h>

//_____________________________________________________________________________
class AliTPCcluster : public TObject {
public:
  AliTPCcluster();
  AliTPCcluster(Float_t *hits, Int_t *lab);
  void Use();
  void SetLabel(Int_t lab, Int_t i);
  void SetQ(Float_t q);
  void SetY(Float_t y);
  void SetZ(Float_t z);
  void SetSigmaY2(Float_t sy2);
  void SetSigmaZ2(Float_t sz2);

  Int_t IsUsed() const;
  Int_t GetLabel(Int_t i) const;
  Float_t GetQ() const;
  Float_t GetY() const;
  Float_t GetZ() const;
  Float_t GetSigmaY2() const;
  Float_t GetSigmaZ2() const;

private:
  Int_t     fTracks[3];//labels of overlapped tracks
  Float_t   fQ ;       //Q of cluster (in ADC counts)
  Float_t   fY ;       //Y of cluster
  Float_t   fZ ;       //Z of cluster
  Float_t   fSigmaY2;  //Sigma Y square of cluster
  Float_t   fSigmaZ2;  //Sigma Z square of cluster
  
  ClassDef(AliTPCcluster,1)  // Time Projection Chamber clusters
};

inline void AliTPCcluster::Use() {
  //if fQ<0 cluster is already associated with a track
  fQ=-fQ;
}

inline Int_t AliTPCcluster::IsUsed() const {
  //is this cluster already associated with any track ?
  return (fQ<0) ? 1 : 0;
}

inline Int_t AliTPCcluster::GetLabel(Int_t i) const {
  //return track label
  return fTracks[i];
}

inline Float_t AliTPCcluster::GetQ() const {
  //just to calm down our rule checker
  return fQ;
}

inline Float_t AliTPCcluster::GetY() const {
  //just to calm down our rule checker
  return fY;
}

inline Float_t AliTPCcluster::GetZ() const {
  //just to calm down our rule checker
  return fZ;
}

inline Float_t AliTPCcluster::GetSigmaY2() const {
  //just to calm down our rule checker
  return fSigmaY2;
}

inline Float_t AliTPCcluster::GetSigmaZ2() const {
  //just to calm down our rule checker
  return fSigmaZ2;
}

inline void AliTPCcluster::SetLabel(Int_t lab, Int_t i) {
  //just to calm down our rule checker
  fTracks[i]=lab;
}

inline void AliTPCcluster::SetQ(Float_t q) {
  //just to calm down our rule checker
  fQ=q;
}

inline void AliTPCcluster::SetY(Float_t y) {
  //just to calm down our rule checker
  fY=y;
}

inline void AliTPCcluster::SetZ(Float_t z) {
  //just to calm down our rule checker
  fZ=z;
}

inline void AliTPCcluster::SetSigmaY2(Float_t sy2) {
  //just to calm down our rule checker
  fSigmaY2=sy2;
}

inline void AliTPCcluster::SetSigmaZ2(Float_t sz2) {
  //just to calm down our rule checker
  fSigmaZ2=sz2;
}

#endif


