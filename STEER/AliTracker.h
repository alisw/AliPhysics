#ifndef ALITRACKER_H
#define ALITRACKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//                          class AliTracker
//   that is the base for AliTPCtracker, AliITStrackerV2 and AliTRDtracker
//       Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------------------------
#include <Rtypes.h>

class AliKalmanTrack;
class AliCluster;
class TFile;

class AliTracker {
public:
  AliTracker() { fX=fY=fZ=0.; fEventN=0; }
  virtual ~AliTracker(){}
  virtual Int_t Clusters2Tracks(const TFile *in, TFile *out)=0;
  virtual Int_t PropagateBack(const TFile *in, TFile *out)=0;
  void SetVertex(Double_t *xyz) { fX=xyz[0]; fY=xyz[1]; fZ=xyz[2]; }
  void SetEventNumber(Int_t ev) { fEventN=ev; }

//protected:
  virtual AliCluster *GetCluster(Int_t index) const=0;
  virtual void  UseClusters(const AliKalmanTrack *t, Int_t from=0) const;
  virtual void  CookLabel(AliKalmanTrack *t,Float_t wrong) const; 
  Double_t GetX() const {return fX;}
  Double_t GetY() const {return fY;}
  Double_t GetZ() const {return fZ;}
  Int_t GetEventNumber() const {return fEventN;}

  static Int_t SetFieldFactor(Char_t* fileName, Bool_t closeFile = kTRUE);
  static Int_t SetFieldFactor(TFile* file, Bool_t deletegAlice = kTRUE);
  static Int_t SetFieldFactor();
  
private:
  Int_t fEventN;//event number

  Double_t fX;  //X-coordinate of the primary vertex
  Double_t fY;  //Y-coordinate of the primary vertex
  Double_t fZ;  //Z-coordinate of the primary vertex

  ClassDef(AliTracker,1) //abstract tracker
};

#endif


