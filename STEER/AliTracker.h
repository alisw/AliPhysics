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
#include <TObject.h>

class AliMagF;
class AliCluster;
class TTree;
class AliKalmanTrack;
class AliESD;
class AliTrackPoint;

class AliTracker : public TObject {
public:
  AliTracker();
  virtual ~AliTracker(){}
  virtual Int_t Clusters2Tracks(AliESD *event)=0;
  virtual Int_t PropagateBack(AliESD *event)=0;
  virtual Int_t RefitInward(AliESD *event)=0;
  void SetVertex(const Double_t *xyz, const Double_t *ers=0) { 
     fX=xyz[0]; fY=xyz[1]; fZ=xyz[2];
     if (ers) { fSigmaX=ers[0]; fSigmaY=ers[1]; fSigmaZ=ers[2]; } 
  }

//protected:
  virtual Int_t LoadClusters(TTree *)=0;
  virtual void UnloadClusters()=0;
  virtual AliCluster *GetCluster(Int_t index) const=0;
  virtual Bool_t GetTrackPoint(Int_t /* index */ , AliTrackPoint& /* p */) const { return kFALSE;}
  virtual void  UseClusters(const AliKalmanTrack *t, Int_t from=0) const;
  virtual void  CookLabel(AliKalmanTrack *t,Float_t wrong) const; 
  Double_t GetX() const {return fX;}
  Double_t GetY() const {return fY;}
  Double_t GetZ() const {return fZ;}
  Double_t GetSigmaX() const {return fSigmaX;}
  Double_t GetSigmaY() const {return fSigmaY;}
  Double_t GetSigmaZ() const {return fSigmaZ;}

  static void SetFieldMap(const AliMagF* map, Bool_t uni);
  static const AliMagF *GetFieldMap() {return fgkFieldMap;}
  static Double_t GetBz(Float_t *r); 
  static Double_t GetBz(Double_t *r) {
    Float_t rr[]={r[0],r[1],r[2]};
    return GetBz(rr);
  }
  static Double_t GetBz() {return fgBz;}
  static Bool_t UniformField() {return fgUniformField;}

protected:
  AliTracker(const AliTracker &atr);
private:
  AliTracker & operator=(const AliTracker & atr);

  static Bool_t fgUniformField;       // uniform field flag
  static const AliMagF *fgkFieldMap;  // field map
  static Double_t fgBz;               // Nominal Bz (kG)

  Double_t fX;  //X-coordinate of the primary vertex
  Double_t fY;  //Y-coordinate of the primary vertex
  Double_t fZ;  //Z-coordinate of the primary vertex
 
  Double_t fSigmaX; // error of the primary vertex position in X
  Double_t fSigmaY; // error of the primary vertex position in Y
  Double_t fSigmaZ; // error of the primary vertex position in Z

  ClassDef(AliTracker,3) //abstract tracker
};

#endif


