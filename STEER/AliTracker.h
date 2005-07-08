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

class AliCluster;
class TTree;
class AliKalmanTrack;
class AliESD;
class AliMagF;

class AliTracker : public TObject {
public:

  enum {kTrackInward, kTrackBack, kTrackRefit} Propagation_t; //type of propagation
  
  AliTracker();
  AliTracker(const AliTracker &atr): TObject(atr)
    {Fatal("Copy ctor","Not Implemented!\n");}
  // AliTracker & operator=(const AliTracker &)
  //  {Fatal("= op","Not Implemented\n");return *this;}
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

private:

  AliTracker & operator=(const AliTracker & atr);

  static const AliMagF *fgkFieldMap; //field map

  Double_t fX;  //X-coordinate of the primary vertex
  Double_t fY;  //Y-coordinate of the primary vertex
  Double_t fZ;  //Z-coordinate of the primary vertex

  Double_t fSigmaX; // error of the primary vertex position in X
  Double_t fSigmaY; // error of the primary vertex position in Y
  Double_t fSigmaZ; // error of the primary vertex position in Z

  ClassDef(AliTracker,2) //abstract tracker
};

#endif


