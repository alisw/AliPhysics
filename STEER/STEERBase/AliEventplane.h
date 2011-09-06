#ifndef ALIEVENTPLANE_H
#define ALIEVENTPLANE_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*****************************************************
//   Class AliEventplane
//   author: Alberica Toia, Johanna Gramling
//*****************************************************

#include "TNamed.h"

class TVector2;
class AliVTrack;
class TObjArray;
class TArrayF;

class AliEventplane : public TNamed
{
 public:

  AliEventplane();  /// constructor
  ~AliEventplane();  /// destructor
  AliEventplane(const AliEventplane& ep); /// copy constructor
  AliEventplane& operator=(const AliEventplane& ep);   /// assignment operator
  virtual void CopyEP(AliEventplane& ep) const;

  /// set event plane result
  void SetQVector(TVector2* qvector) {fQVector = qvector;}
  void SetEventplaneQ(Double_t evp) {fEventplaneQ = evp;} 
  void SetQsub(TVector2* qs1, TVector2* qs2) {fQsub1 = qs1;fQsub2 = qs2;}
  void SetQsubRes (Double_t qsr) {fQsubRes = qsr;}

  /// get event plane result
  TVector2* GetQVector(); 
  Double_t  GetQContributionX(AliVTrack* track);
  Double_t  GetQContributionY(AliVTrack* track);
  Double_t  GetQContributionXsub1(AliVTrack* track);
  Double_t  GetQContributionYsub1(AliVTrack* track);
  Double_t  GetQContributionXsub2(AliVTrack* track);
  Double_t  GetQContributionYsub2(AliVTrack* track);
  TArrayF*  GetQContributionXArray() { return fQContributionX; }
  TArrayF*  GetQContributionYArray() { return fQContributionY; }
  TArrayF*  GetQContributionXArraysub1() { return fQContributionXsub1; }
  TArrayF*  GetQContributionYArraysub1() { return fQContributionYsub1; }
  TArrayF*  GetQContributionXArraysub2() { return fQContributionXsub2; }
  TArrayF*  GetQContributionYArraysub2() { return fQContributionYsub2; }
  Double_t  GetEventplane(const char *method);
  TVector2* GetQsub1();
  TVector2* GetQsub2();
  Double_t  GetQsubRes();
  Bool_t    IsEventInEventplaneClass(Double_t a, Double_t b, const char *method);

  void Reset();

 private:
   TVector2* fQVector;		 // Q-Vector of event
   TArrayF* fQContributionX;	 // array of the tracks' contributions to X component of Q-Vector - index = track ID
   TArrayF* fQContributionY;	 // array of the tracks' contributions to Y component of Q-Vector - index = track ID
   TArrayF* fQContributionXsub1; // array of the tracks' contributions to X component of Q-Vectorsub1 - index = track ID
   TArrayF* fQContributionYsub1; // array of the tracks' contributions to Y component of Q-Vectorsub1 - index = track ID
   TArrayF* fQContributionXsub2; // array of the tracks' contributions to X component of Q-Vectorsub2 - index = track ID
   TArrayF* fQContributionYsub2; // array of the tracks' contributions to Y component of Q-Vectorsub2 - index = track ID
   Double_t fEventplaneQ;	 // Event plane angle from Q-Vector
   TVector2* fQsub1;		 // Q-Vector of subevent 1
   TVector2* fQsub2;		 // Q-Vector of subevent 2
   Double_t fQsubRes;		 // Difference of EP angles of subevents
 
  ClassDef(AliEventplane, 1)
};
#endif //ALIEVENTPLANE_H
