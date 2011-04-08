#ifndef ALIEventplane_H
#define ALIEventplane_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*****************************************************
//   Class AliEventplane
//   author: Alberica Toia, Johanna Gramling
//*****************************************************

#include "TNamed.h"

class TVector2;
class AliESDtrack;
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
  Double_t GetQContributionX(AliESDtrack* track);
  Double_t GetQContributionY(AliESDtrack* track);
  TArrayF* GetQContributionXArray() { return fQContributionX; }
  TArrayF* GetQContributionYArray() { return fQContributionY; }
  Double_t GetEventplane(const char *method);
  TVector2* GetQsub1();
  TVector2* GetQsub2();
  Double_t GetQsubRes();
  Bool_t  IsEventInEventplaneClass(Double_t a, Double_t b, const char *method);

 private:
   TVector2* fQVector;		// Q-Vector of event
   TArrayF* fQContributionX;	// array of the tracks' contributions to X component of Q-Vector - index = track ID
   TArrayF* fQContributionY;	// array of the tracks' contributions to Y component of Q-Vector - index = track ID
   Double_t fEventplaneQ;	// Event plane angle from Q-Vector
   TVector2* fQsub1;		// Q-Vector of subevent 1
   TVector2* fQsub2;		// Q-Vector of subevent 2
   Double_t fQsubRes;		// Difference of EP angles of subevents
 
  ClassDef(AliEventplane, 1)
};
#endif //ALIEVENTPLANE_H
