/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
// AliFemtoPairCutMInv - a pair cut which checks if the sum of the                 //
//                       invariant mass of two particles fit within given range 	  //
// Authors: Malgorzata Janik, Warsaw University of Technology, majanik@cern.ch     //
//          Lukasz Graczykowski, Warsaw University of Technology, lgraczyk@cern.ch //
//          Piotr Modzelewski, Warsaw University of Technology, pmodzele@cern.ch   //
//  	       				                                                           //
/////////////////////////////////////////////////////////////////////////////////////


#ifndef ALIFEMTOPAIRCUTMINV_H
#define ALIFEMTOPAIRCUTMINV_H

#include "AliFemtoPairCut.h"

class AliFemtoPairCutMInv : public AliFemtoPairCut{
public:
  AliFemtoPairCutMInv();
  AliFemtoPairCutMInv(double m1, double m2, double minvmin, double minvmax);
  AliFemtoPairCutMInv(const AliFemtoPairCutMInv& c);
  virtual ~AliFemtoPairCutMInv();
  AliFemtoPairCutMInv& operator=(const AliFemtoPairCutMInv& c);

  virtual bool Pass(const AliFemtoPair* pair);
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  AliFemtoPairCut* Clone() const;

 protected:
  Double_t fNPairsFailed;
  Double_t fNPairsPassed;
  Double_t fMInvMin;          // Minimum allowed pair invariant mass
  Double_t fMInvMax;          // Maximum allowed pair invariant mass
  Double_t fM1;               // First particle mass
  Double_t fM2;               // Second particle mass

#ifdef __ROOT__
  ClassDef(AliFemtoPairCutMInv, 0)
#endif
};

inline AliFemtoPairCut* AliFemtoPairCutMInv::Clone() const
  { AliFemtoPairCutMInv* c = new AliFemtoPairCutMInv(*this); return c;}

#endif
