/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
// AliFemtoPairCutPt - a pair cut which checks if the sum of the transverse        //
// momenta of two particles fit within given range 		                   //
// Authors: Malgorzata Janik, Warsaw University of Technology, majanik@cern.ch     //
//          Lukasz Graczykowski, Warsaw University of Technology, lgraczyk@cern.ch //
//  	       				                                           //
/////////////////////////////////////////////////////////////////////////////////////


#ifndef ALIFEMTOPAIRCUTPT_H
#define ALIFEMTOPAIRCUTPT_H

#include "AliFemtoPairCut.h"

class AliFemtoPairCutPt : public AliFemtoPairCut{
public:
  AliFemtoPairCutPt();
  AliFemtoPairCutPt(double lo, double hi);
  AliFemtoPairCutPt(const AliFemtoPairCutPt& c);
  virtual ~AliFemtoPairCutPt();

  virtual bool Pass(const AliFemtoPair* pair);
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  AliFemtoPairCut* Clone();
  void SetMinSumPt(Double_t sumptmin);
  void SetMaxSumPt(Double_t sumptmax);
  
 protected:
  Double_t fSumPtMin;
  Double_t fSumPtMax;
  Double_t fNPairsFailed;
  Double_t fNPairsPassed;

#ifdef __ROOT__
  ClassDef(AliFemtoPairCutPt, 0)
#endif
};

inline AliFemtoPairCut* AliFemtoPairCutPt::Clone() { AliFemtoPairCutPt* c = new AliFemtoPairCutPt(*this); return c;}

#endif
