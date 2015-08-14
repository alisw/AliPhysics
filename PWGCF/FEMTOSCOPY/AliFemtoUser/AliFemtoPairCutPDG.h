/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
// AliFemtoPairCutPDG - a pair cut which checks if the sum of the transverse       //
// momenta of two particles fit within given range and provides cut for the        //
// pair of particles originated from the same mother particle                      //
// Authors: Malgorzata Janik, Warsaw University of Technology, majanik@cern.ch     //
//          Lukasz Graczykowski, Warsaw University of Technology, lgraczyk@cern.ch //
//          Piotr Modzelewski pmodzele@cern.ch                                     //
//  	       				                                                           //
/////////////////////////////////////////////////////////////////////////////////////


#ifndef ALIFEMTOPAIRCUTPDG_H
#define ALIFEMTOPAIRCUTPDG_H

#include "AliFemtoPairCut.h"

class AliFemtoPairCutPDG : public AliFemtoPairCut{
public:
  AliFemtoPairCutPDG();
  AliFemtoPairCutPDG(double lo, double hi, int pdg_1, int pdg_2);
  AliFemtoPairCutPDG(const AliFemtoPairCutPDG& c);
  virtual ~AliFemtoPairCutPDG();
  AliFemtoPairCutPDG& operator=(const AliFemtoPairCutPDG& c);

  virtual bool Pass(const AliFemtoPair* pair);
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  AliFemtoPairCut* Clone();
  void SetMinSumPt(Double_t sumptmin);
  void SetMaxSumPt(Double_t sumptmax);
  void SetPDG1(Double_t pdg1);
  void SetPDG2(Double_t pdg2);


 protected:
  Double_t fSumPtMin;
  Double_t fSumPtMax;
  Double_t fPDG1;
  Double_t fPDG2;
  Double_t fNPairsFailed;
  Double_t fNPairsPassed;

#ifdef __ROOT__
  ClassDef(AliFemtoPairCutPDG, 0)
#endif
};

inline AliFemtoPairCut* AliFemtoPairCutPDG::Clone() { AliFemtoPairCutPDG* c = new AliFemtoPairCutPDG(*this); return c;}

#endif
