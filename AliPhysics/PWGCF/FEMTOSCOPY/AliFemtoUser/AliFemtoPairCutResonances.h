/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoPairCutResonances - a pair cut which checks     //
// for some pair qualities that attempt to identify slit/doubly            //
// reconstructed tracks and also selects pairs based on their separation   //
// at the entrance to the TPC                                              //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////


#ifndef ALIFEMTOPAIRCUTRESONANCES_H
#define ALIFEMTOPAIRCUTRESONANCES_H

#include "AliFemtoPairCut.h"
#include "AliFemtoShareQualityPairCut.h"

class AliFemtoPairCutResonances : public AliFemtoShareQualityPairCut{
public:
  AliFemtoPairCutResonances();
  AliFemtoPairCutResonances(const AliFemtoPairCutResonances& c);
  virtual ~AliFemtoPairCutResonances();
  AliFemtoPairCutResonances& operator=(const AliFemtoPairCutResonances& c);

  virtual bool Pass(const AliFemtoPair* pair);
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoPairCut* Clone();
  void SetMaxEEMinv(Double_t maxeeminv);
  void SetMaxThetaDiff(Double_t maxdtheta);
  void SetDataType(AliFemtoDataType type);
  void SetChooseResonances(bool onlyResonances);

 protected:
  Double_t fMaxEEMinv; // Maximum allowed ee Minv
  Double_t fMaxDTheta; // Maximum polar angle difference
  AliFemtoDataType fDataType; //Use ESD / AOD / Kinematics.
  bool fSwitchPassFail; // cut resonances (false), choose resonances (true)

#ifdef __ROOT__
  ClassDef(AliFemtoPairCutResonances, 0)
#endif
};

inline AliFemtoPairCut* AliFemtoPairCutResonances::Clone() { AliFemtoPairCutResonances* c = new AliFemtoPairCutResonances(*this); return c;}

#endif
