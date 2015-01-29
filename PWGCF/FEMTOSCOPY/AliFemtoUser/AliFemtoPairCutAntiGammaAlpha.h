/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// AliFemtoPairCutAntiGamma - a pair cut which checks     //
// for some pair qualities that attempt to identify slit/doubly            //
// reconstructed tracks and also selects pairs based on their separation   //
// at the entrance to the TPC                                              //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
/***************************************************************************
 *
 * $Id: AliFemtoPairCutAntiGamma.h,v 1.1.2.1 2007/10/19 13:35:33 akisiel Exp $
 *
 * Author: Adam Kisiel, Ohio State University, kisiel@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *   a cut to remove "shared" and "split" pairs
 *
 ***************************************************************************
 *
 *
 **************************************************************************/


#ifndef ALIFEMTOPAIRCUTANTIGAMMAALPHA_H
#define ALIFEMTOPAIRCUTANTIGAMMAALPHA_H

#include "AliFemtoPairCut.h"
#include "AliFemtoShareQualityPairCut.h"

class AliFemtoPairCutAntiGammaAlpha : public AliFemtoShareQualityPairCut{
public:
  AliFemtoPairCutAntiGammaAlpha();
  AliFemtoPairCutAntiGammaAlpha(const AliFemtoPairCutAntiGammaAlpha& c);
  virtual ~AliFemtoPairCutAntiGammaAlpha();
  AliFemtoPairCutAntiGammaAlpha& operator=(const AliFemtoPairCutAntiGammaAlpha& c);

  virtual bool Pass(const AliFemtoPair* pair);
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoPairCut* Clone();
  void SetMaxEEMinv(Double_t maxeeminv);
  void SetMaxAlphaDiff(Double_t maxdalpha);
  void SetTPCEntranceSepMinimum(double dtpc);
  /* void SetTPCExitSepMinimum(double dtpc); */
  void SetDataType(AliFemtoDataType type);

 protected:
  Double_t fMaxEEMinv; // Maximum allowed ee Minv
  Double_t fMaxDAlpha; // Maximum polar angle difference
  Double_t fDTPCMin;          // Minimum allowed pair nominal separation at the entrance to the TPC
  AliFemtoDataType fDataType; //Use ESD / AOD / Kinematics.

#ifdef __ROOT__
  ClassDef(AliFemtoPairCutAntiGammaAlpha, 0)
#endif
};

inline AliFemtoPairCut* AliFemtoPairCutAntiGammaAlpha::Clone() { AliFemtoPairCutAntiGammaAlpha* c = new AliFemtoPairCutAntiGammaAlpha(*this); return c;}

#endif
