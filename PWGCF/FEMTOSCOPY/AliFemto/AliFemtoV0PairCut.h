///
/// \file AliFemtoV0PairCut.h
/// \date Mar 10, 2008
/// \author Adam Kisiel <kisiel@mps.ohio-state.edu>
///
/// \class AliFemtoV0PairCut
/// \brief A pair cut used for testing two V0 particles
///

#ifndef ALIFEMTOV0PAIRCUT_H
#define ALIFEMTOV0PAIRCUT_H

#include "AliFemtoPairCut.h"

class AliFemtoV0PairCut : public AliFemtoPairCut {
public:
  AliFemtoV0PairCut();
  AliFemtoV0PairCut(const AliFemtoV0PairCut &cut);
  virtual ~AliFemtoV0PairCut();
  AliFemtoV0PairCut &operator=(const AliFemtoV0PairCut &cut);

  virtual bool Pass(const AliFemtoPair *pair);  ///< Checks for pairs that are possibly shared/double reconstruction
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoPairCut *Clone(); ///< Creates a new object with ALL the same attributes as the original
  void SetV0Max(Double_t aAliFemtoV0Max);
  Double_t GetAliFemtoV0Max() const;
  void SetRemoveSameLabel(Bool_t aRemove);
  void SetTPCEntranceSepMinimum(double dtpc);
  void SetTPCExitSepMinimum(double dtpc);
  void SetDataType(AliFemtoDataType type);
  void SetMinAvgSeparation(int type, double minSep);
  void SetNanoAODAnalysis(bool aNanoAOD);

protected:
  long fNPairsPassed;          ///< Number of pairs consideered that passed the cut
  long fNPairsFailed;          ///< Number of pairs consideered that failed the cut
  Double_t fV0Max;             ///< Maximum allowed pair quality
  Double_t fShareFractionMax;  ///< Maximum allowed share fraction
  Bool_t   fRemoveSameLabel;   ///< If 1 pairs with two tracks with the same label will be removed

  AliFemtoDataType fDataType;  ///< Use ESD / AOD / Kinematics.
  Double_t fDTPCMin;           ///< Minimum allowed pair nominal separation at the entrance to the TPC
  Double_t fDTPCExitMin;       ///< Minimum allowed pair nominal separation at the exit of the TPC
  double   fMinAvgSepPosPos;
  double   fMinAvgSepPosNeg;
  double   fMinAvgSepNegPos;
  double   fMinAvgSepNegNeg;
  bool fNanoAODAnalysis;

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoV0PairCut, 0);
  /// \endcond
#endif
};

inline AliFemtoV0PairCut::AliFemtoV0PairCut(const AliFemtoV0PairCut &c):
  AliFemtoPairCut(c),
  fNPairsPassed(c.fNPairsPassed),
  fNPairsFailed(c.fNPairsFailed),
  fV0Max(c.fV0Max),
  fShareFractionMax(c.fShareFractionMax),
  fRemoveSameLabel(c.fRemoveSameLabel),
  fDataType(c.fDataType),
  fDTPCMin(c.fDTPCMin),
  fDTPCExitMin(c.fDTPCExitMin),
  fMinAvgSepPosPos(c.fMinAvgSepPosPos),
  fMinAvgSepPosNeg(c.fMinAvgSepPosNeg),
  fMinAvgSepNegPos(c.fMinAvgSepPosNeg),
  fMinAvgSepNegNeg(c.fMinAvgSepNegNeg),
  fNanoAODAnalysis(c.fNanoAODAnalysis)
{
  /// no-op
}

inline AliFemtoPairCut *AliFemtoV0PairCut::Clone()
{
  AliFemtoV0PairCut *c = new AliFemtoV0PairCut(*this);
  // Should we set fNPairsPassed & fNPairsFailed to 0 here? How will Clone() be used?
  return c;
}
inline void AliFemtoV0PairCut::SetNanoAODAnalysis(bool aNanoAOD) {fNanoAODAnalysis = aNanoAOD;}

#endif
