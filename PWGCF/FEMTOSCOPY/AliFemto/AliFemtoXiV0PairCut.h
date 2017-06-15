#ifndef ALIFEMTOXIV0PAIRCUT_H
#define ALIFEMTOXIV0PAIRCUT_H

#include "AliFemtoPairCut.h"
#include "AliFemtoV0PairCut.h"

class AliFemtoXiV0PairCut : public AliFemtoPairCut {
public:
  AliFemtoXiV0PairCut();
  AliFemtoXiV0PairCut(const AliFemtoXiV0PairCut &cut);
  virtual ~AliFemtoXiV0PairCut();
  AliFemtoXiV0PairCut &operator=(const AliFemtoXiV0PairCut &cut);

  virtual bool Pass(const AliFemtoPair *pair);  ///< Checks for pairs that are possibly shared/double reconstruction
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoPairCut *Clone(); ///< Creates a new object with ALL the same attributes as the original
  void SetDataType(AliFemtoDataType type);

  AliFemtoV0PairCut* GetV0PairCut(); //allows one to set fV0PairCut attributes, so no need to explicitly state here
  void SetMinAvgSepBacPos(double aMin);
  void SetMinAvgSepBacNeg(double aMin);

protected:
  AliFemtoV0PairCut* fV0PairCut;

  long fNPairsPassed;          ///< Number of pairs consideered that passed the cut
  long fNPairsFailed;          ///< Number of pairs consideered that failed the cut

  AliFemtoDataType fDataType;  ///< Use ESD / AOD / Kinematics.

  double fMinAvgSepBacPos, fMinAvgSepBacNeg;

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoXiV0PairCut, 0);
  /// \endcond
#endif
};

inline AliFemtoXiV0PairCut::AliFemtoXiV0PairCut(const AliFemtoXiV0PairCut &c):
  AliFemtoPairCut(c),
  fV0PairCut(c.fV0PairCut ? new AliFemtoV0PairCut(*c.fV0PairCut) : nullptr),
  fNPairsPassed(c.fNPairsPassed),
  fNPairsFailed(c.fNPairsFailed),
  fDataType(c.fDataType),
  fMinAvgSepBacPos(c.fMinAvgSepBacPos),
  fMinAvgSepBacNeg(c.fMinAvgSepBacNeg)

{
  /// no-op
}

inline AliFemtoPairCut *AliFemtoXiV0PairCut::Clone()
{
  AliFemtoXiV0PairCut *c = new AliFemtoXiV0PairCut(*this);
  // Should we set fNPairsPassed & fNPairsFailed to 0 here? How will Clone() be used?
  return c;
}

inline AliFemtoV0PairCut* AliFemtoXiV0PairCut::GetV0PairCut() {return fV0PairCut;}
inline void AliFemtoXiV0PairCut::SetMinAvgSepBacPos(double aMin) {fMinAvgSepBacPos = aMin;}
inline void AliFemtoXiV0PairCut::SetMinAvgSepBacNeg(double aMin) {fMinAvgSepBacNeg = aMin;}

#endif
