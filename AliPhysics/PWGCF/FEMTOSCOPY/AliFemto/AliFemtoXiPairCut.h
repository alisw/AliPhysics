#ifndef ALIFEMTOXIPAIRCUT_H
#define ALIFEMTOXIPAIRCUT_H

#include "AliFemtoPairCut.h"

class AliFemtoXiPairCut : public AliFemtoPairCut {
public:
  AliFemtoXiPairCut();
  AliFemtoXiPairCut(const AliFemtoXiPairCut &cut);
  virtual ~AliFemtoXiPairCut();
  AliFemtoXiPairCut &operator=(const AliFemtoXiPairCut &cut);

  virtual bool Pass(const AliFemtoPair *pair);  ///< Checks for pairs that are possibly shared/double reconstruction
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoPairCut *Clone(); ///< Creates a new object with ALL the same attributes as the original
  void SetDataType(AliFemtoDataType type);

protected:
  long fNPairsPassed;          ///< Number of pairs consideered that passed the cut
  long fNPairsFailed;          ///< Number of pairs consideered that failed the cut

  AliFemtoDataType fDataType;  ///< Use ESD / AOD / Kinematics.

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoXiPairCut, 0);
  /// \endcond
#endif
};

inline AliFemtoXiPairCut::AliFemtoXiPairCut(const AliFemtoXiPairCut &c):
  AliFemtoPairCut(c),
  fNPairsPassed(c.fNPairsPassed),
  fNPairsFailed(c.fNPairsFailed),
  fDataType(c.fDataType)
{
  // no-op
}

inline AliFemtoPairCut *AliFemtoXiPairCut::Clone()
{
  AliFemtoXiPairCut *c = new AliFemtoXiPairCut(*this);
  // Should we set fNPairsPassed & fNPairsFailed to 0 here? How will Clone() be used?
  return c;
}

#endif
