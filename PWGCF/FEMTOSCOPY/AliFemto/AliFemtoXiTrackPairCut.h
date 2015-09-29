#ifndef ALIFEMTOXITRACKPAIRCUT_H
#define ALIFEMTOXITRACKPAIRCUT_H

#include "AliFemtoPairCut.h"

class AliFemtoXiTrackPairCut : public AliFemtoPairCut {
public:
  AliFemtoXiTrackPairCut();
  AliFemtoXiTrackPairCut(const AliFemtoXiTrackPairCut &cut);
  virtual ~AliFemtoXiTrackPairCut();
  AliFemtoXiTrackPairCut &operator=(const AliFemtoXiTrackPairCut &cut);

  virtual bool Pass(const AliFemtoPair *pair);  ///< Checks for pairs that are possibly shared/double reconstruction
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoPairCut *Clone(); ///< Creates a new object with ALL the same attributes as the original
  void SetDataType(AliFemtoDataType type);
  void SetTPCOnly(Bool_t tpconly);

protected:
  long fNPairsPassed;          ///< Number of pairs consideered that passed the cut
  long fNPairsFailed;          ///< Number of pairs consideered that failed the cut

  AliFemtoDataType fDataType;  ///< Use ESD / AOD / Kinematics.
  Bool_t fTrackTPCOnly;       ///< Track ids will be converted to appropriate TPC/Global value 

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoXiTrackPairCut, 0)
  /// \endcond
#endif
};
 
inline AliFemtoXiTrackPairCut::AliFemtoXiTrackPairCut(const AliFemtoXiTrackPairCut &c):
  AliFemtoPairCut(c),
  fNPairsPassed(c.fNPairsPassed),
  fNPairsFailed(c.fNPairsFailed),
  fDataType(c.fDataType),
  fTrackTPCOnly(c.fTrackTPCOnly)
{
  /// no-op
}

inline AliFemtoPairCut *AliFemtoXiTrackPairCut::Clone()
{
  AliFemtoXiTrackPairCut *c = new AliFemtoXiTrackPairCut(*this);
  // Should we set fNPairsPassed & fNPairsFailed to 0 here? How will Clone() be used?
  return c;
}

#endif
