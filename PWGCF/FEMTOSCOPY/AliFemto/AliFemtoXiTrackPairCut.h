#ifndef ALIFEMTOXITRACKPAIRCUT_H
#define ALIFEMTOXITRACKPAIRCUT_H

#include "AliFemtoPairCut.h"
#include "AliFemtoV0TrackPairCut.h"

/// \class AliFemtoXiTrackPairCut
/// \brief A pair cut used for testing a track and a reconstructed Xi particle
/// Note:  It will be useful for this class to have access to the methods available
///        in AliFemtoV0TrackPair cut.  However, at this point, this class will continue
///        to derive from AliFemtoPair cut, and will contain a AliFemtoV0TrackPairCut object.
///        In the future, this idea may want to be revisited, and instead having this class
///        derive directly from AliFemtoV0TrackPairCut may be ideal.

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

  AliFemtoV0TrackPairCut* GetV0TrackPairCut(); //allows one to set fV0TrackPairCut attributes, so no need to explicitly state here
  void SetMinAvgSepTrackBacPion(double aMin);

protected:
  AliFemtoV0TrackPairCut* fV0TrackPairCut;

  long fNPairsPassed;          ///< Number of pairs consideered that passed the cut
  long fNPairsFailed;          ///< Number of pairs consideered that failed the cut

  double fMinAvgSepTrackBacPion;

  AliFemtoDataType fDataType;  ///< Use ESD / AOD / Kinematics.
  Bool_t fTrackTPCOnly;       ///< Track ids will be converted to appropriate TPC/Global value

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoXiTrackPairCut, 0);
  /// \endcond
#endif
};

inline AliFemtoXiTrackPairCut::AliFemtoXiTrackPairCut(const AliFemtoXiTrackPairCut &c):
  AliFemtoPairCut(c),
  fV0TrackPairCut(c.fV0TrackPairCut ? new AliFemtoV0TrackPairCut(*c.fV0TrackPairCut) : nullptr),
  fNPairsPassed(c.fNPairsPassed),
  fNPairsFailed(c.fNPairsFailed),
  fMinAvgSepTrackBacPion(c.fMinAvgSepTrackBacPion),
  fDataType(c.fDataType),
  fTrackTPCOnly(c.fTrackTPCOnly)
{

}

inline AliFemtoPairCut *AliFemtoXiTrackPairCut::Clone()
{
  AliFemtoXiTrackPairCut *c = new AliFemtoXiTrackPairCut(*this);
  // Should we set fNPairsPassed & fNPairsFailed to 0 here? How will Clone() be used?
  return c;
}

inline AliFemtoV0TrackPairCut* AliFemtoXiTrackPairCut::GetV0TrackPairCut() {return fV0TrackPairCut;}
inline void AliFemtoXiTrackPairCut::SetMinAvgSepTrackBacPion(double aMin) {fMinAvgSepTrackBacPion = aMin;}

#endif
