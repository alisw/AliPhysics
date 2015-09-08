///
/// \file AliFemtoV0TrackPairCut.h
///
/// \class AliFemtoV0TrackPairCut
/// \brief A pair cut used for testing a track and a reconstructed V0 particle
///

#ifndef ALIFEMTOV0TRACKPAIRCUT_H
#define ALIFEMTOV0TRACKPAIRCUT_H

#include "AliFemtoPairCut.h"

class AliFemtoV0TrackPairCut : public AliFemtoPairCut {
public:
  enum ParticleType {kLambda = 0, kAntiLambda = 1, kProton = 2, kAntiProton = 3};
  typedef enum ParticleType AliFemtoParticleType;
  AliFemtoV0TrackPairCut();
  AliFemtoV0TrackPairCut(const AliFemtoV0TrackPairCut &cut);
  virtual ~AliFemtoV0TrackPairCut();
  AliFemtoV0TrackPairCut &operator=(const AliFemtoV0TrackPairCut &cut);

  virtual bool Pass(const AliFemtoPair *pair);
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoPairCut *Clone();
  void SetV0Max(Double_t aAliFemtoV0Max);
  Double_t GetAliFemtoV0Max() const;
  void SetRemoveSameLabel(Bool_t aRemove);
  void SetTPCOnly(Bool_t tpconly);
  void SetShareQualityMax(Double_t aShareQualityMax);
  void SetShareFractionMax(Double_t aShareFractionMax);
  void SetTPCEntranceSepMinimum(double dtpc);
  void SetTPCExitSepMinimum(double dtpc);
  void SetDataType(AliFemtoDataType type);
  void SetKstarCut(double kstar, AliFemtoParticleType firstParticle, AliFemtoParticleType secondParticle);
  void SetMinAvgSeparation(int type, double minSep);

protected:
  long fNPairsPassed;  ///< Number of pairs considered that passed the cut
  long fNPairsFailed;  ///< Number of pairs considered that failed the cut

private:
  Double_t fV0Max;            ///< Maximum allowed pair quality
  Double_t fShareQualityMax;  ///< Maximum allowed share quality
  Double_t fShareFractionMax; ///< Maximum allowed share fraction
  Bool_t fRemoveSameLabel;    ///< If 1 pairs with two tracks with the same label will be removed
  Bool_t fTrackTPCOnly;       ///< Track ids will be converted to appropriate TPC/Global value 

  AliFemtoDataType fDataType; ///< Use ESD / AOD / Kinematics.
  Double_t fDTPCMin;          ///< Minimum allowed pair nominal separation at the entrance to the TPC
  Double_t fDTPCExitMin;      ///< Minimum allowed pair nominal separation at the exit of the TPC

  double fKstarCut;                         ///< do we want the K star cut, if yes (>0) then it is the minimum value of k*
  AliFemtoParticleType fFirstParticleType;  ///< for kstar - first particle type (V0 type)
  AliFemtoParticleType fSecondParticleType; ///< for kstar - second particle type (primary track)
  double fMinAvgSepTrackPos;
  double fMinAvgSepTrackNeg;

#ifdef __ROOT__
  ClassDef(AliFemtoV0TrackPairCut, 0)
#endif
};

inline AliFemtoV0TrackPairCut::AliFemtoV0TrackPairCut(const AliFemtoV0TrackPairCut &c):
  AliFemtoPairCut(c),
  fNPairsPassed(c.fNPairsPassed),
  fNPairsFailed(c.fNPairsFailed),
  fV0Max(c.fV0Max),
  fShareQualityMax(c.fShareQualityMax),
  fShareFractionMax(c.fShareFractionMax),
  fRemoveSameLabel(c.fRemoveSameLabel),
  fTrackTPCOnly(c.fTrackTPCOnly),
  fDataType(c.fDataType),
  fDTPCMin(c.fDTPCMin),
  fDTPCExitMin(c.fDTPCExitMin),
  fKstarCut(c.fKstarCut),
  fFirstParticleType(c.fFirstParticleType),
  fSecondParticleType(c.fSecondParticleType),
  fMinAvgSepTrackPos(c.fMinAvgSepTrackPos),
  fMinAvgSepTrackNeg(c.fMinAvgSepTrackNeg)
{
  /* no-op */
}

inline AliFemtoPairCut *AliFemtoV0TrackPairCut::Clone()
{
  AliFemtoV0TrackPairCut *c = new AliFemtoV0TrackPairCut(*this);
  return c;
}

#endif

