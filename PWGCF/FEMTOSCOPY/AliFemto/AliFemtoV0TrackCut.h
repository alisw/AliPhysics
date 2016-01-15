/// \class AliFemtoV0TrackCut

#ifndef ALIFEMTOV0TRACKCUT_H
#define ALIFEMTOV0TRACKCUT_H

#include "AliFemtoTrackCut.h"

/**
 * \class AliFemtoV0TrackCut
 * \brief A track cut designed to cut on V0 (i.e. Lambda) particles.
 *
 * A particle cut object which tests AliFemtoV0 obects for a number of
 * conditions for the reconstructed V0 and the positive and negative
 * (measured) daughter tracks.
 */
class AliFemtoV0TrackCut : public AliFemtoParticleCut {
public:
  /// Enumerated type to easily select correct algorithm for each particle type
  enum V0Type {
    kLambda = 0,
    kAntiLambda = 1,
    kK0s = 2,
    kAll = 99,
    kLambdaMC = 101,
    kAntiLambdaMC = 102,
    kK0sMC = 3
  };
  typedef enum V0Type AliFemtoV0Type;

  AliFemtoV0TrackCut();
  virtual ~AliFemtoV0TrackCut();

  virtual bool Pass(const AliFemtoV0* aV0);

  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoParticleType Type(){return hbtV0;}

  void SetInvariantMassLambda(double,double);
  void SetInvariantMassK0s(double,double);
  void SetMinDaughtersToPrimVertex(double,double);
  void SetMaxDcaV0Daughters(double);
  void SetMaxDcaV0(double);
  void SetMinDcaV0(double);
  void SetMaxCosPointingAngle(double);
  void SetMinCosPointingAngle(double);
  void SetMaxV0DecayLength(double);
  void SetParticleType(short);
  void SetEta(double);
  void SetPt(double,double);
  void SetEtaDaughters(float);
  void SetTPCnclsDaughters(int);
  void SetNdofDaughters(int);
  void SetStatusDaughters(unsigned long);
  void SetPtPosDaughter(float,float);
  void SetPtNegDaughter(float,float);
  void SetOnFlyStatus(bool);
  void SetMinAvgSeparation(double);

  //----n sigma----
  bool IsKaonTPCdEdxNSigma(float mom, float nsigmaK);
  bool IsKaonTOFNSigma(float mom, float nsigmaK);
  bool IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK);
  bool IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi);
  bool IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP);

 protected:   // here are the quantities I want to cut on...

  double fInvMassLambdaMin;        ///< invariant mass lambda min
  double fInvMassLambdaMax;        ///< invariant mass lambda max
  double fInvMassK0sMin;           ///< invariant mass lambda min
  double fInvMassK0sMax;           ///< invariant mass lambda max
  double fMinDcaDaughterPosToVert; ///< DCA of positive daughter to primary vertex
  double fMinDcaDaughterNegToVert; ///< DCA of negative daughter to primary vertex
  double fMaxDcaV0Daughters;       ///< Max DCA of v0 daughters at Decay vertex
  double fMaxDcaV0;
  double fMinDcaV0;
  double fMaxDecayLength;

  double fMaxCosPointingAngle;    //obsolete
  double fMinCosPointingAngle;    //correct
  short fParticleType;             ///< 0-lambda
  double fEta;
  double fPtMin;
  double fPtMax;
  bool fOnFlyStatus;

  float fMaxEtaDaughters;         ///< Eta of positive daughter
  int fTPCNclsDaughters;          ///< No. of cls of pos daughter
  int fNdofDaughters;             ///< No. of degrees of freedom of the pos. daughter track
  unsigned long fStatusDaughters; ///< Status (tpc refit, its refit...)
  float fPtMinPosDaughter;
  float fPtMaxPosDaughter;
  float fPtMinNegDaughter;
  float fPtMaxNegDaughter;
  double fMinAvgSepDaughters;

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoV0TrackCut, 1);
  /// \endcond
#endif

};


#endif
