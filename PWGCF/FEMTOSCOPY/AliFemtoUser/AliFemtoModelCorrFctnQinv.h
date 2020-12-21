///
/// \file AliFemtoModelCorrFctnQinv.h
/// \author Andrew Kubera
///

#pragma once

#ifndef ALIFEMTOMODELCORRFCTNQVAL_H
#define ALIFEMTOMODELCORRFCTNQVAL_H

class TH1F;
class TH2F;

#include "AliFemtoModelCorrFctn.h"
#include "AliFemtoAvgSepCorrFctn.h"
#include "AliFemtoPair.h"
#include "AliFemtoThreeVector.h"
#include "AliFemtoLorentzVector.h"
#include "AliFemtoModelHiddenInfo.h"

#include <cmath>

/// \class AliFemtoModelCorrFctnQinv
/// \brief The correlation function which creates numerator and denominator
///        plots from the true Monte Carlo tracks
///
/// \author: Andrew Kubera, andrew.kubera@cern.ch
///
///
class AliFemtoModelCorrFctnQinv : public AliFemtoModelCorrFctn {
public:

  /// Default constructor
  ///
  /// Histograms are created with default titles/ranges. No other members are
  /// initialized - user **MUST** set pair type and PDG codes.
  ///
  AliFemtoModelCorrFctnQinv();

  /// Construct with histogram parameters
  ///
  /// KStar histograms use these parameters when constructing histograms. No
  /// other memebers are set - user **MUST** set pair type and PDG codes.
  ///
  AliFemtoModelCorrFctnQinv(const char *suffix,
                            const Int_t aNbins,
                            const Float_t aQinvLo,
                            const Float_t aQinvHi);

  /// Copy Constructor
  ///
  /// Copies pair type & PDG codes. Clones histograms.
  ///
  AliFemtoModelCorrFctnQinv(const AliFemtoModelCorrFctnQinv &);

  /// assignment operator
  AliFemtoModelCorrFctnQinv& operator=(const AliFemtoModelCorrFctnQinv &);

  /// Return a pointer to a clone of this correlation function.
  ///
  /// This copies the structure of the internal histograms, but clears
  /// all stored data.
  /// No pointers are shared between this object and the clone.
  ///
  virtual AliFemtoModelCorrFctn* Clone() const;


  /// Destructor
  virtual ~AliFemtoModelCorrFctnQinv();

  /// Return information about the run of the correlation function
  virtual AliFemtoString Report();

  virtual void AddRealPair(AliFemtoPair*);
  virtual void AddMixedPair(AliFemtoPair*);

  virtual TList* GetOutputList();

  virtual TList* AppendOutputList(TList*) const;

  /// Sets whether to expect pairs of Track+Track, Track+V0, or V0+V0 particles
  virtual void SetPairType(const AliFemtoAvgSepCorrFctn::PairType pairtype);

  /// Sets the PDG particle codes that this correlation function will check
  /// particle true.
  virtual void SetExpectedPDGCodes(const Int_t track1_pdg_code,
                                   const Int_t track2_pdg_code);

  /// Returns true if the pair's particles match the ExpectedPDGCodes provided
  virtual bool PairContainsExpectedTypes(const AliFemtoPair*);
  virtual bool PairContainsExpectedTypes(const AliFemtoModelHiddenInfo*,
                                         const AliFemtoModelHiddenInfo*);

  /// Static method which calculates kstar based on two lortenz vectors
  static double CalcQinv(const AliFemtoLorentzVector&,
                         const AliFemtoLorentzVector&);


  /// Calculate kstar of the pair - implementation chosen by the fPairType
  /// member (Track+Track, Track+V0, etc)
  virtual double CalcTrueQinv(const AliFemtoPair *);

  /// Calculate q_{inv} for an AliFemtoPair containing a track and a V0.
  /// The order does not matter
  // virtual double CalcTrueKStarTrackV0(const AliFemtoPair *);

  /// Returns the Track and V0 objects in the pair
  static void GetTrackV0(const AliFemtoPair*,
                         AliFemtoTrack*&,
                         AliFemtoV0*&);

protected:

  /// Adds pair information to appropriate histogram points
  void AddPair(const AliFemtoPair&, TH1F*, TH2F*);

  /// The particle class combination types (track, V0, etc...)
  AliFemtoAvgSepCorrFctn::PairType fPairType;

  Int_t fExpectedTrack1Code;
  Int_t fExpectedTrack2Code;

  /// Numerator with axis for non-expected particle PID
  TH2F *fNumPid;

  /// Denominator with axis for non-expected particle PID
  TH2F *fDenPid;
};


inline void
AliFemtoModelCorrFctnQinv::SetPairType(const AliFemtoAvgSepCorrFctn::PairType pairtype)
{
  fPairType = pairtype;
}

inline void
AliFemtoModelCorrFctnQinv::SetExpectedPDGCodes(const Int_t track1_pdg_code,
                                               const Int_t track2_pdg_code)
{
  fExpectedTrack1Code = track1_pdg_code;
  fExpectedTrack2Code = track2_pdg_code;
}

inline void
AliFemtoModelCorrFctnQinv::GetTrackV0(const AliFemtoPair* pair,
                                      AliFemtoTrack*& track_return,
                                      AliFemtoV0*& v0_return)
{
  v0_return = pair->Track1()->V0();

  // this assumption was correct, get track from track2
  if (v0_return != NULL) {
    track_return = pair->Track2()->Track();

  // V0 must be in track 2
  } else {
    track_return = pair->Track1()->Track();
    v0_return = pair->Track2()->V0();
  }
}

inline double
AliFemtoModelCorrFctnQinv::CalcTrueQinv(const AliFemtoPair *pair)
{
  const AliFemtoModelHiddenInfo *info1 = NULL,
                                *info2 = NULL;
  switch (fPairType) {
  case AliFemtoAvgSepCorrFctn::kTracks:
    info1 = (AliFemtoModelHiddenInfo*)pair->Track1()->Track()->GetHiddenInfo();
    info2 = (AliFemtoModelHiddenInfo*)pair->Track2()->Track()->GetHiddenInfo();
    // info1 = dynamic_cast<const AliFemtoModelHiddenInfo*>(pair->Track1()->Track()->GetHiddenInfo());
    // info2 = dynamic_cast<const AliFemtoModelHiddenInfo*>(pair->Track2()->Track()->GetHiddenInfo());
    break;

  case AliFemtoAvgSepCorrFctn::kV0s:
    info1 = dynamic_cast<const AliFemtoModelHiddenInfo*>(pair->Track1()->V0()->GetHiddenInfo());
    info2 = dynamic_cast<const AliFemtoModelHiddenInfo*>(pair->Track2()->V0()->GetHiddenInfo());
    break;

  case AliFemtoAvgSepCorrFctn::kTrackV0:
    {
    AliFemtoTrack *track = NULL;
    AliFemtoV0 *v0 = NULL;
    GetTrackV0(pair, track, v0);
    info1 = dynamic_cast<const AliFemtoModelHiddenInfo*>(track->GetHiddenInfo());
    info2 = dynamic_cast<const AliFemtoModelHiddenInfo*>(v0->GetHiddenInfo());
    }
    break;

  default:
    return -9999.0;
  }

  const Float_t mass1 = info1->GetMass(),
                mass2 = info2->GetMass();

  const AliFemtoThreeVector *momentum1 = info1->GetTrueMomentum(),
                            *momentum2 = info2->GetTrueMomentum();

  const Float_t e1 = momentum1->MassHypothesis(mass1),
                e2 = momentum2->MassHypothesis(mass2);

  return CalcQinv(AliFemtoLorentzVector(e1, *momentum1),
                  AliFemtoLorentzVector(e2, *momentum2));
}


inline double AliFemtoModelCorrFctnQinv::CalcQinv(
  const AliFemtoLorentzVector &p1,
  const AliFemtoLorentzVector &p2)
{
  return std::fabs((p1 - p2).m());


  const double p_inv = (p1 + p2).m2(),
               q_inv = (p1 - p2).m2(),
           mass_diff = p1.m2() - p2.m2();

  const double tQ = (p_inv == 0.0 ? 0.0 : (mass_diff * mass_diff) / p_inv) - q_inv;
  return ::sqrt(std::fabs(tQ));
}


#endif
