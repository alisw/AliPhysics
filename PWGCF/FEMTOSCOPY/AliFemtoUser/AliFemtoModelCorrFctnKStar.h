///
/// \file AliFemtoModelCorrFctnKStar.h
/// \author Andrew Kubera
///

#ifndef ALIFEMTOMODELCORRFCTNKSTAR_H
#define ALIFEMTOMODELCORRFCTNKSTAR_H

class TH1F;
class TH2F;

#include "AliFemtoModelCorrFctn.h"
#include "AliFemtoAvgSepCorrFctn.h"
#include "AliFemtoPair.h"
#include "AliFemtoThreeVector.h"
#include "AliFemtoLorentzVector.h"
#include "AliFemtoModelHiddenInfo.h"

/// \class AliFemtoModelCorrFctnKStar
/// \brief The correlation function which plots numerators and denominator
///        plots from the real monte caro
///
/// \author: Andrew Kubera, andrew.kubera@cern.ch
///
///
class AliFemtoModelCorrFctnKStar : public AliFemtoModelCorrFctn {
public:

  /// Default constructor
  ///
  /// Histograms are created with default titles/ranges.
  /// No other members are initialized - user **MUST** set pair type
  /// and PDG codes.
  ///
  AliFemtoModelCorrFctnKStar();

  /// Construct with histogram parameters
  ///
  /// KStar histograms use these parameters when constructing
  /// histograms. No other memebers are set - user **MUST** set pair
  /// type and PDG codes.
  ///
  AliFemtoModelCorrFctnKStar(const char *name,
                             const Int_t aNbins,
                             const Float_t aQinvLo,
                             const Float_t aQinvHi);

  /// Copy Constructor
  ///
  /// Copies pair type & PDG codes. Clones histograms.
  ///
  AliFemtoModelCorrFctnKStar(const AliFemtoModelCorrFctnKStar &);

  /// Destructor
  virtual ~AliFemtoModelCorrFctnKStar();

  /// Return information about the run of the correlation function
  virtual AliFemtoString Report();

  virtual void AddRealPair(AliFemtoPair*);
  virtual void AddMixedPair(AliFemtoPair*);

  virtual TList* GetOutputList();

  virtual TList* AppendOutputList(TList*);

  virtual void AddOutputObjectsTo(TCollection &);

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
  static double CalcKStar(const AliFemtoLorentzVector&,
                          const AliFemtoLorentzVector&);


  /// Calculate kstar of the pair - implementation chosen by the fPairType
  /// member (Track+Track, Track+V0, etc)
  virtual double CalcTrueKStar(const AliFemtoPair *);

  /// Calculate k* for an AliFemtoPair containing a track and a V0.
  /// The order does not matter
  // virtual double CalcTrueKStarTrackV0(const AliFemtoPair *);

  /// Returns the Track and V0 objects in the pair
  virtual void GetTrackV0(const AliFemtoPair*,
                          AliFemtoTrack*&,
                          AliFemtoV0*&);

private:
  // AliFemtoModelCorrFctnKStar(const AliFemtoModelCorrFctnKStar&);
  AliFemtoModelCorrFctnKStar& operator=(const AliFemtoModelCorrFctnKStar&);
protected:

  /// The particle class combination types (track, V0, etc...)
  AliFemtoAvgSepCorrFctn::PairType fPairType;

  Int_t fExpectedTrack1Code;
  Int_t fExpectedTrack2Code;

  /// KStar numerator of reconstructed pairs
  TH1F *fResNum;

  /// KStar denominator of reconstructed pairs
  TH1F *fResDen;

  /// Numerator of only true pairs
  TH2F *fTrueNum;

  /// Denominator of only true pairs
  TH2F *fTrueDen;

private:

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoModelCorrFctnKStar, 0);
  /// \endcond
#endif

};

inline void
AliFemtoModelCorrFctnKStar::SetPairType(const AliFemtoAvgSepCorrFctn::PairType pairtype)
{
  fPairType = pairtype;
}

inline
void AliFemtoModelCorrFctnKStar::SetExpectedPDGCodes(const Int_t track1_pdg_code,
                                                     const Int_t track2_pdg_code)
{
  fExpectedTrack1Code = track1_pdg_code;
  fExpectedTrack2Code = track2_pdg_code;
}

inline void AliFemtoModelCorrFctnKStar::GetTrackV0(const AliFemtoPair* pair,
                                                   AliFemtoTrack *&track_return,
                                                   AliFemtoV0 *&v0_return)
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

inline double AliFemtoModelCorrFctnKStar::CalcTrueKStar(const AliFemtoPair *pair)
{
  const AliFemtoModelHiddenInfo *info1 = NULL,
                                *info2 = NULL;
  switch (fPairType) {
  case AliFemtoAvgSepCorrFctn::kTracks:
    info1 = dynamic_cast<const AliFemtoModelHiddenInfo*>(pair->Track1()->Track()->GetHiddenInfo());
    info2 = dynamic_cast<const AliFemtoModelHiddenInfo*>(pair->Track2()->Track()->GetHiddenInfo());
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

  const Float_t e1 = sqrt(mass1 * mass1 + momentum1->Mag2()),
                e2 = sqrt(mass2 * mass2 + momentum2->Mag2());

  return CalcKStar(AliFemtoLorentzVector(e1, *momentum1),
                   AliFemtoLorentzVector(e2, *momentum2));
}


inline double AliFemtoModelCorrFctnKStar::CalcKStar(
  const AliFemtoLorentzVector &p1,
  const AliFemtoLorentzVector &p2)
{
  const double p_inv = (p1 + p2).m2(),
               q_inv = (p1 - p2).m2(),
           mass_diff = p1.m2() - p2.m2();

  const double tQ = ::pow(mass_diff, 2) / p_inv - q_inv;
  return ::sqrt(tQ) / 2.0;
}


#endif
