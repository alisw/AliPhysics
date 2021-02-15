/**
 * \file AliReducedMatchedTrack.h
 * \brief Definition of class AliReducedTrack, a structure with reduced track information
 * at reconstruction level
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Mar 17, 2015
 */
#ifndef ALIREDUCEDMATCHEDTRACK_H
#define ALIREDUCEDMATCHEDTRACK_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

class TVector3;

/**
 * \namespace HighPtTracks
 * \brief Namespace for classes creating trees of events with jets
 *
 * This namespace contains classes descibing reduced events with high
 * jets. A jet event consists of the following classes.
 *    - AliReducedJetEvent
 *    - AliReducedJetParticle
 *    - AliReducedJetConstituent
 *    - AliReducedMatchedTrack
 *  The task AliHighPtReconstructionEfficiency produces the reduced jet events
 *  and is included in the namespace as well. Also the helper class AliParticleMap
 *  is part of this namespace.
 */
namespace HighPtTracks {

/**
 * \class AliReducedMatchedTrack
 * \brief Class with reduced track information at reconstruction level
 *
 * Class with reduced track information for tracks matched to a true particle in a jet cone.
 * This class is part of the reduced jet event. Information stored are the three momentum at
 * reconstruction level, the number of clusters in the TPC, the sign bit of the track label,
 * and a list of track cuts, containing \f$ R_{AA} \f$ standard cuts and hybrid track cuts.
 */
class AliReducedMatchedTrack : public TObject {
public:
  /**
   * \enum ETrackCutsType_t
   * \brief Type of track cuts applied in the track selection
   */
  enum ETrackCutsType_t{
    kHighPtStandard,      ///< Standard track cuts
    kHighPtHybrid         ///< Hybrid track cuts
  };
  AliReducedMatchedTrack();
  AliReducedMatchedTrack(double px, double py, double pz);
  virtual ~AliReducedMatchedTrack() {}

  double Pt() const;
  double Eta() const;
  double Phi() const;
  void FillVector(TVector3 &vec);

  /**
   * Set the 3-momentum vector of the reconstructed track
   *
   * \param px x-component of the momentum vector
   * \param py y-component of the momentum vector
   * \param pz x-component of the momentum vector
   */
  void SetPxPyPz(double px, double py, double pz){
    fPx = px;
    fPy = py;
    fPz = pz;
  }

  /**
   * Set the number of clusters in the TPC.
   *
   * \param ncls Number of clusters in the TPC
   */
  void SetNumberOfClustersTPC(unsigned char ncls) { fNclustersTPC = ncls; }

  /**
   * Define track as good track via the MC Label
   *
   * \param goodtrack True if track has a positive MC Label
   */
  void SetGoodTrackLabel(bool goodtrack) { fGoodTrackLabel = goodtrack; }

  /**
   * Mark track as being selected by given track cuts
   *
   * \param cuts Track cuts the track was selected by
   */
  void SetSurvivedTrackCuts(ETrackCutsType_t cuts) { SETBIT(fTrackCuts, cuts); }

  /**
   * Access to number of clusters in the TPC
   *
   * \return The number of clusters in the TPC
   */
	unsigned char GetNumberOfClustersTPC() const { return fNclustersTPC; }

	/**
	 * Check if the track is a good reconstructed tracks, defined via the sign of the Monte-Carlo label.
	 * Good tracks have a positive Monte-Carlo label.
	 *
	 * \return True if track is a good track, false otherwise
	 */
	bool HasGoodTrackLabel() { return fGoodTrackLabel; }

	/**
	 * Check if track was selected under given track cuts. Currently implemented are the standard
	 * track cuts and the hybrid track cuts
	 *
	 * \param cuts Type of cuts to check
	 * \return True if track was selected, false otherwise
	 */
	bool HasSurvivedTrackCuts(ETrackCutsType_t cuts) { return TESTBIT(fTrackCuts, cuts); }

private:
  double fPx;                     ///< x-component of the 3-momentum
  double fPy;                     ///< y-component of the 3-momentum
  double fPz;                     ///< z-component of the 3-momentum
  bool fGoodTrackLabel;           ///< Flag for tracks with positive track label
  unsigned char fNclustersTPC;    ///< Number of clusters in the TPC
  unsigned int fTrackCuts;        ///< bitmap for track selection bits

  /// \cond CLASSIMP
  ClassDef(AliReducedMatchedTrack, 1);
  /// \endcond
};

} /* namespace HighPtTracks */

#endif /* ALIREDUCEDMATCHEDTRACK_H */
