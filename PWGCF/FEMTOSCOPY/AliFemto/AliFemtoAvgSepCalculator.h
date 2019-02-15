///
/// \file AliFemtoAvgSepCalculator.h
///

#pragma once

#ifndef ALIFEMTOAVGSEPCALCULATOR_H_
#define ALIFEMTOAVGSEPCALCULATOR_H_


#include "AliFemtoV0.h"
#include "AliFemtoThreeVector.h"


/// \class AliFemtoAvgSepCalculator
/// \brief A helper class for calculating the average separation of track and
///        v0 particles.
///
/// This class is designed to be the common location of algorithms which
/// calculate the average separation between two particles. This currently
/// works with AliFemtoTracks and AliFemtoV0s.
///
/// The goal of this class is to provide a single location that various track
/// cuts may use to accurately calculate the average separation of two tracks
/// in the TPC. Previous to this, each cut had its own implementation in its
/// Pass method; when a bug was found, every class had to be checked and fixed.
/// This class minimizes code repetition, and provides a one-stop-shop for
/// anybody interested in understanding the algorithm.
///
/// The current use-case is to generate a 'calculator' object from the two
/// particles, then calling one of the 'passes' methods (i.e.
/// track_v0_passes(...)) with minimum average separation to see if this pair
/// passes the cut. This may change, as having a single object with the minimum
///
///
class AliFemtoAvgSepCalculator {
public:

  /**
   * Static method used to determine if a vector was actually set or just
   * created using default values (-10000.).
   */
  static bool is_valid_point(const AliFemtoThreeVector &v) {
    return v.x() > -9999.0 && v.y() > -9999.0 && v.z() > -9999.0;
  }


  /**
   * Construct out of two tracks
   *
   * This is the simplest case, where each track has one set of TPC points.
   */
  AliFemtoAvgSepCalculator(const AliFemtoTrack*, const AliFemtoTrack*);


  /**
   * Construct from track and V0
   *
   * This keeps a sigle list of track points and two lists (pos+neg) of V0
   * points
   */
  AliFemtoAvgSepCalculator(const AliFemtoTrack*, const AliFemtoV0*);

  /**
   * Construct from two V0s
   *
   * This keeps 4 vectors of track points for pos+neg daugthers of both V0
   * particles
   */
  AliFemtoAvgSepCalculator(const AliFemtoV0*, const AliFemtoV0*);


  /**
   * Calculates the separation between each point in the two vectors and
   * returns the average value.
   *
   * If the two vectors are different sizes, it stops after the shorter is
   * exhausted.
   */
  static float CalcAvgSep(const std::vector<AliFemtoThreeVector>& coll_0,
                          const std::vector<AliFemtoThreeVector>& coll_1);



  bool track_track_passes(const float avg_sep_pos);
  bool track_v0_passes(const float avg_sep_pos, const float avg_sep_neg);
  bool v0_v0_passes(const float avsep_pospos,
                    const float avsep_posneg,
                    const float avsep_negpos,
                    const float avsep_negneg);

  /*
   * Collection of average separations - currently not used
   */
  // union {
  //   float avg_sep;
  //   struct {
  //     float avg_sep_pos;
  //     float avg_sep_neg;
  //   };
  //   struct {
  //     float avg_sep_pospos;
  //     float avg_sep_posneg;
  //     float avg_sep_negpos;
  //     float avg_sep_negneg;
  //   };
  // } min;

  typedef std::vector<std::vector<AliFemtoThreeVector> > PointCollection_t;

protected:

  /// Collection of vectors
  PointCollection_t t1;
  PointCollection_t t2;

  static void push_tpc_point(const AliFemtoThreeVector &v,
                             std::vector<AliFemtoThreeVector> &storage,
                             bool& save)
  {
    if (is_valid_point(v)) {
      storage.push_back(v);
    } else {
      save = false;
    }
  };
};


inline
AliFemtoAvgSepCalculator::AliFemtoAvgSepCalculator(const AliFemtoTrack *track_1,
                                                   const AliFemtoTrack *track_2):
  t1(1),
  t2(1)
{
  bool save_1 = true,
       save_2 = true;

  for (int i = 0; i < 8 && (save_1 || save_2); ++i) {
    if (save_1) {
      push_tpc_point(track_1->NominalTpcPoint(i), t1[0], save_1);
    }
    if (save_2) {
      push_tpc_point(track_2->NominalTpcPoint(i), t2[0], save_2);
    }
  }
}


inline
AliFemtoAvgSepCalculator::AliFemtoAvgSepCalculator(const AliFemtoTrack *track,
                                                   const AliFemtoV0 *v0):
  t1(1),
  t2(2)
{
  bool save_track = true,
       save_v0_pos = true,
       save_v0_neg = true;

  for (int i = 0; i < 8 && (save_track || save_v0_pos || save_v0_neg); ++i) {
    if (save_track) {
      push_tpc_point(track->NominalTpcPoint(i), t1[0], save_track);
    }
    if (save_v0_pos) {
      push_tpc_point(v0->NominalTpcPointPos(i), t2[0], save_v0_pos);
    }
    if (save_v0_neg) {
      push_tpc_point(v0->NominalTpcPointNeg(i), t2[1], save_v0_neg);
    }
  }
}


inline
AliFemtoAvgSepCalculator::AliFemtoAvgSepCalculator(const AliFemtoV0 *v0_1,
                                                   const AliFemtoV0 *v0_2):
  t1(2),
  t2(2)
{
  bool save_v0_1_pos = true,
       save_v0_1_neg = true,
       save_v0_2_pos = true,
       save_v0_2_neg = true;

  for (int i = 0; i < 8; ++i) {
    if (save_v0_1_pos) {
      push_tpc_point(v0_1->NominalTpcPointPos(i), t1[0], save_v0_1_pos);
    }
    if (save_v0_1_neg) {
      push_tpc_point(v0_1->NominalTpcPointNeg(i), t1[1], save_v0_1_neg);
    }
    if (save_v0_2_pos) {
      push_tpc_point(v0_2->NominalTpcPointPos(i), t2[0], save_v0_2_pos);
    }
    if (save_v0_2_neg) {
      push_tpc_point(v0_2->NominalTpcPointNeg(i), t2[1], save_v0_2_neg);
    }
  }
}


inline float
AliFemtoAvgSepCalculator
  ::CalcAvgSep(const std::vector<AliFemtoThreeVector>& coll_0,
               const std::vector<AliFemtoThreeVector>& coll_1)
{
  const UInt_t size = TMath::Min(coll_0.size(), coll_1.size());
  float avg_sep = 0.0;
  for (UInt_t i = 0; i < size; ++i) {
    avg_sep += (coll_0[i] - coll_1[i]).Mag();
  }
  return avg_sep / size;
}


inline bool
AliFemtoAvgSepCalculator::track_track_passes(const float avg_sep)
{
  return CalcAvgSep(t1[0], t2[0]) >= avg_sep;
}


inline bool
AliFemtoAvgSepCalculator::track_v0_passes(const float avg_sep_pos,
                                          const float avg_sep_neg)
{
  return CalcAvgSep(t1[0], t2[0]) >= avg_sep_pos
      && CalcAvgSep(t1[0], t2[1]) >= avg_sep_neg;
}


inline bool
AliFemtoAvgSepCalculator::v0_v0_passes(const float avsep_pospos,
                                       const float avsep_posneg,
                                       const float avsep_negpos,
                                       const float avsep_negneg)
{
  return CalcAvgSep(t1[0], t2[0]) >= avsep_pospos
      && CalcAvgSep(t1[0], t2[1]) >= avsep_posneg
      && CalcAvgSep(t1[1], t2[0]) >= avsep_negpos
      && CalcAvgSep(t1[1], t2[1]) >= avsep_negneg;
}

#endif
