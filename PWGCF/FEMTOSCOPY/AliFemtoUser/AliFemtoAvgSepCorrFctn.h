///
/// \file AliFemtoAvgSepCorrFctn.h
/// \author M. Janik, L. Graczykowski, Warsaw University of Technology
///

#ifndef ALIFEMTOAVGSEPCORRFCTN_H
#define ALIFEMTOAVGSEPCORRFCTN_H

#include "TH1D.h"
#include "AliFemtoCorrFctn.h"

/**
 * \class AliFemtoAvgSepCorrFctn
 * \brief A correlation function built from comparing the separation between
 *        tracks from nominal points in the TPC.
 *
 * The AliFemtoAvgSepCorrFctn stores numerator and denomiator histograms
 * counting the pair frequency of the average separation of the tracks in the
 * pair, measured in centimeters. This average separation is calculated by
 * progressing through the "NominalTpcPoints" as reported by the AliFemtoTrack
 * method NominalTpcPoint(int). The distance between the two tracks at each
 * point is averaged and this quantity is stored.
 *
 * This class enables finding the average separation of multiple "types" of
 * particle pairs. By using the SetPairType method, you specify which algorithm
 * to use when calculating average separation and which histograms will be
 * written upon completion. The options are specified by the enum 'PairType'.
 * The default value is kTracks, which simply accepts pairs of AliFemtoTracks
 * and has a simple numerator and denominator histogram. Setting the PairType
 * to kTrackV0 accepts pairs where the first track is a AliFemtoV0 and the
 * second is an AliFemtoTrack; there are two sets of histograms produced with
 * this configuration, between track and the positive & negative daughters of
 * the V0. The case kV0 accepts pairs of V0 tracks, and produces 4 sets of
 * histograms, for each combination of positive and negative daughters for both
 * of the V0s.
 *
 * Note: If the track does not pass through the entire TPC, some of the points
 * are left at their default value (-9999). The averaging algorithm prior to
 * September 2015 did not check for points which were unset; since then it will
 * average only the points which have been set. This has the potential of
 * calculating the average separation from only a few points, which may be
 * undesireable. The number of misses is very dependent on pT, and thus it is
 * recommended to take the average separation into account when determining
 * your pT cut.
 */
class AliFemtoAvgSepCorrFctn : public AliFemtoCorrFctn {
public:
  /**
   * Enumeration of the possible pair configuration of AliFemtoParticle types
   * to expect in the AliFemtoPair pairs given to the AddRealPair and
   * AddMixedPair methods.
   *  - kTracks: Both AliFemtoTracks
   *  - kTrackV0: Track1=AliFemtoV0, Track2=AliFemtoTrack
   *  - kV0s: Both AliFemtoV0
   */
  enum PairType {
    kTracks = 0,
    kTrackV0 = 1,
    kV0s = 2,
    kTrackXi = 3,
    kXiV0 = 4
  };
  typedef enum PairType AliFemtoPairType;

  /**
   * Construct with histogram parameters
   */
  AliFemtoAvgSepCorrFctn(const char *title, const int &nbins, const float &Low, const float &High);

  /// Copy Constructor
  AliFemtoAvgSepCorrFctn(const AliFemtoAvgSepCorrFctn &aCorrFctn);

  /// Destructor
  virtual ~AliFemtoAvgSepCorrFctn();

  /// Assignment Operator
  AliFemtoAvgSepCorrFctn &operator=(const AliFemtoAvgSepCorrFctn &aCorrFctn);

  virtual AliFemtoString Report();

  /**
   * Extract the average separation between the two particles in the
   * AliFemtoPair and store in appropriate histograms.
   */
  virtual void AddRealPair(AliFemtoPair *aPair);
  virtual void AddMixedPair(AliFemtoPair *aPair);
  virtual void Finish();

  TH1D *Numerator();
  TH1D *Denominator();
  TH1D *Ratio();

  virtual TList *GetOutputList();
  void Write();
  void SetPairType(AliFemtoPairType pairtype);

private:
  //2 tracks
  TH1D *fNumerator;          ///< numerator - real pairs
  TH1D *fDenominator;        ///< denominator - mixed pairs

  //track + V0
  TH1D *fNumeratorPos;       ///< numerator - real pairs
  TH1D *fDenominatorPos;     ///< denominator - mixed pairs
  TH1D *fNumeratorNeg;       ///< numerator - real pairs
  TH1D *fDenominatorNeg;     ///< denominator - mixed pairs

  //2 V0s
  TH1D *fNumeratorPosPos;    ///< numerator - real pairs
  TH1D *fDenominatorPosPos;  ///< denominator - mixed pairs
  TH1D *fNumeratorPosNeg;    ///< numerator - real pairs
  TH1D *fDenominatorPosNeg;  ///< denominator - mixed pairs
  TH1D *fNumeratorNegPos;    ///< numerator - real pairs
  TH1D *fDenominatorNegPos;  ///< denominator - mixed pairs
  TH1D *fNumeratorNegNeg;    ///< numerator - real pairs
  TH1D *fDenominatorNegNeg;  ///< denominator - mixed pairs

  //track + Xi
  //In addition to fNumeratorPos, fDenominatorPos, fNumeratorNeg, and
  //fDenominatorNeg, the kTrackXi case also needs histograms for average
  //separation between bachelor pion and track
  TH1D *fNumeratorBacTrack;       ///< numerator - real pairs
  TH1D *fDenominatorBacTrack;     ///< denominator - mixed pairs

  //Xi + V0
  //In addition to fNumeratorPosPos, fDenominatorPosPos,
  //fNumeratorPosNeg, fDenominatorPosNeg, fNumeratorNegPos, fDenominatorNegPos,
  //fNumeratorNegNeg, fDenominatorNegNeg, the kXiV0 case also needs histograms 
  //for average separation between bachelor pion and Xi's V0 daughters
  TH1D *fNumeratorBacPos;       ///< numerator - real pairs
  TH1D *fDenominatorBacPos;     ///< denominator - mixed pairs
  TH1D *fNumeratorBacNeg;       ///< numerator - real pairs
  TH1D *fDenominatorBacNeg;     ///< denominator - mixed pairs

  TH1D *fRatio;              ///< ratio - correlation function
  AliFemtoPairType fPairType;


#ifdef __ROOT__
  ClassDef(AliFemtoAvgSepCorrFctn, 1)
#endif
};

inline  TH1D *AliFemtoAvgSepCorrFctn::Numerator()
{
  return fNumerator;
}
inline  TH1D *AliFemtoAvgSepCorrFctn::Denominator()
{
  return fDenominator;
}
inline  TH1D *AliFemtoAvgSepCorrFctn::Ratio()
{
  return fRatio;
}


#endif
