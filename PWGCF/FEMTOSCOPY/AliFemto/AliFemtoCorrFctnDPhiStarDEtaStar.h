///
/// \file AliFemtoCorrFctnDPhiStarDEtaStar.h
/// \author M. Szymanski, maszyman@cern.ch Warsaw University of Technology
///

#ifndef ALIFEMTOCORRFCTNDPHISTARDETASTAR_H
#define ALIFEMTOCORRFCTNDPHISTARDETASTAR_H


/**
 * \class AliFemtoCorrFctnDPhiStarDEtaStar
 * \brief A 2D correlation function showing the generalized angular distance between tracks

 * The AliFemtoCorrFctnDPhiStarDEtaStar stores numerator and denominator histograms
 * counting the pair frequency of the generalized angular separation of the tracks
 * in the pair. It is based on the code from the AliAnalysisTaskProtonLambda2d class
 * written by Hans Beck. It is applicable also for secondary particles. The algorithm
 * propagates each track taking into account its spatial position.
 * AliFemtoEventReaderAOD::SetShiftedPositions method sets the spatial position of the
 * track at the selected radius in the shifted coordinate system. From this values, the
 * distance between two tracks (primary or secondary) is calculated in
 * StoreDPhiStarDEtaStarBetween* functions. One should set the radius at which the
 * generalized angular distance is calculated by calling SetRadius method or in
 * constructor. It should be the same value which is given to AliFemtoEventReaderAOD
 * object in SetDoShiftPosition method.
 *
 * More info: http://cds.cern.ch/record/2047435/files/CERN-THESIS-2015-123.pdf
 *
 * This class enables finding the generalized angular separation of multiple "types" of
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

 */

#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "AliFemtoCorrFctn.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"

class AliFemtoCorrFctnDPhiStarDEtaStar : public AliFemtoCorrFctn {
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
    kV0s = 2
  };
  typedef enum PairType AliFemtoPairType;

  /**
   * Construct with histogram parameters
   */
  AliFemtoCorrFctnDPhiStarDEtaStar(const char* title,
                                   double radius,
                                   const int& aEtaBins,
                                   double aEtaRangeLow,
                                   double aEtaRangeUp,
                                   const int& aPhiStarBins,
                                   double aPhiStarRangeLow,
                                   double aPhiStarRangeUp);

  /// Copy Constructor
  AliFemtoCorrFctnDPhiStarDEtaStar(const AliFemtoCorrFctnDPhiStarDEtaStar& aCorrFctn);

  /// Destructor
  virtual ~AliFemtoCorrFctnDPhiStarDEtaStar();

  /// Assignment Operator
  AliFemtoCorrFctnDPhiStarDEtaStar& operator=(const AliFemtoCorrFctnDPhiStarDEtaStar& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();
  void WriteHistos();
  virtual TList* GetOutputList();
  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoCorrFctnDPhiStarDEtaStar(*this); }

  void SetRadius(double minrad);
  void SetPairType(AliFemtoPairType pairtype);

 private:

  //2 tracks
  TH2D *fDPhiStarDEtaStarNumerator;      ///< Numerator of dPhiStar dEtaStar function
  TH2D *fDPhiStarDEtaStarDenominator;    ///< Denominator of dPhiStar dEtaStar function

  //track + V0
  TH2D *fDPhiStarDEtaStarNumeratorTrackPos;      ///< Numerator of dPhiStar dEtaStar function
  TH2D *fDPhiStarDEtaStarDenominatorTrackPos;    ///< Denominator of dPhiStar dEtaStar function
  TH2D *fDPhiStarDEtaStarNumeratorTrackNeg;      ///< Numerator of dPhiStar dEtaStar function
  TH2D *fDPhiStarDEtaStarDenominatorTrackNeg;    ///< Denominator of dPhiStar dEtaStar function

  //2 V0s
  TH2D *fDPhiStarDEtaStarNumeratorPosPos;      ///< Numerator of dPhiStar dEtaStar function
  TH2D *fDPhiStarDEtaStarDenominatorPosPos;    ///< Denominator of dPhiStar dEtaStar function
  TH2D *fDPhiStarDEtaStarNumeratorPosNeg;      ///< Numerator of dPhiStar dEtaStar function
  TH2D *fDPhiStarDEtaStarDenominatorPosNeg;    ///< Denominator of dPhiStar dEtaStar function
  TH2D *fDPhiStarDEtaStarNumeratorNegPos;      ///< Numerator of dPhiStar dEtaStar function
  TH2D *fDPhiStarDEtaStarDenominatorNegPos;    ///< Denominator of dPhiStar dEtaStar function
  TH2D *fDPhiStarDEtaStarNumeratorNegNeg;      ///< Numerator of dPhiStar dEtaStar function
  TH2D *fDPhiStarDEtaStarDenominatorNegNeg;    ///< Denominator of dPhiStar dEtaStar function

  double fEtaStarRangeLow;               ///< Lower range of EtaStar
  double fEtaStarRangeUp;                ///< Upper range of EtaStar

  double fPhiStarRangeLow;           ///< Lower range of PhiStar
  double fPhiStarRangeUp;            ///< Upper range of PhiStar

  Double_t fMinRad;                  ///< Radial distance at which EtaStar and PhiStar are calculated
  AliFemtoPairType fPairType;

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoCorrFctnDPhiStarDEtaStar, 1);
  /// \endcond
#endif
    };


#endif
