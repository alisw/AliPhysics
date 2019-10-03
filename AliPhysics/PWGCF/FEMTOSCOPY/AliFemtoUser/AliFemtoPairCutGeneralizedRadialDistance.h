///
/// \file AliFemtoPairCutGeneralizedRadialDistance.h
/// \author M. Szymanski, maszyman@cern.ch Warsaw University of Technology
///

#ifndef ALIFEMTOPAIRCUTGENERALIZEDRADIALDISTANCE_H
#define ALIFEMTOPAIRCUTGENERALIZEDRADIALDISTANCE_H

/**
 * \class AliFemtoPairCutGeneralizedRadialDistance
 * \brief a pair cut which checks selects pairs based on generalized angular separation between tracks

 * The AliFemtoPairCutGeneralizedRadialDistance selects pairs by checking the generalized 
 * angular separation of the tracksin the pair. It is based on the code from the 
 * AliAnalysisTaskProtonLambda2d classwritten by Hans Beck. It is applicable also for 
 * secondary particles. The algorithmpropagates each track taking into account its spatial position.
 * AliFemtoEventReaderAOD::SetShiftedPositions method sets the spatial position of the
 * track at the selected radius in the shifted coordinate system. From this values, the
 * distance between two tracks is calculated in
 * StoreDPhiStarDEtaStar function. One should set the radius at which the
 * generalized angular distance is calculated by calling SetRadius method or in
 * constructor. It should be the same value which is given to AliFemtoEventReaderAOD
 * object in SetDoShiftPosition method.
 *
 * More info: http://cds.cern.ch/record/2047435/files/CERN-THESIS-2015-123.pdf
 *
 */

#include "AliFemtoPairCut.h"
#include "AliFemtoShareQualityPairCut.h"
#include "AliFemtoPairCutAntiGamma.h"

#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"

class AliFemtoPairCutGeneralizedRadialDistance : public AliFemtoPairCutAntiGamma {
public:
  AliFemtoPairCutGeneralizedRadialDistance();
  AliFemtoPairCutGeneralizedRadialDistance(const AliFemtoPairCutGeneralizedRadialDistance& c);
  virtual ~AliFemtoPairCutGeneralizedRadialDistance();
  AliFemtoPairCutGeneralizedRadialDistance& operator=(const AliFemtoPairCutGeneralizedRadialDistance& c);

  virtual bool Pass(const AliFemtoPair* pair);
  virtual AliFemtoString Report();
  virtual TList *ListSettings();
  virtual AliFemtoPairCut* Clone();

  void SetPhiStarDifferenceMinimum(Double_t);
  void SetEtaStarDifferenceMinimum(Double_t);
  void SetRadius(Double_t);

 protected:
  Double_t fDPhiStarMin;             ///< Minimum allowed pair separation in PhiStar
  Double_t fDEtaStarMin;             ///< Minimum allowed pair separation in EtaStar
  Double_t fRadius;                  ///< Radial distance at which EtaStar and PhiStar are calculated

#ifdef __ROOT__
  ClassDef(AliFemtoPairCutGeneralizedRadialDistance, 0)
#endif
};

inline AliFemtoPairCut* AliFemtoPairCutGeneralizedRadialDistance::Clone() { AliFemtoPairCutGeneralizedRadialDistance* c = new AliFemtoPairCutGeneralizedRadialDistance(*this); return c;}

#endif
