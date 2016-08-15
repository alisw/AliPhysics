#ifndef ALIQNCORRECTIONS_DATAVECTORSCHAN_H
#define ALIQNCORRECTIONS_DATAVECTORSCHAN_H

/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file AliQnCorrectionsDataVectorChannelized.h
/// \brief Class that models data vectors from channelized detectors within the Q vector correction framework
///

#include "AliQnCorrectionsDataVector.h"

/// \class AliQnCorrectionsDataVectorChannelized
/// \brief Data vector class from a channelized detector
///
/// The class expands the data vector one to incorporate channel id and
/// two set of weights to support channel equalization procedures.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Jan 27, 2016
class AliQnCorrectionsDataVectorChannelized : public AliQnCorrectionsDataVector{
public:

  AliQnCorrectionsDataVectorChannelized();
  AliQnCorrectionsDataVectorChannelized(Int_t channelId, Float_t phi, Float_t weight);
  virtual ~AliQnCorrectionsDataVectorChannelized();

  /// Sets the equalized weight
  /// \param weight equalized weight after channel equalization
  void SetEqualizedWeight( Float_t weight) { fEqualizedWeight = weight; }

  /// Gets the weight for the data vector
  /// \return the raw weight
  virtual Float_t Weight() { return fWeight; }
  /// Gets the equalized weight for the data vector
  /// \return the equalized weight
  virtual Float_t EqualizedWeight() { return fEqualizedWeight; }

private:
  Float_t fEqualizedWeight;       //!<! equalized weight after channel equalization

/// \cond CLASSIMP
  ClassDef(AliQnCorrectionsDataVectorChannelized, 1);
/// \endcond
};

#endif /* ALIQNCORRECTIONS_DATAVECTORSCHAN_H */
