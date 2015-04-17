/**
 * \file AliReducedMCHeader.h
 * \brief Declaration of a reduced MC event header
 *
 * In this file the class AliReducedMCHeader, a reduced MC event header for
 * the reduced event structure, is declared.
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Apr 14, 2015
 */
#ifndef ALIREDUCEDMCHEADER_H
#define ALIREDUCEDMCHEADER_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

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
 * \class AliReducedMCHeader
 * \brief A reduced event header with MC information for the reduced event structure
 *
 * This class contains a reduced event header structure storing basic Monte-Carlo information. The
 * event header contains
 * - The event cross section
 * - The number of trials
 * - The \f$ p_{t} \f$ of the hard interaction
 * This is an optional member of the reduced event structure
 */
class AliReducedMCHeader : public TObject {
public:
  AliReducedMCHeader();
  AliReducedMCHeader(Double_t crosssection, Double_t numberOfTrials, Double_t pthard);
  virtual ~AliReducedMCHeader();

  /**
   * Get the event cross section
   * \return Cross section of the event
   */
  Double_t              GetCrossSection() const { return fCrossSection; }
  /**
   * Get the number of trials of the event
   * \return Number of trials
   */
  Int_t                 GetNumberOfTrials() const { return fNumberOfTrials; }
  /**
   * Get the \f$ p_{t} \f$ hard of the event
   * \return the \f$ p_{t} \f$ hard of the event
   */
  Double_t              GetPtHard() const { return fPtHard; }

  /**
   * Set the event cross section
   * \param crosssection The cross section of the event
   */
  void                  SetCrossSection(Double_t crosssection) { fCrossSection = crosssection; }
  /**
   * Set the number of trials of the event
   * @param ntrials Number of trials
   */
  void                  SetNumberOfTrials(Int_t ntrials) { fNumberOfTrials = ntrials; }
  /**
   * Set the \f$ p_{t} \f$ hard of the event
   * @param pthard the \f$ p_{t} \f$ hard of the event
   */
  void                  SetPtHard(Double_t pthard) { fPtHard = pthard; }

protected:
  Double_t                       fCrossSection;                ///< cross section
  Int_t                          fNumberOfTrials;              ///< number of trials
  Double_t                       fPtHard;                      ///< generated \f$ p_{t} \f$ hard

  /// \cond CLASSIMP
  ClassDef(AliReducedMCHeader,1)
  /// \endcond
};

} /* namespace HighPtTracks */

#endif /* ALIREDUCEDMCHEADER_H */
