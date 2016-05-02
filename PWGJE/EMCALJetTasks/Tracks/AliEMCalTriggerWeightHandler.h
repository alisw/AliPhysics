/**
 * \file AliEMCalTriggerWeightHandler.h
 * \brief Weight handler for the analysis of high-\f$ p_{t} \f$ tracks in EMCAL-triggered events
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Mar 15, 2015
 */
#ifndef ALIEMCALTRIGGERWEIGHTHANDLER_H
#define ALIEMCALTRIGGERWEIGHTHANDLER_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

class TF1;
class TObjArray;

class AliGenPythiaEventHeader;
class AliMCEvent;

/**
 * \namespace EMCalTriggerPtAnalysis
 * \brief Analysis of high-\f$ p_{t} \f$ tracks in triggered events
 *
 * This namespace contains classes for the analysis of high-\f$ p_{t} \f$ tracks in
 * triggered events.
 */
namespace EMCalTriggerPtAnalysis {

class AliEMCalTriggerPtHardWeight : public TObject {
public:
  AliEMCalTriggerPtHardWeight():  TObject(), fPtMin(-1), fPtMax(-1), fWeight(-1) {}
  AliEMCalTriggerPtHardWeight(double ptmin, double ptmax, double weight) :  TObject(), fPtMin(ptmin), fPtMax(ptmax), fWeight(weight) {}
  virtual ~AliEMCalTriggerPtHardWeight() {}

  void SetPtRange(Double_t ptmin, Double_t ptmax) { fPtMin = ptmin; fPtMax = ptmax; }
  void SetWeight(double weight) { fWeight = weight; }

  bool IsSelected(Double_t pthard) const { return pthard >= fPtMin && pthard < fPtMax; }
  Double_t GetWeight() const { return fWeight; }
  Double_t GetPtMin() const { return fPtMin; }
  Double_t GetPtMax() const { return fPtMax; }
  void GetPtLimits(Double_t &ptmin, Double_t &ptmax) const { ptmin = fPtMin; ptmax = fPtMax; }

private:
  double                  fPtMin;             ///< Min of the pthard bin
  double                  fPtMax;             ///< Max of the pthard bin
  double                  fWeight;            ///< Weight being applied to the pt-hard bin

  /// \cond CLASSIMP
  ClassDef(AliEMCalTriggerPtHardWeight, 1)
  /// \endcond
};

/**
 * \class AliEMCalTriggerWeightHandler
 * \brief Weight handler
 *
 * Weight handler, assigning an event-dependent weight. The weight is coming from a weight model,
 * which is based on an analytic description. For the moment it is assumed that the event depends
 * on the \f$ p_{t} \f$ of the hard interaction.
 */
class AliEMCalTriggerWeightHandler : public TObject {
public:
  AliEMCalTriggerWeightHandler();
  virtual ~AliEMCalTriggerWeightHandler();

  /**
   * Defines whether we use the cross section as weight.
   * \param useCrossSection Define whether to use the cross section as event weight
   */
  void SetUseCrossSection(bool useCrossSection) { fUseCrossSection = useCrossSection; }

  /**
   * Set the weight model
   * \param model The weight model
   */
  void SetWeightModel(const TF1 *model) { fWeightModel = model; }

  void SetWeightForBin(double ptmin, double ptmax, double weight);

  double GetEventWeight(const AliMCEvent *const event) const;
  double GetEventWeight(const AliGenPythiaEventHeader * const header) const;

protected:
  const AliEMCalTriggerPtHardWeight *FindWeight(Double_t pthard) const;

private:
  const TF1         *fWeightModel;          ///< Weight model
  TObjArray         *fBinWeights;           ///< Container for weights in a given pt-hard bin
  bool               fUseCrossSection;      ///< Calculate weight using pt-hard

  /// \cond CLASSIMP
  ClassDef(AliEMCalTriggerWeightHandler, 1);
  /// \endcond
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERWEIGHTHANDLER_H */
