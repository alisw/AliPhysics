#ifndef ALIEMCALTRIGGERWEIGHTHANDLER_H
#define ALIEMCALTRIGGERWEIGHTHANDLER_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

class TF1;
class TObjArray;

class AliGenPythiaEventHeader;
class AliMCEvent;

namespace PWGJE{ 
  
namespace EMCALJetTasks {

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

  ClassDef(AliEMCalTriggerPtHardWeight, 1)
};

/**
 * @class AliEMCalTriggerWeightHandler
 * @brief Weight handler
 * @ingroup PWGJETASKS
 * @author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * @date Mar 15, 2015
 *
 * Weight handler, assigning an event-dependent weight. The weight is coming from a weight model,
 * which is based on an analytic description. For the moment it is assumed that the event depends
 * on the \f$ p_{t} \f$ of the hard interaction.
 */
class AliEMCalTriggerWeightHandler : public TObject {
public:
  /**
   * Constructor
   */
  AliEMCalTriggerWeightHandler();

  /**
   * Copy constructor
   * @param ref Reference for the copy
   */
  AliEMCalTriggerWeightHandler(const AliEMCalTriggerWeightHandler &ref);

  /**
   * Assignment operator
   * @param ref Reference for the assignment
   * @return Object after assignment
   */
  AliEMCalTriggerWeightHandler &operator=(const AliEMCalTriggerWeightHandler &ref);

  /**
   * Destructor, cleanup memory assigned
   */
  virtual ~AliEMCalTriggerWeightHandler();

  /**
   * Defines whether we use the cross section as weight.
   * @param[in] useCrossSection Define whether to use the cross section as event weight
   */
  void SetUseCrossSection(bool useCrossSection) { fUseCrossSection = useCrossSection; }

  /**
   * Set the weight model
   * @param[in] model The weight model
   */
  void SetWeightModel(const TF1 *model) { fWeightModel = model; }

  /**
   * Set weight for a given pt-hard bin to the list of weights. Creates the
   * container if not yet existing.
   * @param[in] ptmin Min. \f$ p_{t} \f$ of the \f$ p_{t} \f$-hard bin
   * @param[in] ptmax Max. \f$ p_{t} \f$ of the \f$ p_{t} \f$-hard bin
   * @param[in] weight Bin weight
   */
  void SetWeightForBin(double ptmin, double ptmax, double weight);

  /**
   * Get weight for event
   * @param[in] event Input event
   * @return the weight calculated for the event
   */
  double GetEventWeight(const AliMCEvent *const event) const;
  /**
   * Get weight for event using a given pythia event header
   * @param[in] header Pythia Event Header
   * @return the weight calculated for the event
   */
  double GetEventWeight(const AliGenPythiaEventHeader * const header) const;

protected:
  /**
   * Find weihgt for pt-hard value in the list of weights
   * @param[in] pthard Pt-hard value to find a bin for
   * @return weight for the pthard bin (if found), NULL otherwise
   */
  const AliEMCalTriggerPtHardWeight *FindWeight(Double_t pthard) const;

private:
  const TF1         *fWeightModel;          ///< Weight model
  TObjArray         *fBinWeights;           ///< Container for weights in a given pt-hard bin
  bool               fUseCrossSection;      ///< Calculate weight using pt-hard

  ClassDef(AliEMCalTriggerWeightHandler, 1);
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIEMCALTRIGGERWEIGHTHANDLER_H */
