#ifndef ALIEMCALCORRECTIONCELLENERGYVARIATION_H
#define ALIEMCALCORRECTIONCELLENERGYVARIATION_H

#include "AliEmcalCorrectionComponent.h"

/**
 * @class AliEmcalCorrectionCellEnergyVariation
 * @ingroup EMCALCOREFW
 * @brief Cell energy variation component in the EMCal correction framework.
 *
 * This component allows the EMCal cell energy to be scaled by a constant factor.
 * The shifted energy is written in place.
 *
 * This component is necessary when embedding MC into data, in order that the energy scale in MC
 * matches that in data (without relying on the non-linearity correction to correct the energy scale).
 *
 * @author James Mulligan <james.mulligan@yale.edu>, Yale University
 * @date Dec 19 2018
 */

class AliEmcalCorrectionCellEnergyVariation : public AliEmcalCorrectionComponent {
 public:
  AliEmcalCorrectionCellEnergyVariation();
  virtual ~AliEmcalCorrectionCellEnergyVariation();

  // Sets up and runs the task
  Bool_t Initialize();
  void UserCreateOutputObjects();
  void ExecOnce();
  Bool_t Run();
  
protected:
  Double_t               fEnergyScaleShift;               ///< Fraction of cell energy to shift (positive means upward shift)

 private:
  AliEmcalCorrectionCellEnergyVariation(const AliEmcalCorrectionCellEnergyVariation &);               // Not implemented
  AliEmcalCorrectionCellEnergyVariation &operator=(const AliEmcalCorrectionCellEnergyVariation &);    // Not implemented

  // Allows the registration of the class so that it is available to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionCellEnergyVariation> reg;
  
  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionCellEnergyVariation, 1); // EMCal cell energy variation component
  /// \endcond
};

#endif /* ALIEMCALCORRECTIONCELLENERGYVARIATION_H */
