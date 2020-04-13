#ifndef ALIEMCALCORRECTIONCELLENERGYVARIATION_H
#define ALIEMCALCORRECTIONCELLENERGYVARIATION_H

#include "AliEmcalCorrectionComponent.h"

#include "TF1.h"

/**
 * @class AliEmcalCorrectionCellEnergyVariation
 * @ingroup EMCALCORRECTIONFW
 * @brief Cell energy variation component in the EMCal correction framework.
 *
 * This component allows the EMCal cell energy to be scaled by a constant factor or as a function of cell E.
 * The shifted energy is written in place.
 *
 * This component is necessary when embedding MC into data, in order that the energy scale in MC
 * matches that in data (without relying on the MC non-linearity correction to correct the energy scale).
 *
 * @author James Mulligan <james.mulligan@berkeley.edu>, UC Berkeley
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
  
  // Load cell energy scale function TF1 into fEnergyScaleFunction
  void LoadEnergyScaleFunction(const std::string & path, const std::string & name);
  
  Double_t               fMinCellE;                       ///< Min cell E to perform scaling on
  Double_t               fMaxCellE;                       ///< Max cell E to perform scaling on
  Double_t               fEnergyScaleFactorConstant;      ///< Constant factor to scale cell energy by
  TF1*                   fEnergyScaleFunction;            ///< Function describing factor to scale cell energy by, as a function of cell E

 private:
  AliEmcalCorrectionCellEnergyVariation(const AliEmcalCorrectionCellEnergyVariation &);               // Not implemented
  AliEmcalCorrectionCellEnergyVariation &operator=(const AliEmcalCorrectionCellEnergyVariation &);    // Not implemented

  // Allows the registration of the class so that it is available to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionCellEnergyVariation> reg;
  
  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionCellEnergyVariation, 2); // EMCal cell energy variation component
  /// \endcond
};

#endif /* ALIEMCALCORRECTIONCELLENERGYVARIATION_H */
