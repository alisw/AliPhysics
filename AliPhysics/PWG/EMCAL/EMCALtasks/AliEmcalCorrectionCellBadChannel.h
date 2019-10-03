#ifndef ALIEMCALCORRECTIONCELLBADCHANNEL_H
#define ALIEMCALCORRECTIONCELLBADCHANNEL_H

#include "AliEmcalCorrectionComponent.h"

/**
 * @class AliEmcalCorrectionCellBadChannel
 * @ingroup EMCALCORRECTIONFW
 * @brief Bad channel correction component in the EMCal correction framework.
 *
 * Sets cells marked as bad to E = 0, using OADB bad channel map. The original cell information in the event **will be overwritten**.
 *
 * Based on code in AliEMCALTenderSupply.
 *
 * @author Deepa Thomas (Utrecht University), AliEMCALTenderSupply
 * @author Jiri Kral (University of Jyvaskyla), AliEMCALTenderSupply mods/rewrite
 * @author Salvatore Aiola, make AliEMCALTenderSupply work for AODs
 * @author C. Loizides, make AliEMCALTenderSupply work for AODs
 * @author Gustavo Conesa, LPSC-Grenoble, AliEMCALTenderSupply several mods.
 * @author Raymond Ehlers <raymond.ehlers@yale.edu>, Yale University, centralize EMCal corrections using components
 * @author James Mulligan <james.mulligan@yale.edu>, Yale University, centralize EMCal corrections using components
 * @date Jul 8, 2016
 */

class AliEmcalCorrectionCellBadChannel : public AliEmcalCorrectionComponent {
 public:
  AliEmcalCorrectionCellBadChannel();
  virtual ~AliEmcalCorrectionCellBadChannel();

  // Sets up and runs the task
  Bool_t Initialize();
  void UserCreateOutputObjects();
  Bool_t Run();
  Bool_t CheckIfRunChanged();
  
protected:
  TH1F* fCellEnergyDistBefore;              //!<! cell energy distribution, before bad channel correction
  TH1F* fCellEnergyDistAfter;               //!<! cell energy distribution, after bad channel correction
  
private:

  AliEmcalCorrectionCellBadChannel(const AliEmcalCorrectionCellBadChannel &);             // Not implemented
  AliEmcalCorrectionCellBadChannel &operator=(const AliEmcalCorrectionCellBadChannel &);   // Not implemented

  // Allows the registration of the class so that it is availble to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionCellBadChannel> reg;

  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionCellBadChannel, 1); // EMCal cell bad channel correction component
  /// \endcond
};

#endif /* ALIEMCALCORRECTIONCELLBADCHANNEL_H */
