#ifndef ALIEMCALCORRECTIONPHOSCORRECTIONS_H
#define ALIEMCALCORRECTIONPHOSCORRECTIONS_H

#include "AliEmcalCorrectionComponent.h"

#include "AliPHOSTenderSupply.h"

/**
 * @class AliEmcalCorrectionPHOSCorrections
 * @ingroup EMCALCOREFW
 * @brief Wrapper for PHOS Tender to be executed through the EMCal correction framework.
 *
 * @author James Mulligan <james.mulligan@yale.edu>, Yale University
 * @author Raymond Ehlers <raymond.ehlers@yale.edu>, Yale University
 * @date June 27, 2017
 */

class AliEmcalCorrectionPHOSCorrections : public AliEmcalCorrectionComponent {
 public:
  AliEmcalCorrectionPHOSCorrections();
  virtual ~AliEmcalCorrectionPHOSCorrections();

  // Sets up and runs the task
  Bool_t Initialize();
  void UserCreateOutputObjects();
  Bool_t Run();
  Bool_t CheckIfRunChanged();
  
  AliPHOSTenderSupply*          GetPHOSTenderSupply() {return fPHOSTender;}

 protected:
  AliPHOSTenderSupply*          fPHOSTender;                            //!<! Pointer to PHOS Tender Supply
  Bool_t                        fIsMC;                                  ///< Flag for MC event
  TString                       fOptions;                               ///< MC decalibration option

 private:
  AliEmcalCorrectionPHOSCorrections(const AliEmcalCorrectionPHOSCorrections &);               // Not implemented
  AliEmcalCorrectionPHOSCorrections &operator=(const AliEmcalCorrectionPHOSCorrections &);    // Not implemented

  // Allows the registration of the class so that it is availble to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionPHOSCorrections> reg;

  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionPHOSCorrections, 1); // PHOS correction component
  /// \endcond
};

#endif /* ALIEMCALCORRECTIONPHOSCORRECTIONS_H */
