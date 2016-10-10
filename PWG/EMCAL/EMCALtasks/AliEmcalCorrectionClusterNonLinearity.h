#ifndef ALIEMCALCORRECTIONCLUSTERNONLINEARITY_H
#define ALIEMCALCORRECTIONCLUSTERNONLINEARITY_H

#include "AliEmcalCorrectionComponent.h"

#include "AliEMCALRecoUtils.h"

class AliEmcalCorrectionClusterNonLinearity : public AliEmcalCorrectionComponent {
 public:

#if !(defined(__CINT__) || defined(__MAKECINT__))
  std::map <std::string, AliEMCALRecoUtils::NonlinearityFunctions> nonlinearityFunctionMap = {
    { "kPi0MC", AliEMCALRecoUtils::kPi0MC },
    { "kPi0GammaGamma", AliEMCALRecoUtils::kPi0GammaGamma },
    { "kPi0GammaConversion", AliEMCALRecoUtils::kPi0GammaConversion },
    { "kNoCorrection", AliEMCALRecoUtils::kNoCorrection },
    { "kBeamTest", AliEMCALRecoUtils::kBeamTest },
    { "kBeamTestCorrected", AliEMCALRecoUtils::kBeamTestCorrected },
    { "kPi0MCv2", AliEMCALRecoUtils::kPi0MCv2 },
    { "kPi0MCv3", AliEMCALRecoUtils::kPi0MCv3 },
    { "kBeamTestCorrectedv2", AliEMCALRecoUtils::kBeamTestCorrectedv2 },
    { "kSDMv5", AliEMCALRecoUtils::kSDMv5 },
    { "kPi0MCv5", AliEMCALRecoUtils::kPi0MCv5 },
    { "kSDMv6", AliEMCALRecoUtils::kSDMv6 },
    { "kPi0MCv6", AliEMCALRecoUtils::kPi0MCv6 },
    { "kBeamTestCorrectedv3", AliEMCALRecoUtils::kBeamTestCorrectedv3 }
  };
#endif

  AliEmcalCorrectionClusterNonLinearity();
  virtual ~AliEmcalCorrectionClusterNonLinearity();

  // Sets up and runs the task
  Bool_t Initialize();
  Bool_t Run();
  
protected:
  TH1F                  *fEnergyDistBefore;          //!<!energy distribution before
  TH2F                  *fEnergyTimeHistBefore;      //!<!energy/time distribution before
  TH1F                  *fEnergyDistAfter;           //!<!energy distribution after
  TH2F                  *fEnergyTimeHistAfter;       //!<!energy/time distribution after
  
 private:
  AliEmcalCorrectionClusterNonLinearity(const AliEmcalCorrectionClusterNonLinearity &);               // Not implemented
  AliEmcalCorrectionClusterNonLinearity &operator=(const AliEmcalCorrectionClusterNonLinearity &);    // Not implemented

  // Allows the registration of the class so that it is availble to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionClusterNonLinearity> reg;

  ClassDef(AliEmcalCorrectionClusterNonLinearity, 1) // EMCal cluster non-linearity correction component
};

#endif /* ALIEMCALCORRECTIONCLUSTERNONLINEARITY_H */
