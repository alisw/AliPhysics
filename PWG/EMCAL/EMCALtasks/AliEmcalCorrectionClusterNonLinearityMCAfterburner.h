#ifndef ALIEMCALCORRECTIONCLUSTERNONLINEARITYMCAFTERBURNER_H
#define ALIEMCALCORRECTIONCLUSTERNONLINEARITYMCAFTERBURNER_H

#include "AliEmcalCorrectionComponent.h"

/**
 * @class AliEmcalCorrectionClusterNonLinearityMCAfterburner
 * @ingroup EMCALCORRECTIONFW
 * @brief Cluster energy non-linearity correction MC afterburner in the EMCal correction framework.
 *
 * Non-linearity correction to the MC cluster energy is necessary because the current non linearity is done for
 * 50MeV cell threshold and all analyses use a 100MeV cell threshold for the clusterization.
 * Thus, the MC need to be corrected in a second step
 
 The energy of the cluster **after** the MC non-linearity correction can be retrieved using the method `cluster->GetNonLinCorrEnergy()`.
 *
 * Based on code in AliCaloPhotonCuts.cxx.
 *
 * @author 15 active people originally in AliCaloPhotonCuts.cxx
 * @author Eliane Epple, Yale University, include the PCM aproach into the EMCal correction FW using correction components
 * @date Nov 13, 2018
 */


class AliEmcalCorrectionClusterNonLinearityMCAfterburner : public AliEmcalCorrectionComponent {
 public:

  AliEmcalCorrectionClusterNonLinearityMCAfterburner();
  virtual ~AliEmcalCorrectionClusterNonLinearityMCAfterburner();

  // Sets up and runs the task
  Bool_t Initialize();
  void UserCreateOutputObjects();
  Bool_t Run();

  /// MC Nonlinearity afterburner enum list of four possible ways of determining the parameters.
  /// Standard is kPCM_EMCal
  enum  EMCAfterburnerMethod_t
  {
    kPCM_EMCal                 = 0,                    //!<!This is determined by fitting the ratio of MC and Data pi0 mass position
    kEMCal_EMCal               = 1,                    //!<!This is determined by fitting the ratio of MC and Data pi0 mass position
    kPCM_EMCalFunctionRatio    = 2,                    //!<!This is determined by fitting MC and Data pi0 mass position respectivley and dividing the fits
    kEMCal_EMCalFunctionRatio  = 3,                    //!<!This is determined by fitting MC and Data pi0 mass position respectivley and dividing the fits
    kNoCorrection              = 4                     //!<!Default
  };
  /// Relates string to the MC non-linearity afterburner method enumeration for %YAML configuration
  static const std::map <std::string, AliEmcalCorrectionClusterNonLinearityMCAfterburner::EMCAfterburnerMethod_t> fgkMCNonlinearityAfterburnerMethodMap; //!<!


protected:
  TString                SummarizeMCProductions(TString namePeriod);
  void                   InitNonLinearityParam(TString namePeriod, EMCAfterburnerMethod_t method);

  // Task configuration
  TH1F                  *fEnergyDistBefore;          		//!<!energy distribution before
  TH2F                  *fEnergyTimeHistBefore;      		//!<!energy/time distribution before
  TH1F                  *fEnergyDistAfter;          	 		//!<!energy distribution after
  TH2F                  *fEnergyTimeHistAfter;      			//!<!energy/time distribution after
  Float_t                fNLAfterburnerPara[9];     			///< Parameters for the non linearity function
  Int_t                  fNonLinearityAfterburnerFunction; ///< Type of function used for the non linearity afterburner correction

  EMCAfterburnerMethod_t fAfterburnerMethod;        			 ///< The version of the non-linearity afterburner correction (0-3, see enum NonlinearityMCAfterburnerMethod), 4=no correction
  TString                fMCPeriod;                  			///< The MC production name
  Bool_t                 fSetForceClusterE;          			///< Only for backwards compatibility, force cluster->E() to be set to the cluster non-linearity corrected energy. Off by default. For the standard methods, see: http://alidoc.cern.ch/AliPhysics/master/READMEcontainers.html#emcalContainerClusterEnergyCorrections

 private:
  AliEmcalCorrectionClusterNonLinearityMCAfterburner(const AliEmcalCorrectionClusterNonLinearityMCAfterburner &);               // Not implemented
  AliEmcalCorrectionClusterNonLinearityMCAfterburner &operator=(const AliEmcalCorrectionClusterNonLinearityMCAfterburner &);    // Not implemented

  // Allows the registration of the class so that it is availble to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionClusterNonLinearityMCAfterburner> reg;

  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionClusterNonLinearityMCAfterburner, 1); // EMCal cluster non-linearity MC correction afterburner component
  /// \endcond
};

#endif /* ALIEMCALCORRECTIONCLUSTERNONLINEARITYMCAFTERBURNER_H */
