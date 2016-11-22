#ifndef ALIEMCALDOWNSCALEFACTORSOCDB_H
#define ALIEMCALDOWNSCALEFACTORSOCDB_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <map>
#include <vector>
#include <TObject.h>
#include <TString.h>

/**
 * @class AliEmcalDownscaleFactorsOCDB
 * @brief Handler for downscale factors for various triggers obtained from the OCDB
 * @ingroup  EMCALCOREFW
 * @author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
 * @since Nov 22, 2016
 *
 * Handler class for downscale factors read from the OCDB. The class is used as singleton
 * class shared among several wagons. In order to access the cdb connect handler the Instance
 * function is used
 *
 * ~~~{.cxx}
 * AliEmcalDownscaleFactorsOCDB *downscalehandler = AliEmcalDownscaleFactorsOCDB::Instance();
 * downscalehandler->SetRunNumber(196310);
 * ~~~
 *
 * Accessing the downscale factor requires a full name of the trigger class, including trigger
 * cluster:
 *
 * ~~~{.cxx}
 * double ds = downscalehandler->GetDownscaleFactorForTriggerClass("CINT7-B-NOPF-ALLNOTRD");
 * ~~~
 *
 * Attention: The class does not manage OCDB access. When used in analysis, the CDB
 * connect wagon is expected to run before.
 */
class AliEmcalDownscaleFactorsOCDB : public TObject {
public:

  /**
   * Get instance of the downscale OCDB handler. If called for the
   * first time a new object is created
   * @return Downscale OCDB handler
   */
  static AliEmcalDownscaleFactorsOCDB *Instance();

  /**
   * Destructor
   */
  virtual ~AliEmcalDownscaleFactorsOCDB() {}

  /**
   * Set current run number. If the run numbers differ then
   * new downscale factors are loaded from the OCDB.
   * @param[in] runnumber New run number
   */
  void SetRun(int runnumber);

  /**
   * Get the downscale factor for a given trigger class active in the current run
   * @param[in] trigger Trigger class for which to get the downscale factor
   * @return Downscale factor for the trigger (1. if not found)
   */
  Double_t GetDownscaleFactorForTriggerClass(const TString &trigger) const;

  /**
   * Get list of trigger classes found for the given run.
   * @return
   */
  std::vector<TString> GetTriggerClasses() const;

  /**
   * Get the current run number
   * @return Current run number
   */
  Int_t GetCurrentRun() const { return fCurrentRun; }

private:
  Int_t                                       fCurrentRun;                        ///< Current run number (for which downscale factors are loaded)
  std::map<TString, Double_t>                 fDownscaleFactors;                  ///< Downscale factors for the various trigger classes for the current run
  static AliEmcalDownscaleFactorsOCDB         *fgDownscaleFactors;                ///< Singleton object

  AliEmcalDownscaleFactorsOCDB();
  AliEmcalDownscaleFactorsOCDB(const AliEmcalDownscaleFactorsOCDB &);
  AliEmcalDownscaleFactorsOCDB &operator=(const AliEmcalDownscaleFactorsOCDB &);

  /// \cond CLASSIMP
  ClassDef(AliEmcalDownscaleFactorsOCDB, 1);
  /// \endcond
};

#endif /* ALIEMCALDOWNSCALEFACTORSOCDB_H */
