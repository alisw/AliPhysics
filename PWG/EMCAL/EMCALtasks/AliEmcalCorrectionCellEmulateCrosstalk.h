#ifndef ALIEMCALCORRECTIONCELLEMULATECROSSTALK_H
#define ALIEMCALCORRECTIONCELLEMULATECROSSTALK_H

#include <TRandom3.h>

#include "AliEmcalCorrectionComponent.h"

/**
 * @class AliEmcalCorrectionCellEmulateCrosstalk
 * @ingroup EMCALCORRECTIONFW
 * @brief Correction component to emulate cell-level crosstalk in the EMCal correction framework.
 *
 * Performs energy smearing of cells to mimic T-card crosstalk. The original cell information in the event **will be overwritten**.
 *
 * Based on code in AliAnalysisTaskEMCALClusterize.
 *
 * @author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-Grenoble, AliAnalysisTaskEMCALClusterize
 * @author Raymond Ehlers <raymond.ehlers@yale.edu>, Yale University, centralize EMCal corrections using components
 * @author James Mulligan <james.mulligan@yale.edu>, Yale University, centralize EMCal corrections using components
 * @date Feb 11 2018
 */

class AliEmcalCorrectionCellEmulateCrosstalk : public AliEmcalCorrectionComponent {
 public:
  AliEmcalCorrectionCellEmulateCrosstalk();
  virtual ~AliEmcalCorrectionCellEmulateCrosstalk();

  // Sets up and runs the task
  Bool_t Initialize();
  void UserCreateOutputObjects();
  Bool_t Run();
  
  // Print correction paramters
  void PrintTCardParam();
  
 protected:

  // Useful constants
  static const Int_t fgkNEMCalCells = 17664;       ///< Total number of cells in the calorimeter, 10*48*24 (EMCal) + 4*48*8 (EMCal/DCal 1/3) + 6*32*24 (DCal)
  static const Int_t fgkNsm = 22;                  ///< Total number of super-modules

  // T-Card correlation emulation, do on MC
  void           MakeCellTCardCorrelation();
  void           CalculateInducedEnergyInTCardCell(Int_t absId, Int_t absIdRef, 
                                                   Int_t sm, Float_t ampRef, 
                                                   Int_t cellCase) ;

  void           AddInducedEnergiesToExistingCells();
  void           AddInducedEnergiesToNewCells();
  
  virtual void   ResetArrays();
  Bool_t         AcceptCell(Int_t absID);

  // Initialize array parameters.
  #if !(defined(__CINT__) || defined(__MAKECINT__))
  template <typename T>
  void RetrieveAndSetProperties(const T & properties);
  #endif
  void SetProperty(Float_t val[][fgkNsm], std::vector<double> & property, unsigned int iSM, const std::string & name);
  void SetProperty(Float_t val[fgkNsm], std::vector<double> & property, unsigned int iSM, const std::string & name);
  
  TH1F* fCellEnergyDistBefore;        //!<! cell energy distribution, before energy smearing
  TH1F* fCellEnergyDistAfter;         //!<! cell energy distribution, after energy smearing
  
  // T-Card correlation emulation, do on MC
  Bool_t                fTCardCorrClusEnerConserv; ///< When making correlation, subtract from the reference cell the induced energy on the neighbour cells
  Float_t               fTCardCorrCellsEner[fgkNEMCalCells]; ///<  Array with induced cell energy in T-Card neighbour cells
  Bool_t                fTCardCorrCellsNew [fgkNEMCalCells]; ///<  Array with induced cell energy in T-Card neighbour cells, that before had no signal
  
  Float_t               fTCardCorrInduceEner         [4][fgkNsm]; ///< Induced energy loss gauss constant on 0-same row, diff col, 1-up/down cells left/right col 2-left/righ col, and 2nd row cells, param 0
  Float_t               fTCardCorrInduceEnerFrac     [4][fgkNsm]; ///< Induced energy loss gauss fraction param0 on 0-same row, diff col, 1-up/down cells left/right col 2-left/righ col, and 2nd row cells, param 0
  Float_t               fTCardCorrInduceEnerFracP1   [4][fgkNsm]; ///< Induced energy loss gauss fraction param1 on 0-same row, diff col, 1-up/down cells left/right col 2-left/righ col, and 2nd row cells, param1
  Float_t               fTCardCorrInduceEnerFracWidth[4][fgkNsm]; ///< Induced energy loss gauss witdth on 0-same row, diff col, 1-up/down cells left/right col 2-left/righ col, and 2nd row cells
  Float_t               fTCardCorrInduceEnerFracMax[fgkNsm];   ///< In case fTCardCorrInduceEnerFracP1  is non null, restrict the maximum fraction of induced energy per SM
  Float_t               fTCardCorrInduceEnerFracMin[fgkNsm];   ///< In case fTCardCorrInduceEnerFracP1  is non null, restrict the minimum fraction of induced energy per SM
  Float_t               fTCardCorrInduceEnerProb[fgkNsm];      ///< Probability to induce energy loss per SM
  
  TRandom3              fRandom;                   ///<  Random generator
  Bool_t                fRandomizeTCard;           ///<  Use random induced energy
  
  Float_t               fTCardCorrMinAmp;          ///<  Minimum cell energy to induce signal on adjacent cells
  Float_t               fTCardCorrMinInduced;      ///<  Minimum induced energy signal on adjacent cells, sum of induced plus original energy, use same as cell energy clusterization cut
  Float_t               fTCardCorrMaxInducedELeak; ///<  Maximum value of induced energy signal that is always leaked, ~5-10 MeV
  Float_t               fTCardCorrMaxInduced;      ///<  Maximum induced energy signal on adjacent cells
  
  Bool_t                fPrintOnce;                ///< Print once analysis parameters
  
  AliAODCaloCells      *fAODCellsTmp;              //!<! Temporal array of cells copy

  
private:
  
  AliEmcalCorrectionCellEmulateCrosstalk(const AliEmcalCorrectionCellEmulateCrosstalk &);               // Not implemented
  AliEmcalCorrectionCellEmulateCrosstalk &operator=(const AliEmcalCorrectionCellEmulateCrosstalk &);    // Not implemented

  // Allows the registration of the class so that it is availble to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionCellEmulateCrosstalk> reg;

  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionCellEmulateCrosstalk, 3); // EMCal cell crosstalk emulation correction component
  /// \endcond
};

#if !(defined(__CINT__) || defined(__MAKECINT__))
/**
 * Assigns array properties defined in the YAML config to initialize array class members.
 * Can handle any map for which there is an appropirately defined SetProperty(...).
 *
 * NOTE: This function needs to be hidden from CINT because it uses templates.
 *
 * @param[in] properties YAML to class members map
 */
template <typename T>
void AliEmcalCorrectionCellEmulateCrosstalk::RetrieveAndSetProperties(const T & properties)
{
  std::string taskName = GetName();
  std::map <std::string, std::vector<double>> values;

  for (auto& val : properties) {
    AliDebugStream(1) << "Processing value " << val.first << "\n";
    bool enabled = false;
    fYAMLConfig.GetProperty({taskName, val.first, "enabled"}, enabled, true);
    if (enabled) {
      AliDebugStream(1) << val.first << " enabled.\n";
      fYAMLConfig.GetProperty({taskName, val.first, "values"}, values, true);
      // "all"" is prioritized over any individual SM values
      auto property = values.find("all");
      if (property != values.end()) {
        // Set same value in all SMs
        AliDebugStream(1) << "Retrieving all SM settings for property " << val.first << "\n";
        for (unsigned int iSM = 0; iSM < fgkNsm; iSM++) {
          SetProperty(val.second, property->second, iSM, val.first);
        }
      }
      else {
        // Handle per SM values
        AliDebugStream(1) << "Retrieving per SM settings for property " << val.first << "\n";
        for (auto && property : values) {
          unsigned int iSM = std::stoul(property.first);
          // iSM must be >= 0 because it is unsigned
          if (iSM < fgkNsm) {
            SetProperty(val.second, property.second, iSM, val.first);
          }
          else {
            AliWarningStream() << "SM " << iSM << " requested for property " << val.first << " is out of range. Please check your configuration!\n";
          }
        }
      }
    }
  }
}
#endif

#endif /* ALIEMCALCORRECTIONCellEmulateCrosstalk_H */
