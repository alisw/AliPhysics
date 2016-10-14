#ifndef ALIEMCALCORRECTIONCLUSTERIZER_H
#define ALIEMCALCORRECTIONCLUSTERIZER_H

#include "AliEmcalCorrectionComponent.h"

#include "AliEMCALRecParam.h"

/**
 * @class AliEmcalCorrectionClusterizer
 * @brief EMCal clusterizer component in the EMCal correction framework
 * @ingroup EMCALCOREFW
 * 
 * Based on code in AliAnalysisTaskEMCALClusterizeFast, in turn based on code by Deepa Thomas.
 *
 * @author Constantin Loizides, LBNL, AliAnalysisTaskEMCALClusterizeFast
 * @author Salvatore Aiola, LBNL, AliAnalysisTaskEMCALClusterizeFast
 * @author James Mulligan <james.mulligan@yale.edu>, Yale University, centralize EMCal corrections using components
 * @author Raymond Ehlers <raymond.ehlers@yale.edu>, Yale University, centralize EMCal corrections using components
 * @date Jul 8, 2016
 */

class AliEmcalCorrectionClusterizer : public AliEmcalCorrectionComponent {
 public:
  enum InputCellType {
    kFEEData = 0,
    kFEEDataMCOnly,
    kFEEDataExcludeMC,
    kPattern,
    kL0FastORs,
    kL0FastORsTC,
    kL1FastORs
  };

#if !(defined(__CINT__) || defined(__MAKECINT__))
  std::map <std::string, AliEMCALRecParam::AliEMCALClusterizerFlag> clusterizerTypeMap = {
    {"kClusterizerv1", AliEMCALRecParam::kClusterizerv1 },
    {"kClusterizerNxN", AliEMCALRecParam::kClusterizerNxN },
    {"kClusterizerv2", AliEMCALRecParam::kClusterizerv2 },
    {"kClusterizerFW", AliEMCALRecParam::kClusterizerFW }
  };
#endif

  AliEmcalCorrectionClusterizer();
  virtual ~AliEmcalCorrectionClusterizer();

  // Sets up and runs the task
  Bool_t Initialize();
  void UserCreateOutputObjects();
  Bool_t Run();
  
protected:
  void           Clusterize();
  void           FillDigitsArray();
  void           Init();
  void           RecPoints2Clusters(TClonesArray *clus);
  void           UpdateClusters();
  void           CalibrateClusters();
  
  void           RemapMCLabelForAODs(Int_t &label);
  void           SetClustersMCLabelFromOriginalClusters();
  
  TClonesArray          *fDigitsArr;                      //!<!digits array
  TObjArray             *fClusterArr;                     //!<!recpoints array
  AliEMCALRecParam      *fRecParam;                       ///< reconstruction parameters container
  AliEMCALClusterizer   *fClusterizer;                    //!<!clusterizer
  AliEMCALAfterBurnerUF *fUnfolder;                       //!<!unfolding procedure
  Bool_t                 fJustUnfold;                     ///< just unfold, do not recluster
  TString                fGeomName;                       ///< name of geometry to use.
  Bool_t                 fGeomMatrixSet;                  ///< set geometry matrices only once, for the first event.
  Bool_t                 fLoadGeomMatrices;               ///< matrices from configuration, not geometry.root nor ESDs/AODs
  TGeoHMatrix           *fGeomMatrix[AliEMCALGeoParams::fgkEMCALModules]; ///< geometry matrices with alignments
  TString                fOCDBpath;                       ///< path with OCDB location
  AliEMCALCalibData     *fCalibData;                      ///< EMCAL calib data
  AliCaloCalibPedestal  *fPedestalData;                   ///< EMCAL pedestal
  Bool_t                 fLoadCalib;                      ///< access calib object from OCDB (def=off)
  Bool_t                 fLoadPed;                        ///< access ped object from OCDB (def=off)
  Bool_t                 fSubBackground;                  ///< subtract background if true (def=off)
  Int_t                  fNPhi;                           ///< nPhi (for FixedWindowsClusterizer)
  Int_t                  fNEta;                           ///< nEta (for FixedWinoswsClusterizer)
  Int_t                  fShiftPhi;                       ///< shift in phi (for FixedWindowsClusterizer)
  Int_t                  fShiftEta;                       ///< shift in eta (for FixedWindowsClusterizer)
  Bool_t                 fTRUShift;                       ///< shifting inside a TRU (true) or through the whole calorimeter (false) (for FixedWindowsClusterizer)
  InputCellType          fInputCellType;                  ///< input cells type to make clusters
  
  // MC labels
  static const Int_t     fgkTotalCellNumber = 17664 ;     ///< Maximum number of cells in EMCAL/DCAL: (48*24)*(10+4/3.+6*2/3.)
  
  Int_t                  fOrgClusterCellId[fgkTotalCellNumber]; ///< Array ID of cluster to which the cell belongs in unmodified clusters.
  Int_t                  fCellLabels      [fgkTotalCellNumber]; ///< Array with MC label/map
  
  Int_t                  fSetCellMCLabelFromCluster;      ///< Use cluster MC label as cell label:
  // 0 - get the MC label stored in cells (not available for productions done with aliroot < v5-02-Rev09)
  // 1 - assign to the cell the MC label of the cluster
  
  Bool_t                 fSetCellMCLabelFromEdepFrac;     ///< For MC generated with aliroot > v5-07-21, check the EDep information
  // stored in ESDs/AODs to set the cell MC labels
  
  Bool_t                 fRemapMCLabelForAODs;            ///< Remap AOD cells MC label
  
  Bool_t                 fRecalDistToBadChannels;         ///< recalculate distance to bad channel
  Bool_t                 fRecalShowerShape;               ///< switch for recalculation of the shower shape
  
  TClonesArray          *fCaloClusters;                   //!<!calo clusters array
  AliESDEvent           *fEsd;                            //!<!esd event
  AliAODEvent           *fAod;                            //!<!aod event

 private:
  AliEmcalCorrectionClusterizer(const AliEmcalCorrectionClusterizer &);               // Not implemented
  AliEmcalCorrectionClusterizer &operator=(const AliEmcalCorrectionClusterizer &);    // Not implemented

  // Allows the registration of the class so that it is availble to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionClusterizer> reg;

  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionClusterizer, 1); // EMCal correction clusterizer component
  /// \endcond
};

#endif /* ALIEMCALCORRECTIONCLUSTERIZER_H */
