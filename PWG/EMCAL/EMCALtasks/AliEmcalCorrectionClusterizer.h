#ifndef ALIEMCALCORRECTIONCLUSTERIZER_H
#define ALIEMCALCORRECTIONCLUSTERIZER_H

#include "AliEmcalCorrectionComponent.h"

#include "AliEMCALRecParam.h"

class TStopwatch;

/**
 * @class AliEmcalCorrectionClusterizer
 * @ingroup EMCALCOREFW
 * @brief EMCal clusterizer component in the EMCal correction framework.
 * 
 * Clusterizes a collection of cells into a collection of clusters.
 *
 * For some datasets one needs to re-run the clusterizer from the cells. This is dataset-specific. Analyzers that are unfamiliar with a specific dataset should communicate with EMCal/DCal detector experts to determine whether the basic corrections and the bad channel map were already available and used in a specific ESD/AOD production. 
 *
 * The clusterizer type and the time cuts are to be chosen appropriately for each dataset. Usually the v1 clusterizer is used for pp and the v2 clusterizer is used for PbPb. Sometimes for pp reference runs with the same collision energy as PbPb the v2 clusterizer is employed, but this is analysis dependent. EMCal detector experts are to be contacted for the time cuts.
 *
 * The clusterizer will use as input the cell branch specified in the YAML config, and as output will rewrite the cluster branch specified in the YAML config.
 *
 * At this point the energy of the cluster will be available through `cluster->E()` where cluster is the pointer to the AliAODCaloCluster or AliESDCaloCluster object.
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
  /// Relates string to the clusterizer type enumeration for YAML configuration
  static const std::map <std::string, AliEMCALRecParam::AliEMCALClusterizerFlag> fgkClusterizerTypeMap; //!<!

  AliEmcalCorrectionClusterizer();
  virtual ~AliEmcalCorrectionClusterizer();

  // Sets up and runs the task
  Bool_t Initialize();
  void UserCreateOutputObjects();
  Bool_t Run();
  Bool_t CheckIfRunChanged();
  
protected:
  void           Clusterize();
  void           FillDigitsArray();
  void           Init();
  void           RecPoints2Clusters(TClonesArray *clus);
  void           UpdateClusters();
  void           CalibrateClusters();
  
  void           RemapMCLabelForAODs(Int_t &label);
  void           SetClustersMCLabelFromOriginalClusters();
  void           ClearEMCalClusters();
  
  TH1F* fHistCPUTime;                                     //!<! CPU time for the Run() function (event loop)
  TH1F* fHistRealTime;                                    //!<! Real time for the Run() function (event loop)
  TStopwatch * fTimer;                                    //!<! Timer for the Run() function (event loop)

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
  Bool_t                 fTestPatternInput;               ///< Use test pattern as input instead of cells
  
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
  ClassDef(AliEmcalCorrectionClusterizer, 4); // EMCal correction clusterizer component
  /// \endcond
};

#endif /* ALIEMCALCORRECTIONCLUSTERIZER_H */
