#ifndef ALIEMCALCORRECTIONTASK_H
#define ALIEMCALCORRECTIONTASK_H

class AliEmcalCorrectionCellContainer;
class AliEmcalCorrectionComponent;
class AliEMCALGeometry;
class AliVEvent;

#include <AliAnalysisTaskSE.h>
#include <AliVCluster.h>

#include "AliYAMLConfiguration.h"
#include "AliEmcalContainerUtils.h"
#include "AliParticleContainer.h"
#include "AliMCParticleContainer.h"
#include "AliTrackContainer.h"
#include "AliClusterContainer.h"
#include "AliEmcalTrackSelection.h"

/**
 * @class AliEmcalCorrectionTask
 * @ingroup EMCALCORRECTIONFW
 * @brief Steering task for the EMCal correction framework
 *
 * This class is the steering class for the cell and cluster level corrections
 * for the EMCal. A %YAML configuration file is utilized to determine which
 * corrections should be run and how they should be configured. The corrections
 * are initialized by calling their Initialize() function. Similar to
 * AliAnalysisTaskEmcal, the relevant event information is loaded, and then
 * the Run() function of each correction is called.
 *
 * In general, this steering class handles all of the configuration of the
 * corrections, including passing the relevant EMCal containers and event objects.
 *
 * Note: %YAML does not play nicely with CINT and dictionary generation, so it is
 * hidden using conditional inclusion.
 *
 * @author Raymond Ehlers <raymond.ehlers@yale.edu>, Yale University
 * @author James Mulligan <james.mulligan@yale.edu>, Yale University
 * @date Jul 8, 2016
 */

class AliEmcalCorrectionTask : public AliAnalysisTaskSE {
 public:
  /**
   * @enum BeamType
   * @brief Switch for the beam type
   */
  enum BeamType {
    kNA       = -1,//!<! Undefined
    kpp       = 0, //!<! Proton-Proton
    kAA       = 1, //!<! Nucleus-Nucleus
    kpA       = 2  //!<! Proton-Nucleus
  };

  AliEmcalCorrectionTask();
  AliEmcalCorrectionTask(const char * name);
  // Implemented using copy-and-swap mechanism
  AliEmcalCorrectionTask(const AliEmcalCorrectionTask & task);
  AliEmcalCorrectionTask &operator=(AliEmcalCorrectionTask other);
  friend void swap(AliEmcalCorrectionTask & first, AliEmcalCorrectionTask & second);
  AliEmcalCorrectionTask(AliEmcalCorrectionTask && other);
  virtual ~AliEmcalCorrectionTask();

  /**
   * Initializes the Correction Task by initializing the %YAML configuration and selected correction components,
   * including setting up the input objects (cells, clusters, and tracks).
   *
   * This function is the main function for initialization and should be called from a run macro!
   * Once called, most of the configuration of the correction task and the correction components is locked in,
   * so be certain to change any additional configuration options before calling this function!
   *
   * @param[in] removeDummyTask If true, the dummy task created with the configure function is removed.
   */
  void Initialize(bool removeDummyTask = false);

  /** @{
   * @name Functions related to the %YAML configuration files.
   */
  /// Set the path to the user configuration filename
  void SetUserConfigurationFilename(std::string name) { fUserConfigurationFilename = name; }
  /// Set the path to the default configuration filename (Expert use only! The user should set the user configuration!)
  void SetDefaultConfigurationFilename(std::string name) { fDefaultConfigurationFilename = name; }
  // Print the actual configuration string
  std::ostream & PrintConfiguration(std::ostream & in, bool userConfig = false) const;
  // Write configuration to file
  bool WriteConfigurationFile(std::string filename, bool userConfig = false) const;
  // Compare configurations
  bool CompareToStoredConfiguration(std::string filename, bool userConfig = false);
  /// Retrieve the %YAML configurations for direct access
  PWG::Tools::AliYAMLConfiguration & GetYAMLConfiguration() { return fYAMLConfig; }
  /** @} */

  // Options
  // Printing
  friend std::ostream & operator<<(std::ostream &in, const AliEmcalCorrectionTask &myTask);
  void Print(Option_t* opt = "") const;
  std::ostream & Print(std::ostream &in) const;
  std::string toString(bool includeYAMLConfigurationInfo = false) const;
  // Set
  void                        SetForceBeamType(BeamType f)                          { fForceBeamType     = f                              ; }
  void                        SetNeedEmcalGeometry(Bool_t b)                        { fNeedEmcalGeom     = b                              ; }
  // Centrality options
  void                        SetUseNewCentralityEstimation(Bool_t b)               { fUseNewCentralityEstimation = b                     ; }
  void                        SetCentralityEstimator(const char * c)                { fCentEst           = c                              ; }
  virtual void                SetNCentBins(Int_t n)                                 { fNcentBins         = n                              ; }
  void                        SetCentRange(Double_t min, Double_t max)              { fMinCent           = min  ; fMaxCent = max          ; }

  /**
   * Direct access to the correction components.
   *
   * \return std::vector of the correction components. Using this vector, the components and their settings
   * can be modified as desired. However, keep in mind that whatever changes are made here will _NOT_ be
   * reflected in the stored %YAML configuration.
   */
  const std::vector<AliEmcalCorrectionComponent *> & CorrectionComponents() { return fCorrectionComponents; }
  AliEmcalCorrectionComponent * GetCorrectionComponent(const std::string & name) const;

  // Containers and cells
  AliParticleContainer       *AddParticleContainer(const char *n)                   { return AliEmcalContainerUtils::AddContainer<AliParticleContainer>(n, fParticleCollArray); }
  AliTrackContainer          *AddTrackContainer(const char *n)                      { return AliEmcalContainerUtils::AddContainer<AliTrackContainer>(n, fParticleCollArray); }
  AliMCParticleContainer     *AddMCParticleContainer(const char *n)                 { return AliEmcalContainerUtils::AddContainer<AliMCParticleContainer>(n, fParticleCollArray); }
  AliClusterContainer        *AddClusterContainer(const char *n)                    { return AliEmcalContainerUtils::AddContainer<AliClusterContainer>(n, fClusterCollArray); }
  void                        AdoptParticleContainer(AliParticleContainer* cont)    { fParticleCollArray.Add(cont)                        ; }
  void                        AdoptTrackContainer(AliTrackContainer* cont)          { AdoptParticleContainer(cont)                        ; }
  void                        AdoptMCParticleContainer(AliMCParticleContainer* cont){ AdoptParticleContainer(cont)                        ; }
  void                        AdoptClusterContainer(AliClusterContainer* cont)      { fClusterCollArray.Add(cont)                         ; }
  AliParticleContainer       *GetParticleContainer(Int_t i=0)                 const { return AliEmcalContainerUtils::GetContainer<AliParticleContainer>(i, fParticleCollArray); }
  AliParticleContainer       *GetParticleContainer(const char* name)          const { return AliEmcalContainerUtils::GetContainer<AliParticleContainer>(name, fParticleCollArray); }
  AliClusterContainer        *GetClusterContainer(Int_t i=0)                  const { return AliEmcalContainerUtils::GetContainer<AliClusterContainer>(i, fClusterCollArray); }
  AliClusterContainer        *GetClusterContainer(const char* name)           const { return AliEmcalContainerUtils::GetContainer<AliClusterContainer>(name, fClusterCollArray); }
  AliMCParticleContainer     *GetMCParticleContainer(Int_t i=0)               const { return dynamic_cast<AliMCParticleContainer*>(GetParticleContainer(i))   ; }
  AliMCParticleContainer     *GetMCParticleContainer(const char* name)        const { return dynamic_cast<AliMCParticleContainer*>(GetParticleContainer(name)); }
  AliTrackContainer          *GetTrackContainer(Int_t i=0)                    const { return dynamic_cast<AliTrackContainer*>(GetParticleContainer(i))        ; }
  AliTrackContainer          *GetTrackContainer(const char* name)             const { return dynamic_cast<AliTrackContainer*>(GetParticleContainer(name))     ; }
  void                        RemoveParticleContainer(Int_t i=0)                    { fParticleCollArray.RemoveAt(i)                      ; }
  void                        RemoveClusterContainer(Int_t i=0)                     { fClusterCollArray.RemoveAt(i)                       ; }
  // Cells
  AliEmcalCorrectionCellContainer *GetCellContainer(const std::string & cellsContainerName) const;

  // Methods from AliAnalysisTaskSE
  void UserCreateOutputObjects();
  void UserExec(Option_t * option);
  Bool_t UserNotify();

  // Aditional steering functions
  virtual void ExecOnce();
  virtual Bool_t Run();

  /**
   * EMCal Correction Task AddTask. Should be used by most users, except for those on the LEGO train
   * (see below).
   *
   * @param[in] suffix Suffix string used to select components in a %YAML configuration
   *
   * @return A new EMCal Correction Task added to the analysis manager and ready to configure.
   */
  static AliEmcalCorrectionTask* AddTaskEmcalCorrectionTask(TString suffix = "");
  /**
   * Retrieve an existing correction task by name to perform further configuration. This should
   * _ONLY_ be used on the LEGO train. The suffix passed here must be unique to identify a user.
   *
   * To achieve this, a dummy task is created when the configure task is called because AliAnalysisTaskCfg
   * requires that all wagons add a task. Then, when Initialize(true) is called on the correction task, the
   * dummy task is removed. This is a hack, but is required to work around constraints in AliAnalysisTaskCfg.
   *
   * @param[in] suffix Suffix string used to uniquely identify a user and find the corresponding correction task. If using a suffix with the correction task, the suffixes must match.
   *
   * @return An existing (usually unconfigured) EMCal Correction Task which was retrieved from the analysis manager.
   */
  static AliEmcalCorrectionTask* ConfigureEmcalCorrectionTaskOnLEGOTrain(TString suffix);

 private:
  // Utility functions
  // Cell utilities
  void SetCellsObjectInCellContainerBasedOnProperties(AliEmcalCorrectionCellContainer * cellContainer);
  // Container utilities
  void CheckForContainerArray(AliEmcalContainer * cont, AliEmcalContainerUtils::InputObject_t objectType);
  // %YAML configuration utilities
  std::string GetInputFieldNameFromInputObjectType(AliEmcalContainerUtils::InputObject_t inputObjectType);
  bool CheckPossibleNamesForComponentName(std::string & name, std::set <std::string> & possibleComponents);
  // General utilities
  BeamType GetBeamType() const;
  void PrintRequestedContainersInformation(AliEmcalContainerUtils::InputObject_t inputObjectType, std::ostream & stream) const;
  // LEGO Train utilities
  void RemoveDummyTask() const;

  // Retrieve objects in event
  Bool_t RetrieveEventObjects();

  // Execute component functions
  void UserCreateOutputObjectsComponents();
  void ExecOnceComponents();

  // Initialization functions
  void InitializeConfiguration();
  void DetermineComponentsToExecute(std::vector <std::string> & componentsToExecute);
  void CheckForUnmatchedUserSettings();
  void InitializeComponents();

  // Input objects (Cells, Clusters, Tracks) functions
  void CreateInputObjects(AliEmcalContainerUtils::InputObject_t inputObjectType);
  void AddContainersToComponent(AliEmcalCorrectionComponent * component, AliEmcalContainerUtils::InputObject_t inputObjectType, bool checkObjectExists = false);
  
  // Hidden from CINT since it cannot handle %YAML objects well
  // Input objects 
  void SetupContainersFromInputNodes(AliEmcalContainerUtils::InputObject_t inputObjectType, std::set <std::string> & requestedContainers);
  // Cells
  void SetupCellsInfo(std::string containerName);
  // Containers
  void SetupContainer(const AliEmcalContainerUtils::InputObject_t inputObjectType, const std::string containerName);
  AliEmcalContainer * AddContainer(const AliEmcalContainerUtils::InputObject_t contType, const std::string containerName);

  // Utilities
  // %YAML node dependent initialization utlitiles
  void GetPropertyNamesFromNode(const std::string configurationName, const std::string componentName, std::set <std::string> & propertyNames, const bool nodeRequired);

  PWG::Tools::AliYAMLConfiguration fYAMLConfig;            ///< Handles configuration from YAML.

  std::string                 fSuffix;                     ///< Suffix of the Correction Task (used to select specialized components)

  std::string                 fUserConfigurationFilename;  //!<! User %YAML configruation filename
  std::string                 fDefaultConfigurationFilename; //!<! Default %YAML configuration filename

  std::vector <std::string>   fOrderedComponentsToExecute; ///< Ordered set of components to execute
  std::vector <AliEmcalCorrectionComponent *> fCorrectionComponents; ///< Contains the correction components
  bool                        fConfigurationInitialized;   ///< True if the %YAML configuration files are initialized

  bool                        fIsEsd;                      ///< File type
  bool                        fEventInitialized;           ///< If the event is initialized properly
  Double_t                    fCent;                       //!<! Event centrality
  Int_t                       fCentBin;                    //!<! Event centrality bin
  Double_t                    fMinCent;                    ///< min centrality for event selection
  Double_t                    fMaxCent;                    ///< max centrality for event selection
  Int_t                       fNcentBins;                  ///< how many centrality bins
  TString                     fCentEst;                    ///< name of V0 centrality estimator
  Bool_t                      fUseNewCentralityEstimation; ///< Use new centrality estimation (for 2015 data)
  Double_t                    fVertex[3];                  //!<! Event vertex
  Int_t                       fNVertCont;                  //!<! Event vertex number of contributors
  BeamType                    fBeamType;                   //!<! Event beam type
  BeamType                    fForceBeamType;              ///< forced beam type
  Bool_t                      fNeedEmcalGeom;              ///< whether or not the task needs the emcal geometry
  AliEMCALGeometry           *fGeom;                       //!<! Emcal geometry

  TObjArray                   fParticleCollArray;          ///< Particle/track collection array
  TObjArray                   fClusterCollArray;           ///< Cluster collection array
  std::vector <AliEmcalCorrectionCellContainer *> fCellCollArray; ///< Cells collection array
  
  TList *                     fOutput;                     //!<! Output for histograms

  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionTask, 9); // EMCal correction task
  /// \endcond
};

/**
 * @class AliEmcalCorrectionCellContainer
 * @ingroup EMCALCOREFW
 * @brief Wrapper around cells objects for the EMCal Correction Task
 *
 * This class is a container around the cells collections that are used in the
 * EMCal Correction Task. In particular, it keeps track of the cells, the associated
 * names, as well as whether the cells are embedded.
 *
 * Note that because the CaloCells objects are techincally wrappers around the
 * cells collections, this is something of a wrapper around a wrapper. However,
 * due to the difficultly of modifying the base classes, this very simple class
 * will suffice for the EMCal Correction Task.
 *
 * @author Raymond Ehlers <raymond.ehlers@yale.edu>, Yale University
 * @author James Mulligan <james.mulligan@yale.edu>, Yale University
 * @date Nov 8, 2016
 */
class AliEmcalCorrectionCellContainer {
 public:
  AliEmcalCorrectionCellContainer():
    fBranchName(""),
    fName(""),
    fIsEmbedding(""),
    fCells(0)
  {}
  AliEmcalCorrectionCellContainer(std::string branchName, std::string name, std::string branchToCopyName, bool isEmbedded):
    fBranchName(branchName),
    fName(name),
    fIsEmbedding(isEmbedded),
    fCells(0)
  {}
  virtual ~AliEmcalCorrectionCellContainer() {}

  /// Get the name of the cells branch (NOT the same as the name!)
  std::string GetBranchName() const { return fBranchName; }
  /// Get the name of the cells object (NOT the same as the branch!)
  std::string GetName() const { return fName; }
  /// True if the cells are located in the event that is being embedded
  bool GetIsEmbedding() const { return fIsEmbedding; }
  /// Pointer to the actual CaloCells object
  AliVCaloCells * GetCells() const { return fCells; }

  /// Set the name of the cells branch (NOT the same as the name!)
  void SetBranchName(std::string branchName) { fBranchName = branchName; }
  /// Set the name of the cells object (NOT the same as the branch!)
  void SetName(std::string name ) { fName = name; }
  /// Set to true if the cells are located in the event that is being embedded
  void SetIsEmbedding(bool isEmbedded) { fIsEmbedding = isEmbedded; }
  /// Sets the Pointer to the actual CaloCells object
  void SetCells(AliVCaloCells * cells) { fCells = cells; }

 protected:
  std::string fBranchName;                        ///< Name of the cells branch
  std::string fName;                              ///< Name of the cells object
  bool fIsEmbedding;                              ///< Whether the cells should be taken from an external file (embedded)
  AliVCaloCells *fCells;                          //!<! The actual cells object associated with this infomration

  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionCellContainer, 1); // EMCal correction cell container
  /// \endcond
};

#endif /* ALIEMCALCORRECTIONTASK_H */ 
