#ifndef ALIEMCALCORRECTIONTASK_H
#define ALIEMCALCORRECTIONTASK_H

// CINT can't handle the yaml header!
#if !(defined(__CINT__) || defined(__MAKECINT__))
#include <yaml-cpp/yaml.h>
#endif

class AliEmcalCorrectionCellContainer;
class AliEmcalCorrectionComponent;
class AliEMCALGeometry;
class AliVEvent;

#include "AliAnalysisTaskSE.h"
#include "AliParticleContainer.h"
#include "AliMCParticleContainer.h"
#include "AliTrackContainer.h"
#include "AliClusterContainer.h"
#include "AliVCluster.h"
#include "AliEmcalTrackSelection.h"

/**
 * @class AliEmcalCorrectionTask
 * @brief Steering task for the EMCal correction framework
 * @ingroup EMCALCOREFW
 *
 * This class is the steering class for the cell and cluster level corrections
 * for the EMCal. A YAML configuration file is utilized to determine which
 * corrections should be run and how they should be configured. The corrections
 * are initialized by calling their Initialize() function. Similar to
 * AliAnalysisTaskEmcal, the relevant event information is loaded, and then
 * the Run() function of each correction is called.
 *
 * In general, this steering class handles all of the configuration of the
 * corrections, including passing the relevant EMCal containers and event objects.
 *
 * Note: YAML does not play nicely with CINT and dictionary generation, so it is
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

  /** 
   * @enum InputObject_t
   * @brief Type of input object to be created
   */
  enum InputObject_t {
    kNoDefinedInputObject = -1,    //!< Not initialied type
    kCaloCells = 0,                //!< Calo cells
    kCluster = 1,                  //!< Cluster container
    kTrack = 2,                    //!< Track container
  };

#if !(defined(__CINT__) || defined(__MAKECINT__))
  std::map <std::string, AliVCluster::VCluUserDefEnergy_t> clusterEnergyTypeMap = {
    {"kNonLinCorr", AliVCluster::kNonLinCorr },
    {"kHadCorr", AliVCluster::kHadCorr },
    {"kUserDefEnergy1", AliVCluster::kUserDefEnergy1 },
    {"kUserDefEnergy2", AliVCluster::kUserDefEnergy2 }
  };

  std::map <std::string, AliEmcalTrackSelection::ETrackFilterType_t> trackFilterTypeMap = {
    {"kNoTrackFilter", AliEmcalTrackSelection::kNoTrackFilter },
    {"kCustomTrackFilter", AliEmcalTrackSelection::kCustomTrackFilter },
    {"kHybridTracks",  AliEmcalTrackSelection::kHybridTracks },
    {"kTPCOnlyTracks", AliEmcalTrackSelection::kTPCOnlyTracks }
  };
#endif

  AliEmcalCorrectionTask();
  AliEmcalCorrectionTask(const char * name);
  virtual ~AliEmcalCorrectionTask();

  // Methods from AliAnalysisTaskSE
  void UserCreateOutputObjects();
  void UserExec(Option_t * option);
  Bool_t UserNotify();

  // Functions that can be overloaded by users
  /**
   * Executed once per event to initialize objects. The user probably does not need to overload it.
   */
  virtual void ExecOnce();
  /**
   * Executed each event once all objects are available
   */
  virtual Bool_t Run();

  // Options
  void SetUserConfigurationFilename(std::string name) { fUserConfigurationFilename = name; }
  void SetDefaultConfigurationFilename(std::string name) { fDefaultConfigurationFilename = name; }

  void InitializeConfiguration();

  // Containers and cells
  AliParticleContainer       *AddParticleContainer(const char *n);
  AliTrackContainer          *AddTrackContainer(const char *n);
  AliMCParticleContainer     *AddMCParticleContainer(const char *n);
  AliClusterContainer        *AddClusterContainer(const char *n);
  void                        AdoptParticleContainer(AliParticleContainer* cont)    { fParticleCollArray.Add(cont)                        ; }
  void                        AdoptTrackContainer(AliTrackContainer* cont)          { AdoptParticleContainer(cont)                        ; }
  void                        AdoptMCParticleContainer(AliMCParticleContainer* cont){ AdoptParticleContainer(cont)                        ; }
  void                        AdoptClusterContainer(AliClusterContainer* cont)      { fClusterCollArray.Add(cont)                         ; }
  AliParticleContainer       *GetParticleContainer(Int_t i=0)         const;
  AliParticleContainer       *GetParticleContainer(const char* name)  const;
  AliClusterContainer        *GetClusterContainer(Int_t i=0)          const;
  AliClusterContainer        *GetClusterContainer(const char* name)   const;
  AliMCParticleContainer     *GetMCParticleContainer(Int_t i=0)               const { return dynamic_cast<AliMCParticleContainer*>(GetParticleContainer(i))   ; }
  AliMCParticleContainer     *GetMCParticleContainer(const char* name)        const { return dynamic_cast<AliMCParticleContainer*>(GetParticleContainer(name)); }
  AliTrackContainer          *GetTrackContainer(Int_t i=0)                    const { return dynamic_cast<AliTrackContainer*>(GetParticleContainer(i))        ; }
  AliTrackContainer          *GetTrackContainer(const char* name)             const { return dynamic_cast<AliTrackContainer*>(GetParticleContainer(name))     ; }
  void                        RemoveParticleContainer(Int_t i=0)                    { fParticleCollArray.RemoveAt(i)                      ; } 
  void                        RemoveClusterContainer(Int_t i=0)                     { fClusterCollArray.RemoveAt(i)                       ; } 
  AliVCaloCells              *GetCellsFromContainerArray(const std::string & cellsContainerName) const;

  void                        SetForceBeamType(BeamType f)                          { fForceBeamType     = f                              ; }
  void                        SetRunPeriod(const char* runPeriod)                   { fRunPeriod = runPeriod; fRunPeriod.ToLower(); }
  const TString &             GetRunPeriod()                                  const { return fRunPeriod; }

  void                        SetUseNewCentralityEstimation(Bool_t b)               { fUseNewCentralityEstimation = b                     ; }

  const std::vector<AliEmcalCorrectionComponent *> & CorrectionComponents() { return fCorrectionComponents; }

  bool WriteConfigurationFile(std::string filename, bool userConfig = false);

  // Determine branch name using the "usedefault" pattern
  static std::string DetermineUseDefaultName(InputObject_t contType, bool esdMode, bool returnObjectType = false);

  // Get the proper event based on whether embedding is enabled or not
  static AliVEvent * GetEvent(AliVEvent * inputEvent, bool isEmbedding = false);

 private:
  AliEmcalCorrectionTask(const AliEmcalCorrectionTask &);             // Not implemented
  AliEmcalCorrectionTask &operator=(const AliEmcalCorrectionTask &);  // Not implemented

  static inline bool DoesFileExist(const std::string & filename);
  void SetupConfigurationFilePath(std::string & filename, bool userFile = false);

  void RetrieveExecutionOrder(std::vector <std::string> & componentsToAdd);
  void InitializeComponents();

  void SetCellsObjectInCellContainerBasedOnProperties(AliEmcalCorrectionCellContainer * cellContainer);
  void CheckForContainerArray(AliEmcalContainer * cont, InputObject_t objectType);

  std::string GetInputFieldNameFromInputObjectType(InputObject_t inputObjectType);
  void CreateInputObjects(InputObject_t inputObjectType);
  void AddContainersToComponent(AliEmcalCorrectionComponent * component, InputObject_t inputObjectType);
#if !(defined(__CINT__) || defined(__MAKECINT__))
  void SetupContainersFromInputNodes(InputObject_t inputObjectType, YAML::Node & userInputObjectNode, YAML::Node & defaultInputObjectNode, std::set <std::string> & requestedContainers);
  void GetNodeForInputObjects(YAML::Node & inputNode, YAML::Node & nodeToRetrieveFrom, std::string & inputObjectName, bool requiredProperty);

  void SetupCellsInfo(std::string containerName, YAML::Node & userNode, YAML::Node & defaultNode);
  void SetupContainer(InputObject_t inputObjectType, std::string containerName, YAML::Node & userNode, YAML::Node & defaultNode);

  AliEmcalContainer * AddContainer(InputObject_t contType, std::string & containerName, YAML::Node & userNode, YAML::Node & defaultNode);
#endif

  Bool_t                      RetrieveEventObjects();
  void                        ExecOnceComponents();

  BeamType                    GetBeamType();

#if !(defined(__CINT__) || defined(__MAKECINT__))
  YAML::Node              fUserConfiguration;             /// User YAML Configuration
  YAML::Node              fDefaultConfiguration;          /// Default YAML Configuration
#endif

  std::string             fUserConfigurationString;       ///< Store the User YAML configuration as a string so that it can be streamed
  std::string             fDefaultConfigurationString;    ///< Store the default YAML configuration as a string so that it can be streamed

  std::string fUserConfigurationFilename;                 //!<! User YAML configruation filename
  std::string fDefaultConfigurationFilename;              //!<! Default YAML configuration filename

  std::vector <AliEmcalCorrectionComponent *> fCorrectionComponents; //!<! Contains the correction components
  bool fConfigurationInitialized;                         ///< True if the YAML files are initialized

  bool                        fIsEsd;                      ///< File type
  TString                     fRunPeriod;                  ///< Run period (passed by user)
  bool                        fEventInitialized;           ///< If the event is initialized properly
  Double_t                    fCent;                       //!<! Event centrality
  Int_t                       fCentBin;                    //!<! Event centrality bin
  Double_t                    fMinCent;                    ///< min centrality for event selection
  Double_t                    fMaxCent;                    ///< max centrality for event selection
  Int_t                       fNcentBins;                  ///< how many centrality bins
  TString                     fCentEst;                    ///< name of V0 centrality estimator
  Bool_t                      fUseNewCentralityEstimation; ///< Use new centrality estimation (for 2015 data)
  Double_t                    fVertex[3];                  //!<!event vertex
  Int_t                       fNVertCont;                  //!<!event vertex number of contributors
  BeamType                    fBeamType;                   //!<!event beam type
  BeamType                    fForceBeamType;              ///< forced beam type
  Bool_t                      fNeedEmcalGeom;              ///< whether or not the task needs the emcal geometry
  AliEMCALGeometry           *fGeom;                       //!<!emcal geometry

  TObjArray                   fParticleCollArray;         ///< particle/track collection array
  TObjArray                   fClusterCollArray;          ///< cluster collection array
  std::vector <AliEmcalCorrectionCellContainer *> fCellCollArray; ///< Cells collection array
  
  TList * fOutput;                                        //!<! Output for histograms

  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionTask, 1); // EMCal correction task
  /// \endcond
};

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

  std::string GetBranchName() const { return fBranchName; }
  std::string GetName() const { return fName; }
  bool GetIsEmbedding() const { return fIsEmbedding; }
  AliVCaloCells * GetCells() const { return fCells; }

  void SetBranchName(std::string branchName) { fBranchName = branchName; }
  void SetName(std::string name ) { fName = name; }
  void SetIsEmbedding(bool isEmbedded) { fIsEmbedding = isEmbedded; }
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
