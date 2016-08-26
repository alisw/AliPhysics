#ifndef ALIEMCALCORRECTIONTASK_H
#define ALIEMCALCORRECTIONTASK_H



#if !(defined(__CINT__) || defined(__MAKECINT__))
#include <yaml-cpp/yaml.h>
#endif

#include "AliAnalysisTaskSE.h"
#include <Rtypes.h>

// Base component
class AliEmcalCorrectionComponent;
// Geometry
class AliEMCALGeometry;

#include "AliParticleContainer.h"
#include "AliMCParticleContainer.h"
#include "AliTrackContainer.h"
#include "AliClusterContainer.h"

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
    kNA       = -1,//!< Undefined
    kpp       = 0, //!< Proton-Proton
    kAA       = 1, //!< Nucleus-Nucleus
    kpA       = 2  //!< Proton-Nucleus
  };

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
  void                        SetCaloCellsName(const char* name)                    { fCaloCellsName = name; }

  void                        SetForceBeamType(BeamType f)                          { fForceBeamType     = f                              ; }
  void                        SetRunPeriod(const char* runPeriod)                   { fRunPeriod = runPeriod; fRunPeriod.ToLower(); }
  const TString &             GetRunPeriod()                                  const { return fRunPeriod; }

  void                        SetCreateNewObjectBranches(bool flag)                 { fCreateNewObjectBranches = flag; }
  void                        SetUseNewCentralityEstimation(Bool_t b)               { fUseNewCentralityEstimation = b                     ; }

  const std::vector<AliEmcalCorrectionComponent *> & CorrectionComponents() { return fCorrectionComponents; }

  bool WriteConfigurationFile(std::string filename, bool userConfig = false);

 private:
  AliEmcalCorrectionTask(const AliEmcalCorrectionTask &);             // Not implemented
  AliEmcalCorrectionTask &operator=(const AliEmcalCorrectionTask &);  // Not implemented

  static inline bool doesFileExist(const std::string & filename);
  void SetupConfigurationFilePath(std::string & filename, bool userFile = false);

  void RetrieveExecutionOrder(std::vector <std::string> & componentsToAdd);
  void InitializeComponents();

  void                        CreateNewObjectBranches();
  void                        CopyBranchesToNewObjects();
  void                        CopyClusters(TClonesArray *orig, TClonesArray *dest);
  void                        CleanupCreatedBranches();

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
  bool                        fCreateNewObjectBranches;    ///< Create new branches for cells and clusters
  std::string                 fCreatedClusterBranchName;   ///< Name of created cluster branch
  std::string                 fCreatedTrackBranchName;     ///< Name of created track branch
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
  TString                     fCaloCellsName;              ///< name of calo cell collection
  Bool_t                      fNeedEmcalGeom;              ///< whether or not the task needs the emcal geometry
  AliEMCALGeometry           *fGeom;                       //!<!emcal geometry

  TObjArray                   fParticleCollArray;         ///< particle/track collection array
  TObjArray                   fClusterCollArray;          ///< cluster collection array
  AliVCaloCells               *fCaloCells;                //!<! pointer to calo cells
  AliVCaloCells               *fCaloCellsFromInputEvent;  //!<! pointer to calo cells from the input event
  
  TList * fOutput;                                        //!<! Output for histograms

  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionTask, 1); // EMCal correction task
  /// \endcond
};

#endif /* ALIEMCALCORRECTIONTASK_H */ 
