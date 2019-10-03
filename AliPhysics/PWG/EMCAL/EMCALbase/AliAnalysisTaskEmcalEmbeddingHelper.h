#ifndef ALIANALYSISTASKEMCALEMBEDDINGHELPER_H
#define ALIANALYSISTASKEMCALEMBEDDINGHELPER_H
/**
 * \file AliAnalysisTaskEmcalEmbeddingHelper.h
 * \brief Declaration of class AliAnalysisTaskEmcalEmbeddingHelper
 *
 * In this header file the class AliAnalysisTaskEmcalEmbeddingHelper is declared.
 * This class derives from AliAnalysisTaskSE and allows to open an external
 * ESD or AOD file, providing access to the events.
 *
 * \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * \author Raymond Ehlers <raymond.ehlers@cern.ch>, Yale University
 * \date Apr 28, 2016
 */

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class TString;
class TChain;
class TFile;
class AliVEvent;
class AliVHeader;
class AliGenPythiaEventHeader;
class AliEmcalList;

#include <iosfwd>
#include <vector>
#include <string>

#include <TStopwatch.h>
#include <TRandom3.h>
#include <AliAnalysisTaskSE.h>
#include "AliEventCuts.h"
#include "AliYAMLConfiguration.h"
#include "THistManager.h"

/**
 * \class AliAnalysisTaskEmcalEmbeddingHelper
 * \brief Implementation of task to embed external events.
 *
 * This class derives from AliAnalysisTaskSE and allows the user to open an external
 * ESD or AOD file, providing access to the events.
 *
 * The capabilities of the task are as follows:
 * - Open an ESD/AOD file according to certain customizable options, such as
 *   production tag, pt hard bin, run number, file pattern, etc.
 * - It selects an event according to customizable criteria, such as vertex, centrality,
 *   high pt track, pt hard bin, etc.
 * - Load the event in memory: this is the "external" event, as opposed to
 *   the "internal" event provided by the analysis manager
 * - Provide a public method GetExternalEvent() that allows to retrieve a pointer to
 *   the external event.
 *
 * Note that only one instance of this class is allowed in each train (singleton class).
 *
 * For the user, most of these details are handled by AliEmcalContainer derived tasks.
 * To access the embedded input objects, the user simply needs to set
 * AliEmcalContainer::SetIsEmbedding(Bool_t). This design ensures that usage is nearly
 * seamless. For instance, it will "just work" in an analysis task, and it very nearly
 * "just works" in the (new) EMCal correction framework (it only requires one small change
 * in the normal correction procedure beyond noting that the input objects are embedded).
 *
 * For further information on usage (including with the EMCal corrections), see the
 * [EMCal Embedding Documentation](\ref READMEemcEmbedding).
 *
 * \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * \author Raymond Ehlers <raymond.ehlers@cern.ch>, Yale University
 * \date Apr 28, 2016
 */
class AliAnalysisTaskEmcalEmbeddingHelper : public AliAnalysisTaskSE {
 public:

  AliAnalysisTaskEmcalEmbeddingHelper()                          ;
  AliAnalysisTaskEmcalEmbeddingHelper(const char *name)          ;
  virtual ~AliAnalysisTaskEmcalEmbeddingHelper()                 ;

  /**
   * @{
   * @name Inherited from AliAnalysisTaskSE
   */
  void      UserExec(Option_t *option)                           ;
  void      UserCreateOutputObjects()                            ;
  void      Terminate(Option_t *option)                          ;
  /* @} */

  static const AliAnalysisTaskEmcalEmbeddingHelper* GetInstance() { return fgInstance       ; }

  /**
   * @brief Retrieve the embedded event from the embedding helper
   * @return The embedded event
   */
  AliVEvent* GetExternalEvent()                             const { return fExternalEvent   ; }

  /**
   * @{
   * @name Properties of the embedding helper
   */
  /**
   * Initialize the Embedding Helper task. *Must* be called setting configuration for the task,
   * either during the run macro or wagon configuration. Once called, most of the configuration is locked in,
   * so be certain to change any configuration options before calling it.
   *
   * @param[in] removeDummyTask If true, the dummy task created with the configure wagon is removed.
   */
  bool Initialize(bool removeDummyTask = false);

  // Get
  Int_t GetPtHardBin()                                      const { return fPtHardBin; }
  Int_t GetNPtHardBins()                                    const { return fNPtHardBins; }
  TString GetTreeName()                                     const { return fTreeName; }
  Bool_t GetRandomEventNumberAccess()                       const { return fRandomEventNumberAccess; }
  Bool_t GetRandomFileAccess()                              const { return fRandomFileAccess; }
  TString GetFilePattern()                                  const { return fFilePattern; }
  TString GetInputFilename()                                const { return fInputFilename; }
  Int_t GetStartingFileIndex()                              const { return fFilenameIndex; }
  TString GetFileListFilename()                             const { return fFileListFilename; }
  bool GetCreateHistos()                                    const { return fCreateHisto; }

  // Set
  /// Set the pt hard bin which will be added into the file pattern. Can also be omitted and set directly in the pattern.
  void SetPtHardBin(Int_t n)                                      { fPtHardBin           = n; }
  /// Set the number of pt hard bins in the production to properly format the histograms
  void SetNPtHardBins(Int_t n)                                    { fNPtHardBins         = n; }
  /// Set to embed from ESD
  void SetESD(const char * treeName = "esdTree")                  { fTreeName     = treeName; }
  /// Set to embed from AOD
  void SetAOD(const char * treeName = "aodTree")                  { fTreeName     = treeName; }
  /// Set whether to print and plot execution time of InitTree()
  void SetPrintTimingInfoToLog(bool b)                            { fPrintTimingInfoToLog = b;}
  /**
   * Enable to begin embedding at a random entry in each embedded file. Will then loop around in order
   * so that all entries are made available.
   */
  void SetRandomEventNumberAccess(Bool_t b)                       { fRandomEventNumberAccess = b; }
  /// Randomly select the first file to embed from the file list. Continues sequentially afterwards
  void SetRandomFileAccess(Bool_t b)                              { fRandomFileAccess = b; }
  /// Sets the file pattern to select AliEn files. This pattern is used as input to the alien_find command.
  void SetFilePattern(const char * pattern)                       { fFilePattern = pattern; }
  /**
   * Sets the input filename used to select and open files. Note that this is just the filename, not the path!
   * This filename is also used as input to the alien_find command.
   */
  void SetInputFilename(const char * filename)                    { fInputFilename = filename; }
  /// Select the file ID to start embedding from.
  void SetStartingFileIndex(Int_t n)                              { fFilenameIndex = n; }
  /// Set the path to a file containing the list of files to embed
  void SetFileListFilename(const char * filename)                 { fFileListFilename = filename; }
  /// Create QA histograms. These are necessary for proper scaling, so be careful disabling them!
  void SetCreateHistos(bool b)                                    { fCreateHisto = b; }
  /// Set path to %YAML configuration file
  void SetConfigurationPath(const char * path)                    { fConfigurationPath = path; }
  /* @} */

  /**
   * @{
   * @name Internal event selection
   */
  /// Whether internal event selection is enabled.
  bool GetUseInternalEventSelection()                       const { return fUseInternalEventSelection; }
  /// Enable internal event selection. Can also be enabled through the %YAML configuration.
  void SetUseInternalEventSelection(bool b = true)                { fUseInternalEventSelection = b; }
  /**
   * If true, it indicates that the embedded event was used by the embedding helper and that it is
   * available for use. If false, this internal event should be ignored.
   *
   * This value can return false if internal event selection is enabled and the internal event was
   * rejected.
   *
   * @return true if the embedded event was actually embedded and should be used.
   */
  bool EmbeddedEventUsed()                                  const { return fEmbeddedEventUsed; }
  /// Whether to use manual cuts for AliEventCuts.
  bool GetUseManualInternalEventSelection()                 const { return fUseManualInternalEventCuts; }
  /// Enable manual internal event cuts. Can also be enabled through the %YAML configuration.
  /// Must be configured through the retrieving the AliEventCuts object and configuring the cuts.
  /// It is not included via %YAML because it is rather difficult to fully map, especially for an infrequently used mode.
  void SetUseManualInternalEventCuts(bool b = true)               { fUseManualInternalEventCuts = b; }
  /// Event cuts object for accessing centrality, etc from another task if so inclined.
  const AliEventCuts * GetInternalEventCuts()               const { return (fUseInternalEventSelection ? &fInternalEventCuts : nullptr); }
  /// Event cuts object for configuring for setting manual cuts
  /// To access this object, manual event cuts must be enabled!
  AliEventCuts * GetInternalEventCuts();
  /// Set internal event centrality selection
  void SetCentralityRange(double min, double max)                 { fCentMin = min; fCentMax = max; }
  void SetRandomRejectionFactor(Double_t configure = 1.)          { fRandomRejectionFactor = configure; }
  /* @} */

  /**
   * @{
   * @name Options for the embedded event
   */
  UInt_t GetTriggerMask()                                   const { return fTriggerMask; }
  bool GetMCRejectOutliers()                                const { return fMCRejectOutliers; }
  Double_t GetPtHardJetPtRejectionFactor()                  const { return fPtHardJetPtRejectionFactor; }
  Double_t GetZVertexCut()                                  const { return fZVertexCut; }
  Double_t GetMaxVertexDistance()                           const { return fMaxVertexDist; }

  void SetTriggerMask(UInt_t triggerMask)                         { fTriggerMask = triggerMask; }
  void SetMCRejectOutliers(bool reject = true)                    { fMCRejectOutliers = reject; }
  void SetPtHardJetPtRejectionFactor(double factor)               { fPtHardJetPtRejectionFactor = factor; }
  void SetZVertexCut(Double_t zVertex)                            { fZVertexCut = zVertex; }
  void SetMaxVertexDistance(Double_t distance)                    { fMaxVertexDist = distance; }
  /* @} */

  /**
   * @{
   * @name Properties of the embedded event
   */
  AliVHeader * GetEventHeader()                             const { return fExternalHeader; }
  AliGenPythiaEventHeader * GetPythiaHeader()               const { return fPythiaHeader; }
  double GetPythiaXSection()                                const { return fPythiaCrossSection; }
  int GetPythiaTrials()                                     const { return fPythiaTrials; }
  double GetPythiaPtHard()                                  const { return fPythiaPtHard; }
  /* @} */

  /**
   * @{
   * @name pT hard bin auto configuration
   * @brief %Setup pt hard bin auto configuration to be used on the LEGO train. See AutoConfigurePtHardBins() and the variable definitions for the purpose of each variable.
   */
  bool GetAutoConfigurePtHardBins()                         const { return fAutoConfigurePtHardBins; }
  std::string GetAutoConfigureBasePath()                    const { return fAutoConfigureBasePath; }
  std::string GetAutoConfigureTrainTypePath()               const { return fAutoConfigureTrainTypePath; }
  std::string GetAutoConfigureIdentifier()                  const { return fAutoConfigureIdentifier; }

  void SetAutoConfigurePtHardBins(bool configure = true)          { fAutoConfigurePtHardBins = configure; }
  void SetAutoConfigureBasePath(std::string path)                 { fAutoConfigureBasePath = path; }
  void SetAutoConfigureTrainTypePath(std::string path)            { fAutoConfigureTrainTypePath = path; }
  void SetAutoConfigureIdentifier(std::string path)               { fAutoConfigureIdentifier = path; }
  /* @} */

  /**
   * @{
   * @name Utility functions
   */
  /**
   * Add task function. This contains the normal AddTask functionality, except in compiled code, making errors
   * easier to spot than in CINT. The AddTask macro still exists for use on the LEGO train, but simply wraps this
   * function.
   *
   * @return An properly instance of AliAnalysisTaskEmcalEmbeddingHelper, added to the current analysis manager.
   */
  static AliAnalysisTaskEmcalEmbeddingHelper * AddTaskEmcalEmbeddingHelper();
  /**
   * Retrieve an existing embedding helper to perform further configuration. This should
   * _ONLY_ be used on the LEGO train.
   *
   * To achieve this, a dummy task is created when the configure task is called because AliAnalysisTaskCfg
   * requires that all wagons add a task. Then, when Initialize(true) is called on the embedding helper task, the
   * dummy task is removed. This is a hack, but is required to work around constraints in AliAnalysisTaskCfg.
   *
   * @return An existing (usually unconfigured) EMCal Embedding Helper.
   */
  static AliAnalysisTaskEmcalEmbeddingHelper* ConfigureEmcalEmbeddingHelperOnLEGOTrain();

  // Printing
  friend std::ostream & operator<<(std::ostream &in, const AliAnalysisTaskEmcalEmbeddingHelper &myTask);
  void Print(Option_t* opt = "") const;
  std::ostream & Print(std::ostream &in) const;
  std::string toString(bool includeFileList = false) const;
  /* @} */

  /**
   * @{
   * @name To be ignored
   */
  /**
   * @brief **SHOULD NOT BE USED! Use GetExternalEvent()!**
   * **SHOULD NOT BE USED! Use GetExternalEvent()!**
   * Returns the external event by overloading InputEvent() defined in AliAnalysisTaskSE.
   * This is used by AliEmcalCorrectionEventManager, but it should not be used by the user!
   * Instead, use GetExternalEvent().
   *
   * Note that it cannot be protected, because we need to call it externally.
   *
   * @return The external event
   */
  AliVEvent* InputEvent()                                   const { return GetExternalEvent(); }
  /* @} */

 protected:
  virtual void    RetrieveTaskPropertiesFromYAMLConfig();
  bool            GetFilenames()        ;
  void            DeterminePythiaXSecFilename();
  bool            IsRunInRunlist(const std::string & path) const;
  bool            InitializeYamlConfig();
  bool            AutoConfigurePtHardBins();
  std::string     GenerateUniqueFileListFilename() const;
  std::string     RemoveTrailingSlashes(std::string filename) const;
  void            DetermineFirstFileToEmbed();
  void            SetupEmbedding()      ;
  Bool_t          SetupInputFiles()     ;
  std::string     ConstructFullPythiaXSecFilename(std::string inputFilename, const std::string & pythiaFilename, bool testIfExists) const;
  Bool_t          GetNextEntry()        ;
  void            SetEmbeddedEventProperties();
  void            RecordEmbeddedEventProperties();
  Bool_t          IsEventSelected()     ;
  virtual Bool_t  CheckIsEmbeddedEventSelected();
  Bool_t          InitEvent()           ;
  void            InitTree()            ;
  bool            PythiaInfoFromCrossSectionFile(std::string filename);
  // Validation helper
  void            ValidatePhysicsSelectionForInternalEventSelection();
  // Helper functions
  bool            IsFileAccessible() const;
  void            ConnectToAliEn() const;
  // LEGO Train utility
  void            RemoveDummyTask() const;

  UInt_t                                        fTriggerMask;       ///<  Trigger selection mask
  bool                                          fMCRejectOutliers;  ///<  If true, MC outliers will be rejected
  Double_t                                      fPtHardJetPtRejectionFactor; ///<  Factor which the pt hard bin is multiplied by to compare against pythia header jets pt
  Double_t                                      fZVertexCut;        ///<  Z vertex cut on embedded event
  Double_t                                      fMaxVertexDist;     ///<  Max distance between Z vertex of internal and embedded event

  bool                                          fInitializedConfiguration; ///< Notes if the configuration has been initialized
  bool                                          fInitializedNewFile; //!<! Notes where the entry indices have been initialized for a new tree in the chain
  bool                                          fInitializedEmbedding; //!<! Notes where the TChain has been initialized for embedding
  bool                                          fWrappedAroundTree; //!<! Notes whether we have wrapped around the tree, which is important if the offset into the tree is non-zero

  TString                                       fTreeName         ; ///<  Name of the ESD/AOD tree where the events are to be found
  Int_t                                         fNPtHardBins      ; ///<  Total number of pt hard bins
  Int_t                                         fPtHardBin        ; ///<  ptHard bin for the given pythia production
  Bool_t                                        fRandomEventNumberAccess; ///<  If true, it will start embedding from a random entry in the file rather than from the first
  Bool_t                                        fRandomFileAccess ; ///<  If true, it will start embedding from a random file in the input files list
  bool                                          fCreateHisto      ; ///<  If true, create QA histograms
  PWG::Tools::AliYAMLConfiguration              fYAMLConfig       ; ///<  Hanldes configuration from YAML

  bool                                  fUseInternalEventSelection; ///<  If true, apply internal event selection though AliEventCuts
  bool                                 fUseManualInternalEventCuts; ///<  If true, manual event cuts mode will be used for AliEventCuts
  AliEventCuts                                  fInternalEventCuts; ///<  If enabled, Handles internal event selection
  bool                                          fEmbeddedEventUsed; //!<! If true, the internal event was selected, so the embedded event is used. Defaults to true so other tasks are not disrupted if internal event selection is disabled.
  bool                                  fValidatedPhysicsSelection; ///<  Validate that the physics selection is set appropriately.
  UInt_t                                 fInternalEventTriggerMask; ///<  Internal event physics selection (trigger mask) to be used with AliEventCuts.
  double                                        fCentMin          ; ///<  Minimum centrality for internal event selection
  double                                        fCentMax          ; ///<  Maximum centrality for internal event selection
  Double_t                                      fRandomRejectionFactor; ///< factor by which to reject events
  TRandom3                                      fRandom           ; ///< for random rejection of events

  bool                                    fAutoConfigurePtHardBins; ///<  If true, attempt to auto configure pt hard bins. Only works on the LEGO train.
  std::string                               fAutoConfigureBasePath; ///<  The base path to the auto configuration (for example, "/alice/cern.ch/user/a/alitrain/")
  std::string                          fAutoConfigureTrainTypePath; ///<  The path associated with the train type (for example, "PWGJE/Jets_EMC_PbPb/")
  std::string                             fAutoConfigureIdentifier; ///<  How the auto configuration %YAML file should be identified. (for example, "rehlersTrain")

  TString                                       fFilePattern      ; ///<  File pattern to select AliEn files using alien_find
  TString                                       fInputFilename    ; ///<  Filename of input root files
  TString                                       fFileListFilename ; ///<  Name of the file list containing paths to files to embed
  Int_t                                         fFilenameIndex    ; ///<  Index of vector containing paths to files to embed
  std::vector <std::string>                     fFilenames        ; ///<  Paths to the files to embed
  std::string                                   fConfigurationPath; ///<  Path to %YAML configuration
  std::vector <std::string>                     fEmbeddedRunlist  ; ///<  Good runlist for files to embed
  std::string                                  fPythiaXSecFilename; ///<  Name of the pythia x sec filename (either "pyxsec.root" or "pyxsec_hists.root")
  std::vector <std::string>                     fPythiaCrossSectionFilenames; ///< Paths to the pythia xsection files
  TFile                                        *fExternalFile     ; //!<! External file used for embedding
  TChain                                       *fChain            ; //!<! External TChain (tree) containing the events available for embedding
  Int_t                                         fCurrentEntry     ; //!<! Current entry in the current tree
  Int_t                                         fLowerEntry       ; //!<! First entry of the current tree to be used for embedding
  Int_t                                         fUpperEntry       ; //!<! Last entry of the current tree to be used for embedding
  Int_t                                         fOffset           ; //!<! Offset from fLowerEntry where the loop over the tree should start
  UInt_t                                        fMaxNumberOfFiles ; //!<! Max number of files that are in the TChain
  UInt_t                                        fFileNumber       ; //!<! File number corresponding to the current tree
  THistManager                                  fHistManager      ; ///< Manages access to all histograms
  AliEmcalList                                 *fOutput           ; //!<! List which owns the output histograms to be saved
  AliVEvent                                    *fExternalEvent    ; //!<! Current external event available for embedding
  AliVHeader                                   *fExternalHeader   ; //!<! Header of the current external event
  AliGenPythiaEventHeader                      *fPythiaHeader     ; //!<! Pythia header of the current external event

  int                                           fPythiaTrials     ; //!<! Number of pythia trials for the current event (extracted from the pythia header).
  int                                           fPythiaTrialsFromFile; //!<! Average number of trials extracted from a xsec file.
  double                                        fPythiaCrossSection; //!<! Pythia cross section for the current event (extracted from the pythia header).
  double                                        fPythiaCrossSectionFromFile; //!<! Average pythia cross section extracted from a xsec file.
  double                                        fPythiaPtHard     ; //!<! Pt hard of the current event (extracted from the pythia header).
  
  bool                                          fPrintTimingInfoToLog; ///< Flag to print time to execute InitTree(), for logging purposes
  TStopwatch                                    fTimer            ;    //!<! Timer for the InitTree() function

  static AliAnalysisTaskEmcalEmbeddingHelper   *fgInstance        ; //!<! Global instance of this class

 private:
  AliAnalysisTaskEmcalEmbeddingHelper(const AliAnalysisTaskEmcalEmbeddingHelper&)           ; // not implemented
  AliAnalysisTaskEmcalEmbeddingHelper &operator=(const AliAnalysisTaskEmcalEmbeddingHelper&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalEmbeddingHelper, 12);
  /// \endcond
};
#endif
