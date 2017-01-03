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

#include <AliAnalysisTaskSE.h>

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

  void      UserExec(Option_t *option)                           ;
  void      UserCreateOutputObjects()                            ;
  void      SetPtHardBin(Int_t r)                                 { fPtHardBin           = r; }
  void      SetAnchorRun(Int_t r)                                 { fAnchorRun           = r; }
  void      Terminate(Option_t *option)                          ;

  static const AliAnalysisTaskEmcalEmbeddingHelper* GetInstance() { return fgInstance       ; }

  AliVEvent* GetExternalEvent()                             const { return fExternalEvent   ; }

  TString GetTreeName()                                     const { return fTreeName; }
  /// Randomly start from an entry in the file. Will then loop around so that all entries are made available.
  Bool_t GetRandomEventNumberAccess()                       const { return fRandomEventNumberAccess; }
  /// Randomly select the first file. Continues sequentially afterwards
  Bool_t GetRandomFileAccess()                              const { return fRandomFileAccess; }
  TString GetFilePattern()                                  const { return fFilePattern; }
  Int_t GetStartingFileIndex()                              const { return fFilenameIndex; }
  TString GetFileListFilename()                             const { return fFileListFilename; }

  void SetESD(const char * treeName = "esdTree")                  { fTreeName     = treeName; }
  void SetAOD(const char * treeName = "aodTree")                  { fTreeName     = treeName; }
  void SetRandomEventNumberAccess(Bool_t b)                       { fRandomEventNumberAccess = b; }
  void SetRandomFileAccess(Bool_t b)                              { fRandomFileAccess = b; }
  void SetFilePattern(const char * pattern)                       { fFilePattern = pattern; }
  void SetStartingFileIndex(Int_t n)                              { fFilenameIndex = n; }
  void SetFileListFilename(const char * filename)                 { fFileListFilename = filename; }

  UInt_t GetTriggerMask()                                   const { return fTriggerMask; }
  Double_t GetZVertexCut()                                  const { return fZVertexCut; }
  Double_t GetMaxVertexDistance()                           const { return fMaxVertexDist; }

  void SetTriggerMask(UInt_t triggerMask)                         { fTriggerMask = triggerMask; }
  void SetZVertexCut(Double_t zVertex)                            { fZVertexCut = zVertex; }
  void SetMaxVertexDistance(Double_t distance)                    { fMaxVertexDist = distance; }

  static AliAnalysisTaskEmcalEmbeddingHelper * AddTaskEmcalEmbeddingHelper();

 protected:
  void            GetFilenames()        ;
  std::string     GenerateUniqueFileListFilename();
  void            SetupEmbedding()      ;
  Bool_t          SetupInputFiles()     ;
  Bool_t          GetNextEntry()        ;
  Bool_t          IsEventSelected()     ;
  Bool_t          InitEvent()           ;
  void            InitTree()            ;

  UInt_t                                        fTriggerMask;       ///<  Trigger selection mask
  Double_t                                      fZVertexCut;        ///<  Z vertex cut on embedded event
  Double_t                                      fMaxVertexDist;     ///<  Max distance between Z vertex of internal and embedded event

  bool                                          fInitializedNewFile; //!<! Notes where the entry indices have been initialized for a new tree in the chain
  bool                                          fInitializedEmbedding; //!<! Notes where the TChain has been initialized for embedding
  bool                                          fWrappedAroundTree; //!<! Notes whether we have wrapped around the tree, which is important if the offset into the tree is non-zero

  TString                                       fTreeName         ; ///<  Name of the ESD/AOD tree where the events are to be found
  Int_t                                         fAnchorRun        ; ///<  Anchor run for the given pythia production
  Int_t                                         fPtHardBin        ; ///<  ptHard bin for the given pythia production
  Bool_t                                        fRandomEventNumberAccess; ///<  If true, it will start embedding from a random entry in the file rather than from the first
  Bool_t                                        fRandomFileAccess ; ///< If true, it will start embedding from a random file in the input files list

  TString                                       fFilePattern      ; ///<  File pattern to select AliEn files
  TString                                       fFileListFilename ; ///<  Name of the file list containing paths to files to embed
  Int_t                                         fFilenameIndex    ; ///<  Index of vector containing paths to files to embed
  std::vector <std::string>                     fFilenames        ; //!<! Paths to the files to embed
  TFile                                        *fExternalFile     ; //!<! External file used for embedding
  TChain                                       *fChain            ; //!<! External TChain (tree) containing the events available for embedding
  Int_t                                         fCurrentEntry     ; //!<! Current entry in the current tree
  Int_t                                         fLowerEntry       ; //!<! First entry of the current tree to be used for embedding
  Int_t                                         fUpperEntry       ; //!<! Last entry of the current tree to be used for embedding
  Int_t                                         fOffset           ; //!<! Offset from fLowerEntry where the loop over the tree should start
  Int_t                                         fMaxNumberOfFiles ; //!<! Max number of files that are in the TChain
  Int_t                                         fFileNumber       ; //!<! File number corresponding to the current tree
  AliVEvent                                    *fExternalEvent    ; //!<! Current external event available for embedding

  static AliAnalysisTaskEmcalEmbeddingHelper   *fgInstance        ; //!<! Global instance of this class

 private:
  AliAnalysisTaskEmcalEmbeddingHelper(const AliAnalysisTaskEmcalEmbeddingHelper&)           ; // not implemented
  AliAnalysisTaskEmcalEmbeddingHelper &operator=(const AliAnalysisTaskEmcalEmbeddingHelper&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalEmbeddingHelper, 1);
  /// \endcond
};
#endif
