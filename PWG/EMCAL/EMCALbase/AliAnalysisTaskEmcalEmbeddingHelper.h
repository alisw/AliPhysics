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
 * \date Apr 28, 2016
 */

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TString.h>
#include <AliAnalysisTaskSE.h>

class AliVEvent;
class TTree;
class TFile;

/**
 * \class AliAnalysisTaskEmcalEmbeddingHelper
 * \brief Implementation of task to embed external events.
 *
 * This class derives from AliAnalysisTaskSE and allows to open an external
 * ESD or AOD file, providing access to the events.
 * 1) Open an ESD/AOD file according to certain customizable options, such as
 *    production tag, pt hard bin, run number, file pattern, etc.
 * 2) Select an event according to customizable criteria, such as vertex, centrality,
 *    high pt track, pt hard bin, etc.
 * 3) Load the event in memory: this is the "external" event, as opposed to
 *    the "internal" event provided by the analysis manager
 * 4) Provide a public static method GetExternalEvent() that allows
 *    to retrieve a pointer to the external event
 * 5) Only one instance of this class is allowed in each train (singleton class).
 */
class AliAnalysisTaskEmcalEmbeddingHelper : public AliAnalysisTaskSE {
 public:

  AliAnalysisTaskEmcalEmbeddingHelper()                          ;
  AliAnalysisTaskEmcalEmbeddingHelper(const char *name)          ;
  virtual ~AliAnalysisTaskEmcalEmbeddingHelper()                 ;

  void      UserExec(Option_t *option)                           ;
  void      UserCreateOutputObjects()                            ;
  void           SetPYTHIAPath(const char* p)                      { fPYTHIAPath                                = p ; }
  void           SetPtHardBin(Int_t r)                               {fPtHardBin                                   =r;}
  void           SetAnchorRun(Int_t r)                             { fAnchorRun                                 = r ; }
  void      Terminate(Option_t *option)                          ;

  static const AliAnalysisTaskEmcalEmbeddingHelper* GetInstance() { return fgInstance        ; }

  AliVEvent* GetExternalEvent()                             const { return fExternalEvent    ; }

  void SetESD()                  { fTreeName     = "esdTree"; }
  void SetAOD()                  { fTreeName     = "aodTree"; }
  void SetAttempts(Int_t n)      { fAttempts     = n        ; }
  void SetRandomAccess(Bool_t b) { fRandomAccess = b        ; }

 protected:
  TFile          *GetNextFile()         ;
  TString         GetNextFileName()     ;
  Bool_t          OpenNextFile()        ;
  Bool_t          GetNextEntry()        ;
  Bool_t          IsEventSelected()     ;
  Bool_t          InitEvent()           ;

  TString                                       fTreeName         ; ///<  Name of the ESD/AOD tree where the events are to be found
  TString                                       fPYTHIAPath        ; ///< Path to the given Pythia AOD file 
  Int_t                                           fAnchorRun         ; ///<anchor run for the given pythia production
  Int_t                                           fPtHardBin ;///<pthard bin for the given pythia production
  Int_t                                           fAttempts         ; ///<  Number of attempts to open the external file before giving up
  Bool_t                                        fRandomAccess     ; ///<  If true, it will start embedding from a random entry in the file rather than from the first

  TFile                                        *fExternalFile     ; //!<! External file used for embedding
  TTree                                        *fExternalTree     ; //!<! External tree containing the events available for embedding
  Int_t                                         fCurrentEntry     ; //!<! Current entry
  Int_t                                         fFirstEntry       ; //!<! Entry of the current file from which embedding started
  Int_t                                         fLastEntry        ; //!<! Last entry to be used for embedding from the current file
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
