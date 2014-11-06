#ifndef ALIANALYSISTASKJETSFINDER_H
#define ALIANALYSISTASKJETSFINDER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//----------------------------------------------------------------
// Analysis task for interfacing the jet finders with the analysis framework
//
// Author: magali.estienne@subatech.in2p3.fr 
//	   alexandre.shabetai@cern.ch
//----------------------------------------------------------------
 
#include "AliAnalysisTaskSE.h"
#include "AliJetCalTrk.h"

class AliJetFinder;
class AliJetHistos;
class AliAODExtension;
class TTree;
class TString;

class AliAnalysisTaskJetsFinder : public AliAnalysisTaskSE
{
 public:
  AliAnalysisTaskJetsFinder();
  AliAnalysisTaskJetsFinder(const char* name);
  virtual             ~AliAnalysisTaskJetsFinder();
  // Implementation of interface methods
  virtual void        UserCreateOutputObjects();
  virtual void        Init();
  virtual Bool_t      Notify();
  virtual void        LocalInit() {Init();}
  virtual void        ConnectInputData(Option_t *);
  virtual void        UserExec(Option_t *option);
  virtual void        SetConfigFile(const char *c) {fConfigFile = c;}
  virtual void        SetJetFinder(AliJetFinder *finder) {fJetFinder = finder;}
  virtual void        SetNonStdBranch(const char *c){fNonStdBranch = c;}
  virtual const char* GetNonStdBranch(){return fNonStdBranch.Data();}
  virtual void        SetNonStdOutputFile(const char *c){fNonStdFile = c;}
  virtual const char* GetNonStdOutputFile() {return fNonStdFile.Data();}
  virtual void        SetBookAODBackground(Bool_t b){fUseAODBackground = b;}
  virtual void        Terminate(Option_t *option);
  virtual void        SetFilterPt(Float_t f){fFilterPt = f;}
    
  AliJetFinder*       GetJetFinder() {return fJetFinder;}

 private:
  AliAnalysisTaskJetsFinder(const AliAnalysisTaskJetsFinder &det);
  AliAnalysisTaskJetsFinder &operator=(const AliAnalysisTaskJetsFinder &det);
    
  TString             fConfigFile;        //  The name of the ConfigFile
  TString             fNonStdBranch;      //  The name of the non-std branch name
  TString             fNonStdFile;        //  The optional name of the output file the non-std brnach is written to
  AliJetFinder*       fJetFinder;         //  Pointer to the jet finder 
  AliJetHistos*       fHistos;            //! Histogram manager class
  AliAODExtension*    fAODExtension;      //  AOD extension in case we write a non-sdt brnach to a separate file and the aod is standard
  TList*              fListOfHistos;      //! Output list of histograms
  TTree*              fTreeI;             //! Input Tree 
  AliJetCalTrkEvent*  fEvent;             //! Pointer to jet input objects
  Bool_t              fUseAODBackground;  //  Decide wether we book the backround branch or not
  Float_t             fFilterPt;          //  Use this as a switch for writing the AOD, minium p_T of leading jet

  ClassDef(AliAnalysisTaskJetsFinder, 1)  // Jet Finder Analysis task

};
 
#endif
