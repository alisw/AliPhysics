#ifndef ALIANALYSISTASKJETS_H
#define ALIANALYSISTASKJETS_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
#include "AliAnalysisTaskSE.h"


class AliJetFinder;
class AliESDEvent;
class AliAODEvent;
class AliJetHistos;
class AliAODExtension;
class TTree;
class TChain;
class TString;



class AliAnalysisTaskJets : public AliAnalysisTaskSE
{
 public:
    AliAnalysisTaskJets();
    AliAnalysisTaskJets(const char* name);
    AliAnalysisTaskJets(const char* name, TChain* chain);
    virtual ~AliAnalysisTaskJets() {;}
    // Implementation of interface methods
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void UserExec(Option_t *option);
    virtual void SetConfigFile(const char *c) {fConfigFile = c;}
    virtual void SetJetFinder(AliJetFinder *finder) {fJetFinder = finder;}
    virtual void SetNonStdBranch(const char *c){fNonStdBranch = c;}
    virtual void SetNonStdOutputFile(const char *c){fNonStdFile = c;}
    virtual void Terminate(Option_t *option);
    virtual void ReadAODFromOutput() {fReadAODFromOutput = kTRUE;}

    
 private:
  AliAnalysisTaskJets(const AliAnalysisTaskJets &det);
  AliAnalysisTaskJets &operator=(const AliAnalysisTaskJets &det);
    
 private:
  TString       fConfigFile;        // the name of the ConfigFile
  TString       fNonStdBranch;      // the name of the non-std branch name
  TString       fNonStdFile;        // The optional name of the output file the non-std brnach is written to
  AliJetFinder* fJetFinder;         //  Pointer to the jet finder 
  AliJetHistos* fHistos;            //  Histogram manager class
  AliAODExtension* fAODExtension;   //  AOD extension we in case we write a non-sdt brnach to a separate file and the aod is standard
  TList*        fListOfHistos;      //  Output list of histograms
  TChain*       fChain;             //  Chain 
  Int_t         fOpt;               //  Detector configuration used
  Bool_t        fReadAODFromOutput; //  Force reading of the AOD from the output
  
  ClassDef(AliAnalysisTaskJets, 4); // Analysis task for standard jet analysis
};
 
#endif
