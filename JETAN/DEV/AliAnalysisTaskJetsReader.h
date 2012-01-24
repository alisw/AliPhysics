#ifndef ALIANALYSISTASKJETSREADER_H
#define ALIANALYSISTASKJETSREADER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//----------------------------------------------------------------
// Analysis task for interfacing the jet reader with the analysis framework
//
// Authors: magali.estienne@subatech.in2p3.fr 
//	    alexandre.shabetai@cern.ch
//----------------------------------------------------------------
 
#include "AliAnalysisTaskSE.h"

class AliJetReader;
class TTree;
class TString;
class AliJetCalTrkEvent;

class AliAnalysisTaskJetsReader : public AliAnalysisTaskSE
{
 public:
  AliAnalysisTaskJetsReader();
  AliAnalysisTaskJetsReader(const char* name);
  virtual               ~AliAnalysisTaskJetsReader();
  // Implementation of interface methods
  virtual void          UserCreateOutputObjects();
  virtual void          Init();
  virtual void          LocalInit() {Init();}
  virtual void          UserExec(Option_t *option);
  virtual void          SetConfigFile(const char *c) {fConfigFile = c;}
  virtual void          SetJetReader(AliJetReader *reader) {fJetReader = reader;}
  virtual void          Terminate(Option_t *option);
  virtual void          ReadAODFromOutput() {fReadAODFromOutput = kTRUE;}
  virtual AliJetReader* GetJetReader() {return fJetReader;}
    
 private:
  AliAnalysisTaskJetsReader(const AliAnalysisTaskJetsReader& rd);
  AliAnalysisTaskJetsReader& operator = (const AliAnalysisTaskJetsReader& rd);
  TString               fConfigFile;        //  The name of the ConfigFile
  AliJetReader*         fJetReader;         //  Pointer to the jet finder
  Bool_t                fReadAODFromOutput; //  Force reading of the AOD from the output 
  AliJetCalTrkEvent*    fReaderEvent;       //! Pointer to AliJetCalTrkEvent
  TTree*	        fExchangeTree;      //! Tree of AliJetCalTrkEvent
 
  ClassDef(AliAnalysisTaskJetsReader, 1)    // Jet reader Analysis task

};
#endif
