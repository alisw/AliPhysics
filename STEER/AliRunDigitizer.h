#ifndef ALIRUNDIGITIZER_H
#define ALIRUNDIGITIZER_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//
//  Manager Class for Merging/Digitization   
//                  
//  Author: Jiri Chudoba (CERN)
//
////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TNamed.h"
#include "TObjString.h"
#include "TArrayI.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TParticle.h"
#include "TTask.h"

#define MAXDETECTORS 20
#define MAXSTREAMSTOMERGE  4

// --- AliRoot header files ---

class AliDigitizer;
class AliMergeCombi;

class AliRunDigitizer: public TTask {

public:
  AliRunDigitizer(Int_t nInputStream=1, Int_t sperb=1);
  virtual ~AliRunDigitizer();
  void      AddDigitizer(AliDigitizer *digitizer);
  void      SetOutputFile(TString fn);
  TString   GetOutputFile() {return fOutputFileName;}
  void      SetOutputDir(TString dn) {fOutputDirName = dn;}
  TString   GetOutputDir() {return fOutputDirName;}
  void      SetInputStream(Int_t stream, char *inputName);
  void      SetFirstOutputEventNr(Int_t i) {fEvent = i;}
  void      SetNrOfEventsToWrite(Int_t i) {fNrOfEventsToWrite = i;}
  void      SetCopyTreesFromInput(Int_t i) {fCopyTreesFromInput = i;}
  Int_t     GetCopyTreesFromInput() {return fCopyTreesFromInput;}
  Int_t     GetOutputEventNr() {return fEvent;}
  void      SetCombinationFileName(TString fn) {fCombinationFileName = fn;} 
  TString   GetCombinationFileName() {return fCombinationFileName;}
  Int_t     GetNinputs() const {return fNinputs;}
  Int_t     GetMask(Int_t i) const {return fkMASK[i];}
  TTree*    GetInputTreeS(Int_t i) const {return fArrayTreeS[i];}
  TTree*    GetInputTreeH(Int_t i) const {return fArrayTreeH[i];}
  TTree*    GetInputTreeTPCS(Int_t i) const {return fArrayTreeTPCS[i];}
  TTree*    GetInputTreeTRDS(Int_t i) const {return fArrayTreeTRDS[i];}
  TTree*    GetTreeD() const {return fTreeD;}
  TTree*    GetTreeDTPC() const {return fTreeDTPC;} 
  TTree*    GetTreeDTRD() const {return fTreeDTRD;} 
  void      Digitize(Option_t* option = 0);
  void      Exec(Option_t *option) {this->Digitize();}
  void      ExecuteTask(Option_t* option = 0);

  
// Nr of particles in all input files for a given event
//     (as numbered in the output file)
  Int_t GetNParticles(Int_t event);

// Nr of particles in input file input for a given event
//     (as numbered in this input file)
  Int_t GetNParticles(Int_t event, Int_t input);

// return pointer to an int array with input event numbers which were
// merged in the output event event
  Int_t* GetInputEventNumbers(Int_t event);

// return an event number of an eventInput from input file input
// which was merged to create output event event
  Int_t GetInputEventNumber(Int_t event, Int_t input);
  
// return pointer to particle with index i (index with mask)
  TParticle* GetParticle(Int_t i, Int_t event);

// return pointer to particle with index i in the input file input
// (index without mask)
  TParticle* GetParticle(Int_t i, Int_t input, Int_t event);

  
  Int_t     GetDebug() const {return fDebug;}
  void      SetDebug(Int_t level) {fDebug = level;}
  
private:
  Int_t             fkMASK[MAXSTREAMSTOMERGE];  //! masks for track ids from
                                              //  different source files
  Int_t             fkMASKSTEP;           // step to increase MASK for
                                          // each input file
  TString           fOutputFileName;      // output file name
  TString           fOutputDirName;       // output dir name
  TFile *           fOutput;              //! pointer to the output file
  Int_t             fEvent;               // output event nr.
  Int_t             fNrOfEventsToWrite;   // Nr of events to write
  Int_t             fNrOfEventsWritten;   // Nr of events written
  Int_t             fCopyTreesFromInput;  // from which input file the trees
                                          // should be copied, -1 for no copies
  TTree *           fTreeD;               //! output TreeD
  TTree *           fTreeDTPC;            //! output TreeD for TPC
  TTree *           fTreeDTRD;            //! output TreeD for TRD
  Int_t             fNinputs;             // nr of input streams - can be taken from the TClonesArray dimension
  Int_t             fNinputsGiven;        // nr of input streams given by user
  TClonesArray *    fInputStreams;        // input streams
  TFile *           fInputFiles[MAXSTREAMSTOMERGE];   //! p. to current input files
  TTree *           fArrayTreeS[MAXSTREAMSTOMERGE];   //! array with p. to TreeS
  TTree *           fArrayTreeTPCS[MAXSTREAMSTOMERGE];   //! array with p. to TreeD_75x40_100x60_x (TPC Sdigits)
  TTree *           fArrayTreeTRDS[MAXSTREAMSTOMERGE];   //! array with p. to TreeSx_TRD (TRD Sdigits)
  TTree *           fArrayTreeH[MAXSTREAMSTOMERGE];   //! array with p. to TreeH
  AliMergeCombi *   fCombi;               // pointer to the combination object
  TArrayI           fCombination;         //! combination of events from
  TString           fCombinationFileName; // fn with combinations (used
                                          // with type 2 of comb.)
  Bool_t            ConnectInputTrees();
  Bool_t            InitGlobal();
  Bool_t            InitOutputGlobal();
  void              InitEvent();
  void              FinishEvent();
  void              FinishGlobal();
  Int_t             fDebug;                //! specifies debug level, 0 is min
  
  ClassDef(AliRunDigitizer,3)
};

#endif // ALIRUNDIGITIZER_H
