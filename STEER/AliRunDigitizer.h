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

const Int_t kMaxStreamsToMerge = 4;

// --- AliRoot header files ---

class AliDigitizer;
class AliMergeCombi;

class AliRunDigitizer: public TTask {

public:
  AliRunDigitizer();
  AliRunDigitizer(Int_t nInputStreams, Int_t sperb=1);
  virtual ~AliRunDigitizer();
  void      AddDigitizer(AliDigitizer *digitizer);
  void      SetOutputFile(TString fn);
  TString   GetOutputFile() const {return fOutputFileName;}
  void      SetOutputDir(TString dn) {fOutputDirName = dn;}
  TString   GetOutputDir() const {return fOutputDirName;}
  void      SetInputStream(Int_t stream, const char *inputName);
  void      SetFirstOutputEventNr(Int_t i) {fEvent = i;}
  void      SetNrOfEventsToWrite(Int_t i) {fNrOfEventsToWrite = i;}
  void      SetCopyTreesFromInput(Int_t i) {fCopyTreesFromInput = i;}
  Int_t     GetCopyTreesFromInput() {return fCopyTreesFromInput;}
  Int_t     GetOutputEventNr() const {return fEvent;}
  void      SetCombinationFileName(TString fn) {fCombinationFileName = fn;} 
  TString   GetCombinationFileName() const {return fCombinationFileName;}
  Int_t     GetNinputs() const {return fNinputs;}
  Int_t     GetMask(Int_t i) const {return fkMASK[i];}
  TTree*    GetInputTreeS(Int_t i) const {return fArrayTreeS[i];}
  TTree*    GetInputTreeH(Int_t i) const {return fArrayTreeH[i];}
  void      SetInputTreeTPCSBaseName(char * name) {
    fTreeTPCSBaseName = name;}
  TTree*    GetInputTreeTPCS(Int_t i) const {return fArrayTreeTPCS[i];}
  TTree*    GetInputTreeTRDS(Int_t i) const {return fArrayTreeTRDS[i];}
  TTree*    GetTreeD() const {return fTreeD;}
  void      SetTreeDTPCBaseName(char * name) {
    fTreeDTPCBaseName = name;}
  TTree*    GetTreeDTPC() const {return fTreeDTPC;} 
  TTree*    GetTreeDTRD() const {return fTreeDTRD;} 
  TTree*    GetTreeR() const {return fTreeR;} 
  void      Digitize(Option_t* option = 0);
  void      Exec(Option_t *option) {this->Digitize(option);}
  void      ExecuteTask(Option_t* option = 0);

  
// Nr of particles in all input files for a given event
//     (as numbered in the output file)
  Int_t GetNParticles(Int_t event) const;

// Nr of particles in input file input for a given event
//     (as numbered in this input file)
  Int_t GetNParticles(Int_t event, Int_t input) const;

// return pointer to an int array with input event numbers which were
// merged in the output event event
  Int_t* GetInputEventNumbers(Int_t event) const;

// return an event number of an eventInput from input file input
// which was merged to create output event event
  Int_t GetInputEventNumber(Int_t event, Int_t input) const;
  
// return pointer to particle with index i (index with mask)
  TParticle* GetParticle(Int_t i, Int_t event) const;

// return pointer to particle with index i in the input file input
// (index without mask)
  TParticle* GetParticle(Int_t i, Int_t input, Int_t event) const;

// return TString with input file name  
  TString GetInputFileName(const Int_t input, const Int_t order) const;
  
  Int_t     GetDebug() const {return fDebug;}
  void      SetDebug(Int_t level) {fDebug = level;}
  
private:
  Int_t             fkMASK[kMaxStreamsToMerge];  //! masks for track ids from
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
  TTree *           fTreeR;               //! output TreeR for ITS fast points
  Int_t             fNinputs;             // nr of input streams - can be taken from the TClonesArray dimension
  Int_t             fNinputsGiven;        // nr of input streams given by user
  TClonesArray *    fInputStreams;        // input streams
  TFile *           fInputFiles[kMaxStreamsToMerge];   //! p. to current input files
  TTree *           fArrayTreeS[kMaxStreamsToMerge];   //! array with p. to TreeS
  TTree *           fArrayTreeTPCS[kMaxStreamsToMerge];   //! array with p. to TreeD_75x40_100x60_x (TPC Sdigits)
  TTree *           fArrayTreeTRDS[kMaxStreamsToMerge];   //! array with p. to TreeSx_TRD (TRD Sdigits)
  TTree *           fArrayTreeH[kMaxStreamsToMerge];   //! array with p. to TreeH
  char *            fTreeDTPCBaseName;    //! basename of output TreeD for TPC
  char *            fTreeTPCSBaseName;    //! basename of output TreeS for TPC
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
  
  ClassDef(AliRunDigitizer,4)
};

#endif // ALIRUNDIGITIZER_H
