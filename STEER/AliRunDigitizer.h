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

#define MAXDETECTORS 20
#define MAXFILESTOMERGE  4

// --- AliRoot header files ---

class AliDigitizer;

class AliRunDigitizer: public TNamed {

public:
  AliRunDigitizer();
  virtual ~AliRunDigitizer();
  void      AddDigitizer(AliDigitizer *digitizer);
  void      SetOutputFile(TString fn) {fOutputFileName = fn;}
  TString   GetOutputFile() {return fOutputFileName;}
  void      SetOutputDir(TString dn) {fOutputDirName = dn;}
  TString   GetOutputDir() {return fOutputDirName;}
  void      SetFirstOutputEventNr(Int_t i) {fEvent = i;}
  void      SetNrOfEventsToWrite(Int_t i) {fNrOfEventsToWrite = i;}
  void      SetCopyTreesFromInput(Int_t i) {fCopyTreesFromInput = i;}
  Int_t     GetCopyTreesFromInput() {return fCopyTreesFromInput;}
  Int_t     GetOutputEventNr() {return fEvent;}
  Bool_t    SetInput(char* inputFileString);
  void      SetCombinationFileName(TString fn) {fCombinationFileName = fn;} 
  TString   GetCombinationFileName() {return fCombinationFileName;}
  void      SetCombinationType(Int_t i) {fCombinationType = i;}
  Int_t     GetCombinationType() const {return fCombinationType;}
  void      PrintCombination();   // debug method
  Int_t     GetNinputs() const {return fNinputs;}
  Int_t     GetMask(Int_t i) const {return fkMASK[i];}
  TTree*    GetInputTreeS(Int_t i) const {return fArrayTreeS[i];}
  TTree*    GetInputTreeH(Int_t i) const {return fArrayTreeH[i];}
  TTree*    GetInputTreeTPCS(Int_t i) const {return fArrayTreeTPCS[i];}
  Int_t     GetNextMask(Int_t i) const;
  TTree*    GetNextTreeH(TTree *) const;
  TTree*    GetNextTreeS(TTree *) const;
  TTree*    GetNextTreeTPCS(TTree *) const;
  TTree*    GetTreeD() const {return fTreeD;}
  void      Digitize();

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
  AliDigitizer *    fDigitizers[MAXDETECTORS];  //! pointers to registered
                                                //  digitizers
// the constant 20 corresponds to  MAXDETECTORS = 20 - could be done better
  Int_t             fNDigitizers;         //! nr. of registered digitizers
  Int_t             fkMASK[MAXFILESTOMERGE];  //! masks for track ids from
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
  Int_t             fNinputs;             // nr of input files
  TClonesArray *    fInputFileNames;      // input files names
  TClonesArray *    fInputFiles;          //! pointers to file handles
  TTree *           fArrayTreeS[MAXFILESTOMERGE];   //! array with p. to TreeS
  TTree *           fArrayTreeTPCS[MAXFILESTOMERGE];   //! array with p. to TreeD_75x40_100x60_x (TPC Sdigits)
  TTree *           fArrayTreeH[MAXFILESTOMERGE];   //! array with p. to TreeH
  Int_t             fMinNEvents;          // min nr. of events found in inputs
  TArrayI           fCombination;         //! combination of events from
  Int_t             fCombinationType;     // defines the algorithm how
                                          // to generated combination
  TString           fCombinationFileName; // fn with combinations (used
                                          // with type 2 of comb.)
  Bool_t            MakeCombination();    // create next combination
  void              ConnectInputTrees();
  Bool_t            InitGlobal();
  Bool_t            InitOutputGlobal();
  void              InitPerEvent();
  void              FinishPerEvent();
  void              FinishGlobal();
  Int_t             fDebug;                //! specifies debug level, 0 is min
  
  TFile* ConnectInputFile(Int_t input);    // open input file with index input

  ClassDef(AliRunDigitizer,1)
};

#endif // ALIRUNDIGITIZER_H
