#ifndef ALIRUNDIGITIZER_H
#define ALIRUNDIGITIZER_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//
//  Manager Class for Merging/Digitization   
//  This handles Merging and Digitisation of AliRoot events                
//  Author: Jiri Chudoba (CERN)
//
////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

#include "TArrayI.h"
#include "TTask.h"
#include "TClonesArray.h"
class TFile;
class TParticle;
class TTree;

// --- AliRoot header files ---

#include "AliStream.h" 
class AliDigitizer;
class AliMergeCombi;
class AliRunLoader;

#define MAXSTREAMSTOMERGE 4

class AliRunDigitizer: public TTask {

public:
  AliRunDigitizer();
  AliRunDigitizer(Int_t nInputStreams, Int_t sperb=1);

  virtual ~AliRunDigitizer();

  void      ExecuteTask(Option_t* option = 0);
  void      Exec(Option_t *option) {this->Digitize(option);}
  void      Digitize(Option_t* option = 0);
  void      AddDigitizer(AliDigitizer *digitizer);

  void      SetOutputFile(TString fn);
  TString   GetOutputFile() const {return fOutputFileName;}
  
  void      SetOutputDir(TString dn) {fOutputDirName = dn;}
  TString   GetOutputDir() const {return fOutputDirName;}
  
  void      SetInputStream(Int_t stream, const char *inputName, TString foldername = "");
  
  void      SetFirstOutputEventNr(Int_t i) {fEvent = i;}
  void      SetNrOfEventsToWrite(Int_t i) {fNrOfEventsToWrite = i;}
  void      SetCopyTreesFromInput(Int_t i) {fCopyTreesFromInput = i;}
  Int_t     GetCopyTreesFromInput() const {return fCopyTreesFromInput;}
  Int_t     GetOutputEventNr() const {return fEvent;}
  void      SetCombinationFileName(TString fn) {fCombinationFileName = fn;} 
  TString   GetCombinationFileName() const {return fCombinationFileName;}
  Int_t     GetMask(Int_t i) const {return fkMASK[i];}


  Int_t     GetNinputs() const {return fNinputs;}
  const TString& GetInputFolderName(Int_t i) const;
  const char* GetOutputFolderName();


    
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
  
  AliStream * GetInputStream(Int_t index) const { return dynamic_cast<AliStream *>(fInputStreams->At(index)) ; }
// return pointer to particle with index i (index with mask)
  TParticle* GetParticle(Int_t i, Int_t event) const;

// return pointer to particle with index i in the input file input
// (index without mask)
  TParticle* GetParticle(Int_t i, Int_t input, Int_t event) const;

// return TString with input file name  
  TString GetInputFileName(Int_t input, Int_t order) const;
  
private:
  AliRunDigitizer(const AliRunDigitizer& dig); // not implemented
  AliRunDigitizer& operator=(const AliRunDigitizer& dig); // not implemented
  void Copy(TObject& dig) const;
  Bool_t            ConnectInputTrees();
  Bool_t            InitGlobal();
  Bool_t            InitOutputGlobal();
  void              InitEvent();
  void              FinishEvent();
  void              FinishGlobal();

  Int_t             fkMASK[MAXSTREAMSTOMERGE];  //! masks for track ids from
                                              //  different source files
  Int_t             fkMASKSTEP;           // step to increase MASK for
                                          // each input file
  TString           fOutputFileName;      // output file name
  TString           fOutputDirName;       // output dir name

  Int_t             fEvent;               // output event nr.
  Int_t             fNrOfEventsToWrite;   // Nr of events to write
  Int_t             fNrOfEventsWritten;   // Nr of events written
  Int_t             fCopyTreesFromInput;  // from which input file the trees
                                          // should be copied, -1 for no copies
  Int_t             fNinputs;             // nr of input streams - can be taken from the TClonesArray dimension
  Int_t             fNinputsGiven;        // nr of input streams given by user
  TClonesArray *    fInputStreams;        // input signal streams

//  AliStream*        fOutputStream;
  AliRunLoader*     fOutRunLoader;        //!
  Bool_t            fOutputInitialized;   //indicates if outout was initialized 
                                          // 
  AliMergeCombi *   fCombi;               // pointer to the combination object
  TArrayI           fCombination;         //! combination of events from
  TString           fCombinationFileName; // fn with combinations (used
                                          // with type 2 of comb.)

  AliRunLoader*     GetOutRunLoader();
  
  static const TString fgkDefOutFolderName;//default name for output foler 
  static const TString fgkBaseInFolderName;//default name for input foler
  ClassDef(AliRunDigitizer,6)
};

#endif // ALIRUNDIGITIZER_H
