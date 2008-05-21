#ifndef ALIFMDBASEDA_H
#define ALIFMDBASEDA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
//
//This class provides a base interface for the Detector Algorithms (DA) of the FMD. 
//At least three implementations are needed: AliFMDPedestalDA, AliFMDGainDA and AliFMDPhysicsDA . 
//These classes will provide the calibration data for the AliFMDPreprocessor to be used in the shuttle.
//The input for this class are raw data (AliRawReader) and the output is a comma-separated file
//(std::ofstream) that contains the values defined in the implementations of this class.
//
//Author: Hans Hjersing Dalsgaard, hans.dalsgaard@cern.ch
//

#include "TNamed.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "iostream"
#include "fstream"
#include "TString.h"
#include "AliRawReader.h"
#include "AliFMDDigit.h"
#include "AliFMDParameters.h"


class AliFMDBaseDA: public TNamed {
  
 public:
  AliFMDBaseDA() ;
  AliFMDBaseDA(const AliFMDBaseDA & baseDA) ;
  //  AliFMDBaseDA& operator = (const AliFMDBaseDA & baseDA) ; 
  
  ~AliFMDBaseDA() ;
  
  void Run(AliRawReader* fmdReader);
  void SetSaveDiagnostics(Bool_t save) {fSaveHistograms = save;}
  void SetRequiredEvents(Int_t nEvents) {fRequiredEvents = nEvents;}
  Int_t GetRequiredEvents() {return fRequiredEvents ;}
 protected:
  
  virtual void Init()  {};
  virtual void FillChannels(AliFMDDigit* )  {};
  virtual void Analyse(UShort_t, Char_t, UShort_t, UShort_t )  {};
  virtual void WriteHeaderToFile()  {};
  virtual void AddChannelContainer(TObjArray*, UShort_t, Char_t, UShort_t, UShort_t )  {};
  virtual void FinishEvent()  {};
  
  
  Int_t GetCurrentEvent() {return fCurrentEvent;}
  
  static const UInt_t kBaseDDL = 3072;
  //Char_t* fDiagnosticsFilename;
  TString fDiagnosticsFilename;
  std::ofstream fOutputFile;
  std::ofstream fConditionsFile;
  Bool_t fSaveHistograms;
  TObjArray fDetectorArray;
  
  
 private:
 
  void WriteConditionsData();
  void SetCurrentEvent(Int_t currentEvent) {fCurrentEvent = currentEvent; }
  void InitContainer();
  Int_t fRequiredEvents;
  Int_t fCurrentEvent;   
   
  ClassDef(AliFMDBaseDA,0)

};
#endif

