#ifndef ALIANAPARTCORRMAKER_H
#define ALIANAPARTCORRMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//_________________________________________________________________________
// Steering class for particle (gamma, hadron) identification and correlation analysis
// It is called by the task class AliAnalysisTaskParticleCorrelation and it connects the input 
// (ESD/AOD/MonteCarlo) got with AliCaloTrackReader (produces TClonesArrays of AODs 
// (TParticles in MC case if requested)), with the 
// analysis classes that derive from AliAnaPartCorrBaseClass
//
// -- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
class TList; 
class TClonesArray;
#include<TObject.h>
class TString;

// --- Analysis system ---
class AliCaloTrackReader ;

class AliAnaPartCorrMaker : public TObject {

 public: 
  
  AliAnaPartCorrMaker() ; // default ctor
  AliAnaPartCorrMaker(const AliAnaPartCorrMaker & g) ; // cpy ctor
  AliAnaPartCorrMaker & operator = (const AliAnaPartCorrMaker & g) ;//cpy assignment
  virtual ~AliAnaPartCorrMaker() ; //virtual dtor
  
  //Setter and getters
  TList * GetAODBranchList() ;
  TList * GetOutputContainer() ;

  Int_t GetAnaDebug() const  { return fAnaDebug ; }
  void SetAnaDebug(Int_t d)   { fAnaDebug = d ; }

  Bool_t AreHistogramsMade() const  { return fMakeHisto ; }
  void SwitchOnHistogramsMaker()   { fMakeHisto = kTRUE ; }
  void SwitchOffHistogramsMaker()   { fMakeHisto = kFALSE ; }
 
  Bool_t AreAODsMade() const  { return fMakeAOD ; }
  void SwitchOnAODsMaker()   { fMakeAOD = kTRUE ; }
  void SwitchOffAODsMaker()   { fMakeAOD = kFALSE ; }
  
  void Terminate(TList * outputList);

  void AddAnalysis(TObject* ana, Int_t n) {
    if ( fAnalysisContainer) fAnalysisContainer->AddAt(ana,n); 
    else { printf("AliAnaPartCorrMaker::AddAnalysis() - AnalysisContainer not initialized\n");
      abort();}
  }
  
  AliCaloTrackReader * GetReader() const {return fReader ; }
  void SetReader(AliCaloTrackReader * reader) { fReader = reader ; }
  
  //Others
  void Init();
  void InitParameters();
  
  void Print(const Option_t * opt) const;
  
  Bool_t ProcessEvent(Int_t iEntry) ;
  
 private:
  
  //General Data members
  
  TList * fOutputContainer ; // output histograms container
  TList * fAnalysisContainer ; // List with analysis pointers
  Bool_t  fMakeHisto ; // If true makes final analysis with histograms as output
  Bool_t  fMakeAOD ; // If true makes analysis generating AODs
  Int_t   fAnaDebug; // Debugging info.
 
  AliCaloTrackReader *  fReader ;   // Pointer to reader 
  TList * fAODBranchList ;          //! List with AOD branches created and needed in analysis
  
  ClassDef(AliAnaPartCorrMaker,2)
} ;
 

#endif //ALIANAPARTCORRMAKER_H



