#ifndef ALIANAPARTCORRMAKER_H
#define ALIANAPARTCORRMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//_________________________________________________________________________
// Steering class for particle (gamma, hadron) identification and correlation analysis
// It is called by the task class AliAnalysisTaskParticleCorrelation and it connects the input 
// (ESD/AOD/MonteCarlo) got with AliCaloTrackReader (produces TObjArrays of AODs 
// (TParticles in MC case if requested)), with the 
// analysis classes that derive from AliAnaPartCorrBaseClass
//
// -- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system ---
class TList; 
class TClonesArray;
#include<TObject.h>
class TString;
class TH1I;

// --- Analysis system ---
#include "AliCaloTrackReader.h" 
#include "AliCalorimeterUtils.h"

class AliAnaPartCorrMaker : public TObject {

 public: 
  AliAnaPartCorrMaker() ; // default ctor
  virtual ~AliAnaPartCorrMaker() ; //virtual dtor
  AliAnaPartCorrMaker(const AliAnaPartCorrMaker & maker) ; // cpy ctor

 private:
  AliAnaPartCorrMaker & operator = (const AliAnaPartCorrMaker & ) ;//cpy assignment

 public:
	
  //Setter and getters
  TList * GetAODBranchList() ;
  TList * GetListOfAnalysisCuts();
  TList * GetOutputContainer() ;

  Int_t GetAnaDebug() const  { return fAnaDebug ; }
  void SetAnaDebug(Int_t d)  { fAnaDebug = d ; }
	
  Bool_t AreHistogramsMade() const { return fMakeHisto ; }
  void SwitchOnHistogramsMaker()   { fMakeHisto = kTRUE ; }
  void SwitchOffHistogramsMaker()  { fMakeHisto = kFALSE ; }
 
  Bool_t AreAODsMade() const { return fMakeAOD ; }
  void SwitchOnAODsMaker()   { fMakeAOD = kTRUE ; }
  void SwitchOffAODsMaker()  { fMakeAOD = kFALSE ; }
  
  void SwitchOnMixingAnalysis() {fMakeMixing = kTRUE;} //Called by the task, no need to be set by user.
	
  void Terminate(TList * outputList);

  void AddAnalysis(TObject* ana, Int_t n) {
    if ( fAnalysisContainer) fAnalysisContainer->AddAt(ana,n); 
    else { printf("AliAnaPartCorrMaker::AddAnalysis() - AnalysisContainer not initialized\n");
      abort();}
  }
  
  AliCaloTrackReader * GetReader() {if(!fReader) fReader = new AliCaloTrackReader ();return fReader ; }
  void SetReader(AliCaloTrackReader * reader) { fReader = reader ; }
  	
  AliCalorimeterUtils * GetCaloUtils() {if(!fCaloUtils) fCaloUtils = new AliCalorimeterUtils(); return fCaloUtils ; }
  void SetCaloUtils(AliCalorimeterUtils * caloutils) { fCaloUtils = caloutils ; }
	
  //Others
  void Init();
  void InitParameters();
  
  void Print(const Option_t * opt) const;
  
  void ProcessEvent(const Int_t iEntry, const char * currentFileName) ;
  
 private:
  
  //General Data members
  
  TList * fOutputContainer ;   //! Output histograms container
  TList * fAnalysisContainer ; // List with analysis pointers
  Bool_t  fMakeHisto ;         // If true makes final analysis with histograms as output
  Bool_t  fMakeAOD ;           // If true makes analysis generating AODs
  Bool_t  fMakeMixing;         // If true it makes mixing analysis
  Int_t   fAnaDebug;           // Debugging info.
	
  AliCaloTrackReader  *  fReader ;     //  Pointer to reader 
  AliCalorimeterUtils *  fCaloUtils ;  //  Pointer to CalorimeterUtils

  TList * fAODBranchList ;     //! List with AOD branches created and needed in analysis  
  TList * fCuts ;	           //! List with analysis cuts

  TH1I  * fhNEvents;           //! Number of events counter histogram
	
  ClassDef(AliAnaPartCorrMaker,5)
} ;
 

#endif //ALIANAPARTCORRMAKER_H



