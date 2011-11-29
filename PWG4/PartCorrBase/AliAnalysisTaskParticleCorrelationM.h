#ifndef ALIANALYSISTASKPARTICLECORRELATIONM_H
#define ALIANALYSISTASKPARTICLECORRELATIONM_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//_________________________________________________________________________
// Analysis task that executes the analysis classes
// that depend on the PartCorr frame, frame for Particle identification and correlations.
// Specially designed for calorimeters but also can be used for charged tracks
// Input of this task is a configuration file that contains all the settings of the analyis
//
// -- Author: Gustavo Conesa (INFN-LNF)

//--- Root ---
class TList;

//--- AliRoot ---
#include "AliAnalysisTaskME.h"
class AliAnaPartCorrMaker;
class AliMixedEvent; 
class AliMCEvent; 

class AliAnalysisTaskParticleCorrelationM : public AliAnalysisTaskME
{
 public:
  AliAnalysisTaskParticleCorrelationM();
  AliAnalysisTaskParticleCorrelationM(const char* name);
  virtual ~AliAnalysisTaskParticleCorrelationM() ;// virtual dtor
  
  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() ;
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  
  void SetConfigFileName(TString & name ) {fConfigName = name ; }
  TString GetConfigFileName() const {return fConfigName ; }
	
  void SetAnalysisMaker(AliAnaPartCorrMaker * const maker) {fAna = maker;} 
  AliMixedEvent * InputEvent(){ return fInputEvent ; }
  AliMCEvent*   MCEvent()     {return NULL;} // do something about MC event 

 private:
  AliAnalysisTaskParticleCorrelationM(const AliAnalysisTaskParticleCorrelationM&); // Not implemented
  AliAnalysisTaskParticleCorrelationM& operator=(const AliAnalysisTaskParticleCorrelationM&); // Not implemented
  
  AliAnaPartCorrMaker* fAna;  //  Pointer to the manager class 
  TList * fOutputContainer ;  //! Histogram container
  TString fConfigName ;       // Configuration file name
  TList * fCuts ;             //! List with analysis cuts
  
  AliMixedEvent * fInputEvent;
  	
  ClassDef(AliAnalysisTaskParticleCorrelationM, 3); // Analysis task for standard gamma correlation analysis
};

#endif //ALIANALYSISTASKPARTICLECORRELATIONM_H
