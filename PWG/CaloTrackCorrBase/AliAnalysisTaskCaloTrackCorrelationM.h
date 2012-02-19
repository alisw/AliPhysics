#ifndef ALIANALYSISTASKCALOTRACKCORRELATIONM_H
#define ALIANALYSISTASKCALOTRACKCORRELATIONM_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//_________________________________________________________________________
// Analysis task that executes the analysis classes
// that depend on the CaloTrackCorr frame, frame for Particle identification 
// with calorimeters and tracks and correlations.
// Specially designed for calorimeters but also can be used for charged tracks
// Input of this task is a configuration file that contains all the settings 
// of the analysis
//
// -- Author: Gustavo Conesa (INFN-LNF, LPSC-Grenoble)

//--- Root ---
class TList;

//--- AliRoot ---
#include "AliAnalysisTaskME.h"
class AliAnaCaloTrackCorrMaker;
class AliMixedEvent; 
class AliMCEvent; 

class AliAnalysisTaskCaloTrackCorrelationM : public AliAnalysisTaskME
{
 public:
  
  AliAnalysisTaskCaloTrackCorrelationM();
  AliAnalysisTaskCaloTrackCorrelationM(const char* name);
  virtual ~AliAnalysisTaskCaloTrackCorrelationM() ; // virtual dtor
  
  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() ;
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  
  void         SetConfigFileName(TString & name ) { fConfigName = name ; }
  TString      GetConfigFileName()          const { return fConfigName ; }
	
  void         SetAnalysisMaker(AliAnaCaloTrackCorrMaker * const maker) { fAna = maker ; } 
  AliMixedEvent * InputEvent()                    { return fInputEvent ; }
  AliMCEvent*  MCEvent()                    const { return NULL        ; } // do something about MC event 

 private:
  
  AliAnalysisTaskCaloTrackCorrelationM(           const AliAnalysisTaskCaloTrackCorrelationM&); // Not implemented
  AliAnalysisTaskCaloTrackCorrelationM& operator=(const AliAnalysisTaskCaloTrackCorrelationM&); // Not implemented
  
  AliAnaCaloTrackCorrMaker* fAna;  //  Pointer to the manager class 
  TList * fOutputContainer ;       //! Histogram container
  TString fConfigName ;            //  Configuration file name
  TList * fCuts ;                  //! List with analysis cuts
  AliMixedEvent * fInputEvent;     //  Mixed event access pointer
  	
  ClassDef(AliAnalysisTaskCaloTrackCorrelationM, 3); // Analysis task for standard gamma correlation analysis
};

#endif //ALIANALYSISTASKCALOTRACKCORRELATIONM_H
