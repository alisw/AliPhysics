#ifndef ALIANALYSISTASKCALOTRACKCORRELATION_H
#define ALIANALYSISTASKCALOTRACKCORRELATION_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//_________________________________________________________________________
/// \class AliAnalysisTaskCaloTrackCorrelation
/// \brief Main class conecting the CaloTrackCorrelations package and Analysis Frame
///
/// Analysis task that executes the analysis classes
/// that depend on the CaloTrackCorr frame, frame for Particle identification 
/// with calorimeters and tracks and correlations.
/// Specially designed for calorimeters but also can be used for charged tracks
/// Input of this task is a configuration file that contains all the settings 
/// of the analysis.
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
///_________________________________________________________________________


//--- Root ---
class TList;

//--- AliRoot ---
#include "AliAnalysisTaskSE.h"
class AliAnaCaloTrackCorrMaker;
class AliESDEvent;
class AliAODEvent;

class AliAnalysisTaskCaloTrackCorrelation : public AliAnalysisTaskSE
{
 public:
  
  AliAnalysisTaskCaloTrackCorrelation() ;
  
  AliAnalysisTaskCaloTrackCorrelation(const char* name) ;
  
  virtual ~AliAnalysisTaskCaloTrackCorrelation() ;
  
  // Implementation of interface methods
  
  virtual void UserCreateOutputObjects();
  
  virtual void Init();
  
  virtual void LocalInit() ;
  
  virtual void UserExec (Option_t *option);
  
  virtual void Terminate(Option_t *option);
  
  virtual void FinishTaskOutput();
  
  // Setters/Getters
  
  AliAnaCaloTrackCorrMaker* GetAnalysisMaker()    { return fAna          ; }
  void         SetAnalysisMaker(AliAnaCaloTrackCorrMaker * const maker) 
                                                  { fAna = maker         ; } 
  
  void         SetConfigFileName(TString & name ) { fConfigName = name   ; }
  TString      GetConfigFileName()          const { return fConfigName   ; }

  void         SetLastEvent(Int_t event )         { fLastEvent = event   ; }
  Int_t        GetLastEvent()               const { return fLastEvent    ; }

  
 private:
  
  /// Copy constructor not implemented.
  AliAnalysisTaskCaloTrackCorrelation           (const AliAnalysisTaskCaloTrackCorrelation&) ;
  
  /// Assignment operator not implemented.
  AliAnalysisTaskCaloTrackCorrelation& operator=(const AliAnalysisTaskCaloTrackCorrelation&) ; 
  
  AliAnaCaloTrackCorrMaker* fAna;  ///<  Pointer to the manager class. 
  
  TList * fOutputContainer ;       //!<! Histogram container.
  
  TString fConfigName ;            ///<  Configuration file name.
  
  TList * fCuts ;                  //!<! List with analysis cuts.
  
  Int_t   fLastEvent;              //!<! Analyze all the events until this one, for testing.    
  
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskCaloTrackCorrelation, 4) ; 
  /// \endcond

};

#endif //ALIANALYSISTASKCALOTRACKCORRELATION_H
