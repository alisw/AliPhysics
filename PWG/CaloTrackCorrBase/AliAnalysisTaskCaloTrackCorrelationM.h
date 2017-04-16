#ifndef ALIANALYSISTASKCALOTRACKCORRELATIONM_H
#define ALIANALYSISTASKCALOTRACKCORRELATIONM_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//_________________________________________________________________________
/// \class AliAnalysisTaskCaloTrackCorrelationM
/// \ingroup CaloTrackCorrelationsBase 
/// \brief Main class conecting the CaloTrackCorrelations package and Analysis Mixing Frame
///
/// Analysis task that executes the analysis classes
/// that depend on the CaloTrackCorr frame, frame for Particle identification 
/// with calorimeters and tracks and correlations.
/// Specially designed for calorimeters but also can be used for charged tracks
/// Input of this task is a configuration file that contains all the settings 
/// of the analysis.
///
/// **This task was developped for mixing studies (not used since 2011) ...**
///
/// More information can be found in this [twiki](https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PhotonHadronCorrelations).
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
///_________________________________________________________________________



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
  
  AliAnalysisTaskCaloTrackCorrelationM() ;
  
  AliAnalysisTaskCaloTrackCorrelationM(const char* name) ;
  
  virtual ~AliAnalysisTaskCaloTrackCorrelationM() ; 
  
  // Implementation of interface methods
  
  virtual void UserCreateOutputObjects();
  
  virtual void Init();
  
  virtual void LocalInit() ;
  
  virtual void UserExec(Option_t *option);
  
  virtual void Terminate(Option_t *option);
  
  // Setters/Getters
  
  void           SetConfigFileName(TString & name ) { fConfigName = name ; }
  TString        GetConfigFileName()          const { return fConfigName ; }
	
  AliAnaCaloTrackCorrMaker* GetAnalysisMaker()      { return fAna        ; }
  void           SetAnalysisMaker(AliAnaCaloTrackCorrMaker * const maker) 
                                                    { fAna = maker       ; } 
  
  AliMixedEvent* InputEvent()                       { return fInputEvent ; }
  AliMCEvent*    MCEvent()                    const { return NULL        ; }  

 private:
  
  /// Copy constructor not implemented.
  AliAnalysisTaskCaloTrackCorrelationM(           const AliAnalysisTaskCaloTrackCorrelationM&) ; 
  
  /// Assignment operator not implemented.
  AliAnalysisTaskCaloTrackCorrelationM& operator=(const AliAnalysisTaskCaloTrackCorrelationM&) ; 
  
  AliAnaCaloTrackCorrMaker* fAna;  ///<  Pointer to the manager class 
  
  TList * fOutputContainer ;       //!<! Histogram container
  
  TString fConfigName ;            ///<  Configuration file name
  
  TList * fCuts ;                  //!<! List with analysis cuts
  
  AliMixedEvent * fInputEvent;     ///<  Mixed event access pointer
  	
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskCaloTrackCorrelationM, 3) ;   
  /// \endcond

};

#endif //ALIANALYSISTASKCALOTRACKCORRELATIONM_H
