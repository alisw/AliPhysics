#ifndef ALIANALYSISTASKCOUNTER_H
#define ALIANALYSISTASKCOUNTER_H

//_________________________________________________________________________
/// \class AliAnalysisTaskCounter
/// \ingroup CaloTrackCorrelationsBase
/// \brief Count events with different selection criteria
///
/// It produces a histogram, fhNEvents, with the number of events with 9 bins
/// representing different selection criteria:
/// * 1: all events (that passed the physics selection if it was on) 
/// * 2: same but cross check that event pointer did exist (not really necessary) and if the trigger cluster is not FAST if *fAcceptFastCluster=1*
/// * 3: passes z vertex cut, settable *fZVertexCut* (10 cm) 
/// * 4: passes track multiplicity cut, at least one track in: 
///   * eta < 0.8
///   * passing *AliESDtrackCuts fESDtrackCuts* cuts for ESDs or Hybrid tracks for AODs
/// * 5: 3 && 4
/// * 6: pass V0AND
/// * 7: 6 && 3
/// * 8: 6 && 4
/// * 9: 6 && 5
/// * 10: not pileup from SPD *event->IsPileupFromSPD(3, 0.8, 3., 2., 5.)*
/// * 11: Good (primary, not 0,0,0) vertex
/// * 12: 11 && 3
/// * 13: 11 && 4
/// * 14: 11 && 6
/// * 15: 11 && 9
/// * 16: 11 && 10
/// * 17: 10 && 6
/// * 18: special case of EMCal bad events, cluster with too large energy and cells
/// * 19: 18 & 3
/// * 20: special case of EMCal bad events, !18 & too many cells with E avobe a cut
/// * 21: 20 & 3
///
/// Other histograms:
/// * event vertex X,Y,Z before and after event selection
/// * event plane, centrality
///
/// This class also recovers the cross section and number of trials
/// in case of MC PYTHIA productions done in pT-hard bins, and stores
/// them in a histogram *fh1Xsec* and *fh1Trials*
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
///         
//_________________________________________________________________________


class TH1F;
class TList;
class AliESDtrackCuts;
class AliGenPythiaEventHeader;
//class AliTriggerAnalysis;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCounter : public AliAnalysisTaskSE { 
  
 public:
  
  AliAnalysisTaskCounter();  
  
  AliAnalysisTaskCounter(const char *name);  
  
  virtual       ~AliAnalysisTaskCounter();
  
  virtual void   UserCreateOutputObjects();
  
  virtual void   UserExec(Option_t *option);
  
  virtual void   FinishTaskOutput();
  
  virtual Bool_t Notify();
  
  static  Bool_t PythiaInfoFromFile(TString currFile, Float_t & xsec, Float_t & trials) ;

  void           SetTrackMultiplicityEtaCut(Float_t eta) { fTrackMultEtaCut   = eta    ; }
  void           SetZVertexCut(Float_t vcut)             { fZVertexCut        = vcut   ; }
  
  void           AcceptFastCluster()                     { fAcceptFastCluster = kTRUE  ; }
  void           RejectFastCluster()                     { fAcceptFastCluster = kFALSE ; }
  Bool_t         IsFastClusterAccepted()       const     { return fAcceptFastCluster   ; }
  
  Bool_t         CheckForPrimaryVertex() ;

  void           SwitchOnMCCrossSectionCalculation()     { fCheckMCCrossSection = kTRUE  ; }
  void           SwitchOffMCCrossSectionCalculation()    { fCheckMCCrossSection = kFALSE ; }
  
  void           SwitchOnAliCentrality ()                { fUseAliCentrality    = kTRUE  ; }
  void           SwitchOffAliCentrality()                { fUseAliCentrality    = kFALSE ; }
  
  void           SetCentralityClass(TString name)        { fCentralityClass   = name     ; }
  TString        GetCentralityClass()              const { return fCentralityClass       ; }

 private:
  
  Bool_t               fAcceptFastCluster;   ///< Accept events from fast cluster, exclude these events for LHC11a.
  Float_t              fZVertexCut;          ///< Z vertex cut.  
  Float_t              fTrackMultEtaCut;     ///< Track multiplicity eta cut.  
  Float_t              fAvgTrials;           ///< Average number of event trials.
  TList*               fOutputContainer;     //!<! Histogram container.  
  AliESDtrackCuts    * fESDtrackCuts;        ///< Track cut.    
//AliTriggerAnalysis * fTriggerAnalysis;     ///< Trigger algorithm.
  TString              fCurrFileName;        ///< Current file path name.
  Bool_t               fCheckMCCrossSection; ///< Retrieve from the pyxsec.root file only if requested.
  Bool_t               fUseAliCentrality;    ///< Use the centrality estimator from AliCentrality or AliMultSelection
  TString              fCentralityClass;     ///< Multiplicity percentile/centrality estimator, for ex. V0M
  
  //
  // Histograms
  //
  TH1I *  fhNEvents;         //!<! Events that delivers the analysis frame after different assumptions.  
  TH1F *  fhXVertex;         //!<! X Vertex distribution.
  TH1F *  fhYVertex;         //!<! Y Vertex distribution.
  TH1F *  fhZVertex;         //!<! Z Vertex distribution.
  TH1F *  fhXGoodVertex;     //!<! X Vertex distribution, after event selection.
  TH1F *  fhYGoodVertex;     //!<! Y Vertex distribution, after event selection.
  TH1F *  fhZGoodVertex;     //!<! Z Vertex distribution, after event selection.  
  TH1F *  fhCentrality;      //!<! Centrality.
  TH1F *  fhEventPlaneAngle; //!<! Event plane angle.

  TH1F *  fh1Xsec ;          //!<! Cross section in PYTHIA.
  TH1F *  fh1Trials ;        //!<! Number of event trials in PYTHIA.
  
  /// Copy constructor not implemented.
  AliAnalysisTaskCounter(           const AliAnalysisTaskCounter&); 
  
  /// Assignment operator not implemented.
  AliAnalysisTaskCounter& operator=(const AliAnalysisTaskCounter&); 
  
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskCounter, 7) ;
  /// \endcond
  
};

#endif //ALIANALYSISTASKCOUNTER_H
