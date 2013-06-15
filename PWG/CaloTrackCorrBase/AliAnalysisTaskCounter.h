#ifndef ALIANALYSISTASKCOUNTER_H
#define ALIANALYSISTASKCOUNTER_H

//_________________________________________________________________________
//
// Count events with different selections
//
// Author: Gustavo Conesa Balbastre (LPSC)
//
//_________________________________________________________________________

class TH1F;
class TList;
class AliESDtrackCuts;
class AliTriggerAnalysis;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCounter : public AliAnalysisTaskSE { 
  
 public:  
  AliAnalysisTaskCounter();  
  AliAnalysisTaskCounter(const char *name);  
  virtual ~AliAnalysisTaskCounter() ;
  
  virtual void   UserCreateOutputObjects();  
  virtual void   UserExec(Option_t *option);
  virtual void   FinishTaskOutput();
  
  static  Bool_t PythiaInfoFromFile(TString currFile, Float_t & xsec, Float_t & trials) ;
  virtual Bool_t Notify();

  void    SetTrackMultiplicityEtaCut(Float_t eta) { fTrackMultEtaCut   = eta    ; }  
  void    SetZVertexCut(Float_t vcut)             { fZVertexCut        = vcut   ; }  

  void    SwitchOnCaloFilterPatch()               { fCaloFilterPatch   = kTRUE  ; } 
  void    SwitchOffCaloFilterPatch()              { fCaloFilterPatch   = kFALSE ; }  
  Bool_t  IsCaloFilterPatchOn()            const  { return fCaloFilterPatch     ; }   
  
  void    AcceptFastCluster()                     { fAcceptFastCluster = kTRUE  ; } 
  void    RejectFastCluster()                     { fAcceptFastCluster = kFALSE ; }  
  Bool_t  IsFastClusterAccepted()       const     { return fAcceptFastCluster   ; }   
  
  Bool_t  CheckForPrimaryVertex() ;

 private: 
  Bool_t               fAcceptFastCluster; // Accept events from fast cluster, exclude thiese events for LHC11a
  Float_t              fZVertexCut;        // Z vertex cut  
  Float_t              fTrackMultEtaCut;   // Track multiplicity eta cut  
  Float_t              fAvgTrials;         // avg trials
  Bool_t               fCaloFilterPatch;   // CaloFilter patch  
  TList*               fOutputContainer;   //! Histogram container  
  AliESDtrackCuts    * fESDtrackCuts;      // Track cut    
  AliTriggerAnalysis * fTriggerAnalysis;   // Trigger algorithm 
  TString              fCurrFileName;      // current file path name
  
  //Histograms
  TH1I *  fhNEvents;      //! Events that delivers the analysis frame after different assumptions  
  TH1F *  fhXVertex;      //! X Vertex distribution
  TH1F *  fhYVertex;      //! Y Vertex distribution
  TH1F *  fhZVertex;      //! Z Vertex distribution
  TH1F *  fhXGoodVertex;  //! X Vertex good distribution
  TH1F *  fhYGoodVertex;  //! Y Vertex good distribution
  TH1F *  fhZGoodVertex;  //! Z Vertex good distribution  
  TH1F *  fhCentrality;   //! centrality
  TH1F *  fhEventPlaneAngle; //! Histogram with Event plane angle

  TH1F *  fh1Xsec ;                       //! Xsec pythia
  TH1F *  fh1Trials ;                     //! trials pythia
  
  AliAnalysisTaskCounter(           const AliAnalysisTaskCounter&); // not implemented  
  AliAnalysisTaskCounter& operator=(const AliAnalysisTaskCounter&); // not implemented
  
  ClassDef(AliAnalysisTaskCounter, 4);

};

#endif //ALIANALYSISTASKCOUNTER_H
