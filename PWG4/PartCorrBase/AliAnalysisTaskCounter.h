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
  
 private:  
  AliAnalysisTaskCounter(const AliAnalysisTaskCounter&); // not implemented  
  AliAnalysisTaskCounter& operator=(const AliAnalysisTaskCounter&); // not implemented
  
 public: 
  
  virtual void UserCreateOutputObjects();  
  virtual void UserExec(Option_t *option);  
  virtual void FinishTaskOutput();  
  
  void SetTrackMultiplicityEtaCut(Float_t eta) { fTrackMultEtaCut   = eta    ; }  
  void SetZVertexCut(Float_t vcut)             { fZVertexCut        = vcut   ; }  

  void SwitchOnCaloFilterPatch()               { fCaloFilterPatch   = kTRUE  ; } 
  void SwitchOffCaloFilterPatch()              { fCaloFilterPatch   = kFALSE ; }  
  Bool_t IsCaloFilterPatchOn()                 { return fCaloFilterPatch     ; }   
  
  void AcceptFastCluster()                     { fAcceptFastCluster = kTRUE  ; } 
  void RejectFastCluster()                     { fAcceptFastCluster = kFALSE ; }  
  Bool_t IsFastClusterAccepted()               { return fAcceptFastCluster   ; }   
  
  Bool_t CheckForPrimaryVertex() ;
   
 private: 
  Bool_t               fAcceptFastCluster; // Accept events from fast cluster, exclude thiese events for LHC11a
  Float_t              fZVertexCut;        // Z vertex cut  
  Float_t              fTrackMultEtaCut;   // Track multiplicity eta cut  
  Bool_t               fCaloFilterPatch;   // CaloFilter patch  
  TList*               fOutputContainer;   //! Histogram container  
  AliESDtrackCuts    * fESDtrackCuts;      // Track cut    
  AliTriggerAnalysis * fTriggerAnalysis;   // Trigger algorithm 
  
  //Histograms
  TH1I *  fhNEvents;      //! Events that delivers the analysis frame after different assumptions  
  TH1F *  fhXVertex;      //! X Vertex distribution
  TH1F *  fhYVertex;      //! Y Vertex distribution
  TH1F *  fhZVertex;      //! Z Vertex distribution
  TH1F *  fhXGoodVertex;  //! X Vertex good distribution
  TH1F *  fhYGoodVertex;  //! Y Vertex good distribution
  TH1F *  fhZGoodVertex;  //! Z Vertex good distribution  

  ClassDef(AliAnalysisTaskCounter, 1);

};

#endif //ALIANALYSISTASKCOUNTER_H
