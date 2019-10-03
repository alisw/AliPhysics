#ifndef AliAnalysisTaskGammaFlow_h
#define AliAnalysisTaskGammaFlow_h

class TObjArray;
class THashList ;
class TH1F;
class TH1D ;
class TH2I;
class TH2F;
class TH3F;
class TF1 ;
class AliStack ;
class AliESDtrackCuts;
class AliPHOSGeometry;
class AliAODEvent ;
class AliPHOSCalibData;
class AliAODtrack ;
class AliAODCaloCluster ;
class AliEPFlattener ;
class AliFlowTrackCuts ;
class AliFlowEvent ;
class AliFlowVector ;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskGammaFlow : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskGammaFlow(const char *name = "AliAnalysisTaskGammaFlow");
  virtual ~AliAnalysisTaskGammaFlow() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void SetHarmonics(Int_t n=2){fHarmonics=n ;}
  void SetDistCut(Bool_t on=kTRUE){fDistCut=on ;}
  
private:
  AliAnalysisTaskGammaFlow(const AliAnalysisTaskGammaFlow&); // not implemented
  AliAnalysisTaskGammaFlow& operator=(const AliAnalysisTaskGammaFlow&); // not implemented

  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z, Double_t w) const ; //Fill 3D histogram witn name key

  Bool_t TestPHOSEvent(AliAODEvent * event) ; //Checks if PHOS event complete
  Bool_t TestLambda(Double_t pt,Double_t l1,Double_t l2) ;  //Evaluate Dispersion cuts for photons
  Bool_t TestLambda2(Double_t pt,Double_t l1,Double_t l2) ;  //Evaluate Dispersion cuts for photons
  Double_t TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge);
  Bool_t TestTOF(Double_t t, Double_t e) ;
  Int_t ConvertRunNumber(Int_t run) ; 

  void EvalV0ReactionPlane(AliAODEvent * event) ;
//  Double_t ApplyFlattening(Double_t phi, Double_t c) ; //Apply centrality-dependent flattening
//  Double_t ApplyFlatteningV0A(Double_t phi, Double_t c) ; //Apply centrality-dependent flattening
//  Double_t ApplyFlatteningV0C(Double_t phi, Double_t c) ; //Apply centrality-dependent flattening
  void ApplyFinalFlattening() ;//apply final fine flattening
  void ApplyFinalQFlattening() ;//apply final fine flattening
  Bool_t GetTPCEventPlane(Double_t &epAngle, Double_t &qsubRes) ;
  TObjArray* GetEventPlaneTracks(Int_t &maxID) ;
  Double_t GetWeight(TObject* track1) ;
  Double_t GetPhiWeight(TObject* track1) ;

  void EvalResolution() ;
  void EvalQResolution() ;
  Double_t PHOSMultiplicity();
  Double_t CentralityWeight(Double_t c) ;
  
private:
  THashList * fOutputContainer;    //final histogram container
  TList * fPHOSEvents[1][10][11] ; // Containers for events with PHOS photons
  TClonesArray * fPHOSEvent ;      // PHOS photons in current event
 
  //Vector value
  Double_t fQV0A ;   //Lengths of flow vectors
  Double_t fQV0C ;
  Double_t fQTPC ;
 
  //Reaction plain for v2
  Float_t fRP ;        //! Reaction plane calculated with full TPC 
  Float_t fRPV0A ;     //! Reaction plain calculated with A-side TPC: eta>0.15 
  Float_t fRPV0C ;     //! Reaction plain calculated with C-side TPC: eta<-0.15
  Float_t fRPQ ;       //! Reaction plane calculated with full TPC 
  Float_t fRPQV0A ;    //! Reaction plain calculated with A-side TPC: eta>0.15 
  Float_t fRPQV0C ;    //! Reaction plain calculated with C-side TPC: eta<-0.15
  Double_t fV0Ares ;   //! resolution of EP
  Double_t fV0Cres ;   //! resolution of EP
  Double_t fTPCres ;   //! resolution of EP
  Double_t fV0AQres ;  //! resolution of EP
  Double_t fV0CQres ;  //! resolution of EP
  Double_t fTPCQres ;  //! resolution of EP
  Bool_t fHaveTPCRP ;  //! Is TPC RP defined?
  TH1F * fPhiDist ;    //!
  AliEPFlattener * fV0AFlat ; //!
  AliEPFlattener * fV0CFlat ; //!
  AliEPFlattener * fTPCFlat ; //!
  AliEPFlattener * fV0AQFlat ; //!
  AliEPFlattener * fV0CQFlat ; //!
  AliEPFlattener * fTPCQFlat ; //!
  
  Int_t fHarmonics ; //Harminic in use
  Bool_t fDistCut ;

  Int_t fRunNumber ;    //Current run number
  Float_t fCentrality ; //!Centrality of the currecnt event

  Int_t fCenBin ;       //! Current centrality bin

  AliPHOSGeometry  *fPHOSGeo;  //! PHOS geometry
  Int_t fEventCounter;         // number of analyzed events

  TH1D*  fTPCfinalC2 ;   //! HIstos with flattening parameters
  TH1D*  fTPCfinalS2 ;   //!
  TH1D*  fTPCfinalC4 ;   //!
  TH1D*  fTPCfinalS4 ;   //!

  TH1D*  fV0AfinalC2 ;   //!
  TH1D*  fV0AfinalS2 ;   //!
  TH1D*  fV0AfinalC4 ;   //!
  TH1D*  fV0AfinalS4 ;   //!
    
  TH1D*  fV0CfinalC2 ;   //!
  TH1D*  fV0CfinalS2 ;   //!
  TH1D*  fV0CfinalC4 ;   //!
  TH1D*  fV0CfinalS4 ;   //!

  TH1D*  fTPCfinalQC2 ;   //! HIstos with flattening parameters
  TH1D*  fTPCfinalQS2 ;   //!
  TH1D*  fTPCfinalQC4 ;   //!
  TH1D*  fTPCfinalQS4 ;   //!

  TH1D*  fV0AfinalQC2 ;   //!
  TH1D*  fV0AfinalQS2 ;   //!
  TH1D*  fV0AfinalQC4 ;   //!
  TH1D*  fV0AfinalQS4 ;   //!
    
  TH1D*  fV0CfinalQC2 ;   //!
  TH1D*  fV0CfinalQS2 ;   //!
  TH1D*  fV0CfinalQC4 ;   //!
  TH1D*  fV0CfinalQS4 ;    //!
  AliFlowTrackCuts * fCutsV0 ; //! 
  AliFlowTrackCuts * fCutsTPC ; //!
  AliFlowEvent     * fFlowEvent ; //!
  
  ClassDef(AliAnalysisTaskGammaFlow, 2); // PHOS analysis task
};

#endif
