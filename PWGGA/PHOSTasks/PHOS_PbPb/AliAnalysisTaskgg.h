#ifndef AliAnalysisTaskgg_cxx
#define AliAnalysisTaskgg_cxx

// Class for analysis of photon dEta-dPhi correlations
// Authors: D.Peresunko

class THashList ;
class AliPHOSGeometry;
class AliCaloPhoton ;
class AliAODTrack ;
class AliAODEvent ;
class AliEPFlattener ;
class AliPIDResponse ;
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskgg : public AliAnalysisTaskSE {
public:
    
  
  AliAnalysisTaskgg(const char *name = "AliAnalysisTaskgg");
  virtual ~AliAnalysisTaskgg() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
 
  void SelectPbPb(Bool_t isPbPb=kFALSE){fIsPbPb=isPbPb;} 

protected:  
  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z, Double_t w) const ; //Fill 3D histogram witn name key

  
  Int_t  ConvertRunNumber(Int_t run) ; 
  Int_t  FindTrackMatching(Int_t mod,TVector3 *locpos); 
  Int_t  JetRejection(Int_t module) const; //Looks is there is a jet around
  Bool_t PairCut(const AliCaloPhoton * ph1, const AliCaloPhoton * ph2, Int_t cut) const ; 
  Bool_t PHOSCut(const AliCaloPhoton * ph1, Int_t cut) const ;   
  void   ReclusterizeCPV();
  Bool_t TestCPV(Double_t emcX, Double_t emcZ, Double_t e) ;
  Bool_t TestCPVCluster(Double_t cpvX, Double_t cpvZ, Double_t emcX, Double_t emcZ, Double_t e) ; //return true if neutral


  Double_t EtaPhiWeight(Int_t kTbin, Double_t x) const ;
  
  
private:
  AliAnalysisTaskgg(const AliAnalysisTaskgg&); // not implemented
  AliAnalysisTaskgg& operator=(const AliAnalysisTaskgg&); // not implemented

protected:

  THashList *   fOutputContainer;        //final histogram container
  AliAODEvent * fEvent ;                 //!
  AliPIDResponse *fPIDResponse;     //! PID response object
  TList *       fPHOSEvents[10][10][11] ; //Containers for events with PHOS photons
  TClonesArray* fPHOSEvent ;      //! PHOS photons in current event
  TClonesArray* fCPVEvent ;       //! CPV event
 
  //Reaction plain for v2
  AliEPFlattener * fV0AFlat ; //!
  AliEPFlattener * fV0CFlat ; //!  

  Int_t fRunNumber ;    //Current run number
  Float_t fCentrality ; //!Centrality of the currecnt event
  Int_t fCenBin ;       //! Current centrality bin
  Double_t fRP ;        //! Reaction plane

  AliPHOSGeometry  *fPHOSGeo;  //! PHOS geometry
  Int_t fEventCounter;         // number of analyzed events
  Bool_t fIsPbPb ;             // Switch to PbPb or pp/pPb mode
  
  Int_t fNCuts ;   //Number of cuts
  char fCuts[120][20] ;  //Cut names
  Int_t fJetStatus[5] ; //Presence of jets around PHOS

  ClassDef(AliAnalysisTaskgg, 1); // PHOS analysis task
};

#endif
