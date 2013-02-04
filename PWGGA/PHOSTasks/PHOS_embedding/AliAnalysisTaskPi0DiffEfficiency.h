#ifndef ALIANALYSISTASKPI0DIFFEFFICIENCY_H
#define ALIANALYSISTASKPI0DIFFEFFICIENCY_H

// Task calculating effciency of pi0 registration 
// as a difference between spectrum after and before 
// embedding

class TClonesArray;
class AliAODCaloCluster ;
class AliPHOSAodCluster ;
class AliPHOSCalibData ;
class AliPHOSGeometry ;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskPi0DiffEfficiency : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskPi0DiffEfficiency(const char *name = "AliAnalysisTaskPi0DiffEfficiency");
  virtual ~AliAnalysisTaskPi0DiffEfficiency() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);

  void SetPHOSBadMap(Int_t mod,TH2I * h)
  {
    if(fPHOSBadMap[mod]) delete fPHOSBadMap[mod] ;
    fPHOSBadMap[mod]=new TH2I(*h) ;
    printf("Set %s \n",fPHOSBadMap[mod]->GetName());
  }

  protected:
  AliAnalysisTaskPi0DiffEfficiency(const AliAnalysisTaskPi0DiffEfficiency& a) ; // not implemented
  AliAnalysisTaskPi0DiffEfficiency& operator=(const AliAnalysisTaskPi0DiffEfficiency& a ) ; // not implemented
  Bool_t IsGoodChannel(const char * det, Int_t mod,Int_t ix, Int_t iz); //Use addisional bad map for PHOS
  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z, Double_t weight) const ; //Fill 3D histogram witn name key
  Bool_t TestLambda(Double_t e,Double_t l1,Double_t l2) ;  //Evaluate Dispersion cuts for photons
  Bool_t TestLambda2(Double_t e,Double_t l1,Double_t l2) ;  //Evaluate Dispersion cuts for photons
  Double_t TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge);
  Bool_t TestTOF(Double_t t, Double_t e) ;
  
  void ProcessMC() ;
  Double_t CoreEnergy(AliPHOSAodCluster * clu); 
   
private:
  Bool_t IsSameCluster(AliAODCaloCluster * c1,AliAODCaloCluster * c2)const ;
  void EvalLambdas(AliAODCaloCluster * clu, Int_t iR,Double_t &m02, Double_t &m20) ;
  Double_t PrimaryWeight(Double_t x) ;
 
private:
  AliStack * fStack ;               //! Stack of primary
  THashList * fOutputContainer;     //! final histogram container
  TList * fPHOSEvents[1][10][11] ;  //! Containers for events with PHOS photons
  TClonesArray * fPHOSEvent ;       //!PHOS photons in current event
  TClonesArray * fPHOSEvent1 ;      //!PHOS photons in current event
  TClonesArray * fPHOSEvent2 ;      //!PHOS photons in current event
  AliPHOSCalibData *fPHOSCalibData; //! PHOS calibration object
  TF1 *fNonLinCorr;                 // Non-linearity correction
 
  //Reaction plane for v2
  Float_t fRPfull ; //!Reaction plain calculated with full TPC 
  Float_t fRPA ;    //!Reaction plain calculated with A-side TPC: eta>0.15 
  Float_t fRPC ;    //!Reaction plain calculated with C-side TPC: eta<-0.15
  Float_t fRPFar ;  //!Reaction plain calculated with TPC: eta>0.6 
  Float_t fRPAFar ; //!Reaction plain calculated with A-side TPC: eta>0.6 
  Float_t fRPCFar ; //!Reaction plain calculated with C-side TPC: eta<-0.6

  Float_t fCentrality ; //!Centrality of the currecnt event

  Int_t fCenBin ;       //! Current centrality bin

  TH2I *fPHOSBadMap[6] ;    //!Container for PHOS bad channels map

  AliPHOSGeometry  *fPHOSGeo;  //! PHOS geometry
  Int_t fEventCounter;         //! number of analyzed events
 
  ClassDef(AliAnalysisTaskPi0DiffEfficiency, 2); // PHOS analysis task
};

#endif
