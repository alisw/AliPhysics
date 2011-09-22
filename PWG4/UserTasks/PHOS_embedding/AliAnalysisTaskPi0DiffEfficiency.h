#ifndef ALIANALYSISTASKPI0DIFFEFFICIENCY_H
#define ALIANALYSISTASKPI0DIFFEFFICIENCY_H

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

class TObjArray;
class TH1F;
class TH2I;
class TH2F;
class TH3F;
class TF1 ;
class AliStack ;
class AliESDtrackCuts;
class AliPHOSGeometry;
class AliAODEvent ;
class AliPHOSCalibData;
class AliAODTrack ;
class AliAODCaloCluster ;


#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskPi0DiffEfficiency : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskPi0DiffEfficiency(const char *name = "AliAnalysisTaskPi0DiffEfficiency");
  virtual ~AliAnalysisTaskPi0DiffEfficiency() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  void SetPHOSBadMap(Int_t mod,TH2I * h)
  {
    if(fPHOSBadMap[mod]) delete fPHOSBadMap[mod] ;
    fPHOSBadMap[mod]=new TH2I(*h) ;
    printf("Set %s \n",fPHOSBadMap[mod]->GetName());
  }
  
private:
  AliAnalysisTaskPi0DiffEfficiency(const AliAnalysisTaskPi0DiffEfficiency&); // not implemented
  AliAnalysisTaskPi0DiffEfficiency& operator=(const AliAnalysisTaskPi0DiffEfficiency&); // not implemented
  Bool_t IsSameCluster(AliAODCaloCluster * c1,AliAODCaloCluster * c2)const ;
  Bool_t IsGoodChannel(const char * det, Int_t mod,Int_t ix, Int_t iz); //Use addisional bad map for PHOS
  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key
  Bool_t TestLambda(Double_t l1,Double_t l2) ;  //Evaluate Dispersion cuts for photons
  void ProcessMC() ;
 
private:
  AliStack * fStack ;              //stack of MC tracks
  TList * fOutputContainer;        //final histogram container
  TList * fPHOSEvents[1][10][11] ; //Containers for events with PHOS photons
  TClonesArray * fPHOSEvent1 ;      //PHOS photons in current event
  TClonesArray * fPHOSEvent2 ;      //PHOS photons in current event
  AliPHOSCalibData *fPHOSCalibData; // PHOS calibration object
  TF1 *fNonLinCorr;          // Non-linearity correction
 
  //Reaction plain for v2
  Float_t fRPfull ; //!Reaction plain calculated with full TPC 
  Float_t fRPA ;    //!Reaction plain calculated with A-side TPC: eta>0.15 
  Float_t fRPC ;    //!Reaction plain calculated with C-side TPC: eta<-0.15
  Float_t fRPFar ;  //!Reaction plain calculated with TPC: eta>0.6 
  Float_t fRPAFar ; //!Reaction plain calculated with A-side TPC: eta>0.6 
  Float_t fRPCFar ; //!Reaction plain calculated with C-side TPC: eta<-0.6

  Float_t fCentrality ; //!Centrality of the currecnt event

  Int_t fCenBin ;       //! Current centrality bin

  TH2I *fPHOSBadMap[6] ;    //Container for PHOS bad channels map

  AliPHOSGeometry  *fPHOSGeo;  //! PHOS geometry
  Int_t fEventCounter;         // number of analyzed events

  ClassDef(AliAnalysisTaskPi0DiffEfficiency, 1); // PHOS analysis task
};

#endif
