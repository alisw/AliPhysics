#ifndef AliPHOSHijingEfficiency_cxx
#define AliPHOSHijingEfficiency_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

class TObjArray;
class THashList ;
class TH1F;
class TH1D ;
class TH2I;
class TH2F;
class TH3F;
class TF1 ;
class TParticle;
class AliStack ;
class AliESDtrackCuts;
class AliPHOSGeometry;
class AliESDEvent ;
class AliPHOSCalibData;
class AliESDtrack ;
class AliESDCaloCluster ;
class AliCaloPhoton ;

#include "AliAnalysisTaskSE.h"

class AliPHOSHijingEfficiency : public AliAnalysisTaskSE {
public:
  AliPHOSHijingEfficiency(const char *name = "AliPHOSHijingEfficiency");
  virtual ~AliPHOSHijingEfficiency() {}
  
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
  AliPHOSHijingEfficiency(const AliPHOSHijingEfficiency&); // not implemented
  AliPHOSHijingEfficiency& operator=(const AliPHOSHijingEfficiency&); // not implemented
  Bool_t IsGoodChannel(const char * det, Int_t mod,Int_t ix, Int_t iz); //Use addisional bad map for PHOS

  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z, Double_t w) const ; //Fill 3D histogram witn name key
  void FillAllHistograms(const char * key,AliCaloPhoton * ph)const ; //Fill all possible PID combinations for a given photon

  Bool_t TestLambda(Double_t pt,Double_t l1,Double_t l2) ;  //Evaluate Dispersion cuts for photons
  Bool_t TestLambda2(Double_t pt,Double_t l1,Double_t l2) ;  //Evaluate Dispersion cuts for photons
  Double_t TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge);
  Bool_t TestTOF(Double_t t, Double_t e) ;
  Int_t ConvertRunNumber(Int_t run) ; 
  void ProcessMC() ;
  Int_t FindPrimary(AliESDCaloCluster* clu, Bool_t& sure) ;

  void  EvalLambdas(AliESDCaloCluster * clu, Int_t iR,Double_t &m02, Double_t &m20);
  Double_t CoreEnergy(AliESDCaloCluster * clu); 
  void    Reclusterize(AliESDCaloCluster * clu) ;
  Bool_t  AreNeibors(Int_t id1,Int_t id2) ;
  Double_t PrimaryWeight(Int_t primary);
  Double_t PrimaryParticleWeight(TParticle * particle);
  void FillSecondaries() ;
  Int_t FindCommonParent(Int_t iPart, Int_t jPart) ;
  Bool_t HaveParent(Int_t iPart, Int_t pdgParent);
  Bool_t InPi0mass(Double_t m, Double_t pt);
  
private:
  AliStack * fStack ;             //
  THashList * fOutputContainer;   //final histogram container
  TList * fPHOSEvents[1][10][11] ; //Containers for events with PHOS photons
  TClonesArray * fPHOSEvent ;      //PHOS photons in current event
  AliPHOSCalibData *fPHOSCalibData; // PHOS calibration object
 
  //Reaction plain for v2
  Float_t fRP ;       //Reaction plane calculated with full TPC 
  Float_t fRPV0A ;    //Reaction plain calculated with A-side TPC: eta>0.15 
  Float_t fRPV0C ;    //Reaction plain calculated with C-side TPC: eta<-0.15
  Bool_t fHaveTPCRP ; // Is TPC RP defined?
  
  Int_t fRunNumber ;    //Current run number
  Float_t fCentrality ; //Centrality of the currecnt event

  Int_t fCenBin ;       // Current centrality bin

  TH2I *fPHOSBadMap[6] ;    //Container for PHOS bad channels map

  AliPHOSGeometry  *fPHOSGeo;  // PHOS geometry
  Int_t fEventCounter;         // number of analyzed events

  ClassDef(AliPHOSHijingEfficiency, 1); // PHOS analysis task
};

#endif
