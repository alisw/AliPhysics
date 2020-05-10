#ifndef ALIANALYSISTASKTAGGEDPHOTONSLOCAL_H
#define ALIANALYSISTASKTAGGEDPHOTONSLOCAL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Analysis for PHOS Tagged Photons 
// marks photons making pi0 with any other photon
// and calculates necessary corrections for fake pairs and
// decay partners escaped acceptance. If MC info is present 
// fills set of controll histograms.
//
//*-- Dmitry Peresunko
//////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskSE.h"  
class AliAnalysisUtils ;
class AliAODEvent ; 
class THashList ; 
class TH2I ;
class AliPHOSGeometry;
class AliCaloPhoton;
class AliAODMCParticle ;
class AliVCluster ;
class AliTriggerAnalysis ;
class TParticle ;
class AliPHOSTriggerUtils ;

class AliAnalysisTaskTaggedPhotons : public AliAnalysisTaskSE {

public:
  
  enum mcType{kFullMC, kSingleGamma, kSinglePi0, kSingleEta } ; 
  enum cutType{kDefCut, kLowECut} ;
  enum phosTriggerType{kPHOSAny,kPHOSL0,kPHOSL1low,kPHOSL1med,kPHOSL1high} ;
  enum trackSelections{kLHC13x,kFAST,kCENTwoSSD,kCENTwSSD} ;
    
  AliAnalysisTaskTaggedPhotons() ;
  AliAnalysisTaskTaggedPhotons(const char *name) ;
  AliAnalysisTaskTaggedPhotons(const AliAnalysisTaskTaggedPhotons& ap) ;   
  AliAnalysisTaskTaggedPhotons& operator = (const AliAnalysisTaskTaggedPhotons& ap) ;
  virtual ~AliAnalysisTaskTaggedPhotons() ;
   
  virtual void UserCreateOutputObjects(); 
  virtual void Init() ; 
  virtual void LocalInit() { Init() ; }
  virtual void UserExec(Option_t * opt = "") ;
  virtual void Terminate(Option_t * opt = "") ;

  void SetTrigger(Bool_t isPHOSTrig){fIsMB=isPHOSTrig;}   //Analyse MinBias of PHOS triggered events
  void SetPHOSTrigger(phosTriggerType t=kPHOSL0){fPHOSTrigger=t;}
  void UseCaloFast(Bool_t use=kTRUE){fUseCaloFastTr=use ;} //Use PHOS trigger in CaloFast cluster
  void SetMC(Bool_t isMC=kTRUE){fIsMC=isMC;}              //Is is MC or real data
  void SetMCType(mcType a){ fIsMC = kTRUE; fMCType = a; } 
  void SetFastMC(void){fIsFastMC=kTRUE;fIsMC=kTRUE; }     //Only for MC, bypass event quality checks (e.g. for single pi0 MC)
  void SetRunNumber(int run){fRunNumber=run;fForseRun=kTRUE; }  //Use given run number, don't read from data
  void SetPi0WeightParameters(TArrayD * ar) ;             //Introduce weight for primary pi0
  void SetDistanceToBad(Float_t cut=2.5){fMinBCDistance=cut;}     //Distance to bad module
  void SetCluCutType(cutType cut = kDefCut){fCutType = cut ; }      //To use default or lowE (Daiki's) cut
  void SetTimeCut(Float_t cut=25.e-9){fTimeCut=cut;}              //Time selection, use only for data
  void SetCentralityEstimator(Int_t est=1){fCentEstimator=est;}   //Centrality estimator, pPb: 1: V0A/C, 2: V0M, 3: ZNA/C,  4: CL1
  void SetCentralityWeights(TString filename="MBCentralityWeights.root") ;  //for pp: 
  void SetMultiplicityBins(TArrayI *ar){fNCenBin=ar->GetSize() ; fCenBinEdges.Set(ar->GetSize(),ar->GetArray());} //for pp: 
  void SetNonLinearity(Double_t a=1., Double_t b=0., Double_t c=1){ fNonlinA=a; fNonlinB=b; fNonlinC=c;}
  void SetTrackSelection(trackSelections s=kCENTwSSD){fTrackSelection=s;}
protected:
  void    FillMCHistos() ;
  void    FillTaggingHistos() ;
  Int_t   GetFiducialArea(const Float_t * pos)const ;                           //what kind of fiducial area hit the photon
  Int_t   IsSameParent(const AliCaloPhoton *p1, const AliCaloPhoton *p2) const; //Check MC genealogy; return PDG of parent
  Bool_t  IsGoodChannel(Int_t mod, Int_t ix, Int_t iz) ;
  Double_t  InPi0Band(Double_t m, Double_t pt)const; //Check if invariant mass is within pi0 peak
  Bool_t  TestDisp(Double_t l0, Double_t l1, Double_t e)const  ;
  Bool_t  TestTOF(Double_t /*t*/,Double_t /*en*/)const{return kTRUE;} 
  Bool_t  TestCharged(Double_t dr,Double_t en)const ;
  void    InitGeometry() ;  //read geometry from TENDER, if not, get rotation matrixes from AODB
  Int_t   EvalIsolation(TLorentzVector * ph,Bool_t isPhoton) ;
  Bool_t  TestLambda(Double_t pt,Double_t l1,Double_t l2) ;
  Bool_t  TestPID(Int_t iPID, AliCaloPhoton* part) ;
  Bool_t  TestPID(Int_t iPID, AliCaloPhoton* p1, AliCaloPhoton* p2) ;
  Double_t PrimaryParticleWeight(AliAODMCParticle * particle) ;
  Int_t   FindPrimary(AliVCluster*, Bool_t&);
  Double_t TOFCutEff(Double_t x );
  Double_t NonLinearity(Double_t e){ return e*fNonlinA*(1.- fNonlinB*TMath::Exp(-e*fNonlinC)); } 
  
  Bool_t   SelectCentrality(AliVEvent * event) ;
  Double_t CalculateSphericity() ;
  Double_t CalculateSpherocity() ;
  
  Double_t TrigCentralityWeight(Double_t x); //Correction for PHOS trigger centrality bias
  Double_t MBCentralityWeight(Double_t x);   //Correction for Pileup cut centrality bias
  
  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key
  void FillPIDHistograms(const char * name,  AliCaloPhoton * p) const ;
  void FillPIDHistograms(const char * name,  AliCaloPhoton * p ,Double_t y) const ;
  void FillPIDHistograms(const char * name,  AliCaloPhoton * p ,  AliCaloPhoton * p2,Double_t y, Bool_t isReal) const ;

private:

  AliPHOSGeometry  *fPHOSgeom;      //!PHOS geometry
  THashList *   fOutputContainer ;  //!List of output histograms
  TClonesArray *fStack ;            //!Pointer to MC stack (AOD)
  TClonesArray * fTrackEvent ;      //!List of tracks in the event
  TClonesArray * fPHOSEvent ;       //!List of tracks in the event
  TList   * fPHOSEvents[20][20] ;    //!Previous events for mixing
  TList   * fCurrentMixedList;      //! list of previous evetns for given centrality
  AliTriggerAnalysis * fTriggerAnalysis ; //!
  AliAnalysisUtils * fUtils ;       //!
  AliPHOSTriggerUtils * fPHOSTrigUtils ; //! utils to analyze PHOS trigger
 
  Int_t   fCentEstimator;       //Centrality estimator: 1: V0A/C, 2: V0M, 3: ZNA/C,  4: CL1
  Int_t   fNCenBin ;            //NUmber of centrality bins
  TArrayI fCenBinEdges;         //Centrality binning
  Double_t fCentrality;   //!
  Double_t fCentWeight ;  //! Weight to correct bias in PHOS triigeered events
  TH1F * fCentralityWeights[6]; //!Weights to correct centrality non-flatness
  Int_t  fCentBin ;       //! current centrality bin
  Int_t  fRunNumber ;     //! current run number
  Bool_t fForseRun ;      // use fixed run number, dont read from data
  Bool_t fIsMB ;          //which trigger to use
  Bool_t fUseCaloFastTr ; //use also 
  Bool_t fIsMC ;          //Is this is MC
  Bool_t fIsFastMC;       //This is fast MC, bypass event checks
  Double_t fRP;           //! Reaction plane orientation
  
  //Fiducial area parameters
  Float_t fZmax ;               //Rectangular
  Float_t fZmin ;               //area
  Float_t fPhimax ;             //covered by
  Float_t fPhimin ;             //full calorimeter
  Float_t fMinBCDistance;       //minimal distance to bad channel
  Float_t fTimeCut ;            //Time cut
  Double_t fWeightParamPi0[7] ; //!Parameters to calculate weights in MC
  Double_t fNonlinA;
  Double_t fNonlinB;
  Double_t fNonlinC;
  Int_t   fNPID ;               // Number of PID cuts
  mcType  fMCType ;             // Type of MC production: full, single g,pi0,eta,
  cutType fCutType;             // Type of cluster cuts used in analysis
  phosTriggerType fPHOSTrigger; // Kind of PHOS trigger: L0,L1
  trackSelections fTrackSelection; 
  
  //
  TH2I * fPHOSBadMap[6] ;        //! 
  TH2F * fhReMod[5];             //! Real per module
  TH2F * fhMiMod[5];             //! Mixed per module
  TH2F * fhRe[3][10][8];         //! Real: (Emin cut, Centrality, PID cut)
  TH2F * fhMi[3][10][8];         //! Mixed: (Emin cut, Centrality, PID cut) 
  TH2F * fhReSingle[3][10][8];   //!
  TH2F * fhMiSingle[3][10][8];    //!
  TH2F * fhReSingleIso[3][10][8]; //!
  TH2F * fhMiSingleIso[3][10][8]; //!
  TH1F * fhPiIsolation[20][10] ;  //!
  TH2F * fhReTruePi0[3][10][8];   //! Real, true pi0: (Emin cut, Centrality, PID cut)
  TH2F * fhReTrueEta[3][10][8];   //! Real, true eta: (Emin cut, Centrality, PID cut)
  
  TH2F * fhQAAllEpartn ;    //!
  TH2F * fhQAAllzpartn ;    //!
  TH2F * fhQAAllxpartn ;    //!
  TH2F * fhQAIsoEpartn ;    //!
  TH2F * fhQAIsozpartn ;    //!
  TH2F * fhQAIsoxpartn ;    //!
  TH2F * fhQAAllEpartnBg ;  //!
  TH2F * fhQAAllzpartnBg ;  //!
  TH2F * fhQAAllxpartnBg ;  //!
  TH2F * fhQAIsoEpartnBg ;  //!
  TH2F * fhQAIsozpartnBg ;  //!
  TH2F * fhQAIsoxpartnBg ;  //!
      
  ClassDef(AliAnalysisTaskTaggedPhotons, 6);   // a PHOS photon analysis task 
};
#endif // ALIANALYSISTASKTAGGEDPHOTONSLOCAL_H
