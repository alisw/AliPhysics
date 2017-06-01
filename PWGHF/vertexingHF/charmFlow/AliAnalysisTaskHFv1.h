#ifndef AliAnalysisTaskHFv1_H
#define AliAnalysisTaskHFv1_H

/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//**************************************************************************
// AliAnalysisTaskHFv1 directed flow of D mesons with scalar
// product method (modified from AliAnalysisTaskHFv1)
// Authors: Andrea Dubla, Jacopo Margutti
//
// AliAnalysisTaskHFv1 gives the needed tools for the D
// mesons vn analysis with event plane method
// Authors: Chiara Bianchin, Robert Grajcarek, Giacomo Ortona,
//          Carlos Perez Lara, Francesco Prino, Anastasia Barbano,
//          Fabrizio Grosa, Andrea Festanti
//**************************************************************************

/* $Id$ */
#include "THnSparse.h"
#include "AliAnalysisTaskSE.h"
#include "AliFlowVector.h"
#include "AliAnalysisVertexingHF.h"
#include "AliHFAfterBurner.h"
#include "AliQnCorrectionsQnVector.h"
#include "AliAnalysisTaskZDCEP.h"

class TH1F;
class TH2F;
class TH2D;
class AliFlowVector;
class AliMultiDimVector;
class AliRDHFCuts;
class TVector2;
class AliAnalysisTaskZDCEP;

class AliAnalysisTaskHFv1 : public AliAnalysisTaskSE
{
    
 public:
    
  enum DecChannel{kDplustoKpipi,kD0toKpi,kDstartoKpipi,kDstoKKpi}; //more particles can be added
  enum EventPlaneMeth{kTPC,kTPCVZERO,kVZERO,kVZEROA,kVZEROC,kPosTPCVZERO,kNegTPCVZERO,kZDC}; //Event plane to be calculated in the task
  enum FlowMethod{kEP,kSP,kEvShape}; // Event Plane, Scalar Product or Event Shape Engeneering methods
  enum q2Method{kq2TPC,kq2PosTPC,kq2NegTPC,kq2VZERO,kq2VZEROA,kq2VZEROC}; // q2 for Event Shape to be calculated in the task
  enum EventPlaneDet{kNone=-1,kFullTPC,kPosTPC,kNegTPC,kFullV0,kV0A,kV0C,kZDCA,kZDCC};
  //  enum SubEvents{kFullTPC,kPosTPC,kNegTPC,kSingleV0Side}; //Sub-events for V0 EP
    
  AliAnalysisTaskHFv1();
  AliAnalysisTaskHFv1(const char *name, AliRDHFCuts *rdCuts, Int_t decaychannel);
    
  virtual ~AliAnalysisTaskHFv1();

  void SetEventPlaneDetector(Int_t det){
    fEvPlaneDet=det;
  }
  void SetSubEventDetectors(Int_t detsubA, Int_t detsubB){
    fSubEvDetA=detsubA; fSubEvDetB=detsubB;
  }

  void SetAODMismatchProtection(Int_t opt=1) {fAODProtection=opt;}
  void SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}
  void SetMassLimits(Float_t range,Int_t pdg);
  void SetMassLimits(Float_t lowlimit, Float_t uplimit);
  void SetNMassBins(Int_t nbins){fNMassBins=nbins;}
  void SetMinCentrality(Int_t mincentr){fMinCentr=mincentr;}
  void SetMaxCentrality(Int_t maxcentr){fMaxCentr=maxcentr;}
  void SetUseAfterBurner(Bool_t ab){fUseAfterBurner=ab;}
  void SetAfterBurner(AliHFAfterBurner *ab){fAfterBurner=ab;}
  void SetEtaGapFeatureForEventplaneFromTracks (Bool_t etaGap) {fEtaGap = etaGap;}
  void SetCentralityBinWidthPerMil(Int_t w){fCentBinSizePerMil=w;}
  void SetFlowMethod(AliAnalysisTaskHFv1::FlowMethod meth){fFlowMethod=meth;}
  void SetNormMethod(TString normmethod="QoverQlength") {fNormMethod=normmethod;}

  void SetHarmonic(Int_t n){fHarmonic=n;}
  void SetEventPlaneMethod(Int_t epmethod);
  void SetTPCEPOnly(){SetEventPlaneMethod(kTPC);}
  void SetVZEROEP(){SetEventPlaneMethod(kVZERO);}
  void SetVZEROAEP(){SetEventPlaneMethod(kVZEROA);}
  void SetVZEROCEP(){SetEventPlaneMethod(kVZEROC);}
  void SetTPCEP(){SetEventPlaneMethod(kTPCVZERO);}
  void SetEventPlanesCompatibility(Float_t comp) {fEventPlanesComp=comp;}
  void SetUseNewQnCorrFw(Bool_t flag) {fUseNewQnCorrFw=flag;}
  void SetRecomputeTPCEventPlane(Bool_t opt, Bool_t usePtWei, Double_t etagap=-1.){
    fOnTheFlyTPCEP=opt;
    fUsePtWeights=usePtWei;
    fEtaGapInTPCHalves=etagap;
  }
  void SetRecomputeTPCq2(Bool_t opt, Double_t fracKeep=1.1){
    fOnTheFlyTPCq2=opt;
    fFractionOfTracksForTPCq2=fracKeep;
  }
  Float_t GetEventPlanesCompatibility()const {return fEventPlanesComp;}
  Float_t GetUpperMassLimit()const {return fUpmasslimit;}
  Float_t GetLowerMassLimit()const {return fLowmasslimit;}
  Int_t GetNMassBins()const {return fNMassBins;}
  //Float_t GetPhi02Pi(Float_t phi);
  Float_t GetPhiInRange(Float_t phi);
  AliHFAfterBurner *GetAfterBurner()const {return fAfterBurner;}
  const AliQnCorrectionsQnVector *GetQnVectorFromList( const TList *list,
						       const char *subdetector,
						       const char *expectedstep,
						       const char *altstep);

  void SetSeparateD0D0bar(Bool_t separate) {fSeparateD0D0bar=separate;}
  void Setq2Method(Int_t q2method) {fq2Meth=q2method;}
  void Setq2Smearing(TString smearingfilepath, TString histoname, Int_t smearingaxis);
  void SetScalProdLimit(Double_t limit=0.3) {
    fScalProdLimit=limit;
  }
  
  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void LocalInit();// {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
    
 private:
    
  AliAnalysisTaskHFv1(const AliAnalysisTaskHFv1 &source);
  AliAnalysisTaskHFv1& operator=(const AliAnalysisTaskHFv1& source);
    
  void CalculateInvMasses(AliAODRecoDecayHF* d,Float_t* &masses,Int_t& nmasses);
    
  void FillDplus(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin, Float_t dphi, const Float_t* masses,Int_t isSel,Int_t icentr, Double_t phiD, Double_t etaD, Double_t ptD, Double_t QA[2], Double_t QB[2]);
  void FillD02p(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin, Float_t dphi, const Float_t* masses, Int_t isSel,Int_t icentr, Double_t phiD, Double_t etaD, Double_t ptD, Double_t QA[2], Double_t QB[2]);
  void FillDstar(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin, Float_t dphi, const Float_t* masses,Int_t isSel,Int_t icentr, Double_t phiD, Double_t etaD, Double_t ptD, Double_t QA[2], Double_t QB[2]);
  void FillDs(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin, Float_t dphi, const Float_t* masses,Int_t isSel,Int_t icentr, Double_t phiD, Double_t etaD, Double_t ptD, Double_t QA[2], Double_t QB[2]);
  Float_t GetEventPlane(AliAODEvent* aod, AliEventplane *pl, Double_t eventplaneqncorrTPC[3], Double_t eventplaneqncorrVZERO[3], Double_t q2);
  Float_t GetEventPlaneForCandidate(AliAODRecoDecayHF* d, AliEventplane *pl);
  Float_t GetEventPlaneForCandidateNewQnFw(AliAODRecoDecayHF* d, const TList *list);
  //  Float_t GetEventPlaneFromV0(AliAODEvent *aodEvent);
  void ComputeTPCEventPlane(AliAODEvent* aod, Double_t &rpangleTPC, Double_t &rpangleTPCpos,Double_t &rpangleTPCneg) const;
  Double_t ComputeTPCq2(AliAODEvent* aod, Double_t &q2TPCfull, Double_t &q2TPCpos,Double_t &q2TPCneg) const;
  void CreateSparseForEvShapeAnalysis();
  Double_t Getq2(TList* qnlist, Int_t q2meth);
  
  TH1F* fHistEvPlaneQncorrTPC[3];   //! histogram for EP
  TH1F* fHistEvPlaneQncorrVZERO[3]; //! histogram for EP
  TH1F* fhEventsInfo;           //! histogram send on output slot 1
  TH1F *fHistCentrality[3];     //!<!hist. for cent distr (all,sel ev,out of cent)
  TList   *fOutput;             //! list send on output slot 2
  AliRDHFCuts *fRDCuts;         //cut values (saved in slot 3)
  Float_t fLowmasslimit;        //lower inv mass limit for histos
  Float_t fUpmasslimit;         //upper inv mass limit for histos
  Float_t fLowEtaLimit;         //lower limit in eta
  Float_t fUpEtaLimit;          //upper limit in eta
  Int_t fNPtBins;               //number of pt bins
  Int_t fNEtaBins;              //number of eta bins
  Int_t fNMassBins;             //number of bins in the mass histograms
  Bool_t fReadMC;               //flag for access to MC
  Bool_t fUseAfterBurner;      //enable afterburning
  Int_t fDecChannel;            //decay channel identifier
  AliHFAfterBurner *fAfterBurner;//Afterburner options
  Float_t fEventPlanesComp;     // Maximum distance between TPC/VZERO event planes
  Int_t fHarmonic;              //harmonic
  Int_t fMinCentr;              //minimum centrality
  Int_t fMaxCentr;              //maximum centrality
  Bool_t fEtaGap;               // Eta gap feature for Eventplane from tracks; be careful that you do the correct settings in AddTaskEventPlane.C !!!!
  Int_t fEvPlaneDet;            // detector for event plane
  Int_t fSubEvDetA;             // detector for 1st subevent
  Int_t fSubEvDetB;             // detector for 2nd subevent
  Int_t fCentBinSizePerMil;     // width of centrality bins
  Int_t fAODProtection;         /// flag to activate protection against AOD-dAOD mismatch.
  /// -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names
  Bool_t fUseNewQnCorrFw;       //flag to use the new Qn correction framework
  TString fCentrBinName;        //centrality bin string
  TString fDetTPCConfName[3];
  TString fDetV0ConfName[3];
  TString fNormMethod;
  THnSparseF* fHistMassPtPhiq2Centr; //THnSparse for the analysis of v2 as a function of q2
  Int_t fq2Meth;                //flag to select q2 method
  Bool_t fSeparateD0D0bar;      //flag to activate the separation of D0 from D0bar in the THnSparse
  Bool_t fOnTheFlyTPCEP;        // flag to compute the TPC EP in the task
  Double_t fEtaGapInTPCHalves;  // eta gap between two halves of TPC (only if fOnTheFlyTPCEP)
  Bool_t fUsePtWeights;         // use pt weights for TPC EP if fOnTheFlyTPCEP is activated
  Bool_t  fOnTheFlyTPCq2; /// flag to compute the TPC q2 in the task
  Double_t  fFractionOfTracksForTPCq2; /// downscaling factor for tracks used in q2
  TH2F* fq2SmearingHisto;       //-> 2D histo for q2smearing
  Bool_t fq2Smearing;           // flag to activate q2 smearing
  Int_t fq2SmearingAxis;        // axis of the smearing histogram corresponding to the q2 used for the analysis
  Double_t fScalProdLimit;      // max value for the scalar product histograms
  
  AliAnalysisTaskHFv1::FlowMethod fFlowMethod;
    
  ClassDef(AliAnalysisTaskHFv1,4); // AliAnalysisTaskSE for the HF v2 analysis
};

#endif
