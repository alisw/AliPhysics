#ifndef AliAnalysisTaskSEHFv2_H
#define AliAnalysisTaskSEHFv2_H

/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// AliAnalysisTaskSEHFv2 gives the needed tools for the D
// mesons v2 analysis
// Authors: Chiara Bianchin, cbianchi@pd.infn.it,
//          Robert Grajcarek, grajcarek@physi.uni-heidelberg.de
//          Giacomo Ortona, ortona@to.infn.it,
//          Carlos Perez Lara, carlos.eugenio.perez.lara@cern.ch
//          Francesco Prino, prino@to.infn.it
//
//*************************************************************************

/* $Id$ */
#include "THnSparse.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"
#include "AliHFAfterBurner.h"
#include "AliQnCorrectionsQnVector.h"

class TH1F;
class TH2F;
class TH2D;
class AliMultiDimVector;
class AliRDHFCuts;
class TVector2;

class AliAnalysisTaskSEHFv2 : public AliAnalysisTaskSE
{
    
 public:
    
  enum DecChannel{kDplustoKpipi,kD0toKpi,kDstartoKpipi,kDstoKKpi}; //more particles can be added
  enum EventPlaneMeth{kTPC,kTPCVZERO,kVZERO,kVZEROA,kVZEROC}; //Event plane to be calculated in the task
  enum FlowMethod{kEP,kSP,kEvShape}; // Event Plane, Scalar Product or Event Shape Engeneering methods
  enum q2Method{kq2TPC,kq2PosTPC,kq2NegTPC,kq2VZERO,kq2VZEROA,kq2VZEROC}; // q2 for Event Shape to be calculated in the task
  enum EventPlaneDet{kNone=-1,kFullTPC,kPosTPC,kNegTPC,kFullV0,kV0A,kV0C};
  //  enum SubEvents{kFullTPC,kPosTPC,kNegTPC,kSingleV0Side}; //Sub-events for V0 EP
    
  AliAnalysisTaskSEHFv2();
  AliAnalysisTaskSEHFv2(const char *name, AliRDHFCuts *rdCuts, Int_t decaychannel);
    
  virtual ~AliAnalysisTaskSEHFv2();

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
  void SetV0EventPlaneOrder(Int_t n){fV0EPorder=n;}
  void SetMinCentrality(Int_t mincentr){fMinCentr=mincentr;}
  void SetMaxCentrality(Int_t maxcentr){fMaxCentr=maxcentr;}
  void SetUseAfterBurner(Bool_t ab){fUseAfterBurner=ab;}
  void SetAfterBurner(AliHFAfterBurner *ab){fAfterBurner=ab;}
  void SetEtaGapFeatureForEventplaneFromTracks (Bool_t etaGap) {fEtaGap = etaGap;}
  void SetCentralityBinWidthPerMil(Int_t w){fCentBinSizePerMil=w;}
  void SetFlowMethod(AliAnalysisTaskSEHFv2::FlowMethod meth){fFlowMethod=meth;}
  void SetNormMethod(TString normmethod="QoverQlength") {fNormMethod=normmethod;}

  void SetEventPlaneMethod(Int_t epmethod);
  void SetTPCEPOnly(){SetEventPlaneMethod(kTPC);}
  void SetVZEROEP(){SetEventPlaneMethod(kVZERO);}
  void SetVZEROAEP(){SetEventPlaneMethod(kVZEROA);}
  void SetVZEROCEP(){SetEventPlaneMethod(kVZEROC);}
  void SetTPCEP(){SetEventPlaneMethod(kTPCVZERO);}
  void SetEventPlanesCompatibility(Float_t comp) {fEventPlanesComp=comp;}
  void SetUseNewQnCorrFw(Bool_t flag) {fUseNewQnCorrFw=flag;}
  
  Float_t GetEventPlanesCompatibility()const {return fEventPlanesComp;}
  Float_t GetUpperMassLimit()const {return fUpmasslimit;}
  Float_t GetLowerMassLimit()const {return fLowmasslimit;}
  Int_t GetNMassBins()const {return fNMassBins;}
  //Float_t GetPhi02Pi(Float_t phi);
  Float_t GetPhi0Pi(Float_t phi);
  AliHFAfterBurner *GetAfterBurner()const {return fAfterBurner;}
  const AliQnCorrectionsQnVector *GetQnVectorFromList( const TList *list,
						       const char *subdetector,
						       const char *expectedstep,
						       const char *altstep);

  void SetSeparateD0D0bar(Bool_t separate) {fSeparateD0D0bar=separate;}
  void Setq2Method(Int_t q2method) {fq2Meth=q2method;}
  
  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void LocalInit();// {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
    
 private:
    
  AliAnalysisTaskSEHFv2(const AliAnalysisTaskSEHFv2 &source);
  AliAnalysisTaskSEHFv2& operator=(const AliAnalysisTaskSEHFv2& source);
    
  void CalculateInvMasses(AliAODRecoDecayHF* d,Float_t* &masses,Int_t& nmasses);
    
  void FillDplus(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin, Float_t dphi, const Float_t* masses,Int_t isSel,Int_t icentr, Double_t phiD, Float_t dphi2);
  void FillD02p(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin, Float_t dphi, const Float_t* masses, Int_t isSel,Int_t icentr, Double_t phiD, Float_t dphi2);
  void FillDstar(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin, Float_t dphi, const Float_t* masses,Int_t isSel,Int_t icentr, Double_t phiD, Float_t dphi2);
  void FillDs(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin, Float_t dphi, const Float_t* masses,Int_t isSel,Int_t icentr, Double_t phiD, Float_t dphi2);
  Float_t GetEventPlane(AliAODEvent* aod, AliEventplane *pl, Double_t eventplaneqncorrTPC[3], Double_t eventplaneqncorrVZERO[3], Double_t q2);
  Float_t GetEventPlaneForCandidate(AliAODRecoDecayHF* d, AliEventplane *pl);
  Float_t GetEventPlaneForCandidateNewQnFw(AliAODRecoDecayHF* d, const TList *list);
  //  Float_t GetEventPlaneFromV0(AliAODEvent *aodEvent);
  
  void CreateSparseForEvShapeAnalysis();
  Double_t Getq2(TList* qnlist);
  
  TH1F* fHistEvPlaneQncorrTPC[3];   //! histogram for EP
  TH1F* fHistEvPlaneQncorrVZERO[3]; //! histogram for EP
  TH1F* fhEventsInfo;           //! histogram send on output slot 1
  TH1F *fHistCentrality[3];     //!<!hist. for cent distr (all,sel ev,out of cent)
  TList   *fOutput;             //! list send on output slot 2
  AliRDHFCuts *fRDCuts;         //cut values (saved in slot 3)
  Float_t fLowmasslimit;        //lower inv mass limit for histos
  Float_t fUpmasslimit;         //upper inv mass limit for histos
  Int_t fNPtBins;               //number of pt bins
  Int_t fNMassBins;             //number of bins in the mass histograms
  Bool_t fReadMC;               //flag for access to MC
  Bool_t fUseAfterBurner;      //enable afterburning
  Int_t fDecChannel;            //decay channel identifier
  AliHFAfterBurner *fAfterBurner;//Afterburner options
  Float_t fEventPlanesComp;     // Maximum distance between TPC/VZERO event planes
  Int_t  fV0EPorder;            //harmonic for VZERO event plane
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

  AliAnalysisTaskSEHFv2::FlowMethod fFlowMethod;
    
  ClassDef(AliAnalysisTaskSEHFv2,7); // AliAnalysisTaskSE for the HF v2 analysis
};

#endif
