#ifndef AliAnalysisTaskSEHFvn_H
#define AliAnalysisTaskSEHFvn_H

/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// AliAnalysisTaskSEHFvn gives the needed tools for the D
// mesons vn analysis
// Authors: Chiara Bianchin, Robert Grajcarek, Giacomo Ortona,
//          Carlos Perez Lara, Francesco Prino, Anastasia Barbano,
//          Fabrizio Grosa, Andrea Festanti
//
//*************************************************************************

/* $Id$ */
#include "THnSparse.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"
#include "AliHFAfterBurner.h"
#include "AliQnCorrectionsQnVector.h"
#include "AliEventCuts.h"
#include "AliAODTrack.h"
#include "AliESDtrackCuts.h"

class TH1F;
class TH2F;
class TH2D;
class AliMultiDimVector;
class AliRDHFCuts;
class TVector2;

class AliAnalysisTaskSEHFvn : public AliAnalysisTaskSE
{

 public:

  enum DecChannel{kDplustoKpipi,kD0toKpi,kDstartoKpipi,kDstoKKpi,kD0toKpiFromDstar}; //more particles can be added
  enum EventPlaneMeth{kTPC,kTPCVZERO,kVZERO,kVZEROA,kVZEROC,kPosTPCVZERO,kNegTPCVZERO}; //Event plane to be calculated in the task
  enum FlowMethod{kEP,kSP,kEvShape}; // Event Plane, Scalar Product or Event Shape Engeneering methods
  enum q2Method{kq2TPC,kq2PosTPC,kq2NegTPC,kq2VZERO,kq2VZEROA,kq2VZEROC}; // q2 for Event Shape to be calculated in the task
  enum EventPlaneDet{kNone=-1,kFullTPC,kPosTPC,kNegTPC,kFullV0,kV0A,kV0C};
  //  enum SubEvents{kFullTPC,kPosTPC,kNegTPC,kSingleV0Side}; //Sub-events for V0 EP

  AliAnalysisTaskSEHFvn();
  AliAnalysisTaskSEHFvn(const char *name, AliRDHFCuts *rdCuts, Int_t decaychannel);

  virtual ~AliAnalysisTaskSEHFvn();

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
  void SetFlowMethod(AliAnalysisTaskSEHFvn::FlowMethod meth){fFlowMethod=meth;}
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
  void SetRecomputeTPCq2(Bool_t opt, Double_t fracKeep=1.1, Int_t removeDau=0, Bool_t removeNdaurandtracks=kFALSE, Bool_t requiremass=kFALSE, Double_t deltaeta=0., Bool_t removesoftpionfromq2=kFALSE){
    fOnTheFlyTPCq2=opt;
    fFractionOfTracksForTPCq2=fracKeep;
    fRemoveDauFromq2=removeDau;
    fRequireMassForDauRemFromq2=requiremass;
    fRemoverSoftPionFromq2=removesoftpionfromq2;
    if(fracKeep<1. && removeNdaurandtracks) {
      AliWarning("AliAnalysisTaskSEHFvn::Impossible to set fFractionOfTracksForTPCq2<1 and fRemoveNdauRandomTracks at the same time! fRemoveNdauRandomTracks setted kFALSE.\n");
      fRemoveNdauRandomTracks=kFALSE;
    }
    else {
      fRemoveNdauRandomTracks=removeNdaurandtracks;
    }
    if((fracKeep<1. || fRemoveNdauRandomTracks || removeDau==2) && deltaeta>0.) {
      AliWarning("AliAnalysisTaskSEHFvn::Impossible to set fDeltaEtaDmesonq2>0 and fRemoveNdauRandomTracks or fFractionOfTracksForTPCq2<1 or fRemoveDauFromq2=2 at the same time! fDeltaEtaDmesonq2 setted to 0.\n");
      fDeltaEtaDmesonq2=0.;
    }
    else {
      fDeltaEtaDmesonq2=deltaeta;
    }
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
    if(limit<1) fScalProdLimit=limit;
    else fScalProdLimit=1;
  }
  void SetOnTheFlyTPCEtaLimits(Double_t etamin, Double_t etamax) {
    fTPCEtaMin=etamin;
    fTPCEtaMax=etamax;
  }
  void SetRemoveDaughtersFromq2(Int_t removeDau, Bool_t requiremass) {fRemoveDauFromq2=removeDau; fRequireMassForDauRemFromq2=requiremass;}
  void SetEnableQnFrameworkCorrForq2(Bool_t usecorr) {fUseQnFrameworkCorrq2=usecorr;}
  void SetEnableEPVsq2VsCentHistos(Bool_t enablehistos=kTRUE) {fEPVsq2VsCent=enablehistos;}
  void SetEnableNtrklVsq2VsCentHistos(Bool_t enablehistos=kTRUE) {fEnableNtrklHistos=enablehistos;}

  void Setq2PercentileSelection(TString splinesfilepath);
  
  //additional event cuts for Pb-Pb 2015 
  void SetUseCentralityMultiplicityCorrelationCut(Bool_t strongcuts=kFALSE) {
    
    fEnableCentralityCorrCuts=kTRUE;
    fEnableCentralityMultiplicityCorrStrongCuts=strongcuts;
    
    fEventCuts.SetupLHC15o();
    fEventCuts.SetManualMode();
    if(strongcuts) {
      fEventCuts.fUseVariablesCorrelationCuts=true;
      fEventCuts.fUseStrongVarCorrelationCut=true;
    }
  }
  
  AliEventCuts& GetAliEventCut() {
    return fEventCuts;
  }
  
  //options for kD0toKpiFromDstar
  void SetOptD0FromDstar(Int_t option=0, Bool_t useFilterBit4softPion=kFALSE) {
    fOptD0FromDstar=option;
    fUseFiltBit4SoftPion=useFilterBit4softPion;
  }
  Bool_t IsSoftPionSelected(AliAODTrack* track);

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void LocalInit();// {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

 private:

  AliAnalysisTaskSEHFvn(const AliAnalysisTaskSEHFvn &source);
  AliAnalysisTaskSEHFvn& operator=(const AliAnalysisTaskSEHFvn& source);

  void CalculateInvMasses(AliAODRecoDecayHF* d,Float_t* &masses,Int_t& nmasses);

  void FillDplus(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin, Float_t dphi, const Float_t* masses,Int_t isSel,Int_t icentr, Double_t phiD, Double_t etaD, Double_t QA[2], Double_t QB[2]);
  void FillD02p(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin, Float_t dphi, const Float_t* masses, Int_t isSel,Int_t icentr, Double_t phiD, Double_t etaD, Double_t QA[2], Double_t QB[2]);
  void FillDstar(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin, Float_t dphi, const Float_t* masses,Int_t isSel,Int_t icentr, Double_t phiD, Double_t etaD, Double_t QA[2], Double_t QB[2]);
  void FillDs(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin, Float_t dphi, const Float_t* masses,Int_t isSel,Int_t icentr, Double_t phiD, Double_t etaD, Double_t QA[2], Double_t QB[2]);
  Float_t GetEventPlane(AliAODEvent* aod, AliEventplane *pl, Double_t eventplaneqncorrTPC[3], Double_t eventplaneqncorrVZERO[3], Float_t &planereso, Float_t &deltaSubAC, Float_t &deltaSubBC, Int_t &nSubEvents);
  Float_t GetEventPlaneForCandidate(AliAODRecoDecayHF* d, AliEventplane *pl);
  Float_t GetEventPlaneForCandidateNewQnFw(AliAODRecoDecayHF* d, const TList *list);
  //  Float_t GetEventPlaneFromV0(AliAODEvent *aodEvent);
  void ComputeTPCEventPlane(AliAODEvent* aod, Double_t &rpangleTPC, Double_t &rpangleTPCpos,Double_t &rpangleTPCneg) const;
  Double_t ComputeTPCq2(AliAODEvent* aod, Double_t &q2TPCfull, Double_t &q2TPCpos, Double_t &q2TPCneg, Double_t q2VecFullTPC[2], Double_t q2VecPosTPC[2], Double_t q2VecNegTPC[2], Double_t multQvecTPC[3], std::vector<Int_t>& labrejtracks) const;
  void CreateSparseForEvShapeAnalysis();
  Double_t Getq2(TList* qnlist, Int_t q2meth, Double_t &mult);
  Bool_t isInMassRange(Double_t massCand, Double_t pt);
  Double_t GetTPCq2DauSubQnFramework(Double_t qVectWOcorr[2], Double_t multQvec, Int_t nDauRemoved, Double_t qVecDau[2], Double_t corrRec[2], Double_t LbTwist[2], Bool_t isTwistApplied);
  void RemoveTracksInDeltaEtaFromOnTheFlyTPCq2(AliAODEvent* aod, Double_t etaD, Double_t etaLims[2], Double_t qVec[2], Double_t &M, std::vector<Int_t> daulab);
  
  TH1F* fHistEvPlaneQncorrTPC[3];   //! histogram for EP
  TH1F* fHistEvPlaneQncorrVZERO[3]; //! histogram for EP
  TH1F* fhEventsInfo;           //! histogram send on output slot 1
  TH1F *fHistCentrality[3];     //!<!hist. for cent distr (all,sel ev,out of cent)
  TH2F *fHistCandVsCent;        //!<!hist. for number of selected candidates vs. cent
  TH2F *fHistCandMassRangeVsCent; //!<!hist. for number of selected candidates in mass range vs. cent
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
  Int_t fRemoveDauFromq2;      // flag to activate removal of daughter tracks from q2 computation (0->not removed, 1->remove single cand, 2->remove all candidates)
  Double_t fTPCEtaMin;          // min eta for the Q-vector computed on the fly with TPC tracks (both EP and q2)
  Double_t fTPCEtaMax;          // max eta for the Q-vector computed on the fly with TPC tracks (both EP and q2)
  Bool_t fRemoveNdauRandomTracks; // flag to activate the removal of nDau random tracks in the q2TPC computed on the fly
  Bool_t fUseQnFrameworkCorrq2; // flag to activate the Qn-framework corrections for the q2
  Bool_t fRequireMassForDauRemFromq2; // flag to activate mass range when removing daughter tracks from q2
  Double_t fDeltaEtaDmesonq2; //eta gap between q2 and D mesons
  Bool_t fEPVsq2VsCent; //flag to enable EP vs. q2 vs. centrality TH3F in case of kEvShape
  Bool_t fEnableNtrklHistos; //flag to enable Ntrklts vs. q2 vs. centrality TH3F in case of kEvShape
  Bool_t fRemoverSoftPionFromq2; //flag to enable also the removal of the soft pions from q2 for D*
  Bool_t fPercentileq2; //flag to replace q2 with its percentile in the histograms
  TList* fq2SplinesList[6]; //lists of splines used to compute the q2 percentile
  Bool_t fEnableCentralityCorrCuts; //enable V0M - CL0 centrality correlation cuts
  Bool_t fEnableCentralityMultiplicityCorrStrongCuts; //enable centrality vs. multiplicity correlation cuts
  AliEventCuts fEventCuts; //Event cut object for centrality correlation event cuts
  Int_t fOptD0FromDstar; //option for D0 from Dstar analysis (0: starting from Dstar, 1: strarting from Dzero)
  Bool_t fUseFiltBit4SoftPion; //flag to enable filterbit 4 for soft pion in D0 from Dstar analysis
  AliESDtrackCuts *fCutsSoftPion; //track cuts for soft pions used in case of D0 from Dstar analysis (option 1)
  
  AliAnalysisTaskSEHFvn::FlowMethod fFlowMethod;

  ClassDef(AliAnalysisTaskSEHFvn,18); // AliAnalysisTaskSE for the HF v2 analysis
};

#endif
