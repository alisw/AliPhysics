#ifndef ALIANALYSISTASKSEHFV2_H
#define ALIANALYSISTASKSEHFV2_H

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

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"
#include "AliHFAfterBurner.h"


class TH1F;
class TH2D;
class AliMultiDimVector;
class AliRDHFCuts;
class TVector2;

class AliAnalysisTaskSEHFv2 : public AliAnalysisTaskSE
{

 public:

  enum DecChannel{kDplustoKpipi,kD0toKpi,kDstartoKpipi}; //more particles can be added

  AliAnalysisTaskSEHFv2();
  AliAnalysisTaskSEHFv2(const char *name, AliRDHFCuts *rdCuts, Int_t decaychannel,Int_t nbinsphi, Float_t *phibinlimits);
 
  virtual ~AliAnalysisTaskSEHFv2();

  void SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}
  void SetMassLimits(Float_t range,Int_t pdg);
  void SetMassLimits(Float_t lowlimit, Float_t uplimit);
  void SetNMassBins(Int_t nbins){fNMassBins=nbins;}
  //void SetUseV0EP(Bool_t flagV0EP){fUseV0EP=flagV0EP;}
  void SetV0EventPlaneOrder(Int_t n){fV0EPorder=n;}
  void SetMinCentrality(Int_t mincentr){fMinCentr=mincentr;}
  void SetMaxCentrality(Int_t maxcentr){fMaxCentr=maxcentr;}
  void SetUseAfterBurner(Bool_t ab){fUseAfterBurner=ab;}
  void SetAfterBurner(AliHFAfterBurner *ab){fAfterBurner=ab;}
  void SetEtaGapFeatureForEventplaneFromTracks (Bool_t etaGap) {fEtaGap = etaGap;}
  void SetVZEROParHist(TH2D** h);

  Float_t GetUpperMassLimit()const {return fUpmasslimit;}
  Float_t GetLowerMassLimit()const {return fLowmasslimit;}
  Int_t GetNMassBins()const {return fNMassBins;}
  Int_t GetPhiBin(Float_t deltaphi);
  //Float_t GetPhi02Pi(Float_t phi);
  Float_t GetPhi0Pi(Float_t phi);
  AliHFAfterBurner *GetAfterBurner()const {return fAfterBurner;}
  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void LocalInit();// {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
    
 private:

  AliAnalysisTaskSEHFv2(const AliAnalysisTaskSEHFv2 &source);
  AliAnalysisTaskSEHFv2& operator=(const AliAnalysisTaskSEHFv2& source); 

  void CalculateInvMasses(AliAODRecoDecayHF* d,Float_t* &masses,Int_t& nmasses);

  void FillDplus(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin, Float_t dphi,Float_t* masses,Int_t isSel,Int_t icentr);
  void FillD02p(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin, Float_t dphi,Float_t* masses, Int_t isSel,Int_t icentr);
  void FillDstar(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin, Float_t dphi,Float_t* masses,Int_t isSel,Int_t icentr);
  Float_t GetEventPlaneForCandidate(AliAODRecoDecayHF* d, TVector2* q,AliEventplane *pl,TVector2* qsub1,TVector2* qsub2);
  Float_t GetEventPlaneFromV0(AliAODEvent *aodEvent);


  TH1F* fhEventsInfo;           //! histogram send on output slot 1
  TList   *fOutput;             //! list send on output slot 2
  AliRDHFCuts *fRDCuts;         //cut values (saved in slot 3)
  TList *fParHist;               //list for VZERO EP parameters (slot 4)
  TH2D *fHistvzero[6];            //histograms for VZERO EP parameters
  Float_t fLowmasslimit;        //lower inv mass limit for histos
  Float_t fUpmasslimit;         //upper inv mass limit for histos
  Int_t fNPtBins;               //number of pt bins
  Int_t fNPhiBinLims;           //number of delta phi bins limits (= number of bins +1)
  Float_t *fPhiBins;            //[fNPhiBinLims] limits of each phi bin
  Int_t fNMassBins;             //number of bins in the mass histograms
  Bool_t fReadMC;               //flag for access to MC
  Bool_t fUseAfterBurner;      //enable afterburning
  Int_t fDecChannel;            //decay channel identifier
  AliHFAfterBurner *fAfterBurner;//Afterburner options
  Bool_t fUseV0EP;              //flag to select EP method
  Int_t  fV0EPorder;            //harmonic for VZERO event plane
  Int_t fMinCentr;              //minimum centrality
  Int_t fMaxCentr;              //maximum centrality
  Bool_t fEtaGap;               // Eta gap feature for Eventplane from tracks; be careful that you do the correct settings in AddTaskEventPlane.C !!!!

  ClassDef(AliAnalysisTaskSEHFv2,1); // AliAnalysisTaskSE for the HF v2 analysis
};

#endif
