#ifndef ALIRDHFCUTSXICPLUSTOXIPIPIFROMAODTRACKS_H
#define ALIRDHFCUTSXICPLUSTOXIPIPIFROMAODTRACKS_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//***********************************************************
// Class AliRDHFCutsXictoPLUSXiPiPifromAODtracks
// class for cuts on AOD reconstructed Xic-> pi Xi pi
//***********************************************************

#include "AliRDHFCuts.h"

class AliRDHFCutsXicPlustoXiPiPifromAODtracks : public AliRDHFCuts
{
 public:

  enum EPIDStrategy{
    kNSigmaCuts,
    kCombinedCuts
  };

  AliRDHFCutsXicPlustoXiPiPifromAODtracks(const char* name="CutsXicPlustoXiPiPi");
  virtual ~AliRDHFCutsXicPlustoXiPiPifromAODtracks();
  AliRDHFCutsXicPlustoXiPiPifromAODtracks(const AliRDHFCutsXicPlustoXiPiPifromAODtracks& source);
  AliRDHFCutsXicPlustoXiPiPifromAODtracks& operator=(const AliRDHFCutsXicPlustoXiPiPifromAODtracks& source);

  using AliRDHFCuts::GetCutVarsForOpt;
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters);

  using AliRDHFCuts::IsSelected;
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel);
  using AliRDHFCuts::IsSelectedPID;
  virtual Int_t IsSelectedPID(AliAODRecoDecayHF* obj);
  Int_t IsSelectedCombinedPID(AliAODRecoDecayHF* obj);
  Double_t GetPionProbabilityTPCTOF(AliAODTrack *trk);

  void SetPIDStrategy(EPIDStrategy pidStrategy){fPIDStrategy=pidStrategy;}
  EPIDStrategy GetPIDStrategy() const {return fPIDStrategy;}
  void SetCombinedPIDThreshold(Double_t a){fCombinedPIDThreshold=a;}
  Double_t GetCombinedPIDThreshold(){return fCombinedPIDThreshold;}


  Bool_t SingleTrkCuts(AliAODTrack *trk);
  Bool_t SingleCascadeCuts(AliAODcascade *casc);
  Bool_t SelectWithRoughCuts(AliAODcascade *casc, AliAODTrack *trk1, AliAODTrack *trk2);

  void SetProdTrackPtMin(Double_t a){fProdTrackPtMin=a;}
  void SetProdTrackEtaRange(Double_t a){fProdTrackEtaRange=a;}
  void SetProdUseAODFilterBit(Bool_t a){fProdUseAODFilterBit=a;}
  void SetProdMassTolLambda(Double_t a){fProdMassTolLambda=a;}
  void SetProdMassTolXi(Double_t a){fProdMassTolXi=a;}
  void SetProdRfidMinV0(Double_t a){fProdRfidMinV0=a;}
  void SetProdRfidMaxV0(Double_t a){fProdRfidMaxV0=a;}
  void SetProdRfidMinXi(Double_t a){fProdRfidMinXi=a;}
  void SetProdRfidMaxXi(Double_t a){fProdRfidMaxXi=a;}
  void SetProdRoughMassTol(Double_t a){fProdRoughMassTol=a;}
  void SetProdRoughPtMin(Double_t a){fProdRoughPtMin=a;}
  void SetProdLikeSignDcaMax(Double_t a){fProdLikeSignDcaMax=a;}

  Double_t GetProdTrackPtMin(){return fProdTrackPtMin;}
  Double_t GetProdTrackEtaRange(){return fProdTrackEtaRange;}
  Bool_t   GetProdUseAODFilterBit(){return fProdUseAODFilterBit;}
  Double_t GetProdMassTolLambda(){return fProdMassTolLambda;}
  Double_t GetProdMassTolXi(){return fProdMassTolXi;}
  Double_t GetProdRfidMinV0(){return fProdRfidMinV0;}
  Double_t GetProdRfidMaxV0(){return fProdRfidMaxV0;}
  Double_t GetProdRfidMinXi(){return fProdRfidMinXi;}
  Double_t GetProdRfidMaxXi(){return fProdRfidMaxXi;}
  Double_t GetProdRoughMassTol(){return fProdRoughMassTol;}
  Double_t GetProdRoughPtMin(){return fProdRoughPtMin;}
  Double_t GetProdLikeSignDcaMax(){return fProdLikeSignDcaMax;}

  void  SetNCuts(Int_t ncuts){fnCuts=ncuts;}
  Int_t GetNCutsArray(){return fnCuts;}
  void  SetCutsArray(Int_t nCuts, Int_t nVars,Int_t nPtBins,Float_t ***cutsRD);
  void  SetCutsArray(Int_t nTotBins,Float_t *cutsRD);
  void  SetCutsFromArray(Int_t nCuts);
  Int_t GetCutArrayID(Int_t ic,Int_t iv,Int_t ip);

 protected:
	
 private:

  EPIDStrategy fPIDStrategy;        //PID Strategy
  Double_t fCombinedPIDThreshold;   //PID threshold used in IsSelectedCombinedPID

  Double_t fProdTrackPtMin;         //Minimum Bachelor pT 
  Double_t fProdTrackEtaRange;      //Bachelor Eta range
  Bool_t   fProdUseAODFilterBit;    //Use AODfilterBit or not
  Double_t fProdMassTolLambda;      //Tolerance of Lambda mass from PDG value
  Double_t fProdMassTolXi;          //Tolerance of Xi mass from PDG value
  Double_t fProdRfidMinV0;          //Minimum Decay vertex of V0
  Double_t fProdRfidMaxV0;          //Max Decay vertex of V0
  Double_t fProdRfidMinXi;          //Minimum Decay vertex of Xi
  Double_t fProdRfidMaxXi;          //Max Decay vertex of Xi
  Double_t fProdLikeSignDcaMax;     //Maximum DCA of pions
  Double_t fProdRoughMassTol;       //Tolerance of Xic mass from PDG value 
  Double_t fProdRoughPtMin;         //Minimum pT of Xic

  Int_t fnCuts;                    //Number of Cuts
  Int_t fnTotalCutBins;            //fnCuts * fnvars * fnPtBins
  Float_t *fCutsArray;             //[fnTotalCutBins]

  ClassDef(AliRDHFCutsXicPlustoXiPiPifromAODtracks,1); 
};

#endif
