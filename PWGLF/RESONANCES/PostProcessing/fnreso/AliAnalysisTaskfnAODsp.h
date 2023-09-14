

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#ifndef AliAnalysisTaskfnAODsp_H
#define AliAnalysisTaskfnAODsp_H

#include "AliAnalysisTaskSE.h"
#include "TString.h"
#include "TTree.h"
#include "TH1F.h"
#include "AliEventCuts.h"
#include "AliVParticle.h"
#include "TObject.h"
#include "AliInputEventHandler.h"
#include "AliEventPoolManager.h" 
#include "THnSparse.h"
#include "AliVVertex.h"

class TList;

class AliESDEvent;
class AliAODEvent;
class AliVEvent;
class AliESDtrackCuts;
class AliAODTrack;
class AliMCEvent;

class TH1F;
class TH2F;
class TH3F;
class TProfile;
class AliPIDResponse;
class AliMultSelection;
class AliPPVsMultUtils;
class AliAnalysisUtils;


class AliAnalysisTaskfnAODsp : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskfnAODsp();
  AliAnalysisTaskfnAODsp(const char *name);
  virtual ~AliAnalysisTaskfnAODsp();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  Bool_t         GoodEvent(const AliVVertex *vertex);
  Bool_t         AcceptAODtracks(AliAODTrack *pTrack, AliAODTrack *nTrack);    
  void           EventMixing();
  Bool_t         TrackPassesOOBPileupCut(AliAODTrack* t, Double_t b);
  

//---------------------------------------------------------------------------------------
  Bool_t      IsPion(AliVTrack *aodtrack);
  Bool_t      IsKaon(AliVTrack *aodtrack);
  Bool_t      HasTOF(AliAODTrack *track);
  Bool_t      IsV0(AliAODv0 *v1, AliAODEvent *esd);
  Double_t CosThetaStarHel0(TLorentzVector mother, TLorentzVector daughter0);
  Double_t CosThetaStarHel1(TLorentzVector mother, TLorentzVector daughter1);
  Double_t CosThetaStarHel2(TLorentzVector mother, TLorentzVector daughter2);


  void SetFilterBit(Int_t fb)                    {this->fFilterBit   =  fb;}
  void Setkkshmasscut(Double_t mc)               {this->kkshmasscut = mc;}
  void SetPIDnsigtpcpion(Double_t nsigtpcpi)      {this->nsigtpcpion = nsigtpcpi;}
  void SetPIDnsigtofpion(Double_t nsigtofpi)      {this->nsigtofpion = nsigtofpi;}
  void SetPIDnsigtpckaon(Double_t nsigtpcka)      {this->nsigtpckaon = nsigtpcka;}
  void SetPIDnsigtofkaon(Double_t nsigtofka)      {this->nsigtofkaon = nsigtofka;}
  void Setdcaxyposneg(Double_t xypos, Double_t xyneg) {this->dcaxypos=xypos; this->dcaxyneg=xyneg;}
  void Setdcav0daugh (Double_t dca1)             {this->dcav0daugh = dca1;}
  void Setdcav0pv(Double_t dca2)             {this->dcav0pv = dca2;}
  void SetCosPA(Double_t cosPA)              {this->cospa = cosPA;}
  void SetLowradius(Double_t lr)              {this->lowrad = lr;}
  void SetLT(Double_t lt)              {this->lifetime = lt;}
  void SetPIDpion(Double_t pidpi)              {this->pidpion = pidpi;}
  void SetPTC(Float_t cr, Float_t crfc, Double_t chi2global, Double_t chi2ITS)     {this->nCRcut=cr; this->ratiocrfccut=crfc; this->chi2globalcut=chi2global; this->chi2cut=chi2ITS;}




  struct AlikkshPair{
    Int_t charge;
    Int_t trkid;
    Int_t trkidv0;
    TLorentzVector particle;
  };

  struct Alipi{
    Int_t charge;
    Int_t trkid;
    TLorentzVector particle;
  };

  struct Alikaonsp{
    Int_t charge;
    Int_t trkid;
    Int_t trkidv0;
    TLorentzVector particle;
  };

  struct Alikshortsp{
    Int_t charge;
    Int_t trkid;
    Int_t trkidv0;
    TLorentzVector particle;
  };


 private:
  enum
  {
    kMaxTrack=200
  };
  AliEventPoolManager*fPoolMgr; //!
  TObjArray* fpionreduced; //!
  TList  *fOutput;//!
  AliPIDResponse   *fPIDResponse;//!
  AliVEvent        *fVevent;//!                                                                                 
  AliESDEvent       *lESDevent;//!
  AliAODEvent      *lAODevent;//! aod Event                                                                                 
  AliEventCuts fEventCuts; //!
  AliESDtrackCuts  *fESDtrackCuts;//!                                                                                                        
  TH1F    *fHistVz;//!
  TH1F    *fHistCentrality;//!
  TH1F    *fHisteventsummary;//!
  THnSparseD    *f1Unlike0;//!
  THnSparseD    *f1Unlike1;//!
  THnSparseD    *f1Unlike2;//!
  THnSparseD    *f1Like0;//!
  THnSparseD    *f1Like1;//!
  THnSparseD    *f1Like2;//!
  THnSparseD    *f1Mix0;//!
  THnSparseD    *f1Mix1;//!
  THnSparseD    *f1Mix2;//!
  /* TH1D            *hist1;
  TH1D            *hist2;
  TH1D            *hist3;
  TH1D            *hist4;
  TH1D            *hist5;
  TH1D            *hist6;
  TH1D            *hist7;*/


  //variables
  Int_t fFilterBit;
  Double_t kkshmasscut;
  Double_t nsigtpcpion;
  Double_t nsigtofpion;
  Double_t nsigtpckaon;
  Double_t nsigtofkaon;
  Double_t dcaxypos;
  Double_t dcaxyneg;
  Double_t dcav0daugh;
  Double_t dcav0pv;
  Double_t cospa;
  Double_t lowrad;
  Double_t lifetime;
  Double_t pidpion;
  Float_t nCRcut;
  Float_t ratiocrfccut;
  Double_t chi2globalcut;
  Double_t chi2cut;

  
  AliAnalysisTaskfnAODsp(const AliAnalysisTaskfnAODsp&);
  AliAnalysisTaskfnAODsp& operator=(const AliAnalysisTaskfnAODsp&);  
  ClassDef(AliAnalysisTaskfnAODsp, 1);
};

//taken from https://github.com/alisw/AliPhysics/blob/master/PWGCF/Correlations/DPhi/PidPid/AliAnalysisTaskPidPidCorrelations.h
class AliReducedspTrack : public TObject
{
 public:
  AliReducedspTrack(Double_t px, Double_t py, Double_t pz, Short_t charge)
    : fPxReducedsp(px), fPyReducedsp(py), fPzReducedsp(pz), fChargeReducedsp(charge)
  {
  }
  ~AliReducedspTrack() {}
   

  virtual Double_t Px() const { return fPxReducedsp; }
  virtual Double_t Py() const { return fPyReducedsp; }
  virtual Double_t Pz() const { return fPzReducedsp; }
  virtual Short_t Charge() const{ return fChargeReducedsp; }
    
 private:
    
  Double_t fPxReducedsp;    
  Double_t fPyReducedsp;    
  Double_t fPzReducedsp;    
  Short_t fChargeReducedsp; 

  ClassDef(AliReducedspTrack, 1);
};

#endif

