

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
#ifndef AliAnalysisTaskSAphi_H
#define AliAnalysisTaskSAphi_H

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


class AliAnalysisTaskSAphi : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskSAphi();
  AliAnalysisTaskSAphi(const char *name);
  virtual ~AliAnalysisTaskSAphi();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  Bool_t         GoodEvent(const AliVVertex *vertex);
  void           EventMixing();
  void  SetListForTrkCorr(TList *flist)      {this->fListTRKCorr = (TList *) flist->Clone(); }
  void  SetListForNUACorr(TList *flist)      {this->fListNUACorr = (TList *) flist->Clone(); }
  void  SetListForV0MCorr(TList *flist)      {this->fListV0MCorr = (TList *) flist->Clone(); }
  

//---------------------------------------------------------------------------------------
  Bool_t      IsPion(AliVTrack *aodtrack);
  Bool_t      IsKaon(AliVTrack *aodtrack);
  Bool_t      HasTOF(AliAODTrack *track);
  Double_t CosThetaStar(TLorentzVector mother, TLorentzVector daughter0, TLorentzVector daughter1);
  Double_t CosThetaStarHel(TLorentzVector mother, TLorentzVector daughter0, TLorentzVector daughter1);
  
  void  GetNUACorrectionHist(Int_t run=0,Int_t kParticleID=0);
  void  GetEVNTWGTCorrectionHist(Int_t run=0,Int_t kParticleID=0);
  //void  GetV0MCorrectionHist(Int_t run=0);
  void  GetV0MCorrectionHist(Int_t run=0,Int_t kParticleID=0);
  //void  GetMCCorrectionHist(Int_t run=0);
  void  GetMCCorrectionHist(Int_t run=0,Float_t centr=0);

  //--------------------------------------------------------------------

  void  Setframe(Int_t fr)      {this->frame = fr; }



  struct Alikaon{
    Int_t charge;
    Int_t trkid;
    TLorentzVector particle;
  };

  struct Alipion{
    Int_t charge;
    Int_t trkid;
    TLorentzVector particle;
  };

 private:
  enum
  {
    kMaxTrack=2000
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
  TH1F    *fHisteventmult;//!
  THnSparseD    *kstarUnlike;//!
  THnSparseD    *kstarLike;//!
  //THnSparseD    *kstarposLike;//!
  //THnSparseD    *kstarnegLike;//!
  THnSparseD    *kstarMix;//!
  //TH1D    *fHistpionpt;//!
  //TH1D    *fHistkaonpt;//!
  /*TH1D    *fHistnsigtpcpion;
  TH1D    *fHistnsigtpckaon;
  TH1D    *fHistnsigtofpion;
  TH1D    *fHistnsigtofkaon;*/
  Int_t frame;
  TList                 *fListTRKCorr;        //  Supplied from Task
  TList                 *fListNUACorr;        //  Supplied from Task
  TList                 *fListV0MCorr;        //  Supplied from Task  

  ///Used For Corrections:
  TH1D          *fHCorrectMCposChrg;    //! 
  TH1D          *fHCorrectMCnegChrg;    //! 
  TH3F          *fHCorrectNUAposChrg;   //! 
  TH3F          *fHCorrectNUAnegChrg;   //! 

  TH1F    *fHCorrectEVNTWGTChrg;   //!   //eventwgt for charge


  AliAnalysisTaskSAphi(const AliAnalysisTaskSAphi&);
  AliAnalysisTaskSAphi& operator=(const AliAnalysisTaskSAphi&);  
  ClassDef(AliAnalysisTaskSAphi, 1);
};

//taken from https://github.com/alisw/AliPhysics/blob/master/PWGCF/Correlations/DPhi/PidPid/AliAnalysisTaskPidPidCorrelations.h
class AliCompphiTrack : public TObject
{
 public:
  AliCompphiTrack(Double_t px, Double_t py, Double_t pz, Short_t charge)
    : fPxCompphi(px), fPyCompphi(py), fPzCompphi(pz), fChargeCompphi(charge)
  {
  }
  ~AliCompphiTrack() {}
   

  virtual Double_t Px() const { return fPxCompphi; }
  virtual Double_t Py() const { return fPyCompphi; }
  virtual Double_t Pz() const { return fPzCompphi; }
  virtual Short_t Charge() const{ return fChargeCompphi; }
    
 private:
    
  Double_t fPxCompphi;    
  Double_t fPyCompphi;    
  Double_t fPzCompphi;    
  Short_t fChargeCompphi; 

  ClassDef(AliCompphiTrack, 1);
};

#endif

