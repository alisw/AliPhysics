
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
#ifndef AliAnalysisTaskSAEP_H
#define AliAnalysisTaskSAEP_H

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
class TProfile2D;
class AliPIDResponse;
class AliMultSelection;
class AliPPVsMultUtils;
class AliAnalysisUtils;


class AliAnalysisTaskSAEP : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskSAEP();
  AliAnalysisTaskSAEP(const char *name);
  virtual ~AliAnalysisTaskSAEP();
  
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
  Double_t CosThetaStarEP(TLorentzVector mother, TLorentzVector daughter0, TLorentzVector daughter1, TVector2& Qvect);
  
  void  GetNUACorrectionHist(Int_t run=0,Int_t kParticleID=0);
  void  GetV0MCorrectionHist(Int_t run=0,Int_t kParticleID=0);
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
  THnSparseD    *kstarUnlikeA;//!
  THnSparseD    *kstarUnlikeC;//!
  THnSparseD    *kstarLikeA;//!
  THnSparseD    *kstarLikeC;//!
  //THnSparseD    *kstarposLike;//!
  //THnSparseD    *kstarnegLike;//!
  THnSparseD    *kstarMixA;//!
  THnSparseD    *kstarMixC;//!
  //TH1D    *fHistpionpt;//!
  //TH1D    *fHistkaonpt;//!
  /*TH1D    *fHistnsigtpcpion;
  TH1D    *fHistnsigtpckaon;
  TH1D    *fHistnsigtofpion;
  TH1D    *fHistnsigtofkaon;*/
  Int_t frame;
  Int_t         gPsiNSet;   //
  TH2F   *fHCorrectV0ChWgt;  //!
  TH1D    *fHCorrectQNxV0C;  //!
  TH1D    *fHCorrectQNyV0C;  //!
  TH1D    *fHCorrectQNxV0A;  //!
  TH1D    *fHCorrectQNyV0A;  //!
  TList                 *fListTRKCorr;        //  Supplied from Task
  TList                 *fListNUACorr;        //  Supplied from Task
  TList                 *fListV0MCorr;        //  Supplied from Task  
  TProfile2D *hAvgV0ChannelsvsVz;  //!
  TProfile     *hAvgQNXvsCentV0C;  //!
  TProfile     *hAvgQNYvsCentV0C;  //!
  TProfile     *hAvgQNXvsCentV0A;  //!
  TProfile     *hAvgQNYvsCentV0A;  //!
  
  TH2F     *fHistV0CPsiNEventPlane;  //!
  TH2F     *fHistV0APsiNEventPlane;  //!
  TProfile *hV0CV0APsiNCorrelation;  //!
  TProfile *hV0CTPCPsiNCorrelation;  //!
  TProfile *hV0ATPCPsiNCorrelation;  //!
  ///Used For Corrections:
  TH1D          *fHCorrectMCposChrg;    //! 
  TH1D          *fHCorrectMCnegChrg;    //! 
  TH3F          *fHCorrectNUAposChrg;   //! 
  TH3F          *fHCorrectNUAnegChrg;   //! 

  TH1F    *fHCorrectEVNTWGTChrg;   //!   //eventwgt for charge

  ///Rihan:Functions for V0EP:
  void    ApplyV0XqVectRecenter(Float_t fCent,Int_t gPsiN,Double_t &qnxV0C, Double_t &qnyV0C, Double_t &qnxV0A, Double_t &qnyV0A);
  Bool_t  GetGainCorrectedV0Qvector(AliAODEvent *faod,Double_t fVtxZ,Int_t gPsiN,Double_t &qnxV0C,Double_t &qnyV0C,Double_t &qnxV0A,Double_t &qnyV0A); 



  AliAnalysisTaskSAEP(const AliAnalysisTaskSAEP&);
  AliAnalysisTaskSAEP& operator=(const AliAnalysisTaskSAEP&);  
  ClassDef(AliAnalysisTaskSAEP, 1);
};

//taken from https://github.com/alisw/AliPhysics/blob/master/PWGCF/Correlations/DPhi/PidPid/AliAnalysisTaskPidPidCorrelations.h
class AliCompSATrack : public TObject
{
 public:
  AliCompSATrack(Double_t px, Double_t py, Double_t pz, Short_t charge)
    : fPxCompSA(px), fPyCompSA(py), fPzCompSA(pz), fChargeCompSA(charge)
  {
  }
  ~AliCompSATrack() {}
   

  virtual Double_t Px() const { return fPxCompSA; }
  virtual Double_t Py() const { return fPyCompSA; }
  virtual Double_t Pz() const { return fPzCompSA; }
  virtual Short_t Charge() const{ return fChargeCompSA; }
    
 private:
    
  Double_t fPxCompSA;    
  Double_t fPyCompSA;    
  Double_t fPzCompSA;    
  Short_t fChargeCompSA; 

  ClassDef(AliCompSATrack, 1);
};

#endif

