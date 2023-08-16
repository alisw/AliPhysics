

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
#ifndef AliAnalysisTaskfn_H
#define AliAnalysisTaskfn_H

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


class AliAnalysisTaskfn : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskfn();
  AliAnalysisTaskfn(const char *name);
  virtual ~AliAnalysisTaskfn();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  Bool_t         GoodEvent(const AliVVertex *vertex);
  void           EventMixing();
  Bool_t         TrackPassesOOBPileupCut(AliESDtrack* t, Double_t b);
  

//---------------------------------------------------------------------------------------
  Bool_t      IsPion(AliVTrack *esdtrack);
  Bool_t      IsKaon(AliVTrack *esdtrack);
  Bool_t      IsV0(AliESDv0 *v1, AliESDEvent *esd);

  struct AlikkshPair{
    Int_t charge;
    Int_t trkid;
    TLorentzVector particle;
  };

  struct Alipi{
    Int_t charge;
    Int_t trkid;
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
  AliEventCuts fEventCuts; //!
  AliESDtrackCuts  *fESDtrackCuts;//!                                                                                                        
  TH1F    *fHistVz;//!
  TH1F    *fHistCentrality;//!
  TH1F    *fHisteventsummary;//!
  THnSparseD    *f1Unlike;//!
  THnSparseD    *f1Like;//!
  THnSparseD    *f1Mix;//!


  AliAnalysisTaskfn(const AliAnalysisTaskfn&);
  AliAnalysisTaskfn& operator=(const AliAnalysisTaskfn&);  
  ClassDef(AliAnalysisTaskfn, 1);
};

//taken from https://github.com/alisw/AliPhysics/blob/master/PWGCF/Correlations/DPhi/PidPid/AliAnalysisTaskPidPidCorrelations.h
class AliCompactTrack : public TObject
{
 public:
  AliCompactTrack(Double_t px, Double_t py, Double_t pz, Short_t charge)
    : fPxCompact(px), fPyCompact(py), fPzCompact(pz), fChargeCompact(charge)
  {
  }
  ~AliCompactTrack() {}
   

  virtual Double_t Px() const { return fPxCompact; }
  virtual Double_t Py() const { return fPyCompact; }
  virtual Double_t Pz() const { return fPzCompact; }
  virtual Short_t Charge() const{ return fChargeCompact; }
    
 private:
    
  Double_t fPxCompact;    
  Double_t fPyCompact;    
  Double_t fPzCompact;    
  Short_t fChargeCompact; 

  ClassDef(AliCompactTrack, 1);
};

#endif

