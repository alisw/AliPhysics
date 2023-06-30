

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
#ifndef AliAnalysisTaskfnAOD_H
#define AliAnalysisTaskfnAOD_H

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


class AliAnalysisTaskfnAOD : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskfnAOD();
  AliAnalysisTaskfnAOD(const char *name);
  virtual ~AliAnalysisTaskfnAOD();
  
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
  Bool_t      IsV0(AliAODv0 *v1, AliAODEvent *esd);

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
  AliAODEvent      *lAODevent;//! aod Event                                                                                 
  AliEventCuts fEventCuts; //!
  AliESDtrackCuts  *fESDtrackCuts;//!                                                                                                        
  TH1F    *fHistVz;//!
  TH1F    *fHistCentrality;//!
  TH1F    *fHisteventsummary;//!
  THnSparseD    *f1Unlike;//!
  THnSparseD    *f1Like;//!
  THnSparseD    *f1Mix;//!




  AliAnalysisTaskfnAOD(const AliAnalysisTaskfnAOD&);
  AliAnalysisTaskfnAOD& operator=(const AliAnalysisTaskfnAOD&);  
  ClassDef(AliAnalysisTaskfnAOD, 1);
};

//taken from https://github.com/alisw/AliPhysics/blob/master/PWGCF/Correlations/DPhi/PidPid/AliAnalysisTaskPidPidCorrelations.h
class AliReducedTrack : public TObject
{
 public:
  AliReducedTrack(Double_t px, Double_t py, Double_t pz, Short_t charge)
    : fPxReduced(px), fPyReduced(py), fPzReduced(pz), fChargeReduced(charge)
  {
  }
  ~AliReducedTrack() {}
   

  virtual Double_t Px() const { return fPxReduced; }
  virtual Double_t Py() const { return fPyReduced; }
  virtual Double_t Pz() const { return fPzReduced; }
  virtual Short_t Charge() const{ return fChargeReduced; }
    
 private:
    
  Double_t fPxReduced;    
  Double_t fPyReduced;    
  Double_t fPzReduced;    
  Short_t fChargeReduced; 

  ClassDef(AliReducedTrack, 1);
};

#endif

