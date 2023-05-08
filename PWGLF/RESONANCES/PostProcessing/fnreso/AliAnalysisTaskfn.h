

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
#include "AliInputEventHandler.h" // event mixing
#include "AliEventPoolManager.h"  // event mixing
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
  Bool_t         AcceptESDtracks(AliESDtrack *pTrack, AliESDtrack *nTrack);    
  void           SetupForMixing();
  Bool_t         TrackPassesOOBPileupCut(AliESDtrack* t, Double_t b);
  

//---------------------------------------------------------------------------------------
  Bool_t      IsPion(AliVTrack *esdtrack);
  Bool_t      IsKaon(AliVTrack *esdtrack);
  Bool_t      CheckESDV0(AliESDv0 *v1, AliESDEvent *esd);

  struct Alikks0Container{
    Int_t charge;
    Int_t tracknumber;
    TLorentzVector particle;
  };

  struct AlipionContainer{
    Int_t charge;
    Int_t tracknumber;
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
  TH1F    *fHistZVertex;//!
  TH1F    *fHistCentralityEvtCount;//!
  TH1F    *fHisteventsummary;//!
  THnSparseD    *f1Unlike;//!
  THnSparseD    *f1Like;//!
  THnSparseD    *f1Mix;//!


  AliAnalysisTaskfn(const AliAnalysisTaskfn&);
  AliAnalysisTaskfn& operator=(const AliAnalysisTaskfn&);  
  ClassDef(AliAnalysisTaskfn, 1);
};

//_____ Reduced Tracks -- contains only quantities requires for this analysis to reduce memory consumption for event mixing
class AliCompactTrack : public TObject // TObject
{
 public:
  AliCompactTrack(Double_t px, Double_t py, Double_t pz, Short_t charge)
    : fPxReduced(px), fPyReduced(py), fPzReduced(pz), fChargeReduced(charge)
  {
  }
  ~AliCompactTrack() {}
   
  // AliVParticle functions

  virtual Double_t Px() const { return fPxReduced; }
  virtual Double_t Py() const { return fPyReduced; }
  virtual Double_t Pz() const { return fPzReduced; }
  virtual Short_t Charge() const{ return fChargeReduced; }
    
 private:
    
  Double_t fPxReduced;    // eta
  Double_t fPyReduced;     // phi
  Double_t fPzReduced;      // pT
  Short_t fChargeReduced;  // charge

  ClassDef(AliCompactTrack, 1); // reduced track class which contains only quantities requires 
};

#endif

