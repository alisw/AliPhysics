
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
#ifndef Aliresonance_H
#define Aliresonance_H

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


class Aliresonance : public AliAnalysisTaskSE {
 public:
  Aliresonance();
  Aliresonance(const char *name);
  virtual ~Aliresonance();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  Bool_t         GoodEvent(const AliVVertex *vertex);
  Bool_t         AcceptAODtracks(AliAODTrack *pTrack, AliAODTrack *nTrack);    
  void           MixingEvents();
  Bool_t         TrackPassesOOBPileupCut(AliAODTrack* t, Double_t b);
  Bool_t         IsPion(AliVTrack *aodtrack);
  Bool_t         IsKaon(AliVTrack *aodtrack);
  Bool_t         IsV0(AliAODv0 *v1, AliAODEvent *aod);

  struct AlikkshContainer{
    Int_t charge; //for distinguishing like and unlike sign
    Int_t trackid; //to prevent same track counting
    TLorentzVector particle; //particle invmass
  };

  struct AlipiContainer{
    Int_t charge;
    Int_t trackid;
    TLorentzVector particle;
  };

 private:
  enum
  {
    kMaxTrack=500
  };
  AliEventPoolManager *fPoolMgr; //!
  TList  *fOutput;//!
  AliPIDResponse   *fPIDResponse;//!
  AliVEvent        *fVevent;//!VEvent                                                                                 
  AliAODEvent      *lAODevent;//! aod Event                                                                                 
  AliEventCuts fEventCuts; //!
  TH1F    *fHistZVertex;//!
  TH1F    *fHistCentralityEvtCount;//!
  TH1F    *fHisteventsummary;//!
  THnSparseD    *f1Unlike;//!
  THnSparseD    *f1Like;//!
  THnSparseD    *f1Mix;//!


  Aliresonance(const Aliresonance&);
  Aliresonance& operator=(const Aliresonance&);  
  ClassDef(Aliresonance, 1);
};


class AliCompactTrack : public TObject // TObject
{
 public:
  AliCompactTrack(Double_t px, Double_t py, Double_t pz, Short_t charge)
    : PX(px), PY(py), PZ(pz), Chrg(charge)
  {
  }
  ~AliCompactTrack() {}
   

  virtual Double_t Px() const { return PX; }
  virtual Double_t Py() const { return PY; }
  virtual Double_t Pz() const { return PZ; }
  virtual Short_t Charge() const{ return Chrg; }
    
 private:
    
  Double_t PX;    
  Double_t PY;    
  Double_t PZ;    
  Short_t Chrg; 

  ClassDef(AliCompactTrack, 1); 
};

#endif

