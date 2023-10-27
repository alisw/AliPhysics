

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
#ifndef AliAnalysisTaskEP_H
#define AliAnalysisTaskEP_H

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


class AliAnalysisTaskEP : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskEP();
  AliAnalysisTaskEP(const char *name);
  virtual ~AliAnalysisTaskEP();
  
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
  
  void  GetNUACorrectionHist(Int_t run=0,Int_t kParticleID=0);
  void  GetV0MCorrectionHist(Int_t run=0,Int_t kParticleID=0);
  void  GetMCCorrectionHist(Int_t run=0,Float_t centr=0);

  //--------------------------------------------------------------------

  void  Setframe(Int_t fr)      {this->frame = fr; }



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
  Int_t frame;
  Int_t         gPsiNSet;   //
  TList         *fListTRKCorr;        //  Supplied from Task
  TList         *fListNUACorr;        //  Supplied from Task
  TList         *fListV0MCorr;        //  Supplied from Task  
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
  TH2F   *fHCorrectV0ChWgt;  //!
  TH1D    *fHCorrectQNxV0C;  //!
  TH1D    *fHCorrectQNyV0C;  //!
  TH1D    *fHCorrectQNxV0A;  //!
  TH1D    *fHCorrectQNyV0A;  //!
  ///Used For Corrections:
  TH1D          *fHCorrectMCposChrg;    //! 
  TH1D          *fHCorrectMCnegChrg;    //! 
  TH3F          *fHCorrectNUAposChrg;   //! 
  TH3F          *fHCorrectNUAnegChrg;   //! 
  TH1F    *fHCorrectEVNTWGTChrg;   //!   //eventwgt for charge


  ///Rihan:Functions for V0EP:
  void    ApplyV0XqVectRecenter(Float_t fCent,Int_t gPsiN,Double_t &qnxV0C, Double_t &qnyV0C, Double_t &qnxV0A, Double_t &qnyV0A);
  Bool_t  GetGainCorrectedV0Qvector(AliAODEvent *faod,Double_t fVtxZ,Int_t gPsiN,Double_t &qnxV0C,Double_t &qnyV0C,Double_t &qnxV0A,Double_t &qnyV0A); 


  AliAnalysisTaskEP(const AliAnalysisTaskEP&);
  AliAnalysisTaskEP& operator=(const AliAnalysisTaskEP&);  
  ClassDef(AliAnalysisTaskEP, 1);
};


#endif

