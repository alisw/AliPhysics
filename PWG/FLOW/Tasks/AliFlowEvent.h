/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

/*****************************************************************
  AliFlowEvent: Event container for flow analysis                  
                                     
  origin:   Mikolaj Krzewicki  (mikolaj.krzewicki@cern.ch)
  mods:     Redmer A. Bertens (rbertens@cern.ch)
*****************************************************************/

#ifndef ALIFLOWEVENT_H
#define ALIFLOWEVENT_H

class AliFlowTrackCuts;
class AliFlowTrack;
class AliCFManager;
class AliVEvent;
class AliMCEvent;
class AliESDEvent;
class AliAODEvent;
class AliMultiplicity;
class AliESDPmdTrack;
class AliFlowVector;
class TH1F;
class TH1;
class TH2F;
class TArrayD;

#include "AliFlowEventSimple.h"

class AliFlowEvent: public AliFlowEventSimple {
public:

  enum KineSource { kNoKine, kESDkine, kMCkine };

  AliFlowEvent();
  AliFlowEvent(Int_t n);
  AliFlowEvent(const AliFlowEvent& event);
  AliFlowEvent& operator=(const AliFlowEvent& event);
  virtual  ~AliFlowEvent() {}

  //deprecated
  AliFlowEvent( const AliMCEvent* anInput,
                const AliCFManager* rpCFManager=NULL,
                const AliCFManager* poiCFManager=NULL ); //use CF(2x)
  AliFlowEvent( const AliESDEvent* anInput,  
                const AliCFManager* rpCFManager=NULL, 
                const AliCFManager* poiCFManager=NULL ); //use CF(2x)
  AliFlowEvent( const AliAODEvent* anInput, 
                const AliCFManager* rpCFManager=NULL, 
                const AliCFManager* poiCFManager=NULL );  //use CF(2x)
  AliFlowEvent( const AliESDEvent* anInput,
                const AliCFManager* poiCFManager,
                Bool_t hybrid);
  AliFlowEvent( const AliESDEvent* anInput, 
                const AliMCEvent* anInputMc, 
                KineSource anOption=kNoKine,
                const AliCFManager* rpCFManager=NULL, 
                const AliCFManager* poiCFManager=NULL );  //use CF(2x)
  AliFlowEvent( const AliESDEvent* anInput,
                const AliMultiplicity* anInputTracklets,
                const AliCFManager* poiCFManager );
  AliFlowEvent( const AliESDEvent* anInput,
                const TH2F* anInputFMDhist,
                const AliCFManager* poiCFManager );
  //pmd
  AliFlowEvent( const AliESDEvent* anInput,
                const AliESDPmdTrack *pmdtracks,
                const AliCFManager* poiCFManager );
  //pmd
  //end of deprecated

  AliFlowEvent( AliFlowTrackCuts* rpCuts,
                AliFlowTrackCuts* poiCuts );
  
  void Fill( AliFlowTrackCuts* rpCuts,
             AliFlowTrackCuts* poiCuts );

  void FindDaughters(Bool_t keepDaughtersInRPselection=kFALSE);

  void SetMCReactionPlaneAngle(const AliMCEvent* mcEvent);
  using AliFlowEventSimple::SetMCReactionPlaneAngle;

  AliFlowTrack* GetTrack( Int_t i );

  void InsertTrack(AliFlowTrack*);

  virtual AliFlowVector GetQ(Int_t n=2, TList *weightsList=NULL, Bool_t usePhiWeights=kFALSE, Bool_t usePtWeights=kFALSE, Bool_t useEtaWeights=kFALSE);
  virtual void Get2Qsub(AliFlowVector* Qarray, Int_t n = 2, TList *weightsList = 0x0, Bool_t usePhiWeights = 0x0, Bool_t usePtWeights = 0x0, Bool_t useEtaWeights = 0x0);
  void SetVZEROCalibrationForTrackCuts(AliFlowTrackCuts* cuts);
  void SetBetaVZEROCalibrationForTrackCuts(AliFlowTrackCuts* cuts);
  void SetDeltaVZEROCalibrationForTrackCuts(AliFlowTrackCuts* cuts);
  void SetKappaVZEROCalibrationForTrackCuts(AliFlowTrackCuts* cuts);
  void SetHotfixVZEROCalibrationForTrackCuts(AliFlowTrackCuts* cuts);

  virtual void ClearFast();
  virtual void ClearCachedRun();

protected:
  AliFlowTrack* ReuseTrack( Int_t i);

private:
  Int_t         fApplyRecentering;      // apply recentering of q-vectors? 2010 is 10h style, 2011 is 11h style
  Bool_t        fDivSigma;              // divide by st.dev. after recentering?
  Int_t         fCachedRun;             //! cached calibration info for vzero
  Int_t         fVZEROcentralityBin;    //! centrality bin for the current event
  Float_t       fMeanQ[9][2][2];        //! recentering
  Float_t       fWidthQ[9][2][2];       //! recentering
  Float_t       fMeanQv3[9][2][2];      //! recentering
  Float_t       fWidthQv3[9][2][2];     //! recentering
  // BETA testing of new VZERO calibration
  TH1*         fQxavsV0[5];            //! recentering
  TH1*         fQyavsV0[5];            //! recentering
  TH1*         fQxcvsV0[5];            //! recentering
  TH1*         fQycvsV0[5];            //! recentering
  // END OF BETA TESTING
  AliVEvent*    fEvent;                 //! current event
  TArrayD*      fChi2A;                 //! chi vs cent for vzero A ep_2
  TArrayD*      fChi2C;                 //! chi vs cent for vzero C ep_2
  TArrayD*      fChi3A;                 //! chi vs cent for vzero A ep_3
  TArrayD*      fChi3C;                 //! chi vs cent for vzero C ep_3

  ClassDef(AliFlowEvent,5)
};

#endif


