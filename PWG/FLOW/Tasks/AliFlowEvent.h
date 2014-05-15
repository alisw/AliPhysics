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
class TH2F;

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

  virtual void Get2Qsub(AliFlowVector* Qarray, Int_t n = 2, TList *weightsList = 0x0, Bool_t usePhiWeights = 0x0, Bool_t usePtWeights = 0x0, Bool_t useEtaWeights = 0x0);
  void SetVZEROCalibrationForTrackCuts(AliFlowTrackCuts* cuts);
  void SetVZEROCalibrationForTrackCuts2011(AliFlowTrackCuts* cuts);

  virtual void ClearFast();

protected:
  AliFlowTrack* ReuseTrack( Int_t i);

private:
  Int_t         fApplyRecentering;      // apply recentering of q-vectors? 2010 is 10h style, 2011 is 11h style
  Int_t         fCachedRun;             //! cached calibration info for vzero
  Int_t         fVZEROcentralityBin;     //! centrality bin for the current event 
  Float_t       fMeanQ[9][2][2];        //! recentering
  Float_t       fWidthQ[9][2][2];       //! recentering
  Float_t       fMeanQv3[9][2][2];      //! recentering
  Float_t       fWidthQv3[9][2][2];     //! recentering


  ClassDef(AliFlowEvent,3)
};

#endif


