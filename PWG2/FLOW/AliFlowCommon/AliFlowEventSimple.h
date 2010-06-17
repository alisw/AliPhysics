/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

/*****************************************************************
  AliFlowEventSimple: A simple event 
  for flow analysis                  
                                     
  origin: Naomi van der Kolk (kolk@nikhef.nl)           
          Ante Bilandzic     (anteb@nikhef.nl)         
          Raimond Snellings  (Raimond.Snellings@nikhef.nl)    
  mods:   Mikolaj Krzewicki  (mikolaj.krzewicki@cern.ch)
*****************************************************************/

#ifndef ALIFLOWEVENTSIMPLE_H
#define ALIFLOWEVENTSIMPLE_H

#include "TObject.h"
#include "TParameter.h"
#include "TMath.h"
class TTree;
class TF1;
class AliFlowVector;
class AliFlowTrackSimple;
class AliFlowTrackSimpleCuts;

class AliFlowEventSimple: public TObject {

 public:
  AliFlowEventSimple();
  AliFlowEventSimple(Int_t aLenght);
  AliFlowEventSimple( Int_t nParticles,
                      TF1* ptDist,
                      Double_t phiMin=0,
                      Double_t phiMax=TMath::TwoPi(),
                      Double_t etaMin=-1.0,
                      Double_t etaMax= 1.0 );
  AliFlowEventSimple(TTree* anInput, const AliFlowTrackSimpleCuts* rpCuts, const AliFlowTrackSimpleCuts* poiCuts);
  AliFlowEventSimple(const AliFlowEventSimple& anEvent);
  AliFlowEventSimple& operator=(const AliFlowEventSimple& anEvent);
  virtual  ~AliFlowEventSimple();

  virtual void Generate( Int_t nParticles,
                         TF1* ptDist,
                         Double_t phiMin=0,
                         Double_t phiMax=TMath::TwoPi(),
                         Double_t etaMin=-1.0,
                         Double_t etaMax= 1.0 );

  Bool_t  IsFolder() const {return kTRUE;};
  void    Browse(TBrowser *b); 
  void    Print(Option_t* option = "") const;      //method to print stats
  
  Int_t    NumberOfTracks() const                   { return fNumberOfTracks; }
  Int_t    GetReferenceMultiplicity() const               { return fReferenceMultiplicity; }
  void     SetReferenceMultiplicity( Int_t m )            { fReferenceMultiplicity = m; }
  Int_t    GetEventNSelTracksRP() const             { return fNumberOfRPs; } 
  void     SetEventNSelTracksRP(Int_t nr)           { fNumberOfRPs = nr; }  
  Double_t GetMCReactionPlaneAngle() const          { return fMCReactionPlaneAngle; }
  void     SetMCReactionPlaneAngle(Double_t fPhiRP) { fMCReactionPlaneAngle=fPhiRP; fMCReactionPlaneAngleIsSet=kTRUE; }
  Bool_t   IsSetMCReactionPlaneAngle() const        { return fMCReactionPlaneAngleIsSet; }
  void     SetAfterBurnerPrecision(Double_t p)      { fAfterBurnerPrecision=p; }
  Double_t GetAfterBurnerPrecision() const          { return fAfterBurnerPrecision; }

  void ResolutionPt(Double_t res);
  void TagSubeventsInEta(Double_t etaMinA, Double_t etaMaxA, Double_t etaMinB, Double_t etaMaxB );
  void TagRP(AliFlowTrackSimpleCuts* cuts );
  void TagPOI(AliFlowTrackSimpleCuts* cuts );
  void CloneTracks(Int_t n);
  void AddV1( Double_t v1 );
  void AddV2( Double_t v2 );
  void AddV4( Double_t v4 );
  void AddFlow( Double_t v1, Double_t v2, Double_t v4 );
  void DefineDeadZone( Double_t etaMin, Double_t etaMax, Double_t phiMin, Double_t phiMax );
  Int_t CleanUpDeadTracks();
 
  AliFlowTrackSimple* GetTrack(Int_t i);
  void AddTrack( AliFlowTrackSimple* track ); 
 
  AliFlowVector GetQ(Int_t n=2, TList *weightsList=NULL, Bool_t usePhiWeights=kFALSE, Bool_t usePtWeights=kFALSE, Bool_t useEtaWeights=kFALSE);
  void Get2Qsub(AliFlowVector* Qarray, Int_t n=2, TList *weightsList=NULL, Bool_t usePhiWeights=kFALSE, Bool_t usePtWeights=kFALSE, Bool_t useEtaWeights=kFALSE);  

 protected:
  TObjArray*              fTrackCollection;           //-> collection of tracks
  Int_t                   fReferenceMultiplicity;           // reference multiplicity
  Int_t                   fNumberOfTracks;            // number of tracks
  Int_t                   fNumberOfRPs;               // number of tracks that have passed the RP selection
  Double_t                fMCReactionPlaneAngle;      // the angle of the reaction plane from the MC truth
  Bool_t                  fMCReactionPlaneAngleIsSet; // did we set it from MC?
  Double_t                fAfterBurnerPrecision;      // iteration precision in afterburner
  TParameter<Int_t>*      fNumberOfTracksWrap;        //! number of tracks in TBrowser
  TParameter<Int_t>*      fNumberOfRPsWrap;           //! number of tracks that have passed the RP selection in TBrowser
  TParameter<Double_t>*   fMCReactionPlaneAngleWrap;  //! the angle of the reaction plane from the MC truth in TBrowser

  ClassDef(AliFlowEventSimple,1)
};

#endif


