/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliFlowEventSimple_H
#define AliFlowEventSimple_H

#include "AliFlowVector.h" //needed as include
#include "TH1F.h"
#include "TH1D.h"
#include "TFile.h"
#include "TParameter.h"
class AliFlowTrackSimple;

/**************************************
 * AliFlowEventSimple: A simple event *
 *  for flow analysis                 * 
 *                                    * 
 * authors: Naomi van der Kolk        *
 *           (kolk@nikhef.nl)         *  
 *          Ante Bilandzic            *
 *           (anteb@nikhef.nl)        * 
 *          Raimond Snellings         *
 *    (Raimond.Snellings@nikhef.nl)   * 
 * ***********************************/

class AliFlowEventSimple: public TObject {

 public:
  AliFlowEventSimple();
  AliFlowEventSimple(Int_t aLenght);
  AliFlowEventSimple(const AliFlowEventSimple& anEvent);
  AliFlowEventSimple& operator=(const AliFlowEventSimple& anEvent);
  virtual  ~AliFlowEventSimple();

  Bool_t  IsFolder() const {return kTRUE;};
  void Browse(TBrowser *b); 
  void Print(Option_t* option = "") const;      //method to print stats
  
  Int_t NumberOfTracks() const              { return this->fNumberOfTracks; }
  void SetNumberOfTracks(Int_t nr)          { this->fNumberOfTracks = nr; }
  Int_t GetEventNSelTracksIntFlow() const   { return this->fEventNSelTracksIntFlow; }
  void SetEventNSelTracksIntFlow(Int_t nr)  { this->fEventNSelTracksIntFlow = nr; }
  Double_t GetMCReactionPlaneAngle() const   { return this->fMCReactionPlaneAngle; }
  void SetMCReactionPlaneAngle(Double_t fPhiRP)  { this->fMCReactionPlaneAngle = fPhiRP; }
 
  AliFlowTrackSimple* GetTrack(Int_t i);
  TObjArray* TrackCollection() const        { return this->fTrackCollection; }
 
  AliFlowVector GetQ(Int_t n=2, TList *weightsList=NULL, Bool_t usePhiWeights=kFALSE, Bool_t usePtWeights=kFALSE, Bool_t useEtaWeights=kFALSE);
  AliFlowVector GetQsub(Double_t etaMin=-0.9, Double_t etaMax=0.9, Int_t n=2, TList *weightsList=NULL, Bool_t usePhiWeights=kFALSE, Bool_t usePtWeights=kFALSE, Bool_t useEtaWeights=kFALSE);  //for eta subevents
  //add here more subevent options   

 private:
  TObjArray*              fTrackCollection;        // collection of tracks
  Int_t                   fNumberOfTracks;         // number of tracks
  Int_t                   fEventNSelTracksIntFlow; // number of tracks selected for integrated flow calculation
  Double_t                fMCReactionPlaneAngle;   // the angle of the reaction plane from the MC truth
  TParameter<Int_t>*      fNumberOfTracksWrap;         //! number of tracks in TBrowser
  TParameter<Int_t>*      fEventNSelTracksIntFlowWrap; //! number of tracks selected for integrated flow calculation in TBrowser
  TParameter<Double_t>*   fMCReactionPlaneAngleWrap;   //! the angle of the reaction plane from the MC truth in TBrowser
  ClassDef(AliFlowEventSimple,1)                   // simplified event used in flow analysis 
};

#endif


