/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliFlowEventSimple_H
#define AliFlowEventSimple_H

#include "AliFlowVector.h" //needed as include
#include "TH1F.h"
#include "TH1D.h"
#include "TFile.h"
class AliFlowTrackSimple;

/**************************************
 * AliFlowEventSimple: A simple event *
 *  for flow analysis                 * 
 *                                    * 
 * authors: Naomi van der Kolk        *
 *           (kolk@nikhef.nl)         *  
 *          Ante Bilandzic            *
 *           (anteb@nikhef.nl)        * 
 * ***********************************/

class AliFlowEventSimple: public TObject {

 public:
  AliFlowEventSimple(Int_t aLenght);
  AliFlowEventSimple(const AliFlowEventSimple& anEvent);
  AliFlowEventSimple& operator=(const AliFlowEventSimple& anEvent);
  virtual  ~AliFlowEventSimple();
  
  Int_t NumberOfTracks() const              { return this->fNumberOfTracks; }
  Int_t GetEventNSelTracksIntFlow() const   { return this->fEventNSelTracksIntFlow; }
  void SetNumberOfTracks(Int_t nr)          { this->fNumberOfTracks = nr; }
  void SetEventNSelTracksIntFlow(Int_t nr)  { this->fEventNSelTracksIntFlow = nr; }
 
  //setters and getters
  void SetUseWeightsPhi(Bool_t const uwPhi) { this->fUseWeightsPhi = uwPhi;};
  Bool_t GetUseWeightsPhi() const { return this->fUseWeightsPhi;};
  
  void SetUseWeightsPt(Bool_t const uwPt) { this->fUseWeightsPt = uwPt;};
  Bool_t GetUseWeightsPt() const { return this->fUseWeightsPt;};
  
  void SetUseWeightsEta(Bool_t const uwEta) { this->fUseWeightsEta = uwEta;};
  Bool_t GetUseWeightsEta() const { return this->fUseWeightsEta;};
  
  void SetPhiWeights(TH1F* const fphiw)     { this->fPhiWeights = fphiw;};
  void SetPtWeights(TH1D* const fptw)       { this->fPtWeights = fptw;};
  void SetEtaWeights(TH1D* const fetaw)     { this->fEtaWeights = fetaw;};
  
  AliFlowTrackSimple* GetTrack(Int_t i);
  TObjArray* TrackCollection() const        { return this->fTrackCollection; }
  AliFlowVector GetQ(Int_t n=2);                //get Q-vector evaluated in harmonic n (default harmonic n=2) 
  
 private:
  TObjArray*           fTrackCollection;        // collection of tracks
  Int_t                fNumberOfTracks;         // number of tracks
  Int_t                fEventNSelTracksIntFlow; // number of tracks selected for integrated flow calculation
  
  Bool_t fUseWeightsPhi;       // phi weights
  Bool_t fUseWeightsPt;        // v_2(pt) weights
  Bool_t fUseWeightsEta;       // v_2(eta) weights   
  TH1F *fPhiWeights;           // phi weights
  TH1D *fPtWeights;            // v_2(pt) weights   
  TH1D *fEtaWeights;           // v_2(eta) weights
  
  ClassDef(AliFlowEventSimple,1)                // macro for rootcint
};

#endif


