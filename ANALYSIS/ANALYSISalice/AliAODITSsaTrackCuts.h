#ifndef ALIAODITSSATRACKCUTS_H
#define ALIAODITSSATRACKCUTS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// This class applies the ITSsa cuts at the AOD level.
// Needed for MuonCalo pass where the FilterBit information was not properly saved.
// It contains also some quality cuts which can be modifed by user.
//
// Author: Igor Lakomov <Igor.Lakomov@cern.ch>
//

#include <TMath.h>
#include <TFormula.h>
#include "AliLog.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliESDtrackCuts.h"

#include "AliVEvent.h"
#include "AliVVertex.h"

class AliAODITSsaTrackCuts : public AliVCuts
{
 public:
  AliAODITSsaTrackCuts();
  virtual ~AliAODITSsaTrackCuts();

  Bool_t AcceptTrack(const AliAODTrack* aodTrack);				//returns kTRUE if aodTrack passed all the cuts
  virtual Bool_t IsSelected(TObject* obj) {return AcceptTrack((AliAODTrack*)obj);}
  Double_t CalculateDCAXY(const AliAODTrack* aodTrack);				//calculates the DCAZ using coordinates of aodTrack and fPrimaryVertex
  Double_t CalculateDCAZ(const AliAODTrack* aodTrack);				//calculates the DCAZ using coordinates of aodTrack and fPrimaryVertex

//Standard cuts definitions
  static AliAODITSsaTrackCuts* GetStandardAODITSsaTrackCuts2015();		//definition of the standard cuts for the MuonCalo pass2 of 2015 data

//Setters
  void SetMinNClustersITS(Int_t min=-1) {fMinNClustersITS=min;}			//sets minimum number of ITS clusters
  void SetMaxChi2PerClustersITS(Double_t max=-1.) {fMaxChi2PerClustersITS=max;}	//sets max chi2 per ITS cluster
  void SetDefaultDCAXYptdepCut2015();						//defines default pt-dependent cut on DCAXY: (0.0231+0.0315/x^1.3) // 7*(0.0033+0.0045/pt^1.3)
  void SetDefaultDCAZptdepCut2015();						//defines default pt-dependent cut on DCAZ: (const 1)
  void SetUserDCAXYptdepCut(const char *formula);				//defines the user's pt-dependent cut on DCAXY 
  void SetUserDCAZptdepCut(const char *formula);				//defines the user's pt-dependent cut on DCAZ

//Getters
  Double_t GetMinNClustersITS() {return fMinNClustersITS;}			//returns minimum number of ITS clusters
  Double_t GetMaxChi2PerClustersITS() {return fMaxChi2PerClustersITS;}		//returns max chi2 per ITS cluster
  TFormula* GetDCAXYCut() {return fdcaxycut;}					//returns TFormula defining the pt-dependent cut on DCAXY
  TFormula* GetDCAZCut() {return fdcazcut;}					//returns TFormula defining the pt-dependent cut on DCAZ

 private:
  Double_t fMinNClustersITS;							//minimum number of ITS clusters
  Double_t fMaxChi2PerClustersITS;						//max chi2 per ITS cluster
  TFormula *fdcaxycut;								//TFormula defining the pt-dependent cut on DCAXY
  TFormula *fdcazcut;								//TFormula defining the pt-dependent cut on DCAZ
  const AliVVertex* fPrimaryVertex;						//! Primary vertex

  void ExtractAndSetPrimaryVertex(AliVEvent *event) {fPrimaryVertex = event->GetPrimaryVertex();}	//Extracts and sets primary vertex from the AOD event

  ClassDef(AliAODITSsaTrackCuts,1)
};
#endif
