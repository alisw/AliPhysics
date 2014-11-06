#ifndef ALIJETFILLCALTRKTRACKMC_H
#define ALIJETFILLCALTRKTRACKMC_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
/* $Id$ */

//--------------------------------------------------
// Filling of CalTrkTrack objects in the MC reader task
//
// Author: magali.estienne@subatech.in2p3.fr
//         alexandre.shabetai@cern.ch
//-------------------------------------------------

#include "AliJetFillCalTrkEvent.h"

class AliVEvent;
class AliMCEvent;

class AliJetFillCalTrkTrackMC : public AliJetFillCalTrkEvent
{
 public: 
  AliJetFillCalTrkTrackMC();
  AliJetFillCalTrkTrackMC(AliVEvent *fVEvt);
  virtual ~AliJetFillCalTrkTrackMC();
  
  // Setter
  void    SetHadCorrector(AliJetHadronCorrection* const corr) {fHadCorr = corr;}
  void    SetApplyMIPCorrection(Bool_t const val)             {fApplyMIPCorrection = val;}
  void    SetVEvent(AliVEvent* const evt)                     {fVEvt = evt;} 
  void    SetMCEvent(AliMCEvent* const mc)                    {fMCEvent = mc ;}

  // Getter
  Int_t   GetHadCorrection() const                            {return fApplyMIPCorrection;}

  // Other
  void    Exec(Option_t const * option);
  Bool_t  FillKine();
  // Fast Simulation
  Float_t SmearMomentum(Int_t ind, Float_t p);
  Bool_t  Efficiency(Float_t pt, Float_t eta, Float_t phi);

  // we have different cases
  // AOD reading -> MC from AOD
  // ESD reading -> MC from Kinematics
  // this has to match with our selection of input events
  enum {kTrackUndef = 0, kTrackKineAll, kTrackKineCharged, kTrackAODMCAll, kTrackAODMCCharged, kTrackAODMCChargedAcceptance};

 protected:
  AliJetHadronCorrection* fHadCorr;            // Pointer to Hadron Correction Object
  Bool_t                  fApplyMIPCorrection; // Apply MIP or not ? Exclusive with fApplyFractionHadronicCorrection
  AliVEvent*              fVEvt;               // Pointer to AliVEvent object
  AliMCEvent*             fMCEvent;            // Pointer to AliMCEvent object
          
 private:
  AliJetFillCalTrkTrackMC(const AliJetFillCalTrkTrackMC &det);
  AliJetFillCalTrkTrackMC &operator=(const AliJetFillCalTrkTrackMC &det);

  ClassDef(AliJetFillCalTrkTrackMC,1)          // Fill AliJetCalTrkTrack/TrackKine with MC track information
};

#endif
