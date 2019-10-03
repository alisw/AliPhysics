#ifndef ALIJETFILLCALTRKTRACK_H
#define ALIJETFILLCALTRKTRACK_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
/* $Id$ */

//--------------------------------------------------
// Filling of CalTrkTrack objects in the reader task
//
// Author: magali.estienne@subatech.in2p3.fr
//         alexandre.shabetai@cern.ch
//-------------------------------------------------

#include "AliJetFillCalTrkEvent.h"

class AliJetFillCalTrkTrack : public AliJetFillCalTrkEvent
{
 public: 
  AliJetFillCalTrkTrack();
  AliJetFillCalTrkTrack(AliVEvent *fVEvt);
  virtual ~AliJetFillCalTrkTrack();
  
  // Setter
  void  SetHadCorrector(AliJetHadronCorrection* const corr) {fHadCorr = corr;}
  void  SetApplyMIPCorrection(Bool_t const val)             {fApplyMIPCorrection = val;}
  void  SetVEvent(AliVEvent* const evt)                     {fVEvt = evt;} 

  // Getter
  Int_t GetHadCorrection()  const {return fApplyMIPCorrection;}

  // Other
  void  Exec(Option_t const * option);

  // we have different cases
  // AOD reading -> MC from AOD
  // ESD reading -> MC from Kinematics
  // this has to match with our selection of input events
  enum {kTrackUndef = 0, kTrackESD, kTrackAOD, kTrackAODextra, kTrackAODextraonly};

 protected:
  AliJetHadronCorrection* fHadCorr;            // Pointer to Hadron Correction Object
  Bool_t                  fApplyMIPCorrection; // Apply MIP or not ? Exclusive with fApplyFractionHadronicCorrection
  AliVEvent*              fVEvt;               // Pointer to AliVEvent object

 private:
  AliJetFillCalTrkTrack(const AliJetFillCalTrkTrack &det);
  AliJetFillCalTrkTrack &operator=(const AliJetFillCalTrkTrack &det);

  ClassDef(AliJetFillCalTrkTrack,1)            // Fill AliJetCalTrkTrack with track information

};

#endif
