#ifndef ALIJETFILLCALTRKEVENT_H
#define ALIJETFILLCALTRKEVENT_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
/* $Id$ */

//--------------------------------------------------
// Filling of CalTrkEvent objects in the reader task
//
// Author: magali.estienne@subatech.in2p3.fr
//         alexandre.shabetai@cern.ch
//-------------------------------------------------

class AliJetReaderHeader;
class AliEMCALGeometry;
//class AliEMCALRecoUtils;
class AliVEvent;
class AliMCEvent;
class AliJetHadronCorrection;

#include <TMath.h>
#include "AliJetCalTrk.h"

class AliJetFillCalTrkEvent 
{
 public: 
  AliJetFillCalTrkEvent();
  virtual ~AliJetFillCalTrkEvent();
  AliJetFillCalTrkEvent(const AliJetFillCalTrkEvent &det);
  AliJetFillCalTrkEvent &operator=(const AliJetFillCalTrkEvent &det);
  
  // Setter
  virtual void SetReaderHeader(AliJetReaderHeader* const readerHeader) {fReaderHeader = readerHeader;}
  virtual void SetGeom(AliEMCALGeometry* const geom)                   {fGeom = geom;}
  virtual void SetCalTrkEvent(AliJetCalTrkEvent* caltrkevt) {fCalTrkEvent = caltrkevt;}
  virtual void SetHadCorrector(AliJetHadronCorrection* /*corr*/)  {;}
  virtual void SetApplyMIPCorrection(Bool_t /*val*/)              {;}
  virtual void SetVEvent(AliVEvent */*aod*/)                       {;}
  virtual void SetMCEvent(AliMCEvent */*MC*/)                      {;}
  //virtual void SetEMCALRecoUtils(AliEMCALRecoUtils */*ru*/)       {;}
  virtual void SetApplyElectronCorrection(Int_t /*flag*/)         {;}
  virtual void SetApplyFractionHadronicCorrection(Bool_t /*val*/) {;}
  virtual void SetFractionHadronicCorrection(Double_t /*val*/)    {;}

  // Getter
  virtual AliJetCalTrkEvent* GetCalTrkEvent() const {return fCalTrkEvent;}

  // Other
  virtual void          Exec(const Option_t  */*option*/) {;}
  virtual Float_t       EtaToTheta(Float_t arg);

 protected:
  Int_t                 fOpt;             // Detector to be used for jet reconstruction
  Int_t                 fDebug;           // Debug option
  AliJetReaderHeader   *fReaderHeader;    // ReaderHeader
  AliJetCalTrkEvent    *fCalTrkEvent;     // CalTrk event

  AliEMCALGeometry     *fGeom;            // Define EMCal geometry

 private:

  ClassDef(AliJetFillCalTrkEvent,1) // Fill AliJetFillCalTrkEvent with tpc and/or emcal information
};

#endif
