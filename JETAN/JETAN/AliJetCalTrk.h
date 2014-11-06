#ifndef ALIJETCALTRK_H
#define ALIJETCALTRK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *  *  * See cxx source for full Copyright notice     */

/* $Id$ */

//------------------------------------------------------
// CalTrk is used to sore Tracks and CaloCells information
//                    
// Author: alexandre.shabetai@cern.ch & magali.estienne@subatech.in2p3.fr 
//------------------------------------------------------

#include <Riostream.h> 
#include <TObject.h>
#include <TRef.h>
#include <TClonesArray.h>

#include "AliVParticle.h"
#include "AliVTrack.h" 

class AliVCaloCells;

class AliJetCalTrkTrack : public TObject
{
 public:
  AliJetCalTrkTrack();
  AliJetCalTrkTrack(AliVParticle* track, Bool_t CutFlag, Bool_t SignalFlag, Float_t PtCorr = 1.);
  AliJetCalTrkTrack(AliVTrack* track, Bool_t CutFlag, Bool_t SignalFlag, Float_t PtCorr = 1.);
  virtual  ~AliJetCalTrkTrack() {;}
  
  // Setter
  void            SetCutFlag(Bool_t cutFlag)       {fCalTrkCutFlag = cutFlag;}
  void            SetSignalFlag(Bool_t signalFlag) {fCalTrkSignalFlag = signalFlag;}
  void            SetPtCorr(Float_t ptcorr)        {fCalTrkPtCorr = ptcorr;}

  // Getter
  AliVParticle*   GetParticle() const    {return ((AliVParticle*)fCalTrkTrackRef.GetObject());}
  AliVTrack*      GetTrack() const       {return ((AliVTrack*)fCalTrkTrackRef.GetObject());}
  TRef            GetTrackRef() const    {return fCalTrkTrackRef;}
  TObject*        GetTrackObject() const {return fCalTrkTrackRef.GetObject();}
  Bool_t          GetCutFlag() const     {return fCalTrkCutFlag;}
  Bool_t          GetSignalFlag() const  {return fCalTrkSignalFlag;}
  Float_t         GetPtCorr() const      {return fCalTrkPtCorr;}
  Float_t         GetEta() const         {return GetParticle()->Eta();}
  Float_t         GetPhi() const         {return GetParticle()->Phi();}
  virtual Float_t GetE() const           {return GetParticle()->E();}
  virtual Float_t GetPt() const          {return GetParticle()->Pt();}
  virtual Float_t GetPx() const          {return GetParticle()->Px();}
  virtual Float_t GetPy() const          {return GetParticle()->Py();}
  virtual Float_t GetPz() const          {return GetParticle()->Pz();}
  virtual Float_t GetP() const           {return GetParticle()->P();}
  Int_t           GetID() const          {return GetTrack()->GetID();}
  Float_t         GetM() const           {return GetParticle()->M();}
  void            Print(const Option_t* /*option*/) const;

  virtual void    Clear(Option_t* /*option = ""*/);

 private:
  AliJetCalTrkTrack(const AliJetCalTrkTrack& rCalTrk);
  AliJetCalTrkTrack& operator = (const AliJetCalTrkTrack& rhs);
 
  TRef            fCalTrkTrackRef;     //! Reference to a track
  Bool_t          fCalTrkCutFlag;      //  Cut flag of the track in the tpc
  Bool_t          fCalTrkSignalFlag;   //  Signal flag of the track
  Float_t         fCalTrkPtCorr;       //  Correction factor

  ClassDef(AliJetCalTrkTrack,1) // Implementation of AliJetCalTrkTrack

};

//-------------------------------------------------------------
class AliJetCalTrkTrackKine : public AliJetCalTrkTrack
{
 public:
  AliJetCalTrkTrackKine();
  AliJetCalTrkTrackKine(AliVParticle* track, Bool_t CutFlag, Bool_t SignalFlag, Float_t ptReso);
  virtual  ~AliJetCalTrkTrackKine() {;}

  Float_t        GetPtReso() const {return fCalTrkPtReso;}
  Float_t        GetE()  const {return fCalTrkTrackE;}
  Float_t        GetPt() const {return fCalTrkTrackPt;}
  Float_t        GetPx() const {return fCalTrkTrackPx;}
  Float_t        GetPy() const {return fCalTrkTrackPy;}
  Float_t        GetPz() const {return fCalTrkTrackPz;}
  Float_t        GetP()  const {return fCalTrkTrackP;}
  Float_t        CalcPx() {fCalTrkTrackPx = GetPtReso()*GetParticle()->Px(); return fCalTrkTrackPx;}
  Float_t        CalcPy() {fCalTrkTrackPy = GetPtReso()*GetParticle()->Py(); return fCalTrkTrackPy;}
  Float_t        CalcPz() {fCalTrkTrackPz = GetPtReso()*GetParticle()->Pz(); return fCalTrkTrackPz;}
  Float_t        CalcP();
  Float_t        CalcPt();
  Float_t        CalcE();
  void           Clear(Option_t* option = "");

 private:
  Float_t        fCalTrkPtReso;       // Pt resolution rescaling (AliJetFillCalTrkTrackKine::FillKine())
  Float_t        fCalTrkTrackE;       // Particle energy
  Float_t        fCalTrkTrackPt;      // Particle Pt
  Float_t        fCalTrkTrackP;       // Particle P
  Float_t        fCalTrkTrackPx;      // Particle Px
  Float_t        fCalTrkTrackPy;      // Particle Py
  Float_t        fCalTrkTrackPz;      // Particle Pz

  ClassDef(AliJetCalTrkTrackKine,1) // Implementation of AliJetCalTrkTrackKine
};

//-------------------------------------------------------------
class AliJetCalTrkEvent : public TObject
{  
 public:
  AliJetCalTrkEvent();                        //default constructor
  AliJetCalTrkEvent(Short_t opt,Bool_t kine,Bool_t kIsHighMult=kFALSE); // constructor 2
  virtual                ~AliJetCalTrkEvent();
  AliJetCalTrkEvent(const AliJetCalTrkEvent& rCalTrkEvent);
  AliJetCalTrkTrack*     AddCalTrkTrack(AliVTrack* track, Bool_t cutFlag, Bool_t signalFlag, Float_t ptCorr = 1.);
  AliJetCalTrkTrack*     AddCalTrkTrack(AliVParticle* track, Bool_t cutFlag, Bool_t signalFlag, Float_t ptCorr = 1.);
  AliJetCalTrkTrackKine* AddCalTrkTrackKine(AliVParticle* track, Bool_t cutFlag, Bool_t signalFlag, Float_t ptReso = 1.);

  AliJetCalTrkTrack*     GetCalTrkTrack(Int_t i);
  Int_t                  GetNCalTrkTracks() const {return fNJetCalTrkTrack;}

  void                   Clear(Option_t* option = ""); 
  void                   Print(const Option_t* = "") const;
  
 private:
  AliJetCalTrkEvent& operator = (const AliJetCalTrkEvent& rhs);
  TClonesArray*          fJetCalTrkTrack;  //! Array of Tracks
  TClonesArray*          fJetCalTrkCell;   //! Array of Cells
  Int_t	                 fNJetCalTrkTrack; //  Number of tracks 

  ClassDef(AliJetCalTrkEvent,1) // Implementation of AliJetCalTrkEvent

};

#endif

