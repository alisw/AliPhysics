#ifndef ALIEBYEPIDRATIOBASE_H
#define ALIEBYEPIDRATIOBASE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//=========================================================================//
//             AliEbyE Analysis for Particle Ratio Fluctuation             //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
//                  Date: Wed Jul  9 18:38:30 CEST 2014                    //
//          New approch to find particle ratio to reduce memory            //
//                             (Test Only)                                 //
//        Copied from NetParticle Classes
//        Origin: Authors: Jochen Thaeder <jochen@thaeder.de>
//                         Michael Weber <m.weber@cern.ch>
//=========================================================================//

#include "THnSparse.h"
#include "AliEbyEPidRatioHelper.h"

class AliESDEvent;
class AliESDInputHandler;
class AliMCEvent;
class AliStack;
class AliAODInputHandler;

class AliEbyEPidRatioBase : public TNamed {

 public:

  AliEbyEPidRatioBase();
  AliEbyEPidRatioBase(const Char_t* name, const Char_t* title);
  virtual ~AliEbyEPidRatioBase();

  void Initialize(AliEbyEPidRatioHelper* helper, AliESDtrackCuts* cuts = NULL);
  Int_t SetupEvent();
  void ResetEvent();
  virtual void Process() = 0;


 private:

  AliEbyEPidRatioBase(const AliEbyEPidRatioBase&); // not implemented
  AliEbyEPidRatioBase& operator=(const AliEbyEPidRatioBase&); // not implemented

 protected:
  virtual void Init() {};
  virtual void CreateHistograms() {};
  virtual void Reset() {};
  virtual Int_t Setup() { return 0;};
  AliEbyEPidRatioHelper *fHelper;             //! Ptr to helper class
  AliESDEvent        *fESD;                   //! ESD object
  AliESDtrackCuts    *fESDTrackCuts;          //! ESD cuts  
  AliAODEvent        *fAOD;                   //! AOD object
  TClonesArray       *fArrayMC;               //! array of MC particles
  Int_t               fAODtrackCutBit;        //  Track filter bit for AOD tracks
  Bool_t              fIsMC;                  //  Is MC event
  AliMCEvent         *fMCEvent;               //! Ptr to MC event
  AliStack           *fStack;                 //! Ptr to stack
  Float_t             fCentralityBin;         //  Centrality of current event  
  Int_t               fNTracks;               //  N Tracks in the current event
  Bool_t              fIsRatio;               //  Is MC event
  Bool_t              fIsPtBin;               //  Is Pt Bin event
  Bool_t              fIsDetectorWise;        //  Is Detector Wise
  ClassDef(AliEbyEPidRatioBase, 1);
};

#endif
