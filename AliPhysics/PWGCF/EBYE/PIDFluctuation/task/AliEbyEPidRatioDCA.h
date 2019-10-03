#ifndef ALIEBYEPIDRSTIODCA_H
#define ALIEBYEPIDRSTIODCA_H

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
#include "AliEbyEPidRatioBase.h"


class AliEbyEPidRatioDCA : public AliEbyEPidRatioBase {

 public:

  AliEbyEPidRatioDCA();
  virtual ~AliEbyEPidRatioDCA();

  virtual void Process();
  void SetESDTrackCutsBkg(AliESDtrackCuts *p) {fESDTrackCutsBkg = p;}
  THnSparseD* GetHnDCA() {return fHnDCA;}
  
 private:

  AliEbyEPidRatioDCA(const AliEbyEPidRatioDCA&); // not implemented
  AliEbyEPidRatioDCA& operator=(const AliEbyEPidRatioDCA&); // not implemented

  virtual void CreateHistograms();
  Int_t GetContIdxTrack(Int_t label, Int_t sign, Int_t gPdgCode);
  AliESDtrackCuts    *fESDTrackCutsBkg;       //! ESD cuts  
  THnSparseD         *fHnDCA;                 //  THnSparseD contamination DCA
  
  ClassDef(AliEbyEPidRatioDCA, 1);
};

#endif
