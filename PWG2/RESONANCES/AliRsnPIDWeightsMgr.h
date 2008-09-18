/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

/* $Id: AliRsnPIDWeightsMgr.h,v 1.5 2007/02/21 14:33:25 pulvir Exp $ */

//-------------------------------------------------------------------------
//                      Class AliRsnPIDWeightsMgr
//  Simple collection of reconstructed tracks, selected from an ESD event
//
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#ifndef AliRsnPIDWeightsMgr_H
#define AliRsnPIDWeightsMgr_H

#include "AliRsnPID.h"

class AliRsnPIDWeightsMgr : public TObject
{
  public:

    // detectors for customizing PID weights
    enum EDetector
    {
      kITS,
      kTPC,
      kTRD,
      kTOF,
      kHMPID,
      kDetectors
    };

    AliRsnPIDWeightsMgr();
    virtual ~AliRsnPIDWeightsMgr() {}

    void      UseDetector(EDetector det, Bool_t use) {if (CheckBounds(det)) fUseDet[det] = use;}
    void      SetDetectorWeights(EDetector det, Double_t *weights);
    void      SetAcceptanceRange(EDetector det, Double_t ptmin, Double_t ptmax);
    Double_t  GetWeight(AliRsnPID::EType type, Double_t pt);
    Double_t* GetWeightArray(EDetector det) {if (CheckBounds(det)) return fWeights[det]; else return 0x0;}

  private:

    Bool_t   CheckBounds(EDetector det);

    Bool_t   fUseDet[kDetectors];                        // flag to switch off info from a detector
    Double_t fDetPtMin[kDetectors];                      // min value for detector weight acceptance
    Double_t fDetPtMax[kDetectors];                      // max value for detector weight acceptance
    Double_t fWeights[kDetectors][AliRsnPID::kSpecies];  // PID weights of a single detector

    ClassDef(AliRsnPIDWeightsMgr,1);
};

#endif
