/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

/* $Id: AliRsnPIDScheme.h,v 1.5 2007/02/21 14:33:25 pulvir Exp $ */

//-------------------------------------------------------------------------
//                      Class AliRsnPIDDef
//  Simple collection of reconstructed tracks, selected from an ESD event
//
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
//-------------------------------------------------------------------------

#ifndef ALIRSNPIDDEF_H
#define ALIRSNPIDDEF_H

class AliRsnPIDDef : public TObject
{
  public:

    enum EDetector
    {
      kITS,
      kTPC,
      kTRD,
      kTOF,
      kHMPID,
      kDetectors
    };

    AliRsnPIDDef();
    virtual ~AliRsnPIDDef() {}

    void     IncludeDet(EDetector det) {if (CheckBounds(det)) fUseDet[det] = kTRUE;}
    void     ExcludeDet(EDetector det) {if (CheckBounds(det)) fUseDet[det] = kTRUE;}
    void     SetAcceptanceRange(EDetector det, Double_t ptmin, Double_t ptmax);

  private:

    Bool_t   CheckBounds(EDetector det) { return (det >= kITS && det < kDetectors); }

    Bool_t   fUseESD;                 // with this flag, ESD weights are used
    Bool_t   fUseDet[kDetectors];     // flag to switch off info from a detector
    Double_t fDetPtMin[kDetectors];   // min value for detector weight acceptance
    Double_t fDetPtMax[kDetectors];   // max value for detector weight acceptance

    ClassDef(AliRsnPIDScheme,1);
};

#endif
