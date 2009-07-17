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

#ifndef ALIRSNPIDDEFESD_H
#define ALIRSNPIDDEFESD_H

class AliESDtrack;

class AliRsnPIDDefESD : public TObject
{
  public:

    enum EDetector {
      kITS = 0,
      kTPC,
      kTRD,
      kTOF,
      kHMPID,
      kDetectors
    };

    enum EScheme {
      kSchemeESD = 0,
      kSchemeITS,
      kSchemeTPC,
      kSchemeTOF,
      kSchemeITSandTPC,
      kSchemeITSandTOF,
      kSchemeTPCandTOF,
      kSchemeITSandTPCandTOF,
      kSchemeITSandTPCandTOFwithSP,
      kSchemeITSandTPCorTOFwithSP,
      kSchemeLastPIDType
    };

    AliRsnPIDDefESD();
    AliRsnPIDDefESD(const AliRsnPIDDefESD& copy);
    virtual ~AliRsnPIDDefESD() {}

    void        UseESDWeights() {fUseESDWeights = kTRUE;};
    void        NoESDWeights() {fUseESDWeights = kFALSE;}
    void        SetScheme(EScheme scheme, Double_t divValue = 1.0);
    void        IncludeDet(EDetector det) {if (CheckBounds(det)) fUseDet[det] = kTRUE;}
    void        ExcludeDet(EDetector det) {if (CheckBounds(det)) fUseDet[det] = kFALSE;}
    void        IncludeAll() { Int_t det; for (det = 0; det < kDetectors; det++) fUseDet[det] = kTRUE; }
    void        ExcludeAll() { Int_t det; for (det = 0; det < kDetectors; det++) fUseDet[det] = kFALSE; }
    void        SetDivValue(EDetector det, Double_t value,Bool_t userHigher=kTRUE);
    void        ComputeWeights(AliESDtrack *track, Double_t *weights);
    void        PrintStatus();
    const char* DetName(EDetector det) const;
    const char* SchemeName();

  private:

    Bool_t   CheckBounds(EDetector det) const { return (det >= kITS && det < kDetectors); }
    Bool_t   CheckDivValue(EDetector det,Double_t value);

    Bool_t   fUseESDWeights;          // with this flag, ESD weights are used
    Bool_t   fUseDet[kDetectors];     // flag to switch off info from a detector
    Double_t fDivValue[kDetectors];   // division value for detector weight acceptance
    Double_t fUseHigher[kDetectors];  // accepted higher or lower then div value

    ClassDef(AliRsnPIDDefESD,1);
};

#endif
