/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
// PID Response class for the TRD detector
// Based on 1D Likelihood approach
// For further information see implementation file
//
#ifndef ALITRDPIDRESPONSE_H
#define ALITRDPIDRESPONSE_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#ifndef ALIPID_H
#include "AliPID.h"
#endif

class TObjArray;
class AliVTrack;
class AliTRDPIDResponse : public TObject {
public:
 enum ETRDPIDResponseStatus {
   kIsOwner = BIT(14)
 };
 enum ETRDPIDResponseDef {
   kNlayer = 6
  ,kNPBins = 6
 };
 AliTRDPIDResponse(const Char_t *filename = NULL);
 AliTRDPIDResponse(const AliTRDPIDResponse &ref);
 AliTRDPIDResponse& operator=(const AliTRDPIDResponse &);
 ~AliTRDPIDResponse();

 Int_t     GetPIDmethod() const { return fPIDmethod;}
 Bool_t    GetResponse(Int_t n, Double_t *dedx, Float_t *p, Double_t prob[AliPID::kSPECIES], Bool_t kNorm=kTRUE) const;

 Bool_t    IsOwner() const {return TestBit(kIsOwner);}

 void      SetOwner();
 void      SetPIDmethod(Int_t m) {fPIDmethod=m;}
 void      SetGainNormalisationFactor(Double_t gainFactor) { fGainNormalisationFactor = gainFactor; }

private:
 Bool_t    CookdEdx(Double_t *in, Double_t *out) const;
 Int_t     GetLowerMomentumBin(Double_t p) const;
 Double_t  GetProbabilitySingleLayer(Int_t species, Double_t plocal, Double_t dEdx) const;
 Bool_t    Load(const Char_t *filename = NULL);

 static const Double_t fgkPBins[kNPBins];
 TObjArray *fReferences; // Container for reference distributions
 Int_t     fMapRefHists[AliPID::kSPECIES][kNPBins];     
                         // Map for the position of a given historgam in the container 
 Double_t  fGainNormalisationFactor;  // Gain normalisation factor
 UChar_t   fPIDmethod;   // PID method selector  

 ClassDef(AliTRDPIDResponse, 2)    // Tool for TRD PID
};
#endif

