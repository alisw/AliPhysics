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
class AliTRDPIDParams;
class AliTRDPIDReference;
class AliTRDPIDResponse : public TObject {
  public:
    enum ETRDPIDResponseStatus {
      kIsOwner = BIT(14)
    };
    enum ETRDPIDResponseDef {
      kNlayer = 6
      ,kNPBins = 6
    };
    enum ETRDPIDMethod {
      kNN   = 0,
      kLQ2D = 1,
      kLQ1D = 2
    };
    enum ETRDNslices {
      kNslicesLQ1D = 1,
      kNslicesLQ2D = 2,
      kNslicesNN = 7
    };
    AliTRDPIDResponse();
    AliTRDPIDResponse(const AliTRDPIDResponse &ref);
    AliTRDPIDResponse& operator=(const AliTRDPIDResponse &ref);
    ~AliTRDPIDResponse();
    
    ETRDPIDMethod     GetPIDmethod() const { return fPIDmethod;}
    Bool_t    GetResponse(Int_t n, const Double_t * const dedx, const Float_t * const p, Double_t prob[AliPID::kSPECIES], Bool_t kNorm=kTRUE) const;
    inline ETRDNslices  GetNumberOfSlices() const;
    
    Bool_t    IsOwner() const {return TestBit(kIsOwner);}
    
    void      SetOwner();
    void      SetPIDmethod(ETRDPIDMethod m) {fPIDmethod=m;}
    void      SetGainNormalisationFactor(Double_t gainFactor) { fGainNormalisationFactor = gainFactor; }
    void      SetPIDParams(const AliTRDPIDParams * params) { fkPIDParams = params; }

    Bool_t    Load(const Char_t *filename = NULL, const Char_t *refName = "RefTRDLQ1D");
    Bool_t    Load(const AliTRDPIDReference *ref) { fkPIDReference = ref; return kTRUE; }
  
    Bool_t    IdentifiedAsElectron(Int_t nTracklets, const Double_t *like, Double_t p, Double_t level) const;
  
  private:
    Bool_t    CookdEdx(Int_t nSlice, const Double_t * const in, Double_t *out) const;
    Double_t  GetProbabilitySingleLayer(Int_t species, Double_t plocal, Double_t dEdx) const;
    
    const AliTRDPIDReference *fkPIDReference;   // PID References
    const AliTRDPIDParams *fkPIDParams;         // PID Params
    Double_t  fGainNormalisationFactor;         // Gain normalisation factor
    ETRDPIDMethod   fPIDmethod;                 // PID method selector
      
  
  ClassDef(AliTRDPIDResponse, 3)    // Tool for TRD PID
};

AliTRDPIDResponse::ETRDNslices AliTRDPIDResponse::GetNumberOfSlices() const {
  // Get the current number of slices
  ETRDNslices slices = kNslicesLQ1D;
  switch(fPIDmethod){
    case kLQ1D: slices = kNslicesLQ1D; break;
    case kLQ2D: slices = kNslicesLQ2D; break;
    case kNN:   slices = kNslicesNN; break;
  };
  return slices;
}
#endif
