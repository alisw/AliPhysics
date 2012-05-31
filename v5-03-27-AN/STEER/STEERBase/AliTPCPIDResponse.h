#ifndef ALITPCPIDRESPONSE_H
#define ALITPCPIDRESPONSE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------
//                    TPC PID class
// A very naive design... Should be made better by the detector experts...
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
// With many additions and modifications suggested by
//      Alexander Kalweit, GSI, alexander.philipp.kalweit@cern.ch
//      Dariusz Miskowiec, GSI, D.Miskowiec@gsi.de
//-------------------------------------------------------
#include <Rtypes.h>
#include <TObjArray.h>

#include "AliPID.h"

class AliTPCPIDResponse {
public:
  AliTPCPIDResponse();
  AliTPCPIDResponse(const Double_t *param);
  virtual ~AliTPCPIDResponse() {}
  void SetSigma(Float_t res0, Float_t resN2);
  void SetBetheBlochParameters(Double_t kp1,
                               Double_t kp2,
                               Double_t kp3,
                               Double_t kp4,
                               Double_t kp5
                               );
  void SetMip(Float_t mip) { fMIP = mip; } // Set overall normalisation; mean dE/dx for MIP
  Double_t Bethe(Double_t bg) const;
  void SetUseDatabase(Bool_t useDatabase) { fUseDatabase = useDatabase;}
  
  void SetResponseFunction(AliPID::EParticleType type, TObject * const o) { fResponseFunctions.AddAt(o,(Int_t)type); }
  const TObject * GetResponseFunction(AliPID::EParticleType type) { return fResponseFunctions.At((Int_t)type); }
  
  Double_t GetExpectedSignal(const Float_t mom,
                     AliPID::EParticleType n=AliPID::kKaon) const;
  Double_t GetExpectedSigma(const Float_t mom, const Int_t nPoints,
                     AliPID::EParticleType n=AliPID::kKaon) const;
  Float_t  GetNumberOfSigmas(const Float_t mom, const Float_t dEdx, 
			     const Int_t nPoints,
                     AliPID::EParticleType n=AliPID::kKaon) const {

    Double_t bethe=GetExpectedSignal(mom,n);
    Double_t sigma=GetExpectedSigma(mom,nPoints,n);
    return (dEdx-bethe)/sigma;
  }

  Double_t GetMIP() const { return fMIP;} 
  Float_t  GetRes0()  const { return fRes0;  }
  Float_t  GetResN2() const { return fResN2; }

private:
  Float_t fMIP;          // dEdx for MIP
  Float_t fRes0;         // relative dEdx resolution  rel sigma = fRes0*sqrt(1+fResN2/npoint)
  Float_t fResN2;        // relative Npoint dependence rel  sigma = fRes0*sqrt(1+fResN2/npoint)

  Double_t fKp1;   // Parameters
  Double_t fKp2;   //    of
  Double_t fKp3;   // the ALEPH
  Double_t fKp4;   // Bethe-Bloch
  Double_t fKp5;   // formula

  Bool_t fUseDatabase; // flag if fine-tuned database-response or simple ALEPH BB should be used
  TObjArray fResponseFunctions; //! ObjArray of response functions individually for each particle

  ClassDef(AliTPCPIDResponse,3)   // TPC PID class
};

#endif


