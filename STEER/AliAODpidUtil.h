#ifndef ALIAODPIDUTIL_H
#define ALIAODPIDUTIL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAODpidUtil.h 38493 2010-01-26 16:33:03Z hristov $ */

//-------------------------------------------------------
//                    Combined PID class
//                    for the AOD class
//   Origin: Rosa Romita, GSI, r.romita@gsi.de 
//-------------------------------------------------------
#include <Rtypes.h>
#include <TMatrixD.h>
#include "AliAODTrack.h" // Needed for inline functions
#include "AliAODPid.h" // Needed for inline functions
#include "AliTPCPIDResponse.h"
#include "AliITSPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliTRDPIDResponse.h"
//#include "HMPID/AliHMPID.h"

class AliAODEvent;

class AliAODpidUtil {
public:

  AliAODpidUtil(Bool_t isMC = kFALSE): fRange(5.), fTPCResponse(), fITSResponse(isMC), fTOFResponse(), fTRDResponse() {;}
  virtual ~AliAODpidUtil() {;}


  Int_t MakePID(AliAODTrack *track,Double_t *p) const;
  void MakeTPCPID(AliAODTrack *track,Double_t *p) const;
  void MakeITSPID(AliAODTrack *track,Double_t *p) const;
  void MakeTOFPID(AliAODTrack *track,Double_t *p) const;
  //  void MakeHMPIDPID(AliESDtrack *track);
  void MakeTRDPID(AliAODTrack *track,Double_t *p) const;

  Float_t NumberOfSigmasTPC(const AliAODTrack *track, AliPID::EParticleType type) const;
  Float_t NumberOfSigmasTOF(const AliAODTrack *track, AliPID::EParticleType type) const;
  Float_t NumberOfSigmasITS(const AliAODTrack *track, AliPID::EParticleType type) const;

  AliITSPIDResponse &GetITSResponse() {return fITSResponse;}
  AliTPCPIDResponse &GetTPCResponse() {return fTPCResponse;}
  AliTOFPIDResponse &GetTOFResponse() {return fTOFResponse;}

private:
  Float_t           fRange;          // nSigma max in likelihood
  AliTPCPIDResponse fTPCResponse;
  AliITSPIDResponse fITSResponse;
  AliTOFPIDResponse fTOFResponse;
  AliTRDPIDResponse fTRDResponse;
  
  ClassDef(AliAODpidUtil,1)  // PID calculation class
};

#endif


