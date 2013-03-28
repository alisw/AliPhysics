#ifndef ALIAODPIDUTIL_H
#define ALIAODPIDUTIL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAODpidUtil.h 38493 2010-01-26 16:33:03Z hristov $ */

//-------------------------------------------------------
//                    Combined PID class
//                    for the AOD class
//   Origin: Rosa Romita, GSI, r.romita@gsi.de 
//   Modified: Jens Wiechula, Uni Tuebingen, jens.wiechula@cern.ch
//   Modified: Pietro Antonioli, INFN BO, pietro.antonioli@bo.infn.it
//-------------------------------------------------------
#include <Rtypes.h>
#include <TMatrixD.h>
#include <AliLog.h>
#include "AliAODEvent.h" // Needed for inline functions
#include "AliAODTrack.h" // Needed for inline functions
#include "AliAODPid.h" // Needed for inline functions
#include "AliTOFHeader.h" //Needed for inline functions
//#include "HMPID/AliHMPID.h"

#include "AliPIDResponse.h"

class AliAODEvent;
class AliVParticle;

class AliAODpidUtil : public AliPIDResponse  {
public:
  //TODO: isMC???
  AliAODpidUtil(Bool_t isMC = kFALSE): AliPIDResponse(isMC) {;}
  virtual ~AliAODpidUtil() {;}


  Float_t GetTPCsignalTunedOnData(const AliVTrack *t) const;
  Float_t GetTOFsignalTunedOnData(const AliVTrack *t) const;

protected:
  virtual Float_t GetSignalDeltaTOFold(const AliVParticle *track, AliPID::EParticleType type, Bool_t ratio=kFALSE) const;
  virtual Float_t GetNumberOfSigmasTOFold(const AliVParticle *vtrack, AliPID::EParticleType type) const;
  
private:
  
  ClassDef(AliAODpidUtil,3)  // PID calculation class
};


#endif


