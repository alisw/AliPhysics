#ifndef AliAODPid_H
#define AliAODPid_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD Pid object for additional pid information
//     Author: Annalisa Mastroserio, CERN
//-------------------------------------------------------------------------

#include <TObject.h>
class AliESDtrack;

class AliAODPid : public TObject {

 public:
  AliAODPid();
  virtual ~AliAODPid();
  AliAODPid(const AliAODPid& pid); 
  AliAODPid& operator=(const AliAODPid& pid);
  
  enum{kSPECIES=5, kTRDnPlanes=6};

  void SetDetectorRawSignals(AliESDtrack *track, Double_t timezero); 

  Double_t  GetITSsignal()       {return  fITSsignal;}
  Double_t  GetTPCsignal()       {return  fTPCsignal;}
  Int_t     GetTRDnSlices()      {return  fTRDnSlices;}
  Double_t* GetTRDsignal()       {return  fTRDslices;}
  Double_t  GetTOFsignal()       {return  fTOFesdsignal;} 
  void      GetIntegratedTimes(Double_t timeint[5]); 
  Double_t  GetHMPIDsignal()     {return  fHMPIDsignal;}

 private :
  Double32_t fITSsignal;      //[0.,0.,10] detector raw signal
  Double32_t fTPCsignal;      //[0.,0.,10] detector raw signal
  Int_t      fTRDnSlices;     //N slices used for PID in the TRD
  Double32_t* fTRDslices;     //[fTRDnSlices]
  Double32_t fTOFesdsignal;   //TOF signal - t0 (T0 interaction time)
  Double32_t fIntTime[5];     //track time hypothesis
  Double32_t fHMPIDsignal;    //detector raw signal

  ClassDef(AliAODPid,1);
};

#endif
