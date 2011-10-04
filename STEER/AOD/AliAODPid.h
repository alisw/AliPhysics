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

class AliAODPid : public TObject {

 public:
  AliAODPid();
  virtual ~AliAODPid();
  AliAODPid(const AliAODPid& pid); 
  AliAODPid& operator=(const AliAODPid& pid);
  
  enum{kSPECIES=5, kTRDnPlanes=6};

 //setters
  void      SetITSsignal(Double_t its)                         {fITSsignal=its;}
  void      SetITSdEdxSamples(const Double_t s[4]);
  void      SetTPCsignal(Double_t tpc)                         {fTPCsignal=tpc;}
  void      SetTPCsignalN(UShort_t tpcN)                       {fTPCsignalN=tpcN;}
  void      SetTPCmomentum(Double_t tpcMom)                    {fTPCmomentum=tpcMom;}
  inline void  SetTRDsignal(Int_t nslices, const Double_t * const trdslices);  
  void      SetTRDmomentum(Int_t nplane, Float_t trdMom)       {fTRDmomentum[nplane]=trdMom;}
  inline void  SetTRDncls(UChar_t ncls, Int_t layer = -1);
  void      SetTRDntrackletsPID(UChar_t ntls) {fTRDntls = ntls;}
  void      SetTOFsignal(Double_t tof)                         {fTOFesdsignal=tof;}
  void      SetTOFpidResolution(Double_t tofPIDres[5]);
  void      SetIntegratedTimes(Double_t timeint[5]);
  void      SetHMPIDsignal(Double_t hmpid)                     {fHMPIDsignal=hmpid;}
  void      SetHMPIDprobs(Double_t hmpPid[5]);
  void      SetEMCALPosition(Double_t emcalpos[3]);
  void      SetEMCALMomentum(Double_t emcalmom[3]);

  Double_t  GetITSsignal()       const {return  fITSsignal;}
  void      GetITSdEdxSamples(Double_t *s) const;
  Double_t  GetITSdEdxSample(Int_t i) const {
    if(i>=0 && i<4) return fITSdEdxSamples[i];
    else return 0.;
  }
  Double_t  GetTPCsignal()       const {return  fTPCsignal;}
  UShort_t  GetTPCsignalN()      const {return  fTPCsignalN;}
  Double_t  GetTPCmomentum()     const {return  fTPCmomentum;}
  Int_t     GetTRDnSlices()      const {return  fTRDnSlices;}
  Double_t* GetTRDsignal()       const {return  fTRDslices;}
  const Float_t*  GetTRDmomentum() const {return  fTRDmomentum;}
  UChar_t   GetTRDncls(UChar_t layer) const { if(layer > 5) return 0; return fTRDncls[layer];}
  inline UChar_t GetTRDncls() const;
  UChar_t   GetTRDntrackletsPID() const {return fTRDntls;}
  Double_t  GetTOFsignal()       const {return  fTOFesdsignal;}
  Double_t  GetHMPIDsignal()     const {return  fHMPIDsignal;}
  void      GetHMPIDprobs(Double_t *p) const;

  void      GetIntegratedTimes(Double_t timeint[5])  const; 
  void      GetEMCALPosition  (Double_t emcalpos[3]) const;
  void      GetEMCALMomentum  (Double_t emcalmom[3]) const;
  void      GetTOFpidResolution (Double_t tofRes[5]) const;

 private :
  Double32_t  fITSsignal;        //[0.,0.,10] detector raw signal
  Double32_t  fITSdEdxSamples[4];//[0.,0.,10] ITS dE/dx samples
  Double32_t  fTPCsignal;        //[0.,0.,10] detector raw signal
  UShort_t    fTPCsignalN;       // number of points used for TPC dE/dx
  Double_t    fTPCmomentum;      // momentum at the inner wall of TPC;
  Int_t       fTRDnSlices;       // N slices used for PID in the TRD
  UChar_t     fTRDntls;          // number of tracklets used for PID calculation
  UChar_t     fTRDncls[6];       // number of clusters used for dE/dx calculation
  Double32_t* fTRDslices;        //[fTRDnSlices]
  Float_t     fTRDmomentum[6];   // momentum at the TRD layers
  Double32_t  fTOFesdsignal;     // TOF signal - t0 (T0 interaction time)
  Double32_t  fTOFpidResolution[5]; // TOF pid resolution for each mass hypotesys 
  Double32_t  fIntTime[5];       // track time hypothesis
  Double32_t  fHMPIDsignal;      // detector raw signal
  Double32_t  fHMPIDprobs[5];    // detector pid probabilities
  Double32_t  fEMCALPosition[3]; // global position of track
				 // extrapolated to EMCAL surface
  Double32_t  fEMCALMomentum[3]; // momentum of track
				 // extrapolated to EMCAL surface

  ClassDef(AliAODPid, 8);
};

//_____________________________________________________________
void AliAODPid::SetTRDsignal(Int_t nslices, const Double_t * const trdslices) {
  //
  // Set TRD dE/dx slices and the number of dE/dx slices per track
  //
  if(fTRDslices && fTRDnSlices != nslices) {
    delete [] fTRDslices; fTRDslices = NULL;
  };
  if(!fTRDslices) fTRDslices = new Double32_t[nslices];
  fTRDnSlices = nslices; 
  for(Int_t is = 0; is < fTRDnSlices; is++) fTRDslices[is] = trdslices[is];
}

//_____________________________________________________________
void AliAODPid::SetTRDncls(UChar_t ncls, Int_t layer) { 
  //
  // Set the number of clusters / tracklet
  // If no layer is specified the full number of clusters will be put in layer 0
  //
  if(layer > 5) return; 
  if(layer < 0) fTRDncls[0] = ncls;
  else fTRDncls[layer] = ncls;
}

//_____________________________________________________________
UChar_t AliAODPid::GetTRDncls() const{
  //
  // Get number of clusters per track
  // Calculated as sum of the number of clusters per tracklet
  //
  UChar_t ncls = 0;
  for(Int_t ily = 0; ily < 6; ily++) ncls += fTRDncls[ily];
  return ncls;
}
#endif
