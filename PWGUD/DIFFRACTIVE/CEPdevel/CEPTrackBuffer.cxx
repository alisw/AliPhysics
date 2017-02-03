// ////////////////////////////////////////////////////////////////////////////
//
// CEP track buffer
//
// structure to hold track information
//
// Authors:
// P. Buehler, paul.buehler@oeaw.ac.at       27.06.2016
//
//
// ----------------------------------------------------------------------------
#include "CEPTrackBuffer.h"

ClassImp(CEPTrackBuffer)

// ----------------------------------------------------------------------------
CEPTrackBuffer::CEPTrackBuffer()
  : TObject()
  , fTrackStatus(AliCEPBase::kTTUnknown)
  , fChargeSign(AliCEPBase::kdumval)
  , fTPCnclsS(AliCEPBase::kdumval)
  , fZv(-999.9)
  , fMomentum(TVector3(-999.9,-999.9,-999.9))
  , fPID(AliCEPBase::kdumval)
  , fPIDTPCStatus(AliCEPBase::kdumval)
  , fPIDTPCSignal(-999.9)
  , fPIDTOFStatus(AliCEPBase::kdumval)
  , fPIDTOFSignal(-999.9)
  , fPIDBayesStatus(AliCEPBase::kdumval)
  , fMCPID(AliCEPBase::kdumval)
  , fMCMass(-999.9)
  , fMCMomentum(TVector3(-999.9,-999.9,-999.9))
{

  this->Reset();

}

// ----------------------------------------------------------------------------
CEPTrackBuffer::~CEPTrackBuffer()
{
  
}

// ----------------------------------------------------------------------------
void CEPTrackBuffer::Reset()
{

  // general information
  fTrackStatus = AliCEPBase::kTTUnknown;
  fChargeSign  = AliCEPBase::kdumval;
  fTPCnclsS    = AliCEPBase::kdumval;
  fZv          = -999.9;
  fMomentum    = TVector3(-999.9,-999.9,-999.9);
  
  // PID information
  fPID = AliCEPBase::kdumval;
  
  // ... from TPC
  fPIDTPCStatus = AliCEPBase::kdumval;
  fPIDTPCSignal = AliCEPBase::kdumval;
  for(Int_t ii=0;ii<AliPID::kSPECIES;ii++) {
    fPIDTPCnSigma[ii]     = AliCEPBase::kdumval;
    fPIDTPCnSigmaProb[ii] = AliCEPBase::kdumval;
  }
  
  // ... from TOF
  fPIDTOFStatus = AliCEPBase::kdumval;
  fPIDTOFSignal = AliCEPBase::kdumval;
  for(Int_t ii=0;ii<AliPID::kSPECIES;ii++) {
    fPIDTOFnSigma[ii]     = AliCEPBase::kdumval;
    fPIDTOFnSigmaProb[ii] = AliCEPBase::kdumval;
  }
  
  // ... Bayes
  fPIDBayesStatus = AliCEPBase::kdumval;
  for(Int_t ii=0;ii<AliPID::kSPECIES;ii++) {
    fPIDBayesProb[ii]     = AliCEPBase::kdumval;
  }
  
  // MC truth
  fMCPID = AliCEPBase::kdumval;
  fMCMass = -999.9;
  fMCMomentum = TVector3(0,0,0);

}

// ----------------------------------------------------------------------------
void CEPTrackBuffer::SetPIDTPCnSigma(Int_t part, Float_t nsig)
{
  
  if (part < AliPID::kSPECIES) {
    fPIDTPCnSigma[part] = nsig;
  } else {
    printf("Wrong particle index! %i is larger than allowed (%i)\n",
      part,AliPID::kSPECIES-1);
  }

}

// ----------------------------------------------------------------------------
void CEPTrackBuffer::SetPIDTPCProbability(Int_t part, Float_t prob)
{

  if (part < AliPID::kSPECIES) {
    fPIDTPCnSigmaProb[part] = prob;
  } else {
    printf("Wrong particle index! %i is larger than allowed (%i)\n",
      part,AliPID::kSPECIES-1);
  }


}

// ----------------------------------------------------------------------------
void CEPTrackBuffer::SetPIDTOFnSigma(Int_t part, Float_t nsig)
{
  
  if (part < AliPID::kSPECIES) {
    fPIDTOFnSigma[part] = nsig;
  } else {
    printf("Wrong particle index! %i is larger than allowed (%i)\n",
      part,AliPID::kSPECIES-1);
  }

}

// ----------------------------------------------------------------------------
void CEPTrackBuffer::SetPIDTOFProbability(Int_t part, Float_t prob)
{

  if (part < AliPID::kSPECIES) {
    fPIDTOFnSigmaProb[part] = prob;
  } else {
    printf("Wrong particle index! %i is larger than allowed (%i)\n",
      part,AliPID::kSPECIES-1);
  }


}

// ----------------------------------------------------------------------------
void CEPTrackBuffer::SetPIDBayesProbability(Int_t part, Float_t prob)
{

  if (part < AliPID::kSPECIES) {
    fPIDBayesProb[part] = prob;
  } else {
    printf("Wrong particle index! %i is larger than allowed (%i)\n",
      part,AliPID::kSPECIES-1);
  }


}

// ----------------------------------------------------------------------------
Float_t CEPTrackBuffer::GetPIDTPCnSigma(Int_t part)
{

  Float_t nsig = AliCEPBase::kdumval;
  
  if (part < AliPID::kSPECIES) {
    nsig = fPIDTPCnSigma[part];
  } else {
    printf("Wrong particle index! %i is larger than allowed (%i)\n",
      part,AliPID::kSPECIES-1);
  }
  
  return nsig;

}

// ----------------------------------------------------------------------------
Float_t CEPTrackBuffer::GetPIDTPCProbability(Int_t part)
{

  Float_t prob = AliCEPBase::kdumval;
  
  if (part < AliPID::kSPECIES) {
    prob = fPIDTPCnSigmaProb[part];
  } else {
    printf("Wrong particle index! %i is larger than allowed (%i)\n",
      part,AliPID::kSPECIES-1);
  }

  return prob;

}

// ----------------------------------------------------------------------------
Float_t CEPTrackBuffer::GetPIDTOFnSigma(Int_t part)
{

  Float_t nsig = AliCEPBase::kdumval;
  
  if (part < AliPID::kSPECIES) {
    nsig = fPIDTOFnSigma[part];
  } else {
    printf("Wrong particle index! %i is larger than allowed (%i)\n",
      part,AliPID::kSPECIES-1);
  }
  
  return nsig;

}

// ----------------------------------------------------------------------------
Float_t CEPTrackBuffer::GetPIDTOFProbability(Int_t part)
{

  Float_t prob = AliCEPBase::kdumval;
  
  if (part < AliPID::kSPECIES) {
    prob = fPIDTOFnSigmaProb[part];
  } else {
    printf("Wrong particle index! %i is larger than allowed (%i)\n",
      part,AliPID::kSPECIES-1);
  }

  return prob;

}

// ----------------------------------------------------------------------------
Float_t CEPTrackBuffer::GetPIDBayesProbability(Int_t part)
{

  Float_t prob = AliCEPBase::kdumval;
  
  if (part < AliPID::kSPECIES) {
    prob = fPIDBayesProb[part];
  } else {
    printf("Wrong particle index! %i is larger than allowed (%i)\n",
      part,AliPID::kSPECIES-1);
  }

  return prob;

}

// ----------------------------------------------------------------------------
