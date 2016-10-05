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
  , fisSoft(-1)
  , fChargeSign(0)
  , fMomentum(TVector3(0,0,0))
  , fPID(0)
  , fPIDTPCStatus(CEPTrackBuffer::kdumval)
  , fPIDTPCSignal(0)
  , fPIDTOFStatus(CEPTrackBuffer::kdumval)
  , fPIDTOFSignal(0)
  , fPIDBayesStatus(CEPTrackBuffer::kdumval)
  , fMCPID(0)
  , fMCMass(0)
  , fMCMomentum(TVector3(0,0,0))
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
  fisSoft = 0;   // 0: is not, 1: is soft, -1: undefined
  fChargeSign = CEPTrackBuffer::kdumval;
  fMomentum = TVector3(0,0,0);
  
  // PID information
  fPID = CEPTrackBuffer::kdumval;
  
  // ... from TPC
  fPIDTPCStatus = CEPTrackBuffer::kdumval;
  fPIDTPCSignal = CEPTrackBuffer::kdumval;
  for(Int_t ii=0;ii<AliPID::kSPECIES;ii++) {
    fPIDTPCnSigma[ii]     = CEPTrackBuffer::kdumval;
    fPIDTPCnSigmaProb[ii] = CEPTrackBuffer::kdumval;
  }
  
  // ... from TOF
  fPIDTOFStatus = CEPTrackBuffer::kdumval;
  fPIDTOFSignal = CEPTrackBuffer::kdumval;
  for(Int_t ii=0;ii<AliPID::kSPECIES;ii++) {
    fPIDTOFnSigma[ii]     = CEPTrackBuffer::kdumval;
    fPIDTOFnSigmaProb[ii] = CEPTrackBuffer::kdumval;
  }
  
  // ... Bayes
  fPIDBayesStatus = CEPTrackBuffer::kdumval;
  for(Int_t ii=0;ii<AliPID::kSPECIES;ii++) {
    fPIDBayesProb[ii]     = CEPTrackBuffer::kdumval;
  }
  
  // MC truth
  fMCPID = CEPTrackBuffer::kdumval;
  fMCMass = CEPTrackBuffer::kdumval;
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
void CEPTrackBuffer::CEPTrackBuffer::SetPIDTPCProbability(Int_t part, Float_t prob)
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
void CEPTrackBuffer::CEPTrackBuffer::SetPIDBayesProbability(Int_t part, Float_t prob)
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

  Float_t nsig = CEPTrackBuffer::kdumval;
  
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

  Float_t prob = CEPTrackBuffer::kdumval;
  
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

  Float_t nsig = CEPTrackBuffer::kdumval;
  
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

  Float_t prob = CEPTrackBuffer::kdumval;
  
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

  Float_t prob = CEPTrackBuffer::kdumval;
  
  if (part < AliPID::kSPECIES) {
    prob = fPIDBayesProb[part];
  } else {
    printf("Wrong particle index! %i is larger than allowed (%i)\n",
      part,AliPID::kSPECIES-1);
  }

  return prob;

}

// ----------------------------------------------------------------------------
