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
  , fChargeSign(0x0)
  , fMomentum()
  , fPID(0)
  , fMCPID(-1)
  , fMCMass(-999.9)
  , fMCMomentum()
{

  fChargeSign = CEPTrackBuffer::kdumval;
  fMomentum   = TVector3(0,0,0);
  
}

// ----------------------------------------------------------------------------
CEPTrackBuffer::~CEPTrackBuffer()
{
  
  this->Delete();
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
