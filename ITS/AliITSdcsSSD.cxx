
#include "AliITSdcsSSD.h"
#include "AliITSresponseSSD.h"
#include "AliITSsegmentationSSD.h"

ClassImp(AliITSdcsSSD)


//_____________________________________________________________________
//
// Constructor and Destructor
//_____________________________________________________________________


AliITSdcsSSD::AliITSdcsSSD(AliITSsegmentation *seg, AliITSresponse *resp) 
{
  // constructor

    fRandom = new TRandom();

    fNstrips = seg->Npx();
    
    fInvalidP = new TArrayS();
    fInvalidN = new TArrayS();

    Int_t npar=resp->NDetParam();
    if (npar < 6) {
       Warning("AliITSdcsSSD","I need 6 parameters ");
       npar=6;
    }

    Float_t *detpar = new Float_t [npar];
    resp->GetDetParam(detpar);

    fNInvalid = detpar[0];
    fISigma = detpar[1];

    fCouplingPR = detpar[2];
    fCouplingPL = detpar[3];
    fCouplingNR = detpar[4];
    fCouplingNL = detpar[5]; 


    Option_t *opt,*dummy;
    resp->ParamOptions(opt,dummy);
    if (strstr(opt,"SetInvalid")) SetInvalidMC(fNInvalid,fISigma);
	 
	 delete [] detpar;
	 delete detpar;

}

//_____________________________________________________________________


AliITSdcsSSD::~AliITSdcsSSD() {
  // destructor
    delete fInvalidP;
    delete fInvalidN;
}

//__________________________________________________________________________
AliITSdcsSSD::AliITSdcsSSD(const AliITSdcsSSD &source){
  //     Copy Constructor 
  if(&source == this) return;
  this->fCouplingPR = source.fCouplingPR;
  this->fCouplingPL = source.fCouplingPL;
  this->fCouplingNR = source.fCouplingNR;
  this->fCouplingNL = source.fCouplingNL;
  this->fNstrips = source.fNstrips;
  this->fNInvalid = source.fNInvalid;
  this->fISigma = source.fISigma;
  this->fInvalidP = source.fInvalidP;
  this->fInvalidN = source.fInvalidN;
  this->fRandom = source.fRandom;
  return;
}

//_________________________________________________________________________
AliITSdcsSSD& 
  AliITSdcsSSD::operator=(const AliITSdcsSSD &source) {
  //    Assignment operator
  if(&source == this) return *this;
  this->fCouplingPR = source.fCouplingPR;
  this->fCouplingPL = source.fCouplingPL;
  this->fCouplingNR = source.fCouplingNR;
  this->fCouplingNL = source.fCouplingNL;
  this->fNstrips = source.fNstrips;
  this->fNInvalid = source.fNInvalid;
  this->fISigma = source.fISigma;
  this->fInvalidP = source.fInvalidP;
  this->fInvalidN = source.fInvalidN;
  this->fRandom = source.fRandom;
  return *this;
}

//_____________________________________________________________________
//
//  Methods for creating invalid strips
//_____________________________________________________________________
//


void AliITSdcsSSD::SetInvalidMC(Float_t mean, Float_t sigma) {
  // set invalid MC
    SetInvalidParam(mean, sigma);
    SetInvalidMC();
}

//_____________________________________________________________________

void AliITSdcsSSD::SetInvalidMC() {
  // set invalid MC
    Int_t pside;
    Int_t nside;
    Int_t i;
    Int_t strip;

    pside = (Int_t)fRandom->Gaus(fNInvalid, fISigma);
    nside = (Int_t)fRandom->Gaus(fNInvalid, fISigma);
    
    fInvalidP->Set(pside);
    fInvalidN->Set(nside);
     
    for(i=0 ;i<pside; i++) {
       strip = (Int_t)(fRandom->Rndm() * fNstrips);
       fInvalidP->AddAt(strip, i); 
    }
    
    for(i=0 ;i<nside; i++) {
       strip = (Int_t)(fRandom->Rndm() * fNstrips);
       fInvalidN->AddAt(strip, i); 
    }  
}

//_____________________________________________________________________


void AliITSdcsSSD::SetInvalidParam(Float_t mean, Float_t sigma) {
  // set invalid param
    fNInvalid = mean;
    fISigma = sigma;

    fNInvalid = (fNInvalid<0)? 0 : fNInvalid;
    fNInvalid = (fNInvalid>fNstrips)? fNstrips: fNInvalid;
    
    fISigma = (fISigma < 0)? 0 : fISigma;
    fISigma = (fISigma > fNstrips/10) ? fNstrips/10 : fISigma;
}


//_____________________________________________________________________


void AliITSdcsSSD::GetInvalidParam(Float_t &mean, Float_t &sigma) {
  // get invalid param
    mean = fNInvalid;
    sigma = fISigma;

}


//_____________________________________________________________________
//
//  Methods for accessing to invalid strips
//_____________________________________________________________________
//


Bool_t AliITSdcsSSD::IsValidP(Int_t strip) {
  // isvalidP
    Int_t nElem = fInvalidP->GetSize();
    Int_t i;
    for(i = 0; i<nElem; i++)
       if(fInvalidP->At(i) == strip) return kFALSE;
    return kTRUE;
}

//_____________________________________________________________________

Bool_t AliITSdcsSSD::IsValidN(Int_t strip) {
  // is valid N
    Int_t nElem = fInvalidN->GetSize();
    Int_t i;
    for(i = 0; i<nElem; i++)
       if(fInvalidN->At(i) == strip) return kFALSE;
    return kTRUE;
}

//_____________________________________________________________________


TArrayS* AliITSdcsSSD::GetInvalidP() {
  // get invalid P
    return fInvalidP;
}

TArrayS* AliITSdcsSSD::GetInvalidN() {
  // get invalid N
    return fInvalidN;
}

Int_t AliITSdcsSSD::GetNInvalidP(){
  // get numeber of invalid P    
    return fInvalidP->GetSize();
}

Int_t AliITSdcsSSD::GetNInvalidN() {
  // // get number of invalid N
    return fInvalidN->GetSize();
}

//_____________________________________________________________________







