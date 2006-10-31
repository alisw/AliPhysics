#include <TRandom.h>
#include <TArrayS.h>

#include "AliITSdcsSSD.h"
#include "AliITSCalibrationSSD.h"
#include "AliITSresponseSSD.h"
#include "AliITSsegmentationSSD.h"
///////////////////////////////////////////////////////////////////////////
//                                                                       //
//  Class AliITSdcsSSD                                                   //
//  describes Detector Control System parameters for one SSD module.     //
//                                                                       //
//  This class stores parametrers such as gain, threshold                //
//  capacitive coupling.                                                 //
//                                                                       //
//  Class takes care of invalid strip menagement during                  //
//  simulation and runtime                                               //
//                                                                       //
//                                                                       //
//  created at: Warsaw University of Technology                          //
//  ver. 1.0    WARSAW, 23.12.1999                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

ClassImp(AliITSdcsSSD)

// Constructor and Destructor
//______________________________________________________________________
  AliITSdcsSSD::AliITSdcsSSD():
fCouplingPR(0),
fCouplingPL(0),
fCouplingNR(0),
fCouplingNL(0),
fNstrips(0),
fNInvalid(0),
fISigma(0),
fInvalidP(0),
fInvalidN(0){
    // Default Constructor

}
//______________________________________________________________________
AliITSdcsSSD::AliITSdcsSSD(AliITSsegmentation *seg, AliITSCalibration *resp):
fCouplingPR(0),
fCouplingPL(0),
fCouplingNR(0),
fCouplingNL(0),
fNstrips(0),
fNInvalid(0),
fISigma(0),
fInvalidP(0),
fInvalidN(0){
    // Standard constructor

    fNstrips =(Float_t) (((AliITSsegmentationSSD*)seg)->Npx());
    
    fInvalidP = new TArrayS();
    fInvalidN = new TArrayS();

    Int_t npar=((AliITSCalibrationSSD*)resp)->NDetParam();
    if (npar < 6) {
	Warning("AliITSdcsSSD","I need 6 parameters ");
	npar=6;
    } // end if

    Double_t *detpar= new Double_t[npar];
    resp->GetDetParam(detpar);

    fNInvalid = detpar[0];
    fISigma   = detpar[1];

    fCouplingPR = detpar[2];
    fCouplingPL = detpar[3];
    fCouplingNR = detpar[4];
    fCouplingNL = detpar[5];

    char opt[30],dummy[20];
    ((AliITSCalibrationSSD*)resp)->GetParamOptions(opt,dummy);
    if (strstr(opt,"SetInvalid")) SetInvalidMC(fNInvalid,fISigma);

    delete [] detpar;
}
//______________________________________________________________________
AliITSdcsSSD::~AliITSdcsSSD() {
    // destructor

    delete fInvalidP;
    delete fInvalidN;
}
//______________________________________________________________________
AliITSdcsSSD::AliITSdcsSSD(const AliITSdcsSSD &source) : TObject(source),
fCouplingPR(source.fCouplingPR),
fCouplingPL(source.fCouplingPL),
fCouplingNR(source.fCouplingNR),
fCouplingNL(source.fCouplingNL),
fNstrips(source.fNstrips),
fNInvalid(source.fNInvalid),
fISigma(source.fISigma),
fInvalidP(source.fInvalidP),
fInvalidN(source.fInvalidN){
    //     Copy Constructor 

}
//_________________________________________________________________________
AliITSdcsSSD& AliITSdcsSSD::operator=(const AliITSdcsSSD &source) {
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

    return *this;
}
//_____________________________________________________________________
//
//  Methods for creating invalid strips
//_____________________________________________________________________
void AliITSdcsSSD::SetInvalidMC(Float_t mean, Float_t sigma) {
    // set invalid MC

    SetInvalidParam(mean, sigma);
    SetInvalidMC();
}
//______________________________________________________________________
void AliITSdcsSSD::SetInvalidMC() {
    // set invalid MC
    Int_t pside;
    Int_t nside;
    Int_t i;
    Int_t strip;

    pside = (Int_t)gRandom->Gaus(fNInvalid, fISigma);
    nside = (Int_t)gRandom->Gaus(fNInvalid, fISigma);
    
    fInvalidP->Set(pside);
    fInvalidN->Set(nside);
     
    for(i=0 ;i<pside; i++) {
       strip = (Int_t)(gRandom->Rndm() * fNstrips);
       fInvalidP->AddAt(strip, i); 
    } // end for i
    
    for(i=0 ;i<nside; i++) {
       strip = (Int_t)(gRandom->Rndm() * fNstrips);
       fInvalidN->AddAt(strip, i); 
    } // end for i
}
//______________________________________________________________________
void AliITSdcsSSD::SetInvalidParam(Float_t mean, Float_t sigma) {
    // set invalid param

    fNInvalid = mean;
    fISigma = sigma;

    fNInvalid = (fNInvalid<0)? 0 : fNInvalid;
    fNInvalid = (fNInvalid>fNstrips)? fNstrips: fNInvalid;
    
    fISigma = (fISigma < 0)? 0 : fISigma;
    fISigma = (fISigma > fNstrips/10) ? fNstrips/10 : fISigma;
}
//______________________________________________________________________
void AliITSdcsSSD::GetInvalidParam(Float_t &mean, Float_t &sigma) const {
    // get invalid param

    mean = fNInvalid;
    sigma = fISigma;
}
//_____________________________________________________________________
//
//  Methods for accessing to invalid strips
//_____________________________________________________________________
Bool_t AliITSdcsSSD::IsValidP(Int_t strip) {
    // isvalidP
    Int_t nElem = fInvalidP->GetSize();

    for(Int_t i = 0; i<nElem; i++) if(fInvalidP->At(i) == strip) return kFALSE;
    return kTRUE;
}
//______________________________________________________________________
Bool_t AliITSdcsSSD::IsValidN(Int_t strip) {
    // is valid N
    Int_t nElem = fInvalidN->GetSize();

    for(Int_t i = 0; i<nElem; i++) if(fInvalidN->At(i) == strip) return kFALSE;
    return kTRUE;
}
//______________________________________________________________________
TArrayS* AliITSdcsSSD::GetInvalidP() {
    // get invalid P

    return fInvalidP;
}
//______________________________________________________________________
TArrayS* AliITSdcsSSD::GetInvalidN() {
    // get invalid N

    return fInvalidN;
}
//______________________________________________________________________
Int_t AliITSdcsSSD::GetNInvalidP(){
    // get numeber of invalid P

    return fInvalidP->GetSize();
}
//______________________________________________________________________
Int_t AliITSdcsSSD::GetNInvalidN() {
    // get number of invalid N

    return fInvalidN->GetSize();
}
//_____________________________________________________________________





