/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliPIDResponse.cxx 46193 2010-12-21 09:00:14Z wiechula $ */

//-----------------------------------------------------------------
//        Base class for handling the pid response               //
//        functions of all detectors                             //
//        and give access to the nsigmas                         //
//                                                               //
//   Origin: Jens Wiechula, Uni Tuebingen, jens.wiechula@cern.ch //
//-----------------------------------------------------------------

#include <TList.h>
#include <TObjArray.h>
#include <TPRegexp.h>
#include <TF1.h>
#include <TSpline.h>
#include <TFile.h>
#include <TArrayI.h>
#include <TArrayF.h>

#include <AliVEvent.h>
#include <AliVTrack.h>
#include <AliLog.h>
#include <AliPID.h>
#include <AliOADBContainer.h>
#include <AliTRDPIDResponseObject.h>
#include <AliTOFPIDParams.h>

#include "AliPIDResponse.h"
#include "AliDetectorPID.h"

#include "AliCentrality.h"

ClassImp(AliPIDResponse);

AliPIDResponse::AliPIDResponse(Bool_t isMC/*=kFALSE*/) :
TNamed("PIDResponse","PIDResponse"),
fITSResponse(isMC),
fTPCResponse(),
fTRDResponse(),
fTOFResponse(),
fEMCALResponse(),
fRange(5.),
fITSPIDmethod(kITSTruncMean),
fIsMC(isMC),
fOADBPath(),
fCustomTPCpidResponse(),
fBeamType("PP"),
fLHCperiod(),
fMCperiodTPC(),
fMCperiodUser(),
fCurrentFile(),
fRecoPass(0),
fRecoPassUser(-1),
fRun(0),
fOldRun(0),
fArrPidResponseMaster(NULL),
fResolutionCorrection(NULL),
fOADBvoltageMaps(NULL),
fTRDPIDResponseObject(NULL),
fTOFtail(1.1),
fTOFPIDParams(NULL),
fEMCALPIDParams(NULL),
fCurrentEvent(NULL),
fCurrCentrality(0.0),
fTuneMConData(kFALSE)
{
  //
  // default ctor
  //
  AliLog::SetClassDebugLevel("AliPIDResponse",0);
  AliLog::SetClassDebugLevel("AliESDpid",0);
  AliLog::SetClassDebugLevel("AliAODpidUtil",0);

}

//______________________________________________________________________________
AliPIDResponse::~AliPIDResponse()
{
  //
  // dtor
  //
  delete fArrPidResponseMaster;
  delete fTRDPIDResponseObject;
  delete fTOFPIDParams;
}

//______________________________________________________________________________
AliPIDResponse::AliPIDResponse(const AliPIDResponse &other) :
TNamed(other),
fITSResponse(other.fITSResponse),
fTPCResponse(other.fTPCResponse),
fTRDResponse(other.fTRDResponse),
fTOFResponse(other.fTOFResponse),
fEMCALResponse(other.fEMCALResponse),
fRange(other.fRange),
fITSPIDmethod(other.fITSPIDmethod),
fIsMC(other.fIsMC),
fOADBPath(other.fOADBPath),
fCustomTPCpidResponse(other.fCustomTPCpidResponse),
fBeamType("PP"),
fLHCperiod(),
fMCperiodTPC(),
fMCperiodUser(other.fMCperiodUser),
fCurrentFile(),
fRecoPass(0),
fRecoPassUser(other.fRecoPassUser),
fRun(0),
fOldRun(0),
fArrPidResponseMaster(NULL),
fResolutionCorrection(NULL),
fOADBvoltageMaps(NULL),
fTRDPIDResponseObject(NULL),
fTOFtail(1.1),
fTOFPIDParams(NULL),
fEMCALPIDParams(NULL),
fCurrentEvent(NULL),
fCurrCentrality(0.0),
fTuneMConData(kFALSE)
{
  //
  // copy ctor
  //
}

//______________________________________________________________________________
AliPIDResponse& AliPIDResponse::operator=(const AliPIDResponse &other)
{
  //
  // copy ctor
  //
  if(this!=&other) {
    delete fArrPidResponseMaster;
    TNamed::operator=(other);
    fITSResponse=other.fITSResponse;
    fTPCResponse=other.fTPCResponse;
    fTRDResponse=other.fTRDResponse;
    fTOFResponse=other.fTOFResponse;
    fEMCALResponse=other.fEMCALResponse;
    fRange=other.fRange;
    fITSPIDmethod=other.fITSPIDmethod;
    fOADBPath=other.fOADBPath;
    fCustomTPCpidResponse=other.fCustomTPCpidResponse;
    fIsMC=other.fIsMC;
    fBeamType="PP";
    fLHCperiod="";
    fMCperiodTPC="";
    fMCperiodUser=other.fMCperiodUser;
    fCurrentFile="";
    fRecoPass=0;
    fRecoPassUser=other.fRecoPassUser;
    fRun=0;
    fOldRun=0;
    fArrPidResponseMaster=NULL;
    fResolutionCorrection=NULL;
    fOADBvoltageMaps=NULL;
    fTRDPIDResponseObject=NULL;
    fEMCALPIDParams=NULL;
    fTOFtail=1.1;
    fTOFPIDParams=NULL;
    fCurrentEvent=other.fCurrentEvent;

  }
  return *this;
}

//______________________________________________________________________________
Float_t AliPIDResponse::NumberOfSigmas(EDetCode detCode, const AliVParticle *track, AliPID::EParticleType type) const
{
  //
  // NumberOfSigmas for 'detCode'
  //

  switch (detCode){
    case kDetITS: return NumberOfSigmasITS(track, type); break;
    case kDetTPC: return NumberOfSigmasTPC(track, type); break;
    case kDetTOF: return NumberOfSigmasTOF(track, type); break;
    case kDetEMCAL: return NumberOfSigmasEMCAL(track, type); break;
    default: return -999.;
  }

}

//______________________________________________________________________________
Float_t AliPIDResponse::NumberOfSigmas(EDetector detCode, const AliVParticle *track, AliPID::EParticleType type) const
{
  //
  // NumberOfSigmas for 'detCode'
  //
  return NumberOfSigmas((EDetCode)(1<<detCode), track, type);
}

//______________________________________________________________________________
Float_t AliPIDResponse::NumberOfSigmasITS(const AliVParticle *vtrack, AliPID::EParticleType type) const
{
  //
  // Calculate the number of sigmas in the ITS
  //
  
  AliVTrack *track=(AliVTrack*)vtrack;
  
  // look for cached value first
  // only the non SA tracks are cached
  if ( track->GetDetectorPID() ){
    return track->GetDetectorPID()->GetNumberOfSigmas(kITS, type);
  }
  
  Float_t dEdx=track->GetITSsignal();
  if (dEdx<=0) return -999.;
  
  UChar_t clumap=track->GetITSClusterMap();
  Int_t nPointsForPid=0;
  for(Int_t i=2; i<6; i++){
    if(clumap&(1<<i)) ++nPointsForPid;
  }
  Float_t mom=track->P();

  //check for ITS standalone tracks
  Bool_t isSA=kTRUE;
  if( track->GetStatus() & AliVTrack::kTPCin ) isSA=kFALSE;
  
  //TODO: in case of the electron, use the SA parametrisation,
  //      this needs to be changed if ITS provides a parametrisation
  //      for electrons also for ITS+TPC tracks
  return fITSResponse.GetNumberOfSigmas(mom,dEdx,type,nPointsForPid,isSA || (type==AliPID::kElectron));
}

//______________________________________________________________________________
Float_t AliPIDResponse::NumberOfSigmasTPC(const AliVParticle *vtrack, AliPID::EParticleType type) const
{
  //
  // Calculate the number of sigmas in the TPC
  //
  
  AliVTrack *track=(AliVTrack*)vtrack;
  
  // look for cached value first
  if (track->GetDetectorPID()){
    return track->GetDetectorPID()->GetNumberOfSigmas(kTPC, type);
  }
  
  Double_t mom  = track->GetTPCmomentum();
  Double_t sig  = track->GetTPCsignal();
  if(fTuneMConData) sig = this->GetTPCsignalTunedOnData(track);
  UInt_t   sigN = track->GetTPCsignalN();
  
  Double_t nSigma = -999.;
  if (sigN>0) nSigma=fTPCResponse.GetNumberOfSigmas(mom,sig,sigN,type);
  
  return nSigma;
}

//______________________________________________________________________________
Float_t AliPIDResponse::NumberOfSigmasTPC( const AliVParticle *vtrack, 
                                           AliPID::EParticleType type,
                                           AliTPCPIDResponse::ETPCdEdxSource dedxSource) 
{
  //get number of sigmas according the selected TPC gain configuration scenario
  const AliVTrack *track=static_cast<const AliVTrack*>(vtrack);

  Float_t nSigma=fTPCResponse.GetNumberOfSigmas(track, type, dedxSource);

  return nSigma;
}

//______________________________________________________________________________
Float_t AliPIDResponse::NumberOfSigmasEMCAL(const AliVParticle *vtrack, AliPID::EParticleType type) const
{
  //
  // Calculate the number of sigmas in the EMCAL
  //
  
  AliVTrack *track=(AliVTrack*)vtrack;

  // look for cached value first
  if (track->GetDetectorPID()){
    return track->GetDetectorPID()->GetNumberOfSigmas(kEMCAL, type);
  }
  
  AliVCluster *matchedClus = NULL;

  Double_t mom     = -1.; 
  Double_t pt      = -1.; 
  Double_t EovP    = -1.;
  Double_t fClsE   = -1.;
  
  Int_t nMatchClus = -1;
  Int_t charge     = 0;
  
  // Track matching
  nMatchClus = track->GetEMCALcluster();
  if(nMatchClus > -1){
    
    mom    = track->P();
    pt     = track->Pt();
    charge = track->Charge();
    
    matchedClus = (AliVCluster*)fCurrentEvent->GetCaloCluster(nMatchClus);
    
    if(matchedClus){
      
      // matched cluster is EMCAL
      if(matchedClus->IsEMCAL()){
	
	fClsE       = matchedClus->E();
	EovP        = fClsE/mom;
	
	
	// NSigma value really meaningful only for electrons!
	return fEMCALResponse.GetNumberOfSigmas(pt,EovP,type,charge); 
      }
    }
  }
  
  return -999;
  
}

//______________________________________________________________________________
Float_t  AliPIDResponse::NumberOfSigmasEMCAL(const AliVParticle *vtrack, AliPID::EParticleType type, Double_t &eop, Double_t showershape[4]) const {

  AliVTrack *track=(AliVTrack*)vtrack;
  
  AliVCluster *matchedClus = NULL;

  Double_t mom     = -1.; 
  Double_t pt      = -1.; 
  Double_t EovP    = -1.;
  Double_t fClsE   = -1.;

  // initialize eop and shower shape parameters
  eop = -1.;
  for(Int_t i = 0; i < 4; i++){
    showershape[i] = -1.;
  }
  
  Int_t nMatchClus = -1;
  Int_t charge     = 0;
  
  // Track matching
  nMatchClus = track->GetEMCALcluster();
  if(nMatchClus > -1){

    mom    = track->P();
    pt     = track->Pt();
    charge = track->Charge();
    
    matchedClus = (AliVCluster*)fCurrentEvent->GetCaloCluster(nMatchClus);
    
    if(matchedClus){
      
      // matched cluster is EMCAL
      if(matchedClus->IsEMCAL()){
	
	fClsE       = matchedClus->E();
	EovP        = fClsE/mom;
	
	// fill used EMCAL variables here
	eop            = EovP; // E/p
	showershape[0] = matchedClus->GetNCells(); // number of cells in cluster
	showershape[1] = matchedClus->GetM02(); // long axis
	showershape[2] = matchedClus->GetM20(); // short axis
	showershape[3] = matchedClus->GetDispersion(); // dispersion
	
	// NSigma value really meaningful only for electrons!
	return fEMCALResponse.GetNumberOfSigmas(pt,EovP,type,charge); 
      }
    }
  }
  return -999;
}


//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::ComputePIDProbability  (EDetCode detCode,  const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  //
  // Compute PID response of 'detCode'
  //

  switch (detCode){
    case kDetITS: return ComputeITSProbability(track, nSpecies, p); break;
    case kDetTPC: return ComputeTPCProbability(track, nSpecies, p); break;
    case kDetTRD: return ComputeTRDProbability(track, nSpecies, p); break;
    case kDetTOF: return ComputeTOFProbability(track, nSpecies, p); break;
    case kDetPHOS: return ComputePHOSProbability(track, nSpecies, p); break;
    case kDetEMCAL: return ComputeEMCALProbability(track, nSpecies, p); break;
    case kDetHMPID: return ComputeHMPIDProbability(track, nSpecies, p); break;
    default: return kDetNoSignal;
  }
}

//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::ComputePIDProbability  (EDetector detCode,  const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  //
  // Compute PID response of 'detCode'
  //

  return ComputePIDProbability((EDetCode)(1<<detCode),track,nSpecies,p);
}

//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::ComputeITSProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  //
  // Compute PID response for the ITS
  //

  // look for cached value first
  // only the non SA tracks are cached
  if (track->GetDetectorPID()){
    return track->GetDetectorPID()->GetRawProbability(kITS, p, nSpecies);
  }

  // set flat distribution (no decision)
  for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;

  if ((track->GetStatus()&AliVTrack::kITSin)==0 &&
    (track->GetStatus()&AliVTrack::kITSout)==0) return kDetNoSignal;

  //check for ITS standalone tracks
  Bool_t isSA=kTRUE;
  if( track->GetStatus() & AliVTrack::kTPCin ) isSA=kFALSE;
  
  Double_t mom=track->P();
  Double_t dedx=track->GetITSsignal();
  Double_t momITS=mom;
  UChar_t clumap=track->GetITSClusterMap();
  Int_t nPointsForPid=0;
  for(Int_t i=2; i<6; i++){
    if(clumap&(1<<i)) ++nPointsForPid;
  }
  
  if(nPointsForPid<3) { // track not to be used for combined PID purposes
    //       track->ResetStatus(AliVTrack::kITSpid);
    return kDetNoSignal;
  }

  Bool_t mismatch=kTRUE/*, heavy=kTRUE*/;
  for (Int_t j=0; j<AliPID::kSPECIES; j++) {
    Double_t mass=AliPID::ParticleMassZ(j);//GeV/c^2
    const Double_t chargeFactor = TMath::Power(AliPID::ParticleCharge(j),2.);
    Double_t bethe=fITSResponse.Bethe(momITS,mass)*chargeFactor;
    //TODO: in case of the electron, use the SA parametrisation,
    //      this needs to be changed if ITS provides a parametrisation
    //      for electrons also for ITS+TPC tracks
    Double_t sigma=fITSResponse.GetResolution(bethe,nPointsForPid,isSA || (j==(Int_t)AliPID::kElectron));
    if (TMath::Abs(dedx-bethe) > fRange*sigma) {
      p[j]=TMath::Exp(-0.5*fRange*fRange)/sigma;
    } else {
      p[j]=TMath::Exp(-0.5*(dedx-bethe)*(dedx-bethe)/(sigma*sigma))/sigma;
      mismatch=kFALSE;
    }

    // Check for particles heavier than (AliPID::kSPECIES - 1)
    //       if (dedx < (bethe + fRange*sigma)) heavy=kFALSE;

  }

  if (mismatch){
    for (Int_t j=0; j<AliPID::kSPECIES; j++) p[j]=1./AliPID::kSPECIES;
    return kDetNoSignal;
  }

    
  return kDetPidOk;
}
//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::ComputeTPCProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  //
  // Compute PID response for the TPC
  //
  
  // look for cached value first
  if (track->GetDetectorPID()){
    return track->GetDetectorPID()->GetRawProbability(kTPC, p, nSpecies);
  }
  
  // set flat distribution (no decision)
  for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;

  // check quality of the track
  if ( (track->GetStatus()&AliVTrack::kTPCin )==0 && (track->GetStatus()&AliVTrack::kTPCout)==0 ) return kDetNoSignal;

  Double_t mom = track->GetTPCmomentum();

  Double_t dedx=track->GetTPCsignal();
  Bool_t mismatch=kTRUE/*, heavy=kTRUE*/;

  if(fTuneMConData) dedx = this->GetTPCsignalTunedOnData(track);

  for (Int_t j=0; j<AliPID::kSPECIESC; j++) {
    AliPID::EParticleType type=AliPID::EParticleType(j);
    Double_t bethe=fTPCResponse.GetExpectedSignal(mom,type);
    Double_t sigma=fTPCResponse.GetExpectedSigma(mom,track->GetTPCsignalN(),type);
    if (TMath::Abs(dedx-bethe) > fRange*sigma) {
      p[j]=TMath::Exp(-0.5*fRange*fRange)/sigma;
    } else {
      p[j]=TMath::Exp(-0.5*(dedx-bethe)*(dedx-bethe)/(sigma*sigma))/sigma;
      mismatch=kFALSE;
    }
  }

  if (mismatch){
    for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;
    return kDetNoSignal;
  }

  return kDetPidOk;
}
//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::ComputeTOFProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  //
  // Compute PID response for the
  //
  
  // look for cached value first
  if (track->GetDetectorPID()){
    return track->GetDetectorPID()->GetRawProbability(kTOF, p, nSpecies);
  }
  
  Double_t meanCorrFactor = 0.11/fTOFtail; // Correction factor on the mean because of the tail (should be ~ 0.1 with tail = 1.1)

  // set flat distribution (no decision)
  for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;
  
  if ((track->GetStatus()&AliVTrack::kTOFout)==0) return kDetNoSignal;
  if ((track->GetStatus()&AliVTrack::kTIME)==0) return kDetNoSignal;
  
  Bool_t mismatch = kTRUE/*, heavy = kTRUE*/;
  for (Int_t j=0; j<AliPID::kSPECIESC; j++) {
    AliPID::EParticleType type=AliPID::EParticleType(j);
    Double_t nsigmas=NumberOfSigmasTOF(track,type) + meanCorrFactor;

    Double_t expTime = fTOFResponse.GetExpectedSignal(track,type);
    Double_t sig = fTOFResponse.GetExpectedSigma(track->P(),expTime,AliPID::ParticleMassZ(type));
    if (TMath::Abs(nsigmas) > (fRange+2)) {
      if(nsigmas < fTOFtail)
	p[j] = TMath::Exp(-0.5*(fRange+2)*(fRange+2))/sig;
      else
	p[j] = TMath::Exp(-(fRange+2 - fTOFtail*0.5)*fTOFtail)/sig;
    } else{
      if(nsigmas < fTOFtail)
	p[j] = TMath::Exp(-0.5*nsigmas*nsigmas)/sig;
      else
	p[j] = TMath::Exp(-(nsigmas - fTOFtail*0.5)*fTOFtail)/sig;
    }

    if (TMath::Abs(nsigmas)<5.){
      Double_t nsigmasTPC=NumberOfSigmasTPC(track,type);
      if (TMath::Abs(nsigmasTPC)<5.) mismatch=kFALSE;
    }
  }

  if (mismatch){
    return kDetMismatch;    
  }

    // TODO: Light nuclei
    
  return kDetPidOk;
}
//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::ComputeTRDProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[],AliTRDPIDResponse::ETRDPIDMethod PIDmethod) const
{
  //
  // Compute PID response for the
    //
  // look for cached value first
    if (track->GetDetectorPID()&&PIDmethod==AliTRDPIDResponse::kLQ1D){
      AliDebug(3,"Return Cached Value");
      return track->GetDetectorPID()->GetRawProbability(kTRD, p, nSpecies);
  }
  
    UInt_t TRDslicesForPID[2];
    SetTRDSlices(TRDslicesForPID,PIDmethod);
  // set flat distribution (no decision)
  for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;
  if((track->GetStatus()&AliVTrack::kTRDout)==0) return kDetNoSignal;

  Float_t mom[6]={0.};
  Double_t dedx[48]={0.};  // Allocate space for the maximum number of TRD slices
  Int_t nslices = TRDslicesForPID[1] - TRDslicesForPID[0] + 1;
  AliDebug(1, Form("First Slice: %d, Last Slice: %d, Number of slices: %d",  TRDslicesForPID[0], TRDslicesForPID[1], nslices));
  for(UInt_t ilayer = 0; ilayer < 6; ilayer++){
    mom[ilayer] = track->GetTRDmomentum(ilayer);
    for(UInt_t islice = TRDslicesForPID[0]; islice <= TRDslicesForPID[1]; islice++){
      dedx[ilayer*nslices+islice-TRDslicesForPID[0]] = track->GetTRDslice(ilayer, islice);
    }
  }
  fTRDResponse.GetResponse(nslices, dedx, mom, p,PIDmethod);
  return kDetPidOk;
}
//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::ComputeEMCALProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  //
  // Compute PID response for the EMCAL
  //
  
  // look for cached value first
  if (track->GetDetectorPID()){
    return track->GetDetectorPID()->GetRawProbability(kEMCAL, p, nSpecies);
  }

  for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;

  AliVCluster *matchedClus = NULL;

  Double_t mom     = -1.; 
  Double_t pt      = -1.; 
  Double_t EovP    = -1.;
  Double_t fClsE   = -1.;
  
  Int_t nMatchClus = -1;
  Int_t charge     = 0;
  
  // Track matching
  nMatchClus = track->GetEMCALcluster();

  if(nMatchClus > -1){
    
    mom    = track->P();
    pt     = track->Pt();
    charge = track->Charge();
    
    matchedClus = (AliVCluster*)fCurrentEvent->GetCaloCluster(nMatchClus);
    
    if(matchedClus){    
      
      // matched cluster is EMCAL
      if(matchedClus->IsEMCAL()){
      
      fClsE       = matchedClus->E();
      EovP        = fClsE/mom;
      
      
      // compute the probabilities 
        if(fEMCALResponse.ComputeEMCALProbability(nSpecies,pt,EovP,charge,p)){	
          
          // in case everything is OK
	      return kDetPidOk;
        }
      }
    }
  }
  
  // in all other cases set flat distribution (no decision)
  for (Int_t j=0; j<nSpecies; j++) p[j] = 1./nSpecies;
  return kDetNoSignal;
  
}
//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::ComputePHOSProbability (const AliVTrack */*track*/, Int_t nSpecies, Double_t p[]) const
{
  //
  // Compute PID response for the PHOS
  //
  
  // look for cached value first
//   if (track->GetDetectorPID()){
//     return track->GetDetectorPID()->GetRawProbability(kPHOS, p, nSpecies);
//   }
  
  // set flat distribution (no decision)
  for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;
  return kDetNoSignal;
}
//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::ComputeHMPIDProbability(const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  //
  // Compute PID response for the HMPID
  //


  // look for cached value first
  if (track->GetDetectorPID()){
    return track->GetDetectorPID()->GetRawProbability(kHMPID, p, nSpecies);
  }
  
  // set flat distribution (no decision)
  for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;
  if((track->GetStatus()&AliVTrack::kHMPIDpid)==0) return kDetNoSignal;

  track->GetHMPIDpid(p);
  
  return kDetPidOk;
}

//______________________________________________________________________________
void AliPIDResponse::InitialiseEvent(AliVEvent *event, Int_t pass, Int_t run)
{
  //
  // Apply settings for the current event
  //
  fRecoPass=pass;
  
  fCurrentEvent=NULL;
  if (!event) return;
  fCurrentEvent=event;
  if (run>0) fRun=run;
  else fRun=event->GetRunNumber();
  
  if (fRun!=fOldRun){
    ExecNewRun();
    fOldRun=fRun;
  }
  
  //TPC resolution parametrisation PbPb
  if ( fResolutionCorrection ){
    Double_t corrSigma=fResolutionCorrection->Eval(GetTPCMultiplicityBin(event));
    fTPCResponse.SetSigma(3.79301e-03*corrSigma, 2.21280e+04);
  }
  
  //TOF resolution
  SetTOFResponse(event, (AliPIDResponse::EStartTimeType_t)fTOFPIDParams->GetStartTimeMethod());


  // Get and set centrality
  AliCentrality *centrality = event->GetCentrality();
  if(centrality){
    fCurrCentrality = centrality->GetCentralityPercentile("V0M");
  }
  else{
    fCurrCentrality = -1;
  }
}

//______________________________________________________________________________
void AliPIDResponse::ExecNewRun()
{
  //
  // Things to Execute upon a new run
  //
  SetRecoInfo();
  
  SetITSParametrisation();
  
  SetTPCPidResponseMaster();
  SetTPCParametrisation();

  SetTRDPidResponseMaster(); 
  InitializeTRDResponse();

  SetEMCALPidResponseMaster(); 
  InitializeEMCALResponse();
  
  SetTOFPidResponseMaster();
  InitializeTOFResponse();

  if (fCurrentEvent) fTPCResponse.SetMagField(fCurrentEvent->GetMagneticField());
}

//_____________________________________________________
Double_t AliPIDResponse::GetTPCMultiplicityBin(const AliVEvent * const event)
{
  //
  // Get TPC multiplicity in bins of 150
  //
  
  const AliVVertex* vertexTPC = event->GetPrimaryVertex();
  Double_t tpcMulti=0.;
  if(vertexTPC){
    Double_t vertexContribTPC=vertexTPC->GetNContributors();
    tpcMulti=vertexContribTPC/150.;
    if (tpcMulti>20.) tpcMulti=20.;
  }
  
  return tpcMulti;
}

//______________________________________________________________________________
void AliPIDResponse::SetRecoInfo()
{
  //
  // Set reconstruction information
  //
  
  //reset information
  fLHCperiod="";
  fMCperiodTPC="";
  
  fBeamType="";
    
  fBeamType="PP";
  
  TPRegexp reg(".*(LHC1[1-2][a-z]+[0-9]+[a-z_]*)/.*");
  TPRegexp reg12a17(".*(LHC12a17[a-z]+)/.*");

  //find the period by run number (UGLY, but not stored in ESD and AOD... )
  if (fRun>=114737&&fRun<=117223)      { fLHCperiod="LHC10B"; fMCperiodTPC="LHC10D1";  }
  else if (fRun>=118503&&fRun<=121040) { fLHCperiod="LHC10C"; fMCperiodTPC="LHC10D1";  }
  else if (fRun>=122195&&fRun<=126437) { fLHCperiod="LHC10D"; fMCperiodTPC="LHC10F6A"; }
  else if (fRun>=127710&&fRun<=130850) { fLHCperiod="LHC10E"; fMCperiodTPC="LHC10F6A"; }
  else if (fRun>=133004&&fRun<=135029) { fLHCperiod="LHC10F"; fMCperiodTPC="LHC10F6A"; }
  else if (fRun>=135654&&fRun<=136377) { fLHCperiod="LHC10G"; fMCperiodTPC="LHC10F6A"; }
  else if (fRun>=136851&&fRun<=139517) {
    fLHCperiod="LHC10H";
    fMCperiodTPC="LHC10H8";
    if (reg.MatchB(fCurrentFile)) fMCperiodTPC="LHC11A10";
    fBeamType="PBPB";
  }
  else if (fRun>=139699&&fRun<=146860) { fLHCperiod="LHC11A"; fMCperiodTPC="LHC10F6A"; }
  //TODO: periods 11B, 11C are not yet treated assume 11d for the moment
  else if (fRun>=148531&&fRun<=155384) { fLHCperiod="LHC11D"; fMCperiodTPC="LHC10F6A"; }
  else if (fRun>=156477&&fRun<=159635) { fLHCperiod="LHC11D"; fMCperiodTPC="LHC10F6A"; }
  // also for 11e,f use 11d
  else if (fRun>=160676&&fRun<=162740) { fLHCperiod="LHC11D"; fMCperiodTPC="LHC10F6A"; }
  else if (fRun>=162933&&fRun<=165746) { fLHCperiod="LHC11D"; fMCperiodTPC="LHC10F6A"; }
  
  else if (fRun>=166529 && fRun<=170718) {
    fLHCperiod="LHC11H";
    fMCperiodTPC="LHC11A10";
    fBeamType="PBPB";
    if (reg12a17.MatchB(fCurrentFile)) fMCperiodTPC="LHC12A17";
  }
  if (fRun>=170719 && fRun<=177311) { fLHCperiod="LHC12A"; fBeamType="PP"; /*fMCperiodTPC="";*/ }
  // for the moment use LHC12b parameters up to LHC12e
  if (fRun>=177312 /*&& fRun<=179356*/) { fLHCperiod="LHC12B"; fBeamType="PP"; /*fMCperiodTPC="";*/ }
//   if (fRun>=179357 && fRun<=183173) { fLHCperiod="LHC12C"; fBeamType="PP"; /*fMCperiodTPC="";*/ }
//   if (fRun>=183174 && fRun<=186345) { fLHCperiod="LHC12D"; fBeamType="PP"; /*fMCperiodTPC="";*/ }
//   if (fRun>=186346 && fRun<=186635) { fLHCperiod="LHC12E"; fBeamType="PP"; /*fMCperiodTPC="";*/ }

//   if (fRun>=186636 && fRun<=188166) { fLHCperiod="LHC12F"; fBeamType="PP"; /*fMCperiodTPC="";*/ }
//   if (fRun >= 188167 && fRun <= 188355 ) { fLHCperiod="LHC12G"; fBeamType="PP"; /*fMCperiodTPC="";*/ }
//   if (fRun >= 188356 && fRun <= 188503 ) { fLHCperiod="LHC12G"; fBeamType="PPB"; /*fMCperiodTPC="";*/ }
// for the moment use 12g parametrisation for all full gain runs (LHC12f+)
  if (fRun >= 186636  ) { fLHCperiod="LHC12G"; fBeamType="PPB"; /*fMCperiodTPC="";*/ }

  //exception new pp MC productions from 2011
  if (fBeamType=="PP" && reg.MatchB(fCurrentFile)) fMCperiodTPC="LHC11B2";
  // exception for 11f1
  if (fCurrentFile.Contains("LHC11f1/")) fMCperiodTPC="LHC11F1";
}

//______________________________________________________________________________
void AliPIDResponse::SetITSParametrisation()
{
  //
  // Set the ITS parametrisation
  //
}

//______________________________________________________________________________
void AliPIDResponse::SetTPCPidResponseMaster()
{
  //
  // Load the TPC pid response functions from the OADB
  // Load the TPC voltage maps from OADB
  //
  //don't load twice for the moment
   if (fArrPidResponseMaster) return;
 

  //reset the PID response functions
  delete fArrPidResponseMaster;
  fArrPidResponseMaster=NULL;
  
  TString fileName(Form("%s/COMMON/PID/data/TPCPIDResponse.root", fOADBPath.Data()));
  TFile *f=NULL;
  if (!fCustomTPCpidResponse.IsNull()) fileName=fCustomTPCpidResponse;
  
  TString fileNamePIDresponse(Form("%s/COMMON/PID/data/TPCPIDResponse.root", fOADBPath.Data()));
  f=TFile::Open(fileNamePIDresponse.Data());
  if (f && f->IsOpen() && !f->IsZombie()){
    fArrPidResponseMaster=dynamic_cast<TObjArray*>(f->Get("TPCPIDResponse"));
  }
  delete f;

  TString fileNameVoltageMaps(Form("%s/COMMON/PID/data/TPCvoltageSettings.root", fOADBPath.Data()));
  f=TFile::Open(fileNameVoltageMaps.Data());
  if (f && f->IsOpen() && !f->IsZombie()){
    fOADBvoltageMaps=dynamic_cast<AliOADBContainer*>(f->Get("TPCvoltageSettings"));
  }
  delete f;
  
  if (!fArrPidResponseMaster){
    AliFatal(Form("Could not retrieve the TPC pid response from: %s",fileNamePIDresponse.Data()));
    return;
  }
  fArrPidResponseMaster->SetOwner();

  if (!fOADBvoltageMaps)
  {
    AliFatal(Form("Could not retrieve the TPC voltage maps from: %s",fileNameVoltageMaps.Data()));
  }
  fArrPidResponseMaster->SetOwner();
}

//______________________________________________________________________________
void AliPIDResponse::SetTPCParametrisation()
{
  //
  // Change BB parametrisation for current run
  //
  
  if (fLHCperiod.IsNull()) {
    AliFatal("No period set, not changing parametrisation");
    return;
  }
  
  //
  // Set default parametrisations for data and MC
  //
  
  //data type
  TString datatype="DATA";
  //in case of mc fRecoPass is per default 1
  if (fIsMC) {
      if(!fTuneMConData) datatype="MC";
      fRecoPass=1;
  }
  
  //
  //reset old splines
  //
  fTPCResponse.ResetSplines();

  // period
  TString period=fLHCperiod;
  if (fIsMC && !fTuneMConData) period=fMCperiodTPC;

  AliInfo(Form("Searching splines for: %s %s PASS%d %s",datatype.Data(),period.Data(),fRecoPass,fBeamType.Data()));
  Bool_t found=kFALSE;
  //
  //set the new PID splines
  //
  if (fArrPidResponseMaster){
    Int_t recopass = fRecoPass;
    if(fTuneMConData) recopass = fRecoPassUser;
    //for MC don't use period information
    //if (fIsMC) period="[A-Z0-9]*";
    //for MC use MC period information
    //pattern for the default entry (valid for all particles)
    TPRegexp reg(Form("TSPLINE3_%s_([A-Z]*)_%s_PASS%d_%s_MEAN(_*)([A-Z1-9]*)",datatype.Data(),period.Data(),recopass,fBeamType.Data()));

    //find particle id ang gain scenario
    for (Int_t igainScenario=0; igainScenario<AliTPCPIDResponse::fgkNumberOfGainScenarios; igainScenario++)
    {
      TObject *grAll=NULL;
      TString gainScenario = AliTPCPIDResponse::GainScenarioName(igainScenario);
      gainScenario.ToUpper();
      //loop over entries and filter them
      for (Int_t iresp=0; iresp<fArrPidResponseMaster->GetEntriesFast();++iresp)
      {
        TObject *responseFunction=fArrPidResponseMaster->At(iresp);
        if (responseFunction==NULL) continue;
        TString responseName=responseFunction->GetName();
         
        if (!reg.MatchB(responseName)) continue;

        TObjArray *arr=reg.MatchS(responseName); if (!arr) continue;
        TObject* tmp=NULL;
        tmp=arr->At(1); if (!tmp) continue;
        TString particleName=tmp->GetName();
        tmp=arr->At(3); if (!tmp) continue;
        TString gainScenarioName=tmp->GetName();
        delete arr;
        if (particleName.IsNull()) continue;
        if (!grAll && particleName=="ALL" && gainScenarioName==gainScenario) grAll=responseFunction;
        else 
        {
          for (Int_t ispec=0; ispec<(AliTPCPIDResponse::fgkNumberOfParticleSpecies); ++ispec)
          {
            TString particle=AliPID::ParticleName(ispec);
            particle.ToUpper();
            //std::cout<<responseName<<" "<<particle<<" "<<particleName<<" "<<gainScenario<<" "<<gainScenarioName<<std::endl;
            if ( particle == particleName && gainScenario == gainScenarioName )
            {
              fTPCResponse.SetResponseFunction( responseFunction,
                                                (AliPID::EParticleType)ispec,
                                                (AliTPCPIDResponse::ETPCgainScenario)igainScenario );
              fTPCResponse.SetUseDatabase(kTRUE);
              AliInfo(Form("Adding graph: %d %d - %s",ispec,igainScenario,responseFunction->GetName()));
              found=kTRUE;
              // overwrite default with proton spline (for light nuclei)
              if (ispec==AliPID::kProton) grAll=responseFunction;
              break;
            }
          }
        }
      }
      if (grAll)
      {
        for (Int_t ispec=0; ispec<(AliTPCPIDResponse::fgkNumberOfParticleSpecies); ++ispec)
        {
          if (!fTPCResponse.GetResponseFunction( (AliPID::EParticleType)ispec,
                                                 (AliTPCPIDResponse::ETPCgainScenario)igainScenario))
          {
              fTPCResponse.SetResponseFunction( grAll,
                                                (AliPID::EParticleType)ispec,
                                                (AliTPCPIDResponse::ETPCgainScenario)igainScenario );
              fTPCResponse.SetUseDatabase(kTRUE);
              AliInfo(Form("Adding graph: %d %d - %s",ispec,igainScenario,grAll->GetName()));
              found=kTRUE;
          }
        }
      }
    }
  }
  else AliInfo("no fArrPidResponseMaster");

  if (!found){
    AliError(Form("No splines found for: %s %s PASS%d %s",datatype.Data(),period.Data(),fRecoPass,fBeamType.Data()));
  }

  //
  // Setup resolution parametrisation
  //
  
  //default
  fTPCResponse.SetSigma(3.79301e-03, 2.21280e+04);
  
  if (fRun>=122195){
    fTPCResponse.SetSigma(2.30176e-02, 5.60422e+02);
  }

  if (fRun>=186636){
//   if (fRun>=188356){
    fTPCResponse.SetSigma(8.62022e-04, 9.08156e+05);
  }
  
  if (fArrPidResponseMaster)
  fResolutionCorrection=(TF1*)fArrPidResponseMaster->FindObject(Form("TF1_%s_ALL_%s_PASS%d_%s_SIGMA",datatype.Data(),period.Data(),fRecoPass,fBeamType.Data()));
  
  if (fResolutionCorrection) AliInfo(Form("Setting multiplicity correction function: %s",fResolutionCorrection->GetName()));

  //read in the voltage map
  TVectorF* gsm = dynamic_cast<TVectorF*>(fOADBvoltageMaps->GetObject(fRun));
  if (gsm) 
  {
    fTPCResponse.SetVoltageMap(*gsm);
    TString vals;
    AliInfo(Form("Reading the voltage map for run %d\n",fRun));
    vals="IROC A: "; for (Int_t i=0; i<18; i++){vals+=Form("%.2f ",(*gsm)[i]);}
    AliInfo(vals.Data());
    vals="IROC C: "; for (Int_t i=18; i<36; i++){vals+=Form("%.2f ",(*gsm)[i]);}
    AliInfo(vals.Data());
    vals="OROC A: "; for (Int_t i=36; i<54; i++){vals+=Form("%.2f ",(*gsm)[i]);}
    AliInfo(vals.Data());
    vals="OROC C: "; for (Int_t i=54; i<72; i++){vals+=Form("%.2f ",(*gsm)[i]);}
    AliInfo(vals.Data());
  }
  else AliInfo("no voltage map, ideal default assumed");
}

//______________________________________________________________________________
void AliPIDResponse::SetTRDPidResponseMaster()
{
  //
  // Load the TRD pid params and references from the OADB
  //
  if(fTRDPIDResponseObject) return;
  AliOADBContainer contParams("contParams"); 

  Int_t statusResponse = contParams.InitFromFile(Form("%s/COMMON/PID/data/TRDPIDResponse.root", fOADBPath.Data()), "AliTRDPIDResponseObject");
  if(statusResponse){
    AliError("Failed initializing PID Response Object from OADB");
  } else {
    AliInfo(Form("Loading TRD Response from %s/COMMON/PID/data/TRDPIDResponse.root", fOADBPath.Data()));
    fTRDPIDResponseObject = dynamic_cast<AliTRDPIDResponseObject *>(contParams.GetObject(fRun));
    if(!fTRDPIDResponseObject){
      AliError(Form("TRD Response not found in run %d", fRun));
    }
  }
  /*
  AliOADBContainer contRefs("contRefs");
  Int_t statusRefs = contRefs.InitFromFile(Form("%s/COMMON/PID/data/TRDPIDReferenceLQ1D.root", fOADBPath.Data()), "AliTRDPIDReference");
  if(statusRefs){
    AliInfo("Failed Loading References for TRD");
  } else {
    AliInfo(Form("Loading TRD References from %s/COMMON/PID/data/TRDPIDReferenceLQ1D.root", fOADBPath.Data()));
    fTRDPIDReference = dynamic_cast<AliTRDPIDReference *>(contRefs.GetObject(fRun));
    if(!fTRDPIDReference){
      AliError(Form("TRD References not found in OADB Container for run %d", fRun));
    }
    }
    */
}

//______________________________________________________________________________
void AliPIDResponse::InitializeTRDResponse(){
  //
  // Set PID Params and references to the TRD PID response
  // 
    fTRDResponse.SetPIDResponseObject(fTRDPIDResponseObject);
}

//______________________________________________________________________________
void AliPIDResponse::SetTRDSlices(UInt_t TRDslicesForPID[2],AliTRDPIDResponse::ETRDPIDMethod method) const{

    if(fLHCperiod == "LHC10d" || fLHCperiod == "LHC10e"){
	// backward compatibility for setting with 8 slices
	TRDslicesForPID[0] = 0;
	TRDslicesForPID[1] = 7;
    }
    else{
	if(method==AliTRDPIDResponse::kLQ1D){
	    TRDslicesForPID[0] = 0; // first Slice contains normalized dEdx
	    TRDslicesForPID[1] = 0;
	}
	if(method==AliTRDPIDResponse::kLQ2D){
	    TRDslicesForPID[0] = 1;
	    TRDslicesForPID[1] = 7;
	}
    }
    AliDebug(1,Form("Slice Range set to %d - %d",TRDslicesForPID[0],TRDslicesForPID[1]));
}

//______________________________________________________________________________
void AliPIDResponse::SetTOFPidResponseMaster()
{
  //
  // Load the TOF pid params from the OADB
  //

  if (fTOFPIDParams) delete fTOFPIDParams;
  fTOFPIDParams=NULL;

  TFile *oadbf = new TFile(Form("%s/COMMON/PID/data/TOFPIDParams.root",fOADBPath.Data()));
  if (oadbf && oadbf->IsOpen()) {
    AliInfo(Form("Loading TOF Params from %s/COMMON/PID/data/TOFPIDParams.root", fOADBPath.Data()));
    AliOADBContainer *oadbc = (AliOADBContainer *)oadbf->Get("TOFoadb");
    if (oadbc) fTOFPIDParams = dynamic_cast<AliTOFPIDParams *>(oadbc->GetObject(fRun,"TOFparams"));
    oadbf->Close();
    delete oadbc;
  }
  delete oadbf;

  if (!fTOFPIDParams) AliFatal("TOFPIDParams could not be retrieved");
}

//______________________________________________________________________________
void AliPIDResponse::InitializeTOFResponse(){
  //
  // Set PID Params to the TOF PID response
  //

  AliInfo("TOF PID Params loaded from OADB");
  AliInfo(Form("  TOF resolution %5.2f [ps]",fTOFPIDParams->GetTOFresolution()));
  AliInfo(Form("  StartTime method %d",fTOFPIDParams->GetStartTimeMethod()));
  AliInfo(Form("  TOF res. mom. params: %5.2f %5.2f %5.2f %5.2f",
               fTOFPIDParams->GetSigParams(0),fTOFPIDParams->GetSigParams(1),fTOFPIDParams->GetSigParams(2),fTOFPIDParams->GetSigParams(3)));
  
  for (Int_t i=0;i<4;i++) {
    fTOFResponse.SetTrackParameter(i,fTOFPIDParams->GetSigParams(i));
  }
  fTOFResponse.SetTimeResolution(fTOFPIDParams->GetTOFresolution());

}


//_________________________________________________________________________
Bool_t AliPIDResponse::IdentifiedAsElectronTRD(const AliVTrack *vtrack, Double_t efficiencyLevel,Double_t centrality,AliTRDPIDResponse::ETRDPIDMethod PIDmethod) const {
  //
  // Check whether track is identified as electron under a given electron efficiency hypothesis
    //

  Double_t probs[AliPID::kSPECIES];
  ComputeTRDProbability(vtrack, AliPID::kSPECIES, probs,PIDmethod);

  Int_t ntracklets = vtrack->GetTRDntrackletsPID();
  // Take mean of the TRD momenta in the given tracklets
  Float_t p = 0, trdmomenta[AliVTrack::kTRDnPlanes];
  Int_t nmomenta = 0;
  for(Int_t iPl=0;iPl<AliVTrack::kTRDnPlanes;iPl++){
    if(vtrack->GetTRDmomentum(iPl) > 0.){
      trdmomenta[nmomenta++] = vtrack->GetTRDmomentum(iPl); 
    }
  }
  p = TMath::Mean(nmomenta, trdmomenta);

  return fTRDResponse.IdentifiedAsElectron(ntracklets, probs, p, efficiencyLevel,centrality,PIDmethod);
}

//______________________________________________________________________________
void AliPIDResponse::SetEMCALPidResponseMaster()
{
  //
  // Load the EMCAL pid response functions from the OADB
  //
  TObjArray* fEMCALPIDParamsRun      = NULL;
  TObjArray* fEMCALPIDParamsPass     = NULL;

  if(fEMCALPIDParams) return;
  AliOADBContainer contParams("contParams"); 

  Int_t statusPars = contParams.InitFromFile(Form("%s/COMMON/PID/data/EMCALPIDParams.root", fOADBPath.Data()), "AliEMCALPIDParams");
  if(statusPars){
    AliError("Failed initializing PID Params from OADB");
  } 
  else {
    AliInfo(Form("Loading EMCAL Params from %s/COMMON/PID/data/EMCALPIDParams.root", fOADBPath.Data()));

    fEMCALPIDParamsRun = dynamic_cast<TObjArray *>(contParams.GetObject(fRun));
    if(fEMCALPIDParamsRun)  fEMCALPIDParamsPass = dynamic_cast<TObjArray *>(fEMCALPIDParamsRun->FindObject(Form("pass%d",fRecoPass)));
    if(fEMCALPIDParamsPass) fEMCALPIDParams     = dynamic_cast<TObjArray *>(fEMCALPIDParamsPass->FindObject(Form("EMCALPIDParams_Particles")));

    if(!fEMCALPIDParams){
      AliInfo(Form("EMCAL Params not found in run %d pass %d", fRun, fRecoPass));
      AliInfo("Will take the standard LHC11d instead ...");

      fEMCALPIDParamsRun = dynamic_cast<TObjArray *>(contParams.GetObject(156477));
      if(fEMCALPIDParamsRun)  fEMCALPIDParamsPass = dynamic_cast<TObjArray *>(fEMCALPIDParamsRun->FindObject(Form("pass%d",1)));
      if(fEMCALPIDParamsPass) fEMCALPIDParams     = dynamic_cast<TObjArray *>(fEMCALPIDParamsPass->FindObject(Form("EMCALPIDParams_Particles")));

      if(!fEMCALPIDParams){
	AliError(Form("DEFAULT EMCAL Params (LHC11d) not found in file %s/COMMON/PID/data/EMCALPIDParams.root", fOADBPath.Data()));	
      }
    }
  }
}

//______________________________________________________________________________
void AliPIDResponse::InitializeEMCALResponse(){
  //
  // Set PID Params to the EMCAL PID response
  // 
  fEMCALResponse.SetPIDParams(fEMCALPIDParams);

}

//_____________________________________________________
void AliPIDResponse::FillTrackDetectorPID()
{
  //
  // create detector PID information and setup the transient pointer in the track
  //

  if (!fCurrentEvent) return;
  
  //TODO: which particles to include? See also the loops below...
  Double_t values[AliPID::kSPECIESC]={0};
  
  for (Int_t itrack=0; itrack<fCurrentEvent->GetNumberOfTracks(); ++itrack){
    AliVTrack *track=dynamic_cast<AliVTrack*>(fCurrentEvent->GetTrack(itrack));
    if (!track) continue;

    AliDetectorPID *detPID=new AliDetectorPID;
    for (Int_t idet=0; idet<kNdetectors; ++idet){

      //nsigmas
      for (Int_t ipart=0; ipart<AliPID::kSPECIESC; ++ipart)
        values[ipart]=NumberOfSigmas((EDetector)idet,track,(AliPID::EParticleType)ipart);
      detPID->SetNumberOfSigmas((EDetector)idet, values, (Int_t)AliPID::kSPECIESC);

      //probabilities
      EDetPidStatus status=ComputePIDProbability((EDetector)idet,track,AliPID::kSPECIESC,values);
      detPID->SetRawProbability((EDetector)idet, values, (Int_t)AliPID::kSPECIESC, status);
    }

    track->SetDetectorPID(detPID);
  }
}

//_________________________________________________________________________
void AliPIDResponse::SetTOFResponse(AliVEvent *vevent,EStartTimeType_t option){
  //
  // Set TOF response function
  // Input option for event_time used
  //
  
    Float_t t0spread = 0.; //vevent->GetEventTimeSpread();
    if(t0spread < 10) t0spread = 80;

    // T0 from TOF algorithm

    Bool_t flagT0TOF=kFALSE;
    Bool_t flagT0T0=kFALSE;
    Float_t *startTime = new Float_t[fTOFResponse.GetNmomBins()];
    Float_t *startTimeRes = new Float_t[fTOFResponse.GetNmomBins()];
    Int_t *startTimeMask = new Int_t[fTOFResponse.GetNmomBins()];

    // T0-TOF arrays
    Float_t *estimatedT0event = new Float_t[fTOFResponse.GetNmomBins()];
    Float_t *estimatedT0resolution = new Float_t[fTOFResponse.GetNmomBins()];
    for(Int_t i=0;i<fTOFResponse.GetNmomBins();i++){
      estimatedT0event[i]=0.0;
      estimatedT0resolution[i]=0.0;
      startTimeMask[i] = 0;
    }

    Float_t resT0A=75,resT0C=65,resT0AC=55;
    if(vevent->GetT0TOF()){ // check if T0 detector information is available
	flagT0T0=kTRUE;
    }


    AliTOFHeader *tofHeader = (AliTOFHeader*)vevent->GetTOFHeader();

    if (tofHeader) { // read global info and T0-TOF
      fTOFResponse.SetTimeResolution(tofHeader->GetTOFResolution());
      t0spread = tofHeader->GetT0spread(); // read t0 sprad
      if(t0spread < 10) t0spread = 80;

      flagT0TOF=kTRUE;
      for(Int_t i=0;i<fTOFResponse.GetNmomBins();i++){ // read T0-TOF default value
	startTime[i]=tofHeader->GetDefaultEventTimeVal();
	startTimeRes[i]=tofHeader->GetDefaultEventTimeRes();
	if(startTimeRes[i] < 1.e-5) startTimeRes[i] = t0spread;
      }

      TArrayI *ibin=(TArrayI*)tofHeader->GetNvalues();
      TArrayF *t0Bin=(TArrayF*)tofHeader->GetEventTimeValues();
      TArrayF *t0ResBin=(TArrayF*)tofHeader->GetEventTimeRes();
      for(Int_t j=0;j < tofHeader->GetNbins();j++){ // fill T0-TOF in p-bins
	Int_t icurrent = (Int_t)ibin->GetAt(j);
	startTime[icurrent]=t0Bin->GetAt(j);
	startTimeRes[icurrent]=t0ResBin->GetAt(j);
	if(startTimeRes[icurrent] < 1.e-5) startTimeRes[icurrent] = t0spread;
      }
    }

    // for cut of 3 sigma on t0 spread
    Float_t t0cut = 3 * t0spread;
    if(t0cut < 500) t0cut = 500;

    if(option == kFILL_T0){ // T0-FILL is used
	for(Int_t i=0;i<fTOFResponse.GetNmomBins();i++){
	  estimatedT0event[i]=0.0;
	  estimatedT0resolution[i]=t0spread;
	}
	fTOFResponse.SetT0event(estimatedT0event);
	fTOFResponse.SetT0resolution(estimatedT0resolution);
    }

    if(option == kTOF_T0){ // T0-TOF is used when available (T0-FILL otherwise) from ESD
	if(flagT0TOF){
	    fTOFResponse.SetT0event(startTime);
	    fTOFResponse.SetT0resolution(startTimeRes);
	    for(Int_t i=0;i<fTOFResponse.GetNmomBins();i++){
	      if(startTimeRes[i]<t0spread) startTimeMask[i]=1;
	      fTOFResponse.SetT0binMask(i,startTimeMask[i]);
	    }
	}
	else{
	    for(Int_t i=0;i<fTOFResponse.GetNmomBins();i++){
	      estimatedT0event[i]=0.0;
	      estimatedT0resolution[i]=t0spread;
	      fTOFResponse.SetT0binMask(i,startTimeMask[i]);
	    }
	    fTOFResponse.SetT0event(estimatedT0event);
	    fTOFResponse.SetT0resolution(estimatedT0resolution);
	}
    }
    else if(option == kBest_T0){ // T0-T0 or T0-TOF are used when available (T0-FILL otherwise) from ESD
	Float_t t0AC=-10000;
	Float_t t0A=-10000;
	Float_t t0C=-10000;
	if(flagT0T0){
	    t0AC= vevent->GetT0TOF()[0];
	    t0A= vevent->GetT0TOF()[1];
	    t0C= vevent->GetT0TOF()[2];
	}

	Float_t t0t0Best = 0;
	Float_t t0t0BestRes = 9999;
	Int_t t0used=0;
	if(TMath::Abs(t0A) < t0cut && TMath::Abs(t0C) < t0cut && TMath::Abs(t0C-t0A) < 500){
	    t0t0Best = t0AC;
	    t0t0BestRes = resT0AC;
	    t0used=6;
	}
	else if(TMath::Abs(t0C) < t0cut){
	    t0t0Best = t0C;
	    t0t0BestRes = resT0C;
	    t0used=4;
	}
	else if(TMath::Abs(t0A) < t0cut){
	    t0t0Best = t0A;
	    t0t0BestRes = resT0A;
	    t0used=2;
	}

	if(flagT0TOF){ // if T0-TOF info is available
	    for(Int_t i=0;i<fTOFResponse.GetNmomBins();i++){
		if(t0t0BestRes < 999){
		  if(startTimeRes[i] < t0spread){
		    Double_t wtot = 1./startTimeRes[i]/startTimeRes[i] + 1./t0t0BestRes/t0t0BestRes;
		    Double_t t0best = startTime[i]/startTimeRes[i]/startTimeRes[i] + t0t0Best/t0t0BestRes/t0t0BestRes;
		    estimatedT0event[i]=t0best / wtot;
		    estimatedT0resolution[i]=1./TMath::Sqrt(wtot);
		    startTimeMask[i] = t0used+1;
		  }
		  else {
		    estimatedT0event[i]=t0t0Best;
		    estimatedT0resolution[i]=t0t0BestRes;
		    startTimeMask[i] = t0used;
		  }
		}
		else{
		  estimatedT0event[i]=startTime[i];
		  estimatedT0resolution[i]=startTimeRes[i];
		  if(startTimeRes[i]<t0spread) startTimeMask[i]=1;
		}
		fTOFResponse.SetT0binMask(i,startTimeMask[i]);
	    }
	    fTOFResponse.SetT0event(estimatedT0event);
	    fTOFResponse.SetT0resolution(estimatedT0resolution);
	}
	else{ // if no T0-TOF info is available
	    for(Int_t i=0;i<fTOFResponse.GetNmomBins();i++){
	      fTOFResponse.SetT0binMask(i,t0used);
	      if(t0t0BestRes < 999){
		estimatedT0event[i]=t0t0Best;
		estimatedT0resolution[i]=t0t0BestRes;
	      }
	      else{
		estimatedT0event[i]=0.0;
		estimatedT0resolution[i]=t0spread;
	      }
	    }
	    fTOFResponse.SetT0event(estimatedT0event);
	    fTOFResponse.SetT0resolution(estimatedT0resolution);
	}
    }

    else if(option == kT0_T0){ // T0-T0 is used when available (T0-FILL otherwise)
	Float_t t0AC=-10000;
	Float_t t0A=-10000;
	Float_t t0C=-10000;
	if(flagT0T0){
	    t0AC= vevent->GetT0TOF()[0];
	    t0A= vevent->GetT0TOF()[1];
	    t0C= vevent->GetT0TOF()[2];
	}

	if(TMath::Abs(t0A) < t0cut && TMath::Abs(t0C) < t0cut && TMath::Abs(t0C-t0A) < 500){
	    for(Int_t i=0;i<fTOFResponse.GetNmomBins();i++){
	      estimatedT0event[i]=t0AC;
	      estimatedT0resolution[i]=resT0AC;
	      fTOFResponse.SetT0binMask(i,6);
	    }
	}
	else if(TMath::Abs(t0C) < t0cut){
	    for(Int_t i=0;i<fTOFResponse.GetNmomBins();i++){
	      estimatedT0event[i]=t0C;
	      estimatedT0resolution[i]=resT0C;
	      fTOFResponse.SetT0binMask(i,4);
	    }
	}
	else if(TMath::Abs(t0A) < t0cut){
	    for(Int_t i=0;i<fTOFResponse.GetNmomBins();i++){
	      estimatedT0event[i]=t0A;
	      estimatedT0resolution[i]=resT0A;
	      fTOFResponse.SetT0binMask(i,2);
	    }
	}
	else{
	    for(Int_t i=0;i<fTOFResponse.GetNmomBins();i++){
	      estimatedT0event[i]=0.0;
	      estimatedT0resolution[i]=t0spread;
	      fTOFResponse.SetT0binMask(i,0);
	    }
	}
	fTOFResponse.SetT0event(estimatedT0event);
	fTOFResponse.SetT0resolution(estimatedT0resolution);
    }
    delete [] startTime;
    delete [] startTimeRes;
    delete [] startTimeMask;
    delete [] estimatedT0event;
    delete [] estimatedT0resolution;
}
