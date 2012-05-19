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

#include <AliVEvent.h>
#include <AliVTrack.h>
#include <AliLog.h>
#include <AliPID.h>
#include <AliOADBContainer.h>
#include <AliTRDPIDParams.h>
#include <AliTRDPIDReference.h>
#include <AliTOFPIDParams.h>

#include "AliPIDResponse.h"

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
fBeamType("PP"),
fLHCperiod(),
fMCperiodTPC(),
fMCperiodUser(),
fCurrentFile(),
fRecoPass(0),
fRecoPassUser(-1),
fRun(0),
fOldRun(0),
fArrPidResponseMaster(0x0),
fResolutionCorrection(0x0),
fTRDPIDParams(0x0),
fTRDPIDReference(0x0),
fTOFtail(1.1),
fTOFPIDParams(0x0),
fEMCALPIDParams(0x0),
fCurrentEvent(0x0),
fCurrCentrality(0.0)
{
  //
  // default ctor
  //
  AliLog::SetClassDebugLevel("AliPIDResponse",10);
  AliLog::SetClassDebugLevel("AliESDpid",10);
  AliLog::SetClassDebugLevel("AliAODpidUtil",10);

  memset(fTRDslicesForPID,0,sizeof(UInt_t)*2);
}

//______________________________________________________________________________
AliPIDResponse::~AliPIDResponse()
{
  //
  // dtor
  //
  delete fArrPidResponseMaster;
  delete fTRDPIDParams;
  delete fTRDPIDReference;
  if (fTOFPIDParams) delete fTOFPIDParams;
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
fBeamType("PP"),
fLHCperiod(),
fMCperiodTPC(),
fMCperiodUser(other.fMCperiodUser),
fCurrentFile(),
fRecoPass(0),
fRecoPassUser(other.fRecoPassUser),
fRun(0),
fOldRun(0),
fArrPidResponseMaster(0x0),
fResolutionCorrection(0x0),
fTRDPIDParams(0x0),
fTRDPIDReference(0x0),
fTOFtail(1.1),
fTOFPIDParams(0x0),
fEMCALPIDParams(0x0),
fCurrentEvent(0x0),
fCurrCentrality(0.0)
{
  //
  // copy ctor
  //
  memset(fTRDslicesForPID,0,sizeof(UInt_t)*2);
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
    fArrPidResponseMaster=0x0;
    fResolutionCorrection=0x0;
    fTRDPIDParams=0x0;
    fTRDPIDReference=0x0;
    fEMCALPIDParams=0x0;
    memset(fTRDslicesForPID,0,sizeof(UInt_t)*2);
    fTOFtail=1.1;
    fTOFPIDParams=0x0;
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
//     case kDetTRD: return ComputeTRDProbability(track, type); break;
//     case kDetPHOS: return ComputePHOSProbability(track, type); break;
//     case kDetEMCAL: return NumberOfSigmasEMCAL(track, type); break;
//     case kDetHMPID: return ComputeHMPIDProbability(track, type); break;
    default: return -999.;
  }

}

//______________________________________________________________________________
Float_t AliPIDResponse::NumberOfSigmasEMCAL(const AliVTrack *track, AliPID::EParticleType type) const {

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
Float_t  AliPIDResponse::NumberOfSigmasEMCAL(const AliVTrack *track, AliPID::EParticleType type, Double_t &eop, Double_t showershape[4]) const {

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
    case kDetTOF: return ComputeTOFProbability(track, nSpecies, p); break;
    case kDetTRD: return ComputeTRDProbability(track, nSpecies, p); break;
    case kDetPHOS: return ComputePHOSProbability(track, nSpecies, p); break;
    case kDetEMCAL: return ComputeEMCALProbability(track, nSpecies, p); break;
    case kDetHMPID: return ComputeHMPIDProbability(track, nSpecies, p); break;
    default: return kDetNoSignal;
  }
}

//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::ComputeITSProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  //
  // Compute PID response for the ITS
  //

  // set flat distribution (no decision)
  for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;

  if ((track->GetStatus()&AliVTrack::kITSin)==0 &&
    (track->GetStatus()&AliVTrack::kITSout)==0) return kDetNoSignal;
  
  Double_t mom=track->P();
  Double_t dedx=track->GetITSsignal();
  Bool_t isSA=kTRUE;
  Double_t momITS=mom;
  ULong_t trStatus=track->GetStatus();
  if(trStatus&AliVTrack::kTPCin) isSA=kFALSE;
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
    Double_t mass=AliPID::ParticleMass(j);//GeV/c^2
    Double_t bethe=fITSResponse.Bethe(momITS,mass);
    Double_t sigma=fITSResponse.GetResolution(bethe,nPointsForPid,isSA);
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

  // set flat distribution (no decision)
  for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;

  // check quality of the track
  if ( (track->GetStatus()&AliVTrack::kTPCin )==0 && (track->GetStatus()&AliVTrack::kTPCout)==0 ) return kDetNoSignal;

  Double_t mom = track->GetTPCmomentum();

  Double_t dedx=track->GetTPCsignal();
  Bool_t mismatch=kTRUE/*, heavy=kTRUE*/;

  for (Int_t j=0; j<AliPID::kSPECIES; j++) {
    AliPID::EParticleType type=AliPID::EParticleType(j);
    Double_t bethe=fTPCResponse.GetExpectedSignal(mom,type);
    Double_t sigma=fTPCResponse.GetExpectedSigma(mom,track->GetTPCsignalN(),type);
    if (TMath::Abs(dedx-bethe) > fRange*sigma) {
      p[j]=TMath::Exp(-0.5*fRange*fRange)/sigma;
    } else {
      p[j]=TMath::Exp(-0.5*(dedx-bethe)*(dedx-bethe)/(sigma*sigma))/sigma;
      mismatch=kFALSE;
    }

    // TODO: Light nuclei, also in TPC pid response
    
    // Check for particles heavier than (AliPID::kSPECIES - 1)
//     if (dedx < (bethe + fRange*sigma)) heavy=kFALSE;

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

  Double_t meanCorrFactor = 0.11/fTOFtail; // Correction factor on the mean because of the tail (should be ~ 0.1 with tail = 1.1)

  // set flat distribution (no decision)
  for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;
  
  if ((track->GetStatus()&AliVTrack::kTOFout)==0) return kDetNoSignal;
  if ((track->GetStatus()&AliVTrack::kTIME)==0) return kDetNoSignal;
  
  Double_t time[AliPID::kSPECIESN];
  track->GetIntegratedTimes(time);
  
  Double_t sigma[AliPID::kSPECIES];
  for (Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++) {
    sigma[iPart] = fTOFResponse.GetExpectedSigma(track->P(),time[iPart],AliPID::ParticleMass(iPart));
  }
  
  Bool_t mismatch = kTRUE/*, heavy = kTRUE*/;
  for (Int_t j=0; j<AliPID::kSPECIES; j++) {
    AliPID::EParticleType type=AliPID::EParticleType(j);
    Double_t nsigmas=NumberOfSigmasTOF(track,type) + meanCorrFactor;

    Double_t sig = sigma[j];
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

    /* OLD Gaussian shape
    if (TMath::Abs(nsigmas) > (fRange+2)) {
      p[j] = TMath::Exp(-0.5*(fRange+2)*(fRange+2))/sig;
    } else
      p[j] = TMath::Exp(-0.5*nsigmas*nsigmas)/sig;
    */

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
AliPIDResponse::EDetPidStatus AliPIDResponse::ComputeTRDProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  //
  // Compute PID response for the
  //

  // set flat distribution (no decision)
  for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;
  if((track->GetStatus()&AliVTrack::kTRDout)==0) return kDetNoSignal;

  Float_t mom[6];
  Double_t dedx[48];  // Allocate space for the maximum number of TRD slices
  Int_t nslices = fTRDslicesForPID[1] - fTRDslicesForPID[0] + 1;
  AliDebug(1, Form("First Slice: %d, Last Slice: %d, Number of slices: %d",  fTRDslicesForPID[0], fTRDslicesForPID[1], nslices));
  for(UInt_t ilayer = 0; ilayer < 6; ilayer++){
    mom[ilayer] = track->GetTRDmomentum(ilayer);
    for(UInt_t islice = fTRDslicesForPID[0]; islice <= fTRDslicesForPID[1]; islice++){
      dedx[ilayer*nslices+islice-fTRDslicesForPID[0]] = track->GetTRDslice(ilayer, islice);
    }
  }
  fTRDResponse.GetResponse(nslices, dedx, mom, p);
  return kDetPidOk;
}
//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::ComputeEMCALProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  //
  // Compute PID response for the EMCAL
  //

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
      if( 999 != fEMCALResponse.ComputeEMCALProbability(pt,EovP,charge,p)){	

  	// in case everything is OK
  	return kDetPidOk;
	
      }
    }
  }
  }
  
  // in all other cases set flat distribution (no decision)
  for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;
  return kDetNoSignal;
  
}
//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::ComputePHOSProbability (const AliVTrack */*track*/, Int_t nSpecies, Double_t p[]) const
{
  //
  // Compute PID response for the PHOS
  //

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

  // set flat distribution (no decision)
  for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;
  if((track->GetStatus()&AliVTrack::kHMPIDpid)==0) return kDetNoSignal;

  track->GetHMPIDpid(p);
  
  return kDetPidOk;
}

//______________________________________________________________________________
void AliPIDResponse::InitialiseEvent(AliVEvent *event, Int_t pass)
{
  //
  // Apply settings for the current event
  //
  fRecoPass=pass;
  
  fCurrentEvent=0x0;
  if (!event) return;
  fCurrentEvent=event;
  fRun=event->GetRunNumber();
  
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
  
  TPRegexp reg(".*(LHC11[a-z]+[0-9]+[a-z_]*)/.*");
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
  else if (fRun>=166529) {
    fLHCperiod="LHC11H";
    fMCperiodTPC="LHC11A10";
    fBeamType="PBPB";
  }


  //exception new pp MC productions from 2011
  if (fBeamType=="PP" && reg.MatchB(fCurrentFile)) fMCperiodTPC="LHC11B2";
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
  //
  //don't load twice for the moment
   if (fArrPidResponseMaster) return;
 

  //reset the PID response functions
  delete fArrPidResponseMaster;
  fArrPidResponseMaster=0x0;
  
  TString fileName(Form("%s/COMMON/PID/data/TPCPIDResponse.root", fOADBPath.Data()));
  
  TFile *f=TFile::Open(fileName.Data());
  if (f && f->IsOpen() && !f->IsZombie()){
    fArrPidResponseMaster=dynamic_cast<TObjArray*>(f->Get("TPCPIDResponse"));
  }
  delete f;
  
  if (!fArrPidResponseMaster){
    AliFatal(Form("Could not retrieve the TPC pid response from: %s",fileName.Data()));
    return;
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
    datatype="MC";
    fRecoPass=1;
  }
  
  //
  //reset old splines
  //
  for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec){
    fTPCResponse.SetResponseFunction((AliPID::EParticleType)ispec,0x0);
  }
  
  //
  //set the new PID splines
  //
  TString period=fLHCperiod;
  if (fArrPidResponseMaster){
    TObject *grAll=0x0;
    //for MC don't use period information
//     if (fIsMC) period="[A-Z0-9]*";
    //for MC use MC period information
    if (fIsMC) period=fMCperiodTPC;
//pattern for the default entry (valid for all particles)
    TPRegexp reg(Form("TSPLINE3_%s_([A-Z]*)_%s_PASS%d_%s_MEAN",datatype.Data(),period.Data(),fRecoPass,fBeamType.Data()));
    
    //loop over entries and filter them
    for (Int_t iresp=0; iresp<fArrPidResponseMaster->GetEntriesFast();++iresp){
      TObject *responseFunction=fArrPidResponseMaster->At(iresp);
      if (responseFunction==0x0) continue;
      TString responseName=responseFunction->GetName();
      
      if (!reg.MatchB(responseName)) continue;
      
      TObjArray *arr=reg.MatchS(responseName);
      TString particleName=arr->At(1)->GetName();
      delete arr;
      if (particleName.IsNull()) continue;
      if (particleName=="ALL") grAll=responseFunction;
      else {
        //find particle id
        for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec){
          TString particle=AliPID::ParticleName(ispec);
          particle.ToUpper();
          if ( particle == particleName ){
            fTPCResponse.SetResponseFunction((AliPID::EParticleType)ispec,responseFunction);
            fTPCResponse.SetUseDatabase(kTRUE);
            AliInfo(Form("Adding graph: %d - %s",ispec,responseFunction->GetName()));
            break;
          }
        }
      }
    }
    
    //set default response function to all particles which don't have a specific one
    if (grAll){
      for (Int_t ispec=0; ispec<AliPID::kSPECIES; ++ispec){
        if (!fTPCResponse.GetResponseFunction((AliPID::EParticleType)ispec)){
          fTPCResponse.SetResponseFunction((AliPID::EParticleType)ispec,grAll);
          AliInfo(Form("Adding graph: %d - %s",ispec,grAll->GetName()));
        }
      }
    }
  }
  
  //
  // Setup resolution parametrisation
  //
  
  //default
  fTPCResponse.SetSigma(3.79301e-03, 2.21280e+04);
  
  if (fRun>=122195){
    fTPCResponse.SetSigma(2.30176e-02, 5.60422e+02);
  }
  if (fArrPidResponseMaster)
  fResolutionCorrection=(TF1*)fArrPidResponseMaster->FindObject(Form("TF1_%s_ALL_%s_PASS%d_%s_SIGMA",datatype.Data(),period.Data(),fRecoPass,fBeamType.Data()));
  
  if (fResolutionCorrection) AliInfo(Form("Setting multiplicity correction function: %s",fResolutionCorrection->GetName()));
}

//______________________________________________________________________________
void AliPIDResponse::SetTRDPidResponseMaster()
{
  //
  // Load the TRD pid params and references from the OADB
  //
  if(fTRDPIDParams) return;
  AliOADBContainer contParams("contParams"); 

  Int_t statusPars = contParams.InitFromFile(Form("%s/COMMON/PID/data/TRDPIDParams.root", fOADBPath.Data()), "AliTRDPIDParams");
  if(statusPars){
    AliError("Failed initializing PID Params from OADB");
  } else {
    AliInfo(Form("Loading TRD Params from %s/COMMON/PID/data/TRDPIDParams.root", fOADBPath.Data()));
    fTRDPIDParams = dynamic_cast<AliTRDPIDParams *>(contParams.GetObject(fRun));
    if(!fTRDPIDParams){
      AliError(Form("TRD Params not found in run %d", fRun));
    }
  }

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
}

//______________________________________________________________________________
void AliPIDResponse::InitializeTRDResponse(){
  //
  // Set PID Params and references to the TRD PID response
  // 
  fTRDResponse.SetPIDParams(fTRDPIDParams);
  fTRDResponse.Load(fTRDPIDReference);
  if(fLHCperiod == "LHC10b" || fLHCperiod == "LHC10c" || fLHCperiod == "LHC10d" || fLHCperiod == "LHC10e"){
    fTRDslicesForPID[0] = 0;
    fTRDslicesForPID[1] = 7;
  }
}

//______________________________________________________________________________
void AliPIDResponse::SetTOFPidResponseMaster()
{
  //
  // Load the TOF pid params from the OADB
  //
  TFile *oadbf = new TFile(Form("%s/COMMON/PID/data/TOFPIDParams.root",fOADBPath.Data()));
  if (oadbf->IsOpen()) {
    AliInfo(Form("Loading TOF Params from %s/COMMON/PID/data/TOFPIDParams.root", fOADBPath.Data()));
    AliOADBContainer *oadbc = (AliOADBContainer *)oadbf->Get("TOFoadb");
    if (fTOFPIDParams) delete fTOFPIDParams;
    fTOFPIDParams = dynamic_cast<AliTOFPIDParams *>(oadbc->GetObject(fRun));
    oadbf->Close();
    delete oadbc;
  } else {
    AliError(Form("TOFPIDParams.root not found in %s/COMMON/PID/data !!",fOADBPath.Data()));
  }
  delete oadbf;

  }

//______________________________________________________________________________
void AliPIDResponse::InitializeTOFResponse(){
  //
  // Set PID Params to the TOF PID response
  // 
  for (Int_t i=0;i<4;i++) {
    fTOFResponse.SetTrackParameter(i,fTOFPIDParams->GetSigParams(i));
  }
  fTOFResponse.SetTimeResolution(fTOFPIDParams->GetTOFresolution());

}


//_________________________________________________________________________
Bool_t AliPIDResponse::IdentifiedAsElectronTRD(const AliVTrack *vtrack, Double_t efficiencyLevel) const {
  //
  // Check whether track is identified as electron under a given electron efficiency hypothesis
  //
  Double_t probs[AliPID::kSPECIES];
  ComputeTRDProbability(vtrack, AliPID::kSPECIES, probs);

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

  return fTRDResponse.IdentifiedAsElectron(ntracklets, probs, p, efficiencyLevel);
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
