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

/* $Id$ */

//-----------------------------------------------------------------
//           Implementation of the combined PID class
//           For the Event Summary Data Class
//           produced by the reconstruction process
//           and containing information on the particle identification
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------

#include "TArrayI.h"
#include "TArrayF.h"

#include "AliLog.h"
#include "AliPID.h"
#include "AliTOFHeader.h"
#include "AliESDpid.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"

ClassImp(AliESDpid)

Int_t AliESDpid::MakePID(AliESDEvent *event, Bool_t TPConly, Float_t timeZeroTOF) const {
  //
  //  Calculate probabilities for all detectors, except if TPConly==kTRUE
  //  and combine PID
  //  
  //   Option TPConly==kTRUE is used during reconstruction, 
  //  because ITS tracking uses TPC pid
  //  HMPID and TRD pid are done in detector reconstructors
  //

  /*
  Float_t timeZeroTOF = 0;
  if (subtractT0) 
    timeZeroTOF = event->GetT0();
  */
  Int_t nTrk=event->GetNumberOfTracks();
  for (Int_t iTrk=0; iTrk<nTrk; iTrk++) {  
    AliESDtrack *track=event->GetTrack(iTrk);
    MakeTPCPID(track);
    if (!TPConly) {
      MakeITSPID(track);
      MakeTOFPID(track, timeZeroTOF);
      //MakeHMPIDPID(track);
      //MakeTRDPID(track);
    }
    CombinePID(track);
  }
  return 0;
}
//_________________________________________________________________________
void AliESDpid::MakeTPCPID(AliESDtrack *track) const
{
  //
  //  TPC pid using bethe-bloch and gaussian response
  //
  if ((track->GetStatus()&AliESDtrack::kTPCin )==0)
    if ((track->GetStatus()&AliESDtrack::kTPCout)==0) return;

    Double_t mom = track->GetP();
    const AliExternalTrackParam *in=track->GetInnerParam();
    if (in) mom = in->GetP();

    Double_t p[AliPID::kSPECIES];
    Double_t dedx=track->GetTPCsignal(); 
    Bool_t mismatch=kTRUE, heavy=kTRUE;

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

      // Check for particles heavier than (AliPID::kSPECIES - 1)
      if (dedx < (bethe + fRange*sigma)) heavy=kFALSE;

    }

    if (mismatch)
       for (Int_t j=0; j<AliPID::kSPECIES; j++) p[j]=1./AliPID::kSPECIES;

    track->SetTPCpid(p);

    if (heavy) track->ResetStatus(AliESDtrack::kTPCpid);

}
//_________________________________________________________________________
void AliESDpid::MakeITSPID(AliESDtrack *track) const
{
  //
  // ITS PID
  // Two options, depending on fITSPIDmethod:
  //  1) Truncated mean method
  //  2) Likelihood, using charges measured in all 4 layers and 
  //     Landau+gaus response functions
  //

  if ((track->GetStatus()&AliESDtrack::kITSin)==0 &&
      (track->GetStatus()&AliESDtrack::kITSout)==0) return;

  Double_t mom=track->GetP();  
  if (fITSPIDmethod == kITSTruncMean) {
    Double_t dedx=track->GetITSsignal();
    Bool_t isSA=kTRUE;
    Double_t momITS=mom;
    ULong_t trStatus=track->GetStatus();
    if(trStatus&AliESDtrack::kTPCin) isSA=kFALSE;
    UChar_t clumap=track->GetITSClusterMap();
    Int_t nPointsForPid=0;
    for(Int_t i=2; i<6; i++){
      if(clumap&(1<<i)) ++nPointsForPid;
    }

    if(nPointsForPid<3) { // track not to be used for combined PID purposes
      track->ResetStatus(AliESDtrack::kITSpid);
      return;
    }

    Double_t p[10];

    Bool_t mismatch=kTRUE, heavy=kTRUE;
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
      if (dedx < (bethe + fRange*sigma)) heavy=kFALSE;

    }

    if (mismatch)
       for (Int_t j=0; j<AliPID::kSPECIES; j++) p[j]=1./AliPID::kSPECIES;

    track->SetITSpid(p);

    if (heavy) track->ResetStatus(AliESDtrack::kITSpid);
  }
  else {  // Likelihood method
    Double_t condprobfun[AliPID::kSPECIES];
    Double_t qclu[4];
    track->GetITSdEdxSamples(qclu);
    fITSResponse.GetITSProbabilities(mom,qclu,condprobfun);
    track->SetITSpid(condprobfun);
  }

}
//_________________________________________________________________________
void AliESDpid::MakeTOFPID(AliESDtrack *track, Float_t /*timeZeroTOF*/) const
{
  //
  //   TOF PID using gaussian response
  //

  if ((track->GetStatus()&AliESDtrack::kTOFout)==0) return;
  if ((track->GetStatus()&AliESDtrack::kTIME)==0) return;

  Int_t ibin = fTOFResponse.GetMomBin(track->GetP());
  Float_t timezero = fTOFResponse.GetT0bin(ibin);

  Double_t time[AliPID::kSPECIESN];
  track->GetIntegratedTimes(time);

  Double_t sigma[AliPID::kSPECIES];
  for (Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++) {
    sigma[iPart] = fTOFResponse.GetExpectedSigma(track->GetP(),time[iPart],AliPID::ParticleMass(iPart));
  }

  AliDebugGeneral("AliESDpid::MakeTOFPID",2,
	   Form("Expected TOF signals [ps]: %f %f %f %f %f",
		  time[AliPID::kElectron],
		  time[AliPID::kMuon],
		  time[AliPID::kPion],
		  time[AliPID::kKaon],
		  time[AliPID::kProton]));

  AliDebugGeneral("AliESDpid::MakeTOFPID",2,
	   Form("Expected TOF std deviations [ps]: %f %f %f %f %f",
		  sigma[AliPID::kElectron],
		  sigma[AliPID::kMuon],
		  sigma[AliPID::kPion],
		  sigma[AliPID::kKaon],
		  sigma[AliPID::kProton]
		  ));

  Double_t tof = track->GetTOFsignal() - timezero;

  Double_t p[AliPID::kSPECIES];
  Bool_t mismatch = kTRUE, heavy = kTRUE;
  for (Int_t j=0; j<AliPID::kSPECIES; j++) {
    Double_t sig = sigma[j];
    if (TMath::Abs(tof-time[j]) > (fRange+2)*sig) {
	p[j] = TMath::Exp(-0.5*(fRange+2)*(fRange+2))/sig;
    } else
      p[j] = TMath::Exp(-0.5*(tof-time[j])*(tof-time[j])/(sig*sig))/sig;

    // Check the mismatching
    Double_t mass = AliPID::ParticleMass(j);
    Double_t pm = fTOFResponse.GetMismatchProbability(track->GetP(),mass);
    if (p[j]>pm) mismatch = kFALSE;

    // Check for particles heavier than (AliPID::kSPECIES - 1)
    if (tof < (time[j] + fRange*sig)) heavy=kFALSE;

  }

  if (mismatch)
    for (Int_t j=0; j<AliPID::kSPECIES; j++) p[j]=1./AliPID::kSPECIES;

  track->SetTOFpid(p);

  if (heavy) track->ResetStatus(AliESDtrack::kTOFpid);    
  if (!CheckTOFMatching(track)) track->SetStatus(AliESDtrack::kTOFmismatch);
  
}
//_________________________________________________________________________
void AliESDpid::MakeTRDPID(AliESDtrack *track) const
{
  //
  // Method to recalculate the TRD PID probabilities
  //
  if((track->GetStatus()&AliESDtrack::kTRDout)==0) return;
  Double_t prob[AliPID::kSPECIES]; Float_t mom[6];
  Double_t dedx[48];  // Allocate space for the maximum number of TRD slices
  for(Int_t ilayer = 0; ilayer < 6; ilayer++){
    mom[ilayer] = track->GetTRDmomentum(ilayer);
    for(Int_t islice = 0; islice < track->GetNumberOfTRDslices(); islice++){
      dedx[ilayer*track->GetNumberOfTRDslices()+islice] = track->GetTRDslice(ilayer, islice);
    }
  }
  fTRDResponse.GetResponse(track->GetNumberOfTRDslices(), dedx, mom, prob);
  track->SetTRDpid(prob);
}
//_________________________________________________________________________
void AliESDpid::CombinePID(AliESDtrack *track) const
{
  //
  // Combine the information of various detectors
  // to determine the Particle Identification
  //
  Int_t ns=AliPID::kSPECIES;
  Double_t p[10]={1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};

  if (track->IsOn(AliESDtrack::kITSpid)) {
    Double_t d[10];
    track->GetITSpid(d);
    for (Int_t j=0; j<ns; j++) p[j]*=d[j];
  }

  if (track->IsOn(AliESDtrack::kTPCpid)) {
    Double_t d[10];
    track->GetTPCpid(d);
    for (Int_t j=0; j<ns; j++) p[j]*=d[j];
  }

  if (track->IsOn(AliESDtrack::kTRDpid)) {
    Double_t d[10];
    track->GetTRDpid(d);
    for (Int_t j=0; j<ns; j++) p[j]*=d[j];
  }

  if (track->IsOn(AliESDtrack::kTOFpid)) {
    Double_t d[10];
    track->GetTOFpid(d);
    for (Int_t j=0; j<ns; j++) p[j]*=d[j];
  }

  if (track->IsOn(AliESDtrack::kHMPIDpid)) {
    Double_t d[10];
    track->GetHMPIDpid(d);
    for (Int_t j=0; j<ns; j++) p[j]*=d[j];
  }

  track->SetESDpid(p);
}
//_________________________________________________________________________
Bool_t AliESDpid::CheckTOFMatching(AliESDtrack *track) const{
    Bool_t status = kFALSE;
    
    Double_t exptimes[5];
    track->GetIntegratedTimes(exptimes);
    
    Float_t p = track->P();
    
    Float_t dedx = track->GetTPCsignal();
    Float_t time = track->GetTOFsignal() - fTOFResponse.GetStartTime(p);
    
    Double_t ptpc[3];
    track->GetInnerPxPyPz(ptpc);
    Float_t momtpc=TMath::Sqrt(ptpc[0]*ptpc[0] + ptpc[1]*ptpc[1] + ptpc[2]*ptpc[2]);
    
    for(Int_t i=0;i < 5;i++){
	AliPID::EParticleType type=AliPID::EParticleType(i);
	
	Float_t resolutionTOF = fTOFResponse.GetExpectedSigma(p, exptimes[i], AliPID::ParticleMass(i));
	if(TMath::Abs(exptimes[i] - time) < fRange * resolutionTOF){
	    Float_t dedxExp = fTPCResponse.GetExpectedSignal(momtpc,type);
	    Float_t resolutionTPC = fTPCResponse.GetExpectedSigma(momtpc,track->GetTPCsignalN(),type);
	    
	    if(TMath::Abs(dedx - dedxExp) < fRangeTOFMismatch * resolutionTPC){
		status = kTRUE;
	    }
	}
    }
    
    // for nuclei
    Float_t resolutionTOFpr = fTOFResponse.GetExpectedSigma(p, exptimes[4], AliPID::ParticleMass(4));
    if(!status && (exptimes[4] + fRange*resolutionTOFpr < time)) status = kTRUE;
    
    
    return status;
}
//_________________________________________________________________________
void AliESDpid::SetTOFResponse(AliVEvent *vevent,EStartTimeType_t option){
  //
  // Set TOF response function
  // Input option for event_time used
  //

    AliESDEvent *event=(AliESDEvent*)vevent;
  
    Float_t t0spread = 0.; //event->GetEventTimeSpread();
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
    if(event->GetT0TOF()){ // check if T0 detector information is available
	flagT0T0=kTRUE;
    }


    AliTOFHeader *tofHeader =(AliTOFHeader*)event->GetTOFHeader();

    if(tofHeader){ // read global info and T0-TOF info from ESD
      fTOFResponse.SetTimeResolution(tofHeader->GetTOFResolution());
      t0spread = tofHeader->GetT0spread(); // read t0 sprad
      if(t0spread < 10) t0spread = 80;

      flagT0TOF=kTRUE;
      for(Int_t i=0;i<fTOFResponse.GetNmomBins();i++){ // read T0-TOF default value
	startTime[i]=tofHeader->GetDefaultEventTimeVal();
	startTimeRes[i]=tofHeader->GetDefaultEventTimeRes();
      }

      TArrayI *ibin=tofHeader->GetNvalues();
      TArrayF *t0Bin=tofHeader->GetEventTimeValues();
      TArrayF *t0ResBin=tofHeader->GetEventTimeRes();
      for(Int_t j=0;j < tofHeader->GetNbins();j++){ // fill T0-TOF in p-bins
	Int_t icurrent = (Int_t)ibin->GetAt(j);
	startTime[icurrent]=t0Bin->GetAt(j);
	startTimeRes[icurrent]=t0ResBin->GetAt(j);
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
	    t0AC= event->GetT0TOF()[0];
	    t0A= event->GetT0TOF()[1];
	    t0C= event->GetT0TOF()[2];
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

    else if(option == kT0_T0){ // T0-T0 is used when available (T0-FILL otherwise) from ESD
	Float_t t0AC=-10000;
	Float_t t0A=-10000;
	Float_t t0C=-10000;
	if(flagT0T0){
	    t0AC= event->GetT0TOF()[0];
	    t0A= event->GetT0TOF()[1];
	    t0C= event->GetT0TOF()[2];
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
