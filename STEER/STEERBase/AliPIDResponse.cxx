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
#include <TH2D.h>
#include <TSpline.h>
#include <TFile.h>
#include <TArrayI.h>
#include <TArrayF.h>
#include <TLinearFitter.h>

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
fCachePID(kTRUE),
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
fResT0A(75.),
fResT0C(65.),
fResT0AC(55.),
fArrPidResponseMaster(NULL),
fResolutionCorrection(NULL),
fOADBvoltageMaps(NULL),
fUseTPCEtaCorrection(kFALSE),//TODO: In future, default kTRUE
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
fCachePID(other.fCachePID),
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
fResT0A(75.),
fResT0C(65.),
fResT0AC(55.),
fArrPidResponseMaster(NULL),
fResolutionCorrection(NULL),
fOADBvoltageMaps(NULL),
fUseTPCEtaCorrection(other.fUseTPCEtaCorrection),
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
    fCachePID=other.fCachePID;
    fBeamType="PP";
    fLHCperiod="";
    fMCperiodTPC="";
    fMCperiodUser=other.fMCperiodUser;
    fCurrentFile="";
    fRecoPass=0;
    fRecoPassUser=other.fRecoPassUser;
    fRun=0;
    fOldRun=0;
    fResT0A=75.;
    fResT0C=65.;
    fResT0AC=55.;
    fArrPidResponseMaster=NULL;
    fResolutionCorrection=NULL;
    fOADBvoltageMaps=NULL;
	fUseTPCEtaCorrection=other.fUseTPCEtaCorrection;
    fTRDPIDResponseObject=NULL;
    fEMCALPIDParams=NULL;
    fTOFtail=1.1;
    fTOFPIDParams=NULL;
    fCurrentEvent=other.fCurrentEvent;

  }
  return *this;
}

//______________________________________________________________________________
Float_t AliPIDResponse::NumberOfSigmas(EDetector detector, const AliVParticle *vtrack, AliPID::EParticleType type) const
{
  //
  // NumberOfSigmas for 'detCode'
  //
  
  const AliVTrack *track=static_cast<const AliVTrack*>(vtrack);
  // look for cached value first
  const AliDetectorPID *detPID=track->GetDetectorPID();
  
  if ( detPID && detPID->HasNumberOfSigmas(detector)){
    return detPID->GetNumberOfSigmas(detector, type);
  } else if (fCachePID) {
    FillTrackDetectorPID(track, detector);
    detPID=track->GetDetectorPID();
    return detPID->GetNumberOfSigmas(detector, type);
  }
  
  return GetNumberOfSigmas(detector, track, type);
}

//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::NumberOfSigmas(EDetector detCode, const AliVParticle *track,
                                                             AliPID::EParticleType type, Double_t &val) const
{
  //
  // NumberOfSigmas with detector status as return value
  //
  
  val=NumberOfSigmas(detCode, track, type);
  return CheckPIDStatus(detCode, (AliVTrack*)track);
}

//______________________________________________________________________________
// public buffered versions of the PID calculation
//

//______________________________________________________________________________
Float_t AliPIDResponse::NumberOfSigmasITS(const AliVParticle *vtrack, AliPID::EParticleType type) const
{
  //
  // Calculate the number of sigmas in the ITS
  //
  
  return NumberOfSigmas(kITS, vtrack, type);
}

//______________________________________________________________________________
Float_t AliPIDResponse::NumberOfSigmasTPC(const AliVParticle *vtrack, AliPID::EParticleType type) const
{
  //
  // Calculate the number of sigmas in the TPC
  //
  
  return NumberOfSigmas(kTPC, vtrack, type);
}

//______________________________________________________________________________
Float_t AliPIDResponse::NumberOfSigmasTPC( const AliVParticle *vtrack, 
                                           AliPID::EParticleType type,
                                           AliTPCPIDResponse::ETPCdEdxSource dedxSource) const
{
  //get number of sigmas according the selected TPC gain configuration scenario
  const AliVTrack *track=static_cast<const AliVTrack*>(vtrack);

//   return 0.;
  Float_t nSigma=fTPCResponse.GetNumberOfSigmas(track, type, dedxSource, fUseTPCEtaCorrection);

  return nSigma;
}

//______________________________________________________________________________
Float_t AliPIDResponse::NumberOfSigmasTOF(const AliVParticle *vtrack, AliPID::EParticleType type) const
{
  //
  // Calculate the number of sigmas in the TOF
  //
  
  return NumberOfSigmas(kTOF, vtrack, type);
}

//______________________________________________________________________________
Float_t AliPIDResponse::NumberOfSigmasEMCAL(const AliVParticle *vtrack, AliPID::EParticleType type) const
{
  //
  // Calculate the number of sigmas in the EMCAL
  //
  
  return NumberOfSigmas(kEMCAL, vtrack, type);
}

//______________________________________________________________________________
Float_t  AliPIDResponse::NumberOfSigmasEMCAL(const AliVParticle *vtrack, AliPID::EParticleType type, Double_t &eop, Double_t showershape[4])  const
{
  //
  // emcal nsigma with eop and showershape
  //
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

        // look for cached value first
        const AliDetectorPID *detPID=track->GetDetectorPID();
        const EDetector detector=kEMCAL;
        
        if ( detPID && detPID->HasNumberOfSigmas(detector)){
          return detPID->GetNumberOfSigmas(detector, type);
        } else if (fCachePID) {
          FillTrackDetectorPID(track, detector);
          detPID=track->GetDetectorPID();
          return detPID->GetNumberOfSigmas(detector, type);
        }
        
        // NSigma value really meaningful only for electrons!
        return fEMCALResponse.GetNumberOfSigmas(pt,EovP,type,charge);
      }
    }
  }
  return -999;
}

//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::ComputePIDProbability  (EDetCode  detCode, const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  // Compute PID response of 'detCode'
  
  // find detector code from detector bit mask
  Int_t detector=-1;
  for (Int_t idet=0; idet<kNdetectors; ++idet) if ( (detCode&(1<<idet)) ) { detector=idet; break; }
  if (detector==-1) return kDetNoSignal;

  return ComputePIDProbability((EDetector)detector, track, nSpecies, p);
}

//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::ComputePIDProbability  (EDetector detector,  const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  //
  // Compute PID response of 'detector'
  //

  const AliDetectorPID *detPID=track->GetDetectorPID();

  if ( detPID && detPID->HasRawProbability(detector)){
    return detPID->GetRawProbability(detector, p, nSpecies);
  } else if (fCachePID) {
    FillTrackDetectorPID(track, detector);
    detPID=track->GetDetectorPID();
    return detPID->GetRawProbability(detector, p, nSpecies);
  }
  
  //if no caching return values calculated from scratch
  return GetComputePIDProbability(detector, track, nSpecies, p);
}

//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::ComputeITSProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  // Compute PID response for the ITS
  return ComputePIDProbability(kITS, track, nSpecies, p);
}

//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::ComputeTPCProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  // Compute PID response for the TPC
  return ComputePIDProbability(kTPC, track, nSpecies, p);
}

//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::ComputeTOFProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  // Compute PID response for the
  return ComputePIDProbability(kTOF, track, nSpecies, p);
}

//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::ComputeTRDProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  // Compute PID response for the
  return ComputePIDProbability(kTRD, track, nSpecies, p);
}

//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::ComputeEMCALProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  // Compute PID response for the EMCAL
  return ComputePIDProbability(kEMCAL, track, nSpecies, p);
}
//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::ComputePHOSProbability (const AliVTrack */*track*/, Int_t nSpecies, Double_t p[]) const
{
  // Compute PID response for the PHOS
  
  // set flat distribution (no decision)
  for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;
  return kDetNoSignal;
}

//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::ComputeHMPIDProbability(const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  // Compute PID response for the HMPID
  return ComputePIDProbability(kHMPID, track, nSpecies, p);
}

//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::ComputeTRDProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[],AliTRDPIDResponse::ETRDPIDMethod PIDmethod) const
{
  // Compute PID response for the
  return GetComputeTRDProbability(track, nSpecies, p, PIDmethod);
}

//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::CheckPIDStatus(EDetector detector, const AliVTrack *track) const
{
  // calculate detector pid status
  
  const Int_t iDetCode=(Int_t)detector;
  if (iDetCode<0||iDetCode>=kNdetectors) return kDetNoSignal;
  const AliDetectorPID *detPID=track->GetDetectorPID();
  
  if ( detPID ){
    return detPID->GetPIDStatus(detector);
  } else if (fCachePID) {
    FillTrackDetectorPID(track, detector);
    detPID=track->GetDetectorPID();
    return detPID->GetPIDStatus(detector);
  }
  
  // if not buffered and no buffering is requested
  return GetPIDStatus(detector, track);
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
  SetTPCEtaMaps();

  SetTRDPidResponseMaster(); 
  InitializeTRDResponse();

  SetEMCALPidResponseMaster(); 
  InitializeEMCALResponse();
  
  SetTOFPidResponseMaster();
  InitializeTOFResponse();

  if (fCurrentEvent) fTPCResponse.SetMagField(fCurrentEvent->GetMagneticField());
}

//______________________________________________________________________________
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
  TPRegexp reg12a17("LHC1[2-3][a-z]");

  //find the period by run number (UGLY, but not stored in ESD and AOD... )
  if (fRun>=114737&&fRun<=117223)      { fLHCperiod="LHC10B"; fMCperiodTPC="LHC10D1";  }
  else if (fRun>=118503&&fRun<=121040) { fLHCperiod="LHC10C"; fMCperiodTPC="LHC10D1";  }
  else if (fRun>=122195&&fRun<=126437) { fLHCperiod="LHC10D"; fMCperiodTPC="LHC10F6A"; }
  else if (fRun>=127710&&fRun<=130850) { fLHCperiod="LHC10E"; fMCperiodTPC="LHC10F6A"; }
  else if (fRun>=133004&&fRun<=135029) { fLHCperiod="LHC10F"; fMCperiodTPC="LHC10F6A"; }
  else if (fRun>=135654&&fRun<=136377) { fLHCperiod="LHC10G"; fMCperiodTPC="LHC10F6A"; }
  else if (fRun>=136851&&fRun<=139846) {
    fLHCperiod="LHC10H";
    fMCperiodTPC="LHC10H8";
    if (reg.MatchB(fCurrentFile)) fMCperiodTPC="LHC11A10";
    fBeamType="PBPB";
  }
  else if (fRun>=139847&&fRun<=146974) { fLHCperiod="LHC11A"; fMCperiodTPC="LHC10F6A"; }
  //TODO: periods 11B (146975-150721), 11C (150722-155837) are not yet treated assume 11d for the moment
  else if (fRun>=146975&&fRun<=155837) { fLHCperiod="LHC11D"; fMCperiodTPC="LHC10F6A"; }
  else if (fRun>=155838&&fRun<=159649) { fLHCperiod="LHC11D"; fMCperiodTPC="LHC10F6A"; }
  // also for 11e (159650-162750),f(162751-165771) use 11d
  else if (fRun>=159650&&fRun<=162750) { fLHCperiod="LHC11D"; fMCperiodTPC="LHC10F6A"; }
  else if (fRun>=162751&&fRun<=165771) { fLHCperiod="LHC11D"; fMCperiodTPC="LHC10F6A"; }
  
  else if (fRun>=165772 && fRun<=170718) {
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
  if (fRun >= 186636  ) { fLHCperiod="LHC12G"; fBeamType="PPB"; fMCperiodTPC="LHC12G"; }

  //exception new pp MC productions from 2011
  if (fBeamType=="PP" && reg.MatchB(fCurrentFile)) { fMCperiodTPC="LHC11B2"; fBeamType="PP"; }
  // exception for 11f1
  if (fCurrentFile.Contains("LHC11f1/")) fMCperiodTPC="LHC11F1";
  // exception for 12f1a, 12f1b and 12i3
  if (fCurrentFile.Contains("LHC12f1a/") || fCurrentFile.Contains("LHC12f1b/")
      || fCurrentFile.Contains("LHC12i3/")) fMCperiodTPC="LHC12F1";
  // exception for 12c4
  if (fCurrentFile.Contains("LHC12c4/")) fMCperiodTPC="LHC12C4";
}

//______________________________________________________________________________
void AliPIDResponse::SetITSParametrisation()
{
  //
  // Set the ITS parametrisation
  //
}

 
//______________________________________________________________________________
void AliPIDResponse::AddPointToHyperplane(TH2D* h, TLinearFitter* linExtrapolation, Int_t binX, Int_t binY)
{
  if (h->GetBinContent(binX, binY) <= 1e-4)
    return; // Reject bins without content (within some numerical precision) or with strange content
    
  Double_t coord[2] = {0, 0};
  coord[0] = h->GetXaxis()->GetBinCenter(binX);
  coord[1] = h->GetYaxis()->GetBinCenter(binY);
  Double_t binError = h->GetBinError(binX, binY);
  if (binError <= 0) {
    binError = 1000; // Should not happen because bins without content are rejected for the map (TH2D* h)
    printf("ERROR: This should never happen: Trying to add bin in addPointToHyperplane with error not set....\n");
  }
  linExtrapolation->AddPoint(coord, h->GetBinContent(binX, binY, binError));
}


//______________________________________________________________________________
TH2D* AliPIDResponse::RefineHistoViaLinearInterpolation(TH2D* h, Double_t refineFactorX, Double_t refineFactorY)
{
  if (!h)
    return 0x0;
  
  // Interpolate to finer map
  TLinearFitter* linExtrapolation = new TLinearFitter(2, "hyp2", "");
  
  Double_t upperMapBoundY = h->GetYaxis()->GetBinUpEdge(h->GetYaxis()->GetNbins());
  Double_t lowerMapBoundY = h->GetYaxis()->GetBinLowEdge(1);
  Int_t nBinsX = 30;
  // Binning was find to yield good results, if 40 bins are chosen for the range 0.0016 to 0.02. For the new variable range,
  // scale the number of bins correspondingly
  Int_t nBinsY = TMath::Nint((upperMapBoundY - lowerMapBoundY) / (0.02 - 0.0016) * 40);
  Int_t nBinsXrefined = nBinsX * refineFactorX;
  Int_t nBinsYrefined = nBinsY * refineFactorY; 
  
  TH2D* hRefined = new TH2D(Form("%s_refined", h->GetName()),  Form("%s (refined)", h->GetTitle()),
                            nBinsXrefined, h->GetXaxis()->GetBinLowEdge(1), h->GetXaxis()->GetBinUpEdge(h->GetXaxis()->GetNbins()),
                            nBinsYrefined, lowerMapBoundY, upperMapBoundY);
  
  for (Int_t binX = 1; binX <= nBinsXrefined; binX++)  {
    for (Int_t binY = 1; binY <= nBinsYrefined; binY++)  {
      
      hRefined->SetBinContent(binX, binY, 1); // Default value is 1
      
      Double_t centerX = hRefined->GetXaxis()->GetBinCenter(binX);
      Double_t centerY = hRefined->GetYaxis()->GetBinCenter(binY);
      
      /*OLD
      linExtrapolation->ClearPoints();
      
      // For interpolation: Just take the corresponding bin from the old histo.
      // For extrapolation: take the last available bin from the old histo.
      // If the boundaries are to be skipped, also skip the corresponding bins
      Int_t oldBinX = h->GetXaxis()->FindBin(centerX);
      if (oldBinX < 1)  
        oldBinX = 1;
      if (oldBinX > nBinsX)
        oldBinX = nBinsX;
      
      Int_t oldBinY = h->GetYaxis()->FindBin(centerY);
      if (oldBinY < 1)  
        oldBinY = 1;
      if (oldBinY > nBinsY)
        oldBinY = nBinsY;
      
      // Neighbours left column
      if (oldBinX >= 2) {
        if (oldBinY >= 2) {
          AddPointToHyperplane(h, linExtrapolation, oldBinX - 1, oldBinY - 1);
        }
        
        AddPointToHyperplane(h, linExtrapolation, oldBinX - 1, oldBinY);
        
        if (oldBinY < nBinsY) {
          AddPointToHyperplane(h, linExtrapolation, oldBinX - 1, oldBinY + 1);
        }
      }
      
      // Neighbours (and point itself) same column
      if (oldBinY >= 2) {
        AddPointToHyperplane(h, linExtrapolation, oldBinX, oldBinY - 1);
      }
        
      AddPointToHyperplane(h, linExtrapolation, oldBinX, oldBinY);
        
      if (oldBinY < nBinsY) {
        AddPointToHyperplane(h, linExtrapolation, oldBinX, oldBinY + 1);
      }
      
      // Neighbours right column
      if (oldBinX < nBinsX) {
        if (oldBinY >= 2) {
          AddPointToHyperplane(h, linExtrapolation, oldBinX + 1, oldBinY - 1);
        }
        
        AddPointToHyperplane(h, linExtrapolation, oldBinX + 1, oldBinY);
        
        if (oldBinY < nBinsY) {
          AddPointToHyperplane(h, linExtrapolation, oldBinX + 1, oldBinY + 1);
        }
      }
      
      
      // Fit 2D-hyperplane
      if (linExtrapolation->GetNpoints() <= 0)
        continue;
        
      if (linExtrapolation->Eval() != 0)// EvalRobust -> Takes much, much, [...], much more time (~hours instead of seconds)
        continue;
      
      // Fill the bin of the refined histogram with the extrapolated value
      Double_t interpolatedValue = linExtrapolation->GetParameter(0) + linExtrapolation->GetParameter(1) * centerX
                                 + linExtrapolation->GetParameter(2) * centerY;
      */
      Double_t interpolatedValue = h->Interpolate(centerX, centerY) ;
      hRefined->SetBinContent(binX, binY, interpolatedValue);      
    }
  } 
  
  
  // Problem: Interpolation does not work before/beyond center of first/last bin (as the name suggests).
  // Therefore, for each row in dEdx: Take last bin from old map and interpolate values from center and edge.
  // Assume line through these points and extropolate to last bin of refined map
  const Double_t firstOldXbinUpEdge = h->GetXaxis()->GetBinUpEdge(1);
  const Double_t firstOldXbinCenter = h->GetXaxis()->GetBinCenter(1);
  
  const Double_t oldXbinHalfWidth = firstOldXbinUpEdge - firstOldXbinCenter;
  
  const Double_t lastOldXbinLowEdge = h->GetXaxis()->GetBinLowEdge(h->GetNbinsX());
  const Double_t lastOldXbinCenter = h->GetXaxis()->GetBinCenter(h->GetNbinsX());
  
  for (Int_t binY = 1; binY <= nBinsYrefined; binY++)  {
    Double_t centerY = hRefined->GetYaxis()->GetBinCenter(binY);
    
    const Double_t interpolatedCenterFirstXbin = h->Interpolate(firstOldXbinCenter, centerY);
    const Double_t interpolatedUpEdgeFirstXbin = h->Interpolate(firstOldXbinUpEdge, centerY);
    
    const Double_t extrapolationSlopeFirstXbin = (interpolatedUpEdgeFirstXbin - interpolatedCenterFirstXbin) / oldXbinHalfWidth;
    const Double_t extrapolationOffsetFirstXbin = interpolatedCenterFirstXbin;
    
    
    const Double_t interpolatedCenterLastXbin = h->Interpolate(lastOldXbinCenter, centerY);
    const Double_t interpolatedLowEdgeLastXbin = h->Interpolate(lastOldXbinLowEdge, centerY);
    
    const Double_t extrapolationSlopeLastXbin = (interpolatedCenterLastXbin - interpolatedLowEdgeLastXbin) / oldXbinHalfWidth;
    const Double_t extrapolationOffsetLastXbin = interpolatedCenterLastXbin;

    for (Int_t binX = 1; binX <= nBinsXrefined; binX++)  {
      Double_t centerX = hRefined->GetXaxis()->GetBinCenter(binX);
     
      if (centerX < firstOldXbinCenter) {
        Double_t extrapolatedValue = extrapolationOffsetFirstXbin + (centerX - firstOldXbinCenter) * extrapolationSlopeFirstXbin;
        hRefined->SetBinContent(binX, binY, extrapolatedValue);      
      }
      else if (centerX <= lastOldXbinCenter) {
        continue;
      }
      else {
        Double_t extrapolatedValue = extrapolationOffsetLastXbin + (centerX - lastOldXbinCenter) * extrapolationSlopeLastXbin;
        hRefined->SetBinContent(binX, binY, extrapolatedValue);     
      }
    }
  } 
  
  delete linExtrapolation;
  
  return hRefined;
}

//______________________________________________________________________________
void AliPIDResponse::SetTPCEtaMaps(Double_t refineFactorMapX, Double_t refineFactorMapY,
                                   Double_t refineFactorSigmaMapX, Double_t refineFactorSigmaMapY)
{
  //
  // Load the TPC eta correction maps from the OADB
  //
  
  if (fUseTPCEtaCorrection == kFALSE) {
    // Disable eta correction via setting no maps
    if (!fTPCResponse.SetEtaCorrMap(0x0))
      AliInfo("Request to disable TPC eta correction -> Eta correction has been disabled"); 
    else
      AliError("Request to disable TPC eta correction -> Some error occured when unloading the correction maps");
    
    if (!fTPCResponse.SetSigmaParams(0x0, 0))
      AliInfo("Request to disable TPC eta correction -> Using old parametrisation for sigma"); 
    else
      AliError("Request to disable TPC eta correction -> Some error occured when unloading the sigma maps");
    
    return;
  }
  
  TString dataType = "DATA";
  TString period = fLHCperiod.IsNull() ? "No period information" : fLHCperiod;
  
  if (fIsMC)  {
    if (!fTuneMConData) {
      period=fMCperiodTPC;
      dataType="MC";
    }
    fRecoPass = 1;
    
    if (!fTuneMConData && fMCperiodTPC.IsNull()) {
      AliFatal("MC detected, but no MC period set -> Not changing eta maps!");
      return;
    }
  }

  Int_t recopass = fRecoPass;
  if (fTuneMConData)
    recopass = fRecoPassUser;
  
  TString defaultObj = Form("Default_%s_pass%d", dataType.Data(), recopass);
  
  AliInfo(Form("Current period and reco pass: %s.pass%d", period.Data(), recopass));
  
  // Invalidate old maps
  fTPCResponse.SetEtaCorrMap(0x0);
  fTPCResponse.SetSigmaParams(0x0, 0);
  
  // Load the eta correction maps
  AliOADBContainer etaMapsCont(Form("TPCetaMaps_%s_pass%d", dataType.Data(), recopass)); 
  
  Int_t statusCont = etaMapsCont.InitFromFile(Form("%s/COMMON/PID/data/TPCetaMaps.root", fOADBPath.Data()),
                                              Form("TPCetaMaps_%s_pass%d", dataType.Data(), recopass));
  if (statusCont) {
    AliError("Failed initializing TPC eta correction maps from OADB -> Disabled eta correction");
  }
  else {
    AliInfo(Form("Loading TPC eta correction map from %s/COMMON/PID/data/TPCetaMaps.root", fOADBPath.Data()));
    
    TH2D* etaMap = 0x0;
    
    if (fIsMC && !fTuneMConData) {
      TString searchMap = Form("TPCetaMaps_%s_%s_pass%d", dataType.Data(), period.Data(), recopass);
      etaMap = dynamic_cast<TH2D *>(etaMapsCont.GetDefaultObject(searchMap.Data()));
      if (!etaMap) {
        // Try default object
        etaMap = dynamic_cast<TH2D *>(etaMapsCont.GetDefaultObject(defaultObj.Data()));
      }
    }
    else {
      etaMap = dynamic_cast<TH2D *>(etaMapsCont.GetObject(fRun, defaultObj.Data()));
    }
    
        
    if (!etaMap) {
      AliError(Form("TPC eta correction map not found for run %d and also no default map found -> Disabled eta correction!!!", fRun));
    }
    else {
      TH2D* etaMapRefined = RefineHistoViaLinearInterpolation(etaMap, refineFactorMapX, refineFactorMapY);
      
      if (etaMapRefined) {
        if (!fTPCResponse.SetEtaCorrMap(etaMapRefined)) {
          AliError(Form("Failed to set TPC eta correction map for run %d -> Disabled eta correction!!!", fRun));
          fTPCResponse.SetEtaCorrMap(0x0);
        }
        else {
          AliInfo(Form("Loaded TPC eta correction map (refine factors %.2f/%.2f) from %s/COMMON/PID/data/TPCetaMaps.root: %s", 
                       refineFactorMapX, refineFactorMapY, fOADBPath.Data(), fTPCResponse.GetEtaCorrMap()->GetTitle()));
        }
        
        delete etaMapRefined;
      }
      else {
        AliError(Form("Failed to set TPC eta correction map for run %d (map was loaded, but couldn't be refined) -> Disabled eta correction!!!", fRun));
      }
    }
  }
  
  // Load the sigma parametrisation (1/dEdx vs tanTheta_local (~eta))
  AliOADBContainer etaSigmaMapsCont(Form("TPCetaSigmaMaps_%s_pass%d", dataType.Data(), recopass)); 
  
  statusCont = etaSigmaMapsCont.InitFromFile(Form("%s/COMMON/PID/data/TPCetaMaps.root", fOADBPath.Data()),
                                             Form("TPCetaSigmaMaps_%s_pass%d", dataType.Data(), recopass));
  if (statusCont) {
    AliError("Failed initializing TPC eta sigma maps from OADB -> Using old sigma parametrisation");
  }
  else {
    AliInfo(Form("Loading TPC eta sigma map from %s/COMMON/PID/data/TPCetaMaps.root", fOADBPath.Data()));
    
    TObjArray* etaSigmaPars = 0x0;
    
    if (fIsMC && !fTuneMConData) {
      TString searchMap = Form("TPCetaSigmaMaps_%s_%s_pass%d", dataType.Data(), period.Data(), recopass);
      etaSigmaPars = dynamic_cast<TObjArray *>(etaSigmaMapsCont.GetDefaultObject(searchMap.Data()));
      if (!etaSigmaPars) {
        // Try default object
        etaSigmaPars = dynamic_cast<TObjArray *>(etaSigmaMapsCont.GetDefaultObject(defaultObj.Data()));
      }
    }
    else {
      etaSigmaPars = dynamic_cast<TObjArray *>(etaSigmaMapsCont.GetObject(fRun, defaultObj.Data()));
    }
    
    if (!etaSigmaPars) {
      AliError(Form("TPC eta sigma parametrisation not found for run %d -> Using old sigma parametrisation!!!", fRun));
    }
    else {
      TH2D* etaSigmaPar1Map = dynamic_cast<TH2D *>(etaSigmaPars->FindObject("sigmaPar1Map"));
      TNamed* sigmaPar0Info = dynamic_cast<TNamed *>(etaSigmaPars->FindObject("sigmaPar0"));
      Double_t sigmaPar0 = 0.0;
      
      if (sigmaPar0Info) {
        TString sigmaPar0String = sigmaPar0Info->GetTitle();
        sigmaPar0 = sigmaPar0String.Atof();
      }
      else {
        // Something is weired because the object for parameter 0 could not be loaded -> New sigma parametrisation can not be used!
        etaSigmaPar1Map = 0x0;
      }
      
      TH2D* etaSigmaPar1MapRefined = RefineHistoViaLinearInterpolation(etaSigmaPar1Map, refineFactorSigmaMapX, refineFactorSigmaMapY);
      
      
      if (etaSigmaPar1MapRefined) {
        if (!fTPCResponse.SetSigmaParams(etaSigmaPar1MapRefined, sigmaPar0)) {
          AliError(Form("Failed to set TPC eta sigma map for run %d -> Using old sigma parametrisation!!!", fRun));
          fTPCResponse.SetSigmaParams(0x0, 0);
        }
        else {
          AliInfo(Form("Loaded TPC sigma correction map (refine factors %.2f/%.2f) from %s/COMMON/PID/data/TPCetaMaps.root: %s", 
                       refineFactorSigmaMapX, refineFactorSigmaMapY, fOADBPath.Data(), fTPCResponse.GetSigmaPar1Map()->GetTitle()));
        }
        
        delete etaSigmaPar1MapRefined;
      }
      else {
        AliError(Form("Failed to set TPC eta sigma map for run %d (map was loaded, but couldn't be refined) -> Using old sigma parametrisation!!!",
                      fRun));
      }
    }
  }
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
  
  //
  //reset old splines
  //
  fTPCResponse.ResetSplines();
  
  if (fLHCperiod.IsNull()) {
    AliError("No period set, not changing parametrisation");
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

  // period
  TString period=fLHCperiod;
  if (fIsMC && !fTuneMConData) period=fMCperiodTPC;

  Int_t recopass = fRecoPass;
  if(fTuneMConData) recopass = fRecoPassUser;
    
  AliInfo(Form("Searching splines for: %s %s PASS%d %s",datatype.Data(),period.Data(),recopass,fBeamType.Data()));
  Bool_t found=kFALSE;
  //
  //set the new PID splines
  //
  if (fArrPidResponseMaster){
    //for MC don't use period information
    //if (fIsMC) period="[A-Z0-9]*";
    //for MC use MC period information
    //pattern for the default entry (valid for all particles)
    TPRegexp reg(Form("TSPLINE3_%s_([A-Z]*)_%s_PASS%d_%s_MEAN(_*)([A-Z1-9]*)",datatype.Data(),period.Data(),recopass,fBeamType.Data()));

    //find particle id and gain scenario
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
              break;
            }
          }
        }
      }
      
      // Retrieve responsefunction for pions - will (if available) be used for muons if there are no dedicated muon splines.
      // For light nuclei, try to set the proton spline, if no dedicated splines are available.
      // In both cases: Use default splines, if no dedicated splines and no pion/proton splines are available.
      TObject* responseFunctionPion = fTPCResponse.GetResponseFunction( (AliPID::EParticleType)AliPID::kPion,                             
                                                                        (AliTPCPIDResponse::ETPCgainScenario)igainScenario);
      TObject* responseFunctionProton = fTPCResponse.GetResponseFunction( (AliPID::EParticleType)AliPID::kProton,                             
                                                                          (AliTPCPIDResponse::ETPCgainScenario)igainScenario);
      
      for (Int_t ispec=0; ispec<(AliTPCPIDResponse::fgkNumberOfParticleSpecies); ++ispec)
      {
        if (!fTPCResponse.GetResponseFunction( (AliPID::EParticleType)ispec,
          (AliTPCPIDResponse::ETPCgainScenario)igainScenario))
        {
          if (ispec == AliPID::kMuon) { // Muons
            if (responseFunctionPion) {
              fTPCResponse.SetResponseFunction( responseFunctionPion,
                                                (AliPID::EParticleType)ispec,
                                                (AliTPCPIDResponse::ETPCgainScenario)igainScenario );
              fTPCResponse.SetUseDatabase(kTRUE);
              AliInfo(Form("Adding graph: %d %d - %s",ispec,igainScenario,responseFunctionPion->GetName()));
              found=kTRUE;  
            }
            else if (grAll) {
              fTPCResponse.SetResponseFunction( grAll,
                                                (AliPID::EParticleType)ispec,
                                                (AliTPCPIDResponse::ETPCgainScenario)igainScenario );
              fTPCResponse.SetUseDatabase(kTRUE);
              AliInfo(Form("Adding graph: %d %d - %s",ispec,igainScenario,grAll->GetName()));
              found=kTRUE;
            }
            //else
            //  AliError(Form("No splines found for muons (also no pion splines and no default splines) for gain scenario %d!", igainScenario));
          }
          else if (ispec >= AliPID::kSPECIES) { // Light nuclei
            if (responseFunctionProton) {
              fTPCResponse.SetResponseFunction( responseFunctionProton,
                                                (AliPID::EParticleType)ispec,
                                                (AliTPCPIDResponse::ETPCgainScenario)igainScenario );
              fTPCResponse.SetUseDatabase(kTRUE);
              AliInfo(Form("Adding graph: %d %d - %s",ispec,igainScenario,responseFunctionProton->GetName()));
              found=kTRUE;  
            }
            else if (grAll) {
              fTPCResponse.SetResponseFunction( grAll,
                                                (AliPID::EParticleType)ispec,
                                                (AliTPCPIDResponse::ETPCgainScenario)igainScenario );
              fTPCResponse.SetUseDatabase(kTRUE);
              AliInfo(Form("Adding graph: %d %d - %s",ispec,igainScenario,grAll->GetName()));
              found=kTRUE;
            }
            //else
            //  AliError(Form("No splines found for species %d (also no proton splines and no default splines) for gain scenario %d!",
            //                ispec, igainScenario));
          }
        }
      }
    }
  }
  else AliInfo("no fArrPidResponseMaster");

  if (!found){
    AliError(Form("No splines found for: %s %s PASS%d %s",datatype.Data(),period.Data(),recopass,fBeamType.Data()));
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
  fResolutionCorrection=(TF1*)fArrPidResponseMaster->FindObject(Form("TF1_%s_ALL_%s_PASS%d_%s_SIGMA",datatype.Data(),period.Data(),recopass,fBeamType.Data()));
  
  if (fResolutionCorrection) AliInfo(Form("Setting multiplicity correction function: %s",fResolutionCorrection->GetName()));

  //read in the voltage map
  TVectorF* gsm = 0x0;
  if (fOADBvoltageMaps) gsm=dynamic_cast<TVectorF*>(fOADBvoltageMaps->GetObject(fRun));
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

  AliInfo("TZERO resolution loaded from ESDrun/AODheader");
  Float_t t0Spread[4];
  for (Int_t i=0;i<4;i++) t0Spread[i]=fCurrentEvent->GetT0spread(i);
  AliInfo(Form("  TZERO spreads from data: (A+C)/2 %f A %f C %f (A'-C')/2: %f",t0Spread[0],t0Spread[1],t0Spread[2],t0Spread[3]));
  Float_t a = t0Spread[1]*t0Spread[1]-t0Spread[0]*t0Spread[0]+t0Spread[3]*t0Spread[3];
  Float_t c = t0Spread[2]*t0Spread[2]-t0Spread[0]*t0Spread[0]+t0Spread[3]*t0Spread[3];
  if ( (t0Spread[0] > 50. && t0Spread[0] < 400.) && (a > 0.) && (c>0.)) {
    fResT0AC=t0Spread[3];
    fResT0A=TMath::Sqrt(a);
    fResT0C=TMath::Sqrt(c);
  } else {
    AliInfo("  TZERO spreads not present or inconsistent, loading default");
    fResT0A=75.;
    fResT0C=65.;
    fResT0AC=55.;
  }
  AliInfo(Form("  TZERO resolution set to: T0A: %f [ps] T0C: %f [ps] T0AC %f [ps]",fResT0A,fResT0C,fResT0AC));

}


//______________________________________________________________________________
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

//______________________________________________________________________________
void AliPIDResponse::FillTrackDetectorPID(const AliVTrack *track, EDetector detector) const
{
  //
  // create detector PID information and setup the transient pointer in the track
  //
  
  // check if detector number is inside accepted range
  if (detector == kNdetectors) return;
  
  // get detector pid
  AliDetectorPID *detPID=const_cast<AliDetectorPID*>(track->GetDetectorPID());
  if (!detPID) {
    detPID=new AliDetectorPID;
    (const_cast<AliVTrack*>(track))->SetDetectorPID(detPID);
  }
  
  //check if values exist
  if (detPID->HasRawProbability(detector) && detPID->HasNumberOfSigmas(detector)) return;
  
  //TODO: which particles to include? See also the loops below...
  Double_t values[AliPID::kSPECIESC]={0};

  //probabilities
  EDetPidStatus status=GetComputePIDProbability(detector,track,AliPID::kSPECIESC,values);
  detPID->SetRawProbability(detector, values, (Int_t)AliPID::kSPECIESC, status);
  
  //nsigmas
  for (Int_t ipart=0; ipart<AliPID::kSPECIESC; ++ipart)
    values[ipart]=GetNumberOfSigmas(detector,track,(AliPID::EParticleType)ipart);
  // the pid status is the same for probabilities and nSigmas, so it is
  // fine to use the one from the probabilities also here
  detPID->SetNumberOfSigmas(detector, values, (Int_t)AliPID::kSPECIESC, status);
  
}

//______________________________________________________________________________
void AliPIDResponse::FillTrackDetectorPID()
{
  //
  // create detector PID information and setup the transient pointer in the track
  //

  if (!fCurrentEvent) return;
  
  for (Int_t itrack=0; itrack<fCurrentEvent->GetNumberOfTracks(); ++itrack){
    AliVTrack *track=dynamic_cast<AliVTrack*>(fCurrentEvent->GetTrack(itrack));
    if (!track) continue;

    for (Int_t idet=0; idet<kNdetectors; ++idet){
      FillTrackDetectorPID(track, (EDetector)idet);
    }
  }
}

//______________________________________________________________________________
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

    Float_t resT0A=fResT0A;
    Float_t resT0C=fResT0C;
    Float_t resT0AC=fResT0AC;
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
	    t0A= vevent->GetT0TOF()[1];
	    t0C= vevent->GetT0TOF()[2];
        //      t0AC= vevent->GetT0TOF()[0];
        t0AC= t0A/resT0A/resT0A + t0C/resT0C/resT0C;
        resT0AC= TMath::Sqrt(1./resT0A/resT0A + 1./resT0C/resT0C);
        t0AC /= resT0AC*resT0AC;
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
	    t0A= vevent->GetT0TOF()[1];
	    t0C= vevent->GetT0TOF()[2];
        //      t0AC= vevent->GetT0TOF()[0];
        t0AC= t0A/resT0A/resT0A + t0C/resT0C/resT0C;
        resT0AC= TMath::Sqrt(1./resT0A/resT0A + 1./resT0C/resT0C);
        t0AC /= resT0AC*resT0AC;
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

//______________________________________________________________________________
// private non cached versions of the PID calculation
//


//______________________________________________________________________________
Float_t AliPIDResponse::GetNumberOfSigmas(EDetector detector, const AliVParticle *vtrack, AliPID::EParticleType type) const
{
  //
  // NumberOfSigmas for 'detCode'
  //

  const AliVTrack *track=static_cast<const AliVTrack*>(vtrack);
  
  switch (detector){
    case kITS:   return GetNumberOfSigmasITS(track, type); break;
    case kTPC:   return GetNumberOfSigmasTPC(track, type); break;
    case kTOF:   return GetNumberOfSigmasTOF(track, type); break;
    case kEMCAL: return GetNumberOfSigmasEMCAL(track, type); break;
    default: return -999.;
  }

  return -999.;
}

//______________________________________________________________________________
Float_t AliPIDResponse::GetNumberOfSigmasITS(const AliVParticle *vtrack, AliPID::EParticleType type) const
{
  //
  // Calculate the number of sigmas in the ITS
  //
  
  AliVTrack *track=(AliVTrack*)vtrack;

  const EDetPidStatus pidStatus=GetITSPIDStatus(track);
  if (pidStatus!=kDetPidOk) return -999.;
    
  UChar_t clumap=track->GetITSClusterMap();
  Int_t nPointsForPid=0;
  for(Int_t i=2; i<6; i++){
    if(clumap&(1<<i)) ++nPointsForPid;
  }
  Float_t mom=track->P();
  
  //check for ITS standalone tracks
  Bool_t isSA=kTRUE;
  if( track->GetStatus() & AliVTrack::kTPCin ) isSA=kFALSE;
  
  const Float_t dEdx=track->GetITSsignal();

  //TODO: in case of the electron, use the SA parametrisation,
  //      this needs to be changed if ITS provides a parametrisation
  //      for electrons also for ITS+TPC tracks
  return fITSResponse.GetNumberOfSigmas(mom,dEdx,type,nPointsForPid,isSA || (type==AliPID::kElectron));
}

//______________________________________________________________________________
Float_t AliPIDResponse::GetNumberOfSigmasTPC(const AliVParticle *vtrack, AliPID::EParticleType type) const
{
  //
  // Calculate the number of sigmas in the TPC
  //
  
  AliVTrack *track=(AliVTrack*)vtrack;

  const EDetPidStatus pidStatus=GetTPCPIDStatus(track);
  if (pidStatus!=kDetPidOk) return -999.;
  
  Double_t nSigma = -999.;
  
  if (fTuneMConData)
    this->GetTPCsignalTunedOnData(track);
  
  nSigma = fTPCResponse.GetNumberOfSigmas(track, type, AliTPCPIDResponse::kdEdxDefault, fUseTPCEtaCorrection);
  
  return nSigma;
}

//______________________________________________________________________________
Float_t AliPIDResponse::GetNumberOfSigmasTOF(const AliVParticle *vtrack, AliPID::EParticleType type) const
{
  //
  // Calculate the number of sigmas in the TOF
  //
  
  AliVTrack *track=(AliVTrack*)vtrack;

  const EDetPidStatus pidStatus=GetTOFPIDStatus(track);
  if (pidStatus!=kDetPidOk) return -999.;

  
  return GetNumberOfSigmasTOFold(vtrack, type);
}

//______________________________________________________________________________
Float_t AliPIDResponse::GetNumberOfSigmasEMCAL(const AliVParticle *vtrack, AliPID::EParticleType type) const
{
  //
  // Calculate the number of sigmas in the EMCAL
  //
  
  AliVTrack *track=(AliVTrack*)vtrack;

  const EDetPidStatus pidStatus=GetEMCALPIDStatus(track);
  if (pidStatus!=kDetPidOk) return -999.;

  const Int_t nMatchClus = track->GetEMCALcluster();
  AliVCluster *matchedClus = (AliVCluster*)fCurrentEvent->GetCaloCluster(nMatchClus);
  
  const Double_t mom    = track->P();
  const Double_t pt     = track->Pt();
  const Int_t    charge = track->Charge();
  const Double_t fClsE  = matchedClus->E();
  const Double_t EovP   = fClsE/mom;
  
  return fEMCALResponse.GetNumberOfSigmas(pt,EovP,type,charge);
}


//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::GetComputePIDProbability  (EDetector detCode,  const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  //
  // Compute PID response of 'detCode'
  //

  switch (detCode){
    case kITS: return GetComputeITSProbability(track, nSpecies, p); break;
    case kTPC: return GetComputeTPCProbability(track, nSpecies, p); break;
    case kTRD: return GetComputeTRDProbability(track, nSpecies, p); break;
    case kTOF: return GetComputeTOFProbability(track, nSpecies, p); break;
    case kPHOS: return GetComputePHOSProbability(track, nSpecies, p); break;
    case kEMCAL: return GetComputeEMCALProbability(track, nSpecies, p); break;
    case kHMPID: return GetComputeHMPIDProbability(track, nSpecies, p); break;
    default: return kDetNoSignal;
  }
}

//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::GetComputeITSProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  //
  // Compute PID response for the ITS
  //
  
  // set flat distribution (no decision)
  for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;
  
  const EDetPidStatus pidStatus=GetITSPIDStatus(track);
  if (pidStatus!=kDetPidOk) return pidStatus;
  
  if (track->GetDetectorPID()){
    return track->GetDetectorPID()->GetRawProbability(kITS, p, nSpecies);
  }
  
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

  Bool_t mismatch=kTRUE/*, heavy=kTRUE*/;
  for (Int_t j=0; j<nSpecies; j++) {
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
  }

  if (mismatch){
    for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;
  }

  return kDetPidOk;
}
//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::GetComputeTPCProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  //
  // Compute PID response for the TPC
  //
  
  // set flat distribution (no decision)
  for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;
  
  const EDetPidStatus pidStatus=GetTPCPIDStatus(track);
  if (pidStatus!=kDetPidOk) return pidStatus;
  
  Double_t dedx=track->GetTPCsignal();
  Bool_t mismatch=kTRUE/*, heavy=kTRUE*/;
  
  if(fTuneMConData) dedx = this->GetTPCsignalTunedOnData(track);
  
  Double_t bethe = 0.;
  Double_t sigma = 0.;
  
  for (Int_t j=0; j<nSpecies; j++) {
    AliPID::EParticleType type=AliPID::EParticleType(j);
    
    bethe=fTPCResponse.GetExpectedSignal(track, type, AliTPCPIDResponse::kdEdxDefault, fUseTPCEtaCorrection);
    sigma=fTPCResponse.GetExpectedSigma(track, type, AliTPCPIDResponse::kdEdxDefault, fUseTPCEtaCorrection);
    
    if (TMath::Abs(dedx-bethe) > fRange*sigma) {
      p[j]=TMath::Exp(-0.5*fRange*fRange)/sigma;
    } else {
      p[j]=TMath::Exp(-0.5*(dedx-bethe)*(dedx-bethe)/(sigma*sigma))/sigma;
      mismatch=kFALSE;
    }
  }
  
  if (mismatch){
    for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;
  }
  
  return kDetPidOk;
}
//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::GetComputeTOFProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  //
  // Compute PID probabilities for TOF
  //
  
  // set flat distribution (no decision)
  for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;
  
  const EDetPidStatus pidStatus=GetTOFPIDStatus(track);
  if (pidStatus!=kDetPidOk) return pidStatus;

  const Double_t meanCorrFactor = 0.11/fTOFtail; // Correction factor on the mean because of the tail (should be ~ 0.1 with tail = 1.1)
  
  for (Int_t j=0; j<nSpecies; j++) {
    AliPID::EParticleType type=AliPID::EParticleType(j);
    const Double_t nsigmas=GetNumberOfSigmasTOFold(track,type) + meanCorrFactor;
    
    const Double_t expTime = fTOFResponse.GetExpectedSignal(track,type);
    const Double_t sig     = fTOFResponse.GetExpectedSigma(track->P(),expTime,AliPID::ParticleMassZ(type));
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
  }
  
  return kDetPidOk;
}
//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::GetComputeTRDProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[],AliTRDPIDResponse::ETRDPIDMethod PIDmethod/*=AliTRDPIDResponse::kLQ1D*/) const
{
  //
  // Compute PID probabilities for the TRD
  //
  
  // set flat distribution (no decision)
  for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;
  
  const EDetPidStatus pidStatus=GetTRDPIDStatus(track);
  if (pidStatus!=kDetPidOk) return pidStatus;

  UInt_t TRDslicesForPID[2];
  SetTRDSlices(TRDslicesForPID,PIDmethod);
  
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
AliPIDResponse::EDetPidStatus AliPIDResponse::GetComputeEMCALProbability  (const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  //
  // Compute PID response for the EMCAL
  //
  
  for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;

  const EDetPidStatus pidStatus=GetEMCALPIDStatus(track);
  if (pidStatus!=kDetPidOk) return pidStatus;

  const Int_t nMatchClus = track->GetEMCALcluster();
  AliVCluster *matchedClus = (AliVCluster*)fCurrentEvent->GetCaloCluster(nMatchClus);
  
  const Double_t mom    = track->P();
  const Double_t pt     = track->Pt();
  const Int_t    charge = track->Charge();
  const Double_t fClsE  = matchedClus->E();
  const Double_t EovP   = fClsE/mom;
  
  // compute the probabilities
  fEMCALResponse.ComputeEMCALProbability(nSpecies,pt,EovP,charge,p);
  return kDetPidOk;
}

//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::GetComputePHOSProbability (const AliVTrack */*track*/, Int_t nSpecies, Double_t p[]) const
{
  //
  // Compute PID response for the PHOS
  //
  
  // set flat distribution (no decision)
  for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;
  return kDetNoSignal;
}

//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::GetComputeHMPIDProbability(const AliVTrack *track, Int_t nSpecies, Double_t p[]) const
{
  //
  // Compute PID response for the HMPID
  //
  
  // set flat distribution (no decision)
  for (Int_t j=0; j<nSpecies; j++) p[j]=1./nSpecies;
  
  const EDetPidStatus pidStatus=GetHMPIDPIDStatus(track);
  if (pidStatus!=kDetPidOk) return pidStatus;
  
  track->GetHMPIDpid(p);
  
  return kDetPidOk;
}

//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::GetITSPIDStatus(const AliVTrack *track) const
{
  // compute ITS pid status

  // check status bits
  if ((track->GetStatus()&AliVTrack::kITSin)==0 &&
    (track->GetStatus()&AliVTrack::kITSout)==0) return kDetNoSignal;

  const Float_t dEdx=track->GetITSsignal();
  if (dEdx<=0) return kDetNoSignal;
  
  // requite at least 3 pid clusters
  const UChar_t clumap=track->GetITSClusterMap();
  Int_t nPointsForPid=0;
  for(Int_t i=2; i<6; i++){
    if(clumap&(1<<i)) ++nPointsForPid;
  }
  
  if(nPointsForPid<3) { 
    return kDetNoSignal;
  }
  
  return kDetPidOk;
}

//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse:: GetTPCPIDStatus(const AliVTrack *track) const
{
  // compute TPC pid status
  
  // check quality of the track
  if ( (track->GetStatus()&AliVTrack::kTPCin )==0 && (track->GetStatus()&AliVTrack::kTPCout)==0 ) return kDetNoSignal;

  // check pid values
  const Double_t dedx=track->GetTPCsignal();
  const UShort_t signalN=track->GetTPCsignalN();
  if (signalN<10 || dedx<10) return kDetNoSignal;

  if (!(fArrPidResponseMaster && fArrPidResponseMaster->At(AliPID::kPion))) return kDetNoParams;
  
  return kDetPidOk;
}

//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::GetTRDPIDStatus(const AliVTrack *track) const
{
  // compute TRD pid status

  if((track->GetStatus()&AliVTrack::kTRDout)==0) return kDetNoSignal;
  return kDetPidOk;
}

//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::GetTOFPIDStatus(const AliVTrack *track) const
{
  // compute TOF pid status

  if ((track->GetStatus()&AliVTrack::kTOFout)==0) return kDetNoSignal;
  if ((track->GetStatus()&AliVTrack::kTIME)==0) return kDetNoSignal;

  return kDetPidOk;
}

//______________________________________________________________________________
Float_t AliPIDResponse::GetTOFMismatchProbability(const AliVTrack *track) const
{
  // compute mismatch probability cross-checking at 5 sigmas with TPC
  // currently just implemented as a 5 sigma compatibility cut

  // check pid status
  const EDetPidStatus tofStatus=GetTOFPIDStatus(track);
  if (tofStatus!=kDetPidOk) return 0.;

  //mismatch
  const EDetPidStatus tpcStatus=GetTPCPIDStatus(track);
  if (tpcStatus!=kDetPidOk) return 0.;
  
  const Double_t meanCorrFactor = 0.11/fTOFtail; // Correction factor on the mean because of the tail (should be ~ 0.1 with tail = 1.1)
  Bool_t mismatch = kTRUE/*, heavy = kTRUE*/;
  for (Int_t j=0; j<AliPID::kSPECIESC; j++) {
    AliPID::EParticleType type=AliPID::EParticleType(j);
    const Double_t nsigmas=GetNumberOfSigmasTOFold(track,type) + meanCorrFactor;
    
    if (TMath::Abs(nsigmas)<5.){
      const Double_t nsigmasTPC=GetNumberOfSigmasTPC(track,type);
      if (TMath::Abs(nsigmasTPC)<5.) mismatch=kFALSE;
    }
  }
  
  if (mismatch){
    return 1.;
  }
  
  return 0.;
}



//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse:: GetHMPIDPIDStatus(const AliVTrack *track) const
{
  // compute HMPID pid status
  if((track->GetStatus()&AliVTrack::kHMPIDpid)==0) return kDetNoSignal;
  return kDetPidOk;
}

//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse:: GetPHOSPIDStatus(const AliVTrack */*track*/) const
{
  // compute PHOS pid status
  return kDetNoSignal;  
}

//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse:: GetEMCALPIDStatus(const AliVTrack *track) const
{
  // compute EMCAL pid status


  // Track matching
  const Int_t nMatchClus = track->GetEMCALcluster();
  if (nMatchClus<0) return kDetNoSignal;

  AliVCluster *matchedClus = (AliVCluster*)fCurrentEvent->GetCaloCluster(nMatchClus);

  if (!(matchedClus && matchedClus->IsEMCAL())) return kDetNoSignal;

  const Int_t charge = track->Charge();
  if (TMath::Abs(charge)!=1) return kDetNoSignal;

  if (!(fEMCALPIDParams && fEMCALPIDParams->At(AliPID::kElectron))) return kDetNoParams;
  
  return kDetPidOk;

}

//______________________________________________________________________________
AliPIDResponse::EDetPidStatus AliPIDResponse::GetPIDStatus(EDetector detector, const AliVTrack *track) const
{
  //
  // check pid status for a track
  //

  switch (detector){
    case kITS:   return GetITSPIDStatus(track);   break;
    case kTPC:   return GetTPCPIDStatus(track);   break;
    case kTRD:   return GetTRDPIDStatus(track);   break;
    case kTOF:   return GetTOFPIDStatus(track);   break;
    case kPHOS:  return GetPHOSPIDStatus(track);  break;
    case kEMCAL: return GetEMCALPIDStatus(track); break;
    case kHMPID: return GetHMPIDPIDStatus(track); break;
    default: return kDetNoSignal;
  }
  return kDetNoSignal;
  
}
