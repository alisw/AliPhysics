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

//-----------------------------------------------------------------
//           Implementation of the TPC PID class
// Very naive one... Should be made better by the detector experts...
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
// With many additions and modifications suggested by
//      Alexander Kalweit, GSI, alexander.philipp.kalweit@cern.ch
//      Dariusz Miskowiec, GSI, D.Miskowiec@gsi.de
// ...and some modifications by
//      Mikolaj Krzewicki, GSI, mikolaj.krzewicki@cern.ch
// ...and some modifications plus eta correction functions by
//      Benjamin Hess, University of Tuebingen, bhess@cern.ch
//-----------------------------------------------------------------

#include <TGraph.h>
#include <TObjArray.h>
#include <TSpline.h>
#include <TBits.h>
#include <TMath.h>
#include <TH2D.h>

#include <AliLog.h>
#include "AliExternalTrackParam.h"
#include "AliVTrack.h"
#include "AliTPCPIDResponse.h"
#include "AliTPCdEdxInfo.h"

ClassImp(AliTPCPIDResponse);

const char* AliTPCPIDResponse::fgkGainScenarioName[fgkNumberOfGainScenarios+1]=
{
  "", //default - no name
  "1", //all high
  "2", //OROC only
  "unknown"
};

//_________________________________________________________________________
AliTPCPIDResponse::AliTPCPIDResponse():
  TNamed(),
  fMIP(50.),
  fRes0(),
  fResN2(),
  fKp1(0.0283086),
  fKp2(2.63394e+01),
  fKp3(5.04114e-11),
  fKp4(2.12543),
  fKp5(4.88663),
  fUseDatabase(kFALSE),
  fResponseFunctions(fgkNumberOfParticleSpecies*fgkNumberOfGainScenarios),
  fVoltageMap(72),
  fLowGainIROCthreshold(-40),
  fBadIROCthreshhold(-70),
  fLowGainOROCthreshold(-40),
  fBadOROCthreshhold(-40),
  fMaxBadLengthFraction(0.5),
  fMagField(0.),
  fhEtaCorr(0x0),
  fhEtaSigmaPar1(0x0),
  fSigmaPar0(0.0)
{
  //
  //  The default constructor
  //
  for (Int_t i=0; i<fgkNumberOfGainScenarios; i++) {fRes0[i]=0.07;fResN2[i]=0.0;}
}
/*TODO remove?
//_________________________________________________________________________
AliTPCPIDResponse::AliTPCPIDResponse(const Double_t *param):
  TNamed(),
  fMIP(param[0]),
  fRes0(),
  fResN2(),
  fKp1(0.0283086),
  fKp2(2.63394e+01),
  fKp3(5.04114e-11),
  fKp4(2.12543),
  fKp5(4.88663),
  fUseDatabase(kFALSE),
  fResponseFunctions(fgkNumberOfParticleSpecies*fgkNumberOfGainScenarios),
  fVoltageMap(72),
  fLowGainIROCthreshold(-40),
  fBadIROCthreshhold(-70),
  fLowGainOROCthreshold(-40),
  fBadOROCthreshhold(-40),
  fMaxBadLengthFraction(0.5),
  fMagField(0.),
  fhEtaCorr(0x0),
  fhEtaSigmaPar1(0x0),
  fSigmaPar0(0.0)
{
  //
  //  The main constructor
  //
  for (Int_t i=0; i<fgkNumberOfGainScenarios; i++) {fRes0[i]=param[1];fResN2[i]=param[2];}
}
*/

//_________________________________________________________________________
AliTPCPIDResponse::~AliTPCPIDResponse()
{
  //
  // Destructor
  //
  
  delete fhEtaCorr;
  fhEtaCorr = 0x0;
  
  delete fhEtaSigmaPar1;
  fhEtaSigmaPar1 = 0x0;
}


//_________________________________________________________________________
AliTPCPIDResponse::AliTPCPIDResponse(const AliTPCPIDResponse& that):
  TNamed(that),
  fMIP(that.fMIP),
  fRes0(),
  fResN2(),
  fKp1(that.fKp1),
  fKp2(that.fKp2),
  fKp3(that.fKp3),
  fKp4(that.fKp4),
  fKp5(that.fKp5),
  fUseDatabase(that.fUseDatabase),
  fResponseFunctions(that.fResponseFunctions),
  fVoltageMap(that.fVoltageMap),
  fLowGainIROCthreshold(that.fLowGainIROCthreshold),
  fBadIROCthreshhold(that.fBadIROCthreshhold),
  fLowGainOROCthreshold(that.fLowGainOROCthreshold),
  fBadOROCthreshhold(that.fBadOROCthreshhold),
  fMaxBadLengthFraction(that.fMaxBadLengthFraction),
  fMagField(that.fMagField),
  fhEtaCorr(0x0),
  fhEtaSigmaPar1(0x0),
  fSigmaPar0(that.fSigmaPar0)
{
  //copy ctor
  for (Int_t i=0; i<fgkNumberOfGainScenarios; i++) {fRes0[i]=that.fRes0[i];fResN2[i]=that.fResN2[i];}
 
  // Copy eta maps
  if (that.fhEtaCorr) {
    fhEtaCorr = new TH2D(*(that.fhEtaCorr));
    fhEtaCorr->SetDirectory(0);
  }
  
  if (that.fhEtaSigmaPar1) {
    fhEtaSigmaPar1 = new TH2D(*(that.fhEtaSigmaPar1));
    fhEtaSigmaPar1->SetDirectory(0);
  }
}

//_________________________________________________________________________
AliTPCPIDResponse& AliTPCPIDResponse::operator=(const AliTPCPIDResponse& that)
{
  //assignment
  if (&that==this) return *this;
  TNamed::operator=(that);
  fMIP=that.fMIP;
  fKp1=that.fKp1;
  fKp2=that.fKp2;
  fKp3=that.fKp3;
  fKp4=that.fKp4;
  fKp5=that.fKp5;
  fUseDatabase=that.fUseDatabase;
  fResponseFunctions=that.fResponseFunctions;
  fVoltageMap=that.fVoltageMap;
  fLowGainIROCthreshold=that.fLowGainIROCthreshold;
  fBadIROCthreshhold=that.fBadIROCthreshhold;
  fLowGainOROCthreshold=that.fLowGainOROCthreshold;
  fBadOROCthreshhold=that.fBadOROCthreshhold;
  fMaxBadLengthFraction=that.fMaxBadLengthFraction;
  fMagField=that.fMagField;
  for (Int_t i=0; i<fgkNumberOfGainScenarios; i++) {fRes0[i]=that.fRes0[i];fResN2[i]=that.fResN2[i];}

  delete fhEtaCorr;
  fhEtaCorr=0x0;
  if (that.fhEtaCorr) {
    fhEtaCorr = new TH2D(*(that.fhEtaCorr));
    fhEtaCorr->SetDirectory(0);
  }
  
  delete fhEtaSigmaPar1;
  fhEtaSigmaPar1=0x0;
  if (that.fhEtaSigmaPar1) {
    fhEtaSigmaPar1 = new TH2D(*(that.fhEtaSigmaPar1));
    fhEtaSigmaPar1->SetDirectory(0);
  }
  
  fSigmaPar0 = that.fSigmaPar0;

  return *this;
}

//_________________________________________________________________________
Double_t AliTPCPIDResponse::Bethe(Double_t betaGamma) const {
  //
  // This is the Bethe-Bloch function normalised to 1 at the minimum
  // WARNING
  // Simulated and reconstructed Bethe-Bloch differs
  //           Simulated  curve is the dNprim/dx
  //           Reconstructed is proportianal dNtot/dx
  // Temporary fix for production -  Simple linear correction function
  // Future    2 Bethe Bloch formulas needed
  //           1. for simulation
  //           2. for reconstructed PID
  //
  
//   const Float_t kmeanCorrection =0.1;
  Double_t bb=
    AliExternalTrackParam::BetheBlochAleph(betaGamma,fKp1,fKp2,fKp3,fKp4,fKp5);
  return bb*fMIP;
}

//_________________________________________________________________________
void AliTPCPIDResponse::SetBetheBlochParameters(Double_t kp1,
                             Double_t kp2,
                             Double_t kp3,
                             Double_t kp4,
                             Double_t kp5) {
  //
  // Set the parameters of the ALEPH Bethe-Bloch formula
  //
  fKp1=kp1;
  fKp2=kp2;
  fKp3=kp3;
  fKp4=kp4;
  fKp5=kp5;
}

//_________________________________________________________________________
void AliTPCPIDResponse::SetSigma(Float_t res0, Float_t resN2) {
  //
  // Set the relative resolution  sigma_rel = res0 * sqrt(1+resN2/npoint)
  //
  for (Int_t i=0; i<fgkNumberOfGainScenarios; i++) {fRes0[i]=res0;fResN2[i]=resN2;}
}

//_________________________________________________________________________
Double_t AliTPCPIDResponse::GetExpectedSignal(const Float_t mom,
					      AliPID::EParticleType n) const {
  //
  // Deprecated function (for backward compatibility). Please use 
  // GetExpectedSignal(const AliVTrack* track, AliPID::EParticleType species, ETPCdEdxSource dedxSource,
  //                   Bool_t correctEta = kTRUE);
  // instead!
  //
  //
  // Calculates the expected PID signal as the function of 
  // the information stored in the track, for the specified particle type 
  //  
  // At the moment, these signals are just the results of calling the 
  // Bethe-Bloch formula. 
  // This can be improved. By taking into account the number of
  // assigned clusters and/or the track dip angle, for example.  
  //
  
  //charge factor. BB goes with z^2, however in reality it is slightly larger (calibration, threshold effects, ...)
  // !!! Splines for light nuclei need to be normalised to this factor !!!
  const Double_t chargeFactor = TMath::Power(AliPID::ParticleCharge(n),2.3);
  
  Double_t mass=AliPID::ParticleMassZ(n);
  if (!fUseDatabase) return Bethe(mom/mass) * chargeFactor;
  //
  const TSpline3 * responseFunction = (TSpline3 *) fResponseFunctions.UncheckedAt(n);

  if (!responseFunction) return Bethe(mom/mass) * chargeFactor;
  
  return fMIP*responseFunction->Eval(mom/mass)*chargeFactor;

}

//_________________________________________________________________________
Double_t AliTPCPIDResponse::GetExpectedSigma(const Float_t mom, 
                                             const Int_t nPoints,
                                             AliPID::EParticleType n) const {
  //
  // Deprecated function (for backward compatibility). Please use 
  // GetExpectedSigma(onst AliVTrack* track, AliPID::EParticleType species, 
  // ETPCdEdxSource dedxSource, Bool_t correctEta) instead!
  //
  //
  // Calculates the expected sigma of the PID signal as the function of 
  // the information stored in the track, for the specified particle type 
  //  
  
  if (nPoints != 0) 
    return GetExpectedSignal(mom,n)*fRes0[0]*sqrt(1. + fResN2[0]/nPoints);
  else
    return GetExpectedSignal(mom,n)*fRes0[0];
}

////////////////////////////////////////////////////NEW//////////////////////////////

//_________________________________________________________________________
void AliTPCPIDResponse::SetSigma(Float_t res0, Float_t resN2, ETPCgainScenario gainScenario) {
  //
  // Set the relative resolution  sigma_rel = res0 * sqrt(1+resN2/npoint)
  //
  if ((Int_t)gainScenario>(Int_t)fgkNumberOfGainScenarios) return; //TODO: better handling!
  fRes0[gainScenario]=res0;
  fResN2[gainScenario]=resN2;
}


//_________________________________________________________________________
Double_t AliTPCPIDResponse::GetExpectedSignal(const AliVTrack* track,
                                              AliPID::EParticleType species,
                                              Double_t /*dEdx*/,
                                              const TSpline3* responseFunction,
                                              Bool_t correctEta) const 
{
  // Calculates the expected PID signal as the function of 
  // the information stored in the track and the given parameters,
  // for the specified particle type 
  //  
  // At the moment, these signals are just the results of calling the 
  // Bethe-Bloch formula plus, if desired, taking into account the eta dependence. 
  // This can be improved. By taking into account the number of
  // assigned clusters and/or the track dip angle, for example.  
  //
  
  Double_t mom=track->GetTPCmomentum();
  Double_t mass=AliPID::ParticleMassZ(species);
  
  //charge factor. BB goes with z^2, however in reality it is slightly larger (calibration, threshold effects, ...)
  // !!! Splines for light nuclei need to be normalised to this factor !!!
  const Double_t chargeFactor = TMath::Power(AliPID::ParticleCharge(species),2.3);
  
  if (!responseFunction)
    return Bethe(mom/mass) * chargeFactor;
  
  Double_t dEdxSplines = fMIP*responseFunction->Eval(mom/mass) * chargeFactor;
  
  if (!correctEta)
    return dEdxSplines;
  
  //TODO Alternatively take current track dEdx
  //return dEdxSplines * GetEtaCorrection(track, dEdx);
  return dEdxSplines * GetEtaCorrection(track, dEdxSplines);
}


//_________________________________________________________________________
Double_t AliTPCPIDResponse::GetExpectedSignal(const AliVTrack* track,
                                              AliPID::EParticleType species,
                                              ETPCdEdxSource dedxSource,
                                              Bool_t correctEta) const
{
  // Calculates the expected PID signal as the function of 
  // the information stored in the track, for the specified particle type 
  //  
  // At the moment, these signals are just the results of calling the 
  // Bethe-Bloch formula plus, if desired, taking into account the eta dependence. 
  // This can be improved. By taking into account the number of
  // assigned clusters and/or the track dip angle, for example.  
  //
  
  if (!fUseDatabase) {
    //charge factor. BB goes with z^2, however in reality it is slightly larger (calibration, threshold effects, ...)
    // !!! Splines for light nuclei need to be normalised to this factor !!!
    const Double_t chargeFactor = TMath::Power(AliPID::ParticleCharge(species),2.3);
  
    return Bethe(track->GetTPCmomentum() / AliPID::ParticleMassZ(species)) * chargeFactor;
  }
  
  Double_t dEdx = -1;
  Int_t nPoints = -1;
  ETPCgainScenario gainScenario = kGainScenarioInvalid;
  TSpline3* responseFunction = 0x0;
    
  if (!ResponseFunctiondEdxN(track, species, dedxSource, dEdx, nPoints, gainScenario, &responseFunction)) {
    // Something is wrong with the track -> Return obviously invalid value
    return -999;
  }
  
  // Charge factor already taken into account inside the following function call
  return GetExpectedSignal(track, species, dEdx, responseFunction, correctEta);
}
  
//_________________________________________________________________________
TSpline3* AliTPCPIDResponse::GetResponseFunction( AliPID::EParticleType type,
                                                  AliTPCPIDResponse::ETPCgainScenario gainScenario ) const
{
  //get response function
  return dynamic_cast<TSpline3*>(fResponseFunctions.At(ResponseFunctionIndex(type,gainScenario)));
}

//_________________________________________________________________________
TSpline3* AliTPCPIDResponse::GetResponseFunction( const AliVTrack* track,
                               AliPID::EParticleType species,
                               ETPCdEdxSource dedxSource) const 
{
  //the splines are stored in an array, different scenarios

  Double_t dEdx = -1;
  Int_t nPoints = -1;
  ETPCgainScenario gainScenario = kGainScenarioInvalid;
  TSpline3* responseFunction = 0x0;
  
  if (ResponseFunctiondEdxN(track, species, dedxSource, dEdx, nPoints, gainScenario, &responseFunction))
    return responseFunction;
  
  return NULL;
}

//_________________________________________________________________________
void AliTPCPIDResponse::ResetSplines()
{
  //reset the array with splines
  for (Int_t i=0;i<fResponseFunctions.GetEntriesFast();i++)
  {
    fResponseFunctions.AddAt(NULL,i);
  }
}
//_________________________________________________________________________
Int_t AliTPCPIDResponse::ResponseFunctionIndex( AliPID::EParticleType species,
                                                ETPCgainScenario gainScenario ) const
{
  //get the index in fResponseFunctions given type and scenario
  return Int_t(species)+Int_t(gainScenario)*fgkNumberOfParticleSpecies;
}

//_________________________________________________________________________
void AliTPCPIDResponse::SetResponseFunction( TObject* o,
                                             AliPID::EParticleType species,
                                             ETPCgainScenario gainScenario )
{
  fResponseFunctions.AddAtAndExpand(o,ResponseFunctionIndex(species,gainScenario));
}


//_________________________________________________________________________
Double_t AliTPCPIDResponse::GetExpectedSigma(const AliVTrack* track, 
                                             AliPID::EParticleType species,
                                             ETPCgainScenario gainScenario,
                                             Double_t dEdx,
                                             Int_t nPoints,
                                             const TSpline3* responseFunction,
                                             Bool_t correctEta) const 
{
  // Calculates the expected sigma of the PID signal as the function of 
  // the information stored in the track and the given parameters,
  // for the specified particle type 
  //
  
  if (!responseFunction)
    return 999;
  
  // If no sigma map is available or if no eta correction is requested (sigma maps only for corrected eta!), use the old parametrisation
  if (!fhEtaSigmaPar1 || !correctEta) {  
    if (nPoints != 0) 
      return GetExpectedSignal(track, species, dEdx, responseFunction, kFALSE) *
               fRes0[gainScenario] * sqrt(1. + fResN2[gainScenario]/nPoints);
    else
      return GetExpectedSignal(track, species, dEdx, responseFunction, kFALSE)*fRes0[gainScenario];
  }
    
  if (nPoints > 0) {
    Double_t sigmaPar1 = GetSigmaPar1(track, species, dEdx, responseFunction);
    
    return GetExpectedSignal(track, species, dEdx, responseFunction, kTRUE)*sqrt(fSigmaPar0 * fSigmaPar0 + sigmaPar1 * sigmaPar1 / nPoints);
  }
  else { 
    // One should never have/take tracks with 0 dEdx clusters!
    return 999;
  }
}


//_________________________________________________________________________
Double_t AliTPCPIDResponse::GetExpectedSigma(const AliVTrack* track, 
                                             AliPID::EParticleType species,
                                             ETPCdEdxSource dedxSource,
                                             Bool_t correctEta) const 
{
  // Calculates the expected sigma of the PID signal as the function of 
  // the information stored in the track, for the specified particle type 
  // and dedx scenario
  //
  
  Double_t dEdx = -1;
  Int_t nPoints = -1;
  ETPCgainScenario gainScenario = kGainScenarioInvalid;
  TSpline3* responseFunction = 0x0;
  
  if (!ResponseFunctiondEdxN(track, species, dedxSource, dEdx, nPoints, gainScenario, &responseFunction))
    return 999; //TODO: better handling!
  
  return GetExpectedSigma(track, species, gainScenario, dEdx, nPoints, responseFunction, correctEta);
}


//_________________________________________________________________________
Float_t AliTPCPIDResponse::GetNumberOfSigmas(const AliVTrack* track, 
                             AliPID::EParticleType species,
                             ETPCdEdxSource dedxSource,
                             Bool_t correctEta) const
{
  //Calculates the number of sigmas of the PID signal from the expected value
  //for a given particle species in the presence of multiple gain scenarios
  //inside the TPC
  
  Double_t dEdx = -1;
  Int_t nPoints = -1;
  ETPCgainScenario gainScenario = kGainScenarioInvalid;
  TSpline3* responseFunction = 0x0;
  
  if (!ResponseFunctiondEdxN(track, species, dedxSource, dEdx, nPoints, gainScenario, &responseFunction))
    return -999; //TODO: Better handling!
    
  Double_t bethe = GetExpectedSignal(track, species, dEdx, responseFunction, correctEta);
  Double_t sigma = GetExpectedSigma(track, species, gainScenario, dEdx, nPoints, responseFunction, correctEta);
  // 999 will be returned by GetExpectedSigma e.g. in case of 0 dEdx clusters
  if (sigma >= 998) 
    return -999;
  else
    return (dEdx-bethe)/sigma;
}

//_________________________________________________________________________
Float_t AliTPCPIDResponse::GetSignalDelta(const AliVTrack* track,
                                          AliPID::EParticleType species,
                                          ETPCdEdxSource dedxSource,
                                          Bool_t correctEta) const
{
  //Calculates the number of sigmas of the PID signal from the expected value
  //for a given particle species in the presence of multiple gain scenarios
  //inside the TPC

  Double_t dEdx = -1;
  Int_t nPoints = -1;
  ETPCgainScenario gainScenario = kGainScenarioInvalid;
  TSpline3* responseFunction = 0x0;

  if (!ResponseFunctiondEdxN(track, species, dedxSource, dEdx, nPoints, gainScenario, &responseFunction))
    return -9999.; //TODO: Better handling!

  Double_t bethe = GetExpectedSignal(track, species, dEdx, responseFunction, correctEta);
  // 999 will be returned by GetExpectedSigma e.g. in case of 0 dEdx clusters
  return dEdx-bethe;
}

//_________________________________________________________________________
Bool_t AliTPCPIDResponse::ResponseFunctiondEdxN( const AliVTrack* track, 
                                                 AliPID::EParticleType species,
                                                 ETPCdEdxSource dedxSource,
                                                 Double_t& dEdx,
                                                 Int_t& nPoints,
                                                 ETPCgainScenario& gainScenario,
                                                 TSpline3** responseFunction) const 
{
  // Calculates the right parameters for PID
  //   dEdx parametrization for the proper gain scenario, dEdx 
  //   and NPoints used for dEdx
  // based on the track geometry (which chambers it crosses) for the specified particle type 
  // and preferred source of dedx.
  // returns true on success
  
  
  if (dedxSource == kdEdxDefault) {
    // Fast handling for default case. In addition: Keep it simple (don't call additional functions) to
    // avoid possible bugs
    
    // GetTPCsignalTunedOnData will be non-positive, if it has not been set (i.e. in case of MC NOT tuned to data).
    // If this is the case, just take the normal signal
    dEdx = track->GetTPCsignalTunedOnData();
    if (dEdx <= 0) {
      dEdx = track->GetTPCsignal();
    }
    
    nPoints = track->GetTPCsignalN();
    gainScenario = kDefault;
    
    TObject* obj = fResponseFunctions.UncheckedAt(ResponseFunctionIndex(species,gainScenario));
    *responseFunction = dynamic_cast<TSpline3*>(obj); //TODO:maybe static cast?
  
    return kTRUE;
  }
  
  //TODO Proper handle of tuneMConData for other dEdx sources
  
  Double32_t signal[4]; //0: IROC, 1: OROC medium, 2:OROC long, 3: OROC all (def. truncation used)
  Char_t ncl[3];        //same
  Char_t nrows[3];      //same
  const AliTPCdEdxInfo* dEdxInfo = track->GetTPCdEdxInfo();
  
  if (!dEdxInfo && dedxSource!=kdEdxDefault)  //in one case its ok if we dont have the info
  {
    AliError("AliTPCdEdxInfo not available");
    return kFALSE;
  }

  if (dEdxInfo) dEdxInfo->GetTPCSignalRegionInfo(signal,ncl,nrows);

  //check if we cross a bad OROC in which case we reject
  EChamberStatus trackOROCStatus = TrackStatus(track,2);
  if (trackOROCStatus==kChamberOff || trackOROCStatus==kChamberLowGain)
  {
    return kFALSE;
  }

  switch (dedxSource)
  {
    case kdEdxOROC:
      {
        if (trackOROCStatus==kChamberInvalid) return kFALSE; //never reached OROC
        dEdx = signal[3];
        nPoints = ncl[2]+ncl[1];
        gainScenario = kOROChigh;
        break;
      }
    case kdEdxHybrid:
      {
        //if we cross a bad IROC we use OROC dedx, if we dont we use combined
        EChamberStatus status = TrackStatus(track,1);
        if (status!=kChamberHighGain)
        {
          dEdx = signal[3];
          nPoints = ncl[2]+ncl[1];
          gainScenario = kOROChigh;
        }
        else
        {
          dEdx = track->GetTPCsignal();
          nPoints = track->GetTPCsignalN();
          gainScenario = kALLhigh;
        }
        break;
      }
    default:
      {
         dEdx = 0.;
         nPoints = 0;
         gainScenario = kGainScenarioInvalid;
         return kFALSE;
      }
  }
  TObject* obj = fResponseFunctions.UncheckedAt(ResponseFunctionIndex(species,gainScenario));
  *responseFunction = dynamic_cast<TSpline3*>(obj); //TODO:maybe static cast?
  
  return kTRUE;
}


//_________________________________________________________________________
Double_t AliTPCPIDResponse::GetEtaCorrection(const AliVTrack *track, Double_t dEdxSplines) const
{
  //
  // Get eta correction for the given parameters.
  //
  
  if (!fhEtaCorr)
    return 1.;
  
  Double_t tpcSignal = dEdxSplines;
  
  if (tpcSignal < 1.)
    return 1.;
  
  // For ESD tracks, the local tanTheta could be used (esdTrack->GetInnerParam()->GetTgl()).
  // However, this value is not available for AODs and, thus, not for AliVTrack.
  // Fortunately, the following formula allows to approximate the local tanTheta with the 
  // global theta angle -> This is for by far most of the tracks the same, but gives at
  // maybe the percent level differences within +- 0.2 in tanTheta -> Which is still ok.
  Double_t tanTheta = TMath::Tan(-track->Theta() + TMath::Pi() / 2.0);
  Int_t binX = fhEtaCorr->GetXaxis()->FindBin(tanTheta);
  Int_t binY = fhEtaCorr->GetYaxis()->FindBin(1. / tpcSignal);
  
  if (binX == 0) 
    binX = 1;
  if (binX > fhEtaCorr->GetXaxis()->GetNbins())
    binX = fhEtaCorr->GetXaxis()->GetNbins();
  
  if (binY == 0)
    binY = 1;
  if (binY > fhEtaCorr->GetYaxis()->GetNbins())
    binY = fhEtaCorr->GetYaxis()->GetNbins();
  
  return fhEtaCorr->GetBinContent(binX, binY);
}


//_________________________________________________________________________
Double_t AliTPCPIDResponse::GetEtaCorrection(const AliVTrack *track, AliPID::EParticleType species, ETPCdEdxSource dedxSource) const
{
  //
  // Get eta correction for the given track.
  //
  
  if (!fhEtaCorr)
    return 1.;
  
  Double_t dEdx = -1;
  Int_t nPoints = -1;
  ETPCgainScenario gainScenario = kGainScenarioInvalid;
  TSpline3* responseFunction = 0x0;
  
  if (!ResponseFunctiondEdxN(track, species, dedxSource, dEdx, nPoints, gainScenario, &responseFunction))
    return 1.; 
  
  Double_t dEdxSplines = GetExpectedSignal(track, species, dEdx, responseFunction, kFALSE);
  
  //TODO Alternatively take current track dEdx
  //return GetEtaCorrection(track, dEdx);
  
  return GetEtaCorrection(track, dEdxSplines);
}


//_________________________________________________________________________
Double_t AliTPCPIDResponse::GetEtaCorrectedTrackdEdx(const AliVTrack *track, AliPID::EParticleType species, ETPCdEdxSource dedxSource) const
{
  //
  // Get eta corrected dEdx for the given track. For the correction, the expected dEdx of
  // the specified species will be used. If the species is set to AliPID::kUnknown, the
  // dEdx of the track is used instead.
  // WARNING: In the latter case, the eta correction might not be as good as if the
  // expected dEdx is used, which is the way the correction factor is designed
  // for.
  // In any case, one has to decide carefully to which expected signal one wants to
  // compare the corrected value - to the corrected or uncorrected.
  // Anyhow, a safer way of looking e.g. at the n-sigma is to call
  // the corresponding function GetNumberOfSigmas!
  //
  
  Double_t dEdx = -1;
  Int_t nPoints = -1;
  ETPCgainScenario gainScenario = kGainScenarioInvalid;
  TSpline3* responseFunction = 0x0;
  
  // Note: In case of species == AliPID::kUnknown, the responseFunction might not be set. However, in this case
  // it is not used anyway, so this causes no trouble.
  if (!ResponseFunctiondEdxN(track, species, dedxSource, dEdx, nPoints, gainScenario, &responseFunction))
    return -1.;
  
  Double_t etaCorr = 0.;
  
  if (species < AliPID::kUnknown) {
    Double_t dEdxSplines = GetExpectedSignal(track, species, dEdx, responseFunction, kFALSE);
    etaCorr = GetEtaCorrection(track, dEdxSplines);
  }
  else {
    etaCorr = GetEtaCorrection(track, dEdx);
  }
    
  if (etaCorr <= 0)
    return -1.;
  
  return dEdx / etaCorr; 
}



//_________________________________________________________________________
Double_t AliTPCPIDResponse::GetSigmaPar1(const AliVTrack *track, AliPID::EParticleType species, Double_t dEdx, const TSpline3* responseFunction) const
{
  //
  // Get parameter 1 of sigma parametrisation of TPC dEdx from the histogram for the given track.
  //
  
  if (!fhEtaSigmaPar1)
    return 999;
  
  // The sigma maps are created with data sets that are already eta corrected and for which the 
  // splines have been re-created. Therefore, the value for the lookup needs to be the value of
  // the splines without any additional eta correction.
  // NOTE: This is due to the method the maps are created. The track dEdx (not the expected one!)
  // is corrected to uniquely related a momemtum bin with an expected dEdx, where the expected dEdx
  // equals the track dEdx for all eta (thanks to the correction and the re-fit of the splines).
  // Consequently, looking up the uncorrected expected dEdx at a given tanTheta yields the correct
  // sigma parameter!
  // Also: It has to be the spline dEdx, since one wants to get the sigma for the assumption(!)
  // of such a particle, which by assumption then has this dEdx value
    
  Double_t dEdxExpected = GetExpectedSignal(track, species, dEdx, responseFunction, kFALSE);
  
  if (dEdxExpected < 1.)
    return 999;
  
  // For ESD tracks, the local tanTheta could be used (esdTrack->GetInnerParam()->GetTgl()).
  // However, this value is not available for AODs and, thus, not or AliVTrack.
  // Fortunately, the following formula allows to approximate the local tanTheta with the 
  // global theta angle -> This is for by far most of the tracks the same, but gives at
  // maybe the percent level differences within +- 0.2 in tanTheta -> Which is still ok.
  Double_t tanTheta = TMath::Tan(-track->Theta() + TMath::Pi() / 2.0);
  Int_t binX = fhEtaSigmaPar1->GetXaxis()->FindBin(tanTheta);
  Int_t binY = fhEtaSigmaPar1->GetYaxis()->FindBin(1. / dEdxExpected);
    
  if (binX == 0) 
    binX = 1;
  if (binX > fhEtaSigmaPar1->GetXaxis()->GetNbins())
    binX = fhEtaSigmaPar1->GetXaxis()->GetNbins();
    
  if (binY == 0)
    binY = 1;
  if (binY > fhEtaSigmaPar1->GetYaxis()->GetNbins())
    binY = fhEtaSigmaPar1->GetYaxis()->GetNbins();
    
  return fhEtaSigmaPar1->GetBinContent(binX, binY);
}


//_________________________________________________________________________
Double_t AliTPCPIDResponse::GetSigmaPar1(const AliVTrack *track, AliPID::EParticleType species, ETPCdEdxSource dedxSource) const
{
  //
  // Get eta correction for the given track.
  //
  
  if (!fhEtaSigmaPar1)
    return 999;
  
  Double_t dEdx = -1;
  Int_t nPoints = -1;
  ETPCgainScenario gainScenario = kGainScenarioInvalid;
  TSpline3* responseFunction = 0x0;
  
  if (!ResponseFunctiondEdxN(track, species, dedxSource, dEdx, nPoints, gainScenario, &responseFunction))
    return 999; 
  
  return GetSigmaPar1(track, species, dEdx, responseFunction);
}


//_________________________________________________________________________
Bool_t AliTPCPIDResponse::SetEtaCorrMap(TH2D* hMap)
{
  //
  // Load map for TPC eta correction (a copy is stored and will be deleted automatically).
  // If hMap is 0x0,the eta correction will be disabled and kFALSE is returned. 
  // If the map can be set, kTRUE is returned.
  //
  
  delete fhEtaCorr;
  
  if (!hMap) {
    fhEtaCorr = 0x0;
    
    return kFALSE;
  }
  
  fhEtaCorr = (TH2D*)(hMap->Clone());
  fhEtaCorr->SetDirectory(0);
      
  return kTRUE;
}

//_________________________________________________________________________
Bool_t AliTPCPIDResponse::SetSigmaParams(TH2D* hSigmaPar1Map, Double_t sigmaPar0)
{
  //
  // Load map for TPC sigma map (a copy is stored and will be deleted automatically):
  // Parameter 1 is stored as a 2D map (1/dEdx vs. tanTheta_local) and parameter 0 is
  // a constant. If hSigmaPar1Map is 0x0, the old sigma parametrisation will be used
  // (and sigmaPar0 is ignored!) and kFALSE is returned. 
  // If the map can be set, sigmaPar0 is also set and kTRUE will be returned.
  //
  
  delete fhEtaSigmaPar1;
  
  if (!hSigmaPar1Map) {
    fhEtaSigmaPar1 = 0x0;
    fSigmaPar0 = 0.0;
    
    return kFALSE;
  }
  
  fhEtaSigmaPar1 = (TH2D*)(hSigmaPar1Map->Clone());
  fhEtaSigmaPar1->SetDirectory(0);
  fSigmaPar0 = sigmaPar0;
  
  return kTRUE;
}


//_________________________________________________________________________
Bool_t AliTPCPIDResponse::sectorNumbersInOut(Double_t* trackPositionInner,
                                             Double_t* trackPositionOuter,
                                             Float_t& inphi,
                                             Float_t& outphi,
                                             Int_t& in,
                                             Int_t& out ) const
{
  //calculate the sector numbers (equivalent to IROC chamber numbers) a track crosses
  //for OROC chamber numbers add 36
  //returned angles are between (0,2pi)

  inphi = TMath::ATan2(trackPositionInner[1],trackPositionInner[0]);
  outphi = TMath::ATan2(trackPositionOuter[1], trackPositionOuter[0]);

  if (inphi<0) {inphi+=TMath::TwoPi();} //because ATan2 gives angles between -pi,pi
  if (outphi<0) {outphi+=TMath::TwoPi();} //because ATan2 gives angles between -pi,pi

  in = sectorNumber(inphi);
  out = sectorNumber(outphi);

  //for the C side (positive z) offset by 18
  if (trackPositionInner[2]>0.0)
  {
    in+=18;
    out+=18;
  }
  return kTRUE;
}


//_____________________________________________________________________________
Int_t AliTPCPIDResponse::sectorNumber(Double_t phi) const
{
  //calculate sector number
  const Float_t width=TMath::TwoPi()/18.;
  return TMath::Floor(phi/width);
}


//_____________________________________________________________________________
void AliTPCPIDResponse::Print(Option_t* /*option*/) const
{
  //Print info
  fResponseFunctions.Print();
}


//_____________________________________________________________________________
AliTPCPIDResponse::EChamberStatus AliTPCPIDResponse::TrackStatus(const AliVTrack* track, Int_t layer) const
{
  //status of the track: if it crosses any bad chambers return kChamberOff;
  //IROC:layer=1, OROC:layer=2
  if (layer<1 || layer>2) layer=1;
  Int_t in=0;
  Int_t out=0;
  Float_t inphi=0.;
  Float_t outphi=0.;
  Float_t innerRadius = (layer==1)?83.0:133.7;
  Float_t outerRadius = (layer==1)?133.5:247.7;

  /////////////////////////////////////////////////////////////////////////////
  //find out where track enters and leaves the layer.
  //
  Double_t trackPositionInner[3];
  Double_t trackPositionOuter[3];
 
  //if there is no inner param this could mean we're using the AOD track,
  //we still can extrapolate from the vertex - so use those params.
  const AliExternalTrackParam* ip = track->GetInnerParam();
  if (ip) track=ip;

  Bool_t trackAtInner = track->GetXYZAt(innerRadius, fMagField, trackPositionInner);
  Bool_t trackAtOuter = track->GetXYZAt(outerRadius, fMagField, trackPositionOuter);

  if (!trackAtInner)
  {
    //if we dont even enter inner radius we do nothing and return invalid
    inphi=0.0;
    outphi=0.0;
    in=0;
    out=0;
    return kChamberInvalid;
  }

  if (!trackAtOuter)
  {
    //if we don't reach the outer radius check that the apex is indeed within the outer radius and use apex position
    Bool_t haveApex = TrackApex(track, fMagField, trackPositionOuter);
    Float_t apexRadius = TMath::Sqrt(trackPositionOuter[0]*trackPositionOuter[0]+trackPositionOuter[1]*trackPositionOuter[1]);
    if ( haveApex && apexRadius<=outerRadius && apexRadius>innerRadius)
    {
      //printf("pt: %.2f, apexRadius: %.2f(%s), x: %.2f, y: %.2f\n",track->Pt(),apexRadius,(haveApex)?"OK":"BAD",trackPositionOuter[0],trackPositionOuter[1]);
    }
    else
    {
      inphi=0.0;
      outphi=0.0;
      in=0;
      out=0;
      return kChamberInvalid;
    }
  }


  if (!sectorNumbersInOut(trackPositionInner,
                          trackPositionOuter,
                          inphi,
                          outphi,
                          in,
                          out)) return kChamberInvalid;

  /////////////////////////////////////////////////////////////////////////////
  //now we have the location of the track we can check
  //if it is in a good/bad chamber
  //
  Bool_t sideA = kTRUE;
 
  if (((in/18)%2==1) && ((out/18)%2==1)) sideA=kFALSE;

  in=in%18;
  out=out%18;

  if (TMath::Abs(in-out)>9)
  {
    if (TMath::Max(in,out)==out)
    {
      Int_t tmp=out;
      out = in;
      in = tmp;
      Float_t tmpphi=outphi;
      outphi=inphi;
      inphi=tmpphi;
    }
    in-=18;
    inphi-=TMath::TwoPi();
  }
  else
  {
    if (TMath::Max(in,out)==in)
    {
      Int_t tmp=out;
      out=in;
      in=tmp;
      Float_t tmpphi=outphi;
      outphi=inphi;
      inphi=tmpphi;
    }
  }

  Float_t trackLengthInBad = 0.;
  Float_t trackLengthInLowGain = 0.;
  Float_t trackLengthTotal = TMath::Abs(outphi-inphi);
  Float_t lengthFractionInBadSectors = 0.;

  const Float_t sectorWidth = TMath::TwoPi()/18.;  
 
  for (Int_t i=in; i<=out; i++)
  {
    int j=i;
    if (i<0) j+=18;    //correct for the negative values
    if (!sideA) j+=18; //move to the correct side
   
    Float_t deltaPhi = 0.;
    Float_t phiEdge=sectorWidth*i;
    if (inphi>phiEdge) {deltaPhi=phiEdge+sectorWidth-inphi;}
    else if ((outphi>=phiEdge) && (outphi<(phiEdge+sectorWidth))) {deltaPhi=outphi-phiEdge;}
    else {deltaPhi=sectorWidth;}
   
    Float_t v = fVoltageMap[(layer==1)?(j):(j+36)];
    if (v<=fBadIROCthreshhold)
    {
      trackLengthInBad+=deltaPhi;
      lengthFractionInBadSectors=1.;
    }
    if (v<=fLowGainIROCthreshold && v>fBadIROCthreshhold)
    {
      trackLengthInLowGain+=deltaPhi;
      lengthFractionInBadSectors=1.;
    }
  }

  //for now low gain and bad (off) chambers are treated equally
  if (trackLengthTotal>0)
    lengthFractionInBadSectors = (trackLengthInLowGain+trackLengthInBad)/trackLengthTotal;

  //printf("### side: %s, pt: %.2f, pz: %.2f, in: %i, out: %i, phiIN: %.2f, phiOUT: %.2f, rIN: %.2f, rOUT: %.2f\n",(sideA)?"A":"C",track->Pt(),track->Pz(),in,out,inphi,outphi,innerRadius,outerRadius);
 
  if (lengthFractionInBadSectors>fMaxBadLengthFraction)
  {
    //printf("%%%%%%%% %s kChamberLowGain\n",(layer==1)?"IROC":"OROC");
    return kChamberLowGain;
  }
 
  return kChamberHighGain;
}


//_____________________________________________________________________________
Float_t AliTPCPIDResponse::MaxClusterRadius(const AliVTrack* track) const
{
  //return the radius of the outermost padrow containing a cluster in TPC
  //for the track
  const TBits* clusterMap=track->GetTPCClusterMapPtr();
  if (!clusterMap) return 0.;

  //from AliTPCParam, radius of first IROC row
  const Float_t rfirstIROC = 8.52249984741210938e+01;
  const Float_t padrowHeightIROC = 0.75;
  const Float_t rfirstOROC0 = 1.35100006103515625e+02;
  const Float_t padrowHeightOROC0 = 1.0;
  const Float_t rfirstOROC1 = 1.99350006103515625e+02;
  const Float_t padrowHeightOROC1 = 1.5;

  Int_t maxPadRow=160;
  while ((--maxPadRow)>0 && !clusterMap->TestBitNumber(maxPadRow)){}
  if (maxPadRow>126) return rfirstOROC1+(maxPadRow-126-1)*padrowHeightOROC1;
  if (maxPadRow>62) return rfirstOROC0+(maxPadRow-62-1)*padrowHeightOROC0;
  if (maxPadRow>0) return rfirstIROC+(maxPadRow-1)*padrowHeightIROC;
  return 0.0;
}


//_____________________________________________________________________________
Bool_t AliTPCPIDResponse::TrackApex(const AliVTrack* track, Float_t magField, Double_t position[3]) const
{
  //calculate the coordinates of the apex of the track
  Double_t x[3];
  track->GetXYZ(x);
  Double_t p[3];
  track->GetPxPyPz(p);
  Double_t r = 1./track->OneOverPt()/0.0299792458/magField; //signed - will determine the direction of b
  //printf("b: %.2f, x:%.2f, y:%.2f, pt: %.2f, px:%.2f, py%.2f, r: %.2f\n",magField, x[0],x[1],track->Pt(), p[0],p[1],r);
  //find orthogonal vector (Gram-Schmidt)
  Double_t alpha = (p[0]*x[0] + p[1]*x[1])/(p[0]*p[0] + p[1]*p[1]);
  Double_t b[2];
  b[0]=x[0]-alpha*p[0];
  b[1]=x[1]-alpha*p[1];
 
  Double_t norm = TMath::Sqrt(b[0]*b[0]+b[1]*b[1]);
  if (TMath::AreEqualAbs(norm,0.0,1e-10)) return kFALSE;
  b[0]/=norm;
  b[1]/=norm;
  b[0]*=r;
  b[1]*=r;
  b[0]+=x[0];
  b[1]+=x[1];
  //printf("center: x:%.2f, y:%.2f\n",b[0],b[1]);
 
  norm = TMath::Sqrt(b[0]*b[0]+b[1]*b[1]);
  if (TMath::AreEqualAbs(norm,0.0,1e-10)) return kFALSE;
 
  position[0]=b[0]+b[0]*TMath::Abs(r)/norm;
  position[1]=b[1]+b[1]*TMath::Abs(r)/norm;
  position[2]=0.;
  return kTRUE;
}
