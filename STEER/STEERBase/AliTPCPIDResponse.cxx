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
//-----------------------------------------------------------------

#include <TGraph.h>
#include <TObjArray.h>
#include <TSpline.h>
#include <TBits.h>
#include <TMath.h>

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
  fCurrentResponseFunction(NULL),
  fCurrentdEdx(0.0),
  fCurrentNPoints(0),
  fCurrentGainScenario(kGainScenarioInvalid),
  fMagField(0.)
{
  //
  //  The default constructor
  //
  for (Int_t i=0; i<fgkNumberOfGainScenarios; i++) {fRes0[i]=0.07;fResN2[i]=0.0;}
}

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
  fCurrentResponseFunction(NULL),
  fCurrentdEdx(0.0),
  fCurrentNPoints(0),
  fCurrentGainScenario(kGainScenarioInvalid),
  fMagField(0.)
{
  //
  //  The main constructor
  //
  for (Int_t i=0; i<fgkNumberOfGainScenarios; i++) {fRes0[i]=param[1];fResN2[i]=param[2];}
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
  fCurrentResponseFunction(NULL),
  fCurrentdEdx(0.0),
  fCurrentNPoints(0),
  fCurrentGainScenario(kGainScenarioInvalid),
  fMagField(that.fMagField)
{
  //copy ctor
  for (Int_t i=0; i<fgkNumberOfGainScenarios; i++) {fRes0[i]=that.fRes0[i];fResN2[i]=that.fResN2[i];}
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
  // Calculates the expected PID signal as the function of 
  // the information stored in the track, for the specified particle type 
  //  
  // At the moment, these signals are just the results of calling the 
  // Bethe-Bloch formula. 
  // This can be improved. By taking into account the number of
  // assigned clusters and/or the track dip angle, for example.  
  //
  
  Double_t mass=AliPID::ParticleMassZ(n);
  if (!fUseDatabase) return Bethe(mom/mass);
  //
  const TSpline3 * responseFunction = (TSpline3 *) fResponseFunctions.UncheckedAt(n);

  if (!responseFunction) return Bethe(mom/mass);
  //charge factor. BB goes with z^2, however in reality it is slightly larger (calibration, threshold effects, ...)
  // !!! Splines for light nuclei need to be normalised to this factor !!!
  const Double_t chargeFactor = TMath::Power(AliPID::ParticleCharge(n),2.3);
  return fMIP*responseFunction->Eval(mom/mass)*chargeFactor;

}

//_________________________________________________________________________
Double_t AliTPCPIDResponse::GetExpectedSigma(const Float_t mom, 
			                     const Int_t nPoints,
                                             AliPID::EParticleType n) const {
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
Double_t AliTPCPIDResponse::GetExpectedSignal(Double_t mom,
                                              AliPID::EParticleType species,
                                              const TSpline3* responseFunction) const 
{
  // Calculates the expected PID signal as the function of 
  // the information stored in the track, for the specified particle type 
  //  
  // At the moment, these signals are just the results of calling the 
  // Bethe-Bloch formula. 
  // This can be improved. By taking into account the number of
  // assigned clusters and/or the track dip angle, for example.  
  //
  
  
  Double_t mass=AliPID::ParticleMass(species);
  
  if (!fUseDatabase) return Bethe(mom/mass);
  
  if (!responseFunction) return Bethe(mom/mass);
  return fMIP*responseFunction->Eval(mom/mass);
}

//_________________________________________________________________________
Double_t AliTPCPIDResponse::GetExpectedSignal(const AliVTrack* track,
                                              AliPID::EParticleType species,
                                              ETPCdEdxSource dedxSource) 
{
  // Calculates the expected PID signal as the function of 
  // the information stored in the track, for the specified particle type 
  //  
  // At the moment, these signals are just the results of calling the 
  // Bethe-Bloch formula. 
  // This can be improved. By taking into account the number of
  // assigned clusters and/or the track dip angle, for example.  
  //
  
  Double_t mom = track->GetTPCmomentum();
  Double_t mass=AliPID::ParticleMass(species);
  
  if (!fUseDatabase) return Bethe(mom/mass);
  
  const TSpline3* responseFunction = GetResponseFunction(track,species,dedxSource);
  if (!responseFunction) return Bethe(mom/mass);
  return fMIP*responseFunction->Eval(mom/mass);

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
                               ETPCdEdxSource dedxSource ) 
{
  //the splines are stored in an array, different scenarios

  if (ResponseFunctiondEdxN(track,
                            species,
                            dedxSource))
    return fCurrentResponseFunction;
  
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
Double_t AliTPCPIDResponse::GetExpectedSigma( Double_t mom,
                                              Int_t nPoints,
                                              AliPID::EParticleType species,
                                              ETPCgainScenario gainScenario,
                                              const TSpline3* responseFunction) const
{
  // Calculates the expected sigma of the PID signal as the function of 
  // the information stored in the track, for the specified particle type 
  // and dedx scenario
  //
  
  if (!responseFunction) return 999;
  if (nPoints != 0) 
    return GetExpectedSignal(mom,species,responseFunction) *
           fRes0[gainScenario] * sqrt(1. + fResN2[gainScenario]/nPoints);
  else
    return GetExpectedSignal(mom,species,responseFunction)*fRes0[gainScenario];
}

//_________________________________________________________________________
Double_t AliTPCPIDResponse::GetExpectedSigma(const AliVTrack* track, 
                                             AliPID::EParticleType species,
                                             ETPCdEdxSource dedxSource) 
{
  // Calculates the expected sigma of the PID signal as the function of 
  // the information stored in the track, for the specified particle type 
  // and dedx scenario
  //
  
  if (!ResponseFunctiondEdxN(track,
                             species,
                             dedxSource)) return 998; //TODO: better handling!

  return GetExpectedSigma( track->GetTPCmomentum(),
                           fCurrentNPoints,
                           species,
                           fCurrentGainScenario,
                           fCurrentResponseFunction );
}

//_________________________________________________________________________
Float_t AliTPCPIDResponse::GetNumberOfSigmas(const AliVTrack* track, 
                             AliPID::EParticleType species,
                             ETPCdEdxSource dedxSource) 
{
  //calculates the number of sigmas of the PID signal from the expected value
  //for a gicven particle species in the presence of multiple gain scenarios
  //inside the TPC

  if (!ResponseFunctiondEdxN(track,
                             species,
                             dedxSource)) return 997; //TODO: better handling!

  Double_t mom = track->GetTPCmomentum();
  Double_t bethe=GetExpectedSignal(mom,species,fCurrentResponseFunction);
  Double_t sigma=GetExpectedSigma( mom,
                                   fCurrentNPoints,
                                   species,
                                   fCurrentGainScenario,
                                   fCurrentResponseFunction );
  return (fCurrentdEdx-bethe)/sigma;
}

//_________________________________________________________________________
Bool_t AliTPCPIDResponse::ResponseFunctiondEdxN( const AliVTrack* track, 
                                                 AliPID::EParticleType species,
                                                 ETPCdEdxSource dedxSource ) 
{
  // Calculates the right parameters for PID
  //   dEdx parametrization for the proper gain scenario, dEdx 
  //   and NPoints used for dEdx
  // based on the track geometry (which chambers it crosses) for the specified particle type 
  // and preferred source of dedx.
  // returns true on success
  
  Double32_t signal[4]; //0: IROC, 1: OROC medium, 2:OROC long, 3: OROC all (def. truncation used)
  Char_t ncl[3];        //same
  Char_t nrows[3];      //same
  AliTPCdEdxInfo* dEdxInfo = track->GetTPCdEdxInfo();
  
  if (!dEdxInfo && dedxSource!=kdEdxDefault)  //in one case its ok if we dont have the info
  {
    AliError("AliTPCdEdxInfo not available");
    InvalidateCurrentValues();
    return kFALSE;
  }

  if (dEdxInfo) dEdxInfo->GetTPCSignalRegionInfo(signal,ncl,nrows);

  //check if we cross a bad OROC in which case we reject
  EChamberStatus trackStatus = TrackStatus(track,2);
  if (trackStatus==kChamberOff || trackStatus==kChamberLowGain)
  {
    InvalidateCurrentValues();
    return kFALSE;
  }

  switch (dedxSource)
  {
    case kdEdxDefault:
      {
        fCurrentdEdx = track->GetTPCsignal();
        fCurrentNPoints = track->GetTPCsignalN();
        fCurrentGainScenario = kDefault;
        break;
      }
    case kdEdxOROC:
      {
        fCurrentdEdx = signal[3];
        fCurrentNPoints = ncl[2]+ncl[1];
        fCurrentGainScenario = kOROChigh;
        break;
      }
    case kdEdxHybrid:
      {
        //if we cross a bad IROC we use OROC dedx, if we dont we use combined
        EChamberStatus status = TrackStatus(track,1);
        if (status!=kChamberHighGain)
        {
          fCurrentdEdx = signal[3];
          fCurrentNPoints = ncl[2]+ncl[1];
          fCurrentGainScenario = kOROChigh;
        }
        else
        {
          fCurrentdEdx = track->GetTPCsignal();
          fCurrentNPoints = track->GetTPCsignalN();
          fCurrentGainScenario = kALLhigh;
        }
        break;
      }
    default:
      {
         fCurrentdEdx = 0.;
         fCurrentNPoints = 0;
         fCurrentGainScenario = kGainScenarioInvalid;
         return kFALSE;
      }
  }
  TObject* obj = fResponseFunctions.UncheckedAt(ResponseFunctionIndex(species,fCurrentGainScenario));
  fCurrentResponseFunction = dynamic_cast<TSpline3*>(obj); //TODO:maybe static cast?
  return kTRUE;
}

//_________________________________________________________________________
Bool_t AliTPCPIDResponse::sectorNumbersInOut(const AliVTrack* track,
                                             Double_t innerRadius,
                                             Double_t outerRadius,
                                             Float_t& inphi, 
                                             Float_t& outphi,
                                             Int_t& in, 
                                             Int_t& out ) const
{
  //calculate the sector numbers (equivalent to IROC chamber numbers) a track crosses
  //for OROC chamber numbers add 36
  //returned angles are between (0,2pi)
  
  Double_t trackPositionInner[3]; 
  Double_t trackPositionOuter[3]; 

  Bool_t trackAtInner=kTRUE;
  Bool_t trackAtOuter=kTRUE;
  const AliExternalTrackParam* ip = track->GetInnerParam();
  
  //if there is no inner param this could mean we're using the AOD track,
  //we still can extrapolate from the vertex - so use those params.
  if (ip) track=ip;

  if (!track->GetXYZAt(innerRadius, fMagField, trackPositionInner))
    trackAtInner=kFALSE;

  if (!track->GetXYZAt(outerRadius, fMagField, trackPositionOuter))
    trackAtOuter=kFALSE;

  if (!trackAtInner)
  {
    //if we dont even enter inner radius we do nothing
    inphi=0.0;
    outphi=0.0;
    in=0;
    out=0;
    return kFALSE;
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
      return kFALSE;
    }
  }

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
  if (!sectorNumbersInOut(track, 
                          innerRadius, 
                          outerRadius, 
                          inphi, 
                          outphi, 
                          in, 
                          out)) return kChamberOff;

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
    if (v<=fBadIROCthreshhold) trackLengthInBad+=deltaPhi;
    if (v<=fLowGainIROCthreshold && v>fBadIROCthreshhold) trackLengthInLowGain+=deltaPhi;
  }

  //for now low gain and bad (off) chambers are treated equally
  Float_t lengthFractionInBadSectors = (trackLengthInLowGain+trackLengthInBad)/trackLengthTotal;

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
void AliTPCPIDResponse::InvalidateCurrentValues()
{
  //make the "current" stored values from the last processed track invalid
  fCurrentResponseFunction=NULL;
  fCurrentdEdx=0.;
  fCurrentNPoints=0;
  fCurrentGainScenario=kGainScenarioInvalid;
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

