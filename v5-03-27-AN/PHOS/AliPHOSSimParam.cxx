/**************************************************************************
 * Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved.      *
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

// Base class for the PHOS simulation parameters.
// Do not use in the simulation; use derivative classes instead.
// Author: Dmitri Peressounko, RRC KI

// --- AliRoot header files ---
#include "AliPHOSSimParam.h"
#include "AliLog.h"

ClassImp(AliPHOSSimParam)

AliPHOSSimParam  * AliPHOSSimParam::fgSimParam = 0 ;
//-----------------------------------------------------------------------------
AliPHOSSimParam::AliPHOSSimParam() :
  TNamed(),
  fLightYieldMean(0.),fIntrinsicAPDEfficiency(0.),
  fLightFactor(0.),fAPDFactor(0.),         
  fAPDNoise(0.),fEMCDigitThreshold(0.),
  fEMCADCchannel(0.),fTOFa(0.),fTOFb(0.),
  fCellNonLineaityA(0.),fCellNonLineaityB(1.),fCellNonLineaityC(1.),
  fEMCSubtractPedestals(kFALSE),
  fGlobalAltroOffset(0),fGlobalAltroThreshold(0),fEMCSampleQualityCut(0),
  fADCpedestalCpv(0.),fADCchanelCpv(0.),
  fCPVNoise(0.),fCPVDigitThreshold(0.),fNADCcpv(0),
  fDigitizeE(0),fCellNonLineaityOn(1)
{
  //Default constructor.
  for(Int_t i=0; i<10; i++) fDStream[i] = 0 ;
}

//-----------------------------------------------------------------------------
AliPHOSSimParam::AliPHOSSimParam(Int_t) :
  TNamed(),
  fLightYieldMean(0.),fIntrinsicAPDEfficiency(0.),
  fLightFactor(0.),fAPDFactor(0.),         
  fAPDNoise(0.),fEMCDigitThreshold(0.),
  fEMCADCchannel(0.),fTOFa(0.),fTOFb(0.),
  fCellNonLineaityA(0.),fCellNonLineaityB(1.),fCellNonLineaityC(1.),
  fEMCSubtractPedestals(kFALSE),
  fGlobalAltroOffset(0),fGlobalAltroThreshold(0),fEMCSampleQualityCut(0),
  fADCpedestalCpv(0.),fADCchanelCpv(0.),
  fCPVNoise(0.),fCPVDigitThreshold(0.),
  fNADCcpv(0),
  fDigitizeE(0),fCellNonLineaityOn(1)
{
  //Real (private) constructor 
  //Set default parameters

  //Parameters describing energy deposition and light collection by APD, used in AliPHOSv1
  //Photoelectron statistics:
  // The light yield is a poissonian distribution of the number of
  // photons created in the PbWo4 crystal, calculated using following formula
  // NumberOfPhotons = EnergyLost * LightYieldMean* APDEfficiency 
  // LightYieldMean is parameter calculated to be over 47000 photons per GeV
  // APDEfficiency is 0.02655
  // k_0 is 0.0045 from Valery Antonenko
  // The number of electrons created in the APD is
  // NumberOfElectrons = APDGain * LightYield
  // The APD Gain is 300
  fLightYieldMean = 47000;            //Average number of photoelectrons per GeV
  fIntrinsicAPDEfficiency = 0.02655 ; //APD efficiency including geometric coverage
//  fLightYieldAttenuation  = 0.0045 ;  //light attenuation in PWO. Last analysis shows no z-position dependence
//                                      //so we removed this dependence from simulations 
  fLightFactor            = fLightYieldMean * fIntrinsicAPDEfficiency ; //Average number of photons collected by 
                            //APD per GeV deposited energy
  fAPDFactor              = (13.418/fLightYieldMean/100.) * 300. ; //factor relating light yield and APD response
                            //evaluated as (13.418/fLightYieldMean/100) * APDGain ;


  //Parameters defining electronic noise calculation and Digits noise thresholds
  //used in AliPHOSDigitizer
  fAPDNoise           = 0.004 ;  // [GeV]
  fEMCDigitThreshold  = 2.5   ;  // [ADC counts]
  fEMCADCchannel      = 0.005 ;  // [GeV]
  fTOFa               = 0.5e-9 ; // [sec] constant term
  fTOFb               = 1.e-9 ;  // [sec/sqrt(GeV)]] stohastic term
  fCellNonLineaityA   = 0.18 ;   //Amp of non-linearity of cell responce
  fCellNonLineaityB   = 0.109;   //Scale of non-linearity of cell responce
  fCellNonLineaityC   = 0.976;   //Overall calibration

  fADCpedestalCpv     = 0.012 ;  // [aux units]
  fADCchanelCpv       = 0.0012;  // [aux units]    
  fCPVNoise           = 0.01;    // [aux units]
  fCPVDigitThreshold  = 0.09 ;   // [aux units]
  fNADCcpv  =  (Int_t)TMath::Power(2,12) ;

  fGlobalAltroOffset = 10;
  fGlobalAltroThreshold = 5;
  fEMCSampleQualityCut = 4.;

  //Imput streams for merging. If true => this stream contains digits (and thus noise) and not SDigits.
  for(Int_t i=0; i<10; i++){
    fDStream[i] = 0 ;
  }
  fgSimParam = this ;
}

//-----------------------------------------------------------------------------
AliPHOSSimParam::AliPHOSSimParam(const AliPHOSSimParam& ):
  TNamed(),
  fLightYieldMean(0.),fIntrinsicAPDEfficiency(0.),
  fLightFactor(0.),fAPDFactor(0.),         
  fAPDNoise(0.),fEMCDigitThreshold(0.),
  fEMCADCchannel(0.),fTOFa(0.),fTOFb(0.),
  fCellNonLineaityA(0.),fCellNonLineaityB(1.),fCellNonLineaityC(1.),
  fEMCSubtractPedestals(kFALSE),
  fGlobalAltroOffset(0),fGlobalAltroThreshold(0),fEMCSampleQualityCut(1.),
  fADCpedestalCpv(0.),fADCchanelCpv(0.),
  fCPVNoise(0.),fCPVDigitThreshold(0.),fNADCcpv(0),
  fDigitizeE(0),fCellNonLineaityOn(1)
{
  //Copy constructor.
  AliError("Should not use copy constructor for singleton") ;
  for(Int_t  i=0; i<10; i++){
    fDStream[i] = 0 ;
  }
  fgSimParam = this ;
}
//-----------------------------------------------------------------------------                                                            
AliPHOSSimParam * AliPHOSSimParam::GetInstance(){

  if(!fgSimParam)
    new AliPHOSSimParam(0) ;
  return fgSimParam ;
}
//-----------------------------------------------------------------------------
AliPHOSSimParam& AliPHOSSimParam::operator = (const AliPHOSSimParam& simParam)
{
  //Assignment operator.

  if(this != &simParam) {
    AliError("Should not use operator= for singleton\n") ;
  }

  return *this;
}

