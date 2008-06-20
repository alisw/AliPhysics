/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id:$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Implementation of the class to store the parameters used in   //
// the simulation of SPD, SDD and SSD detectors                  //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "AliITSSimuParam.h"
#include <TMath.h>

const Float_t  AliITSSimuParam::fgkSPDBiasVoltageDefault = 18.182;
const Double_t AliITSSimuParam::fgkSPDThreshDefault = 3000.;
const Double_t AliITSSimuParam::fgkSPDSigmaDefault = 250.;
const TString  AliITSSimuParam::fgkSPDCouplingOptDefault = "old";
const Double_t AliITSSimuParam::fgkSPDCouplColDefault = 0.;
const Double_t AliITSSimuParam::fgkSPDCouplRowDefault = 0.055;
const Float_t  AliITSSimuParam::fgkSPDEccDiffDefault = 0.85;
const Float_t  AliITSSimuParam::fgkSDDDiffCoeffDefault = 3.23;
const Float_t  AliITSSimuParam::fgkSDDDiffCoeff1Default = 30.;
const Float_t  AliITSSimuParam::fgkSDDJitterErrorDefault = 20.; // 20 um from beam test 2001
const Float_t  AliITSSimuParam::fgkSDDDynamicRangeDefault = 132.;
const Int_t    AliITSSimuParam::fgkSDDMaxAdcDefault = 1024;
const Float_t  AliITSSimuParam::fgkSDDChargeLossDefault = 0.;
const Double_t AliITSSimuParam::fgkSSDCouplingPRDefault = 0.01;
const Double_t AliITSSimuParam::fgkSSDCouplingPLDefault = 0.01;
const Double_t AliITSSimuParam::fgkSSDCouplingNRDefault = 0.01;
const Double_t AliITSSimuParam::fgkSSDCouplingNLDefault = 0.01;
const Int_t    AliITSSimuParam::fgkSSDZSThresholdDefault = 3;

const Float_t AliITSSimuParam::fgkNsigmasDefault = 3.;
const Int_t AliITSSimuParam::fgkNcompsDefault = 121;

ClassImp(AliITSSimuParam)

//______________________________________________________________________
AliITSSimuParam::AliITSSimuParam():
  TObject(),
fGeVcharge(0.),
fDOverV(0.),
fSPDBiasVoltage(fgkSPDBiasVoltageDefault),
fSPDThresh(fgkSPDThreshDefault),
fSPDSigma(fgkSPDSigmaDefault),
fSPDCouplOpt(0),
fSPDCouplCol(fgkSPDCouplColDefault),
fSPDCouplRow(fgkSPDCouplRowDefault),
fSPDEccDiff(0.),
fSDDElectronics(0),
fSDDDiffCoeff(0.),
fSDDDiffCoeff1(0.),
fSDDJitterError(fgkSDDJitterErrorDefault),
fSDDDynamicRange(fgkSDDDynamicRangeDefault),
fSDDMaxAdc(0.),
fSDDChargeLoss(fgkSDDChargeLossDefault),
fSSDADCpereV(0.),
fSSDCouplingPR(0),
fSSDCouplingPL(0),
fSSDCouplingNR(0),
fSSDCouplingNL(0),
fSSDZSThreshold(fgkSSDZSThresholdDefault),
fNsigmas(fgkNsigmasDefault),
fNcomps(fgkNcompsDefault),
fGaus()
{  
  // default constructor
  SetGeVToCharge();
  SetDistanceOverVoltage();
  SetSPDCouplingOption(fgkSPDCouplingOptDefault);
  SetSPDSigmaDiffusionAsymmetry(fgkSPDEccDiffDefault);
  SetSDDElectronics();
  SetSDDDiffCoeff(fgkSDDDiffCoeffDefault,fgkSDDDiffCoeff1Default);
  SetSDDMaxAdc((Double_t)fgkSDDMaxAdcDefault);
  SetSSDADCpereV();
  SetSSDCouplings(fgkSSDCouplingPRDefault,fgkSSDCouplingPLDefault,fgkSSDCouplingNRDefault,fgkSSDCouplingNLDefault);
}
//______________________________________________________________________
AliITSSimuParam::AliITSSimuParam(const AliITSSimuParam &simpar):
TObject(),
fGeVcharge(simpar.fGeVcharge),
fDOverV(simpar.fDOverV),
fSPDBiasVoltage(simpar.fSPDBiasVoltage),
fSPDThresh(simpar.fSPDThresh),
fSPDSigma(simpar.fSPDSigma),
fSPDCouplOpt(simpar.fSPDCouplOpt),
fSPDCouplCol(simpar.fSPDCouplCol),
fSPDCouplRow(simpar.fSPDCouplRow),
fSPDEccDiff(simpar.fSPDEccDiff),
fSDDElectronics(simpar.fSDDElectronics),
fSDDDiffCoeff(simpar.fSDDDiffCoeff),
fSDDDiffCoeff1(simpar.fSDDDiffCoeff1),
fSDDJitterError(simpar.fSDDJitterError),
fSDDDynamicRange(simpar.fSDDDynamicRange),
fSDDMaxAdc(simpar.fSDDMaxAdc),
fSDDChargeLoss(simpar.fSDDChargeLoss),
fSSDADCpereV(simpar.fSSDADCpereV),
fSSDCouplingPR(simpar.fSSDCouplingPR),
fSSDCouplingPL(simpar.fSSDCouplingPL),
fSSDCouplingNR(simpar.fSSDCouplingNR),
fSSDCouplingNL(simpar.fSSDCouplingNL),
fSSDZSThreshold(simpar.fSSDZSThreshold),
fNsigmas(simpar.fNsigmas),
fNcomps(simpar.fNcomps),
fGaus(){
  // copy constructor
}

//______________________________________________________________________
AliITSSimuParam& AliITSSimuParam::operator=(const AliITSSimuParam& source){
    // Assignment operator. 
    this->~AliITSSimuParam();
    new(this) AliITSSimuParam(source);
    return *this;

}


//______________________________________________________________________
AliITSSimuParam::~AliITSSimuParam() {
  // destructor
  if(fGaus) delete fGaus;
}
//________________________________________________________________________
void AliITSSimuParam::SetNLookUp(Int_t p1){
  // Set number of sigmas over which cluster disintegration is performed
  fNcomps=p1;
  if (fGaus) delete fGaus;
  fGaus = new TArrayF(fNcomps+1);
  for(Int_t i=0; i<=fNcomps; i++) {
    Float_t x = -fNsigmas + (2.*i*fNsigmas)/(fNcomps-1);
    (*fGaus)[i] = exp(-((x*x)/2));
  }
}
//________________________________________________________________________
void AliITSSimuParam::PrintParameters() const{
  printf("GeVToCharge               = %G\n",fGeVcharge);
  printf("DistanveOverVoltage       = %f \n",fDOverV);
  printf("\n");
  printf("=====  SPD parameters  =====\n");
  printf("Bias Voltage              = %f \n",fSPDBiasVoltage);
  printf("Threshold and sigma       = %f %f\n",fSPDThresh,fSPDSigma);
  printf("Coupling Option           = %s\n",fSPDCouplOpt.Data());
  printf("Coupling value (column)   = %f\n",fSPDCouplCol);
  printf("Coupling value (row)      = %f\n",fSPDCouplRow);
  printf("Eccentricity in diffusion = %f\n",fSPDEccDiff);
  printf("\n");
  printf("=====  SDD parameters  =====\n");
  printf("Electronic chips          = %d\n",fSDDElectronics);
  printf("Diffusion Coefficients    = %f %f\n",fSDDDiffCoeff,fSDDDiffCoeff1);
  printf("Jitter Error              = %f um\n",fSDDJitterError);
  printf("Dynamic Range             = %f\n",fSDDDynamicRange);
  printf("Max. ADC                  = %f\n",fSDDMaxAdc);
  printf("Charge Loss               = %f\n",fSDDChargeLoss);  
  printf("\n");
  printf("=====  SSD parameters  =====\n");
  printf("ADC per eV                = %f\n",fSSDADCpereV);
  printf("Coupling PR               = %f\n",fSSDCouplingPR);
  printf("Coupling PL               = %f\n",fSSDCouplingPL);
  printf("Coupling NR               = %f\n",fSSDCouplingNR);
  printf("Coupling NL               = %f\n",fSSDCouplingNL);
  printf("Zero Supp threshold       = %d\n",fSSDZSThreshold);
}
