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

/*
$Log$
Revision 1.4  2000/10/06 16:49:46  cblume
Made Getters const

Revision 1.3.2.1  2000/09/18 13:45:30  cblume
New class AliTRDsim that simulates TR photons

Revision 1.2  1999/09/29 09:24:35  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD simulation - multimodule (regular rad.)                              //
//  after: M. CASTELLANO et al., COMP. PHYS. COMM. 51 (1988) 431             //
//                             + COMP. PHYS. COMM. 61 (1990) 395             //
//                                                                           //
//   17.07.1998 - A.Andronic                                                 //
//   08.12.1998 - simplified version                                         //
//   11.07.2000 - Adapted code to aliroot environment (C.Blume)              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>

#include <TH1.h>
#include <TRandom.h>
#include <TMath.h>
#include <TParticle.h>

#include "AliModule.h"

#include "AliTRDsim.h"

ClassImp(AliTRDsim)

//_____________________________________________________________________________
AliTRDsim::AliTRDsim():TObject()
{
  //
  // AliTRDsim default constructor
  // 

  Init();

}

//_____________________________________________________________________________
AliTRDsim::AliTRDsim(AliModule *mod, Int_t foil, Int_t gap)
{
  //
  // AliTRDsim constructor. Takes the material properties of the radiator
  // foils and the gas in the gaps from AliModule <mod>.
  // The default number of foils is 100 with a thickness of 20 mu. The 
  // thickness of the gaps is 500 mu.
  //

  Float_t aFoil, zFoil, rhoFoil;
  Float_t aGap,  zGap,  rhoGap;
  Float_t rad, abs;
  Char_t  name[21];

  Init();

  mod->AliGetMaterial(foil,name,aFoil,zFoil,rhoFoil,rad,abs);
  mod->AliGetMaterial(gap ,name,aGap ,zGap ,rhoGap ,rad,abs);

  fFoilDens  = rhoFoil;
  fFoilA     = aFoil;
  fFoilZ     = zFoil;
  fFoilOmega = Omega(fFoilDens,fFoilZ,fFoilA);

  fGapDens   = rhoGap;
  fGapA      = aGap;
  fGapZ      = zGap;
  fGapOmega  = Omega(fGapDens ,fGapZ ,fGapA );

}

//_____________________________________________________________________________
AliTRDsim::AliTRDsim(const AliTRDsim &s)
{
  //
  // AliTRDsim copy constructor
  //

  ((AliTRDsim &) s).Copy(*this);

}

//_____________________________________________________________________________
AliTRDsim::~AliTRDsim() 
{
  //
  // AliTRDsim destructor
  //

  if (fSpectrum) delete fSpectrum;
  if (fSigma)    delete fSigma;

}

//_____________________________________________________________________________
AliTRDsim &AliTRDsim::operator=(const AliTRDsim &s)
{
  //
  // Assignment operator
  //

  if (this != &s) ((AliTRDsim &) s).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDsim::Copy(TObject &s)
{
  //
  // Copy function
  //

  ((AliTRDsim &) s).fNFoils     = fNFoils;
  ((AliTRDsim &) s).fFoilThick  = fFoilThick;
  ((AliTRDsim &) s).fFoilDens   = fFoilDens;
  ((AliTRDsim &) s).fFoilOmega  = fFoilOmega;
  ((AliTRDsim &) s).fFoilZ      = fFoilZ;
  ((AliTRDsim &) s).fFoilA      = fFoilA;
  ((AliTRDsim &) s).fGapThick   = fGapThick;
  ((AliTRDsim &) s).fGapDens    = fGapDens;
  ((AliTRDsim &) s).fGapOmega   = fGapOmega;
  ((AliTRDsim &) s).fGapZ       = fGapZ;
  ((AliTRDsim &) s).fGapA       = fGapA;
  ((AliTRDsim &) s).fTemp       = fTemp;
  ((AliTRDsim &) s).fSpNBins    = fSpNBins;
  ((AliTRDsim &) s).fSpRange    = fSpRange;
  ((AliTRDsim &) s).fSpBinWidth = fSpBinWidth;
  ((AliTRDsim &) s).fSpLower    = fSpLower;
  ((AliTRDsim &) s).fSpUpper    = fSpUpper;

  if (((AliTRDsim &) s).fSigma) delete ((AliTRDsim &) s).fSigma;
  ((AliTRDsim &) s).fSigma = new Double_t[fSpNBins];
  for (Int_t iBin = 0; iBin < fSpNBins; iBin++) {
    ((AliTRDsim &) s).fSigma[iBin] = fSigma[iBin];
  }  

  fSpectrum->Copy(*((AliTRDsim &) s).fSpectrum);

}

//_____________________________________________________________________________
void AliTRDsim::Init()
{
  //
  // Initialization 
  // The default radiator are 100 prolypropilene foils of 20 mu thickness
  // with gaps of 500 mu filled with CO2.
  //      
  // 

  fNFoils     = 100;

  fFoilThick  = 0.0020;
  fFoilDens   = 0.92;   
  fFoilZ      = 5.28571;
  fFoilA      = 10.4286;
  fFoilOmega  = Omega(fFoilDens,fFoilZ,fFoilA);

  fGapThick   = 0.0500;
  fGapDens    = 0.001977;  
  fGapZ       = 7.45455;
  fGapA       = 14.9091;
  fGapOmega   = Omega(fGapDens ,fGapZ ,fGapA );

  fTemp       = 293.16;

  fSpNBins    = 200;
  fSpRange    = 100;
  fSpBinWidth = fSpRange / fSpNBins;
  fSpLower    = 1.0 - 0.5 * fSpBinWidth;
  fSpUpper    = fSpLower + fSpRange;

  if (fSpectrum) delete fSpectrum;
  fSpectrum   = new TH1D("TRspectrum","TR spectrum",fSpNBins,fSpLower,fSpUpper);

  // Set the sigma values 
  SetSigma();

}

//_____________________________________________________________________________
Int_t AliTRDsim::CreatePhotons(Int_t pdg, Float_t p
                             , Int_t &nPhoton, Float_t *ePhoton)
{
  //
  // Create TRD photons for a charged particle of type <pdg> with the total 
  // momentum <p>. 
  // Number of produced TR photons:       <nPhoton>
  // Energies of the produced TR photons: <ePhoton>
  //

  // PDG codes
  const Int_t kPdgEle  =  11;
  const Int_t kPdgMuon =  13;
  const Int_t kPdgPion = 211;
  const Int_t kPdgKaon = 321;

  Float_t  mass        = 0;
  switch (TMath::Abs(pdg)) {
  case kPdgEle:
    mass      =  5.11e-4;
    break;
  case kPdgMuon:
    mass      =  0.10566;
    break;
  case kPdgPion:
    mass      =  0.13957;
    break;
  case kPdgKaon:
    mass      =  0.4937;
    break;
  default:
    return 0;
    break;
  };

  // Calculate gamma
  Double_t gamma = TMath::Sqrt(p*p + mass*mass) / mass;

  // Calculate the TR photons
  return TrPhotons(gamma, nPhoton, ePhoton);

}

//_____________________________________________________________________________
Int_t AliTRDsim::TrPhotons(Double_t gamma, Int_t &nPhoton, Float_t *ePhoton)
{
  //
  // Produces TR photons.
  //

  const Double_t kAlpha  = 0.0072973;
  const Int_t    kSumMax = 10;

  Double_t kappa = fGapThick / fFoilThick;

  fSpectrum->Reset();

  // The TR spectrum
  Double_t stemp = 0;
  for (Int_t iBin = 0; iBin < fSpNBins; iBin++) {

    // keV -> eV
    Double_t energyeV = (fSpBinWidth * iBin + 1.0) * 1e3;

    Double_t csFoil   = fFoilOmega / energyeV;
    Double_t csGap    = fGapOmega  / energyeV;

    Double_t rho1     = energyeV * fFoilThick * 1e4 * 2.5 
                                 * (1.0 / (gamma*gamma) + csFoil*csFoil);
    Double_t rho2     = energyeV * fFoilThick * 1e4 * 2.5 
                                 * (1.0 / (gamma*gamma) + csGap *csGap);

    // Calculate the sum
    Double_t sum = 0;
    for (Int_t iSum = 0; iSum < kSumMax; iSum++) {
      Double_t tetan = (TMath::Pi() * 2.0 * (iSum+1) - (rho1 + kappa * rho2)) 
                     / (kappa + 1.0);
      if (tetan < 0.0) tetan = 0.0;
      Double_t aux   = 1.0 / (rho1 + tetan) - 1.0 / (rho2 + tetan);
               sum  += tetan * (aux*aux) * (1.0 - TMath::Cos(rho1 + tetan));
    }

    // Absorbtion
    Double_t conv      = 1.0 - TMath::Exp(-fNFoils * fSigma[iBin]);

    // eV -> keV
    Float_t  energykeV = energyeV * 0.001;

    // dN / domega
    Double_t wn        = kAlpha * 4.0 / (fSigma[iBin] * (kappa + 1.0)) 
                                * conv * sum / energykeV;
    fSpectrum->SetBinContent(iBin,wn);

    stemp += wn;

  }

  // <nTR> (binsize corr.)
  Float_t ntr = stemp * fSpBinWidth;
  // Number of TR photons from Poisson distribution with mean <ntr>
  nPhoton = gRandom->Poisson(ntr);
  // Energy of the TR photons
  for (Int_t iPhoton = 0; iPhoton < nPhoton; iPhoton++) {
    ePhoton[iPhoton] = fSpectrum->GetRandom();
  }

  return 1;

}

//_____________________________________________________________________________
void AliTRDsim::SetSigma() 
{
  //
  // Sets the absorbtion crosssection for the energies of the TR spectrum
  //

  if (fSigma) delete fSigma;
  fSigma = new Double_t[fSpNBins];
  for (Int_t iBin = 0; iBin < fSpNBins; iBin++) {
    Double_t energykeV = iBin * fSpBinWidth + 1.0;
    fSigma[iBin]       = Sigma(energykeV);
    //printf("SetSigma(): iBin = %d fSigma %g\n",iBin,fSigma[iBin]);
  }

}

//_____________________________________________________________________________
Double_t AliTRDsim::Sigma(Double_t energykeV)
{
  //
  // Calculates the absorbtion crosssection for a one-foil-one-gap-radiator
  //

  // Gas at 0 C
  const Double_t kTemp0 = 273.16;

  // keV -> MeV
  Double_t energyMeV = energykeV * 0.001;
  if (energyMeV >= 0.001) {
    return(GetMuPo(energyMeV) * fFoilDens * fFoilThick + 
           GetMuCO(energyMeV) * fGapDens  * fGapThick  * fTemp/kTemp0);
  }
  else {
    return 1e6;
  }

}

//_____________________________________________________________________________
Double_t AliTRDsim::GetMuPo(Double_t energyMeV)
{
  //
  // Returns the photon absorbtion cross section for polypropylene
  //

  const Int_t kN = 36;

  Double_t mu[kN] = { 1.894E+03, 5.999E+02, 2.593E+02
                    , 7.743E+01, 3.242E+01, 1.643E+01
                    , 9.432E+00, 3.975E+00, 2.088E+00
                    , 7.452E-01, 4.315E-01, 2.706E-01
                    , 2.275E-01, 2.084E-01, 1.970E-01
                    , 1.823E-01, 1.719E-01, 1.534E-01
                    , 1.402E-01, 1.217E-01, 1.089E-01
                    , 9.947E-02, 9.198E-02, 8.078E-02
                    , 7.262E-02, 6.495E-02, 5.910E-02   
                    , 5.064E-02, 4.045E-02, 3.444E-02
                    , 3.045E-02, 2.760E-02, 2.383E-02
		    , 2.145E-02, 1.819E-02, 1.658E-02 };

  Double_t en[kN] = { 1.000E-03, 1.500E-03, 2.000E-03
                    , 3.000E-03, 4.000E-03, 5.000E-03
                    , 6.000E-03, 8.000E-03, 1.000E-02
                    , 1.500E-02, 2.000E-02, 3.000E-02
                    , 4.000E-02, 5.000E-02, 6.000E-02
                    , 8.000E-02, 1.000E-01, 1.500E-01
                    , 2.000E-01, 3.000E-01, 4.000E-01
                    , 5.000E-01, 6.000E-01, 8.000E-01
                    , 1.000E+00, 1.250E+00, 1.500E+00
                    , 2.000E+00, 3.000E+00, 4.000E+00
                    , 5.000E+00, 6.000E+00, 8.000E+00
		    , 1.000E+01, 1.500E+01, 2.000E+01 };

  return Interpolate(energyMeV,en,mu,kN);

}

//_____________________________________________________________________________
Double_t AliTRDsim::GetMuCO(Double_t energyMeV)
{
  //
  // Returns the photon absorbtion cross section for CO2
  //

  const Int_t kN = 36;

  Double_t mu[kN] = { 0.39383E+04, 0.13166E+04, 0.58750E+03
                    , 0.18240E+03, 0.77996E+02, 0.40024E+02
                    , 0.23116E+02, 0.96997E+01, 0.49726E+01
                    , 0.15543E+01, 0.74915E+00, 0.34442E+00
                    , 0.24440E+00, 0.20589E+00, 0.18632E+00
                    , 0.16578E+00, 0.15394E+00, 0.13558E+00
                    , 0.12336E+00, 0.10678E+00, 0.95510E-01
                    , 0.87165E-01, 0.80587E-01, 0.70769E-01
                    , 0.63626E-01, 0.56894E-01, 0.51782E-01
                    , 0.44499E-01, 0.35839E-01, 0.30825E-01
                    , 0.27555E-01, 0.25269E-01, 0.22311E-01
		    , 0.20516E-01, 0.18184E-01, 0.17152E-01 };

  Double_t en[kN] = { 0.10000E-02, 0.15000E-02, 0.20000E-02
                    , 0.30000E-02, 0.40000E-02, 0.50000E-02
                    , 0.60000E-02, 0.80000E-02, 0.10000E-01
                    , 0.15000E-01, 0.20000E-01, 0.30000E-01
                    , 0.40000E-01, 0.50000E-01, 0.60000E-01
                    , 0.80000E-01, 0.10000E+00, 0.15000E+00
                    , 0.20000E+00, 0.30000E+00, 0.40000E+00
                    , 0.50000E+00, 0.60000E+00, 0.80000E+00
                    , 0.10000E+01, 0.12500E+01, 0.15000E+01
                    , 0.20000E+01, 0.30000E+01, 0.40000E+01
                    , 0.50000E+01, 0.60000E+01, 0.80000E+01
		    , 0.10000E+02, 0.15000E+02, 0.20000E+02 };

  return Interpolate(energyMeV,en,mu,kN);

}

//_____________________________________________________________________________
Double_t AliTRDsim::GetMuXe(Double_t energyMeV)
{
  //
  // Returns the photon absorbtion cross section for xenon
  //

  const Int_t kN = 48;

  Double_t mu[kN] = { 9.413E+03, 8.151E+03, 7.035E+03
                    , 7.338E+03, 4.085E+03, 2.088E+03
                    , 7.780E+02, 3.787E+02, 2.408E+02
                    , 6.941E+02, 6.392E+02, 6.044E+02
                    , 8.181E+02, 7.579E+02, 6.991E+02
                    , 8.064E+02, 6.376E+02, 3.032E+02
                    , 1.690E+02, 5.743E+01, 2.652E+01
                    , 8.930E+00, 6.129E+00, 3.316E+01
                    , 2.270E+01, 1.272E+01, 7.825E+00
                    , 3.633E+00, 2.011E+00, 7.202E-01
                    , 3.760E-01, 1.797E-01, 1.223E-01
                    , 9.699E-02, 8.281E-02, 6.696E-02
                    , 5.785E-02, 5.054E-02, 4.594E-02
                    , 4.078E-02, 3.681E-02, 3.577E-02
                    , 3.583E-02, 3.634E-02, 3.797E-02
		    , 3.987E-02, 4.445E-02, 4.815E-02 };

  Double_t en[kN] = { 1.00000E-03, 1.07191E-03, 1.14900E-03
                    , 1.14900E-03, 1.50000E-03, 2.00000E-03
                    , 3.00000E-03, 4.00000E-03, 4.78220E-03
                    , 4.78220E-03, 5.00000E-03, 5.10370E-03
                    , 5.10370E-03, 5.27536E-03, 5.45280E-03
                    , 5.45280E-03, 6.00000E-03, 8.00000E-03
                    , 1.00000E-02, 1.50000E-02, 2.00000E-02
                    , 3.00000E-02, 3.45614E-02, 3.45614E-02
                    , 4.00000E-02, 5.00000E-02, 6.00000E-02
                    , 8.00000E-02, 1.00000E-01, 1.50000E-01
                    , 2.00000E-01, 3.00000E-01, 4.00000E-01
                    , 5.00000E-01, 6.00000E-01, 8.00000E-01
                    , 1.00000E+00, 1.25000E+00, 1.50000E+00
                    , 2.00000E+00, 3.00000E+00, 4.00000E+00
                    , 5.00000E+00, 6.00000E+00, 8.00000E+00
		    , 1.00000E+01, 1.50000E+01, 2.00000E+01 };

  return Interpolate(energyMeV,en,mu,kN);

}

//_____________________________________________________________________________
Double_t AliTRDsim::GetMuBu(Double_t energyMeV)
{
  //
  // Returns the photon absorbtion cross section for isobutane
  //

  const Int_t kN = 36;

  Double_t mu[kN] = { 0.38846E+03, 0.12291E+03, 0.53225E+02
                    , 0.16091E+02, 0.69114E+01, 0.36541E+01
                    , 0.22282E+01, 0.11149E+01, 0.72887E+00
                    , 0.45053E+00, 0.38167E+00, 0.33920E+00
                    , 0.32155E+00, 0.30949E+00, 0.29960E+00
                    , 0.28317E+00, 0.26937E+00, 0.24228E+00
                    , 0.22190E+00, 0.19289E+00, 0.17288E+00
                    , 0.15789E+00, 0.14602E+00, 0.12829E+00
                    , 0.11533E+00, 0.10310E+00, 0.93790E-01
                    , 0.80117E-01, 0.63330E-01, 0.53229E-01
                    , 0.46390E-01, 0.41425E-01, 0.34668E-01
		    , 0.30267E-01, 0.23910E-01, 0.20509E-01 };

  Double_t en[kN] = { 0.10000E-02, 0.15000E-02, 0.20000E-02
                    , 0.30000E-02, 0.40000E-02, 0.50000E-02
                    , 0.60000E-02, 0.80000E-02, 0.10000E-01
                    , 0.15000E-01, 0.20000E-01, 0.30000E-01
                    , 0.40000E-01, 0.50000E-01, 0.60000E-01
                    , 0.80000E-01, 0.10000E+00, 0.15000E+00
                    , 0.20000E+00, 0.30000E+00, 0.40000E+00
                    , 0.50000E+00, 0.60000E+00, 0.80000E+00
                    , 0.10000E+01, 0.12500E+01, 0.15000E+01
                    , 0.20000E+01, 0.30000E+01, 0.40000E+01
                    , 0.50000E+01, 0.60000E+01, 0.80000E+01
		    , 0.10000E+02, 0.15000E+02, 0.20000E+02 };

  return Interpolate(energyMeV,en,mu,kN);

}

//_____________________________________________________________________________
Double_t AliTRDsim::GetMuMy(Double_t energyMeV)
{
  //
  // Returns the photon absorbtion cross section for mylar
  //

  const Int_t kN = 36;

  Double_t mu[kN] = { 2.911E+03, 9.536E+02, 4.206E+02
                    , 1.288E+02, 5.466E+01, 2.792E+01
                    , 1.608E+01, 6.750E+00, 3.481E+00
                    , 1.132E+00, 5.798E-01, 3.009E-01
                    , 2.304E-01, 2.020E-01, 1.868E-01
                    , 1.695E-01, 1.586E-01, 1.406E-01
                    , 1.282E-01, 1.111E-01, 9.947E-02
                    , 9.079E-02, 8.395E-02, 7.372E-02
                    , 6.628E-02, 5.927E-02, 5.395E-02
                    , 4.630E-02, 3.715E-02, 3.181E-02
                    , 2.829E-02, 2.582E-02, 2.257E-02
                    , 2.057E-02, 1.789E-02, 1.664E-02 };

  Double_t en[kN] = { 1.00000E-03, 1.50000E-03, 2.00000E-03
                    , 3.00000E-03, 4.00000E-03, 5.00000E-03
                    , 6.00000E-03, 8.00000E-03, 1.00000E-02
                    , 1.50000E-02, 2.00000E-02, 3.00000E-02
                    , 4.00000E-02, 5.00000E-02, 6.00000E-02
                    , 8.00000E-02, 1.00000E-01, 1.50000E-01
                    , 2.00000E-01, 3.00000E-01, 4.00000E-01
                    , 5.00000E-01, 6.00000E-01, 8.00000E-01
                    , 1.00000E+00, 1.25000E+00, 1.50000E+00
                    , 2.00000E+00, 3.00000E+00, 4.00000E+00
                    , 5.00000E+00, 6.00000E+00, 8.00000E+00
                    , 1.00000E+01, 1.50000E+01, 2.00000E+01 };

  return Interpolate(energyMeV,en,mu,kN);

}

//_____________________________________________________________________________
Double_t AliTRDsim::GetMuN2(Double_t energyMeV)
{
  //
  // Returns the photon absorbtion cross section for nitrogen
  //

  const Int_t kN = 36;

  Double_t mu[kN] = { 3.311E+03, 1.083E+03, 4.769E+02
                    , 1.456E+02, 6.166E+01, 3.144E+01
                    , 1.809E+01, 7.562E+00, 3.879E+00
                    , 1.236E+00, 6.178E-01, 3.066E-01
                    , 2.288E-01, 1.980E-01, 1.817E-01
                    , 1.639E-01, 1.529E-01, 1.353E-01
                    , 1.233E-01, 1.068E-01, 9.557E-02
                    , 8.719E-02, 8.063E-02, 7.081E-02
                    , 6.364E-02, 5.693E-02, 5.180E-02
                    , 4.450E-02, 3.579E-02, 3.073E-02
                    , 2.742E-02, 2.511E-02, 2.209E-02
                    , 2.024E-02, 1.782E-02, 1.673E-02 };

  Double_t en[kN] = { 1.00000E-03, 1.50000E-03, 2.00000E-03
                    , 3.00000E-03, 4.00000E-03, 5.00000E-03
                    , 6.00000E-03, 8.00000E-03, 1.00000E-02
                    , 1.50000E-02, 2.00000E-02, 3.00000E-02
                    , 4.00000E-02, 5.00000E-02, 6.00000E-02
                    , 8.00000E-02, 1.00000E-01, 1.50000E-01
                    , 2.00000E-01, 3.00000E-01, 4.00000E-01
                    , 5.00000E-01, 6.00000E-01, 8.00000E-01
                    , 1.00000E+00, 1.25000E+00, 1.50000E+00
                    , 2.00000E+00, 3.00000E+00, 4.00000E+00
                    , 5.00000E+00, 6.00000E+00, 8.00000E+00
                    , 1.00000E+01, 1.50000E+01, 2.00000E+01 };

  return Interpolate(energyMeV,en,mu,kN);

}

//_____________________________________________________________________________
Double_t AliTRDsim::GetMuO2(Double_t energyMeV)
{
  //
  // Returns the photon absorbtion cross section for oxygen
  //

  const Int_t kN = 36;

  Double_t mu[kN] = { 4.590E+03, 1.549E+03, 6.949E+02
                    , 2.171E+02, 9.315E+01, 4.790E+01
                    , 2.770E+01, 1.163E+01, 5.952E+00
                    , 1.836E+00, 8.651E-01, 3.779E-01
                    , 2.585E-01, 2.132E-01, 1.907E-01
                    , 1.678E-01, 1.551E-01, 1.361E-01
                    , 1.237E-01, 1.070E-01, 9.566E-02
                    , 8.729E-02, 8.070E-02, 7.087E-02
                    , 6.372E-02, 5.697E-02, 5.185E-02
                    , 4.459E-02, 3.597E-02, 3.100E-02
                    , 2.777E-02, 2.552E-02, 2.263E-02
                    , 2.089E-02, 1.866E-02, 1.770E-02 };

  Double_t en[kN] = { 1.00000E-03, 1.50000E-03, 2.00000E-03
                    , 3.00000E-03, 4.00000E-03, 5.00000E-03
                    , 6.00000E-03, 8.00000E-03, 1.00000E-02
                    , 1.50000E-02, 2.00000E-02, 3.00000E-02
                    , 4.00000E-02, 5.00000E-02, 6.00000E-02
                    , 8.00000E-02, 1.00000E-01, 1.50000E-01
                    , 2.00000E-01, 3.00000E-01, 4.00000E-01
                    , 5.00000E-01, 6.00000E-01, 8.00000E-01
                    , 1.00000E+00, 1.25000E+00, 1.50000E+00
                    , 2.00000E+00, 3.00000E+00, 4.00000E+00
                    , 5.00000E+00, 6.00000E+00, 8.00000E+00
                    , 1.00000E+01, 1.50000E+01, 2.00000E+01 };

  return Interpolate(energyMeV,en,mu,kN);

}

//_____________________________________________________________________________
Double_t AliTRDsim::GetMuHe(Double_t energyMeV)
{
  //
  // Returns the photon absorbtion cross section for helium
  //

  const Int_t kN = 36;

  Double_t mu[kN] = { 6.084E+01, 1.676E+01, 6.863E+00
                    , 2.007E+00, 9.329E-01, 5.766E-01
                    , 4.195E-01, 2.933E-01, 2.476E-01
                    , 2.092E-01, 1.960E-01, 1.838E-01
                    , 1.763E-01, 1.703E-01, 1.651E-01
                    , 1.562E-01, 1.486E-01, 1.336E-01
                    , 1.224E-01, 1.064E-01, 9.535E-02
                    , 8.707E-02, 8.054E-02, 7.076E-02
                    , 6.362E-02, 5.688E-02, 5.173E-02
                    , 4.422E-02, 3.503E-02, 2.949E-02
                    , 2.577E-02, 2.307E-02, 1.940E-02
                    , 1.703E-02, 1.363E-02, 1.183E-02 };

  Double_t en[kN] = { 1.00000E-03, 1.50000E-03, 2.00000E-03
                    , 3.00000E-03, 4.00000E-03, 5.00000E-03
                    , 6.00000E-03, 8.00000E-03, 1.00000E-02
                    , 1.50000E-02, 2.00000E-02, 3.00000E-02
                    , 4.00000E-02, 5.00000E-02, 6.00000E-02
                    , 8.00000E-02, 1.00000E-01, 1.50000E-01
                    , 2.00000E-01, 3.00000E-01, 4.00000E-01
                    , 5.00000E-01, 6.00000E-01, 8.00000E-01
                    , 1.00000E+00, 1.25000E+00, 1.50000E+00
                    , 2.00000E+00, 3.00000E+00, 4.00000E+00
                    , 5.00000E+00, 6.00000E+00, 8.00000E+00
                    , 1.00000E+01, 1.50000E+01, 2.00000E+01 };

  return Interpolate(energyMeV,en,mu,kN);

}

//_____________________________________________________________________________
Double_t AliTRDsim::Interpolate(Double_t energyMeV
                              , Double_t *en, Double_t *mu, Int_t n)
{
  //
  // Interpolates the photon absorbtion cross section 
  // for a given energy <energyMeV>.
  //

  Double_t de    = 0;
  Int_t    index = 0;
  Int_t    istat = Locate(en,n,energyMeV,index,de);
  if (istat == 0) {
    return (mu[index] - de * (mu[index]   - mu[index+1]) 
                           / (en[index+1] - en[index]  ));
  }
  else {
    return 0.0; 
  }

}

//_____________________________________________________________________________
Int_t AliTRDsim::Locate(Double_t *xv, Int_t n, Double_t xval
                      , Int_t &kl, Double_t &dx) 
{
  //
  // Locates a point (xval) in a 1-dim grid (xv(n))
  //

  if (xval >= xv[n-1]) return  1;
  if (xval <  xv[0])   return -1;

  Int_t km;
  Int_t kh = n - 1;

  kl = 0;
  while (kh - kl > 1) {
    if (xval < xv[km = (kl+kh)/2]) kh = km; 
    else                           kl = km;
  }
  if (xval < xv[kl] || xval > xv[kl+1] || kl >= n-1) {
    printf("Locate failed xv[%d] %f xval %f xv[%d] %f!!!\n"
          ,kl,xv[kl],xval,kl+1,xv[kl+1]);
    exit(1);
  }

  dx = xval - xv[kl];

  return 0;

}

//_____________________________________________________________________________
void AliTRDsim::Streamer(TBuffer &R__b)
{
  //
  // Stream an object of class AliTRDsim.
  //

  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); if (R__v) { }
    TObject::Streamer(R__b);
    R__b >> fNFoils;
    R__b >> fFoilThick;
    R__b >> fGapThick;
    R__b >> fFoilDens;
    R__b >> fGapDens;
    R__b >> fFoilOmega;
    R__b >> fGapOmega;
    R__b >> fFoilZ;
    R__b >> fGapZ;
    R__b >> fFoilA;
    R__b >> fGapA;
    R__b >> fTemp;
    R__b >> fSpNBins;
    R__b >> fSpRange;
    R__b >> fSpBinWidth;
    R__b >> fSpLower;
    R__b >> fSpUpper;
    R__b.ReadArray(fSigma);
    R__b >> fSpectrum;
  } 
  else {
    R__b.WriteVersion(AliTRDsim::IsA());
    TObject::Streamer(R__b);
    R__b << fNFoils;
    R__b << fFoilThick;
    R__b << fGapThick;
    R__b << fFoilDens;
    R__b << fGapDens;
    R__b << fFoilOmega;
    R__b << fGapOmega;
    R__b << fFoilZ;
    R__b << fGapZ;
    R__b << fFoilA;
    R__b << fGapA;
    R__b << fTemp;
    R__b << fSpNBins;
    R__b << fSpRange;
    R__b << fSpBinWidth;
    R__b << fSpLower;
    R__b << fSpUpper;
    R__b.WriteArray(fSigma, fSpNBins);
    R__b << (TObject*) fSpectrum;
  }

}
