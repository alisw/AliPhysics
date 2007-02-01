
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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD simulation - multimodule (regular rad.)                              //
//  after: M. CASTELLANO et al., COMP. PHYS. COMM. 51 (1988) 431             //
//                             + COMP. PHYS. COMM. 61 (1990) 395             //
//                                                                           //
//   17.07.1998 - A.Andronic                                                 //
//   08.12.1998 - simplified version                                         //
//   11.07.2000 - Adapted code to aliroot environment (C.Blume)              //
//   04.06.2004 - Momentum dependent parameters implemented (CBL)            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>

#include <TH1.h>
#include <TRandom.h>
#include <TMath.h>
#include <TParticle.h>
#include <TVirtualMC.h>
#include <TVirtualMCStack.h>

#include "AliModule.h"
#include "AliLog.h"
#include "AliMC.h"

#include "AliTRDsim.h"

ClassImp(AliTRDsim)

//_____________________________________________________________________________
AliTRDsim::AliTRDsim()
  :TObject()
  ,fNFoilsDim(0)
  ,fNFoils(0)
  ,fNFoilsUp(0)
  ,fFoilThick(0)
  ,fGapThick(0)
  ,fFoilDens(0)
  ,fGapDens(0)
  ,fFoilOmega(0)
  ,fGapOmega()
  ,fFoilZ(0)
  ,fGapZ(0)
  ,fFoilA(0)
  ,fGapA(0)
  ,fTemp(0)
  ,fSpNBins(0)
  ,fSpRange(0)
  ,fSpBinWidth(0)
  ,fSpLower(0)
  ,fSpUpper(0)
  ,fSigma(0)
  ,fSpectrum(0)
{
  //
  // AliTRDsim default constructor
  // 

  Init();

}

//_____________________________________________________________________________
AliTRDsim::AliTRDsim(AliModule *mod, Int_t foil, Int_t gap)
  :TObject()
  ,fNFoilsDim(0)
  ,fNFoils(0)
  ,fNFoilsUp(0)
  ,fFoilThick(0)
  ,fGapThick(0)
  ,fFoilDens(0)
  ,fGapDens(0)
  ,fFoilOmega(0)
  ,fGapOmega()
  ,fFoilZ(0)
  ,fGapZ(0)
  ,fFoilA(0)
  ,fGapA(0)
  ,fTemp(0)
  ,fSpNBins(0)
  ,fSpRange(0)
  ,fSpBinWidth(0)
  ,fSpLower(0)
  ,fSpUpper(0)
  ,fSigma(0)
  ,fSpectrum(0)
{
  //
  // AliTRDsim constructor. Takes the material properties of the radiator
  // foils and the gas in the gaps from AliModule <mod>.
  // The default number of foils is 100 with a thickness of 20 mu. The 
  // thickness of the gaps is 500 mu.
  //

  Float_t aFoil;
  Float_t zFoil;
  Float_t rhoFoil;

  Float_t aGap;
  Float_t zGap;
  Float_t rhoGap;

  Float_t rad;
  Float_t abs;

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
  :TObject(s)
  ,fNFoilsDim(s.fNFoilsDim)
  ,fNFoils(0)
  ,fNFoilsUp(0)
  ,fFoilThick(s.fFoilThick)
  ,fGapThick(s.fGapThick)
  ,fFoilDens(s.fFoilDens)
  ,fGapDens(s.fGapDens)
  ,fFoilOmega(s.fFoilOmega)
  ,fGapOmega(s.fGapOmega)
  ,fFoilZ(s.fFoilZ)
  ,fGapZ(s.fGapZ)
  ,fFoilA(s.fFoilA)
  ,fGapA(s.fGapA)
  ,fTemp(s.fTemp)
  ,fSpNBins(s.fSpNBins)
  ,fSpRange(s.fSpRange)
  ,fSpBinWidth(s.fSpBinWidth)
  ,fSpLower(s.fSpLower)
  ,fSpUpper(s.fSpUpper)
  ,fSigma(0)
  ,fSpectrum(0)
{
  //
  // AliTRDsim copy constructor
  //

  if (((AliTRDsim &) s).fNFoils) {
    delete [] ((AliTRDsim &) s).fNFoils;
  }
  ((AliTRDsim &) s).fNFoils   = new Int_t[fNFoilsDim];
  for (Int_t iFoil = 0; iFoil < fNFoilsDim; iFoil++) {
    ((AliTRDsim &) s).fNFoils[iFoil]   = fNFoils[iFoil];
  }  

  if (((AliTRDsim &) s).fNFoilsUp) {
    delete [] ((AliTRDsim &) s).fNFoilsUp;
  }
  ((AliTRDsim &) s).fNFoilsUp = new Double_t[fNFoilsDim];
  for (Int_t iFoil = 0; iFoil < fNFoilsDim; iFoil++) {
    ((AliTRDsim &) s).fNFoilsUp[iFoil] = fNFoilsUp[iFoil];
  }  

  if (((AliTRDsim &) s).fSigma) {
    delete [] ((AliTRDsim &) s).fSigma;
  }
  ((AliTRDsim &) s).fSigma    = new Double_t[fSpNBins];
  for (Int_t iBin = 0; iBin < fSpNBins; iBin++) {
    ((AliTRDsim &) s).fSigma[iBin]     = fSigma[iBin];
  }  

  fSpectrum->Copy(*((AliTRDsim &) s).fSpectrum);

}

//_____________________________________________________________________________
AliTRDsim::~AliTRDsim() 
{
  //
  // AliTRDsim destructor
  //

  if (fSigma) {
    delete [] fSigma;
    fSigma    = 0;
  }

  if (fNFoils) {
    delete [] fNFoils;
    fNFoils   = 0;
  }

  if (fNFoilsUp) {
    delete [] fNFoilsUp;
    fNFoilsUp = 0;
  }

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
void AliTRDsim::Copy(TObject &s) const
{
  //
  // Copy function
  //

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

  if (((AliTRDsim &) s).fNFoils) {
    delete [] ((AliTRDsim &) s).fNFoils;
  }
  ((AliTRDsim &) s).fNFoils   = new Int_t[fNFoilsDim];
  for (Int_t iFoil = 0; iFoil < fNFoilsDim; iFoil++) {
    ((AliTRDsim &) s).fNFoils[iFoil]   = fNFoils[iFoil];
  }  

  if (((AliTRDsim &) s).fNFoilsUp) {
    delete [] ((AliTRDsim &) s).fNFoilsUp;
  }
  ((AliTRDsim &) s).fNFoilsUp = new Double_t[fNFoilsDim];
  for (Int_t iFoil = 0; iFoil < fNFoilsDim; iFoil++) {
    ((AliTRDsim &) s).fNFoilsUp[iFoil] = fNFoilsUp[iFoil];
  }  

  if (((AliTRDsim &) s).fSigma) {
    delete [] ((AliTRDsim &) s).fSigma;
  }
  ((AliTRDsim &) s).fSigma    = new Double_t[fSpNBins];
  for (Int_t iBin = 0; iBin < fSpNBins; iBin++) {
    ((AliTRDsim &) s).fSigma[iBin]     = fSigma[iBin];
  }  

  fSpectrum->Copy(*((AliTRDsim &) s).fSpectrum);

}

//_____________________________________________________________________________
void AliTRDsim::Init()
{
  //
  // Initialization 
  // The default radiator are prolypropilene foils of 10 mu thickness
  // with gaps of 80 mu filled with N2.
  // 

  fNFoilsDim   = 7;

  if (fNFoils) {
    delete [] fNFoils;
  }
  fNFoils      = new Int_t[fNFoilsDim];
  fNFoils[0]   = 170;
  fNFoils[1]   = 225;
  fNFoils[2]   = 275;
  fNFoils[3]   = 305;
  fNFoils[4]   = 325;
  fNFoils[5]   = 340;
  fNFoils[6]   = 350;

  if (fNFoilsUp) {
    delete [] fNFoilsUp;
  }
  fNFoilsUp    = new Double_t[fNFoilsDim];
  fNFoilsUp[0] = 1.25;
  fNFoilsUp[1] = 1.75;
  fNFoilsUp[2] = 2.50;
  fNFoilsUp[3] = 3.50;
  fNFoilsUp[4] = 4.50;
  fNFoilsUp[5] = 5.50;
  fNFoilsUp[6] = 10000.0;

  fFoilThick  = 0.0013;
  fFoilDens   = 0.92;   
  fFoilZ      = 5.28571;
  fFoilA      = 10.4286;
  fFoilOmega  = Omega(fFoilDens,fFoilZ,fFoilA);

  fGapThick   = 0.0060;
  fGapDens    = 0.00125;  
  fGapZ       = 7.0;
  fGapA       = 14.00674;
  fGapOmega   = Omega(fGapDens ,fGapZ ,fGapA );

  fTemp       = 293.16;

  fSpNBins    = 200;
  fSpRange    = 100;
  fSpBinWidth = fSpRange / fSpNBins;
  fSpLower    = 1.0 - 0.5 * fSpBinWidth;
  fSpUpper    = fSpLower + fSpRange;

  if (fSpectrum) delete fSpectrum;
  fSpectrum   = new TH1D("TRspectrum","TR spectrum",fSpNBins,fSpLower,fSpUpper);
  fSpectrum->SetDirectory(0);

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

  // Calculate the TR photons
  return TrPhotons(p, mass, nPhoton, ePhoton);

}

//_____________________________________________________________________________
Int_t AliTRDsim::TrPhotons(Float_t p, Float_t mass
                         , Int_t &nPhoton, Float_t *ePhoton)
{
  //
  // Produces TR photons using a parametric model for regular radiator. Photons
  // with energy larger than 15 keV are included in the MC stack and tracked by VMC
  // machinary.
  //
  // Input parameters:
  // p    - parent momentum [GeV/c]
  // mass - parent mass
  //
  // Output :
  // nPhoton - number of photons which have to be processed by custom code
  // ePhoton - energy of this photons in keV.
  //   

  const Double_t kAlpha  = 0.0072973;
  const Int_t    kSumMax = 30;
	
  Double_t tau   = fGapThick / fFoilThick;

  // Calculate gamma
  Double_t gamma = TMath::Sqrt(p*p + mass*mass) / mass;

  // Select the number of foils corresponding to momentum
  Int_t    foils = SelectNFoils(p);

  fSpectrum->Reset();

  // The TR spectrum
  Double_t csi1;
  Double_t csi2;
  Double_t rho1;
  Double_t rho2;
  Double_t sigma;
  Double_t sum;
  Double_t nEqu;
  Double_t thetaN;
  Double_t aux;
  Double_t energyeV;
  Double_t energykeV;
  for (Int_t iBin = 1; iBin <= fSpNBins; iBin++) {

    energykeV = fSpectrum->GetBinCenter(iBin);
    energyeV  = energykeV * 1.0e3;

    sigma    = Sigma(energykeV);

    csi1      = fFoilOmega / energyeV;
    csi2      = fGapOmega  / energyeV;

    rho1      = 2.5 * energyeV * fFoilThick * 1.0e4 
                               * (1.0 / (gamma*gamma) + csi1*csi1);
    rho2      = 2.5 * energyeV * fFoilThick * 1.0e4 
                               * (1.0 / (gamma*gamma) + csi2 *csi2);

    // Calculate the sum
    sum = 0.0;
    for (Int_t n = 1; n <= kSumMax; n++) {
      thetaN = (TMath::Pi() * 2.0 * n - (rho1 + tau * rho2)) / (1.0 + tau);
      if (thetaN < 0.0) {
        thetaN = 0.0;
      }
      aux   = 1.0 / (rho1 + thetaN) - 1.0 / (rho2 + thetaN);
      sum  += thetaN * (aux*aux) * (1.0 - TMath::Cos(rho1 + thetaN));
    }

    // Equivalent number of foils
    nEqu = (1.0 - TMath::Exp(-foils * sigma)) / (1.0 - TMath::Exp(-sigma));

    // dN / domega
    fSpectrum->SetBinContent(iBin,4.0 * kAlpha * nEqu * sum /  (energykeV * (1.0 + tau)));

  }

  // <nTR> (binsize corr.)
  Float_t nTr     = fSpBinWidth * fSpectrum->Integral();
  // Number of TR photons from Poisson distribution with mean <nTr>
  Int_t   nPhCand = gRandom->Poisson(nTr);
  
  // Link the MC stack and get info about parent electron
  TVirtualMCStack *stack       = gMC->GetStack();
  TParticle       *trGenerator = stack->GetCurrentTrack();
  Int_t    track = stack->GetCurrentTrackNumber();
  Double_t px    = trGenerator->Px() / trGenerator->P();
  Double_t py    = trGenerator->Py() / trGenerator->P();
  Double_t pz    = trGenerator->Pz() / trGenerator->P();

  // Current position of electron
  Double_t x;
  Double_t y;
  Double_t z; 
  gMC->TrackPosition(x,y,z);
  
  // Counter for TR analysed in custom code (e < 15keV)
  nPhoton = 0;  

  for (Int_t iPhoton = 0; iPhoton < nPhCand; iPhoton++) {

    // Energy of the TR photon
    Double_t e = fSpectrum->GetRandom();

    // Put TR photon on particle stack
    if (e > 15.0 ) { 

      e *= 1.0e-6; // Convert it to GeV

      Int_t phtrack;
      stack-> PushTrack(1                 // Must be 1
		       ,track             // Identifier of the parent track, -1 for a primary
		       ,22                // Particle code.
		       ,px*e              // 4 momentum (The photon is generated on the same  
                       ,py*e              // direction as the parent. For irregular radiator one
                       ,pz*e              // can calculate also the angle but this is a secondary
                       ,e                 // order effect)
    	               ,x,y,z,0.0         // 4 vertex	
		       ,0.0,0.0,0.0       // Polarisation
		       ,kPFeedBackPhoton  // Production mechanism (there is no TR in G3 so one
                                          // has to make some convention)
		       ,phtrack           // On output the number of the track stored
		       ,1.0
                       ,1);

    }
    // Custom treatment of TR photons
    else {
  
      ePhoton[nPhoton++] = e;

    }

  }

  return 1;

}

//_____________________________________________________________________________
void AliTRDsim::SetSigma() 
{
  //
  // Sets the absorbtion crosssection for the energies of the TR spectrum
  //

  if (fSigma) {
    delete [] fSigma;
  }
  fSigma = new Double_t[fSpNBins];

  for (Int_t iBin = 0; iBin < fSpNBins; iBin++) {
    Double_t energykeV = iBin * fSpBinWidth + 1.0;
    fSigma[iBin]       = Sigma(energykeV);
  }

}

//_____________________________________________________________________________
Double_t AliTRDsim::Sigma(Double_t energykeV)
{
  //
  // Calculates the absorbtion crosssection for a one-foil-one-gap-radiator
  //

  // keV -> MeV
  Double_t energyMeV = energykeV * 0.001;
  if (energyMeV >= 0.001) {
    return(GetMuPo(energyMeV) * fFoilDens * fFoilThick +
           GetMuAi(energyMeV) * fGapDens  * fGapThick  * GetTemp());
  }
  else {
    return 1.0e6;
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
Double_t AliTRDsim::GetMuAi(Double_t energyMeV)
{
  //
  // Returns the photon absorbtion cross section for air
  // Implemented by Oliver Busch
  //

  const Int_t kN = 38;

  Double_t mu[kN] = { 0.35854E+04, 0.11841E+04, 0.52458E+03,
                      0.16143E+03, 0.14250E+03, 0.15722E+03,
                      0.77538E+02, 0.40099E+02, 0.23313E+02,
                      0.98816E+01, 0.51000E+01, 0.16079E+01,
                      0.77536E+00, 0.35282E+00, 0.24790E+00,
                      0.20750E+00, 0.18703E+00, 0.16589E+00,
                      0.15375E+00, 0.13530E+00, 0.12311E+00,
                      0.10654E+00, 0.95297E-01, 0.86939E-01,
                      0.80390E-01, 0.70596E-01, 0.63452E-01,
                      0.56754E-01, 0.51644E-01, 0.44382E-01,
                      0.35733E-01, 0.30721E-01, 0.27450E-01,
                      0.25171E-01, 0.22205E-01, 0.20399E-01,
                      0.18053E-01, 0.18057E-01 };



  Double_t en[kN] = { 0.10000E-02, 0.15000E-02, 0.20000E-02,
                      0.30000E-02, 0.32029E-02, 0.32029E-02,
                      0.40000E-02, 0.50000E-02, 0.60000E-02,
                      0.80000E-02, 0.10000E-01, 0.15000E-01,
                      0.20000E-01, 0.30000E-01, 0.40000E-01,
                      0.50000E-01, 0.60000E-01, 0.80000E-01,
                      0.10000E+00, 0.15000E+00, 0.20000E+00,
                      0.30000E+00, 0.40000E+00, 0.50000E+00,
                      0.60000E+00, 0.80000E+00, 0.10000E+01,
                      0.12500E+01, 0.15000E+01, 0.20000E+01,
                      0.30000E+01, 0.40000E+01, 0.50000E+01,
                      0.60000E+01, 0.80000E+01, 0.10000E+02,
                      0.15000E+02, 0.20000E+02 };

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

  if (xval >= xv[n-1]) {
    return  1;
  }
  if (xval <  xv[0]) {
    return -1;
  }

  Int_t km;
  Int_t kh = n - 1;

  kl = 0;
  while (kh - kl > 1) {
    if (xval < xv[km = (kl+kh)/2]) {
      kh = km; 
    }
    else {
      kl = km;
    }
  }
  if ((xval <  xv[kl])   || 
      (xval >  xv[kl+1]) || 
      (kl   >= n-1)) {
    AliFatal(Form("Locate failed xv[%d] %f xval %f xv[%d] %f!!!\n"
                 ,kl,xv[kl],xval,kl+1,xv[kl+1]));
    exit(1);
  }

  dx = xval - xv[kl];

  return 0;

}

//_____________________________________________________________________________
Int_t AliTRDsim::SelectNFoils(Float_t p)
{
  //
  // Selects the number of foils corresponding to the momentum
  //

  Int_t foils = fNFoils[fNFoilsDim-1];

  for (Int_t iFoil = 0; iFoil < fNFoilsDim; iFoil++) {
    if (p < fNFoilsUp[iFoil]) {
      foils = fNFoils[iFoil];
      break;
    }
  }

  return foils;

}
