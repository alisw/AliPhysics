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


#include "AliDecayerExodus.h"
#include <Riostream.h>
#include <TMath.h>
#include <AliLog.h>
#include <TH1.h>
#include <TRandom.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TPDGCode.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TF1.h>


ClassImp(AliDecayerExodus)

//---------------------------------------------------------------------------------------------------
//                                 
// Generate electron-pair mass distributions for Dalitz decays according
// to the Kroll-Wada parametrization: N. Kroll, W. Wada: Phys. Rev 98(1955)1355
// and generate electron-pair mass distributions for resonances according
// to the Gounaris-Sakurai parametrization: G.J. Gounaris, J.J. Sakurai: Phys.Rev.Lett. 21(1968)244 
//
// For the electromagnetic form factor the parameterization from
// Lepton-G is used: L.G. Landsberg et al.: Phys. Rep. 128(1985)301
//
// Ralf Averbeck (R.Averbeck@gsi.de) 
// Irem Erdemir  (irem.erdemir@cern.ch)
//
//---------------------------------------------------------------------------------------------------


AliDecayerExodus::AliDecayerExodus():
    AliDecayer(),
    fEPMassPion(0),
    fEPMassEta(0),
    fEPMassEtaPrime(0),
    fEPMassRho(0),
    fEPMassOmega(0),
    fEPMassOmegaDalitz(0),
    fEPMassPhi(0),
    fEPMassPhiDalitz(0),
    fEPMassJPsi(0),
    fPol(new TF1("dsigdcostheta","1.+[0]*x*x",-1.,1.)), /* Polarization Function for resonances */
    fInit(0)

{
// Constructor
}


void AliDecayerExodus::Init()
{
 
// Initialisation
//
   Int_t ibin, nbins;
   Double_t min, maxpion, maxeta, maxomega, maxetaprime, maxphi, binwidth_pion, binwidth_eta, binwidth_omega, binwidth_etaprime, binwidth_phi;
   Double_t pionmass, etamass, omegamass, etaprimemass, phimass, emass, proton_mass, omasspion, omasseta, omassgamma;
   Double_t epsilon_pion, epsilon_eta, epsilon_omega, epsilon_etaprime, epsilon_phi;
   Double_t delta_pion, delta_eta, delta_omega, delta_etaprime, delta_phi;
   Double_t mLL_pion, mLL_eta, mLL_omega, mLL_etaprime, mLL_phi;
   Double_t q_pion, q_eta, q_omega, q_etaprime, q_phi;
   Double_t kwHelp_pion, kwHelp_eta, kwHelp_omega, kwHelp_etaprime, kwHelp_phi;
   Double_t krollWada_pion, krollWada_eta, krollWada_omega, krollWada_etaprime, krollWada_phi;
   Double_t formFactor_pion, formFactor_eta, formFactor_omega, formFactor_etaprime, formFactor_phi;
   Double_t weight_pion, weight_eta, weight_omega_dalitz, weight_etaprime, weight_phi_dalitz;

   Float_t  binwidth;
   Float_t  mass_bin, mass_min, mass_max;
   Double_t vmass_rho, vmass_omega, vmass_phi, vmass_jpsi, vwidth_rho, vwidth_omega, vwidth_phi, vwidth_jpsi;
   Double_t weight_rho, weight_omega, weight_phi, weight_jpsi;

//================================================================================//
//          Create electron pair mass histograms from dalitz decays               //
//================================================================================//

    // Get the particle masses
    // parent
    nbins = 1000;
    mass_min = 0.;
    mass_max = 0.;
    pionmass     = (TDatabasePDG::Instance()->GetParticle(111))->Mass();
    etamass      = (TDatabasePDG::Instance()->GetParticle(221))->Mass();  
    omegamass    = (TDatabasePDG::Instance()->GetParticle(223))->Mass();  
    etaprimemass = (TDatabasePDG::Instance()->GetParticle(331))->Mass();  
    phimass      = (TDatabasePDG::Instance()->GetParticle(333))->Mass();
    // child - electron
    emass = (TDatabasePDG::Instance()->GetParticle(11))->Mass();
    // child - other : third childs from Dalitz decays   
    omasspion  = pionmass;
    omasseta   = etamass;
    omassgamma = 0.;
       
    min         = 2.0 * emass;
    maxpion     = pionmass - omassgamma;
    maxeta      = etamass - omassgamma;
    maxomega    = omegamass - pionmass;
    maxetaprime = etaprimemass - omassgamma;
    maxphi      = phimass - omasseta; 

    binwidth_pion     = (maxpion - min) / (Double_t)nbins;
    binwidth_eta      = (maxeta - min) / (Double_t)nbins;
    binwidth_omega    = (maxomega - min) / (Double_t)nbins;
    binwidth_etaprime = (maxetaprime - min) / (Double_t)nbins;
    binwidth_phi      = (maxphi - min) / (Double_t)nbins;


    epsilon_pion     = (emass / pionmass) * (emass / pionmass);
    epsilon_eta      = (emass / etamass) * (emass / etamass);
    epsilon_omega    = (emass / omegamass) * (emass / omegamass);
    epsilon_etaprime = (emass / etaprimemass) * (emass / etaprimemass);
    epsilon_phi      = (emass / phimass) * (emass / phimass);    


    delta_pion       = (omassgamma / pionmass) * (omassgamma / pionmass);
    delta_eta        = (omassgamma / etamass) * (omassgamma / etamass);
    delta_omega      = (omasspion / omegamass) * (omasspion / omegamass);
    delta_etaprime   = (omassgamma / etaprimemass) * (omassgamma / etaprimemass);
    delta_phi        = (omasseta / phimass) * (omasseta / phimass);    

    // create pair mass histograms for Dalitz decays of pi0, eta, omega, eta' and phi
    if (!fEPMassPion)          {delete fEPMassPion;        fEPMassPion          = new TH1F("fEPMassPion", "Dalitz electron pair from pion", nbins, min, maxpion); }
    if (!fEPMassEta)           {delete fEPMassEta;         fEPMassEta           = new TH1F("fEPMassEta", "Dalitz electron pair from eta", nbins, min, maxeta);}
    if (!fEPMassOmegaDalitz)   {delete fEPMassOmegaDalitz; fEPMassOmegaDalitz   = new TH1F("fEPMassOmegaDalitz", "Dalitz electron pair from omega ", nbins, min, maxomega);}
    if (!fEPMassEtaPrime)      {delete fEPMassEtaPrime;    fEPMassEtaPrime      = new TH1F("fEPMassEtaPrime", "Dalitz electron pair from etaprime", nbins, min, maxetaprime);}
    if (!fEPMassPhiDalitz)     {delete fEPMassPhiDalitz;   fEPMassPhiDalitz     = new TH1F("fEPMassPhiDalitz", "Dalitz electron pair from phi ", nbins, min, maxphi);}


    mLL_pion =  mLL_eta = mLL_omega = mLL_etaprime = mLL_phi = 0.;

    for (ibin = 1; ibin <= nbins; ibin++ )
        {
         mLL_pion     = min + (Double_t)(ibin - 1) * binwidth_pion + binwidth_pion / 2.0;
         mLL_eta      = min + (Double_t)(ibin - 1) * binwidth_eta + binwidth_eta / 2.0;
         mLL_omega    = min + (Double_t)(ibin - 1) * binwidth_omega + binwidth_omega / 2.0;
         mLL_etaprime = min + (Double_t)(ibin - 1) * binwidth_etaprime + binwidth_etaprime / 2.0;
         mLL_phi      = min + (Double_t)(ibin - 1) * binwidth_phi + binwidth_phi / 2.0;


         q_pion        = (mLL_pion / pionmass) * (mLL_pion / pionmass);
         q_eta         = (mLL_eta / etamass) * (mLL_eta / etamass);
         q_omega       = (mLL_omega / omegamass)*(mLL_omega / omegamass);
         q_etaprime    = (mLL_etaprime / etaprimemass) * (mLL_etaprime / etaprimemass);
         q_phi         = (mLL_phi / phimass) * (mLL_phi / phimass);

    if ( q_pion <= 4.0 * epsilon_pion || q_eta <= 4.0 * epsilon_eta || q_omega <= 4.0 * epsilon_omega || q_etaprime <= 4.0 * epsilon_etaprime || q_phi <= 4.0 * epsilon_phi )
       {
         AliFatal("Error in calculating Dalitz mass histogram binning!");
       }
  

    kwHelp_pion     = (1.0 + q_pion /  (1.0 - delta_pion)) * (1.0 + q_pion / (1.0 - delta_pion))
                                     - 4.0 * q_pion / ((1.0 - delta_pion) * (1.0 - delta_pion));

    kwHelp_eta      = (1.0 + q_eta /  (1.0 - delta_eta)) * (1.0 + q_eta / (1.0 - delta_eta))
                                    - 4.0 * q_eta / ((1.0 - delta_eta) * (1.0 - delta_eta));

    kwHelp_omega    = (1.0 + q_omega /  (1.0 - delta_omega)) * (1.0 + q_omega / (1.0 - delta_omega))
                                      - 4.0 * q_omega / ((1.0 - delta_omega) * (1.0 - delta_omega));

    kwHelp_etaprime = (1.0 + q_etaprime /  (1.0 - delta_etaprime)) * (1.0 + q_etaprime / (1.0 - delta_etaprime))
                                         - 4.0 * q_etaprime / ((1.0 - delta_etaprime) * (1.0 - delta_etaprime));

    kwHelp_phi      = (1.0 + q_phi /  (1.0 - delta_phi)) * (1.0 + q_phi / (1.0 - delta_phi))
                                    - 4.0 * q_phi / ((1.0 - delta_phi) * (1.0 - delta_phi));




    if ( kwHelp_pion <= 0.0 || kwHelp_eta <= 0.0 || kwHelp_omega <= 0.0 || kwHelp_etaprime <= 0.0 || kwHelp_phi <= 0.0 )
       {
         AliFatal("Error in calculating Dalitz mass histogram binning!");
    
       }


 // Invariant mass distributions of electron pairs from Dalitz decays
 // using Kroll-Wada function   

      krollWada_pion     = (2.0 / mLL_pion) * TMath::Exp(1.5 * TMath::Log(kwHelp_pion))
                                            * TMath::Sqrt(1.0 - 4.0 * epsilon_pion / q_pion)
                                            * (1.0 + 2.0 * epsilon_pion / q_pion);
    
     
      krollWada_eta      = (2.0 / mLL_eta) * TMath::Exp(1.5 * TMath::Log(kwHelp_eta))
                                           * TMath::Sqrt(1.0 - 4.0 * epsilon_eta / q_eta)
                                           * (1.0 + 2.0 * epsilon_eta / q_eta);
   
   
      krollWada_omega    = (2.0 / mLL_omega) * TMath::Exp(1.5 * TMath::Log(kwHelp_omega))
                                             * TMath::Sqrt(1.0 - 4.0 * epsilon_omega / q_omega)
                                             * (1.0 + 2.0 * epsilon_omega / q_omega);
   
   
      krollWada_etaprime = (2.0 / mLL_etaprime) * TMath::Exp(1.5 * TMath::Log(kwHelp_etaprime))
                                                * TMath::Sqrt(1.0 - 4.0 * epsilon_etaprime / q_etaprime)
                                                * (1.0 + 2.0 * epsilon_etaprime / q_etaprime);
   
      krollWada_phi      = (2.0 / mLL_phi) * TMath::Exp(1.5 * TMath::Log(kwHelp_phi))
                                           * TMath::Sqrt(1.0 - 4.0 * epsilon_phi / q_phi)
                                           * (1.0 + 2.0 * epsilon_phi / q_phi);   



    // Form factors from Lepton-G  
    formFactor_pion     = 1.0/(1.0-5.5*mLL_pion*mLL_pion);
    formFactor_eta      = 1.0/(1.0-1.9*mLL_eta*mLL_eta);
    formFactor_omega    = (TMath::Power(TMath::Power(0.6519,2),2))
                          /(TMath::Power(TMath::Power(0.6519,2)-TMath::Power(mLL_omega, 2), 2)
                          + TMath::Power(0.04198, 2)*TMath::Power(0.6519, 2));
    formFactor_etaprime = (TMath::Power(TMath::Power(0.764,2),2))
                          /(TMath::Power(TMath::Power(0.764,2)-TMath::Power(mLL_etaprime, 2), 2)
                          + TMath::Power(0.1020, 2)*TMath::Power(0.764, 2));
    formFactor_phi      = 1.0; 




    weight_pion         = krollWada_pion * formFactor_pion * formFactor_pion;
    weight_eta          = krollWada_eta * formFactor_eta * formFactor_eta;
    weight_omega_dalitz = krollWada_omega * formFactor_omega;
    weight_etaprime     = krollWada_etaprime * formFactor_etaprime;
    weight_phi_dalitz   = krollWada_phi * formFactor_phi * formFactor_phi;


    // Fill histograms of electron pair masses from dalitz decays
    fEPMassPion       ->AddBinContent(ibin, weight_pion);
    fEPMassEta        ->AddBinContent(ibin, weight_eta);
    fEPMassOmegaDalitz->AddBinContent(ibin, weight_omega_dalitz);
    fEPMassEtaPrime   ->AddBinContent(ibin, weight_etaprime);
    fEPMassPhiDalitz  ->AddBinContent(ibin, weight_phi_dalitz);
    }


   

//===================================================================================//
//         Create electron pair mass histograms from resonance decays                //
//===================================================================================//

   Double_t pimass = 0.13956995;

   // Get the particle masses
   // parent
   vmass_rho   = (TDatabasePDG::Instance()->GetParticle(113))->Mass();  
   vmass_omega = (TDatabasePDG::Instance()->GetParticle(223))->Mass();  
   vmass_phi   = (TDatabasePDG::Instance()->GetParticle(333))->Mass();  
   vmass_jpsi  = (TDatabasePDG::Instance()->GetParticle(443))->Mass();
   // Get the particle widths
   // parent
   vwidth_rho   = (TDatabasePDG::Instance()->GetParticle(113))->Width();   
   vwidth_omega = (TDatabasePDG::Instance()->GetParticle(223))->Width();  
   vwidth_phi   = (TDatabasePDG::Instance()->GetParticle(333))->Width();  
   vwidth_jpsi  = (TDatabasePDG::Instance()->GetParticle(443))->Width();


       if ( mass_min == 0. && mass_max == 0. )
          {
           mass_min  = 2.*pimass;
           mass_max  = 5.;
          }

     binwidth  = (mass_max-mass_min)/(Double_t)nbins;

     // create pair mass histograms for resonances of rho, omega, phi and jpsi
     if (!fEPMassRho)   {delete fEPMassRho;   fEPMassRho    = new TH1F("fEPMassRho","mass rho",nbins,mass_min,mass_max);}
     if (!fEPMassOmega) {delete fEPMassOmega; fEPMassOmega  = new TH1F("fEPMassOmega","mass omega",nbins,mass_min,mass_max);}
     if (!fEPMassPhi)   {delete fEPMassPhi;   fEPMassPhi    = new TH1F("fEPMassPhi","mass phi",nbins,mass_min,mass_max);}
     if (!fEPMassJPsi)  {delete fEPMassJPsi;  fEPMassJPsi   = new TH1F("fEPMassJPsi","mass jpsi",nbins,mass_min,mass_max);}

     for (ibin=1; ibin<=nbins; ibin++ )
     {
     mass_bin = mass_min+(Double_t)(ibin-1)*binwidth+binwidth/2.0;

     weight_rho     = (Float_t)GounarisSakurai(mass_bin,vmass_rho,vwidth_rho,emass);
     weight_omega   = (Float_t)GounarisSakurai(mass_bin,vmass_omega,vwidth_omega,emass);
     weight_phi     = (Float_t)GounarisSakurai(mass_bin,vmass_phi,vwidth_phi,emass); 
     weight_jpsi    = (Float_t)Lorentz(mass_bin,vmass_jpsi,vwidth_jpsi);

     // Fill histograms of electron pair masses from resonance decays
     fEPMassRho  ->AddBinContent(ibin,weight_rho);
     fEPMassOmega->AddBinContent(ibin,weight_omega);
     fEPMassPhi  ->AddBinContent(ibin,weight_phi);
     fEPMassJPsi ->AddBinContent(ibin,weight_jpsi);
    }  

}

Double_t AliDecayerExodus::GounarisSakurai(Float_t mass, Double_t vmass, Double_t vwidth, Double_t emass)
{
// Invariant mass distributions of electron pairs from resonance decays
// of rho, omega and phi
// using Gounaris-Sakurai function

  Double_t corr = 0.;
  Double_t epsilon = 0.;
  Double_t weight = 0.;

  Double_t pimass = 0.13956995;
 
  if(mass>pimass){
  corr = vwidth*(vmass/mass)*exp(1.5*log((mass*mass/4.0-pimass*pimass)
         /(vmass*vmass/4.0-pimass*pimass)));
  }

  epsilon = (emass/mass)*(emass/mass);
       
  if ( 1.0-4.0*epsilon>=0.0 )
  {
   weight = sqrt(1.0-4.0*epsilon)*(1.0+2.0*epsilon)/
                ((vmass*vmass-mass*mass)*(vmass*vmass-mass*mass)+
                (vmass*corr)*(vmass*corr));
  }
  return weight;  
}


Double_t AliDecayerExodus::Lorentz(Float_t mass, Double_t vmass, Double_t vwidth)
{
// Invariant mass distributions of electron pairs from resonance decay
// of jpsi (and it can also be used for other particles except rho, omega and phi) 
// using Lorentz function

  Double_t weight;
  
  weight = (vwidth*vwidth/4.0)/(vwidth*vwidth/4.0+(vmass-mass)*(vmass-mass));

  return weight;

}

void AliDecayerExodus::Decay(Int_t idpart, TLorentzVector* pparent)
{
 
    if (!fInit) {
        Init();
        fInit=1;  
    }

   //local variables for dalitz/2-body decay:
   Double_t pmass, epmass, realp_mass, e1, p1, e3, p3;
   Double_t wp_res, mp_res, md_res, epmass_res, Ed_res, pd_res;
   Double_t PolPar;
   TLorentzVector fProducts_res[2], fProducts_dalitz[3];
   Int_t idRho=113;
   Int_t idOmega=223;
   Int_t idPhi=333;
   Int_t idJPsi=443;
   Int_t idPi0=111;
   Int_t idEta=221;
   Int_t idEtaPrime=331;

   // Get the particle masses of daughters
   Double_t emass, proton_mass, omass_pion, omass_eta, omass_gamma;
   emass       = (TDatabasePDG::Instance()->GetParticle(11)) ->Mass();  
   proton_mass = (TDatabasePDG::Instance()->GetParticle(2212)) ->Mass();  
   omass_pion  = (TDatabasePDG::Instance()->GetParticle(111))->Mass();
   omass_eta   = (TDatabasePDG::Instance()->GetParticle(221))->Mass();  
   omass_gamma = (TDatabasePDG::Instance()->GetParticle(22)) ->Mass();   

   //flat angular distributions
   Double_t costheta, sintheta, cosphi, sinphi, phi;
   Double_t beta_square, lambda;
   costheta = (2.0 * gRandom->Rndm()) - 1.;
   sintheta = TMath::Sqrt((1. + costheta) * (1. - costheta));
   phi      = 2.0 * TMath::ACos(-1.) * gRandom->Rndm();
   sinphi   = TMath::Sin(phi);
   cosphi   = TMath::Cos(phi); 


//-----------------------------------------------------------------------------//
//             Generate Dalitz decays: Pi0/Eta/Omega/EtaPrime/Phi              //
//-----------------------------------------------------------------------------//

  if(idpart==idPi0||idpart==idEta||idpart==idOmega||idpart==idEtaPrime||idpart==idPhi){

   //get the parent mass
   pmass = pparent->M();

   // Sample the electron pair mass from a histogram
   for(;;){
        if(idpart==idPi0){
         epmass = fEPMassPion->GetRandom();
         realp_mass=omass_gamma;
        }else if(idpart==idEta){
         epmass = fEPMassEta->GetRandom();
         realp_mass=omass_gamma;
        }else if(idpart==idOmega){
         epmass = fEPMassOmegaDalitz->GetRandom();
         realp_mass=omass_pion;
        }else if(idpart==idEtaPrime){
         epmass = fEPMassEtaPrime->GetRandom();
         realp_mass=omass_gamma;
        }else if(idpart==idPhi){
         epmass = fEPMassPhiDalitz->GetRandom();
         realp_mass=omass_eta;
        }else{ printf(" Exodus ERROR: Dalitz mass parametrization not found \n");
               return;
        }
        if(pmass-realp_mass>epmass && epmass/2.>emass) break;
   }

   // electron pair kinematics in virtual photon rest frame
   e1 = epmass / 2.;
   p1 = TMath::Sqrt((e1 + emass) * (e1 - emass));


   //Polarization parameters (lambda) for Dalitz:
   if ( realp_mass<0.01 ){
    beta_square = 1.0 - 4.0*(emass*emass)/(epmass*epmass);
    lambda      = beta_square/(2.0-beta_square);
    do{
     costheta = (2.0*gRandom->Rndm())-1.;
    }
    while ( (1.0+lambda*costheta*costheta)<(2.0*gRandom->Rndm()) );
    sintheta = TMath::Sqrt((1. + costheta) * (1. - costheta));
    phi      = 2.0 * TMath::ACos(-1.) * gRandom->Rndm();
    sinphi   = TMath::Sin(phi);
    cosphi   = TMath::Cos(phi); 
   }

   // momentum vectors of electrons in virtual photon rest frame
   Double_t pProd1[3] = {p1 * sintheta * cosphi,
                         p1 * sintheta * sinphi,
                         p1 * costheta};
   Double_t pProd2[3] = {-1.0 * p1 * sintheta * cosphi,
                         -1.0 * p1 * sintheta * sinphi,
                         -1.0 * p1 * costheta};
   fProducts_dalitz[0].SetPx(pProd1[0]);
   fProducts_dalitz[0].SetPy(pProd1[1]);
   fProducts_dalitz[0].SetPz(pProd1[2]);
   fProducts_dalitz[0].SetE(e1);
   fProducts_dalitz[1].SetPx(pProd2[0]);
   fProducts_dalitz[1].SetPy(pProd2[1]);
   fProducts_dalitz[1].SetPz(pProd2[2]);
   fProducts_dalitz[1].SetE(e1);

   // third child kinematics in parent meson rest frame
   e3 = (pmass*pmass + realp_mass*realp_mass - epmass*epmass)/(2. * pmass);
   p3 = TMath::Sqrt((e3+realp_mass) * (e3-realp_mass));
   
   // third child 4-vector in parent meson rest frame
   costheta = (2.0 * gRandom->Rndm()) - 1.;
   sintheta = TMath::Sqrt((1. + costheta) * (1. - costheta));
   phi      = 2.0 * TMath::ACos(-1.) * gRandom->Rndm();
   sinphi   = TMath::Sin(phi);
   cosphi   = TMath::Cos(phi); 
   fProducts_dalitz[2].SetPx(p3 * sintheta * cosphi);
   fProducts_dalitz[2].SetPy(p3 * sintheta * sinphi);
   fProducts_dalitz[2].SetPz(p3 * costheta);
   fProducts_dalitz[2].SetE(e3);

   // boost the dielectron into the parent meson's rest frame
   Double_t eLPparent = TMath::Sqrt(p3*p3 + epmass*epmass);
   TVector3 boostPair( -1.0 * fProducts_dalitz[2].Px() / eLPparent,
                       -1.0 * fProducts_dalitz[2].Py() / eLPparent,
                       -1.0 * fProducts_dalitz[2].Pz() / eLPparent);
   fProducts_dalitz[0].Boost(boostPair);
   fProducts_dalitz[1].Boost(boostPair);

   // boost all decay products into the lab frame
   TVector3 boostLab(pparent->Px() / pparent->E(),
                     pparent->Py() / pparent->E(),
                     pparent->Pz() / pparent->E());
   fProducts_dalitz[0].Boost(boostLab);
   fProducts_dalitz[1].Boost(boostLab);
   fProducts_dalitz[2].Boost(boostLab);

   if(idpart==idPi0) {
     fProducts_pion[0]=fProducts_dalitz[0];
     fProducts_pion[1]=fProducts_dalitz[1];
     fProducts_pion[2]=fProducts_dalitz[2];
   }else if(idpart==idEta){
     fProducts_eta[0]=fProducts_dalitz[0];
     fProducts_eta[1]=fProducts_dalitz[1];
     fProducts_eta[2]=fProducts_dalitz[2];
   }else if(idpart==idOmega){
     fProducts_omega_dalitz[0]=fProducts_dalitz[0];
     fProducts_omega_dalitz[1]=fProducts_dalitz[1];
     fProducts_omega_dalitz[2]=fProducts_dalitz[2];
   }else if(idpart==idEtaPrime){
     fProducts_etaprime[0]=fProducts_dalitz[0];
     fProducts_etaprime[1]=fProducts_dalitz[1];
     fProducts_etaprime[2]=fProducts_dalitz[2];
   }else if(idpart==idPhi){
     fProducts_phi_dalitz[0]=fProducts_dalitz[0];
     fProducts_phi_dalitz[1]=fProducts_dalitz[1];
     fProducts_phi_dalitz[2]=fProducts_dalitz[2];
   }

  }


//-----------------------------------------------------------------------------//
//             Generate 2-body resonance decays: Rho/Omega/Phi/JPsi            //
//-----------------------------------------------------------------------------//
   
  if(idpart==idRho||idpart==idOmega||idpart==idPhi||idpart==idJPsi){


   //get the parent mass
   mp_res = pparent->M();

   //check daughters mass
   md_res=emass;
   if ( mp_res < 2.*md_res ){
        printf("res into ee Decay kinematically impossible! \n");
        return;
   }

   // Sample the electron pair mass from a histogram and set Polarization
   for( ;; ) {
        if(idpart==idRho){
         epmass_res = fEPMassRho->GetRandom();
	 PolPar=0.;
        }else if(idpart==idOmega){
	 epmass_res = fEPMassOmega->GetRandom();
	 PolPar=0.;
        }else if(idpart==idPhi){
	 epmass_res = fEPMassPhi->GetRandom();
	 PolPar=0.;
        }else if(idpart==idPhi){
	 epmass_res = fEPMassPhi->GetRandom();
	 PolPar=0.;
        }else if(idpart==idJPsi){
	 epmass_res = fEPMassJPsi->GetRandom();
	 PolPar=0.;
        }else{ printf(" Exodus ERROR: Resonance mass G-S parametrization not found \n");
               return;
        }
        if ( mp_res < 2.*epmass_res ) break;
   }

   // electron pair kinematics in virtual photon rest frame
   Ed_res = epmass_res/2.;
   pd_res = TMath::Sqrt((Ed_res+md_res)*(Ed_res-md_res));

   // momentum vectors of electrons in virtual photon rest frame 
   fPol->SetParameter(0,PolPar);
   costheta = fPol->GetRandom();
   sintheta = TMath::Sqrt((1. + costheta)*(1. - costheta));
   fProducts_res[0].SetPx(pd_res * sintheta * cosphi);
   fProducts_res[0].SetPy(pd_res * sintheta * sinphi);
   fProducts_res[0].SetPz(pd_res * costheta);
   fProducts_res[0].SetE(Ed_res);
   fProducts_res[1].SetPx(-1.0 * pd_res * sintheta * cosphi);
   fProducts_res[1].SetPy(-1.0 * pd_res * sintheta * sinphi);
   fProducts_res[1].SetPz(-1.0 * pd_res * costheta);
   fProducts_res[1].SetE(Ed_res);

   // Beam parameters in LAB frame
   TLorentzVector pProj, pTarg; 
   Double_t BeamE=3500.;
   pProj.SetPxPyPzE(0.,0.,-1.*BeamE,TMath::Sqrt(BeamE*BeamE+proton_mass*proton_mass)); // Beam 1
   pTarg.SetPxPyPzE(0.,0.,BeamE,TMath::Sqrt(BeamE*BeamE+proton_mass*proton_mass)); // Beam 2

   //re-build parent with G-S mass
   TLorentzVector pparent_corr;
   pparent_corr.SetPx(pparent->Px());
   pparent_corr.SetPy(pparent->Py());
   pparent_corr.SetPz(pparent->Pz());
   pparent_corr.SetE(sqrt(pow(pparent->P(),2)+pow(epmass_res,2)));

   //Boost Beam from CM to Resonance rest frame
   TVector3 betaResCM;
   betaResCM = (-1./pparent_corr.E()*pparent_corr.Vect());
   pProj.Boost(betaResCM);   
   pTarg.Boost(betaResCM);

   //Define Zaxis in C-S frame and rotate legs to it
   TVector3 zaxisCS;
   zaxisCS=(((pProj.Vect()).Unit())-((pTarg.Vect()).Unit())).Unit();
   fProducts_res[0].RotateUz(zaxisCS);
   fProducts_res[1].RotateUz(zaxisCS);

   // boost decay products into the lab frame 
   TVector3 boostLab_res_corr(pparent_corr.Px() / pparent_corr.E(),
                         pparent_corr.Py() / pparent_corr.E(),
                         pparent_corr.Pz() / pparent_corr.E());
   fProducts_res[0].Boost(boostLab_res_corr);
   fProducts_res[1].Boost(boostLab_res_corr);

   if(idpart==idRho) {
    fProducts_rho[0]=fProducts_res[0];
    fProducts_rho[1]=fProducts_res[1];
   }else if(idpart==idOmega){
    fProducts_omega[0]=fProducts_res[0];
    fProducts_omega[1]=fProducts_res[1];
   }else  if(idpart==idPhi){
    fProducts_phi[0]=fProducts_res[0];
    fProducts_phi[1]=fProducts_res[1];
   }else  if(idpart==idJPsi){
    fProducts_jpsi[0]=fProducts_res[0];
    fProducts_jpsi[1]=fProducts_res[1];
   }

 }
  
 return;
}

void AliDecayerExodus::Rot(Double_t pin[3], Double_t pout[3], Double_t costheta, Double_t sintheta,
                           Double_t cosphi, Double_t sinphi) const
{
// Perform rotation
   pout[0] = pin[0]*costheta*cosphi-pin[1]*sinphi+pin[2]*sintheta*cosphi;
   pout[1] = pin[0]*costheta*sinphi+pin[1]*cosphi+pin[2]*sintheta*sinphi;
   pout[2] = -1.0  * pin[0] * sintheta + pin[2] * costheta;
   return;
}


Int_t AliDecayerExodus::ImportParticles(TClonesArray *particles)
{
//
//   Import particles for Dalitz and resonance decays
//

  TClonesArray &clonesParticles = *particles;

  Int_t i, k;
  Double_t px, py, pz, e;

  Int_t pdgD  [3][3] = { {kElectron, -kElectron, 22},     // pizero, eta, etaprime
                         {kElectron, -kElectron, 111},    // omega dalitz
                         {kElectron, -kElectron, 221} };  // phi dalitz

  Int_t pdgR [2] = {kElectron, -kElectron}; // rho, omega, phi, jpsi



    Int_t parentD[3] = { 0,  0, -1};
    Int_t dauD1  [3] = {-1, -1,  1};
    Int_t dauD2  [3] = {-1, -1,  2};

    Int_t parentR[2] = {  0,  0};
    Int_t dauR1  [2] = { -1, -1};
    Int_t dauR2  [2] = { -1, -1};

    for (Int_t j = 0; j < 9; j++){

    // pizero
    if(j==0){
        for (i = 2; i > -1; i--) {
        px = fProducts_pion[i].Px();
        py = fProducts_pion[i].Py();
        pz = fProducts_pion[i].Pz();
        e  = fProducts_pion[i].E();
      new(clonesParticles[2 - i]) TParticle(pdgD[0][i], 1, parentD[i], -1, dauD1[i], dauD2[i], px, py, pz, e, 0., 0., 0., 0.);
      }
      return (3);
      }

    // rho
    else if(j==1){
        for (k = 1; k > -1; k--) {
        px = fProducts_rho[k].Px();
        py = fProducts_rho[k].Py();
        pz = fProducts_rho[k].Pz();
        e  = fProducts_rho[k].E();
      new(clonesParticles[1 - k]) TParticle(pdgR[k], 1, parentR[k], -1, dauR1[k], dauR2[k], px, py, pz, e, 0., 0., 0., 0.);
      }
      return (2);
      }

    // eta
    else if(j==2){
        for (i = 2; i > -1; i--) {
        px = fProducts_eta[i].Px();
        py = fProducts_eta[i].Py();
        pz = fProducts_eta[i].Pz();
        e  = fProducts_eta[i].E();
      new(clonesParticles[2 - i]) TParticle(pdgD[0][i], 1, parentD[i], -1, dauD1[i], dauD2[i], px, py, pz, e, 0., 0., 0., 0.);
      }
      return (3);
      }

    // omega dalitz
    else if(j==3){
        for (i = 2; i > -1; i--) {
        px = fProducts_omega_dalitz[i].Px();
        py = fProducts_omega_dalitz[i].Py();
        pz = fProducts_omega_dalitz[i].Pz();
        e  = fProducts_omega_dalitz[i].E();
      new(clonesParticles[2 - i]) TParticle(pdgD[1][i], 1, parentD[i], -1, dauD1[i], dauD2[i], px, py, pz, e, 0., 0., 0., 0.);
      }
      return (3);
      }

    // omega direct
    else if(j==4){
         for (k = 1; k > -1; k--) {
         px = fProducts_rho[k].Px();
         py = fProducts_rho[k].Py();
         pz = fProducts_rho[k].Pz();
         e  = fProducts_rho[k].E();
       new(clonesParticles[1 - k]) TParticle(pdgR[k], 1, parentR[k], -1, dauR1[k], dauR2[k], px, py, pz, e, 0., 0., 0., 0.);
       }
       return (2);
       }

    // etaprime
    else if(j==5){
        for (i = 2; i > -1; i--) {
        px = fProducts_etaprime[i].Px();
        py = fProducts_etaprime[i].Py();
        pz = fProducts_etaprime[i].Pz();
        e  = fProducts_etaprime[i].E();
      new(clonesParticles[2 - i]) TParticle(pdgD[0][i], 1, parentD[i], -1, dauD1[i], dauD2[i], px, py, pz, e, 0., 0., 0., 0.);
      }
      return (3);
     }

    // phi dalitz
     else if(j==6){
         for (i = 2; i > -1; i--) {
         px = fProducts_phi_dalitz[i].Px();
         py = fProducts_phi_dalitz[i].Py();
         pz = fProducts_phi_dalitz[i].Pz();
         e  = fProducts_phi_dalitz[i].E();
       new(clonesParticles[2 - i]) TParticle(pdgD[2][i], 1, parentD[i], -1, dauD1[i], dauD2[i], px, py, pz, e, 0., 0., 0., 0.);
       }
       return (3);
      }


    // phi direct
    else if(j==7){
        for (k = 1; k > -1; k--) {
        px = fProducts_phi[k].Px();
        py = fProducts_phi[k].Py();
        pz = fProducts_phi[k].Pz();
        e  = fProducts_phi[k].E();
      new(clonesParticles[1 - k]) TParticle(pdgR[k], 1, parentR[k], -1, dauR1[k], dauR2[k], px, py, pz, e, 0., 0., 0., 0.);
      }
      return (2);
    }

    // jpsi direct
    else if(j==8){
          for (k = 1; k > -1; k--) {
          px = fProducts_jpsi[k].Px();
          py = fProducts_jpsi[k].Py();
          pz = fProducts_jpsi[k].Pz();
          e  = fProducts_jpsi[k].E();
       new(clonesParticles[1 - k]) TParticle(pdgR[k], 1, parentR[k], -1, dauR1[k], dauR2[k], px, py, pz, e, 0., 0., 0., 0.);
       }
       return (2);
    }

   }

   return particles->GetEntries();

}


void AliDecayerExodus::Decay(TClonesArray *array)
{
  // Replace all Dalitz(pi0,eta,omega,eta',phi) and resonance(rho,omega,phi,jpsi) decays with the correct matrix element decays
  // for di-electron cocktail calculations


  Int_t nt = array->GetEntriesFast();
  TParticle* dp3[3];
  TParticle* dp2[2];
  Int_t fd3, ld3, fd2, ld2, fd, ld;
  Int_t j, k;

  for (Int_t i = 0; i < nt; i++) {
  TParticle* part = (TParticle*) (array->At(i));
  if (part->GetPdgCode() != 111 || part->GetPdgCode() != 221 || part->GetPdgCode() != 223 || part->GetPdgCode() != 331 || part->GetPdgCode() != 333 || part->GetPdgCode() != 443 ) continue;

  //
  // Pizero Dalitz
  //
  if(part->GetPdgCode() == 111){
  
  fd3 = part->GetFirstDaughter() - 1;
  ld3 = part->GetLastDaughter()  - 1;
  
  if (fd3 < 0)                           continue;
  if ((ld3 - fd3) != 2)                  continue;
  
  for (j = 0; j < 3; j++) dp3[j] = (TParticle*) (array->At(fd3+j));

  if((dp3[0]->GetPdgCode() != 22) && (TMath::Abs(dp3[1]->GetPdgCode()) != 11))   continue;

  TLorentzVector Pizero(part->Px(), part->Py(), part->Pz(), part->Energy());
  Decay(111, &Pizero);
  for (j = 0; j < 3; j++) dp3[j]->SetMomentum(fProducts_pion[2-j]);
  }


  //
  // Eta Dalitz
  //

  if(part->GetPdgCode() == 221){
      
  fd3 = part->GetFirstDaughter() - 1;
  ld3 = part->GetLastDaughter()  - 1;

  if (fd3 < 0)                           continue;
  if ((ld3 - fd3) != 2)                  continue;
                      
  for (j = 0; j < 3; j++) dp3[j] = (TParticle*) (array->At(fd3+j));

  if((dp3[0]->GetPdgCode() != 22) && ((TMath::Abs(dp3[1]->GetPdgCode()) != 11)))   continue;

  TLorentzVector Eta(part->Px(), part->Py(), part->Pz(), part->Energy());
  Decay(221, &Eta);
  for (j = 0; j < 3; j++) dp3[j]->SetMomentum(fProducts_eta[2-j]);
  }

  //
  // Rho
  //

  if(part->GetPdgCode() == 113){

  fd2 = part->GetFirstDaughter() - 1;
  ld2 = part->GetLastDaughter()  - 1;

  if (fd2 < 0)                           continue;
  if ((ld2 - fd2) != 1)                  continue;

  for (k = 0; k < 2; k++) dp2[k] = (TParticle*) (array->At(fd2+k));

  if((dp2[0]->GetPdgCode() != 11) && ((TMath::Abs(dp2[1]->GetPdgCode()) != 11)))   continue;

  TLorentzVector Rho(part->Px(), part->Py(), part->Pz(), part->Energy());
  Decay(113, &Rho);
  for (k = 0; k < 2; k++) dp2[k]->SetMomentum(fProducts_rho[1-k]);
  }

  //
  // Omega dalitz and direct
  //

  if(part->GetPdgCode() == 223){

  fd = part->GetFirstDaughter() - 1;
  ld = part->GetLastDaughter()  - 1;

  if (fd < 0)               continue;

  if ((ld - fd) == 2){

  for (j = 0; j < 3; j++) dp3[j] = (TParticle*) (array->At(fd+j));
  if( dp3[0]->GetPdgCode() != 111 && (TMath::Abs(dp3[1]->GetPdgCode()) != 11)) continue;

  TLorentzVector Omegadalitz(part->Px(), part->Py(), part->Pz(), part->Energy());
  Decay(223, &Omegadalitz);
  for (j = 0; j < 3; j++) dp3[j]->SetMomentum(fProducts_omega_dalitz[2-j]);
  }

  else if ((ld - fd) == 1) {

  for (k = 0; k < 2; k++) dp2[k] = (TParticle*) (array->At(fd+k));
  if( dp2[0]->GetPdgCode() != 11 && (TMath::Abs(dp2[1]->GetPdgCode()) != 11))   continue;

  TLorentzVector Omega(part->Px(), part->Py(), part->Pz(), part->Energy());
  Decay(223, &Omega);
  for (k = 0; k < 2; k++) dp2[k]->SetMomentum(fProducts_omega[1-k]);
  }
 }
 
  //
  // Etaprime dalitz
  //

  if(part->GetPdgCode() == 331){

  fd3 = part->GetFirstDaughter() - 1;
  ld3 = part->GetLastDaughter()  - 1;

  if (fd3 < 0)                           continue;
  if ((ld3 - fd3) != 2)                  continue;

  for (j = 0; j < 3; j++) dp3[j] = (TParticle*) (array->At(fd3+j));

  if((dp3[0]->GetPdgCode() != 22) && ((TMath::Abs(dp3[1]->GetPdgCode()) != 11)))   continue;

  TLorentzVector Etaprime(part->Px(), part->Py(), part->Pz(), part->Energy());
  Decay(331, &Etaprime);
  for (j = 0; j < 3; j++) dp3[j]->SetMomentum(fProducts_etaprime[2-j]);
  }

  //
  // Phi dalitz and direct
  //
  if(part->GetPdgCode() == 333){

  fd = part->GetFirstDaughter() - 1;
  ld = part->GetLastDaughter()  - 1;

  if (fd < 0)               continue;
  if ((ld - fd) == 2){
  for (j = 0; j < 3; j++) dp3[j] = (TParticle*) (array->At(fd+j));
  if( dp3[0]->GetPdgCode() != 221 && (TMath::Abs(dp3[1]->GetPdgCode()) != 11)) continue;

  TLorentzVector Phidalitz(part->Px(), part->Py(), part->Pz(), part->Energy());
  Decay(333, &Phidalitz);
  for (j = 0; j < 3; j++) dp3[j]->SetMomentum(fProducts_phi_dalitz[2-j]);
  } 

  else if ((ld - fd) == 1) {
  for (k = 0; k < 2; k++) dp2[k] = (TParticle*) (array->At(fd+k));
  if( dp2[0]->GetPdgCode() != 11 && (TMath::Abs(dp2[1]->GetPdgCode()) != 11))   continue;

  TLorentzVector Phi(part->Px(), part->Py(), part->Pz(), part->Energy());
  Decay(333, &Phi);
  for (k = 0; k < 2; k++) dp2[k]->SetMomentum(fProducts_phi[1-k]);
   }
  } 

  //
  // JPsi
  //

  if(part->GetPdgCode() == 443){

  fd2 = part->GetFirstDaughter() - 1;
  ld2 = part->GetLastDaughter()  - 1;

  if (fd2 < 0)                           continue;
  if ((ld2 - fd2) != 1)                  continue;

  for (k = 0; k < 2; k++) dp2[k] = (TParticle*) (array->At(fd2+k));

  if((dp2[0]->GetPdgCode() != 11) && ((TMath::Abs(dp2[1]->GetPdgCode()) != 11)))   continue;

  TLorentzVector JPsi(part->Px(), part->Py(), part->Pz(), part->Energy());
  Decay(443, &JPsi);
  for (k = 0; k < 2; k++) dp2[k]->SetMomentum(fProducts_jpsi[1-k]);
  }

 }
}


AliDecayerExodus& AliDecayerExodus::operator=(const AliDecayerExodus& rhs)
{
  // Assignment operator
  rhs.Copy(*this);
  return *this;
}

void AliDecayerExodus::Copy(TObject&) const
{
  //
  // Copy 
  //
  Fatal("Copy","Not implemented!\n");
}


AliDecayerExodus::AliDecayerExodus(const AliDecayerExodus &decayer)
     : AliDecayer(), 
      fEPMassPion(0),
      fEPMassEta(0),
      fEPMassEtaPrime(0),
      fEPMassRho(0),
      fEPMassOmega(0),
      fEPMassOmegaDalitz(0),
      fEPMassPhi(0),
      fEPMassPhiDalitz(0), 
      fEPMassJPsi(0),
      fInit(0)
{
 // Copy Constructor
    decayer.Copy(*this);
}


   
