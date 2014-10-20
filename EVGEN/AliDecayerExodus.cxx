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
   Double_t pionmass, etamass, omegamass, etaprimemass, phimass, emass, omasspion, omasseta, omassgamma;
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


   Double_t pmass_pion, pmass_eta, pmass_omega_dalitz, pmass_etaprime, pmass_phi_dalitz;
   Double_t emass, omass_pion, omass_eta, omass_gamma, epmass_pion, epmass_eta, epmass_omega_dalitz, epmass_etaprime, epmass_phi_dalitz;
   Double_t e1_pion, e1_eta, e1_omega, e1_etaprime, e1_phi;
   Double_t p1_pion, p1_eta, p1_omega, p1_etaprime, p1_phi;
   Double_t e3_gamma_pion, e3_gamma_eta, e3_pion, e3_gamma_etaprime, e3_eta; 
   Double_t p3_gamma_pion, p3_gamma_eta, p3_pion, p3_gamma_etaprime, p3_eta; 

   Double_t wp_rho, wp_omega, wp_phi, wp_jpsi, epmass_rho, epmass_omega, epmass_phi, epmass_jpsi;
   Double_t mp_rho, mp_omega, mp_phi, mp_jpsi, md_rho, md_omega, md_phi, md_jpsi;
   Double_t Ed_rho, Ed_omega, Ed_phi, Ed_jpsi, pd_rho, pd_omega, pd_phi, pd_jpsi;


    md_rho =  md_omega =  md_phi =  md_jpsi = 0.;


   Double_t costheta, sintheta, cosphi, sinphi, phi;

   // Get the particle masses of daughters
   emass       = (TDatabasePDG::Instance()->GetParticle(11)) ->Mass();  
   omass_pion  = (TDatabasePDG::Instance()->GetParticle(111))->Mass();
   omass_eta   = (TDatabasePDG::Instance()->GetParticle(221))->Mass();  
   omass_gamma = (TDatabasePDG::Instance()->GetParticle(22)) ->Mass();   

   // Get the particle widths of mothers for resonances
   wp_rho   = (TDatabasePDG::Instance()->GetParticle(113))->Width();
   wp_omega = (TDatabasePDG::Instance()->GetParticle(223))->Width();
   wp_phi   = (TDatabasePDG::Instance()->GetParticle(333))->Width();
   wp_jpsi  = (TDatabasePDG::Instance()->GetParticle(443))->Width();

   costheta = (2.0 * gRandom->Rndm()) - 1.;
   sintheta = TMath::Sqrt((1. + costheta) * (1. - costheta));
   phi      = 2.0 * TMath::ACos(-1.) * gRandom->Rndm();
   sinphi   = TMath::Sin(phi);
   cosphi   = TMath::Cos(phi); 


//-----------------------------------------------------------------------------//
//                        Generate Pizero Dalitz decay                         //
//-----------------------------------------------------------------------------//
   
   if(idpart==111){
   pmass_pion = pparent->M();

   for(;;){
   // Sample the electron pair mass from a histogram
   epmass_pion = fEPMassPion->GetRandom();
   if (pmass_pion-omass_gamma>epmass_pion && epmass_pion/2.>emass) break;
   }

   // electron pair kinematics in virtual photon rest frame
   e1_pion = epmass_pion / 2.;
   p1_pion = TMath::Sqrt((e1_pion + emass) * (e1_pion - emass));

   // momentum vectors of electrons in virtual photon rest frame
   Double_t pProd1_pion[3] = {p1_pion * sintheta * cosphi,
                              p1_pion * sintheta * sinphi,
                              p1_pion * costheta};

   Double_t pProd2_pion[3] = {-1.0 * p1_pion * sintheta * cosphi,
                              -1.0 * p1_pion * sintheta * sinphi,
                              -1.0 * p1_pion * costheta};


   // third child kinematics in parent meson rest frame
   e3_gamma_pion       = (pmass_pion * pmass_pion - epmass_pion * epmass_pion)/(2. * pmass_pion);
   p3_gamma_pion       = TMath::Sqrt((e3_gamma_pion  * e3_gamma_pion));


   // third child 4-vector in parent meson rest frame
   fProducts_pion[2].SetPx(p3_gamma_pion * sintheta * cosphi);
   fProducts_pion[2].SetPy(p3_gamma_pion * sintheta * sinphi);
   fProducts_pion[2].SetPz(p3_gamma_pion * costheta);
   fProducts_pion[2].SetE(e3_gamma_pion);


   // electron 4-vectors in properly rotated virtual photon rest frame
   Double_t pRot1_pion[3] = {0.};
   Rot(pProd1_pion, pRot1_pion, costheta, -sintheta, -cosphi, -sinphi);
   Double_t pRot2_pion[3] = {0.};
   Rot(pProd2_pion, pRot2_pion, costheta, -sintheta, -cosphi, -sinphi);
   fProducts_pion[0].SetPx(pRot1_pion[0]);
   fProducts_pion[0].SetPy(pRot1_pion[1]);
   fProducts_pion[0].SetPz(pRot1_pion[2]);
   fProducts_pion[0].SetE(e1_pion);
   fProducts_pion[1].SetPx(pRot2_pion[0]);
   fProducts_pion[1].SetPy(pRot2_pion[1]);
   fProducts_pion[1].SetPz(pRot2_pion[2]);
   fProducts_pion[1].SetE(e1_pion);

   // boost the dielectron into the parent meson's rest frame
   Double_t eLPparent_pion = TMath::Sqrt(p3_gamma_pion * p3_gamma_pion + epmass_pion * epmass_pion);
   TVector3 boostPair_pion( -1.0 * fProducts_pion[2].Px() / eLPparent_pion,
                       -1.0 * fProducts_pion[2].Py() / eLPparent_pion,
                       -1.0 * fProducts_pion[2].Pz() / eLPparent_pion);
   fProducts_pion[0].Boost(boostPair_pion);
   fProducts_pion[1].Boost(boostPair_pion);

   // boost all decay products into the lab frame
   TVector3 boostLab_pion(pparent->Px() / pparent->E(),
                     pparent->Py() / pparent->E(),
                     pparent->Pz() / pparent->E());

   fProducts_pion[0].Boost(boostLab_pion);
   fProducts_pion[1].Boost(boostLab_pion);
   fProducts_pion[2].Boost(boostLab_pion);

   } 


//-----------------------------------------------------------------------------//
//                        Generate Rho resonance decay                         //
//-----------------------------------------------------------------------------//

   else if(idpart==113){
   // calculate rho mass
        if(wp_rho!=0.0){
          mp_rho = pparent->M();
          }
        else{
        Double_t x_rho=pparent->Px(); Double_t y_rho=pparent->Py(); Double_t z_rho=pparent->Pz();
        Double_t t_rho=pparent->E();
        Double_t p_rho=x_rho*x_rho+y_rho*y_rho+z_rho*z_rho;
        Double_t Q2_rho=abs((t_rho*t_rho)-(p_rho*p_rho));
        mp_rho = sqrt(Q2_rho);
        }
   // daughter
       if ( mp_rho < 2.*md_rho )
          {
           printf("Rho into ee Decay kinematically impossible! \n");
           return;
           }
  
   for( ;; ) {
   // Sample the electron pair mass from a histogram  
   epmass_rho = fEPMassRho->GetRandom();
   if ( mp_rho < 2.*epmass_rho ) break;
   }

   // electron pair kinematics in virtual photon rest frame
   Ed_rho = epmass_rho/2.;
   pd_rho = TMath::Sqrt((Ed_rho+md_rho)*(Ed_rho-md_rho));
  
   // momentum vectors of electrons in virtual photon rest frame 
   Double_t pProd1_rho[3] = {pd_rho * sintheta * cosphi,
                             pd_rho * sintheta * sinphi,
                             pd_rho * costheta};
  
   Double_t pProd2_rho[3] = {-1.0 * pd_rho * sintheta * cosphi,
                             -1.0 * pd_rho * sintheta * sinphi,
                             -1.0 * pd_rho * costheta};
                                                                                                                                                                        

   // electron 4 vectors in properly rotated virtual photon rest frame
   Double_t pRot1_rho[3] = {0.};
   Rot(pProd1_rho, pRot1_rho, costheta, -sintheta, -cosphi, -sinphi);
   Double_t pRot2_rho[3] = {0.};
   Rot(pProd2_rho, pRot2_rho, costheta, -sintheta, -cosphi, -sinphi);
   fProducts_rho[0].SetPx(pRot1_rho[0]);
   fProducts_rho[0].SetPy(pRot1_rho[1]);
   fProducts_rho[0].SetPz(pRot1_rho[2]);
   fProducts_rho[0].SetE(Ed_rho);
   fProducts_rho[1].SetPx(pRot2_rho[0]);
   fProducts_rho[1].SetPy(pRot2_rho[1]);
   fProducts_rho[1].SetPz(pRot2_rho[2]);
   fProducts_rho[1].SetE(Ed_rho);
   
   
   // boost decay products into the lab frame 
   TVector3 boostLab_rho(pparent->Px() / pparent->E(),
                         pparent->Py() / pparent->E(),
                         pparent->Pz() / pparent->E());
   
   fProducts_rho[0].Boost(boostLab_rho);
   fProducts_rho[1].Boost(boostLab_rho);
   }


//-----------------------------------------------------------------------------//
//                        Generate Eta Dalitz decay                            //
//-----------------------------------------------------------------------------// 
  
   else if(idpart==221){
   pmass_eta = pparent->M();

   for(;;){
   // Sample the electron pair mass from a histogram
   epmass_eta = fEPMassEta->GetRandom();
   if(pmass_eta-omass_gamma>epmass_eta && epmass_eta/2.>emass)
   break;
   }
   
   // electron pair kinematics in virtual photon rest frame
   e1_eta = epmass_eta / 2.;
   p1_eta = TMath::Sqrt((e1_eta + emass) * (e1_eta - emass));

   // momentum vectors of electrons in virtual photon rest frame
   Double_t pProd1_eta[3] = {p1_eta * sintheta * cosphi,
                             p1_eta * sintheta * sinphi,
                             p1_eta * costheta};
   Double_t pProd2_eta[3] = {-1.0 * p1_eta * sintheta * cosphi,
                             -1.0 * p1_eta * sintheta * sinphi,
                             -1.0 * p1_eta * costheta};

   // third child kinematics in parent meson rest frame
   e3_gamma_eta       = (pmass_eta * pmass_eta - epmass_eta * epmass_eta)/(2. * pmass_eta);
   p3_gamma_eta       = TMath::Sqrt((e3_gamma_eta * e3_gamma_eta));


   // third child 4-vector in parent meson rest frame
   fProducts_eta[2].SetPx(p3_gamma_eta * sintheta * cosphi);
   fProducts_eta[2].SetPy(p3_gamma_eta * sintheta * sinphi);
   fProducts_eta[2].SetPz(p3_gamma_eta * costheta);
   fProducts_eta[2].SetE(e3_gamma_eta); 

   // electron 4-vectors in properly rotated virtual photon rest frame
   Double_t pRot1_eta[3] = {0.};
   Rot(pProd1_eta, pRot1_eta, costheta, -sintheta, -cosphi, -sinphi);
   Double_t pRot2_eta[3] = {0.};
   Rot(pProd2_eta, pRot2_eta, costheta, -sintheta, -cosphi, -sinphi);
   fProducts_eta[0].SetPx(pRot1_eta[0]);
   fProducts_eta[0].SetPy(pRot1_eta[1]);
   fProducts_eta[0].SetPz(pRot1_eta[2]);
   fProducts_eta[0].SetE(e1_eta);
   fProducts_eta[1].SetPx(pRot2_eta[0]);
   fProducts_eta[1].SetPy(pRot2_eta[1]);
   fProducts_eta[1].SetPz(pRot2_eta[2]);
   fProducts_eta[1].SetE(e1_eta);

   // boost the dielectron into the parent meson's rest frame
   Double_t eLPparent_eta = TMath::Sqrt(p3_gamma_eta * p3_gamma_eta + epmass_eta * epmass_eta);
   TVector3 boostPair_eta( -1.0 * fProducts_eta[2].Px() / eLPparent_eta,
                       -1.0 * fProducts_eta[2].Py() / eLPparent_eta,
                       -1.0 * fProducts_eta[2].Pz() / eLPparent_eta);
   fProducts_eta[0].Boost(boostPair_eta);
   fProducts_eta[1].Boost(boostPair_eta);

   // boost all decay products into the lab frame
   TVector3 boostLab_eta(pparent->Px() / pparent->E(),
                     pparent->Py() / pparent->E(),
                     pparent->Pz() / pparent->E());

   fProducts_eta[0].Boost(boostLab_eta);
   fProducts_eta[1].Boost(boostLab_eta);
   fProducts_eta[2].Boost(boostLab_eta);

    }
 
   
//-----------------------------------------------------------------------------//
//                        Generate Omega Dalitz decay                          //
//-----------------------------------------------------------------------------//

   else if(idpart==223){
   pmass_omega_dalitz = pparent->M();
   for(;;){
   // Sample the electron pair mass from a histogram
   epmass_omega_dalitz = fEPMassOmegaDalitz->GetRandom();
   if(pmass_omega_dalitz-omass_pion>epmass_omega_dalitz && epmass_omega_dalitz/2.>emass)
   break;}

   // electron pair kinematics in virtual photon rest frame
   e1_omega = epmass_omega_dalitz / 2.;
   p1_omega = TMath::Sqrt((e1_omega + emass) * (e1_omega - emass)); 

   // momentum vectors of electrons in virtual photon rest frame
   Double_t pProd1_omega_dalitz[3] = {p1_omega * sintheta * cosphi,
                               p1_omega * sintheta * sinphi,
                               p1_omega * costheta};
   Double_t pProd2_omega_dalitz[3] = {-1.0 * p1_omega * sintheta * cosphi,
                               -1.0 * p1_omega * sintheta * sinphi,
                               -1.0 * p1_omega * costheta};

   // third child kinematics in parent meson rest frame
   e3_pion       = (pmass_omega_dalitz * pmass_omega_dalitz + omass_pion * omass_pion - epmass_omega_dalitz * epmass_omega_dalitz)/(2. * pmass_omega_dalitz);
   p3_pion       = TMath::Sqrt((e3_pion + omass_pion)  * (e3_pion - omass_pion));

   // third child 4-vector in parent meson rest frame
   fProducts_omega_dalitz[2].SetPx(p3_pion * sintheta * cosphi);
   fProducts_omega_dalitz[2].SetPy(p3_pion * sintheta * sinphi);
   fProducts_omega_dalitz[2].SetPz(p3_pion * costheta);
   fProducts_omega_dalitz[2].SetE(e3_pion);

   // lepton 4-vectors in properly rotated virtual photon rest frame
   Double_t pRot1_omega_dalitz[3] = {0.};
   Rot(pProd1_omega_dalitz, pRot1_omega_dalitz, costheta, -sintheta, -cosphi, -sinphi);
   Double_t pRot2_omega_dalitz[3] = {0.};
   Rot(pProd2_omega_dalitz, pRot2_omega_dalitz, costheta, -sintheta, -cosphi, -sinphi);
   fProducts_omega_dalitz[0].SetPx(pRot1_omega_dalitz[0]);
   fProducts_omega_dalitz[0].SetPy(pRot1_omega_dalitz[1]);
   fProducts_omega_dalitz[0].SetPz(pRot1_omega_dalitz[2]);
   fProducts_omega_dalitz[0].SetE(e1_omega);
   fProducts_omega_dalitz[1].SetPx(pRot2_omega_dalitz[0]);
   fProducts_omega_dalitz[1].SetPy(pRot2_omega_dalitz[1]);
   fProducts_omega_dalitz[1].SetPz(pRot2_omega_dalitz[2]);
   fProducts_omega_dalitz[1].SetE(e1_omega); 

   // boost the dielectron into the parent meson's rest frame
   Double_t eLPparent_omega = TMath::Sqrt(p3_pion * p3_pion + epmass_omega_dalitz * epmass_omega_dalitz);
   TVector3 boostPair_omega( -1.0 * fProducts_omega_dalitz[2].Px() / eLPparent_omega,
                       -1.0 * fProducts_omega_dalitz[2].Py() / eLPparent_omega,
                       -1.0 * fProducts_omega_dalitz[2].Pz() / eLPparent_omega);
   fProducts_omega_dalitz[0].Boost(boostPair_omega);
   fProducts_omega_dalitz[1].Boost(boostPair_omega);

   // boost all decay products into the lab frame
   TVector3 boostLab_omega_dalitz(pparent->Px() / pparent->E(),
                     pparent->Py() / pparent->E(),
                     pparent->Pz() / pparent->E());

   fProducts_omega_dalitz[0].Boost(boostLab_omega_dalitz);
   fProducts_omega_dalitz[1].Boost(boostLab_omega_dalitz);
   fProducts_omega_dalitz[2].Boost(boostLab_omega_dalitz);
    

//-----------------------------------------------------------------------------//
//                       Generate Omega resonance decay                        //
//-----------------------------------------------------------------------------//

      if(wp_omega!=0.0){
      // calculate omega mass  
         mp_omega = pparent->M();
         }
      else{
           Double_t x_omega=pparent->Px(); Double_t y_omega=pparent->Py(); Double_t z_omega=pparent->Pz();
           Double_t t_omega=pparent->E();
           Double_t p_omega=x_omega*x_omega+y_omega*y_omega+z_omega*z_omega;
           Double_t Q2_omega= abs((t_omega*t_omega)-(p_omega*p_omega));
           mp_omega = sqrt(Q2_omega);
           }

   // daughter
   if ( mp_omega< 2.*md_omega )
      {
       printf("Omega into ee Decay kinematically impossible! \n");
       return;
       }

   for( ;; ) {
   // Sample the electron pair mass from a histogram 
   epmass_omega = fEPMassOmega->GetRandom();
   if( mp_omega < 2.*epmass_omega ) break;
   }

   // electron pair kinematics in virtual photon rest frame
   Ed_omega = epmass_omega/2.;
   pd_omega = TMath::Sqrt((Ed_omega+md_omega)*(Ed_omega-md_omega));

   // momentum vectors of electrons in virtual photon rest frame
   Double_t pProd1_omega[3] = {pd_omega * sintheta * cosphi,
                               pd_omega * sintheta * sinphi,
                               pd_omega * costheta};

   Double_t pProd2_omega[3] = {-1.0 * pd_omega * sintheta * cosphi,
                               -1.0 * pd_omega * sintheta * sinphi,
                               -1.0 * pd_omega * costheta}; 


   // lepton 4 vectors in properly rotated virtual photon rest frame
   Double_t pRot1_omega[3] = {0.};
   Rot(pProd1_omega, pRot1_omega, costheta, -sintheta, -cosphi, -sinphi);
   Double_t pRot2_omega[3] = {0.};
   Rot(pProd2_omega, pRot2_omega, costheta, -sintheta, -cosphi, -sinphi);
   fProducts_omega[0].SetPx(pRot1_omega[0]);
   fProducts_omega[0].SetPy(pRot1_omega[1]);
   fProducts_omega[0].SetPz(pRot1_omega[2]);
   fProducts_omega[0].SetE(Ed_omega);
   fProducts_omega[1].SetPx(pRot2_omega[0]);
   fProducts_omega[1].SetPy(pRot2_omega[1]);
   fProducts_omega[1].SetPz(pRot2_omega[2]);
   fProducts_omega[1].SetE(Ed_omega);

   // boost decay products into the lab frame 
   TVector3 boostLab_omega(pparent->Px() / pparent->E(),
                           pparent->Py() / pparent->E(),
                           pparent->Pz() / pparent->E());

   fProducts_omega[0].Boost(boostLab_omega);
   fProducts_omega[1].Boost(boostLab_omega);

   }

//-----------------------------------------------------------------------------//
//                      Generate Etaprime Dalitz decay                         //
//-----------------------------------------------------------------------------//  

   else if(idpart==331){
   pmass_etaprime = pparent->M();
   for(;;){
   // Sample the electron pair mass from a histogram
   epmass_etaprime = fEPMassEtaPrime->GetRandom();
   if(pmass_etaprime-omass_gamma>epmass_etaprime && epmass_etaprime/2.>emass)
   break;}
  
   // electron pair kinematics in virtual photon rest frame
   e1_etaprime = epmass_etaprime / 2.;
   p1_etaprime = TMath::Sqrt((e1_etaprime + emass) * (e1_etaprime - emass));

   // momentum vectors of electrons in virtual photon rest frame
   Double_t pProd1_etaprime[3] = {p1_etaprime * sintheta * cosphi,
                                  p1_etaprime * sintheta * sinphi,
                                  p1_etaprime * costheta};
   Double_t pProd2_etaprime[3] = {-1.0 * p1_etaprime * sintheta * cosphi,
                                  -1.0 * p1_etaprime * sintheta * sinphi,
                                  -1.0 * p1_etaprime * costheta};

   // third child kinematics in parent meson rest frame
   e3_gamma_etaprime       = (pmass_etaprime * pmass_etaprime + omass_gamma * omass_gamma - epmass_etaprime * epmass_etaprime)/(2. * pmass_etaprime);
   p3_gamma_etaprime       = TMath::Sqrt((e3_gamma_etaprime + omass_gamma)  * (e3_gamma_etaprime - omass_gamma));

   // third child 4-vector in parent meson rest frame
   fProducts_etaprime[2].SetPx(p3_gamma_etaprime * sintheta * cosphi);
   fProducts_etaprime[2].SetPy(p3_gamma_etaprime * sintheta * sinphi);
   fProducts_etaprime[2].SetPz(p3_gamma_etaprime * costheta);
   fProducts_etaprime[2].SetE(e3_gamma_etaprime);

   // electron 4-vectors in properly rotated virtual photon rest frame
   Double_t pRot1_etaprime[3] = {0.};
   Rot(pProd1_etaprime, pRot1_etaprime, costheta, -sintheta, -cosphi, -sinphi);
   Double_t pRot2_etaprime[3] = {0.};
   Rot(pProd2_etaprime, pRot2_etaprime, costheta, -sintheta, -cosphi, -sinphi);
   fProducts_etaprime[0].SetPx(pRot1_etaprime[0]);
   fProducts_etaprime[0].SetPy(pRot1_etaprime[1]);
   fProducts_etaprime[0].SetPz(pRot1_etaprime[2]);
   fProducts_etaprime[0].SetE(e1_etaprime);
   fProducts_etaprime[1].SetPx(pRot2_etaprime[0]);
   fProducts_etaprime[1].SetPy(pRot2_etaprime[1]);
   fProducts_etaprime[1].SetPz(pRot2_etaprime[2]);
   fProducts_etaprime[1].SetE(e1_etaprime);

   // boost the dielectron into the parent meson's rest frame 
   Double_t eLPparent_etaprime = TMath::Sqrt(p3_gamma_etaprime * p3_gamma_etaprime + epmass_etaprime * epmass_etaprime);
   TVector3 boostPair_etaprime( -1.0 * fProducts_etaprime[2].Px() / eLPparent_etaprime,
                       -1.0 * fProducts_etaprime[2].Py() / eLPparent_etaprime,
                       -1.0 * fProducts_etaprime[2].Pz() / eLPparent_etaprime);
   fProducts_etaprime[0].Boost(boostPair_etaprime);
   fProducts_etaprime[1].Boost(boostPair_etaprime);

   // boost all decay products into the lab frame
   TVector3 boostLab_etaprime(pparent->Px() / pparent->E(),
                     pparent->Py() / pparent->E(),
                     pparent->Pz() / pparent->E());

   fProducts_etaprime[0].Boost(boostLab_etaprime);
   fProducts_etaprime[1].Boost(boostLab_etaprime);
   fProducts_etaprime[2].Boost(boostLab_etaprime);

   }

//-----------------------------------------------------------------------------//
//                        Generate Phi Dalitz decay                            //
//-----------------------------------------------------------------------------//   

   else if(idpart==333){
   pmass_phi_dalitz = pparent->M();
   for(;;){
   // Sample the electron pair mass from a histogram
   epmass_phi_dalitz = fEPMassPhiDalitz->GetRandom();
   if(pmass_phi_dalitz-omass_eta>epmass_phi_dalitz && epmass_phi_dalitz/2.>emass)
   break;}

   // electron pair kinematics in virtual photon rest frame
   e1_phi = epmass_phi_dalitz / 2.;
   p1_phi = TMath::Sqrt((e1_phi + emass) * (e1_phi - emass));

   // momentum vectors of electrons in virtual photon rest frame
   Double_t pProd1_phi_dalitz[3] = {p1_phi * sintheta * cosphi,
                                    p1_phi * sintheta * sinphi,
                                    p1_phi * costheta};
   Double_t pProd2_phi_dalitz[3] = {-1.0 * p1_phi * sintheta * cosphi,
                                    -1.0 * p1_phi * sintheta * sinphi,
                                    -1.0 * p1_phi * costheta};

   // third child kinematics in parent meson rest frame
   e3_eta       = (pmass_phi_dalitz * pmass_phi_dalitz + omass_eta * omass_eta - epmass_phi_dalitz * epmass_phi_dalitz)/(2. * pmass_phi_dalitz);
   p3_eta       = TMath::Sqrt((e3_eta + omass_eta)  * (e3_eta - omass_eta));

   // third child 4-vector in parent meson rest frame
   fProducts_phi_dalitz[2].SetPx(p3_eta * sintheta * cosphi);
   fProducts_phi_dalitz[2].SetPy(p3_eta * sintheta * sinphi);
   fProducts_phi_dalitz[2].SetPz(p3_eta * costheta);
   fProducts_phi_dalitz[2].SetE(e3_eta);

   // electron 4-vectors in properly rotated virtual photon rest frame
   Double_t pRot1_phi_dalitz[3] = {0.};
   Rot(pProd1_phi_dalitz, pRot1_phi_dalitz, costheta, -sintheta, -cosphi, -sinphi);
   Double_t pRot2_phi_dalitz[3] = {0.};
   Rot(pProd2_phi_dalitz, pRot2_phi_dalitz, costheta, -sintheta, -cosphi, -sinphi);
   fProducts_phi_dalitz[0].SetPx(pRot1_phi_dalitz[0]);
   fProducts_phi_dalitz[0].SetPy(pRot1_phi_dalitz[1]);
   fProducts_phi_dalitz[0].SetPz(pRot1_phi_dalitz[2]);
   fProducts_phi_dalitz[0].SetE(e1_phi);
   fProducts_phi_dalitz[1].SetPx(pRot2_phi_dalitz[0]);
   fProducts_phi_dalitz[1].SetPy(pRot2_phi_dalitz[1]);
   fProducts_phi_dalitz[1].SetPz(pRot2_phi_dalitz[2]);
   fProducts_phi_dalitz[1].SetE(e1_phi);

   // boost the dielectron into the parent meson's rest frame
   Double_t eLPparent_phi = TMath::Sqrt(p3_eta * p3_eta + epmass_phi_dalitz * epmass_phi_dalitz);
   TVector3 boostPair_phi( -1.0 * fProducts_phi_dalitz[2].Px() / eLPparent_phi,
                           -1.0 * fProducts_phi_dalitz[2].Py() / eLPparent_phi,
                           -1.0 * fProducts_phi_dalitz[2].Pz() / eLPparent_phi);
   fProducts_phi_dalitz[0].Boost(boostPair_phi);
   fProducts_phi_dalitz[1].Boost(boostPair_phi);

   // boost all decay products into the lab frame
   TVector3 boostLab_phi_dalitz(pparent->Px() / pparent->E(),
                                pparent->Py() / pparent->E(),
                                pparent->Pz() / pparent->E());

   fProducts_phi_dalitz[0].Boost(boostLab_phi_dalitz);
   fProducts_phi_dalitz[1].Boost(boostLab_phi_dalitz);
   fProducts_phi_dalitz[2].Boost(boostLab_phi_dalitz);


//-----------------------------------------------------------------------------//
//                        Generate Phi resonance decay                         //
//-----------------------------------------------------------------------------//

      if(wp_phi!=0.0){
     // calculate phi mass   
         mp_phi = pparent->M();
         }
      else{
           Double_t x_phi=pparent->Px(); Double_t y_phi=pparent->Py(); Double_t z_phi=pparent->Pz();
           Double_t t_phi=pparent->E();
           Double_t p_phi=x_phi*x_phi+y_phi*y_phi+z_phi*z_phi;
           Double_t Q2_phi= abs((t_phi*t_phi)-(p_phi*p_phi));
           mp_phi = sqrt(Q2_phi);
          }
    
   if ( mp_phi< 2.*md_phi )
   {
    printf("Phi into ee Decay kinematically impossible! \n");
    return;
   }

   for( ;; ) {
   // Sample the electron pair mass from a histogram
   epmass_phi = fEPMassPhi->GetRandom();
   if(mp_phi < 2.*epmass_phi) break;
   }

   // electron pair kinematics in virtual photon rest frame
   Ed_phi = epmass_phi/2.;
   pd_phi = TMath::Sqrt((Ed_phi+md_phi)*(Ed_phi-md_phi));

   // momentum vectors of electrons in virtual photon rest frame
   Double_t pProd1_phi[3] = {pd_phi * sintheta * cosphi,
                             pd_phi * sintheta * sinphi,
                             pd_phi * costheta};
   Double_t pProd2_phi[3] = {-1.0 * pd_phi * sintheta * cosphi,
                             -1.0 * pd_phi * sintheta * sinphi,
                             -1.0 * pd_phi * costheta};

   // electron 4 vectors in properly rotated virtual photon rest frame
   Double_t pRot1_phi[3] = {0.};
   Rot(pProd1_phi, pRot1_phi, costheta, -sintheta, -cosphi, -sinphi);
   Double_t pRot2_phi[3] = {0.};
   Rot(pProd2_phi, pRot2_phi, costheta, -sintheta, -cosphi, -sinphi);
   fProducts_phi[0].SetPx(pRot1_phi[0]);
   fProducts_phi[0].SetPy(pRot1_phi[1]);
   fProducts_phi[0].SetPz(pRot1_phi[2]);
   fProducts_phi[0].SetE(Ed_phi);
   fProducts_phi[1].SetPx(pRot2_phi[0]);
   fProducts_phi[1].SetPy(pRot2_phi[1]);
   fProducts_phi[1].SetPz(pRot2_phi[2]);
   fProducts_phi[1].SetE(Ed_phi);

   // boost decay products into the lab frame
   TVector3 boostLab_phi(pparent->Px() / pparent->E(),
                     pparent->Py() / pparent->E(),
                     pparent->Pz() / pparent->E());

   fProducts_phi[0].Boost(boostLab_phi);
   fProducts_phi[1].Boost(boostLab_phi);

   }

//-----------------------------------------------------------------------------//
//                        Generate Jpsi resonance decay                        //
//-----------------------------------------------------------------------------//
   
   else if(idpart==443){
   // calculate jpsi mass
     if(wp_jpsi!=0.0){
        mp_jpsi = pparent->M();
        }
     else{
      /*Double_t x_jpsi=pparent->Px(); 
      Double_t y_jpsi=pparent->Py(); 
      Double_t z_jpsi=pparent->Pz();
      Double_t t_jpsi=pparent->E();
      Double_t p_jpsi=x_jpsi*x_jpsi+y_jpsi*y_jpsi+z_jpsi*z_jpsi;
      Double_t Q2_jpsi= abs((t_jpsi*t_jpsi)-(p_jpsi*p_jpsi));
      mp_jpsi = sqrt(Q2_jpsi);*/
       
      mp_jpsi = 3.096;

     }
    
     // daughter  
     if ( mp_jpsi < 2.*md_jpsi )
        {
         printf("JPsi into ee Decay kinematically impossible! \n");
         return;
        }

  for( ;; ) {
  // Sample the electron pair mass from a histogram 
  epmass_jpsi = fEPMassJPsi->GetRandom();
  if ( mp_jpsi < 2.*epmass_jpsi ) break;
  } 
  // electron pair kinematics in virtual photon rest frame
  Ed_jpsi = epmass_jpsi/2.;
  pd_jpsi = TMath::Sqrt((Ed_jpsi+md_jpsi)*(Ed_jpsi-md_jpsi));

  // momentum vectors of electrons in virtual photon rest frame 
  Double_t pProd1_jpsi[3] = {pd_jpsi * sintheta * cosphi,
                             pd_jpsi * sintheta * sinphi,
                             pd_jpsi * costheta};

  Double_t pProd2_jpsi[3] = {-1.0 * pd_jpsi * sintheta * cosphi,
                             -1.0 * pd_jpsi * sintheta * sinphi,
                             -1.0 * pd_jpsi * costheta};

  
  // electron 4 vectors in properly rotated virtual photon rest frame
  Double_t pRot1_jpsi[3] = {0.};
  Rot(pProd1_jpsi, pRot1_jpsi, costheta, -sintheta, -cosphi, -sinphi);
  Double_t pRot2_jpsi[3] = {0.};
  Rot(pProd2_jpsi, pRot2_jpsi, costheta, -sintheta, -cosphi, -sinphi);
  fProducts_jpsi[0].SetPx(pRot1_jpsi[0]);
  fProducts_jpsi[0].SetPy(pRot1_jpsi[1]);
  fProducts_jpsi[0].SetPz(pRot1_jpsi[2]);
  fProducts_jpsi[0].SetE(Ed_jpsi);
  fProducts_jpsi[1].SetPx(pRot2_jpsi[0]);
  fProducts_jpsi[1].SetPy(pRot2_jpsi[1]);
  fProducts_jpsi[1].SetPz(pRot2_jpsi[2]);
  fProducts_jpsi[1].SetE(Ed_jpsi);


  // boost decay products into the lab frame
  TVector3 boostLab_jpsi(pparent->Px() / pparent->E(),
                         pparent->Py() / pparent->E(),
                         pparent->Pz() / pparent->E());

  fProducts_jpsi[0].Boost(boostLab_jpsi);
  fProducts_jpsi[1].Boost(boostLab_jpsi);
           
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


   
