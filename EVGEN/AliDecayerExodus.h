#ifndef ALIDECAYEREXODUS_H
#define ALIDECAYEREXODUS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *  * See cxx source for full Copyright notice                               */

/* $Id$ */

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



#include "AliDecayer.h"
#include <TLorentzVector.h>
#include <TF1.h>

class TH1F;
class TClonesArray;

class AliDecayerExodus : public AliDecayer
{
 public:
    AliDecayerExodus();
    virtual void    Init();
    virtual void    Decay(Int_t idpart,TLorentzVector* pparent);
    virtual Int_t   ImportParticles(TClonesArray *particles);
    virtual void    SetForceDecay(Int_t)                      {;}
    virtual void    ForceDecay()                              {;}
    virtual Float_t GetPartialBranchingRatio(Int_t /*ipart*/) {return -1;}
    virtual Float_t GetLifetime(Int_t /*kf*/)                 {return -1;}
    virtual void    ReadDecayTable()                          {;}
    
    virtual TH1F*   ElectronPairMassHistoPion()          {return  fEPMassPion;}
    virtual TH1F*   ElectronPairMassHistoEta()           {return  fEPMassEta;}
    virtual TH1F*   ElectronPairMassHistoEtaPrime()      {return  fEPMassEtaPrime;}
    virtual TH1F*   ElectronPairMassHistoRho()           {return  fEPMassRho;}
    virtual TH1F*   ElectronPairMassHistoOmega()         {return  fEPMassOmega;}
    virtual TH1F*   ElectronPairMassHistoOmegaDalitz()   {return  fEPMassOmegaDalitz;}
    virtual TH1F*   ElectronPairMassHistoPhi()           {return  fEPMassPhi;}
    virtual TH1F*   ElectronPairMassHistoPhiDalitz()     {return  fEPMassPhiDalitz;}
    virtual TH1F*   ElectronPairMassHistoJPsi()          {return  fEPMassJPsi;}

    virtual void    Decay(TClonesArray* array);

    virtual const   TLorentzVector* Products_pion()         const {return fProducts_pion;}
    virtual const   TLorentzVector* Products_eta()          const {return fProducts_eta;}
    virtual const   TLorentzVector* Products_etaprime()     const {return fProducts_etaprime;}
    virtual const   TLorentzVector* Products_rho()          const {return fProducts_rho;}
    virtual const   TLorentzVector* Products_omega()        const {return fProducts_omega;}
    virtual const   TLorentzVector* Products_omega_dalitz() const {return fProducts_omega_dalitz;}
    virtual const   TLorentzVector* Products_phi()          const {return fProducts_phi;}
    virtual const   TLorentzVector* Products_phi_dalitz()   const {return fProducts_phi_dalitz;}
    virtual const   TLorentzVector* Products_jpsi()         const {return fProducts_jpsi;}

    virtual void    Copy(TObject&) const;

 protected:
    // Histograms for electron pair mass
    TH1F*         fEPMassPion;          
    TH1F*         fEPMassEta;       
    TH1F*         fEPMassEtaPrime;
    TH1F*         fEPMassRho;
    TH1F*         fEPMassOmega;
    TH1F*         fEPMassOmegaDalitz;
    TH1F*         fEPMassPhi;
    TH1F*         fEPMassPhiDalitz;
    TH1F*         fEPMassJPsi;

    TF1* fPol;

    // Decay products
    TLorentzVector  fProducts_pion[3];  
    TLorentzVector  fProducts_eta[3];  
    TLorentzVector  fProducts_etaprime[3];
    TLorentzVector  fProducts_rho[2];
    TLorentzVector  fProducts_omega[2];
    TLorentzVector  fProducts_omega_dalitz[3];
    TLorentzVector  fProducts_phi[2];
    TLorentzVector  fProducts_phi_dalitz[3];
    TLorentzVector  fProducts_jpsi[2];

    Bool_t fInit;

 private:
    Double_t GounarisSakurai(Float_t mass, Double_t vmass, Double_t vwidth, Double_t emass);
    Double_t Lorentz(Float_t mass, Double_t vmass, Double_t vwidth); 
    virtual void    Rot(Double_t pin[3], Double_t pout[3],
                        Double_t costheta, Double_t sintheta,
                        Double_t cosphi, Double_t sinphi) const;
    AliDecayerExodus(const AliDecayerExodus &decayer);
    AliDecayerExodus & operator=(const AliDecayerExodus & rhs);


    ClassDef(AliDecayerExodus, 1) // AliDecayer implementation using Exodus  
};
#endif







