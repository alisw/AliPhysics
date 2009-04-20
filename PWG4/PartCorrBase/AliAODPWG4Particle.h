#ifndef ALIAODPWG4PARTICLE_H
#define ALIAODPWG4PARTICLE_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAODPWG4Particle.h  $ */

//-------------------------------------------------------------------------
//     Copy of AOD photon class, adapted for particle identification
//     and correlations analysis
//     Author: Yves Schutz, CERN, Gustavo Conesa, INFN
//-------------------------------------------------------------------------

//-- ROOT system --
#include <TLorentzVector.h>
class TString;

//-- Analysis system
#include "AliVParticle.h"

class AliAODPWG4Particle : public AliVParticle {
  
 public:
  AliAODPWG4Particle();
  AliAODPWG4Particle(Double_t px, Double_t py, Double_t pz, Double_t e);
  AliAODPWG4Particle(TLorentzVector & p);  
  virtual ~AliAODPWG4Particle();
  AliAODPWG4Particle(const AliAODPWG4Particle& photon); 
  AliAODPWG4Particle& operator=(const AliAODPWG4Particle& photon);
  
  // AliVParticle methods
  virtual Double_t Px()         const { return fMomentum->Px();      }
  virtual Double_t Py()         const { return fMomentum->Py();      }
  virtual Double_t Pz()         const { return fMomentum->Pz();      }
  virtual Double_t Pt()         const { return fMomentum->Pt();      }
  virtual Double_t P()          const { return fMomentum->P();       }
  virtual Bool_t   PxPyPz(Double_t p[3]) const { p[0] = Px(); p[1] = Py(); p[2] = Pz(); return kTRUE; }
  virtual Double_t OneOverPt()  const { return 1. / fMomentum->Pt(); }
  virtual Double_t Phi()        const;
  virtual Double_t Theta()      const { return fMomentum->Theta();   }
  virtual Double_t E()          const { return fMomentum->E();       }
  virtual Double_t M()          const { return fMomentum->M();       }
  virtual Double_t Eta()        const { return fMomentum->Eta();     }
  virtual Double_t Y()          const { return fMomentum->Rapidity();}
  virtual Double_t Xv()         const {return -999.;} // put reasonable values here
  virtual Double_t Yv()         const {return -999.;} //
  virtual Double_t Zv()         const {return -999.;} //
  virtual Bool_t   XvYvZv(Double_t x[3]) const { x[0] = Xv(); x[1] = Yv(); x[2] = Zv(); return kTRUE; }  
  virtual void     Print(Option_t* /*option*/) const;
  
  //
  //Dummy
  virtual Short_t Charge()      const { return 0;}
  virtual const Double_t* PID() const { return NULL;}
  //
  
  virtual Int_t GetPdg()   const {return fPdg ; }
  virtual Int_t GetTag()   const {return fTag ; }
  virtual Int_t GetLabel()   const {return fLabel ; }
  virtual Int_t GetCaloLabel (Int_t i) const {return fCaloLabel[i];}
  virtual Int_t GetTrackLabel(Int_t i) const {return fTrackLabel[i];}
  virtual TString GetDetector()   const {return fDetector ; }

  virtual Bool_t   GetDispBit(void) const {return fDisp;}
  virtual Bool_t   GetTOFBit(void) const {return fTof;}
  virtual Bool_t   GetChargedBit(void) const {return fCharged;}
  virtual Int_t    DistToBad() const  {return fBadDist ;} 
  
  virtual void SetPdg(Int_t pdg)   {fPdg = pdg ; }
  virtual void SetTag(Int_t tag)   {fTag = tag ; }
  virtual void SetLabel(Int_t l)   {fLabel = l ; }
  virtual void SetCaloLabel (Int_t a, Int_t b)   {fCaloLabel [0] = a; fCaloLabel [1] = b  ; }
  virtual void SetTrackLabel(Int_t a, Int_t b)   {fTrackLabel[0] = a; fTrackLabel[1] = b  ; }
  virtual void SetDetector(TString d)   {fDetector = d ; }
  
  virtual void SetDispBit(Bool_t chi2){fDisp = chi2 ;} 
  virtual void SetTOFBit(Bool_t tof){fTof = tof ;} 
  virtual void SetChargedBit(Bool_t ch){fCharged = ch; }
  virtual void SetDistToBad(Int_t dist){fBadDist=dist;} 
  
  TLorentzVector	* Momentum() const {return fMomentum ; }
  virtual void SetMomentum(TLorentzVector *lv) {fMomentum = lv ; }
  
  Bool_t IsPIDOK(const Int_t ipid, const Int_t pdgwanted) const;
  
 private:
  TLorentzVector* fMomentum;  // Photon 4-momentum vector
  Int_t           fPdg;       // id of particle
  Int_t           fTag;       // tag of particle (decay, fragment, prompt photon)
  Int_t           fLabel;     // MC label
  Int_t		  fCaloLabel[2];  // CaloCluster index, 1 for photons, 2 for pi0.
  Int_t           fTrackLabel[2]; // Track lable, 1 for pions, 2 for conversion photons 
  TString         fDetector;  // Detector where particle was measured.
  
  Bool_t	  fDisp ;     //Dispersion bit
  Bool_t	  fTof ;      //TOF bit
  Bool_t	  fCharged ;  //Charged bit
  Int_t           fBadDist ;  //Distance to bad module in module units
  
  ClassDef(AliAODPWG4Particle,1);
};

inline Double_t AliAODPWG4Particle::Phi() const
{
  // Return phi
  Double_t phi = fMomentum->Phi();
  if (phi < 0.) phi += 2. * TMath::Pi();
  return phi;
}

#endif //ALIAODPWG4PARTICLE_H
