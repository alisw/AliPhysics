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

#include <TMath.h>
#include <TF1.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TParticle.h>
#include <TVector3.h>
#include <TRandom.h>
#include <TObjectTable.h>

#include "AliDecayerPolarized.h"
#include "AliLog.h"

ClassImp(AliDecayerPolarized)



//____________________________________________________________
AliDecayerPolarized::AliDecayerPolarized():
    fAlpha(0),
    fSystRef(kHelicity),
    fDecProd(kMuon),
    fPol(0),
    fMother(0),
    fDaughter1(0),
    fDaughter2(0)
{
// Default constructor

}

//____________________________________________________________
AliDecayerPolarized::AliDecayerPolarized(Double_t alpha, Polar_t systref, FinState_t decprod):
    fAlpha(alpha),
    fSystRef(systref),
    fDecProd(decprod),
    fPol(new TF1("dsigdcostheta","1.+[0]*x*x",-1.,1.)),
    fMother(0),
    fDaughter1(0),
    fDaughter2(0)
{
// Another constructor
    fPol->SetParameter(0,fAlpha);
    if(fDecProd!=kMuon && fDecProd!=kElectron)
	AliFatal("Only polarized decay into muons or electrons is implemented!");
}

AliDecayerPolarized::AliDecayerPolarized(const AliDecayerPolarized &decayer):
    AliDecayer(decayer),
    fAlpha(0),
    fSystRef(kHelicity),
    fDecProd(kMuon),
    fPol(new TF1("dsigdcostheta","1.+[0]*x*x",-1.,1.)),
    fMother(0),
    fDaughter1(0),
    fDaughter2(0)   
{
    // Copy Constructor
    decayer.Copy(*this);
}

//____________________________________________________________
AliDecayerPolarized::~AliDecayerPolarized()
{
// Destructor
    if(fPol) delete fPol; fPol=0;
    if(fMother) delete fMother; fMother=0;
    if(fDaughter1) delete fDaughter1; fDaughter1=0;
    if(fDaughter2) delete fDaughter2; fDaughter2=0;
}
//____________________________________________________________
void AliDecayerPolarized::Decay(Int_t ipart, TLorentzVector *p)
{
// Polarized 2- body decay 
    TDatabasePDG *pDataBase = TDatabasePDG::Instance();
    if(ipart!= (pDataBase->GetParticle("J/psi")->PdgCode()) &&
       ipart!= (pDataBase->GetParticle("psi'")->PdgCode())  &&
       ipart!= (pDataBase->GetParticle("Upsilon")->PdgCode())  &&
       ipart!= (pDataBase->GetParticle("Upsilon'")->PdgCode()))  
       AliFatal("Polarized decay only implemented for J/psi, psi', Upsilon and Upsilon' !");

    TParticlePDG *d1 = 0, *d2 = 0;
    if(fDecProd==kMuon){
      d1 = pDataBase->GetParticle("mu-");
      d2 = pDataBase->GetParticle("mu+");
    }
    else if(fDecProd==kElectron){
      d1 = pDataBase->GetParticle("e-");
      d2 = pDataBase->GetParticle("e+");
    }
// energies and momenta in lab frame 
    Double_t e1 = p->M() / 2.;
    Double_t e2 = e1; 
    Double_t p1 = TMath::Sqrt((e1 + d1->Mass())*(e1 - d1->Mass())); 

    Double_t costheta = fPol->GetRandom();
    
    Double_t sintheta = TMath::Sqrt((1. + costheta)*(1. - costheta));
    Double_t phi      = 2. * TMath::Pi() * gRandom->Rndm(); 
    Double_t px1 = p1 * sintheta * TMath::Cos(phi); 
    Double_t py1 = p1 * sintheta * TMath::Sin(phi); 
    Double_t pz1 = p1 * costheta; 

    TLorentzVector v1,v2, boosted1, boosted2;	
    v1.SetPxPyPzE(-px1,-py1,-pz1,e1);   //in the dimuon rest frame
    v2.SetPxPyPzE(px1,py1,pz1,e2); 
	
    TLorentzVector PProj, PTarg; // In the CM frame
    Float_t mp = 0.938;
    PProj.SetPxPyPzE(0.,0.,-7000.,TMath::Sqrt(7000.*7000.+mp*mp)); // projectile
    PTarg.SetPxPyPzE(0.,0.,7000.,TMath::Sqrt(7000.*7000.+mp*mp)); // target
 
    TVector3 betajpsicm;
    betajpsicm = (-1./p->E()*p->Vect());
      
    if(fSystRef == kHelicity) {
      //  polarization axis: direction of the JPsi in the CM 
      TVector3 v3jpsi = (p->Vect()).Unit();  
      v1.RotateUz(v3jpsi);  
      v2.RotateUz(v3jpsi); 	
    } else if (fSystRef == kColSop){
      //  polarization axis: bisector of proj and target in the dimuon rest frame
      PProj.Boost(betajpsicm);   //boost proj and targ from CM to DIMU rest frame
      PTarg.Boost(betajpsicm);

      TVector3 zaxisCS;
      zaxisCS=(((PProj.Vect()).Unit())-((PTarg.Vect()).Unit())).Unit();
      v1.RotateUz(zaxisCS);
      v2.RotateUz(zaxisCS);
    }
 
//    printf("v1 components (mu1 with jpsi at rest)%f %f %f %f\n",v1.Px(),v1.Py(),v1.Pz(),v1.E());
//    printf("v2 components (mu2 with jpsi at rest)%f %f %f %f\n",v2.Px(),v2.Py(),v2.Pz(),v2.E());

    v1.Boost(-betajpsicm);  //boost muons from DIMUON rest frame to CM
    v2.Boost(-betajpsicm); 

//    printf("v1 components (mu1 in polar. ref. syst.) %f %f %f %f\n",v1.Px(),v1.Py(),v1.Pz(),v1.E());
//    printf("v2 components (mu2 in polar. ref. syst.) %f %f %f %f\n",v2.Px(),v2.Py(),v2.Pz(),v2.E());

    Int_t status_decayed=11;
    Int_t status_undecayed=1;
        
   
    fMother = new TParticle(ipart,status_decayed,0,-1,2,3,p->Px(),p->Py(),p->Pz(),p->E(),0.,0.,0.,0);   
    fDaughter1 = new TParticle(d1->PdgCode(),status_undecayed,1,-1,0,0,v1.Px(),v1.Py(),v1.Pz(),v1.E(),0.,0.,0.,0);   
    fDaughter2 = new TParticle(d2->PdgCode(),status_undecayed,1,-1,0,0,v2.Px(),v2.Py(),v2.Pz(),v2.E(),0.,0.,0.,0);

}

//____________________________________________________________
Int_t AliDecayerPolarized::ImportParticles(TClonesArray *part)
{
// Return array of particles  
  TClonesArray &cloneparticles = *part;
  cloneparticles.Clear();
  
  new(cloneparticles[0])TParticle(*fMother);
  new(cloneparticles[1])TParticle(*fDaughter1);
  new(cloneparticles[2])TParticle(*fDaughter2);

  return part->GetEntries();
}

void  AliDecayerPolarized::Copy(TObject &) const
{
    //
    // Copy *this onto AliDecayerPolarized -- not implemented
    //
    Fatal("Copy","Not implemented!\n");
}

void AliDecayerPolarized::SetForceDecay(Int_t)
{
    // This method is dummy
}

void AliDecayerPolarized::ForceDecay()
{
    // This method is dummy
    AliWarning("Method not implemented for this class !\n");
}

Float_t AliDecayerPolarized::GetPartialBranchingRatio(Int_t)
{
    // This method is dummy
    return  1.;
}

Float_t AliDecayerPolarized::GetLifetime(Int_t)
{
    // This method is dummy
    AliWarning("Method not implemented for this class !\n");
    return -1.;
}

void AliDecayerPolarized::ReadDecayTable()
{
    // This method is dummy
    AliWarning("Method not implemented for this class !\n");
}

