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
 * about the suitability of this software for any purpose. It is         *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 *
 */

//_________________________________________________________________________
// Class for reading data (Kinematics) in order to do prompt gamma correlations
//  Class created from old AliPHOSGammaJet 
//  (see AliRoot versions previous Release 4-09)
//
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---

#include <TFile.h>
#include <TParticle.h>
#include <TH2.h>
#include <TChain.h>
#include <TRandom.h>

//---- ANALYSIS system ----
#include "AliGammaMCReader.h" 
#include "Riostream.h"
#include "AliLog.h"

ClassImp(AliGammaMCReader)

//____________________________________________________________________________
AliGammaMCReader::AliGammaMCReader() : 
  AliGammaReader(), fDecayPi0(0), fCheckOverlapping(kFALSE)
{
  //Ctor
  
  //Initialize parameters
  InitParameters();
  fDataType = kMC;  
  
}

//____________________________________________________________________________
AliGammaMCReader::AliGammaMCReader(const AliGammaMCReader & g) :   
  AliGammaReader(g), fDecayPi0(g.fDecayPi0), fCheckOverlapping(g.fCheckOverlapping)
{
  // cpy ctor
}

//_________________________________________________________________________
AliGammaMCReader & AliGammaMCReader::operator = (const AliGammaMCReader & source)
{
  // assignment operator

  if(&source == this) return *this;

  fDecayPi0 = source.fDecayPi0; 
  fCheckOverlapping = source.fCheckOverlapping;

  return *this;

}

//___________________________________________________________________________
void AliGammaMCReader::CaseDecayGamma(Int_t index, TParticle * particle, AliStack * sta,
				  TClonesArray * plEMCAL, Int_t &indexEMCAL,
				  TClonesArray * plPHOS, Int_t &indexPHOS){
  //In case pi0 are decayed by pythia, check if mother is pi0 and in such case look if 
  //there is overlapping. Send particle=pi0 if there is overlapping. 
  
  TParticle * pmother =sta->Particle(particle->GetFirstMother());
  if(pmother->GetPdgCode() == 111 && pmother->GetNDaughters() == 2) {//Do not consider 3 particle decay case
    Int_t idaug0 = pmother->GetDaughter(0);
    Int_t idaug1 = pmother->GetDaughter(1);
    TParticle * pdaug0 = sta -> Particle(idaug0);
    TParticle * pdaug1 = sta -> Particle(idaug1);
    
    if((index ==  idaug0 &&  pdaug0->Pt() > fNeutralPtCut) ||
       (index ==  idaug1 &&  pdaug0->Pt() <= fNeutralPtCut))//Check decay when first daughter arrives, do nothing with second.
      FillListWithDecayGammaOrPi0(pmother, pdaug0, pdaug1, plEMCAL, indexEMCAL, plPHOS, indexPHOS);
    
  }//mother is a pi0 with 2 daughters
  else{
    if(IsInPHOS(particle->Phi(),particle->Eta()))
      new((*plPHOS)[indexPHOS++])  TParticle(*particle) ;
    else if(IsInEMCAL(particle->Phi(),particle->Eta()))
      new((*plEMCAL)[indexEMCAL++])  TParticle(*particle) ;
  }

}

//___________________________________________________________________________
 void AliGammaMCReader::CaseGeantDecay(TParticle * particle, AliStack * stack,
				   TClonesArray * plEMCAL, Int_t &indexEMCAL,
				   TClonesArray * plPHOS, Int_t &indexPHOS){
   //Find decay gamma from pi0, decayed by GEANT and put them in the list.
   
   Int_t ndaug = particle->GetNDaughters() ;
   TParticle * d1 = new TParticle();
   TParticle * d2 = new TParticle();
   
   if(ndaug > 0){//At least 1 daugther
     d1 = stack->Particle(particle->GetDaughter(0));
     if (ndaug > 1 ) 
       d2 = stack->Particle(particle->GetDaughter(1));
     
     FillListWithDecayGammaOrPi0(particle, d1, d2, plEMCAL, indexEMCAL, plPHOS, indexPHOS);
     
   }// ndaugh > 0
 }

//___________________________________________________________________________
void AliGammaMCReader::CasePi0Decay(TParticle * pPi0, 
				    TClonesArray * plEMCAL, Int_t &indexEMCAL,
				    TClonesArray * plPHOS, Int_t &indexPHOS){
  
  //Decays pi0, see if aperture angle is small and then add the pi0 or the 2 gamma
  
  TLorentzVector lvPi0, lvGamma1, lvGamma2 ;
  //Double_t angle = 0;
  
  lvPi0.SetPxPyPzE(pPi0->Px(),pPi0->Py(),pPi0->Pz(),pPi0->Energy());

  //Decay
  MakePi0Decay(lvPi0,lvGamma1,lvGamma2);//,angle);
  
  //Check if Pi0 is in the acceptance of the calorimeters, if aperture angle is small, keep it
  TParticle * pPhoton1 = new TParticle(22,1,0,0,0,0,lvGamma1.Px(),lvGamma1.Py(),
				      lvGamma1.Pz(),lvGamma1.E(),0,0,0,0);   
  TParticle * pPhoton2 = new TParticle(22,1,0,0,0,0,lvGamma2.Px(),lvGamma2.Py(),
				      lvGamma2.Pz(),lvGamma2.E(),0,0,0,0);

  FillListWithDecayGammaOrPi0(pPi0, pPhoton1, pPhoton2, plEMCAL, indexEMCAL, plPHOS, indexPHOS);
 
}

//_______________________________________________________________
void AliGammaMCReader::InitParameters()
{
  
  //Initialize the parameters of the analysis.

  fDecayPi0 = kGeantDecay;
  fCheckOverlapping = kTRUE ;
}

//____________________________________________________________________________
void AliGammaMCReader::CreateParticleList(TObject * data, TObject *, 
					  TClonesArray * plCh, 
					  TClonesArray * plEMCAL,  
					  TClonesArray * plPHOS,
					  TClonesArray *plParton,TClonesArray *,TClonesArray *)
{
  //Create list of particles from EMCAL, PHOS and CTS. 
  AliStack * stack = (AliStack *) data ;
  Int_t indexCh     = plCh->GetEntries() ;
  Int_t indexEMCAL = plEMCAL->GetEntries() ;
  Int_t indexPHOS = plPHOS->GetEntries() ;
  Int_t indexParton = plParton->GetEntries() ;
  Int_t iParticle = 0 ;
  Double_t charge = 0.;
    
  for (iParticle=0 ; iParticle <  stack->GetNprimary() ; iParticle++) {
    TParticle * particle = stack->Particle(iParticle); 
    

    //Keep partons
    if(particle->GetStatusCode() == 21 && iParticle>=2){//All partons, not nucleus
      new((*plParton)[indexParton++])  TParticle(*particle) ;
    }

    //Keep Stable particles 
    if((particle->GetStatusCode() == 1) && (particle->Pt() > 0)){
      
      charge = TDatabasePDG::Instance()->GetParticle(particle->GetPdgCode())->Charge();
      
      //---------- Charged particles ----------------------
      if((charge != 0) && (particle->Pt() > fChargedPtCut)){
	//Particles in CTS acceptance
	if(TMath::Abs(particle->Eta())<fCTSEtaCut){  
	  //Fill lists
	  new((*plCh)[indexCh++])       TParticle(*particle) ;
	}
      }
      //-------------Neutral particles ----------------------
      else if((charge == 0) && particle->Pt() > fNeutralPtCut &&  
	      TMath::Abs(particle->GetPdgCode())>16){//Avoid neutrinos
	
	if(particle->GetPdgCode()!=111){
	  //In case that we force PYTHIA to decay pi0, and we want to check the overlapping of 
	  // the decay gamma.
	  if(particle->GetPdgCode() == 22 && fDecayPi0==kDecayGamma){
	    CaseDecayGamma(iParticle,particle, stack,plEMCAL, indexEMCAL, plPHOS, indexPHOS); //If pythia decays pi0
	  }
	  else{
	    if(IsInPHOS(particle->Phi(),particle->Eta()))
	      new((*plPHOS)[indexPHOS++])  TParticle(*particle) ;
	    else if(IsInEMCAL(particle->Phi(),particle->Eta()))
	      new((*plEMCAL)[indexEMCAL++])  TParticle(*particle) ;
	  }
	}//no pi0
	else{
	  if(fDecayPi0 == kNoDecay){//keep the pi0 do not decay
	    if(IsInPHOS(particle->Phi(),particle->Eta()))
	      new((*plPHOS)[indexPHOS++])  TParticle(*particle) ;
	    else if(IsInEMCAL(particle->Phi(),particle->Eta()))
	      new((*plEMCAL)[indexEMCAL++])  TParticle(*particle) ;
	  }
	  else if(fDecayPi0 == kDecay)
	    CasePi0Decay(particle,plEMCAL,indexEMCAL,plPHOS, indexPHOS);
	  else if(fDecayPi0 == kGeantDecay)
	    CaseGeantDecay(particle, stack,plEMCAL, indexEMCAL, plPHOS, indexPHOS);
	}//pi0  
      }//neutral particle
    }//stable particle
  }//particle loop
}


//___________________________________________________________________________
void AliGammaMCReader::FillListWithDecayGammaOrPi0(TParticle * pPi0, TParticle * pdaug0, TParticle * pdaug1,
				   TClonesArray * plEMCAL, Int_t &indexEMCAL,
				   TClonesArray * plPHOS, Int_t &indexPHOS){

  //Check if decay gamma overlapp in calorimeter, in such case keep the pi0, if not keep both photons.
  
  Bool_t  overlap = kFALSE ;
  TLorentzVector lv1 , lv2 ;
  pdaug0->Momentum(lv1);
  pdaug1->Momentum(lv2);
  Double_t angle = lv1.Angle(lv2.Vect());
  
  if(fCheckOverlapping){//Check if decay products overlapp
    if(IsInEMCAL(pPi0->Phi(), pPi0->Eta())){
      if (angle < fEMCALMinAngle){
	new((*plEMCAL)[indexEMCAL++])       TParticle(*pPi0) ;
	overlap = kTRUE;
      }
    }
    else if(IsInPHOS(pPi0->Phi(), pPi0->Eta())){
      if (angle < fPHOSMinAngle){
	new((*plPHOS)[indexPHOS++])       TParticle(*pPi0) ;
	overlap = kTRUE;
      }
    }
  }//fCheckOverlapping
  
  //Fill with gammas if not overlapp
  if(!overlap){
    if(pdaug0->GetPdgCode() == 22 || TMath::Abs(pdaug0->GetPdgCode() ) == 11 ){
      if(IsInEMCAL(pdaug0->Phi(),pdaug0->Eta()) &&  pdaug0->Pt() > fNeutralPtCut)
	new((*plEMCAL)[indexEMCAL++])       TParticle(*pdaug0) ;
      else if(IsInPHOS(pdaug0->Phi(),pdaug0->Eta()) &&  pdaug0->Pt() > fNeutralPtCut)
	new((*plPHOS)[indexPHOS++])       TParticle(*pdaug0) ;
    }
    
    if(pdaug1->GetPdgCode() == 22 || TMath::Abs(pdaug1->GetPdgCode() ) == 11 ){
      if(IsInEMCAL(pdaug1->Phi(),pdaug1->Eta()) &&  pdaug1->Pt() > fNeutralPtCut)
	new((*plEMCAL)[indexEMCAL++])       TParticle(*pdaug1) ;
      else if(IsInPHOS(pdaug1->Phi(),pdaug1->Eta()) &&  pdaug1->Pt() > fNeutralPtCut)
	new((*plPHOS)[indexPHOS++])       TParticle(*pdaug1) ;
    }
  }// overlap?
 }


 
//___________________________________________________________________________
Bool_t  AliGammaMCReader::IsInEMCAL(Double_t phi, Double_t eta){
  //Check if particle is in EMCAL acceptance
  if(phi<0)
     phi+=TMath::TwoPi();
     if( phi > fPhiEMCALCut[0] && phi < fPhiEMCALCut[1] && 
       TMath::Abs(eta)<fEMCALEtaCut) return kTRUE ;
  else  return kFALSE;     
  
  return kFALSE ;
}

//___________________________________________________________________________
Bool_t  AliGammaMCReader::IsInPHOS(Double_t phi, Double_t eta){
  //Check if particle is in EMCAL acceptance
  if(phi<0)
    phi+=TMath::TwoPi();
  if( phi > fPhiPHOSCut[0] && phi < fPhiPHOSCut[1] && 
      TMath::Abs(eta)<fPHOSEtaCut) return kTRUE ;
  else  return kFALSE;
  
  return kFALSE ;
}

//____________________________________________________________________________
void AliGammaMCReader::MakePi0Decay(TLorentzVector &p0, TLorentzVector &p1, 
				TLorentzVector &p2){//, Double_t &angle) {
  // Perform isotropic decay pi0 -> 2 photons
  // p0 is pi0 4-momentum (inut)
  // p1 and p2 are photon 4-momenta (output)
  //  cout<<"Boost vector"<<endl;
  Double_t mPi0 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
  TVector3 b = p0.BoostVector();
  //cout<<"Parameters"<<endl;
  //Double_t mPi0   = p0.M();
  Double_t phi    = TMath::TwoPi() * gRandom->Rndm();
  Double_t cosThe = 2 * gRandom->Rndm() - 1;
  Double_t cosPhi = TMath::Cos(phi);
  Double_t sinPhi = TMath::Sin(phi);
  Double_t sinThe = TMath::Sqrt(1-cosThe*cosThe);
  Double_t ePi0   = mPi0/2.;
  //cout<<"ePi0 "<<ePi0<<endl;
  //cout<<"Components"<<endl;
  p1.SetPx(+ePi0*cosPhi*sinThe);
  p1.SetPy(+ePi0*sinPhi*sinThe);
  p1.SetPz(+ePi0*cosThe);
  p1.SetE(ePi0);
  //cout<<"p1: "<<p1.Px()<<" "<<p1.Py()<<" "<<p1.Pz()<<" "<<p1.E()<<endl;
  //cout<<"p1 Mass: "<<p1.Px()*p1.Px()+p1.Py()*p1.Py()+p1.Pz()*p1.Pz()-p1.E()*p1.E()<<endl;
  p2.SetPx(-ePi0*cosPhi*sinThe);
  p2.SetPy(-ePi0*sinPhi*sinThe);
  p2.SetPz(-ePi0*cosThe);
  p2.SetE(ePi0);
  //cout<<"p2: "<<p2.Px()<<" "<<p2.Py()<<" "<<p2.Pz()<<" "<<p2.E()<<endl;
  //cout<<"p2 Mass: "<<p2.Px()*p2.Px()+p2.Py()*p2.Py()+p2.Pz()*p2.Pz()-p2.E()*p2.E()<<endl;
  //cout<<"Boost "<<b.X()<<" "<<b.Y()<<" "<<b.Z()<<endl;
  p1.Boost(b);
  //cout<<"p1: "<<p1.Px()<<" "<<p1.Py()<<" "<<p1.Pz()<<" "<<p1.E()<<endl;
  p2.Boost(b);
  //cout<<"p2: "<<p2.Px()<<" "<<p2.Py()<<" "<<p2.Pz()<<" "<<p2.E()<<endl;
  //cout<<"angle"<<endl;
  //angle = p1.Angle(p2.Vect());
  //cout<<angle<<endl;
}

//________________________________________________________________
void AliGammaMCReader::Print(const Option_t * opt) const
{
  
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  Info("Print", "%s %s", GetName(), GetTitle() ) ;
  
  printf("Decay Pi0?          : %d\n", fDecayPi0) ;
  printf("Check Overlapping?          : %d\n", fCheckOverlapping) ;

}




