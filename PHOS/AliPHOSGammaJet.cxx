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

//_________________________________________________________________________
// Class for the analysis of gamma-jet correlations 
//   
// 
//*-- Author: Gustavo Conesa & Yves Schutz (IFIC, CERN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include "TRandom.h"

#include "TLorentzVector.h"
#include "TList.h"
#include "TParticle.h"
#include "AliPHOSGammaJet.h" 
#include "AliPHOSGetter.h" 

ClassImp(AliPHOSGammaJet)

//____________________________________________________________________________
AliPHOSGammaJet::AliPHOSGammaJet() {
  // ctor
}

//____________________________________________________________________________
AliPHOSGammaJet::AliPHOSGammaJet(const TString inputfilename) {
  // ctor 
  AliPHOSGetter::Instance(inputfilename) ; 
}

//____________________________________________________________________________
AliPHOSGammaJet::AliPHOSGammaJet(const AliPHOSGammaJet & gj) : TTask(gj) {
}

//____________________________________________________________________________
AliPHOSGammaJet::~AliPHOSGammaJet() {
}

//____________________________________________________________________________
void AliPHOSGammaJet::Exec(Option_t * opt) 
{
  // does the job
  
  if (! opt) 
    return ; 

  AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
  Int_t nEvent1 = gime->MaxEvent() ;   
  Int_t iEvent = 0 ; 
  for ( iEvent = 0 ; iEvent < nEvent1 ; iEvent++) {
    if (iEvent <= 100 || iEvent%10 == 0)
      Info("Exec", "Event %d", iEvent) ;     
    TParticle * particle = 0 ;
    //-----------------Fill list with particles--------------------
    TLorentzVector pPi0, pGamma1, pGamma2 ;
    Double_t angle = 0 ;
    TList particleList ;
    Int_t n = -1; 
    gime->Event(iEvent, "X") ; 
    Int_t  nparticles = gime->NPrimaries() ; 
    Int_t iParticle=0 ;
    for (iParticle=0 ; iParticle < nparticles ; iParticle++) {
      //Keep original partons
      particle = gime->Primary(iParticle) ; 
      if((particle->GetStatusCode()== 21)){
	particleList.Add(particle);
	n++;
      }
      //Keep Stable particles within eta range
      Float_t etacut = 0. ; 
      if((particle->GetStatusCode() == 1)&&
	 (TMath::Abs(particle->Eta())<etacut)){
	// Keep particles different from Pi0
	if(particle->GetPdgCode() != 111){
	  particleList.Add(particle);
	  n++;
	} 
	//Decay Pi0 and keep it with different status name
	//Keep decay photons
	if(particle->GetPdgCode() == 111) {
	  //cout<<"Pi0 "<<endl;
	  n += 3 ; 
	  particle->Momentum(pPi0);
	  Pi0Decay(particle->GetMass(),pPi0,pGamma1,pGamma2,angle);	
	  TParticle * photon1 = new TParticle(22,1,0,0,0,0,pGamma1.Px(),pGamma1.Py(),
					      pGamma1.Pz(),pGamma1.E(),0,0,0,0);
	  TParticle * photon2 = new TParticle(22,1,0,0,0,0,pGamma2.Px(),pGamma2.Py(),
					      pGamma2.Pz(),pGamma2.E(),0,0,0,0);
	  photon1->SetWeight(1);
	  photon2->SetWeight(2);
	  particle->SetStatusCode(2);
	  particleList.Add(particle);
	  particleList.Add(photon1);
	  particleList.Add(photon2);
	  //photon1->Print();
	  //photon2->Print();  
	}//if pi0
      }//final particle etacut
    }//for (iParticle<nParticle)
    TLorentzVector gamma,charge,pi0, gammapair;
    Int_t idg = -1;
    GetGammaJet(particleList,gamma, idg);
    GetLeadingCharge(particleList,charge, idg);
    GetLeadingPi0(particleList,pi0);
    Info("Pi0Decay", "Gamma: %f %d", gamma.Energy(), idg) ;
    Info("Pi0Decay", "Charge: %f", charge.Energy()) ;
    Info("Pi0Decay", "Pi0: %f", pi0.Energy()) ;
    //    GetLeadingGammaPair(particleList, gammapair, idg, 
    //			0.04,0.2, 1.0,0.13,0.14);
    Info("Pi0Decay", "Pair: %f", gammapair.Energy()) ;
  }//loop: events
}    

//____________________________________________________________________________
void AliPHOSGammaJet::Pi0Decay(Double_t mPi0, TLorentzVector &p0, 
			       TLorentzVector &p1, TLorentzVector &p2, Double_t &angle) {
  // Perform isotropic decay pi0 -> 2 photons
  // p0 is pi0 4-momentum (input)
  // p1 and p2 are photon 4-momenta (output)
  //  cout<<"Boost vector"<<endl;
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
  angle = p1.Angle(p2.Vect());
  //cout<<angle<<endl;
}
//____________________________________________________________________________
void AliPHOSGammaJet::GetGammaJet(TList &particleList, TLorentzVector &gamma, Int_t & id) 
{
  // Get the lists of jet particles and gamma
  TParticle *particle = 0x0;
  
  Int_t iPrimary=-1, id_motherg = -1;
  
  TIter next(&particleList);
  while ( (particle = (TParticle*)next()) ) {
    iPrimary++;  
    Int_t ksCode = particle->GetStatusCode();
    Int_t iMother= particle->GetMother(0);
    
    if (ksCode == 21 && iMother == -1)
      if(particle->GetPdgCode()==22){
	id_motherg = iPrimary;
	//cout<<"idmother "<<id_motherg<<endl;
      }
    if(ksCode == 1){
      
      if(id_motherg == particle->GetMother(0)){
 	particle->Momentum(gamma);
 	id = iPrimary;
 	break;
      }
    }// kscode == 1
  }// while
}
 
//____________________________________________________________________________
void AliPHOSGammaJet::GetLeadingCharge(TList &particleList, TLorentzVector &charge, Int_t & id) 
{
  // Gets the leading particle from the list of charged particles
  TParticle *particle = 0x0;
  
  Int_t iPrimary=-1;
  Double_t ptmax = 0, pti = 0;
  TIter next(&particleList);
  while ( (particle = (TParticle*)next()) ) {
    iPrimary++;  
    Int_t ksCode = particle->GetStatusCode();
    
    if((ksCode == 1)&&(id != iPrimary)
       &&(particle->GetPDG(0)->Charge()!=0)){
      pti = particle->Pt(); 
      if(pti> ptmax){
	ptmax = pti;
 	particle->Momentum(charge);
      }//ptmax   
    }// kscode == 1
  }// while
}

//____________________________________________________________________________
void AliPHOSGammaJet::GetLeadingPi0(TList &particleList, TLorentzVector &pi0) 
{
  // Gets the leading pi0 from the list of particles
  TParticle *particle = 0x0;
  
  Int_t iPrimary=-1;
  Double_t ptmax = 0, pti = 0;
  TIter next(&particleList);
  while ( (particle = (TParticle*)next()) ) {
    iPrimary++;  
    Int_t ksCode = particle->GetStatusCode();
    
    if((ksCode == 2))
      {
	pti = particle->Pt(); 
	if(pti> ptmax){
	  ptmax = pti;
	  particle->Momentum(pi0);
	}//ptmax   
      }// kscode == 1
  }// while
}

//____________________________________________________________________________
//  void AliPHOSGammaJet::GetLeadingGammaPair(TList &particleList, TLorentzVector &gammapair, Int_t & id, 
//  			 Double_t & thetacut,Double_t & ratiocut1, Double_t & ratiocut2,
//  			 Double_t & invmasscut1,Double_t & invmasscut2) 
//  {
//   TParticle *particle = 0x0;
  
//   Int_t  iPrimary=-1;
//   Double_t El = 0, E12 = 0;
//   TLorentzVector gamma_i,gamma_j;
//   TIter next(&particleList);
//   while ( (particle = (TParticle*)next()) ) {
//     iPrimary++;	  
//     Int_t ksCode = particle->GetStatusCode();
//     Int_t ksPdg = particle->GetPdgCode();
//     if((ksCode == 1) && (iPrimary != id) && (ksPdg == 22)){
//       particle->Momentum(gamma_i);
//       Int_t jPrimary=-1;
//       TIter next2(&particleList);
//       while ( (particle = (TParticle*)next2()) ) {
// 	jPrimary++;
// 	if(jPrimary>iPrimary){
// 	  ksCode = particle->GetStatusCode();
// 	  ksPdg = particle->GetPdgCode();
// 	  if((ksCode == 1) && (iPrimary != id) && (ksPdg == 22)){
// 	    particle->Momentum(gamma_j);
// 	    if(gamma_j.Angle(gamma_i.Vect())<thetacut){
// 	      Float_t invmass = (gamma_i+gamma_j).M();
// 	      h_invmass->Fill(Eg,invmass);
// 	      if((invmass>invmasscut1) && (invmass<invmasscut2)){
// 		E12 =  (gamma_i+gamma_j).Energy(); 
// 		if(E12>El && (E12/Eg>ratiocut1) && (E12/Eg<ratiocut2)){
// 		  //cout<<E12<<" "<<E12/Eg<<endl;
// 		  El = E12;
// 		  id_i = iPrimary;
// 		  id_j = jPrimary;
// 		  gammapair = gamma_i+gamma_j;
// 		}//E12>El && (E12/Eg>0.2 && E12/Eg<1.)
// 	      }//(invmass>0.125) && (invmass<0.145)
// 	    }//gamma_j.Angle(gamma_i.Vect())<0.04
// 	  }//(ksCode == 1)
// 	}
//       }//while
//       //	    cout<<"jPrimary "<<jPrimary<<endl;
//     }// if kscode 1
//   }//while
//   //	cout<<"iPrimary "<<iPrimary<<endl;
//  }
