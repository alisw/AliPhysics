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

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.1.2.1  2007/07/26 10:32:09  schutz
 * new analysis classes in the the new analysis framework
 *
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
  AliGammaReader(),
  fEMCALIPDistance(0.),   fPHOSIPDistance(0.), 
  fEMCALMinDistance(0.),   fPHOSMinDistance(0.), 
  fDecayPi0(kFALSE)
{
  //Ctor
  
  //Initialize parameters
  fDataType = kMC;
  InitParameters();
  
  
}

//____________________________________________________________________________
AliGammaMCReader::AliGammaMCReader(const AliGammaMCReader & g) :   
  AliGammaReader(g),
  fEMCALIPDistance(g.fEMCALIPDistance),  fPHOSIPDistance(g.fPHOSIPDistance),  
  fEMCALMinDistance(g.fEMCALMinDistance),  fPHOSMinDistance(g.fPHOSMinDistance), 
  fDecayPi0(g.fDecayPi0)
{
  // cpy ctor
}

//_________________________________________________________________________
AliGammaMCReader & AliGammaMCReader::operator = (const AliGammaMCReader & source)
{
  // assignment operator

  if(&source == this) return *this;

  fEMCALIPDistance = source.fEMCALIPDistance; 
  fPHOSIPDistance = source.fPHOSIPDistance; 
  fEMCALMinDistance = source.fEMCALMinDistance; 
  fPHOSMinDistance = source.fPHOSMinDistance; 
  fDecayPi0 = source.fDecayPi0;
 

  return *this;

}


//____________________________________________________________________________
void AliGammaMCReader::CreateParticleList(TObject * data, TObject *, 
					 TClonesArray * plCh, 
					 TClonesArray * plEMCAL,  
					  TClonesArray * plPHOS,
					  TClonesArray *plParton)
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
    if((particle->GetStatusCode() == 0) && (particle->Pt() > 0)){
      
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
	  if(IsInPHOS(particle->Phi(),particle->Eta()))
	    new((*plPHOS)[indexPHOS++])  TParticle(*particle) ;
	  else if(IsInEMCAL(particle->Phi(),particle->Eta()))
	    new((*plEMCAL)[indexEMCAL++])  TParticle(*particle) ;
	}//no pi0
	else{
	  if(fDecayPi0 == kNoDecay){//keep the pi0 do not decay
	    if(IsInPHOS(particle->Phi(),particle->Eta()))
	      new((*plPHOS)[indexPHOS++])  TParticle(*particle) ;
	    else if(IsInEMCAL(particle->Phi(),particle->Eta()))
	      new((*plEMCAL)[indexEMCAL++])  TParticle(*particle) ;
	  }
	  else if(fDecayPi0 == kDecay)
	    MakePi0Decay(particle,plEMCAL,indexEMCAL,plPHOS, indexPHOS);
	  else if(fDecayPi0 == kGeantDecay)
	    SetGeantDecay(particle, stack,plEMCAL, indexEMCAL, plPHOS, indexPHOS);
	}//pi0  
      }//neutral particle
    }//stable particle
  }//particle loop
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

//___________________________________________________________________________
void AliGammaMCReader::SetGeantDecay(TParticle * particle, AliStack * stack,
				     TClonesArray * plEMCAL, Int_t &indexEMCAL,
				     TClonesArray * plPHOS, Int_t &indexPHOS){
  //Find decay gamma from pi0 and put them in the list.
  
  Int_t ndaug = particle->GetNDaughters() ;
  if(ndaug<=2 && ndaug >0){//At least 1 daugther
    TParticle * d1 = stack->Particle(particle->GetDaughter(0));
    if(d1->GetPdgCode()==22){
      if(IsInEMCAL(d1->Phi(),d1->Eta()))
	new((*plEMCAL)[indexEMCAL++])       TParticle(*d1) ;
      else if(IsInPHOS(d1->Phi(),d1->Eta()))
	new((*plPHOS)[indexPHOS++])       TParticle(*d1) ;
    }
    
    if(ndaug>1){//second daugther if present
      TParticle * d2 = stack->Particle(particle->GetDaughter(1));
      if(IsInEMCAL(d2->Phi(),d2->Eta()))
	new((*plEMCAL)[indexEMCAL++])       TParticle(*d2) ;
      else if(IsInPHOS(d2->Phi(),d2->Eta()))
	new((*plPHOS)[indexPHOS++])       TParticle(*d2) ;
    }
  }
}

//___________________________________________________________________________
void AliGammaMCReader::MakePi0Decay(TParticle * particle, 
				    TClonesArray * plEMCAL, Int_t &indexEMCAL,
				    TClonesArray * plPHOS, Int_t &indexPHOS){
  
  //Decays pi0, see if aperture angle is small and then add the pi0 or the 2 gamma
  
  TLorentzVector pPi0, pGamma1, pGamma2 ;
  Double_t angle = 0, cellDistance = 0.;
  Bool_t checkPhoton = kTRUE;
  
  pPi0.SetPxPyPzE(particle->Px(),particle->Py(),particle->Pz(),particle->Energy());
  
  //Decay
  Pi0Decay(pPi0,pGamma1,pGamma2,angle);
  
  //Check if Pi0 is in the acceptance of the calorimeters, if aperture angle is small, keep it
  if(IsInPHOS(particle->Phi(), particle->Eta())){
    cellDistance = angle*fPHOSIPDistance;
    if (cellDistance < fPHOSMinDistance){
      new((*plPHOS)[indexPHOS++])       TParticle(*particle) ;
      checkPhoton = kFALSE;
    }
  } 
  else if(IsInEMCAL(particle->Phi(), particle->Eta())){
    cellDistance = angle*fEMCALIPDistance;
    if (cellDistance < fEMCALMinDistance) {
      new((*plEMCAL)[indexEMCAL++])       TParticle(*particle) ;
      checkPhoton = kFALSE;
    }
  } 
  else checkPhoton = kTRUE ;
  
  if (checkPhoton) {
    //Gamma Not overlapped
    TParticle * photon1 = new TParticle(22,1,0,0,0,0,pGamma1.Px(),pGamma1.Py(),
					pGamma1.Pz(),pGamma1.E(),0,0,0,0);    
    if(photon1->Pt() > fNeutralPtCut) {
      
      if(IsInPHOS(photon1->Phi(), photon1->Eta()))
	new((*plPHOS)[indexPHOS++])       TParticle(*photon1) ;
      
      else if(IsInEMCAL(photon1->Phi(), photon1->Eta()))
	new((*plEMCAL)[indexEMCAL++])       TParticle(*photon1) ;
      
    }// photon 1 of pi0 in acceptance
    
    
    TParticle * photon2 = new TParticle(22,1,0,0,0,0,pGamma2.Px(),pGamma2.Py(),
					pGamma2.Pz(),pGamma2.E(),0,0,0,0);
    
    if(photon2->Pt() > fNeutralPtCut) {
      
      if(IsInPHOS(photon2->Phi(), photon2->Eta()))
	new((*plPHOS)[indexPHOS++])       TParticle(*photon2) ;
      
      else if(IsInEMCAL(photon2->Phi(), photon2->Eta()))
	new((*plEMCAL)[indexEMCAL++])       TParticle(*photon2) ;
      
    }// photon 2 of pi0 in acceptance
  }//Not overlapped gamma
}


//_______________________________________________________________
void AliGammaMCReader::InitParameters()
{
  
  //Initialize the parameters of the analysis.
  
  //Fill particle lists when PID is ok
  fEMCALMinDistance    = 10. ;  
  fPHOSMinDistance    = 3.6 ;
  fEMCALIPDistance    = 450. ;//cm  
  fPHOSIPDistance    = 460. ;//cm
  fDecayPi0 = kGeantDecay;
}

//____________________________________________________________________________
void AliGammaMCReader::Pi0Decay(TLorentzVector &p0, TLorentzVector &p1, 
				TLorentzVector &p2, Double_t &angle) {
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
  angle = p1.Angle(p2.Vect());
  //cout<<angle<<endl;
}

//________________________________________________________________
void AliGammaMCReader::Print(const Option_t * opt) const
{
  
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  Info("Print", "%s %s", GetName(), GetTitle() ) ;
  
  printf("IP distance to PHOS         : %f\n", fPHOSIPDistance) ;
  printf("IP distance to EMCAL         : %f\n", fEMCALIPDistance) ;
  printf("Min gamma decay distance in PHOS         : %f\n", fPHOSMinDistance) ;
  printf("Min gamma decay distance in EMCAL         : %f\n", fEMCALMinDistance) ;
  printf("Decay Pi0?          : %d\n", fDecayPi0) ;
  
}
