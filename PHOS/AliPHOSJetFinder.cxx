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

//_________________________________________________________________________
// C++ version of UA2 and/or Lund jet finding algorithm
// UA1 jet algorithm from LUND JETSET (LUCELL)
// Find jets at the level of no detector and Digits.
// Needs modifications. 
//*-- Author : D.Peressounko after UA1 coll. etc
//////////////////////////////////////////////////////////////////////////////

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.8  2005/05/28 14:19:04  schutz
 * Compilation warnings fixed by T.P.
 *
 */

// --- ROOT system ---
#include "TClonesArray.h"
//      #include "TIter.h"
#include "TParticle.h"
#include "TTask.h"
// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliPHOSJet.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSDigit.h"
#include "AliPHOSJetFinder.h"
#include "AliPHOSDigitizer.h"

ClassImp(AliPHOSJetFinder)


//____________________________________________________________________________ 
AliPHOSJetFinder::AliPHOSJetFinder():
  TNamed("AliPHOSJetFinder",""),
  fNJets(0),
  fStatusCode(-999),
  fMode(0),
  fConeRad(1.),
  fMaxConeMove(0.15),
  fMinConeMove(0.05),
  fEtSeed(4.),
  fEtMin(5.),
  fPrecBg(0.00035),
  fSimGain(0.),
  fSimPedestal(0.),
  fParticles(0),
  fJets(0)
{
  //Initialize jet parameters
}

//____________________________________________________________________________ 
AliPHOSJetFinder::AliPHOSJetFinder(const AliPHOSJetFinder & jet) : 
  TNamed(jet),
  fNJets(0),
  fStatusCode(-999),
  fMode(0),
  fConeRad(1.),
  fMaxConeMove(0.15),
  fMinConeMove(0.05),
  fEtSeed(4.),
  fEtMin(5.),
  fPrecBg(0.00035),
  fSimGain(0.),
  fSimPedestal(0.),
  fParticles(0),
  fJets(0)
{
  // copy ctor: no implementation yet
  Fatal("cpy ctor", "not implemented");
}

//____________________________________________________________________________ 
  AliPHOSJetFinder::~AliPHOSJetFinder()
{
  //dtor
  if(fParticles){
    delete fParticles ;
    fParticles = 0 ;
  }
  if(fJets){
    delete fJets ;
    fJets = 0 ;
  } 
}

//____________________________________________________________________________ 
void  AliPHOSJetFinder::FindJetsFromParticles(const TClonesArray * plist,TObjArray * jetslist) 
{
  //Find jets in the case without detector.
  TIter next(plist) ;

  TIter nextJet(jetslist) ;

  fNJets = 0 ;
  TParticle * p ;
  AliPHOSJet * jet ;
  //In this cicle we find number of jets and define approx. their directions
  //note, that we do not really add particles to jet (index =-1)
  while((p=static_cast<TParticle*>(next()))){
    if(fStatusCode==-999 || p->GetStatusCode()==fStatusCode){
    if(p->Energy() >= fEtSeed){ //Energetic enough
      //cout << "p " << p->Energy() << endl ;
      //cout << "Status "<<fStatusCode<<" "<<p->GetName()<<" " << p->Energy() << " "<<p->Eta()<< " "<<p->Phi()<<endl ;
      Bool_t startnew = kTRUE ;
      //Do not start new jet if part of older jet
      nextJet.Reset() ;
      while((jet=static_cast<AliPHOSJet*>(nextJet()))){
	//jet->Print() ;
	if(jet->AcceptConeDeviation(p)){
	  startnew = kFALSE ;
	  //cout << "false" << endl ;
	  break ;
	}
      }
      if(startnew){
	//cout << "new " << endl ;
	jet = new AliPHOSJet() ;
	jetslist->Add(jet) ;
	//	jet = static_cast<AliPHOSJet*>(jetslist->Last()) ;
	jet->SetConeRadius(fConeRad) ;
	jet->SetMaxConeMove(fMaxConeMove) ;
	//	jet->SetMinConeMove(fMinConeMove) ;
	jet->AddParticle(p,-1) ;
	fNJets++;
      }
    }
    while((jet=static_cast<AliPHOSJet*>(nextJet()))){
      if(jet->AcceptConeDeviation(p))
	jet->AddParticle(p,-1) ; //Just recalculate direction of jet
    }
    }
  }
    
  //now calculate directions of jets using collected information
  nextJet.Reset() ;
  while((jet=static_cast<AliPHOSJet*>(nextJet()))){
    jet->CalculateAll() ;
    if(jet->Energy() < fEtMin){
      jetslist->Remove(jet) ;
      delete jet ;
    }
  }

  jetslist->Compress() ;
  //And finally, really add particles to jets
  for(Int_t iPart=0; iPart<plist->GetEntries();iPart++){
    p=static_cast<TParticle*>(plist->At(iPart)) ;
    if(fStatusCode == -999 || p->GetStatusCode()==fStatusCode){
    Double_t dist = 999999. ; //big distance
    Int_t iJet = -1 ;
    for(Int_t i=0; i<jetslist->GetEntriesFast();i++){
      jet=static_cast<AliPHOSJet*>(jetslist->At(i)) ;
      if(jet->IsInCone(p)){
	Double_t cdist = jet->DistanceToJet(p); 
	if(cdist < dist){
	  dist = cdist ;
	  iJet = i ;
	}
      }
    }
    if(iJet>-1)
      (static_cast<AliPHOSJet*>(jetslist->At(iJet)))->AddParticle(p,iPart); //assign particle to closest jet
  }
  }
  
  //Calculate jet parameters 
  nextJet.Reset() ;
  while((jet=static_cast<AliPHOSJet*>(nextJet()))){
    jet->CalculateAll() ;
  }

}
//____________________________________________________________________________ 
void AliPHOSJetFinder::FindJetsFromDigits(const TClonesArray * digits, TObjArray * jets){
  //Find jets in the case witht detector at the level of digits.
  if(digits->GetEntries()==0){
    AliError(Form("No entries in digits list \n")) ;
    return ;
  }


  TClonesArray * copyDigits = new TClonesArray(*digits) ;

  //Remove CPV digits if any
  AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance("GPS2","") ;
  Int_t iDigit ;
  AliPHOSDigit * digit ;
  for(iDigit=copyDigits->GetEntries()-1;iDigit>=0;iDigit--){
    digit=static_cast<AliPHOSDigit *>(copyDigits->At(iDigit)) ;
    if(!geom->IsInEMC(digit->GetId()))
      copyDigits->RemoveAt(iDigit) ;
    else
      break ;
  }
  copyDigits->Compress() ;

  Double_t totalEnergy = 0 ;
  Float_t * energy = new Float_t[copyDigits->GetEntries()] ;
  //calculate average energy of digits
  //fill array of energies
  for(iDigit=0;iDigit<copyDigits->GetEntries();iDigit++){
    digit=static_cast<AliPHOSDigit *>(copyDigits->At(iDigit)) ;
    energy[iDigit] = Calibrate(digit) ;
    totalEnergy+=energy[iDigit] ;
  }
  
  //Sort digits in decreasing energy.
  Int_t * index = new Int_t[copyDigits->GetEntries()] ;
  TMath::Sort(copyDigits->GetEntries(),energy,index) ;
  
  Double_t eAverage = totalEnergy/copyDigits->GetEntries()  ;
  //remove digits below average energy
  for(iDigit=copyDigits->GetEntries()-1;iDigit>=0;iDigit--){
    digit=static_cast<AliPHOSDigit *>(copyDigits->At(index[iDigit])) ;
    if(energy[index[iDigit]] < eAverage)
      copyDigits->RemoveAt(iDigit) ;
    else
      break ;
  }
  
  
  AliPHOSJet * jet ;
  Int_t iIter = 0 ;
  while(iIter < 10){//less than 10 iterations
    
    //while digits above seed
    for(Int_t ind=0;ind<copyDigits->GetEntriesFast();ind++){
      digit=static_cast<AliPHOSDigit*>(copyDigits->At(index[ind])) ;
      if(energy[index[ind]] > fEtSeed && digit){ //start new jet      
 	jet = new AliPHOSJet() ;
	Double_t e,eta,phi ;
	CalculateEEtaPhi(digit,e,eta,phi) ;
	jet->AddDigit(e,eta,phi,-1) ;  
	//loop over left digits
	for(iDigit = 0 ; iDigit < copyDigits->GetEntries() ; iDigit++){
	  if(iDigit!= ind){ //first digit already in jet
	    digit = static_cast<AliPHOSDigit *>(copyDigits->At(iDigit));
	    CalculateEEtaPhi(digit,e,eta,phi) ;	    
	    if(jet->IsInCone(eta,phi) && //is cell in cone
	       jet->AcceptConeDeviation(e,eta,phi)){//if cone does not move too much	      
	      jet->AddDigit(e,eta,phi,-1) ;  //accept new direction
	    }
	  }
	}//end of loop over cells
	
	//accept  all anused cells incide cone
	//note, that digits might be returned as anused later
	for(Int_t icell = 0 ; icell < copyDigits->GetEntries() ; icell++){
	  digit = static_cast<AliPHOSDigit *>(copyDigits->At(icell));
	  if(jet->IsInCone(eta,phi)){ //is cell in cone
	    CalculateEEtaPhi(digit,e,eta,phi) ;
	    jet->AddDigit(e,eta,phi,digit->GetIndexInList()) ; 
	  }
	}
	
	//Accept Jet with Et > Et_min and remove all belonging digits
	if(jet->Energy()/TMath::CosH(jet->Eta()) > fEtMin){
	  Int_t nIndxs ;
	  const Int_t * indxs = jet->Indexs(nIndxs) ;
	  for(Int_t i=0;i<nIndxs;i++){
	    copyDigits->RemoveAt(indxs[i]) ;
	  }
	  jet->CalculateAll() ;
	  jets->AddAt(jet,fNJets++);
	}
	else{ //remove jet and do not touch digits
	  delete jet ;
	}
      }
      else{ 
	if(energy[index[ind]] < fEtSeed){ // no more digits above threshold left, return from loop
	  break ;
	}
      }
      
      iIter++ ;
      //calculate new energy of backrgound
      Double_t oldTotalEnergy = totalEnergy ;
      totalEnergy = 0 ;
      for(Int_t i=0 ; i<copyDigits->GetEntriesFast() ; i++){
	digit=static_cast<AliPHOSDigit*>(copyDigits->At(index[ind])) ;
	if(digit)
	  totalEnergy+=energy[i] ;
      }
      if(!fMode || ((oldTotalEnergy != 0) && 
		    (TMath::Abs(oldTotalEnergy - totalEnergy)/oldTotalEnergy < fPrecBg)))
	break ;	
    }
  }
  delete [] energy;
  delete [] index;
  copyDigits->Delete() ;
  
}
//____________________________________________________________________________ 
Double_t AliPHOSJetFinder::Calibrate(const AliPHOSDigit * digit){ 
  //Both simulated and raw digits are already calibrated
  
  return digit->GetEnergy() ;        
  
  //  }  
}
//____________________________________________________________________________ 
void AliPHOSJetFinder::CalculateEEtaPhi(const AliPHOSDigit * d,Double_t &e, Double_t &eta, Double_t &phi){
  //Calculate direction of the jet
  e=Calibrate(d) ;
  AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance("GPS2","") ;
  TVector3 pos ;
  geom->RelPosInAlice(d->GetId(), pos) ;
  eta = pos.Eta() ;
  phi = pos.Phi() ;
}
//____________________________________________________________________________ 
void AliPHOSJetFinder::Print(const Option_t *) const {	
  //Print parameters of the found jet
  printf("\n --------------- AliPHOSJetFinder --------------- \n") ;
  printf(" Jets found .........%d \n",fNJets) ;
  printf(" Seed energy cut ....%f \n",fEtSeed) ;
  printf(" Cone radius ........%f \n",fConeRad) ;
  printf(" Minimal cone move ..%f \n",fMinConeMove) ;
  printf(" Maximal cone move ..%f \n",fMaxConeMove) ;
  printf("------------------------------------------------- \n") ;
}
