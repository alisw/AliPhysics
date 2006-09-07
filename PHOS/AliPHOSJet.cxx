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
// Class to calculate jet chararacteristics
//
// This class implements for PHOS a jet finder for PHOS. It depends on a 
// energy seed
// minimum energy, cone radius and movement of the cone.
//*-- Author :  D.Peressounko
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TParticle.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliPHOSJet.h"
//  #include "AliPHOSGetter.h"

ClassImp(AliPHOSJet)
  
//____________________________________________________________________________ 
AliPHOSJet::AliPHOSJet():
  fNpart(0),
  fList(0),
  fConeRad(0),
  fMaxConeMove(0),
  fMinConeMove(0),
  fSumEnergy(0),
  fSumEta(0),
  fSumPhi(0),
  fEnergy(0),
  fEta(0),
  fPhi(0),
  fLEnergy(0),
  fLEta(0),
  fLPhi(0) 
{
  //Initialize members
}

//____________________________________________________________________________ 
AliPHOSJet::AliPHOSJet(const AliPHOSJet & jet) : 
  TObject(jet),
  fNpart(0),
  fList(0),
  fConeRad(0),
  fMaxConeMove(0),
  fMinConeMove(0),
  fSumEnergy(0),
  fSumEta(0),
  fSumPhi(0),
  fEnergy(0),
  fEta(0),
  fPhi(0),
  fLEnergy(0),
  fLEta(0),
  fLPhi(0) 
{
  // copy ctor: no implementation yet
  Fatal("cpy ctor", "not implemented") ;
}


//____________________________________________________________________________ 
AliPHOSJet::~AliPHOSJet(){
  //dtor
  if(fList){
    delete fList ;
    fList = 0 ;
  }
  
}
//____________________________________________________________________________ 
void AliPHOSJet::AddParticle(const TParticle * p, Int_t index){
  //adds particle to jet. Calculates change in jet direction, 
  //due to addition of this particle and if it is smaller, than fMaxDev, 
  //add particle, axcept new direction and return true.
  
  fSumEnergy+=p->Energy() ;
  if(p->Pt()/p->Energy() > 0.001)
    fSumEta+=p->Eta()*p->Energy() ;
  else
    fSumEta+=100.*p->Pz() ; 
  fSumPhi+=p->Phi()*p->Energy() ;    
  
  //check if this a leading particle?
  if(fLEnergy < p->Energy()){
    fLEnergy = p->Energy() ;
    if(p->Pt()/p->Energy() > 0.001)
      fLEta = p->Eta() ;
    else 
      fLEta = 100.*p->Pz()/p->Energy() ;
    fLPhi = p->Phi() ;
  }

  if(index >=0){ //add index to list of indexes
    if(!fList){
      fList = new TArrayI(100) ;
    }
    if(fList->GetSize()<=fNpart+1){
      TArrayI * tmp = new TArrayI(fNpart*2) ;
      tmp->Adopt(fList->GetSize(),fList->GetArray()) ; //note, we should not delete old array!
      fList=tmp ;
    }
    fList->AddAt(index,fNpart++) ;
  }
}
//____________________________________________________________________________ 
void AliPHOSJet::AddDigit(Double_t e, Double_t eta, Double_t phi, Int_t index){
  //adds particle to jet. Calculates change in jet direction, 
  //due to addition of this particle and if it is smaller, than fMaxDev, 
  //add particle, axcept new direction and return true.
  
  fSumEnergy+=e ;
  fSumEta+=eta*e ;
  fSumPhi+=phi*e ;    
  
  //check if this a leading particle?
  if(fLEnergy < e){
    fLEnergy = e;
    fLEta = eta ;
    fLPhi = phi ;
  }

  if(index >=0){ //add index to list of indexes
    if(!fList){
      fList = new TArrayI(100) ;
    }
    if(fList->GetSize()<=fNpart+1){
      TArrayI * tmp = new TArrayI(fNpart*2) ;
      tmp->Adopt(fList->GetSize(),fList->GetArray()) ; //note, we should not delete old array!
      fList=tmp ;
    }
    fList->AddAt(index,fNpart++) ;
  }
}
// //____________________________________________________________________________ 
// Bool_t AliPHOSJet::IsInJet(TParticle * p,Int_t mode)
// {
//   Double_t dEta ;
//   Double_t dPhi ;
//   Double_t energy ;
//   if(!fEnergy){ //Final values not calculated yet, use intermediate
//     if(fSumEnergy==0)
//       return kTRUE ; //First particle
// note p->Eta() causes fpe! 
//     dEta=(p->Eta() - fSumEta/fSumEnergy) ;
//     dPhi=(p->Phi() - fSumPhi/fSumEnergy) ;
//     energy = fSumEnergy ;
//   }
//   else{ //calculated with respect to final 
//     dEta=(p->Eta() - fEta) ;
//     dPhi=(p->Phi() - fPhi) ;
//     energy = fEnergy ;    
//   }
  
//   switch (mode){
//   case 0: 
//     return IsInCone(eta,phi) ; //pure geometrical comparison
//   case 1: 
//     return AcceptConeDeviation(dEta,dPhi,p->Energy() );
//   default:
//     AliError(Form("Unknown mode of cone calculation %d \n",mode ));
//   }
//   return kFALSE ;
//}
//____________________________________________________________________________ 
Bool_t AliPHOSJet::AcceptConeDeviation(const TParticle * p)const
{ //Calculate cone deviation in case of inclusion of the given
  //particle to jet. 

  Double_t tmpEnergy = fSumEnergy + p->Energy() ;  
  Double_t tmpEta = fSumEta ;
  if(p->Pt()/p->Energy() >0.001)
    tmpEta+= p->Eta()*p->Energy() ;
  else
    tmpEta+= 100.*p->Pz() ;
  tmpEta = tmpEta/tmpEnergy / (fSumEta/fSumEnergy) - 1 ;
  Double_t tmpPhi = fSumPhi + p->Phi()*p->Energy() ;
  tmpPhi = tmpPhi/tmpEnergy /(fSumPhi/fSumEnergy) - 1 ;
  Double_t dev = TMath::Sqrt(tmpEta*tmpEta + tmpPhi*tmpPhi) ;
  if(dev < fMaxConeMove)
    return kTRUE ;
  else
    return kFALSE ;
}
//____________________________________________________________________________ 
Bool_t AliPHOSJet::AcceptConeDeviation(Double_t e, Double_t eta, Double_t phi)const
{ //Calculate cone deviation in case of inclusion of the given
  //particle to jet. 

  Double_t tmpEnergy = fSumEnergy + e ;
  Double_t tmpEta = fSumEta + eta*e ;
  tmpEta = tmpEta/tmpEnergy / (fSumEta/fSumEnergy) ;
  Double_t tmpPhi = fSumPhi + phi*e ;
  tmpPhi = tmpPhi/tmpEnergy /(fSumPhi/fSumEnergy) ;
  Double_t dev = TMath::Sqrt(tmpEta*tmpEta + tmpPhi*tmpPhi) ;
  if(dev<fMaxConeMove && dev > fMinConeMove)
    return kTRUE ;
  else
    return kFALSE ;
}
//____________________________________________________________________________ 
Bool_t AliPHOSJet::IsInCone(const TParticle * p)const
{
  //Say if  particle is inside the defined cone
  Double_t dEta ;
  Double_t dPhi ;
  if(!fEnergy){ //Final values not calculated yet, use intermediate
    if(fSumEnergy==0)
      return kTRUE ; //First particle    
    if(p->Pt()/p->Energy() > 0.001)
      dEta=(p->Eta() - fSumEta/fSumEnergy) ;
    else
      dEta=(100.*p->Pz()/p->Energy() - fSumEta/fSumEnergy) ;
    dPhi=(p->Phi() - fSumPhi/fSumEnergy) ;
  }
  else{ //calculated with respect to final 
    if(p->Pt()/p->Energy() > 0.001)
      dEta=(p->Eta() - fEta) ;
    else
      dEta=(100.*p->Pz()/p->Energy() - fEta) ;
    dPhi=(p->Phi() - fPhi) ;
  }
  if(TMath::Sqrt(dEta*dEta + dPhi*dPhi) < fConeRad)
    return kTRUE ;
  else
    return kFALSE ;
}
//____________________________________________________________________________ 
Bool_t AliPHOSJet::IsInCone(Double_t eta, Double_t phi)const
{
  //Says if particle is inside the defined cone
  Double_t dEta ;
  Double_t dPhi ;
  if(!fEnergy){ //Final values not calculated yet, use intermediate
    if(fSumEnergy==0)
      return kTRUE ; //First particle    
    dEta=(eta - fSumEta/fSumEnergy) ;
    dPhi=(phi - fSumPhi/fSumEnergy) ;
  }
  else{ //calculated with respect to final 
    dEta=(eta - fEta) ;
    dPhi=(phi - fPhi) ;
  }
  if(TMath::Sqrt(dEta*dEta + dPhi*dPhi) < fConeRad)
    return kTRUE ;
  else
    return kFALSE ;
}
//____________________________________________________________________________ 
Double_t AliPHOSJet::DistanceToJet(const TParticle *p)const{
  //Calculate radius 
  Double_t dEta ;
  Double_t dPhi ;
  if(!fEnergy){ //Final values not calculated yet, use intermediate
    if(fSumEnergy==0)
      return kTRUE ; //First particle    
    if(p->Pt()/p->Energy() > 0.001)
      dEta=(p->Eta() - fSumEta/fSumEnergy) ;
    else
      dEta=(100.*p->Pz()/p->Energy() - fSumEta/fSumEnergy) ;      
    dPhi=(p->Phi() - fSumPhi/fSumEnergy) ;
  }
  else{ //calculated with respect to final 
    if(p->Pt()/p->Energy() > 0.001)
      dEta=(p->Eta() - fEta) ;
    else
      dEta=(100*p->Pz()/p->Energy() - fEta) ;    
    dPhi=(p->Phi() - fPhi) ;
  }
  return TMath::Sqrt(dEta*dEta + dPhi*dPhi) ;

}
//____________________________________________________________________________ 
void AliPHOSJet::CalculateAll(void){
  //Calculate all jet parameters
  if(fSumEnergy==0)
    return  ; //Nothing to calculate    
  
  fEta = fSumEta/fSumEnergy ;
  fPhi = fSumPhi/fSumEnergy ;
  fEnergy = fSumEnergy ;
  
  fSumEnergy = 0. ;
  fSumEta = 0. ;
  fSumPhi = 0. ;
}
//____________________________________________________________________________ 
void AliPHOSJet::Print(const Option_t *) const {
  //Print jet parameters
  printf("-------------- AliPHOSJet ------------\n") ;
  printf(" Energy............. %f \n",fEnergy) ;
  printf(" Eta................ %f \n",fEta ) ;
  printf(" Phi................ %f \n",fPhi ) ;
  printf(" Leading Energy..... %f \n",fLEnergy) ;
  printf(" Leading Eta........ %f \n",fLEta) ;
  printf(" Leading Phi........ %f \n",fLPhi) ;
  printf(" N particles in jet  %d \n",fNpart) ;
  printf("----------------------------------\n") ;
}



