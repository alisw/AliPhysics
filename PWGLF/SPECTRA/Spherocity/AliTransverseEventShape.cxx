/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/*
Class AliTransverseEventShape  
It allow to calculate Sphericity and Spherocity variables

*/
//           Please report any Bugs to:
//*************************************************************************
//   author: Antonio Ortiz Velasquez  ( antonio.ortiz@nucleares.unam.mx )
//           Eleazar Cuautle Flores   ( ecuautle@nucleares.unam.mx )
/**************************************************************************


  First version, June 4 2015
***************************************************************************/


#include "AliStack.h"

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVVertex.h"
#include "AliLog.h"
#include "AliAODVertex.h"
#include "AliVTrack.h"
#include "AliVEvent.h"
#include <TMatrixDSym.h>
#include <TMath.h>
#include <TParticlePDG.h>
#include <TParticle.h>
#include "AliESDUtils.h"
#include "AliESDtrackCuts.h"
#include "AliTransverseEventShape.h"
#include <TFile.h>
#include "AliAODHeader.h"
// STL includes
#include <iostream>
using namespace std;


ClassImp(AliTransverseEventShape)

	//______________________________________________________________________
AliTransverseEventShape::AliTransverseEventShape():TObject(),
  fUseHybrid(0),
  fTrackFilterHybrid1(0),
  fTrackFilterHybrid2(0),
  fTrackFilterESA(0),
  fMinMultESA(0),
  fSizeStepESA(0),
  fIsAbsEtaESA(0),
  fEtaMaxCutESA(0),
  fEtaMinCutESA(0),
  fPtMaxCutESA(0),
  fPtMinCutESA(0),
  fRunNumber(0),
  fheta(0),
  fhphi(0),
  fhpt(0),
  fhetaMC(0),
  fhphiMC(0),
  fhptMC(0),
  fAverageAmplitudes(0)
{
  // Default contructor
}


//_____________________________________________________________________________
AliTransverseEventShape::AliTransverseEventShape(const AliTransverseEventShape &c) : 
  TObject(c),
  fUseHybrid(0),
  fTrackFilterHybrid1(0),
  fTrackFilterHybrid2(0),
  fTrackFilterESA(0),
  fMinMultESA(0),
  fSizeStepESA(0),
  fIsAbsEtaESA(0),
  fEtaMaxCutESA(0),
  fEtaMinCutESA(0),
  fPtMaxCutESA(0),
  fPtMinCutESA(0),
  fRunNumber(0),
  fheta(0),
  fhphi(0),
  fhpt(0),
  fhetaMC(0),
  fhphiMC(0),
  fhptMC(0),
  fAverageAmplitudes(0)
  
{
  //
  // copy constructor - untested
  //
  ((AliTransverseEventShape &) c).Copy(*this);
}//_____________________________________________________________________________
AliTransverseEventShape &AliTransverseEventShape::operator=(const AliTransverseEventShape &c)
{
  //
  // Assignment operator - untested
  //
  
  if (this != &c) ((AliTransverseEventShape &) c).Copy(*this);
  return *this;
}
//__________________________________________________________________________
void AliTransverseEventShape::Init(){
  fheta = new TH1D("fheta","ESD;#eta; entries",30,-1.5,1.5);
  fhphi = new TH1D("fhphi","ESD;#phi (rad); entries", 64, 0, 2*TMath::Pi());
  fhpt  = new TH1D("fhpt","ESD;#it{p}_{T} (GeV/#it{c}); entries", 1000, 0, 100);
  
  fhetaMC = new TH1D("fhetaMC","MC;#eta; entries",30,-1.5,1.5);
  fhphiMC = new TH1D("fhphiMC","MC;#phi (rad); entries", 64, 0, 2*TMath::Pi());
  fhptMC  = new TH1D("fhptMC","MC;#it{p}_{T} (GeV/#it{c}); entries", 1000, 0, 100);
    
} 


//_____________________________________________________________________________
Float_t AliTransverseEventShape::MinVal( Float_t A, Float_t B ) {
  if( A < B ) {
    return A;
  }
  else {
    return B;
  }
}

//______________________________________________________________________
Float_t AliTransverseEventShape::GetEventShape( AliVEvent *event, TString lMethod, Bool_t fillHist )
{
  //  cout<<"executing GetEventShape"<<endl;
  
  Int_t lRequestedRunNumber = event->GetRunNumber();
  Double_t lreturnval = 1.0;
  
  if ( lMethod == "SO" ) lreturnval = GetSpherocity(event, fillHist);
  if ( lMethod == "ST" ) lreturnval = GetSphericity(event, fillHist);
  
  return lreturnval;
}
//______________________________________________________________________
Float_t AliTransverseEventShape::GetEventShapeTrue( AliStack *event, TString lMethod, Bool_t fillHist )
{
  
  //  cout<<"executing GetEventShapeTrue"<<endl;
  
  //Int_t lRequestedRunNumber = event->GetRunNumber();
  Double_t lreturnval = 1.0;
  
  if ( lMethod == "SO" ) lreturnval = GetSpherocityMC(event, fillHist);
  if ( lMethod == "ST" ) lreturnval = GetSphericityMC(event, fillHist);
  
  return lreturnval;
}
//_____________________________________________________________________
Float_t AliTransverseEventShape::GetSpherocity( AliVEvent *event, Bool_t fillHist )
{
  
  //  cout<<"Executig  GetSpherocity   fEtaMaxCutESA="<<fEtaMaxCutESA<<endl;
  //  cout<<"fEtaMinCutESA="<<fEtaMinCutESA<<endl;
  Float_t spherocity = -10.0;
  
  if (event->InheritsFrom("AliESDEvent")) {
    AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(event);
    if (!esdevent) return kFALSE;
    
    Int_t nTracks = esdevent->GetNumberOfTracks();
    Int_t nTracksForESA = 0;
    
    
    for(Int_t iT = 0; iT < nTracks; iT++) {            
      AliESDtrack* Track = 0;
      Track = esdevent->GetTrack(iT);
      if(!Track)
	continue;
            

      if(fIsAbsEtaESA){  //cuts in pseudorapidity
	if( TMath::Abs(Track->Eta()) > fEtaMaxCutESA || TMath::Abs(Track->Eta()) < fEtaMinCutESA )					
	  continue;
      }
      else{
	
	if( Track->Eta() > fEtaMaxCutESA || Track->Eta() < fEtaMinCutESA )  
	  continue;
      }
      //cuts in pt
      if( Track->Pt() > fPtMaxCutESA || Track->Pt() <  fPtMinCutESA )
	continue;
      
      //quality cuts
      if(!fUseHybrid){//golden track cuts
	if(!fTrackFilterESA->IsSelected(Track))
	  continue;
      }
      else{//hybrid track cuts
	Bool_t cutset1 = kFALSE;
	Bool_t cutset2 = kFALSE;
	
	cutset1 = fTrackFilterHybrid1->IsSelected(Track);
	cutset2 = (!fTrackFilterESA->IsSelected(Track)) && fTrackFilterHybrid2->IsSelected(Track);
	if(!(cutset1 || cutset2))
	  continue;
      }
      
      if(fillHist){
	fheta->Fill(Track->Eta());
	fhphi->Fill(Track->Phi()); 
	fhpt->Fill(Track->Pt());
      }
      nTracksForESA++;
      
    } //close first loop on nTracks
    
    cout<<"nTracksForESA="<<nTracksForESA<<endl;
    
    if( nTracksForESA < fMinMultESA )
      return -0.5;
    else{
      
      
      Double_t *pxA=new Double_t[nTracksForESA];
      Double_t *pyA=new Double_t[nTracksForESA];
      Double_t sumapt=0;
      Int_t counter=0;
      
      for(Int_t iT = 0; iT < nTracks; iT++) { // second loop over nTracks
	AliESDtrack* Track = 0;
	Track = esdevent->GetTrack(iT);
	if(!Track)
	  continue;
	
	if(fIsAbsEtaESA){  	//cuts in pseudorapidity
	  if( TMath::Abs(Track->Eta()) > fEtaMaxCutESA || TMath::Abs(Track->Eta()) < fEtaMinCutESA )					
	    continue;
	}
	else{
	  if( Track->Eta() > fEtaMaxCutESA || Track->Eta() < fEtaMinCutESA )  
	    continue;
	}
	//cuts in pt
	if( Track->Pt() > fPtMaxCutESA || Track->Pt() <  fPtMinCutESA )
	  continue;
	
	//quality cuts
	if(!fUseHybrid){//golden cuts
	  if(!fTrackFilterESA->IsSelected(Track))
	    continue;
	}
	else{//hybrid track cuts
	  Bool_t cutset1 = kFALSE;
	  Bool_t cutset2 = kFALSE;
	  
	  cutset1 = fTrackFilterHybrid1->IsSelected(Track);
	  cutset2 = (!fTrackFilterESA->IsSelected(Track)) && fTrackFilterHybrid2->IsSelected(Track);
	  if(!(cutset1 || cutset2))
	    continue;
	}
	
		
	pxA[counter] = Track->Pt() * TMath::Cos( Track->Phi() );
	pyA[counter] = Track->Pt() * TMath::Sin( Track->Phi() );
	sumapt += Track->Pt();
	counter++;
	
      }//close second loop
      
      Double_t pFull = 0;
      Double_t Spherocity = 2;
      //Getting thrust
      for(Int_t i = 0; i < 360/(fSizeStepESA); ++i){
	Double_t numerador = 0;
	Double_t phiparam  = 0;
	Double_t nx = 0;
	Double_t ny = 0;
	phiparam=( (TMath::Pi()) * i * fSizeStepESA ) / 180; // parametrization of the angle
	nx = TMath::Cos(phiparam);            // x component of an unitary vector n
	ny = TMath::Sin(phiparam);            // y component of an unitary vector n
	for(Int_t i1 = 0; i1 < nTracksForESA; ++i1){
	  numerador += TMath::Abs(ny * pxA[i1] - nx * pyA[i1]);//product between p  proyection in XY plane and the unitary vector
	}
	pFull=TMath::Power( (numerador / sumapt),2 );
	if(pFull < Spherocity)//maximization of pFull
	  {
	    Spherocity = pFull;
	  }
      }
      
      spherocity=((Spherocity)*TMath::Pi()*TMath::Pi())/4.0;
      
      if(pxA){// clean up array memory used for TMath::Sort
	delete[] pxA;
	pxA=0;
      }
      if(pyA){// clean up array memory used for TMath::Sort
	delete[] pyA;
	pyA=0;
      }
      
      
    }//close case when multiplicity is ok!
        
  }//close ESD case
  //Redo equivalent test
  else if (event->InheritsFrom("AliAODEvent")) {
    AliAODEvent *aodevent = dynamic_cast<AliAODEvent *>(event);
    if (!aodevent) return kFALSE;
    
    cout<<"AOD ana"<<endl;
    
  }
  
  return spherocity;
  
}
//____________________________________________________________________
Float_t AliTransverseEventShape::GetSphericity( AliVEvent *event, Bool_t fillHist  )
{
  
  //  cout<<"Executig  GetSpherIcity   fEtaMaxCutESA="<<fEtaMaxCutESA<<endl;
  //  cout<<"fEtaMinCutESA="<<fEtaMinCutESA<<endl;
  
  
  Float_t sphericity = -10.0;
  
  if (event->InheritsFrom("AliESDEvent")) {
    
    AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(event);
    if (!esdevent) return kFALSE;
    
    Int_t nTracks = esdevent->GetNumberOfTracks();
    Double_t s00=0;
    Double_t s01=0;
    Double_t s11=0;
    Double_t totalpt=0;
    Int_t nTracksForESA = 0;
    
    for(Int_t iT = 0; iT < nTracks; iT++) {
      
      AliESDtrack* Track = 0;
      Track = esdevent->GetTrack(iT);
      if(!Track)
	continue;
         

      if(fIsAbsEtaESA){       //cuts in pseudorapidity
	if( TMath::Abs(Track->Eta()) > fEtaMaxCutESA || TMath::Abs(Track->Eta()) < fEtaMinCutESA )
	  continue;
      }
      else{
	if( Track->Eta() > fEtaMaxCutESA || Track->Eta() < fEtaMinCutESA )  
	  continue;
      }
      //cuts in pt
      if( Track->Pt() > fPtMaxCutESA || Track->Pt() <  fPtMinCutESA )
	continue;
      
      //quality cuts
      if(!fUseHybrid){//golden track cuts
	if(!fTrackFilterESA->IsSelected(Track))
	  continue;
      }
      else{//hybrid track cuts
	Bool_t cutset1 = kFALSE;
	Bool_t cutset2 = kFALSE;
	
	cutset1 = fTrackFilterHybrid1->IsSelected(Track);
	cutset2 = (!fTrackFilterESA->IsSelected(Track)) && fTrackFilterHybrid2->IsSelected(Track);
	if(!(cutset1 || cutset2))
	  continue;
      }
                  
      Double_t px = Track->Pt() * TMath::Cos( Track->Phi() );
      Double_t py = Track->Pt() * TMath::Sin( Track->Phi() );
      Double_t pt = Track->Pt();
      
      s00 += (px * px)/pt;
      s01 += (py * px)/pt;
      s11 += (py * py)/pt;
      totalpt += pt;
      
      if(fillHist){
	fheta->Fill(Track->Eta());
	fhphi->Fill(Track->Phi());
	fhpt->Fill(Track->Pt());
      }      
      nTracksForESA++;      
    }//close first loop
        
    if( nTracksForESA < fMinMultESA )
      return -0.5;
    else{
      
      Double_t S00=s00/totalpt;
      Double_t S01=s01/totalpt;
      Double_t S11=s11/totalpt;
      
      Float_t lambda1=((S00+S11)+TMath::Sqrt((S00+S11)*(S00+S11)-4*(S00*S11-S01*S01)))/2;
      Float_t lambda2=((S00+S11)-TMath::Sqrt((S00+S11)*(S00+S11)-4*(S00*S11-S01*S01)))/2;
      if((lambda2==0)&&(lambda1==0))
	sphericity=0;
      if(lambda1+lambda2!=0)
	sphericity=2*TMath::Min( lambda1,lambda2 )/( lambda1+lambda2 );
      
    }//close case when multiplicity is ok!
        
  }//close ESD case
  else if (event->InheritsFrom("AliAODEvent")) {
    AliAODEvent *aodevent = dynamic_cast<AliAODEvent *>(event);
    if (!aodevent) return kFALSE;
    
    cout<<"AOD ana"<<endl;
    
  }
    
  return sphericity;
  
}
//____________________________________________________________________
Float_t AliTransverseEventShape::GetSphericityMC( AliStack *event, Bool_t fillHist )
{
  
  Float_t sphericity = -10.0;
  
  Int_t nTracks = event->GetNtrack();
  
  Double_t s00=0;
  Double_t s01=0;
  Double_t s11=0;
  Double_t totalpt=0;
  Int_t nTracksForESA = 0;
  
  for(Int_t iT = 0; iT < nTracks; iT++) {       
    //Cuts
    if(!(event->IsPhysicalPrimary(iT)))
      continue;
    
    TParticle* Track = event->Particle(iT);
    if(!Track)
      continue;
    
    TParticlePDG* pdgPart = Track->GetPDG();
    Double_t chargeMC = pdgPart->Charge();
    
    if( TMath::Abs(chargeMC) > 0.1 )
      continue;
    

    if(fIsAbsEtaESA){       //cuts in pseudorapidity
      if( TMath::Abs(Track->Eta()) > fEtaMaxCutESA || TMath::Abs(Track->Eta()) < fEtaMinCutESA )
	continue;
    }
    else{
      if( Track->Eta() > fEtaMaxCutESA || Track->Eta() < fEtaMinCutESA )  
	continue;
    }
    //cuts in pt
    if( Track->Pt() > fPtMaxCutESA || Track->Pt() <  fPtMinCutESA )
      continue;
    
    Double_t px = Track->Px();
    Double_t py = Track->Py();
    Double_t pt = Track->Pt();
    
    s00 += (px * px)/pt;
    s01 += (py * px)/pt;
    s11 += (py * py)/pt;
    totalpt += pt;
    
    if(fillHist){
      fhetaMC->Fill(Track->Eta());
      fhphiMC->Fill(Track->Phi());
      fhptMC->Fill(Track->Pt());
    }
    
    nTracksForESA++;
    
  }//close first loop
  
  
  if( nTracksForESA < fMinMultESA )
    return -0.5;
  else{
    
    Double_t S00=s00/totalpt;
    Double_t S01=s01/totalpt;
    Double_t S11=s11/totalpt;
    
    Float_t lambda1=((S00+S11)+TMath::Sqrt((S00+S11)*(S00+S11)-4*(S00*S11-S01*S01)))/2;
    Float_t lambda2=((S00+S11)-TMath::Sqrt((S00+S11)*(S00+S11)-4*(S00*S11-S01*S01)))/2;
    if((lambda2==0)&&(lambda1==0))
      sphericity=0;
    if(lambda1+lambda2!=0)
      sphericity=2*TMath::Min( lambda1,lambda2 )/( lambda1+lambda2 );
    
  }//close case when multiplicity is ok!
    
  return sphericity;
  
}
//_____________________________________________________________________
Float_t AliTransverseEventShape::GetSpherocityMC( AliStack *event, Bool_t fillHist )
{
  
  //  cout<<"Executig  GetSpherocityMC   fEtaMaxCutESA="<<fEtaMaxCutESA<<endl;
  //  cout<<"fEtaMinCutESA="<<fEtaMinCutESA<<endl;
  
  Float_t spherocity = -10.0;
  Int_t nTracks = event->GetNtrack();
  Int_t nTracksForESA = 0;
  
  for(Int_t iT = 0; iT < nTracks; iT++) {  // first loop over nTracks    
    //Cuts
    if(!(event->IsPhysicalPrimary(iT)))
      continue;
    
    TParticle* Track = event->Particle(iT);
    if(!Track)
      continue;
    
    TParticlePDG* pdgPart = Track->GetPDG();
    Double_t chargeMC = pdgPart->Charge();
    
    if( TMath::Abs(chargeMC) < 0.1 )
      continue;        

    if(fIsAbsEtaESA){       //cuts in pseudorapidity
      if( TMath::Abs(Track->Eta()) > fEtaMaxCutESA || TMath::Abs(Track->Eta()) < fEtaMinCutESA ) 
	continue;
    }
    else{
      if( Track->Eta() > fEtaMaxCutESA || Track->Eta() < fEtaMinCutESA )  
	continue;
    }
    //cuts in pt
    if( Track->Pt() > fPtMaxCutESA || Track->Pt() <  fPtMinCutESA )
      continue;
        
    if(fillHist){
      fhetaMC->Fill(Track->Eta());
      fhphiMC->Fill(Track->Phi()); 
      fhptMC->Fill(Track->Pt());
    }
    
    nTracksForESA++;
    
  }//close first loop
  
  
  if( nTracksForESA < fMinMultESA )
    return -0.5;
  else{    
    Double_t *pxA=new Double_t[nTracksForESA];
    Double_t *pyA=new Double_t[nTracksForESA];
    Double_t sumapt=0;
    Int_t counter=0;
    
    for(Int_t iT = 0; iT < nTracks; iT++) {
      
      //Cuts
      if(!(event->IsPhysicalPrimary(iT)))
	continue;
      
      TParticle* Track = event->Particle(iT);
      if(!Track)
	continue;
      
      TParticlePDG* pdgPart = Track->GetPDG();
      Double_t chargeMC = pdgPart->Charge();
      
      if( TMath::Abs(chargeMC) < 0.1 )
	continue;

      if(fIsAbsEtaESA){         //cuts in pseudorapidity
	if( TMath::Abs(Track->Eta()) > fEtaMaxCutESA || TMath::Abs(Track->Eta()) < fEtaMinCutESA ) 
	  continue;
      }
      else{
	if( Track->Eta() > fEtaMaxCutESA || Track->Eta() < fEtaMinCutESA )  
	  continue;
      }
      //cuts in pt
      if( Track->Pt() > fPtMaxCutESA || Track->Pt() <  fPtMinCutESA )
	continue;
        
      pxA[counter] = Track->Px();
      pyA[counter] = Track->Py();
      sumapt += Track->Pt();
      counter++;
      
    }//close second loop
    
    Double_t pFull = 0;
    Double_t Spherocity = 2;
    //Getting thrust
    for(Int_t i = 0; i < 360/(fSizeStepESA); ++i){
      Double_t numerador = 0;
      Double_t phiparam  = 0;
      Double_t nx = 0;
      Double_t ny = 0;
      phiparam=( (TMath::Pi()) * i * fSizeStepESA ) / 180; // parametrization of the angle
      nx = TMath::Cos(phiparam);            // x component of an unitary vector n
      ny = TMath::Sin(phiparam);            // y component of an unitary vector n
      for(Int_t i1 = 0; i1 < nTracksForESA; ++i1){
	numerador += TMath::Abs(ny * pxA[i1] - nx * pyA[i1]);//product between p proyection in XY plane and the unitary vector
      }
      pFull=TMath::Power( (numerador / sumapt),2 );
      if(pFull < Spherocity)//maximization of pFull
	{
	  Spherocity = pFull;
	}
    }
    
    spherocity=((Spherocity)*TMath::Pi()*TMath::Pi())/4.0;
    
    if(pxA){// clean up array memory used for TMath::Sort
      delete[] pxA;
      pxA=0;
    }
    if(pyA){// clean up array memory used for TMath::Sort
      delete[] pyA;
      pyA=0;
    }    
    
  }//close case when multiplicity is ok!
  
  return spherocity;
  
}
//________________________________________________________
TH1D * AliTransverseEventShape::GetHistData( Int_t bin_histo ){
  
  switch(bin_histo){
  case 0:
    return fheta;
    break;
  case 1:
    return fhphi;
    break;
  case 2:
    return fhpt;
    break;
  default:
    return 0;
  }
  
  
}
//__________________________________________________________
void AliTransverseEventShape::SaveHistos( const char* folder ){
  
  if(!fheta)
    return;
  
  if (folder)
    {
      gDirectory->mkdir(folder);
      gDirectory->cd(folder);
    }
  
  
  fheta->Write();
  fhphi->Write();
  fhpt->Write();
  
  fhetaMC->Write();
  fhphiMC->Write();
  fhptMC->Write();
  
  if (folder)
    gDirectory->cd("..");
  
}


