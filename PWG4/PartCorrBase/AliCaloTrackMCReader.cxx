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
/* $Id:  $ */

//_________________________________________________________________________
// Class for reading data (Kinematics) in order to do prompt gamma 
// or other particle identification and correlations
// Separates generated particles into charged (CTS) 
// and neutral (PHOS or EMCAL acceptance)
//
//*-- Author: Gustavo Conesa (LNF-INFN) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---

#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TRandom.h>
#include <TArrayI.h>
#include "TParticle.h"
//#include "Riostream.h"

//---- ANALYSIS system ----
#include "AliCaloTrackMCReader.h" 
#include "AliGenEventHeader.h"
#include "AliStack.h"
#include "AliAODCaloCluster.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliFiducialCut.h"
#include "AliMCAnalysisUtils.h"

  ClassImp(AliCaloTrackMCReader)

//____________________________________________________________________________
AliCaloTrackMCReader::AliCaloTrackMCReader() : 
  AliCaloTrackReader(), fDecayPi0(0), 
  fNeutralParticlesArray(0x0), fChargedParticlesArray(0x0), 
  fStatusArray(0x0), fKeepAllStatus(0), fCheckOverlap(0),  
  fEMCALOverlapAngle(0),fPHOSOverlapAngle(0), fIndex2ndPhoton(0)
{
  //Ctor
  
  //Initialize parameters
  InitParameters();
  fDataType = kMC;  
  fReadStack          = kTRUE;
  fReadAODMCParticles = kFALSE;
  
}
/*
//____________________________________________________________________________
AliCaloTrackMCReader::AliCaloTrackMCReader(const AliCaloTrackMCReader & g) :   
  AliCaloTrackReader(g), fDecayPi0(g.fDecayPi0), 
  fNeutralParticlesArray(g.fNeutralParticlesArray?new TArrayI(*g.fNeutralParticlesArray):0x0),
  fChargedParticlesArray(g.fChargedParticlesArray?new TArrayI(*g.fChargedParticlesArray):0x0),
  fStatusArray(g.fStatusArray?new TArrayI(*g.fStatusArray):0x0),
  fKeepAllStatus(g.fKeepAllStatus), fCheckOverlap(g.fCheckOverlap),
  fEMCALOverlapAngle( g.fEMCALOverlapAngle), fPHOSOverlapAngle(g.fPHOSOverlapAngle),
  fIndex2ndPhoton(g.fIndex2ndPhoton)
{
  // cpy ctor
}
*/
//_________________________________________________________________________
//AliCaloTrackMCReader & AliCaloTrackMCReader::operator = (const AliCaloTrackMCReader & source)
//{
//  // assignment operator
//
//  if(&source == this) return *this;
//
//  fDecayPi0 = source.fDecayPi0; 
//
//  delete fChargedParticlesArray;
//  fChargedParticlesArray = source.fChargedParticlesArray?new TArrayI(*source.fChargedParticlesArray):0x0;
//
//  delete fNeutralParticlesArray;
//  fNeutralParticlesArray = source.fNeutralParticlesArray?new TArrayI(*source.fNeutralParticlesArray):0x0;
//
//  delete fStatusArray;
//  fStatusArray = source.fStatusArray?new TArrayI(*source.fStatusArray):0x0;
// 
//  fKeepAllStatus = source.fKeepAllStatus ;
//
//  return *this;
//
//}
//
//_________________________________
AliCaloTrackMCReader::~AliCaloTrackMCReader() {
  //Dtor

  if(fChargedParticlesArray) delete fChargedParticlesArray ;
  if(fNeutralParticlesArray) delete fNeutralParticlesArray ;
  if(fStatusArray) delete fStatusArray ;

}

//____________________________________________________________________________
void AliCaloTrackMCReader::GetVertex(Double_t  v[3]) const {
  //Return vertex position

  TArrayF pv;
  GetGenEventHeader()->PrimaryVertex(pv);
  v[0]=pv.At(0);
  v[1]=pv.At(1);
  v[2]=pv.At(2);

}


//_______________________________________________________________
void AliCaloTrackMCReader::InitParameters()
{
  
  //Initialize the parameters of the analysis.

  fDecayPi0 = kFALSE;

  fChargedParticlesArray = new TArrayI(1);
  fChargedParticlesArray->SetAt(11,0);  
  //Int_t pdgarray[]={12,14,16};// skip neutrinos
  //fNeutralParticlesArray = new TArrayI(3, pdgarray);
  fNeutralParticlesArray = new TArrayI(3);
  fNeutralParticlesArray->SetAt(12,0); fNeutralParticlesArray->SetAt(14,1); fNeutralParticlesArray->SetAt(16,2); 
  fStatusArray = new TArrayI(1);
  fStatusArray->SetAt(1,0); 
 
  fKeepAllStatus = kTRUE;

  fCheckOverlap = kFALSE;
  fEMCALOverlapAngle = 2.5 * TMath::DegToRad();
  fPHOSOverlapAngle = 0.5 * TMath::DegToRad();
  fIndex2ndPhoton = -1;
}
//____________________________________________________________________________
void  AliCaloTrackMCReader::CheckOverlap(const Float_t anglethres, const Int_t imom, Int_t & iPrimary, Int_t & index, TLorentzVector & mom, Int_t & pdg) {
  //Check overlap of decay photons
  if( fIndex2ndPhoton==iPrimary ){
    fIndex2ndPhoton=-1;
    return;
  }
  else fIndex2ndPhoton=-1;
  

  if(pdg!=22) return;
  
  TLorentzVector ph1, ph2;
  TParticle *meson = GetStack()->Particle(imom);
  Int_t mepdg = meson->GetPdgCode();
  Int_t idaug1 = meson->GetFirstDaughter();
  if((mepdg == 111 || mepdg == 221 ) && meson->GetNDaughters() == 2){ //Check only decay in 2 photons
    TParticle * d1 = GetStack()->Particle(idaug1);
    TParticle  *d2 = GetStack()->Particle(idaug1+1);
    if(d1->GetPdgCode() == 22 && d2->GetPdgCode() == 22 ){
      d1->Momentum(ph1);
      d2->Momentum(ph2);
      //printf("angle %2.2f\n",ph1.Angle(ph2.Vect()));
      
      if(anglethres >  ph1.Angle(ph2.Vect())){ 	  
	//Keep the meson
	pdg=mepdg;
	index=imom;
	meson->Momentum(mom);
	//printf("Overlap:: pt %2.2f, phi %2.2f, eta %2.2f\n",mom.Pt(),mom.Phi(),mom.Eta());
	if(iPrimary == idaug1) iPrimary++; //skip next photon in list
      }
      else{
	//Do not check overlapping for next decay photon from same meson
	if(iPrimary == idaug1) {fIndex2ndPhoton = idaug1+1;
	}

      }
    }
  }//Meson Decay with 2 photon daughters
}

//____________________________________________________________________________
void  AliCaloTrackMCReader::FillCalorimeters(Int_t & iParticle, TParticle* particle, TLorentzVector momentum,
					     Int_t &ncalo) {
  //Fill AODCaloClusters or TParticles lists of PHOS or EMCAL
  //In PHOS
  if(fFillPHOS && momentum.Pt() > fPHOSPtMin){
	  
	if(!fFiducialCut->IsInFiducialCut(momentum,"PHOS")) return;
	  
    Int_t index = iParticle ;
    Int_t pdg = TMath::Abs(particle->GetPdgCode());
    if(fCheckOverlap) 
      CheckOverlap(fPHOSOverlapAngle,particle->GetFirstMother(),index, iParticle, momentum, pdg);
    
    Char_t ttype= AliAODCluster::kPHOSNeutral;	
    Int_t labels[] = {index};
    Float_t x[] = {momentum.X(), momentum.Y(), momentum.Z()};
    //Create object and write it to file
    AliAODCaloCluster *calo = new((*(fOutputEvent->GetCaloClusters()))[ncalo++]) 
      AliAODCaloCluster(index,1,labels,momentum.E(), x, NULL, ttype, 0);
    
    SetCaloClusterPID(pdg,calo) ;
    if(fDebug > 3 && momentum.Pt() > 0.2)
      printf("AliCaloTrackMCReader::FillCalorimeters() - PHOS : Selected cluster %s E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
	     particle->GetName(),momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());			
    fAODPHOS->Add(calo);//reference the selected object to the list
  }
  
  //In EMCAL
  if(fFillEMCAL  && momentum.Pt() > fEMCALPtMin){
	  
	if(!fFiducialCut->IsInFiducialCut(momentum,"EMCAL")) return;
	  
    Int_t index = iParticle ;
    Int_t pdg = TMath::Abs(particle->GetPdgCode());
    //Int_t pdgorg=pdg;
    if(fCheckOverlap) 
      CheckOverlap(fEMCALOverlapAngle,particle->GetFirstMother(),iParticle, index, momentum, pdg);
    
    Char_t ttype= AliAODCluster::kEMCALClusterv1;
    Int_t labels[] = {index};
    Float_t x[] = {momentum.X(), momentum.Y(), momentum.Z()};
    //Create object and write it to file
    AliAODCaloCluster *calo = new((*(fOutputEvent->GetCaloClusters()))[ncalo++]) 
      AliAODCaloCluster(iParticle,1,labels,momentum.E(), x, NULL, ttype, 0);
    
    SetCaloClusterPID(pdg,calo) ;
    if(fDebug > 3 && momentum.Pt() > 0.2)
      printf("AliCaloTrackMCReader::FillCalorimeters() - EMCAL : Selected cluster %s E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
	     particle->GetName(),momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());	
    fAODEMCAL->Add(calo);//reference the selected object to the list
  }
}

//____________________________________________________________________________
Bool_t AliCaloTrackMCReader::FillInputEvent(const Int_t iEntry, const char * currentFileName){
  //Fill the event counter and input lists that are needed, called by the analysis maker.
  
  fEventNumber = iEntry;
  fCurrentFileName = TString(currentFileName);
	
  //In case of analysis of events with jets, skip those with jet pt > 5 pt hard	
  if(fComparePtHardAndJetPt && GetStack()) {
	if(!ComparePtHardAndJetPt()) return kFALSE ;
  }
	
  Int_t iParticle = 0 ;
  Double_t charge = 0.;
  Int_t ncalo  = (fOutputEvent->GetCaloClusters())->GetEntriesFast();
  Int_t ntrack = (fOutputEvent->GetTracks())->GetEntriesFast();
	
  for (iParticle = 0 ; iParticle <  GetStack()->GetNtrack() ; iParticle++) {
    TParticle * particle = GetStack()->Particle(iParticle);
    TLorentzVector momentum;
    Float_t p[3];
    Float_t x[3];
    Int_t pdg = particle->GetPdgCode();						
    
    //Keep particles with a given status 
    if(KeepParticleWithStatus(particle->GetStatusCode()) && (particle->Pt() > 0) ){
      
      //Skip bizarre particles, they crash when charge is calculated
      //	if(TMath::Abs(pdg) == 3124 || TMath::Abs(pdg) > 10000000) continue ;
      
      charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
      particle->Momentum(momentum);
      //---------- Charged particles ----------------------
      if(charge != 0){
		 if(fFillCTS && (momentum.Pt() > fCTSPtMin)){
	  //Particles in CTS acceptance
		
	  if(!fFiducialCut->IsInFiducialCut(momentum,"CTS")) continue;
		
	  if(fDebug > 3 && momentum.Pt() > 0.2)
	    printf("AliCaloTrackMCReader::FillInputEvent() - CTS : Selected tracks E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f\n",
		   momentum.E(),momentum.Pt(),momentum.Phi()*TMath::RadToDeg(),momentum.Eta());
	  
	  x[0] = particle->Vx(); x[1] = particle->Vy(); x[2] = particle->Vz();
	  p[0] = particle->Px(); p[1] = particle->Py(); p[2] = particle->Pz();
	  //Create object and write it to file
	  AliAODTrack *aodTrack = new((*(fOutputEvent->GetTracks()))[ntrack++]) 
	    AliAODTrack(0, iParticle, p, kTRUE, x, kFALSE,NULL, 0, 0, 
			NULL,
			0x0,//primary,
			kFALSE, // No fit performed
			kFALSE, // No fit performed
			AliAODTrack::kPrimary, 
			0);
	  SetTrackChargeAndPID(pdg, aodTrack);
	  fAODCTS->Add(aodTrack);//reference the selected object to the list
	}
	//Keep some charged particles in calorimeters lists
	if((fFillPHOS || fFillEMCAL) && KeepChargedParticles(pdg)) FillCalorimeters(iParticle, particle, momentum, ncalo);
	
      }//Charged
      
      //-------------Neutral particles ----------------------
      else if(charge == 0 && (fFillPHOS || fFillEMCAL)){
	//Skip neutrinos or other neutral particles
	//if(SkipNeutralParticles(pdg) || particle->IsPrimary()) continue ; // protection added (MG)
	if(SkipNeutralParticles(pdg)) continue ;
	//Fill particle/calocluster arrays
	if(!fDecayPi0) {
	  FillCalorimeters(iParticle, particle, momentum, ncalo);
	}
	else {
	  //Sometimes pi0 are stable for the generator, if needed decay it by hand
	  if(pdg == 111 ){
	    if(momentum.Pt() >  fPHOSPtMin || momentum.Pt() >  fEMCALPtMin){
	      TLorentzVector lvGamma1, lvGamma2 ;
	      //Double_t angle = 0;
	      
	      //Decay
	      MakePi0Decay(momentum,lvGamma1,lvGamma2);//,angle);
	      
	      //Check if Pi0 is in the acceptance of the calorimeters, if aperture angle is small, keep it
	      TParticle * pPhoton1 = new TParticle(22,1,iParticle,0,0,0,lvGamma1.Px(),lvGamma1.Py(),
						   lvGamma1.Pz(),lvGamma1.E(),0,0,0,0);   
	      TParticle * pPhoton2 = new TParticle(22,1,iParticle,0,0,0,lvGamma2.Px(),lvGamma2.Py(),
						   lvGamma2.Pz(),lvGamma2.E(),0,0,0,0);
	      //Fill particle/calocluster arrays
	      FillCalorimeters(iParticle,pPhoton1,lvGamma1, ncalo);
	      FillCalorimeters(iParticle,pPhoton2,lvGamma2, ncalo);
	    }//pt cut
	  }//pi0
	  else FillCalorimeters(iParticle,particle, momentum, ncalo); //Add the rest
	}
      }//neutral particles
    } //particle with correct status
  }//particle loop
  
  fIndex2ndPhoton = -1; //In case of overlapping studies, reset for each event	
 	
  return kTRUE;	

}

//________________________________________________________________
Bool_t AliCaloTrackMCReader::KeepParticleWithStatus(Int_t status) const {
  //Check if status is equal to one of the  list
  //These particles will be used in analysis.
  if(!fKeepAllStatus){
    for(Int_t i= 0; i < fStatusArray->GetSize(); i++)
      if(status ==  fStatusArray->At(i)) return kTRUE ;

    return kFALSE; 
    
  }
  else
    return kTRUE ;  
}

//________________________________________________________________
Bool_t AliCaloTrackMCReader::KeepChargedParticles(Int_t pdg) const {
  //Check if pdg is equal to one of the charged particles list
  //These particles will be added to the calorimeters lists.

  for(Int_t i= 0; i < fChargedParticlesArray->GetSize(); i++)
    if(TMath::Abs(pdg) ==  fChargedParticlesArray->At(i)) return kTRUE ;
  
  return kFALSE; 
  
}

//________________________________________________________________
void AliCaloTrackMCReader::Print(const Option_t * opt) const
{
  //Print some relevant parameters set for the analysis
  if(! opt)
    return;
  
  printf("**** Print **** %s %s ****\n", GetName(), GetTitle() ) ;
  
  printf("Decay Pi0?          : %d\n", fDecayPi0) ;
  printf("Check Overlap in Calo?    : %d\n", fCheckOverlap) ;
  printf("Keep all status?    : %d\n", fKeepAllStatus) ;
  
  if(!fKeepAllStatus) printf("Keep particles with status : ");
  for(Int_t i= 0; i < fStatusArray->GetSize(); i++)
    printf(" %d ; ", fStatusArray->At(i));
  printf("\n");
  
  printf("Skip neutral particles in calo : ");
  for(Int_t i= 0; i < fNeutralParticlesArray->GetSize(); i++)
    printf(" %d ; ", fNeutralParticlesArray->At(i));
  printf("\n");
  
  printf("Keep charged particles in calo : ");
  for(Int_t i= 0; i < fChargedParticlesArray->GetSize(); i++)
    printf(" %d ; ", fChargedParticlesArray->At(i));
  printf("\n");

}

//____________________________________________________________________________
void AliCaloTrackMCReader::MakePi0Decay(TLorentzVector &p0, TLorentzVector &p1, 
				TLorentzVector &p2) const {//, Double_t &angle) {
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

//____________________________________________________________________________
void AliCaloTrackMCReader::SetInputOutputMCEvent(AliVEvent* /*esd*/, AliAODEvent* aod, AliMCEvent* mc) {
  // Connect the data pointer
  SetMC(mc);
  SetOutputEvent(aod);
}


//________________________________________________________________
Bool_t AliCaloTrackMCReader::SkipNeutralParticles(Int_t pdg) const {
  //Check if pdg is equal to one of the neutral particles list
  //These particles will be skipped from analysis.

  for(Int_t i= 0; i < fNeutralParticlesArray->GetSize(); i++)
    if(TMath::Abs(pdg) ==  fNeutralParticlesArray->At(i)) return kTRUE ;
  
  return kFALSE; 
  
}


//____________________________________________________________________
void AliCaloTrackMCReader::SetTrackChargeAndPID(const Int_t pdgCode, AliAODTrack *track) const {
//Give a PID weight for tracks equal to 1 depending on the particle type

  Float_t pid[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

  switch (pdgCode) {

  case 22: // gamma
    track->SetCharge(0);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;

  case 11: // e- 
    track->SetCharge(-1);
    pid[AliAODTrack::kElectron] = 1.;
    track->SetPID(pid);
    break;
    
  case -11: // e+
    track->SetCharge(+1);
    pid[AliAODTrack::kElectron] = 1.;
    track->SetPID(pid);
    break;
    
  case 13: // mu- 
    track->SetCharge(-1);
    pid[AliAODTrack::kMuon] = 1.;
    track->SetPID(pid);
    break;
    
  case -13: // mu+
    track->SetCharge(+1);
    pid[AliAODTrack::kMuon] = 1.;
    track->SetPID(pid);
    break;
    
  case 111: // pi0
    track->SetCharge(0);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;
    
  case 211: // pi+
    track->SetCharge(+1);
    pid[AliAODTrack::kPion] = 1.;
    track->SetPID(pid);
    break;
    
  case -211: // pi-
    track->SetCharge(-1);
    pid[AliAODTrack::kPion] = 1.;
    track->SetPID(pid);
    break;
    
  case 130: // K0L
    track->SetCharge(0);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;
    
  case 321: // K+
    track->SetCharge(+1);
    pid[AliAODTrack::kKaon] = 1.;
    track->SetPID(pid);
    break;
    
  case -321: // K- 
    track->SetCharge(-1);
    pid[AliAODTrack::kKaon] = 1.;
    track->SetPID(pid);
    break;
    
  case 2112: // n
    track->SetCharge(0);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;
    
  case 2212: // p
    track->SetCharge(+1);
    pid[AliAODTrack::kProton] = 1.;
    track->SetPID(pid);
    break;
    
  case -2212: // anti-p
    track->SetCharge(-1);
    pid[AliAODTrack::kProton] = 1.;
    track->SetPID(pid);
    break;

  case 310: // K0S
    track->SetCharge(0);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;
    
  case 311: // K0
    track->SetCharge(0);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;
    
  case -311: // anti-K0
    track->SetCharge(0);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;
    
  case 221: // eta
    track->SetCharge(0);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;

  case 3122: // lambda
    track->SetCharge(0);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;

  case 3222: // Sigma+
    track->SetCharge(+1);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;

  case 3212: // Sigma0
    track->SetCharge(-1);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;

  case 3112: // Sigma-
    track->SetCharge(-1);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;

  case 3322: // Xi0
    track->SetCharge(0);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;

  case 3312: // Xi-
    track->SetCharge(-1);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;

  case 3334: // Omega-
    track->SetCharge(-1);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;

  case -2112: // n-bar
    track->SetCharge(0);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;

  case -3122: // anti-Lambda
    track->SetCharge(0);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;

  case -3222: // anti-Sigma-
    track->SetCharge(-1);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;

  case -3212: // anti-Sigma0
    track->SetCharge(0);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;

  case -3112: // anti-Sigma+
    track->SetCharge(+1);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;

  case -3322: // anti-Xi0
    track->SetCharge(0);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;

  case -3312: // anti-Xi+
    track->SetCharge(+1);
    break;

  case -3334: // anti-Omega+
    track->SetCharge(+1);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;

  case 411: // D+
    track->SetCharge(+1);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;

  case -411: // D- 
    track->SetCharge(-1);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;

  case 421: // D0
    track->SetCharge(0);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;

  case -421: // anti-D0
    track->SetCharge(0);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
    break;

  default : // unknown
    track->SetCharge(-99);
    pid[AliAODTrack::kUnknown] = 1.;
    track->SetPID(pid);
 }

  track->SetPID(pid);

  return;
}

//____________________________________________________________________
void AliCaloTrackMCReader::SetCaloClusterPID(const Int_t pdgCode, AliAODCaloCluster *calo) const {
//Give a PID weight for CaloClusters equal to 1 depending on the particle type

  Float_t pid[13] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

  switch (pdgCode) {

  case 22: // gamma
    pid[AliAODCaloCluster::kPhoton] = 1.;
    calo->SetPID(pid);
    break;

  case 11: // e- 
    pid[AliAODCaloCluster::kElectron] = 1.;
    calo->SetPID(pid);
    break;
    
  case -11: // e+
    pid[AliAODCaloCluster::kElectron] = 1.;
    calo->SetPID(pid);
    break;
    
  case 13: // mu- 
    pid[AliAODCaloCluster::kCharged] = 1.;
    calo->SetPID(pid);
    break;
    
  case -13: // mu+
    pid[AliAODCaloCluster::kCharged] = 1.;
    calo->SetPID(pid);
    break;
    
  case 111: // pi0
    pid[AliAODCaloCluster::kPi0] = 1.;
    calo->SetPID(pid);
    break;
    
  case 211: // pi+
    pid[AliAODCaloCluster::kCharged] = 1.;
    calo->SetPID(pid);
    break;
    
  case -211: // pi-
    pid[AliAODCaloCluster::kCharged] = 1.;
    calo->SetPID(pid);
    break;
    
  case 130: // K0L
    pid[AliAODCaloCluster::kKaon0] = 1.;
    pid[AliAODCaloCluster::kNeutral] = 1;
    calo->SetPID(pid);
    break;
    
  case 321: // K+
    pid[AliAODCaloCluster::kCharged] = 1.;
    calo->SetPID(pid);
    break;
    
  case -321: // K- 
    pid[AliAODCaloCluster::kCharged] = 1.;
    calo->SetPID(pid);
    break;
    
  case 2112: // n
    pid[AliAODCaloCluster::kNeutron] = 1.;
    pid[AliAODCaloCluster::kNeutral] = 1.;
    calo->SetPID(pid);
    break;
    
  case 2212: // p
    pid[AliAODCaloCluster::kCharged] = 1.;
    calo->SetPID(pid);
    break;
    
  case -2212: // anti-p
    pid[AliAODCaloCluster::kCharged] = 1.;
    calo->SetPID(pid);
    break;

  case 310: // K0S
    pid[AliAODCaloCluster::kKaon0] = 1.;
    pid[AliAODCaloCluster::kNeutral] = 1.;
    calo->SetPID(pid);
    break;
    
  case 311: // K0
    pid[AliAODCaloCluster::kKaon0] = 1.;
    pid[AliAODCaloCluster::kNeutral] = 1.;
    calo->SetPID(pid);
    break;
    
  case -311: // anti-K0
    pid[AliAODCaloCluster::kKaon0] = 1.;
    pid[AliAODCaloCluster::kNeutral] = 1.;
    calo->SetPID(pid);
    break;
    
  case 221: // eta
    pid[AliAODCaloCluster::kNeutral] = 1.;
    calo->SetPID(pid);
    break;

  case 3122: // lambda
    pid[AliAODCaloCluster::kUnknown] = 1.;
    calo->SetPID(pid);
    break;

  case 3222: // Sigma+
    pid[AliAODCaloCluster::kUnknown] = 1.;
    calo->SetPID(pid);
    break;

  case 3212: // Sigma0
    pid[AliAODCaloCluster::kUnknown] = 1.;
    calo->SetPID(pid);
    break;

  case 3112: // Sigma-
    pid[AliAODCaloCluster::kUnknown] = 1.;
    calo->SetPID(pid);
    break;

  case 3322: // Xi0
    pid[AliAODCaloCluster::kUnknown] = 1.;
    calo->SetPID(pid);
    break;

  case 3312: // Xi-
    pid[AliAODCaloCluster::kUnknown] = 1.;
    calo->SetPID(pid);
    break;

  case 3334: // Omega-
    pid[AliAODCaloCluster::kUnknown] = 1.;
    calo->SetPID(pid);
    break;

  case -2112: // n-bar
    pid[AliAODCaloCluster::kNeutron] = 1.;
    pid[AliAODCaloCluster::kNeutral] = 1.;
    calo->SetPID(pid);
    break;

  case -3122: // anti-Lambda
    pid[AliAODCaloCluster::kUnknown] = 1.;
    calo->SetPID(pid);
    break;

  case -3222: // anti-Sigma-
    pid[AliAODCaloCluster::kUnknown] = 1.;
    calo->SetPID(pid);
    break;

  case -3212: // anti-Sigma0
    pid[AliAODCaloCluster::kUnknown] = 1.;
    calo->SetPID(pid);
    break;

  case -3112: // anti-Sigma+
    pid[AliAODCaloCluster::kUnknown] = 1.;
    calo->SetPID(pid);
    break;

  case -3322: // anti-Xi0
    pid[AliAODCaloCluster::kUnknown] = 1.;
    calo->SetPID(pid);
    break;

  case -3312: // anti-Xi+
    pid[AliAODCaloCluster::kUnknown] = 1.;
    calo->SetPID(pid);
    break;

  case -3334: // anti-Omega+
    pid[AliAODCaloCluster::kUnknown] = 1.;
    calo->SetPID(pid);
    break;

  case 411: // D+
    pid[AliAODCaloCluster::kUnknown] = 1.;
    calo->SetPID(pid);
    break;

  case -411: // D- 
    pid[AliAODCaloCluster::kUnknown] = 1.;
    calo->SetPID(pid);
    break;

  case 421: // D0
    pid[AliAODCaloCluster::kUnknown] = 1.;
    calo->SetPID(pid);
    break;

  case -421: // anti-D0
    pid[AliAODCaloCluster::kUnknown] = 1.;
    calo->SetPID(pid);
    break;

  default : // unknown
    pid[AliAODCaloCluster::kUnknown] = 1.;
    calo->SetPID(pid);
 }

 
  return;
}
