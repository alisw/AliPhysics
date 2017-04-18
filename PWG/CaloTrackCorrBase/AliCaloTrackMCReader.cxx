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

// --- ROOT system ---
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TRandom.h>
#include <TArrayI.h>
#include "TParticle.h"

//---- ANALYSIS system ----
#include "AliCaloTrackMCReader.h" 
#include "AliGenEventHeader.h"
#include "AliStack.h"
#include "AliVCluster.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliFiducialCut.h"
#include "AliMCAnalysisUtils.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliCaloTrackMCReader) ;
/// \endcond 

//____________________________________________
/// Default constructor. Initialize parameters
//____________________________________________
AliCaloTrackMCReader::AliCaloTrackMCReader() : 
AliCaloTrackReader(),        fDecayPi0(0), 
fNeutralParticlesArray(0x0), fChargedParticlesArray(0x0), 
fStatusArray(0x0),           fKeepAllStatus(0), 
fCheckOverlap(0),            fEMCALOverlapAngle(0),fPHOSOverlapAngle(0), 
fIndex2ndPhoton(0),          fOnlyGeneratorParticles(kTRUE),
fMomentum(),                 fPi0Momentum(),
fGamDecayMom1(),             fGamDecayMom2()
{
  InitParameters();
}

//___________________________________________
/// Destructor.
//___________________________________________
AliCaloTrackMCReader::~AliCaloTrackMCReader()
{
  AliCaloTrackReader::DeletePointers();

  if(fChargedParticlesArray) delete fChargedParticlesArray ;
  if(fNeutralParticlesArray) delete fNeutralParticlesArray ;
  if(fStatusArray)           delete fStatusArray ;
}

//________________________________________________________
/// \return Vertex position from header.
//________________________________________________________
void AliCaloTrackMCReader::GetVertex(Double_t  v[3]) const 
{  
  TArrayF pv;
  GetGenEventHeader()->PrimaryVertex(pv);
  v[0]=pv.At(0);
  v[1]=pv.At(1);
  v[2]=pv.At(2);
}

//_________________________________________
/// Initialize the parameters of the analysis.
//_________________________________________
void AliCaloTrackMCReader::InitParameters()
{
  fDecayPi0 = kFALSE;
  
  fChargedParticlesArray = new TArrayI(1);
  fChargedParticlesArray->SetAt(11,0);  
  //Int_t pdgarray[]={12,14,16};// skip neutrinos
  //fNeutralParticlesArray = new TArrayI(3, pdgarray);
  fNeutralParticlesArray = new TArrayI(3);
  fNeutralParticlesArray->SetAt(12,0); fNeutralParticlesArray->SetAt(14,1); fNeutralParticlesArray->SetAt(16,2); 
  fStatusArray = new TArrayI(1);
  fStatusArray->SetAt(1,0); 
  
  fOnlyGeneratorParticles = kTRUE;
  fKeepAllStatus          = kTRUE;
  
  fCheckOverlap       = kFALSE;
  fEMCALOverlapAngle  = 2.5 * TMath::DegToRad();
  fPHOSOverlapAngle   = 0.5 * TMath::DegToRad();
  fIndex2ndPhoton     = -1;
  
  fDataType           = kMC;  
  fReadStack          = kTRUE;
  fReadAODMCParticles = kFALSE; //This class only works with Kinematics.root input.
  
  //For this reader we own the objects of the arrays
  fCTSTracks    ->SetOwner(kTRUE); 
  fEMCALClusters->SetOwner(kTRUE); 
  fPHOSClusters ->SetOwner(kTRUE); 
}

//____________________________________________________________________________________
/// Check overlap of decay photons
//____________________________________________________________________________________
void  AliCaloTrackMCReader::CheckOverlap(Float_t anglethres, Int_t imom,
                                         Int_t & iPrimary, Int_t & index, Int_t & pdg)
{
  if( fIndex2ndPhoton==iPrimary )
  {
    fIndex2ndPhoton=-1;
    return;
  }
  else 
    fIndex2ndPhoton=-1;
  
  if(pdg!=22) return;
  
  TParticle *meson = GetStack()->Particle(imom);
  Int_t mepdg  = meson->GetPdgCode();
  Int_t idaug1 = meson->GetFirstDaughter();
  if((mepdg == 111 || mepdg == 221 ) && meson->GetNDaughters() == 2)
  { //Check only decay in 2 photons
    TParticle * d1 = GetStack()->Particle(idaug1);
    TParticle * d2 = GetStack()->Particle(idaug1+1);
    if(d1->GetPdgCode() == 22 && d2->GetPdgCode() == 22 )
    {
      d1->Momentum(fGamDecayMom1);
      d2->Momentum(fGamDecayMom2);
      //printf("angle %2.2f\n",ph1.Angle(ph2.Vect()));
      
      if(anglethres >  fGamDecayMom1.Angle(fGamDecayMom2.Vect()))
      { 	  
        //Keep the meson
        pdg=mepdg;
        index=imom;
        meson->Momentum(fMomentum);
        //printf("Overlap:: pt %2.2f, phi %2.2f, eta %2.2f\n",mom.Pt(),mom.Phi(),mom.Eta());
        if(iPrimary == idaug1) iPrimary++; //skip next photon in list
      }
      else
      {
        //Do not check overlapping for next decay photon from same meson
        if(iPrimary == idaug1)
        {
          fIndex2ndPhoton = idaug1+1;
        }
        
      }
    }
  } // Meson Decay with 2 photon daughters
}

//__________________________________________________________________________________
/// Fill CaloClusters or TParticles lists of PHOS or EMCAL.
//__________________________________________________________________________________
void  AliCaloTrackMCReader::FillCalorimeters(Int_t & iParticle, TParticle* particle)
{
  Char_t  ttype          = 0;
  Float_t overAngleLimit = 100;
  
  if     (fFillPHOS)
  {
    if( particle->Pt() < fPHOSPtMin ) return;
    if( fCheckFidCut && !fFiducialCut->IsInFiducialCut(particle->Eta(),particle->Phi(),kPHOS )) return;
    ttype = AliVCluster::kPHOSNeutral;
    overAngleLimit = fPHOSOverlapAngle;
  }
  else if(fFillEMCAL)
  {
    if( particle->Pt() < fEMCALPtMin ) return;
    if( fCheckFidCut && !fFiducialCut->IsInFiducialCut(particle->Eta(),particle->Phi(),kEMCAL)) return;
    ttype= AliVCluster::kEMCALClusterv1;
    overAngleLimit = fEMCALOverlapAngle;
  }
  
  particle->Momentum(fMomentum);
  
  Int_t index = iParticle ;
  Int_t pdg = TMath::Abs(particle->GetPdgCode());
  if(fCheckOverlap)
    CheckOverlap(overAngleLimit,particle->GetFirstMother(),index, iParticle, pdg);
  
  Int_t labels[] = {index};
  Double_t x[] = {fMomentum.X(), fMomentum.Y(), fMomentum.Z()};
  
  // Create object and write it to file
  AliAODCaloCluster *calo = new AliAODCaloCluster(index,1,labels,fMomentum.E(), x, NULL, ttype, 0);
  
  SetCaloClusterPID(pdg,calo) ;
 
  AliDebug(3,Form("PHOS %d?, EMCAL %d? : Selected cluster pdg %d, E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f",
                  fFillPHOS, fFillEMCAL, pdg,fMomentum.E(), fMomentum.Pt(), fMomentum.Phi()*TMath::RadToDeg(),fMomentum.Eta()));
  
  // Reference the selected object to the list
  if(fFillPHOS)fPHOSClusters ->Add(calo);
  else         fEMCALClusters->Add(calo);
}

//___________________________________________________________________________
/// Fill the event counter and input lists that are needed, called by AliAnaCaloTrackCorrMaker.
//___________________________________________________________________________
Bool_t AliCaloTrackMCReader::FillInputEvent(Int_t iEntry,
                                            const char * /*currentFileName*/)
{  
  fEventNumber     = iEntry;
  //fCurrentFileName = TString(currentFileName);
  
  for(Int_t iptCut = 0; iptCut < fTrackMultNPtCut; iptCut++ )
  {
    fTrackMult [iptCut] = 0;
    fTrackSumPt[iptCut] = 0;
  }  
  
  // In case of analysis of events with jets, skip those with jet pt > 5 pt hard	
  if(fComparePtHardAndJetPt && GetStack()) 
  {
    if(!ComparePtHardAndJetPt()) return kFALSE ;
  }
	
  // Fill Vertex array
  FillVertexArray();
  
  Int_t iParticle  = 0 ;
  Double_t charge  = 0.;
  Int_t nparticles = GetStack()->GetNtrack() ;
  
  if(fOnlyGeneratorParticles) nparticles=GetStack()->GetNprimary();
  
  for (iParticle = 0 ; iParticle <  nparticles ; iParticle++) 
  {
    TParticle * particle = GetStack()->Particle(iParticle);
    Float_t p[3];
    Float_t x[3];
    Int_t pdg = particle->GetPdgCode();						
    
    //Keep particles with a given status 
    if(KeepParticleWithStatus(particle->GetStatusCode()) && (particle->Pt() > 0) )
    {
      //Skip bizarre particles, they crash when charge is calculated
      //	if(TMath::Abs(pdg) == 3124 || TMath::Abs(pdg) > 10000000) continue ;
      
      charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
      
      Float_t en  = particle->Energy();
      Float_t pt  = particle->Pt();
      Float_t eta = particle->Eta();
      Float_t phi = particle->Phi();
      //---------- Charged particles ----------------------
      if(charge != 0)
      {
//      if(TMath::Abs(eta)< fTrackMultEtaCut) fTrackMult++;
        if(TMath::Abs(eta)< fTrackMultEtaCut) 
        {
          for(Int_t iptCut = 0; iptCut < fTrackMultNPtCut; iptCut++ )
          {
            if(pt > fTrackMultPtCut[iptCut]) 
            {
              fTrackMult [iptCut]++;
              fTrackSumPt[iptCut]+=pt;
            }
          }
        }

//      if(fFillCTS && (pt > fCTSPtMin))
        if(fFillCTS && (fCTSPtMin > pt || fCTSPtMax < pt)) continue ;
        {
          // Particles in CTS acceptance
          
          if(fCheckFidCut && !fFiducialCut->IsInFiducialCut(eta,phi,kCTS)) continue;
          
          if(TMath::Abs(pdg) == 11 && GetStack()->Particle(particle->GetFirstMother())->GetPdgCode()==22) continue ;
          
          AliDebug(2,Form("CTS : Selected tracks E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f",
                          en,pt,phi*TMath::RadToDeg(),eta));
          
          x[0] = particle->Vx(); x[1] = particle->Vy(); x[2] = particle->Vz();
          p[0] = particle->Px(); p[1] = particle->Py(); p[2] = particle->Pz();
          //Create object and write it to file
          AliAODTrack *aodTrack = new AliAODTrack(0, iParticle, p, kTRUE, x, kFALSE,NULL, 0, 0, 
                                                  // NULL,
                                                  0x0,//primary,
                                                  kFALSE, // No fit performed
                                                  kFALSE, // No fit performed
                                                  AliAODTrack::kPrimary, 
                                                  0);
          SetTrackChargeAndPID(pdg, aodTrack);
          fCTSTracks->Add(aodTrack);//reference the selected object to the list
        }
        //Keep some charged particles in calorimeters lists
        if((fFillPHOS || fFillEMCAL) && KeepChargedParticles(pdg)) FillCalorimeters(iParticle, particle);
        
      }//Charged
      
      //-------------Neutral particles ----------------------
      else if(charge == 0 && (fFillPHOS || fFillEMCAL))
      {
        //Skip neutrinos or other neutral particles
        //if(SkipNeutralParticles(pdg) || particle->IsPrimary()) continue ; // protection added (MG)
        if(SkipNeutralParticles(pdg)) continue ;
        //Fill particle/calocluster arrays
        if(!fDecayPi0) 
        {
          AliDebug(2,Form("Calo : Selected tracks E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f",
                          en,pt,phi*TMath::RadToDeg(),eta));
          FillCalorimeters(iParticle, particle);
        }
        else 
        {
          //Sometimes pi0 are stable for the generator, if needed decay it by hand
          if(pdg == 111 )
          {
            if(pt >  fPHOSPtMin || pt >  fEMCALPtMin)
            {
              //Decay
              //Double_t angle = 0;
              particle->Momentum(fPi0Momentum);

              MakePi0Decay();//,angle);
              
              //Check if Pi0 is in the acceptance of the calorimeters, if aperture angle is small, keep it
              TParticle * pPhoton1 = new TParticle(22,1,iParticle,0,0,0,fGamDecayMom1.Px(),fGamDecayMom1.Py(),
                                                   fGamDecayMom1.Pz(),fGamDecayMom1.E(),0,0,0,0);
              TParticle * pPhoton2 = new TParticle(22,1,iParticle,0,0,0,fGamDecayMom2.Px(),fGamDecayMom2.Py(),
                                                   fGamDecayMom2.Pz(),fGamDecayMom2.E(),0,0,0,0);
              //Fill particle/calocluster arrays
              FillCalorimeters(iParticle,pPhoton1);
              FillCalorimeters(iParticle,pPhoton2);
            }//pt cut
          }//pi0
          else
          {
            AliDebug(2,Form("Calo : Selected tracks E %3.2f, pt %3.2f, phi %3.2f, eta %3.2f",
                            en,pt,phi*TMath::RadToDeg(),eta));
            FillCalorimeters(iParticle,particle); //Add the rest
          }
        }
      }//neutral particles
    } //particle with correct status
  }//particle loop
  
  fIndex2ndPhoton = -1; //In case of overlapping studies, reset for each event	
 	
  return kTRUE;	
}

//_____________________________________________________________________
/// Check if status requiered for the particles is equal to one of the list.
/// These particles will be used in analysis.
//_____________________________________________________________________
Bool_t AliCaloTrackMCReader::KeepParticleWithStatus(Int_t status) const 
{
  if(!fKeepAllStatus)
  {
    for(Int_t i= 0; i < fStatusArray->GetSize(); i++)
      if(status ==  fStatusArray->At(i)) return kTRUE ;
    
    return kFALSE; 
  }
  else
    return kTRUE ;  
}

//________________________________________________________________
/// Check if pdg is equal to one of the charged particles list
/// These particles will be added to the calorimeters lists.
//________________________________________________________________
Bool_t AliCaloTrackMCReader::KeepChargedParticles(Int_t pdg) const 
{
  for(Int_t i= 0; i < fChargedParticlesArray->GetSize(); i++)
  {
    if(TMath::Abs(pdg) ==  fChargedParticlesArray->At(i)) return kTRUE ;
  }
  
  return kFALSE; 
}

//__________________________________________________________
/// Print some relevant parameters set for the analysis.
//__________________________________________________________
void AliCaloTrackMCReader::Print(const Option_t * opt) const
{
  if(! opt)
    return;
  
  AliCaloTrackReader::Print(opt);
  
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

//_______________________________________
/// Perform isotropic decay pi0 -> 2 photons.
/// fPi0Momentum is pi0 4-momentum (ipnut).
/// fGamDecayMom1 and fGamDecayMom2 are photon 4-momenta (output).
//_______________________________________
void AliCaloTrackMCReader::MakePi0Decay() //, Double_t &angle)
{   
  //  cout<<"Boost vector"<<endl;
  Double_t mPi0 = TDatabasePDG::Instance()->GetParticle(111)->Mass();
  TVector3 b = fPi0Momentum.BoostVector();
  //cout<<"Parameters"<<endl;
  //Double_t mPi0   = fPi0Momentum.M();
  Double_t phi    = TMath::TwoPi() * gRandom->Rndm();
  Double_t cosThe = 2 * gRandom->Rndm() - 1;
  Double_t cosPhi = TMath::Cos(phi);
  Double_t sinPhi = TMath::Sin(phi);
  Double_t sinThe = TMath::Sqrt(1-cosThe*cosThe);
  Double_t ePi0   = mPi0/2.;
  //cout<<"ePi0 "<<ePi0<<endl;
  //cout<<"Components"<<endl;
  fGamDecayMom1.SetPx(+ePi0*cosPhi*sinThe);
  fGamDecayMom1.SetPy(+ePi0*sinPhi*sinThe);
  fGamDecayMom1.SetPz(+ePi0*cosThe);
  fGamDecayMom1.SetE(ePi0);
  //cout<<"fGamDecayMom1: "<<fGamDecayMom1.Px()<<" "<<fGamDecayMom1.Py()<<" "<<fGamDecayMom1.Pz()<<" "<<fGamDecayMom1.E()<<endl;
  //cout<<"fGamDecayMom1 Mass: "<<fGamDecayMom1.Px()*fGamDecayMom1.Px()+fGamDecayMom1.Py()*fGamDecayMom1.Py()+fGamDecayMom1.Pz()*fGamDecayMom1.Pz()-fGamDecayMom1.E()*fGamDecayMom1.E()<<endl;
  fGamDecayMom2.SetPx(-ePi0*cosPhi*sinThe);
  fGamDecayMom2.SetPy(-ePi0*sinPhi*sinThe);
  fGamDecayMom2.SetPz(-ePi0*cosThe);
  fGamDecayMom2.SetE(ePi0);
  //cout<<"fGamDecayMom2: "<<fGamDecayMom2.Px()<<" "<<fGamDecayMom2.Py()<<" "<<fGamDecayMom2.Pz()<<" "<<fGamDecayMom2.E()<<endl;
  //cout<<"fGamDecayMom2 Mass: "<<fGamDecayMom2.Px()*fGamDecayMom2.Px()+fGamDecayMom2.Py()*fGamDecayMom2.Py()+fGamDecayMom2.Pz()*fGamDecayMom2.Pz()-fGamDecayMom2.E()*fGamDecayMom2.E()<<endl;
  //cout<<"Boost "<<b.X()<<" "<<b.Y()<<" "<<b.Z()<<endl;
  fGamDecayMom1.Boost(b);
  //cout<<"fGamDecayMom1: "<<fGamDecayMom1.Px()<<" "<<fGamDecayMom1.Py()<<" "<<fGamDecayMom1.Pz()<<" "<<fGamDecayMom1.E()<<endl;
  fGamDecayMom2.Boost(b);
  //cout<<"fGamDecayMom2: "<<fGamDecayMom2.Px()<<" "<<fGamDecayMom2.Py()<<" "<<fGamDecayMom2.Pz()<<" "<<fGamDecayMom2.E()<<endl;
  //cout<<"angle"<<endl;
  //angle = fGamDecayMom1.Angle(fGamDecayMom2.Vect());
  //cout<<angle<<endl;
}

//__________________________________________________________________
/// Connect the input data pointer.
//__________________________________________________________________
void AliCaloTrackMCReader::SetInputOutputMCEvent(AliVEvent* /*esd*/, 
                                                 AliAODEvent* aod, 
                                                 AliMCEvent* mc) 
{
  SetMC(mc);
  SetOutputEvent(aod);
}

//________________________________________________________________
/// Check if pdg is equal to one of the neutral particles list
/// These particles will be skipped from analysis.
//________________________________________________________________
Bool_t AliCaloTrackMCReader::SkipNeutralParticles(Int_t pdg) const 
{
  for(Int_t i= 0; i < fNeutralParticlesArray->GetSize(); i++)
  {
    if(TMath::Abs(pdg) ==  fNeutralParticlesArray->At(i)) return kTRUE ;
  }
  
  return kFALSE; 
}

//_______________________________________________________________________
/// Give a PID weight for tracks equal to 1 depending on the particle type.
//_______________________________________________________________________
void AliCaloTrackMCReader::SetTrackChargeAndPID(Int_t pdgCode,
                                                AliAODTrack *track) const 
{  
  Float_t pid[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  
  switch (pdgCode) 
  {
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
/// Give a PID weight for CaloClusters equal to 1 depending on the particle type.
//____________________________________________________________________
void AliCaloTrackMCReader::SetCaloClusterPID(const Int_t pdgCode, 
                                             AliVCluster *calo) const 
{  
  Float_t pid[13] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  
  switch (pdgCode) 
  {
    case 22: // gamma
      pid[AliVCluster::kPhoton] = 1.;
      calo->SetPID(pid);
      break;
      
    case 11: // e- 
      pid[AliVCluster::kElectron] = 1.;
      calo->SetPID(pid);
      break;
      
    case -11: // e+
      pid[AliVCluster::kElectron] = 1.;
      calo->SetPID(pid);
      break;
      
    case 13: // mu- 
      pid[AliVCluster::kCharged] = 1.;
      calo->SetPID(pid);
      break;
      
    case -13: // mu+
      pid[AliVCluster::kCharged] = 1.;
      calo->SetPID(pid);
      break;
      
    case 111: // pi0
      pid[AliVCluster::kPi0] = 1.;
      calo->SetPID(pid);
      break;
      
    case 211: // pi+
      pid[AliVCluster::kCharged] = 1.;
      calo->SetPID(pid);
      break;
      
    case -211: // pi-
      pid[AliVCluster::kCharged] = 1.;
      calo->SetPID(pid);
      break;
      
    case 130: // K0L
      pid[AliVCluster::kKaon0] = 1.;
      pid[AliVCluster::kNeutral] = 1;
      calo->SetPID(pid);
      break;
      
    case 321: // K+
      pid[AliVCluster::kCharged] = 1.;
      calo->SetPID(pid);
      break;
      
    case -321: // K- 
      pid[AliVCluster::kCharged] = 1.;
      calo->SetPID(pid);
      break;
      
    case 2112: // n
      pid[AliVCluster::kNeutron] = 1.;
      pid[AliVCluster::kNeutral] = 1.;
      calo->SetPID(pid);
      break;
      
    case 2212: // p
      pid[AliVCluster::kCharged] = 1.;
      calo->SetPID(pid);
      break;
      
    case -2212: // anti-p
      pid[AliVCluster::kCharged] = 1.;
      calo->SetPID(pid);
      break;
      
    case 310: // K0S
      pid[AliVCluster::kKaon0] = 1.;
      pid[AliVCluster::kNeutral] = 1.;
      calo->SetPID(pid);
      break;
      
    case 311: // K0
      pid[AliVCluster::kKaon0] = 1.;
      pid[AliVCluster::kNeutral] = 1.;
      calo->SetPID(pid);
      break;
      
    case -311: // anti-K0
      pid[AliVCluster::kKaon0] = 1.;
      pid[AliVCluster::kNeutral] = 1.;
      calo->SetPID(pid);
      break;
      
    case 221: // eta
      pid[AliVCluster::kNeutral] = 1.;
      calo->SetPID(pid);
      break;
      
    case 3122: // lambda
      pid[AliVCluster::kUnknown] = 1.;
      calo->SetPID(pid);
      break;
      
    case 3222: // Sigma+
      pid[AliVCluster::kUnknown] = 1.;
      calo->SetPID(pid);
      break;
      
    case 3212: // Sigma0
      pid[AliVCluster::kUnknown] = 1.;
      calo->SetPID(pid);
      break;
      
    case 3112: // Sigma-
      pid[AliVCluster::kUnknown] = 1.;
      calo->SetPID(pid);
      break;
      
    case 3322: // Xi0
      pid[AliVCluster::kUnknown] = 1.;
      calo->SetPID(pid);
      break;
      
    case 3312: // Xi-
      pid[AliVCluster::kUnknown] = 1.;
      calo->SetPID(pid);
      break;
      
    case 3334: // Omega-
      pid[AliVCluster::kUnknown] = 1.;
      calo->SetPID(pid);
      break;
      
    case -2112: // n-bar
      pid[AliVCluster::kNeutron] = 1.;
      pid[AliVCluster::kNeutral] = 1.;
      calo->SetPID(pid);
      break;
      
    case -3122: // anti-Lambda
      pid[AliVCluster::kUnknown] = 1.;
      calo->SetPID(pid);
      break;
      
    case -3222: // anti-Sigma-
      pid[AliVCluster::kUnknown] = 1.;
      calo->SetPID(pid);
      break;
      
    case -3212: // anti-Sigma0
      pid[AliVCluster::kUnknown] = 1.;
      calo->SetPID(pid);
      break;
      
    case -3112: // anti-Sigma+
      pid[AliVCluster::kUnknown] = 1.;
      calo->SetPID(pid);
      break;
      
    case -3322: // anti-Xi0
      pid[AliVCluster::kUnknown] = 1.;
      calo->SetPID(pid);
      break;
      
    case -3312: // anti-Xi+
      pid[AliVCluster::kUnknown] = 1.;
      calo->SetPID(pid);
      break;
      
    case -3334: // anti-Omega+
      pid[AliVCluster::kUnknown] = 1.;
      calo->SetPID(pid);
      break;
      
    case 411: // D+
      pid[AliVCluster::kUnknown] = 1.;
      calo->SetPID(pid);
      break;
      
    case -411: // D- 
      pid[AliVCluster::kUnknown] = 1.;
      calo->SetPID(pid);
      break;
      
    case 421: // D0
      pid[AliVCluster::kUnknown] = 1.;
      calo->SetPID(pid);
      break;
      
    case -421: // anti-D0
      pid[AliVCluster::kUnknown] = 1.;
      calo->SetPID(pid);
      break;
      
    default : // unknown
      pid[AliVCluster::kUnknown] = 1.;
      calo->SetPID(pid);
  }
  
  return;
}

