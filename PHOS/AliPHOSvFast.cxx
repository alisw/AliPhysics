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
 * Revision 1.30  2006/09/13 07:31:01  kharlov
 * Effective C++ corrections (T.Pocheptsov)
 *
 * Revision 1.29  2005/05/28 14:19:05  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
// Implementation of the PHOS manager class for fast simulations     
// Tracks particles until the reach a grossly designed PHOS module
// Modify the particles property (momentum, energy, type) according to
//  the PHOS response function. The result is called a virtual reconstructed
//  particle.                
//
//*-- Author: Yves Schutz (SUBATECH)

// --- ROOT system ---
 
#include <TBRIK.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TParticle.h>
#include "TClonesArray.h" 
#include <TVirtualMC.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliPHOSFastRecParticle.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSLoader.h"
#include "AliPHOSvFast.h"
#include "AliRun.h"

ClassImp(AliPHOSvFast)

AliPHOSvFast::AliPHOSvFast() :
  fBigBoxX(0.),
  fBigBoxY(0.),
  fBigBoxZ(0.),
  fFastRecParticles(0),
  fNRecParticles(0),
  fRan(0),
  fResPara1(0.),
  fResPara2(0.),
  fResPara3(0.),
  fPosParaA0(0.),
  fPosParaA1(0.),
  fPosParaB0(0.),
  fPosParaB1(0.),
  fPosParaB2(0.)    
{
  // default ctor : initialize data member
}

//____________________________________________________________________________
AliPHOSvFast::AliPHOSvFast(const char *name, const char *title):
  AliPHOS(name,title),
  fBigBoxX(0.),
  fBigBoxY(0.),
  fBigBoxZ(0.),
  fFastRecParticles(new AliPHOSFastRecParticle::FastRecParticlesList("AliPHOSFastRecParticle", 100)),
  fNRecParticles(0),
  fRan(0),
  fResPara1(0.030), // GeV
  fResPara2(0.00003),
  fResPara3(0.00001),
  fPosParaA0(2.87),    // mm
  fPosParaA1(-0.0975),
  fPosParaB0(0.257),
  fPosParaB1(0.137),
  fPosParaB2(0.00619)
{
  // ctor
  // create the Loader 
  SetBigBox(0, GetGeometry()->GetOuterBoxSize(0) ) ;
  SetBigBox(1, GetGeometry()->GetOuterBoxSize(3) + GetGeometry()->GetCPVBoxSize(1) ) ; 
  SetBigBox(2, GetGeometry()->GetOuterBoxSize(2) ); 
}

//____________________________________________________________________________
AliPHOSvFast::~AliPHOSvFast()
{
  // dtor
 
  fFastRecParticles->Delete() ; 
  delete fFastRecParticles ;
  fFastRecParticles = 0 ; 

}

//____________________________________________________________________________
void AliPHOSvFast::AddRecParticle(const AliPHOSFastRecParticle & rp)
{  
  // Add a virtually reconstructed particle to the list 

  new( (*fFastRecParticles)[fNRecParticles] ) AliPHOSFastRecParticle(rp) ;
  fNRecParticles++ ; 
}

//____________________________________________________________________________
void AliPHOSvFast::BuildGeometry()
{
  // Build the PHOS geometry for the ROOT display
   //BEGIN_HTML
  /*
    <H2>
     PHOS FAST in ALICE displayed by root
    </H2>
    <H4> All Views </H4>
    <P>
    <CENTER>
    <IMG Align=BOTTOM ALT="Fast All Views" SRC="../images/AliPHOSvFastAllViews.gif"> 
    </CENTER></P>
    <H4> Front View </H4>
    <P>
    <CENTER>
    <IMG Align=BOTTOM ALT="Fast Front View" SRC="../images/AliPHOSvFastFrontView.gif"> 
    </CENTER></P>
  */
  //END_HTML  
 
  const Int_t kColorPHOS = kRed ;
  
  Double_t const kRADDEG = 180.0 / TMath::Pi() ;
  
  new TBRIK( "BigBox", "PHOS box", "void", GetBigBox(0)/2, 
	     GetBigBox(1)/2, 
	     GetBigBox(2)/2 );
  
  // position PHOS into ALICE

  Float_t r = GetGeometry()->GetIPtoCrystalSurface() + GetBigBox(1) / 2.0 ;
  Int_t number = 988 ; 
  Float_t pphi =  TMath::ATan( GetBigBox(0)  / ( 2.0 * GetGeometry()->GetIPtoCrystalSurface() ) ) ;
  pphi *= kRADDEG ;
  TNode * top = gAlice->GetGeometry()->GetNode("alice") ;
 
  char * nodename = new char[20] ;  
  char * rotname  = new char[20] ; 

  for( Int_t i = 1; i <= GetGeometry()->GetNModules(); i++ ) { 
   Float_t angle = pphi * 2 * ( i - GetGeometry()->GetNModules() / 2.0 - 0.5 ) ;
   sprintf(rotname, "%s%d", "rot", number++) ;
   new TRotMatrix(rotname, rotname, 90, angle, 90, 90 + angle, 0, 0);
   top->cd();
   sprintf(nodename,"%s%d", "Module", i) ;    
   Float_t x =  r * TMath::Sin( angle / kRADDEG ) ;
   Float_t y = -r * TMath::Cos( angle / kRADDEG ) ;
   TNode * bigboxnode = new TNode(nodename, nodename, "BigBox", x, y, 0, rotname ) ;
   bigboxnode->SetLineColor(kColorPHOS) ;
   fNodes->Add(bigboxnode) ;
  }
  delete[] nodename ; 
  delete[] rotname ; 
}

//____________________________________________________________________________
void AliPHOSvFast::CreateGeometry()
{
  // Create the geometry for GEANT
  
  AliPHOSvFast *phostmp = (AliPHOSvFast*)gAlice->GetModule("PHOS") ;
  
  if ( phostmp == NULL ) {
    
    fprintf(stderr, "PHOS detector not found!\n") ;
    return ;
    
  }

  // Get pointer to the array containing media indeces
  Int_t *idtmed = fIdtmed->GetArray() - 699 ;
  
  Float_t bigbox[3] ; 
  bigbox[0] =   GetBigBox(0) / 2.0 ;
  bigbox[1] =   GetBigBox(1) / 2.0 ;
  bigbox[2] =   GetBigBox(2) / 2.0 ;
  
  gMC->Gsvolu("PHOS", "BOX ", idtmed[798], bigbox, 3) ;
  
  // --- Position  PHOS mdules in ALICE setup ---
  
  Int_t idrotm[99] ;
  Double_t const kRADDEG = 180.0 / TMath::Pi() ;
  
  for( Int_t i = 1; i <= GetGeometry()->GetNModules(); i++ ) {
    
    Float_t angle = GetGeometry()->GetPHOSAngle(i) ;
    AliMatrix(idrotm[i-1], 90.0, angle, 90.0, 90.0+angle, 0.0, 0.0) ;
 
    Float_t r = GetGeometry()->GetIPtoCrystalSurface() + GetBigBox(1) / 2.0 ;

    Float_t xP1 = r * TMath::Sin( angle / kRADDEG ) ;
    Float_t yP1 = -r * TMath::Cos( angle / kRADDEG ) ;
    gMC->Gspos("PHOS", i, "ALIC", xP1, yP1, 0.0, idrotm[i-1], "ONLY") ;
 
  } // for GetNModules

}


//____________________________________________________________________________
void AliPHOSvFast::Init(void)
{
  // Prints out an information message
  
  Int_t i;

  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" FAST PHOS_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");

  // Here the PHOS initialisation code (if any!)

  for(i=0;i<80;i++) printf("*");
  printf("\n");
  
}

//___________________________________________________________________________
Float_t AliPHOSvFast::GetBigBox(Int_t index) const
{
  // Get the X, Y or Z dimension of the box describing a PHOS module
  
  Float_t rv = 0 ; 

  switch (index) {
  case 0:
    rv = fBigBoxX ; 
    break ; 
  case 1:
     rv = fBigBoxY ; 
    break ; 
  case 2:
     rv = fBigBoxZ ; 
    break ; 
 }
  return rv ; 
}
//___________________________________________________________________________

void AliPHOSvFast::MakeBranch(Option_t* opt)
{  
  // Create new branch in the current reconstructed Root Tree
  AliDetector::MakeBranch(opt);
  const char *cd = strstr(opt,"R");
  
  if (fFastRecParticles && fLoader->TreeR() && cd) {
    MakeBranchInTree(fLoader->TreeR(), GetName(), &fFastRecParticles, fBufferSize, 0);
  }
}
//____________________________________________________________________________

Double_t AliPHOSvFast::MakeEnergy(Double_t energy)
{  
  // Smears the energy according to the energy dependent energy resolution.
  // A gaussian distribution is assumed

  Double_t sigma  = SigmaE(energy) ; 
  return  fRan.Gaus(energy, sigma) ;   
}
//____________________________________________________________________________

TVector3 AliPHOSvFast::MakePosition(Double_t energy, TVector3 pos, Double_t theta, Double_t phi)
{
  // Smears the impact position according to the energy dependent position resolution
  // A gaussian position distribution is assumed

  TVector3 newpos ;
  Double_t sigma = SigmaP( energy, theta*180./TMath::Pi() ) ;
  Double_t x = fRan.Gaus( pos.X(), sigma ) ;
  sigma = SigmaP( energy, phi*180./TMath::Pi() ) ;
  Double_t z = fRan.Gaus( pos.Z(), sigma ) ;
  Double_t y = pos.Y() ; 
  
  newpos.SetX(x) ; 
  newpos.SetY(y) ; 
  newpos.SetZ(z) ; 
	      
  return newpos ; 
}

//____________________________________________________________________________
void AliPHOSvFast::MakeRecParticle(Int_t modid, TVector3 pos, AliPHOSFastRecParticle & rp)
{
  // Modify the primary particle properties according
  //  1. the response function of PHOS
  //  2. the performance of the EMC+PPSD setup
  
  Int_t type = MakeType( rp ) ;
  rp.SetType(type) ;

  
  // get the detected energy

  TLorentzVector momentum ;  
  rp.Momentum(momentum) ; 
  Double_t kineticenergy = TMath::Sqrt( TMath::Power(momentum.E(), 2) - TMath::Power(rp.GetMass(), 2) ) ; 
  Double_t modifiedkineticenergy = MakeEnergy(kineticenergy ) ;
  Double_t modifiedenergy = TMath::Sqrt( TMath::Power(modifiedkineticenergy, 2)  
					 + TMath::Power( rp.GetMass(), 2) ) ;
 
  // get the angle of incidence 
  
  Double_t incidencetheta = 90. * TMath::Pi() /180 - rp.Theta() ; 
  Double_t incidencephi   = ( 270 + GetGeometry()->GetPHOSAngle(modid) ) * TMath::Pi() / 180. - rp.Phi() ;   

  // get the detected direction
  
  TVector3 modifiedposition = MakePosition(kineticenergy, pos, incidencetheta, incidencephi) ; 
  modifiedposition *= modifiedkineticenergy / modifiedposition.Mag() ; 

  // Set the modified 4-momentum of the reconstructed particle

  rp.SetMomentum(modifiedposition.X(), modifiedposition.Y(), modifiedposition.Z(), modifiedenergy) ; 

 }

//____________________________________________________________________________
Int_t AliPHOSvFast::MakeType(AliPHOSFastRecParticle & rp )
{
  // Generate a particle type using the performance of the EMC+PPSD setup

  Int_t rv =   AliPHOSFastRecParticle::kUNDEFINED ;
  Int_t charge = (Int_t)rp.GetPDG()->Charge() ;
  Int_t test ; 
  Float_t ran ; 
  if ( charge != 0 && ( TMath::Abs(rp.GetPdgCode()) != 11 ) ) 
    test = - 1 ;
  else
    test = rp.GetPdgCode() ; 

  Fatal("MakeType", "SHOULD NOT BE USED until values of probabilities are properly set ") ;
   // NB: ALL VALUES SHOULD BE CHECKED !!!!
  switch (test) { 

  case 22:    // it's a photon              // NB: ALL VALUES SHOLD BE CHECKED !!!!
    ran = fRan.Rndm() ; 
    if( ran <= 0.9498 )
      rv =  AliPHOSFastRecParticle::kNEUTRALHAFAST ; 
    else
      rv =  AliPHOSFastRecParticle::kNEUTRALEMFAST ;     
    break ; 

  case 2112:  // it's a neutron
    ran = fRan.Rndm() ; 
    if ( ran <= 0.9998 )
      rv =  AliPHOSFastRecParticle::kNEUTRALHASLOW ; 
    else 
      rv = AliPHOSFastRecParticle::kNEUTRALEMSLOW ; 
    break ; 
    
  case -2112: // it's a anti-neutron
    ran = fRan.Rndm() ; 
    if ( ran <= 0.9984 )
      rv =  AliPHOSFastRecParticle::kNEUTRALHASLOW ; 
    else 
      rv =  AliPHOSFastRecParticle::kNEUTRALEMSLOW ; 
    break ; 
    
  case 11:    // it's a electron
    ran = fRan.Rndm() ; 
    if ( ran <= 0.9996 )
      rv =  AliPHOSFastRecParticle::kCHARGEDEMFAST ; 
    else 
      rv =  AliPHOSFastRecParticle::kCHARGEDHAFAST ; 
    break; 

  case -11:   // it's a positon
    ran = fRan.Rndm() ; 
    if ( ran <= 0.9996 )
      rv =  AliPHOSFastRecParticle::kCHARGEDEMFAST ; 
    else 
      rv =  AliPHOSFastRecParticle::kCHARGEDHAFAST ; 
    break; 

  case -1:    // it's a charged
    ran = fRan.Rndm() ; 
    if ( ran <= 0.9996 )
      rv =  AliPHOSFastRecParticle::kCHARGEDHAFAST ; 
    else 
      rv =  AliPHOSFastRecParticle::kNEUTRALHAFAST ; 

    break ; 
  }
    
  
  return rv ;
}

//___________________________________________________________________________
void AliPHOSvFast::ResetPoints()
{
  // This overloads the method in AliDetector
  
  ResetFastRecParticles() ; 
}

//___________________________________________________________________________
void AliPHOSvFast::ResetFastRecParticles()
{
  // Resets the list of virtual reconstructed particles
 
  if (fFastRecParticles) 
    fFastRecParticles->Clear() ;
  fNRecParticles = 0 ; 
}

//___________________________________________________________________________
void AliPHOSvFast::SetBigBox(Int_t index, Float_t value)
{
  // Set the size of the Box describing a PHOS module
  
  switch (index) {
  case 0:
    fBigBoxX = value ; 
    break ; 
  case 1:
    fBigBoxY = value ; 
    break ; 
  case 2:
    fBigBoxZ = value ; 
    break ; 
 }

}

//____________________________________________________________________________
Double_t AliPHOSvFast::SigmaE(Double_t energy)
{
  // Calculates the energy dependent energy resolution
  
  Double_t rv = -1 ; 
  
  rv = TMath::Sqrt( TMath::Power(fResPara1/energy, 2) 
	       + TMath::Power(fResPara2/TMath::Sqrt(energy), 2) 
	       + TMath::Power(fResPara3, 2) ) ;  

  return rv * energy ; 
}

//____________________________________________________________________________
Double_t AliPHOSvFast::SigmaP(Double_t energy, Double_t incidence)
{
  // Calculates the energy dependent position resolution 

  Double_t paraA = fPosParaA0 + fPosParaA1 * incidence ; 
  Double_t paraB = fPosParaB0 + fPosParaB1 * incidence + fPosParaB2 * incidence * incidence ; 

  return ( paraA / TMath::Sqrt(energy) + paraB ) * 0.1   ; // in cm  
}

//____________________________________________________________________________
void AliPHOSvFast::StepManager(void)
{
  // Only verifies if the particle reaches PHOS and stops the tracking 

  TLorentzVector lv ; 
  gMC->TrackPosition(lv) ;
  TVector3 pos = lv.Vect() ; 
  Int_t modid  ; 
  gMC->CurrentVolID(modid);
  
  Float_t energy = gMC->Etot() ; //Total energy of current track

  //Calculating mass of current particle
  TDatabasePDG * pdg = TDatabasePDG::Instance() ;
  TParticlePDG * partPDG = pdg->GetParticle(gMC->TrackPid()) ;
  Float_t mass = partPDG->Mass() ;

  if(energy > mass){
    pos.SetMag(TMath::Sqrt(energy*energy-mass*mass)) ;
    TLorentzVector pTrack(pos, energy) ;  

    TParticle * part = new TParticle(gMC->TrackPid(), 0,-1,-1,-1,-1, pTrack, lv)  ;
        
    AliPHOSFastRecParticle rp(*part) ;

    // Adds the response of PHOS to the particle
    MakeRecParticle(modid, pos, rp) ;
    
    // add the `track' particle to the FastRecParticles list
  
    AddRecParticle(rp) ;

    part->Delete() ;
  }
  // stop the track as soon PHOS is reached
  
  gMC->StopTrack() ; 

}

