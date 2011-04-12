/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id:   AliAODPWG4Particle.h $ */

//-------------------------------------------------------------------------
//     AOD class for photon and other particles storage and 
//     correlation studies
//     Author: Yves Schutz, CERN, Gustavo Conesa, INFN
//-------------------------------------------------------------------------

//-- ROOT system --

//-- Analysis system
#include "AliAODPWG4Particle.h"

ClassImp(AliAODPWG4Particle)


//______________________________________________________________________________
AliAODPWG4Particle::AliAODPWG4Particle() :
AliVParticle(),
fMomentum(0),fPdg(-1), fTag(0), fBtag(-1), fLabel(-1), fCaloLabel(), fTrackLabel(),
fDetector(""), fDisp(0), fTof(0), fCharged(0), fTagged(0), fBadDist(0), fFidArea(0), fInputFileIndex(0)
{
  // constructor
  fCaloLabel [0] = -1;
  fCaloLabel [1] = -1;
  fTrackLabel[0] = -1;
  fTrackLabel[1] = -1;
  fTrackLabel[2] = -1;
  fTrackLabel[3] = -1;
}

//______________________________________________________________________________
AliAODPWG4Particle::AliAODPWG4Particle(Double_t px, Double_t py, Double_t pz, Double_t e):
  AliVParticle(),
  fMomentum(0),fPdg(-1), fTag(0), fBtag(-1), fLabel(-1),fCaloLabel(), fTrackLabel(),
  fDetector(""), fDisp(0), fTof(0),fCharged(0), fTagged(0), fBadDist(0), fFidArea(0), fInputFileIndex(0)
{
  // constructor
  fMomentum = new TLorentzVector(px, py, pz, e);
  
  fCaloLabel [0] = -1;
  fCaloLabel [1] = -1;
  fTrackLabel[0] = -1;
  fTrackLabel[1] = -1;	
  fTrackLabel[2] = -1;
  fTrackLabel[3] = -1;	
}

//______________________________________________________________________________
AliAODPWG4Particle::AliAODPWG4Particle(TLorentzVector & p):
  AliVParticle(),
  fMomentum(0),fPdg(-1), fTag(0), fBtag(-1), fLabel(-1),fCaloLabel(), fTrackLabel(),
  fDetector(""), fDisp(0), fTof(0), fCharged(0), fTagged(0), fBadDist(0), fFidArea(0), fInputFileIndex(0)
{
  // constructor
  fMomentum = new TLorentzVector(p);
  
  fCaloLabel [0] = -1;
  fCaloLabel [1] = -1;
  fTrackLabel[0] = -1;
  fTrackLabel[1] = -1;
  fTrackLabel[2] = -1;
  fTrackLabel[3] = -1;
}


//______________________________________________________________________________
AliAODPWG4Particle::~AliAODPWG4Particle() 
{
  // destructor
    delete fMomentum;
}

//______________________________________________________________________________
void AliAODPWG4Particle::Clear(const Option_t* /*opt*/) 
{
  //clear
  delete fMomentum;
}

//______________________________________________________________________________
AliAODPWG4Particle::AliAODPWG4Particle(const AliAODPWG4Particle& part) :
  AliVParticle(part),
  fMomentum(0), fPdg(part.fPdg), fTag(part.fTag), fBtag(part.fBtag), fLabel(part.fLabel), 
  fCaloLabel(), fTrackLabel(), fDetector(part.fDetector),fDisp(part.fDisp), 
  fTof(part.fTof), fCharged(part.fCharged), fTagged(part.fTagged), fBadDist(part.fBadDist), 
  fFidArea(part.fFidArea), fInputFileIndex(part.fInputFileIndex)
{
  // Copy constructor
  fMomentum = new TLorentzVector(*part.fMomentum);
  
  fCaloLabel [0] = part.fCaloLabel[0];
  fCaloLabel [1] = part.fCaloLabel[1];
  fTrackLabel[0] = part.fTrackLabel[0];
  fTrackLabel[1] = part.fTrackLabel[1];
  fTrackLabel[2] = part.fTrackLabel[2];
  fTrackLabel[3] = part.fTrackLabel[3];
}

//______________________________________________________________________________
AliAODPWG4Particle& AliAODPWG4Particle::operator=(const AliAODPWG4Particle& part)
{
  // Assignment operator
  if(this!=&part) {
    
    fPdg   = part.fPdg;
    fTag   = part.fTag;
    fBtag  = part.fBtag;
    fLabel = part.fLabel;
	
    fCaloLabel [0] = part.fCaloLabel[0];
    fCaloLabel [1] = part.fCaloLabel[1];
    fTrackLabel[0] = part.fTrackLabel[0];
    fTrackLabel[1] = part.fTrackLabel[1];
	
    fDetector = part.fDetector;
    fDisp     = part.fDisp; 
    fTof      = part.fTof; 
    fCharged  = part.fCharged; 
	fTagged   = part.fTagged;
    fBadDist  = part.fBadDist;
	fFidArea  = part.fFidArea;
	fInputFileIndex =  part.fInputFileIndex;
	  
    if (fMomentum ) delete fMomentum;	
    fMomentum = new TLorentzVector(*part.fMomentum);
  }
  
  return *this;
}


//_______________________________________________________________
Bool_t AliAODPWG4Particle::IsPIDOK(const Int_t ipid, const Int_t pdgwanted) const{
  // returns true if particle satisfies given PID criterium
	switch(ipid){
	case 0: return kTRUE ; //No PID at all
	case 1: 
	  {
	    if (fPdg == pdgwanted) return kTRUE; 
	    else return kFALSE; //Overall PID calculated with bayesian methods.
	  }
	case 2: return fDisp ;   //only dispersion cut
	case 3: return fTof ;    //Only TOF cut
	case 4: return fCharged ;    //Only Charged cut
	case 5: return fDisp && fTof ;  //Dispersion and TOF
	case 6: return fDisp && fCharged ;  //Dispersion and Charged
	case 7: return fTof  && fCharged ;  //TOF and Charged
	case 8: return fDisp && fTof && fCharged ; // all 3 cuts
	default: return kFALSE ; //Not known combination
	}
}

//______________________________________________________________________________
void AliAODPWG4Particle::Print(Option_t* /*option*/) const 
{
  // Print information of all data members
  printf("Particle 4-vector:\n");
  printf("     E  = %13.3f", E() );
  printf("     Px = %13.3f", Px());
  printf("     Py = %13.3f", Py());
  printf("     Pz = %13.3f\n", Pz());
  printf("PID bits :\n");
  printf("     TOF        : %d",fTof);
  printf("     Charged    : %d",fCharged);
  printf("     Dispersion : %d\n",fDisp);
  printf("PDG       : %d\n",fPdg);
  printf("Tag       : %d\n",fTag); 
  printf("Btag      : %d\n",fBtag);  
  printf("Pi0 Tag   : %d\n",fTagged);  
  printf("Dist. to bad channel : %d\n",fBadDist);  
  printf("Fid Area  : %d\n",fFidArea);  
  printf("Input File Index : %d\n",fInputFileIndex);  
  printf("Detector  : %s\n",fDetector.Data());
  
}
