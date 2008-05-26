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

/* $Id:   AliAODParticleCorrelation.h $ */

//-------------------------------------------------------------------------
//     AOD class for photon and other particles storage and 
//     correlation studies
//     Author: Yves Schutz, CERN, Gustavo Conesa, INFN
//-------------------------------------------------------------------------

//-- ROOT system --

//-- Analysis system
#include "AliAODParticleCorrelation.h"

ClassImp(AliAODParticleCorrelation)


//______________________________________________________________________________
AliAODParticleCorrelation::AliAODParticleCorrelation() :
    AliVParticle(),
    fMomentum(0),fPdg(-1), fTag(-1),fLabel(-1),fDetector(""),fR(0),
    fRefTracks(new TRefArray()), fRefClusters(new TRefArray()),
    fRefIsolationConeTracks(new TRefArray()), fRefIsolationConeClusters(new TRefArray()),
    fRefBackgroundTracks(new TRefArray()), fRefBackgroundClusters(new TRefArray()),  
    fLeadingDetector(""), fLeading(), fCorrJet(),  fCorrBkg(), fRefJet(0)
{
  // constructor
}

//______________________________________________________________________________
AliAODParticleCorrelation::AliAODParticleCorrelation(Double_t px, Double_t py, Double_t pz, Double_t e):
    AliVParticle(),
    fMomentum(0),fPdg(-1), fTag(-1),fLabel(-1),fDetector(""),fR(0),
    fRefTracks(new TRefArray()), fRefClusters(new TRefArray()),
    fRefIsolationConeTracks(new TRefArray()), fRefIsolationConeClusters(new TRefArray()),
    fRefBackgroundTracks(new TRefArray()), fRefBackgroundClusters(new TRefArray()),
    fLeadingDetector(""),  fLeading(), fCorrJet(),
    fCorrBkg(), fRefJet(0)
{
  // constructor
    fMomentum = new TLorentzVector(px, py, pz, e);
}

//______________________________________________________________________________
AliAODParticleCorrelation::AliAODParticleCorrelation(TLorentzVector & p):
    AliVParticle(),
    fMomentum(0),fPdg(-1), fTag(-1),fLabel(-1),fDetector(""),fR(0),
    fRefTracks(new TRefArray()), fRefClusters(new TRefArray()),
    fRefIsolationConeTracks(new TRefArray()), fRefIsolationConeClusters(new TRefArray()),
    fRefBackgroundTracks(new TRefArray()), fRefBackgroundClusters(new TRefArray()),  
    fLeadingDetector(""),  fLeading(), fCorrJet(), fCorrBkg(),fRefJet(0)
{
  // constructor
    fMomentum = new TLorentzVector(p);
}


//______________________________________________________________________________
AliAODParticleCorrelation::~AliAODParticleCorrelation() 
{
  // destructor
    delete fMomentum;
    delete fRefTracks;
    delete fRefClusters;
    delete fRefIsolationConeTracks;
    delete fRefIsolationConeClusters;
    delete fRefBackgroundTracks;
    delete fRefBackgroundClusters;

}

//______________________________________________________________________________
AliAODParticleCorrelation::AliAODParticleCorrelation(const AliAODParticleCorrelation& part) :
    AliVParticle(part),
    fMomentum(0) ,fPdg(part.fPdg), fTag(part.fTag),fLabel(part.fLabel),
    fDetector(part.fDetector),
    fR(part.fR),
    fRefTracks(), fRefClusters(),
    fRefIsolationConeTracks(), fRefIsolationConeClusters(),
    fRefBackgroundTracks(), fRefBackgroundClusters(),   
    fLeadingDetector(part.fLeadingDetector), fLeading(part.fLeading),  
    fCorrJet(part.fCorrJet), fCorrBkg(part.fCorrBkg), fRefJet(part.fRefJet)
{
  // Copy constructor
  fMomentum = new TLorentzVector(*part.fMomentum);
  fRefTracks = new TRefArray(*part.fRefTracks);
  fRefClusters = new TRefArray(*part.fRefClusters);
  fRefIsolationConeTracks = new TRefArray(*part.fRefIsolationConeTracks);
  fRefIsolationConeClusters = new TRefArray(*part.fRefIsolationConeClusters);
  fRefBackgroundTracks = new TRefArray(*part.fRefBackgroundTracks);
  fRefBackgroundClusters = new TRefArray(*part.fRefBackgroundClusters);
}

//______________________________________________________________________________
AliAODParticleCorrelation& AliAODParticleCorrelation::operator=(const AliAODParticleCorrelation& part)
{
  // Assignment operator
  if(this!=&part) {
  }
  
  fPdg = part.fPdg;
  fTag = part.fTag;
  fLabel = part.fLabel;
  fR = part.fR;
  fRefJet = part.fRefJet ;
  fDetector =part.fDetector;
  fLeadingDetector =part.fLeadingDetector;
  fLeading = part.fLeading;
  fCorrJet = part.fCorrJet ;
  fCorrBkg = part.fCorrBkg;

  if (fMomentum ) delete fMomentum;
  if( fRefTracks ) delete fRefTracks ;
  if( fRefClusters) delete fRefClusters ;
  if( fRefIsolationConeTracks ) delete fRefIsolationConeTracks ;
  if( fRefIsolationConeClusters) delete fRefIsolationConeClusters ;
  if( fRefBackgroundTracks ) delete fRefBackgroundTracks ;
  if( fRefBackgroundClusters ) delete fRefBackgroundClusters ;

  fMomentum = new TLorentzVector(*part.fMomentum);
  fRefTracks = new TRefArray(*part.fRefTracks);
  fRefClusters = new TRefArray(*part.fRefClusters);
  fRefIsolationConeTracks = new TRefArray(*part.fRefIsolationConeTracks);
  fRefIsolationConeClusters = new TRefArray(*part.fRefIsolationConeClusters);
  fRefBackgroundTracks = new TRefArray(*part.fRefBackgroundTracks);
  fRefBackgroundClusters = new TRefArray(*part.fRefBackgroundClusters);  
  
  return *this;
}

void AliAODParticleCorrelation::Print(Option_t* /*option*/) const 
{
  // Print information of all data members
  printf("Particle 4-vector:\n");
  printf("     E  = %13.3f\n", E() );
  printf("     Px = %13.3f\n", Px());
  printf("     Py = %13.3f\n", Py());
  printf("     Pz = %13.3f\n", Pz());
  printf("pdg : %d\n",fPdg);
  printf("tag : %d\n",fTag);
  printf("R : %2.2f\n",fR);
  printf("Trigger Detector : %s\n",fDetector.Data());
  printf("Leading Detector : %s\n",fLeadingDetector.Data());
  // if(fRefJet) fRefJet.Print();

}
