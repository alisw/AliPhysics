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

/* $Id:   AliAODPWG4ParticleCorrelation.h $ */

//-------------------------------------------------------------------------
//     AOD class for photon and other particles storage and 
//     correlation studies
//     Author: Yves Schutz, CERN, Gustavo Conesa, INFN
//-------------------------------------------------------------------------

//-- ROOT system --

//-- Analysis system
#include "AliAODPWG4ParticleCorrelation.h"
#include "AliAODJet.h"

ClassImp(AliAODPWG4ParticleCorrelation)


//______________________________________________________________________________
AliAODPWG4ParticleCorrelation::AliAODPWG4ParticleCorrelation() :
    AliAODPWG4Particle(), fIsolated(kFALSE),
    fRefTracks(new TRefArray()), fRefClusters(new TRefArray()),
    fRefIsolationConeTracks(new TRefArray()), fRefIsolationConeClusters(new TRefArray()),
    fRefBackgroundTracks(new TRefArray()), fRefBackgroundClusters(new TRefArray()),  
    fLeadingDetector(""), fLeading(), fCorrJet(),  fCorrBkg(), fRefJet(0)
{
  // constructor
}

//______________________________________________________________________________
AliAODPWG4ParticleCorrelation::AliAODPWG4ParticleCorrelation(Double_t px, Double_t py, Double_t pz, Double_t e):
    AliAODPWG4Particle(), fIsolated(kFALSE),
    fRefTracks(new TRefArray()), fRefClusters(new TRefArray()),
    fRefIsolationConeTracks(new TRefArray()), fRefIsolationConeClusters(new TRefArray()),
    fRefBackgroundTracks(new TRefArray()), fRefBackgroundClusters(new TRefArray()),
    fLeadingDetector(""),  fLeading(), fCorrJet(),
    fCorrBkg(), fRefJet(0)
{
  // constructor
    SetMomentum(new TLorentzVector(px, py, pz, e));
}

//______________________________________________________________________________
AliAODPWG4ParticleCorrelation::AliAODPWG4ParticleCorrelation(TLorentzVector & p):
    AliAODPWG4Particle(p), fIsolated(kFALSE),
    fRefTracks(new TRefArray()), fRefClusters(new TRefArray()),
    fRefIsolationConeTracks(new TRefArray()), fRefIsolationConeClusters(new TRefArray()),
    fRefBackgroundTracks(new TRefArray()), fRefBackgroundClusters(new TRefArray()),  
    fLeadingDetector(""),  fLeading(), fCorrJet(), fCorrBkg(),fRefJet(0)
{
  // constructor
}

//______________________________________________________________________________
AliAODPWG4ParticleCorrelation::AliAODPWG4ParticleCorrelation(AliAODPWG4Particle & p):
    AliAODPWG4Particle(p), fIsolated(kFALSE),
    fRefTracks(new TRefArray()), fRefClusters(new TRefArray()),
    fRefIsolationConeTracks(new TRefArray()), fRefIsolationConeClusters(new TRefArray()),
    fRefBackgroundTracks(new TRefArray()), fRefBackgroundClusters(new TRefArray()),  
    fLeadingDetector(""),  fLeading(), fCorrJet(), fCorrBkg(),fRefJet(0)
{
  // constructor

}

//______________________________________________________________________________
AliAODPWG4ParticleCorrelation::~AliAODPWG4ParticleCorrelation() 
{
  // destructor
    delete fRefTracks;
    delete fRefClusters;
    delete fRefIsolationConeTracks;
    delete fRefIsolationConeClusters;
    delete fRefBackgroundTracks;
    delete fRefBackgroundClusters;

}

//______________________________________________________________________________
AliAODPWG4ParticleCorrelation::AliAODPWG4ParticleCorrelation(const AliAODPWG4ParticleCorrelation& part) :
    AliAODPWG4Particle(part), fIsolated(part.fIsolated),
    fRefTracks(), fRefClusters(),
    fRefIsolationConeTracks(), fRefIsolationConeClusters(),
    fRefBackgroundTracks(), fRefBackgroundClusters(),   
    fLeadingDetector(part.fLeadingDetector), fLeading(part.fLeading),  
    fCorrJet(part.fCorrJet), fCorrBkg(part.fCorrBkg), fRefJet(part.fRefJet)
{
  // Copy constructor
  fRefTracks                = new TRefArray(*part.fRefTracks);
  fRefClusters              = new TRefArray(*part.fRefClusters);
  fRefIsolationConeTracks   = new TRefArray(*part.fRefIsolationConeTracks);
  fRefIsolationConeClusters = new TRefArray(*part.fRefIsolationConeClusters);
  fRefBackgroundTracks      = new TRefArray(*part.fRefBackgroundTracks);
  fRefBackgroundClusters    = new TRefArray(*part.fRefBackgroundClusters);
}

//______________________________________________________________________________
AliAODPWG4ParticleCorrelation& AliAODPWG4ParticleCorrelation::operator=(const AliAODPWG4ParticleCorrelation& part)
{
  // Assignment operator
  if(this!=&part) {
  
	fIsolated = part.fIsolated;
	fRefJet   = part.fRefJet ;
	fLeadingDetector =part.fLeadingDetector;
	fLeading  = part.fLeading;
	fCorrJet  = part.fCorrJet ;
	fCorrBkg  = part.fCorrBkg;

	if( fRefTracks )               delete fRefTracks ;
	if( fRefClusters)              delete fRefClusters ;
	if( fRefIsolationConeTracks )  delete fRefIsolationConeTracks ;
	if( fRefIsolationConeClusters) delete fRefIsolationConeClusters ;
	if( fRefBackgroundTracks )     delete fRefBackgroundTracks ;
	if( fRefBackgroundClusters )   delete fRefBackgroundClusters ;

	fRefTracks                = new TRefArray(*part.fRefTracks);
	fRefClusters              = new TRefArray(*part.fRefClusters);
	fRefIsolationConeTracks   = new TRefArray(*part.fRefIsolationConeTracks);
	fRefIsolationConeClusters = new TRefArray(*part.fRefIsolationConeClusters);
	fRefBackgroundTracks      = new TRefArray(*part.fRefBackgroundTracks);
	fRefBackgroundClusters    = new TRefArray(*part.fRefBackgroundClusters);  
  
  }
  
  return *this;
}

//______________________________________________________________________________
void AliAODPWG4ParticleCorrelation::Print(Option_t* /*option*/) const 
{
  // Print information of all data members
  AliAODPWG4Particle::Print("");
  if(fIsolated) printf("Isolated! \n");
  printf("Leading Detector : %s\n",fLeadingDetector.Data());
  // if(fRefJet) fRefJet.Print();

}
