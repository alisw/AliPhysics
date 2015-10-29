// $Id$
//
// Jet embedding task.
//
// Author: S.Aiola, C.Loizides

#include <TRandom3.h>
#include <TFile.h>
#include <AliLog.h>

#include "AliJetEmbeddingTask.h"

ClassImp(AliJetEmbeddingTask)

//________________________________________________________________________
AliJetEmbeddingTask::AliJetEmbeddingTask() : 
  AliJetModelBaseTask("AliJetEmbeddingTask"),
  fMassless(kFALSE),
  fMassFromDistr(kFALSE),
  fNeutralFraction(0),
  fNeutralMass(0.135),
  fMass(0.1396),
  fHMassDistrib(0)
{
  // Default constructor.
  SetSuffix("Embedded");
}

//________________________________________________________________________
AliJetEmbeddingTask::AliJetEmbeddingTask(const char *name) : 
  AliJetModelBaseTask(name),
  fMassless(kFALSE),
  fMassFromDistr(kFALSE),
  fNeutralFraction(0),
  fNeutralMass(0.135),
  fMass(0.1396),
  fHMassDistrib(0)
{
  // Standard constructor.
  SetSuffix("Embedded");
}

//________________________________________________________________________
AliJetEmbeddingTask::~AliJetEmbeddingTask()
{
  // Destructor
}

//________________________________________________________________________
void AliJetEmbeddingTask::Run() 
{
  // Embed particles.
  
  if (fNClusters > 0 && fOutClusters) {
    if (fCopyArray) 
      CopyClusters();
    for (Int_t i = 0; i < fNClusters; ++i) {
      AddCluster();
    }
  }
 
  if (fNTracks > 0 && fOutTracks) {
    if (fCopyArray) 
      CopyTracks();
    for (Int_t i = 0; i < fNTracks; ++i) {
      Double_t mass = fMass;
      Short_t charge = 1;
      if(fNeutralFraction>0.) {
	Double_t rnd = gRandom->Rndm();
	if(rnd<fNeutralFraction) {
	  charge = 0;
	  mass = fNeutralMass;
	}
      }
      if(fMassless) mass = 0.;
      if(fMassFromDistr) {
      	 if(fHMassDistrib)
      	 mass = fHMassDistrib->GetRandom();
      	 else {
      	    AliError(Form("Template distribution for mass of track embedding not found, use %f", fMass));
      	    mass = fMass;
      	 }
      }
      AddTrack(-1,-999,-1,0,0,0,0,kFALSE,0,charge,mass);
    }
  }
}

//________________________________________________________________________
void AliJetEmbeddingTask::SetMassDistributionFromFile(TString filename, TString histoname){
   TFile *f = TFile::Open(filename);
   if(!f){
      AliError(Form("File %s not found, cannot SetMassDistribution", filename.Data()));
      return;
   }
   
   TH1D* h = (TH1D*) f->Get(histoname);
   SetMassDistribution(h);
   return;

}

//________________________________________________________________________

void AliJetEmbeddingTask::SetMassAndPtDistributionFromFile(TString filenameM, TString filenamepT, TString histonameM, TString histonamepT){
   SetMassDistributionFromFile(filenameM, histonameM);
   SetpTDistributionFromFile(filenamepT, histonamepT);
   return;
}

//________________________________________________________________________
void AliJetEmbeddingTask::SetpTDistributionFromFile(TString filename, TString histoname){
   TFile *f = TFile::Open(filename);
   if(!f){
      AliError(Form("File %s not found, cannot SetpTDistribution", filename.Data()));
      return;
   }
   
   TH1D* h = (TH1D*) f->Get(histoname);
   AliJetModelBaseTask::SetPtSpectrum((TH1F*)h);
   return;

}
