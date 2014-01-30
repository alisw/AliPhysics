// $Id$
//
// Copy tags from particle level constituent to detector level
//
// Author: S. Aiola

#include "AliJetConstituentTagCopier.h"

#include <TClonesArray.h>
#include <TMath.h>
#include <TLorentzVector.h>

#include "AliNamedArrayI.h"
#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliLog.h"

ClassImp(AliJetConstituentTagCopier)

//________________________________________________________________________
AliJetConstituentTagCopier::AliJetConstituentTagCopier() : 
  AliAnalysisTaskEmcal("AliJetConstituentTagCopier", kFALSE),
  fCleanBeforeCopy(kFALSE),
  fMCLabelShift(0),
  fMCParticleContainer(0)
{
  // Default constructor.
}

//________________________________________________________________________
AliJetConstituentTagCopier::AliJetConstituentTagCopier(const char *name) : 
  AliAnalysisTaskEmcal(name, kFALSE),
  fCleanBeforeCopy(kFALSE),
  fMCLabelShift(0),
  fMCParticleContainer(0)
{
  // Standard constructor.
}

//________________________________________________________________________
AliJetConstituentTagCopier::~AliJetConstituentTagCopier()
{
  // Destructor
}

//________________________________________________________________________
Bool_t AliJetConstituentTagCopier::Run()
{
  for (Int_t i = 0; i < fParticleCollArray.GetEntriesFast(); i++) {
    AliParticleContainer *cont = static_cast<AliParticleContainer*>(fParticleCollArray.At(i));
    if (!cont) continue;
    if (cont == fMCParticleContainer) continue;
    DoParticleLoop(cont);
  }

  for (Int_t i = 0; i < fClusterCollArray.GetEntriesFast(); i++) {
    AliClusterContainer *cont = static_cast<AliClusterContainer*>(fClusterCollArray.At(i));
    if (!cont) continue;
    DoClusterLoop(cont);
  }

  return kTRUE;
}

//________________________________________________________________________
void AliJetConstituentTagCopier::DoClusterLoop(AliClusterContainer *cont)
{
  AliVCluster *cluster = 0;

  if (fCleanBeforeCopy) {
    cont->ResetCurrentID();
    while ((cluster = static_cast<AliVCluster*>(cont->GetNextAcceptCluster()))) {
      Int_t mcLabel = cluster->GetLabel();
      if (mcLabel > 0) cluster->SetBit(TObject::kBitMask, kFALSE);
    }
  }

  if (!fMCParticleContainer) return;
  
  Double_t totalEnergy = 0;
  cont->ResetCurrentID();
  while ((cluster = static_cast<AliVCluster*>(cont->GetNextAcceptCluster()))) {
    Int_t mcLabel = cluster->GetLabel();
    if (mcLabel > fMCLabelShift) mcLabel -= fMCLabelShift;
    if (mcLabel > 0) {
      TLorentzVector vect;
      cluster->GetMomentum(vect, fVertex);
      AliDebug(2, Form("Cluster %d, pt = %f, eta = %f, phi = %f, label = %d",
		       cont->GetCurrentID(), cluster->E(), vect.Eta(), vect.Phi(), mcLabel));
      totalEnergy += cluster->E();
      Int_t index = fMCParticleContainer->GetIndexFromLabel(mcLabel);
      if (index < 0) continue;
      AliVParticle *part = fMCParticleContainer->GetParticle(index);
      if (!part) {
	AliError(Form("%s: Could not get MC particle %d", GetName(), index));
	continue;
      }      
      AliDebug(2, Form("Matched with particle %d, pt = %f, eta = %f, phi = %f", 
		       index, part->E(), part->Eta(), part->Phi()));
      UInt_t bits = (UInt_t)part->TestBits(TObject::kBitMask);
      cluster->SetBit(bits);
    }
  }

  AliDebug(2, Form("Total energy of MC clusters = %f", totalEnergy));
}

//________________________________________________________________________
void AliJetConstituentTagCopier::DoParticleLoop(AliParticleContainer *cont)
{
  AliVParticle *track = 0;

  if (fCleanBeforeCopy) {
    cont->ResetCurrentID();
    while ((track = static_cast<AliVParticle*>(cont->GetNextAcceptParticle()))) {
      Int_t mcLabel = TMath::Abs(track->GetLabel());
      if (mcLabel > 0) track->SetBit(TObject::kBitMask, kFALSE);
    }
  }

  if (!fMCParticleContainer) return;

  cont->ResetCurrentID();
  while ((track = static_cast<AliVParticle*>(cont->GetNextAcceptParticle()))) {
    Int_t mcLabel = TMath::Abs(track->GetLabel());
    if (mcLabel > fMCLabelShift) mcLabel -= fMCLabelShift;
    if (mcLabel > 0) {
      Int_t index = fMCParticleContainer->GetIndexFromLabel(mcLabel);
      if (index < 0) continue;
      AliVParticle *part = fMCParticleContainer->GetParticle(index);
      if (!part) {
	AliError(Form("%s: Could not get MC particle %d", GetName(), index));
	continue;
      }
      AliDebug(3, Form("Track %d, pt = %f, eta = %f, phi = %f, label = %d is matched with particle %d, pt = %f, eta = %f, phi = %f", 
		       cont->GetCurrentID(), track->Pt(), track->Eta(), track->Phi(), mcLabel, index, part->Pt(), part->Eta(), part->Phi()));
      UInt_t bits = (UInt_t)part->TestBits(TObject::kBitMask);
      track->SetBit(bits);
    }
  }
}
