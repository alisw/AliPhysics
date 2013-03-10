// $Id: AliJetConstituentTagCopier.cxx  $
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
#include "AliEmcalParticle.h"
#include "AliLog.h"

ClassImp(AliJetConstituentTagCopier)

//________________________________________________________________________
AliJetConstituentTagCopier::AliJetConstituentTagCopier() : 
  AliAnalysisTaskEmcal("AliJetConstituentTagCopier", kFALSE),
  fMCParticlesName(),
  fMCParticles(0),
  fMCParticlesMap(0)
{
  // Default constructor.
}

//________________________________________________________________________
AliJetConstituentTagCopier::AliJetConstituentTagCopier(const char *name) : 
  AliAnalysisTaskEmcal(name, kFALSE),
  fMCParticlesName("MCParticles"),
  fMCParticles(0),
  fMCParticlesMap(0)
{
  // Standard constructor.
}

//________________________________________________________________________
AliJetConstituentTagCopier::~AliJetConstituentTagCopier()
{
  // Destructor
}

//________________________________________________________________________
void AliJetConstituentTagCopier::ExecOnce()
{
  // Execute once.

  if (!fMCParticles) {
    fMCParticles = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fMCParticlesName));
    if (!fMCParticles) {
      AliError(Form("%s: Could not retrieve MC particles %s!", GetName(), fMCParticlesName.Data()));
      return;
    }
    else if (!fMCParticles->GetClass()->GetBaseClass("AliVParticle")) {
      AliError(Form("%s: Collection %s does not contain AliVParticle objects!", GetName(), fMCParticlesName.Data())); 
      fMCParticles = 0;
      return;
    }
  }

  if (!fMCParticlesMap) {
    fMCParticlesMap = dynamic_cast<AliNamedArrayI*>(InputEvent()->FindListObject(fMCParticlesName + "_Map"));
    // this is needed to map the MC labels with the indexes of the MC particle collection
    // if teh map is not given, the MC labels are assumed to be consistent with the indexes (which is not the case if AliEmcalMCTrackSelector is used)
    if (!fMCParticlesMap) {
      AliWarning(Form("%s: Could not retrieve map for MC particles %s! Will assume MC labels consistent with indexes...", GetName(), fMCParticlesName.Data())); 
      fMCParticlesMap = new AliNamedArrayI("tracksMap",9999);
      for (Int_t i = 0; i < 9999; i++) {
	fMCParticlesMap->AddAt(i,i);
      }
    }
  }

  AliAnalysisTaskEmcal::ExecOnce();
}

//________________________________________________________________________
Bool_t AliJetConstituentTagCopier::Run()
{
  if (fTracks) {
    if (fTracks->GetClass()->GetBaseClass("AliVParticle"))
      DoTrackLoop(fTracks);
    else if (fTracks->GetClass()->GetBaseClass("AliEmcalParticle"))
      DoEmcalParticleLoop(fTracks);
    else 
      AliError(Form("%s: Object type not recognized in collection %s. Nothing will be done.", GetName(), fTracks->GetName()));
  }

  if (fCaloClusters) {
    if (fCaloClusters->GetClass()->GetBaseClass("AliVCluster"))
      DoClusterLoop(fCaloClusters);
    else if (fCaloClusters->GetClass()->GetBaseClass("AliEmcalParticle"))
      DoEmcalParticleLoop(fCaloClusters);
    else 
      AliError(Form("%s: Object type not recognized in collection %s. Nothing will be done.", GetName(), fCaloClusters->GetName()));
  }

  return kTRUE;
}

//________________________________________________________________________
void AliJetConstituentTagCopier::DoClusterLoop(TClonesArray *array)
{
  Double_t totalEnergy = 0;
  for (Int_t i = 0; i < array->GetEntries(); i++) {
    AliVCluster *cluster = static_cast<AliVCluster*>(array->At(i));
    if (!cluster) {
      AliError(Form("%s: Could not get cluster %d", GetName(), i));
      continue;
    }
    if (!AcceptCluster(cluster))
      continue;
    Int_t mcLabel = cluster->GetLabel();
    if (mcLabel > 0) {
      TLorentzVector vect;
      cluster->GetMomentum(vect, fVertex);
      AliDebug(2, Form("Cluster %d, pt = %f, eta = %f, phi = %f, label = %d",
		       i, cluster->E(), vect.Eta(), vect.Phi(), mcLabel));
      totalEnergy += cluster->E();
      Int_t index = -1;
      if (mcLabel < fMCParticlesMap->GetSize())
	index = fMCParticlesMap->At(mcLabel);
      if (index < 0)
	continue;
      AliVParticle *part = static_cast<AliVParticle*>(fMCParticles->At(index));
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

  AliDebug(2, Form("Total energy of MC clusters = %f", 
		   totalEnergy));
}

//________________________________________________________________________
void AliJetConstituentTagCopier::DoTrackLoop(TClonesArray *array)
{
  for (Int_t i = 0; i < array->GetEntries(); i++) {
    AliVParticle *track = static_cast<AliVParticle*>(array->At(i));
    if (!track) {
      AliError(Form("%s: Could not get track %d", GetName(), i));
      continue;
    }
    if (!AcceptTrack(track))
      continue;
    Int_t mcLabel = TMath::Abs(track->GetLabel());
    if (mcLabel != 0) {
      Int_t index = -1;
      if (mcLabel < fMCParticlesMap->GetSize())
	index = fMCParticlesMap->At(mcLabel);
      if (index < 0)
	continue;
      AliVParticle *part = static_cast<AliVParticle*>(fMCParticles->At(index));
      if (!part) {
	AliError(Form("%s: Could not get MC particle %d", GetName(), index));
	continue;
      }
      AliDebug(3, Form("Track %d, pt = %f, eta = %f, phi = %f, label = %d is matched with particle %d, pt = %f, eta = %f, phi = %f", 
		       i, track->Pt(), track->Eta(), track->Phi(), mcLabel, index, part->Pt(), part->Eta(), part->Phi()));
      UInt_t bits = (UInt_t)part->TestBits(TObject::kBitMask);
      track->SetBit(bits);
    }
  }
}

//________________________________________________________________________
void AliJetConstituentTagCopier::DoEmcalParticleLoop(TClonesArray *array)
{
  for (Int_t i = 0; i < array->GetEntries(); i++) {
    AliEmcalParticle *emcpart = static_cast<AliEmcalParticle*>(array->At(i));
    if (!emcpart) {
      AliError(Form("%s: Could not get EmcalParticle %d", GetName(), i));
      continue;
    }
    if (!AcceptEmcalPart(emcpart))
      continue;
    AliVCluster *cluster = emcpart->GetCluster();
    AliVParticle *track = emcpart->GetTrack();
    Int_t mcLabel = 0;
    if (cluster)
      mcLabel = cluster->GetLabel();
    else if (track)
      mcLabel = TMath::Abs(track->GetLabel());
    if (mcLabel != 0) {
      Int_t index = -1;
      if (mcLabel < fMCParticlesMap->GetSize())
	index = fMCParticlesMap->At(mcLabel);
      if (index < 0)
	continue;
      AliVParticle *part = static_cast<AliVParticle*>(fMCParticles->At(index));
      if (!part) {
	AliError(Form("%s: Could not get MC particle %d", GetName(), index));
	continue;
      }
      UInt_t bits = (UInt_t)part->TestBits(TObject::kBitMask);
      emcpart->SetBit(bits);
    }
  }
}
