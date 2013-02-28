// $Id: AliJetConstituentTagCopier.cxx  $
//
// Copy tags from particle level constituent to detector level
//
// Author: S. Aiola

#include "AliJetConstituentTagCopier.h"

#include <TClonesArray.h>
#include <TH1I.h>

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
    fMCParticlesMap = dynamic_cast<TH1I*>(InputEvent()->FindListObject(fMCParticlesName + "_Map"));
    // this is needed to map the MC labels with the indexes of the MC particle collection
      // if teh map is not given, the MC labels are assumed to be consistent with the indexes (which is not the case if AliEmcalMCTrackSelector is used)
    if (!fMCParticlesMap) {
      AliWarning(Form("%s: Could not retrieve map for MC particles %s! Will assume MC labels consistent with indexes...", GetName(), fMCParticlesName.Data())); 
      fMCParticlesMap = new TH1I("tracksMap","tracksMap",9999,0,1);
      for (Int_t i = 0; i < 9999; i++) {
	fMCParticlesMap->SetBinContent(i,i);
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
      Int_t index = fMCParticlesMap->GetBinContent(mcLabel);
      AliVParticle *part = static_cast<AliVParticle*>(fMCParticles->At(index));
      if (!part) {
	AliError(Form("%s: Could not get MC particle %d", GetName(), index));
	continue;
      }
      UInt_t bits = (UInt_t)part->TestBits(TObject::kBitMask);
      cluster->SetBit(bits);
    }
  }
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
    Int_t mcLabel = track->GetLabel();
    if (mcLabel != 0) {
      Int_t index = fMCParticlesMap->GetBinContent(mcLabel);
      AliVParticle *part = static_cast<AliVParticle*>(fMCParticles->At(index));
      if (!part) {
	AliError(Form("%s: Could not get MC particle %d", GetName(), index));
	continue;
      }
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
      mcLabel = track->GetLabel();
    if (mcLabel != 0) {
      Int_t index = fMCParticlesMap->GetBinContent(mcLabel);
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
