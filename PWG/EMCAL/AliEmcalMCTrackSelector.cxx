//
// Class to select particles in MC events.
//
// Author: S. Aiola

#include "AliEmcalMCTrackSelector.h"

#include <TClonesArray.h>

#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliAODMCParticle.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliNamedArrayI.h"

#include "AliLog.h"

ClassImp(AliEmcalMCTrackSelector)

//________________________________________________________________________
AliEmcalMCTrackSelector::AliEmcalMCTrackSelector() : 
  AliAnalysisTaskSE("AliEmcalMCTrackSelector"),
  fParticlesOutName("MCParticlesSelected"),
  fOnlyPhysPrim(kTRUE),
  fRejectNK(kFALSE),
  fChargedMC(kFALSE),
  fOnlyHIJING(kFALSE),
  fEtaMax(1),
  fParticlesMapName(""),
  fInit(kFALSE),
  fParticlesIn(0),
  fParticlesOut(0),
  fParticlesMap(0),
  fEvent(0),
  fMC(0),
  fIsESD(kFALSE),
  fDisabled(kFALSE)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalMCTrackSelector::AliEmcalMCTrackSelector(const char *name) : 
  AliAnalysisTaskSE(name),
  fParticlesOutName("MCParticlesSelected"),
  fOnlyPhysPrim(kTRUE),
  fRejectNK(kFALSE),
  fChargedMC(kFALSE),
  fOnlyHIJING(kFALSE),
  fEtaMax(1),
  fParticlesMapName(""),
  fInit(kFALSE),
  fParticlesIn(0),
  fParticlesOut(0),
  fParticlesMap(0),
  fEvent(0),
  fMC(0),
  fIsESD(kFALSE),
  fDisabled(kFALSE)
{
  // Constructor.
}

//________________________________________________________________________
AliEmcalMCTrackSelector::~AliEmcalMCTrackSelector()
{
  // Destructor.
}

//________________________________________________________________________
void AliEmcalMCTrackSelector::UserCreateOutputObjects()
{
  // Create my user objects.
}

//________________________________________________________________________
void AliEmcalMCTrackSelector::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  if (fDisabled) return;

  if (!fInit) {
    fEvent = InputEvent();
    if (!fEvent) {
      AliError("Could not retrieve event! Returning");
      return;
    }

    if (fEvent->InheritsFrom("AliESDEvent")) fIsESD = kTRUE;
    else fIsESD = kFALSE;

    TObject *obj = fEvent->FindListObject(fParticlesOutName);
    if (obj) { // the output array is already present in the array!
      AliError(Form("The output array %s is already present in the event! Task will be disabled.", fParticlesOutName.Data()));
      fDisabled = kTRUE;
      return;
    }
    else {  // copy the array from the standard ESD/AOD collections, and filter if requested      

      fParticlesOut = new TClonesArray("AliAODMCParticle");  // the output will always be of AliAODMCParticle, regardless of the input
      fParticlesOut->SetName(fParticlesOutName);
      fEvent->AddObject(fParticlesOut);

      fParticlesMapName = fParticlesOutName;
      fParticlesMapName += "_Map";

      if (fEvent->FindListObject(fParticlesMapName)) {
	AliError(Form("The output array map %s is already present in the event! Task will be disabled.", fParticlesMapName.Data()));
	fDisabled = kTRUE;
	return;
      }
      else {
	fParticlesMap = new AliNamedArrayI(fParticlesMapName, 99999);
	fEvent->AddObject(fParticlesMap);
      }

      if (!fIsESD) {
	fParticlesIn = static_cast<TClonesArray*>(InputEvent()->FindListObject(AliAODMCParticle::StdBranchName()));
	if (!fParticlesIn) {
	  AliError("Could not retrieve AOD MC particles! Task will be disabled.");
	  fDisabled = kTRUE;
	  return;
	}
	TClass *cl = fParticlesIn->GetClass();
	if (!cl->GetBaseClass("AliAODMCParticle")) {
	  AliError(Form("%s: Collection %s does not contain AliAODMCParticle! Task will be disabled.", GetName(), AliAODMCParticle::StdBranchName())); 
	  fDisabled = kTRUE;
	  fParticlesIn = 0;
	  return;
	}
      }
    }

    fMC = MCEvent();
    if (!fMC) {
      AliError("Could not retrieve MC event! Returning");
      fDisabled = kTRUE;
      return;
    }

    fInit = kTRUE;
  }

  if (fIsESD) ConvertMCParticles();
  else CopyMCParticles();
}

//________________________________________________________________________
void AliEmcalMCTrackSelector::ConvertMCParticles() 
{
  // Convert MC particles in MC AOD articles.

  // clear container (normally a null operation as the event should clean it already)
  fParticlesOut->Delete();

  // clear particles map
  fParticlesMap->Clear();

  const Int_t Nparticles = fMC->GetNumberOfTracks();
  const Int_t nprim = fMC->GetNumberOfPrimaries();

  if (fParticlesMap->GetSize() <= Nparticles) fParticlesMap->Set(Nparticles*2);

  // loop over particles
  for (Int_t iPart = 0, nacc = 0; iPart < Nparticles; iPart++) {

    fParticlesMap->AddAt(-1, iPart);

    AliMCParticle* part = static_cast<AliMCParticle*>(fMC->GetTrack(iPart));
    
    if (!part) continue;

    Bool_t isPhysPrim = fMC->IsPhysicalPrimary(iPart);

    Int_t flag = 0;
    if (iPart < nprim) flag |= AliAODMCParticle::kPrimary;
    if (isPhysPrim) flag |= AliAODMCParticle::kPhysicalPrim;
    if (fMC->IsSecondaryFromWeakDecay(iPart)) flag |= AliAODMCParticle::kSecondaryFromWeakDecay;
    if (fMC->IsSecondaryFromMaterial(iPart)) flag |= AliAODMCParticle::kSecondaryFromMaterial;

    AliAODMCParticle *aodPart = new ((*fParticlesOut)[nacc]) AliAODMCParticle(part, iPart, flag);
    aodPart->SetGeneratorIndex(part->GetGeneratorIndex());    
    aodPart->SetStatus(part->Particle()->GetStatusCode());
    aodPart->SetMCProcessCode(part->Particle()->GetUniqueID());

    if (!AcceptParticle(aodPart)) continue;    

    fParticlesMap->AddAt(nacc, iPart);
    nacc++;
  }
}

//________________________________________________________________________
void AliEmcalMCTrackSelector::CopyMCParticles() 
{
  // Convert standard MC AOD particles in a new array, and filter if requested.

  if (!fParticlesIn) return;

  // clear container (normally a null operation as the event should clean it already)
  fParticlesOut->Delete();

  // clear particles map
  fParticlesMap->Clear();

  const Int_t Nparticles = fParticlesIn->GetEntriesFast();
  
  if (fParticlesMap->GetSize() <= Nparticles) fParticlesMap->Set(Nparticles*2);

  AliDebug(2, Form("Total number of particles = %d", Nparticles));

  Int_t nacc = 0;
  
  // loop over particles
  for (Int_t iPart = 0; iPart < Nparticles; iPart++) {
    fParticlesMap->AddAt(-1, iPart);
    
    AliAODMCParticle* part = static_cast<AliAODMCParticle*>(fParticlesIn->At(iPart));
    
    if (!AcceptParticle(part)) continue;
    
    fParticlesMap->AddAt(nacc, iPart);

    AliAODMCParticle *newPart = new ((*fParticlesOut)[nacc]) AliAODMCParticle(*part);
    newPart->SetGeneratorIndex(part->GetGeneratorIndex());

    nacc++;
  }
}

//________________________________________________________________________
Bool_t AliEmcalMCTrackSelector::AcceptParticle(AliAODMCParticle* part) const
{
  // Determine whether the MC particle is accepted.

  if (!part) return kFALSE;
    
  Int_t partPdgCode = TMath::Abs(part->PdgCode());

  if (fOnlyHIJING && (part->GetGeneratorIndex() != 0)) return kFALSE;

  if (fEtaMax > 0. && TMath::Abs(part->Eta()) > fEtaMax) return kFALSE;
    
  if (fRejectNK && (partPdgCode == 130 || partPdgCode == 2112)) return kFALSE;
    
  if (fChargedMC && part->Charge() == 0) return kFALSE;

  if (fOnlyPhysPrim && !part->IsPhysicalPrimary()) return kFALSE;

  return kTRUE;
}
