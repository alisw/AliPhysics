#include "AliEmcalJetUtilityEventSubtractor.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliEmcalJetTask.h"

ClassImp(AliEmcalJetUtilityEventSubtractor)

//______________________________________________________________________________
AliEmcalJetUtilityEventSubtractor::AliEmcalJetUtilityEventSubtractor() :
  AliEmcalJetUtility(),
  fJetsSubName(""),
  fParticlesSubName(""),
  fUseExternalBkg(kFALSE),
  fRhoName(""),
  fRhomName(""),
  fRho(1e-6),
  fRhom(1e-6),
  fJetsSub(0x0),
  fParticlesSub(0x0),
  fRhoParam(0),
  fRhomParam(0)
{
  // Dummy constructor.

}

//______________________________________________________________________________
AliEmcalJetUtilityEventSubtractor::AliEmcalJetUtilityEventSubtractor(const char* name) :
  AliEmcalJetUtility(name),
  fJetsSubName(""),
  fParticlesSubName(""),
  fUseExternalBkg(kFALSE),
  fRhoName(""),
  fRhomName(""),
  fRho(1e-6),
  fRhom(1e-6),
  fJetsSub(0x0),
  fParticlesSub(0x0),
  fRhoParam(0),
  fRhomParam(0)
{
  // Default constructor.
}

//______________________________________________________________________________
AliEmcalJetUtilityEventSubtractor::AliEmcalJetUtilityEventSubtractor(const AliEmcalJetUtilityEventSubtractor &other) :
  AliEmcalJetUtility(other),
  fJetsSubName(other.fJetsSubName),
  fParticlesSubName(other.fParticlesSubName),
  fUseExternalBkg(other.fUseExternalBkg),
  fRhoName(other.fRhoName),
  fRhomName(other.fRhomName),
  fRho(other.fRho),
  fRhom(other.fRhom),
  fJetsSub(other.fJetsSub),
  fParticlesSub(other.fParticlesSub),
  fRhoParam(other.fRhoParam),
  fRhomParam(other.fRhomParam)
{
  // Copy constructor.
}

//______________________________________________________________________________
AliEmcalJetUtilityEventSubtractor& AliEmcalJetUtilityEventSubtractor::operator=(const AliEmcalJetUtilityEventSubtractor &other)
{
  // Assignment.

  if (&other == this) return *this;
  AliEmcalJetUtility::operator=(other);
  fJetsSubName = other.fJetsSubName;
  fParticlesSubName = other.fParticlesSubName;
  fUseExternalBkg = other.fUseExternalBkg;
  fRhoName = other.fRhoName;
  fRhomName = other.fRhomName;
  fRho = other.fRho;
  fRhom = other.fRhom;
  fJetsSub = other.fJetsSub;
  fParticlesSub = other.fParticlesSub;
  fRhoParam = other.fRhoParam;
  fRhomParam = other.fRhomParam;
  return *this;
}

//______________________________________________________________________________
void AliEmcalJetUtilityEventSubtractor::Init()
{
  // Initialize the utility.
  // Add constituent subtracted jets to event
  if (!fJetsSubName.IsNull()) {
    if (!(fJetTask->InputEvent()->FindListObject(fJetsSubName)) ) {
      fJetsSub = new TClonesArray("AliEmcalJet");
      fJetsSub->SetName(fJetsSubName);
      fJetTask->InputEvent()->AddObject(fJetsSub);
    } 
    else {
      AliError(Form("%s: Object for subtracted jet branch with name %s already in event! Returning", GetName(), fJetsSubName.Data()));
      return;
    }
  }

  // Add tracks from constituent subtracted jets to event
  if (!fParticlesSubName.IsNull()) {
    if (!(fJetTask->InputEvent()->FindListObject(fParticlesSubName))) {
      fParticlesSub = new TClonesArray("AliEmcalParticle");
      fParticlesSub->SetName(fParticlesSubName);
      fJetTask->InputEvent()->AddObject(fParticlesSub);
    } else {
      AliError(Form("%s: Object with name %s already in event! Returning", GetName(), fParticlesSubName.Data()));
      return;
    }
  }

  if (!fRhoName.IsNull() && !fRhoParam) { // get rho from the event
    fRhoParam = dynamic_cast<AliRhoParameter*>(fJetTask->InputEvent()->FindListObject(fRhoName));
    if (!fRhoParam) {
      AliError(Form("%s: Could not retrieve rho %s!", GetName(), fRhoName.Data()));
      return;
    }
  }

  if (!fRhomName.IsNull() && !fRhomParam) { // get rhom from the event
    fRhomParam = dynamic_cast<AliRhoParameter*>(fJetTask->InputEvent()->FindListObject(fRhomName));
    if (!fRhomParam) {
      AliError(Form("%s: Could not retrieve rho_m %s!", GetName(), fRhomName.Data()));
      return;
    }
  }

  // Create particle container in the main jet task so the indices of the constituents
  // added later are properly mapped by the jet task.
  // It is safe to create it before the actual jet finding because the underlying array
  // will be empty until after jet finding.
  fJetTask->AddParticleContainer(fParticlesSubName);

  fInit = kTRUE;
}

//______________________________________________________________________________
void AliEmcalJetUtilityEventSubtractor::InitEvent(AliFJWrapper& fjw)
{
  // Prepare the utility.
  if (!fInit) return;
  if (fRhoParam) fRho = fRhoParam->GetVal();
  if (fRhomParam) fRhom = fRhomParam->GetVal();
  if(fRho < 1e-6) {
    fRho = 1e-6;
  }
  if(fRhom < 1e-6) {
    fRhom = 1e-6;
  }

  if (fJetsSub) fJetsSub->Delete();
  fjw.SetEventSub(kTRUE);
  fjw.SetMaxDelR(fMaxDelR);
  fjw.SetUseExternalBkg(fUseExternalBkg, fRho, fRhom);
}

//______________________________________________________________________________
void AliEmcalJetUtilityEventSubtractor::Prepare(AliFJWrapper& /*fjw*/)
{
  // Prepare the utility.

}

//______________________________________________________________________________
void AliEmcalJetUtilityEventSubtractor::ProcessJet(AliEmcalJet* /*jet*/, Int_t /*ij*/, AliFJWrapper& /*fjw*/)
{
  // Process each jet.
}

//______________________________________________________________________________
void AliEmcalJetUtilityEventSubtractor::Terminate(AliFJWrapper& fjw)
{
  // Run termination of the utility (after each event).

  if (!fInit) return;

  if (!fJetsSub) {
    AliWarning(Form("No jet branch to write to for subtracted jets. fJetsSubName: %s", fJetsSubName.Data()));
    return;
  }

#ifdef FASTJET_VERSION
  std::vector<fastjet::PseudoJet> jets_event_sub;
  jets_event_sub = fjw.GetEventSubJets();
  AliDebug(1,Form("%d event constituent subtracted jets found", (Int_t)jets_event_sub.size()));
  for (UInt_t ijet = 0, jetCount = 0; ijet < jets_event_sub.size(); ++ijet) {
    //printf("Jet pt = %f, area = %f", jets_event_sub[ijet].perp(), fjw.GetEventSubJetArea(ijet));
    //Only storing 4-vector and jet area of unsubtracted jet
    if(jets_event_sub[ijet].E()>0.) {
      AliEmcalJet *jet_event_sub = new ((*fJetsSub)[ijet])
        AliEmcalJet(jets_event_sub[ijet].perp(), jets_event_sub[ijet].eta(), jets_event_sub[ijet].phi(), jets_event_sub[ijet].m());
      jet_event_sub->SetLabel(ijet);
      jet_event_sub->SetJetAcceptanceType(fJetTask->FindJetAcceptanceType(jet_event_sub->Eta(), jet_event_sub->Phi_0_2pi(), fJetTask->GetRadius()));

      fastjet::PseudoJet area(fjw.GetEventSubJetAreaVector(ijet));
      jet_event_sub->SetArea(area.perp()); //changed from area.perp()
      jet_event_sub->SetAreaEta(area.eta());
      jet_event_sub->SetAreaPhi(area.phi());
      jet_event_sub->SetAreaEmc(area.perp());
      
      // Fill constituent info
      std::vector<fastjet::PseudoJet> constituents_sub = jets_event_sub[ijet].constituents();
      fJetTask->FillJetConstituents(jet_event_sub, constituents_sub, constituents_sub, 1, fParticlesSubName);
      jetCount++;
    }
  }

#endif

}

