#include "AliEmcalJetUtilitySoftDrop.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliEmcalJetTask.h"

ClassImp(AliEmcalJetUtilitySoftDrop)

//______________________________________________________________________________
AliEmcalJetUtilitySoftDrop::AliEmcalJetUtilitySoftDrop() :
  AliEmcalJetUtility(),
  fGroomedJetsName(""),
  fGroomedJetParticlesName(""),
  fUseExternalBkg(kFALSE),
  fRhoName(""),
  fRhomName(""),
  fRho(1e-6),
  fRhom(1e-6),
  fGroomedJets(0x0),
  fGroomedJetParticles(0x0),
  fRhoParam(0),
  fRhomParam(0)
{
  // Dummy constructor.

}

//______________________________________________________________________________
AliEmcalJetUtilitySoftDrop::AliEmcalJetUtilitySoftDrop(const char* name) :
  AliEmcalJetUtility(name),
  fGroomedJetsName(""),
  fGroomedJetParticlesName(""),
  fUseExternalBkg(kFALSE),
  fRhoName(""),
  fRhomName(""),
  fRho(1e-6),
  fRhom(1e-6),
  fGroomedJets(0x0),
  fGroomedJetParticles(0x0),
  fRhoParam(0),
  fRhomParam(0)
{
  // Default constructor.
}

//______________________________________________________________________________
AliEmcalJetUtilitySoftDrop::AliEmcalJetUtilitySoftDrop(const AliEmcalJetUtilitySoftDrop &other) :
  AliEmcalJetUtility(other),
  fGroomedJetsName(other.fGroomedJetsName),
  fGroomedJetParticlesName(other.fGroomedJetParticlesName),
  fUseExternalBkg(other.fUseExternalBkg),
  fRhoName(other.fRhoName),
  fRhomName(other.fRhomName),
  fRho(other.fRho),
  fRhom(other.fRhom),
  fGroomedJets(other.fGroomedJets),
  fGroomedJetParticles(other.fGroomedJetParticles),
  fRhoParam(other.fRhoParam),
  fRhomParam(other.fRhomParam)
{
  // Copy constructor.
}

//______________________________________________________________________________
AliEmcalJetUtilitySoftDrop& AliEmcalJetUtilitySoftDrop::operator=(const AliEmcalJetUtilitySoftDrop &other)
{
  // Assignment.

  if (&other == this) return *this;
  AliEmcalJetUtility::operator=(other);
  fGroomedJetsName = other.fGroomedJetsName;
  fGroomedJetParticlesName = other.fGroomedJetParticlesName;
  fUseExternalBkg = other.fUseExternalBkg;
  fRhoName = other.fRhoName;
  fRhomName = other.fRhomName;
  fRho = other.fRho;
  fRhom = other.fRhom;
  fGroomedJets = other.fGroomedJets;
  fGroomedJetParticles = other.fGroomedJetParticles;
  fRhoParam = other.fRhoParam;
  fRhomParam = other.fRhomParam;
  return *this;
}

//______________________________________________________________________________
void AliEmcalJetUtilitySoftDrop::Init()
{
  // Initialize the utility.

  // Add constituent subtracted jets to event
  if (!fGroomedJetsName.IsNull()) {
    if (!(fJetTask->InputEvent()->FindListObject(fGroomedJetsName)) ) {
      fGroomedJets = new TClonesArray("AliEmcalJet");
      fGroomedJets->SetName(fGroomedJetsName);
      fJetTask->InputEvent()->AddObject(fGroomedJets);
    } 
    else {
      AliError(Form("%s: Object for subtracted jet branch with name %s already in event! Returning", GetName(), fGroomedJetsName.Data()));
      return;
    }
  }

  // Add tracks from constituent subtracted jets to event
  if (!fGroomedJetParticlesName.IsNull()) {
    if (!(fJetTask->InputEvent()->FindListObject(fGroomedJetParticlesName))) {
      fGroomedJetParticles = new TClonesArray("AliEmcalParticle");
      fGroomedJetParticles->SetName(fGroomedJetParticlesName);
      fJetTask->InputEvent()->AddObject(fGroomedJetParticles);
    } else {
      AliError(Form("%s: Object with name %s already in event! Returning", GetName(), fGroomedJetParticlesName.Data()));
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

  fInit = kTRUE;
}


//______________________________________________________________________________
void AliEmcalJetUtilitySoftDrop::InitEvent(AliFJWrapper& /*fjw*/)
{
  // Prepare the utility.

}

//______________________________________________________________________________
void AliEmcalJetUtilitySoftDrop::Prepare(AliFJWrapper& fjw)
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
 
  if (fGroomedJets) fGroomedJets->Delete();

  fjw.SetUseExternalBkg(fUseExternalBkg, fRho, fRhom);
  fjw.DoSoftDrop();
  
}

void AliEmcalJetUtilitySoftDrop::ProcessJet(AliEmcalJet* jet, Int_t ij, AliFJWrapper& fjw)
{
  // Process each jet.
  if (!fInit) return;

  #ifdef FASTJET_VERSION

  std::vector<fastjet::PseudoJet> jets_inclusive;
  std::vector<fastjet::PseudoJet> jets_groomed;
  jets_inclusive = fjw.GetInclusiveJets();
  Int_t ninc = (Int_t)jets_inclusive.size();
  jets_groomed = fjw.GetGroomedJets();
  Int_t ngrmd = (Int_t)jets_groomed.size();
  if( (ngrmd > 0) && (ij<ngrmd) ) {

  	//do we really need this area stuff here?
    fastjet::PseudoJet area(fjw.GetJetAreaVector(ij));
    jet->SetArea(area.perp());
    jet->SetAreaEta(area.eta());
    jet->SetAreaPhi(area.phi());
    jet->SetAreaEmc(area.perp());

    jet->GetShapeProperties()->SetSoftDropZg(jets_groomed[ij].structure_of<fastjet::contrib::SoftDrop>().symmetry());
    jet->GetShapeProperties()->SetSoftDropdR(jets_groomed[ij].structure_of<fastjet::contrib::SoftDrop>().delta_R());

    //getting ungroomed pt
    unsigned k = jets_groomed[ij].user_index();
    if ( (k>0) && (k<ninc) ) jet->GetShapeProperties()->SetSoftDropPtfrac( jets_groomed[ij].perp() / jets_inclusive[k].perp() );

    jet->GetShapeProperties()->SetSoftDropDropCount(jets_groomed[ij].structure_of<fastjet::contrib::SoftDrop>().dropped_count());

  }

  #endif
}

//______________________________________________________________________________
void AliEmcalJetUtilitySoftDrop::Terminate(AliFJWrapper& fjw)
{
  // Run termination of the utility (after each event).

  if (!fInit) return;

  if (!fGroomedJets) {
    AliWarning(Form("No jet branch to write to for subtracted jets. fGroomedJetsName: %s", fGroomedJetsName.Data()));
    return;
  }

#ifdef FASTJET_VERSION
  

#endif

}

