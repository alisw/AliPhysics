#include <iostream>
#include <vector>

#include "THn.h"
#include "TMath.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"

#include "AliAnalysisTaskMCSmearing.h"
#include "AliAnalysisC2Utils.h"

using std::cout;
using std::endl;

AliAnalysisTaskMCSmearing::AliAnalysisTaskMCSmearing()
  : AliAnalysisTaskSE(),
    fOutputList(0),
    fNprimaries(0),
    fNsecondaries(0),
    fPiCheck(0)
{
}

AliAnalysisTaskMCSmearing::AliAnalysisTaskMCSmearing(const char* name)
  : AliAnalysisTaskSE(name),
    fOutputList(0),
    fNprimaries(0),
    fNsecondaries(0),
    fPiCheck(0)
{
  DefineOutput(1, TList::Class());
}

void AliAnalysisTaskMCSmearing::UserCreateOutputObjects()
{
  this->fOutputList = new TList();
  this->fOutputList->SetOwner();
  {
    // Dimensions: chaintype, eta_prim, eta_sec, phi_prim, phi_sec, p_prim
    const Int_t ndims = 6;
    Int_t dims[ndims] = {17, 17, 20, 20, 1, cPrimaryType::kNSPECIES};
    Double_t mins[ndims] = {-3.5, -3.5, 0, 0, 0, 0};
    Double_t maxs[ndims] = {5.0, 5.0, 2*TMath::Pi(), 2*TMath::Pi(), 10, cPrimaryType::kNSPECIES};
    this->fNsecondaries = new THnF("Nsecondaries",
				   ("N_{sec};"
				    "#eta_{prim};#eta_{sec};"
				    "#varphi_{prim};#varphi_{sec};"
				    "p_{prim};"
				    "type;"),
				   ndims,
				   dims,
				   mins,
				   maxs
				   );
    TAxis *ax = this->fNsecondaries->GetAxis(ndims-1);
    ax->SetBinLabel(cPrimaryType::kPI0 + 1, "#pi_{0}");
    ax->SetBinLabel(cPrimaryType::kPICHARGED + 1, "#pi_{ch}");
    ax->SetBinLabel(cPrimaryType::kOTHERS + 1, "others");
    this->fOutputList->Add(this->fNsecondaries);
  }
  {
    Int_t ndims = 3;
    Int_t dims[] = {cPrimaryType::kNSPECIES, 50, 20};
    Double_t mins[] = {0, -4., 0};
    Double_t maxs[] = {cPrimaryType::kNSPECIES, 6., 5};
    this->fNprimaries = new THnF("Nprimaries", "N_{prim};species;#eta;p;",
				 ndims,
				 dims,
				 mins,
				 maxs);
    TAxis *ax = this->fNprimaries->GetAxis(0);
    ax->SetBinLabel(cPrimaryType::kPI0 + 1, "#pi_{0}");
    ax->SetBinLabel(cPrimaryType::kPICHARGED + 1, "#pi_{ch}");
    ax->SetBinLabel(cPrimaryType::kOTHERS + 1, "others");
    this->fOutputList->Add(this->fNprimaries);
  }
  this->fPiCheck = new TH1F("pi0check", "pi0check", cPrimaryType::kNSPECIES, 0, cPrimaryType::kNSPECIES);
  this->fPiCheck->GetXaxis()->SetBinLabel(cPrimaryType::kPI0 + 1, "#pi_{0}");
  this->fPiCheck->GetXaxis()->SetBinLabel(cPrimaryType::kPICHARGED + 1, "#pi_{ch}");
  this->fPiCheck->GetXaxis()->SetBinLabel(cPrimaryType::kOTHERS + 1, "others");
  this->fOutputList->Add(this->fPiCheck);
  PostData(1, fOutputList);
}

AliMCParticle* AliAnalysisTaskMCSmearing::GetMother(AliMCParticle* p) {
  // Recurses until the mother IsPhysicalPrimary
  // Return NULL if no mother was found
  AliMCEvent* event = this->MCEvent();
  // GetLabel() is the index on the Stack!
  // event->Stack()->IsPhysicalPrimary(p->GetLabel());
  Bool_t isPP = this->IsRedefinedPhysicalPrimary(p);
  // Return this particle if it is stable
  if (isPP) {
    return p;
  }
  else {
    // No stable particle found and no mother left !?
    if (p->GetMother() < 0) {
      return 0x0;
    }
    AliMCParticle* ancestor = dynamic_cast< AliMCParticle* >(event->GetTrack(p->GetMother()));
    return GetMother(ancestor);
  }
}

AliMCParticle* AliAnalysisTaskMCSmearing::GetChargedMother(AliMCParticle* p) {
  AliMCParticle *mother = this->GetMother(p);
  if (!mother || mother->Charge() == 0) {
    return 0x0;
  }
  return mother;
}

AliMCParticle* AliAnalysisTaskMCSmearing::GetNeutralMother(AliMCParticle* p) {
  AliMCParticle *mother = this->GetMother(p);
  if (!mother || mother->Charge() != 0) {
    return 0x0;
  }
  return mother;
}

Bool_t AliAnalysisTaskMCSmearing::IsRedefinedPhysicalPrimary(AliMCParticle* p) {
  AliMCEvent* event = this->MCEvent();
  // Is this a pi0 which was produced as a primary particle?
  if (TMath::Abs(p->PdgCode()) == 111 /*pi0*/ &&
      p->GetLabel() < event->Stack()->GetNprimary()) {
    return true;
  }
  // Is it a Physical Primary by the standard definition?
  Bool_t isPPStandardDef = event->Stack()->IsPhysicalPrimary(p->GetLabel());
  AliMCParticle *pi0Candidate = dynamic_cast< AliMCParticle* >(event->GetTrack(p->GetMother()));
  // Check if this is a primary originating from a pi0
  if (isPPStandardDef && pi0Candidate) {
    if (TMath::Abs(pi0Candidate->PdgCode()) == 111/*pi0*/) {
      return false; // Don't allow stable particles stemming from pi0!
    }
  }
  return isPPStandardDef;
}

std::vector< AliMCParticle* > AliAnalysisTaskMCSmearing::GetDaughters(AliMCParticle* p) {
  AliMCEvent* event = this->MCEvent();
  std::vector< AliMCParticle* > daughters;
  // Find the decays ("edges") leading downstream from this particle ("vertex")
  AliMCParticle* daughterFirst =
    static_cast< AliMCParticle* >(event->GetTrack(p->GetFirstDaughter()));
  // p's mother does not have daughters (p == mother)
  if (!daughterFirst) {
    return daughters;
  }
  AliMCParticle* daughterLast =
    static_cast< AliMCParticle* >(event->GetTrack(p->GetLastDaughter()));
  // We only have one daughter
  if (!daughterLast) {
    daughterLast = daughterFirst;
  }
  // Perform depth-first-search in decay chain for hits on FMD
  for (Int_t iDaughter = daughterFirst->GetLabel(); iDaughter <= daughterLast->GetLabel(); iDaughter++){
    AliMCParticle* daughter  = static_cast< AliMCParticle* >(event->GetTrack(iDaughter));
    daughters.push_back(daughter);
  }
  return daughters;
}

Int_t AliAnalysisTaskMCSmearing::ParticleProducedNHitsOnFMD(AliMCParticle* p) {
  AliMCEvent* event = this->MCEvent();
  // "Explore" the current particle (Graph theory wise)
  Int_t counter = this->IsHitFMD(p) ? 1 : 0;

  // Find the decays ("edges") leading downstream from this particle ("vertex")
  // Perform depth-first-search in decay chain for hits on FMD
  for (auto daughter : this->GetDaughters(p)){
    counter += ParticleProducedNHitsOnFMD(daughter);
  }
  return counter;
}


AliMCParticle* AliAnalysisTaskMCSmearing::GetIncidentParticleFromFirstMaterialInteraction(AliMCParticle* p) {
  AliMCEvent* event = this->MCEvent();
  // Is this particle from material, but there is no mother?!
  Bool_t pIsFromMat = event->Stack()->IsSecondaryFromMaterial(p->GetLabel());
  // If `p` is not from material, we don't need to look further
  if (!pIsFromMat) {
    return 0x0;
  }
  // `p` is from material, but has no mother; This should not happen
  if (pIsFromMat && (p->GetMother() < 0)) {
    return 0x0;
  }
  AliMCParticle* ancestor = dynamic_cast< AliMCParticle* >(event->GetTrack(p->GetMother()));
  // Return the ancestor if `p` is from material but `ancestor` is not
  if (pIsFromMat &&
      !event->Stack()->IsSecondaryFromMaterial(ancestor->GetLabel())) {
    return ancestor;
  }
  // Recurse if non of the above patterns mached
  else {
    return GetIncidentParticleFromFirstMaterialInteraction(ancestor);
  }
}

Bool_t AliAnalysisTaskMCSmearing::IsHitFMD(AliMCParticle* p) {
  for (Int_t iTrRef = 0; iTrRef < p->GetNumberOfTrackReferences(); iTrRef++) { 
    const AliTrackReference* ref = p->GetTrackReference(iTrRef);
    // Check hit on FMD
    if (!ref || AliTrackReference::kFMD != ref->DetectorId()) {
      continue;
    }
    else {
      return true;
    }
  }
  return false;
}

void AliAnalysisTaskMCSmearing::GetTrackRefEtaPhi(AliMCParticle* p, Double_t* etaPhi) {
  AliTrackReference* ref = 0x0;
  for (Int_t iTrRef = 0; iTrRef < p->GetNumberOfTrackReferences(); iTrRef++) { 
    ref = p->GetTrackReference(iTrRef);
    // Check hit on FMD
    if (ref && AliTrackReference::kFMD == ref->DetectorId()) {
      break;
    }
    else {
      ref = 0x0;
    }
  }
  if (!ref) {
    etaPhi = 0x0;
    return;
  }
  const AliVVertex* vertex = this->MCEvent()->GetPrimaryVertex();
  // Calculate the vector pointing from the vertex to the track reference on the detector
  Double_t x      = ref->X() - vertex->GetX();
  Double_t y      = ref->Y() - vertex->GetY();
  Double_t z      = ref->Z() - vertex->GetZ();
  Double_t rr     = TMath::Sqrt(x * x + y * y);
  Double_t thetaR = TMath::ATan2(rr, z);
  Double_t phiR   = TMath::ATan2(y,x);
  // Correct angles
  if (thetaR < 0) {
    thetaR += 2*TMath::Pi();
  }
  if (phiR < 0) {
    phiR += 2*TMath::Pi();
  }
  etaPhi[0] = -TMath::Log(TMath::Tan(thetaR / 2));
  etaPhi[1] = phiR;
  // cout << x << " " << y << " " << z << endl << endl;
  // cout << etaPhi[0] - p->Eta() << " " << etaPhi[1] - p->Phi() << " "
  //      << this->GetDaughters(p).size() << " "
  //      << p->PdgCode() << " "
  //      << endl;
}

void AliAnalysisTaskMCSmearing::UserExec(Option_t *option)
{
  AliMCEvent* event = this->MCEvent();
  AliStack* stack = event->Stack();
  if (!stack) {
    return;
  }
  
  Int_t nTracks   = stack->GetNtrack();
  Int_t nPrim     = stack->GetNprimary();
  for (Int_t iTr = 0; iTr < nTracks; iTr++) { // 
    AliMCParticle* p = static_cast< AliMCParticle* >(this->MCEvent()->GetTrack(iTr));
    AliMCParticle* mom = this->GetMother(p);

    // Fill the pi0 check histogram if we are among still in the generator-primary section
    if (iTr < nPrim) {
      switch (TMath::Abs(p->PdgCode())) {
      case 111:
	this->fPiCheck->Fill(cPrimaryType::kPI0);
	break;
      case 211:
	this->fPiCheck->Fill(cPrimaryType::kPICHARGED);
	break;
      default:
	this->fPiCheck->Fill(cPrimaryType::kOTHERS);
	break;
      }
    }

    if (!mom) {
      continue;
    }
    Int_t primaryMotherEnum;
    switch (TMath::Abs(mom->PdgCode())) {
    case 111:
      primaryMotherEnum = cPrimaryType::kPI0;
      break;
    case 211:
      primaryMotherEnum = cPrimaryType::kPICHARGED;
      break;
    default:
      primaryMotherEnum = cPrimaryType::kOTHERS;
      break;
    }
    
    // Fill the secondary distribution
    if (// p == mom //disregard self correlations where p == mom
	p->Charge() != 0
	&& this->IsHitFMD(p)
	// Only deal with chains that have at least two hits on the FMD
	// && this->ParticleProducedNHitsOnFMD(mom) >= 2
	) {
      // Direction of this particle
      Double_t *etaPhi = new Double_t[2];
      // etaPhi = {p->Eta(), p->Phi()};
      // its the mother hitting the fmd, but it might be
      // deflected from its original direction. If we are looking at a primary hitting
      // the FMD, we should compare its impact on the FMD with the true direction it had
      if (true) {
	this->GetTrackRefEtaPhi(p, etaPhi);
	if (!etaPhi)
	  cout << "NASTY ERROR!" << endl;
      }

      Double_t stuffing[] = {
	mom->Eta(),
	etaPhi[0],
	AliAnalysisC2Utils::WrapAngle(mom->Phi(),
				      this->fNsecondaries->GetAxis(3)),
	AliAnalysisC2Utils::WrapAngle(etaPhi[1],
				      this->fNsecondaries->GetAxis(4)),
	mom->P(),
	Double_t(primaryMotherEnum)
      };
      this->fNsecondaries->Fill(stuffing);
    }
    // Fill the primary distribution; Only fill with primary particles
    if (this->IsRedefinedPhysicalPrimary(p)) {
      Double_t stuffing[] = {
	Double_t(primaryMotherEnum),
	mom->Eta(),
	mom->P()
      };
      this->fNprimaries->Fill(stuffing);
    }
  }
}
