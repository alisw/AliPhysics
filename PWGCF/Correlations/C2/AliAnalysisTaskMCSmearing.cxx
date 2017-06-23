#include <iostream>
#include <vector>
#include <numeric>
#include <assert.h>

#include "THn.h"
#include "TMath.h"
#include "TCutG.h"
#include "TParticle.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"

#include "AliAnalysisTaskMCSmearing.h"
#include "AliAnalysisC2Utils.h"

using std::cout;
using std::endl;

AliAnalysisTaskMCSmearing::AliAnalysisTaskMCSmearing()
  : AliAnalysisTaskSE(),
    fdNdeta(0),
    fEventCounter(0),
    fPiCheck(0),
    fdNdetaOrigin(0),
    fxray(0),
    fNsecondaries(0),    
    fNprimaries(0),
    fITS(0),
    fFMD1(0),
    fFMD2(0),
    fFMD3(0),
    fPipe(0),
    fEarlyDecay(0),
    fOutputList(0)
{
}

AliAnalysisTaskMCSmearing::AliAnalysisTaskMCSmearing(const char* name)
  : AliAnalysisTaskSE(name),
    fdNdeta(0),
    fEventCounter(0),
    fPiCheck(0),
    fdNdetaOrigin(0),
    fxray(0),
    fNsecondaries(0),    
    fNprimaries(0),
    fITS(0),
    fFMD1(0),
    fFMD2(0),
    fFMD3(0),
    fPipe(0),
    fEarlyDecay(0),
    fOutputList(0)
{
  DefineOutput(1, TList::Class());
}

void AliAnalysisTaskMCSmearing::UserCreateOutputObjects()
{
  this->fOutputList = new TList();
  this->fOutputList->SetOwner();
  this->SetupCuts();
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
  // dN/deta with origin
  {
    this->fdNdetaOrigin = new TH2F("dNdetaOrigin", "dN/d#eta Origins",
				   180, -4, 5,
				   cOriginType::kNORIGINTYPES, 0, cOriginType::kNORIGINTYPES);
    TAxis *ax = this->fdNdetaOrigin->GetYaxis();
    ax->SetBinLabel(cOriginType::kPRIMARY + 1, "Primary");
    ax->SetBinLabel(cOriginType::kEARLYDECAY + 1, "Early decay");
    ax->SetBinLabel(cOriginType::kPIPE + 1, "Beam pipe");
    ax->SetBinLabel(cOriginType::kITS + 1, "ITS & support structure");
    ax->SetBinLabel(cOriginType::kFMD + 1, "FMD & support structure");
    ax->SetBinLabel(cOriginType::kOTHER + 1, "Others");
    this->fOutputList->Add(this->fdNdetaOrigin);
  }
  this->fxray = new TH2F("xray", "xray;z;r",
			 900, -100, 350, // z
			 400, 0, 100);    // r
  this->fOutputList->Add(this->fxray);

  this->fPiCheck = new TH1F("pi0check", "pi0check", cPrimaryType::kNSPECIES, 0, cPrimaryType::kNSPECIES);
  this->fPiCheck->GetXaxis()->SetBinLabel(cPrimaryType::kPI0 + 1, "#pi_{0}");
  this->fPiCheck->GetXaxis()->SetBinLabel(cPrimaryType::kPICHARGED + 1, "#pi_{ch}");
  this->fPiCheck->GetXaxis()->SetBinLabel(cPrimaryType::kOTHERS + 1, "others");
  this->fOutputList->Add(this->fPiCheck);

  this->fdNdeta = new TH1F("dNdeta", "dNdeta;#eta;counts",
			   80, -3.1, 4.9);
  this->fOutputList->Add(this->fdNdeta);
  this->fEventCounter = new TH1F("EventCounter", "events vs z;z;counts",
				 20, -20, 20);
  this->fOutputList->Add(this->fEventCounter);

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
  // Recurse if non of the above patterns matched
  else {
    return GetIncidentParticleFromFirstMaterialInteraction(ancestor);
  }
}

AliMCParticle* AliAnalysisTaskMCSmearing::GetFirstNonPrimaryMother(AliMCParticle* p) {
  AliMCEvent* event = this->MCEvent();
  // If `p` is not from material, we don't need to look further
  if (event->Stack()->IsPhysicalPrimary(p->GetLabel())) {
    return 0x0;
  }
  // Are there no more mothers left?
  if (p->GetMother() < 0) {
    return 0x0;
  }
  AliMCParticle* ancestor = dynamic_cast< AliMCParticle* >(event->GetTrack(p->GetMother()));
  // Return p if its ancestor is a primary particle; else recurse
  if (event->Stack()->IsPhysicalPrimary(ancestor->GetLabel())) {
    return p;
  }
  // Recurse if non of the above patterns matched
  else {
    return GetIncidentParticleFromFirstMaterialInteraction(ancestor);
  }
}


AliTrackReference* AliAnalysisTaskMCSmearing::IsHitFMD(AliMCParticle* p) {
  for (Int_t iTrRef = 0; iTrRef < p->GetNumberOfTrackReferences(); iTrRef++) { 
    AliTrackReference* ref = p->GetTrackReference(iTrRef);
    // Check hit on FMD
    if (!ref || AliTrackReference::kFMD != ref->DetectorId()) {
      continue;
    }
    else {
      return ref;
    }
  }
  return 0x0;
}

AliTrackReference* AliAnalysisTaskMCSmearing::IsHitITS(AliMCParticle* p) {
  for (Int_t iTrRef = 0; iTrRef < p->GetNumberOfTrackReferences(); iTrRef++) { 
    AliTrackReference* ref = p->GetTrackReference(iTrRef);
    // Check hit on ITS
    if (!ref || AliTrackReference::kITS != ref->DetectorId()) {
      continue;
    }
    // We are interested if it produced a signal, not only a hit in the support structure.
    // This is an envelop around the active area
    if (ref->R() > 3.5 && ref->R() < 4.5 && TMath::Abs(ref->Z()) < 14.1) {
      return ref;
    }
  }
  return 0x0;
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

void AliAnalysisTaskMCSmearing::UserExec(Option_t* /*option*/)
{
  AliMCEvent* event = this->MCEvent();
  AliStack* stack = event->Stack();
  if (!stack) {
    return;
  }
  // Disregard events without reconstructed vertex
  Float_t event_vtx_z = event->GetPrimaryVertex()->GetZ();
  if (!(TMath::Abs(event_vtx_z) > 0)) {
    return;
  }
  this->fEventCounter->Fill(event_vtx_z);

  // Small helper function to get the eta value of a hit
  auto get_ref_eta = [event_vtx_z](AliTrackReference *ref) {
		       Double_t new_ref_z = ref->Z() - event_vtx_z;
		       Double_t ref_r = TMath::Sqrt(ref->X()*ref->X() + ref->Y()*ref->Y());
		       Double_t theta = TMath::ATan2(ref_r, new_ref_z);
		       if (theta < 0){
			 theta += TMath::TwoPi();
		       }
		       Double_t ref_eta = -TMath::Log(TMath::Tan(theta/2.));
		       return ref_eta;
		     };

  Int_t nTracks   = stack->GetNtrack();
  Int_t nPrim     = stack->GetNprimary();

  // Analysis part concerned with origin of secondaries
  for (Int_t iTr = 0; iTr < nTracks; iTr++) { // 
    AliMCParticle* p = static_cast< AliMCParticle* >(this->MCEvent()->GetTrack(iTr));
    if (p->Charge() == 0) {
      continue;
    }
    if (TMath::Abs(event->GetPrimaryVertex()->GetZ()) > 2) {
      continue;
    }
    // Ignore things that do not make a signal in the FMD or ITS
    if (!(this->IsHitFMD(p) || this->IsHitITS(p))) {
      continue;
    }
    this->fdNdeta->Fill(p->Eta());
    // Fill xray plot with particle's creation vertex position
    this->fxray->Fill(p->Zv(), TMath::Sqrt(p->Xv()*p->Xv() + p->Yv()*p->Yv()));
    if (stack->IsPhysicalPrimary(iTr)) {
      this->fdNdetaOrigin->Fill(p->Eta(), cOriginType::kPRIMARY);
      continue;
    }
    // Fill at most one reference per particle from an FMD hit into the origin hist
    if (AliTrackReference *ref = this->IsHitFMD(p)) {
      this->fdNdetaOrigin->Fill(get_ref_eta(ref), this->GetOriginType(p));
    }
    // Fill at most one reference per particle from an ITS hit into the origin hist
    if (AliTrackReference *ref = this->IsHitITS(p)) {
      this->fdNdetaOrigin->Fill(get_ref_eta(ref), this->GetOriginType(p));
    }
  }

  // Hadron chemstry studies
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
	AliAnalysisC2Utils::Wrap02pi(mom->Phi()),
	AliAnalysisC2Utils::Wrap02pi(etaPhi[1]),
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
  PostData(1, this->fOutputList);
}

Int_t AliAnalysisTaskMCSmearing::GetOriginType(AliMCParticle *p) {
  AliMCEvent* event = this->MCEvent();
  AliStack *stack = event->Stack();
  if (stack->IsPhysicalPrimary(p->GetLabel())) {
    return cOriginType::kPRIMARY;
  }
  Double_t r = TMath::Sqrt(p->Yv() * p->Yv() + p->Xv() * p->Xv());
  if (this->fITS->IsInside(p->Zv(), r)) {
    return cOriginType::kITS;
  }
  if (this->fFMD1->IsInside(p->Zv(), r) ||
      this->fFMD2->IsInside(p->Zv(), r) ||
      this->fFMD3->IsInside(p->Zv(), r)) {
    return cOriginType::kFMD;
  }
  if (this->fPipe->IsInside(p->Zv(), r)) {
    return cOriginType::kPIPE;
  }
  if (this->fEarlyDecay->IsInside(p->Zv(), r)) {
    return cOriginType::kEARLYDECAY;
  }
  return cOriginType::kOTHER;
}

void AliAnalysisTaskMCSmearing::SetupCuts() {
  {
    const Int_t npoints = 19;
    const Double_t xs[npoints] = {
      -73.46064836419426, 74.04822648521639, 75.21315343368164, 99.09415587721901,
      98.80292414010268, 76.37808038214686, 76.52369625070503, 46.52682732772516,
      35.60563718586357, 27.742380283723207, -30.358351270980478, -36.91106535609744,
      -46.08486507526118, -75.06242291833397, -79.13966723796229, -102.43820620726703,
      -101.5645109959181, -74.62557531265949, -73.46064836419426
    };
    const Double_t ys[npoints] = {
      53.309457291078104, 53.236891746925636, 50.697097701589115, 49.898876715911925,
      44.81928862523888, 41.84410131498753, 36.11142332694225, 6.43211576858122,
      6.3595502244287445, 3.3843629141774016, 3.4569284583298696, 6.577246856886163,
      6.649812401038631, 35.09550570880764, 41.989232403292476, 44.96441971354383,
      50.987359878199, 51.78558086387619, 53.309457291078104
    };
    fITS = new TCutG("ITS", npoints, xs, ys);
    this->fOutputList->Add(fITS);
  }
  {
    const Int_t npoints = 7;
    const Double_t xs[npoints] = {
      318.96463952403485, 324.6855425320387, 324.7579590258109, 321.8612992749229,
      318.7473900427183, 318.81980653649043, 318.96463952403485};
    const Double_t ys[npoints] = {
      3.964887267397174, 4.037452811549642, 22.68679765873494, 22.68679765873494,
      20.87265905492314, 8.46395100485043, 3.964887267397174
    };
    fFMD1 = new TCutG("FMD1", npoints, xs, ys);
    this->fOutputList->Add(fFMD1);
  }
  {
    const Int_t npoints = 10;
    const Double_t xs[npoints] = {
      73.46576301098378, 73.61137887954195, 77.39739146205395, 77.54300733061211,
      91.08528310652048, 91.37651484363681, 89.62912442093895, 89.33789268382264,
      81.91148338735675, 73.46576301098378
    };
    const Double_t ys[npoints] = {
      4.037452811549642, 32.84597384008101, 36.7645132243145, 41.69897022668259,
      43.51310883049439, 33.281367104995844, 31.176966324574153, 4.25514944400706,
      3.964887267397174, 4.037452811549642
    };
    fFMD2 = new TCutG("FMD2", npoints, xs, ys);
    this->fOutputList->Add(fFMD2);
  }
  {
    const Int_t npoints = 27;
    const Double_t xs[npoints] = {
      -65.74300733061206, -62.539458222332655, -61.81137887954188, -61.665763010983724,
      -49.434030052098734, -46.812944418051956, -44.91993812679594, -48.85156657786612,
      -66.32547080484468, -65.88862319917021, -73.16941662707794, -77.39227681526444,
      -77.53789268382258, -75.64488639256658, -62.97630582800713, -74.47995944410134,
      -74.47995944410134, -78.1203561580552, -78.99405136940413, -78.99405136940413,
      -75.35365465545027, -74.18872770698502, -74.33434357554319, -65.5973914620539,
      -66.76231841051914, -66.61670254196099, -65.74300733061206
    };
    const Double_t ys[npoints] = {
      4.182583899854585, 3.8923217232446987, 4.182583899854585, 17.389512935604486,
      5.416198150446611, 4.690542708921889, 6.0692880478188584, 7.448033386715828,
      24.791198439156624, 25.44428833652887, 32.84597384008101, 36.90964431261944,
      31.249531868726628, 31.176966324574153, 18.840823818653924, 18.840823818653924,
      28.63717227923764, 29.653089897372247, 24.71863289500415, 18.33286500958662,
      13.9789323604383, 15.212546611030326, 18.477996097891562, 18.477996097891562,
      14.849718890267965, 4.908239341379307, 4.182583899854585
    };
    fFMD3 = new TCutG("FMD3", npoints, xs, ys);
    this->fOutputList->Add(fFMD3);
  }
  {
    const Int_t npoints = 28;
    const Double_t xs[npoints] = {
      -304.0187936742833, -77.6357005305573, -77.6357005305573, -46.320459980416445,
      -44.36325744603266, 403.8361229278572, 394.70251110073286, 393.39770941114364,
      385.5688992736084, 386.22130011840306, 361.4300680162082, 352.9488570338784,
      93.94572165042246, 91.98851911603867, 86.76931235768188, 80.89770475453042,
      -33.92484392931908, -45.01565829082733, -45.01565829082733, -51.539666738773235,
      -50.88726589397868, -79.59290306494108, -80.24530390973564, -95.25052334001145,
      -95.25052334001145, -292.9279793127751, -320.3288147941483, -304.0187936742833
    };
    const Double_t ys[npoints] = {
      2.8823548793052556, 2.8823548793052556, 2.769878285868924, 2.7721278177376507,
      2.891353006780162, 2.891353006780162, 3.116306193652825, 3.5639630355294245,
      3.5324695893672513, 3.0645669606721126, 3.116306193652825, 3.0353230463786662,
      3.030823982641213, 4.072357237861643, 3.9688787719002177, 3.042071641984846,
      3.0353230463786662, 3.22653325522043, 4.052111451043103, 4.000372218062391,
      3.6854377564406624, 3.6854377564406624, 4.387291699483371, 4.346800125846292,
      3.5639630355294245, 3.5639630355294245, 3.0038296002164935, 2.8823548793052556
    };
    fPipe = new TCutG("Pipe", npoints, xs, ys);
    this->fOutputList->Add(fPipe);
  }
  {
    const Int_t npoints = 11;
    const Double_t xs[npoints] = {
      -302.71399198469413, -80.89770475453031, -78.94050222014641, -43.710856601238106,
      -41.101253222059654, 404.4885237726518, 409.7077305310087, -330.1148274660673,
      -448.1993803738898, -358.82046463702966, -302.71399198469413
    };
    const Double_t ys[npoints] = {
      2.8735599450546316, 2.8697865367587028, 2.756584287880847, 2.7641311044727037,
      2.888653578238346, 2.87733335335056, -0.047058075994059756, -0.03951125940220246,
      1.6698426986534267, 2.8320524537994176, 2.8735599450546316
    };
    fEarlyDecay = new TCutG("EarlyDecay", npoints, xs, ys);
    this->fOutputList->Add(fEarlyDecay);
  }
}
