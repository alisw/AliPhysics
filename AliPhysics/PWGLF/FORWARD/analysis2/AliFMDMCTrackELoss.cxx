#include "AliFMDMCTrackELoss.h"
#include "AliESDFMD.h"
#include "AliTrackReference.h"
#include <TMath.h>
#include "AliFMDStripIndex.h"
#include "AliFMDEncodedEdx.h"
#include "AliLog.h"
#include "AliMCParticle.h"
#include <TH2D.h>
#include <TH1D.h>
#include <TList.h>
#include <TROOT.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include <iostream>

//====================================================================
void
AliFMDMCTrackELoss::State::Clear(Bool_t alsoCount)
{
  angle       = 0;
  oldDetector = 0;
  oldRing     = '\0';
  oldSector   = 1024;
  oldStrip    = 1024;
  startStrip  = 1024;
  nRefs       = 0;
  nStrips     = 0;
  longest     = 0x0;
  de          = 0;
  dx          = 0;      
  if (alsoCount) count = 0;
}

//____________________________________________________________________
AliFMDMCTrackELoss::State&
AliFMDMCTrackELoss::State::operator=(const State& o)
{
  if (&o == this) return *this;
  angle          = o.angle;
  oldDetector    = o.oldDetector;
  oldRing        = o.oldRing;
  oldSector      = o.oldSector;
  oldStrip       = o.oldStrip;
  startStrip     = o.startStrip;
  nRefs          = o.nRefs;
  nStrips        = o.nStrips;
  count          = o.count;
  longest        = o.longest;
  de             = o.de;
  dx             = o.dx;      
  return *this;
}

//====================================================================
AliFMDMCTrackELoss::Hit::Hit() 
  : fGamma(0), fBeta(0), fEta(0), fDe(0), fDx(0), fDetId(0), fPdg(0), 
    fPrimary(true), fDetector(0), fRing('\0'), 
    fSector(0xFFFF), fStrip(0xFFFF)
{}
//____________________________________________________________________
void 
AliFMDMCTrackELoss::Hit::Decode() const
{ 
  if (fDetector > 0) return;
  AliFMDStripIndex::Unpack(fDetId, fDetector, fRing, fSector, fStrip);
}
//____________________________________________________________________
UInt_t
AliFMDMCTrackELoss::Hit::AbsPdg() const 
{ 
  return TMath::Abs(fPdg); 
}

//====================================================================
AliFMDMCTrackELoss::AliFMDMCTrackELoss()
  : AliBaseMCTrackDensity(), 
    fState(),
    fEvent(),
    fMaxConsequtiveStrips(3), 
    fUseTree(false), 
    fHits(0),
    fTree(0),
    fNr(0), 
    fNt(0), 
    fNc(0),
    fNcr(0),
    fBetaGammadEdx(0),
    fBetaGammaEta(0),
    fDEdxEta(0),
    fPrimaries(0),
    fSecondaries(0),
    fAll(0),
    fEta(0)
{
  // Default constructor 
}

//____________________________________________________________________
AliFMDMCTrackELoss::AliFMDMCTrackELoss(const char*)
  : AliBaseMCTrackDensity("fmdMCTrackELoss"), 
    fState(),
    fEvent(),
    fMaxConsequtiveStrips(3), 
    fUseTree(false), 
    fHits(0),
    fTree(0),
    fNr(0), 
    fNt(0), 
    fNc(0),
    fNcr(0),
    fBetaGammadEdx(0),
    fBetaGammaEta(0),
    fDEdxEta(0),
    fPrimaries(0),
    fSecondaries(0),
    fAll(0),
    fEta(0)
{
  // Normal constructor constructor 
  SetTitle("mcTrackELoss");
  fHits = new TClonesArray("AliFMDMCTrackELoss::Hit");
}

#if 0
//____________________________________________________________________
AliFMDMCTrackELoss::AliFMDMCTrackELoss(const AliFMDMCTrackELoss& o)
  : AliBaseMCTrackDensity(o),
    fState(o.fState),
    fMaxConsequtiveStrips(o.fMaxConsequtiveStrips), 
    fNr(o.fNr), 
    fNt(o.fNt), 
    fNc(o.fNc),
    fNcr(o.fNcr),
    fOutput(o.fOutput)
{
  // Normal constructor constructor 
}

//____________________________________________________________________
AliFMDMCTrackELoss&
AliFMDMCTrackELoss::operator=(const AliFMDMCTrackELoss& o)
{
  // Assignment operator 
  if (&o == this) return *this; 
  AliBaseMCTrackDensity::operator=(o);
  fMaxConsequtiveStrips = o.fMaxConsequtiveStrips;
  fUseTree             = o.fUseTree;
  fNr                   = o.fNr;
  fNt                   = o.fNt;
  fNc                   = o.fNc;
  fNcr                  = o.fNcr;
  fState                = o.fState;

  if (o.fTree) 
    fTree               = static_cast<TTree*>(o.fTree->Clone());

  return *this;
}
#endif

//____________________________________________________________________
void
AliFMDMCTrackELoss::CreateOutputObjects(TList* l)
{
  AliBaseMCTrackDensity::CreateOutputObjects(l);
  TList* ll = static_cast<TList*>(l->FindObject(GetTitle()));
  if (!ll) ll = l;

  fNr = new TH1D("clusterRefs", "# track references per cluster",
		 21, -.5, 20.5);
  fNr->SetXTitle("# of track references/cluster");
  fNr->SetFillColor(kRed+1);
  fNr->SetFillStyle(3001);
  fNr->SetDirectory(0);
  ll->Add(fNr);

  fNt = new TH1D("clusterSize", "cluster length in strips", 21, -.5, 20.5);
  fNt->SetXTitle("Cluster size [strips]");
  fNt->SetFillColor(kBlue+1);
  fNt->SetFillStyle(3001);
  fNt->SetDirectory(0);
  ll->Add(fNt);

  fNc = new TH1D("nClusters", "# clusters per track", 21, -.5, 20.5);
  fNc->SetXTitle("# clusters");
  fNc->SetFillColor(kGreen+1);
  fNc->SetFillStyle(3001);
  fNc->SetDirectory(0);
  ll->Add(fNc);

  fNcr = new TH2D("clusterVsRefs", "# clusters vs. # refs", 
		  21, -.5, 20.5, 21, -.5, 20.5);
  fNcr->SetXTitle("# References");
  fNcr->SetYTitle("# Clusters");
  fNcr->SetOption("COLZ");
  fNcr->SetDirectory(0);
  ll->Add(fNcr);
  

  TArrayD edepArray;
  TArrayD bgArray(601);
  AliFMDEncodedEdx::GetdEAxis().FillBinArray(edepArray);
  AliForwardUtil::MakeLogScale(600, -2, 5, bgArray);
  for (Int_t i = 0; i < edepArray.GetSize(); i++) 
    edepArray[i] *= 10;

  fBetaGammadEdx = new TH2D("betaGammadEdx", "Energy loss", 
			    bgArray.GetSize()-1, bgArray.GetArray(),
			    edepArray.GetSize()-1, edepArray.GetArray());
  fBetaGammadEdx->SetDirectory(0);
  fBetaGammadEdx->SetXTitle("#it{#beta#gamma}");
  fBetaGammadEdx->SetYTitle("d#it{#Delta}/d#it{x}");
  fBetaGammadEdx->Sumw2();
  ll->Add(fBetaGammadEdx);
  
  fBetaGammaEta = new TH2D("betaGammaEta", "#beta#gamma", 
			   200, -4, 6, bgArray.GetSize()-1, bgArray.GetArray());
  fBetaGammaEta->SetXTitle("#eta");
  fBetaGammaEta->SetYTitle("#it{#beta#gamma}");
  fBetaGammaEta->Sumw2();
  ll->Add(fBetaGammaEta);

  fDEdxEta = new TH2D("dEdxEta", "d#it{#Delta}/d#it{x}", 
		      200, -4, 6, edepArray.GetSize()-1, edepArray.GetArray());
  fDEdxEta->SetXTitle("#eta");
  fDEdxEta->SetYTitle("d#it{#Delta}/d#it{x}");
  fDEdxEta->Sumw2();
  ll->Add(fDEdxEta);

  if (!fUseTree) return;

  fTree  = new TTree("tree", "Tree of hits");
  fTree->Branch("event", &fEvent, "ipZ/D:cent");
  fTree->Branch("hits", &fHits);
}

//____________________________________________________________________
void
AliFMDMCTrackELoss::Clear(Option_t*)
{
  fPrimaries.Reset(0);
  fSecondaries.Reset(0);
  fAll.Reset(0);
  fEta.Reset(AliESDFMD::kInvalidEta);
  fHits->Clear();
}

//____________________________________________________________________
Int_t
AliFMDMCTrackELoss::GetDetectorId() const
{
  return AliTrackReference::kFMD;
}

//____________________________________________________________________
void
AliFMDMCTrackELoss::BeginTrackRefs()
{
  fState.Clear(true);
}

//____________________________________________________________________
void
AliFMDMCTrackELoss::EndTrackRefs(Int_t nRefs)
{
  fNc->Fill(fState.count);
  fNcr->Fill(nRefs, fState.count);
  fState.Clear(true);
}
  
//____________________________________________________________________
AliTrackReference*
AliFMDMCTrackELoss::ProcessRef(AliMCParticle*       particle,
			       const AliMCParticle* mother,
			       AliTrackReference*   ref)
{
  // Get the detector coordinates 
  UShort_t d, s, t;
  Char_t r;
  AliFMDStripIndex::Unpack(ref->UserId(), d, r, s, t);
  Double_t edep, length, dEdep, dLength;
  AliFMDEncodedEdx::Decode((ref->UserId() >> 19), edep, length, dEdep, dLength);

		 
		 
		 
  // Calculate distance of previous reference to base of cluster 
  UShort_t nT = TMath::Abs(t - fState.startStrip) + 1;

  // Now check if we should flush to output 
  Bool_t used = false;

  // If this is a new detector/ring, then reset the other one 
  // Check if we have a valid old detectorm ring, and sector 
  if (fState.oldDetector >  0 && 
      fState.oldRing     != '\0' && 
      fState.oldSector   != 1024) {
    // New detector, new ring, or new sector 
    if (d != fState.oldDetector   || 
	r != fState.oldRing       || 
	s != fState.oldSector) {
      if (fDebug) Info("Process", "New because new sector");
      used = true;
    }
    else if (nT > fMaxConsequtiveStrips) {
      if (fDebug) Info("Process", "New because too long: %d (%d,%d,%d)", 
		       fState.nStrips, t, fState.oldStrip, fState.startStrip);
      used = true;
    }
  }
  if (used) {
    if (fDebug) 
      Info("Process", "I=%p L=%p D=%d (was %d), R=%c (was %c), "
	   "S=%2d (was %2d) t=%3d (was %3d) nT=%3d/%4d",
	   ref, fState.longest, 
	   d, fState.oldDetector, 
	   r, fState.oldRing, 
	   s, fState.oldSector, 
	   t, fState.oldStrip, 
	   fState.nStrips, fMaxConsequtiveStrips);
    // Int_t nnT   = TMath::Abs(fState.oldStrip - fState.startStrip) + 1;
    StoreParticle(particle, mother, fState.longest);
    fState.Clear(false);
  }

  // If base of cluster not set, set it here. 
  if (fState.startStrip == 1024) fState.startStrip = t;
  
  // Calculate distance of previous reference to base of cluster 
  fState.nStrips = TMath::Abs(t - fState.startStrip) + 1;

  // Count number of track refs in this sector 
  fState.nRefs++;

  fState.oldDetector = d;
  fState.oldRing     = r;
  fState.oldSector   = s;
  fState.oldStrip    = t;
  fState.de          += edep;
  fState.dx          += length;

  // Debug output 
  if (fDebug) {
    if (t == fState.startStrip) 
      Info("Process", "New cluster starting at FMD%d%c[%2d,%3d]", 
	   d, r, s, t);
    else 
      Info("Process", "Adding to cluster starting at FMD%d%c[%2d,%3d], "
	   "length=%3d (now in %3d, previous %3d)", 
	   d, r, s, fState.startStrip, fState.nStrips, t, fState.oldStrip);
  }
    
  // The longest passage is determined through the angle 
  Double_t ang  = GetTrackRefTheta(ref);
  if (ang > fState.angle) {
    fState.longest = ref;
    fState.angle   = ang;
  }

  return fState.longest;
}

//____________________________________________________________________
Double_t
AliFMDMCTrackELoss::StoreParticle(AliMCParticle*       particle, 
				  const AliMCParticle* mother, 
				  AliTrackReference*   ref) const
{
  Double_t w = 
    AliBaseMCTrackDensity::StoreParticle(particle, mother, ref);
  if (w <= 0) return w;

  // Get the detector coordinates 
  UShort_t d, s, t;
  Char_t r;
  AliFMDStripIndex::Unpack(ref->UserId(), d, r, s, t);
  
  Double_t de   = fState.de;
  Double_t dx   = fState.dx;
  Bool_t   prim = particle == mother;
  fAll(d,r,s,t) += de;
  if (prim) fPrimaries(d,r,s,t)   += de;
  else      fSecondaries(d,r,s,t) += de;

  // Fill histograms 
  fNr->Fill(fState.nRefs);
  fNt->Fill(fState.nStrips);

  fState.count++;

  if (de <= 0 || dx <= 0) return w;
  
  TLorentzVector v;  
  particle->Particle()->Momentum(v);
  if (v.E() <= 0) return w;
  if (v.Beta() < 0 || v.Beta() > 1) return w;

  Double_t eta   = fEta(d,r,s,t);
  Double_t beta  = v.Beta();
  Double_t gamma = v.Gamma();
  Double_t dEdx  = de/dx;


  fBetaGammadEdx->Fill(beta*gamma, dEdx);
  fBetaGammaEta->Fill(eta, beta*gamma);
  fDEdxEta->Fill(eta, dEdx);

  Int_t nHits = fHits->GetEntries();
  Hit* hit = new((*fHits)[nHits]) Hit;
  hit->fPrimary   = prim;
  hit->fEta       = eta;
  hit->fDe        = de;
  hit->fDx        = dx;
  hit->fDetId     = ref->UserId() & AliFMDStripIndex::kIdMask;
  hit->fBeta      = beta;
  hit->fGamma     = gamma;
  hit->fPdg       = particle->PdgCode();

  return w;
}  

//____________________________________________________________________
Bool_t
AliFMDMCTrackELoss::Calculate(const AliESDFMD&  input, 
			      const AliMCEvent& event,
			      const TVector3&   ip,
			      Double_t          cent)
{
  // 
  // Filter the input kinematics and track references, using 
  // some of the ESD information
  // 
  // Parameters:
  //    input   Input ESD event
  //    event   Input MC event
  //    vz      Vertex position 
  //    output  Output ESD-like object
  //    primary Per-event histogram of primaries 
  //
  // Return:
  //    True on succes, false otherwise 
  //
  Clear();
  fEvent.fCent = cent;
  fEvent.fIpZ  = ip.Z();

  // Copy eta values to output 
  for (UShort_t ed = 1; ed <= 3; ed++) { 
    UShort_t nq = (ed == 1 ? 1 : 2);
    for (UShort_t eq = 0; eq < nq; eq++) {
      Char_t   er = (eq == 0 ? 'I' : 'O');
      UShort_t ns = (eq == 0 ?  20 :  40);
      UShort_t nt = (eq == 0 ? 512 : 256);
      for (UShort_t es = 0; es < ns; es++) 
	for (UShort_t et = 0; et < nt; et++) 
	  fEta(ed, er, es, et) = input.Eta(ed, er, es, et);
    }
  }

  Bool_t ret = ProcessTracks(event, ip, 0);

  if (fTree) fTree->Fill();

  return ret;
}

#define PFV(N,VALUE)					\
  do {							\
    AliForwardUtil::PrintName(N);			\
    std::cout << (VALUE) << std::endl; } while(false)
//____________________________________________________________________
void
AliFMDMCTrackELoss::Print(Option_t* option) const 
{
  AliBaseMCTrackDensity::Print(option);
  gROOT->IncreaseDirLevel();
  PFV("Max cluster size", fMaxConsequtiveStrips);
  gROOT->DecreaseDirLevel();
}

//____________________________________________________________________
//
// EOF
//
