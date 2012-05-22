#include "AliFMDMCTrackDensity.h"
#include "AliESDFMD.h"
#include "AliTrackReference.h"
#include <TMath.h>
#include "AliFMDStripIndex.h"
#include "AliLog.h"
#include <TH2D.h>
#include <TH1D.h>
#include <TList.h>
#include <TROOT.h>
#include <iostream>

//____________________________________________________________________
void
AliFMDMCTrackDensity::State::Clear(Bool_t alsoCount)
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
  if (alsoCount) count = 0;
}

//____________________________________________________________________
AliFMDMCTrackDensity::State&
AliFMDMCTrackDensity::State::operator=(const State& o)
{
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
  return *this;
}

//____________________________________________________________________
AliFMDMCTrackDensity::AliFMDMCTrackDensity()
  : AliBaseMCTrackDensity(), 
    fState(),
    fMaxConsequtiveStrips(3), 
    fNr(0), 
    fNt(0), 
    fNc(0),
    fNcr(0),
    fOutput(0)
{
  // Default constructor 
}

//____________________________________________________________________
AliFMDMCTrackDensity::AliFMDMCTrackDensity(const char*)
  : AliBaseMCTrackDensity("fmdMCTrackDensity"), 
    fState(),
    fMaxConsequtiveStrips(3), 
    fNr(0), 
    fNt(0), 
    fNc(0),
    fNcr(0),
    fOutput(0)
{
  // Normal constructor constructor 
}

//____________________________________________________________________
AliFMDMCTrackDensity::AliFMDMCTrackDensity(const AliFMDMCTrackDensity& o)
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
AliFMDMCTrackDensity&
AliFMDMCTrackDensity::operator=(const AliFMDMCTrackDensity& o)
{
  // Assignment operator 
  if (&o == this) return *this; 
  AliBaseMCTrackDensity::operator=(o);
  fMaxConsequtiveStrips = o.fMaxConsequtiveStrips;
  fNr                   = o.fNr;
  fNt                   = o.fNt;
  fNc                   = o.fNc;
  fNcr                  = o.fNcr;
  fState                = o.fState;
  fOutput               = o.fOutput;

  return *this;
}

//____________________________________________________________________
void
AliFMDMCTrackDensity::DefineOutput(TList* l)
{
  AliBaseMCTrackDensity::DefineOutput(l);
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
  
		  
}
//____________________________________________________________________
Int_t
AliFMDMCTrackDensity::GetDetectorId() const
{
  return AliTrackReference::kFMD;
}

//____________________________________________________________________
void
AliFMDMCTrackDensity::BeginTrackRefs()
{
  fState.Clear(true);
}

//____________________________________________________________________
void
AliFMDMCTrackDensity::EndTrackRefs(Int_t nRefs)
{
  fNc->Fill(fState.count);
  fNcr->Fill(nRefs, fState.count);
  fState.Clear(true);
}
  
//____________________________________________________________________
AliTrackReference*
AliFMDMCTrackDensity::ProcessRef(AliMCParticle*       particle,
				 const AliMCParticle* mother,
				 AliTrackReference*   ref)
{
  // Get the detector coordinates 
  UShort_t d, s, t;
  Char_t r;
  AliFMDStripIndex::Unpack(ref->UserId(), d, r, s, t);
    
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
AliFMDMCTrackDensity::StoreParticle(AliMCParticle*       particle, 
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

  // Check if we have value already 
  Double_t old = fOutput->Multiplicity(d,r,s,t);

  // If invalid, force it valid 
  if (old == AliESDFMD::kInvalidMult) old = 0;

  // Increment count 
  fOutput->SetMultiplicity(d,r,s,t,old+w);

  // Fill histograms 
  fNr->Fill(fState.nRefs);
  fNt->Fill(fState.nStrips);

  fState.count++;

  return w;
}  

//____________________________________________________________________
Bool_t
AliFMDMCTrackDensity::Calculate(const AliESDFMD&  input, 
				const AliMCEvent& event,
				Double_t          vz,
				AliESDFMD&        output, 
				TH2D*             primary)
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
  output.Clear();
  fOutput = &output;

  // Copy eta values to output 
  for (UShort_t ed = 1; ed <= 3; ed++) { 
    UShort_t nq = (ed == 1 ? 1 : 2);
    for (UShort_t eq = 0; eq < nq; eq++) {
      Char_t   er = (eq == 0 ? 'I' : 'O');
      UShort_t ns = (eq == 0 ?  20 :  40);
      UShort_t nt = (eq == 0 ? 512 : 256);
      for (UShort_t es = 0; es < ns; es++) 
	for (UShort_t et = 0; et < nt; et++) 
	  output.SetEta(ed, er, es, et, input.Eta(ed, er, es, et));
    }
  }

  return ProcessTracks(event, vz, primary);
}
//____________________________________________________________________
void
AliFMDMCTrackDensity::Print(Option_t* option) const 
{
  AliBaseMCTrackDensity::Print(option);
  char ind[gROOT->GetDirLevel()+1];
  for (Int_t i = 0; i < gROOT->GetDirLevel(); i++) ind[i] = ' ';
  ind[gROOT->GetDirLevel()] = '\0';
  std::cout << ind << " Max cluster size:       " << fMaxConsequtiveStrips
	    << std::endl;
  
}

//____________________________________________________________________
//
// EOF
//
