#include "AliFMDMCTrackDensity.h"
#include <AliESDFMD.h>
#include <AliMCEvent.h>
#include <AliTrackReference.h>
#include <AliStack.h>
#include <TMath.h>
#include "AliFMDStripIndex.h"
#include "AliFMDFloatMap.h"
#include <AliLog.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TList.h>
#include <TROOT.h>
#include <iostream>

//____________________________________________________________________
AliFMDMCTrackDensity::AliFMDMCTrackDensity()
  : TNamed(), 
    fUseOnlyPrimary(false), 
    fMaxConsequtiveStrips(3), 
    fNr(0), 
    fNt(0),
    fBinFlow(0), 
    fEtaBinFlow(0),
    fPhiBinFlow(0),
    fDebug(false)
{
  // Default constructor 
}

//____________________________________________________________________
AliFMDMCTrackDensity::AliFMDMCTrackDensity(const char*)
  : TNamed("fmdMCTrackDensity","fmdMCTrackDensity"), 
    fUseOnlyPrimary(false), 
    fMaxConsequtiveStrips(3), 
    fNr(0), 
    fNt(0),
    fBinFlow(0), 
    fEtaBinFlow(0),
    fPhiBinFlow(0),
    fDebug(false)
{
  // Normal constructor constructor 
}

//____________________________________________________________________
AliFMDMCTrackDensity::AliFMDMCTrackDensity(const AliFMDMCTrackDensity& o)
  : TNamed(o),
    fUseOnlyPrimary(o.fUseOnlyPrimary), 
    fMaxConsequtiveStrips(o.fMaxConsequtiveStrips), 
    fNr(o.fNr), 
    fNt(o.fNt),
    fBinFlow(o.fBinFlow), 
    fEtaBinFlow(o.fEtaBinFlow),
    fPhiBinFlow(o.fPhiBinFlow),
    fDebug(o.fDebug)
{
  // Normal constructor constructor 
}

//____________________________________________________________________
AliFMDMCTrackDensity&
AliFMDMCTrackDensity::operator=(const AliFMDMCTrackDensity& o)
{
  // Assignment operator 
  if (&o == this) return *this; 
  TNamed::operator=(o);
  fUseOnlyPrimary       = o.fUseOnlyPrimary;
  fMaxConsequtiveStrips = o.fMaxConsequtiveStrips;
  fNr                   = o.fNr;
  fNt                   = o.fNt;
  fBinFlow              = o.fBinFlow;
  fEtaBinFlow           = o.fEtaBinFlow;
  fPhiBinFlow           = o.fPhiBinFlow;
  fDebug                = o.fDebug;
  return *this;
}

//____________________________________________________________________
void
AliFMDMCTrackDensity::DefineOutput(TList* l)
{
  fNr = new TH1D("clusterRefs", "# track references per cluster",
		 21, -.5, 20.5);
  fNr->SetXTitle("# of track references/cluster");
  fNr->SetDirectory(0);
  l->Add(fNr);

  fNt = new TH1D("clusterSize", "cluster length in strips", 21, -.5, 20.5);
  fNt->SetXTitle("Cluster size [strips]");
  fNt->SetDirectory(0);
  l->Add(fNt);

  fBinFlow = new TH2D("binFlow", "#eta and #varphi bin flow", 
		      200, -5, 5, 40, -180, 180);
  fBinFlow->SetXTitle("#Delta#eta");
  fBinFlow->SetYTitle("#Delta#varphi");
  fBinFlow->SetDirectory(0);
  l->Add(fBinFlow);

  fEtaBinFlow = new TH2D("binFlowEta", "#eta bin flow vs #eta", 
			 200, -4, 6, 200, -5, 5);
  fEtaBinFlow->SetXTitle("#eta of strip");
  fEtaBinFlow->SetYTitle("#Delta#eta");
  fEtaBinFlow->SetDirectory(0);
  l->Add(fEtaBinFlow);

  fPhiBinFlow = new TH2D("binFlowPhi", "#phi bin flow vs #phi", 
			 40, 0, 360, 40, -180, 180);
  fPhiBinFlow->SetXTitle("#varphi of strip");
  fPhiBinFlow->SetYTitle("#Delta#varphi");
  fPhiBinFlow->SetDirectory(0);
  l->Add(fPhiBinFlow);
}

//____________________________________________________________________
void
AliFMDMCTrackDensity::StoreParticle(AliMCParticle*       particle, 
				    const AliMCParticle* mother, 
				    Int_t                longest,
				    Double_t             vz,
				    UShort_t             nC, 
				    UShort_t             nT,
				    AliESDFMD&           output) const
{
  // Store a particle. 
  if (longest < 0) return;

  AliTrackReference* ref = particle->GetTrackReference(longest);
  if (!ref) return;
    
  // Get the detector coordinates 
  UShort_t d, s, t;
  Char_t r;
  AliFMDStripIndex::Unpack(ref->UserId(), d, r, s, t);

  // Check if we have value already 
  Double_t old = output.Multiplicity(d,r,s,t);

  // If invalid, force it valid 
  if (old == AliESDFMD::kInvalidMult) old = 0;

  // Increment count 
  output.SetMultiplicity(d,r,s,t,old+1);

  // Get track-reference stuff 
  Double_t x      = ref->X();
  Double_t y      = ref->Y();
  Double_t z      = ref->Z()-vz;
  Double_t rr     = TMath::Sqrt(x*x+y*y);
  Double_t thetaR = TMath::ATan2(rr,z);
  Double_t phiR   = TMath::ATan2(y,x);
  Double_t etaR   = -TMath::Log(TMath::Tan(thetaR/2));
  
  // Correct angle and convert to degrees 
  if (thetaR < 0) thetaR += 2*TMath::Pi();
  thetaR *= 180. / TMath::Pi();
  if (phiR < 0) phiR += 2*TMath::Pi();
  phiR *= 180. / TMath::Pi();

  // Fill histograms 
  fNr->Fill(nC);
  fNt->Fill(nT);

  const AliMCParticle* mp = (mother ? mother : particle);
  Double_t dEta = mp->Eta() - etaR;
  Double_t dPhi = mp->Phi() * 180 / TMath::Pi() - phiR;
  if (dPhi >  180) dPhi -= 360;
  if (dPhi < -180) dPhi += 360;
  fBinFlow->Fill(dEta, dPhi);
  fEtaBinFlow->Fill(etaR, dEta);
  fPhiBinFlow->Fill(phiR, dPhi);

  if (fDebug) 
    Info("StoreParticle", "Store @ FMD%d%c[%2d,%3d] nStrips=%3d, "
	 "dEta=%7.4f, dPhi=%d",
	 d, r, s, t, nT, dEta, int(dPhi+.5));
}

//____________________________________________________________________
Double_t
AliFMDMCTrackDensity::GetTrackRefTheta(const AliTrackReference* ref,
				       Double_t vz) const
{
  // Get the incidient angle of the track reference. 
  Double_t x    = ref->X();
  Double_t y    = ref->Y();
  Double_t z    = ref->Z()-vz;
  Double_t rr   = TMath::Sqrt(x*x+y*y);
  Double_t theta= TMath::ATan2(rr,z);
  Double_t ang  = TMath::Abs(TMath::Pi()-theta);
  return ang;
}
				    
//____________________________________________________________________
const AliMCParticle*
AliFMDMCTrackDensity::GetMother(Int_t     iTr,
				const AliMCEvent& event) const
{
  // 
  // Track down primary mother 
  // 
  Int_t i  = iTr;
  do { 
    const AliMCParticle* p = static_cast<AliMCParticle*>(event.GetTrack(i));
    if (const_cast<AliMCEvent&>(event).Stack()->IsPhysicalPrimary(i)) return p;
    
    i = p->GetMother();
  } while (i > 0);

  return 0;
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
  AliStack* stack = const_cast<AliMCEvent&>(event).Stack();
  Int_t nTracks   = stack->GetNtrack();//event.GetNumberOfTracks();
  Int_t nPrim     = stack->GetNprimary();//event.GetNumberOfPrimary();
  for (Int_t iTr = 0; iTr < nTracks; iTr++) { 
    AliMCParticle* particle = 
      static_cast<AliMCParticle*>(event.GetTrack(iTr));
    
    // Check the returned particle 
    if (!particle) continue;
    
    // Check if this charged and a primary 
    Bool_t isCharged = particle->Charge() != 0;
    if (!isCharged) continue;
    Bool_t isPrimary = stack->IsPhysicalPrimary(iTr);

    // Pseudo rapidity and azimuthal angle 
    Double_t eta = particle->Eta();
    Double_t phi = particle->Phi();
    
    // Fill 'dn/deta' histogram 
    if (isPrimary && iTr < nPrim) {
      if (primary) primary->Fill(eta, phi);
    }

    // Bail out if we're only processing primaries - perhaps we should
    // track back to the original primary?
    if (fUseOnlyPrimary && !isPrimary) continue;

    Int_t    nTrRef   = particle->GetNumberOfTrackReferences();
    Int_t    longest  = -1;
    Double_t angle    = 0;
    UShort_t oD       = 0, oS = 1024, oT = 1024, ooT=1024;
    Char_t   oR       = '\0';
    UShort_t nC       = 0;
    UShort_t nT       = 0;
    // Double_t oTheta= 0;
    for (Int_t iTrRef = 0; iTrRef < nTrRef; iTrRef++) { 
      AliTrackReference* ref = particle->GetTrackReference(iTrRef);
      
      // Check existence 
      if (!ref) continue;

      // Check that we hit an FMD element 
      if (ref->DetectorId() != AliTrackReference::kFMD) 
	continue;

      // Get the detector coordinates 
      UShort_t d, s, t;
      Char_t r;
      AliFMDStripIndex::Unpack(ref->UserId(), d, r, s, t);

      // Calculate length of last and second to last strips. 

      // IF base of cluster not set, set it here. 
      if (ooT == 1024) ooT = t;

      // Calculate distance of previous reference to base of cluster 
      nT = TMath::Abs(t - ooT) + 1;

      // Count number of track refs in this sector 
      nC++;

      Bool_t used = false;
      // If this is a new detector/ring, then reset the other one 
      // Check if we have a valid old detectorm ring, and sector 
      if (oD >  0 && oR != '\0' && oS != 1024) {
	// New detector, new ring, or new sector 
	if (d != oD   || r != oR   || s != oS) {
	  if (fDebug) Info("Process", "New because new sector");
	  used = true;
	}
      }
      if (nT > fMaxConsequtiveStrips) {
	if (fDebug) Info("Process", "New because too long: %d (%d,%d,%d)", 
			 nT, t, oT, ooT);
	used = true;
      }
      if (used) {
	if (fDebug) 
	  Info("Process", "I=%3d L=%3d D=%d (was %d), R=%c (was %c), "
	       "S=%2d (was %2d) t=%3d (was %3d) nT=%3d/%4d",
	       iTr, longest, d, oD, r, oR, s, oS, t, oT, nT, 
	       fMaxConsequtiveStrips);
	Int_t nnT   = TMath::Abs(oT - ooT) + 1;
	const AliMCParticle* mother = GetMother(iTr, event);
	StoreParticle(particle, mother, longest, vz, nC, nnT, output);
	longest = -1;
	angle   = 0;
	nC  = 1;    // Reset track-ref counter - we have this->1
	nT  = 0;    // Reset cluster size 
	oD  = 0;    // Reset detector
	oR  = '\0'; // Reset ring 
	oS  = 1024; // Reset sector 
	oT  = 1024; // Reset old strip
	ooT = t;    // Reset base 
      }
      else if (fDebug) {
	if (ooT == t) 
	  Info("Process", "New cluster starting at FMD%d%c[%2d,%3d]", 
	       d, r, s, t);
	else 
	  Info("Process", "Adding to cluster starting at FMD%d%c[%2d,%3d], "
	       "length=%3d (now in %3d, previous %3d)", 
	       d, r, s, ooT, nT, t, oT);
      }
	
      oD = d;
      oR = r;
      oS = s;
      oT = t;
      nT = TMath::Abs(t-ooT)+1;

      // The longest passage is determined through the angle 
      Double_t ang  = GetTrackRefTheta(ref, vz);
      if (ang > angle) {
	longest = iTrRef;
	angle   = ang;
      }
      // oTheta = ang;
    } // Loop over track references
    if (longest < 0) continue; // Nothing found

    // Get the reference corresponding to the longest path through the detector
    if (fDebug) 
      Info("Process", "I=%3d L=%3d nT=%3d (out of %3d)", 
	   iTr, longest, nT, fMaxConsequtiveStrips);
    const AliMCParticle* mother = GetMother(iTr, event);
    StoreParticle(particle, mother, longest, vz, nC, nT, output);
  } // Loop over tracks
  return kTRUE;
}
//____________________________________________________________________
void
AliFMDMCTrackDensity::Print(Option_t* /*option*/) const 
{
  char ind[gROOT->GetDirLevel()+1];
  for (Int_t i = 0; i < gROOT->GetDirLevel(); i++) ind[i] = ' ';
  ind[gROOT->GetDirLevel()] = '\0';
  std::cout << ind << ClassName() << ": " << GetName() << '\n'
	    << std::boolalpha 
	    << ind << " Only primary tracks:    " << fUseOnlyPrimary << '\n'
	    << ind << " Max cluster size:       " << fMaxConsequtiveStrips
	    << std::noboolalpha << std::endl;
  
}

//____________________________________________________________________
//
// EOF
//
