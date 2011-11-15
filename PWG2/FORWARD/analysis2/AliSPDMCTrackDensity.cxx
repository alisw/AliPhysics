#include "AliSPDMCTrackDensity.h"
#include <AliMCEvent.h>
#include <AliTrackReference.h>
#include <AliStack.h>
#include <TMath.h>
#include <AliLog.h>
#include <TH3D.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TList.h>
#include <TROOT.h>
#include <iostream>

//____________________________________________________________________
AliSPDMCTrackDensity::AliSPDMCTrackDensity()
  : TNamed(), 
    fUseOnlyPrimary(false), 
    fMinR(3.5), 
    fMaxR(4.5),
    fMinZ(-15), // -14.1), 
    fMaxZ(+15), // +14.1),
    fRZ(0), 
    fXYZ(0),
    fNRefs(0),
    fBinFlow(0)
{
  // Default constructor 
}

//____________________________________________________________________
AliSPDMCTrackDensity::AliSPDMCTrackDensity(const char*)
  : TNamed("spdMCTrackDensity","spdMCTrackDensity"), 
    fUseOnlyPrimary(false), 
    fMinR(3.5), 
    fMaxR(4.5),
    fMinZ(-14.1), 
    fMaxZ(+14.1),
    fRZ(0), 
    fXYZ(0),
    fNRefs(0),
    fBinFlow(0)
{
  // Normal constructor constructor 
}

//____________________________________________________________________
AliSPDMCTrackDensity::AliSPDMCTrackDensity(const AliSPDMCTrackDensity& o)
  : TNamed(o),
    fUseOnlyPrimary(o.fUseOnlyPrimary), 
    fMinR(o.fMinR), 
    fMaxR(o.fMaxR),
    fMinZ(o.fMinZ), 
    fMaxZ(o.fMaxZ),
    fRZ(o.fRZ), 
    fXYZ(o.fXYZ),
    fNRefs(o.fNRefs),
    fBinFlow(o.fBinFlow)
{
  // Normal constructor constructor 
}

//____________________________________________________________________
AliSPDMCTrackDensity&
AliSPDMCTrackDensity::operator=(const AliSPDMCTrackDensity& o)
{
  // Assignment operator 
  if (&o == this) return *this;
  TNamed::operator=(o);
  fUseOnlyPrimary       = o.fUseOnlyPrimary;
  fMinR             = o.fMinR;
  fMaxR             = o.fMaxR;
  fMinZ             = o.fMinZ;
  fMaxZ             = o.fMaxZ;
  fRZ               = o.fRZ;
  fXYZ              = o.fXYZ;
  fNRefs            = o.fNRefs;
  fBinFlow          = o.fBinFlow;
  return *this;
}

//____________________________________________________________________
void
AliSPDMCTrackDensity::DefineOutput(TList* l)
{
  fRZ = new TH2D("rz", "(r,z) of used track references", 
		 30, -15, 15, 50, 0, 5);
  fRZ->SetXTitle("z [cm]");
  fRZ->SetYTitle("r [cm]");
  fRZ->SetDirectory(0);
  l->Add(fRZ);

  fXYZ = new TH3D("xyz", "(x,y,z) of used track references", 
		  100, -5., 5., 100, -5., 5., 30, -15., 15.);
  fXYZ->SetXTitle("x [cm]");
  fXYZ->SetYTitle("y [cm]");
  fXYZ->SetZTitle("z [cm]");
  fXYZ->SetDirectory(0);
  l->Add(fXYZ);

  fNRefs = new TH1D("nrefs", "Number of references used per track", 
		    11, -.5, 10.5); 
  fNRefs->SetXTitle("# of references per track");
  fNRefs->SetDirectory(0);
  l->Add(fNRefs);

  fBinFlow = new TH2D("binFlow", "#eta and #varphi bin flow", 
		      120, -3, 3, 40, -180, 180);
  fBinFlow->SetXTitle("#Delta#eta");
  fBinFlow->SetYTitle("#Delta#varphi");
  fBinFlow->SetDirectory(0);
  l->Add(fBinFlow);
}

//____________________________________________________________________
void
AliSPDMCTrackDensity::StoreParticle(AliMCParticle* particle, 
				    const AliMCParticle* mother, 
				    Int_t          refNo,
				    Double_t       vz, 
				    TH2D&     output) const
{
  // Store a particle. 
  if (refNo < 0) return;

  AliTrackReference* ref = particle->GetTrackReference(refNo);
  if (!ref) return;
    
  Double_t r = ref->R();
  Double_t x = ref->X();
  Double_t y = ref->Y();
  Double_t z = ref->Z();

  Double_t zr = z-vz;
  Double_t th = TMath::ATan2(r,zr);
  if (th < 0) th += 2*TMath::Pi();
  Double_t et = -TMath::Log(TMath::Tan(th/2));
  Double_t ph = TMath::ATan2(y,x);
  if (ph < 0) ph += 2*TMath::Pi();
  output.Fill(et,ph);

  const AliMCParticle* mp = (mother ? mother : particle);
  Double_t dEta = mp->Eta() - et;
  Double_t dPhi = (mp->Phi() - ph) * 180 / TMath::Pi();
  if (dPhi >  180) dPhi -= 360;
  if (dPhi < -180) dPhi += 360;
  fBinFlow->Fill(dEta, dPhi);
}


//____________________________________________________________________
const AliMCParticle*
AliSPDMCTrackDensity::GetMother(Int_t     iTr,
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
AliSPDMCTrackDensity::Calculate(const AliMCEvent& event, 
				Double_t          vz,
				TH2D&             output, 
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

  AliStack* stack = const_cast<AliMCEvent&>(event).Stack();
  Int_t nTracks   = stack->GetNtrack();
  Int_t nPrim     = stack->GetNtrack();
  for (Int_t iTr = 0; iTr < nTracks; iTr++) { 
    AliMCParticle* particle = 
      static_cast<AliMCParticle*>(event.GetTrack(iTr));
    
    // Check the returned particle 
    if (!particle) continue;
    
    // Check if this charged and a primary 
    Bool_t isCharged = particle->Charge() != 0;
    if (!isCharged) continue;

    Bool_t isPrimary = stack->IsPhysicalPrimary(iTr);

    // Fill 'dn/deta' histogram 
    if (isPrimary && iTr < nPrim) {
      if (primary) {
	Double_t eta = particle->Eta();
	Double_t phi = particle->Phi();
	primary->Fill(eta, phi);
      }
    }

    // Bail out if we're only processing primaries - perhaps we should
    // track back to the original primary?
    if (fUseOnlyPrimary && !isPrimary) continue;

    Int_t nTrRef  = particle->GetNumberOfTrackReferences();
    Int_t nRef    = 0;
    for (Int_t iTrRef = 0; iTrRef < nTrRef; iTrRef++) { 
      AliTrackReference* ref = particle->GetTrackReference(iTrRef);
      
      // Check existence 
      if (!ref) continue;

      // Check that we hit an FMD element 
      if (ref->DetectorId() != AliTrackReference::kITS) 
	continue;
      
      // Get radius and z where the track reference was made 
      Double_t r = ref->R();
      Double_t x = ref->X();
      Double_t y = ref->Y();
      Double_t z = ref->Z();
      if (r > fMaxR || r < fMinR) continue;
      if (z > fMaxZ || z < fMinZ) continue;

      fRZ->Fill(z,r);
      fXYZ->Fill(x, y, z);
      nRef++;
      // Only fill first reference 
      if (nRef == 1) { 
	const AliMCParticle* mother = GetMother(iTr, event);
	StoreParticle(particle, mother, iTrRef, vz, output);
      }
    }
    fNRefs->Fill(nRef);
  }
  return kTRUE;
}
//____________________________________________________________________
void
AliSPDMCTrackDensity::Print(Option_t* /*option*/) const 
{
  char ind[gROOT->GetDirLevel()+1];
  for (Int_t i = 0; i < gROOT->GetDirLevel(); i++) ind[i] = ' ';
  ind[gROOT->GetDirLevel()] = '\0';
  std::cout << ind << ClassName() << ": " << GetName() << '\n'
	    << std::boolalpha 
	    << ind << " Only primary tracks:    " << fUseOnlyPrimary << '\n'
	    << ind << " R range:                [" << fMinR << ',' << fMaxR
	    << "]\n"
	    << ind << " Z range:                [" << fMinZ << ',' << fMaxZ
	    << "]" << std::noboolalpha << std::endl;
  
}

//____________________________________________________________________
//
// EOF
//
