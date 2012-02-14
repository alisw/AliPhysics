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
#include "AliGenHijingEventHeader.h"
#include <TF1.h>
#include <TGraph.h>

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
				    TH2D&     output,
                                    Double_t  w) const
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
  output.Fill(et,ph,w);

  const AliMCParticle* mp = (mother ? mother : particle);
  Double_t dEta = mp->Eta() - et;
  Double_t dPhi = (mp->Phi() - ph) * 180 / TMath::Pi();
  if (dPhi >  180) dPhi -= 360;
  if (dPhi < -180) dPhi += 360;
  fBinFlow->Fill(dEta, dPhi,w);
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

// #define USE_FLOW_WEIGHTS
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

#ifdef USE_FLOW_WEIGHTS
  AliGenHijingEventHeader* hd = dynamic_cast<AliGenHijingEventHeader*>
                (event.GenEventHeader());
  Double_t rp = (hd ? hd->ReactionPlaneAngle() : 0.);
  Double_t b = (hd ? hd->ImpactParameter() : -1 );
#endif

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
#ifdef USE_FLOW_WEIGHTS
        Double_t phi = (mother ? mother->Phi() : particle->Phi());
        Double_t eta = (mother ? mother->Eta() : particle->Eta());
        Double_t pt  = (mother ? mother->Pt() : particle->Pt());
        Int_t    id  = (mother ? mother->PdgCode() : 2212);
	Double_t weight = CalculateWeight(eta, pt, b, phi, rp, id);
#else
        Double_t weight = 1.; 
#endif
	StoreParticle(particle, mother, iTrRef, vz, output, weight);
      }
    }
    fNRefs->Fill(nRef);
  }
  return kTRUE;
}

//____________________________________________________________________
Double_t
AliSPDMCTrackDensity::CalculateWeight(Double_t eta, Double_t pt, Double_t b, 
				      Double_t phi, Double_t rp, Int_t id) const
{
  static TF1 gaus = TF1("gaus", "gaus", -6, 6);
  gaus.SetParameters(0.1, 0., 9);
  //  gaus.SetParameters(0.1, 0., 3);
  //  gaus.SetParameters(0.1, 0., 15);
  
  const Double_t xCumulant2nd4050ALICE[] = {0.00, 0.25, 0.350,
					    0.45, 0.55, 0.650, 
					    0.75, 0.85, 0.950,
					    1.10, 1.30, 1.500,
					    1.70, 1.90, 2.250,
					    2.75, 3.25, 3.750,
					    4.50};
  const Double_t yCumulant2nd4050ALICE[] = {0.00000, 0.043400,
					    0.059911,0.073516,
					    0.089756,0.105486,
					    0.117391,0.128199,
					    0.138013,0.158271,
					    0.177726,0.196383,
					    0.208277,0.216648,
					    0.242954,0.249961,
					    0.240131,0.269006,
					    0.207796};
  const Int_t nPointsCumulant2nd4050ALICE = 
    sizeof(xCumulant2nd4050ALICE)/sizeof(Double_t);                                      
  static TGraph alicePointsPt2(nPointsCumulant2nd4050ALICE,xCumulant2nd4050ALICE,yCumulant2nd4050ALICE);
#if 0
  const Double_t xCumulant4th3040ALICE[] = {0.00,0.250,0.35,
					    0.45,0.550,0.65,
					    0.75,0.850,0.95,
					    1.10,1.300,1.50,
					    1.70,1.900,2.25,
					    2.75,3.250,3.75,
					    4.50,5.500,7.00,
					    9.000000};
  const Double_t yCumulant4th3040ALICE[] = {0.000000,0.037071,
					    0.048566,0.061083,
					    0.070910,0.078831,
					    0.091396,0.102026,
					    0.109691,0.124449,
					    0.139819,0.155561,
					    0.165701,0.173678,
					    0.191149,0.202015,
					    0.204540,0.212560,
					    0.195885,0.000000,
					    0.000000,0.000000};
#endif
  const Double_t xCumulant4th4050ALICE[] = {0.00,0.25,0.350,
					    0.45,0.55,0.650,
					    0.75,0.85,0.950,
					    1.10,1.30,1.500,
					    1.70,1.90,2.250,
					    2.75,3.25,3.750,
					    4.50};
  const Double_t yCumulant4th4050ALICE[] = {0.000000,0.038646,
					    0.049824,0.066662,
					    0.075856,0.081583,
					    0.099778,0.104674,
					    0.118545,0.131874,
					    0.152959,0.155348,
					    0.169751,0.179052,
					    0.178532,0.198851,
					    0.185737,0.239901,
					    0.186098};
  const Int_t nPointsCumulant4th4050ALICE = 
    sizeof(xCumulant4th4050ALICE)/sizeof(Double_t);   
  static TGraph alicePointsPt4(nPointsCumulant4th4050ALICE, 
			       xCumulant4th4050ALICE, 
			       yCumulant4th4050ALICE);

  const Double_t xCumulant4thTPCrefMultTPConlyAll[] = {1.75,
						       4.225,
						       5.965,
						       7.765,
						       9.215,
						       10.46,
						       11.565,
						       12.575};
  const Double_t yCumulant4thTPCrefMultTPConlyAll[] = {0.017855,0.032440,
						       0.055818,0.073137,
						       0.083898,0.086690,
						       0.082040,0.077777};
  const Int_t nPointsCumulant4thTPCrefMultTPConlyAll = 
    sizeof(xCumulant4thTPCrefMultTPConlyAll)/sizeof(Double_t);
  TGraph aliceCent(nPointsCumulant4thTPCrefMultTPConlyAll,
		   xCumulant4thTPCrefMultTPConlyAll,
		   yCumulant4thTPCrefMultTPConlyAll);


  Double_t weight = (20. * gaus.Eval(eta) * (alicePointsPt2.Eval(pt) * 0.5 + 
					     alicePointsPt4.Eval(pt) * 0.5) 
		     * (aliceCent.Eval(b) / aliceCent.Eval(10.46)) 
		     * 2. * TMath::Cos(2. * (phi - rp)));
  if      (TMath::Abs(id) == 211)  weight *= 1.3; //pion flow
  else if (TMath::Abs(id) == 2212) weight *= 1.0;  //proton flow
  else                             weight *= 0.7;
  
  return weight;
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
